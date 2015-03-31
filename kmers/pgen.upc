#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_hash.h"

const size_t PACKED_KMER_SIZE = KMER_PACKED_LENGTH_WITH_EXT * sizeof(char);


///////////////////////////////////////////////////////////////////////////////
// begin packed_kmer_t
///////////////////////////////////////////////////////////////////////////////

typedef struct packed_kmer_t packed_kmer_t {
  const unsigned char[KMER_PACKED_LENGTH_WITH_EXT] kmer;
}

// bool isBackwardStart(packed_kmer_t *kmer) {
//   return 
// }

///////////////////////////////////////////////////////////////////////////////
// end packed_kmer_t
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// begin kmer_list_t
///////////////////////////////////////////////////////////////////////////////

typedef struct kmer_list_t kmer_list_t {
  const unsigned char *kmer;
  const kmer_list_t *next;
}

int listSize(kmer_list_t *list) {
  int size = 0;
  while (list != NULL) {
    size++;
    list = list->next;
  }
}

int makeFlatCopy(unsigned char **dst, const kmer_list_t *list) {
  const int size = listSize(list);
  *dst = (unsigned char *) malloc(size*PACKED_KMER_SIZE);
  int i = 0;
  while (list != NULL) {
    memcpy(&dst[i*PACKED_KMER_SIZE], list->kmer, PACKED_KMER_SIZE);
    list = list->next;;
    i++;
  }
  return size;
}

kmer_list_t *addFront(const kmer_list_t *list, const unsigned char *kmer) {
  kmer_list_t* newFront = (kmer_list_t *) malloc(sizeof(kmer_list_t));
  newFront->kmer = kmer;
  newFront->next = list;
  return newFront
}

void freeKmerList(kmer_list_t* list) {
  while (true) {
    const kmer_list_t *next = list.next;
    free(list);
    if (next == NULL) {
      break;
    } else {
      list = next;
    }
  }
}

typedef struct kmer_memory_heap_t kmer_memory_heap_t {
  kmer_list_t *heap;
  int64_t nextHeapLocation;
}

kmer_list_t *addFrontWithHeap(const kmer_list_t *list, const unsigned char *kmer, kmer_memory_heap_t *heap) {
  kmer_list_t *newFront = &heap->heap[heap->nextHeapLocation];
  newFront->kmer = kmer;
  newFront->next = list;
  heap->nextHeapLocation++;
  return newFront;
}

void setFrontWithHeap(kmer_list_t **list, const unsigned char *kmer, kmer_memory_heap_t *heap) {
  *list = addFrontWithHeap(*list, kmer, heap);
}

kmer_memory_heap_t *makeMemoryHeap(int64_t size) {
  kmer_list_t *kmers = (kmer_list_t *) malloc(sizeof(kmer_list_t)*size);
  kmer_memory_heap_t *heap = (kmer_memory_heap_t *) malloc(sizeof(kmer_memory_heap_t));
  heap->heap = kmers;
  heap->nextHeapLocation = 0;
}

void freeMemoryHeap(kmer_memory_heap_t *heap) {
  free(heap.heap);
}

///////////////////////////////////////////////////////////////////////////////
// end kmer_list_t
///////////////////////////////////////////////////////////////////////////////

/** Read kmers for this thread from @inputFileName. */
unsigned char *readLocalKmers(const char *inputFileName, int nKmersReadTypically, int nKmersToRead, int threadIdx) {
  const int localCharsToRead = nKmersToRead * LINE_SIZE;
  FILE *inputFile = fopen(inputFileName, "r");
  fseek(inputFile, threadIdx*nKmersReadTypically*LINE_SIZE, SEEK_SET);
  unsigned char *localReadBuffer = (unsigned char*) malloc(localCharsToRead * sizeof(unsigned char));
  fread(localReadBuffer, sizeof(unsigned char), localCharsToRead, inputFile);
  fclose(inputFile);
  return localReadBuffer;
}


///////////////////////////////////////////////////////////////////////////////
// begin shared data and implementation for hashShuffleKmers
///////////////////////////////////////////////////////////////////////////////

/** NUM_KMERS_PER_SOURCE[nThreads*i+j] will be the
  * number of kmers originally read by thread i that are to be sent
  * to thread j.  It is populated below. */
shared int *NUM_KMERS_PER_SOURCE;

typedef shared [] unsigned char *kmer_send_buffer_shrpt_t;
/** Buffers for sending packed kmers.  The buffer for data sent by thread i
  * to thread j is SEND_BUFFERS[i*nThreads+j].  That buffer is of size
  * NUM_KMERS_PER_SOURCE[i*nThreads+j].  The sending thread is responsible
  * for allocating (to shared memory with affinity to itself) the buffer.
  */
shared kmer_send_buffer_shrpt_t *SEND_BUFFERS;

/** Group the kmers in @localRawKmers by their destination thread.
  * kmersByDestination[threadIdx] becomes a list (of length
  * numKmersByDestination[threadIdx]) of kmers to be handled by thread
  * threadIdx.  Each char is a kmer in packed form.
  * @param localRawKmers consists of @nKmersReadLocally kmers,
  *   each taking up LINE_SIZE chars.
  */
void groupKmersByDestination(const unsigned char *localRawKmers, unsigned char **kmersByDestination, int *numKmersByDestination, const int nKmersReadLocally, const int nThreads, const int globalHashTableSize) {
  unsigned char *packedKmers = (unsigned char *) malloc(nKmersReadLocally*PACKED_KMER_SIZE);
  kmer_memory_heap_t *heap = makeMemoryHeap(nKmersReadLocally);
  kmer_list_t **lists = (kmer_list_t **) malloc(nThreads*sizeof(kmer_list_t *));
  for (int rawKmerIdx = 0; rawKmerIdx < nKmersReadLocally; rawKmerIdx++) {
    const unsigned char *rawKmer = &localRawKmers[rawKmerIdx*LINE_SIZE];
    unsigned char *packedKmer = &packedKmers[rawKmerIdx*PACKED_KMER_SIZE];
    packSequence(rawKmer, packedKmer, KMER_LENGTH);
    int64_t destinationThread = coarseHash(hashkmer(globalHashTableSize, packedKmer), nThreads, globalHashTableSize);
    setFrontWithHeap(&lists[destinationThread], packedKmer, heap);
  }
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    int size = makeFlatCopy(&kmersByDestination[otherThreadIdx], lists[otherThreadIdx]);
    numKmersByDestination[otherThreadIdx] = size;
  }
}

/** Hash-shuffle @localRawKmers, returning the (packed) kmers that should be
  * handled by this thread.
  * 
  * Leaves SEND_BUFFERS and NUM_KMERS_PER_SOURCE unfreed.  The caller should
  * free them when all threads are done.  We do things this way so that this
  * thread doesn't have to wait for all other threads to finish before it
  * starts building its local hashtable.  Note that the buffers themselves in
  * SEND_BUFFERS are freed.
  */
unsigned char* hashShuffleKmers(const unsigned char *localRawKmers, const int nKmersReadLocally, const int threadIdx, const int nThreads) {
  /* Figure out the destination thread for each local kmer.  This also 
   * compresses the raw kmers into packed form. */
  unsigned char **kmersByDestination = (unsigned char **) malloc(nThreads*sizeof(char *));
  int *numKmersByDestination = (int *) malloc(nThreads*sizeof(int));
  groupKmersByDestination(localRawKmersBuffer, kmersByDestination, numKmersByDestination, nKmersReadLocally, nThreads);
  
  /* Tell everyone how many kmers to expect from us. */
  NUM_KMERS_PER_SOURCE = (shared int *) upc_all_alloc(nThreads, nThreads*sizeof(int));
  upc_memput(numKmersByDestination, NUM_KMERS_PER_SOURCE + threadIdx*nThreads, nThreads);
  
  /* Place our local data in the appropriate shared buffer. */
  SEND_BUFFERS = (shared kmer_send_buffer_shrpt_t *) upc_all_alloc(nThreads, nThreads*sizeof(kmer_send_buffer_shrpt_t));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    SEND_BUFFERS[threadIdx*nThreads + otherThreadIdx] = (kmer_send_buffer_shrpt_t) upc_alloc(PACKED_KMER_SIZE * numKmersByDestination[otherThreadIdx]);
    upc_memput(SEND_BUFFERS[threadIdx*nThreads + otherThreadIdx], kmersByDestination[otherThreadIdx], numKmersByDestination[otherThreadIdx]);
  }
  
  /* We need just this one barrier to ensure that readers can read sent data.
   * We could use more fine-grained barriers, since each reader only needs to
   * wait for its writer to finish, not everybody.  That would be more
   * complicated. */
  upc_barrier;
  
  /* Compute how much data we are to receive. */
  int totalKmersReceived = 0;
  int numKmersBySource = malloc(nThreads*sizeof(int));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    // Better to use a prefix sum here, but this will only take O(p) time anyway.
    const int numReceived = NUM_KMERS_PER_SOURCE[otherThreadIdx*nThreads + threadIdx];
    totalKmersReceived += numReceived;
    numKmersBySource[otherThreadIdx] = numReceived;
    if (otherThreadIdx > 0) {
      numKmersBeforeSource[otherThreadIdx] = numKmersBeforeSource[otherThreadIdx-1] + numReceived;
    } else {
      numKmersBeforeSource[otherThreadIdx] = 0;
    }
  }
  
  /* Receive our new hash-matched data from the appropriate shared buffer. */
  unsigned char *receiveBuffer = malloc(PACKED_KMER_SIZE*totalKmersReceived);
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    // It would probably be much better if we could do asynchronous copying
    // here.  Instead we must wait for each copy to finish before starting the
    // next.
    unsigned char *receiveBufferStart = receiveBuffer + numKmersBeforeSource[otherThreadIdx];
    upc_memget(receiveBufferStart, SEND_BUFFERS[otherThreadIdx*nThreads + threadIdx], numKmersBySource[otherThreadIdx];
    upc_free(SEND_BUFFERS[otherThreadIdx*nThreads + threadIdx]);
  }
  
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    free(kmersByDestination[otherThreadIdx]);
  }
  free(kmersByDestination);
  free(numKmersByDestination);
  free(numKmersBySource);
  
  return receiveBuffer;
}

///////////////////////////////////////////////////////////////////////////////
// end shared data and implementation for hashShuffleKmers
///////////////////////////////////////////////////////////////////////////////
  

int main(int argc, char *argv[]) {
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
	
	init_LookupTable();
  const char *inputFileName = argv[1];


	/** Read input **/
	inputTime -= gettime();
	
  /* Extract the number of k-mers in the input file */
  const int nKmers = getNumKmersInUFX(inputFileName);
  const int nThreads = THREADS;
  const int threadIdx = MYTHREAD;
  const bool isLastThread = (threadIdx == (nThreads - 1));
  const int nKmersReadTypically = DIV_ROUND_UP(nKmers, nThreads);
  
  /* The last thread might get fewer kmers. */
  const int nKmersReadLocally = isLastThread ?
    nKmersReadTypically - (nKmers % nThreads) :
    nKmersReadTypically;

  /* Read the kmers from the input file into a local buffer. */
  unsigned char *localRawKmersBuffer = readLocalKmers(inputFileName, nKmersReadTypically, nKmersReadLocally, threadIdx);
  
	inputTime += gettime();


	/** Graph construction **/
	constrTime -= gettime();
	
	// On each machine:
	//  - Assume kmers have been read into a private array
	//  - For each kmer on this machine, put it in a byte (wasting a bit of space) in a new shared array
	//  - Sort this machine's sublist by hashed machine (preferably using bucket sort for O(n) runtime, but probably quicksort is fine)
	//  - Find the start indices of each machine's kmers in the sorted sublist
	//  - Compute the number of kmers we will receive from each other machine by reducing and broadcasting
	//  - Assign final indices in a shared sorted list to each kmer on our machine
	// Now, once:
	//  - upc_permute_all using these indices
	// On each machine:
	//  - Insert (serially) into the global hashtable.  No locking or communication is required for this step.
	// To build the start list, on each machine:
	//  - Make a shared list of pointers to local start lists (these will never need to be communicated)
	//  - For each kmer on this machine, if it is a start kmer, add it to the local start list
	//  - Now the machine can use this local start list to start building contigs
	//  - Question: Do we need to load-balance these start lists?  That would be annoying but doable.
  const unsigned char *localKmers = hashShuffleKmers(localRawKmersBuffer);
	const   
	
	
  upc_barrier;
	constrTime += gettime();

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	upc_barrier;
	traversalTime += gettime();

	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
