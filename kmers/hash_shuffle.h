#ifndef HASH_SHUFFLE_H
#define HASH_SHUFFLE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include "shared_hashtable.h"

///////////////////////////////////////////////////////////////////////////////
// begin shared data and implementation for hashShuffleKmers
///////////////////////////////////////////////////////////////////////////////

typedef shared [] packed_kmer_t *kmer_send_buffer_shrpt_t;

/** Group the kmers in @localRawKmers by their destination thread.
  * kmersByDestination[threadIdx] becomes a list (of length
  * numKmersByDestination[threadIdx]) of kmers to be handled by thread
  * threadIdx.  Each char is a kmer in packed form.
  * @param localRawKmers consists of @nKmersReadLocally kmers,
  *   each taking up LINE_SIZE chars.
  */
kmer_memory_heap_t *groupKmersByDestination(const unpacked_kmer_t *localRawKmers, packed_kmer_t **kmersByDestination, int *numKmersByDestination, const int nKmersReadLocally, const int nThreads, const int globalHashTableSize) {
  packed_kmer_t *packedKmers = (packed_kmer_t *) malloc(nKmersReadLocally*sizeof(packed_kmer_t));
  kmer_memory_heap_t *heap = makeMemoryHeap(nKmersReadLocally);
  packed_kmer_list_t **kmersByDestinationLists = (packed_kmer_list_t **) malloc(nThreads*sizeof(packed_kmer_list_t *));
  memset(kmersByDestinationLists, 0, nThreads*sizeof(packed_kmer_list_t *));
  for (int rawKmerIdx = 0; rawKmerIdx < nKmersReadLocally; rawKmerIdx++) {
    const unpacked_kmer_t *rawKmer = &localRawKmers[rawKmerIdx];
    packed_kmer_t *packedKmer = &packedKmers[rawKmerIdx];
    toPacked(rawKmer, packedKmer);
    int64_t hashValue = hashPackedKmer(packedKmer, globalHashTableSize);
    int64_t destinationThread = coarseHash(hashValue, nThreads, globalHashTableSize);
    assert(destinationThread >= 0);
    assert(destinationThread < nThreads);
    setFrontWithHeap(&kmersByDestinationLists[destinationThread], packedKmer, heap);
  }
  
  // Now we have the raw kmers, the packed kmers in a single list, and a linked list
  // of kmers for each destination thread.  For efficiency, we pack the linked
  // lists (whose size we now know) into arrays, and free the other stuff.
  // (The caller is responsible for freeing the raw local kmers, of course.)
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    const int size = makeFlatCopy(&kmersByDestination[otherThreadIdx], kmersByDestinationLists[otherThreadIdx]);
    numKmersByDestination[otherThreadIdx] = size;
    printf("Thread %d sending %d kmers to thread %d\n", MYTHREAD, size, otherThreadIdx);
  }
  
  free(kmersByDestinationLists);
  freeMemoryHeap(heap);
  free(packedKmers);
  
  return heap;
}

/** NUM_KMERS_PER_SOURCE[nThreads*j+i] will be the
  * number of kmers originally read by thread i that are to be sent
  * to thread j.  It is populated below.  The blocking strategy is designed
  * so that NUM_KMERS_PER_SOURCE[nThreads*j+i] is resident on thread i,
  * the thread that read the kmers and knows how many it will send. */
shared [1] int *NUM_KMERS_PER_SOURCE;

/** Buffers for sending packed kmers.  The buffer for data sent by thread i
  * to thread j is SEND_BUFFERS[j*nThreads+i].  That buffer is of size
  * NUM_KMERS_PER_SOURCE[j*nThreads+i].  The sending thread is responsible
  * for allocating (to shared memory with affinity to itself) the buffer.
  * 
  * Note that it would be cleaner to use a directory structure here, with
  * private-pointers-to-shared on each machine.  That is, for example, how we
  * handle a similar situation in shared_hashtable.h.  But each entry of this
  * directory is used only once, so broadcasting the directory would not save
  * any time, and would be slightly more complicated.  Instead we use a
  * shared-pointers-to-shared directory.
  */
shared [1] kmer_send_buffer_shrpt_t *SEND_BUFFERS;

/** Hash-shuffle @localRawKmers, returning in @dst the (packed) kmers that
  * should be handled by this thread.
  * 
  * FIXME: Leaves SEND_BUFFERS and NUM_KMERS_PER_SOURCE unfreed.
  */
int hashShuffleKmers(packed_kmer_t **dst, const unpacked_kmer_t *localRawKmers, const int nKmersReadLocally, const int threadIdx, const int nThreads, const int globalHashTableSize) {
  /* Figure out the destination thread for each local kmer.  This also 
   * compresses the raw kmers into packed form. */
  packed_kmer_t **kmersByDestination = malloc(nThreads*sizeof(packed_kmer_t *));
  int *numKmersByDestination = malloc(nThreads*sizeof(int));
  groupKmersByDestination(localRawKmers, kmersByDestination, numKmersByDestination, nKmersReadLocally, nThreads, globalHashTableSize);

  /* Tell everyone how many kmers to expect from us. */
  NUM_KMERS_PER_SOURCE = upc_all_alloc(nThreads*nThreads, sizeof(int));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    NUM_KMERS_PER_SOURCE[otherThreadIdx*nThreads+threadIdx] = numKmersByDestination[otherThreadIdx];
  }
  
  /* Place our local data in the appropriate shared buffer. */
  SEND_BUFFERS = upc_all_alloc(nThreads*nThreads, sizeof(kmer_send_buffer_shrpt_t));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    kmer_send_buffer_shrpt_t sendBuffer = (kmer_send_buffer_shrpt_t) upc_alloc(numKmersByDestination[otherThreadIdx]*sizeof(packed_kmer_t));
    // Treat the buffer as a private pointer for faster memcpy, and to clarify
    // that it indeed has local affinity.
    packed_kmer_t *sendBufferPrivatePtr = (packed_kmer_t *) sendBuffer;
    memcpy(sendBufferPrivatePtr, kmersByDestination[otherThreadIdx], numKmersByDestination[otherThreadIdx]*sizeof(packed_kmer_t));
    free(kmersByDestination[otherThreadIdx]);
    SEND_BUFFERS[otherThreadIdx*nThreads + threadIdx] = sendBuffer;
  }
  free(kmersByDestination);
  free(numKmersByDestination);
  
  printf("Waiting for other threads to finish sending on thread %d.\n", threadIdx);
  
  /* We need just this one barrier to ensure that readers can read sent data.
   * We could use more fine-grained barriers, since each reader only needs to
   * wait for its writer to finish, not everybody.  That would be more
   * complicated. */
  upc_barrier;
  
  printf("Receiving data on thread %d.\n", threadIdx);
  
  /* Compute how much data we are to receive. */
  int totalKmersReceived = 0;
  int *numKmersBySource = (int *) malloc(nThreads*sizeof(int));
  int *numKmersBeforeSource = (int *) malloc(nThreads*sizeof(int));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    // Better to use a prefix sum here, but this will only take O(p) time anyway.
    const int numReceived = NUM_KMERS_PER_SOURCE[threadIdx*nThreads + otherThreadIdx];
    totalKmersReceived += numReceived;
    numKmersBySource[otherThreadIdx] = numReceived;
    numKmersBeforeSource[otherThreadIdx] = totalKmersReceived - numReceived;
  }
  
  /* Receive our new hash-matched data from the appropriate shared buffer. */
  *dst = malloc(totalKmersReceived*sizeof(packed_kmer_t));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    // It would probably be much better if we could do asynchronous copying
    // here.  Instead we must wait for each copy to finish before starting the
    // next.
    //NOTE: The copy here is contiguous, so potentially we could use a single
    // upc_memget here, simplifying the code and solving the issue noted above.
    packed_kmer_t *receiveBufferStart = *dst + numKmersBeforeSource[otherThreadIdx];
    upc_memget(receiveBufferStart, SEND_BUFFERS[threadIdx*nThreads + otherThreadIdx], numKmersBySource[otherThreadIdx]*sizeof(packed_kmer_t));
    upc_free(SEND_BUFFERS[threadIdx*nThreads + otherThreadIdx]);
  }
  
  printf("Finished receiving %d hash-matched kmers on thread %d.\n", totalKmersReceived, threadIdx);
  
  free(numKmersBySource);
  free(numKmersBeforeSource);
  upc_all_free(SEND_BUFFERS);
  upc_all_free(NUM_KMERS_PER_SOURCE);
  
  return totalKmersReceived;
}

///////////////////////////////////////////////////////////////////////////////
// end shared data and implementation for hashShuffleKmers
///////////////////////////////////////////////////////////////////////////////

#endif // HASH_SHUFFLE_H