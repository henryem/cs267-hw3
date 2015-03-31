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

/** NUM_KMERS_PER_SOURCE[nThreads*i+j] will be the
  * number of kmers originally read by thread i that are to be sent
  * to thread j.  It is populated below. */
shared int *NUM_KMERS_PER_SOURCE;

typedef shared [] packed_kmer_t *kmer_send_buffer_shrpt_t;
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
void groupKmersByDestination(const unpacked_kmer_t *localRawKmers, packed_kmer_t **kmersByDestination, int *numKmersByDestination, const int nKmersReadLocally, const int nThreads, const int globalHashTableSize) {
  packed_kmer_t *packedKmers = (packed_kmer_t *) malloc(nKmersReadLocally*sizeof(packed_kmer_t));
  kmer_memory_heap_t *heap = makeMemoryHeap(nKmersReadLocally);
  packed_kmer_list_t **lists = (packed_kmer_list_t **) malloc(nThreads*sizeof(packed_kmer_list_t *));
  for (int rawKmerIdx = 0; rawKmerIdx < nKmersReadLocally; rawKmerIdx++) {
    const unpacked_kmer_t *rawKmer = &localRawKmers[rawKmerIdx];
    packed_kmer_t *packedKmer = &packedKmers[rawKmerIdx];
    toPacked(rawKmer, packedKmer);
    int64_t hashValue = hashPackedKmer(packedKmerValue(packedKmer), globalHashTableSize);
    int64_t destinationThread = coarseHash(hashValue, nThreads, globalHashTableSize);
    setFrontWithHeap(&lists[destinationThread], packedKmer, heap);
  }
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    int size = makeFlatCopy(&kmersByDestination[otherThreadIdx], lists[otherThreadIdx]);
    numKmersByDestination[otherThreadIdx] = size;
  }
}

/** Hash-shuffle @localRawKmers, returning in @dst the (packed) kmers that
  * should be handled by this thread.
  * 
  * Leaves SEND_BUFFERS and NUM_KMERS_PER_SOURCE unfreed.  The caller should
  * free them when all threads are done.  We do things this way so that this
  * thread doesn't have to wait for all other threads to finish before it
  * starts building its local hashtable.  Note that the buffers themselves in
  * SEND_BUFFERS are freed.
  */
int hashShuffleKmers(packed_kmer_t **dst, const unpacked_kmer_t *localRawKmers, const int nKmersReadLocally, const int threadIdx, const int nThreads, const int globalHashTableSize) {
  /* Figure out the destination thread for each local kmer.  This also 
   * compresses the raw kmers into packed form. */
  packed_kmer_t **kmersByDestination = (packed_kmer_t **) malloc(nThreads*sizeof(packed_kmer_t *));
  int *numKmersByDestination = (int *) malloc(nThreads*sizeof(int));
  groupKmersByDestination(localRawKmersBuffer, kmersByDestination, numKmersByDestination, nKmersReadLocally, nThreads, globalHashTableSize);
  
  /* Tell everyone how many kmers to expect from us. */
  NUM_KMERS_PER_SOURCE = (shared int *) upc_all_alloc(nThreads, nThreads*sizeof(int));
  upc_memput(numKmersByDestination, NUM_KMERS_PER_SOURCE + threadIdx*nThreads, nThreads);
  
  /* Place our local data in the appropriate shared buffer. */
  SEND_BUFFERS = (shared kmer_send_buffer_shrpt_t *) upc_all_alloc(nThreads, nThreads*sizeof(kmer_send_buffer_shrpt_t));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    SEND_BUFFERS[threadIdx*nThreads + otherThreadIdx] = (kmer_send_buffer_shrpt_t) upc_alloc(numKmersByDestination[otherThreadIdx]*sizeof(packed_kmer_t));
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
  *dst = malloc(totalKmersReceived*sizeof(packed_kmer_t));
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    // It would probably be much better if we could do asynchronous copying
    // here.  Instead we must wait for each copy to finish before starting the
    // next.
    packed_kmer_t *receiveBufferStart = *dst + numKmersBeforeSource[otherThreadIdx];
    upc_memget(receiveBufferStart, SEND_BUFFERS[otherThreadIdx*nThreads + threadIdx], numKmersBySource[otherThreadIdx]);
    upc_free(SEND_BUFFERS[otherThreadIdx*nThreads + threadIdx]);
  }
  
  for (int otherThreadIdx = 0; otherThreadIdx < nThreads; otherThreadIdx++) {
    free(kmersByDestination[otherThreadIdx]);
  }
  free(kmersByDestination);
  free(numKmersByDestination);
  free(numKmersBySource);
  
  return totalKmersReceived;
}

///////////////////////////////////////////////////////////////////////////////
// end shared data and implementation for hashShuffleKmers
///////////////////////////////////////////////////////////////////////////////

#endif // HASH_SHUFFLE_H