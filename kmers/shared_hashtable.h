#ifndef SHARED_HASHTABLE_H
#define SHARED_HASHTABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "contig_generation.h"
#include "kmer_packing.h"

typedef shared [] packed_kmer_t *shared_bucket_t;
typedef shared [] int *local_bucket_sizes_t;

typedef struct shared_hash_table_t {
  int64_t size;
  int64_t sizePerThread;
  int nThreads;
} shared_hash_table_t;

int64_t localHashValue(int64_t hashValue, int threadIdx, int nThreads, int64_t maxHashValue) {
  int64_t range = maxHashValue / nThreads;
  return hashValue % range; //NOTE: Assumes maxHashValue > nThreads.
}

void addToTemporaryBucket(packed_kmer_list_t **temporaryBuckets, kmer_memory_heap_t *heap, const packed_kmer_t *kmer, const int threadIdx, const int nThreads, const int tableSize) {
  int64_t hashValue = hashPackedKmer(kmer, tableSize);
  int64_t localHash = localHashValue(hashValue, threadIdx, nThreads, tableSize);
  assert(coarseHash(hashValue, nThreads, tableSize) == threadIdx);
  temporaryBuckets[localHash] = addFrontWithHeap(temporaryBuckets[localHash], kmer, heap);
}

/** DIRECTORY[threadIdx+nThreads*i] is a shared pointer to local memory on
  * thread threadIdx for kmers with localHashValue == i.  That block of memory
  * is of size BUCKET_SIZES[threadIdx][i].
  */
shared [1] shared_bucket_t *DIRECTORY;
shared [1] local_bucket_sizes_t *BUCKET_SIZES;

void buildBroadcastDirectory(const packed_kmer_list_t **localBucketLists, const int numLocalKmers, const int sizePerThread, const int threadIdx, const int nThreads) {
  const int numBuckets = sizePerThread*nThreads;
  DIRECTORY = upc_all_alloc(numBuckets, sizeof(shared_bucket_t));
  BUCKET_SIZES = upc_all_alloc(nThreads, sizeof(local_bucket_sizes_t));
  
  local_bucket_sizes_t localBucketSizes = upc_alloc(sizePerThread*sizeof(int));
  BUCKET_SIZES[threadIdx] = localBucketSizes;
  int *privateBucketSizes = (int *) localBucketSizes;
  
  // Flatten each list into (our private pointer to) shared space.  We use
  // a contiguous block of shared memory for all the kmers on this thread,
  // to avoid the overhead of many calls to upc_alloc().  (In fact, calling
  // upc_alloc() once per bucket seems to cause a segfault in the version of
  // BUPC on Hopper, though not on my laptop.)
  shared_bucket_t localBucketStorage = upc_alloc(numLocalKmers*sizeof(packed_kmer_t));
  packed_kmer_t *privateBucketStorage = (packed_kmer_t *) localBucketStorage;
  int numKmersEncountered = 0;
  for (int localBucketIdx = 0; localBucketIdx < sizePerThread; localBucketIdx++) {
    packed_kmer_list_t *linkedList = localBucketLists[localBucketIdx];
    int bucketSize = packedListSize(linkedList);
    privateBucketSizes[localBucketIdx] = bucketSize;
    // The following write should be local.  We are putting a shared pointer
    // (which points to local data) into a shared array.
    DIRECTORY[threadIdx+nThreads*localBucketIdx] = localBucketStorage + numKmersEncountered;
    packed_kmer_t *privateBucket = privateBucketStorage + numKmersEncountered;
    toFlatCopy(privateBucket, linkedList, bucketSize);
    numKmersEncountered += bucketSize;
    
  }

  //FIXME: May not need a upc_barrier here.
  upc_barrier;
}

shared_hash_table_t* buildSharedHashTable(const packed_kmer_t *localKmers, const int numLocalKmers, const int threadIdx, const int nThreads, const int tableSize) {
  assert(tableSize % nThreads == 0);
  
  shared_hash_table_t *result = (shared_hash_table_t *) malloc(sizeof(shared_hash_table_t));
  
  result->size = tableSize;
  result->sizePerThread = tableSize / nThreads;
  result->nThreads = nThreads;
  
  /* Add all the local kmers to temporary buckets before copying them to shared
   * memory.  This is done so that we can avoid implementing a dynamic array
   * in shared memory. */
  packed_kmer_list_t **temporaryBuckets = malloc(result->sizePerThread*sizeof(packed_kmer_list_t *));
  memset(temporaryBuckets, 0, result->sizePerThread*sizeof(packed_kmer_list_t *));
  kmer_memory_heap_t *heap = makeMemoryHeap(numLocalKmers);
  for (int i = 0; i < numLocalKmers; i++) {
#ifdef CHECKSUM
    checkPacked(&localKmers[i]);
#endif
    addToTemporaryBucket(temporaryBuckets, heap, &localKmers[i], threadIdx, nThreads, tableSize);
  }
  
  printf("Added %d kmers to the hash table on thread %d.\n", numLocalKmers, threadIdx);
  
  buildBroadcastDirectory(temporaryBuckets, numLocalKmers, result->sizePerThread, threadIdx, nThreads);
  
  printf("Finished broadcasting directory on thread %d.\n", threadIdx);
  
  free(temporaryBuckets);
  freeMemoryHeap(heap);
  
  return result;
}

char lookupForwardExtension(const shared_hash_table_t *table, const packed_kmer_t *kmer) {
#ifdef CHECKSUM
  checkPacked(kmer);
#endif
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  int64_t ownerThread = coarseHash(hashValue, table->nThreads, table->size);
  int64_t localHash = localHashValue(hashValue, ownerThread, table->nThreads, table->size);
  
  const int bucketSize = BUCKET_SIZES[ownerThread][localHash];
  const packed_kmer_t *bucket = malloc(bucketSize*sizeof(packed_kmer_t));
  upc_memget(bucket, DIRECTORY[ownerThread+table->nThreads*localHash], bucketSize*sizeof(packed_kmer_t));
  
  //FIXME
  printfDebug("Looking up matches for kmer ");
#ifdef DEBUG
  printPacked(kmer);
#endif
  printfDebug(".  hashValue = %lld, ownerThread = %lld, localHash = %lld, bucketSize = %d, bucket affinity = %d\n",
    hashValue, ownerThread, localHash, bucketSize, upc_threadof(DIRECTORY[ownerThread+table->nThreads*localHash]));
  for (int kmerIdx = 0; kmerIdx < bucketSize; kmerIdx++) {
#ifdef CHECKSUM
    checkPacked(&bucket[kmerIdx]);
#endif
    if (equalsOnKmer(&bucket[kmerIdx], kmer)) {
      return forwardExtensionPacked(&bucket[kmerIdx]);
    }
    //FIXME
    printfDebug("Potential match ");
#ifdef DEBUG
    printPacked(&bucket[kmerIdx]);
#endif
    printfDebug(" does not match.\n");
  }

  printf("No match found for kmer ");
  printPacked(kmer);
  printf("\n");
  assert(0);
  return 0;
}

//FIXME: Make this freeable.

#endif // SHARED_HASHTABLE_H