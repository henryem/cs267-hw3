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
typedef shared [] shared_bucket_t *local_buckets_t;
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

shared [1] local_buckets_t *DIRECTORY;
shared [1] local_bucket_sizes_t *BUCKET_SIZES;

/** Note: The below description is currently false.  The directory returned
  * is itself shared; each thread does not get its own private copy.  This
  * may be very inefficient, adding an extra communication step to every
  * table lookup unless the compiler (or runtime?) is smart about caching.
  * But apparently different threads may have different memory representations
  * of the same shared pointer, so that copying a shared pointer to another
  * thread via memcpy does not work.  There must be some way to translate the
  * pointer correctly, but I cannot find it and it is probably UPC compiler
  * dependent. :-(
  * 
  * Builds, on each thread, a private array of @nThreads pointers-to-shared.
  * The ith pointer in the returned array points to the start of the local
  * bucket-list on thread i.  The (shared, with local affinity) bucket-list for
  * this thread is passed as an argument.
  * 
  * The caller is responsible for freeing both the returned list and the
  * bucket-lists on each thread.
  * 
  * The reason for this bizarre, inefficient directory strategy is that
  * we do not want to fix the number of threads or kmers at compile time, and
  * UPC does not support dynamic block strategies.
  * upc_all_alloc(blockSize, nThreads) might
  * do the right thing, but its value must be cast to some kind of static
  * shared type.  shared [] foo * tells the compiler that everything lives
  * on one thread, while shared [1] foo * tells the compiler that the layout
  * is 1-cyclic.  We cannot say shared [LOCAL_LIST_SIZE] because we do not know
  * the local list size (i.e. nKmers / nThreads) at compile time.
  */
void buildBroadcastDirectory(const packed_kmer_list_t **localBucketLists, const int numLocalKmers, const int sizePerThread, const int threadIdx, const int nThreads) {
  DIRECTORY = upc_all_alloc(nThreads, sizeof(local_buckets_t));
  BUCKET_SIZES = upc_all_alloc(nThreads, sizeof(local_bucket_sizes_t));
  
  local_buckets_t localBuckets = upc_alloc(sizePerThread*sizeof(shared_bucket_t *));
  DIRECTORY[threadIdx] = localBuckets;
  shared_bucket_t *privateBuckets = (shared_bucket_t *) localBuckets;
  local_bucket_sizes_t localBucketSizes = upc_alloc(sizePerThread*sizeof(int));
  BUCKET_SIZES[threadIdx] = localBucketSizes;
  int *privateBucketSizes = (int *) localBucketSizes;
  
  // Flatten each list into (our private pointer to) shared space.
  shared_bucket_t localBucketStorage = upc_alloc(numLocalKmers*sizeof(packed_kmer_t));
  packed_kmer_t *privateBucketStorage = (packed_kmer_t *) localBucketStorage;
  int numKmersEncountered = 0;
  for (int localBucketIdx = 0; localBucketIdx < sizePerThread; localBucketIdx++) {
    packed_kmer_list_t *linkedList = localBucketLists[localBucketIdx];
    int bucketSize = packedListSize(linkedList);
    privateBucketSizes[localBucketIdx] = bucketSize;
    shared_bucket_t localBucket = localBucketStorage + numKmersEncountered;
    privateBuckets[localBucketIdx] = localBucket;
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
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  int64_t ownerThread = coarseHash(hashValue, table->nThreads, table->size);
  int64_t localHash = localHashValue(hashValue, ownerThread, table->nThreads, table->size);
  
  const int bucketSize = BUCKET_SIZES[ownerThread][localHash];
  const packed_kmer_t *bucket = malloc(bucketSize*sizeof(packed_kmer_t));
  upc_memget(bucket, DIRECTORY[ownerThread][localHash], bucketSize*sizeof(packed_kmer_t));
  
  //FIXME
  printf("Looking up matches for kmer ");
  printPacked(kmer);
  printf(".  hashValue = %lld, ownerThread = %lld, localHash = %lld, bucketSize = %d, bucket affinity = %d\n",
    hashValue, ownerThread, localHash, bucketSize, upc_threadof(DIRECTORY[ownerThread][localHash]));
  for (int kmerIdx = 0; kmerIdx < bucketSize; kmerIdx++) {
    if (equalsOnKmer(&bucket[kmerIdx], kmer)) {
      return forwardExtensionPacked(&bucket[kmerIdx]);
    }
    //FIXME
    printf("Potential match ");
    printPacked(&bucket[kmerIdx]);
    printf(" does not match.\n");
  }

  printf("No match found for kmer ");
  printPacked(kmer);
  printf("\n");
  assert(0);
  return 0;
}

//FIXME: Make this freeable.

#endif // SHARED_HASHTABLE_H