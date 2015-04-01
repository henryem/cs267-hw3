#ifndef SHARED_HASHTABLE_H
#define SHARED_HASHTABLE_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include "contig_generation.h"
#include "kmer_packing.h"

/** A linked list designed to hold a short list of packed kmers, optimized for
  * searching when the first element of this list is what we're looking for.
  * Designed for use in a distributed hashtable where we want to minimize
  * pointer chasing.  Can be shared. */
typedef struct shareable_kmer_list_t shareable_kmer_list_t;
typedef shared [] shareable_kmer_list_t *ptr_to_shared_list_t;

struct shareable_kmer_list_t {
  packed_kmer_t kmer;
  ptr_to_shared_list_t next;
  // We may encounter uninitialized lists.  Check this flag to see whether
  // a list is initialized.
  bool isInitialized;
};

/** Add @kmer to the front of the list, making a new one.  The new head
  * replaces the current one at *list, and the old head (if it was ever
  * initialized) is given a new place in memory, with affinity to the calling
  * thread.
  */
void addFront(ptr_to_shared_list_t list, const packed_kmer_t *kmer) {
  shareable_kmer_list_t oldHead = *list;
  ptr_to_shared_list_t newNext;
  if (oldHead.isInitialized) {
    newNext = (ptr_to_shared_list_t) upc_alloc(sizeof(shareable_kmer_list_t));
    upc_memput(newNext, &oldHead, sizeof(shareable_kmer_list_t));
  } else {
    newNext = NULL;
  }
  shareable_kmer_list_t newFrontTmp;
  memcpy(&(newFrontTmp.kmer), kmer, sizeof(packed_kmer_t));
  newFrontTmp.next = newNext;
  newFrontTmp.isInitialized = true;
  upc_memput(list, &newFrontTmp, sizeof(shareable_kmer_list_t));
}

//FIXME: Make this freeable.
//FIXME: Could use a memory heap for the extra allocations above.


typedef struct shared_hash_table_t {
  int64_t size;
  int64_t sizePerThread;
  shared [] shareable_kmer_list_t *localLists;
} shared_hash_table_t;

void add(shared_hash_table_t *table, const packed_kmer_t *kmer) {
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  assert(coarseHash(hashValue, THREADS, table->size) == MYTHREAD);
  //TODO: Also assert that the blocking strategy worked correctly.
  addFront(&(table->localLists[hashValue]), kmer);
}

shared_hash_table_t* buildSharedHashTable(const packed_kmer_t *localKmers, const int numLocalKmers, const int threadIdx, const int nThreads, const int tableSize) {
  assert(tableSize % nThreads == 0);
  
  shared_hash_table_t *result = (shared_hash_table_t *) malloc(sizeof(shared_hash_table_t));
  
  result->size = tableSize;
  result->sizePerThread = tableSize / nThreads;
  size_t blockSizeBytes = result->sizePerThread * sizeof(shareable_kmer_list_t);
  result->localLists = (shared [] shareable_kmer_list_t *) upc_all_alloc(blockSizeBytes, nThreads);
  shared [] shareable_kmer_list_t *myBlockStart = &(result->localLists[threadIdx*result->sizePerThread]);
  upc_memset(myBlockStart, 0, blockSizeBytes);
  
  /* Add all the local kmers. */
  for (int i = 0; i < numLocalKmers; i++) {
    add(result, &localKmers[i]);
  }
  
  return result;
}

char lookupForwardExtension(const shared_hash_table_t *table, const packed_kmer_t *kmer) {
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  
  shareable_kmer_list_t remainingBucket = table->localLists[hashValue];
  // The unusual pointer chasing here is because of our weird scheme to
  // optimize communication in the case when the list has 1 element.  The
  // element at table->localLists[hashValue] might be an uninitialized list-
  // plus-kmer; or else it is initialized and the list is a regular linked list
  // terminating at NULL.
  if (remainingBucket.isInitialized) {
    while (true) {
      if (equalsOnKmer(kmer, &remainingBucket.kmer)) {
        return forwardExtensionPacked(&remainingBucket.kmer);
      }
      if (remainingBucket.next == NULL) {
        break;
      } else {
        upc_memget(&remainingBucket, remainingBucket.next, sizeof(shareable_kmer_list_t));
      }
    }
  }
  
  assert(0);
  return 0;
}

//FIXME: Make this freeable.

#endif // SHARED_HASHTABLE_H