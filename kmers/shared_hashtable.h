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
  const packed_kmer_t kmer;
  const ptr_to_shared_list_t next;
};

/** Add @kmer to the front of the list, making a new one.  Note that this
  * typically be followed by a upc_memput() to place the new head in whatever
  * shared memory it belongs to, and then freeing the head returned by this
  * method.  This is a little less efficient than directly
  * mallocing the shared memory, but probably not too bad. */
shareable_kmer_list_t *addFront(const ptr_to_shared_list_t list, const packed_kmer_t *kmer) {
  shareable_kmer_list_t *newFront = (shareable_kmer_list_t *) malloc(sizeof(shareable_kmer_list_t));
  memcpy(&newFront.kmer, kmer, sizeof(packed_kmer_t));
  newFront->next = list;
  return newFront;
}

//FIXME: Make this freeable.


typedef struct shared_hash_table_t shared_hash_table_t {
  int64_t size;
  int64_t sizePerThread;
  shared [] shareable_kmer_list_t *localLists;
};

shared_hash_table_t* makeSharedHashTable(const packed_kmer_t *localKmers, const int numLocalKmers, const int threadIdx, const int nThreads, const int tableSize) {
  assert(tableSize % nThreads == 0);
  
  shared_hash_table *result = (*shared_hash_table *) malloc(sizeof(shared_hash_table_t));
  
  result->size = tableSize;
  result->sizePerThread = tableSize / nThreads;
  result->localLists = (shared [] shareable_kmer_list_t *) upc_all_alloc(result->sizePerThread*sizeof(shareable_kmer_list_t), nThreads);
  
  /* Add all the local kmers. */
  for (int i = 0; i < numLocalKmers; i++) {
    add(&localKmers[i]);
  }
  
  return result;
}

void add(shared_hash_table_t *table, const packed_kmer_t *kmer) {
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  assert(coarseHash(hashValue, THREADS, table->size) == MYTHREAD);
  //TODO: Also assert that the blocking strategy worked correctly.
  const shareable_kmer_list_t *newFront = addFront(table->localLists[hashValue], kmer);
  upc_memput(&table->localLists[hashValue], newFront, sizeof(shareable_kmer_list_t));
}


char lookupForwardExtension(const shared_hash_table_t *table, const packet_kmer_t *kmer) {
  int64_t hashValue = hashPackedKmer(kmer, table->size);
  
  shareable_kmer_list_t remainingBucket = hashtable->localLists[hashValue];
  while (remainingBucket != NULL) {
     if (equalsKmerValue(kmer, &remainingBucket.kmer)) {
        return forwardExtension(&remainingBucket.kmer);
     }
     remainingBucket = *(remainingBucket.next);
  }
  return 0;
}

//FIXME: Make this freeable.

#endif // SHARED_HASHTABLE_H