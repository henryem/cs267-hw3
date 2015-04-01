#ifndef KMER_PACKING_H
#define KMER_PACKING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <string.h>
#include <stdbool.h>
#include "contig_generation.h"
#include "packingDNAseq.h"

// Forward declarations:

typedef struct unpacked_kmer_t {
  unsigned char data[LINE_SIZE];
} unpacked_kmer_t;

/** A packed kmer.  The structure is as follows:
  *  - The first KMER_PACKED_LENGTH chars are for the middle of the kmer.
  *    Each char contains 4 bases, represented by 2 bits.  If KMER_LENGTH
  *    is not divisible by 4, the last char may have some extra As, which
  *    should be ignored.
  *  - The last char contains the backward extension (in 2 bits), the forward
  *    extension (in the next 2 bits), a bit to indicate whether the forward
  *    extension is F (the start extension), and last a bit to indicate
  *    whether the backward extension is F.  The last two bits are ignored.
  *    If an extension is F, it will be marked with an A.
  * These packed kmers are designed mostly to be convenient and reasonably
  * efficient for transport.
  */
typedef struct packed_kmer_t {
  unsigned char data[KMER_PACKED_LENGTH_WITH_EXT];
} packed_kmer_t;

void toPacked(const unpacked_kmer_t *rawKmer, packed_kmer_t *packedKmer);
int64_t hashPackedKmer(const packed_kmer_t *packedKmer, int64_t hashTableSize);

int64_t hashBytes(int64_t hashtable_size, char *seq, int size) {
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer (without paying attention to its
 * extensions). */
int64_t hashKmerBytes(int64_t hashtable_size, char *seq) {
   return hashBytes(hashtable_size, seq, KMER_PACKED_LENGTH);
}


/* Further coarsen a hash value.  Useful for partitioning a hashtable across
 * @coarseTableSize processors.
 */
int64_t coarseHash(int64_t hash, int64_t coarseTableSize, int64_t hashTableSize) {
  assert(coarseTableSize < hashTableSize);
  return hash % coarseTableSize;
}

///////////////////////////////////////////////////////////////////////////////
// begin unpacked_kmer_t
///////////////////////////////////////////////////////////////////////////////

const int RAW_KMER_BACKWARD_EXT_POS = KMER_LENGTH+1;
const int RAW_KMER_FORWARD_EXT_POS = KMER_LENGTH+2;

const char START_CHAR = 'F';

const unsigned char *kmerValue(const unpacked_kmer_t *kmer) {
  return kmer->data;
}

unsigned char backwardExtensionUnpacked(const unpacked_kmer_t *kmer) {
  return kmer->data[RAW_KMER_BACKWARD_EXT_POS];
}

unsigned char forwardExtensionUnpacked(const unpacked_kmer_t *kmer) {
  return kmer->data[RAW_KMER_FORWARD_EXT_POS];
}

bool isBackwardStartUnpacked(const unpacked_kmer_t *kmer) {
  return backwardExtensionUnpacked(kmer) == START_CHAR;
}

bool isForwardStartUnpacked(const unpacked_kmer_t *kmer) {
  return forwardExtensionUnpacked(kmer) == START_CHAR;
}

int64_t hashUnpackedKmer(const unpacked_kmer_t *unpackedKmer, int64_t hashTableSize) {
  packed_kmer_t packedKmer;
  toPacked(unpackedKmer, &packedKmer);
  return hashPackedKmer(&packedKmer, hashTableSize);
}

///////////////////////////////////////////////////////////////////////////////
// end unpacked_kmer_t
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// begin packed_kmer_t
///////////////////////////////////////////////////////////////////////////////
const int KMER_PACKED_EXT_POS = KMER_PACKED_LENGTH_WITH_EXT - 1;

unsigned char *packedKmerValue(packed_kmer_t *kmer) {
  return kmer->data;
}

/* True if @a and @b agree on their kmer values (not necessarily their
 * extensions). */
bool equalsOnKmer(const packed_kmer_t *a, const packed_kmer_t *b) {
  return memcmp(packedKmerValue(a), packedKmerValue(b), KMER_PACKED_LENGTH * sizeof(unsigned char)) == 0; 
}

bool isBackwardStartPacked(const packed_kmer_t *kmer) {
  return (kmer->data[KMER_PACKED_EXT_POS] & (1 << 3)) != 0;
}

bool isForwardStartPacked(const packed_kmer_t *kmer) {
  return (kmer->data[KMER_PACKED_EXT_POS] & (1 << 2)) != 0;
}

unsigned char backwardExtensionPacked(const packed_kmer_t *kmer) {
  if (isBackwardStartPacked(kmer)) {
    return START_CHAR;
  }
  unsigned char code = (kmer->data[KMER_PACKED_EXT_POS] & (0x3 << 6)) >> 6;
  return convertPackedCodeToMer(code);
}

unsigned char forwardExtensionPacked(const packed_kmer_t *kmer) {
  if (isForwardStartPacked(kmer)) {
    return START_CHAR;
  }
  unsigned char code = (kmer->data[KMER_PACKED_EXT_POS] & (0x3 << 4)) >> 4;
  return convertPackedCodeToMer(code);
}

void toPacked(const unpacked_kmer_t *rawKmer, packed_kmer_t *packedKmer) {
	// See formatting above.
  unsigned char *rawData = kmerValue(rawKmer);
  unsigned char *packedData = packedKmerValue(packedKmer);
  packSequence(rawData, packedData, KMER_LENGTH);
	packedData[KMER_PACKED_EXT_POS] = 0;
	unsigned char backExt = backwardExtensionUnpacked(rawKmer);
	unsigned char forExt = forwardExtensionUnpacked(rawKmer);
	packedData[KMER_PACKED_EXT_POS] |= (
		(backExt == START_CHAR ? 0 : convertMerToPackedCode(backExt))
		<< 6);
	packedData[KMER_PACKED_EXT_POS] |= (
		(forExt == START_CHAR ? 0 : convertMerToPackedCode(forExt))
		<< 4);
	packedData[KMER_PACKED_EXT_POS] |= (
		(backExt == START_CHAR ? 1 : 0)
		<< 3);
	packedData[KMER_PACKED_EXT_POS] |= (
		(forExt == START_CHAR ? 1 : 0)
		<< 2);
}

int64_t hashPackedKmer(const packed_kmer_t *packedKmer, int64_t hashTableSize) {
  return hashKmerBytes(hashTableSize, packedKmerValue(packedKmer));
}

///////////////////////////////////////////////////////////////////////////////
// end packed_kmer_t
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// begin packed_kmer_list_t
///////////////////////////////////////////////////////////////////////////////

/** A linked list designed to hold packed kmers somewhat efficiently. */
typedef struct packed_kmer_list_t packed_kmer_list_t;
struct packed_kmer_list_t {
  packed_kmer_t *kmer;
  packed_kmer_list_t *next;
};

int packedListSize(const packed_kmer_list_t *list) {
  int size = 0;
  while (list != NULL) {
    size++;
    list = list->next;
  }
  return size;
}

/* (*dst) is now a flat array of packed_kmer_t. */
int makeFlatCopy(packed_kmer_t **dst, const packed_kmer_list_t *list) {
  const int size = packedListSize(list);
  packed_kmer_t *flatCopy = (packed_kmer_t *) malloc(size*sizeof(packed_kmer_t));
  int i = 0;
  while (list != NULL) {
    memcpy(&flatCopy[i], list->kmer, sizeof(packed_kmer_t));
    list = list->next;;
    i++;
  }
  *dst = flatCopy;
  return size;
}

typedef struct kmer_memory_heap_t {
  packed_kmer_list_t *heap;
  int64_t nextHeapLocation;
} kmer_memory_heap_t;

packed_kmer_list_t *addFrontWithHeap(const packed_kmer_list_t *list, const packed_kmer_t *kmer, kmer_memory_heap_t *heap) {
  packed_kmer_list_t *newFront = &heap->heap[heap->nextHeapLocation];
  newFront->kmer = kmer;
  newFront->next = list;
  heap->nextHeapLocation++;
  return newFront;
}

void setFrontWithHeap(packed_kmer_list_t **list, const packed_kmer_t *kmer, kmer_memory_heap_t *heap) {
  *list = addFrontWithHeap(*list, kmer, heap);
}

kmer_memory_heap_t *makeMemoryHeap(int64_t size) {
  packed_kmer_list_t *kmers = (packed_kmer_list_t *) malloc(sizeof(packed_kmer_list_t)*size);
  kmer_memory_heap_t *heap = (kmer_memory_heap_t *) malloc(sizeof(kmer_memory_heap_t));
  heap->heap = kmers;
  heap->nextHeapLocation = 0;
  return heap;
}

void freeMemoryHeap(kmer_memory_heap_t *heap) {
  free(heap->heap);
  free(heap);
}

///////////////////////////////////////////////////////////////////////////////
// end packed_kmer_list_t
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// begin unpacked_kmer_list_t
///////////////////////////////////////////////////////////////////////////////

/** A simple linked list with moderately inefficient operations.  Designed to
  * hold a short list of start kmers. */
typedef struct unpacked_kmer_list_t unpacked_kmer_list_t;
struct unpacked_kmer_list_t {
  unpacked_kmer_t *kmer;
  unpacked_kmer_list_t *next;
};

unpacked_kmer_list_t *makeUnpackedList() {
  return NULL;
}

int unpackedListSize(const unpacked_kmer_list_t *list) {
  int size = 0;
  while (list != NULL) {
    size++;
    list = list->next;
  }
  return size;
}

unpacked_kmer_list_t *addUnpackedFront(const unpacked_kmer_list_t *list, const unpacked_kmer_t *kmer) {
  unpacked_kmer_list_t *newFront = (unpacked_kmer_list_t *) malloc(sizeof(unpacked_kmer_list_t));
  newFront->kmer = kmer;
  newFront->next = list;
  return newFront;
}

void freeUnpackedKmerList(unpacked_kmer_list_t *list) {
  while (list != NULL) {
    unpacked_kmer_list_t *next = list->next;
    free(list);
    list = next;
  }
}

///////////////////////////////////////////////////////////////////////////////
// end packed_kmer_list_t
///////////////////////////////////////////////////////////////////////////////


#endif // KMER_PACKING_H