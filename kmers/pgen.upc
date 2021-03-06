#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>

#include "packingDNAseq.h"
#include "kmer_packing.h"
#include "hash_shuffle.h"

/** Read kmers for this thread from @inputFileName. */
unpacked_kmer_t *readLocalKmers(const char *inputFileName, int nKmersReadTypically, int nKmersToRead, int threadIdx) {
  printf("Reading %d kmers on thread %d.\n", nKmersToRead, MYTHREAD);
  const int localCharsToRead = nKmersToRead * LINE_SIZE;
  FILE *inputFile = fopen(inputFileName, "r");
  fseek(inputFile, threadIdx*nKmersReadTypically*LINE_SIZE, SEEK_SET);
  unsigned char *localReadBuffer = (unsigned char*) malloc(localCharsToRead * sizeof(unsigned char));
  fread(localReadBuffer, sizeof(unsigned char), localCharsToRead, inputFile);
  fclose(inputFile);
  unpacked_kmer_t *localReadKmers = (unpacked_kmer_t *) malloc(nKmersToRead*sizeof(unpacked_kmer_t));
  for (int i = 0; i < nKmersToRead; i++) {
    //TODO: Have to copy one by one out of the buffer into the structs, because 
    // the structs might not be byte-identical to the arrays they contain (the
    // compiler may add padding to word-align them).  Perhaps this could be
    // optimized if it is a problem.
    fromRawChars(&localReadKmers[i], &localReadBuffer[i*LINE_SIZE]);
#ifdef CHECKSUM
    checkUnpacked(&localReadKmers[i]);
#endif
  }
  free(localReadBuffer);
  printf("Finished reading %d kmers on thread %d.\n", nKmersToRead, MYTHREAD);
  return localReadKmers;
}

int filterBackwardStartKmers(packed_kmer_t **dst, const unpacked_kmer_t *kmers, int nKmers) {
  unpacked_kmer_list_t *list = makeUnpackedList();
  for (int i = 0; i < nKmers; i++) {
    if (isBackwardStartUnpacked(&kmers[i])) {
      list = addUnpackedFront(list, &kmers[i]);
#ifdef CHECKSUM
      checkUnpacked(list->kmer);
#endif
    }
  }
  int numStartKmers = makeFlatPackedCopy(dst, list);
  freeUnpackedKmerList(list);
  return numStartKmers;
}

/* @param dst should already be allocated with at least maxContigSize chars. 
 *   It will be set to the contig starting at @start. 
 * @return the length of the contig.
 */
int buildContig(unsigned char *dst, const packed_kmer_t *start, const shared_hash_table_t *table) {
  /* Initialize current contig with the seed content */
  unpacked_kmer_t startUnpacked;
  unpack(&startUnpacked, start);
  memcpy(dst, kmerValue(&startUnpacked), KMER_LENGTH * sizeof(unsigned char));
  
  int posInContig = KMER_LENGTH;
  char forwardExt = forwardExtensionPacked(start);
  packed_kmer_t packedStorage;
  unpacked_kmer_t unpackedStorage;
  //HACK
  unpackedStorage.data[RAW_KMER_BACKWARD_EXT_POS] = 'A';
  unpackedStorage.data[RAW_KMER_FORWARD_EXT_POS] = 'A';

  /* Keep adding bases while not finding a terminal node. */
  while (forwardExt != 'F') {
    dst[posInContig] = forwardExt;
    posInContig++;
    /* dst[posInContig-KMER_LENGTH] is the start of the last k-mer in the
    * current contig. */
    //FIXME: This is slow.  Should maintain the packed form only.
    memcpy(&(unpackedStorage.data), &dst[posInContig-KMER_LENGTH], KMER_LENGTH*sizeof(unsigned char));
#ifdef CHECKSUM
    unpackedStorage.checksum = unpackedChecksum(&unpackedStorage);
#endif
    toPacked(&unpackedStorage, &packedStorage);
    forwardExt = lookupForwardExtension(table, &packedStorage);
  }

  dst[posInContig] = '\0';
  return posInContig;
}

int main(int argc, char *argv[]) {
  double inputTime=0.0, constrTime=0.0, traversalTime=0.0;
  
  init_LookupTable();
  const char *inputFileName = argv[1];
  char outputFileName[256];
  snprintf(outputFileName, sizeof outputFileName, "pgen-p%dt%d.out", MYTHREAD, THREADS);

  /** Read input **/
  inputTime -= gettime();
  
  /* Extract the number of k-mers in the input file */
  const int nKmers = getNumKmersInUFX(inputFileName);
  const int nThreads = THREADS;
  const int threadIdx = MYTHREAD;
  const bool isLastThread = (threadIdx == (nThreads - 1));
  const int nKmersReadTypically = divRoundUp(nKmers, nThreads);
  const int64_t nHashBuckets = divRoundUp(nKmers * LOAD_FACTOR, nThreads) * nThreads;
  
  /* The last thread might get fewer kmers. */
  const int nKmersReadLocally = isLastThread ?
    nKmers - nKmersReadTypically*(nThreads-1) :
    nKmersReadTypically;

  /* Read the kmers from the input file into a local buffer. */
  const unpacked_kmer_t *localRawKmers = readLocalKmers(inputFileName, nKmersReadTypically, nKmersReadLocally, threadIdx);
  
  inputTime += gettime();

  /** Graph construction **/
  constrTime -= gettime();
  
  /* Build the local start list. */
  packed_kmer_t *startList;
  int numStartKmers = filterBackwardStartKmers(&startList, localRawKmers, nKmersReadLocally);
  
  /** Get the kmers to handle on this thread. */
  packed_kmer_t *localKmers;
  int numLocalKmers = hashShuffleKmers(&localKmers, localRawKmers, nKmersReadLocally, threadIdx, nThreads, nHashBuckets);
  
  free(localRawKmers);
  
  printf("Read %d/%d kmers on thread %d\n", nKmersReadLocally, nKmers, threadIdx);
  printf("Handling %d/%d kmers on thread %d\n", numLocalKmers, nKmers, threadIdx);
  
  /** Build the distributed hash table. */
  shared_hash_table_t *table = buildSharedHashTable(localKmers, numLocalKmers, threadIdx, nThreads, nHashBuckets);
  
  free(localKmers);
  
  upc_barrier;
  constrTime += gettime();
  
  printf("Finished building the hashtable on thread %d.  Now traversing the graph.\n", threadIdx);

  /** Graph traversal **/
  traversalTime -= gettime();
  unsigned char *contig = (unsigned char *) malloc(MAXIMUM_CONTIG_SIZE*sizeof(unsigned char));
  FILE *outputFile = fopen(outputFileName, "w");
  
  for (int startKmerIdx = 0; startKmerIdx < numStartKmers; startKmerIdx++) {
    buildContig(contig, &startList[startKmerIdx], table);
    fprintf(outputFile, "%s\n", contig);
  }
  fclose(outputFile);
  free(contig);
  free(startList);
  free(table);
  upc_barrier;
  traversalTime += gettime();
  
  //FIXME: Should technically free the stuff in the hash table...

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
