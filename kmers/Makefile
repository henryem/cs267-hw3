## NOTES ##
# I cannot build BUPC on OSX 10.10.  This Makefile is targeted at Ubuntu 14.04.1.

CC = g++
UPCC = upcc

KMER_LENGTH 		= 19
KMER_EXT_LENGTH = 2
KMER_PACKED_LENGTH 	= $(shell echo $$((($(KMER_LENGTH)+3)/4)))
KMER_PACKED_LENGTH_WITH_EXT 	= $(shell echo $$(($(KMER_LENGTH)+1)))

CFLAGS  = -std=gnu99
ifdef DEBUG
	CFLAGS += "-d -g"
else
	CFLAGS += "-O3"
endif

DEFINE 	= -DKMER_LENGTH=$(KMER_LENGTH) -DKMER_PACKED_LENGTH=$(KMER_PACKED_LENGTH) -DKMER_PACKED_LENGTH_WITH_EXT=$(KMER_PACKED_LENGTH_WITH_EXT)
HEADERS	= contig_generation.h kmer_hash.h packingDNAseq.h kmer_packing.h
LIBS	=

TARGETS	= serial pgen

all: 	$(TARGETS)

serial: serial.c $(HEADERS)
		$(CC) $(CFLAGS) -o $@ $< $(DEFINE) $(LIBS)

pgen:	pgen.upc $(HEADERS)
		$(UPCC) $(UPCFLAGS) -Wc,"$(CFLAGS)" -o $@ $< $(DEFINE) $(LIBS)

GRADING_PACKAGE_NAME = derek_henry_hw3.tar.gz
grading-package: $(GRADING_PACKAGE_NAME)
$(GRADING_PACKAGE_NAME): *.c *.upc *.h MAKEFILE members.txt report.pdf
	tar cvzf $(GRADING_PACKAGE_NAME) $^
report.pdf: ../report/report.pdf
	cp $< $@

clean :
	rm -f *.o
	rm -rf $(TARGETS)
