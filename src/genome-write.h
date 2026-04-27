/* $Id: genome-write.h 224922 2022-10-08 13:23:45Z twu $ */
#ifndef GENOME_WRITE_INCLUDED
#define GENOME_WRITE_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "types.h"
#include "genomicpos.h"

extern void
Genome_write_comp32 (char *filename, char *fileroot, FILE *input, 
		     Univ_IIT_T contig_iit, Univ_IIT_T chromosome_iit,
		     bool uncompressedp, bool rawp, bool writefilep,
		     Univcoord_T genomelength, int index1part, int nmessages,
		     bool force_revcomp_p);

#endif
