/* $Id: 848f42806b452454277bc89c2bce25da64384206 $ */
#ifndef TRANSCRIPTOME_SEARCH_INCLUDED
#define TRANSCRIPTOME_SEARCH_INCLUDED

#include "bool.h"
#include "chrnum.h"		/* For Trnum_T */
#include "list.h"
#include "indexdb.h"
#include "ef64.h"

#include "iit-read-univ.h"	/* For transcript_iit */
#include "transcriptome.h"
#include "compress.h"
#include "shortread.h"
#include "genomebits.h"

#include "method.h"
#include "stage1hr.h"
#include "indel.h"
#include "mergeinfo.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "hitlistpool.h"


extern int
Transcriptome_exact1 (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
		      Trcoord_T **_sense_trdiagonals, int *n_sense_trdiagonals,
		      Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
		      Trcoord_T **_antisense_trdiagonals, int *n_antisense_trdiagonals,
		      Stage1_T stage1, int querylength);

extern int
Transcriptome_anypair (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
		       Trcoord_T **_sense_trdiagonals, int **sense_tstarts, int **sense_tends, int *n_sense_trdiagonals,
		       Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
		       Trcoord_T **_antisense_trdiagonals, int **antisense_tstarts, int **antisense_tends,
		       int *n_antisense_trdiagonals, Stage1_T stage1, int querylength);

extern int
Transcriptome_prevalent (Trnum_T **sense_trnums, Trcoord_T **sense_troffsets, Trcoord_T **sense_trhighs,
			 Trcoord_T **_sense_trdiagonals, int **sense_tstarts, int **sense_tends, int *n_sense_trdiagonals,
			 Trnum_T **antisense_trnums, Trcoord_T **antisense_troffsets, Trcoord_T **antisense_trhighs,
			 Trcoord_T **_antisense_trdiagonals, int **antisense_tstarts, int **antisense_tends,
			 int *n_antisense_trdiagonals,
			 Stage1_T stage1, int querylength, Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern void
Transcriptome_search_setup (int index1part_tr_in, int index1interval_tr_in,
			    Indexdb_T tr_indexdb_in, Transcriptome_T transcriptome_in,
			    EF64_T transcript_ef64_in, Genomebits_T genomebits_in,
			    Genomebits_T genomebits_alt_in, Genomebits_T transcriptomebits_in);

#endif

