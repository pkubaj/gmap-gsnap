/* $Id: 3136f868f96702ccf9ac09a95cabe75230eb5d57 $ */
#ifndef TR_EXTENSION_SEARCH_INCLUDED
#define TR_EXTENSION_SEARCH_INCLUDED

#include "types.h"
#include "genomicpos.h"
#include "mode.h"
#include "method.h"

#include "list.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "genomebits.h"
#include "compress.h"
#include "shortread.h"

#include "trdiag.h"

#include "stage1hr.h"
#include "indexdb.h"

#include "trdiagpool.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "hitlistpool.h"


#define T Tr_elt_T
typedef struct T *T;
struct T {
  int min_qstart;
  int max_qend;

  int nmatches;

  Trdiag_T *all_trdiags;
  int n_all_trdiags;

  Trdiag_T *trdiags;	/* Filtered by binary search to generate lowi and highi */
  int ntrdiags;

  int lowi;
  int highi;
};


extern void
Tr_extension_search_setup (Transcriptome_T transcriptome_in, Trcoord_T transcriptomelength_in, EF64_T transcript_ef64_in,
			   Genomebits_T transcriptomebits_in, Indexdb_T tr_indexdb_in,
			   int index1part_tr_in, int max_insertionlen_in, int max_deletionlen_in,
			   int maxpaths_search_in);

extern void
Tr_elt_gc (List_T *set, Listpool_T listpool, Trdiagpool_T trdiagpool);

extern void
Tr_extension_search (int *found_score, int sufficient_score,
		     List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
		     List_T *sense_trpaths, List_T *antisense_trpaths,

		     Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     bool paired_end_p, bool first_read_p,

		     Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
		     Hitlistpool_T hitlistpool,

		     int nmismatches_allowed, int genestrand, Method_T method);

#undef T
#endif


