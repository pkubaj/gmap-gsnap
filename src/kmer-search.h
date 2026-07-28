/* $Id: df710e2c49bc3a730a0cf145144cf3421734a394 $ */
#ifndef KMER_SEARCH_INCLUDED
#define KMER_SEARCH_INCLUDED

#include "method.h"
#include "list.h"
#include "indexdb.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "transcriptome.h"
#include "compress.h"
#include "shortread.h"
#include "genomebits.h"
#include "indel.h"

#include "auxinfo.h"
#include "pathstore.h"

#include "bool.h"
#include "pass.h"
#include "univdiag.h"
#include "stage1hr.h"
#include "mergeinfo.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "univdiagpool.h"
#include "auxinfopool.h"
#include "intlistpool.h"
#include "uintlistpool.h"
#include "univcoord.h"
#include "pathpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "transcriptpool.h"
#include "vectorpool.h"
#include "listpool.h"
#include "hitlistpool.h"


/* Does not take paired_end_p as a parameter.  ? Generates both sense and antisense */
extern int
Kmer_exact1 (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	     Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	     Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	     Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
	     Stage1_T stage1, int querystart, int queryend, int querylength,
	     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool);

extern Auxinfo_T
Kmer_compute_auxinfo_univdiags (Univcoord_T main_univdiagonal, int qstart, int qend, int i,
				Univcoord_T *exhaustive, int *qstarts, int *qends, int nexhaustive,
				Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
				Method_T method);

extern int
Kmer_segment (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	      Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	      Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	      Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
	      Stage1_T stage1, int querylength, EF64_T repetitive_ef64,
	      Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool);

extern void
Kmer_search_setup (int index1part_in, int index1interval_in,
		   Indexdb_T indexdb_in, EF64_T chromosome_ef64_in,
		   Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		   Univcoord_T genomelength_in,
		   int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist,
		   bool splicingp_in);
#endif
