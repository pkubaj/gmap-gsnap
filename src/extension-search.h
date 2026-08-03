#ifndef EXTENSION_SEARCH_INCLUDED
#define EXTENSION_SEARCH_INCLUDED

#include "method.h"
#include "types.h"
#include "genomicpos.h"
#include "mode.h"
#include "pass.h"

#include "auxinfo.h"
#include "pathstore.h"

#include "list.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "genomebits.h"
#include "compress.h"
#include "shortread.h"

#include "stage1hr.h"
#include "indexdb.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "segmentelt.h"

#include "univdiagpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "univcoord.h"
#include "intlistpool.h"
#include "listpool.h"


#define T Elt_T
typedef struct T *T;
struct T {
  int min_qstart;
  int max_qend;

  int nmatches;

  Univdiag_T *all_univdiags;
  int n_all_univdiags;

  Univdiag_T *univdiags;	/* Filtered by binary search to generate lowi and highi */
  int nunivdiags;

  int lowi;
  int highi;
};


extern void
Extension_search_setup (Mode_T mode,
			Univcoord_T genomelength_in, int circular_typeint_in, bool *circularp_in, EF64_T chromosome_ef64_in,
			Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
			Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in,
			int index1part_in, int index1interval_in, int maxpaths_search_in,
			int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist);

extern void
Elt_gc (List_T *set, Listpool_T listpool);

extern void
Extension_search (int ***indices_by_nsplices_gplus, int **nindices_by_nsplices_gplus,
		  int ***indices_by_nsplices_gminus, int **nindices_by_nsplices_gminus,

		  Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
		  Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
		  Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
		  Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,

		  Stage1_T stage1, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  int querylength, int max_nsplices, bool first_read_p, bool retain_terminal_p,
		  bool prune_equal_score_p,

		  Segmenteltpool_T segmenteltpool,
		  Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		  Univcoordlistpool_T univcoordlistpool, Univdiagpool_T univdiagpool,
		  Intlistpool_T intlistpool, Listpool_T listpool);

#undef T
#endif

