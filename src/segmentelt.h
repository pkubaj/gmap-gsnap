#ifndef SEGMENTELT_INCLUDED
#define SEGMENTELT_INCLUDED


typedef struct Segmentelt_T *Segmentelt_T;

#include "univcoord.h"
#include "auxinfo.h"
#include "pathstore.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "intlist.h"
#include "list.h"

#define SEGMENTELT_DEFAULT_TERMINAL_ANCHOR_MIN_SCORE 30
#define SEGMENTELT_DEFAULT_TERMINAL_ANCHOR_SLOP 2


struct Segmentelt_T {
  int qstart;
  int qend;
  Univcoord_T *univdiagonals;
  int ndiagonals;

  /* For dynamic programming. Allows equi-optimal paths */
  List_T *prev_elt_lists;
  Intlist_T *prev_diagonali_lists;
  int *scores;
  int *nsplices;

};


#include "univdiagpool.h"
#include "segmenteltpool.h"
#include "intlistpool.h"
#include "listpool.h"

#define T Segmentelt_T

extern void
Segmentelt_setup (int max_insertionlen, Chrpos_T shortsplicedist, int max_deletionlen_in,
		  int terminal_anchor_min_score_in, int terminal_anchor_slop_in);

extern void
Segmentelt_free (T *old, Segmenteltpool_T segmenteltpool,
		 Intlistpool_T intlistpool, Listpool_T listpool);

extern void
Segmentelt_gc (List_T set, Segmenteltpool_T segmenteltpool,
	       Intlistpool_T intlistpool, Listpool_T listpool);

extern T
Segmentelt_new (int qstart, int qend, Univcoord_T *univdiagonals,
		int ndiagonals, Segmenteltpool_T segmenteltpool);

extern T *
Segmentelt_new_array (int n);

extern void
Segmentelt_print (T this);

extern void
Segmentelt_dump (T this);

extern int
Segmentelt_cmp (const void *a, const void *b);

extern int **
Segmentelt_left_chains (int **nindices_by_nsplices, int max_nsplices,
			Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores, int *nunivdiagonals,
			T *eltarray, int narray, int querylength, bool retain_terminal_p,
			bool prune_equal_score_p,
			Univcoordlistpool_T univcoordlistpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Listpool_T listpool);

extern int **
Segmentelt_right_chains (int **nindices_by_nsplices, int max_nsplices,
			 Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores, int *nunivdiagonals,
			 T *eltarray, int narray, int querylength, bool retain_terminal_p,
			 bool prune_equal_score_p,
			 Univcoordlistpool_T univcoordlistpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			 Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Listpool_T listpool);

#undef T
#endif
