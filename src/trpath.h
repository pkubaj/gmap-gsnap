/* $Id: 96dfa48e21a7a7ecae1ca693e73c4ae1c5eb50d2 $ */
#ifndef TRPATH_INCLUDED
#define TRPATH_INCLUDED

typedef struct Trpath_T *Trpath_T;

#include "bool.h"
#include "types.h"
#include "transcript.h"
#include "list.h"
#include "chrnum.h"
#include "univcoord.h"
#include "intlist.h"
#include "uintlist.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"

#include "compress.h"
#include "indel.h"
#include "splice.h"
#include "knownsplicing.h"

#include "path.h"
#include "method.h"


#define T Trpath_T
struct T {
  bool tplusp;

  Trnum_T trnum;
  Trcoord_T troffset;
  Trcoord_T trhigh;

  Chrnum_T chrnum;
  int transcript_genestrand;

  /* Note: We cannot have ref_nmismatches, because we don't have an
     alt version of transcriptome */
  Uintlist_T trdiagonals;
  Intlist_T endpoints;
  Intlist_T nmismatches;
  List_T junctions;

  int found_score;
  int end_score;
  int total_ninserts;
  Method_T method;
};


extern void
Trpath_free (T *old, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool);

extern void
Trpath_gc (List_T *list, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Trpathpool_T trpathpool,
	   Hitlistpool_T hitlistpool);

extern int
Trpath_trnum_cmp (const void *x, const void *y);

/* Called by combine_leftright_paths in Trpath_solve_from_diagonals */
extern T
Trpath_create (Intlist_T endpoints, Uintlist_T trdiagonals, Intlist_T nmismatches,
	       List_T junctions, bool tplusp,
	       Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
	       Chrnum_T chrnum, int transcript_genestrand,
	       Trpathpool_T trpathpool, int found_score, int total_ninserts, Method_T method);

extern T
Trpath_new_for_tstart_extension (Trcoord_T univdiagonal, int tstart, int tend, int nmismatches,
				 bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
				 Chrnum_T chrnum, int transcript_genestrand,
				 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				 Trpathpool_T trpathpool, Method_T method);

extern T
Trpath_new_for_tend_extension (Trcoord_T univdiagonal, int tstart, int tend, int nmismatches,
			       bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			       Chrnum_T chrnum, int transcript_genestrand,
			       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			       Trpathpool_T trpathpool, Method_T method);

/* Called by approx method */
extern T
Trpath_new_from_ends (Trcoord_T trdiagonal5, int tstart5, int tend5,
		      Trcoord_T trdiagonal3, int tstart3, int tend3,
		      bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
		      Chrnum_T chrnum, int transcript_genestrand,
		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Listpool_T listpool, Trpathpool_T trpathpool,
		      int found_score, int total_ninserts, Method_T method);

/* Called by exact method */
extern T
Trpath_new_exact (Trcoord_T trdiagonal, int tstart, int tend, int nmismatches,
		  bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
		  Chrnum_T chrnum, int transcript_genestrand,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		  Trpathpool_T trpathpool, int found_score, Method_T method);

extern T
Trpath_reverse (T this, bool expect_fwd_p);

extern List_T
Trpath_filter_best (List_T *subopt, List_T trpaths, int found_score);

extern void
Trpath_print (T this);


#undef T
#endif
