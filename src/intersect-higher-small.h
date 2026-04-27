#ifndef INTERSECT_HIGHER_SMALL_INCLUDED
#define INTERSECT_HIGHER_SMALL_INCLUDED

#include "univcoord.h"
#include "types.h"

/* GSNAP uses intersect-higher-small and GSNAPL uses intersect-higher-large for this procedure */

extern int
Intersect_higher (Univcoord_T *diagonals,
		  const UINT4 *set1, const int length1, int diagterm1,
		  const Univcoord_T *set2, const int length2,
		  const Univcoord_T slop, const Univcoord_T insertion_slop);

extern int
Intersect_approx_higher_old (unsigned int *diagonals,
			     unsigned int *positions1, int npositions1, int diagterm1,
			     unsigned int *positions0, int npositions0,
			     const unsigned int maxdistance);

extern void
Intersect_higher_setup ();

#endif
