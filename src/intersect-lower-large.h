#ifndef INTERSECT_LOWER_LARGE_INCLUDED
#define INTERSECT_LOWER_LARGE_INCLUDED

#include "univcoord.h"
#include "types.h"

/* GSNAP uses intersect-lower-small and GSNAPL uses intersect-lower-large for this procedure */

extern int
Intersect_lower (Univcoord_T *diagonals,
		 const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		 const Univcoord_T *set2, const int length2,
		 const Univcoord_T slop, const Univcoord_T insertion_slop);

extern int
Intersect_approx_lower_old (Univcoord_T *diagonals,
			    unsigned char *positions1_high, UINT4 *positions1_low, int npositions1, int diagterm1,
			    Univcoord_T *positions0, int npositions0,
			    const Univcoord_T maxdistance);

extern void
Intersect_lower_setup ();

#endif
