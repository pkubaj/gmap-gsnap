#ifndef INTERSECT_CONCORDANCE_UINT8_INCLUDED
#define INTERSECT_CONCORDANCE_UINT8_INCLUDED

#include "univcoord.h"

extern int
Intersect_concordance (int **result,
		       const Univcoord_T *set1, const int length1,
		       const Univcoord_T *set2, const int length2,
		       const Univcoord_T slop1_neg, const Univcoord_T slop1_pos,
		       const Univcoord_T slop2_neg, const Univcoord_T slop2_pos);

extern void
Intersect_concordance_setup ();

#endif


