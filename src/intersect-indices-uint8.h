#ifndef INTERSECT_INDICES_UINT8_INCLUDED
#define INTERSECT_INDICES_UINT8_INCLUDED

#include "univcoord.h"

/* Segment_search wants the indices of set2 */
extern int
Intersect_indices_uint8 (int *indices,
			 const Univcoord_T *set1, const int length1, int diagterm1,
			 const Univcoord_T *set2, const int length2, int diagterm2);

extern void
Intersect_indices_uint8_setup ();

#endif
