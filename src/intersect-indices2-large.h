#ifndef INTERSECT_INDICES2_LARGE_INCLUDED
#define INTERSECT_INDICES2_LARGE_INCLUDED

#include "univcoord.h"
#include "types.h"

/* Kmer_segment wants the indices of set2 */
extern int
Intersect_indices2_large (int *indices,
			 const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
			 const Univcoord_T *set2, const int length2, int diagterm2);

extern void
Intersect_indices2_large_setup ();

#endif
