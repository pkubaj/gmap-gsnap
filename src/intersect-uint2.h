#ifndef INTERSECT_UINT2_INCLUDED
#define INTERSECT_UINT2_INCLUDED

#include "univcoord.h"

extern int
Intersect_uint2_scalar (Univcoord_T *diagonals, unsigned short **set2_ptr,
			unsigned short *set1, const int length1, int diagterm1,
			unsigned short *set2, const int length2, int diagterm2,
			Univcoord_T region_term);

extern int
Intersect_uint2 (Univcoord_T *diagonals, unsigned short *localdb_alloc,
		 unsigned short *set1, const int length1, int diagterm1,
		 unsigned short *set2, const int length2, int diagterm2,
		 Univcoord_T region_term);

#endif
