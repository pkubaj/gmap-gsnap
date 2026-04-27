#ifndef INTERSECT_LARGE_INCLUDED
#define INTERSECT_LARGE_INCLUDED

#include "types.h"

extern UINT8 *
Intersect_large (int *ndiagonals,
		 const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		 const unsigned char *set2_high, const UINT4 *set2_low, const int length2, int diagterm2);

extern UINT8 *
Intersect_exact_large_old (int *ndiagonals,
			   const unsigned char *positionsa_high, const UINT4 *positionsa_low,
			   const int npositionsa, int diagterma,
			   const unsigned char *positionsb_high, const UINT4 *positionsb_low,
			   const int npositionsb, int diagtermb);
#endif
