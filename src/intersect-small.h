#ifndef INTERSECT_SMALL_INCLUDED
#define INTERSECT_SMALL_INCLUDED

#include "bool.h"

extern unsigned int *
Intersect_small (int *ndiagonals,
		 const unsigned int *set1, const int length1, int diagterm1,
		 const unsigned int *set2, const int length2, int diagterm2,
		 bool alignp);

extern unsigned int *
Intersect_exact_old (int *ndiagonals,
		     const unsigned int *positionsa, const int npositionsa, int diagterma,
		     const unsigned int *positionsb, const int npositionsb, int diagtermb);
#endif
