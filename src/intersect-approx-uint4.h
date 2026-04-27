#ifndef INTERSECT_APPROX_UINT4_INCLUDED
#define INTERSECT_APPROX_UINT4_INCLUDED

extern unsigned int *
Intersect_approx_uint4 (int *ndiagpairs,
			const unsigned int *set1, const int length1, int diagterm1,
			const unsigned int *set2, const int length2, int diagterm2,
			const unsigned int below_slop, const unsigned int above_slop);

extern void
Intersect_approx_uint4_setup ();

#endif
