#ifndef INTERSECT_APPROX_INDICES_UINT4_INCLUDED
#define INTERSECT_APPROX_INDICES_UINT4_INCLUDED

extern int *
Intersect_approx_indices_uint4 (int *ndiagpairs,
				const unsigned int *set1, const int length1, int diagterm1,
				const unsigned int *set2, const int length2, int diagterm2,
				const unsigned int below_slop, const unsigned int above_slop);

extern void
Intersect_approx_indices_uint4_setup ();

#endif
