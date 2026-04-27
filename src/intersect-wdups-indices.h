#ifndef INTERSECT_WDUPS_INDICES_INCLUDED
#define INTERSECT_WDUPS_INDICES_INCLUDED

extern int *
Intersect_wdups_indices (int *ndiagpairs,
			 const unsigned int *set1, const int length1,
			 const unsigned int *set2, const int length2);

extern void
Intersect_wdups_indices_setup ();

#endif
