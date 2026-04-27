#ifndef INTERSECT_INDICES_SMALL_INCLUDED
#define INTERSECT_INDICES_SMALL_INCLUDED

/* Segment_search wants the indices of set2 */
extern int
Intersect_indices_small (int *indices,
			 const unsigned int *set1, const int length1, int diagterm1,
			 const unsigned int *set2, const int length2, int diagterm2);

extern void
Intersect_indices_small_setup ();

#endif
