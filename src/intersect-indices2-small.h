#ifndef INTERSECT_INDICES2_SMALL_INCLUDED
#define INTERSECT_INDICES2_SMALL_INCLUDED

/* Kmer_segment wants the indices of set2 */
extern int
Intersect_indices2_small (int *indices,
			 const unsigned int *set1, const int length1, int diagterm1,
			 const unsigned int *set2, const int length2, int diagterm2);

extern void
Intersect_indices2_small_setup ();

#endif
