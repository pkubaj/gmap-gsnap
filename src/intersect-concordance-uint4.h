#ifndef INTERSECT_CONCORDANCE_UINT4_INCLUDED
#define INTERSECT_CONCORDANCE_UINT4_INCLUDED

extern int
Intersect_concordance (int **result,
		       const unsigned int *set1, const int length1,
		       const unsigned int *set2, const int length2,
		       const unsigned int slop1_neg, const unsigned int slop1_pos,
		       const unsigned int slop2_neg, const unsigned int slop2_pos);

extern void
Intersect_concordance_setup ();

#endif


