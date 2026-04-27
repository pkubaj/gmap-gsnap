/* $Id: sedgesort.h 225220 2022-11-01 22:46:01Z twu $ */
#ifndef SEDGESORT_INCLUDED
#define SEDGESORT_INCLUDED

#include "types.h"

extern void
Sedgesort_uint4 (register unsigned int array[], register int len);
#ifdef LARGE_GENOMES
extern void
Sedgesort_uint8 (register UINT8 array[], register int len);
#endif

extern int *
Sedgesort_order_uint4 (register unsigned int array[], register int len);
extern int *
Sedgesort_order_int (register int array[], register int len);
extern int *
Sedgesort_order_int2 (register int array1[], register int array2[], register int len);
#ifdef HAVE_64_BIT
extern int *
Sedgesort_order_uint8 (register UINT8 array[], register int len);
#endif

#endif

