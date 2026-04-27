/* $Id: merge-uint4.h 223349 2020-10-28 02:49:25Z twu $ */
#ifndef MERGE_UINT4_INCLUDED
#define MERGE_UINT4_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

extern unsigned int *
Merge_uint4 (unsigned int *__restrict__ dest, unsigned int *__restrict__ A,
	     unsigned int *__restrict__ B, int nA, int nB);

#endif


