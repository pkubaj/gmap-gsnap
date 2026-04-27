/* $Id: merge-uint8.h 225220 2022-11-01 22:46:01Z twu $ */
#ifndef MERGE_UINT8_INCLUDED
#define MERGE_UINT8_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "types.h"

extern UINT8 *
Merge_uint8 (UINT8 *__restrict__ dest, UINT8 *__restrict__ A,
	     UINT8 *__restrict__ B, int nA, int nB);

#endif
