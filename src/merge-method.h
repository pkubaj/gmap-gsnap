/* $Id: 5cde2b3ade1ecc30e8194a7ed99b9805531eb6d5 $ */
#ifndef MERGE_METHOD_INCLUDED
#define MERGE_METHOD_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif


#if 1

/* Original algorithm */
#ifdef HAVE_AVX2
#define USE_SIMD_MERGE 1
#elif defined(HAVE_AVX512)
#define USE_SIMD_MERGE 1
#else
#define USE_HEAP_MERGE 1
#endif

#else

/* Might even be faster, since it is O(n log k) */
#define USE_HEAP_MERGE 1

#endif

#endif
