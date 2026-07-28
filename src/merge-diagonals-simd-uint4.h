#ifndef MERGE_DIAGONALS_SIMD_UINT4_INCLUDED
#define MERGE_DIAGONALS_SIMD_UINT4_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "types.h"
#include "univcoord.h"
#include "mergeinfo.h"
#include "merge-method.h"	/* For USE_HEAP_MERGE or USE_SIMD_MERGE */

#ifdef USE_SIMD_MERGE
extern UINT4 *
Merge_diagonals (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		 int *diagterm_array, int nstreams, Mergeinfo_uint4_T mergeinfo);
#endif

#ifdef USE_SIMD_MERGE
extern UINT4 *
Merge_diagonals_uint4 (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint4_T mergeinfo);
#endif

#endif


