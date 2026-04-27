/* $Id: merge-diagonals-simd-uint4.h 224755 2021-12-13 00:42:15Z twu $ */
#ifndef MERGE_DIAGONALS_SIMD_INCLUDED
#define MERGE_DIAGONALS_SIMD_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "types.h"
#include "univcoord.h"
#include "mergeinfo.h"

extern UINT4 *
Merge_diagonals (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		 int *diagterm_array, int nstreams, Mergeinfo_uint4_T mergeinfo);

extern UINT4 *
Merge_diagonals_uint4 (int *nelts1, UINT4 **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint4_T mergeinfo);

#endif


