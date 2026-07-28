#ifndef MERGE_DIAGONALS_HEAP_UINT8_INCLUDED
#define MERGE_DIAGONALS_HEAP_UINT8_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "types.h"
#include "univcoord.h"
#include "mergeinfo.h"
#include "merge-method.h"	/* For USE_HEAP_MERGE or USE_SIMD_MERGE */

/* kmer-search.c calls Merge_diagonals_large with stream_high_list and
   stream_low_list.  localdb.c calls Merge_diagonals_uint4 and
   Merge_diagonals_uint8 */

#ifdef USE_HEAP_MERGE
extern Univcoord_T *
Merge_diagonals_large (int *nelts, unsigned char **stream_high_array, UINT4 **stream_low_array,
		       int *streamsize_array, int *diagterm_array, int nstreams,
		       Mergeinfo_uint8_T mergeinfo);
#endif

#ifdef USE_HEAP_MERGE
extern Univcoord_T *
Merge_diagonals_uint8 (int *nelts, Univcoord_T **stream_array, int *streamsize_array,
		       int nstreams, Mergeinfo_uint8_T mergeinfo);
#endif

#endif

