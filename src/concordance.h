/* $Id: c2d8248301f76c304e492b8326a95cece3910ac2 $ */
#ifndef CONCORDANCE_INCLUDED
#define CONCORDANCE_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "list.h"
#include "types.h"
#include "genomicpos.h"
#include "pass.h"

#include "shortread.h"
#include "compress.h"
#include "knownsplicing.h"
#include "knownindels.h"
#include "stage1hr.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "hitlistpool.h"
#include "pathpool.h"
#include "transcriptpool.h"
#include "vectorpool.h"


extern List_T
Concordance_tr (int *found_score_paired, int *found_score_5, int *found_score_3, List_T pathpairs,

		List_T newtrpaths5, List_T newtrpaths3,
		List_T trpaths5, List_T trpaths3,
		    
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

		Shortread_T queryseq5, Shortread_T queryseq3,
		int querylength5, int querylength3,
		Stage1_T stage1_5, Stage1_T stage1_3, Knownsplicing_T knownsplicing,

		int nmismatches_filter_5, int nmismatches_filter_3,
		int mincoverage_filter_5, int mincoverage_filter_3,

		Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, int sensedir);

extern void
Concordance_setup (int subopt_levels_in, Chrpos_T pairmax_transcriptome_in,
		   Chrpos_T pairmax_linear_in, Chrpos_T pairmax_circular_in,
		   bool *circularp_in, bool merge_samechr_p_in, bool two_pass_p);

#undef T
#endif

