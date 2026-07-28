#ifndef PATHPAIR_EVAL_INCLUDED
#define PATHPAIR_EVAL_INCLUDED

#include "path.h"
#include "pathpair.h"

#include "bool.h"
#include "univcoord.h"
#include "univdiag.h"
#include "shortread.h"
#include "compress.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "univdiagpool.h"
#include "hitlistpool.h"
#include "transcriptpool.h"

#include "stage1hr.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "genomebits.h"
#include "transcriptome.h"


#define T Path_T

extern List_T
Pathpair_uniq (List_T pathpairs, Pathpair_T *pathpairarray, int npaths,
	       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool);

extern Pathpair_T *
Pathpair_eval_and_sort (int *found_score_5, int *found_score_3,
			int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
			Pathpair_T *pathpairarray, int npaths, Stage1_T stage1_5, Stage1_T stage1_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Shortread_T queryseq5, Shortread_T queryseq3,
			char *queryuc_ptr_5, char *queryrc5, char *queryuc_ptr_3, char *queryrc3, 
			char *quality_string_5, char *quality_string_3,
			int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

			int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, int localdb_nmismatches_allowed,

			Knownsplicing_T knownsplicing, Knownindels_T knownindels,

			int nmismatches_filter_5, int nmismatches_filter_3,
			int mincoverage_filter_5, int mincoverage_filter_3,
			int querylength5, int querylength3,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
			Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, bool filterp);

extern void
Pathpair_eval_setup (int expected_pairlength, int pairlength_deviation,
		     Transcriptome_T transcriptome_in,
		     bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		     Outputtype_T output_type_in,
		     bool splicingp_in, bool resolve_inner_p_in, bool want_random_p_in,
		     bool allow_soft_clips_p_in);

#undef T
#endif


