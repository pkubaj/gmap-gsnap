#ifndef PATH_EVAL_INCLUDED
#define PATH_EVAL_INCLUDED

#include "path.h"
#include "pathpair.h"

#include "bool.h"
#include "univcoord.h"
#include "compress.h"

#include "intlistpool.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"

#include "stage1hr.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "genomebits.h"
#include "outputtype.h"


#define T Path_T

extern int
Path_eval_nmatches (int *found_score, int *end_score, T this,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev);

extern int
Path_eval_nmatches_chopped (T this);

extern int
Path_nmatches_cmp (T a, T b);

extern void
Path_mark_alignment (T path, Compress_T query_compress_fwd, char *queryuc_ptr,
		     Compress_T query_compress_rev, char *queryrc, Pathpool_T pathpool);

extern bool
Path_eval_perfect_ends_p (T this, Compress_T query_compress_fwd, char *queryuc_ptr,
			  Compress_T query_compress_rev, char *queryrc,
			  int querystart, int queryend, Pathpool_T pathpool);

extern int
Path_best_cmp (const void *x, const void *y);

extern int
Path_best_pair_cmp (T a, T b);

extern T *
Path_eval_and_sort (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T *patharray, int npaths,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    char *queryuc_ptr, char *queryrc, char *quality_string,
		    int nmismatches_filter, int mincoverage_filter, int max_softclip_filter,
		    bool split_paths_p, int split_path_max_overlap,
		    int split_path_min_segment, int split_path_min_combined_coverage,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Hitlistpool_T hitlistpool, bool filterp);

extern List_T
Path_consolidate (List_T paths, Shortread_T queryseq,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		  Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool,
		  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Path_eval_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		 Transcriptome_T transcriptome_in,
		 bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		 int index1part_in, int index1interval_in,
		 Outputtype_T output_type_in, bool md_report_snps_p_in,
		 bool want_random_p_in, bool allow_soft_clips_start_p_in,
		 bool allow_soft_clips_end_p_in,
		 bool quality_ranker_p_in, double quality_match_reward_in,
		 double quality_mismatch_penalty_in,
		 double quality_indel_penalty_in,
		 double quality_softclip_penalty_in,
		 double quality_error_rate_penalty_in,
		 double quality_rank_margin_in,
		 bool quality_split_ranker_p_in,
		 double quality_split_gap_penalty_in,
		 double quality_split_segment_penalty_in);

#undef T
#endif
