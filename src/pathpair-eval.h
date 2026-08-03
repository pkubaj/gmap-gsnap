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

#define PATHPAIR_DEFAULT_TERMINAL_ANCHOR_RESCUE_MIN_SCORE 15

#define T Path_T

extern bool
Pathpair_near_exact_terminal_anchor_partner_p (Path_T path);

extern double
Path_final_score (Path_T path);

extern bool
Path_final_score_reportable_p (Path_T path);

extern void
Path_set_final_score (Path_T path, int pair_score);

extern List_T
Pathpair_uniq (List_T pathpairs, Pathpair_T *pathpairarray, int npaths,
	       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool);

extern Pathpair_T *
Pathpair_eval_and_sort (int *found_score_5, int *found_score_3, int *end_score_5, int *end_score_3,
			int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
			List_T *singlepaths5_from_rejected, List_T *singlepaths3_from_rejected,
			Pathpair_T *pathpairarray, int npaths, Stage1_T stage1_5, Stage1_T stage1_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Shortread_T queryseq5, Shortread_T queryseq3,
			char *queryuc_ptr_5, char *queryrc5, char *queryuc_ptr_3, char *queryrc3, 
			char *quality_string_5, char *quality_string_3,
			int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			int max_nblocks, int localdb_nmismatches_allowed,

			Knownsplicing_T knownsplicing, Knownindels_T knownindels,

			int nmismatches_filter_5, int nmismatches_filter_3,
			int mincoverage_filter_5, int mincoverage_filter_3,
			int querylength5, int querylength3,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
			Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			List_T *terminal_anchor_pathpairs_deferred,
			List_T *terminal_anchor_paths5_deferred,
			List_T *terminal_anchor_paths3_deferred,
			bool filterp, bool allow_postfallback_terminal_anchor_rescue_p,
			bool allow_near_exact_terminal_anchor_rescue_p);

extern void
Pathpair_eval_setup (int expected_pairlength, int pairlength_deviation,
		     Transcriptome_T transcriptome_in,
		     bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		     Outputtype_T output_type_in,
		     bool splicingp_in, bool resolve_inner_p_in, bool want_random_p_in,
		     bool allow_soft_clips_start_p_in, bool allow_soft_clips_end_p_in,
		     int paired_terminal_anchor_min_score_in,
		     int paired_terminal_anchor_slop_in,
		     bool mincoverage_filter_per_read_p_in,
		     bool pathpair_equal_coverage_diff_tiebreak_p_in,
		     double score_match_reward_in,
		     double score_mismatch_penalty_in,
		     double score_indel_open_penalty_in,
		     double score_indel_extend_penalty_in,
		     double score_splice_unannotated_penalty_in,
		     double score_pair_transcript_discordance_penalty_in,
		     double score_pair_margin_in,
		     double score_min_read_in,
		     double score_min_pair_in);

#undef T
#endif
