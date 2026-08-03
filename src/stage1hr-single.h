#ifndef STAGE1HR_SINGLE_INCLUDED
#define STAGE1HR_SINGLE_INCLUDED

#include <stdint.h>
#include "bool.h"
#include "method.h"
#include "align-cascade.h"
#include "transcriptome.h"

#include "ef64.h"
#include "auxinfo.h"
#include "path.h"

#include "types.h"
#include "compress.h"
#include "shortread.h"
#include "knownsplicing.h"
#include "knownindels.h"
#include "localdb-read.h"

#include "intlistpool.h"
#include "univcoord.h"
#include "trdiagpool.h"
#include "segmenteltpool.h"
#include "univdiagpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "hitlistpool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "transcriptpool.h"
#include "spliceendsgen.h"

#include "pass.h"

#define T Stage1_T

extern Method_T
single_read_next_method_trdiagonal (Method_T last_method, T this, int querylength,
				    Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    bool first_read_p);

extern Method_T
single_read_trdiagonal_method_goal (void);

extern bool
single_read_trdiagonal_method_available_p (Method_T last_method, Method_T method_goal);

extern Method_T
single_read_next_method_tr (int *found_score, int *end_score, int sufficient_score, Method_T last_method,

			    List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
			    List_T *sense_trpaths, List_T *antisense_trpaths,

			    T this, int genestrand, int querylength,
			    int *mismatch_positions_alloc,
			    char *queryuc_ptr, char *queryrc,
			    Compress_T query_compress_fwd, Compress_T query_compress_rev,

			    int read_nmismatches_allowed,
			 
			    Trdiagpool_T trdiagpool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			    Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			    Hitlistpool_T hitlistpool, bool paired_end_p, bool first_read_p,
			    bool appendp);

extern Method_T
single_read_tr_method_goal (void);

extern bool
single_read_tr_method_available_p (Method_T last_method, Method_T method_goal);

extern void
single_read_unsolved_tr_paths (int *found_score, int *end_score, 

			       List_T *sense_paths_gplus, List_T *sense_paths_gminus, 
			       List_T *antisense_paths_gplus, List_T *antisense_paths_gminus, 

			       T this, Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,

			       int max_insertionlen, int max_deletionlen,
			 
			       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
			       Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			       Hitlistpool_T hitlistpool);

extern Path_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T this, Shortread_T queryseq, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Localdb_T localdb,

		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		    unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		    int max_nblocks, int localdb_nmismatches_allowed,

		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
		    Spliceendsgen_T spliceendsgen, bool single_cell_p, bool first_read_p,
		    Pass_T pass);


extern void
Stage1hr_single_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       bool refine_partial_tr_p_in, bool tr_chain_search_p_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       bool split_paths_p_in, int split_path_max_overlap_in,
		       int split_path_min_segment_in, double split_path_min_combined_coverage_float_in,
		       bool splicingp_in, bool se_ext_dedup_exact_chains_p_in,
		       bool se_ext_equiv_prev_prune_p_in,
		       bool se_keep_trpath_candidates_p_in,
		       bool se_segment_stream_chains_p_in,
		       bool se_segment_dedup_exact_chains_p_in,
		       bool se_segment_prefix_skip_p_in,
		       bool se_segment_qstart_anchor_cache_p_in,
		       bool se_path_qend_local_cache_p_in,
		       bool se_emitted_path_dedup_p_in,
		       bool se_segment_equiv_prev_prune_p_in,
		       int se_segment_cache_threshold_in,
		       int se_segment_max_window_in,
		       double se_segment_min_kmer_avg_quality_in,
		       int se_segment_max_kmer_bad_bases_in,
		       int se_segment_bad_base_quality_in,
		       int segment_max_avg_kmer_positions_in,
		       double se_max_mismatches_filter_float_in,
		       double se_max_softclip_filter_float_in,
		       Align_cascade_T align_cascade_in);

#undef T
#endif
