#ifndef STAGE1HR_PAIRED_INCLUDED
#define STAGE1HR_PAIRED_INCLUDED

#include <stdint.h>
#include "bool.h"
#include "ef64.h"
#include "align-cascade.h"
#include "cascade-timing.h"

#include "path.h"
#include "pathpair.h"

#include "types.h"
#include "shortread.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "intlistpool.h"
#include "univcoord.h"
#include "univdiagpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "hitlistpool.h"
#include "transcriptpool.h"
#include "spliceendsgen.h"

#include "pass.h"

#define STAGE1HR_PAIRED_DEFAULT_SEGMENT_EXACT_MAX_MISMATCH_DIVISOR 20

#define T Stage1_T

extern Pathpair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Path_T **patharray5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Path_T **patharray3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    T this5, T this3, Shortread_T queryseq5, Shortread_T queryseq3, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Chrpos_T pairmax_linear,

		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		    int localdb_nmismatches_allowed,

		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		    Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3,
		    Cascade_Timing_T *cascade_timing, Pass_T pass);

extern void
Stage1hr_paired_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       Genomebits_T genomebits_in, EF64_T chromosome_ef64_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       int max_deletionlen, int max_insertlength, Chrpos_T shortsplicedist_in, bool splicingp_in,
		       int maxpaths_search_in, int maxpaths_report_in,
		       bool *circularp_in, int pairmax_linear_in, int pairmax_circular_in,
		       int paired_segment_exact_max_mismatch_divisor_in,
		       bool paired_mate_rescue_gapped_p_in,
		       bool paired_preeval_ungapped_mate_p_in,
		       bool paired_preeval_rescue_skip_p_in,
		       bool paired_segment_stream_chains_p_in,
		       bool paired_solve_all_anchor_first_p_in,
		       bool paired_postsolve_reliable_anchors_first_p_in,
		       bool paired_seg1_exhaustive_concordance_p_in,
		       bool paired_duplicate_auxinfo_frontier_p_in,
		       unsigned int paired_mate_minpure_floor_in,
		       bool paired_transcript_completion_p_in,
		       bool paired_tr_adjacent_indel_repair_p_in,
		       bool tr_chain_search_p_in,
		       double se_segment_min_kmer_avg_quality_in,
		       int se_segment_max_kmer_bad_bases_in,
		       int se_segment_bad_base_quality_in,
		       int segment_max_avg_kmer_positions_in,
		       Align_cascade_T align_cascade_in);

extern void
Stage1hr_paired_pass2_setup (int max_insertlength, Chrpos_T shortsplicedist);

extern void
Stage1hr_paired_localdb_setup (Localdb_T localdb_in);


#undef T
#endif
