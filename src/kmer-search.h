#ifndef KMER_SEARCH_INCLUDED
#define KMER_SEARCH_INCLUDED

typedef enum {
  SEGMENT_DIAGONAL_ALGORITHM_WIDTH_OR_COUNT,
  SEGMENT_DIAGONAL_ALGORITHM_ADAPTIVE_HYBRID
} Segment_diagonal_algorithm_T;

#include "method.h"
#include "list.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "compress.h"
#include "shortread.h"
#include "genomebits.h"
#include "indel.h"

#include "auxinfo.h"
#include "pathstore.h"

#include "bool.h"
#include "pass.h"
#include "univdiag.h"
#include "stage1hr.h"
#include "mergeinfo.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "univdiagpool.h"
#include "univcoord.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "vectorpool.h"
#include "listpool.h"
#include "hitlistpool.h"

/* Does not take paired_end_p as a parameter.  ? Generates both sense and antisense */
extern int
Kmer_exact1 (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	     Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	     Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	     Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
	     Stage1_T stage1, int querystart, int queryend, int querylength,
	     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool);

extern List_T
Kmer_solve_eltarray_left (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			  Univdiagpool_T univdiagpool, Listpool_T listpool);

extern List_T
Kmer_solve_eltarray_right (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			   Univdiagpool_T univdiagpool, Listpool_T listpool);

extern void
Kmer_apply_eltarray_left (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			  Listpool_T listpool,
			  Singleelt_chain_skip_T skip_chain_p,
			  bool equiv_prev_prune_p, int max_chain_window,
			  bool materialized_order_p,
			  Singleelt_chain_apply_T apply_chain, void *data);

extern void
Kmer_apply_eltarray_right (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			   Listpool_T listpool,
			   Singleelt_chain_skip_T skip_chain_p,
			   bool equiv_prev_prune_p, int max_chain_window,
			   bool materialized_order_p,
			   Singleelt_chain_apply_T apply_chain, void *data);

extern void
Kmer_segment (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	      Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	      Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	      Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,

	      Stage1_T stage1, int querylength, int max_nsplices, bool first_read_p,
	      const char *read_name, const char *label,

	      EF64_T repetitive_ef64, Univdiagpool_T univdiagpool,
	      Listpool_T listpool,
	      Univcoordlistpool_T univcoordlistpool,
	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool);

extern void
Kmer_search_setup (int index1part_in, int index1interval_in, EF64_T chromosome_ef64_in,
		   Univ_IIT_T chromosome_iit_in,
		   Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		   Univcoord_T genomelength_in,
		   int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist,
		   bool splicingp_in, Segment_diagonal_algorithm_T segment_diagonal_algorithm_in,
		   int adaptive_hybrid_floor_per_kmer_in,
		   int adaptive_hybrid_budget_per_kmer_in,
		   int adaptive_hybrid_hard_max_in,
		   int adaptive_hybrid_weak_width_per_kmer_in,
		   int adaptive_hybrid_weak_support_divisor_in,
		   double adaptive_hybrid_strong_support_fraction_in,
		   double stats_weak_support_fraction_in,
		   bool segment_repetitive_kmer_filter_p_in,
		   int segment_max_kmer_positions_in,
		   bool minpure_diagonal_debug_p_in);

#include "minpure.h"
extern void
Kmer_search_setup_minpure (Minpure_Table_T genome_minpure_in);

#endif
