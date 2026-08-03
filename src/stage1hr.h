#ifndef STAGE1HR_INCLUDED
#define STAGE1HR_INCLUDED

typedef struct Stage1_T *Stage1_T;

#include "bool.h"
#include "univcoord.h"
#include "types.h"
#include "mode.h"
#include "pass.h"

#include "reader.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "shortread.h"
#include "iit-read-univ.h"
#include "ef64.h"

#include "singleelt.h"

#include "auxinfo.h"
#include "pathstore.h"

#include "genome.h"
#include "transcriptome.h"
#include "genomebits.h"

#include "localdb-read.h"
#include "mergeinfo.h"
#include "splice.h"
#include "indel.h"

#include "knownsplicing.h"
#include "knownindels.h"

#include "trdiagpool.h"
#include "univdiagpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "uintlistpool.h"
#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "trpathpool.h"
#include "hitlistpool.h"

#define STAGE1_MINPURE_SPECIFICITY_SCALE 1000000U


#define T Stage1_T
struct T {
  /* Results */
  bool first_read_p;
  int querylength;

  /* Trdiagonals */
  Trnum_T *sense_trnums;	/* Could have duplicates */
  Trcoord_T *sense_troffsets;
  Trcoord_T *sense_trhighs;
  Trcoord_T *_sense_trdiagonals; /* Aligned for !LARGE_GENOMES, AVX512, AVX2.  Could have multiple trdiagonals for the same trnum */
  int *sense_tstarts;
  int *sense_tends;
  int n_sense_trdiagonals;

  Trnum_T *antisense_trnums;	/* Could have duplicates */
  Trcoord_T *antisense_troffsets;
  Trcoord_T *antisense_trhighs;
  Trcoord_T *_antisense_trdiagonals; /* Aligned for !LARGE_GENOMES, AVX512, AVX2.  Could have multiple trdiagonals for the same trnum */
  int *antisense_tstarts;
  int *antisense_tends;
  int n_antisense_trdiagonals;

  /* Trpaths */
  List_T partial_sense_trpaths;	/* Used only to generate univdiagonals */
  List_T partial_antisense_trpaths;

  List_T sense_trpaths;
  List_T antisense_trpaths;


  /* Univdiagonals */
  bool exact_paths_computed_p;

  /* Result from Kmer_segment */
  Univcoord_T *exhaustive_gplus; /* Memory is aligned */
  int *exhaustive_counts_gplus;
  Singleelt_T *exhaustive_eltarray_gplus;
  unsigned int *exhaustive_minpure_specificity_gplus;
  int *exhaustive_minpure_observations_gplus;
  unsigned int *exhaustive_minpure_specificity_min_gplus;
  unsigned int *exhaustive_minpure_specificity_max_gplus;
  unsigned int *exhaustive_minpure_frequency_min_gplus;
  unsigned int *exhaustive_minpure_frequency_max_gplus;
  unsigned long long *exhaustive_minpure_frequency_sum_gplus;
  Auxinfo_T *exhaustive_auxinfo_gplus;
  int nexhaustive_gplus;

  Univcoord_T *exhaustive_gminus; /* Memory is aligned */
  int *exhaustive_counts_gminus;
  Singleelt_T *exhaustive_eltarray_gminus;
  unsigned int *exhaustive_minpure_specificity_gminus;
  int *exhaustive_minpure_observations_gminus;
  unsigned int *exhaustive_minpure_specificity_min_gminus;
  unsigned int *exhaustive_minpure_specificity_max_gminus;
  unsigned int *exhaustive_minpure_frequency_min_gminus;
  unsigned int *exhaustive_minpure_frequency_max_gminus;
  unsigned long long *exhaustive_minpure_frequency_sum_gminus;
  Auxinfo_T *exhaustive_auxinfo_gminus;
  int nexhaustive_gminus;

  /* From all search methods */
  Univcoord_T *all_univdiagonals_gplus; /* Memory is aligned */
  Univcoord_T *all_univdiagonals_gminus;
  Auxinfo_T *all_auxinfo_gplus;
  Auxinfo_T *all_auxinfo_gminus;
  Pathstore_T *all_pathstores_gplus;
  Pathstore_T *all_pathstores_gminus;
  int all_nunivdiagonals_gplus;
  int all_nunivdiagonals_gminus;

#if 0
  /* Univdiagonal links */
  Univcoordtable_T stored_paths_gplus;	/* For storing paths of each univdiagonal */
  Univcoordtable_T stored_paths_gminus;
  Univcoordtableuint_T stored_methods_gplus;	 /* For storing methods of each univdiagonal */
  Univcoordtableuint_T stored_methods_gminus;
#endif


  /* Global list of paths found so far */
  Path_T *sense_paths_gplus;
  Path_T *sense_paths_gminus;
  Path_T *antisense_paths_gplus;
  Path_T *antisense_paths_gminus;

#if 0
  /* Replaced by all_univdiagonals/all_auxinfo */
  /* These are the univdiagonals corresponding to the paths found so
     far and sent to Concordance_gen */
  /* Need to be aligned memory for merging */
  Univcoord_T *sense_coords_gplus;
  Univcoord_T *sense_coords_gminus;
  Univcoord_T *antisense_coords_gplus;
  Univcoord_T *antisense_coords_gminus;

  int *sense_indices_gplus;
  int *sense_indices_gminus;
  int *antisense_indices_gplus;
  int *antisense_indices_gminus;

  int n_sense_paths_gplus;
  int n_sense_paths_gminus;
  int n_antisense_paths_gplus;
  int n_antisense_paths_gminus;

  int nunique_sense_coords_gplus;
  int nunique_sense_coords_gminus;
  int nunique_antisense_coords_gplus;
  int nunique_antisense_coords_gminus;
#endif


  Reader_T reader;
  Reader_T tr_reader;

  /* Initialized by Stage1_fill_all_positions */
  bool *validp;			/* Need only one, since we allocate
				   forward_oligos[querypos] and
				   revcomp_oligos[querypos] at same
				   time */
  bool *segment_positions_filterp;
  int *segment_positions_filter_plus_npositions;
  int *segment_positions_filter_minus_npositions;
  bool *tr_validp;

  bool all_oligos_gen_filledp;
  bool all_positions_gen_filledp;

  Oligospace_T *forward_oligos;
  Oligospace_T *revcomp_oligos;
  Oligospace_T *tr_forward_oligos;
  Oligospace_T *tr_revcomp_oligos;

  Signature_T *forward_signatures;
  Signature_T *revcomp_signatures;

  /* Need plus_retrievedp and minus_retrievedp because
     Extension_search and Tr_extension_search can retrieve either plus
     or minus */
  bool *retrievedp_allocated;
  bool *plus_retrievedp;	       /* points to above[index1interval-1] */
  bool *minus_retrievedp;	       /* points to above[index1interval-1] */

  bool *tr_retrievedp_allocated;
  bool *tr_plus_retrievedp;
  bool *tr_minus_retrievedp;

#ifdef LARGE_GENOMES
  unsigned char **positions_high_allocated;
  unsigned char **plus_positions_high; /* points to above[index1interval-1] */
  unsigned char **minus_positions_high; /* points to above[index1interval-1] */
#endif

  /* plus positions are the alignments to the plus genome strand of
     the forward oligo, and also the alignments to the minus genome
     strand of the revcomp oligo */
  UINT4 **positions_allocated;
  UINT4 **plus_positions; /* points to above[index1interval-1] */
  UINT4 **minus_positions; /* points to above[index1interval-1] */
  int *plus_diagterms;
  int *minus_diagterms;

  UINT4 **tr_positions_allocated;
  UINT4 **tr_plus_positions;
  UINT4 **tr_minus_positions;

  int *npositions_allocated;
  int *plus_npositions;		/* points to above[index1interval-1] */
  int *minus_npositions;	/* points to above[index1interval-1] */

  int *tr_npositions_allocated;
  int *tr_plus_npositions;
  int *tr_minus_npositions;
  int *tr_plus_diagterms;
  int *tr_minus_diagterms;

  int segment_positions_filter_removed;
  
  int *mismatch_positions_alloc; /* Non-pool, based on querylength */
  int *positions_alloc;		 /* Non-pool, based on querylength */

  /* Memory allocated for indelinfo */
  Indelinfo_T indelinfo;	/* Non-pool, based on querylength */

  /* Memory allocated for spliceinfo in kmer-search.c and
     path-solve.c, used by Splice_resolve_sense and
     Splice_resolve_antisense */
  Spliceinfo_T spliceinfo;	/* Non-pool, based on querylength */

  Mergeinfo_uint4_T mergeinfo_tr; /* Non-pool, based on querylength */
#ifdef LARGE_GENOMES
  Mergeinfo_uint8_T mergeinfo;	/* Non-pool, based on querylength */
#else
  Mergeinfo_uint4_T mergeinfo;	/* Non-pool, based on querylength */
#endif

  /* Memory allocated for Segment_identify in segment-search.c, and
     Merge_diagonals in kmer-search.c (which needs four sets of
     arrays) */
#ifdef LARGE_GENOMES
  /* Non-pool, based on querylength */
  unsigned char **stream_high_alloc, **gplus_stream_high_array_5, **gminus_stream_high_array_5, **gplus_stream_high_array_3, **gminus_stream_high_array_3;
  UINT4 **stream_low_alloc, **gplus_stream_low_array_5, **gminus_stream_low_array_5, **gplus_stream_low_array_3, **gminus_stream_low_array_3;
#endif

  /* int streamspace_max_alloc;	-- Entries in streamspace_alloc; now allocated by worker thread */
  /* Univcoord_T *streamspace_alloc; -- Depends on max_localdb_regions; now allocated by worker thread */

  int max_nstreams;		/* Length of streamptr_alloc */


  /* Non-pool, based on querylength */
  Univcoord_T **streamptr_alloc, **gplus_stream_array_5, **gminus_stream_array_5, **gplus_stream_array_3, **gminus_stream_array_3;
  Trcoord_T **tplus_stream_array, **tminus_stream_array;
  int *querypos_diagterm_alloc, *tplus_diagterm_array, *tminus_diagterm_array; /* Non-pool, based on querylength */

  int *streamsize_alloc, *tplus_streamsize_array, *tminus_streamsize_array;
#if 0
  /* Used for Kmer_exact2, not being used any more */
  
  int *gplus_streamsize_array_5, *gminus_streamsize_array_5, *gplus_streamsize_array_3, *gminus_streamsize_array_3;
  int *gplus_diagterm_array_5, *gminus_diagterm_array_5, *gplus_diagterm_array_3, *gminus_diagterm_array_3;
#endif


  /* Returned from Tr_extension_search */
  List_T tr_queryfwd_plus_set;
  List_T tr_queryfwd_minus_set;
  List_T tr_queryrev_plus_set;
  List_T tr_queryrev_minus_set;

#if 0
  /* Intermediate calculations for transcriptome */
  /* Returned from Transcriptome_search_ends and used for Transcriptome_search_complete */
  Trcoord_T *tplus_positions_5, *tminus_positions_5, *tplus_positions_3, *tminus_positions_3;
  int n_tplus_positions_5, n_tminus_positions_5, n_tplus_positions_3, n_tminus_positions_3;
  int tplus_diagterm_5, tminus_diagterm_5, tplus_diagterm_3, tminus_diagterm_3;
#endif
};


extern void
Stage1_list_trpaths (T this);

extern void
Stage1_list_extension (T this);

extern void
Stage1_list_exhaustive (T this);

extern void
Stage1_list_all_univdiagonals (T this);

extern void
Stage1_trdiagonals_gc (T this);

extern void
Stage1_free (T *old, Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
	     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool, Intlistpool_T intlistpool,
	     Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Trpathpool_T trpathpool,
	     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern T
Stage1_new (bool first_read_p);

extern void
Stage1_set_sequence (T this, Shortread_T queryseq, Localdb_T localdb);

extern void
Stage1_collect_unextended_paths (List_T *unextended_sense_paths_gplus,
				 List_T *unextended_sense_paths_gminus,
				 List_T *unextended_antisense_paths_gplus,
				 List_T *unextended_antisense_paths_gminus,
				 T this, Hitlistpool_T hitlistpool);

extern bool
Stage1_collect_paths (List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,
		      T this, Hitlistpool_T hitlistpool);

extern void
Stage1_dump (T this, int querylength);
extern void
Stage1_dump_tr (T this, int querylength);

extern void
Stage1_init_end_gen (int *querystart, int *queryend, T this,
		     int querylength, int genestrand);
extern void
Stage1_init_end_tr (T this, int querylength);

extern void
Stage1_fill_all_oligos_gen (T this, int querylength, int genestrand);
extern void
Stage1_fill_all_oligos_tr (T this, int querylength);

extern void
Stage1_fill_all_positions_gen (int *total_npositions_plus, int *total_npositions_minus,
			       T this, int querylength, int genestrand);
extern int
Stage1_filter_segment_kmer_positions_gen (T this, int querylength,
					  int max_avg_kmer_positions, const char *label);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
extern void
Stage1_print_segment_kmer_positions_gen (T this, int querylength, const char *label);
#endif
extern void
Stage1_fill_all_positions_tr (T this, int querylength);

extern void
Stage1hr_setup (Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in, Indexdb_T indexdb_tr_in,
		EF64_T repetitive_ef64_in, bool segment_repetitive_kmer_filter_p_in,
		int index1part_in, int index1interval_in, int msize_in,
		int index1part_tr_in, int index1interval_tr_in, 
		int max_deletionlen, Chrpos_T shortsplicedist,
		Transcriptome_T transcriptome_in);

#undef T
#endif
