#ifndef PATH_SOLVE_INCLUDED
#define PATH_SOLVE_INCLUDED

#include "types.h"
#include "genomicpos.h"
#include "univcoord.h"
#include "sense.h"
#include "method.h"
#include "pass.h"

#include "bool.h"
#include "list.h"
#include "univdiag.h"
#include "chrnum.h"
#include "compress.h"
#include "shortread.h"

#include "genomebits.h"
#include "localdb-read.h"
#include "indel.h"
#include "splice.h"
#include "stage1hr.h"
#include "knownsplicing.h"
#include "knownindels.h"


#ifdef LARGE_GENOMES
#include "uint8listpool.h"
#else
#include "uintlistpool.h"
#endif


#include "intlistpool.h"
#include "uintlistpool.h"
#include "univdiagpool.h"
#include "listpool.h"
#include "hitlistpool.h"
#include "pathpool.h"
#include "univdiagpool.h"
#include "vectorpool.h"
#include "spliceendsgen.h"

#include "auxinfo.h"
#include "path.h"


#define T Path_T
typedef struct T *T;

typedef struct Path_segment_cache_T *Path_segment_cache_T;
typedef struct Path_ext_anchor_cache_T *Path_ext_anchor_cache_T;

extern Path_ext_anchor_cache_T
Path_ext_anchor_cache_new (void);

extern void
Path_ext_anchor_cache_free (Path_ext_anchor_cache_T *old,
			    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool,
			    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern Path_segment_cache_T
Path_segment_cache_new (bool emitted_path_dedup_p);

extern void
Path_segment_cache_free (Path_segment_cache_T *old,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool,
			 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern bool
Path_segment_cache_qend_known_p (Path_segment_cache_T cache,
				 int head_qstart, int head_qend, int head_nmismatches,
				 Univcoord_T head_univdiagonal,
				 int next_qstart, int next_qend, int next_nmismatches,
				 Univcoord_T next_univdiagonal,
				 int querylength, int localdb_nmismatches_allowed,
				 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				 bool plusp, int genestrand, bool first_read_p,
				 bool find_splices_p);

extern bool
Path_segment_cache_qstart_known_p (Path_segment_cache_T cache,
				   int head_qstart, int head_qend, int head_nmismatches,
				   Univcoord_T head_univdiagonal,
				   int next_qstart, int next_qend, int next_nmismatches,
				   Univcoord_T next_univdiagonal,
				   int querylength, int localdb_nmismatches_allowed,
				   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				   bool plusp, int genestrand, bool first_read_p,
				   bool find_splices_p);


/* Used by kmer-search, extension-search, and segment-search */
/* Unsolved paths are those that could not be found using indel
   resolve or splice resolve functions.  They potentially represent
   splice plus indels, which require a more complicated solution,
   deferred until later */
extern void
Path_solve_from_diagonals (int *found_score, int *end_score,

			   List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			   List_T *sense_paths, List_T *antisense_paths,

			   Univcoord_T middle_univdiagonal, Auxinfo_T auxinfo,
			   
			   Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			   int localdb_nmismatches_allowed,

			   Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand,
			   bool paired_end_p, bool first_read_p,

			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			   Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern bool
Path_solve_exact (int *found_score, int *end_score,

		  List_T *sense_paths, List_T *antisense_paths,

		  Univcoord_T univdiagonal, Auxinfo_T auxinfo, int querylength,
		  int max_nmismatches,
		  bool plusp, bool first_read_p, int genestrand,
		  Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Shortread_T queryseq, char *queryuc_ptr, char *queryrc,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool, Method_T method);

extern int
Path_solve_count_mismatches_substring (int *ref_nmismatches, Compress_T query_compress,
				       Univcoord_T univdiagonal, int querylength,
				       int pos5, int pos3, bool plusp, int genestrand);

extern void
Path_solve_from_univdiagonal (int *found_score, int *end_score,

			      List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			      List_T *sense_paths, List_T *antisense_paths,

			      Univcoord_T univdiagonal, Auxinfo_T auxinfo,
			      Shortread_T queryseq, char *queryptr, Compress_T query_compress,
			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      bool plusp, int querylength, int *mismatch_positions_alloc,
			      
			      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			      int localdb_nmismatches_allowed,

			      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      bool paired_end_p, bool first_read_p,

			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			      Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			      Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern void
Path_solve_from_right_chain (int *found_score, int *end_score,

			     List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			     List_T *sense_paths, List_T *antisense_paths,

			     List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			     int localdb_nmismatches_allowed,

			     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			     Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern void
Path_solve_from_right_chain_cached (int *found_score, int *end_score,

				    List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
				    List_T *sense_paths, List_T *antisense_paths,

				    List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

				    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
				    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
				    int localdb_nmismatches_allowed,

				    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
				    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
				    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
				    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
				    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
				    Path_ext_anchor_cache_T ext_anchor_cache,
				    Path_segment_cache_T segment_cache,
				    bool segment_prefix_skip_p,
				    bool segment_qstart_anchor_cache_p,
				    bool path_qend_local_cache_p,
				    bool emitted_path_dedup_p);

extern void
Path_solve_from_left_chain (int *found_score, int *end_score,

			    List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			    List_T *sense_paths, List_T *antisense_paths,

			    List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			    int localdb_nmismatches_allowed,

			    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern void
Path_solve_from_left_chain_cached (int *found_score, int *end_score,

				   List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
				   List_T *sense_paths, List_T *antisense_paths,

				   List_T left_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

				   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
				   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
				   int localdb_nmismatches_allowed,

				   Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
				   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				   bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
				   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
				   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
				   Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
				   Path_ext_anchor_cache_T ext_anchor_cache,
				   Path_segment_cache_T segment_cache,
				   bool segment_prefix_skip_p,
				   bool segment_qstart_anchor_cache_p,
				   bool path_qend_local_cache_p,
				   bool emitted_path_dedup_p);

extern List_T
Path_extend (int *found_score, int *end_score, List_T *global_unextended_paths,
	     T original_path, Shortread_T queryseq, char *queryptr, int querylength,
	     int *mismatch_positions_alloc,
	     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
	     int localdb_nmismatches_allowed,

	     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
	     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     int genestrand, bool paired_end_p, bool lowp,
	     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
	     Spliceendsgen_T spliceendsgen, bool extend_qstart_p, bool extend_qend_p);

extern void
Path_qstart_resolve (int *found_score, int *end_score, T path,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal, char *queryptr, int querylength,
		     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		     int localdb_nmismatches_allowed,

		     Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev, int genestrand, 
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool);

extern void
Path_qend_resolve (int *found_score, int *end_score, T path,
		   Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal, char *queryptr, int querylength,
		   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		   int localdb_nmismatches_allowed,

		   Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev, int genestrand,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool);

extern T
Path_solve_junctions (int *found_score, int *end_score, T this, int sensedir, int genestrand,
		      Compress_T query_compress,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      Shortread_T queryseq, char *queryptr, int querylength,

		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Localdb_T localdb, int localdb_nmismatches_allowed,

		      Stage1_T stage1, Knownsplicing_T knownsplicing,
		      Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool);

extern void
Path_solve_setup (bool *circularp_in, Transcriptome_T transcriptome_in, EF64_T chromosome_ef64_in,
		  Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in, Univcoord_T genomelength_in,
		  int index1part_in, int min_intronlength_in,
		  int max_insertionlen_in, int max_deletionlen_in, int shortsplicedist_in,
		  bool novelsplicingp, bool knownsplicingp, int n_splice_paths_in,
		  int max_recursive_ends_in, int min_support_indel_in,
		  int indel_support_mismatch_penalty_in,
		  bool materialize_equiv_splice_alts_p_in);

extern void
Path_solve_localdb_setup (Localdb_T localdb_in);

#undef T
#endif
