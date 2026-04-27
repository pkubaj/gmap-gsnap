/* $Id: deb9b58c9cd5445c68d12e911889616987c71129 $ */
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


/* Used by kmer-search, extension-search, and segment-search */
/* Unsolved paths are those that could not be found using indel
   resolve or splice resolve functions.  They potentially represent
   splice plus indels, which require a more complicated solution,
   deferred until later */
extern void
Path_solve_from_diagonals (int *found_score,

			   List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			   List_T *sense_paths, List_T *antisense_paths,

			   Univcoord_T middle_diagonal_univdiagonal, Auxinfo_T auxinfo,
			   
			   Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,
			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,

			   Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand, int localdb_nmismatches_allowed,
			   bool paired_end_p, bool first_read_p,

			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			   Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern bool
Path_solve_exact (int *found_score,

		  List_T *sense_paths, List_T *antisense_paths,

		  Univcoord_T univdiagonal, Auxinfo_T auxinfo, int querylength,
		  bool plusp, bool first_read_p, int genestrand,
		  Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Shortread_T queryseq, char *queryuc_ptr, char *queryrc,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool, Method_T method);

extern void
Path_solve_from_univdiagonal (int *found_score,

			      List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			      List_T *sense_paths, List_T *antisense_paths,

			      Univcoord_T univdiagonal, Auxinfo_T auxinfo,
			      Shortread_T queryseq, char *queryptr, Compress_T query_compress,
			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      bool plusp, int querylength, int *mismatch_positions_alloc,
			      
			      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      int localdb_nmismatches_allowed, bool paired_end_p, bool first_read_p,

			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			      Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			      Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p);

extern List_T
Path_extend (int *found_score, List_T *global_unextended_paths,
	     T original_path, Shortread_T queryseq, char *queryptr, int querylength,
	     int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
	     Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
	     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     int genestrand, int localdb_nmismatches_allowed, bool paired_end_p, bool lowp,
	     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
	     Spliceendsgen_T spliceendsgen, bool extend_qstart_p, bool extend_qend_p);

extern void
Path_qstart_resolve (int *found_score, T path,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		     char *queryptr, int querylength,
		     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		     Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     int genestrand, int localdb_nmismatches_allowed,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool);

extern void
Path_qend_resolve (int *found_score, T path,
		   Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		   char *queryptr, int querylength,
		   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		   Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		   int genestrand, int localdb_nmismatches_allowed,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool);

extern T
Path_solve_junctions (int *found_score, T this, int sensedir, int genestrand,
		      Compress_T query_compress,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      Shortread_T queryseq, int querylength,
		      Stage1_T stage1, Knownsplicing_T knownsplicing,
		      Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool);

extern void
Path_solve_setup (bool *circularp_in, Transcriptome_T transcriptome_in, EF64_T chromosome_ef64_in,
		  Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in, Univcoord_T genomelength_in,
		  int index1part_in, Localdb_T localdb_in, int min_intronlength_in,
		  int max_insertionlen_in, int max_deletionlen_in,
		  bool novelsplicingp, bool knownsplicingp);

#undef T
#endif

