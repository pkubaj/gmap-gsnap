/* $Id: c7ef43f3b7c2431bf22df66e57b009641660ed05 $ */
#ifndef STAGE1HR_SINGLE_INCLUDED
#define STAGE1HR_SINGLE_INCLUDED

#include "stage1hr-single.h"

#include "bool.h"
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
single_read_next_method_tr (int *found_score, int sufficient_score, Method_T last_method,

			    List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
			    List_T *sense_trpaths, List_T *antisense_trpaths,

			    T this, int genestrand, int querylength,
			    int *mismatch_positions_alloc,
			    Compress_T query_compress_fwd, Compress_T query_compress_rev,

			    int nmismatches_allowed,
			 
			    Trdiagpool_T trdiagpool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			    Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			    Hitlistpool_T hitlistpool, bool paired_end_p, bool first_read_p,
			    bool appendp);

extern void
single_read_unsolved_tr_paths (int *found_score,

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
		    Shortread_T queryseq, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Localdb_T localdb,
		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
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
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       bool splicingp_in);

#undef T
#endif

