/* $Id: 1d78ca259d112cf935cace1976bbcd4325f14352 $ */
#ifndef TRPATH_CONVERT_INCLUDED
#define TRPATH_CONVERT_INCLUDED

#include "bool.h"
#include "transcriptome.h"
#include "ef64.h"

#include "shortread.h"
#include "stage1hr.h"
#include "knownsplicing.h"
#include "compress.h"

#include "list.h"
#include "intlistpool.h"
#include "uintlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"


extern void
Trpath_convert_sense (int *found_score,
		      List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		      List_T sense_trpaths, bool first_read_p,
		      Shortread_T queryseq, int querylength,
		      Stage1_T this, Knownsplicing_T knownsplicing,

		      Compress_T query_compress_fwd, Compress_T query_compress_rev,

		      int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, Localdb_T localdb, int localdb_nmismatches_allowed,

		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		      Hitlistpool_T hitlistpool, bool crossover_sense_p);

extern void
Trpath_convert_antisense (int *found_score,
			  List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,
			  List_T *sense_paths_gplus, List_T *sense_paths_gminus,

			  List_T antisense_trpaths, bool first_read_p,
			  Shortread_T queryseq, int querylength,
			  Stage1_T this, Knownsplicing_T knownsplicing,

			  Compress_T query_compress_fwd, Compress_T query_compress_rev,
			  
			  int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, Localdb_T localdb, int localdb_nmismatches_allowed,

			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			  Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			  Hitlistpool_T hitlistpool, bool crossover_sense_p);

extern void
Trpath_convert_setup (Transcriptome_T transcriptome_in,
		      EF64_T chromosome_ef64_in);


#endif
