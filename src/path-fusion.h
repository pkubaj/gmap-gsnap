/* $Id: a6fe395024af59a332f63d732509155add4924ae $ */
#ifndef PATH_FUSION_INCLUDED
#define PATH_FUSION_INCLUDED

#include "path.h"

#include "bool.h"
#include "univcoord.h"
#include "univdiag.h"
#include "compress.h"
#include "shortread.h"

#include "genomebits.h"
#include "stage1.h"
#include "knownsplicing.h"

#include "ef64.h"
#include "transcriptome.h"
#include "genomebits.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "hitlistpool.h"
#include "univdiagpool.h"
#include "transcriptpool.h"
#include "pass.h"


#define T Path_T

extern List_T
Path_fusion_outer_querystart_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1, 

				   Compress_T query_compress_fwd, Compress_T query_compress_rev,
				   Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				   int genestrand, int nmismatches_allowed,
				   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				   Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern List_T
Path_fusion_outer_queryend_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				 Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				 int genestrand, int nmismatches_allowed,
				 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				 Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern List_T
Path_fusion_outer_querystart_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				    Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    Shortread_T queryseq,  int querylength, Knownsplicing_T knownsplicing,
				    int genestrand, int nmismatches_allowed,
				    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				    Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);


extern List_T
Path_fusion_outer_queryend_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				  int genestrand, int nmismatches_allowed,
				  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				  Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern T
Path_fusion_inner_qend (int *found_score, int *end_score, T fusion5, T anchor3,
			char *queryptr_main, bool main_plusp, int querylength,
			Chrnum_T main_chrnum, Univcoord_T main_chroffset, Univcoord_T main_chrhigh,
			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			Stage1_T stage1, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
			Compress_T query_compress_main, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Shortread_T queryseq, int genestrand, int nmismatches_allowed,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern T
Path_fusion_inner_qstart (int *found_score, int *end_score, T fusion3, T anchor5,
			  char *queryptr_main, bool main_plusp, int querylength,
			  Chrnum_T main_chrnum, Univcoord_T main_chroffset, Univcoord_T main_chrhigh,
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			  Stage1_T stage1, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
			  Compress_T query_compress_main, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			  Shortread_T queryseq, int genestrand, int nmismatches_allowed,
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);


extern List_T
Path_fusion_querystart_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			     int genestrand, int nmismatches_allowed,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);
extern List_T
Path_fusion_querystart_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			      int genestrand, int nmismatches_allowed,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern List_T
Path_fusion_queryend_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			   int genestrand, int nmismatches_allowed,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern List_T
Path_fusion_queryend_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			    Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			    int genestrand,  int nmismatches_allowed,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Path_fusion_setup (bool *circularp_in, EF64_T chromosome_ef64_in, Univcoord_T genomelength_in,
		   Chrpos_T shortsplicedist_in, Transcriptome_T transcriptome_in,
		   Genomebits_T genomebits_in);

#undef T
#endif
