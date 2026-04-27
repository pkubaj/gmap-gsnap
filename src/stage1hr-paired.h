/* $Id: 318372917628cbe9a2673f83142dc565567f7cf6 $ */
#ifndef STAGE1HR_PAIRED_INCLUDED
#define STAGE1HR_PAIRED_INCLUDED

#include "stage1hr-paired.h"

#include "bool.h"
#include "ef64.h"

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

extern Pathpair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Path_T **patharray5, int *nhits5_primary, int *nhits5_altloc, int *first_absmq5, int *second_absmq5,
		    Path_T **patharray3, int *nhits3_primary, int *nhits3_altloc, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Chrpos_T pairmax_linear,

		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		    Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3, Pass_T pass);

extern void
Stage1hr_paired_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       Genomebits_T genomebits_in, Localdb_T localdb_in, EF64_T chromosome_ef64_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       int max_deletionlen, int max_insertlength, Chrpos_T shortsplicedist, bool splicingp_in,
		       int maxpaths_search_in, int maxpaths_report_in,
		       bool *circularp_in, int pairmax_linear_in, int pairmax_circular_in);

extern void
Stage1hr_paired_pass2_setup (int max_insertlength, Chrpos_T shortsplicedist);

#undef T
#endif

