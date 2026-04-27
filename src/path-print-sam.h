/* $Id: 1c1af9f9f951624da3f56d1be3078ac92a7e28f9 $ */
#ifndef PATH_PRINT_SAM_INCLUDED
#define PATH_PRINT_SAM_INCLUDED

#include "path.h"
#include "pathpair.h"

#include "bool.h"
#include "univcoord.h"
#include "filestring.h"
#include "shortread.h"
#include "iit-read-univ.h"
#include "resulthr.h"

#include "listpool.h"


#define T Path_T

void
Path_print_sam_nomapping (Filestring_T fp, char *abbrev, Shortread_T queryseq, Shortread_T mate_queryseq,
			  Shortread_T single_cell_infoseq,
			  char *acc1, char *acc2, Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
			  int pathnum, int npaths_primary, int npaths_altloc, bool artificial_mate_p, int npaths_mate,

			  T mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p);

void
Path_print_sam (Filestring_T fp, Filestring_T *fp_failedinput, char *abbrev,
		T this, char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
		int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit,
		Shortread_T queryseq, Shortread_T queryseq_mate, Shortread_T single_cell_infoseq,
		int pairedlength, int pair_relationship, T mate, Resulttype_T resulttype,
		bool paired_read_p, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		Listpool_T listpool);

void
Path_print_sam_paired (Filestring_T fp, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
		       Result_T result, Resulttype_T resulttype, Univ_IIT_T chromosome_iit,
		       Shortread_T queryseq1, Shortread_T queryseq2,
		       bool invert_first_p, bool invert_second_p,
		       bool nofailsp, bool failsonlyp, int quality_shift, char *sam_read_group_id,
		       Listpool_T listpool);

extern void
Path_print_sam_setup (bool add_paired_nomappers_p_in, bool paired_flag_means_concordant_p_in, bool sam_insert_0M_p_in,
		      bool quiet_if_excessive_p_in, int maxpaths_report_in,
		      char *failedinput_root_in, bool fastq_format_p_in, bool extend_soft_clips_p_in, bool method_print_p_in,
		      bool only_concordant_p_in, bool omit_concordant_uniq_p_in, bool omit_concordant_mult_p_in, 
		      bool only_tr_consistent_p_in,
		      bool *circularp_in, bool clip_overlap_p_in, bool merge_overlap_p_in, bool merge_samechr_p_in,
		      bool sam_multiple_primaries_p_in, bool sam_sparse_secondaries_p_in,
		      Univ_IIT_T chromosome_iit_in, Univ_IIT_T transcript_iit_in,
		      IIT_T snps_iit_in, bool maskedp_in);

#undef T
#endif

