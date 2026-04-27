/* $Id: b4214026f8f7656df0a0b09735914ced5f36ec68 $ */
#ifndef OUTPUT_INCLUDED
#define OUTPUT_INCLUDED

#include "types.h"
#include "outputtype.h"
#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "iit-read.h"
#include "filestring.h"
#include "listpool.h"

#include "request.h"
#include "mem.h"		/* To get MEMUSAGE */

#ifdef GSNAP
#include "resulthr.h"
#else
#include "result.h"
#include "sequence.h"
#endif

extern void
Output_setup (Univ_IIT_T chromosome_iit_in,
	      bool nofailsp_in, bool failsonlyp_in, bool quiet_if_excessive_p_in, int maxpaths_report_in,
	      char *failedinput_root_in, int quality_shift_in,
#ifdef GSNAP
	      Outputtype_T output_type_in, bool invert_first_p_in, bool invert_second_p_in, bool single_cell_p_in,
	      bool only_concordant_p_in, bool only_tr_consistent_p_in,
#else
	      Printtype_T printtype_in, int invertmode_in, int wraplength_in, int ngap_in,
	      bool nointronlenp_in, bool sam_paired_p_in, bool sam_flippedp, int cds_startpos_in,
	      bool fulllengthp_in, bool truncatep_in, bool strictp_in, bool checksump_in,

	      char *dbversion_in, char *chrsubset_name_in,
	      Univ_IIT_T contig_iit_in, IIT_T altstrain_iit_in, bool chimeras_allowed_p_in,
	      IIT_T map_iit_in, int *map_divint_crosstable_in, bool map_exons_p_in,
	      bool map_bothstrands_p_in, int nflanking_in, bool print_comment_p_in,
#endif
	      char *sam_read_group_id_in);


#ifdef GSNAP
extern Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			      Result_T result, Request_T request, Listpool_T listpool);
#else
extern Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput, Result_T result, Request_T request,
			      Sequence_T headerseq);
#endif


#endif

