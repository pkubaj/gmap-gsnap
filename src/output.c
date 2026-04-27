static char rcsid[] = "$Id: 2d64dbfd9aceb5fbf89815579caf5ac6286d52b7 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "output.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef GSNAP
#include "shortread.h"
#include "path-print-alignment.h"
#include "path-print-m8.h"
#include "path-print-sam.h"
#endif

#include "samheader.h"
#include "samflags.h"		/* For output types */


/* For GSNAP, now handling --failsonlyp in sam_sort.  Still
   handling quiet-if-excessive in GMAP/GSNAP, because that changes the
   SAM line to a nomapping format.  */

static Univ_IIT_T chromosome_iit;
static bool nofailsp;
static bool failsonlyp;
static bool quiet_if_excessive_p;

static int maxpaths_report;
static int quality_shift;

#ifdef GSNAP
static Outputtype_T output_type;
static bool invert_first_p;
static bool invert_second_p;
static bool single_cell_p;
static bool only_concordant_p;
static bool only_tr_consistent_p;

#else
static Printtype_T printtype;
static int invertmode;
static int wraplength;
static int ngap;
static bool nointronlenp;
static bool sam_paired_p;
static bool sam_flippedp;

static int cds_startpos;
static bool fulllengthp;
static bool truncatep;
static bool strictp;
static bool checksump;

static char *dbversion;
static char *chrsubset_name;
static Univ_IIT_T contig_iit;
static IIT_T altstrain_iit;
static bool chimeras_allowed_p;

static IIT_T map_iit;
static int *map_divint_crosstable;
static bool map_exons_p;
static bool map_bothstrands_p;
static int nflanking;
static bool print_comment_p;
#endif

static char *failedinput_root;
static char *sam_read_group_id;


void
Output_setup (Univ_IIT_T chromosome_iit_in,
	      bool nofailsp_in, bool failsonlyp_in, bool quiet_if_excessive_p_in, int maxpaths_report_in,
	      char *failedinput_root_in, int quality_shift_in,
#ifdef GSNAP
	      Outputtype_T output_type_in, bool invert_first_p_in, bool invert_second_p_in, bool single_cell_p_in,
	      bool only_concordant_p_in, bool only_tr_consistent_p_in,
#else
	      Printtype_T printtype_in, int invertmode_in, int wraplength_in, int ngap_in,
	      bool nointronlenp_in, bool sam_paired_p_in, bool sam_flippedp_in, int cds_startpos_in,
	      bool fulllengthp_in, bool truncatep_in, bool strictp_in, bool checksump_in,

	      char *dbversion_in, char *chrsubset_name_in,
	      Univ_IIT_T contig_iit_in, IIT_T altstrain_iit_in, bool chimeras_allowed_p_in,
	      IIT_T map_iit_in, int *map_divint_crosstable_in, bool map_exons_p_in,
	      bool map_bothstrands_p_in, int nflanking_in, bool print_comment_p_in,
#endif
	      char *sam_read_group_id_in) {

  chromosome_iit = chromosome_iit_in;

  nofailsp = nofailsp_in;
  failsonlyp = failsonlyp_in;
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  maxpaths_report = maxpaths_report_in;

  failedinput_root = failedinput_root_in;
  quality_shift = quality_shift_in;

#ifdef GSNAP
  output_type = output_type_in;
  invert_first_p = invert_first_p_in;
  invert_second_p = invert_second_p_in;
  single_cell_p = single_cell_p_in;
  only_concordant_p = only_concordant_p_in;
  only_tr_consistent_p = only_tr_consistent_p_in;

#else
  printtype = printtype_in;
  invertmode = invertmode_in;
  wraplength = wraplength_in;
  ngap = ngap_in;
  nointronlenp = nointronlenp_in;
  sam_paired_p = sam_paired_p_in;
  sam_flippedp = sam_flippedp_in;

  cds_startpos = cds_startpos_in;
  fulllengthp = fulllengthp_in;
  truncatep = truncatep_in;
  strictp = strictp_in;
  checksump = checksump_in;

  dbversion = dbversion_in;
  chrsubset_name = chrsubset_name_in;
  contig_iit = contig_iit_in;
  altstrain_iit = altstrain_iit_in;
  chimeras_allowed_p = chimeras_allowed_p_in;

  map_iit = map_iit_in;
  map_divint_crosstable = map_divint_crosstable_in;
  map_exons_p = map_exons_p_in;
  map_bothstrands_p = map_bothstrands_p_in;
  nflanking = nflanking_in;
  print_comment_p = print_comment_p_in;
#endif

  sam_read_group_id = sam_read_group_id_in;

  return;
}


#ifdef GSNAP
/************************************************************************
 *   Print routines and threads for GSNAP
 ************************************************************************/

/* Taken from print_result_sam from old outbuffer.c */
static Filestring_T
filestring_fromresult_sam (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			   Result_T result, Request_T request, Listpool_T listpool) {
  Filestring_T fp;
  Resulttype_T resulttype;
  Shortread_T headerseq, single_cell_infoseq, queryseq;
  Path_T *patharray, path;
  int npaths_primary, npaths_altloc, pathnum, first_absmq, second_absmq;
  char *abbrev;

  fp = Filestring_new();
  *fp_failedinput = *fp_failedinput_1 = *fp_failedinput_2 = (Filestring_T) NULL;

  resulttype = Result_resulttype(result);
  if (resulttype == SINGLEEND_NOMAPPING) {
    if (single_cell_p == true) {
      single_cell_infoseq = headerseq = Request_queryseq1(request);
      queryseq = Request_queryseq2(request);
    } else {
      single_cell_infoseq = (Shortread_T) NULL;
      headerseq = queryseq = Request_queryseq1(request);
    }

    if (nofailsp == true) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NM); /* Needs to go outside of nofailsp */
      Path_print_sam_nomapping(fp,ABBREV_NOMAPPING_1,queryseq,/*mate_queryseq*/NULL,single_cell_infoseq,
			       /*acc1*/Shortread_accession(headerseq),/*acc2*/NULL,chromosome_iit,resulttype,
			       /*first_read_p*/true,/*pathnum*/0,/*npaths_primary*/0,/*npaths_altloc*/0,
			       /*artificial_mate_p*/false,/*npaths_mate*/0,
			       /*mate*/NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq,headerseq);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    if (single_cell_p == true) {
      single_cell_infoseq = headerseq = Request_queryseq1(request);
      queryseq = Request_queryseq2(request);
    } else {
      single_cell_infoseq = (Shortread_T) NULL;
      headerseq = queryseq = Request_queryseq1(request);
    }

    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq,headerseq);
      }

    } else {
      if (single_cell_p == true) {
	single_cell_infoseq = headerseq = Request_queryseq1(request);
	queryseq = Request_queryseq2(request);
      } else {
	single_cell_infoseq = (Shortread_T) NULL;
	headerseq = queryseq = Request_queryseq1(request);
      }

      patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
      path = patharray[0];
#if 0
      if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	chrnum = 0;
	chrpos_low = 0;
      } else {
	chrpos_low = SAM_compute_chrpos(&chrnum,/*hardclip_low*/0,/*hardclip_high*/0,
					stage3,Shortread_fulllength(queryseq),/*first_read_p*/true);
      }
#endif
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
      Path_print_sam(fp,&(*fp_failedinput),abbrev,
		     path,/*acc1*/Shortread_accession(headerseq),/*acc2*/NULL,
		     /*pathnum*/1,npaths_primary,npaths_altloc,
		     path->absmq_score,first_absmq,second_absmq,path->mapq_score,chromosome_iit,
		     queryseq,/*queryseq2*/NULL,single_cell_infoseq,
		     /*pairedlength*/0,/*pair_relationship*/0,/*mate*/NULL,resulttype,
		     /*paired_read_p*/false,/*first_read_p*/true,/*artificial_mate_p*/false,
		     /*npaths_mate*/0,quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		     listpool);
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    if (single_cell_p == true) {
      single_cell_infoseq = headerseq = Request_queryseq1(request);
      queryseq = Request_queryseq2(request);
    } else {
      single_cell_infoseq = (Shortread_T) NULL;
      headerseq = queryseq = Request_queryseq1(request);
    }

    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);
    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && only_tr_consistent_p == false) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq,headerseq);
      }
      
    } else {
      if (single_cell_p == true) {
	single_cell_infoseq = headerseq = Request_queryseq1(request);
	queryseq = Request_queryseq2(request);
      } else {
	single_cell_infoseq = (Shortread_T) NULL;
	headerseq = queryseq = Request_queryseq1(request);
      }

      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
#if 0
	if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	  chrnum = 0;
	  chrpos_low = 0;
	} else {
	  chrpos_low = SAM_compute_chrpos(&chrnum,/*hardclip_low*/0,/*hardclip_high*/0,
					  stage3,Shortread_fulllength(queryseq),/*first_read_p*/true);
	}
#endif
	Path_print_sam(fp,&(*fp_failedinput),ABBREV_UNPAIRED_TRANSLOC,
		       path,/*acc1*/Shortread_accession(headerseq),
		       /*acc2*/NULL,pathnum,npaths_primary,npaths_altloc,
		       path->absmq_score,first_absmq,second_absmq,path->mapq_score,chromosome_iit,
		       queryseq,/*queryseq2*/NULL,single_cell_infoseq,
		       /*pairedlength*/0,/*pair_relationship*/0,/*mate*/NULL,resulttype,
		       /*paired_read_p*/false,/*first_read_p*/true,/*artificial_mate_p*/false,
		       /*npaths_mate*/0,quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		       listpool);
      }
    }

  } else if (resulttype == SINGLEEND_MULT) {
    if (single_cell_p == true) {
      single_cell_infoseq = headerseq = Request_queryseq1(request);
      queryseq = Request_queryseq2(request);
    } else {
      single_cell_infoseq = (Shortread_T) NULL;
      headerseq = queryseq = Request_queryseq1(request);
    }

    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
	/* No output */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	Shortread_print_query_singleend(*fp_failedinput,queryseq,headerseq);
      }

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
      if (nofailsp == true) {
	/* No output */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	Path_print_sam_nomapping(fp,ABBREV_UNPAIRED_MULT_XS,queryseq,/*mate_queryseq*/NULL,single_cell_infoseq,
				 /*acc1*/Shortread_accession(headerseq),/*acc2*/NULL,chromosome_iit,resulttype,
				 /*first_read_p*/true,/*pathnum*/1,npaths_primary,npaths_altloc,
				 /*artificial_mate_p*/false,/*npaths_mate*/0,
				 /*mate*/NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	/* Don't treat XS output as a failed input */
      }

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
#if 0
	if (Stage3end_hittype(stage3) == SAMECHR_SPLICE || Stage3end_hittype(stage3) == TRANSLOC_SPLICE) {
	  chrnum = 0;
	  chrpos_low = 0U;
	} else {
	  chrpos_low = SAM_compute_chrpos(&chrnum,/*hardclip_low*/0,/*hardclip_high*/0,
					  stage3,Shortread_fulllength(queryseq),/*first_read_p*/true);
	}
#endif
	Path_print_sam(fp,&(*fp_failedinput),ABBREV_UNPAIRED_MULT,
		       path,/*acc1*/Shortread_accession(headerseq),
		       /*acc2*/NULL,pathnum,npaths_primary,npaths_altloc,
		       path->absmq_score,first_absmq,second_absmq,path->mapq_score,chromosome_iit,
		       queryseq,/*queryseq2*/NULL,single_cell_infoseq,
		       /*pairedlength*/0,/*pair_relationship*/0,/*mate*/NULL,resulttype,
		       /*paired_read_p*/false,/*first_read_p*/true,/*artificial_mate_p*/false,
		       /*npaths_mate*/0,quality_shift,sam_read_group_id,invert_first_p,invert_second_p,
		       listpool);
      }
    }

  } else {
    /* Paired-end result */
    Path_print_sam_paired(fp,&(*fp_failedinput_1),&(*fp_failedinput_2),
			  result,resulttype,chromosome_iit,
			  Request_queryseq1(request),Request_queryseq2(request),
			  invert_first_p,invert_second_p,nofailsp,failsonlyp,
			  quality_shift,sam_read_group_id,listpool);
  }

  return fp;
}


static void
print_header_singleend (Filestring_T fp, Request_T request, bool translocationp, int npaths_primary, int npaths_altloc) {
  Shortread_T headerseq, queryseq;

  if (output_type == M8_OUTPUT) {
    /* Skip header */
  } else {
    if (single_cell_p == true) {
      /* single_cell_infoseq = */ headerseq = Request_queryseq1(request);
      queryseq = Request_queryseq2(request);
    } else {
      /* single_cell_infoseq = (Shortread_T) NULL; */
      headerseq = queryseq = Request_queryseq1(request);
    }

    FPRINTF(fp,">");
    Shortread_print_oneline(fp,queryseq);
    FPRINTF(fp,"\t%d",npaths_primary + npaths_altloc);
    if (translocationp == true) {
      FPRINTF(fp," (transloc)");
    }

    /* No sequence inversion on single-end reads */
    if (Shortread_quality_string(queryseq) != NULL) {
      FPRINTF(fp,"\t");
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    }

    FPRINTF(fp,"\t");
    Shortread_print_header(fp,headerseq,/*queryseq2*/NULL);
    /* FPRINTF(fp,"\n"); -- included in header */
  }

  return;
}


/* Taken from print_result_gsnap from old outbuffer.c */
static Filestring_T
filestring_fromresult_alignment (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
				 Result_T result, Request_T request, Listpool_T listpool) {
  Filestring_T fp;
  Resulttype_T resulttype;
  Shortread_T queryseq1, queryseq2, headerseq;
  Path_T *patharray, path;
  int npaths_primary, npaths_altloc, pathnum, first_absmq, second_absmq;
  bool print1p = false, print2p = false, single_end_p = false;

  fp = Filestring_new();
  *fp_failedinput = *fp_failedinput_1 = *fp_failedinput_2 = (Filestring_T) NULL;

  resulttype = Result_resulttype(result);
  if (resulttype == SINGLEEND_NOMAPPING) {
    single_end_p = true;
    if (nofailsp == true) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NM);
      print_header_singleend(fp,request,/*translocationp*/false,/*npaths_primary*/0,/*npaths_altloc*/0);
      FPRINTF(fp,"\n");

      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    single_end_p = true;
    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);

      patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
      path = patharray[0];

      print1p = true;
      queryseq1 = Request_queryseq1(request);

      print_header_singleend(fp,request,/*translocationp*/false,npaths_primary,npaths_altloc);
      Path_print_alignment(fp,path,/*pathpair*/NULL,queryseq1,/*invertp*/invert_first_p,listpool);
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    single_end_p = true;
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);

    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && only_tr_consistent_p == false) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else {
      print1p = true;
      queryseq1 = Request_queryseq1(request);

      print_header_singleend(fp,request,/*translocationp*/true,npaths_primary,npaths_altloc);
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
	Path_print_alignment(fp,path,/*pathpair*/NULL,queryseq1,/*invertp*/invert_first_p,listpool);
      }
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == SINGLEEND_MULT) {
    single_end_p = true;
    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* A success that is not wanted => No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      /* Not a success, but a failure => No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
      if (nofailsp == true) {
	/* No output */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	print_header_singleend(fp,request,/*translocationp*/false,npaths_primary,npaths_altloc);
	FPRINTF(fp,"\n");
	/* Don't treat XS output as a failed input */
      }

    } else {
      print1p = true;
      queryseq1 = Request_queryseq1(request);

      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      print_header_singleend(fp,request,/*translocationp*/false,npaths_primary,npaths_altloc);
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
	Path_print_alignment(fp,path,/*pathpair*/NULL,queryseq1,/*invertp*/invert_first_p,listpool);
      }
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (only_concordant_p == true || only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	queryseq1 = Request_queryseq1(request);
	queryseq2 = Request_queryseq2(request);
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }

    } else {
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      /* Stage3pair_print_end will call Filestring_set_split_output(), based on resulttype */

      /* First end */
      print1p = Pathpair_print_end_alignment(fp,result,resulttype,'>',/*firstp*/true,
					     /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
					     maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_first_p,
					     quality_shift,listpool);

      /* Second end */
      print2p = Pathpair_print_end_alignment(fp,result,resulttype,'<',/*firstp*/false,
					     /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
					     maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_second_p,
					     quality_shift,listpool);
    }

  } else {
    /* Paired-end mapping */
    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      
    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	queryseq1 = Request_queryseq1(request);
	queryseq2 = Request_queryseq2(request);
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }
    
    } else {
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      /* Stage3pair_print_end will call Filestring_set_split_output() based on resulttype */

      /* First end */
      print1p = Pathpair_print_end_alignment(fp,result,resulttype,'>',/*firstp*/true,
					     /*queryseq*/queryseq1,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
					     maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_first_p,quality_shift,
					     listpool);

      /* Second end */
      print2p = Pathpair_print_end_alignment(fp,result,resulttype,'<',/*firstp*/false,
					     /*queryseq*/queryseq2,/*headerseq1*/queryseq1,/*headerseq2*/queryseq2,
					     maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_second_p,quality_shift,
					     listpool);
    }
  }

  if (single_end_p == true) {
    if (failedinput_root != NULL && print1p == failsonlyp) {
      *fp_failedinput = Filestring_new();
      headerseq = queryseq1 = Request_queryseq1(request);
      Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
    }

  } else {
    if (failedinput_root != NULL && print1p == failsonlyp && print2p == failsonlyp) {
      *fp_failedinput_1 = Filestring_new();
      *fp_failedinput_2 = Filestring_new();
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
    }
  }

  return fp;
}


static Filestring_T
filestring_fromresult_m8 (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			  Result_T result, Request_T request, Listpool_T listpool) {
  Filestring_T fp;
  Resulttype_T resulttype;
  Shortread_T headerseq, queryseq1, queryseq2;
  Path_T *patharray, path;
  int npaths_primary, npaths_altloc, pathnum, first_absmq, second_absmq;
  bool print1p = false, print2p = false, single_end_p = false;

  fp = Filestring_new();
  *fp_failedinput = *fp_failedinput_1 = *fp_failedinput_2 = (Filestring_T) NULL;

  resulttype = Result_resulttype(result);
  if (resulttype == SINGLEEND_NOMAPPING) {
    single_end_p = true;
    /* Skip: M8 does not print failed alignments */
    if (failedinput_root != NULL) {
      *fp_failedinput = Filestring_new();
      headerseq = queryseq1 = Request_queryseq1(request);
      Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
    }

  } else if (resulttype == SINGLEEND_UNIQ) {
    single_end_p = true;

    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);

      patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
      path = patharray[0];

      print1p = true;
      headerseq = Request_queryseq1(request);
      Path_print_m8(fp,path,/*accession*/Shortread_accession(headerseq),/*acc_suffix*/"",
		    /*invertp*/invert_first_p,listpool);
    }

  } else if (resulttype == SINGLEEND_TRANSLOC) {
    single_end_p = true;
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);

    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && only_tr_consistent_p == false) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else {
      print1p = true;
      headerseq = Request_queryseq1(request);
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
	Path_print_m8(fp,path,/*accession*/Shortread_accession(headerseq),/*acc_suffix*/"",
		      /*invertp*/invert_first_p,listpool);
      }
    }

  } else if (resulttype == SINGLEEND_MULT) {
    single_end_p = true;
    patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
	/* No output */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput = Filestring_new();
	headerseq = queryseq1 = Request_queryseq1(request);
	Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
      }

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
      if (nofailsp == true) {
	/* No output */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	/* Don't treat XS output as a failed input */
      }

    } else {
      print1p = true;
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      headerseq = Request_queryseq1(request);
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	path = patharray[pathnum-1];
	Path_print_m8(fp,path,/*accession*/Shortread_accession(headerseq),/*acc_suffix*/"",
		      /*invertp*/invert_first_p,listpool);
      }
    }

  } else if (resulttype == PAIREDEND_NOMAPPING) {
    /* Skip: M8 does not print failed alignments */
    if (only_concordant_p == true || only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	queryseq1 = Request_queryseq1(request);
	queryseq2 = Request_queryseq2(request);
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }
    }
      
  } else {
    if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	queryseq1 = Request_queryseq1(request);
	queryseq2 = Request_queryseq2(request);
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }
    
    } else {
      /* Stage3pair_print_end will call Filestring_set_split_output() based on resulttype */
      headerseq = Request_queryseq1(request);

      /* First end */
      print1p = Pathpair_print_end_m8(fp,result,resulttype,/*firstp*/true,
				      /*accession*/Shortread_accession(headerseq),
				      maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_first_p,
				      listpool);

      /* Second end */
      print2p = Pathpair_print_end_m8(fp,result,resulttype,/*firstp*/false,
				      /*accession*/Shortread_accession(headerseq),
				      maxpaths_report,quiet_if_excessive_p,/*invertp*/invert_second_p,
				      listpool);
    }
  }

  if (single_end_p == true) {
    if (failedinput_root != NULL && print1p == failsonlyp) {
      *fp_failedinput = Filestring_new();
      headerseq = queryseq1 = Request_queryseq1(request);
      Shortread_print_query_singleend(*fp_failedinput,queryseq1,headerseq);
    }

  } else {
    if (failedinput_root != NULL && print1p == failsonlyp && print2p == failsonlyp) {
      *fp_failedinput_1 = Filestring_new();
      *fp_failedinput_2 = Filestring_new();
      queryseq1 = Request_queryseq1(request);
      queryseq2 = Request_queryseq2(request);
      Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
    }
  }

  return fp;
}


Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
			      Result_T result, Request_T request, Listpool_T listpool) {
  if (output_type == SAM_OUTPUT) {
    return filestring_fromresult_sam(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				     result,request,listpool);
  } else if (output_type == M8_OUTPUT) {
    return filestring_fromresult_m8(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				    result,request,listpool);
  } else {
    return filestring_fromresult_alignment(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
					   result,request,listpool);
  }
}

#else
/************************************************************************
 *   Print routines and threads for GMAP
 ************************************************************************/

static void
print_npaths (Filestring_T fp, int npaths, char *chrsubset_name,
	      Chimera_T chimera, Failure_T failuretype) {

  if (npaths == 0) {
    FPRINTF(fp,"Paths (0):");
#if 0
  } else if (mergedp == true) {
    FPRINTF(fp,"Paths (1):");
#endif
  } else {
    FPRINTF(fp,"Paths (%d):",npaths);
  }
  if (chrsubset_name != NULL) {
    FPRINTF(fp,"  [chrsubset: %s]",chrsubset_name);
  }
  if (failuretype == NO_FAILURE) {
    if (chimera != NULL) {
      Chimera_print(fp,chimera);
    }
  } else if (failuretype == EMPTY_SEQUENCE) {
    FPRINTF(fp," *** Empty sequence ***");
  } else if (failuretype == SHORT_SEQUENCE) {
    FPRINTF(fp," *** Short sequence < index oligo size ***");
  } else if (failuretype == POOR_SEQUENCE) {
    FPRINTF(fp," *** Poor sequence (use -p flag to change pruning behavior) ***");
  } else if (failuretype == REPETITIVE) {
    FPRINTF(fp," *** Repetitive sequence (use -p flag to change pruning behavior) ***");
  }
  FPRINTF(fp,"\n");
  if (npaths == 0) {
    FPRINTF(fp,"\n");
  }
  return;
}


/* Taken from Outbuffer_print_result */
Filestring_T
Output_filestring_fromresult (Filestring_T *fp_failedinput, Result_T result, Request_T request,
			      Sequence_T headerseq) {
  Filestring_T fp;
  char *abbrev;
  Sequence_T queryseq;
  Genome_T genome, genomealt;
  Stage3_T *stage3array;
  int npaths_primary, npaths_altloc, pathnum, effective_maxpaths, first_absmq, second_absmq;
  Chimera_T chimera = NULL;
  /* int chimerapos, chimeraequivpos, chimera_cdna_direction; */
  int querylength;
  /* double donor_prob, acceptor_prob; */
  /* bool mergedp = false; */
  bool printp = true;

  Sequence_T genome_sequence;
  Univcoord_T genomiclength;
  char *chrstring;


  fp = Filestring_new();
  if (failedinput_root == NULL) {
    *fp_failedinput = (Filestring_T) NULL;
  } else {
    *fp_failedinput = Filestring_new();
  }

  queryseq = Request_queryseq(request);
  querylength = Sequence_fulllength_given(queryseq);

  genome = Request_genome(request);
  genomealt = Request_genomealt(request);

  stage3array = Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

  /* chimerapos = chimeraequivpos = -1; */
  /* chimera_cdna_direction = 0; */
  /* donor_prob = acceptor_prob = 0.0; */

  /* Translation */
  if (npaths_primary + npaths_altloc == 0) {
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NM);
    abbrev = ABBREV_NOMAPPING_1;
    effective_maxpaths = 0;
    if (nofailsp == true) {
      printp = false;
    }

    if (Result_failuretype(result) == POOR_SEQUENCE) {
      fprintf(stderr,"Accession %s skipped (poor sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else if (Result_failuretype(result) == REPETITIVE) {
      fprintf(stderr,"Accession %s skipped (repetitive sequence).  Use -p flag to change pruning behavior\n",Sequence_accession(headerseq));
    } else {
      fprintf(stderr,"No paths found for %s\n",Sequence_accession(headerseq));
    }

#if 0
  } else if ((mergedp = Result_mergedp(result)) == true) {
    if (Stage3_circularpos(stage3array[0]) > 0) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UC);
      abbrev = ABBREV_UNPAIRED_CIRCULAR;
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }
    effective_maxpaths = 1;
    if (failsonlyp == true) {
      printp = false;
    } else {
      Stage3_translate(stage3array[0],
#ifdef PMAP
		       queryseq,
#endif
		       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }
#endif

  } else if ((chimera = Result_chimera(result)) != NULL) {
    if (chimeras_allowed_p == true) {
      effective_maxpaths = 2;
    } else {
      effective_maxpaths = 0;
    }
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);
    abbrev = ABBREV_UNPAIRED_TRANSLOC;

    if (failsonlyp == true) {
      printp = false;
    } else {
      /* chimerapos = Chimera_pos(chimera); */
      /* chimeraequivpos = Chimera_equivpos(chimera); */
      /* donor_prob = Chimera_donor_prob(chimera); */
      /* acceptor_prob = Chimera_acceptor_prob(chimera); */
      /* chimera_cdna_direction = Chimera_cdna_direction(chimera); */

      Stage3_translate_chimera(stage3array[0],stage3array[1],
#ifdef PMAP
			       queryseq,
#endif
			       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }

  } else if (maxpaths_report == 0) {
    effective_maxpaths = 1;
    if (npaths_primary + npaths_altloc > 1) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      abbrev = ABBREV_UNPAIRED_MULT;
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (failsonlyp == true) {
      printp = false;
    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
      printp = false;
    } else {
      Stage3_translate(stage3array[0],
#ifdef PMAP
		       queryseq,
#endif
		       querylength,fulllengthp,cds_startpos,truncatep,strictp);
    }

  } else {
    if (npaths_primary + npaths_altloc > 1) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      abbrev = ABBREV_UNPAIRED_MULT;
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;
    }

    if (npaths_primary + npaths_altloc < maxpaths_report) {
      effective_maxpaths = npaths_primary + npaths_altloc;
    } else {
      effective_maxpaths = maxpaths_report;
    }

    if (failsonlyp == true) {
      printp = false;
    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
      printp = false;
    } else {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_translate(stage3array[pathnum-1],
#ifdef PMAP
			 queryseq,
#endif
			 querylength,fulllengthp,cds_startpos,truncatep,strictp);
      }
    }
  }

  /* Printing */
  if (printp == false) {
    /* No output, either because of --nofails or --quiet-if-excessive */

  } else {
    if (*fp_failedinput != NULL &&
	(npaths_primary + npaths_altloc == 0 && quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report)) {
      PUTC('>',*fp_failedinput);
      Sequence_print_header(*fp_failedinput,headerseq,checksump);
      Sequence_print(*fp_failedinput,queryseq,/*uppercasep*/false,wraplength,/*trimmedp*/false);
    }

    if (printtype == SIMPLE || printtype == SUMMARY || printtype == ALIGNMENT) {
      /* Print header, even if no alignment is found */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);

      if (npaths_primary + npaths_altloc == 0) {
	print_npaths(fp,/*npaths*/0,chrsubset_name,
		     /*chimera*/NULL,Result_failuretype(result));


      } else {
	print_npaths(fp,npaths_primary + npaths_altloc,chrsubset_name,chimera,NO_FAILURE);
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_pathsummary(fp,stage3array[pathnum-1],pathnum,
				   chromosome_iit,contig_iit,
				   altstrain_iit,queryseq,dbversion);
	}
      }

      if (printtype != SIMPLE) {
	FPRINTF(fp,"Alignments:\n");
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  FPRINTF(fp,"  Alignment for path %d:\n\n",pathnum);
	  Stage3_print_alignment(fp,stage3array[pathnum-1],
				 genome,genomealt,chromosome_iit,printtype,
				 /*continuousp*/false,/*continuous_by_exon_p*/false,
				 /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
	}
      }

      if (map_iit != NULL) {
	FPRINTF(fp,"Maps:\n");
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_map(fp,stage3array[pathnum-1],map_iit,map_divint_crosstable,
			   chromosome_iit,pathnum,map_exons_p,map_bothstrands_p,
			   nflanking,print_comment_p);
	}
      }

    } else if (printtype == CONTINUOUS) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      if (npaths_primary + npaths_altloc == 0) {
	FPRINTF(fp,"\n\n\n");
      } else {
	Stage3_print_alignment(fp,stage3array[0],genome,genomealt,chromosome_iit,printtype,
			       /*continuousp*/true,/*continuous_by_exon_p*/false,
			       /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
      }

    } else if (printtype == CONTINUOUS_BY_EXON) {
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      print_npaths(fp,npaths_primary + npaths_altloc,chrsubset_name,chimera,NO_FAILURE);
      if (npaths_primary + npaths_altloc == 0) {
	FPRINTF(fp,"\n\n\n");
      } else {
	Stage3_print_pathsummary(fp,stage3array[0],/*pathnum*/1,
				 chromosome_iit,contig_iit,
				 altstrain_iit,queryseq,dbversion);
	FPRINTF(fp,"Alignments:\n");
	FPRINTF(fp,"  Alignment for path %d:\n\n",/*pathnum*/1);
	Stage3_print_alignment(fp,stage3array[0],genome,genomealt,chromosome_iit,printtype,
			       /*continuousp*/false,/*continuous_by_exon_p*/true,
			       /*flipgenomep*/true,invertmode,nointronlenp,wraplength);
      }

    } else if (printtype == EXONS_CDNA || printtype == EXONS_CDNA_WINTRONS) {
      /* Introns are printed with ngap set to infinity, and not printed with ngap == 0 */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,"<path %d>\n",pathnum);
	Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			 wraplength,ngap,/*cdna*/true);
	FPRINTF(fp,"</path>\n");
      }

    } else if (printtype == EXONS_GENOMIC || printtype == EXONS_GENOMIC_WINTRONS) {
      /* Introns are printed with ngap set to infinity, and not printed with ngap == 0 */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,"<path %d>\n",pathnum);
	Pair_print_exons(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			 wraplength,ngap,/*cdna*/false);
	FPRINTF(fp,"</path>\n");
      }

    } else if (printtype == MASK_INTRONS) {
      /* Print only best path */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      Pair_print_mask_introns(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			      Stage3_chrlength(stage3array[0]),wraplength,/*include_utr_p*/false);

    } else if (printtype == MASK_UTR_INTRONS) {
      /* Print only best path */
      PUTC('>',fp);
      Sequence_print_header(fp,headerseq,checksump);
      Pair_print_mask_introns(fp,Stage3_pairarray(stage3array[0]),Stage3_npairs(stage3array[0]),
			      Stage3_chrlength(stage3array[0]),wraplength,/*include_utr_p*/true);

    } else if (printtype == CDNA) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	PUTC('>',fp);
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_cdna(fp,stage3array[pathnum-1],wraplength);
      }

    } else if (printtype == PROTEIN_GENOMIC) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	PUTC('>',fp);
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_protein_genomic(fp,stage3array[pathnum-1],wraplength);
      }

    } else if (printtype == PSL_NT) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pslformat_nt(fp,stage3array[pathnum-1],
				  genome,chromosome_iit,queryseq);
      }

#ifdef PMAP
    } else if (printtype == PSL_PRO) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_pslformat_pro(fp,stage3array[pathnum-1],
				   genome,chromosome_iit,queryseq,strictp);
      }
#endif

    } else if (printtype == GFF3_GENE || printtype == GFF3_MATCH_CDNA ||
	       printtype == GFF3_MATCH_EST) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_gff3(fp,stage3array[pathnum-1],pathnum,
			  genome,chromosome_iit,queryseq,querylength,printtype,
			  /*sourcename*/Genome_accession(genome) != NULL ? Genome_accession(genome) : dbversion);
      }

#ifndef PMAP
    } else if (printtype == BEDPE) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_bedpe(fp,stage3array[pathnum-1],chromosome_iit);
      }

    } else if (printtype == SAM) {
      if (npaths_primary + npaths_altloc == 0) {
	if (sam_flippedp) {
	  genomiclength = Genome_genomelength(genome);
	  genome_sequence = Genome_get_segment(genome,/*genomicstart*/0,genomiclength,chromosome_iit,/*revcomp*/false);
	  chrstring = Genome_accession(genome);

	  Pair_print_sam_nomapping_flipped(fp,abbrev,chrstring,
					   Sequence_fullpointer(genome_sequence),(int) genomiclength,
					   Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);
	  Sequence_free(&genome_sequence);

	} else {
	  Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
				   Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
				   Sequence_fulllength(queryseq),quality_shift,
				   Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);
	}

      } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
	if (sam_flippedp) {
	  genomiclength = Genome_genomelength(genome);
	  genome_sequence = Genome_get_segment(genome,/*genomicstart*/0,genomiclength,chromosome_iit,/*revcomp*/false);
	  chrstring = Genome_accession(genome);

	  Pair_print_sam_nomapping_flipped(fp,abbrev,chrstring,
					   Sequence_fullpointer(genome_sequence),(int) genomiclength,
					   Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);
	  Sequence_free(&genome_sequence);

	} else {
	  Pair_print_sam_nomapping(fp,abbrev,/*acc1*/Sequence_accession(headerseq),/*acc2*/NULL,
				   Sequence_fullpointer(queryseq),Sequence_quality_string(queryseq),
				   Sequence_fulllength(queryseq),quality_shift,
				   Sequence_firstp(queryseq),sam_paired_p,sam_read_group_id);
	}
#if 0
      } else if (mergedp == true) {
	Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
			 Stage3_absmq_score(stage3array[0]),second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 genome,chromosome_iit,queryseq,
			 /*chimera_part*/0,/*chimera*/NULL,quality_shift,sam_paired_p,
			 sam_read_group_id);
#endif

      } else if (chimera != NULL) {
	Stage3_print_sam(fp,abbrev,stage3array[0],/*pathnum*/1,npaths_primary,npaths_altloc,
			 Stage3_absmq_score(stage3array[0]),second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 genome,chromosome_iit,queryseq,
			 /*chimera_part*/-1,chimera,quality_shift,sam_paired_p,
			 sam_read_group_id);
	Stage3_print_sam(fp,abbrev,stage3array[1],/*pathnum*/1,npaths_primary,npaths_altloc,
			 Stage3_absmq_score(stage3array[0]),second_absmq,
			 Stage3_mapq_score(stage3array[0]),
			 genome,chromosome_iit,queryseq,
			 /*chimera_part*/+1,chimera,quality_shift,sam_paired_p,
			 sam_read_group_id);

      } else {
	for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	  Stage3_print_sam(fp,abbrev,stage3array[pathnum-1],pathnum,npaths_primary,npaths_altloc,
			   Stage3_absmq_score(stage3array[pathnum-1]),second_absmq,
			   Stage3_mapq_score(stage3array[pathnum-1]),
			   genome,chromosome_iit,queryseq,
			   /*chimera_part*/0,/*chimera*/NULL,quality_shift,sam_paired_p,
			   sam_read_group_id);
	}
      }
#endif

    } else if (printtype == COORDS) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	FPRINTF(fp,">");
	Sequence_print_header(fp,headerseq,checksump);
	Stage3_print_coordinates(fp,stage3array[pathnum-1],chromosome_iit,invertmode);
      }

    } else if (printtype == SPLICESITES) {
      /* Print only best path */
      if (npaths_primary + npaths_altloc > 0) {
	Stage3_print_splicesites(fp,stage3array[0],chromosome_iit,queryseq);
      }

    } else if (printtype == INTRONS) {
      /* Print only best path */
      if (npaths_primary + npaths_altloc > 0) {
	Stage3_print_introns(fp,stage3array[0],chromosome_iit,queryseq);
      }

    } else if (printtype == MAP_RANGES) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_iit_map(fp,stage3array[pathnum-1],chromosome_iit,queryseq);
      }
      
    } else if (printtype == MAP_EXONS) {
      for (pathnum = 1; pathnum <= effective_maxpaths; pathnum++) {
	Stage3_print_iit_exon_map(fp,stage3array[pathnum-1],chromosome_iit,queryseq);
      }

    } else {
      fprintf(stderr,"Unexpected printtype %d\n",printtype);
      abort();

    }
  }

  return fp;
}

#endif


