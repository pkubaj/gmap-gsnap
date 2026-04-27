static char rcsid[] = "$Id: 5f3588790dcf2b71b0ff218b1afc12f5cadc3009 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-print-sam.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For isupper, toupper */

#include "bool.h"
#include "assert.h"
#include "complement.h"
#include "univcoord.h"
#include "junction.h"
#include "mapq.h"
#include "simplepair.h"
#include "transcript.h"

#include "single-cell.h"


static bool add_paired_nomappers_p;
static bool paired_flag_means_concordant_p;
static bool sam_insert_0M_p;

static bool quiet_if_excessive_p;
static int maxpaths_report;

static char *failedinput_root;
static bool fastq_format_p;
static bool extend_soft_clips_p;
static bool method_print_p;

static bool only_tr_consistent_p;
static bool only_concordant_p;
static bool omit_concordant_uniq_p;
static bool omit_concordant_mult_p;

static bool *circularp;

static bool clip_overlap_p;
static bool merge_overlap_p;
static bool merge_samechr_p;

static bool sam_multiple_primaries_p;
static bool sam_sparse_secondaries_p;

static Univ_IIT_T chromosome_iit;
static Univ_IIT_T transcript_iit;

static IIT_T snps_iit;
static bool maskedp;

static char complCode[128] = COMPLEMENT_LC;


#define PRINT_ALTS_COORDS 1


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Printing of fusions */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif



#define T Path_T


static unsigned int
compute_flag (bool plusp, T mate, Resulttype_T resulttype,
	      bool first_read_p, int pathnum, int npaths, bool artificial_mate_p, int npaths_mate,
	      int absmq_score, int first_absmq, bool invertp, bool invert_mate_p,
	      bool supplementaryp) {
  unsigned int flag = 0U;


  debug(printf("compute_flag: resulttype %s, mate %p\n",
	       Resulttype_string(resulttype),mate));

  if (npaths == 0) {
    debug(printf("npaths = 0, so QUERY_UNMAPPED %d\n",QUERY_UNMAPPED));
    flag |= QUERY_UNMAPPED;
  } else if (plusp == invertp) {
    debug(printf("plusp %d and invertp %d, so QUERY_MINUSP %d\n",
		 plusp,invertp,QUERY_MINUSP));
    flag |= QUERY_MINUSP;
  }

  if (resulttype == SINGLEEND_NOMAPPING || resulttype == SINGLEEND_UNIQ ||
      resulttype == SINGLEEND_TRANSLOC || resulttype == SINGLEEND_MULT) {
    /* No first or second read or mate */
  } else {
    debug(printf("PAIRED_READ %d\n",PAIRED_READ));
    flag |= PAIRED_READ;
    if (first_read_p == true) {
      debug(printf("FIRST_READ %d\n",FIRST_READ_P));
      flag |= FIRST_READ_P;
    } else {
      debug(printf("SECOND_READ %d\n",SECOND_READ_P));
      flag |= SECOND_READ_P;
    }
    if (artificial_mate_p == true) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (npaths_mate == 0) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (quiet_if_excessive_p && npaths_mate > maxpaths_report) {
      debug(printf("MATE_UNMAPPED %d\n",MATE_UNMAPPED));
      flag |= MATE_UNMAPPED;

    } else if (mate == (Path_T) NULL) {
      /* Unpaired; no mate.  Not clear if should be MATE_UNMAPPED. */
      /* Picard says MATE_UNMAPPED flag should not be set for unpaired reads */

    } else {
      if (mate->plusp == invert_mate_p) {
	debug(printf("MATE_MINUSP %d\n",MATE_MINUSP));
	flag |= MATE_MINUSP;
      }

      if (npaths == 0) {
	/* Need to check npaths == 0 in case clipping of overlaps results in a nomapping */

      } else if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
	flag |= PAIRED_MAPPING;

      } else if (resulttype == PAIRED_UNIQ_INV || resulttype == PAIRED_UNIQ_SCR ||
		 resulttype == PAIRED_UNIQ_TOOLONG || resulttype == PAIRED_MULT) {
	/* Note: We are counting PAIRED_UNIQ and PAIRED_MULT as "paired" mappings.
	   However, we are no longer counting UNPAIRED_UNIQ as a "paired" mapping. */
	if (paired_flag_means_concordant_p == true) {
	  /* Not concordant, so don't turn on paired flag */
	} else {
	  debug(printf("PAIRED_MAPPING %d\n",PAIRED_MAPPING));
	  flag |= PAIRED_MAPPING;
	}

      } else {
	/* Not concordant or paired */
      }
    }
  }

  if (pathnum > 1) {
    if (sam_multiple_primaries_p == false) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else if (absmq_score != first_absmq) {
      debug(printf("NOT_PRIMARY %d\n",NOT_PRIMARY));
      flag |= NOT_PRIMARY;
    } else {
      /* Just as good as first alignment, so don't mark as altloc */
    }
  }

  if (supplementaryp == true) {
    flag |= SUPPLEMENTARY;
  }

  debug(printf("Returning flag %d\n",flag));
  return flag;
}


static void
Path_cigar_md (int *hardclip_low, int *hardclip_high,
	       Filestring_T *cigar_fp, Filestring_T *md_fp, T this,
	       Shortread_T queryseq) {
  Intlist_T q;
  List_T j;
  Univcoordlist_T u;
  Univcoord_T deletionpos;
  Junction_T junction;
  Junctiontype_T type;
  char *deletion_string;
  int qpos, nconsecutive, n, i;
  int ninserts, softclip, choplength;
  char *p, c;			/* p is needed to compute MD string */
  bool md_startp = true;

  debug(printf("Entered Path_cigar_md\n"));

  *cigar_fp = Filestring_new();
  *md_fp = Filestring_new();

  if (this == NULL) {
    *hardclip_low = *hardclip_high = 0;
    FPRINTF(*cigar_fp,"*");
    FPRINTF(*md_fp,"*");

  } else {
    ninserts = 0;
    nconsecutive = 0;

    assert(this->genomic_diff != NULL);
    p = this->genomic_diff;
    q = this->endpoints;
    j = this->junctions;
    u = this->univdiagonals;

    /* Compute softclip from chop */
    if (this->plusp == true) {
      choplength = Shortread_left_choplength(queryseq);
    } else {
      choplength = Shortread_right_choplength(queryseq);
    }	

    /* Compute softclip from alignment and hardclip */
    qpos = Intlist_head(q);
    if ((this->plusp == true && this->fusion_querystart_junction != NULL) ||
	(this->plusp == false && this->fusion_queryend_junction != NULL)) {
      *hardclip_low = choplength + qpos;
      softclip = 0;

    } else if (this->circular_endpoints != NULL && this->circular_high_p == false) {
      *hardclip_low = choplength + qpos;
      softclip = 0;

    } else {
      *hardclip_low = 0;
      softclip = choplength + qpos;
    }
    p += qpos;			/* Ignore choplength, which is not in genomic_diff */

    if (softclip > 0) {
      FPRINTF(*cigar_fp,"%dS",softclip);
    }
    
    if (*hardclip_low > 0) {
      FPRINTF(*cigar_fp,"%dH",*hardclip_low);
    }

    while (j != NULL) {
      q = Intlist_next(q);
      debug(printf("n %d = endpoint %d - qpos %d, nconsecutive %d\n",
		   Intlist_head(q) - qpos,Intlist_head(q),qpos,nconsecutive));
      n = Intlist_head(q) - qpos;
      if (n > 0 || sam_insert_0M_p == true) {
	FPRINTF(*cigar_fp,"%dM",n);
      }
      for (i = 0; i < n; i++) {
	debug(printf("isupper %c, nconsecutive %d\n",*p,nconsecutive));
	if (isupper(c = *p++)) {
	  /* Match */
	  nconsecutive++;

	} else if (nconsecutive > 0) {
	  /* Mismatch */
	  FPRINTF(*md_fp,"%d%c",nconsecutive,toupper(c));
	  nconsecutive = 0;
	  md_startp = false;
	  
	} else if (md_startp == false) {
	  /* Consecutive mismatches */
	  FPRINTF(*md_fp,"%c",toupper(c));
	  
	} else {
	  /* Initial mismatch: MD string needs to start with a number */
	  FPRINTF(*md_fp,"0%c",toupper(c));
	  md_startp = false;
	}
      }
      debug(printf("After cigar M, nconsecutive is %d\n",nconsecutive));

      ninserts = 0;
      junction = (Junction_T) List_head(j);
      if (junction == JUNCTION_UNSOLVED) {
	FPRINTF(*cigar_fp,"X");

      } else if ((type = Junction_type(junction)) == DEL_JUNCTION) {
	FPRINTF(*cigar_fp,"%dD",Junction_nindels(junction));

	deletionpos = (Univcoordlist_head(u) - this->querylength) + /*qpos*/Intlist_head(q);
	deletion_string = Junction_deletion_string(junction,deletionpos);
	FPRINTF(*md_fp,"%d^%s",nconsecutive,deletion_string);
	FREE(deletion_string);

	nconsecutive = 0;
	md_startp = false;

      } else if (type == INS_JUNCTION) {
	FPRINTF(*cigar_fp,"%dI",ninserts = Junction_nindels(junction));
	p += ninserts;

      } else if (type == SPLICE_JUNCTION) {
	FPRINTF(*cigar_fp,"%dN",Junction_splice_distance(junction));

      } else {
	fprintf(stderr,"Unknown junction type %d\n",type);
	abort();
      }
      qpos = Intlist_head(q) + ninserts;

      j = List_next(j);
      u = Univcoordlist_next(u);
    }

    q = Intlist_next(q);
    debug(printf("n %d = endpoint %d - qpos %d, nconsecutive %d\n",
		 Intlist_head(q) - qpos,Intlist_head(q),qpos,nconsecutive));
    n = Intlist_head(q) - qpos;
    if (n > 0 || sam_insert_0M_p == true) {
      FPRINTF(*cigar_fp,"%dM",n);
    }
    for (i = 0; i < n; i++) {
      debug(printf("isupper %c, nconsecutive %d\n",*p,nconsecutive));
      if (isupper(c = *p++)) {
	/* Match */
	nconsecutive++;

      } else if (nconsecutive > 0) {
	/* Mismatch */
	FPRINTF(*md_fp,"%d%c",nconsecutive,toupper(c));
	nconsecutive = 0;
	md_startp = false;
	  
      } else if (md_startp == false) {
	/* Consecutive mismatches */
	FPRINTF(*md_fp,"%c",toupper(c));
	  
      } else {
	/* Initial mismatch: MD string needs to start with a number */
	FPRINTF(*md_fp,"0%c",toupper(c));
	md_startp = false;
      }
    }
    debug(printf("After cigar M, nconsecutive is %d\n",nconsecutive));

    /* Previously checked for nconsecutive > 0, but MD string needs to end with a number */
    FPRINTF(*md_fp,"%d",nconsecutive);

    /* Compute softclip from chop */
    if (this->plusp == true) {
      choplength = Shortread_right_choplength(queryseq);
    } else {
      choplength = Shortread_left_choplength(queryseq);
    }	

    /* Compute softclip from alignment and hardclip */
    qpos = Intlist_head(q);
    if ((this->plusp == true && this->fusion_queryend_junction != NULL) ||
	(this->plusp == false && this->fusion_querystart_junction != NULL)) {
      *hardclip_high = (this->querylength - qpos) + choplength;
      softclip = 0;

    } else if (this->circular_endpoints != NULL && this->circular_high_p == true) {
      *hardclip_high = (this->querylength - qpos) + choplength;
      softclip = 0;

    } else {
      *hardclip_high = 0;
      softclip = (this->querylength - qpos) + choplength;
    }

    if (*hardclip_high > 0) {
      FPRINTF(*cigar_fp,"%dH",*hardclip_high);
    }

    if (softclip > 0) {
      FPRINTF(*cigar_fp,"%dS",softclip);
    }
  }

  /* Filestring_stringify required before Filestring_merge */
  Filestring_stringify(*cigar_fp);
  Filestring_stringify(*md_fp);

  return;
}



static void
Path_supplemental_cigar_md (int *hardclip_low, int *hardclip_high,
			    Filestring_T *cigar_fp, Filestring_T *md_fp,
			    char *genomic_diff, Intlist_T suppl_endpoints, List_T suppl_junctions,
			    Univcoordlist_T suppl_univdiagonals, bool segment_plusp, bool query_plusp,
			    Shortread_T queryseq, int querylength,
			    bool supplemental_hardclip_low_p, bool supplemental_hardclip_high_p) {
  Intlist_T q;
  List_T j;
  Univcoordlist_T u;
  Univcoord_T deletionpos;
  Junction_T junction;
  Junctiontype_T type;
  char *deletion_string;
  int qpos, nconsecutive, n, i;
  int ninserts, softclip, choplength;
  char *p, c;
  bool md_startp = true;

  debug(printf("Entered Path_supplemental_cigar_md\n"));

  *cigar_fp = Filestring_new();
  *md_fp = Filestring_new();

  ninserts = 0;
  nconsecutive = 0;

  if (segment_plusp != query_plusp) {
    p = &(genomic_diff[querylength-1]);
  } else {
    p = genomic_diff;
  }
  q = suppl_endpoints;
  j = suppl_junctions;
  u = suppl_univdiagonals;

  /* Compute softclip from chop */
  if (segment_plusp == true) {
    choplength = Shortread_left_choplength(queryseq);
  } else {
    choplength = Shortread_right_choplength(queryseq);
  }	

  /* Compute softclip from alignment and hardclip */
  qpos = Intlist_head(q);
  if (supplemental_hardclip_low_p == true) {
    *hardclip_low = choplength + qpos;
    softclip = 0;

  } else {
    *hardclip_low = 0;
    softclip = choplength + qpos;
  }
  /* Ignore choplength, which is not included in genomic_diff */
  p += (segment_plusp == query_plusp) ? +qpos : -qpos;

  if (softclip > 0) {
    FPRINTF(*cigar_fp,"%dS",softclip);
  }
    
  if (*hardclip_low > 0) {
    FPRINTF(*cigar_fp,"%dH",*hardclip_low);
  }

  while (j != NULL) {
    q = Intlist_next(q);
    n = Intlist_head(q) - qpos;
    if (n > 0 || sam_insert_0M_p == true) {
      FPRINTF(*cigar_fp,"%dM",n);
    }
    for (i = 0; i < n; i++) {
      c = (segment_plusp != query_plusp) ? complCode[(int) *p--] : *p++;
      if (isupper(c)) {
	/* Match */
	nconsecutive++;

      } else if (nconsecutive > 0) {
	/* Mismatch */
	FPRINTF(*md_fp,"%d%c",nconsecutive,toupper(c));
	nconsecutive = 0;
	md_startp = false;
	  
      } else if (md_startp == false) {
	/* Consecutive mismatches */
	FPRINTF(*md_fp,"%c",toupper(c));
	  
      } else {
	/* Initial mismatch: MD string needs to start with a number */
	FPRINTF(*md_fp,"0%c",toupper(c));
	md_startp = false;
      }
    }

    ninserts = 0;
    junction = (Junction_T) List_head(j);
    if (junction == JUNCTION_UNSOLVED) {
      FPRINTF(*cigar_fp,"X");

    } else if ((type = Junction_type(junction)) == DEL_JUNCTION) {
      FPRINTF(*cigar_fp,"%dD",Junction_nindels(junction));

      deletionpos = (Univcoordlist_head(u) - querylength) + /*qpos*/Intlist_head(q);
      deletion_string = Junction_deletion_string(junction,deletionpos);
      FPRINTF(*md_fp,"%d^%s",nconsecutive,deletion_string);
      FREE(deletion_string);

      nconsecutive = 0;
      md_startp = false;

    } else if (type == INS_JUNCTION) {
      FPRINTF(*cigar_fp,"%dI",ninserts = Junction_nindels(junction));
      p += ninserts;

    } else if (type == SPLICE_JUNCTION) {
      FPRINTF(*cigar_fp,"%dN",Junction_splice_distance(junction));

    } else {
      fprintf(stderr,"Unknown junction type %d\n",type);
      abort();
    }
    qpos = Intlist_head(q) + ninserts;

    j = List_next(j);
    u = Univcoordlist_next(u);
  }

  q = Intlist_next(q);
  n = Intlist_head(q) - qpos;
  if (n > 0 || sam_insert_0M_p == true) {
    FPRINTF(*cigar_fp,"%dM",n);
  }
  for (i = 0; i < n; i++) {
    c = (segment_plusp != query_plusp) ? complCode[(int) *p--] : *p++;
    if (isupper(c)) {
      /* Match */
      nconsecutive++;

    } else if (nconsecutive > 0) {
      /* Mismatch */
      FPRINTF(*md_fp,"%d%c",nconsecutive,toupper(c));
      nconsecutive = 0;
      md_startp = false;
	  
    } else if (md_startp == false) {
      /* Consecutive mismatches */
      FPRINTF(*md_fp,"%c",toupper(c));
	  
    } else {
      /* Initial mismatch: MD string needs to start with a number */
      FPRINTF(*md_fp,"0%c",toupper(c));
      md_startp = false;
    }
  }

  /* Previously checked for nconsecutive > 0, but MD string needs to end with a number */
  FPRINTF(*md_fp,"%d",nconsecutive);

  /* Compute softclip from chop */
  if (segment_plusp == true) {
    choplength = Shortread_right_choplength(queryseq);
  } else {
    choplength = Shortread_left_choplength(queryseq);
  }	

  /* Compute softclip from alignment and hardclip */
  qpos = Intlist_head(q);
  if (supplemental_hardclip_high_p == true) {
    *hardclip_high = (querylength - qpos) + choplength;
    softclip = 0;
    
  } else {
    *hardclip_high = 0;
    softclip = (querylength - qpos) + choplength;
  }

  if (*hardclip_high > 0) {
    FPRINTF(*cigar_fp,"%dH",*hardclip_high);
  }

  if (softclip > 0) {
    FPRINTF(*cigar_fp,"%dS",softclip);
  }
    
  /* Filestring_stringify required before Filestring_merge */
  Filestring_stringify(*cigar_fp);
  Filestring_stringify(*md_fp);

  return;
}


/* npaths could be non-zero, if user selected --quiet-if-excessive */
void
Path_print_sam_nomapping (Filestring_T fp, char *abbrev, Shortread_T queryseq, Shortread_T mate_queryseq,
			  Shortread_T single_cell_infoseq,
			  char *acc1, char *acc2, Univ_IIT_T chromosome_iit, Resulttype_T resulttype, bool first_read_p,
			  int pathnum, int npaths_primary, int npaths_altloc, bool artificial_mate_p, int npaths_mate,

			  T mate, int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag;

  Univcoord_T mate_chrpos_low;
  Filestring_T mate_cigar_fp, mate_md_fp;
  int mate_hardclip_low, mate_hardclip_high;

  char *chr;
  bool allocp;


  /* If mate is non-NULL, then paired_read_p must be true, which means
     that Path_print_sam_paired has taken care of calling Path_trim_circular on mate */

  /* mate_first_read_p = (first_read_p == true ? false : true); */
  debug(printf("Calling Cigar_compute_main on mate %p\n",mate));
  Path_cigar_md(&mate_hardclip_low,&mate_hardclip_high,&mate_cigar_fp,&mate_md_fp,mate,
		mate_queryseq);
  debug(printf("Computed mate hardclips %d and %d\n",mate_hardclip_low,mate_hardclip_high));

  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }
  
  /* 2. FLAG */
  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  /* 5. MAPQ: Mapping quality.  Picard says MAPQ should be 0 for an unmapped read */
  /* 6. CIGAR */
  flag = compute_flag(/*plusp (NA)*/true,mate,resulttype,first_read_p,
		      /*pathnum*/0,/*npaths*/0,artificial_mate_p,npaths_mate,
		      /*absmq_score*/0,/*first_absmq*/0,invertp,invert_mate_p,
		      /*supplementaryp*/false);
  FPRINTF(fp,"\t%u\t*\t0\t0\t*",flag);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t*\t0");
  } else if ((mate_chrpos_low = Path_genomiclow_softclipped(mate)) == 0U) {
    FPRINTF(fp,"\t*\t0");
  } else {
    chr = Univ_IIT_label(chromosome_iit,mate->chrnum,&allocp);
    FPRINTF(fp,"\t%s\t%u",chr,mate_chrpos_low - mate->chroffset + 1U);
    if (allocp == true) {
      FREE(chr);
    }
  }


  /* 9. ISIZE: Insert size */
  FPRINTF(fp,"\t0");

  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Since there is no mapping, we print the original query sequence. */
  if (sam_sparse_secondaries_p == true && (flag & NOT_PRIMARY) != 0) {
    /* SAM format specification says that secondary mappings should not print SEQ or QUAL to reduce file size */
    FPRINTF(fp,"\t*\t*");

  } else if (invertp == false) {
#if 0
    /* Intended for primers, not poly-A/T */
    Shortread_print_chopped_sam(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0);
    Shortread_print_quality(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0,
			    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_sam(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0);
    Shortread_print_hardclipped_quality(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0,
					quality_shift);
#endif
  } else {
#if 0
    /* Intended for primers, not poly-A/T */
    Shortread_print_chopped_revcomp_sam(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0);
    Shortread_print_quality_revcomp(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0,
				    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_revcomp_sam(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0);
    Shortread_print_hardclipped_reverse_quality(fp,queryseq,/*sequence_hardclip_low*/0,/*sequence_hardclip_high*/0,
						quality_shift);
#endif
  }

  /* 12. TAGS: XM (mate CIGAR), XD (mate MD string), and XN (mate mismatches) */
  if (mate == NULL || (flag & MATE_UNMAPPED) != 0) {
    /* Previously checked if queryseq_mate == NULL */
    /* Unpaired alignment.  Don't print XM. */
  } else {
    FPRINTF(fp,"\tXM:Z:");
    Filestring_merge(fp,mate_cigar_fp);
    FPRINTF(fp,"\tXD:Z:");
    Filestring_merge(fp,mate_md_fp);
    FPRINTF(fp,"\tXN:i:%d",Path_ndiffs(mate));
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }
  
  /* 12. TAGS: NH */
  if (npaths_primary + npaths_altloc > 0) {
    FPRINTF(fp,"\tNH:i:%d",npaths_primary + npaths_altloc);
    if (add_paired_nomappers_p == true) {
      FPRINTF(fp,"\tHI:i:%d",pathnum);
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

#if 0
  /* 12. TAGS: XP */
  /* Intended for primers, not poly-A/T */
  Shortread_print_chop(fp,queryseq,invertp);
#endif

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: (BC,) CB, CR, CY, UR, UY */
  if (single_cell_infoseq != NULL) {
    Single_cell_print_fields(fp,single_cell_infoseq);
  }

  FPRINTF(fp,"\n");

  Filestring_free(&mate_md_fp,/*free_string_p*/true);
  Filestring_free(&mate_cigar_fp,/*free_string_p*/true);

  return;
}


static void
print_substrings (int *hardclip_low, int *hardclip_high,
		  Filestring_T fp, char *abbrev, T path,
		  
		  char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
		  int absmq_score, int first_absmq, int second_absmq, int mapq_score,
		  Shortread_T queryseq, Shortread_T mate_queryseq, Shortread_T single_cell_infoseq,
		  int pairedlength, int pair_relationship,

		  T mate,
		  Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		  int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  int mate_hardclip_low, mate_hardclip_high;

  Univcoord_T univcoord_low, mate_univcoord_low;
  Filestring_T cigar_fp, mate_cigar_fp, md_fp, mate_md_fp;

  char *chr;
  bool allocp;

  int sensedir;
  int transcript_genestrand;

  int n, i;
  Univcoord_T *alts_univdiagonals, medial_univdiagonal;
#ifdef PRINT_ALTS_COORDS
  Univcoord_T chroffset;
#endif

  
  /* Compute on mate first, so we don't use its MD pointers */
  /* mate_first_read_p = (first_read_p == true ? false : true); */
  debug(printf("print_substrings calling Path_cigar_md on mate %p\n",mate));
  Path_cigar_md(&mate_hardclip_low,&mate_hardclip_high,&mate_cigar_fp,&mate_md_fp,mate,
		mate_queryseq);
  debug(printf("computed mate hardclips %d and %d\n",mate_hardclip_low,mate_hardclip_high));

  debug(printf("print_substrings calling Path_cigar_md on main part\n"));
  Path_cigar_md(&(*hardclip_low),&(*hardclip_high),&cigar_fp,&md_fp,path,queryseq);


  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = compute_flag(path->plusp,mate,resulttype,first_read_p,
		      pathnum,npaths_primary + npaths_altloc,artificial_mate_p,npaths_mate,
		      absmq_score,first_absmq,invertp,invert_mate_p,/*supplementaryp*/false);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chr = Univ_IIT_label(chromosome_iit,path->chrnum,&allocp);
  univcoord_low = Path_genomiclow_softclipped(path);
  FPRINTF(fp,"\t%s\t%u",chr,univcoord_low - path->chroffset + 1U);
  assert(univcoord_low >= path->chroffset);
  if (allocp == true) {
    FREE(chr);
  }

  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d\t",mapq_score);

  /* 6. CIGAR */
  Filestring_merge(fp,cigar_fp);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t*\t0");
  } else {
    mate_univcoord_low = Path_genomiclow_softclipped(mate);
    if (mate->chrnum == path->chrnum) {
      FPRINTF(fp,"\t=\t%u",mate_univcoord_low - mate->chroffset + 1U);
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate->chrnum,&allocp);
      FPRINTF(fp,"\t%s\t%u",chr,mate_univcoord_low - mate->chroffset + 1U);
      if (allocp == true) {
	FREE(chr);
      }
    }
  }

  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (pair_relationship > 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",pairedlength);
      } else {
	FPRINTF(fp,"\t%d",-pairedlength);
      }

    } else if (pair_relationship < 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",-pairedlength);
      } else {
	FPRINTF(fp,"\t%d",pairedlength);
      }

    } else if (path->plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }

  } else if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t0");
  } else if (univcoord_low < mate_univcoord_low) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (univcoord_low > mate_univcoord_low) {
    FPRINTF(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  if (sam_sparse_secondaries_p == true && (flag & NOT_PRIMARY) != 0) {
    /* SAM format specification says that secondary mappings should not print SEQ or QUAL to reduce file size */
    FPRINTF(fp,"\t*\t*");

  } else if (path->plusp == true) {
#if 0
    /* Intended for primers, not poly-A/T */
    Shortread_print_chopped_sam(fp,queryseq,*hardclip_low,*hardclip_high);
    Shortread_print_quality(fp,queryseq,*hardclip_low,*hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_sam(fp,queryseq,*hardclip_low,*hardclip_high);
    Shortread_print_hardclipped_quality(fp,queryseq,*hardclip_low,*hardclip_high,
					quality_shift);
#endif
  } else {
#if 0
    /* Intended for primers, not poly-A/T */
    Shortread_print_chopped_revcomp_sam(fp,queryseq,*hardclip_low,*hardclip_high);
    Shortread_print_quality_revcomp(fp,queryseq,*hardclip_low,*hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_revcomp_sam(fp,queryseq,*hardclip_low,*hardclip_high);
    Shortread_print_hardclipped_reverse_quality(fp,queryseq,*hardclip_low,*hardclip_high,
						quality_shift);
#endif
  }


  /* 12. TAGS: XM (mate CIGAR), XD (mate MD string), and XN (mate mismatches) */
  if (mate == NULL || (flag & MATE_UNMAPPED) != 0) {
    /* Previously checked if queryseq_mate == NULL */
    /* Unpaired alignment.  Don't print XM. */
  } else {
    FPRINTF(fp,"\tXM:Z:");
    Filestring_merge(fp,mate_cigar_fp);
    FPRINTF(fp,"\tXD:Z:");
    Filestring_merge(fp,mate_md_fp);
    FPRINTF(fp,"\tXN:i:%d",Path_ndiffs(mate));
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI */
  if (*hardclip_low > 0 || *hardclip_high > 0) {
    FPRINTF(fp,"\tXH:Z:");
    if (path->plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,*hardclip_low,*hardclip_high);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,*hardclip_low,*hardclip_high);
    }

    if (Shortread_quality_string(queryseq) != NULL) {
      FPRINTF(fp,"\tXI:Z:");
      if (path->plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,*hardclip_low,*hardclip_high,quality_shift);
      } else {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,*hardclip_low,*hardclip_high,quality_shift);
      }
    }
  }

  /* 12. TAGS: XB */
  Shortread_print_barcode(fp,queryseq);

#if 0
  /* 12. TAGS: XP.  Logically should be last in reconstructing a read. */
  /* Intended for primers, not poly-A/T */
  Shortread_print_chop(fp,queryseq,invertp);
#endif

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  Filestring_merge(fp,md_fp);

  /* 12. TAGS: NH */
  /* 12. TAGS: HI */
  FPRINTF(fp,"\tNH:i:%d\tHI:i:%d",npaths_primary + npaths_altloc,pathnum);

  /* 12. TAGS: NM (mismatches) */
  FPRINTF(fp,"\tNM:i:%d",Path_ndiffs(path));
  
  /* 12. TAGS: XE */
#ifdef TO_FIX
  if (maskedp) {
    FPRINTF(fp,"\tXE:i:%d",nmatches_exonic);
  }
#endif

#ifdef TO_FIX
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }
#endif

  /* 12. TAGS: SM */
  /* 12. TAGS: XQ */
  /* 12. TAGS: X2 */
  FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq);

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XA */
  if (path->qstart_alts != NULL || path->qend_alts != NULL) {
    FPRINTF(fp,"\tXA:Z:");

    if (path->qstart_alts != NULL) {
      alts_univdiagonals = path->qstart_alts->univdiagonals;
      n = path->qstart_alts->nunivdiagonals;
#ifdef PRINT_ALTS_COORDS
      FPRINTF(fp,"%u",alts_univdiagonals[0] - path->querylength + path->qstart_alts->splice_qpos - path->chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_univdiagonals[i] - path->querylength + path->qstart_alts->splice_qpos - path->chroffset + 1U);
      }
#else
      medial_univdiagonal = Univcoordlist_head(path->univdiagonals);
      FPRINTF(fp,"%u",medial_univdiagonal - alts_univdiagonals[0]);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",medial_univdiagonal - alts_univdiagonals[i]);
      }
#endif
    }
    FPRINTF(fp,"|");
    if (path->qend_alts != NULL) {
      alts_univdiagonals = path->qend_alts->univdiagonals;
      n = path->qend_alts->nunivdiagonals;
#ifdef PRINT_ALTS_COORDS
      FPRINTF(fp,"%u",alts_univdiagonals[0] - path->querylength + path->qend_alts->splice_qpos - path->chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_univdiagonals[i] - path->querylength + path->qend_alts->splice_qpos - path->chroffset + 1U);
      }
#else
      medial_univdiagonal = Univcoordlist_last_value(path->univdiagonals);
      FPRINTF(fp,"%u",alts_univdiagonals[0] - medial_univdiagonal);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_univdiagonals[i] - medial_univdiagonal);
      }
#endif
    }
  }

#if 0
  /* 12. TAG: XT (possible with transcriptome-guided genomic alignment) */
  /* Now printing only on the supplemental alignment */
  if ((donor = Stage3end_substring_donor(stage3end)) != NULL &&
      (acceptor = Stage3end_substring_acceptor(stage3end)) != NULL) {
    print_xt_info(fp,donor,acceptor,Stage3end_sensedir(stage3end));
  }
#endif

  /* 12. TAGS: XX, XY */
  if (path->transcripts != NULL) {
    Transcript_print_list(fp,path->transcripts,transcript_iit,/*header*/"\tXX:Z:");
  } else if (path->invalid_transcripts != NULL) {
    Transcript_print_list(fp,path->invalid_transcripts,transcript_iit,/*header*/"\tXY:Z:");
  }

  /* 12. TAGS: XS */
  if ((sensedir = Path_effective_sensedir(path)) == SENSE_NULL) {
    transcript_genestrand = 0;

  } else if (sensedir == SENSE_FORWARD) {
    if (path->plusp == true) {
      transcript_genestrand = +1;
    } else {
      transcript_genestrand = -1;
    }

  } else if (sensedir == SENSE_ANTI) {
    if (path->plusp == true) {
      transcript_genestrand = -1;
    } else {
      transcript_genestrand = +1;
    }

  }

  if (transcript_genestrand == 0 && mate != NULL &&
      (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT)) {
    /* Use mate effective sensedir */
    sensedir = Path_effective_sensedir(mate);
    if (sensedir == SENSE_FORWARD) {
      if (path->plusp == true) {
	transcript_genestrand = +1;
      } else {
	transcript_genestrand = -1;
      }
    } else if (sensedir == SENSE_ANTI) {
      if (path->plusp == true) {
	transcript_genestrand = -1;
      } else {
	transcript_genestrand = +1;
      }
    } else {
      transcript_genestrand = 0;
    }
  }
    
  if (transcript_genestrand > 0) {
    FPRINTF(fp,"\tXS:A:+");
  } else if (transcript_genestrand < 0) {
    FPRINTF(fp,"\tXS:A:-");
  }


#if 0
  /* 12. TAGS: XC */
  if (circularp == true) {
    FPRINTF(fp,"\tXC:A:+");
  }
#endif

  /* 12. TAGS: XG */
  if (method_print_p == true) {
    Method_samprint(fp,path->method);
  }

#if 0
  /* 12. TAGS: XE (BLAST E-value) */
  FPRINTF(fp,"\tXE:f:%.2g",Stage3end_min_evalue(stage3end));
#endif


  /* 12. TAGS: (BC,) CB, CR, CY, UR, UY */
  if (single_cell_infoseq != NULL) {
    Single_cell_print_fields(fp,single_cell_infoseq);
  }

  FPRINTF(fp,"\n");

  Filestring_free(&mate_md_fp,/*free_string_p*/true);
  Filestring_free(&md_fp,/*free_string_p*/true);
  Filestring_free(&mate_cigar_fp,/*free_string_p*/true);
  Filestring_free(&cigar_fp,/*free_string_p*/true);

  return;
}


static void
print_circular (Filestring_T fp, char *abbrev, T path,
		char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
		int absmq_score, int first_absmq, int mapq_score,
		Shortread_T queryseq, Shortread_T mate_queryseq, Shortread_T single_cell_infoseq,
		int pairedlength, int pair_relationship,
		T mate, Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {

  unsigned int flag = 0U;

  int hardclip_low, hardclip_high, mate_hardclip_low, mate_hardclip_high;
  bool supplemental_hardclip_low_p, supplemental_hardclip_high_p;

  Univcoord_T univcoord_low, mate_univcoord_low;
  Filestring_T cigar_fp, mate_cigar_fp, md_fp, mate_md_fp;

  char *chr;
  bool allocp;

#ifdef PRINT_ALTS_COORDS
  Univcoord_T chroffset;
#endif

  
  /* Compute on mate first, so we don't use its MD pointers */
  debug(printf("Calling Cigar_compute_main on mate %p\n",mate));
  /* mate_first_read_p = (first_read_p == true ? false : true); */
  Path_cigar_md(&mate_hardclip_low,&mate_hardclip_high,&mate_cigar_fp,&mate_md_fp,mate,
		mate_queryseq);
  debug(printf("computed mate hardclips %d and %d\n",mate_hardclip_low,mate_hardclip_high));

  /* Side of the supplemental hardclip is the opposite of that of the main path */
  if (path->circular_high_p == true) {
    supplemental_hardclip_low_p = true;
    supplemental_hardclip_high_p = false;
  } else {
    supplemental_hardclip_low_p = false;
    supplemental_hardclip_high_p = true;
  }

  /* path->querylength does not include chopped sequence, but queryseq->querylength does */
  Path_supplemental_cigar_md(&hardclip_low,&hardclip_high,&cigar_fp,&md_fp,
			     path->genomic_diff,path->circular_endpoints,
			     path->circular_junctions,path->circular_univdiagonals,
			     /*segment_plusp*/path->plusp,/*query_plusp*/path->plusp,
			     queryseq,path->querylength,
			     supplemental_hardclip_low_p,supplemental_hardclip_high_p);


  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = compute_flag(path->plusp,mate,resulttype,first_read_p,
		      pathnum,npaths_primary + npaths_altloc,artificial_mate_p,npaths_mate,
		      absmq_score,first_absmq,invertp,invert_mate_p,/*supplementaryp*/true);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chr = Univ_IIT_label(chromosome_iit,path->chrnum,&allocp);
  univcoord_low = Path_genomiclow_circular_softclipped(path);
  assert(univcoord_low >= path->chroffset);
  FPRINTF(fp,"\t%s\t%u",chr,univcoord_low - path->chroffset + 1U);
  if (allocp == true) {
    FREE(chr);
  }

  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d\t",mapq_score);

  /* 6. CIGAR */
  Filestring_merge(fp,cigar_fp);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t*\t0");
  } else {
    mate_univcoord_low = Path_genomiclow_softclipped(mate);
    if (mate->chrnum == path->chrnum) {
      FPRINTF(fp,"\t=\t%u",mate_univcoord_low - mate->chroffset + 1U);
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate->chrnum,&allocp);
      FPRINTF(fp,"\t%s\t%u",chr,mate_univcoord_low - mate->chroffset + 1U);
      if (allocp == true) {
	FREE(chr);
      }
    }
  }

  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (pair_relationship > 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",pairedlength);
      } else {
	FPRINTF(fp,"\t%d",-pairedlength);
      }

    } else if (pair_relationship < 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",-pairedlength);
      } else {
	FPRINTF(fp,"\t%d",pairedlength);
      }

    } else if (path->plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }

  } else if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t0");
  } else if (univcoord_low < mate_univcoord_low) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (univcoord_low > mate_univcoord_low) {
    FPRINTF(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  /* Same as XH and XI fields in primary mapping, except on genomic plus strand */
  if (sam_sparse_secondaries_p == true && (flag & NOT_PRIMARY) != 0) {
    /* SAM format specification says that secondary mappings should not print SEQ or QUAL to reduce file size */
    /* We can use the XH field in the primary mapping to reconstruct the original sequence */
    FPRINTF(fp,"\t*\t*");

  } else if (path->plusp == true) {
#if 0
    Shortread_print_chopped_sam(fp,queryseq,hardclip_low,hardclip_high);
    Shortread_print_quality(fp,queryseq,hardclip_low,hardclip_high,
			    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_sam(fp,queryseq,hardclip_low,hardclip_high);
    Shortread_print_hardclipped_quality(fp,queryseq,hardclip_low,hardclip_high,
					quality_shift);
#endif
  } else {
#if 0
    Shortread_print_chopped_revcomp_sam(fp,queryseq,hardclip_low,hardclip_high);
    Shortread_print_quality_revcomp(fp,queryseq,hardclip_low,hardclip_high,
				    quality_shift,/*show_chopped_p*/false);
#else
    Shortread_print_hardclipped_revcomp_sam(fp,queryseq,hardclip_low,hardclip_high);
    Shortread_print_hardclipped_reverse_quality(fp,queryseq,hardclip_low,hardclip_high,
						quality_shift);
#endif
  }


  /* 12. TAGS: XM (mate CIGAR), XD (mate MD string), and XN (mate mismatches) */
  if (mate == NULL || (flag & MATE_UNMAPPED) != 0) {
    /* Previously checked if queryseq_mate == NULL */
    /* Unpaired alignment.  Don't print XM. */
  } else {
    FPRINTF(fp,"\tXM:Z:");
    Filestring_merge(fp,mate_cigar_fp);
    FPRINTF(fp,"\tXD:Z:");
    Filestring_merge(fp,mate_md_fp);
    FPRINTF(fp,"\tXN:i:%d",Path_ndiffs(mate));
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI.  Not in supplemental. */

  /* 12. TAGS: XB.  Not in supplemental */
  /* Shortread_print_barcode(fp,queryseq); */

  /* 12. TAGS: XP.  Not in supplemental */
  /* Shortread_print_chop(fp,queryseq,invertp); */

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  Filestring_merge(fp,md_fp);

  /* 12. TAGS: NH.  Not in supplemental */
  /* 12. TAGS: HI */
  FPRINTF(fp,"\tHI:i:%d",npaths_primary + npaths_altloc,pathnum);

  /* 12. TAGS: NM (mismatches) */
  FPRINTF(fp,"\tNM:i:%d",Path_ndiffs(path));
  
  /* 12. TAGS: XE */
#ifdef TO_FIX
  if (maskedp) {
    FPRINTF(fp,"\tXE:i:%d",nmatches_exonic);
  }
#endif

#ifdef TO_FIX
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }
#endif

  /* 12. TAGS: SM, XQ, X2.  Not in supplemental */
  /* FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq); */

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: (BC,) CB, CR, CY, UR, UY */
  if (single_cell_infoseq != NULL) {
    Single_cell_print_fields(fp,single_cell_infoseq);
  }

  FPRINTF(fp,"\n");

  Filestring_free(&mate_md_fp,/*free_string_p*/true);
  Filestring_free(&md_fp,/*free_string_p*/true);
  Filestring_free(&mate_cigar_fp,/*free_string_p*/true);
  Filestring_free(&cigar_fp,/*free_string_p*/true);

  return;
}


static void
print_fusion (Filestring_T fp, char *abbrev, T path,
	      char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
	      int absmq_score, int first_absmq, int mapq_score,
	      Shortread_T queryseq, Shortread_T mate_queryseq, Shortread_T single_cell_infoseq,
	      int pairedlength, int pair_relationship, int hardclip_low, int hardclip_high,

	      T mate, Resulttype_T resulttype, bool first_read_p, bool artificial_mate_p, int npaths_mate,
	      int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p) {
  unsigned int flag = 0U;

  int hardclip_low_ignore, hardclip_high_ignore;
  int mate_hardclip_low, mate_hardclip_high;
  bool supplemental_hardclip_low_p, supplemental_hardclip_high_p;

  Univcoord_T univcoord_low, mate_univcoord_low;
  Filestring_T cigar_fp, mate_cigar_fp, md_fp, mate_md_fp;

  Junction_T fusion_junction;
  char *chr, *donor_chr, *acceptor_chr, donor_strand, acceptor_strand;
  Univcoord_T donor_position, acceptor_position;
  Chrpos_T donor_chrpos, acceptor_chrpos;
  bool allocp, alloc1p, alloc2p;

#ifdef PRINT_ALTS_COORDS
  Univcoord_T chroffset;
#endif

  debug(printf("Entered print_fusion with hardclip_low %d, hardclip_high %d\n",
	       hardclip_low,hardclip_high));
  
  /* Compute on mate first, so we don't use its MD pointers */
  debug(printf("print_fusion calling Path_cigar_md on mate %p with mate hardclips %d and %d\n",
	       mate,mate_hardclip_low,mate_hardclip_high));
  /* mate_first_read_p = (first_read_p == true ? false : true); */
  Path_cigar_md(&mate_hardclip_low,&mate_hardclip_high,&mate_cigar_fp,&mate_md_fp,mate,
		mate_queryseq);

  /* Rules for supplemental hardclip sides are the same of those
     for the main path, but use fusion_plusp instead of main path
     plusp */
  debug(printf("fusion_plusp %d, fusion_querystart_junction %p, fusion_queryend_junction %p\n",
	       path->fusion_plusp,path->fusion_querystart_junction,path->fusion_queryend_junction));
  if ((path->fusion_plusp == true && path->fusion_querystart_junction != NULL) ||
      (path->fusion_plusp == false && path->fusion_queryend_junction != NULL)) {
    debug(printf("supplemental hardclip is on high side\n"));
    supplemental_hardclip_low_p = false;
    supplemental_hardclip_high_p = true;
  } else {
    debug(printf("supplemental hardclip is on low side\n"));
    supplemental_hardclip_low_p = true;
    supplemental_hardclip_high_p = false;
  }

  debug(printf("print_fusion calling Path_supplemental_cigar_md\n"));
  Path_supplemental_cigar_md(&hardclip_low_ignore,&hardclip_high_ignore,&cigar_fp,&md_fp,
			     path->genomic_diff,path->fusion_endpoints,
			     path->fusion_junctions,path->fusion_univdiagonals,
			     /*segment_plusp*/path->fusion_plusp,/*query_plusp*/path->plusp,queryseq,
			     path->querylength,supplemental_hardclip_low_p,supplemental_hardclip_high_p);

  /* Printing the supplemental (fusion) part */
  /* 1. QNAME */
  if (acc2 == NULL) {
    FPRINTF(fp,"%s",acc1);
  } else {
    FPRINTF(fp,"%s,%s",acc1,acc2);
  }

  /* 2. FLAG */
  flag = compute_flag(path->fusion_plusp,mate,resulttype,first_read_p,
		      pathnum,npaths_primary + npaths_altloc,artificial_mate_p,npaths_mate,
		      absmq_score,first_absmq,invertp,invert_mate_p,/*supplementaryp*/true);
  FPRINTF(fp,"\t%u",flag);

  /* 3. RNAME: chr */
  /* 4. POS: chrpos */
  chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&allocp);
  univcoord_low = Path_genomiclow_fusion_softclipped(path);
  FPRINTF(fp,"\t%s\t%u",chr,univcoord_low - path->fusion_chroffset + 1U);
  assert(univcoord_low >= path->fusion_chroffset);
  if (allocp == true) {
    FREE(chr);
  }

  /* 5. MAPQ: Mapping quality */
  FPRINTF(fp,"\t%d\t",mapq_score);

  /* 6. CIGAR */
  Filestring_merge(fp,cigar_fp);

  /* 7. MRNM: Mate chr */
  /* 8. MPOS: Mate chrpos */
  if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t*\t0");
  } else {
    mate_univcoord_low = Path_genomiclow_softclipped(mate);
    if (mate->chrnum == path->fusion_chrnum) {
      FPRINTF(fp,"\t=\t%u",mate_univcoord_low - mate->chroffset + 1U);
    } else {
      chr = Univ_IIT_label(chromosome_iit,mate->chrnum,&allocp);
      FPRINTF(fp,"\t%s\t%u",chr,mate_univcoord_low - mate->chroffset + 1U);
      if (allocp == true) {
	FREE(chr);
      }
    }
  }

  /* 9. ISIZE: Insert size */
  if (resulttype == CONCORDANT_UNIQ || resulttype == CONCORDANT_TRANSLOC || resulttype == CONCORDANT_MULT) {
    if (pair_relationship > 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",pairedlength);
      } else {
	FPRINTF(fp,"\t%d",-pairedlength);
      }

    } else if (pair_relationship < 0) {
      if (first_read_p == true) {
	FPRINTF(fp,"\t%d",-pairedlength);
      } else {
	FPRINTF(fp,"\t%d",pairedlength);
      }

    } else if (path->plusp == invertp) {
      FPRINTF(fp,"\t%d",-pairedlength);
    } else {
      FPRINTF(fp,"\t%d",pairedlength);
    }

  } else if (mate == (Path_T) NULL) {
    FPRINTF(fp,"\t0");
  } else if (univcoord_low < mate_univcoord_low) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else if (univcoord_low > mate_univcoord_low) {
    FPRINTF(fp,"\t%d",-pairedlength);
  } else if (first_read_p == true) {
    FPRINTF(fp,"\t%d",pairedlength);
  } else {
    FPRINTF(fp,"\t%d",-pairedlength);
  }


  /* 10. SEQ: queryseq and 11. QUAL: quality scores */
  /* Queryseq has already been inverted, so just measure plusp relative to its current state */
  /* Same as XH and XI fields in primary mapping, except on genomic plus strand */
  if (sam_sparse_secondaries_p == true && (flag & NOT_PRIMARY) != 0) {
    /* SAM format specification says that secondary mappings should not print SEQ or QUAL to reduce file size */
    /* We can use the XH field in the primary mapping to reconstruct the original sequence */
    FPRINTF(fp,"\t*\t*");

  } else {
    FPRINTF(fp,"\t");
    if (path->plusp == true && path->fusion_plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_low,hardclip_high);
    } else if (path->plusp == true && path->fusion_plusp == false) {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_high,hardclip_low);
    } else if (path->plusp == false && path->fusion_plusp == true) {
      Shortread_print_chopped_end(fp,queryseq,hardclip_high,hardclip_low);
    } else {
      Shortread_print_chopped_end_revcomp(fp,queryseq,hardclip_low,hardclip_high);
    }

    if (Shortread_quality_string(queryseq) == NULL) {
      FPRINTF(fp,"\t*");
    } else {
      FPRINTF(fp,"\t");
      if (path->plusp == true && path->fusion_plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,hardclip_low,hardclip_high,quality_shift);
      } else if (path->plusp == true && path->fusion_plusp == false) {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,hardclip_high,hardclip_low,quality_shift);
      } else if (path->plusp == false && path->fusion_plusp == true) {
	Shortread_print_chopped_end_quality(fp,queryseq,hardclip_high,hardclip_low,quality_shift);
      } else {
	Shortread_print_chopped_end_quality_reverse(fp,queryseq,hardclip_low,hardclip_high,quality_shift);
      }
    }
  }


  /* 12. TAGS: XM (mate CIGAR), XD (mate MD string), and XN (mate mismatches) */
  if (mate == NULL || (flag & MATE_UNMAPPED) != 0) {
    /* Previously checked if queryseq_mate == NULL */
    /* Unpaired alignment.  Don't print XM. */
  } else {
    FPRINTF(fp,"\tXM:Z:");
    Filestring_merge(fp,mate_cigar_fp);
    FPRINTF(fp,"\tXD:Z:");
    Filestring_merge(fp,mate_md_fp);
    FPRINTF(fp,"\tXN:i:%d",Path_ndiffs(mate));
  }

  /* 12. TAGS: RG */
  if (sam_read_group_id != NULL) {
    FPRINTF(fp,"\tRG:Z:%s",sam_read_group_id);
  }

  /* 12. TAGS: XH and XI.  Not in supplemental */

  /* 12. TAGS: XB.  Not in supplemental */
  /* Shortread_print_barcode(fp,queryseq); */

  /* 12. TAGS: XP.  Not in supplemental */
  /* Shortread_print_chop(fp,queryseq,invertp); */

  /* 12. TAGS: MD */
  FPRINTF(fp,"\tMD:Z:");
  Filestring_merge(fp,md_fp);

  /* 12. TAGS: NH.  Not in supplemental */
  /* 12. TAGS: HI */
  FPRINTF(fp,"\tHI:i:%d",npaths_primary + npaths_altloc,pathnum);

  /* 12. TAGS: NM (mismatches) */
  FPRINTF(fp,"\tNM:i:%d",Path_ndiffs(path));
  
  /* 12. TAGS: XE */
#ifdef TO_FIX
  if (maskedp) {
    FPRINTF(fp,"\tXE:i:%d",nmatches_exonic);
  }
#endif

#ifdef TO_FIX
  if (snps_iit) {
    /* 12. TAGS: XW and XV */
    FPRINTF(fp,"\tXW:i:%d",nmismatches_bothdiff);
    FPRINTF(fp,"\tXV:i:%d",nmismatches_refdiff - nmismatches_bothdiff);
  }
#endif

  /* 12. TAGS: SM, XQ, X2.  Not in supplemental */
  /* FPRINTF(fp,"\tSM:i:%d\tXQ:i:%d\tX2:i:%d",mapq_score,absmq_score,second_absmq); */

  /* 12. TAGS: XO */
  FPRINTF(fp,"\tXO:Z:%s",abbrev);

  /* 12. TAGS: XT */
  /* Add 1 to calls to Intlist_head, but not to Intlist_last_value */
  debug2(Path_print(path));
  if (path->fusion_querystart_junction != NULL) {
    fusion_junction = path->fusion_querystart_junction;

    if (path->plusp == true) {
      /* Case 1, 2 */
      /* 5' is fusion; 3' is main (plus strand) */
      if (path->sensedir == SENSE_FORWARD) {
	/* 5' (fusion) is donor; 3' (main) is acceptor */
	if (path->fusion_plusp == true) {
	  debug2(printf("Case 1, sense forward\n"));
	  donor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	} else {
	  debug2(printf("Case 2, sense forward\n"));
	  donor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	}
	acceptor_position = Univcoordlist_head(path->univdiagonals) - path->querylength +
	  Intlist_head(path->endpoints) + 1;

	donor_strand = path->fusion_plusp != invertp ? '+' : '-';
	acceptor_strand = path->plusp != invertp ? '+' : '-';
	donor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc1p);
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc2p);
	donor_chrpos = donor_position - path->fusion_chroffset;
	acceptor_chrpos = acceptor_position - path->chroffset;

      } else {
	/* 5' (fusion) is acceptor; 3' (main) is donor */
	if (path->fusion_plusp  == true) {
	  debug2(printf("Case 1, sense anti\n"));
	  acceptor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	} else {
	  debug2(printf("Case 2, sense anti\n"));
	  acceptor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	}
	donor_position = Univcoordlist_head(path->univdiagonals) - path->querylength +
	  Intlist_head(path->endpoints) + 1;

	acceptor_strand = path->fusion_plusp != invertp ? '+' : '-';
	donor_strand = path->plusp != invertp ? '-' : '+';
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc2p);
	donor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc1p);
	acceptor_chrpos = acceptor_position - path->fusion_chroffset;
	donor_chrpos = donor_position - path->chroffset;
      }

    } else {
      /* Case 4', 3' */
      /* 5' is fusion; 3' is main (minus strand) */
      if (path->sensedir == SENSE_FORWARD) {
	/* 5' (fusion) is donor; 3' (main) is acceptor */
	acceptor_position = Univcoordlist_last_value(path->univdiagonals) - path->querylength +
	  Intlist_last_value(path->endpoints);

	if (path->fusion_plusp == true) {
	  debug2(printf("Case 4', sense forward\n"));
	  donor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	} else {
	  debug2(printf("Case 3', sense forward\n"));
	  donor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	}

	acceptor_strand = path->plusp != invertp ? '+' : '-';
	donor_strand = path->fusion_plusp != invertp ? '+' : '-';
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc2p);
	donor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc1p);
	acceptor_chrpos = acceptor_position - path->chroffset;
	donor_chrpos = donor_position - path->fusion_chroffset;

      } else {
	/* 5' (fusion) is acceptor; 3' (main) is donor */
	debug2(printf("Case 4', 3', sense anti\n"));
	donor_position = Univcoordlist_last_value(path->univdiagonals) - path->querylength +
	  Intlist_last_value(path->endpoints);

	if (path->fusion_plusp == true) {
	  debug2(printf("Case 4', sense anti\n"));
	  acceptor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	} else {
	  debug2(printf("Case 3', sense anti\n"));
	  acceptor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	}

	donor_strand = path->plusp != invertp ? '-' : '+';
	acceptor_strand = path->fusion_plusp != invertp ? '-' : '+';
	donor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc1p);
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc2p);
	donor_chrpos = donor_position - path->chroffset;
	acceptor_chrpos = acceptor_position - path->fusion_chroffset;
      }
    }

  } else {
    fusion_junction = path->fusion_queryend_junction;

    if (path->plusp == true) {
      /* Case 3, 4 */
      /* 5' is main (plus strand); 3' is fusion */
      if (path->sensedir == SENSE_FORWARD) {
	/* 5' (main) is donor; 3' (fusion) is acceptor */
	debug2(printf("Case 3, 4, sense forward\n"));
	donor_position = Univcoordlist_last_value(path->univdiagonals) - path->querylength +
	  Intlist_last_value(path->endpoints);

	if (path->fusion_plusp == true) {
	  debug2(printf("Case 3, sense forward\n"));
	  acceptor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	} else {
	  debug2(printf("Case 4, sense forward\n"));
	  acceptor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	}

	donor_strand = path->plusp != invertp ? '-' : '+';
	acceptor_strand = path->fusion_plusp != invertp ? '-' : '+';
	donor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc1p);
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc2p);
	donor_chrpos = donor_position - path->chroffset;
	acceptor_chrpos = acceptor_position - path->fusion_chroffset;

      } else {
	/* 5' (main) is acceptor; 3' (fusion) is donor */
	debug2(printf("Case 3, 4, sense anti\n"));
	acceptor_position = Univcoordlist_last_value(path->univdiagonals) - path->querylength +
	  Intlist_last_value(path->endpoints);

	if (path->fusion_plusp == true) {
	  debug2(printf("Case 3, sense anti\n"));
	  donor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	} else {
	  debug2(printf("Case 4, sense anti\n"));
	  donor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	}

	acceptor_strand = path->plusp != invertp ? '+' : '-';
	donor_strand = path->fusion_plusp != invertp ? '+' : '-';
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc2p);
	donor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc1p);
	acceptor_chrpos = acceptor_position - path->chroffset;
	donor_chrpos = donor_position - path->fusion_chroffset;
      }

    } else {
      /* Case 2', 1' */
      /* 5' is main (minus strand); 3' is fusion */
      if (path->sensedir == SENSE_FORWARD) {
	/* 5' (main) is donor; 3' (fusion) is acceptor */
	if (path->fusion_plusp == true) {
	  debug2(printf("Case 2', sense forward\n"));
	  acceptor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	} else {
	  debug2(printf("Case 1', sense forward\n"));
	  acceptor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	}
	donor_position = Univcoordlist_head(path->univdiagonals) - path->querylength +
	  Intlist_head(path->endpoints) + 1;

	acceptor_strand = path->fusion_plusp != invertp ? '+' : '-';
	donor_strand = path->plusp != invertp ? '-' : '+';
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc2p);
	donor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc1p);
	acceptor_chrpos = acceptor_position - path->fusion_chroffset;
	donor_chrpos = donor_position - path->chroffset;

      } else {
	/* 5' (main) is acceptor; 3' (fusion) is donor */
	debug2(printf("Case 2', 1', sense anti\n"));
	if (path->fusion_plusp == true) {
	  debug2(printf("Case 2', sense anti\n"));
	  donor_position = Univcoordlist_head(path->fusion_univdiagonals) - path->querylength +
	    Intlist_head(path->fusion_endpoints) + 1;
	} else {
	  debug2(printf("Case 1', sense anti\n"));
	  donor_position = Univcoordlist_last_value(path->fusion_univdiagonals) - path->querylength +
	    Intlist_last_value(path->fusion_endpoints);
	}
	acceptor_position = Univcoordlist_head(path->univdiagonals) - path->querylength +
	  Intlist_head(path->endpoints) + 1;

	donor_strand = path->fusion_plusp != invertp ? '+' : '-';
	acceptor_strand = path->plusp != invertp ? '+' : '-';
	donor_chr = Univ_IIT_label(chromosome_iit,path->fusion_chrnum,&alloc1p);
	acceptor_chr = Univ_IIT_label(chromosome_iit,path->chrnum,&alloc2p);
	donor_chrpos = donor_position - path->fusion_chroffset;
	acceptor_chrpos = acceptor_position - path->chroffset;
      }
    }
  }
  
  FPRINTF(fp,"\tXT:Z:%c%c-%c%c,%.2f,%.2f",
	  fusion_junction->donor1,fusion_junction->donor2,
	  fusion_junction->acceptor2,fusion_junction->acceptor1,
	  fusion_junction->donor_prob,fusion_junction->acceptor_prob);
  FPRINTF(fp,",%c%s@%u..%c%s@%u",donor_strand,donor_chr,donor_chrpos,
	  acceptor_strand,acceptor_chr,acceptor_chrpos);
  

  /* 12. TAGS: XX, XY */
  if (path->fusion_transcripts != NULL) {
    Transcript_print_list(fp,path->fusion_transcripts,transcript_iit,/*header*/"\tXX:Z:");
  } else if (path->fusion_invalid_transcripts != NULL) {
    Transcript_print_list(fp,path->fusion_invalid_transcripts,transcript_iit,/*header*/"\tXY:Z:");
  }


  if (alloc1p) {
    FREE(donor_chr);
  }
  if (alloc2p) {
    FREE(acceptor_chr);
  }


#if 0
  /* TO_FIX: Use fusion_alts */
  /* 12. TAGS: XA */
  if (path->fusion_alts != NULL) {
    FPRINTF(fp,"\tXA:Z:");

    if (path->plusp == true) {
      alts_coords = path->qstart_alts->coords;
      n = path->qstart_alts->ncoords;
#ifdef PRINT_ALTS_COORDS
      FPRINTF(fp,"%u",alts_coords[0] - path->chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_coords[i] - path->chroffset + 1U);
      }
#else
      splicecoord = Univcoordlist_head(path->univdiagonals) - path->querylength + Intlist_head(path->endpoints);
      FPRINTF(fp,"%u",splicecoord - alts_coords[0]);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",splicecoord - alts_coords[i]);
      }
#endif
    }
    FPRINTF(fp,"|");
    if (path->qend_alts != NULL) {
      alts_coords = path->qend_alts->coords;
      n = path->qend_alts->ncoords;
#ifdef PRINT_ALTS_COORDS
      FPRINTF(fp,"%u",alts_coords[0] - path->chroffset + 1U);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_coords[i] - path->chroffset + 1U);
      }
#else
      splicecoord = Univcoordlist_last_value(path->univdiagonals) - path->querylength + Intlist_last_value(path->endpoints);
      FPRINTF(fp,"%u",alts_coords[0] - splicecoord);
      for (i = 1; i < n; i++) {
	FPRINTF(fp,",%u",alts_coords[i] - splicecoord);
      }
#endif
    }
  }
#endif

  /* 12. TAGS: (BC,) CB, CR, CY, UR, UY */
  if (single_cell_infoseq != NULL) {
    Single_cell_print_fields(fp,single_cell_infoseq);
  }

  FPRINTF(fp,"\n");

  Filestring_free(&mate_md_fp,/*free_string_p*/true);
  Filestring_free(&md_fp,/*free_string_p*/true);
  Filestring_free(&mate_cigar_fp,/*free_string_p*/true);
  Filestring_free(&cigar_fp,/*free_string_p*/true);

  return;
}


void
Path_print_sam (Filestring_T fp, Filestring_T *fp_failedinput, char *abbrev,
		T this, char *acc1, char *acc2, int pathnum, int npaths_primary, int npaths_altloc,
		int absmq_score, int first_absmq, int second_absmq, int mapq_score, Univ_IIT_T chromosome_iit,
		Shortread_T queryseq, Shortread_T mate_queryseq, Shortread_T single_cell_infoseq,
		int pairedlength, int pair_relationship,

		T mate, Resulttype_T resulttype, bool paired_read_p, bool first_read_p,
		bool artificial_mate_p, int npaths_mate,
		int quality_shift, char *sam_read_group_id, bool invertp, bool invert_mate_p,
		Listpool_T listpool) {

  int hardclip_low, hardclip_high;

  debug(printf("Entered Path_print_sam of hit %p and abbrev %s\n",this,abbrev));
  debug(Path_print(this));

  /* Test for nomapping was chrpos == 0, but we can actually align to chrpos 0 */
  /* Also, can use this test here because --quiet-if-excessive cases go directly to Path_print_sam_nomapping */
  assert(npaths_primary + npaths_altloc > 0);
  if (npaths_primary + npaths_altloc == 0) {
    Path_print_sam_nomapping(fp,abbrev,queryseq,mate_queryseq,single_cell_infoseq,
			     acc1,acc2,chromosome_iit,resulttype,first_read_p,
			     /*pathnum*/0,/*npaths_primary*/0,/*npaths_altloc*/0,artificial_mate_p,npaths_mate,
			     mate,quality_shift,sam_read_group_id,invertp,invert_mate_p);

    if (failedinput_root != (char *) NULL && *fp_failedinput == (Filestring_T) NULL) {
      *fp_failedinput = Filestring_new();
#if 0
      if (first_read_p == true) {
	Shortread_print_query_singleend(*fp_failedinput,queryseq,/*headerseq*/queryseq);
      } else {
	Shortread_print_query_singleend(*fp_failedinput,queryseq,/*headerseq*/mate_queryseq);
      }
#else
      Shortread_print_query_singleend(*fp_failedinput,queryseq,/*headerseq*/queryseq);
#endif
    }

  } else {
    print_substrings(&hardclip_low,&hardclip_high,fp,abbrev,this,
		     acc1,acc2,pathnum,npaths_primary,npaths_altloc,
		     absmq_score,first_absmq,second_absmq,mapq_score,
		     queryseq,mate_queryseq,single_cell_infoseq,pairedlength,pair_relationship,
		     mate,resulttype,first_read_p,artificial_mate_p,
		     npaths_mate,quality_shift,sam_read_group_id,
		     invertp,invert_mate_p);
    debug(printf("print_substrings returning hardclips %d and %d\n",hardclip_low,hardclip_high));
  }

  if (this->circular_endpoints != NULL) {
    /* Do not pass hardclip_low and hardclip_high to print_circular */
    print_circular(fp,abbrev,this,
		   acc1,acc2,pathnum,npaths_primary,npaths_altloc,
		   absmq_score,first_absmq,mapq_score,
		   queryseq,mate_queryseq,single_cell_infoseq,pairedlength,pair_relationship,
		   mate,resulttype,first_read_p,artificial_mate_p,
		   npaths_mate,quality_shift,sam_read_group_id,
		   invertp,invert_mate_p);
  }

  if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
    /* Do pass hardclip_low and hardclip_high to print_fusion */
    /* Print fusion part */
    print_fusion(fp,abbrev,this,
		 acc1,acc2,pathnum,npaths_primary,npaths_altloc,
		 absmq_score,first_absmq,mapq_score,
		 queryseq,mate_queryseq,single_cell_infoseq,pairedlength,pair_relationship,
		 hardclip_low,hardclip_high,mate,
		 resulttype,first_read_p,artificial_mate_p,
		 npaths_mate,quality_shift,sam_read_group_id,invertp,invert_mate_p);
  }

  return;
}


void
Path_print_sam_paired (Filestring_T fp, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
		       Result_T result, Resulttype_T resulttype, Univ_IIT_T chromosome_iit,
		       Shortread_T queryseq1, Shortread_T queryseq2, bool invert_first_p, bool invert_second_p,
		       bool nofailsp, bool failsonlyp, int quality_shift, char *sam_read_group_id,
		       Listpool_T listpool) {
  Pathpair_T *pathpairarray, pathpair;
  T *patharray1, *patharray2, path, mate, path5, path3;
  int npaths_primary, npaths_altloc, npaths_max, npaths_primary_max, npaths_altloc_max,
    npaths1_primary, npaths1_altloc, npaths2_primary, npaths2_altloc, pathnum;
  int first_absmq, second_absmq, first_absmq1, second_absmq1, first_absmq2, second_absmq2;
  int clipdir;
  char *acc1, *acc2;
  char *abbrev;
  bool concordant_softclipped_p;


  acc1 = Shortread_accession(queryseq1);
  acc2 = Shortread_accession(queryseq2); /* NULL, unless --allow-pe-name-mismatch is specified */

  debug(printf("Entered Path_print_sam_paired with resulttype %d\n",resulttype));

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (nofailsp == true) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */
      
    } else if (only_concordant_p == true || only_tr_consistent_p == true) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NM);
      Path_print_sam_nomapping(fp,ABBREV_NOMAPPING_1,queryseq1,/*mate_queryseq*/queryseq2,
			       /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
			       /*first_read_p*/true,/*pathnum*/0,/*npaths_primary*/0,/*npaths_altloc*/0,
			       /*artificial_mate_p*/false,/*npaths_mate*/0,
			       /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
      Path_print_sam_nomapping(fp,ABBREV_NOMAPPING_2,queryseq2,/*mate_queryseq*/queryseq1,
			       /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
			       /*first_read_p*/false,/*pathnum*/0,/*npaths_primary*/0,/*npaths_altloc*/0,
			       /*artificial_mate_p*/false,/*npaths_mate*/0,
			       /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);

      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }
    }

  } else {
    if (nofailsp == true && (only_tr_consistent_p == true && Result_tr_consistent_p(result) == false)) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (failsonlyp == true && (only_tr_consistent_p == false || Result_tr_consistent_p(result) == true)) {
      if (failedinput_root != NULL) {
	*fp_failedinput_1 = Filestring_new();
	*fp_failedinput_2 = Filestring_new();
	Shortread_print_query_pairedend(*fp_failedinput_1,*fp_failedinput_2,queryseq1,queryseq2);
      }

    } else if (resulttype == CONCORDANT_UNIQ) {
      pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
      /* Pathpair_eval(pathpairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      pathpair = pathpairarray[0];
      path5 = pathpair->path5;
      path3 = pathpair->path3;
      
      if (omit_concordant_uniq_p == true) {
	/* Skip printing: This means that filestring->string is NULL, so cannot attempt to call Filestring_string */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

      } else if (clip_overlap_p == false && merge_overlap_p == false) {
	clipdir = 0;
	/* hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0; */
	  
      } else {
#ifdef TO_FIX
	clipdir = Stage3pair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,stage3pair);
	debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
		      clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
#endif
      }

      if (Path_softclippedp(path5) == true || Path_softclippedp(path3) == true) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/true,OUTPUT_CU);
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CU);
      }
      abbrev = ABBREV_CONCORDANT_UNIQ;

#if 0
      univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					path5,querylength5,/*first_read_p*/true);
      univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					path3,querylength3,/*first_read_p*/false);
#endif

      if (merge_overlap_p == false || clipdir == 0) {
	/* print first end */
	Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path5,
		       acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		       pathpair->absmq_score,first_absmq,/*second_absmq*/0,
		       pathpair->mapq_score,chromosome_iit,
		       /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
		       pathpair->insertlength,pathpair->pair_relationship,
		       /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
		       /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
		       quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	/* print second end */
	Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path3,
		       acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		       pathpair->absmq_score,first_absmq,/*second_absmq*/0,
		       pathpair->mapq_score,chromosome_iit,
		       /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
		       pathpair->insertlength,pathpair->pair_relationship,
		       /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
		       /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
		       quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);

      } else {
#ifdef TO_FIX
	/* merge_overlap_p == true and overlap was found */
	pairarray = Pathpair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
				   pathpair,queryseq1,queryseq2,querylength5,querylength3,
				   clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	/* printf("queryseq_merged: %s\n",queryseq_merged); */

#if 0
	if (clipdir >= 0) {
	  chrnum = chrnum5;
	  univcoord_low = univcoord_low_5;
	} else {
	  chrnum = chrnum3;
	  univcoord_low = univcoord_low_3;
	}
#else
	univcoord_low = Simplepair_genomicpos_low(/*hardclip_low*/0,/*hardclip_high*/0,pairarray,npairs,
					       querylength_merged,/*plusp*/path5->plusp,
					       extend_soft_clips_p);
#endif

	/* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */

	flag = compute_flag(path5->plusp,/*mate*/NULL,
			    /*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
			    /*pathnum*/1,/*npaths*/1,/*artificial_mate_p*/false,/*npaths_mate*/0,
			    pathpair->absmq_score,first_absmq,/*invertp*/false,
			    /*invert_mate_p*/false,/*supplementaryp*/false);

	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
	Simplepair_overlap_print_sam(fp,/*abbrev*/ABBREV_UNPAIRED_UNIQ,pairarray,npairs,
				     acc1,/*acc2*/NULL,path5->chrnum,chromosome_iit,
				     /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
				     /*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
				     path5->plusp,Path_effective_sensedir(path5),quality_shift,
				     /*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
#if 0
				     pathpair->absmq_score,/*second_absmq*/0,
#else
				     /*absmq_score*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
				     univcoord_low,Path_chrlength(path5),/*queryseq*/NULL,flag,
				     /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
				     path5->method,sam_read_group_id);
	
	if (quality_merged != NULL) {
	  FREE_OUT(quality_merged);
	}
	FREE_OUT(queryseq_merged);
	FREE_OUT(pairarray);
#endif
      }

    } else if (resulttype == CONCORDANT_MULT) {
      pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

      if (omit_concordant_mult_p == true) {
	/* Skip printing: This means that filestring->string is NULL, so cannot attempt to call Filestring_string */
	/* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

      } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
	if (nofailsp == true) {
	  /* No output */
	} else {
	  /* Print as nomapping, but send to fp_*_xs */
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CX);
	  abbrev = ABBREV_CONCORDANT_MULT_XS;

	  Path_print_sam_nomapping(fp,abbrev,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	concordant_softclipped_p = false;

	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	  pathpair = pathpairarray[pathnum-1];
	  path5 = pathpair->path5;
	  path3 = pathpair->path3;
	  if (Path_softclippedp(path5) == true || Path_softclippedp(path3) == true) {
	    concordant_softclipped_p = true;
	  }
	}

	Filestring_set_split_output(fp,concordant_softclipped_p,OUTPUT_CM);
	abbrev = ABBREV_CONCORDANT_MULT;
	      
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {
	  pathpair = pathpairarray[pathnum-1];
	  path5 = pathpair->path5;
	  path3 = pathpair->path3;
#ifdef TO_FIX
	  clipdir = Pathpair_overlap(&hardclip5_low,&hardclip5_high,&hardclip3_low,&hardclip3_high,pathpair);
	  debug3(printf("clipdir %d with hardclip5 = %d..%d, hardclip3 = %d..%d\n",
			clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high));
#endif

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					       path5,querylength5,/*first_read_p*/true);
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					       path3,querylength3,/*first_read_p*/false);
#endif

	  if (1 || merge_overlap_p == false || clipdir == 0) {
	    /* print first end */
	    Path_print_sam(fp,&(*fp_failedinput_1),abbrev,
			   path5,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			   pathpair->absmq_score,first_absmq,second_absmq,
			   pathpair->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			   pathpair->insertlength,pathpair->pair_relationship,
			   /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			   quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);
	    
	    /* print second end */
	    Path_print_sam(fp,&(*fp_failedinput_2),abbrev,
			   path3,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			   pathpair->absmq_score,first_absmq,second_absmq,
			   pathpair->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			   pathpair->insertlength,pathpair->pair_relationship,
			   /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			   quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	  
	  } else {
#ifdef TO_FIX
	    /* merge_overlap_p == true and overlap was found */
	    pairarray = Pathpair_merge(&npairs,&querylength_merged,&queryseq_merged,&quality_merged,
				       pathpair,queryseq1,queryseq2,querylength5,querylength3,
				       clipdir,hardclip5_low,hardclip5_high,hardclip3_low,hardclip3_high);
	    /* printf("queryseq_merged: %s\n",queryseq_merged); */
#endif
	    
#if 0
	    if (clipdir >= 0) {
	      chrnum = chrnum5;
	      univcoord_low = univcoord_low_5;
	    } else {
	      chrnum = chrnum3;
	      univcoord_low = univcoord_low_3;
	    }
#endif

#ifdef TO_FIX
	    univcoord_low = Simplepair_genomicpos_low(/*hardclip_low*/0,/*hardclip_high*/0,pairarray,npairs,
						      querylength_merged,/*plusp*/path5->plusp,
						      extend_soft_clips_p);
	    /* merging changes resulttype from UNPAIRED_UNIQ to SINGLEEND_UNIQ */
	    flag = compute_flag(path5->plusp,/*mate*/NULL,
				/*resulttype*/SINGLEEND_UNIQ,/*first_read_p*/true,
				/*pathnum*/1,/*npaths*/1,/*artificial_mate_p*/false,/*npaths_mate*/0,
				pathpair->absmq_score,first_absmq,/*invertp*/false,
				/*invert_mate_p*/false,/*supplementaryp*/false);
	    
	    Simplepair_overlap_print_sam(fp,abbrev,pairarray,npairs,
					 acc1,/*acc2*/NULL,path5->chrnum,chromosome_iit,
					 /*queryseq_ptr*/queryseq_merged,/*quality_string*/quality_merged,
					 /*hardclip_low*/0,/*hardclip_high*/0,/*querylength*/querylength_merged,
					 path5->plusp,Path_effective_sensedir(path5),
					 quality_shift,pathnum,npaths_primary,npaths_altloc,
#if 0
					 pathpair->absmq_score,/*second_absmq*/0,
#else
					 /*absmq_score*/MAX_QUALITY_SCORE,/*second_absmq*/0,
#endif
					 univcoord_low,Path_chrlength(path5),/*queryseq*/NULL,flag,
					 /*pair_mapq_score*/MAX_QUALITY_SCORE,/*end_mapq_score*/MAX_QUALITY_SCORE,
					 path5->method,sam_read_group_id);
	    
	    if (quality_merged != NULL) {
	      FREE_OUT(quality_merged);
	    }
	    FREE_OUT(queryseq_merged);
	    FREE_OUT(pairarray);
#endif	
	  }
	}
      }

    } else if (only_concordant_p == true) {
      /* No output */
      /* Filestring_set_split_output(fp,OUTPUT_NONE); -- Already initialized to be OUTPUT_NONE */

    } else if (resulttype == CONCORDANT_TRANSLOC) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CT);
      pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
	if (nofailsp == true) {
	  /* No output */
	} else {
	  /* Print as nomapping, but send to fp_concordant_transloc */
	  Path_print_sam_nomapping(fp,ABBREV_CONCORDANT_TRANSLOC,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  Path_print_sam_nomapping(fp,ABBREV_CONCORDANT_TRANSLOC,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* Pathpair_eval(pathpairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	/* Note: We are ignoring merge_overlap for concordant_transloc */
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {

	  pathpair = pathpairarray[pathnum-1];
	  path5 = pathpair->path5;
	  path3 = pathpair->path3;

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					    path5,querylength5,/*first_read_p*/true);
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					    path3,querylength3,/*first_read_p*/false);
#endif

	  /* print first end */
	  Path_print_sam(fp,&(*fp_failedinput_1),ABBREV_CONCORDANT_TRANSLOC,
			 path5,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			 pathpair->absmq_score,first_absmq,second_absmq,
			 pathpair->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 pathpair->insertlength,pathpair->pair_relationship,
			 /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	  /* print second end */
	  Path_print_sam(fp,&(*fp_failedinput_2),ABBREV_CONCORDANT_TRANSLOC,
			 path3,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			 pathpair->absmq_score,first_absmq,second_absmq,
			 pathpair->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 pathpair->insertlength,pathpair->pair_relationship,
			 /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	}
      }
    

    } else if (resulttype == PAIRED_UNIQ_INV || resulttype == PAIRED_UNIQ_SCR || resulttype == PAIRED_UNIQ_TOOLONG) {
      pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
      /* Pathpair_eval(pathpairarray,npaths,maxpaths_report,queryseq1,queryseq2); */

      pathpair = pathpairarray[0];
      path5 = pathpair->path5;
      path3 = pathpair->path3;

      if (resulttype == PAIRED_UNIQ_INV) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PI);
	abbrev = ABBREV_PAIRED_UNIQ_INV;
      } else if (resulttype == PAIRED_UNIQ_SCR) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PS);
	abbrev = ABBREV_PAIRED_UNIQ_SCR;
      } else if (resulttype == PAIRED_UNIQ_TOOLONG) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PL);
	abbrev = ABBREV_PAIRED_UNIQ_LONG;
      } else {
	fprintf(stderr,"Unexpected resulttype %d\n",resulttype);
	abort();
      }

#if 0
      univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					path5,querylength5,/*first_read_p*/true);
      univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					path3,querylength3,/*first_read_p*/false);
#endif

      /* print first end */
      Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path5,
		     acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		     pathpair->absmq_score,first_absmq,/*second_absmq*/0,
		     pathpair->mapq_score,chromosome_iit,
		     /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
		     pathpair->insertlength,pathpair->pair_relationship,
		     /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
		     /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
		     quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

      /* print second end */
      Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path3,
		     acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		     pathpair->absmq_score,first_absmq,/*second_absmq*/0,
		     pathpair->mapq_score,chromosome_iit,
		     /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
		     pathpair->insertlength,pathpair->pair_relationship,
		     /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
		     /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
		     quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);

    } else if (resulttype == PAIRED_MULT) {
      pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths_report) {
	if (nofailsp == true) {
	  /* No output */
	} else {
	  /* Print as nomapping, but send to fp_paired_mult */
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PX);
	  Path_print_sam_nomapping(fp,ABBREV_PAIRED_MULT_XS,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  Path_print_sam_nomapping(fp,ABBREV_PAIRED_MULT_XS,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths_primary,npaths_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
				   /*mate*/(T) NULL,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* Pathpair_eval(pathpairarray,npaths,maxpaths_report,queryseq1,queryseq2); */
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PM);
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths_report; pathnum++) {

	  pathpair = pathpairarray[pathnum-1];
	  /* hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0; */

	  path5 = pathpair->path5;
	  path3 = pathpair->path3;

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					    path5,querylength5,/*first_read_p*/true);
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					    path3,querylength3,/*first_read_p*/false);
#endif

	  /* print first end */
	  Path_print_sam(fp,&(*fp_failedinput_1),ABBREV_PAIRED_MULT,
			 path5,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			 pathpair->absmq_score,first_absmq,second_absmq,
			 pathpair->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 pathpair->insertlength,pathpair->pair_relationship,
			 /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	  /* print second end */
	  Path_print_sam(fp,&(*fp_failedinput_2),ABBREV_PAIRED_MULT,
			 path3,acc1,acc2,pathnum,npaths_primary,npaths_altloc,
			 pathpair->absmq_score,first_absmq,second_absmq,
			 pathpair->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 pathpair->insertlength,pathpair->pair_relationship,
			 /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_primary + npaths_altloc,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	}
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      /* Even though they are not related, we should print mate information in this situation */
      patharray1 = (T *) Result_array(&npaths1_primary,&npaths1_altloc,&first_absmq1,&second_absmq1,result);
      patharray2 = (T *) Result_array2(&npaths2_primary,&npaths2_altloc,&first_absmq2,&second_absmq2,result);

      /* hardclip5_low = hardclip5_high = hardclip3_low = hardclip3_high = 0; */

      path5 = patharray1[0];
      path3 = patharray2[0];

#if 0
      univcoord_low_5 = SAM_compute_chrpos(&chrnum5,hardclip5_low,hardclip5_high,
					path5,querylength5,/*first_read_p*/true);
      univcoord_low_3 = SAM_compute_chrpos(&chrnum3,hardclip3_low,hardclip3_high,
					path3,querylength3,/*first_read_p*/false);
#endif

      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);
      abbrev = ABBREV_UNPAIRED_UNIQ;

      /* print first end */
      /* Path_eval_and_sort(patharray1,npaths1,maxpaths_report,queryseq1); */
      /* Previously set hardclips to be 0.  Not sure why. */
      Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path5,
		     acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		     patharray1[0]->absmq_score,first_absmq1,/*second_absmq*/0,
		     patharray1[0]->mapq_score,chromosome_iit,
		     /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
		     /*pairedlength*/0U,/*pair_relationship*/0,
		     /*mate*/path3,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
		     /*artificial_mate_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		     invert_first_p,invert_second_p,listpool);

      /* Note: Do not act on add_paired_nomappers_p, since the two reads are artificially paired up already */

      /* print second end */
      /* Path_eval_and_sort(patharray2,npaths2,maxpaths_report,queryseq2); */
      /* Previously set hardclips to be 0.  Not sure why. */
      Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path3,
		     acc1,acc2,/*pathnum*/1,/*npaths_primary*/1,/*npaths_altloc*/0,
		     patharray2[0]->absmq_score,first_absmq2,/*second_absmq*/0,
		     patharray2[0]->mapq_score,chromosome_iit,
		     /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
		     /*pairedlength*/0U,/*pair_relationship*/0,
		     /*mate*/path5,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
		     /*artificial_mate_p*/false,/*npaths_mate*/1,quality_shift,sam_read_group_id,
		     invert_second_p,invert_first_p,listpool);

    } else if (resulttype == UNPAIRED_MULT || resulttype == UNPAIRED_TRANSLOC) {
      if (resulttype == UNPAIRED_MULT) {
	if (quiet_if_excessive_p && npaths1_primary + npaths1_altloc > maxpaths_report &&
	    npaths2_primary + npaths2_altloc > maxpaths_report) {
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	} else {
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
	}
	abbrev = ABBREV_UNPAIRED_MULT;
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);
	abbrev = ABBREV_UNPAIRED_TRANSLOC;
      }

      patharray1 = (T *) Result_array(&npaths1_primary,&npaths1_altloc,&first_absmq1,&second_absmq1,result);
      patharray2 = (T *) Result_array2(&npaths2_primary,&npaths2_altloc,&first_absmq2,&second_absmq2,result);

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 1) {
	/* Path_eval_and_sort(patharray1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Path_eval_and_sort(patharray1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 1) {
	/* Path_eval_and_sort(patharray2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Path_eval_and_sort(patharray2,npaths2,maxpaths_report,queryseq2); */
      }
#endif

      if (add_paired_nomappers_p == true) {
	/* Artificially pair up results */
	if (npaths1_primary + npaths1_altloc > npaths2_primary + npaths2_altloc) {
	  npaths_primary_max = npaths1_primary;
	  npaths_altloc_max = npaths1_altloc;
	  npaths_max = npaths1_primary + npaths1_altloc;
	} else {
	  npaths_primary_max = npaths2_primary;
	  npaths_altloc_max = npaths2_altloc;
	  npaths_max = npaths2_primary + npaths2_altloc;
	}
	for (pathnum = 1; pathnum <= npaths1_primary + npaths1_altloc &&
	       pathnum <= npaths2_primary + npaths2_altloc && pathnum <= maxpaths_report; pathnum++) {
#if 0
	  /* hardclip5_low = hardclip5_high = 0; */
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					    /*path*/patharray1[pathnum-1],
					    querylength5,/*first_read_p*/true);

	  /* hardclip3_low = hardclip3_high = 0; */
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					    /*path*/patharray2[pathnum-1],
					    querylength3,/*first_read_p*/false);
#endif

	  path = patharray1[pathnum-1];
	  Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			 acc1,acc2,pathnum,npaths_primary_max,npaths_altloc_max,
			 path->absmq_score,first_absmq1,second_absmq1,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 /*mate*/patharray2[pathnum-1],
			 resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	  path = patharray2[pathnum-1];
	  Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			 acc1,acc2,pathnum,npaths_primary_max,npaths_altloc_max,
			 path->absmq_score,first_absmq2,second_absmq2,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 /*mate*/patharray1[pathnum-1],
			 resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	}

	/* Print remaining results with non-mappers */
	if (npaths1_primary + npaths1_altloc > npaths2_primary + npaths2_altloc) {
	  for ( ; pathnum <= npaths1_primary + npaths1_altloc && pathnum <= maxpaths_report; pathnum++) {
	    path = patharray1[pathnum-1];
#if 0
	    /* hardclip5_low = hardclip5_high = 0; */
	    univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					      path,querylength5,/*first_read_p*/true);
	    chrnum3 = 0;
	    univcoord_low_3 = 0;
#endif

	    Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			   acc1,acc2,pathnum,npaths_primary_max,npaths_altloc_max,
			   path->absmq_score,first_absmq1,second_absmq1,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,/*mate*/NULL,
			   resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			   /*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
			   quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	    /* matching nomapper for second end */
	    Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				     /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				     /*first_read_p*/false,pathnum,npaths_primary_max,npaths_altloc_max,
				     /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				     /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  }

	} else if (npaths2_primary + npaths2_altloc > npaths1_primary + npaths1_altloc) {
	  for ( ; pathnum <= npaths2_primary + npaths2_altloc && pathnum <= maxpaths_report; pathnum++) {
	    path = patharray2[pathnum-1];
#if 0
	    /* hardclip3_low = hardclip3_high = 0; */
	    chrnum5 = 0;
	    univcoord_low_5 = 0;
	    univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					      path,querylength3,/*first_read_p*/false);
#endif

	    /* matching nomapper for first end */
	    Path_print_sam_nomapping(fp,abbrev,queryseq1,/*mate_queryseq*/queryseq2,
				     /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				     /*first_read_p*/true,pathnum,npaths_primary_max,npaths_altloc_max,
				     /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				     /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

	    Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			   acc1,acc2,pathnum,npaths_primary_max,npaths_altloc_max,
			   path->absmq_score,first_absmq2,second_absmq2,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,/*mate*/NULL,
			   resulttype,/*paired_read_p*/true,/*first_read_p*/false,/*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
			   quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	  }
	}

      } else {
	/* print first end results */
	if (npaths2_primary + npaths2_altloc == 0) {
	  mate = (T) NULL;
	  /* chrnum3 = 0; */
	  /* univcoord_low_3 = 0U; */
	} else if (quiet_if_excessive_p && npaths2_primary + npaths2_altloc > maxpaths_report) {
	  mate = (T) NULL;
	  /* chrnum3 = 0; */
	  /* univcoord_low_3 = 0U; */
	} else {
	  mate = patharray2[0];
	  /* hardclip3_low = hardclip3_high = 0; */

#if 0
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					    mate,querylength3,/*first_read_p*/false);
#endif
	}

	if (npaths1_primary + npaths1_altloc == 1) {
	  path = patharray1[0];
	  /* hardclip5_low = hardclip5_high = 0; */

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					    path,querylength5,/*first_read_p*/true);
#endif

	  Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths1_primary,npaths1_altloc,
			 path->absmq_score,first_absmq1,second_absmq1,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,mate,
			 resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);

	} else if (quiet_if_excessive_p && npaths1_primary + npaths1_altloc > maxpaths_report) {
	  /* Just printing one end as nomapping */
	  Path_print_sam_nomapping(fp,abbrev,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths1_primary,npaths1_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
				   mate,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);

	} else {
	  for (pathnum = 1; pathnum <= npaths1_primary + npaths1_altloc && pathnum <= maxpaths_report; pathnum++) {
	    path = patharray1[pathnum-1];
	    /* hardclip5_low = hardclip5_high = 0; */

#if 0
	    univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/hardclip5_high,
					      path,querylength5,/*first_read_p*/true);
#endif
	    
	    Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			   acc1,acc2,pathnum,npaths1_primary,npaths1_altloc,
			   path->absmq_score,first_absmq1,second_absmq1,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,mate,
			   resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
			   quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);
	  }
	}
			  
	/* print second end results */
	if (npaths1_primary + npaths1_altloc == 0) {
	  mate = (T) NULL;
	  /* chrnum5 = 0; */
	  /* univcoord_low_5 = 0U; */
	} else if (quiet_if_excessive_p && npaths1_primary + npaths1_altloc > maxpaths_report) {
	  mate = (T) NULL;
	  /* chrnum5 = 0;*/
	  /* univcoord_low_5 = 0U; */
	} else {
	  mate = patharray1[0];
	  /* hardclip5_low = hardclip5_high = 0; */

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					    mate,querylength5,/*first_read_p*/true);
#endif
	}

	if (npaths2_primary + npaths2_altloc == 1) {
	  path = patharray2[0];
	  /* hardclip3_low = hardclip3_high = 0; */

#if 0
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					    path,querylength3,/*first_read_p*/false);
#endif
	  
	  Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths2_primary,npaths2_altloc,
			 path->absmq_score,first_absmq2,second_absmq2,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,mate,
			 resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	  
	} else if (quiet_if_excessive_p && npaths2_primary + npaths2_altloc > maxpaths_report) {
	  /* Just printing one end as nomapping */
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths2_primary,npaths2_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
				   mate,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	  
	} else {
	  for (pathnum = 1; pathnum <= npaths2_primary + npaths2_altloc && pathnum <= maxpaths_report; pathnum++) {
	    path = patharray2[pathnum-1];
	    /* hardclip3_low = hardclip3_high = 0; */

#if 0
	    univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					      path,querylength3,/*first_read_p*/false);
#endif

	    Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			   acc1,acc2,pathnum,npaths2_primary,npaths2_altloc,
			   path->absmq_score,first_absmq2,second_absmq2,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,
			   mate,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
			   quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	  }
	}
      }

    } else {
      patharray1 = (T *) Result_array(&npaths1_primary,&npaths1_altloc,&first_absmq1,&second_absmq1,result);
      patharray2 = (T *) Result_array2(&npaths2_primary,&npaths2_altloc,&first_absmq2,&second_absmq2,result);

      if (resulttype == HALFMAPPING_UNIQ) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HU);
	abbrev = ABBREV_HALFMAPPING_UNIQ;

      } else if (resulttype == HALFMAPPING_TRANSLOC) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HT);
	abbrev = ABBREV_HALFMAPPING_TRANSLOC;

      } else if (resulttype == HALFMAPPING_MULT) {
	if (quiet_if_excessive_p == true && npaths1_primary + npaths1_altloc > maxpaths_report && npaths2_primary + npaths2_altloc > maxpaths_report) {
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HX);
	  abbrev = ABBREV_HALFMAPPING_MULT_XS;
	} else {
	  Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HM);
	  abbrev = ABBREV_HALFMAPPING_MULT;
	}
      } else {
	abort();
      }

#if 0
      /* Do eval and sorting first */
      if (npaths1 == 0) {
	/* Nothing to sort */
      } else if (npaths1 == 1) {
	/* Path_eval_and_sort(patharray1,npaths1,maxpaths_report,queryseq1); */
      } else if (quiet_if_excessive_p && npaths1 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Path_eval_and_sort(patharray1,npaths1,maxpaths_report,queryseq1); */
      }
#endif

#if 0
      if (npaths2 == 0) {
	/* Nothing to sort */
      } else if (npaths2 == 1) {
	/* Path_eval_and_sort(patharray2,npaths2,maxpaths_report,queryseq2); */
      } else if (quiet_if_excessive_p && npaths2 > maxpaths_report) {
	/* Don't sort */
      } else {
	/* Path_eval_and_sort(patharray2,npaths2,maxpaths_report,queryseq2); */
      }
#endif


      /* print first end results */
      if (npaths2_primary + npaths2_altloc == 0) {
	mate = (T) NULL;
	/* chrnum3 = 0; */
	/* univcoord_low_3 = 0U; */
      } else if (quiet_if_excessive_p && npaths2_primary + npaths2_altloc > maxpaths_report) {
	mate = (T) NULL;
	/* chrnum3 = 0; */
	/* univcoord_low_3 = 0U; */
      } else {
	mate = patharray2[0];
	/* hardclip3_low = hardclip3_high = 0; */

#if 0
	univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					  mate,querylength3,/*first_read_p*/false);
#endif
      }

      if (npaths1_primary + npaths1_altloc == 0) {
	/* just printing one end as nomapping */
	/* mate should be non-NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  Path_print_sam_nomapping(fp,abbrev,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/0,npaths1_primary,npaths1_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
				   mate,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else if (npaths1_primary + npaths1_altloc == 1) {
	/* mate should be NULL here */

	path = patharray1[0];
	/* hardclip5_low = hardclip5_high = 0; */

#if 0
	univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					  path,querylength5,/*first_read_p*/true);
#endif

	if (add_paired_nomappers_p == true) {
	  /* matching nomapper for second end */
	  npaths_max = npaths1_primary + npaths1_altloc; /* since npaths2_primary + npaths2_altloc == 0 */
	  Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths1_primary,npaths1_altloc,
			 path->absmq_score,first_absmq1,/*second_absmq1*/0,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 mate,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths1_primary,npaths1_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				   /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	} else {
	  Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths1_primary,npaths1_altloc,
			 path->absmq_score,first_absmq1,/*second_absmq1*/0,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 mate,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
			 quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);
	}

      } else if (quiet_if_excessive_p && npaths1_primary + npaths1_altloc > maxpaths_report) {
	/* Just printing one end as nomapping */
	/* mate should be NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  Path_print_sam_nomapping(fp,abbrev,queryseq1,/*mate_queryseq*/queryseq2,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths1_primary,npaths1_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
				   mate,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths1_primary + npaths1_altloc && pathnum <= maxpaths_report; pathnum++) {
	  path = patharray1[pathnum-1];
	  /* hardclip5_low = hardclip5_high = 0; */

#if 0
	  univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					    path,querylength5,/*first_read_p*/true);
#endif

	  if (add_paired_nomappers_p == true) {
	    /* matching nomapper for second end */
	    npaths_max = npaths1_primary + npaths1_altloc; /* since npaths2 == 0 */
	    Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				     /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				     /*first_read_p*/false,pathnum,npaths1_primary,npaths1_altloc,
				     /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				     /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	    Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			   acc1,acc2,pathnum,npaths1_primary,npaths1_altloc,
			   path->absmq_score,first_absmq1,second_absmq1,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,
			   mate,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			   /*artificial_mate_p*/true,/*npaths_mate*/npaths_max,quality_shift,sam_read_group_id,
			   invert_first_p,invert_second_p,listpool);

	  } else {
	    Path_print_sam(fp,&(*fp_failedinput_1),abbrev,path,
			   acc1,acc2,pathnum,npaths1_primary,npaths1_altloc,
			   path->absmq_score,first_absmq1,second_absmq1,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq1,/*mate_queryseq*/queryseq2,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,
			   mate,resulttype,/*paired_read_p*/true,/*first_read_p*/true,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths2_primary + npaths2_altloc,
			   quality_shift,sam_read_group_id,invert_first_p,invert_second_p,listpool);
	  }
	}
      }
			  
      /* print second end results */
      if (npaths1_primary + npaths1_altloc == 0) {
	mate = (T) NULL;
	/* chrnum5 = 0; */
	/* univcoord_low_5 = 0U; */
      } else if (quiet_if_excessive_p && npaths1_primary + npaths1_altloc > maxpaths_report) {
	mate = (T) NULL;
	/* chrnum5 = 0; */
	/* univcoord_low_5 = 0U; */
      } else {
	mate = patharray1[0];
	/* hardclip5_low = hardclip5_high = 0; */

#if 0
	univcoord_low_5 = SAM_compute_chrpos(&chrnum5,/*hardclip_low*/0,/*hardclip_high*/0,
					  mate,querylength5,/*first_read_p*/true);
#endif
      }

      if (npaths2_primary + npaths2_altloc == 0) {
	/* Just printing one end as nomapping */
	/* mate should be non-NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/0,npaths2_primary,npaths2_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
				   mate,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else if (npaths2_primary + npaths2_altloc == 1) {
	/* mate should be NULL here */

	path = patharray2[0];
	/* hardclip3_low = hardclip3_high = 0; */

#if 0
	univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					  path,querylength3,/*first_read_p*/false);
#endif

	if (add_paired_nomappers_p == true) {
	  /* matching nomapper for first end */
	  npaths_max = npaths2_primary + npaths2_altloc; /* since npaths1_primary + npaths1_altloc == 0 */
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/true,/*pathnum*/1,npaths2_primary,npaths2_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				   /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	  Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths2_primary,npaths2_altloc,
			 path->absmq_score,first_absmq2,/*second_absmq2*/0,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 mate,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);

	} else {
	  Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			 acc1,acc2,/*pathnum*/1,npaths2_primary,npaths2_altloc,
			 path->absmq_score,first_absmq2,/*second_absmq2*/0,
			 path->mapq_score,chromosome_iit,
			 /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			 /*pairedlength*/0U,/*pair_relationship*/0,
			 mate,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			 /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
			 quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	}

      } else if (quiet_if_excessive_p && npaths2_primary + npaths2_altloc > maxpaths_report) {
	/* Just printing one end as nomapping */
	/* mate should be NULL here */
	if (add_paired_nomappers_p == true) {
	  /* Handle nomappers with each mapped mate */
	} else {
	  Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				   /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				   /*first_read_p*/false,/*pathnum*/1,npaths2_primary,npaths2_altloc,
				   /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
				   mate,quality_shift,sam_read_group_id,invert_second_p,invert_first_p);
	}

      } else {
	/* mate should be NULL here */
	for (pathnum = 1; pathnum <= npaths2_primary + npaths2_altloc && pathnum <= maxpaths_report; pathnum++) {
	  path = patharray2[pathnum-1];
	  /* hardclip3_low = hardclip3_high = 0; */

#if 0
	  univcoord_low_3 = SAM_compute_chrpos(&chrnum3,/*hardclip_low*/0,/*hardclip_high*/0,
					    path,querylength3,/*first_read_p*/false);
#endif

	  if (add_paired_nomappers_p == true) {
	    /* matching nomapper for first end */
	    npaths_max = npaths2_primary + npaths2_altloc; /* since npaths1_primary + npaths1_altloc == 0 */
	    Path_print_sam_nomapping(fp,abbrev,queryseq2,/*mate_queryseq*/queryseq1,
				     /*single_cell_infoseq*/NULL,acc1,acc2,chromosome_iit,resulttype,
				     /*first_read_p*/true,pathnum,npaths2_primary,npaths2_altloc,
				     /*artificial_mate_p*/false,/*npaths_mate*/npaths_max,
				     /*mate*/path,quality_shift,sam_read_group_id,invert_first_p,invert_second_p);
	    Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			   acc1,acc2,pathnum,npaths2_primary,npaths2_altloc,
			   path->absmq_score,first_absmq2,second_absmq2,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,
			   mate,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			   /*artificial_mate_p*/true,/*npaths_mate*/npaths_max,
			   quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);

	  } else {
	    Path_print_sam(fp,&(*fp_failedinput_2),abbrev,path,
			   acc1,acc2,pathnum,npaths2_primary,npaths2_altloc,
			   path->absmq_score,first_absmq2,second_absmq2,
			   path->mapq_score,chromosome_iit,
			   /*queryseq*/queryseq2,/*mate_queryseq*/queryseq1,/*single_cell_infoseq*/NULL,
			   /*pairedlength*/0U,/*pair_relationship*/0,
			   mate,resulttype,/*paired_read_p*/true,/*first_read_p*/false,
			   /*artificial_mate_p*/false,/*npaths_mate*/npaths1_primary + npaths1_altloc,
			   quality_shift,sam_read_group_id,invert_second_p,invert_first_p,listpool);
	  }
	}
      }

    }
  }

  return;
}


void
Path_print_sam_setup (bool add_paired_nomappers_p_in, bool paired_flag_means_concordant_p_in, bool sam_insert_0M_p_in,
		      bool quiet_if_excessive_p_in, int maxpaths_report_in,
		      char *failedinput_root_in, bool fastq_format_p_in, bool extend_soft_clips_p_in, bool method_print_p_in,
		      bool only_concordant_p_in, bool omit_concordant_uniq_p_in, bool omit_concordant_mult_p_in, 
		      bool only_tr_consistent_p_in,
		      bool *circularp_in, bool clip_overlap_p_in, bool merge_overlap_p_in, bool merge_samechr_p_in,
		      bool sam_multiple_primaries_p_in, bool sam_sparse_secondaries_p_in,
		      Univ_IIT_T chromosome_iit_in, Univ_IIT_T transcript_iit_in,
		      IIT_T snps_iit_in, bool maskedp_in) {

  add_paired_nomappers_p = add_paired_nomappers_p_in;
  paired_flag_means_concordant_p = paired_flag_means_concordant_p_in;
  sam_insert_0M_p = sam_insert_0M_p_in;

  quiet_if_excessive_p = quiet_if_excessive_p_in;
  maxpaths_report = maxpaths_report_in;

  failedinput_root = failedinput_root_in;
  fastq_format_p = fastq_format_p_in;
  extend_soft_clips_p = extend_soft_clips_p_in;
  method_print_p = method_print_p_in;

  only_concordant_p = only_concordant_p_in;
  omit_concordant_uniq_p = omit_concordant_uniq_p_in;
  omit_concordant_mult_p = omit_concordant_mult_p_in;
  only_tr_consistent_p = only_tr_consistent_p_in;

  circularp = circularp_in;

  clip_overlap_p = clip_overlap_p_in;
  merge_overlap_p = merge_overlap_p_in;
  merge_samechr_p = merge_samechr_p_in;

  sam_multiple_primaries_p = sam_multiple_primaries_p_in;
  sam_sparse_secondaries_p = sam_sparse_secondaries_p_in;

  /* force_xs_direction_p = force_xs_direction_p_in; */

  /* find_dna_chimeras_p = find_dna_chimeras_p_in; */
  /* splicingp = splicingp_in; */
  /* donor_typeint = donor_typeint_in; */
  /* acceptor_typeint = acceptor_typeint_in; */

  chromosome_iit = chromosome_iit_in;
  transcript_iit = transcript_iit_in;
  snps_iit = snps_iit_in;
  maskedp = maskedp_in;

  /* transcript_iit = transcript_iit_in; */

  return;
}


