#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-print-alignment.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For tolower */

#include "bool.h"
#include "assert.h"
#include "univcoord.h"
#include "junction.h"
#include "altsplice.h"
#include "complement.h"


/* Previously turned off if CHECK_ASSERTIONS was true */
#define PRINT_SUB 1


/* making substrings */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* embellish_genomic */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

static char complCode[128] = COMPLEMENT_LC;
static Univ_IIT_T chromosome_iit;
static bool method_print_p;
static bool print_univdiagonal_p;
static bool print_nsnpdiffs_p;


#define T Path_T


static char *
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return sequence;
}


static void
print_forward (Filestring_T fp, char *string, int n) {
  
  FPRINTF(fp,"%.*s",n,string);
  return;
}

static void
print_lc (Filestring_T fp, char *string, int n) {
  int i;
  
  for (i = 0; i < n; i++) {
    FPRINTF(fp,"%c",(char) tolower(string[i]));
  }
  return;
}

static void
print_revcomp (Filestring_T fp, char *nt, int len) {

  FPRINTF(fp,"%.*R",len,nt);
  return;
}

static void
print_revcomp_lc (Filestring_T fp, char *nt, int len) {
  int i;

  for (i = len-1; i >= 0; --i) {
    FPRINTF(fp,"%c",(char) tolower(complCode[(int) nt[i]]));
  }
  return;
}


#if 0
/* Taken from Shortread_print_subseq_uc */
static void
print_subseq_uc (Filestring_T fp, char *queryptr,
		 int qstart, int qend, int querylength) {
  int i;

  for (i = 0; i < qstart; i++) {
    FPRINTF(fp,"-");
  }

  FPRINTF(fp,"%.*s",qend - qstart,&(queryptr[qstart]));

  for (i = qend; i < querylength; i++) {
    FPRINTF(fp,"-");
  }

  return;
}
#endif


#if 0
/* Taken from Shortread_print_subseq_revcomp_uc */
static void
print_subseq_revcomp_uc (Filestring_T fp, char *queryptr,
			 int qstart, int qend, int querylength) {
  int i;

  for (i = querylength-1; i >= qend; --i) {
    FPRINTF(fp,"-");
  }

  FPRINTF(fp,"%.*R",qend - qstart,&(queryptr[qstart]));

  for (i = qstart-1; i >= 0; --i) {
    FPRINTF(fp,"-");
  }

  return;
}
#endif


static void
fill_w_dashes (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '-';
  }
  return;
}

static void
fill_w_stars (char *string, int start, int end) {
  int i;

  for (i = start; i < end; i++) {
    string[i] = '*';
  }
  return;
}


/* The result variable should be an existing genomic_bothdiff or genomic_refdiff, and then reassigned to that */
static char *
embellish_genomic (char *genomic_diff, int qstart, int qend, int querylength,
		   int chrbound_trim_qstart, int chrbound_trim_qend,
		   Univcoord_T left, int extraleft, int extraright,
		   bool substring_plusp, bool main_plusp) {
  char *result, *gbuffer;
  int i;

  /* Previously used qstart < qend, but a deletion can yield an empty segment */
  assert(qstart <= qend);

  debug1(printf("Entered embellish_genomic with qstart %d, qend %d, querylength %d, chrbound_trim_qstart %d, chrbound_trim_qend %d\n",
		qstart,qend,querylength,chrbound_trim_qstart,chrbound_trim_qend));
  if (extraleft > qstart) {
    extraleft = qstart;
  }
  if (qend + extraright > querylength) {
    extraright = querylength - qend;
  }

  assert(chrbound_trim_qstart >= 0);
  assert(chrbound_trim_qend >= 0);
  assert(chrbound_trim_qstart < querylength);
  assert(chrbound_trim_qend < querylength);

#if 0
  if (plusp == true) {
    chrbound_trim_qstart = chrbound_trim_qstart_in;
    chrbound_trim_qend = chrbound_trim_qend_in;
  } else {
    chrbound_trim_qstart = chrbound_trim_qend_in;
    chrbound_trim_qend = chrbound_trim_qstart_in;
  }
#endif


  assert(genomic_diff != NULL);
  result = (char *) MALLOC_OUT((querylength+1) * sizeof(char));
  strcpy(result,genomic_diff);
  debug1(printf("A g_diff: %s (%d..%d) extraleft:%d extraright:%d\n",result,qstart,qend,extraleft,extraright));
#if 0
  for (i = 0; i < querylength; i++) {
    result[i] = '?';
  }
#endif

  if (extraleft == 0 && extraright == 0) {
    gbuffer = (char *) NULL;
  } else {
    gbuffer = (char *) MALLOC((querylength+1)*sizeof(char));
    Genome_fill_buffer(left,querylength,gbuffer);
    if (substring_plusp != main_plusp) {
      make_complement_inplace(gbuffer,querylength);
    }
    debug1(printf("A buffer: %s\n",gbuffer));
  }

  /* Add aligned region with lower-case diffs, surrounded by dashes */
  fill_w_dashes(result,0,qstart-extraleft);
  fill_w_stars(result,0,chrbound_trim_qstart);

  /* Don't need to know adj anymore, because each substring has its own left */
#if 0
  debug1(printf("Copying from genomic_diff[%d] to result[%d] for a length of %d - %d\n",qstart,qstart,qend,qstart));
#endif
  strncpy(&(result[qstart-extraleft]),&(gbuffer[qstart-extraleft]),extraleft);
  for (i = qstart - extraleft; i < qstart; i++) {
    result[i] = tolower(result[i]);
  }
  strncpy(&(result[qstart]),&(genomic_diff[qstart]),qend-qstart);
  strncpy(&(result[qend]),&(gbuffer[qend]),extraright);
  for (i = qend; i < qend + extraright; i++) {
    result[i] = tolower(result[i]);
  }
  debug1(printf("B extra:  %s\n",result));

  fill_w_dashes(result,qend+extraright,querylength);
  fill_w_stars(result,querylength - chrbound_trim_qend,querylength);
  debug1(printf("C dashes: %s\n",result));

  if (gbuffer != NULL) {
    FREE(gbuffer);
  }

  return result;
}


static void
print_genomic (Filestring_T fp, char *genomic_diff, int querystart, int queryend, int querylength,
	       char *deletion, int deletionlength, Univcoord_T left, int extraleft, int extraright,
	       bool substring_plusp, bool main_plusp, Shortread_T queryseq, bool invertp) {
  char *fragment;
  int qstart, qend;
  int querystart_choplength;
  
  /* Previously handled by calculation of queryleft/queryright by Stage3end_display_prep */
  /* deletion = NULL; */
  /* deletionlength = 0; */

  /* TODO: extraleft and extraright are not printed correctly */
  if (invertp == false) {
    Shortread_print_left_chop_symbols(fp,queryseq);
    querystart_choplength = Shortread_left_choplength(queryseq);
  } else {
    Shortread_print_right_chop_symbols(fp,queryseq);
    querystart_choplength = Shortread_right_choplength(queryseq);
  }

  if (main_plusp != invertp) {
    qstart = querystart;
    qend = queryend;

    fragment =
      embellish_genomic(genomic_diff,qstart,qend,querylength,
			/*chrbound_trim_qstart*/0,/*chrbound_trim_qend*/0,
			left,extraleft,extraright,substring_plusp,main_plusp);

    print_forward(fp,fragment,qend);
    if (deletion != NULL) {
      print_lc(fp,deletion,deletionlength);
    }
    print_forward(fp,&(fragment[qend]),querylength - deletionlength - qend);

  } else {
    qstart = querylength - queryend;
    qend = querylength - querystart;

    fragment =
      embellish_genomic(genomic_diff,qstart,qend,querylength,
			/*chrbound_trim_querystart*/0,/*chrbound_trim_queryend*/0,
			left,extraright,extraleft,substring_plusp,main_plusp);

    print_revcomp(fp,&(fragment[qstart]),querylength - qstart);
    if (deletion != NULL) {
      print_revcomp_lc(fp,deletion,deletionlength);
    }
    print_revcomp(fp,&(fragment[deletionlength]),qstart - deletionlength);
  }

  if (invertp == false) {
    Shortread_print_right_chop_symbols(fp,queryseq);
  } else {
    Shortread_print_left_chop_symbols(fp,queryseq);
  }

  FPRINTF(fp,"\t");
  FPRINTF(fp,"%d..%d",1 + querystart_choplength + querystart,querystart_choplength + queryend);
  FREE_OUT(fragment);

  return;
}


static void
print_coordinates (Filestring_T fp, Univcoord_T univdiagonal, int querystart, int queryend, int querylength,
		   char *chr, Univcoord_T chroffset, bool substring_plusp, bool invertp) {
  int qstart, qend;

  if (substring_plusp != invertp) {
    qstart = querystart;
    qend = queryend;
    if (print_univdiagonal_p == true) {
      FPRINTF(fp,"+%u..%u",univdiagonal - querylength + qstart + 1,
	      univdiagonal - querylength + qend);
    } else {
      FPRINTF(fp,"+%s:%u..%u",chr,univdiagonal - querylength + qstart - chroffset + 1,
	      univdiagonal - querylength + qend - chroffset);
    }
  } else {
    qstart = querylength - queryend;
    qend = querylength - querystart;
    if (print_univdiagonal_p == true) {
      FPRINTF(fp,"-%u..%u",univdiagonal - querylength + qend,
	      univdiagonal - querylength + qstart + 1);
    } else {
      FPRINTF(fp,"-%s:%u..%u",chr,univdiagonal - querylength + qend - chroffset,
	      univdiagonal - querylength + qstart - chroffset + 1);
    }
  }

  return;
}



typedef struct Substring_T *Substring_T;
struct Substring_T {
  Univcoord_T univdiagonal;
  int nmismatches;
  int ref_nmismatches;
  int querystart;
  int queryend;
  Junction_T pre_junction;
  Junction_T post_junction;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  bool plusp;
};


static void
Substring_free (Substring_T *old) {
  FREE(*old);
  return;
}


static Substring_T
Substring_new (Univcoord_T univdiagonal, int nmismatches, int ref_nmismatches,
	       int querystart, int queryend,
	       Junction_T pre_junction, Junction_T post_junction,
	       Chrnum_T chrnum, Univcoord_T chroffset, bool plusp) {
  Substring_T new = (Substring_T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->nmismatches = nmismatches;
  new->ref_nmismatches = ref_nmismatches;
  new->querystart = querystart;
  new->queryend = queryend;
  new->pre_junction = pre_junction;
  new->post_junction = post_junction;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  return new;
}


/* Modified from Substring_print_alignment */
static void
print_substring (Filestring_T fp, Substring_T substring, Path_T path, Shortread_T queryseq,
		 bool invertp) {
  char *chr;
  char *deletion_string;
  int deletion_length;
  int extraleft, extraright;
  Junctiontype_T type1, type2;
  Altsplice_T altsplice;
  Chrpos_T splice_distance_1, splice_distance_2;
  Univcoord_T deletionpos;
  bool sense_forward_p, allocp;

  sense_forward_p = (path->sensedir == SENSE_FORWARD) ? true : false;

  /* assert(substring->querystart < substring->queryend); */

  chr = Univ_IIT_label(chromosome_iit,substring->chrnum,&allocp);

  if (substring->pre_junction == NULL) {
#if 1
    /* Extend to the start */
    extraleft = substring->querystart;
#else
    if (substring->plusp != invertp && path->splicetype5 != NO_SPLICE) {
      extraleft = 2;
    } else if (substring->plusp == invertp && path->splicetype3 != NO_SPLICE) {
      extraleft = 2;
    } else {
      extraleft = 0;
    }
#endif
  } else if (Junction_type(substring->pre_junction) == INS_JUNCTION) {
    extraleft = 0;
    if (substring->plusp != invertp) {
      debug0(printf("pre-junction is INS, so adding %d to querystart %d\n",
		    Junction_ninserts(substring->pre_junction),substring->querystart));
      substring->querystart += Junction_ninserts(substring->pre_junction);
    }
  } else if (Junction_type(substring->pre_junction) == SPLICE_JUNCTION) {
    extraleft = 2;
  } else if (Junction_type(substring->pre_junction) == CHIMERA_JUNCTION) {
    extraleft = 2;
  } else {
    extraleft = 0;
  }

  if (substring->post_junction == NULL) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
#if 1
    /* Extend to the end */
    extraright = path->querylength - substring->queryend;
#else
    if (substring->plusp != invertp && path->splicetype3 != NO_SPLICE) {
      extraright = 2;
    } else if (substring->plusp == invertp && path->splicetype5 != NO_SPLICE) {
      extraright = 2;
    } else {
      extraright = 0;
    }
#endif
  } else if ((type2 = Junction_type(substring->post_junction)) == INS_JUNCTION) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
    extraright = 0;
    if (substring->plusp == invertp) {
      debug0(printf("post-junction is INS, so subtracting %d from queryend %d\n",
		    Junction_ninserts(substring->post_junction),substring->queryend));
      substring->queryend -= Junction_ninserts(substring->post_junction);
    }
  } else if (type2 == DEL_JUNCTION) {
    deletion_length = Junction_nindels(substring->post_junction);
    if (substring->plusp == invertp) {
      deletionpos = substring->univdiagonal - substring->queryend - deletion_length;
    } else {
      deletionpos = substring->univdiagonal - path->querylength + substring->queryend;
    }
    deletion_string = Junction_deletion_string(substring->post_junction,deletionpos);
    extraright = 0;
  } else if (type2 == SPLICE_JUNCTION) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
    extraright = 2;
  } else if (type2 == CHIMERA_JUNCTION) {
    deletion_string = (char *) NULL;
    deletion_length = 0;
    extraright = 2;
  } else {
    deletion_string = (char *) NULL;
    deletion_length = 0;
    extraright = 0;
  }

  print_genomic(fp,path->genomic_diff,substring->querystart,substring->queryend,
		path->querylength,deletion_string,deletion_length,
		/*left*/substring->univdiagonal - path->querylength,
		extraleft,extraright,substring->plusp,path->plusp,queryseq,invertp);

  FREE(deletion_string);
  FPRINTF(fp,"\t");
  print_coordinates(fp,substring->univdiagonal,substring->querystart,substring->queryend,
		    path->querylength,chr,substring->chroffset,substring->plusp,invertp);

  FPRINTF(fp,"\t");

#ifdef DEBUG1
  FPRINTF(fp,"segment_plusp:%d,main_plusp:%d,invertp:%d",substring->plusp,path->plusp,invertp);
  FPRINTF(fp,"\t");
#endif

  if (substring->pre_junction == NULL) {
    type1 = NO_JUNCTION;
    if (sense_forward_p == invertp) {
      /* antisense, so donor is after the intron */
      if (path->plusp != invertp && (path->splicetype5 == DONOR || path->splicetype5 == ANTIDONOR)) {
	FPRINTF(fp,"donor:%.2f",path->ambig_prob_5);
      } else if (path->plusp == invertp && (path->splicetype3 == DONOR || path->splicetype3 == ANTIDONOR)) {
	FPRINTF(fp,"donor:%.2f",path->ambig_prob_3);
      } else if (path->plusp != invertp && (altsplice = path->qstart_alts) != NULL) {
	FPRINTF(fp,"donor:%.2f",altsplice->medial_prob);
      } else if (path->plusp == invertp && (altsplice = path->qend_alts) != NULL) {
	FPRINTF(fp,"donor:%.2f",altsplice->medial_prob);
      } else {
	FPRINTF(fp,"start:%d",substring->querystart);
      }

    } else {
      /* sense, so acceptor is after the intron */
      if (path->plusp != invertp && (path->splicetype5 == ACCEPTOR || path->splicetype5 == ANTIACCEPTOR)) {
	FPRINTF(fp,"acceptor:%.2f",path->ambig_prob_5);
      } else if (path->plusp == invertp && (path->splicetype3 == ACCEPTOR || path->splicetype3 == ANTIACCEPTOR)) {
	FPRINTF(fp,"acceptor:%.2f",path->ambig_prob_3);
      } else if (path->plusp != invertp && (altsplice = path->qstart_alts) != NULL) {
	FPRINTF(fp,"acceptor:%.2f",altsplice->medial_prob);
      } else if (path->plusp == invertp && (altsplice = path->qend_alts) != NULL) {
	FPRINTF(fp,"acceptor:%.2f",altsplice->medial_prob);
      } else {
	FPRINTF(fp,"start:%d",substring->querystart);
      }
    }

  } else if ((type1 = Junction_type(substring->pre_junction)) == INS_JUNCTION) {
    FPRINTF(fp,"ins:%d",Junction_nindels(substring->pre_junction));

  } else if (type1 == DEL_JUNCTION) {
    FPRINTF(fp,"del:%d",Junction_nindels(substring->pre_junction));

  } else if (type1 == SPLICE_JUNCTION || type1 == CHIMERA_JUNCTION) {
    if (sense_forward_p == invertp) {
      /* antisense, so donor is after the intron */
      FPRINTF(fp,"donor:%.2f",Junction_donor_prob(substring->pre_junction));
    } else {
      /* sense, so acceptor is after the intron */
      FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(substring->pre_junction));
    }

  } else {
    abort();
  }

  FPRINTF(fp,"..");

  if (substring->post_junction == NULL) {
    type2 = NO_JUNCTION;
    if (sense_forward_p == invertp) {
      /* antisense, so acceptor is before the intron */
      if (path->plusp != invertp && (path->splicetype3 == ACCEPTOR || path->splicetype3 == ANTIACCEPTOR)) {
	FPRINTF(fp,"acceptor:%.2f",path->ambig_prob_3);
      } else if (path->plusp == invertp && (path->splicetype5 == ACCEPTOR || path->splicetype5 == ANTIACCEPTOR)) {
	FPRINTF(fp,"acceptor:%.2f",path->ambig_prob_5);
      } else if (path->plusp != invertp && (altsplice = path->qend_alts) != NULL) {
	FPRINTF(fp,"acceptor:%.2f",altsplice->medial_prob);
      } else if (path->plusp == invertp && (altsplice = path->qstart_alts) != NULL) {
	FPRINTF(fp,"acceptor:%.2f",altsplice->medial_prob);
      } else {
	FPRINTF(fp,"end:%d",path->querylength - substring->queryend);
      }

    } else {
      /* sense, so donor is before the intron */
      if (path->plusp != invertp && (path->splicetype3 == DONOR || path->splicetype3 == ANTIDONOR)) {
	FPRINTF(fp,"donor:%.2f",path->ambig_prob_3);
      } else if (path->plusp == invertp && (path->splicetype5 == DONOR || path->splicetype5 == ANTIDONOR)) {
	FPRINTF(fp,"donor:%.2f",path->ambig_prob_5);
      } else if (path->plusp != invertp && (altsplice = path->qend_alts) != NULL) {
	FPRINTF(fp,"donor:%.2f",altsplice->medial_prob);
      } else if (path->plusp == invertp && (altsplice = path->qstart_alts) != NULL) {
	FPRINTF(fp,"donor:%.2f",altsplice->medial_prob);
      } else {
	FPRINTF(fp,"end:%d",path->querylength - substring->queryend);
      }
    }

  } else if ((type2 = Junction_type(substring->post_junction)) == INS_JUNCTION) {
    FPRINTF(fp,"ins:%d",Junction_nindels(substring->post_junction));

  } else if (type2 == DEL_JUNCTION) {
    FPRINTF(fp,"del:%d",Junction_nindels(substring->post_junction));

  } else if (type2 == SPLICE_JUNCTION || type2 == CHIMERA_JUNCTION) {
    if (sense_forward_p == invertp) {
      /* antisense, so acceptor is before the intron */
      FPRINTF(fp,"acceptor:%.2f",Junction_acceptor_prob(substring->post_junction));
    } else {
      /* sense, so donor is before the intron */
      FPRINTF(fp,"donor:%.2f",Junction_donor_prob(substring->post_junction));
    }

  } else {
    abort();
  }


#if 0
  FPRINTF(fp,",matches:%d,sub:%d",nmatches,nmismatches);
  if (print_nsnpdiffs_p) {
    FPRINTF(fp,"+%d=%d",substring->nmismatches_refdiff - substring->nmismatches_bothdiff,substring->nmismatches_refdiff);
    if (print_snplabels_p && substring->nmismatches_refdiff > substring->nmismatches_bothdiff) {
      print_snp_labels(fp,substring,queryseq);
    }
  }
#elif defined(PRINT_SUB)
  FPRINTF(fp,",sub:%d",substring->nmismatches);
  if (print_nsnpdiffs_p) {
    FPRINTF(fp,"+%d=%d",substring->ref_nmismatches - substring->nmismatches,substring->ref_nmismatches);
  }
#endif

  if (type1 == SPLICE_JUNCTION && type2 == SPLICE_JUNCTION) {
    if (sense_forward_p == invertp) {
      FPRINTF(fp,",dir:antisense");
    } else {
      FPRINTF(fp,",dir:sense");
    }
    splice_distance_1 = Junction_splice_distance(substring->pre_junction);
    splice_distance_2 = Junction_splice_distance(substring->post_junction);
    if (splice_distance_1 == 0 && splice_distance_2 == 0) {
      /* Skip */
    } else if (splice_distance_1 == 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    } else if (splice_distance_2 == 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
    } else {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    }

  } else if (type1 == SPLICE_JUNCTION) {
    if (sense_forward_p == invertp) {
      FPRINTF(fp,",dir:antisense");
    } else {
      FPRINTF(fp,",dir:sense");
    }
    if ((splice_distance_1 = Junction_splice_distance(substring->pre_junction)) > 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_1:%u",splice_distance_1);
    }

  } else if (type2 == SPLICE_JUNCTION) {
    if (sense_forward_p == invertp) {
      FPRINTF(fp,",dir:antisense");
    } else {
      FPRINTF(fp,",dir:sense");
    }
    if ((splice_distance_2 = Junction_splice_distance(substring->post_junction)) > 0) {
      FPRINTF(fp,",splice_type:consistent");
      FPRINTF(fp,",splice_dist_2:%u",splice_distance_2);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  return;
}


static void
print_pair_info (Filestring_T fp, Pairtype_T pairtype, int insertlength) {

#ifdef TO_FIX
  FPRINTF(fp,"pair_score:%d,",pairscore);
#endif
#ifndef NO_COMPARE
  /* FPRINTF(fp,"insert_length:%d",insertlength); */
#endif
  switch (pairtype) {
  case CONCORDANT: break;
  case CONCORDANT_TRANSLOCATIONS: break;
  case PAIRED_INVERSION: FPRINTF(fp,",pairtype:inversion"); break;
  case PAIRED_SCRAMBLE: FPRINTF(fp,",pairtype:scramble"); break;
  case PAIRED_TOOLONG: FPRINTF(fp,",pairtype:toolong"); break;
  case UNPAIRED: abort();
  case PAIRED_UNSPECIFIED: abort(); /* Used only as an initial pairtype */
  case UNSPECIFIED: abort(); /* Used only as an initial pairtype */
  }

  return;
}


static void
print_substring_list (Filestring_T fp, List_T substrings, Path_T path, Pathpair_T pathpair,
		      Shortread_T queryseq, bool invertp) {
  Substring_T substring;
  List_T p;

  /* First line */
  FPRINTF(fp," ");
  substring = List_head(substrings);
  print_substring(fp,substring,path,queryseq,invertp);

  /* Print alignment info */
  /* FPRINTF(fp,"\tsegs:%d",List_length(substrings)); */
  /* FPRINTF(fp,",align_score:%d,mapq:%d",score,mapq_score); */
  if (method_print_p == true) {
    /* FPRINTF(fp,"\t"); -- Method_print already prints this */
    Method_print(fp,path->method);
  }
    
#if 0
  /* Print pairing info */
  if (pathpair != NULL) {
    FPRINTF(fp,"\t");
    print_pair_info(fp,pathpair->pairtype,pathpair->insertlength);
  }
#endif
  FPRINTF(fp,"\n");


  /* Remaining lines */
  for (p = List_next(substrings); p != NULL; p = List_next(p)) {
    FPRINTF(fp,",");
    substring = List_head(p);
    print_substring(fp,substring,path,queryseq,invertp);
    FPRINTF(fp,"\n");
  }

  return;
}



void
Path_print_alignment (Filestring_T fp, Path_T path, Pathpair_T pathpair,
		      Shortread_T queryseq, bool invertp, Listpool_T listpool) {
  List_T substrings = NULL, fusion_substrings = NULL;
  List_T j;
  Intlist_T q, x, y;
  Univcoordlist_T p;
  Univcoord_T univdiagonal;
  Junction_T pre_junction, post_junction, fusion_junction, temp_junction;
  int querystart, queryend;
  Substring_T substring;
  List_T a;


  debug0(printf("Printing path %p\n",path));
  debug0(Path_print(path));

  /* Process along watson strand */

  if (path->fusion_querystart_junction != NULL) {
    fusion_junction = path->fusion_querystart_junction;
  } else if (path->fusion_queryend_junction != NULL) {
    fusion_junction = path->fusion_queryend_junction;
  } else {
    fusion_junction = (Junction_T) NULL;
  }


  /* 1a.  Make fusion_substrings.  Process along watson strand, but convert
     qstart and qend to querystart and queryend along query
     direction */
  if (fusion_junction != NULL) {
    pre_junction = (Junction_T) NULL;
    for (p = path->fusion_univdiagonals, x = path->fusion_nmismatches, y = path->fusion_ref_nmismatches,
	   q = path->fusion_endpoints, j = path->fusion_junctions;
	 p != NULL; p = Univcoordlist_next(p), x = Intlist_next(x), y = Intlist_next(y),
	   q = Intlist_next(q), j = List_next(j)) {
      univdiagonal = Univcoordlist_head(p);
      post_junction = (j == NULL) ? (Junction_T) NULL : (Junction_T) List_head(j);
      if (path->fusion_plusp != invertp) {
	querystart = Intlist_head(q);
	queryend = Intlist_second_value(q);
      } else {
	querystart = path->querylength - Intlist_second_value(q);
	queryend = path->querylength - Intlist_head(q);
      }
      debug0(printf("Creating pre-fusion substring at %u, %d:%u, %d..%d, plusp:%d\n",
		    univdiagonal,path->fusion_chrnum,univdiagonal - path->fusion_chroffset,
		    querystart,queryend,path->fusion_plusp));
      fusion_substrings = Listpool_push(fusion_substrings,listpool,
					(void *) Substring_new(univdiagonal,/*nmismatches*/Intlist_head(x),
							       /*ref_nmismatches*/Intlist_head(y),querystart,queryend,
							       pre_junction,post_junction,path->fusion_chrnum,
							       path->fusion_chroffset,path->fusion_plusp)
					listpool_trace(__FILE__,__LINE__));
      pre_junction = post_junction;
    }
			   
    /* 1b.  Order fusion substrings */
    if (path->fusion_plusp != invertp) {
      debug0(printf("Putting fusion substrings back into original order\n"));
      fusion_substrings = List_reverse(fusion_substrings);

    } else {
      debug0(printf("Keeping fusion substrings reversed, but swapping all pre- and post-junctions\n"));
      for (a = fusion_substrings; a != NULL; a = List_next(a)) {
	substring = (Substring_T) List_head(a);
	temp_junction = substring->pre_junction;
	substring->pre_junction = substring->post_junction;
	substring->post_junction = temp_junction;
      }
    }
  }


  /* 2a.  Make main substrings.  Process along watson strand, but convert
     qstart and qend to querystart and queryend along query
     direction */
  pre_junction = (Junction_T) NULL;
  for (p = path->univdiagonals, x = path->nmismatches, y = path->ref_nmismatches,
	 q = path->endpoints, j = path->junctions;
       p != NULL; p = Univcoordlist_next(p), x = Intlist_next(x), y = Intlist_next(y),
	 q = Intlist_next(q), j = List_next(j)) {
    univdiagonal = Univcoordlist_head(p);
    post_junction = (j == NULL) ? (Junction_T) NULL : (Junction_T) List_head(j);
    if (path->plusp != invertp) {
      querystart = Intlist_head(q);
      queryend = Intlist_second_value(q);
    } else {
      querystart = path->querylength - Intlist_second_value(q);
      queryend = path->querylength - Intlist_head(q);
    }
    debug0(printf("Creating main substring at %u, %d:%u, %d..%d, plusp:%d\n",
		  univdiagonal,path->chrnum,univdiagonal - path->chroffset,
		  querystart,queryend,path->plusp));
    substrings = Listpool_push(substrings,listpool,
			       (void *) Substring_new(univdiagonal,/*nmismatches*/Intlist_head(x),
						      /*ref_nmismatches*/Intlist_head(y),querystart,queryend,
						      pre_junction,post_junction,path->chrnum,
						      path->chroffset,path->plusp)
			       listpool_trace(__FILE__,__LINE__));
    pre_junction = post_junction;
  }

  /* 2b.  Order main substrings */
  if (path->plusp != invertp) {
    debug0(printf("Putting substrings back into original order\n"));
    substrings = List_reverse(substrings);

  } else {
    debug0(printf("Keeping substrings reversed, but swapping all pre- and post-junctions\n"));
    for (a = substrings; a != NULL; a = List_next(a)) {
      substring = (Substring_T) List_head(a);
      temp_junction = substring->pre_junction;
      substring->pre_junction = substring->post_junction;
      substring->post_junction = temp_junction;
    }
  }


  /* 3. Merge fusion substrings and main substrings */
  if (fusion_substrings == NULL) {
    /* Skip */
  } else if (invertp == false) {
    if (path->fusion_querystart_junction != NULL) {
      /* Put fusion at querystart */
      debug0(printf("Putting fusion at querystart\n"));
      substring = List_last_value(fusion_substrings,NULL);
      substring->post_junction = fusion_junction;
      substring = List_head(substrings);
      substring->pre_junction = fusion_junction;
      substrings = List_append(fusion_substrings,substrings);
      
    } else {
      /* Put fusion at queryend */
      debug0(printf("Putting fusion at queryend\n"));
      substring = List_last_value(substrings,NULL);
      substring->post_junction = fusion_junction;
      substring = List_head(fusion_substrings);
      substring->pre_junction = fusion_junction;
      substrings = List_append(substrings,fusion_substrings);
    }
  } else {
    if (path->fusion_querystart_junction != NULL) {
      /* Put fusion at queryend */
      debug0(printf("Putting fusion at queryend\n"));
      substring = List_last_value(substrings,NULL);
      substring->post_junction = fusion_junction;
      substring = List_head(fusion_substrings);
      substring->pre_junction = fusion_junction;
      substrings = List_append(substrings,fusion_substrings);
    } else {
      /* Put fusion at querystart */
      debug0(printf("Putting fusion at querystart\n"));
      substring = List_last_value(fusion_substrings,NULL);
      substring->post_junction = fusion_junction;
      substring = List_head(substrings);
      substring->pre_junction = fusion_junction;
      substrings = List_append(fusion_substrings,substrings);
    }
  }


  /* 4.  Print the substrings */
  print_substring_list(fp,substrings,path,pathpair,queryseq,invertp);

  for (a = substrings; a != NULL; a = List_next(a)) {
    substring = (Substring_T) List_head(a);
    Substring_free(&substring);
  }

  /* fusion_substrings are appended to substrings */
  Listpool_free_list(&substrings,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */

  return;
}



void
Path_print_alignment_setup (Univ_IIT_T chromosome_iit_in, bool method_print_p_in,
			    bool print_univdiagonal_p_in, bool print_nsnpdiffs_p_in) {

  chromosome_iit = chromosome_iit_in;
  method_print_p = method_print_p_in;
  print_univdiagonal_p = print_univdiagonal_p_in;
  print_nsnpdiffs_p = print_nsnpdiffs_p_in;

  return;
}
