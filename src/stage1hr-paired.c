static char rcsid[] = "$Id: 5d5b09a79cc1f61b88042c59cd2d29a7fa9026f6 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include "stage1hr-single.h"
#include "stage1hr-paired.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For rint */
#include <string.h>		/* For memset */

#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"
#include "univdiag.h"

#include "list.h"
#include "compress.h"
#include "record.h"

#include "genomebits_mismatches.h" /* For MISMATCH_EXTRA */
#include "genomebits_kmer.h" /* For MISMATCH_EXTRA */

#include "stage1hr.h"
#include "kmer-search.h"
#include "extension-search.h"
#include "stage1hr-single.h"
#include "transcriptome-search.h"
#include "tr-extension-search.h"

#include "intersect-wdups-indices.h" /* For taking intersection of trnums */
#ifdef LARGE_GENOMES
#include "merge-uint8.h"
#include "intersect-approx-indices-uint8.h"
#else
#include "merge-uint4.h"
#include "intersect-approx-indices-uint4.h"
#endif

#include "concordance.h"
#include "trpath-solve.h"
#include "trpath-convert.h"
#include "path-solve.h"
#include "path-fusion.h"

#include "transcript-remap.h"
#include "transcript-velocity.h"
#include "path-eval.h"
#include "pathpair-eval.h"

#include "orderstat.h"


/* Mutual concordance is greedy and can miss optimal solutions, but
   needed when there are too many pairs */
#define MAX_OPTIMAL 100

#define LOCALDB_REGION_SIZE 65536
#define QUERYLENGTH_FOR_LOCALDB_MATE 50

/* #define MIN_SIZELIMIT 100 */
#define MAX_HITS_EXACT 100	/* Excessive exact paths indicate a repetitive sequence */
#define MAX_DENSITY_INTERSECTION 100 /* Limit on ndense5 * ndense3 */
#define MAX_NPATHS_FIND_SPLICES 0    /* Merging seems not have much effect */

#define MAX_SINGLE_END_HITS 1000
#define MAX_UNEXTENDED_HITS 1000

#define FILTERING_NHITS 100

#if 0
/* Too slow */
#define USE_ALL_UNIVDIAGONALS 1	/* Needed to obtain correct results, and can speed up alignment */
#endif


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Trdiagonals */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* find_search_mates */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif

/* determine_pairtype */
#ifdef DEBUG14
#define debug14(x) x
#else
#define debug14(x)
#endif

/* consolidate_paired_results */
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif


static Mode_T mode;
static int index1part;
static int index1interval;
static int index1part_tr;

static Transcriptome_T transcriptome;
static bool genome_align_p;
static bool transcriptome_align_p;

static Genomebits_T genomebits;
static Localdb_T localdb;

static EF64_T chromosome_ef64;

/* static double defect_rate = 0.01; */

static double user_nmismatches_filter_float;
static double user_mincoverage_filter_float;

static Chrpos_T positive_gap_distance;
static Chrpos_T concordance_distance;

static bool splicingp;
static int maxpaths_search;	/* Not currently used */
static int maxpaths_report;	/* Not currently used */

static bool *circularp;
static int pairmax_linear;
static int pairmax_circular;


#define T Stage1_T


#if 0
/* Do not need to extend all paths, since fusion procedures will extend fusion candidates */
static void
eval_unextended_paths (int *found_score, T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  List_T p;

  for (p = this->unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),(Path_T) List_head(p),query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),(Path_T) List_head(p),query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),(Path_T) List_head(p),query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),(Path_T) List_head(p),query_compress_fwd,query_compress_rev);
  }

  return;
}
#endif


static List_T
paired_search_trdiagonals (int *found_score_paired, int *found_score_5, int *found_score_3,
			   int sufficient_score_5, int sufficient_score_3,
			   Method_T *last_method_5, Method_T *last_method_3,
		  
			   List_T pathpairs, T this5, T this3, Knownsplicing_T knownsplicing,
	      
			   Shortread_T queryseq5, Shortread_T queryseq3,
			   int querylength5, int querylength3,
	       
			   Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			   Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
	      
			   int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

			   int nmismatches_filter_5, int nmismatches_filter_3,
			   int mincoverage_filter_5, int mincoverage_filter_3,

			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, Method_T method_goal) {
  
  int *sense_indices, *antisense_indices;
  int sense_nindices, antisense_nindices;
  int index1, index2;
  int tstart5, tend5, tstart3, tend3;
  Trpath_T trpath5, trpath3;

  Trnum_T trnum;
  Chrnum_T chrnum;
  int transcript_genestrand;
  bool partialp;

  int i, j, k;

  debug1(printf(">>Entered paired_search_trdiagonals\n"));

  if (*last_method_5 < TR_EXACT1 && *last_method_3 < TR_EXACT1) {
    *last_method_5 = single_read_next_method_trdiagonal(*last_method_5,this5,querylength5,
							query5_compress_fwd,query5_compress_rev,
							/*first_read_p*/true);
    *last_method_3 = single_read_next_method_trdiagonal(*last_method_3,this3,querylength3,
							query3_compress_fwd,query3_compress_rev,
							/*first_read_p*/false);
    
  } else if (*last_method_5 >= method_goal) {
    *last_method_3 = single_read_next_method_trdiagonal(*last_method_3,this3,querylength3,
							query3_compress_fwd,query3_compress_rev,
							/*first_read_p*/false);

  } else if (*last_method_3 >= method_goal) {
    *last_method_5 = single_read_next_method_trdiagonal(*last_method_5,this5,querylength5,
							query5_compress_fwd,query5_compress_rev,
							/*first_read_p*/true);

  } else if ((*found_score_5) >= (*found_score_3)) {
    *last_method_5 = single_read_next_method_trdiagonal(*last_method_5,this5,querylength5,
							query5_compress_fwd,query5_compress_rev,
							/*first_read_p*/true);

  } else {
    *last_method_3 = single_read_next_method_trdiagonal(*last_method_3,this3,querylength3,
							query3_compress_fwd,query3_compress_rev,
							/*first_read_p*/false);
  }

  debug(printf("Have %d and %d sense trdiagonals, %d and %d antisense trdiagonals\n",
	       this5->n_sense_trdiagonals,this3->n_sense_trdiagonals,
	       this5->n_antisense_trdiagonals,this3->n_antisense_trdiagonals));

#ifdef DEBUG1
  for (i = 0; i < this5->n_sense_trdiagonals; i++) {
    printf("%u %u\n",this5->sense_trnums[i],this5->_sense_trdiagonals[i]);
  }
  printf("\n");

  for (i = 0; i < this5->n_antisense_trdiagonals; i++) {
    printf("%u %u\n",this5->antisense_trnums[i],this5->_antisense_trdiagonals[i]);
  }
  printf("\n");

  for (i = 0; i < this3->n_sense_trdiagonals; i++) {
    printf("%u %u\n",this3->sense_trnums[i],this3->_sense_trdiagonals[i]);
  }
  printf("\n");

  for (i = 0; i < this3->n_antisense_trdiagonals; i++) {
    printf("%u %u\n",this3->antisense_trnums[i],this3->_antisense_trdiagonals[i]);
  }
  printf("\n");
#endif


  /* Find concordant trnums */
  if (this5->n_sense_trdiagonals == 0 || this3->n_sense_trdiagonals == 0) {
    sense_indices = (int *) NULL;
    sense_nindices = 0;
  } else {
#if 0
    /* Previously used this because Intersect_indices_uint4 doesn't allow duplicates */
    sense_indices =
      Intersect_approx_indices_uint4(&sense_nindices,
				     this5->sense_trnums,this5->n_sense_trdiagonals,/*diagterm1*/0,
				     this3->sense_trnums,this3->n_sense_trdiagonals,/*diagterm2*/0,
				     /*below_slop*/0,/*above_slop*/0);
#else
    sense_indices =
      Intersect_wdups_indices(&sense_nindices,
			      this5->sense_trnums,this5->n_sense_trdiagonals,
			      this3->sense_trnums,this3->n_sense_trdiagonals);
#endif
  }
  
  if (this5->n_antisense_trdiagonals == 0 || this3->n_antisense_trdiagonals == 0) {
    antisense_indices = (int *) NULL;
    antisense_nindices = 0;
  } else {
#if 0
    /* Previously used this because Intersect_indices_uint4 doesn't allow duplicates */
    antisense_indices =
      Intersect_approx_indices_uint4(&antisense_nindices,
				     this5->antisense_trnums,this5->n_antisense_trdiagonals,/*diagterm1*/0,
				     this3->antisense_trnums,this3->n_antisense_trdiagonals,/*diagterm2*/0,
				     /*below_slop*/0,/*above_slop*/0);
#else
    antisense_indices =
      Intersect_wdups_indices(&antisense_nindices,
			      this5->antisense_trnums,this5->n_antisense_trdiagonals,
			      this3->antisense_trnums,this3->n_antisense_trdiagonals);
#endif
  }

  debug1(printf("Have %d sense concordant and %d antisense concordant\n",
		sense_nindices,antisense_nindices));
    
  tstart5 = 0;
  tend5 = querylength5;
  tstart3 = 0;
  tend3 = querylength3;

  for (i = 0, k = 0; i < sense_nindices; i++, k += 2) {
    index1 = sense_indices[k];
    index2 = sense_indices[k+1];

    /* Previously checked last_method_5 */
    if (this5->sense_tstarts != NULL) {
      tstart5 = this5->sense_tstarts[index1];
      tend5 = this5->sense_tends[index1];
    }

    trnum = this5->sense_trnums[index1];
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						/*trdiagonal*/this5->_sense_trdiagonals[index1],
						tstart5,tend5,trnum,/*troffset*/this5->sense_troffsets[index1],
						/*trhigh*/this5->sense_trhighs[index1],
						/*query_compress_tr*/query5_compress_fwd,/*tplusp*/true,querylength5,
						mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						/*paired_end_p*/true,/*first_read_p*/true,
						this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						*last_method_5)) != NULL) {
      if (partialp == true) {
	this5->partial_sense_trpaths = Hitlist_push(this5->partial_sense_trpaths,hitlistpool,(void *) trpath5
						    hitlistpool_trace(__FILE__,__LINE__));
      } else {
	this5->sense_trpaths = Hitlist_push(this5->sense_trpaths,hitlistpool,(void *) trpath5
					    hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    /* Previously checked last_method_3 */
    if (this3->sense_tstarts != NULL) {
      tstart3 = this3->sense_tstarts[index2];
      tend3 = this3->sense_tends[index2];
    }

    trnum = this3->sense_trnums[index2];
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						/*trdiagonal*/this3->_sense_trdiagonals[index2],
						tstart3,tend3,/*trnum*/this3->sense_trnums[index2],
						/*troffset*/this3->sense_troffsets[index2],
						/*trhigh*/this3->sense_trhighs[index2],
						/*query_compress_tr*/query3_compress_fwd,/*tplusp*/true,querylength3,
						mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						/*paired_end_p*/true,/*first_read_p*/false,
						this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						*last_method_3)) != NULL) {
      if (partialp == true) {
	this3->partial_sense_trpaths = Hitlist_push(this3->partial_sense_trpaths,hitlistpool,(void *) trpath3
						    hitlistpool_trace(__FILE__,__LINE__));
      } else {
	this3->sense_trpaths = Hitlist_push(this3->sense_trpaths,hitlistpool,(void *) trpath3
					    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }
    

  tstart5 = 0;
  tend5 = querylength5;
  tstart3 = 0;
  tend3 = querylength3;

  for (i = 0, k = 0; i < antisense_nindices; i++, k += 2) {
    index1 = antisense_indices[k];
    index2 = antisense_indices[k+1];

    /* Previously checked last_method_5 */
    if (this5->antisense_tstarts != NULL) {
      tstart5 = this5->antisense_tstarts[index1];
      tend5 = this5->antisense_tends[index1];
    }

    trnum = this5->antisense_trnums[index1];
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						/*trdiagonal*/this5->_antisense_trdiagonals[index1],
						tstart5,tend5,trnum,/*troffset*/this5->antisense_troffsets[index1],
						/*trhigh*/this5->antisense_trhighs[index1],
						/*query_compress_tr*/query5_compress_rev,/*tplusp*/false,querylength5,
						mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						/*paired_end_p*/true,/*first_read_p*/true,
						this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						*last_method_5)) != NULL) {
      if (partialp == true) {
	this5->partial_antisense_trpaths = Hitlist_push(this5->partial_antisense_trpaths,hitlistpool,(void *) trpath5
						hitlistpool_trace(__FILE__,__LINE__));
      } else {
	this5->antisense_trpaths = Hitlist_push(this5->antisense_trpaths,hitlistpool,(void *) trpath5
						hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    /* Previously checked last_method_3 */
    if (this3->antisense_tstarts != NULL) {
      tstart3 = this3->antisense_tstarts[index2];
      tend3 = this3->antisense_tends[index2];
    }
    
    trnum = this3->antisense_trnums[index2];
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						/*trdiagonal*/this3->_antisense_trdiagonals[index2],
						tstart3,tend3,trnum,/*troffset*/this3->antisense_troffsets[index2],
						/*trhigh*/this3->antisense_trhighs[index2],
						/*query_compress_tr*/query3_compress_rev,/*tplusp*/false,querylength3,
						mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						/*paired_end_p*/true,/*first_read_p*/false,
						this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						*last_method_3)) != NULL) {
      if (partialp == true) {
	this3->partial_antisense_trpaths = Hitlist_push(this3->partial_antisense_trpaths,hitlistpool,(void *) trpath3
							hitlistpool_trace(__FILE__,__LINE__));
      } else {
	this3->antisense_trpaths = Hitlist_push(this3->antisense_trpaths,hitlistpool,(void *) trpath3
						hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*new*/this5->sense_trpaths,/*new*/this3->sense_trpaths,
			     /*trpaths5*/NULL,/*trpaths3*/NULL,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
			     
			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*new*/this5->antisense_trpaths,/*new*/this3->antisense_trpaths,
			     /*trpaths5*/NULL,/*trpaths3*/NULL,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI);

  /* Comparing found_score_paired with found_score_5 + found_score_3 does not seem to help */
  if (pathpairs != NULL /*&& (*found_score_paired) == (*found_score_5) + (*found_score_3)*/) {
    FREE(antisense_indices);
    FREE(sense_indices);
    return pathpairs;

  } else {
    /* Create trpaths out of remaining trdiagonals */
    /* Pathpair_gc(&pathpairs,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool); */

    /* 5' read, sense */
    tstart5 = 0;
    tend5 = querylength5;

    i = j = 0; k = 0;
    while (i < this5->n_sense_trdiagonals && j < sense_nindices) {
      index1 = sense_indices[k];
      if (i < index1) {
	if (this5->sense_tstarts != NULL) {
	  tstart5 = this5->sense_tstarts[i];
	  tend5 = this5->sense_tends[i];
	}

	trnum = this5->sense_trnums[i];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						    /*trdiagonal*/this5->_sense_trdiagonals[i],
						    tstart5,tend5,trnum,/*troffset*/this5->sense_troffsets[i],
						    /*trhigh*/this5->sense_trhighs[i],
						    /*query_compress_tr*/query5_compress_fwd,/*tplusp*/true,querylength5,
						    mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						    /*paired_end_p*/true,/*first_read_p*/true,
						    this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						    *last_method_5)) != NULL) {
	  if (partialp == true) {
	    this5->partial_sense_trpaths = Hitlist_push(this5->partial_sense_trpaths,hitlistpool,(void *) trpath5
							hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this5->sense_trpaths = Hitlist_push(this5->sense_trpaths,hitlistpool,(void *) trpath5
						hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
	i++;

      } else if (index1 < i) {
	j++; k += 2;

      } else {
	i++; j++; k += 2;
      }
    }

    while (i < this5->n_sense_trdiagonals) {
      if (this5->sense_tstarts != NULL) {
	tstart5 = this5->sense_tstarts[i];
	tend5 = this5->sense_tends[i];
      }

      trnum = this5->sense_trnums[i];
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						  /*trdiagonal*/this5->_sense_trdiagonals[i],
						  tstart5,tend5,trnum,/*troffset*/this5->sense_troffsets[i],
						  /*trhigh*/this5->sense_trhighs[i],
						  /*query_compress_tr*/query5_compress_fwd,/*tplusp*/true,querylength5,
						  mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						  /*paired_end_p*/true,/*first_read_p*/true,
						  this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						  *last_method_5)) != NULL) {
	if (partialp == true) {
	  this5->partial_sense_trpaths = Hitlist_push(this5->partial_sense_trpaths,hitlistpool,(void *) trpath5
					      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  this5->sense_trpaths = Hitlist_push(this5->sense_trpaths,hitlistpool,(void *) trpath5
					      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      i++;
    }
      

    /* 3' read, sense */
    tstart3 = 0;
    tend3 = querylength3;

    i = j = 0; k = 0;
    while (i < this3->n_sense_trdiagonals && j < sense_nindices) {
      index2 = sense_indices[k+1];
      if (i < index2) {
	if (this3->sense_tstarts != NULL) {
	  tstart3 = this3->sense_tstarts[i];
	  tend3 = this3->sense_tends[i];
	}

	trnum = this3->sense_trnums[i];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						    /*trdiagonal*/this3->_sense_trdiagonals[i],
						    tstart3,tend3,trnum,/*troffset*/this3->sense_troffsets[i],
						    /*trhigh*/this3->sense_trhighs[i],
						    /*query_compress_tr*/query3_compress_fwd,/*tplusp*/true,querylength3,
						    mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						    /*paired_end_p*/true,/*first_read_p*/false,
						    this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						    *last_method_3)) != NULL) {
	  if (partialp == true) {
	    this3->partial_sense_trpaths = Hitlist_push(this3->partial_sense_trpaths,hitlistpool,(void *) trpath3
						hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this3->sense_trpaths = Hitlist_push(this3->sense_trpaths,hitlistpool,(void *) trpath3
						hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
	i++;

      } else if (index2 < i) {
	j++; k += 2;

      } else {
	i++; j++; k += 2;
      }
    }

    while (i < this3->n_sense_trdiagonals) {
      if (this3->sense_tstarts != NULL) {
	tstart3 = this3->sense_tstarts[i];
	tend3 = this3->sense_tends[i];
      }

      trnum = this3->sense_trnums[i];
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						  /*trdiagonal*/this3->_sense_trdiagonals[i],
						  tstart3,tend3,trnum,/*troffset*/this3->sense_troffsets[i],
						  /*trhigh*/this3->sense_trhighs[i],
						  /*query_compress_tr*/query3_compress_fwd,/*tplusp*/true,querylength3,
						  mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						  /*paired_end_p*/true,/*first_read_p*/false,
						  this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						  *last_method_3)) != NULL) {
	if (partialp == true) {
	  this3->partial_sense_trpaths = Hitlist_push(this3->partial_sense_trpaths,hitlistpool,(void *) trpath3
						      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  this3->sense_trpaths = Hitlist_push(this3->sense_trpaths,hitlistpool,(void *) trpath3
					      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      i++;
    }


    /* 5' read, antisense */
    tstart5 = 0;
    tend5 = querylength5;

    i = j = 0; k = 0;
    while (i < this5->n_antisense_trdiagonals && j < antisense_nindices) {
      index1 = antisense_indices[k];
      if (i < index1) {
	if (this5->antisense_tstarts != NULL) {
	  tstart5 = this5->antisense_tstarts[i];
	  tend5 = this5->antisense_tends[i];
	}

	trnum = this5->antisense_trnums[i];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						    /*trdiagonal*/this5->_antisense_trdiagonals[i],
						    tstart5,tend5,trnum,/*troffset*/this5->antisense_troffsets[i],
						    /*trhigh*/this5->antisense_trhighs[i],
						    /*query_compress_tr*/query5_compress_rev,/*tplusp*/false,querylength5,
						    mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						    /*paired_end_p*/true,/*first_read_p*/true,
						    this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						    *last_method_5)) != NULL) {
	  if (partialp == true) {
	    this5->partial_antisense_trpaths = Hitlist_push(this5->partial_antisense_trpaths,hitlistpool,(void *) trpath5
						    hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this5->antisense_trpaths = Hitlist_push(this5->antisense_trpaths,hitlistpool,(void *) trpath5
						    hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
	i++;

      } else if (index1 < i) {
	j++; k += 2;

      } else {
	i++; j++; k += 2;
      }
    }

    while (i < this5->n_antisense_trdiagonals) {
      if (this5->antisense_tstarts != NULL) {
	tstart5 = this5->antisense_tstarts[i];
	tend5 = this5->antisense_tends[i];
      }

      trnum = this5->antisense_trnums[i];
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&partialp,sufficient_score_5,
						  /*trdiagonal*/this5->_antisense_trdiagonals[i],
						  tstart5,tend5,trnum,/*troffset*/this5->antisense_troffsets[i],
						  /*trhigh*/this5->antisense_trhighs[i],
						  /*query_compress_tr*/query5_compress_rev,/*tplusp*/false,querylength5,
						  mismatch_positions_alloc_5,chrnum,transcript_genestrand,
						  /*paired_end_p*/true,/*first_read_p*/true,
						  this5->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						  *last_method_5)) != NULL) {
	if (partialp == true) {
	  this5->partial_antisense_trpaths = Hitlist_push(this5->partial_antisense_trpaths,hitlistpool,(void *) trpath5
							  hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  this5->antisense_trpaths = Hitlist_push(this5->antisense_trpaths,hitlistpool,(void *) trpath5
						  hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      i++;
    }
      

    /* 3' read, antisense */
    tstart3 = 0;
    tend3 = querylength3;

    i = j = 0; k = 0;
    while (i < this3->n_antisense_trdiagonals && j < antisense_nindices) {
      index2 = antisense_indices[k+1];
      if (i < index2) {
	if (this3->antisense_tstarts != NULL) {
	  tstart3 = this3->antisense_tstarts[i];
	  tend3 = this3->antisense_tends[i];
	}

	trnum = this3->antisense_trnums[i];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						    /*trdiagonal*/this3->_antisense_trdiagonals[i],
						    tstart3,tend3,trnum,/*troffset*/this3->antisense_troffsets[i],
						    /*trhigh*/this3->antisense_trhighs[i],
						    /*query_compress_tr*/query3_compress_rev,/*tplusp*/false,querylength3,
						    mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						    /*paired_end_p*/true,/*first_read_p*/false,
						    this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						    *last_method_3)) != NULL) {
	  if (partialp == true) {
	    this3->partial_antisense_trpaths = Hitlist_push(this3->partial_antisense_trpaths,hitlistpool,(void *) trpath3
						    hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this3->antisense_trpaths = Hitlist_push(this3->antisense_trpaths,hitlistpool,(void *) trpath3
						    hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
	i++;

      } else if (index2 < i) {
	j++; k += 2;

      } else {
	i++; j++; k += 2;
      }
    }

    while (i < this3->n_antisense_trdiagonals) {
      if (this3->antisense_tstarts != NULL) {
	tstart3 = this3->antisense_tstarts[i];
	tend3 = this3->antisense_tends[i];
      }

      trnum = this3->antisense_trnums[i];
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&partialp,sufficient_score_3,
						  /*trdiagonal*/this3->_antisense_trdiagonals[i],
						  tstart3,tend3,trnum,
						  /*troffset*/this3->antisense_troffsets[i],
						  /*trhigh*/this3->antisense_trhighs[i],
						  /*query_compress_tr*/query3_compress_rev,/*tplusp*/false,querylength3,
						  mismatch_positions_alloc_3,chrnum,transcript_genestrand,
						  /*paired_end_p*/true,/*first_read_p*/false,
						  this3->indelinfo,intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						  *last_method_3)) != NULL) {
	if (partialp == true) {
	  this3->partial_antisense_trpaths = Hitlist_push(this3->partial_antisense_trpaths,hitlistpool,(void *) trpath3
						  hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  this3->antisense_trpaths = Hitlist_push(this3->antisense_trpaths,hitlistpool,(void *) trpath3
						  hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      i++;
    }
    
    FREE(antisense_indices);
    FREE(sense_indices);

    return (List_T) NULL;
  }
}


static List_T
paired_read_next_method_tr_5 (int *found_score_paired, int *found_score_5, int *found_score_3,
			      int sufficient_score_5, Method_T *last_method_5,

			      List_T pathpairs, T this5, T this3, Knownsplicing_T knownsplicing,
			      
			      Shortread_T queryseq5, Shortread_T queryseq3,
			      int querylength5, int querylength3,

			      int *mismatch_positions_alloc_5, int nmismatches_allowed_5,

			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

			      int nmismatches_filter_5, int nmismatches_filter_3,
			      int mincoverage_filter_5, int mincoverage_filter_3,

			      int genestrand, Trdiagpool_T trdiagpool, Intlistpool_T intlistpool,
			      Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Transcriptpool_T transcriptpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T partial_sense_trpaths5, partial_antisense_trpaths5, sense_trpaths5, antisense_trpaths5;

  *last_method_5 = single_read_next_method_tr(&(*found_score_5),sufficient_score_5,*last_method_5,
					     
					      &partial_sense_trpaths5,&partial_antisense_trpaths5,
					      &sense_trpaths5,&antisense_trpaths5,
					     
					      this5,genestrand,querylength5,
					      mismatch_positions_alloc_5,
					      query5_compress_fwd,query5_compress_rev,
					      nmismatches_allowed_5,
					      trdiagpool,intlistpool,uintlistpool,
					      listpool,trpathpool,pathpool,hitlistpool,
					      /*paired_end_p*/true,/*first_read_p*/true,/*appendp*/false);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*new*/sense_trpaths5,/*sense_trpaths3*/NULL,
			     /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD);
  this5->partial_sense_trpaths = List_append(partial_sense_trpaths5,this5->partial_sense_trpaths);
  this5->sense_trpaths = List_append(sense_trpaths5,this5->sense_trpaths);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*new*/antisense_trpaths5,/*antisense_trpaths3*/NULL,
			     /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI);
  this5->partial_antisense_trpaths = List_append(partial_antisense_trpaths5,this5->partial_antisense_trpaths);
  this5->antisense_trpaths = List_append(antisense_trpaths5,this5->antisense_trpaths);

  return pathpairs;
}


static List_T
paired_read_next_method_tr_3 (int *found_score_paired, int *found_score_5, int *found_score_3,
			      int sufficient_score_3, Method_T *last_method_3,

			      List_T pathpairs, T this5, T this3, Knownsplicing_T knownsplicing,
			      
			      Shortread_T queryseq5, Shortread_T queryseq3,
			      int querylength5, int querylength3,

			      int *mismatch_positions_alloc_3, int nmismatches_allowed_3,
			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    
			      int nmismatches_filter_5, int nmismatches_filter_3,
			      int mincoverage_filter_5, int mincoverage_filter_3,

			      int genestrand, Trdiagpool_T trdiagpool, Intlistpool_T intlistpool,
			      Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Transcriptpool_T transcriptpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T partial_sense_trpaths3, partial_antisense_trpaths3, sense_trpaths3, antisense_trpaths3;

  *last_method_3 = single_read_next_method_tr(&(*found_score_3),sufficient_score_3,*last_method_3,
					     
					      &partial_sense_trpaths3,&partial_antisense_trpaths3,
					      &sense_trpaths3,&antisense_trpaths3,
					     
					      this3,genestrand,querylength3,
					      mismatch_positions_alloc_3,
					      query3_compress_fwd,query3_compress_rev,
					      nmismatches_allowed_3,
					      trdiagpool,intlistpool,uintlistpool,
					      listpool,trpathpool,pathpool,hitlistpool,
					      /*paired_end_p*/true,/*first_read_p*/false,/*appendp*/false);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*sense_trpaths5*/NULL,/*new*/sense_trpaths3,
			     /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD);
  this3->partial_sense_trpaths = List_append(partial_sense_trpaths3,this3->partial_sense_trpaths);
  this3->sense_trpaths = List_append(sense_trpaths3,this3->sense_trpaths);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 
			     /*antisense_trpaths5*/NULL,/*new*/antisense_trpaths3,
			     /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI);
  this3->partial_antisense_trpaths = List_append(partial_antisense_trpaths3,this3->partial_antisense_trpaths);
  this3->antisense_trpaths = List_append(antisense_trpaths3,this3->antisense_trpaths);

  return pathpairs;
}


/* set_best_paths_p is set to false when computing from univdiagonals,
   and scanning over entire genome, but true when computing for mates */
static void
solve_univdiagonal_auxinfo (bool *complete_sense_p, bool *complete_antisense_p,
			    int *found_score, Univcoord_T univdiagonal, Auxinfo_T *auxinfo, Pathstore_T pathstore,
		  
			    Shortread_T queryseq,
			    char *queryptr, char *queryuc_ptr, char *queryrc, int querylength,
			    T this, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			    int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
			    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   
			    int localdb_nmismatches_allowed, int genestrand,
			     
			    Auxinfopool_T auxinfopool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			    Spliceendsgen_T spliceendsgen, bool plusp, bool first_read_p, bool lowp,
			    bool set_best_paths_p) {

  if (univdiagonal < (Univcoord_T) querylength) {
    /* Skip */

  } else {
    while (pathstore->solvedp == false && *auxinfo != NULL) {
      debug(printf("Method is %s\n",Method_string((*auxinfo)->method)));
      if ((*auxinfo)->method == KMER_EXACT1) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif
	debug(printf("Solving auxinfo %p using method EXACT1 with %d left univdiags and %d right univdiags\n",
		     *auxinfo,List_length((*auxinfo)->left_univdiags),List_length((*auxinfo)->right_univdiags)));
	Path_solve_exact(&(*found_score),
			 
			 &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,
			 
			 univdiagonal,*auxinfo,querylength,
			 plusp,first_read_p,genestrand,
			 query_compress,query_compress_fwd,query_compress_rev,
			 queryseq,queryuc_ptr,queryrc,
			 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			 /*chrhigh*/pathstore->chrhigh,
			 intlistpool,uintlistpool,univcoordlistpool,
			 listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			 /*method*/KMER_EXACT1);
	if (set_best_paths_p == true) {
	  Pathstore_set_best_paths(pathstore,hitlistpool);
	}
	if (pathstore->complete_sense_paths != NULL) {
	  *complete_sense_p = true;
	  pathstore->solvedp = true;
	}
	if (pathstore->complete_antisense_paths != NULL) {
	  *complete_antisense_p = true;
	  pathstore->solvedp = true;
	}
	
      } else if ((*auxinfo)->method == EXT) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif
	debug(printf("Solving auxinfo %p using method EXT with %d left univdiags and %d right univdiags\n",
		     *auxinfo,List_length((*auxinfo)->left_univdiags),List_length((*auxinfo)->right_univdiags)));
	Path_solve_from_diagonals(&(*found_score),
				  
				  &pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
				  &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,
				  
				  univdiagonal,*auxinfo,queryseq,queryptr,querylength,
				  mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
				  localdb_alloc,this,knownsplicing,knownindels,
				  query_compress,query_compress_fwd,query_compress_rev,
				  /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				  /*chrhigh*/pathstore->chrhigh,
				  plusp,genestrand,localdb_nmismatches_allowed,
				  /*paired_end_p*/true,first_read_p,intlistpool,
				  uintlistpool,univcoordlistpool,listpool,pathpool,
				  transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				  /*method*/EXT,/*find_splices_p*/true);
	if (set_best_paths_p == true) {
	  Pathstore_set_best_paths(pathstore,hitlistpool);
	}
	if (pathstore->complete_sense_paths != NULL) {
	  *complete_sense_p = true;
	  pathstore->solvedp = true;
	}
	if (pathstore->complete_antisense_paths != NULL) {
	  *complete_antisense_p = true;
	  pathstore->solvedp = true;
	}
	
      } else if ((*auxinfo)->method == SEGMENT1) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif
	debug(printf("Solving auxinfo %p using method SEGMENT1 with %d left univdiags and %d right univdiags\n",
		     *auxinfo,List_length((*auxinfo)->left_univdiags),List_length((*auxinfo)->right_univdiags)));
	Path_solve_from_diagonals(&(*found_score),
				  
				  &pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
				  &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,
				  
				  univdiagonal,*auxinfo,queryseq,queryptr,querylength,
				  mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
				  localdb_alloc,this,knownsplicing,knownindels,
				  query_compress,query_compress_fwd,query_compress_rev,
				  /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				  /*chrhigh*/pathstore->chrhigh,
				  plusp,genestrand,localdb_nmismatches_allowed,
				  /*paired_end_p*/true,first_read_p,intlistpool,
				  uintlistpool,univcoordlistpool,listpool,pathpool,
				  transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				  /*method*/SEGMENT1,/*find_splices_p*/true);
	if (set_best_paths_p == true) {
	  Pathstore_set_best_paths(pathstore,hitlistpool);
	}
	if (pathstore->complete_sense_paths != NULL) {
	  *complete_sense_p = true;
	  pathstore->solvedp = true;
	}
	if (pathstore->complete_antisense_paths != NULL) {
	  *complete_antisense_p = true;
	  pathstore->solvedp = true;
	}
	
      } else {
	fprintf(stderr,"Unexpected method %s\n",Method_string((*auxinfo)->method));
	abort();
      }

      *auxinfo = Auxinfo_pop(*auxinfo,univdiagpool,auxinfopool);
    }
  }

#if 0
  /* ? Time-consuming.  Can call Path_extend later */
  /* Extend sense paths */
  if (pathstore->complete_sense_paths == NULL && pathstore->unextended_sense_paths != NULL) {
    /* Not sure if we should iterate on all unextended paths or just one */
    path = (Path_T) List_head(pathstore->unextended_sense_paths);
    if (path->extendedp == true) {
      /* Already called Path_extend and got no complete results */
    } else if ((complete_paths = Path_extend(&(*found_score),&unextended_paths,
					     path,queryseq,queryptr,querylength,
					     mismatch_positions_alloc,novel_univdiags_alloc,
					     novel_diagonals_alloc,localdb_alloc,
					     this,knownsplicing,knownindels,query_compress,
					     query_compress_fwd,query_compress_rev,genestrand,
					     localdb_nmismatches_allowed,/*paired_end_p*/true,lowp,
					     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     univdiagpool,vectorpool,hitlistpool,spliceendsgen,
					     /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      pathstore->complete_sense_paths = List_append(complete_paths,pathstore->complete_sense_paths);
    }
  }

  /* Extend antisense paths */
  if (pathstore->complete_antisense_paths == NULL && pathstore->unextended_antisense_paths != NULL) {
    /* Not sure if we should iterate on all unextended paths or just one */
    path = (Path_T) List_head(pathstore->unextended_antisense_paths);
    if (path->extendedp == true) {
      /* Already called Path_extend and got no complete results */
    } else if ((complete_paths = Path_extend(&(*found_score),&unextended_paths,
					     path,queryseq,queryptr,querylength,
					     mismatch_positions_alloc,novel_univdiags_alloc,
					     novel_diagonals_alloc,localdb_alloc,
					     this,knownsplicing,knownindels,query_compress,
					     query_compress_fwd,query_compress_rev,genestrand,
					     localdb_nmismatches_allowed,/*paired_end_p*/true,lowp,
					     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     univdiagpool,vectorpool,hitlistpool,spliceendsgen,
					     /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      pathstore->complete_antisense_paths = List_append(complete_paths,pathstore->complete_antisense_paths);
    }
  }
#endif

  return;
}


#define GENERATE_ALL_PAIRS 1

static List_T
make_pathpairs (int *found_score_paired, List_T *unresolved_pathpairs, List_T pathpairs,
		Pathstore_T pathstoreL, Pathstore_T pathstoreH,
		Shortread_T queryseqL, Shortread_T queryseqH, bool plusp,

		int nmismatches_filter_5, int nmismatches_filter_3,
		int mincoverage_filter_5, int mincoverage_filter_3,

		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
		bool only_completeL_p, bool only_completeH_p) {

  Pathpair_T pathpair;
  Path_T pathL, pathH;

#ifdef GENERATE_ALL_PAIRS
  List_T pathsL, pathsH, p, q;
#endif  


  debug(printf("Making pathpairs from pathstoreL %d sense, %d antisense, and pathstoreH %d sense, %d antisense\n",
	       List_length(pathstoreL->best_sense_paths),List_length(pathstoreL->best_antisense_paths),
	       List_length(pathstoreH->best_sense_paths),List_length(pathstoreH->best_antisense_paths)));

  /* Make allsense pathpairs */
  if (only_completeL_p == false) {
    pathsL = pathstoreL->best_sense_paths; /* Could be complete or unextended */
  } else if (pathstoreL->complete_sense_p == true) {
    pathsL = pathstoreL->best_sense_paths; /* Complete */
  } else {
    pathsL = (List_T) NULL;
  }
    
  if (only_completeH_p == false) {
    pathsH = pathstoreH->best_sense_paths; /* Could be complete or unextended */
  } else if (pathstoreH->complete_sense_p == true) {
    pathsH = pathstoreH->best_sense_paths; /* Complete */
  } else {
    pathsH = (List_T) NULL;
  }

  debug(printf("Checking %d x %d sense paths for concordance\n",List_length(pathsL),List_length(pathsH)));
  for (p = pathsL; p != NULL; p = List_next(p)) {
    pathL = (Path_T) List_head(p);
    for (q = pathsH; q != NULL; q = List_next(q)) {
      pathH = (Path_T) List_head(q);
      debug(printf("  %u and %u\n",pathL->main_univdiagonal,pathH->main_univdiagonal));

      if (pathL->main_univdiagonal > pathH->main_univdiagonal) {
	/* Skip */
      } else if ((pathpair =
		  Pathpair_new_concordant(&(*unresolved_pathpairs),pathL,pathH,queryseqL,queryseqH,plusp,

					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,

					  intlistpool,univcoordlistpool,listpool,
					  pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*check_inner_p*/true,/*copyLp*/true,/*copyHp*/true)) != NULL) {
	debug0(printf("Pathpair found\n"));
	debug0(Pathpair_print(pathpair));
	if (Pathpair_found_score(pathpair) < *found_score_paired) {
	  *found_score_paired = Pathpair_found_score(pathpair);
	}
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }


  /* Make antisense pathpairs */
  if (only_completeL_p == false) {
    pathsL = pathstoreL->best_antisense_paths; /* Could be complete or unextended */
  } else if (pathstoreL->complete_antisense_p == true) {
    pathsL = pathstoreL->best_antisense_paths; /* Complete */
  } else {
    pathsL = (List_T) NULL;
  }
    
  if (only_completeH_p == false) {
    pathsH = pathstoreH->best_antisense_paths; /* Could be complete or unextended */
  } else if (pathstoreH->complete_antisense_p == true) {
    pathsH = pathstoreH->best_antisense_paths; /* Complete */
  } else {
    pathsH = (List_T) NULL;
  }
    
  debug(printf("Checking %d x %d antisense paths for concordance\n",List_length(pathsL),List_length(pathsH)));
  for (p = pathsL; p != NULL; p = List_next(p)) {
    pathL = (Path_T) List_head(p);
    for (q = pathsH; q != NULL; q = List_next(q)) {
      pathH = (Path_T) List_head(q);

      if (pathL->main_univdiagonal > pathH->main_univdiagonal) {
	/* Skip */
      } else if ((pathpair =
		  Pathpair_new_concordant(&(*unresolved_pathpairs),pathL,pathH,queryseqL,queryseqH,plusp,

					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,

					  intlistpool,univcoordlistpool,listpool,
					  pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*check_inner_p*/true,/*copyLp*/true,/*copyHp*/true)) != NULL) {
	debug0(printf("Pathpair found\n"));
	debug0(Pathpair_print(pathpair));
	if (Pathpair_found_score(pathpair) < *found_score_paired) {
	  *found_score_paired = Pathpair_found_score(pathpair);
	}
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  return pathpairs;
}


#ifndef CHECK_ASSERTIONS
static inline void
check_ascending_list (Intlist_T list) {
  return;
}

static inline void
check_descending_list (Intlist_T list) {
  return;
}

#else

static void
check_ascending_list (Intlist_T list) {
  int prevpos;
  Intlist_T p;

  if (list != NULL) {
    /* printf("Ascending? %s\n",Intlist_to_string(list)); */
    prevpos = Intlist_head(list);
    for (p = Intlist_next(list); p != NULL; p = Intlist_next(p)) {
      if (Intlist_head(p) <= prevpos) {
	printf("Expecting ascending, but got %d <= %d\n",
	       Intlist_head(p),prevpos);
	abort();
      }
      prevpos = Intlist_head(p);
    }
  }
 
  return;
}

static void
check_descending_list (Intlist_T list) {
  int prevpos;
  Intlist_T p;

  if (list != NULL) {
    /* printf("Descending? %s\n",Intlist_to_string(list)); */
    prevpos = Intlist_head(list);
    for (p = Intlist_next(list); p != NULL; p = Intlist_next(p)) {
      if (Intlist_head(p) >= prevpos) {
	printf("Expecting descending, but got %d <= %d\n",
	       Intlist_head(p),prevpos);
	abort();
      }
      prevpos = Intlist_head(p);
    }
  }
 
  return;
}

#endif



static List_T
gather_pathpairs_mutual (bool *completeL_p, bool *completeH_p, int *found_score_paired,
			 List_T *unresolved_pathpairs, List_T pathpairs,

			 int *indices, int nindices,

			 Univcoord_T *univdiagonalsL_array, Pathstore_T *pathstoreL_array,
			 Univcoord_T *univdiagonalsH_array, Pathstore_T *pathstoreH_array,

			 Shortread_T queryseqL, Shortread_T queryseqH,
			 
			 int nmismatches_filter_5, int nmismatches_filter_3,
			 int mincoverage_filter_5, int mincoverage_filter_3,

			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			 Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, bool plusp) {

  int index1, index2, i, k;
  Pathstore_T pathstoreL, pathstoreH;
#ifdef DEBUG
  Univcoord_T univdiagonalL, univdiagonalH;
  Pathstore_T last_pathstoreL, last_pathstoreH;
#endif
  List_T paths;
  int cmp;


  /* Find best partners for pathstoreL.  Go through index pairs in
     reverse genomic order so that partners are ascending */

  debug(printf("\nFinding best partners for L univdiagonals\n"));

  for (i = nindices - 1, k = 2*(nindices - 1); i >= 0; i--, k -= 2) {
    index1 = indices[k];
    index2 = indices[k+1];
    pathstoreL = pathstoreL_array[index1];
    pathstoreH = pathstoreH_array[index2];
    debug(printf("Approximately close univdiagonals are %u and %u\n",
		 univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
    
#if 0
    printf("reverse %d %d univdiagonalL %u, univdiagonalH %u\n",index1,index2,
	   univdiagonalsL_array[index1],univdiagonalsH_array[index2]);
#endif
    
    if ((paths = pathstoreH->best_sense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best sense paths\n",index2,univdiagonalsH_array[index2]));
      
    } else if (pathstoreL->best_sense_partners == NULL) {
      debug(printf("Partner #%d at %u is first\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_sense_partners = Intlistpool_push(NULL,intlistpool,index2
						       intlistpool_trace(__FILE__,__LINE__));
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreH_array[pathstoreL->best_sense_partners->first]->best_sense_paths->first)) < 0) {
      debug(printf("Partner #%d at %u is better\n",index2,univdiagonalsH_array[index2]));
      Intlistpool_free_list(&pathstoreL->best_sense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreL->best_sense_partners = Intlistpool_push(NULL,intlistpool,index2
						       intlistpool_trace(__FILE__,__LINE__));
    } else if (cmp == 0) {
      debug(printf("Partner #%d at %u is the same\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_sense_partners = Intlistpool_push(pathstoreL->best_sense_partners,intlistpool,index2
						       intlistpool_trace(__FILE__,__LINE__));
    }
    
    if ((paths = pathstoreH->best_antisense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best antisense paths\n",index2,univdiagonalsH_array[index2]));
      
    } else if (pathstoreL->best_antisense_partners == NULL) {
      debug(printf("Partner #%d at %u is first\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index2
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreH_array[pathstoreL->best_antisense_partners->first]->best_antisense_paths->first)) < 0) {
      debug(printf("Partner #%d at %u is better\n",index2,univdiagonalsH_array[index2]));
      Intlistpool_free_list(&pathstoreL->best_antisense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreL->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index2
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if (cmp == 0) {
      debug(printf("Partner #%d at %u is the same\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_antisense_partners = Intlistpool_push(pathstoreL->best_antisense_partners,intlistpool,index2
							   intlistpool_trace(__FILE__,__LINE__));
    }
  }
  
#ifdef DEBUG
  for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    univdiagonalL = univdiagonalsL_array[index1];
    pathstoreL = pathstoreL_array[index1];
    
    if (i == 0 || pathstoreL != last_pathstoreL) {
      debug(printf("Best sense partners for L at %u are: %s\n",
		   univdiagonalL,Intlist_to_string(pathstoreL->best_sense_partners)));
      debug(printf("Best antisense partners for L at %u are: %s\n",
		   univdiagonalL,Intlist_to_string(pathstoreL->best_antisense_partners)));
      last_pathstoreL = pathstoreL;
    }
  }	
#endif
  
  
  /* Find best partners for pathstoreH.  Go through index pairs in
     forward genomic order so that partners are descending */
  debug(printf("\nFinding best partners for H univdiagonals\n"));
  
  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    index2 = indices[k+1];
    pathstoreL = pathstoreL_array[index1];
    pathstoreH = pathstoreH_array[index2];
    debug(printf("Approximately close univdiagonals are %u and %u\n",
		 univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
    
#if 0
    printf("forward %d %d univdiagonalL %u, univdiagonalH %u\n",index1,index2,
	   univdiagonalsL_array[index1],univdiagonalsH_array[index2]);
#endif
    
    if ((paths = pathstoreL->best_sense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best sense paths\n",index1,univdiagonalsL_array[index1]));
      
    } else if (pathstoreH->best_sense_partners == NULL) {
      debug(printf("Sense partner #%d at %u is first\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_sense_partners = Intlistpool_push(NULL,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreL_array[pathstoreH->best_sense_partners->first]->best_sense_paths->first)) < 0) {
      debug(printf("Sense partner #%d at %u is better\n",index1,univdiagonalsL_array[index1]));
      Intlistpool_free_list(&pathstoreH->best_sense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreH->best_sense_partners = Intlistpool_push(NULL,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
      
    } else if (cmp == 0) {
      debug(printf("Sense partner #%d at %u is the same\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_sense_partners = Intlistpool_push(pathstoreH->best_sense_partners,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
    }
    
    if ((paths = pathstoreL->best_antisense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best antisense paths\n",index1,univdiagonalsL_array[index1]));
      
    } else if (pathstoreH->best_antisense_partners == NULL) {
      debug(printf("Antisense partner #%d at %u is first\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index1
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreL_array[pathstoreH->best_antisense_partners->first]->best_antisense_paths->first)) < 0) {
      debug(printf("Antisense partner #%d at %u is better\n",index1,univdiagonalsL_array[index1]));
      Intlistpool_free_list(&pathstoreH->best_antisense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreH->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index1
							   intlistpool_trace(__FILE__,__LINE__));
    } else if (cmp == 0) {
      debug(printf("Antisense partner #%d at %u is the same\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_antisense_partners = Intlistpool_push(pathstoreH->best_antisense_partners,intlistpool,index1
							   intlistpool_trace(__FILE__,__LINE__));
    }
  }
  
#ifdef DEBUG
  for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
    index2 = indices[k+1];
    univdiagonalH = univdiagonalsH_array[index2];
    pathstoreH = pathstoreH_array[index2];
    
    if (i == 0 || pathstoreH != last_pathstoreH) {
      debug(printf("Best sense partners for H at %u are: %s\n",
		   univdiagonalH,Intlist_to_string(pathstoreH->best_sense_partners)));
      debug(printf("Best antisense partners for H at %u are: %s\n",
		   univdiagonalH,Intlist_to_string(pathstoreH->best_antisense_partners)));
      last_pathstoreH = pathstoreH;
    }
  }	
#endif
  
  
  /* Mutual solutions.  Faster, but may be too greedy and sometimes does not work */
  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    index2 = indices[k+1];
    pathstoreL = pathstoreL_array[index1];
    pathstoreH = pathstoreH_array[index2];
    
    check_ascending_list(pathstoreL->best_sense_partners);
    check_ascending_list(pathstoreL->best_antisense_partners);
    check_descending_list(pathstoreH->best_sense_partners);
    check_descending_list(pathstoreH->best_antisense_partners);
    
#if 0
    /* Tries all mutually best pairs */
    if ((Intlist_exists_p(pathstoreH->best_sense_partners,index1) == true ||
	 Intlist_exists_p(pathstoreH->best_antisense_partners,index1) == true) &&
	(Intlist_exists_p(pathstoreL->best_sense_partners,index2) == true ||
	 Intlist_exists_p(pathstoreL->best_antisense_partners,index2) == true)) {
      pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				 pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				 transcriptpool,hitlistpool,
				 /*only_completeL_p*/*completeL_p,
				 /*only_completeH_p*/*completeH_p);
    }
#else
    /* Tries closest mutually best pairs */
    if ((Intlist_first_equals_p(pathstoreH->best_sense_partners,index1) == true ||
	 Intlist_first_equals_p(pathstoreH->best_antisense_partners,index1) == true) &&
	(Intlist_first_equals_p(pathstoreL->best_sense_partners,index2) == true ||
	 Intlist_first_equals_p(pathstoreL->best_antisense_partners,index2) == true)) {
      debug(printf("(1) Calling make_pathpairs\n"));
      pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				 pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				 transcriptpool,hitlistpool,
				 /*only_completeL_p*/*completeL_p,
				 /*only_completeH_p*/*completeH_p);
    }
#endif
  }
  
  return pathpairs;
}


static List_T
gather_pathpairs_optimal (bool *completeL_p, bool *completeH_p, int *found_score_paired,
			  List_T *unresolved_pathpairs, List_T pathpairs,

			  int *indices, int nindices,

			  Univcoord_T *univdiagonalsL_array, Pathstore_T *pathstoreL_array,
			  Univcoord_T *univdiagonalsH_array, Pathstore_T *pathstoreH_array,
			  
			  Shortread_T queryseqL, Shortread_T queryseqH,
			  
			  int nmismatches_filter_5, int nmismatches_filter_3,
			  int mincoverage_filter_5, int mincoverage_filter_3,
			  
			  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, bool plusp) {

  int index1, index2, i, k;
  Univcoord_T univdiagonalL, univdiagonalH;
  Pathstore_T pathstoreL, pathstoreH, last_pathstoreL, last_pathstoreH;
  List_T paths;
  int cmp;
  Intlist_T p, q;


  /* Find best partners for pathstoreL.  Go through index pairs in
     forward genomic order so that closest partners are seen first */

  debug(printf("\nFinding best partners for L univdiagonals\n"));

  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    index2 = indices[k+1];
    pathstoreL = pathstoreL_array[index1];
    pathstoreH = pathstoreH_array[index2];
    debug(printf("Approximately close univdiagonals are %u and %u\n",
		 univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
    
#if 0
    printf("reverse %d %d univdiagonalL %u, univdiagonalH %u\n",index1,index2,
	   univdiagonalsL_array[index1],univdiagonalsH_array[index2]);
#endif
    
    if ((paths = pathstoreH->best_sense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best sense paths\n",index2,univdiagonalsH_array[index2]));
      
    } else if (pathstoreL->best_sense_partners == NULL) {
      debug(printf("Partner #%d at %u is first\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_sense_partners = Intlistpool_push(NULL,intlistpool,index2
						       intlistpool_trace(__FILE__,__LINE__));
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreH_array[pathstoreL->best_sense_partners->first]->best_sense_paths->first)) < 0) {
      debug(printf("Partner #%d at %u is better\n",index2,univdiagonalsH_array[index2]));
      Intlistpool_free_list(&pathstoreL->best_sense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreL->best_sense_partners = Intlistpool_push(NULL,intlistpool,index2
						       intlistpool_trace(__FILE__,__LINE__));
    } else if (cmp == 0) {
      debug(printf("Partner #%d at %u is the same\n",index2,univdiagonalsH_array[index2]));
      /* pathstoreL->best_sense_partners = Intlistpool_push(pathstoreL->best_sense_partners,intlistpool,index2
	 intlistpool_trace(__FILE__,__LINE__)); */
    }
    
    if ((paths = pathstoreH->best_antisense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best antisense paths\n",index2,univdiagonalsH_array[index2]));
      
    } else if (pathstoreL->best_antisense_partners == NULL) {
      debug(printf("Partner #%d at %u is first\n",index2,univdiagonalsH_array[index2]));
      pathstoreL->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index2
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreH_array[pathstoreL->best_antisense_partners->first]->best_antisense_paths->first)) < 0) {
      debug(printf("Partner #%d at %u is better\n",index2,univdiagonalsH_array[index2]));
      Intlistpool_free_list(&pathstoreL->best_antisense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreL->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index2
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if (cmp == 0) {
      debug(printf("Partner #%d at %u is the same\n",index2,univdiagonalsH_array[index2]));
      /* pathstoreL->best_antisense_partners = Intlistpool_push(pathstoreL->best_antisense_partners,intlistpool,index2
	 intlistpool_trace(__FILE__,__LINE__)); */
    }
  }
  
#ifdef DEBUG
  for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    univdiagonalL = univdiagonalsL_array[index1];
    pathstoreL = pathstoreL_array[index1];
    
    if (i == 0 || pathstoreL != last_pathstoreL) {
      debug(printf("Best sense partners for L at %u are: %s\n",
		   univdiagonalL,Intlist_to_string(pathstoreL->best_sense_partners)));
      debug(printf("Best antisense partners for L at %u are: %s\n",
		   univdiagonalL,Intlist_to_string(pathstoreL->best_antisense_partners)));
      last_pathstoreL = pathstoreL;
    }
  }	
#endif
  
  
  /* Find best partners for pathstoreH.  Go through index pairs in
     reverse genomic order so that partners are ascending, and we
     see the closest one first */

  debug(printf("\nFinding best partners for H univdiagonals\n"));
  
  for (i = nindices - 1, k = 2*(nindices - 1); i >= 0; i--, k -= 2) {
    index1 = indices[k];
    index2 = indices[k+1];
    pathstoreL = pathstoreL_array[index1];
    pathstoreH = pathstoreH_array[index2];
    debug(printf("Approximately close univdiagonals are %u and %u\n",
		 univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
    
#if 0
    printf("forward %d %d univdiagonalL %u, univdiagonalH %u\n",index1,index2,
	   univdiagonalsL_array[index1],univdiagonalsH_array[index2]);
#endif
    
    if ((paths = pathstoreL->best_sense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best sense paths\n",index1,univdiagonalsL_array[index1]));
      
    } else if (pathstoreH->best_sense_partners == NULL) {
      debug(printf("Sense partner #%d at %u is first\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_sense_partners = Intlistpool_push(NULL,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreL_array[pathstoreH->best_sense_partners->first]->best_sense_paths->first)) < 0) {
      debug(printf("Sense partner #%d at %u is better\n",index1,univdiagonalsL_array[index1]));
      Intlistpool_free_list(&pathstoreH->best_sense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreH->best_sense_partners = Intlistpool_push(NULL,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
      
    } else if (cmp == 0) {
      debug(printf("Sense partner #%d at %u is the same\n",index1,univdiagonalsL_array[index1]));
#ifdef MUTUAL_CONCORDANCE
      pathstoreH->best_sense_partners = Intlistpool_push(pathstoreH->best_sense_partners,intlistpool,index1
						       intlistpool_trace(__FILE__,__LINE__));
#endif
    }
    
    if ((paths = pathstoreL->best_antisense_paths) == NULL) {
      /* Skip */
      debug(printf("Partner #%d at %u has no best antisense paths\n",index1,univdiagonalsL_array[index1]));
      
    } else if (pathstoreH->best_antisense_partners == NULL) {
      debug(printf("Antisense partner #%d at %u is first\n",index1,univdiagonalsL_array[index1]));
      pathstoreH->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index1
							   intlistpool_trace(__FILE__,__LINE__));
      
    } else if ((cmp = Path_best_pair_cmp(paths->first,
					 pathstoreL_array[pathstoreH->best_antisense_partners->first]->best_antisense_paths->first)) < 0) {
      debug(printf("Antisense partner #%d at %u is better\n",index1,univdiagonalsL_array[index1]));
      Intlistpool_free_list(&pathstoreH->best_antisense_partners,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      pathstoreH->best_antisense_partners = Intlistpool_push(NULL,intlistpool,index1
							   intlistpool_trace(__FILE__,__LINE__));
    } else if (cmp == 0) {
      debug(printf("Antisense partner #%d at %u is the same\n",index1,univdiagonalsL_array[index1]));
      /* pathstoreH->best_antisense_partners = Intlistpool_push(pathstoreH->best_antisense_partners,intlistpool,index1
	 intlistpool_trace(__FILE__,__LINE__)); */
    }
  }
  
#ifdef DEBUG
  for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
    index2 = indices[k+1];
    univdiagonalH = univdiagonalsH_array[index2];
    pathstoreH = pathstoreH_array[index2];
    
    if (i == 0 || pathstoreH != last_pathstoreH) {
      debug(printf("Best sense partners for H at %u are: %s\n",
		   univdiagonalH,Intlist_to_string(pathstoreH->best_sense_partners)));
      debug(printf("Best antisense partners for H at %u are: %s\n",
		   univdiagonalH,Intlist_to_string(pathstoreH->best_antisense_partners)));
      last_pathstoreH = pathstoreH;
    }
  }	
#endif
  
  
  /* Optimal solutions for L side */
  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    index1 = indices[k];
    pathstoreL = pathstoreL_array[index1];
    check_ascending_list(pathstoreL->best_sense_partners);
    check_ascending_list(pathstoreL->best_antisense_partners);
    
    if (i == 0 || pathstoreL != last_pathstoreL) {
      p = pathstoreL->best_sense_partners;
      q = pathstoreL->best_antisense_partners;
      
      while (p != NULL && q != NULL) {
	if (Intlist_head(p) < Intlist_head(q)) {
	  index2 = Intlist_head(p);
	} else if (Intlist_head(q) < Intlist_head(p)) {
	  index2 = Intlist_head(q);
	} else {
	  index2 = Intlist_head(p);
	}
	
	pathstoreH = pathstoreH_array[index2];
	debug(printf("(1a) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	
	if (Intlist_head(p) < Intlist_head(q)) {
	  p = Intlist_next(p);
	} else if (Intlist_head(q) < Intlist_head(p)) {
	  q = Intlist_next(q);
	} else {
	  p = Intlist_next(p);
	  q = Intlist_next(q);
	}
      }
      
      while (p != NULL) {
	index2 = Intlist_head(p);
	
	pathstoreH = pathstoreH_array[index2];
	debug(printf("(1b) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	p = Intlist_next(p);
      }
      
      while (q != NULL) {
	index2 = Intlist_head(q);
	
	pathstoreH = pathstoreH_array[index2];
	debug(printf("(1c) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	q = Intlist_next(q);
      }
      
      last_pathstoreL = pathstoreL;
    }
  }
  
  
  /* Optimal solutions for H side */
  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    index2 = indices[k+1];
    pathstoreH = pathstoreH_array[index2];
    check_ascending_list(pathstoreH->best_sense_partners);
    check_ascending_list(pathstoreH->best_antisense_partners);
    
    if (i == 0 || pathstoreH != last_pathstoreH) {
      p = pathstoreH->best_sense_partners;
      q = pathstoreH->best_antisense_partners;
      
      while (p != NULL && q != NULL) {
	if (Intlist_head(p) < Intlist_head(q)) {
	  index1 = Intlist_head(p);
	} else if (Intlist_head(q) < Intlist_head(p)) {
	  index1 = Intlist_head(q);
	} else {
	  index1 = Intlist_head(p);
	}
	
	pathstoreL = pathstoreL_array[index1];
	debug(printf("(1d) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	
	if (Intlist_head(p) < Intlist_head(q)) {
	  p = Intlist_next(p);
	} else if (Intlist_head(q) < Intlist_head(p)) {
	  q = Intlist_next(q);
	} else {
	  p = Intlist_next(p);
	  q = Intlist_next(q);
	}
      }
      
      while (p != NULL) {
	index1 = Intlist_head(p);
	
	pathstoreL = pathstoreL_array[index1];
	debug(printf("(1e) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	p = Intlist_next(p);
      }
      
      while (q != NULL) {
	index1 = Intlist_head(q);
	
	pathstoreL = pathstoreL_array[index1];
	debug(printf("(1f) Calling make_pathpairs for %u and %u\n",
		     univdiagonalsL_array[index1],univdiagonalsH_array[index2]));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstoreL,pathstoreH,queryseqL,queryseqH,plusp,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/*completeL_p,
				   /*only_completeH_p*/*completeH_p);
	q = Intlist_next(q);
      }
      
      last_pathstoreH = pathstoreH;
    }
  }
  
#ifdef DEBUG
  printf("returning %d pathpairs\n",List_length(pathpairs));
  for (List_T p = pathpairs; p != NULL; p = List_next(p)) {
    Pathpair_print((Pathpair_T) List_head(p));
  }
#endif

  return pathpairs;
}


static List_T
concordance_univdiagonals (bool *completeL_p, bool *completeH_p, int *found_score_paired,
			   int *found_score_L, int *found_score_H,
		  
			   List_T *unresolved_pathpairs, List_T pathpairs, T thisL, T thisH,
	      
			   Univcoord_T *univdiagonalsL_array, Auxinfo_T *auxinfoL_array,
			   Pathstore_T *pathstoreL_array, int nunivdiagonalsL_array,
			   Univcoord_T *univdiagonalsH_array, Auxinfo_T *auxinfoH_array,
			   Pathstore_T *pathstoreH_array, int nunivdiagonalsH_array,

			   Shortread_T queryseqL, Shortread_T queryseqH,
			   char *queryptrL, char *queryuc_ptr_L, char *queryrcL, int querylengthL,
			   char *queryptrH, char *queryuc_ptr_H, char *queryrcH, int querylengthH,
			   Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			   int *mismatch_positions_alloc_L, int *mismatch_positions_alloc_H,
			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   Compress_T queryL_compress, Compress_T queryL_compress_fwd, Compress_T queryL_compress_rev,
			   Compress_T queryH_compress, Compress_T queryH_compress_fwd, Compress_T queryH_compress_rev,
			   
			   int localdb_nmismatches_allowed_L, int localdb_nmismatches_allowed_H, int genestrand,
			     
			   int nmismatches_filter_5, int nmismatches_filter_3,
			   int mincoverage_filter_5, int mincoverage_filter_3,

			   Auxinfopool_T auxinfopool, Intlistpool_T intlistpool,
			   Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			   Spliceendsgen_T spliceendsgenL, Spliceendsgen_T spliceendsgenH,
			   bool plusp, bool L_first_read_p, bool H_first_read_p) {

  int *indices, index1, index2;
  int nindices, i, k;
  Univcoord_T univdiagonalL, univdiagonalH;
  Pathstore_T pathstoreL, pathstoreH, last_pathstoreL, last_pathstoreH;
  bool complete_sense_p, complete_antisense_p;
  /* List_T paths; */
  /* int cmp; */


  debug(printf("Entering concordance_univdiagonals with found_scores %d and %d, array sizes %d and %d\n",
	       *found_score_L,*found_score_H,nunivdiagonalsL_array,nunivdiagonalsH_array));

#ifdef DEBUG
  for (i = 0; i < nunivdiagonalsL_array; i++) {
    printf("univdiagonalL %u, auxinfoL ",univdiagonalsL_array[i]);
    Auxinfo_print(auxinfoL_array[i]);
    Pathstore_print(pathstoreL_array[i]);
    printf("\n");
  }
  for (i = 0; i < nunivdiagonalsH_array; i++) {
    printf("univdiagonalH %u, auxinfoH ",univdiagonalsH_array[i]);
    Auxinfo_print(auxinfoH_array[i]);
    Pathstore_print(pathstoreH_array[i]);
    printf("\n");
  }
#endif  

  if (nunivdiagonalsL_array > 0 && nunivdiagonalsH_array > 0) {
#ifdef LARGE_GENOMES
    indices = Intersect_approx_indices_uint8(&nindices,
					     univdiagonalsL_array,nunivdiagonalsL_array,/*diagterm1*/-querylengthL,
					     univdiagonalsH_array,nunivdiagonalsH_array,/*diagterm2*/0,
					     /*below_slop*/0,/*above_slop*/concordance_distance);
#else
    indices = Intersect_approx_indices_uint4(&nindices,
					     univdiagonalsL_array,nunivdiagonalsL_array,/*diagterm1*/-querylengthL,
					     univdiagonalsH_array,nunivdiagonalsH_array,/*diagterm2*/0,
					     /*below_slop*/0,/*above_slop*/concordance_distance);
#endif

    /* Solve only intersecting univdiagonals */
    debug(printf("Before removing multiple links, have %d indices\n",nindices));
  
    /* Solve each side separately to reduce combinatorics */
    /* Solving L paths */
    complete_sense_p = complete_antisense_p = false;
    for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
      index1 = indices[k];
      univdiagonalL = univdiagonalsL_array[index1];
      pathstoreL = pathstoreL_array[index1];

      if (i > 0 && pathstoreL == last_pathstoreL) {
	debug(printf("%d/%d: Already solved L path for univdiagonal %u\n",i,nindices,univdiagonalL));

      } else {
	if (pathstoreL->best_sense_partners != NULL) {
	  Intlistpool_free_list(&pathstoreL->best_sense_partners,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  pathstoreL->best_sense_partners = (Intlist_T) NULL;
	}
	if (pathstoreL->best_antisense_partners != NULL) {
	  Intlistpool_free_list(&pathstoreL->best_antisense_partners,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  pathstoreL->best_antisense_partners = (Intlist_T) NULL;
	}

	debug(printf("%d/%d: Solving L path for univdiagonal %u\n",i,nindices,univdiagonalL));
	solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				   &(*found_score_L),univdiagonalL,&(auxinfoL_array[index1]),pathstoreL,
				   
				   queryseqL,queryptrL,queryuc_ptr_L,queryrcL,querylengthL,
				   thisL,knownsplicing,knownindels,
				   
				   mismatch_positions_alloc_L,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   queryL_compress,queryL_compress_fwd,queryL_compress_rev,
				   
				   localdb_nmismatches_allowed_L,genestrand,
				   auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
				   hitlistpool,spliceendsgenL,plusp,L_first_read_p,/*lowp*/true,
				   /*set_best_paths_p*/false);
	last_pathstoreL = pathstoreL;
      }
    }

    if (complete_sense_p == true && complete_antisense_p == true) {
      *completeL_p = true;
      for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
	index1 = indices[k];
	pathstoreL = pathstoreL_array[index1];
	Pathstore_set_best_sense_paths(pathstoreL,hitlistpool,/*only_complete_p*/true);
	Pathstore_set_best_antisense_paths(pathstoreL,hitlistpool,/*only_complete_p*/true);
      }
    } else if (complete_sense_p == true) {
      *completeL_p = true;
      for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
	index1 = indices[k];
	pathstoreL = pathstoreL_array[index1];
	Pathstore_set_best_sense_paths(pathstoreL,hitlistpool,/*only_complete_p*/true);
      }
    } else if (complete_antisense_p == true) {
      *completeL_p = true;
      for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
	index1 = indices[k];
	pathstoreL = pathstoreL_array[index1];
	Pathstore_set_best_antisense_paths(pathstoreL,hitlistpool,/*only_complete_p*/true);
      }
    } else if (*completeL_p == true) {
      /* Already found complete results, so don't store incomplete ones */
    } else {
      /* *completeL_p = false; -- Keep track over all univdiagonals */
      for (i = 0, /*index1*/k = 0; i < nindices; i++, k += 2) {
	index1 = indices[k];
	pathstoreL = pathstoreL_array[index1];
	Pathstore_set_best_sense_paths(pathstoreL,hitlistpool,/*only_complete_p*/false);
	Pathstore_set_best_antisense_paths(pathstoreL,hitlistpool,/*only_complete_p*/false);
      }
    }


    /* Solving H paths */
    complete_sense_p = complete_antisense_p = false;
    for (i = 0, /*index2*/k = 1; i < nindices; i++, k += 2) {
      index2 = indices[k];
      univdiagonalH = univdiagonalsH_array[index2];
      pathstoreH = pathstoreH_array[index2];

      if (i > 0 && pathstoreH == last_pathstoreH) {
	debug(printf("%d/%d: Already solved H path for univdiagonal %u\n",i,nindices,univdiagonalH));

      } else {
	if (pathstoreH->best_sense_partners != NULL) {
	  Intlistpool_free_list(&pathstoreH->best_sense_partners,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  pathstoreH->best_sense_partners = (Intlist_T) NULL;
	}
	if (pathstoreH->best_antisense_partners != NULL) {
	  Intlistpool_free_list(&pathstoreH->best_antisense_partners,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  pathstoreH->best_antisense_partners = (Intlist_T) NULL;
	}

	debug(printf("%d/%d: Solving H path for univdiagonal %u\n",i,nindices,univdiagonalH));
	solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				   &(*found_score_H),univdiagonalH,&(auxinfoH_array[index2]),pathstoreH,

				   queryseqH,queryptrH,queryuc_ptr_H,queryrcH,querylengthH,
				   thisH,knownsplicing,knownindels,

				   mismatch_positions_alloc_H,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   queryH_compress,queryH_compress_fwd,queryH_compress_rev,

				   localdb_nmismatches_allowed_H,genestrand,
				   auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
				   hitlistpool,spliceendsgenH,plusp,H_first_read_p,/*lowp*/false,
				   /*set_best_paths_p*/false);
	last_pathstoreH = pathstoreH;
      }
    }

    if (complete_sense_p == true && complete_antisense_p == true) {
      *completeH_p = true;
      for (i = 0, /*index2*/k = 1; i < nindices; i++, k += 2) {
	index2 = indices[k];
	pathstoreH = pathstoreH_array[index2];
	Pathstore_set_best_sense_paths(pathstoreH,hitlistpool,/*only_complete_p*/true);
	Pathstore_set_best_antisense_paths(pathstoreH,hitlistpool,/*only_complete_p*/true);
      }
    } else if (complete_sense_p == true) {
      *completeH_p = true;
      for (i = 0, /*index2*/k = 1; i < nindices; i++, k += 2) {
	index2 = indices[k];
	pathstoreH = pathstoreH_array[index2];
	Pathstore_set_best_sense_paths(pathstoreH,hitlistpool,/*only_complete_p*/true);
      }
    } else if (complete_antisense_p == true) {
      *completeH_p = true;
      for (i = 0, /*index2*/k = 1; i < nindices; i++, k += 2) {
	index2 = indices[k];
	pathstoreH = pathstoreH_array[index2];
	Pathstore_set_best_antisense_paths(pathstoreH,hitlistpool,/*only_complete_p*/true);
      }
    } else if (*completeH_p == true) {
      /* Already found complete results, so don't store incomplete ones */
    } else {
      /* *completeH_p = false; -- Keep track over all univdiagonals */
      for (i = 0, /*index2*/k = 1; i < nindices; i++, k += 2) {
	index2 = indices[k];
	pathstoreH = pathstoreH_array[index2];
	Pathstore_set_best_sense_paths(pathstoreH,hitlistpool,/*only_complete_p*/false);
	Pathstore_set_best_antisense_paths(pathstoreH,hitlistpool,/*only_complete_p*/false);
      }
    }

#ifdef DEBUG
    printf("%d approximate pairs\n",nindices);
    for (i = 0, k = 0; i < nindices; i++, k += 2) {
      index1 = indices[k];
      index2 = indices[k+1];
      printf("%d %d\n",index1,index2);
    }
    printf("\n");
#endif

    /* If there are many approximate pairs, use mutual concordance
       algorithm, which avoids computational complexity and excess
       memory usage when there are too many possibilities */
    if (nindices > MAX_OPTIMAL) {
      pathpairs =
	gather_pathpairs_mutual(&(*completeL_p),&(*completeH_p),&(*found_score_paired),
				 &(*unresolved_pathpairs),pathpairs,indices,nindices,
					   
				 univdiagonalsL_array,pathstoreL_array,univdiagonalsH_array,pathstoreH_array,
				 
				 queryseqL,queryseqH,nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 
				 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,plusp);
    } else {
      pathpairs =
	gather_pathpairs_optimal(&(*completeL_p),&(*completeH_p),&(*found_score_paired),
				 &(*unresolved_pathpairs),pathpairs,indices,nindices,
					   
				 univdiagonalsL_array,pathstoreL_array,univdiagonalsH_array,pathstoreH_array,
				 
				 queryseqL,queryseqH,nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 
				 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,plusp);
    }

    FREE(indices);
  }


  return pathpairs;
}


/* _univdiagonals_gplus and _univdiagonals_gminus are aligned */
static void
single_read_gen_univdiagonals (Method_T *last_method, int *querystart, int *queryend,

			       Univcoord_T **_univdiagonals_gplus, Univcoord_T **_univdiagonals_gminus,
			       Auxinfo_T **auxinfo_gplus, Auxinfo_T **auxinfo_gminus,
			       Pathstore_T **pathstores_gplus, Pathstore_T **pathstores_gminus,
			       int *nunivdiagonals_gplus, int *nunivdiagonals_gminus,

			       T this, int genestrand,			       
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       int querylength, EF64_T repetitive_ef64,
			       
			       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			       Univdiagpool_T univdiagpool, Univcoordlistpool_T univcoordlistpool,
			       Listpool_T listpool, bool first_read_p) {

  int total_npositions_plus, total_npositions_minus;

  debug(printf("Entering single_read_gen_univdiagonals with last method %s\n",Method_string(*last_method)));

  if (*last_method < KMER_EXACT1) {
    /* 1. Exact search */
    debug(printf("%s Read: 1.  Running Kmer exact1\n",first_read_p ? "5'" : "3'"));
    Stage1_init_end_gen(&(*querystart),&(*queryend),this,querylength,genestrand);
    Kmer_exact1(&(*_univdiagonals_gplus),&(*auxinfo_gplus),
		&(*pathstores_gplus),&(*nunivdiagonals_gplus),
		&(*_univdiagonals_gminus),&(*auxinfo_gminus),
		&(*pathstores_gminus),&(*nunivdiagonals_gminus),
		this,*querystart,*queryend,querylength,auxinfopool,pathstorepool);

    Pathstore_assign_chrinfo(*_univdiagonals_gplus,*pathstores_gplus,*nunivdiagonals_gplus,querylength);
    Pathstore_assign_chrinfo(*_univdiagonals_gminus,*pathstores_gminus,*nunivdiagonals_gminus,querylength);
    debug(printf("Kmer exact1 search returning %d plus and %d minus univdiagonals\n",
		 *nunivdiagonals_gplus,*nunivdiagonals_gminus));
    *last_method = KMER_EXACT1;
    return;

  } else if (*last_method < EXT) {
    /* 2. Extension search */
    debug(printf("%s Read: 2.  Running Extension search\n",first_read_p ? "5'" : "3'"));
    Stage1_fill_all_oligos_gen(this,querylength,genestrand);
    Extension_search(&(*_univdiagonals_gplus),&(*auxinfo_gplus),
		     &(*pathstores_gplus),&(*nunivdiagonals_gplus),
		     &(*_univdiagonals_gminus),&(*auxinfo_gminus),
		     &(*pathstores_gminus),&(*nunivdiagonals_gminus),
		     
		     this,query_compress_fwd,query_compress_rev,querylength,
		     auxinfopool,pathstorepool,
		     univdiagpool,univcoordlistpool,listpool);

    Pathstore_assign_chrinfo(*_univdiagonals_gplus,*pathstores_gplus,*nunivdiagonals_gplus,querylength);
    Pathstore_assign_chrinfo(*_univdiagonals_gminus,*pathstores_gminus,*nunivdiagonals_gminus,querylength);
    debug(printf("Extension search returning %d plus and %d minus univdiagonals\n",
		 *nunivdiagonals_gplus,*nunivdiagonals_gminus));
    *last_method = EXT;
    return;

  } else if (*last_method < SEGMENT1) {
    /* 3. Segment search */
    debug(printf("%s Read: 3.  Running Segment search\n",first_read_p ? "5'" : "3'"));
    Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				  this,querylength,genestrand);
    Kmer_segment(&(*_univdiagonals_gplus),&(*auxinfo_gplus),
		 &(*pathstores_gplus),&(*nunivdiagonals_gplus),
		 &(*_univdiagonals_gminus),&(*auxinfo_gminus),
		 &(*pathstores_gminus),&(*nunivdiagonals_gminus),
		 this,querylength,repetitive_ef64,univdiagpool,auxinfopool,pathstorepool);

    Pathstore_assign_chrinfo(*_univdiagonals_gplus,*pathstores_gplus,*nunivdiagonals_gplus,querylength);
    Pathstore_assign_chrinfo(*_univdiagonals_gminus,*pathstores_gminus,*nunivdiagonals_gminus,querylength);
    debug(printf("Kmer_segment returning %d plus and %d minus univdiagonals\n",
		 *nunivdiagonals_gplus,*nunivdiagonals_gminus));
    *last_method = SEGMENT1;
    return;

#if 0
  } else if (*last_method < KMER_PREVALENT) {
    /* 3. Prevalent (merging).  Equivalent of segment search */
    debug(printf("%s Read: 3.  Running Kmer prevalent\n",first_read_p ? "5'" : "3'"));
    assert(this->all_oligos_gen_filledp == true); /* From Extension_search */
    Stage1_fill_all_positions_gen(this,querylength,genestrand);

    Kmer_prevalent(&(*_univdiagonals_gplus),&(*auxinfo_gplus),&(*nunivdiagonals_gplus),
		   &(*_univdiagonals_gminus),&(*auxinfo_gminus),&(*nunivdiagonals_gminus),
		   this,querylength);

    debug(printf("Paired: Kmer_prevalent returning %d plus and %d minus prevalent\n",
		 *nunivdiagonals_gplus,*nunivdiagonals_gminus));
    *last_method = KMER_PREVALENT;
    return;
#endif

  } else {
    fprintf(stderr,"No method after SEGMENT1\n");
    abort();
  }
}


/* Follows paired_search_trdiagonals */
static List_T
paired_search_univdiagonals (bool *complete5_p, bool *complete3_p, int *found_score_paired,
			     int *found_score_5, int *found_score_3,
			     Method_T *last_method_5, Method_T *last_method_3,
			     bool *any_imperfect_ends5_p, bool *any_imperfect_ends3_p,
		  
			     List_T *unresolved_pathpairs, List_T pathpairs, T this5, T this3,
			     
			     Shortread_T queryseq5, Shortread_T queryseq3,
			     char *queryuc_ptr_5, char *queryrc5, int querylength5,
			     char *queryuc_ptr_3, char *queryrc3, int querylength3,
			     Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			     int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
			     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
	      
			     int nmismatches_allowed_5, int nmismatches_allowed_3, int genestrand,
			     
			     int nmismatches_filter_5, int nmismatches_filter_3,
			     int mincoverage_filter_5, int mincoverage_filter_3,

			     EF64_T repetitive_ef64,
			     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			     Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
			     Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			     Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {
  
  Method_T method_goal = SEGMENT1;
  int kmer5_querystart, kmer5_queryend, kmer3_querystart, kmer3_queryend;

  Univcoord_T *_univdiagonals5_gplus = NULL, *_univdiagonals5_gminus = NULL, *_univdiagonals3_gplus = NULL, *_univdiagonals3_gminus = NULL;
  Auxinfo_T *auxinfo5_gplus = NULL, *auxinfo5_gminus = NULL, *auxinfo3_gplus = NULL, *auxinfo3_gminus = NULL;
  Pathstore_T *pathstores5_gplus = NULL, *pathstores5_gminus = NULL, *pathstores3_gplus = NULL, *pathstores3_gminus = NULL;
  int nunivdiagonals5_gplus = 0, nunivdiagonals5_gminus = 0, nunivdiagonals3_gplus = 0, nunivdiagonals3_gminus = 0;
  int nnovel5_gplus, nnovel5_gminus, nnovel3_gplus, nnovel3_gminus;
  bool run5p = false, run3p = false;

  /* int ndense5_gplus, ndense3_gplus, ndense5_gminus, ndense3_gminus; */


  /* Actually want to provide pair distance */
  debug(printf("\n>>Entered paired_search_univdiagonals with found scores %d and %d, methods %s and %s, nunivdiagonals %d,%d and %d,%d\n",
	       *found_score_5,*found_score_3,Method_string(*last_method_5),Method_string(*last_method_3),
	       this5->all_nunivdiagonals_gplus,this5->all_nunivdiagonals_gminus,
	       this3->all_nunivdiagonals_gplus,this3->all_nunivdiagonals_gminus));
	       
  /* Since we are now solving paths only when we find concordance, we won't know what found_scores are */
  if (*last_method_5 < KMER_EXACT1 && *last_method_3 < KMER_EXACT1
      /* && this5->all_nunivdiagonals_gplus == 0 && this5->all_nunivdiagonals_gminus == 0 */
      /* && this3->all_nunivdiagonals_gplus == 0 && this3->all_nunivdiagonals_gminus == 0 */
      && (*found_score_5) == querylength5 && (*found_score_3) == querylength3) {
    /* Previously ran both only when there are no old univdiagonals
       (from transcriptome search), so we need only one call to a
       concordance procedure */
    /* Now running both when we could not find any trpath */
    run5p = true;
    run3p = true;
  } else if (*last_method_5 >= method_goal) {
    run3p = true;
  } else if (*last_method_3 >= method_goal) {
    run5p = true;
  } else if ((*found_score_5) < querylength5 || (*found_score_3) < querylength3) {
    /* Must created some paths already from concordance, so use that information */
    if ((*found_score_5) > (*found_score_3)) {
      run5p = true;
    } else if ((*found_score_3) > (*found_score_5)) {
      run3p = true;
    } else if ((*last_method_5) < (*last_method_3)) {
      run5p = true;
    } else if ((*last_method_3) < (*last_method_5)) {
      run3p = true;
    } else if (querylength5 > querylength3) {
      run5p = true;
    } else if (querylength3 > querylength5) {
      run3p = true;
    } else {
      run5p = true;
    }
  } else {
    /* Still looking for concordant univdiagonals */
    if ((*last_method_5) < (*last_method_3)) {
      run5p = true;
    } else if ((*last_method_3) < (*last_method_5)) {
      run3p = true;
    } else if (querylength5 > querylength3) {
      run5p = true;
    } else if (querylength3 > querylength5) {
      run3p = true;
    } else {
      run5p = true;
    }
  }
  debug(printf("run5p %d, run3p %d\n",run5p,run3p));

  if (run5p == true) {
    single_read_gen_univdiagonals(&(*last_method_5),&kmer5_querystart,&kmer5_queryend,
				  &_univdiagonals5_gplus,&_univdiagonals5_gminus,
				  &auxinfo5_gplus,&auxinfo5_gminus,
				  &pathstores5_gplus,&pathstores5_gminus,
				  &nunivdiagonals5_gplus,&nunivdiagonals5_gminus,
				  this5,genestrand,query5_compress_fwd,query5_compress_rev,
				  querylength5,repetitive_ef64,
				  auxinfopool,pathstorepool,
				  univdiagpool,univcoordlistpool,listpool,
				  /*first_read_p*/true);
  }
  if (run3p == true) {
    single_read_gen_univdiagonals(&(*last_method_3),&kmer3_querystart,&kmer3_queryend,
				  &_univdiagonals3_gplus,&_univdiagonals3_gminus,
				  &auxinfo3_gplus,&auxinfo3_gminus,
				  &pathstores3_gplus,&pathstores3_gminus,
				  &nunivdiagonals3_gplus,&nunivdiagonals3_gminus,
				  this3,genestrand,query3_compress_fwd,query3_compress_rev,
				  querylength3,repetitive_ef64,
				  auxinfopool,pathstorepool,
				  univdiagpool,univcoordlistpool,listpool,
				  /*first_read_p*/false);
  }


#ifdef DEBUG
  printf("After methods %s and %s, have %d and %d plus univdiagonals, %d and %d minus univdiagonals\n",
	 Method_string(*last_method_5),Method_string(*last_method_3),
	 nunivdiagonals5_gplus,nunivdiagonals3_gplus,
	 nunivdiagonals5_gminus,nunivdiagonals3_gminus);

  printf("5' plus:");
  for (int i = 0; i < nunivdiagonals5_gplus; i++) {
    printf(" %u %s",_univdiagonals5_gplus[i],Method_string(auxinfo5_gplus[i]->method));
  }
  printf("\n");

  printf("3' plus:");
  for (int i = 0; i < nunivdiagonals3_gplus; i++) {
    printf(" %u %s",_univdiagonals3_gplus[i],Method_string(auxinfo3_gplus[i]->method));
  }
  printf("\n");
    

  printf("5' minus:");
  for (int i = 0; i < nunivdiagonals5_gminus; i++) {
    printf(" %u %s",_univdiagonals5_gminus[i],Method_string(auxinfo5_gminus[i]->method));
  }
  printf("\n");

  printf("3' minus:");
  for (int i = 0; i < nunivdiagonals3_gminus; i++) {
    printf(" %u %s",_univdiagonals3_gminus[i],Method_string(auxinfo3_gminus[i]->method));
  }
  printf("\n");
#endif


  /* Identify novel univdiagonals, so we can reuse results from the old Auxinfo_T objects */
  if (run5p == false) {
    nnovel5_gplus = nnovel5_gminus = 0;
  } else {
    nnovel5_gplus = Path_setdiff_univdiagonals_auxinfo(this5->all_univdiagonals_gplus,this5->all_auxinfo_gplus,
						       this5->all_pathstores_gplus,this5->all_nunivdiagonals_gplus,
						       _univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,
						       nunivdiagonals5_gplus,pathstorepool,intlistpool);
    nnovel5_gminus = Path_setdiff_univdiagonals_auxinfo(this5->all_univdiagonals_gminus,this5->all_auxinfo_gminus,
							this5->all_pathstores_gminus,this5->all_nunivdiagonals_gminus,
							_univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,
							nunivdiagonals5_gminus,pathstorepool,intlistpool);
  }

  if (run3p == false) {
    nnovel3_gplus = nnovel3_gminus = 0;
  } else {
    nnovel3_gplus = Path_setdiff_univdiagonals_auxinfo(this3->all_univdiagonals_gplus,this3->all_auxinfo_gplus,
						       this3->all_pathstores_gplus,this3->all_nunivdiagonals_gplus,
						       _univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,
						       nunivdiagonals3_gplus,pathstorepool,intlistpool);
    nnovel3_gminus = Path_setdiff_univdiagonals_auxinfo(this3->all_univdiagonals_gminus,this3->all_auxinfo_gminus,
							this3->all_pathstores_gminus,this3->all_nunivdiagonals_gminus,
							_univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,
							nunivdiagonals3_gminus,pathstorepool,intlistpool);
  }

#ifdef CHECK_ASSERTIONS
  int ncases = 0;
  if (nnovel5_gplus > 0 && nnovel3_gplus > 0) {
    ncases++;
  } else if (nnovel5_gplus > 0 && this3->all_nunivdiagonals_gplus > 0) {
    ncases++;
  } else if (nnovel3_gplus > 0 && this5->all_nunivdiagonals_gplus > 0) {
    ncases++;
  }
  assert(ncases <= 1);

  ncases = 0;
  if (nnovel5_gminus > 0 && nnovel3_gminus > 0) {
    ncases++;
  } else if (nnovel5_gminus > 0 && this3->all_nunivdiagonals_gminus > 0) {
    ncases++;
  } else if (nnovel3_gminus > 0 && this5->all_nunivdiagonals_gminus > 0) {
    ncases++;
  }
  assert(ncases <= 1);
#endif


  if (run5p == true && run3p == true) {
    /* Novel5 vs Novel3 plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),
					  &(*found_score_paired),&(*found_score_5),&(*found_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/_univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,nnovel5_gplus,
					  /*H*/_univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,nnovel3_gplus,

					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,
					  knownsplicing,knownindels,
			      
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
			   
					  nmismatches_allowed_5,nmismatches_allowed_3,genestrand,
			     
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false);

    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/_univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,nnovel3_gminus,
					  /*H*/_univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,nnovel5_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  
					  nmismatches_allowed_3,nmismatches_allowed_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true);

  } else if (run5p == true) {
    /* Novel5 vs Old3, plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),&(*found_score_paired),
					  &(*found_score_5),&(*found_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/_univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,nnovel5_gplus,
					  /*H*/this3->all_univdiagonals_gplus,this3->all_auxinfo_gplus,
					  this3->all_pathstores_gplus,this3->all_nunivdiagonals_gplus,

					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					  
					  nmismatches_allowed_5,nmismatches_allowed_3,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false);

    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/this3->all_univdiagonals_gminus,this3->all_auxinfo_gminus,
					  this3->all_pathstores_gminus,this3->all_nunivdiagonals_gminus,
					  /*H*/_univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,nnovel5_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  
					  nmismatches_allowed_3,nmismatches_allowed_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true);
  } else if (run3p == true) {
    /* Old5 vs Novel3, plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),
					  &(*found_score_paired),&(*found_score_5),&(*found_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/this5->all_univdiagonals_gplus,this5->all_auxinfo_gplus,
					  this5->all_pathstores_gplus,this5->all_nunivdiagonals_gplus,
					  /*H*/_univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,nnovel3_gplus,
					  
					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					  
					  nmismatches_allowed_5,nmismatches_allowed_3,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false);
    
    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/_univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,nnovel3_gminus,
					  /*H*/this5->all_univdiagonals_gminus,this5->all_auxinfo_gminus,
					  this5->all_pathstores_gminus,this5->all_nunivdiagonals_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  
					  nmismatches_allowed_3,nmismatches_allowed_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true);
  }


  /* Always merge information, since we may return to this procedure even if pathpairs is not NULL */
  Path_merge_univdiagonals_auxinfo(&this5->all_univdiagonals_gplus,&this5->all_auxinfo_gplus,
				   &this5->all_pathstores_gplus,&this5->all_nunivdiagonals_gplus,
				   _univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,nnovel5_gplus,
				   pathstorepool,intlistpool);
  Path_merge_univdiagonals_auxinfo(&this5->all_univdiagonals_gminus,&this5->all_auxinfo_gminus,
				   &this5->all_pathstores_gminus,&this5->all_nunivdiagonals_gminus,
				   _univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,nnovel5_gminus,
				   pathstorepool,intlistpool);
  
  Path_merge_univdiagonals_auxinfo(&this3->all_univdiagonals_gplus,&this3->all_auxinfo_gplus,
				   &this3->all_pathstores_gplus,&this3->all_nunivdiagonals_gplus,
				   _univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,nnovel3_gplus,
				   pathstorepool,intlistpool);
  Path_merge_univdiagonals_auxinfo(&this3->all_univdiagonals_gminus,&this3->all_auxinfo_gminus,
				   &this3->all_pathstores_gminus,&this3->all_nunivdiagonals_gminus,
				   _univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,nnovel3_gminus,
				   pathstorepool,intlistpool);

  return pathpairs;

#if 0
  if (pathpairs != NULL) {
    /* Also check for sufficient scores, where we might want to merge univdiagonals */
    /* TODO: Keep the last set of univdiagonals and auxinfo, in case
       Pathpair_eval_and_sort fails and we want to merge them later */
    FREE_ALIGN(_univdiagonals5_gplus);
    FREE_ALIGN(_univdiagonals5_gminus);
    FREE_ALIGN(_univdiagonals3_gplus);
    FREE_ALIGN(_univdiagonals3_gminus);

    Auxinfo_gc(auxinfo5_gplus,nnovel5_gplus,univdiagpool,auxinfopool);
    Auxinfo_gc(auxinfo5_gminus,nnovel5_gminus,univdiagpool,auxinfopool);
    Auxinfo_gc(auxinfo3_gplus,nnovel3_gplus,univdiagpool,auxinfopool);
    Auxinfo_gc(auxinfo3_gminus,nnovel3_gminus,univdiagpool,auxinfopool);

    Pathstore_gc(pathstores5_gplus,nnovel5_gplus,pathstorepool,
		 intlistpool,univcoordlistpool,listpool,pathpool,
		 transcriptpool,hitlistpool);
    Pathstore_gc(pathstores5_gminus,nnovel5_gminus,pathstorepool,
		 intlistpool,univcoordlistpool,listpool,pathpool,
		 transcriptpool,hitlistpool);
    Pathstore_gc(pathstores3_gplus,nnovel3_gplus,pathstorepool,
		 intlistpool,univcoordlistpool,listpool,pathpool,
		 transcriptpool,hitlistpool);
    Pathstore_gc(pathstores3_gminus,nnovel3_gminus,pathstorepool,
		 intlistpool,univcoordlistpool,listpool,pathpool,
		 transcriptpool,hitlistpool);

    return pathpairs;

  } else {
    /* Stage1_list_all_univdiagonals(this3); */

    return pathpairs;
  }
#endif
}



#if 0
static List_T
find_inner_fusions_path5 (int *found_score_5, List_T pathpairs, Path_T path5, List_T singlepaths3,
			  Stage1_T this5, Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Shortread_T queryseq5, Shortread_T queryseq3, char *queryuc_ptr_5, char *queryrc5,
			  Univdiag_T *novel_univdiags_alloc,Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			  Knownsplicing_T knownsplicing, int nmismatches_allowed_5,
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T p;
  Path_T newpath, path3;
  Pathpair_T pathpair;

  for (p = singlepaths3; p != NULL; p = List_next(p)) {
    path3 = (Path_T) List_head(p); /* anchor */
    /* path3 determines overall plusp, and query_compress and queryptr for main part of path5 */
    if (Path_unextended_querystart_p(path5,/*endtrim_allowed*/8,/*allow_ambig_p*/true) == true) {
      /* Skip.  Want fusion alignments to extend to the end */
    } else if (path3->plusp == true) {
      /* Case 1: fusion5 plus, main5 plus; anchor3 plus */
      /* Case 2: fusion5 minus, main5 plus; anchor3 plus */
      if ((newpath = Path_fusion_inner_qend(&(*found_score_5),/*fusion5*/path5,/*anchor3*/path3,
					    /*queryptr_main*/queryuc_ptr_5,/*main_plusp*/true,path5->querylength,
					    /*main_chrnum*/path3->chrnum,path3->chroffset,path3->chrhigh,
					    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this5,
					    this5->spliceinfo,knownsplicing,
					    /*query_compress_main*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					    queryseq5,path5->genestrand,nmismatches_allowed_5,
					    intlistpool,uintlistpool,univcoordlistpool,listpool,
					    pathpool,vectorpool,transcriptpool,hitlistpool)) != NULL &&
	  (pathpair = Pathpair_new_inner_fusion(/*pathL*/newpath,/*pathH*/path3,
						/*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
						intlistpool,univcoordlistpool,listpool,pathpool,
						vectorpool,transcriptpool,hitlistpool,
						/*copyLp*/false,/*copyHp*/true)) != NULL) {
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    } else {
      /* Case 4': anchor3 minus; main5 minus, fusion5 plus */
      /* Case 3': anchor3 minus; main5 minus, fusion5 minus */
      if ((newpath = Path_fusion_inner_qend(&(*found_score_5),/*fusion5*/path5,/*anchor3*/path3,
					    /*queryptr_main*/queryrc5,/*main_plusp*/false,path5->querylength,
					    /*main_chrnum*/path3->chrnum,path3->chroffset,path3->chrhigh,
					    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this5,
					    this5->spliceinfo,knownsplicing,
					    /*query_compress_main*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					    queryseq5,path5->genestrand,nmismatches_allowed_5,
					    intlistpool,uintlistpool,univcoordlistpool,listpool,
					    pathpool,vectorpool,transcriptpool,hitlistpool)) != NULL &&
	  (pathpair = Pathpair_new_inner_fusion(/*pathL*/path3,/*pathH*/newpath,
						/*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
						intlistpool,univcoordlistpool,listpool,pathpool,
						vectorpool,transcriptpool,hitlistpool,
						/*copyLp*/true,/*copyHp*/false)) != NULL) {
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  return pathpairs;
}
#endif


#if 0
static List_T
find_inner_fusions_path3 (int *found_score_3, List_T pathpairs, Path_T path3, List_T singlepaths5,
			  Stage1_T this3, Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  Shortread_T queryseq5, Shortread_T queryseq3, char *queryuc_ptr_3, char *queryrc3,
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			  Knownsplicing_T knownsplicing, int nmismatches_allowed_3,
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T p;
  Path_T newpath, path5;
  Pathpair_T pathpair;

  for (p = singlepaths5; p != NULL; p = List_next(p)) {
    path5 = (Path_T) List_head(p);
    /* path5 determines query_compress and queryptr for main part of path3 */
    if (Path_unextended_queryend_p(path3,/*endtrim_allowed*/8,/*allow_ambig_p*/true) == true) {
      /* Skip.  Want fusion alignments to extend to the end */
    } else if (path5->plusp == true) {
      /* Case 3: anchor5 plus; main3 plus, fusion3 plus */
      /* Case 4: anchor5 plus; main3 plus, fusion3 minus */
      if ((newpath = Path_fusion_inner_qstart(&(*found_score_3),/*fusion3*/path3,/*anchor5*/path5,
					      /*queryptr_main*/queryuc_ptr_3,/*main_plusp*/true,path3->querylength,
					      /*main_chrnum*/path5->chrnum,path5->chroffset,path5->chrhigh,
					      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this3,
					      this3->spliceinfo,knownsplicing,
					      /*query_compress_main*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					      queryseq3,path3->genestrand,nmismatches_allowed_3,
					      intlistpool,uintlistpool,univcoordlistpool,listpool,
					      pathpool,vectorpool,transcriptpool,hitlistpool)) != NULL &&
	  (pathpair = Pathpair_new_inner_fusion(/*pathL*/path5,/*pathH*/newpath,
						/*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
						intlistpool,univcoordlistpool,listpool,pathpool,
						vectorpool,transcriptpool,hitlistpool,
						/*copyLp*/true,/*copyHp*/false)) != NULL) {
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }

    } else {
      /* Case 2': fusion3 plus, main3 minus; anchor5 minus */
      /* Case 1': fusion3 minus, main3 minus; anchor5 minus */
      if ((newpath = Path_fusion_inner_qstart(&(*found_score_3),/*fusion3*/path3,/*anchor5*/path5,
					      /*queryptr_main*/queryrc3,/*main_plusp*/false,path3->querylength,
					      /*main_chrnum*/path5->chrnum,path5->chroffset,path5->chrhigh,
					      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this3,
					      this3->spliceinfo,knownsplicing,
					      /*query_compress_main*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					      queryseq3,path3->genestrand,nmismatches_allowed_3,
					      intlistpool,uintlistpool,univcoordlistpool,listpool,
					      pathpool,vectorpool,transcriptpool,hitlistpool)) != NULL &&
	  (pathpair = Pathpair_new_inner_fusion(/*pathL*/newpath,/*pathH*/path5,
						/*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
						intlistpool,univcoordlistpool,listpool,pathpool,
						vectorpool,transcriptpool,hitlistpool,
						/*copyLp*/false,/*copyHp*/true)) != NULL) {
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  return pathpairs;
}
#endif


#ifdef DEBUG
static void
print_pathpairs_contents (List_T pathpairs) {
  List_T p;
  Pathpair_T pathpair;

  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (Pathpair_T) List_head(p);
    printf("%p %p pathpair\n",pathpair->path5,pathpair->path3);
  }

  return;
}
#endif


static void
create_univdiagonals_auxinfo_from_tr (int *found_score, Stage1_T this, Knownsplicing_T knownsplicing,
				      Shortread_T queryseq, int querylength,
				      List_T sense_trpaths, List_T antisense_trpaths,
				      Compress_T query_compress_fwd, Compress_T query_compress_rev,

				      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
				      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
				      Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
				      Hitlistpool_T hitlistpool, bool first_read_p) {

  List_T unsolved_sense_pathlist_gplus = NULL, unsolved_sense_pathlist_gminus = NULL,
    unsolved_antisense_pathlist_gplus = NULL, unsolved_antisense_pathlist_gminus = NULL;
  List_T sense_pathlist_gplus = NULL, sense_pathlist_gminus = NULL,
    antisense_pathlist_gplus = NULL, antisense_pathlist_gminus = NULL;
  Path_T *sense_paths_gplus, *sense_paths_gminus, *antisense_paths_gplus, *antisense_paths_gminus, path;
  int sense_npaths_gplus, sense_npaths_gminus, antisense_npaths_gplus, antisense_npaths_gminus, npaths;
  int i, j, ii, jj, k, l;

  Univcoord_T univdiagonal_i, univdiagonal_j;
  Pathstore_T pathstore;


  /* Convert sense_trpaths */
  Trpath_convert_sense(&(*found_score),
		       &sense_pathlist_gplus,&sense_pathlist_gminus,
		       &antisense_pathlist_gplus,&antisense_pathlist_gminus,

		       sense_trpaths,first_read_p,
		       queryseq,querylength,this,knownsplicing,
		       query_compress_fwd,query_compress_rev,
		       intlistpool,uintlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool,
		       /*crossover_sense_p*/true);

  /* Convert antisense_trpaths */
  Trpath_convert_antisense(&(*found_score),
			   &antisense_pathlist_gplus,&antisense_pathlist_gminus,
			   &sense_pathlist_gplus,&sense_pathlist_gminus,

			   antisense_trpaths,first_read_p,
			   queryseq,querylength,this,knownsplicing,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,uintlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool,
			   /*crossover_sense_p*/true);

  sense_pathlist_gplus = List_append(unsolved_sense_pathlist_gplus,sense_pathlist_gplus);
  sense_pathlist_gminus = List_append(unsolved_sense_pathlist_gminus,sense_pathlist_gminus);
  
  sense_paths_gplus = (Path_T *) List_to_array_n(&sense_npaths_gplus,sense_pathlist_gplus);
  sense_paths_gminus = (Path_T *) List_to_array_n(&sense_npaths_gminus,sense_pathlist_gminus);
  qsort(sense_paths_gplus,sense_npaths_gplus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  qsort(sense_paths_gminus,sense_npaths_gminus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  Hitlistpool_free_list(&sense_pathlist_gplus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&sense_pathlist_gminus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  
  antisense_pathlist_gplus = List_append(unsolved_antisense_pathlist_gplus,antisense_pathlist_gplus);
  antisense_pathlist_gminus = List_append(unsolved_antisense_pathlist_gminus,antisense_pathlist_gminus);

  antisense_paths_gplus = (Path_T *) List_to_array_n(&antisense_npaths_gplus,antisense_pathlist_gplus);
  antisense_paths_gminus = (Path_T *) List_to_array_n(&antisense_npaths_gminus,antisense_pathlist_gminus);
  qsort(antisense_paths_gplus,antisense_npaths_gplus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  qsort(antisense_paths_gminus,antisense_npaths_gminus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  Hitlistpool_free_list(&antisense_pathlist_gplus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&antisense_pathlist_gminus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));


  /* Build plus */
  if ((npaths = sense_npaths_gplus + antisense_npaths_gplus) > 0) {
    MALLOC_ALIGN(this->all_univdiagonals_gplus,npaths*sizeof(Univcoord_T));
    this->all_auxinfo_gplus = (Auxinfo_T *) MALLOC(npaths*sizeof(Auxinfo_T));
    this->all_pathstores_gplus = (Pathstore_T *) MALLOC(npaths*sizeof(Pathstore_T));

    i = j = 0;
    k = 0;

    while (i < sense_npaths_gplus && j < antisense_npaths_gplus) {
      if ((univdiagonal_i = sense_paths_gplus[i]->main_univdiagonal) < (univdiagonal_j = antisense_paths_gplus[j]->main_univdiagonal)) {
	this->all_univdiagonals_gplus[k] = univdiagonal_i;
	path = sense_paths_gplus[i];
	this->all_auxinfo_gplus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gplus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);


	ii = i + 1;
	while (ii < sense_npaths_gplus && sense_paths_gplus[ii]->main_univdiagonal == univdiagonal_i) {
	  ii++;
	}
	for (l = i; l < ii; l++) {
	  path = sense_paths_gplus[l];
	  pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
							 hitlistpool_trace(__FILE__,__LINE__));
	}
	pathstore->solvedp = true;
	pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);

	i = ii; k++;

      } else if (univdiagonal_j < univdiagonal_i) {
	this->all_univdiagonals_gplus[k] = univdiagonal_j;
	path = antisense_paths_gplus[j];
	this->all_auxinfo_gplus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gplus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

	jj = j + 1;
	while (jj < antisense_npaths_gplus && antisense_paths_gplus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gplus[l];
	  pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							     hitlistpool_trace(__FILE__,__LINE__));
	}
	pathstore->solvedp = true;
	pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

	j = jj; k++;
      
      } else {
	this->all_univdiagonals_gplus[k] = univdiagonal_i;
	path = sense_paths_gplus[i];
	this->all_auxinfo_gplus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gplus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

	ii = i + 1;
	while (ii < sense_npaths_gplus && sense_paths_gplus[ii]->main_univdiagonal == univdiagonal_i) {
	  ii++;
	}
	for (l = i; l < ii; l++) {
	  path = sense_paths_gplus[l];
	  pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
							 hitlistpool_trace(__FILE__,__LINE__));
	}

	jj = j + 1;
	while (jj < antisense_npaths_gplus && antisense_paths_gplus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gplus[l];
	  pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							     hitlistpool_trace(__FILE__,__LINE__));
	}

	pathstore->solvedp = true;
	pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);
	pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

	i = ii; j = jj; k++;
      }
    }

    while (i < sense_npaths_gplus) {
      univdiagonal_i = sense_paths_gplus[i]->main_univdiagonal;
      this->all_univdiagonals_gplus[k] = univdiagonal_i;
      path = sense_paths_gplus[i];
      this->all_auxinfo_gplus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
      pathstore = this->all_pathstores_gplus[k] =
	Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

      ii = i + 1;
      while (ii < sense_npaths_gplus && sense_paths_gplus[ii]->main_univdiagonal == univdiagonal_i) {
	ii++;
      }
      for (l = i; l < ii; l++) {
	path = sense_paths_gplus[l];
	pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
						       hitlistpool_trace(__FILE__,__LINE__));
      }
      pathstore->solvedp = true;
      pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);

      i = ii; k++;
    }

    while (j < antisense_npaths_gplus) {
      univdiagonal_j = antisense_paths_gplus[j]->main_univdiagonal;
      this->all_univdiagonals_gplus[k] = univdiagonal_j;
      path = antisense_paths_gplus[j];
      this->all_auxinfo_gplus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
      pathstore = this->all_pathstores_gplus[k] =
	Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

      jj = j + 1;
      while (jj < antisense_npaths_gplus && antisense_paths_gplus[jj]->main_univdiagonal == univdiagonal_j) {
	jj++;
      }
      for (l = j; l < jj; l++) {
	path = antisense_paths_gplus[l];
	pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							   hitlistpool_trace(__FILE__,__LINE__));
      }
      pathstore->solvedp = true;
      pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

      j = jj; k++;
    }

    this->all_nunivdiagonals_gplus = k;
  }

  /* Build minus */
  if ((npaths = sense_npaths_gminus + antisense_npaths_gminus) > 0) {
    MALLOC_ALIGN(this->all_univdiagonals_gminus,npaths*sizeof(Univcoord_T));
    this->all_auxinfo_gminus = (Auxinfo_T *) MALLOC(npaths*sizeof(Auxinfo_T));
    this->all_pathstores_gminus = (Pathstore_T *) MALLOC(npaths*sizeof(Pathstore_T));

    i = j = 0;
    k = 0;

    while (i < sense_npaths_gminus && j < antisense_npaths_gminus) {
      if ((univdiagonal_i = sense_paths_gminus[i]->main_univdiagonal) < (univdiagonal_j = antisense_paths_gminus[j]->main_univdiagonal)) {
	this->all_univdiagonals_gminus[k] = univdiagonal_i;
	path = sense_paths_gminus[i];
	this->all_auxinfo_gminus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gminus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

	ii = i + 1;
	while (ii < sense_npaths_gminus && sense_paths_gminus[ii]->main_univdiagonal == univdiagonal_i) {
	  ii++;
	}
	for (l = i; l < ii; l++) {
	  path = sense_paths_gminus[l];
	  pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
							 hitlistpool_trace(__FILE__,__LINE__));
	}
	pathstore->solvedp = true;
	pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);

	i = ii; k++;

      } else if (univdiagonal_j < univdiagonal_i) {
	this->all_univdiagonals_gminus[k] = univdiagonal_j;
	path = antisense_paths_gminus[j];
	this->all_auxinfo_gminus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gminus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

	jj = j + 1;
	while (jj < antisense_npaths_gminus && antisense_paths_gminus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gminus[l];
	  pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							     hitlistpool_trace(__FILE__,__LINE__));
	}
	pathstore->solvedp = true;
	pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

	j = jj; k++;
      
      } else {
	this->all_univdiagonals_gminus[k] = univdiagonal_i;
	path = sense_paths_gminus[i];
	this->all_auxinfo_gminus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
	pathstore = this->all_pathstores_gminus[k] =
	  Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

	ii = i + 1;
	while (ii < sense_npaths_gminus && sense_paths_gminus[ii]->main_univdiagonal == univdiagonal_i) {
	  ii++;
	}
	for (l = i; l < ii; l++) {
	  path = sense_paths_gminus[l];
	  pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
							 hitlistpool_trace(__FILE__,__LINE__));
	}

	jj = j + 1;
	while (jj < antisense_npaths_gminus && antisense_paths_gminus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gminus[l];
	  pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							     hitlistpool_trace(__FILE__,__LINE__));
	}

	pathstore->solvedp = true;
	pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);
	pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

	i = ii; j = jj; k++;
      }
    }

    while (i < sense_npaths_gminus) {
      univdiagonal_i = sense_paths_gminus[i]->main_univdiagonal;
      this->all_univdiagonals_gminus[k] = univdiagonal_i;
      path = sense_paths_gminus[i];
      this->all_auxinfo_gminus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
      pathstore = this->all_pathstores_gminus[k] =
	Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

      ii = i + 1;
      while (ii < sense_npaths_gminus && sense_paths_gminus[ii]->main_univdiagonal == univdiagonal_i) {
	ii++;
      }
      for (l = i; l < ii; l++) {
	path = sense_paths_gminus[l];
	pathstore->complete_sense_paths = Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
						       hitlistpool_trace(__FILE__,__LINE__));
      }
      pathstore->solvedp = true;
      pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);

      i = ii; k++;
    }

    while (j < antisense_npaths_gminus) {
      univdiagonal_j = antisense_paths_gminus[j]->main_univdiagonal;
      this->all_univdiagonals_gminus[k] = univdiagonal_j;
      path = antisense_paths_gminus[j];
      this->all_auxinfo_gminus[k] = (Auxinfo_T) NULL;  /* Auxinfo_new_tr(auxinfopool); */
      pathstore = this->all_pathstores_gminus[k] =
	Pathstore_new_tr(pathstorepool,path->chrnum,path->chroffset,path->chrhigh);

      jj = j + 1;
      while (jj < antisense_npaths_gminus && antisense_paths_gminus[jj]->main_univdiagonal == univdiagonal_j) {
	jj++;
      }
      for (l = j; l < jj; l++) {
	path = antisense_paths_gminus[l];
	pathstore->complete_antisense_paths = Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
							   hitlistpool_trace(__FILE__,__LINE__));
      }
      pathstore->solvedp = true;
      pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);

      j = jj; k++;
    }

    this->all_nunivdiagonals_gminus = k;
  }

  FREE(sense_paths_gplus);
  FREE(sense_paths_gminus);
  FREE(antisense_paths_gplus);
  FREE(antisense_paths_gminus);

  return;
}



static List_T
paired_search_tr (int *found_score_paired, int *found_score_5, int *found_score_3,
		  int sufficient_score_5, int sufficient_score_3,

		  Shortread_T queryseq5, Shortread_T queryseq3,
		  int querylength5, int querylength3, int genestrand,
		  
		  int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
		  T this5, T this3, Knownsplicing_T knownsplicing,
		  
		  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		  
		  int nmismatches_allowed_5, int nmismatches_allowed_3,
		  int nmismatches_filter_5, int nmismatches_filter_3,
		  int mincoverage_filter_5, int mincoverage_filter_3,

		  Trdiagpool_T trdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		  Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {
  
  List_T pathpairs = NULL;
  Method_T last_method_5 = METHOD_INIT, last_method_3 = METHOD_INIT;

  debug(printf(">>Entered paired_search_tr with methods %s and %s\n",
	       Method_string(last_method_5),Method_string(last_method_3)));

#if 0
  /* Take the larger of maxpaths_search and 10*maxpaths_report */
  maxpairedpaths = maxpaths_search;
  if (maxpairedpaths < 10*maxpaths_report) {
    maxpairedpaths = 10*maxpaths_report;
  }
#endif

  /* TODO: Return level from single_read.  Then if we don't have concordant pathpairs, take the min(level5,level3) and start paired_readfrom t
here */


  /* A.  Initial search.  Try each side until we get a result */
  /* Disjunction: Stops when we have found pathpairs or both ends have gone to TR_PREVALENT (conjunction would mean either end) */
  while (last_method_5 < TR_PREVALENT || last_method_3 < TR_PREVALENT) {
    if ((pathpairs = paired_search_trdiagonals(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
					       sufficient_score_5,sufficient_score_3,
					       &last_method_5,&last_method_3,
					       
					       pathpairs,this5,this3,knownsplicing,
					       queryseq5,queryseq3,querylength5,querylength3,
					       query5_compress_fwd,query5_compress_rev,
					       query3_compress_fwd,query3_compress_rev,
					       mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					       
					       nmismatches_filter_5,nmismatches_filter_3,
					       mincoverage_filter_5,mincoverage_filter_3,
					       intlistpool,uintlistpool,univcoordlistpool,
					       listpool,trpathpool,pathpool,transcriptpool,vectorpool,
					       hitlistpool,/*method_goal*/TR_PREVALENT)) != NULL) {
      debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		   Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		   *found_score_paired,*found_score_5,*found_score_3));
      /* Path_print(((Pathpair_T) pathpairs->first)->path5); */
      /* Path_print(((Pathpair_T) pathpairs->first)->path3); */

      return pathpairs;
      
    } else {
      debug1(printf("Continuing with found scores %d and %d\n",*found_score_5,*found_score_3));
    }
  }
  
  
  /* B.  Advance by side separately to TR_EXT */
  /* Helps to avoid genomic search */
  /* Disjunction: Stops when we have found pathpairs or both ends have gone to TR_EXT (conjunction would mean either end) */
  while (last_method_5 < TR_EXT || last_method_3 < TR_EXT) {
    if (last_method_3 >= TR_EXT) {
      if ((pathpairs = paired_read_next_method_tr_5(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    sufficient_score_5,&last_method_5,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_5,nmismatches_allowed_5,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,
						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	return pathpairs;
      }
      
    } else if (last_method_5 >= TR_EXT) {
      if ((pathpairs = paired_read_next_method_tr_3(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    sufficient_score_3,&last_method_3,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_3,nmismatches_allowed_3,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,
						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	return pathpairs;
      }
      
    } else if ((*found_score_5) >= (*found_score_3)) {
      if ((pathpairs = paired_read_next_method_tr_5(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    sufficient_score_5,&last_method_5,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_5,nmismatches_allowed_5,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,
						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	return pathpairs;
      }
      
    } else {
      if ((pathpairs = paired_read_next_method_tr_3(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    sufficient_score_3,&last_method_3,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_3,nmismatches_allowed_3,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,
						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	return pathpairs;
      }
    }
  }
  

  /* D.  Prep for genomic search.  Previously converted all trpaths to
     paths.  Now create univdiagonals and auxinfo in Stage1_T
     objects */
  debug(printf("Preparing for genomic search\n"));

  create_univdiagonals_auxinfo_from_tr(&(*found_score_5),this5,knownsplicing,
				       queryseq5,querylength5,
				       this5->sense_trpaths,this5->antisense_trpaths,
				       query5_compress_fwd,query5_compress_rev,
				       auxinfopool,pathstorepool,intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
				       /*first_read_p*/true);

  create_univdiagonals_auxinfo_from_tr(&(*found_score_3),this3,knownsplicing,
				       queryseq3,querylength3,
				       this3->sense_trpaths,this3->antisense_trpaths,
				       query3_compress_fwd,query3_compress_rev,
				       auxinfopool,pathstorepool,intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
				       /*first_read_p*/false);

  debug(printf("Returning from paired_search_tr with %d pathpairs\n",List_length(pathpairs)));
  return pathpairs;
}


#if 0
static List_T
find_inner_fusions (int *found_score_5, int *found_score_3,

		    List_T local_sense_paths5, List_T local_antisense_paths5,
		    List_T local_sense_paths3, List_T local_antisense_paths3,

		    Shortread_T queryseq5, Shortread_T queryseq3,
		    char *queryuc_ptr_5, char *queryrc5, char *queryuc_ptr_3, char *queryrc3,
		    Knownsplicing_T knownsplicing, Univdiag_T *novel_univdiags_alloc,
		    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc, T this5, T this3,

		    Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		    Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		    int nmismatches_allowed_5, int nmismatches_allowed_3,
		    
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {


  List_T pathpairs = NULL, p;
  Path_T path5, path3;


  for (p = local_sense_paths5; p != NULL; p = List_next(p)) {
    path5 = (Path_T) List_head(p);
    if ((path5->plusp == true && Path_unextended_qend_p(path5,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true) ||
	(path5->plusp == false && Path_unextended_qstart_p(path5,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true)) {
      pathpairs = find_inner_fusions_path5(&(*found_score_5),pathpairs,path5,local_sense_paths3,
					   this5,query5_compress_fwd,query5_compress_rev,queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,knownsplicing,nmismatches_allowed_5,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,vectorpool,hitlistpool,transcriptpool);

    }
  }

  for (p = local_antisense_paths5; p != NULL; p = List_next(p)) {
    path5 = (Path_T) List_head(p);
    if ((path5->plusp == true && Path_unextended_qend_p(path5,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true) ||
	(path5->plusp == false && Path_unextended_qstart_p(path5,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true)) {
      pathpairs = find_inner_fusions_path5(&(*found_score_5),pathpairs,path5,local_antisense_paths3,
					   this5,query5_compress_fwd,query5_compress_rev,queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,knownsplicing,nmismatches_allowed_5,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
    }
  }

  for (p = local_sense_paths3; p != NULL; p = List_next(p)) {
    path3 = (Path_T) List_head(p);
    if ((path3->plusp == true && Path_unextended_qstart_p(path3,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true) ||
	(path3->plusp == false && Path_unextended_qend_p(path3,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true)) {
      pathpairs = find_inner_fusions_path3(&(*found_score_3),pathpairs,path3,local_sense_paths5,
					   this3,query3_compress_fwd,query3_compress_rev,queryseq5,queryseq3,
					   queryuc_ptr_3,queryrc3,novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,knownsplicing,nmismatches_allowed_3,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,vectorpool,
					   hitlistpool,transcriptpool);
    }
  }

  for (p = local_antisense_paths3; p != NULL; p = List_next(p)) {
    path3 = (Path_T) List_head(p);
    if ((path3->plusp == true && Path_unextended_qstart_p(path3,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true) ||
	(path3->plusp == false && Path_unextended_qend_p(path3,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == true)) {
      pathpairs = find_inner_fusions_path3(&(*found_score_3),pathpairs,path3,local_antisense_paths5,
					   this3,query3_compress_fwd,query3_compress_rev,queryseq5,queryseq3,
					   queryuc_ptr_3,queryrc3,novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,knownsplicing,nmismatches_allowed_3,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,vectorpool,
					   hitlistpool,transcriptpool);
    }
  }

  return pathpairs;
}
#endif



/* pairtype can be CONCORDANT, CONCORDANT_TRANSLOCATIONS, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG, UNPAIRED */
static Pairtype_T
determine_pairtype (Path_T path5, Path_T path3) {
  Univcoord_T genomicstart5, genomicend5, genomicstart3, genomicend3;
  Chrpos_T pairmax;

  debug14(printf("Entered determine_pairtype\n"));
  if (path5->chrnum != path3->chrnum) {
    debug14(printf("Returning unpaired because path5 chrnum %d != path3 chrnum %d\n",
		   path5->chrnum,path3->chrnum));
    return UNPAIRED;

  } else if (path5->fusion_querystart_junction != NULL ||
      path5->fusion_queryend_junction != NULL ||
      path3->fusion_querystart_junction != NULL ||
      path3->fusion_queryend_junction != NULL) {
    debug14(printf("Returning translocations\n"));
    /* On the same chromosome, so consider them concordant, but we may need to check their distance */
    return CONCORDANT_TRANSLOCATIONS;

#ifdef TO_FIX
  } else if (Transcript_concordant_p(path5->transcripts_consistent,path3->transcripts_consistent) == true) {
    debug14(printf("Returning concordant based on transcriptome\n"));
    return CONCORDANT;
#endif

  } else if (path5->plusp != path3->plusp) {
    debug14(printf("Returning paired_inversion\n"));
    return PAIRED_INVERSION;

  } else if (path5->plusp == true) {
#if 0
    genomicstart5 = Path_genomicstart(path5);
    genomicend5 = Path_genomicend(path5);
    genomicstart3 = Path_genomicstart(path3);
    genomicend3 = Path_genomicend(path3);
#else
    genomicstart5 = path5->main_univdiagonal;
    genomicend5 = path5->main_univdiagonal;
    genomicstart3 = path3->main_univdiagonal;
    genomicend3 = path3->main_univdiagonal;
#endif

#ifdef TO_FIX
    if (path5->circularalias == 0 && path3->circularalias == 0) {
      /* Keep coordinates as is */
    } else if (path5->circularalias == 0) {
      if (path3->circularalias == -1) {
	debug14(printf("Have to alias path3\n"));
	assert(circularp[path3->chrnum] == true);
	chrlength = (path3->chrhigh - path3->chroffset)/2;
	genomicstart3 += chrlength;
	genomicend3 += chrlength;
      }
    } else if (path3->circularalias == 0) {
      if (path5->circularalias == +1) {
	debug14(printf("Have to unalias path5\n"));
	assert(circularp[path5->chrnum] == true);
	chrlength = (path5->chrhigh - path5->chroffset)/2;
	genomicstart5 -= chrlength;
	genomicend5 -= chrlength;
      }
    }
#endif

    if (genomicend3 < genomicstart5) {
      debug14(printf("(plus) path3 genomicend %u < path5 genomicstart %u => Returning paired_scramble\n",
		     genomicend3,genomicstart5));
      return PAIRED_SCRAMBLE;
    } else {
      if (circularp[path5->chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }

      if (genomicstart3 > genomicend5 + pairmax) {
	debug14(printf("Returning paired_toolong\n"));
	return PAIRED_TOOLONG;
      } else {
	debug14(printf("Returning concordant\n"));
	return CONCORDANT;
      }
    }

  } else {
#if 0
    genomicstart5 = Path_genomicstart(path5);
    genomicend5 = Path_genomicend(path5);
    genomicstart3 = Path_genomicstart(path3);
    genomicend3 = Path_genomicend(path3);
#else
    genomicstart5 = path5->main_univdiagonal;
    genomicend5 = path5->main_univdiagonal;
    genomicstart3 = path3->main_univdiagonal;
    genomicend3 = path3->main_univdiagonal;
#endif

#ifdef TO_FIX
    if (path5->circularalias == 0 && path3->circularalias == 0) {
      /* Keep coordinates as is */
    } else if (path5->circularalias == 0) {
      if (path3->circularalias == +1) {
	debug14(printf("Have to unalias path3\n"));
	assert(circularp[path3->chrnum] == true);
	chrlength = (path3->chrhigh - path3->chroffset)/2;
	genomicstart3 -= chrlength;
	genomicend3 -= chrlength;
      }
    } else if (path3->circularalias == 0) {
      if (path5->circularalias == -1) {
	debug14(printf("Have to alias path5\n"));
	assert(circularp[path5->chrnum] == true);
	chrlength = (path5->chrhigh - path5->chroffset)/2;
	genomicstart5 += chrlength;
	genomicend5 += chrlength;
      }
    }
#endif

    if (genomicend3 > genomicstart5) {
      debug14(printf("(minus) path3 genomicend %u > path5 genomicstart %u => Returning paired_scramble\n",
		     genomicend3,genomicstart5));
      return PAIRED_SCRAMBLE;
    } else {
      if (circularp[path3->chrnum] == true) {
	pairmax = pairmax_circular;
      } else {
	pairmax = pairmax_linear;
      }

      if (genomicstart3 + pairmax < genomicend5) {
	debug14(printf("Returning paired_toolong\n"));
	return PAIRED_TOOLONG;
      } else {
	debug14(printf("Returning concordant\n"));
	return CONCORDANT;
      }
    }
  }
}


/* Have three lists: pathpairs, samechr, and conc_transloc => result */
/* final_pairtype can be CONCORDANT_TRANSLOCATIONS, CONCORDANT, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG, UNPAIRED */


static Pathpair_T *
consolidate_results (int *found_score_5, int *found_score_3, Pairtype_T *final_pairtype,
		     int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		     Path_T **patharray5, int *npaths5_primary, int *npaths5_altloc, int *first_absmq5, int *second_absmq5,
		     Path_T **patharray3, int *npaths3_primary, int *npaths3_altloc, int *first_absmq3, int *second_absmq3,

		     List_T sense_paths5, List_T antisense_paths5,
		     List_T sense_paths3, List_T antisense_paths3,

		     Shortread_T queryseq5, Shortread_T queryseq3,
		     char *queryuc_ptr_5, char *queryrc5, int querylength5,
		     char *queryuc_ptr_3, char *queryrc3, int querylength3,

		     Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		     Univdiag_T *novel_univdiags_alloc,Univcoord_T *novel_diagonals_alloc,
		     unsigned short *localdb_alloc, T this5, T this3,

		     int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

		     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

		     int nmismatches_allowed_5, int nmismatches_allowed_3,
		     int nmismatches_filter_5, int nmismatches_filter_3,
		     int mincoverage_filter_5, int mincoverage_filter_3,
		
		     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Univdiagpool_T univdiagpool,
		     Pathpool_T pathpool, Vectorpool_T vectorpool,
		     Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		     Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {

  Pathpair_T *pathpairarray, pathpair;
  List_T unresolved_pathpairs = NULL;
  List_T pathpairs, complete_paths, unextended_paths5 = NULL, unextended_paths3 = NULL,
    singlepaths5, singlepaths3, paths5, paths3, p;
  Path_T path, path5, path3, copy;
  Compress_T query_compress;
  char *queryptr;
  int sensedir5, sensedir3;
  /* int tr_insertlength; */
  int i, j;


  debug16(printf("Entered consolidate_results with %d + %d paths5 and %d + %d paths3\n",
		 List_length(sense_paths5),List_length(antisense_paths5),
		 List_length(sense_paths3),List_length(antisense_paths3)));

  /* Unpaired ends */
  if (sense_paths5 != NULL && antisense_paths5 == NULL) {
    sensedir5 = SENSE_FORWARD;
  } else if (sense_paths5 == NULL && antisense_paths5 != NULL) {
    sensedir5 = SENSE_ANTI;
  } else {
    sensedir5 = SENSE_NULL;
  }
	
  if (sense_paths3 != NULL && antisense_paths3 == NULL) {
    sensedir3 = SENSE_FORWARD;
  } else if (sense_paths3 == NULL && antisense_paths3 != NULL) {
    sensedir3 = SENSE_ANTI;
  } else {
    sensedir3 = SENSE_NULL;
  }

  if (sensedir5 == SENSE_NULL && sensedir3 == SENSE_NULL) {
    /* No restriction */
    paths5 = List_append(sense_paths5,antisense_paths5);
    paths3 = List_append(sense_paths3,antisense_paths3);
  } else if (sensedir5 == SENSE_FORWARD && sensedir3 == SENSE_ANTI) {
    /* Contradiction */
    paths5 = List_append(sense_paths5,antisense_paths5);
    paths3 = List_append(sense_paths3,antisense_paths3);
  } else if (sensedir5 == SENSE_ANTI && sensedir3 == SENSE_FORWARD) {
    /* Contradiction */
    paths5 = List_append(sense_paths5,antisense_paths5);
    paths3 = List_append(sense_paths3,antisense_paths3);
  } else if (sensedir5 == SENSE_FORWARD || sensedir3 == SENSE_FORWARD) {
    /* Restrict to sense */
    paths5 = sense_paths5;
    paths3 = sense_paths3;
  } else if (sensedir5 == SENSE_ANTI || sensedir3 == SENSE_ANTI) {
    /* Restrict to antisense */
    paths5 = antisense_paths5;
    paths3 = antisense_paths3;
  } else {
    fprintf(stderr,"Unexpected combination of sensedirs\n");
    abort();
  }

  debug16(printf("HAVE %d PATHS5\n",List_length(paths5)));
  debug16(printf("HAVE %d PATHS3\n",List_length(paths3)));

  singlepaths5 = (List_T) NULL;
  for (p = paths5; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->plusp == true) {
      queryptr = queryuc_ptr_5;
      query_compress = query5_compress_fwd;
    } else {
      queryptr = queryrc5;
      query_compress = query5_compress_rev;
    }
    /* Need to call Path_extend primarily for the trimming, which can weed out false alignments past the end of the chromosome,
       but might as well extend also if possible */

    if (path->transcriptome_method_p == true) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      singlepaths5 = Hitlist_push(singlepaths5,hitlistpool,(void *) copy
				  hitlistpool_trace(__FILE__,__LINE__));
    } else if (path->extendedp == true) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      singlepaths5 = Hitlist_push(singlepaths5,hitlistpool,(void *) copy
				  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((complete_paths = Path_extend(&(*found_score_5),&unextended_paths5,
					     /*original_path*/path,queryseq5,queryptr,querylength5,
					     mismatch_positions_alloc_5,novel_univdiags_alloc,
					     novel_diagonals_alloc,localdb_alloc,
					     this5,knownsplicing,knownindels,
					     query_compress,query5_compress_fwd,query5_compress_rev,
					     /*genestrand*/0,nmismatches_allowed_5,
					     /*paired_end_p*/false,/*lowp*/true,
					     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     univdiagpool,vectorpool,hitlistpool,spliceendsgen5,
					     /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      /* Found complete paths */
      singlepaths5 = List_append(complete_paths,singlepaths5);

    } else {
      /* unextended_paths5 accumulate */
    }
  }
  Hitlistpool_free_list(&paths5,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));


  singlepaths3 = (List_T) NULL;
  for (p = paths3; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->plusp == true) {
      queryptr = queryuc_ptr_3;
      query_compress = query3_compress_fwd;
    } else {
      queryptr = queryrc3;
      query_compress = query3_compress_rev;
    }
    /* Need to call Path_extend primarily for the trimming, which can weed out false alignments past the end of the chromosome,
       but might as well extend also if possible */
    if (path->transcriptome_method_p == true) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      singlepaths3 = Hitlist_push(singlepaths3,hitlistpool,(void *) copy
				  hitlistpool_trace(__FILE__,__LINE__));

    } else if (path->extendedp == true) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      singlepaths3 = Hitlist_push(singlepaths3,hitlistpool,(void *) copy
				  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((complete_paths = Path_extend(&(*found_score_3),&unextended_paths3,
					     /*original_path*/path,queryseq3,queryptr,querylength3,
					     mismatch_positions_alloc_3,novel_univdiags_alloc,
					     novel_diagonals_alloc,localdb_alloc,
					     this3,knownsplicing,knownindels,
					     query_compress,query3_compress_fwd,query3_compress_rev,
					     /*genestrand*/0,nmismatches_allowed_3,
					     /*paired_end_p*/false,/*lowp*/true,
					     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     univdiagpool,vectorpool,hitlistpool,spliceendsgen3,
					     /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      /* Found complete paths */
      singlepaths3 = List_append(complete_paths,singlepaths3);

    } else {
      /* unextended paths3 accumulate */
    }
  }
  Hitlistpool_free_list(&paths3,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));

  debug16(printf("HAVE %d SINGLEPATHS5, COMPLETE\n",List_length(singlepaths5)));
  debug16(printf("HAVE %d SINGLEPATHS5, UNEXTENDED\n",List_length(unextended_paths5)));

  debug16(printf("HAVE %d SINGLEPATHS3, COMPLETE\n",List_length(singlepaths3)));
  debug16(printf("HAVE %d SINGLEPATHS3, UNEXTENDED\n",List_length(unextended_paths3)));

  if (singlepaths5 != NULL) {
    Path_gc(&unextended_paths5,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
  } else {
    singlepaths5 = unextended_paths5;
  }

  if (singlepaths3 != NULL) {
    Path_gc(&unextended_paths3,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
  } else {
    singlepaths3 = unextended_paths3;
  }

  if (singlepaths5 == (List_T) NULL) {
    *npaths5_primary = *npaths5_altloc = 0;
    *patharray5 = (Path_T *) NULL;
  } else {
    *patharray5 = (Path_T *) List_to_array_out(singlepaths5,NULL); /* Return value */
    *patharray5 = Path_eval_and_sort(&(*npaths5_primary),&(*npaths5_altloc),&(*first_absmq5),&(*second_absmq5),
				     *patharray5,List_length(singlepaths5),
				     query5_compress_fwd,query5_compress_rev,queryuc_ptr_5,queryrc5,
				     Shortread_quality_string(queryseq5),
				     nmismatches_filter_5,mincoverage_filter_5,
				     intlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,hitlistpool,/*filterp*/true);

    Hitlistpool_free_list(&singlepaths5,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }

  if (singlepaths3 == (List_T) NULL) {
    *npaths3_primary = *npaths3_altloc = 0;
    *patharray3 = (Path_T *) NULL;
  } else {
    *patharray3 = (Path_T *) List_to_array_out(singlepaths3,NULL); /* Return value */
    *patharray3 = Path_eval_and_sort(&(*npaths3_primary),&(*npaths3_altloc),&(*first_absmq3),&(*second_absmq3),
				     *patharray3,List_length(singlepaths3),
				     query3_compress_fwd,query3_compress_rev,queryuc_ptr_3,queryrc3,
				     Shortread_quality_string(queryseq3),
				     nmismatches_filter_3,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,hitlistpool,/*filterp*/true);

    Hitlistpool_free_list(&singlepaths3,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }

#ifdef DEBUG16
  printf("For single paths5, got %d+%d paths\n",*npaths5_primary,*npaths5_altloc);
  for (int i = 0; i < *npaths5_primary; i++) {
    Path_print((*patharray5)[i]);
  }

  printf("For single paths3, got %d+%d paths\n",*npaths3_primary,*npaths3_altloc);
  for (int i = 0; i < *npaths3_primary; i++) {
    Path_print((*patharray3)[i]);
  }
#endif

  /* Check for concordance */
  pathpairs = NULL;
  for (i = 0; i < *npaths5_primary + *npaths5_altloc; i++) {
    path5 = (*patharray5)[i];
    for (j = 0; j < *npaths3_primary + *npaths3_altloc; j++) {
      path3 = (*patharray3)[j];
      if (determine_pairtype(path5,path3) != CONCORDANT) {
	/* Skip */
      } else if (path5->plusp == true) {
	if ((pathpair = Pathpair_new_concordant(&unresolved_pathpairs,/*pathL*/path5,/*pathH*/path3,
						/*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,

						nmismatches_filter_5,nmismatches_filter_3,
						mincoverage_filter_5,mincoverage_filter_3,

						intlistpool,univcoordlistpool,listpool,
						pathpool,vectorpool,transcriptpool,hitlistpool,
						/*check_inner_p*/false,/*copyLp*/true,/*copyHp*/true)) != NULL) {
	  debug16(printf("Found a concordant pair, gplus\n"));
	  pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				   hitlistpool_trace(__FILE__,__LINE__));
	  debug16(Path_print(pathpair->path5));
	  debug16(Path_print(pathpair->path3));
	}
      } else {
	if ((pathpair = Pathpair_new_concordant(&unresolved_pathpairs,/*pathL*/path3,/*pathH*/path5,
						/*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,

						nmismatches_filter_5,nmismatches_filter_3,
						mincoverage_filter_5,mincoverage_filter_3,

						intlistpool,univcoordlistpool,listpool,
						pathpool,vectorpool,transcriptpool,hitlistpool,
						/*check_inner_p*/false,/*copyLp*/true,/*copyHp*/true)) != NULL) {
	  debug16(printf("Found a concordant pair, gminus\n"));
	  pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				   hitlistpool_trace(__FILE__,__LINE__));
	  debug16(Path_print(pathpair->path5));
	  debug16(Path_print(pathpair->path3));
	}
      }
    }
  }

  /* Handle pathpairs computed from single ends */
  if (pathpairs != NULL) {
    debug16(printf("HAVE %d pathpairs\n",List_length(pathpairs)));

    /* Don't want to filter here, since Pathpair_eval_and_sort needs to keep both sensedirs and perform resolve */
    /* pathpairs = Pathpair_filter(pathpairs,
       intlistpool,univcoordlistpool,listpool,pathpool,hitlistpool); */

    *npaths_primary = List_length(pathpairs);
    *npaths_altloc = 0;	/* TODO: Determine whether any paths are on the altloc chromosome */

    pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
    pathpairarray = Pathpair_eval_and_sort(&(*found_score_5),&(*found_score_3),
					   &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					   pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					   /*quality_string_5*/Shortread_quality_string(queryseq5),
					   /*quality_string_3*/Shortread_quality_string(queryseq3),
					   mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					   knownsplicing,knownindels,

					   nmismatches_allowed_5,nmismatches_allowed_3,
					   nmismatches_filter_5,nmismatches_filter_3,
					   mincoverage_filter_5,mincoverage_filter_3,
					   querylength5,querylength3,
					   univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,/*filterp*/true);

    Hitlistpool_free_list(&pathpairs,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));  /* Allocated by hitlistpool */

    debug16(printf("After Pathpair_eval_and_sort, have %d + %d pathpairs\n",
		   *npaths_primary,*npaths_altloc));
  }

  if ((*npaths_primary) + (*npaths_altloc) > 0) {
    *final_pairtype = CONCORDANT;
    debug(printf("Success in finding pathpairs, so clearing single paths\n"));
    for (i = 0; i < *npaths5_primary + *npaths5_altloc; i++) {
      path5 = (*patharray5)[i];
      Path_free(&path5,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    *npaths5_primary = *npaths5_altloc = 0;
    FREE_OUT(*patharray5);
    
    for (i = 0; i < *npaths3_primary + *npaths3_altloc; i++) {
      path3 = (*patharray3)[i];
      Path_free(&path3,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    *npaths3_primary = *npaths3_altloc = 0;
    FREE_OUT(*patharray3);

    *patharray5 = *patharray3 = (Path_T *) NULL;
    debug(printf("(1) Returning %d %d %d\n",*npaths_primary,*npaths5_primary,*npaths3_primary));

  } else if (*npaths5_primary != 1 || *npaths3_primary != 1) {
    *final_pairtype = UNPAIRED;
    *npaths_primary = *npaths_altloc = 0;
    pathpairarray = (Pathpair_T *) NULL;
    debug(printf("(2) Returning %d %d %d\n",*npaths_primary,*npaths5_primary,*npaths3_primary));

  } else {
    *final_pairtype = determine_pairtype((*patharray5)[0],(*patharray3)[0]);
    *npaths_primary = *npaths_altloc = 0;
    pathpairarray = (Pathpair_T *) NULL;
    debug(printf("(3) Returning %d %d %d\n",*npaths_primary,*npaths5_primary,*npaths3_primary));
  }

  return pathpairarray;
}


#ifdef DEBUG
static void
print_pathpairarray_contents (Pathpair_T *pathpairarray, int n) {
  int i;

  for (i = 0; i < n; i++) {
    printf("%p %p pathpairarray\n",pathpairarray[i]->path5,pathpairarray[i]->path3);
  }

  return;
}
#endif



#if 0
static void
list_paths (List_T list) {
  List_T p;
  Path_T path;

  for (p = list; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    Path_print(path);
  }

  return;
}
#endif


static List_T
paired_search_exhaustive (int *found_score_paired, int *found_score_5, int *found_score_3,
			  List_T *unresolved_pathpairs,

			  Shortread_T queryseq5, char *queryuc_ptr_5, char *queryrc5,
			  Shortread_T queryseq3, char *queryuc_ptr_3, char *queryrc3,
			  int querylength5, int querylength3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  
			  int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			  T this5, T this3, Knownsplicing_T knownsplicing, Knownindels_T knownindels,

			  int nmismatches_allowed_5, int nmismatches_allowed_3,

			  int nmismatches_filter_5, int nmismatches_filter_3,
			  int mincoverage_filter_5, int mincoverage_filter_3,

			  Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			  Univdiagpool_T univdiagpool, 
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			  Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			  Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {

  List_T pathpairs = NULL;

  int *indices_gplus, *indices_gminus;
  int nindices_gplus, nindices_gminus, i, k;
  int index1, index2;
  int max_count, max_count_gplus, max_count_gminus, count;
  
  int ignore;
  int qstart5, qend5, qstart3, qend3;
  Univcoord_T univdiagonal5, univdiagonal3;
  Auxinfo_T auxinfo5, auxinfo3;
  Pathstore_T pathstore5, pathstore3;
  
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;


  /* plus */
#ifdef LARGE_GENOMES
  indices_gplus =
    Intersect_approx_indices_uint8(&nindices_gplus,
				   this5->exhaustive_gplus,this5->nexhaustive_gplus,/*diagterm1*/-querylength5,
				   this3->exhaustive_gplus,this3->nexhaustive_gplus,/*diagterm2*/0,
				   /*below_slop*/0,/*above_slop*/concordance_distance);
#else
  indices_gplus =
    Intersect_approx_indices_uint4(&nindices_gplus,
				   this5->exhaustive_gplus,this5->nexhaustive_gplus,/*diagterm1*/-querylength5,
				   this3->exhaustive_gplus,this3->nexhaustive_gplus,/*diagterm2*/0,
				   /*below_slop*/0,/*above_slop*/concordance_distance);
#endif

  max_count_gplus = 0;
  for (i = 0, k = 0; i < nindices_gplus; i++, k += 2) {
    index1 = indices_gplus[k];
    index2 = indices_gplus[k+1];

    if ((count = this5->exhaustive_counts_gplus[index1] + this3->exhaustive_counts_gplus[index2]) > max_count_gplus) {
      max_count_gplus = count;
    }
  }
  debug(printf("plus nindices %d, max count %d\n",nindices_gplus,max_count_gplus));


  /* minus */
#ifdef LARGE_GENOMES
  indices_gminus =
    Intersect_approx_indices_uint8(&nindices_gminus,
				   this3->exhaustive_gminus,this3->nexhaustive_gminus,/*diagterm1*/-querylength3,
				   this5->exhaustive_gminus,this5->nexhaustive_gminus,/*diagterm2*/0,
				   /*below_slop*/0,/*above_slop*/concordance_distance);
#else
  indices_gminus =
    Intersect_approx_indices_uint4(&nindices_gminus,
				   this3->exhaustive_gminus,this3->nexhaustive_gminus,/*diagterm1*/-querylength3,
				   this5->exhaustive_gminus,this5->nexhaustive_gminus,/*diagterm2*/0,
				   /*below_slop*/0,/*above_slop*/concordance_distance);
#endif

  max_count_gminus = 0;
  for (i = 0, k = 0; i < nindices_gminus; i++, k += 2) {
    index1 = indices_gminus[k];
    index2 = indices_gminus[k+1];

    if ((count = this3->exhaustive_counts_gminus[index1] + this5->exhaustive_counts_gminus[index2]) > max_count_gminus) {
      max_count_gminus = count;
    }
  }
  debug(printf("minus nindices %d, max count %d\n",nindices_gminus,max_count_gminus));


  /* Find threshold */
  if (max_count_gplus > max_count_gminus) {
    max_count = max_count_gplus;
  } else {
    max_count = max_count_gminus;
  }
  if (max_count < 4) {
    max_count = 1;
  } else {
    max_count -= 3;		/* to allow for suboptimal solutions */
  }


  /* plus */
  if (max_count_gplus >= max_count) {
    for (i = 0, k = 0; i < nindices_gplus; i++, k += 2) {
      index1 = indices_gplus[k];
      index2 = indices_gplus[k+1];

      if (this5->exhaustive_counts_gplus[index1] + this3->exhaustive_counts_gplus[index2] >= max_count) {
	univdiagonal5 = this5->exhaustive_gplus[index1];
	qstart5 = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					     univdiagonal5,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);
	qend5 = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					    univdiagonal5,querylength5,/*pos5*/0,/*pos3*/querylength5,
					    /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);

	if (univdiagonal5 < (Univcoord_T) querylength5) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,0,univdiagonal5);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,univdiagonal5 - querylength5,univdiagonal5);
	}
	auxinfo5 = Auxinfo_new(EXHAUSTIVE,qstart5,qend5,auxinfopool);
	pathstore5 = Pathstore_new(pathstorepool);
	Path_solve_from_univdiagonal(&(*found_score_5),
				     
				     &pathstore5->unextended_sense_paths,&pathstore5->unextended_antisense_paths,
				     &pathstore5->complete_sense_paths,&pathstore5->complete_antisense_paths,
				     
				     univdiagonal5,auxinfo5,queryseq5,/*queryptr*/queryuc_ptr_5,
				     /*query_compress*/query5_compress_fwd,
				     query5_compress_fwd,query5_compress_rev,
				     /*plusp*/true,querylength5,mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     this5,knownsplicing,knownindels,
				     nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				     chrnum,chroffset,chrhigh,
				     intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen5,/*method*/EXHAUSTIVE,/*find_splices_p*/true);
	pathstore5->solvedp = true;
	Pathstore_set_best_paths(pathstore5,hitlistpool); /* Needed for make_pathpairs to work */
      
	univdiagonal3 = this3->exhaustive_gplus[index2];
	qstart3 = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					     univdiagonal3,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);
	qend3 = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					    univdiagonal3,querylength3,/*pos5*/0,/*pos3*/querylength3,
					    /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);

	if (univdiagonal3 < (Univcoord_T) querylength3) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,0,univdiagonal3);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,univdiagonal3 - querylength3,univdiagonal3);
	}
	auxinfo3 = Auxinfo_new(EXHAUSTIVE,qstart3,qend3,auxinfopool);
	pathstore3 = Pathstore_new(pathstorepool);
	Path_solve_from_univdiagonal(&(*found_score_3),
				   
				     &pathstore3->unextended_sense_paths,&pathstore3->unextended_antisense_paths,
				     &pathstore3->complete_sense_paths,&pathstore3->complete_antisense_paths,
				   
				     univdiagonal3,auxinfo3,queryseq3,/*queryptr*/queryuc_ptr_3,
				     /*query_compress*/query3_compress_fwd,
				     query3_compress_fwd,query3_compress_rev,
				     /*plusp*/true,querylength3,mismatch_positions_alloc_3,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     this3,knownsplicing,knownindels,
				     nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				     chrnum,chroffset,chrhigh,
				     intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen3,/*method*/EXHAUSTIVE,/*find_splices_p*/true);
	pathstore3->solvedp = true;
	Pathstore_set_best_paths(pathstore3,hitlistpool); /* Needed for make_pathpairs to work */
      
	debug(printf("plus max count %d, %u %d..%d to %u %d..%d\n",
		     max_count,univdiagonal5,qstart5,qend5,univdiagonal3,qstart3,qend3));
	debug(printf("(2) Calling make_pathpairs\n"));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstore5,pathstore3,queryseq5,queryseq3,/*plusp*/true,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/false,/*only_completeH_p*/false);

	Auxinfo_free(&auxinfo5,univdiagpool,auxinfopool);
	Auxinfo_free(&auxinfo3,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore5,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
	Pathstore_free(&pathstore3,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }


  /* minus */
  if (max_count_gminus >= max_count) {
    for (i = 0, k = 0; i < nindices_gminus; i++, k += 2) {
      index1 = indices_gminus[k];
      index2 = indices_gminus[k+1];

      if (this3->exhaustive_counts_gminus[index1] + this5->exhaustive_counts_gminus[index2] >= max_count) {
	univdiagonal3 = this3->exhaustive_gminus[index1];
	qstart3 = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					     univdiagonal3,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);
	qend3 = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					    univdiagonal3,querylength3,/*pos5*/0,/*pos3*/querylength3,
					    /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);

	auxinfo3 = Auxinfo_new(EXHAUSTIVE,qstart3,qend3,auxinfopool);
	pathstore3 = Pathstore_new(pathstorepool);
	if (univdiagonal3 < (Univcoord_T) querylength3) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,0,univdiagonal3);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,univdiagonal3 - querylength3,univdiagonal3);
	}
	Path_solve_from_univdiagonal(&(*found_score_3),
				     
				     &pathstore3->unextended_sense_paths,&pathstore3->unextended_antisense_paths,
				     &pathstore3->complete_sense_paths,&pathstore3->complete_antisense_paths,
				     
				     univdiagonal3,auxinfo3,queryseq3,/*queryptr*/queryrc3,
				     /*query_compress*/query3_compress_rev,
				     query3_compress_fwd,query3_compress_rev,
				     /*plusp*/false,querylength3,mismatch_positions_alloc_3,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     this3,knownsplicing,knownindels,
				     nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				     chrnum,chroffset,chrhigh,
				     intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen3,/*method*/EXHAUSTIVE,/*find_splices_p*/true);
	pathstore3->solvedp = true;
	Pathstore_set_best_paths(pathstore3,hitlistpool); /* Needed for make_pathpairs to work */

	univdiagonal5 = this5->exhaustive_gminus[index2];
	qstart5 = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					     univdiagonal5,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);
	qend5 = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					    univdiagonal5,querylength5,/*pos5*/0,/*pos3*/querylength5,
					    /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/index1part);

	if (univdiagonal5 < (Univcoord_T) querylength5) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,0,univdiagonal5);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,univdiagonal5 - querylength5,univdiagonal5);
	}
	auxinfo5 = Auxinfo_new(EXHAUSTIVE,qstart5,qend5,auxinfopool);
	pathstore5 = Pathstore_new(pathstorepool);
	Path_solve_from_univdiagonal(&(*found_score_5),
				     
				     &pathstore5->unextended_sense_paths,&pathstore5->unextended_antisense_paths,
				     &pathstore5->complete_sense_paths,&pathstore5->complete_antisense_paths,

				     univdiagonal5,auxinfo5,queryseq5,/*queryptr*/queryrc5,
				     /*query_compress*/query5_compress_rev,
				     query5_compress_fwd,query5_compress_rev,
				     /*plusp*/false,querylength5,mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     this5,knownsplicing,knownindels,
				     nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				     chrnum,chroffset,chrhigh,
				     intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen5,/*method*/EXHAUSTIVE,/*find_splices_p*/true);
	pathstore5->solvedp = true;
	Pathstore_set_best_paths(pathstore5,hitlistpool); /* Needed for make_pathpairs to work */

	debug(printf("minus max count %d, %u %d..%d to %u %d..%d\n",
		     max_count,univdiagonal3,qstart3,qend3,univdiagonal5,qstart5,qend5));
	debug(printf("(3) Calling make_pathpairs\n"));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstore3,pathstore5,queryseq3,queryseq5,/*plusp*/false,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/false,/*only_completeH_p*/false);

	Auxinfo_free(&auxinfo5,univdiagpool,auxinfopool);
	Auxinfo_free(&auxinfo3,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore5,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
	Pathstore_free(&pathstore3,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }

  FREE(indices_gminus);
  FREE(indices_gplus);

  return pathpairs;
}


static int
find_max_count_upstream (int *bestk, Univcoord_T *exhaustive, int *counts, int nexhaustive,
			 int k, Univcoord_T high_univdiagonal) {
  int max_counts = 0;

  *bestk = -1;
  while (k < nexhaustive && exhaustive[k] < high_univdiagonal) {
    if (counts[k] > max_counts) {
      *bestk = k;
      max_counts = counts[k];
    }
    k++;
  }

  return max_counts;
}

static int
find_max_count_downstream (int *bestk, Univcoord_T *exhaustive, int *counts,
			   int k, Univcoord_T low_univdiagonal) {
  int max_counts = 0;

  *bestk = -1;
  while (k >= 0 && exhaustive[k] > low_univdiagonal) {
    if (counts[k] > max_counts) {
      *bestk = k;
      max_counts = counts[k];
    }
    k--;
  }

  return max_counts;
}


#if 0
static int
find_prevalent_upstream (Univcoord_T *exhaustive, int *counts, int nexhaustive,
			 int k, Univcoord_T high_univdiagonal, int required_count) {
  int bestk = -1;
  int max_counts = 0;

  while (k < nexhaustive && exhaustive[k] < high_univdiagonal) {
    if (counts[k] > max_counts) {
      bestk = k;
      max_counts = counts[k];
    }
    k++;
  }

  if (max_counts < required_count) {
    return -1;
  } else {
    return bestk;
  }
}
#endif


#if 0
static int
find_prevalent_downstream (Univcoord_T *exhaustive, int *counts,
			   int k, Univcoord_T low_univdiagonal, int required_count) {
  int bestk = -1;
  int max_counts = 0;

  while (k >= 0 && exhaustive[k] > low_univdiagonal) {
    if (counts[k] > max_counts) {
      bestk = k;
      max_counts = counts[k];
    }
    k--;
  }

  if (max_counts < required_count) {
    return -1;
  } else {
    return bestk;
  }
}
#endif


static List_T
paths5_mates (int *found_score_paired, int *found_score_3,
	      List_T *unresolved_pathpairs, List_T pathpairs,
		     
	      Shortread_T queryseq5, Shortread_T queryseq3,
	      char *queryuc_ptr_3, char *queryrc3, int querylength3,
	      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	      int *mismatch_positions_alloc_3, Univdiag_T *novel_univdiags_alloc,
	      Univcoord_T *novel_diagonals_alloc,
	      unsigned short *localdb_alloc, T this5, T this3,
		  
	      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
	      int nmismatches_allowed_3,
		  
	      int nmismatches_filter_5, int nmismatches_filter_3,
	      int mincoverage_filter_5, int mincoverage_filter_3,

	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
	      Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
	      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
	      Spliceendsgen_T spliceendsgen3, bool only_complete5_p) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo_mate;
  Pathstore_T pathstore, pathstore_mate;
  int qstart, qend;

  int i;
  int k_gplus3, k_gminus3, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  int ignore;


  /* From solved paths, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered paths5_mates\n"));

  /* First, find threshold for mates */
  k_gplus3 = 0;
  k_gminus3 = this3->nexhaustive_gminus - 1;

  /* Go upstream */
  for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
    univdiagonal = this5->all_univdiagonals_gplus[i];
    while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
      k_gplus3++;
    }

    /* auxinfo = this5->all_auxinfo_gplus[i]; */
    pathstore = this5->all_pathstores_gplus[i];

    if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      /* (1,2) Look upstream for 3' end matches */
      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
      if ((pathstore->mate_count =
	   find_max_count_upstream(&pathstore->mate_bestk,this3->exhaustive_gplus,
				   this3->exhaustive_counts_gplus,this3->nexhaustive_gplus,
				   k_gplus3,high_univdiagonal)) > max_count_gplus) {
	max_count_gplus = pathstore->mate_count;
      }
    }
  }


  /* Go downstream */
  for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
    univdiagonal = this5->all_univdiagonals_gminus[i];
    while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
      k_gminus3--;
    }

    /* auxinfo = this5->all_auxinfo_gminus[i]; */
    pathstore = this5->all_pathstores_gminus[i];

    if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      /* (3,4) Look downstream for 3' end matches */
      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
      if ((pathstore->mate_count =
	   find_max_count_downstream(&pathstore->mate_bestk,
				     this3->exhaustive_gminus,this3->exhaustive_counts_gminus,
				     k_gminus3,low_univdiagonal)) > max_count_gminus) {
	max_count_gminus = pathstore->mate_count;
      }
    }
  }
  

  /* Find threshold */
  debug11(printf("Got max_count_gplus %d, max_count_gminus %d\n",
		 max_count_gplus,max_count_gminus));
  if (max_count_gplus > max_count_gminus) {
    max_count = max_count_gplus;
  } else {
    max_count = max_count_gminus;
  }
  if (max_count < 4) {
    max_count = 1;
  } else {
    max_count -= 3;		/* to allow for suboptimal solutions */
  }


  /* Second, find mates */
  k_gplus3 = 0;
  k_gminus3 = this3->nexhaustive_gminus - 1;

  if (max_count_gplus >= max_count) {
    /* Go upstream */
    for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
      univdiagonal = this5->all_univdiagonals_gplus[i];
      while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
	k_gplus3++;
      }

      /* auxinfo = this5->all_auxinfo_gplus[i]; */
      pathstore = this5->all_pathstores_gplus[i];

      if (only_complete5_p == true && pathstore->complete_sense_p == false) {
	/* Skip */
      } else if (pathstore->best_sense_paths != NULL) {
	/* (1) Look upstream for 3' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this3->exhaustive_gplus[k];
	  qstart = this3->exhaustive_qstart_gplus[k];
	  qend = this3->exhaustive_qend_gplus[k];
	  debug11(printf("(1) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gplus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this3->exhaustive_gplus,this3->exhaustive_qstart_gplus,
					   this3->exhaustive_qend_gplus,this3->nexhaustive_gplus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,this3,
							    /*queryptr*/queryuc_ptr_3,querylength3,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query3_compress_fwd,
							    /*plusp*/true,/*genestrand*/0,
							    genomebits,nmismatches_allowed_3)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(1) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}

	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_3),
				     
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				     
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryuc_ptr_3,
				       /*query_compress*/query3_compress_fwd,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/true,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	
	  debug(printf("(4) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete5_p,/*only_completeH_p*/false);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      if (only_complete5_p == true && pathstore->complete_antisense_p == false) {
	/* Skip */
      } else if (pathstore->best_antisense_paths != NULL) {
	/* (2) Look upstream for 3' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this3->exhaustive_gplus[k];
	  qstart = this3->exhaustive_qstart_gplus[k];
	  qend = this3->exhaustive_qend_gplus[k];
	  debug11(printf("(2) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gplus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this3->exhaustive_gplus,this3->exhaustive_qstart_gplus,
					   this3->exhaustive_qend_gplus,this3->nexhaustive_gplus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,this3,
							    /*queryptr*/queryuc_ptr_3,querylength3,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query3_compress_fwd,
							    /*plusp*/true,/*genestrand*/0,
							    genomebits,nmismatches_allowed_3)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(2) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}

	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_3),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryuc_ptr_3,
				       /*query_compress*/query3_compress_fwd,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/true,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	
	  debug(printf("(5) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete5_p,/*only_completeH_p*/false);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }


  if (max_count_gminus >= max_count) {
    /* Go downstream */
    for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      univdiagonal = this5->all_univdiagonals_gminus[i];
      while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
	k_gminus3--;
      }

      /* auxinfo = this5->all_auxinfo_gminus[i]; */
      pathstore = this5->all_pathstores_gminus[i];

      if (only_complete5_p == true && pathstore->complete_sense_p == false) {
	/* Skip */
      } else if (pathstore->best_sense_paths != NULL) {
	/* (3) Look downstream for 3' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this3->exhaustive_gminus[k];
	  qstart = this3->exhaustive_qstart_gminus[k];
	  qend = this3->exhaustive_qend_gminus[k];
	  debug11(printf("(3) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gminus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this3->exhaustive_gminus,this3->exhaustive_qstart_gminus,
					   this3->exhaustive_qend_gminus,this3->nexhaustive_gminus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,this3,
							     /*queryptr*/queryrc3,querylength3,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query3_compress_rev,
							     /*plusp*/false,/*genestrand*/0,
							     genomebits,nmismatches_allowed_3)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(3) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}
	
	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_3),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryrc3,
				       /*query_compress*/query3_compress_rev,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/false,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(6) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete5_p);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      if (only_complete5_p == true && pathstore->complete_antisense_p == false) {
	/* Skip */
      } else if (pathstore->best_antisense_paths != NULL) {
	/* (4) Look downstream for 3' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this3->exhaustive_gminus[k];
	  qstart = this3->exhaustive_qstart_gminus[k];
	  qend = this3->exhaustive_qend_gminus[k];
	  debug11(printf("(4) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gminus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this3->exhaustive_gminus,this3->exhaustive_qstart_gminus,
					   this3->exhaustive_qend_gminus,this3->nexhaustive_gminus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,this3,
							     /*queryptr*/queryrc3,querylength3,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query3_compress_rev,
							     /*plusp*/false,/*genestrand*/0,
							     genomebits,nmismatches_allowed_3)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(4) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}
	
	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_3),
				     
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				     
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryrc3,
				       /*query_compress*/query3_compress_rev,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/false,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	
	  debug(printf("(7) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete5_p);

	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  return pathpairs;
}



static List_T
paths3_mates (int *found_score_paired, int *found_score_5,
	      List_T *unresolved_pathpairs, List_T pathpairs,
		     
	      Shortread_T queryseq5, Shortread_T queryseq3,
	      char *queryuc_ptr_5, char *queryrc5, int querylength5,
	      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	      int *mismatch_positions_alloc_5, Univdiag_T *novel_univdiags_alloc,
	      Univcoord_T *novel_diagonals_alloc,
	      unsigned short *localdb_alloc, T this5, T this3,
		  
	      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
	      int nmismatches_allowed_5,
		  
	      int nmismatches_filter_5, int nmismatches_filter_3,
	      int mincoverage_filter_5, int mincoverage_filter_3,

	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
	      Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
	      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
	      Spliceendsgen_T spliceendsgen5, bool only_complete3_p) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo_mate;
  Pathstore_T pathstore, pathstore_mate;
  int qstart, qend;

  int i;
  int k_gplus5, k_gminus5, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  int ignore;


  /* From solved paths, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered paths3_mates\n"));

  /* First, find threshold for mates */
  k_gplus5 = this5->nexhaustive_gplus - 1;
  k_gminus5 = 0;

  /* Go downstream */
  for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
    univdiagonal = this3->all_univdiagonals_gplus[i];
    while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
      k_gplus5--;
    }

    /* auxinfo = this3->all_auxinfo_gplus[i]; */
    pathstore = this3->all_pathstores_gplus[i];

    if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      /* (5,6) Look downstream for 5' end matches */
      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
      if ((pathstore->mate_count =
	   find_max_count_downstream(&pathstore->mate_bestk,
				     this5->exhaustive_gplus,this5->exhaustive_counts_gplus,
				     k_gplus5,low_univdiagonal)) > max_count_gplus) {
	max_count_gplus = pathstore->mate_count;
      }
    }
  }

  /* Go upstream */
  for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
    univdiagonal = this3->all_univdiagonals_gminus[i];
    while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
      k_gminus5++;
    }

    /* auxinfo = this3->all_auxinfo_gminus[i]; */
    pathstore = this3->all_pathstores_gminus[i];

    if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      /* (7,8) Look upstream for 5' end matches */
      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
      if ((pathstore->mate_count =
	   find_max_count_upstream(&pathstore->mate_bestk,
				   this5->exhaustive_gminus,this5->exhaustive_counts_gminus,this5->nexhaustive_gminus,
				   k_gminus5,high_univdiagonal)) > max_count_gminus) {
	max_count_gminus = pathstore->mate_count;
      }
    }
  }

  /* Find threshold */
  debug11(printf("Got max_count_gplus %d, max_count_gminus %d\n",
		 max_count_gplus,max_count_gminus));
  if (max_count_gplus > max_count_gminus) {
    max_count = max_count_gplus;
  } else {
    max_count = max_count_gminus;
  }
  if (max_count < 4) {
    max_count = 1;
  } else {
    max_count -= 3;		/* Allows searching for suboptimal solutions */
  }


  /* Second, find mates */
  k_gplus5 = this5->nexhaustive_gplus - 1;
  k_gminus5 = 0;

  if (max_count_gplus >= max_count) {
    /* Go downstream */
    for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
      univdiagonal = this3->all_univdiagonals_gplus[i];
      while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
	k_gplus5--;
      }

      debug11(printf("Finding downstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gplus[i]; */
      pathstore = this3->all_pathstores_gplus[i];

      if (only_complete3_p == true && pathstore->complete_sense_p == false) {
	/* Skip */
      } else if (pathstore->best_sense_paths != NULL) {
	/* (5) Look downstream for 5' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this5->exhaustive_gplus[k];
	  qstart = this5->exhaustive_qstart_gplus[k];
	  qend = this5->exhaustive_qend_gplus[k];
	  debug11(printf("(5) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gplus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this5->exhaustive_gplus,this5->exhaustive_qstart_gplus,
					   this5->exhaustive_qend_gplus,this5->nexhaustive_gplus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,this5,
							     /*queryptr*/queryuc_ptr_5,querylength5,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query5_compress_fwd,
							     /*plusp*/true,/*genestrand*/0,
							     genomebits,nmismatches_allowed_5)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(5) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}
	
	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryuc_ptr_5,
				       /*query_compress*/query5_compress_fwd,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/true,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(8) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete3_p);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      if (only_complete3_p == true && pathstore->complete_antisense_p == false) {
	/* Skip */
      } else if (pathstore->best_antisense_paths != NULL) {
	/* (6) Look downstream for 5' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
	
	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this5->exhaustive_gplus[k];
	  qstart = this5->exhaustive_qstart_gplus[k];
	  qend = this5->exhaustive_qend_gplus[k];
	  debug11(printf("(6) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gplus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this5->exhaustive_gplus,this5->exhaustive_qstart_gplus,
					   this5->exhaustive_qend_gplus,this5->nexhaustive_gplus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,this5,
							     /*queryptr*/queryuc_ptr_5,querylength5,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query5_compress_fwd,
							     /*plusp*/true,/*genestrand*/0,
							     genomebits,nmismatches_allowed_5)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(6) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}

	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryuc_ptr_5,
				       /*query_compress*/query5_compress_fwd,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/true,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(9) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete3_p);

	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }


  if (max_count_gminus >= max_count) {
    /* Go upstream */
    for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
      univdiagonal = this3->all_univdiagonals_gminus[i];
      while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
	k_gminus5++;
      }
      
      debug11(printf("Finding upstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gminus[i]; */
      pathstore = this3->all_pathstores_gminus[i];

      if (only_complete3_p == true && pathstore->complete_sense_p == false) {
	/* Skip */
      } else if (pathstore->best_sense_paths != NULL) {
	/* (7) Look upstream for 5' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this5->exhaustive_gminus[k];
	  qstart = this5->exhaustive_qstart_gminus[k];
	  qend = this5->exhaustive_qend_gminus[k];
	  debug11(printf("(7) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gminus[k]));

	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this5->exhaustive_gminus,this5->exhaustive_qstart_gminus,
					   this5->exhaustive_qend_gminus,this5->nexhaustive_gminus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);

	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,this5,
							    /*queryptr*/queryrc5,querylength5,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query5_compress_rev,
							    /*plusp*/false,/*genestrand*/0,
							    genomebits,nmismatches_allowed_5)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(7) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}
	
	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryrc5,
				       /*query_compress*/query5_compress_rev,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/false,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(10) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete3_p,/*only_completeH_p*/false);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      if (only_complete3_p == true && pathstore->complete_antisense_p == false) {
	/* Skip */
      } else if (pathstore->best_antisense_paths != NULL) {
	/* (8) Look upstream for 5' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
	
	auxinfo_mate = (Auxinfo_T) NULL;
	if (pathstore->mate_count >= max_count) {
	  k = pathstore->mate_bestk;
	  assert(k >= 0);

	  univdiagonal_mate = this5->exhaustive_gminus[k];
	  qstart = this5->exhaustive_qstart_gminus[k];
	  qend = this5->exhaustive_qend_gminus[k];
	  debug11(printf("(8) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gminus[k]));
	  
	  auxinfo_mate =
	    Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					   this5->exhaustive_gminus,this5->exhaustive_qstart_gminus,
					   this5->exhaustive_qend_gminus,this5->nexhaustive_gminus,
					   univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	  pathstore_mate = Pathstore_new(pathstorepool);
	  
	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,this5,
							    /*queryptr*/queryrc5,querylength5,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query5_compress_rev,
							    /*plusp*/false,/*genestrand*/0,
							    genomebits,nmismatches_allowed_5)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  debug11(printf("(8) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	  auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	  pathstore_mate = Pathstore_new(pathstorepool);
	}
	
	if (auxinfo_mate != NULL) {
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryrc5,
				       /*query_compress*/query5_compress_rev,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/false,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(11) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete3_p,/*only_completeH_p*/false);
	  
	  Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  return pathpairs;
}


#if 0
static List_T
univdiagonals5_mates (int *found_score_paired, int *found_score_5, int *found_score_3,
		      List_T *unresolved_pathpairs, List_T pathpairs,
		     
		      Shortread_T queryseq5, Shortread_T queryseq3,
		      char *queryuc_ptr_5, char *queryrc5, int querylength5,
		      char *queryuc_ptr_3, char *queryrc3, int querylength3,
		      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
		      int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, T this5, T this3,
		  
		      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		      int nmismatches_allowed_5, int nmismatches_allowed_3,
		  
		      int nmismatches_filter_5, int nmismatches_filter_3,
		      int mincoverage_filter_5, int mincoverage_filter_3,

		      Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
		      Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo, auxinfo_mate;
  int qstart, qend;

  int i;
  int k_gplus3, k_gminus3, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  bool complete_sense_p, complete_antisense_p; /* ignored */
  int ignore;


  /* From univdiagonals, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered univdiagonals5_mates\n"));

  /* First, find threshold for mates */
  k_gplus3 = 0;
  k_gminus3 = this3->nexhaustive_gminus - 1;

  /* Go upstream */
  for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
    univdiagonal = this5->all_univdiagonals_gplus[i];
    while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
      k_gplus3++;
    }

    /* (1,2) Look upstream for 3' end matches */
    /* auxinfo = this5->all_auxinfo_gplus[i]; */
    pathstore = this5->all_pathstores_gplus[i];

    low_univdiagonal = univdiagonal;
    high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

    if ((pathstore->mate_count =
	 find_max_count_upstream(&pathstore->mate_bestk,
				 this3->exhaustive_gplus,this3->exhaustive_counts_gplus,this3->nexhaustive_gplus,
				 k_gplus3,high_univdiagonal)) > max_count_gplus) {
      max_count_gplus = pathstore->mate_count;
    }
  }

  /* Go downstream */
  for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
    univdiagonal = this5->all_univdiagonals_gminus[i];
    while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
      k_gminus3--;
    }

    /* (3,4) Look downstream for 3' end matches */
    /* auxinfo = this5->all_auxinfo_gminus[i]; */
    pathstore = this5->all_pathstores_gminus[i];

    high_univdiagonal = univdiagonal;
    low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

    if ((pathstore->mate_count =
	 find_max_count_downstream(&pathstore->mate_bestk,
				   this3->exhaustive_gminus,this3->exhaustive_counts_gminus,
				   k_gminus3,low_univdiagonal)) > max_count_gminus) {
      max_count_gminus = pathstore->mate_count;
    }
  }
  

  /* Find threshold */
  debug11(printf("Got max_count_gplus %d, max_count_gminus %d\n",max_count_gplus,max_count_gminus));
  if (max_count_gplus > max_count_gminus) {
    max_count = max_count_gplus;
  } else {
    max_count = max_count_gminus;
  }
  if (max_count < 4) {
    max_count = 1;
  } else {
    max_count -= 3;		/* to allow for suboptimal solutions */
  }


  /* Second, find mates */
  k_gplus3 = 0;
  k_gminus3 = this3->nexhaustive_gminus - 1;

  if (max_count_gplus >= max_count) {
    /* Go upstream */
    for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
      univdiagonal = this5->all_univdiagonals_gplus[i];
      while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
	k_gplus3++;
      }

      /* (1,2) Look upstream for 3' end matches */
      debug11(printf("Finding upstream mates for 5' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this5->all_auxinfo_gplus[i]; */
      pathstore = this5->all_pathstores_gplus[i];

      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

      auxinfo_mate = (Auxinfo_T) NULL;
      if (pathstore->mate_count >= max_count) {
	k = pathstore->mate_bestk;
	assert(k >= 0);

	univdiagonal_mate = this3->exhaustive_gplus[k];
	qstart = this3->exhaustive_qstart_gplus[k];
	qend = this3->exhaustive_qend_gplus[k];
	debug11(printf("(1,2) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gplus[k]));

	auxinfo_mate =
	  Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					 this3->exhaustive_gplus,this3->exhaustive_qstart_gplus,
					 this3->exhaustive_qend_gplus,this3->nexhaustive_gplus,
					 univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	pathstore_mate = Pathstore_new(pathstorepool);

      } else if (localdb != NULL &&
		 querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		 (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							  novel_diagonals_alloc,localdb_alloc,this3,
							  /*queryptr*/queryuc_ptr_3,querylength3,
							  /*low_univdiagonal*/univdiagonal,high_univdiagonal,
							  /*query_compress*/query3_compress_fwd,
							  /*plusp*/true,/*genestrand*/0,
							  genomebits,nmismatches_allowed_3)) != 0) {
	qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					    univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					    /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					   univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					   /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	debug11(printf("(1,2) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - univdiagonal,qstart,qend));
	auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	pathstore_mate = Pathstore_new(pathstorepool);
      }

      if (auxinfo_mate != NULL) {
	if (pathstore->solvedp == false) {
	  /* Solve for anchor if necessary */
	  solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				     &(*found_score_5),univdiagonal,&(this5->all_auxinfo_gplus[i]),pathstore,
				     
				     queryseq5,/*queryptr*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,
				     this5,knownsplicing,knownindels,
				     
				     mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     /*query_compress*/query5_compress_fwd,
				     query5_compress_fwd,query5_compress_rev,
				     
				     nmismatches_allowed_5,/*genestrand*/0,
				     auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen5,/*plusp*/true,/*first_read_p*/true,/*lowp*/true,
				     /*set_best_paths_p*/true);

	}
	
	if (pathstore->complete_sense_paths != NULL || pathstore->complete_antisense_paths != NULL ||
	    pathstore->unextended_sense_paths != NULL || pathstore->unextended_antisense_paths != NULL) {
	  /* Solve for mate */
	  Path_solve_from_univdiagonal(&(*found_score_3),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryuc_ptr_3,
				       /*query_compress*/query3_compress_fwd,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/true,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(12) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/false);
	}

	Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }


  if (max_count_gminus >= max_count) {
    /* Go downstream */
    for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      univdiagonal = this5->all_univdiagonals_gminus[i];
      while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
	k_gminus3--;
      }

      /* (3,4) Look downstream for 3' end matches */
      debug11(printf("Finding downstream mates for 5' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this5->all_auxinfo_gminus[i]; */
      pathstore = this5->all_pathstores_gminus[i];

      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

      auxinfo_mate = (Auxinfo_T) NULL;
      if (pathstore->mate_count >= max_count) {
	k = pathstore->mate_bestk;
	assert(k >= 0);

	univdiagonal_mate = this3->exhaustive_gminus[k];
	qstart = this3->exhaustive_qstart_gminus[k];
	qend = this3->exhaustive_qend_gminus[k];
	debug11(printf("(3,4) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this3->exhaustive_counts_gminus[k]));

	auxinfo_mate =
	  Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					 this3->exhaustive_gminus,this3->exhaustive_qstart_gminus,
					 this3->exhaustive_qend_gminus,this3->nexhaustive_gminus,
					 univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	pathstore_mate = Pathstore_new(pathstorepool);

      } else if (localdb != NULL &&
		 querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		 (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							   novel_diagonals_alloc,localdb_alloc,this3,
							   /*queryptr*/queryrc3,querylength3,
							   low_univdiagonal,/*high_univdiagonal*/univdiagonal,
							   /*query_compress*/query3_compress_rev,
							   /*plusp*/false,/*genestrand*/0,
							   genomebits,nmismatches_allowed_3)) != 0) {
	qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					    univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					    /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					   univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					   /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	debug11(printf("(3,4) highest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal - univdiagonal_mate,qstart,qend));
	auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	pathstore_mate = Pathstore_new(pathstorepool);
      }
	
      if (auxinfo_mate != NULL) {
	if (pathstore->solvedp == false) {
	  /* Solve for anchor if necessary */
	  solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				     &(*found_score_5),univdiagonal,&(this5->all_auxinfo_gminus[i]),pathstore,

				     queryseq5,/*queryptr*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,
				     this5,knownsplicing,knownindels,
				 
				     mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     /*query_compress*/query5_compress_rev,
				     query5_compress_fwd,query5_compress_rev,
				 
				     nmismatches_allowed_5,/*genestrand*/0,
				     auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen5,/*plusp*/false,/*first_read_p*/true,/*lowp*/false,
				     /*set_best_paths_p*/true);
	}

	if (auxinfo->complete_sense_paths != NULL || auxinfo->complete_antisense_paths != NULL ||
	    auxinfo->unextended_sense_paths != NULL || auxinfo->unextended_antisense_paths != NULL) {
	  /* Solve for mate */
	  Path_solve_from_univdiagonal(&(*found_score_3),
				       
				       &auxinfo_mate->unextended_sense_paths,&auxinfo_mate->unextended_antisense_paths,
				       &auxinfo_mate->complete_sense_paths,&auxinfo_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryrc3,
				       /*query_compress*/query3_compress_rev,
				       query3_compress_fwd,query3_compress_rev,
				       /*plusp*/false,querylength3,mismatch_positions_alloc_3,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this3,knownsplicing,knownindels,
				       nmismatches_allowed_3,/*paired_end_p*/true,/*first_read_p*/false,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  auxinfo_mate->solvedp = true;
	  Auxinfo_set_best_paths(auxinfo_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(13) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*auxinfoL*/auxinfo_mate,/*auxinfoH*/auxinfo,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/false);
	}

	Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }

  return pathpairs;
}
#endif


#if 0
static List_T
univdiagonals3_mates (int *found_score_paired, int *found_score_5, int *found_score_3,
		      List_T *unresolved_pathpairs, List_T pathpairs,
		     
		      Shortread_T queryseq5, Shortread_T queryseq3,
		      char *queryuc_ptr_5, char *queryrc5, int querylength5,
		      char *queryuc_ptr_3, char *queryrc3, int querylength3,
		      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
		      int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, T this5, T this3,
		  
		      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		      int nmismatches_allowed_5, int nmismatches_allowed_3,
		  
		      int nmismatches_filter_5, int nmismatches_filter_3,
		      int mincoverage_filter_5, int mincoverage_filter_3,

		      Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
		      Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo, auxinfo_mate;
  int qstart, qend;

  int i;
  int k_gplus5, k_gminus5, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  bool complete_sense_p, complete_antisense_p; /* ignored */
  int ignore;


  /* From univdiagonals, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered univdiagonals3_mates\n"));

  /* First, find threshold for mates */
  k_gplus5 = this5->nexhaustive_gplus - 1;
  k_gminus5 = 0;

  /* Go downstream */
  for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
    univdiagonal = this3->all_univdiagonals_gplus[i];
    while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
      k_gplus5--;
    }

    /* (5,6) Look downstream for 5' end matches */
    /* auxinfo = this3->all_auxinfo_gplus[i]; */
    pathstore = this3->all_pathstores_gplus[i];

    high_univdiagonal = univdiagonal;
    low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

    if ((pathstore->mate_count =
	 find_max_count_downstream(&pathstore->mate_bestk,
				   this5->exhaustive_gplus,this5->exhaustive_counts_gplus,
				   k_gplus5,low_univdiagonal)) > max_count_gplus) {
      max_count_gplus = pathstore->mate_count;
    }
  }

  /* Go upstream */
  for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
    univdiagonal = this3->all_univdiagonals_gminus[i];
    while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
      k_gminus5++;
    }

    /* (7,8) Look upstream for 5' end matches */
    /* auxinfo = this3->all_auxinfo_gminus[i]; */
    pathstore = this3->all_pathstores_gminus[i];

    low_univdiagonal = univdiagonal;
    high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

    if ((pathstore->mate_count =
	 find_max_count_upstream(&pathstore->mate_bestk,
				 this5->exhaustive_gminus,this5->exhaustive_counts_gminus,this5->nexhaustive_gminus,
				 k_gminus5,high_univdiagonal)) > max_count_gminus) {
      max_count_gminus = pathstore->mate_count;
    }
  }


  /* Find threshold */
  debug11(printf("Got max_count_gplus %d, max_count_gminus %d\n",
		 max_count_gplus,max_count_gminus));
  if (max_count_gplus > max_count_gminus) {
    max_count = max_count_gplus;
  } else {
    max_count = max_count_gminus;
  }
  if (max_count < 4) {
    max_count = 1;
  } else {
    max_count -= 3;		/* Allows searching for suboptimal solutions */
  }


  /* Second, find mates */
  k_gplus5 = this5->nexhaustive_gplus - 1;
  k_gminus5 = 0;

  if (max_count_gplus >= max_count) {
    /* Go downstream */
    for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
      univdiagonal = this3->all_univdiagonals_gplus[i];
      while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
	k_gplus5--;
      }

      /* (5,6) Look downstream for 5' end matches */
      debug11(printf("Finding downstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gplus[i]; */
      pathstore = this3->all_pathstores_gplus[i];

      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);

      auxinfo_mate = (Auxinfo_T) NULL;
      if (pathstore->mate_count >= max_count) {
	k = pathstore->mate_bestk;
	assert(k >= 0);

	univdiagonal_mate = this5->exhaustive_gplus[k];
	qstart = this5->exhaustive_qstart_gplus[k];
	qend = this5->exhaustive_qend_gplus[k];
	debug11(printf("(5,6) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gplus[k]));

	auxinfo_mate =
	  Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					 this5->exhaustive_gplus,this5->exhaustive_qstart_gplus,
					 this5->exhaustive_qend_gplus,this5->nexhaustive_gplus,
					 univdiagpool,auxinfopool,/*method*/LOCAL_MATE);

      } else if (localdb != NULL &&
		 querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		 (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							   novel_diagonals_alloc,localdb_alloc,this5,
							   /*queryptr*/queryuc_ptr_5,querylength5,
							   low_univdiagonal,/*high_univdiagonal*/univdiagonal,
							   /*query_compress*/query5_compress_fwd,
							   /*plusp*/true,/*genestrand*/0,
							   genomebits,nmismatches_allowed_5)) != 0) {
	qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					    univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					    /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					   univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					   /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	debug11(printf("(5) highest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal - univdiagonal_mate,qstart,qend));
	auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	pathstore_mate = Pathstore_new(pathstorepool);
      }
	
      if (auxinfo_mate != NULL) {
	if (pathstore->solvedp == false) {
	  /* Solve for anchor if necessary */
	  solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				     &(*found_score_3),univdiagonal,&(this3->all_auxinfo_gplus[i]),pathstore,

				     queryseq3,/*queryptr*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,
				     this3,knownsplicing,knownindels,
				 
				     mismatch_positions_alloc_3,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     /*query_compress*/query3_compress_fwd,
				     query3_compress_fwd,query3_compress_rev,
				 
				     nmismatches_allowed_3,/*genestrand*/0,
				     auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen3,/*plusp*/true,/*first_read_p*/false,/*lowp*/false,
				     /*set_best_paths_p*/true);
	}
	  
	if (auxinfo->complete_sense_paths != NULL || auxinfo->complete_antisense_paths != NULL ||
	    auxinfo->unextended_sense_paths != NULL || auxinfo->unextended_antisense_paths != NULL) {
	  /* Solve for mate */
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &auxinfo_mate->unextended_sense_paths,&auxinfo_mate->unextended_antisense_paths,
				       &auxinfo_mate->complete_sense_paths,&auxinfo_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryuc_ptr_5,
				       /*query_compress*/query5_compress_fwd,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/true,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  auxinfo_mate->solvedp = true;
	  Auxinfo_set_best_paths(auxinfo_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(14) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*auxinfoL*/auxinfo_mate,/*auxinfoH*/auxinfo,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/false);
	}

	Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }


  if (max_count_gminus >= max_count) {
    /* Go upstream */
    for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
      univdiagonal = this3->all_univdiagonals_gminus[i];
      while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
	k_gminus5++;
      }
      
      /* (7,8) Look upstream for 5' end matches */
      debug11(printf("Finding upstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gminus[i]; */
      pathstore = this3->all_pathstores_gminus[i];

      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);

      auxinfo_mate = (Auxinfo_T) NULL;
      if (pathstore->mate_count >= max_count) {
	k = pathstore->mate_bestk;
	assert(k >= 0);

	univdiagonal_mate = this5->exhaustive_gminus[k];
	qstart = this5->exhaustive_qstart_gminus[k];
	qend = this5->exhaustive_qend_gminus[k];
	debug11(printf("(7,8) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,qstart,qend,
		       this5->exhaustive_counts_gminus[k]));

	auxinfo_mate =
	  Kmer_compute_auxinfo_univdiags(univdiagonal_mate,qstart,qend,k,
					 this5->exhaustive_gminus,this5->exhaustive_qstart_gminus,
					 this5->exhaustive_qend_gminus,this5->nexhaustive_gminus,
					 univdiagpool,auxinfopool,/*method*/LOCAL_MATE);
	pathstore_mate = Pathstore_new(pathstorepool);

      } else if (localdb != NULL &&
		 querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		 (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							  novel_diagonals_alloc,localdb_alloc,this5,
							  /*queryptr*/queryrc5,querylength5,
							  /*low_univdiagonal*/univdiagonal,high_univdiagonal,
							  /*query_compress*/query5_compress_rev,
							  /*plusp*/false,/*genestrand*/0,
							  genomebits,nmismatches_allowed_5)) != 0) {
	qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					    univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					    /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					   univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					   /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	debug11(printf("(7,8) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - univdiagonal,qstart,qend));
	auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	pathstore_mate = Pathstore_new(pathstorepool);
      }
	
      if (auxinfo_mate != NULL) {
	if (auxinfo->solvedp == false) {
	  /* Solve for anchor if necessary */
	  solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				     &(*found_score_3),univdiagonal,&(this3->all_auxinfo_gminus[i]),pathstore,

				     queryseq3,/*queryptr*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,
				     this3,knownsplicing,knownindels,
				 
				     mismatch_positions_alloc_3,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     /*query_compress*/query3_compress_rev,
				     query3_compress_fwd,query3_compress_rev,
				 
				     nmismatches_allowed_3,/*genestrand*/0,
				 
				     auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen3,/*plusp*/false,/*first_read_p*/false,/*lowp*/true,
				     /*set_best_paths_p*/true);
	}

	if (pathstore->complete_sense_paths != NULL || pathstore->complete_antisense_paths != NULL ||
	    pathstore->unextended_sense_paths != NULL || pathstore->unextended_antisense_paths != NULL) {
	  /* Solve for mate */
	  Path_solve_from_univdiagonal(&(*found_score_5),
				       
				       &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
				       &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
				       univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryrc5,
				       /*query_compress*/query5_compress_rev,
				       query5_compress_fwd,query5_compress_rev,
				       /*plusp*/false,querylength5,mismatch_positions_alloc_5,
				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       this5,knownsplicing,knownindels,
				       nmismatches_allowed_5,/*paired_end_p*/true,/*first_read_p*/true,
				       /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
				       intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				       spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	  pathstore_mate->solvedp = true;
	  Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  
	  debug(printf("(15) Calling make_pathpairs\n"));
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/false);
	}

	Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool);
	Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }

  return pathpairs;
}
#endif


static void
solve_all (bool *completep, int *found_score, Shortread_T queryseq,
	   char *queryuc_ptr, char *queryrc, int querylength,
	   T this, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	   int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
	   Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	   
	   Compress_T query_compress_fwd, Compress_T query_compress_rev,
	   int nmismatches_allowed, bool first_read_p,

	   Auxinfopool_T auxinfopool,
	   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	   Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
	   Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
	   Hitlistpool_T hitlistpool, Spliceendsgen_T spliceendsgen) {

  Univcoord_T univdiagonal;
  Pathstore_T pathstore;
  bool lowp;
  bool complete_sense_p = false, complete_antisense_p = false;
  int i;

  /* plus */
  lowp = (first_read_p == true) ? true : false;
  for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
    univdiagonal = this->all_univdiagonals_gplus[i];
    pathstore = this->all_pathstores_gplus[i];
    solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
			       &(*found_score),univdiagonal,&(this->all_auxinfo_gplus[i]),pathstore,

			       queryseq,/*queryptr*/queryuc_ptr,queryuc_ptr,queryrc,querylength,
			       this,knownsplicing,knownindels,
				 
			       mismatch_positions_alloc,
			       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			       /*query_compress*/query_compress_fwd,
			       query_compress_fwd,query_compress_rev,
				 
			       nmismatches_allowed,/*genestrand*/0,
			       auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
			       hitlistpool,spliceendsgen,/*plusp*/true,first_read_p,lowp,
			       /*set_best_paths_p*/false);
    if (complete_sense_p == true || complete_antisense_p == true) {
      *completep = true;
    }
  }

  /* minus */
  lowp = (first_read_p == true) ? false : true;
  for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
    univdiagonal = this->all_univdiagonals_gminus[i];
    pathstore = this->all_pathstores_gminus[i];
    solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
			       &(*found_score),univdiagonal,&(this->all_auxinfo_gminus[i]),pathstore,

			       queryseq,/*queryptr*/queryrc,queryuc_ptr,queryrc,querylength,
			       this,knownsplicing,knownindels,

			       mismatch_positions_alloc,
			       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			       /*query_compress*/query_compress_rev,
			       query_compress_fwd,query_compress_rev,
			       
			       nmismatches_allowed,/*genestrand*/0,
			       auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
			       hitlistpool,spliceendsgen,/*plusp*/false,first_read_p,lowp,
			       /*set_best_paths_p*/false);
  }

  if (complete_sense_p == true && complete_antisense_p == true) {
    *completep = true;

    /* plus */
    lowp = (first_read_p == true) ? true : false;
    for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
      pathstore = this->all_pathstores_gplus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }
    
    /* minus */
    lowp = (first_read_p == true) ? false : true;
    for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      pathstore = this->all_pathstores_gminus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }

  } else if (complete_sense_p == true) {
    *completep = true;

    /* plus */
    lowp = (first_read_p == true) ? true : false;
    for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
      pathstore = this->all_pathstores_gplus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }
    
    /* minus */
    lowp = (first_read_p == true) ? false : true;
    for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      pathstore = this->all_pathstores_gminus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }
    
  } else if (complete_antisense_p == true) {
    *completep = true;

    /* plus */
    lowp = (first_read_p == true) ? true : false;
    for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
      pathstore = this->all_pathstores_gplus[i];
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }
    
    /* minus */
    lowp = (first_read_p == true) ? false : true;
    for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      pathstore = this->all_pathstores_gminus[i];
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/true);
    }

  } else if (*completep == true) {
    /* Already found complete solutions */

  } else {
    /* *completep = false; -- Keep track over all univdiagonals */

    /* plus */
    lowp = (first_read_p == true) ? true : false;
    for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
      pathstore = this->all_pathstores_gplus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/false);
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/false);
    }
    
    /* minus */
    lowp = (first_read_p == true) ? false : true;
    for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      pathstore = this->all_pathstores_gminus[i];
      Pathstore_set_best_sense_paths(pathstore,hitlistpool,/*only_complete_p*/false);
      Pathstore_set_best_antisense_paths(pathstore,hitlistpool,/*only_complete_p*/false);
    }
  }


  return;
}


/* final_pairtype can be CONCORDANT_TRANSLOCATIONS, CONCORDANT, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG, UNPAIRED */
Pathpair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Path_T **patharray5, int *npaths5_primary, int *npaths5_altloc, int *first_absmq5, int *second_absmq5,
		    Path_T **patharray3, int *npaths3_primary, int *npaths3_altloc, int *first_absmq3, int *second_absmq3,
		    Shortread_T queryseq5, Shortread_T queryseq3, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Chrpos_T pairmax_linear,
		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		    Transcriptpool_T transcriptpool,
		    Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3, Pass_T pass) {
  Pathpair_T *pathpairarray;
  List_T pathpairs = NULL, unresolved_pathpairs = NULL, unresolved_pathpairs_from_exhaustive = NULL;
  Method_T last_method_5, last_method_3;

  List_T sense_paths5 = NULL, antisense_paths5 = NULL, sense_paths3 = NULL, antisense_paths3 = NULL;

  /* For single-end reads */
  List_T sense_paths5_gplus, sense_paths5_gminus,
    antisense_paths5_gplus, antisense_paths5_gminus,
    sense_paths3_gplus, sense_paths3_gminus,
    antisense_paths3_gplus, antisense_paths3_gminus;

  /* List_T local_sense_paths5 = NULL, local_antisense_paths5 = NULL,
     local_sense_paths3 = NULL, local_antisense_paths3 = NULL; */

  /* List_T geneplus_pathpairs, geneminus_pathpairs; */
  /* List_T geneplus_paths5, geneplus_paths3, geneminus_paths5, geneminus_paths3; */

  /* Maximum number of localdb regions possible, based on region size of 65536 bp */
  int max_localdb_nregions = (positive_gap_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  int max_nintersections = max_localdb_nregions * LOCALDB_REGION_SIZE;
  int sufficient_score_5, sufficient_score_3, querylength5, querylength3;

  bool complete5_p, complete3_p;
  int found_score_paired, found_score_5, found_score_3;
  bool any_imperfect_ends5_p = false, any_imperfect_ends3_p = false;
  bool solved_all_p = false;

  int nmismatches_filter_5, nmismatches_filter_3;
  int mincoverage_filter_5, mincoverage_filter_3;
  int nmismatches_allowed_5, nmismatches_allowed_3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *queryrc5, *queryrc3;
  int *mismatch_positions_alloc_5, *mismatch_positions_alloc_3;
  Univdiag_T *novel_univdiags_alloc;
  Univcoord_T *novel_diagonals_alloc;
  unsigned short *localdb_alloc;
  Compress_T query5_compress_fwd, query5_compress_rev, query3_compress_fwd, query3_compress_rev;
  T this5, this3;


  if ((querylength5 = Shortread_fulllength(queryseq5)) < index1part + index1interval - 1 ||
      (querylength3 = Shortread_fulllength(queryseq3)) < index1part + index1interval - 1) {
    return (Pathpair_T *) NULL;
  }

  queryuc_ptr_5 = Shortread_queryuc_ptr(queryseq5);
  queryuc_ptr_3 = Shortread_queryuc_ptr(queryseq3);
  queryrc5 = Shortread_queryrc(queryseq5);
  queryrc3 = Shortread_queryrc(queryseq3);

  this5 = Stage1_new(queryuc_ptr_5,querylength5,/*first_read_p*/true);
  this3 = Stage1_new(queryuc_ptr_3,querylength3,/*first_read_p*/false);

  found_score_5 = querylength5;
  found_score_3 = querylength3;

  /* nmismatches_allowed means nmismatches_search and is not specified
     by the user.  The user-specified value for -m represents
     nmismatches_filter */
  /* TODO: make this dependent upon the defect rate */
  nmismatches_allowed_5 = querylength5/20; /* was querylength/index1part */
  nmismatches_allowed_3 = querylength3/20; /* was querylength/index1part */

  if (user_nmismatches_filter_float < 0.0) {
    /* Not specified, so don't filter */
    nmismatches_filter_5 = querylength5;
    nmismatches_filter_3 = querylength3;
  } else if (user_nmismatches_filter_float < 1.0) {
    nmismatches_filter_5 = (int) rint(user_nmismatches_filter_float * (double) querylength5);
    nmismatches_filter_3 = (int) rint(user_nmismatches_filter_float * (double) querylength3);
  } else {
    nmismatches_filter_5 = nmismatches_filter_3 = (int) user_nmismatches_filter_float;
  }

  if (user_mincoverage_filter_float <= 0.0) {
    /* Not specified, so don't filter */
    mincoverage_filter_5 = 0;
    mincoverage_filter_3 = 0;
  } else if (user_mincoverage_filter_float <= 1.0) {
    /* Assuming that --min-coverage=1 must mean 1.0 and not a coverage of 1 bp */
    mincoverage_filter_5 = (int) rint(user_mincoverage_filter_float * (double) querylength5);
    mincoverage_filter_3 = (int) rint(user_mincoverage_filter_float * (double) querylength3);
  } else {
    mincoverage_filter_5 = mincoverage_filter_3 = (int) user_mincoverage_filter_float;
  }

  mismatch_positions_alloc_5 = (int *) MALLOC((querylength5+MISMATCH_EXTRA)*sizeof(int));
  mismatch_positions_alloc_3 = (int *) MALLOC((querylength3+MISMATCH_EXTRA)*sizeof(int));

  /* Intersection indices are paired */
  novel_univdiags_alloc = Univdiag_new_array(2 * max_nintersections);
  novel_diagonals_alloc = (Univcoord_T *) MALLOC(2 * max_nintersections * sizeof(Univcoord_T));
  MALLOC_ALIGN(localdb_alloc,max_nintersections * sizeof(unsigned short));

  query5_compress_fwd = Compress_new_fwd(queryuc_ptr_5,querylength5);
  query5_compress_rev = Compress_new_rev(queryuc_ptr_5,querylength5);
  query3_compress_fwd = Compress_new_fwd(queryuc_ptr_3,querylength3);
  query3_compress_rev = Compress_new_rev(queryuc_ptr_3,querylength3);



  if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    /* Not implemented yet */
    fprintf(stderr,"Nonstranded modes not yet implemented\n");
    exit(9);

    pathpairarray = (Pathpair_T *) NULL;
#ifdef TO_FIX
    geneplus_pathpairs = paired_read(&geneplus_abort_pairing_p,&geneplus_hits5,&geneplus_hits3,
				     &geneplus_samechr,&geneplus_conc_transloc,
				     queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				     knownsplicing,knownindels,mismatch_positions_alloc_5,mismatch_positions_alloc_3,
				     novel_diagonals_alloc,localdb_alloc,this5,this3,
				     query5_compress_fwd,query5_compress_rev,
				     query3_compress_fwd,query3_compress_rev,
				     max_insertionlen_5,max_insertionlen_3,max_deletionlen_5,max_deletionlen_3,
				     /*genestrand*/+1,pairmax_linear,overall_max_distance_5,overall_max_distance_3,
				     overall_end_distance_5,overall_end_distance_3,
				     max_mismatches_refalt_5,max_mismatches_refalt_3,
				     max_mismatches_ref_5,max_mismatches_ref_3,
				     min_coverage_5,min_coverage_3,
				     intlistpool,univcoordlistpool,listpool,univdiagpool,
				     hitlistpool,pathpool,vectorpool,
				     transcriptpool,spliceendsgen5,spliceendsgen3,pass);

    geneminus_pathpairs = paired_read(&geneminus_abort_pairing_p,&geneminus_hits5,&geneminus_hits3,
				      &geneminus_samechr,&geneminus_conc_transloc,
				      queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				      knownsplicing,knownindels,mismatch_positions_alloc_5,mismatch_positions_alloc_3,
				      novel_diagonals_alloc,localdb_alloc,this5,this3,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      max_insertionlen_5,max_insertionlen_3,max_deletionlen_5,max_deletionlen_3,
				      /*genestrand*/+2,pairmax_linear,overall_max_distance_5,overall_max_distance_3,
				      overall_end_distance_5,overall_end_distance_3,
				      max_mismatches_refalt_5,max_mismatches_refalt_3,
				      max_mismatches_ref_5,max_mismatches_ref_3,
				      min_coverage_5,min_coverage_3,
				      intlistpool,univcoordlistpool,listpool,univdiagpool,
				      hitlistpool,pathpool,vectorpool,
				      transcriptpool,spliceendsgen5,spliceendsgen3,pass);
#endif

  } else { /*mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED */
    found_score_5 = querylength5;
    found_score_3 = querylength3;
    found_score_paired = querylength5 + querylength3 + 1; /* Add 1 to enter loop */
    sufficient_score_5 = querylength5/20;
    sufficient_score_3 = querylength3/20;
    complete5_p = complete3_p = false;

    if (transcriptome_align_p == true) {
      pathpairs = paired_search_tr(&found_score_paired,&found_score_5,&found_score_3,
				   sufficient_score_5,sufficient_score_3,
				   
				   queryseq5,queryseq3,querylength5,querylength3,/*genestrand*/0,
				   mismatch_positions_alloc_5,mismatch_positions_alloc_3,
				   this5,this3,knownsplicing,
				   query5_compress_fwd,query5_compress_rev,
				   query3_compress_fwd,query3_compress_rev,

				   nmismatches_allowed_5,nmismatches_allowed_3,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,

				   trdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool);
    }

    if (genome_align_p == true) {
      last_method_5 = last_method_3 = METHOD_INIT;
      /* Comparing found_score_paired against found_score_5 +
	 found_score_3 tells us whether we have missed any good
	 ends */
      /* Exit if found_score_5 <= sufficient_score && found_score_3 <= sufficient_score &&
	 found_score_paired <= found_score_5 + found_score_3 */
      while ((found_score_paired > found_score_5 + found_score_3 ||
	      found_score_5 > sufficient_score_5 || found_score_3 > sufficient_score_3) &&
	     (last_method_5 < SEGMENT1 || last_method_3 < SEGMENT1)) {
	pathpairs =
	  paired_search_univdiagonals(&complete5_p,&complete3_p,
				      &found_score_paired,&found_score_5,&found_score_3,
				      &last_method_5,&last_method_3,
				      &any_imperfect_ends5_p,&any_imperfect_ends3_p,
				      
				      &unresolved_pathpairs,pathpairs,this5,this3,queryseq5,queryseq3,
				      queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				      knownsplicing,knownindels,mismatch_positions_alloc_5,mismatch_positions_alloc_3,
				      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      nmismatches_allowed_5,nmismatches_allowed_3,/*genestrand*/0,
				      
				      nmismatches_filter_5,nmismatches_filter_3,
				      mincoverage_filter_5,mincoverage_filter_3,
				      
				      repetitive_ef64,auxinfopool,pathstorepool,univdiagpool,
				      intlistpool,uintlistpool,univcoordlistpool,listpool,
				      pathpool,transcriptpool,vectorpool,hitlistpool,
				      spliceendsgen5,spliceendsgen3);
	debug(printf("found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     found_score_paired,found_score_5,found_score_3));
      }
    }
  }


  /* Want to continue searching until we found a resolved pathpair,
     but can now consider the unresolved ones.  Pathpair_eval_and_sort
     will call Pathpair_resolve */

  if (pathpairs == NULL) {
    /* Univdiagonals solved paths vs Exhaustive */
    /* Does yield results, since we might have been too strict on complete_p for the mates */
    pathpairs = paths5_mates(&found_score_paired,&found_score_3,&unresolved_pathpairs,
			     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_3,queryrc3,querylength3,
			     knownsplicing,knownindels,mismatch_positions_alloc_3,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     this5,this3,query3_compress_fwd,query3_compress_rev,
			     nmismatches_allowed_3,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen3,/*only_complete5_p*/complete5_p);
    
    pathpairs = paths3_mates(&found_score_paired,&found_score_5,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_5,queryrc5,querylength5,
			     knownsplicing,knownindels,mismatch_positions_alloc_5,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     this5,this3,query5_compress_fwd,query5_compress_rev,
			     nmismatches_allowed_5,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen5,/*only_complete3_p*/complete3_p);

    debug(Stage1_list_all_univdiagonals(this5));
    debug(Stage1_list_all_univdiagonals(this3));

    /* debug(Stage1_list_exhaustive(this5)); */
    /* debug(Stage1_list_exhaustive(this3)); */

#if 0
    /* Demonstration that this step yields results */
    if (pathpairs != NULL) {
      for (List_T p = pathpairs; p != NULL; p = List_next(p)) {
	Pathpair_print((Pathpair_T) List_head(p));
      }
      exit(0);
    }
#endif
  }


#if 0
  if (pathpairs == NULL) {
    /* Univdiagonals vs Exhaustive.  This procedure is too
       time-consuming.  Better to go to paired_search_exhaustive,
       which takes counts into account */
    pathpairs =
      univdiagonals5_mates(&found_score_paired,&found_score_5,&found_score_3,&unresolved_pathpairs,
				  
			   pathpairs,queryseq5,queryseq3,
			   queryuc_ptr_5,queryrc5,querylength5,
			   queryuc_ptr_3,queryrc3,querylength3,
			   knownsplicing,knownindels,
			   mismatch_positions_alloc_5,mismatch_positions_alloc_3,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this5,this3,
			   query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_rev,
			   nmismatches_allowed_5,nmismatches_allowed_3,

			   nmismatches_filter_5,nmismatches_filter_3,
			   mincoverage_filter_5,mincoverage_filter_3,

			   univdiagpool,auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			   spliceendsgen5,spliceendsgen3);
    
    pathpairs =
      univdiagonals3_mates(&found_score_paired,&found_score_5,&found_score_3,&unresolved_pathpairs,
				  
			   pathpairs,queryseq5,queryseq3,
			   queryuc_ptr_5,queryrc5,querylength5,
			   queryuc_ptr_3,queryrc3,querylength3,
			   knownsplicing,knownindels,
			   mismatch_positions_alloc_5,mismatch_positions_alloc_3,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,this5,this3,
			   query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_rev,
			   nmismatches_allowed_5,nmismatches_allowed_3,

			   nmismatches_filter_5,nmismatches_filter_3,
			   mincoverage_filter_5,mincoverage_filter_3,

			   univdiagpool,auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			   spliceendsgen5,spliceendsgen3);

    debug(Stage1_list_all_univdiagonals(this5));
    debug(Stage1_list_all_univdiagonals(this3));

    /* debug(Stage1_list_exhaustive(this5)); */
    /* debug(Stage1_list_exhaustive(this3)); */
  }
#endif


  debug(printf("Beginning paired_search_exhaustive\n"));
  if (pathpairs == NULL) {
    /* Exhaustive vs Exhaustive.  Takes counts into account, so should
       be better than univdiagonals_mates procedures */
    pathpairs = paired_search_exhaustive(&found_score_paired,&found_score_5,&found_score_3,
					 &unresolved_pathpairs_from_exhaustive,
					 queryseq5,queryuc_ptr_5,queryrc5,
					 queryseq3,queryuc_ptr_3,queryrc3,
					 querylength5,querylength3,

					 query5_compress_fwd,query5_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 
					 this5,this3,knownsplicing,knownindels,
					 nmismatches_allowed_5,nmismatches_allowed_3,

					 nmismatches_filter_5,nmismatches_filter_3,
					 mincoverage_filter_5,mincoverage_filter_3,

					 auxinfopool,pathstorepool,univdiagpool,intlistpool,uintlistpool,
					 univcoordlistpool,listpool,pathpool,vectorpool,
					 transcriptpool,hitlistpool,
					 spliceendsgen5,spliceendsgen3);
  }
  if (pathpairs == NULL) {
    pathpairs = paired_search_exhaustive(&found_score_paired,&found_score_5,&found_score_3,
					 &unresolved_pathpairs_from_exhaustive,
					 queryseq5,queryuc_ptr_5,queryrc5,
					 queryseq3,queryuc_ptr_3,queryrc3,
					 querylength5,querylength3,

					 query5_compress_fwd,query5_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 
					 this5,this3,knownsplicing,knownindels,
					 nmismatches_allowed_5,nmismatches_allowed_3,

					 nmismatches_filter_5,nmismatches_filter_3,
					 mincoverage_filter_5,mincoverage_filter_3,

					 auxinfopool,pathstorepool,univdiagpool,intlistpool,uintlistpool,
					 univcoordlistpool,listpool,pathpool,vectorpool,
					 transcriptpool,hitlistpool,
					 spliceendsgen5,spliceendsgen3);
  }
  debug(printf("Done with paired_search_exhaustive\n"));

  /* Exhaustive should succeed, unless one querylength < index1part */
  debug(printf("After exhaustive, have %d pathpairs, %d unresolved pathpairs, %d unresolved from exhaustive\n",
	       List_length(pathpairs),List_length(unresolved_pathpairs),
	       List_length(unresolved_pathpairs_from_exhaustive)));

  if (pathpairs != NULL) {
    Pathpair_gc(&unresolved_pathpairs_from_exhaustive,intlistpool,univcoordlistpool,listpool,pathpool,
		transcriptpool,hitlistpool);
  } else if (List_length(unresolved_pathpairs_from_exhaustive) > 100) {
    Pathpair_gc(&unresolved_pathpairs_from_exhaustive,intlistpool,univcoordlistpool,listpool,pathpool,
		transcriptpool,hitlistpool);
  } else {
    pathpairs = unresolved_pathpairs_from_exhaustive;
    /* unresolved_pathpairs_from_exhaustive = (List_T) NULL; */
  }

  if (pathpairs == NULL) {
    solved_all_p = true;

    /* Solve all 5' solutions and find mates locally */
    solve_all(&complete5_p,&found_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
	      this5,knownsplicing,knownindels,mismatch_positions_alloc_5,
	      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		
	      query5_compress_fwd,query5_compress_rev,
	      nmismatches_allowed_5,/*first_read_p*/true,
		
	      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
	      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
	      spliceendsgen5);

    pathpairs = paths5_mates(&found_score_paired,&found_score_3,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_3,queryrc3,querylength3,
			     knownsplicing,knownindels,mismatch_positions_alloc_3,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     this5,this3,query3_compress_fwd,query3_compress_rev,
			     nmismatches_allowed_3,
				     
			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
				     
			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen3,/*only_complete5_p*/complete5_p);

    /* Solve all 3' solutions and find mates locally */
    solve_all(&complete3_p,&found_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
	      this3,knownsplicing,knownindels,mismatch_positions_alloc_3,
	      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		
	      query3_compress_fwd,query3_compress_rev,
	      nmismatches_allowed_3,/*first_read_p*/false,
		
	      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
	      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
	      spliceendsgen3);
    
    pathpairs = paths3_mates(&found_score_paired,&found_score_5,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_5,queryrc5,querylength5,
			     knownsplicing,knownindels,mismatch_positions_alloc_5,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     this5,this3,query5_compress_fwd,query5_compress_rev,
			     nmismatches_allowed_5,
			     
			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
				     
			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen5,/*only_complete3_p*/complete3_p);
    solved_all_p = true;
  }

  debug(Stage1_list_all_univdiagonals(this5));
  debug(Stage1_list_all_univdiagonals(this3));

  /* Previously had a call to paired_search_anchored, but now we have mates procedures */

  debug(printf("Have %d pathpairs, %d unresolved pathpairs, complete5p %d, complete3p %d\n",
	       List_length(pathpairs),List_length(unresolved_pathpairs),
	       complete5_p,complete3_p));

  if (pathpairs != NULL) {
    Pathpair_gc(&unresolved_pathpairs,intlistpool,univcoordlistpool,listpool,pathpool,
		transcriptpool,hitlistpool);
  } else if (List_length(unresolved_pathpairs) > 100) {
    /* Otherwise, Pathpair_eval_and_sort takes too long */
    debug(printf("Deleting unresolved pathpairs because there are too many\n"));
    Pathpair_gc(&unresolved_pathpairs,intlistpool,univcoordlistpool,listpool,pathpool,
		transcriptpool,hitlistpool);
  } else {
    pathpairs = unresolved_pathpairs;
  }
    
  if (pathpairs == NULL) {
    if (solved_all_p == true) {
      debug(printf("Already solved all\n"));
      /* paired_search5_mates and paired_search3_mates were called previously after solve_all */

    } else {
      debug(printf("Calling solve_all on 5' end\n"));
      solve_all(&complete5_p,&found_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
		this5,knownsplicing,knownindels,mismatch_positions_alloc_5,
		novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		
		query5_compress_fwd,query5_compress_rev,
		nmismatches_allowed_5,/*first_read_p*/true,
		
		auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		spliceendsgen5);
      debug(printf("Done with solve_all on 5' end\n"));

      debug(printf("Calling solve_all on 3' end\n"));
      solve_all(&complete3_p,&found_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
		this3,knownsplicing,knownindels,mismatch_positions_alloc_3,
		novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		
		query3_compress_fwd,query3_compress_rev,
		nmismatches_allowed_3,/*first_read_p*/false,
		
		auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		spliceendsgen3);
      debug(printf("Done with solve_all on 3' end\n"));
    }

    Stage1_collect_paths(&sense_paths5_gplus,&sense_paths5_gminus,
			 &antisense_paths5_gplus,&antisense_paths5_gminus,
			 this5,hitlistpool);
    Stage1_collect_paths(&sense_paths3_gplus,&sense_paths3_gminus,
			 &antisense_paths3_gplus,&antisense_paths3_gminus,
			 this3,hitlistpool);

    sense_paths5 = List_append(sense_paths5_gplus,sense_paths5_gminus);
    antisense_paths5 = List_append(antisense_paths5_gplus,antisense_paths5_gminus);
    sense_paths3 = List_append(sense_paths3_gplus,sense_paths3_gminus);
    antisense_paths3 = List_append(antisense_paths3_gplus,antisense_paths3_gminus);

#if 0
    if (splicingp == true) {
      /* Inner path fusions require further work, and is too time-consuming */
      if ((pathpairs = find_inner_fusions(&found_score_5,&found_score_3,
				     
					  sense_paths5,antisense_paths5,
					  sense_paths3,antisense_paths3,
					  queryseq5,queryseq3,queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					  knownsplicing,novel_univdiags_alloc,novel_diagonals_alloc,
					  localdb_alloc,this5,this3,
					  query5_compress_fwd,query5_compress_rev,
					  query3_compress_fwd,query3_compress_rev,
					  nmismatches_allowed_5,nmismatches_allowed_3,
				     
					  intlistpool,uintlistpool,univcoordlistpool,listpool,
					  pathpool,vectorpool,hitlistpool,transcriptpool)) != NULL) {
	Hitlistpool_free_list(&sense_paths5,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	Hitlistpool_free_list(&antisense_paths5,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	Hitlistpool_free_list(&sense_paths5,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	Hitlistpool_free_list(&antisense_paths5,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
      }
    }
#endif
  }

  /* If pathpairs == NULL, then we have collected all of the singlepaths */

#if 0
  /* Previously called Pathpair_filter, but this can lead to poor answers */
  /* This call to Pathpair_filter can rule out correct loci that have
     not been extended, so need to make sure extensions work */
  debug(printf("Have %d pathpairs before filtering\n",List_length(pathpairs)));
  debug(print_pathpairs_contents(pathpairs));
  if (pathpairs != NULL) {
    pathpairs = Pathpair_filter(pathpairs,
				intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  }
  debug(printf("Have %d pathpairs after filtering\n",List_length(pathpairs)));
  debug(print_pathpairs_contents(pathpairs));
#endif
    
    /* Check for coverage */
  *npaths_primary = List_length(pathpairs);
  *npaths_altloc = 0;	/* TODO: Determine whether any paths are on the altloc chromosome */
    
  if (*npaths_primary == 0) {
    pathpairarray = (Pathpair_T *) NULL;
  } else {
    debug(printf("Starting Pathpair_eval_and_sort\n"));
    pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
    pathpairarray = Pathpair_eval_and_sort(&found_score_5,&found_score_3,
					   &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					   pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					   /*quality_string_5*/Shortread_quality_string(queryseq5),
					   /*quality_string_3*/Shortread_quality_string(queryseq3),
					   mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					   knownsplicing,knownindels,
					   
					   nmismatches_allowed_5,nmismatches_allowed_3,
					   nmismatches_filter_5,nmismatches_filter_3,
					   mincoverage_filter_5,mincoverage_filter_3,
					   querylength5,querylength3,
					   univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,/*filterp*/true);
    debug(printf("Done with Pathpair_eval_and_sort\n"));
  }
    
  if ((*npaths_primary) + (*npaths_altloc) > 0) {
    *final_pairtype = CONCORDANT;

  } else {
    if (pathpairs == NULL) {
      /* We never generated pathpairpairarray, and we already collected singlepaths */
      assert(solved_all_p == true);

    } else {
      FREE_OUT(pathpairarray);
      Hitlistpool_free_list(&pathpairs,hitlistpool
			    hitlistpool_trace(__FILE__,__LINE__));

      if (solved_all_p == true) {
	/* paired_search5_mates and paired_search3_mates were called previously after solve_all */
	
      } else {
	solve_all(&complete5_p,&found_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
		  this5,knownsplicing,knownindels,mismatch_positions_alloc_5,
		  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		  
		  query5_compress_fwd,query5_compress_rev,
		  nmismatches_allowed_5,/*first_read_p*/true,
		  
		  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		  pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		  spliceendsgen5);
	
	solve_all(&complete3_p,&found_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
		  this3,knownsplicing,knownindels,mismatch_positions_alloc_3,
		  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		  
		  query3_compress_fwd,query3_compress_rev,
		  nmismatches_allowed_3,/*first_read_p*/false,
		  
		  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		  pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		  spliceendsgen3);
      }

      /* These were not collected before, since we tried to use Pathpair_eval_and_sort */

      Stage1_collect_paths(&sense_paths5_gplus,&sense_paths5_gminus,
			   &antisense_paths5_gplus,&antisense_paths5_gminus,
			   this5,hitlistpool);
      Stage1_collect_paths(&sense_paths3_gplus,&sense_paths3_gminus,
			   &antisense_paths3_gplus,&antisense_paths3_gminus,
			   this3,hitlistpool);

      sense_paths5 = List_append(sense_paths5_gplus,sense_paths5_gminus);
      antisense_paths5 = List_append(antisense_paths5_gplus,antisense_paths5_gminus);
      sense_paths3 = List_append(sense_paths3_gplus,sense_paths3_gminus);
      antisense_paths3 = List_append(antisense_paths3_gplus,antisense_paths3_gminus);
    }

    debug(printf("Calling consolidate_results\n"));

    pathpairarray =
      consolidate_results(&found_score_5,&found_score_3,&(*final_pairtype),
			  &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
			  &(*patharray5),&(*npaths5_primary),&(*npaths5_altloc),&(*first_absmq5),&(*second_absmq5),
			  &(*patharray3),&(*npaths3_primary),&(*npaths3_altloc),&(*first_absmq3),&(*second_absmq3),
			  sense_paths5,antisense_paths5,sense_paths3,antisense_paths3,
			    
			  queryseq5,queryseq3,
			  queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
			    
			  knownsplicing,knownindels,novel_univdiags_alloc,novel_diagonals_alloc,
			  localdb_alloc,this5,this3,
			  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
			  query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
			  
			  nmismatches_allowed_5,nmismatches_allowed_3,
			  nmismatches_filter_5,nmismatches_filter_3,mincoverage_filter_5,mincoverage_filter_3,
			  
			  intlistpool,uintlistpool,univcoordlistpool,listpool,univdiagpool,pathpool,vectorpool,
			  hitlistpool,transcriptpool,spliceendsgen5,spliceendsgen3);

    Hitlistpool_free_list(&sense_paths5,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&sense_paths3,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&antisense_paths5,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&antisense_paths3,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }

  Compress_free(&query3_compress_rev);
  Compress_free(&query3_compress_fwd);
  Compress_free(&query5_compress_rev);
  Compress_free(&query5_compress_fwd);

  FREE_ALIGN(localdb_alloc);
  FREE(novel_diagonals_alloc);
  FREE(novel_univdiags_alloc[0]);
  FREE(novel_univdiags_alloc);

  FREE(mismatch_positions_alloc_3);
  FREE(mismatch_positions_alloc_5);

  Stage1_free(&this3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
	      intlistpool,uintlistpool,univcoordlistpool,
	      listpool,pathpool,trpathpool,transcriptpool,hitlistpool,
	      /*free_paths_p*/true);
  Stage1_free(&this5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
	      intlistpool,uintlistpool,univcoordlistpool,
	      listpool,pathpool,trpathpool,transcriptpool,hitlistpool,
	      /*free_paths_p*/true);

  /* FREE(queryrc3); -- Taken from Shortread */
  /* FREE(queryrc5); -- Taken from Shortread */

  debug(printf("Returning with final_pairtype %s\n",Pairtype_string(*final_pairtype)));
  debug(print_pathpairarray_contents(pathpairarray,/*n*/(*npaths_primary) + (*npaths_altloc)));

  return pathpairarray;
}


void
Stage1hr_paired_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       Genomebits_T genomebits_in, Localdb_T localdb_in, EF64_T chromosome_ef64_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       int max_deletionlen, int max_insertlength, Chrpos_T shortsplicedist, bool splicingp_in,
		       int maxpaths_search_in, int maxpaths_report_in,
		       bool *circularp_in, int pairmax_linear_in, int pairmax_circular_in) {

  mode = mode_in;
  index1part = index1part_in;
  index1interval = index1interval_in;
  index1part_tr = index1part_tr_in;

  transcriptome = transcriptome_in;
  genome_align_p = genome_align_p_in;
  transcriptome_align_p = transcriptome_align_p_in;

  genomebits = genomebits_in;
  localdb = localdb_in;

  chromosome_ef64 = chromosome_ef64_in;

  user_nmismatches_filter_float = user_nmismatches_filter_float_in;
  user_mincoverage_filter_float = user_mincoverage_filter_float_in;

  concordance_distance = (Chrpos_T) max_insertlength + shortsplicedist;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  splicingp = splicingp_in;
  maxpaths_search = maxpaths_search_in;
  maxpaths_report = maxpaths_report_in;

  circularp = circularp_in;
  pairmax_linear = pairmax_linear_in;
  pairmax_circular = pairmax_circular_in;

  return;
}


void
Stage1hr_paired_pass2_setup (int max_insertlength, Chrpos_T shortsplicedist) {

  concordance_distance = (Chrpos_T) max_insertlength + shortsplicedist;

  return;
}

