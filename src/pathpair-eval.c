static char rcsid[] = "$Id: 6f5fe675402d07bced4e52181df60ac4ba436472 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "pathpair-eval.h"
#include "path-eval.h"		/* For Path_mark_alignment */

#include "path-solve.h"
#include "path-fusion.h"
#include "path-trim.h"

#include <stdio.h>
#include <math.h>		/* For rint */
#include "fastlog.h"		/* For fasterexp */

#include "assert.h"
#include "list.h"
#include "genomebits_count.h"
#include "junction.h"
#include "mapq.h"
#include "outputtype.h"
#include "stage1hr.h"


static bool *circularp;
static bool *chrsubsetp;
static bool *altlocp;

static Outputtype_T output_type;
static bool splicingp;
static bool resolve_inner_p;
static bool want_random_p;
static bool allow_soft_clips_p;
static Chrpos_T max_insertlength;

static Transcriptome_T transcriptome;


/* Splice goodness now determined by Splice_accept_p */
/* #define OMIT_BAD_SPLICES 1 */

#define NMATCHES_FACTOR 8
#define INSERTLENGTH_FACTOR 1.5
#define OUTERLENGTH_FACTOR 1.5

#define INNER_NEEDED 15


/* Pathpair_eval_and_sort */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

#define T Path_T

/* For DNA-seq reads, do not want to use insertlength or outerlength */
/* Duplicates with respect to method have already been taken care of */
/* Ignore sensedir, so we keep both sensedirs if they have equivalent matches */
static int
Pathpair_method_cmp (const void *x, const void *y) {
  Pathpair_T a = * (Pathpair_T *) x;
  Pathpair_T b = * (Pathpair_T *) y;

  int cmpL, cmpH;
  int nindels_a, nindels_b;
  int method_a, method_b;


  cmpL = Path_nmatches_cmp(a->pathL,b->pathL);
  cmpH = Path_nmatches_cmp(a->pathH,b->pathH);

  if (cmpL <= 0 && cmpH < 0) {
    return -1;
  } else if (cmpL < 0 && cmpH <= 0) {
    return -1;
  } else if (cmpL >= 0 && cmpH > 0) {
    return +1;
  } else if (cmpL > 0 && cmpH >= 0) {
    return +1;
  } else if (cmpL < 0 && cmpH > 0) {
    /* Conflicting information */
    return 0;
  } else if (cmpL > 0 && cmpH < 0) {
    /* Conflicting information */
    return 0;
  }

  nindels_a = Path_nindels(a->pathL) + Path_nindels(a->pathH);
  nindels_b = Path_nindels(b->pathL) + Path_nindels(b->pathH);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;

  } else {
    method_a = (int) a->pathL->method + (int) a->pathH->method;
    method_b = (int) b->pathL->method + (int) b->pathH->method;

    if (method_a > method_b) {
      return -1;
    } else if (method_b > method_a) {
      return +1;
    } else {
      return 0;
    }
  }
}


static int
Pathpair_local_cmp (const void *x, const void *y) {
  Pathpair_T a = * (Pathpair_T *) x;
  Pathpair_T b = * (Pathpair_T *) y;

  int cmpL, cmpH;
  int chopped_nmatches_a, chopped_nmatches_b, nmatches_a, nmatches_b;

  int coverage_a, coverage_b;
  int nindels_a, nindels_b;
  int nalts_a = 0, nalts_b = 0;
#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif

#if 0
  double junction_splice_prob_a, junction_splice_prob_b,
    total_splice_prob_a, total_splice_prob_b;;
#endif
  bool overlappingp_x, overlappingp_y;
  bool insertlength_knownp_x, insertlength_knownp_y,
    outerlength_knownp_x, outerlength_knownp_y;

  /* Fusions can eliminate good results.  Can judge all by their nmatches */
  bool fusionp_x = false, fusionp_y = false;

#ifdef DEBUG8
  printf("Pathpair_local_cmp comparing\n");
  Pathpair_print(a);
  Pathpair_print(b);
#endif

  if (a->pathL->fusion_querystart_junction != NULL ||
      a->pathL->fusion_queryend_junction != NULL ||
      a->pathH->fusion_querystart_junction != NULL ||
      a->pathH->fusion_queryend_junction != NULL) {
    fusionp_x = true;
  }
  if (b->pathL->fusion_querystart_junction != NULL ||
      b->pathL->fusion_queryend_junction != NULL ||
      b->pathH->fusion_querystart_junction != NULL ||
      b->pathH->fusion_queryend_junction != NULL) {
    fusionp_y = true;
  }

  if (fusionp_x == false && fusionp_y == true) {
    return -1;
  } else if (fusionp_y == false && fusionp_x == true) {
    return +1;
  }

  coverage_a = Path_coverage(a->pathL) + Path_coverage(a->pathH);
  coverage_b = Path_coverage(b->pathL) + Path_coverage(b->pathH);
  
#ifdef OMIT_BAD_SPLICES
  nbadsplices_a = Pathpair_nbadsplices(a);
  nbadsplices_b = Pathpair_nbadsplices(b);
#endif

#if 0
  junction_splice_prob_a = a->pathL->junction_splice_prob + a->pathH->junction_splice_prob;
  junction_splice_prob_b = b->pathL->junction_splice_prob + b->pathH->junction_splice_prob;
  total_splice_prob_a = a->pathL->total_splice_prob + a->pathH->total_splice_prob;
  total_splice_prob_b = b->pathL->total_splice_prob + b->pathH->total_splice_prob;
#endif
    
  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  }

#ifdef OMIT_BAD_SPLICES				/* Allow for local_cmp */
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif
    
#if 0
  /* With better control over long splices, this is no longer useful */
  /* Prioritize outerlength before nmatches, which could vary by a single base at the end */
  outerlength_knownp_x = Pathpair_outerlength_knownp(a);
  outerlength_knownp_y = Pathpair_outerlength_knownp(b);
  if (outerlength_knownp_x == false) {
    /* Fall through */
  } else if (outerlength_knownp_y == false) {
    /* Fall through */
  } else if (a->outerlength * OUTERLENGTH_FACTOR < b->outerlength) {
    return -1;
  } else if (b->outerlength * OUTERLENGTH_FACTOR < a->outerlength) {
    return +1;
  }
#endif

  cmpL = Path_nmatches_cmp(a->pathL,b->pathL);
  cmpH = Path_nmatches_cmp(a->pathH,b->pathH);
  debug8(printf("cmpL %d, cmpH %d\n",cmpL,cmpH));

  if (cmpL <= 0 && cmpH < 0) {
    return -1;
  } else if (cmpL < 0 && cmpH <= 0) {
    return -1;
  } else if (cmpL >= 0 && cmpH > 0) {
    return +1;
  } else if (cmpL > 0 && cmpH >= 0) {
    return +1;
  }

  nmatches_a = a->pathL->nmatches + a->pathH->nmatches;
  nmatches_b = b->pathL->nmatches + b->pathH->nmatches;
  
#if 0
  if (nmatches_a > nmatches_b + NMATCHES_FACTOR) {
    /* Use nmatches instead of chopped */
  } else if (nmatches_b > nmatches_a + NMATCHES_FACTOR) {
    /* Use nmatches instead of chopped */
  } else {
    /* Resolve by adding the values */
    chopped_nmatches_a = Path_eval_nmatches_chopped(a->pathL) + Path_eval_nmatches_chopped(a->pathH);
    chopped_nmatches_b = Path_eval_nmatches_chopped(b->pathL) + Path_eval_nmatches_chopped(b->pathH);
    if (chopped_nmatches_a > chopped_nmatches_b) {
      return -1;
    } else if (chopped_nmatches_b > chopped_nmatches_a) {
      return +1;
    }
  }
#endif

  if (nmatches_a > nmatches_b) {
    return -1;
  } else if (nmatches_b > nmatches_a) {
    return +1;
  } else {
    /* Equal.  Resolve below */
  }

  if (a->transcript_concordant_p > b->transcript_concordant_p) {
    return -1;
  } else if (b->transcript_concordant_p > a->transcript_concordant_p) {
    return +1;
  }

  overlappingp_x = Pathpair_overlappingp(a);
  overlappingp_y = Pathpair_overlappingp(b);
  if (overlappingp_x == false && overlappingp_y == true) {
    return -1;
  } else if (overlappingp_y == false && overlappingp_x == true) {
    return +1;
  }

  nindels_a = Path_nindels(a->pathL) + Path_nindels(a->pathH);
  nindels_b = Path_nindels(b->pathL) + Path_nindels(b->pathH);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  }

  if (a->pathL->qstart_alts != NULL) {
    nalts_a++;
  }
  if (a->pathL->qend_alts != NULL) {
    nalts_a++;
  }
  if (a->pathH->qstart_alts != NULL) {
    nalts_a++;
  }
  if (a->pathH->qend_alts != NULL) {
    nalts_a++;
  }

  if (b->pathL->qstart_alts != NULL) {
    nalts_b++;
  }
  if (b->pathL->qend_alts != NULL) {
    nalts_b++;
  }
  if (b->pathH->qstart_alts != NULL) {
    nalts_b++;
  }
  if (b->pathH->qend_alts != NULL) {
    nalts_b++;
  }

  if (nalts_a < nalts_b) {
    return -1;
  } else if (nalts_b < nalts_a) {
    return +1;
  }

#if 0
    /* Already compared sense/antisense */
    /* We want to keep all possible transcript results at a locus */
  if (junction_splice_prob_a > junction_splice_prob_b) {
    return -1;
  } else if (junction_splice_prob_b > junction_splice_prob_a) {
    return +1;
  } else if (total_splice_prob_a > total_splice_prob_b) {
    return -1;
  } else if (total_splice_prob_b > total_splice_prob_a) {
    return +1;
  }
#endif

  insertlength_knownp_x = Pathpair_insertlength_knownp(a);
  insertlength_knownp_y = Pathpair_insertlength_knownp(b);

  if (insertlength_knownp_x == false) {
    /* Fall through */
  } else if (insertlength_knownp_y == false) {
    /* Fall through */
  } else if (a->insertlength < b->insertlength) {
    return -1;
  } else if (b->insertlength < a->insertlength) {
    return +1;
  }

  return 0;
}


static int
Pathpair_global_cmp (const void *x, const void *y) {
  Pathpair_T a = * (Pathpair_T *) x;
  Pathpair_T b = * (Pathpair_T *) y;

  int cmpL, cmpH;
  int chopped_nmatches_a, chopped_nmatches_b, nmatches_a, nmatches_b;

  int coverage_a, coverage_b;
  /* bool insertlength_knownp_x, insertlength_knownp_y; */
  Chrpos_T insertlength_x, insertlength_y;
  Chrpos_T outerlength_x, outerlength_y;

#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif


  coverage_a = Path_coverage(a->pathL) + Path_coverage(a->pathH);
  coverage_b = Path_coverage(b->pathL) + Path_coverage(b->pathH);

#ifdef OMIT_BAD_SPLICES
  nbadsplices_a = Pathpair_nbadsplices(a);
  nbadsplices_b = Pathpair_nbadsplices(b);
#endif

  /* Fusions can eliminate good results.  Can judge all by their nmatches */
  bool fusionp_x = false, fusionp_y = false;
  if (a->pathL->fusion_querystart_junction != NULL ||
      a->pathL->fusion_queryend_junction != NULL ||
      a->pathH->fusion_querystart_junction != NULL ||
      a->pathH->fusion_queryend_junction != NULL) {
    fusionp_x = true;
  }
  if (b->pathL->fusion_querystart_junction != NULL ||
      b->pathL->fusion_queryend_junction != NULL ||
      b->pathH->fusion_querystart_junction != NULL ||
      b->pathH->fusion_queryend_junction != NULL) {
    fusionp_y = true;
  }

  if (fusionp_x == false && fusionp_y == true) {
    return -1;
  } else if (fusionp_y == false && fusionp_x == true) {
    return +1;
  }

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  }

#ifdef OMIT_BAD_SPLICES
    /* Do not consider bad splices for global_cmp */
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif
    
  cmpL = Path_nmatches_cmp(a->pathL,b->pathL);
  cmpH = Path_nmatches_cmp(a->pathH,b->pathH);

  if (cmpL <= 0 && cmpH < 0) {
    return -1;
  } else if (cmpL < 0 && cmpH <= 0) {
    return -1;
  } else if (cmpL >= 0 && cmpH > 0) {
    return +1;
  } else if (cmpL > 0 && cmpH >= 0) {
    return +1;
  }

  nmatches_a = a->pathL->nmatches + a->pathH->nmatches;
  nmatches_b = b->pathL->nmatches + b->pathH->nmatches;

#if 0
  if (nmatches_a > nmatches_b + NMATCHES_FACTOR) {
    /* Use nmatches instead of chopped */
  } else if (nmatches_b > nmatches_a + NMATCHES_FACTOR) {
    /* Use nmatches instead of chopped */
  } else {
    /* Resolve by adding the values */
    chopped_nmatches_a = Path_eval_nmatches_chopped(a->pathL) + Path_eval_nmatches_chopped(a->pathH);
    chopped_nmatches_b = Path_eval_nmatches_chopped(b->pathL) + Path_eval_nmatches_chopped(b->pathH);
    if (chopped_nmatches_a > chopped_nmatches_b) {
      return -1;
    } else if (chopped_nmatches_b > chopped_nmatches_a) {
      return +1;
    }
  }
#endif

  if (nmatches_a > nmatches_b) {
    return -1;
  } else if (nmatches_b > nmatches_a) {
    return +1;
  } else {
    /* Equal.  Resolve below */
  }


  /* Previously made transcript_concordant_p more important than
     nmatches, but this allows a good solution to be missed if it
     does not correspond to a known transcript */
  if (a->transcript_concordant_p > b->transcript_concordant_p) {
    return -1;
  } else if (b->transcript_concordant_p > a->transcript_concordant_p) {
    return +1;

  } else if (splicingp == false) {
    return 0;

  } else {
    /* We want to keep all possible transcript results across loci */
    /* insertlength_knownp_x = Pathpair_insertlength_knownp(a); */
    /* insertlength_knownp_y = Pathpair_insertlength_knownp(b); */
    insertlength_x = Pathpair_insertlength(a);
    insertlength_y = Pathpair_insertlength(b);

    if (insertlength_x > max_insertlength && insertlength_y > max_insertlength) {
      /* Don't use insertlength */
    } else if (insertlength_x * INSERTLENGTH_FACTOR < insertlength_y) {
      return -1;
    } else if (insertlength_y * INSERTLENGTH_FACTOR < insertlength_x) {
      return +1;
    } else {
      return 0;
    }

    outerlength_x = Path_genomichigh(a->pathH) - Path_genomiclow(a->pathL);
    outerlength_y = Path_genomichigh(b->pathH) - Path_genomiclow(b->pathL);

    if (outerlength_x * OUTERLENGTH_FACTOR < outerlength_y) {
      return -1;
    } else if (outerlength_y * OUTERLENGTH_FACTOR < outerlength_x) {
      return +1;
    } else {
      return 0;
    }
  }
}


/* filtering is creating an issue between pathpairs != NULL, stopping
   path creation, and then an empty pathpairarray */

Pathpair_T *
Pathpair_eval_and_sort (int *found_score_5, int *found_score_3,
			int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
			Pathpair_T *pathpairarray, int npaths, Stage1_T stage1_5, Stage1_T stage1_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Shortread_T queryseq5, Shortread_T queryseq3,
			char *queryuc_ptr_5, char *queryrc5, char *queryuc_ptr_3, char *queryrc3, 
			char *quality_string_5, char *quality_string_3,
			int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, Knownsplicing_T knownsplicing,
			Knownindels_T knownindels,

			int nmismatches_allowed_5, int nmismatches_allowed_3,
			int nmismatches_filter_5, int nmismatches_filter_3,
			int mincoverage_filter_5, int mincoverage_filter_3,
			int querylength5, int querylength3,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
			Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, bool filterp) {
  int pathi;
  float maxlik, loglik;

  float total, qual;
  int mapq_score;

  int randomi, i, j, k, l;
  Pathpair_T temp, pathpair;
  Path_T path5, path3;
  Univcoord_T high, low, distal5, distal3;
  int tr_insertlength;

  bool valid5p, valid3p;
  bool possible_fusion5_p = false, possible_fusion3_p = false;

  List_T paths, fusion_paths_5, fusion_paths_3, p, q;
  List_T unresolved_pathpairs;
  double best_avg_junction_sense_prob, best_avg_junction_antisense_prob,
    best_end_sense_prob, best_end_antisense_prob, prob;
  /* int best_querystart5, best_queryend5, best_querystart3, best_queryend3, querystart, queryend; */
  int chop_qstart_L, chop_qend_L, chop_qstart_H, chop_qend_H, qstart, qend;

  
#ifdef DEBUG8
  printf("Entered Pathpair_eval_and_sort with %d pathpairs\n",npaths);
  for (pathi = 0; pathi < npaths; pathi++) {
    printf("Pathpair (sense %d and %d): %p and %p\n",
	   pathpairarray[pathi]->pathL->sensedir,pathpairarray[pathi]->pathH->sensedir,
	   pathpairarray[pathi]->pathL,pathpairarray[pathi]->pathH);
  }
  printf("\n");
#endif

  if (npaths == 0) {
    /* Skip */
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;
    return (Pathpair_T *) NULL;

  }

  /* Have already called Path_extend */

  /* 0.  Sort by structure to remove duplicates */
  if (npaths > 1) {
    qsort(pathpairarray,npaths,sizeof(T),Pathpair_structure_cmp);
    
    k = 0;
    i = 0;
    while (i < npaths) {
      chop_qstart_L = Intlist_head(pathpairarray[i]->pathL->endpoints);
      chop_qend_L = Intlist_last_value(pathpairarray[i]->pathL->endpoints);
      chop_qstart_H = Intlist_head(pathpairarray[i]->pathH->endpoints);
      chop_qend_H = Intlist_last_value(pathpairarray[i]->pathH->endpoints);

      j = i + 1;
      while (j < npaths && Pathpair_structure_cmp(&(pathpairarray[j]),&(pathpairarray[i])) == 0) {
	if ((qstart = Intlist_head(pathpairarray[j]->pathL->endpoints)) > chop_qstart_L) {
	  chop_qstart_L = qstart;
	}
	if ((qend = Intlist_last_value(pathpairarray[j]->pathL->endpoints)) < chop_qend_L) {
	  chop_qend_L = qend;
	}
	if ((qstart = Intlist_head(pathpairarray[j]->pathH->endpoints)) > chop_qstart_H) {
	  chop_qstart_H = qstart;
	}
	if ((qend = Intlist_last_value(pathpairarray[j]->pathH->endpoints)) < chop_qend_H) {
	  chop_qend_H = qend;
	}
	j++;
      }
      debug8(printf("Found an identical group by structure (except sensedir) of %d paths, chop region %d..%d and %d..%d => Re-sorting by method_cmp\n",
		    j - i,chop_qstart_L,chop_qend_L,chop_qstart_H,chop_qend_H));
      
      for (l = i; l < j; l++) {
	pathpairarray[l]->pathL->chop_qstart = chop_qstart_L;
	pathpairarray[l]->pathL->chop_qend = chop_qend_L;
	pathpairarray[l]->pathH->chop_qstart = chop_qstart_H;
	pathpairarray[l]->pathH->chop_qend = chop_qend_H;
      }

      qsort(&(pathpairarray[i]),j - i,sizeof(Pathpair_T),Pathpair_method_cmp);
      debug8(printf("(0) Keeping by method_cmp\n")); debug8(Pathpair_print(pathpairarray[i]));
      pathpairarray[k++] = pathpairarray[i];
      
      for (l = i + 1; l < j; l++) {
	debug8(printf("(0) Eliminating by method_cmp\n")); debug8(Pathpair_print(pathpairarray[l]));
	pathpair = pathpairarray[l];
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
      
      i = j;
    }
    npaths = k;
  }
  debug8(printf("After removing duplicates (except for sensedir), have %d paths\n",npaths));


  /* 1.  Unalias circular alignments and trim chrbounds */
  for (i = 0; i < npaths; i++) {
    pathpair = pathpairarray[i];
    path5 = pathpair->path5;
    path3 = pathpair->path3;
    if (circularp[pathpair->pathL->chrnum] == true) {
      Path_trim_circular_unalias_pair(pathpair->pathL,pathpair->pathH);
    }
    Path_trim_chrbounds(path5,query5_compress_fwd,query5_compress_rev,
			intlistpool,univcoordlistpool,listpool,pathpool,
			transcriptpool);
    Path_trim_chrbounds(path3,query3_compress_fwd,query3_compress_rev,
			intlistpool,univcoordlistpool,listpool,pathpool,
			transcriptpool);
  }
  

  /* 2.  Trim and resolve inner splices */
  /* Previously sorted, but not necessary to group the pathpairs, and we sort at step 4 anyway */
  /* qsort(pathpairarray,npaths,sizeof(Pathpair_T),Pathpair_interval_cmp); */
  
  if (resolve_inner_p == true && npaths < 100) {
    /* Do not try to resolve, which is time-comsuming, if there are many pathpairs */
    /* Note: Need to handle cases where the reads overlap, so trimming needs to be determined by the distal univdiagonals */
    k = 0;
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      path5 = pathpair->path5;
      path3 = pathpair->path3;
      valid5p = valid3p = true;
      
      debug8(printf("Resolving paths %p and %p\n",path5,path3));
      if (pathpair->plusp == true) {
	distal5 = Path_genomiclow(path5);
	distal3 = Path_genomichigh(path3);
	
	/* Need to shift 5' univdiagonals down by querylength5 */
	if (Univcoordlist_last_value(path5->univdiagonals) < (Univcoord_T) querylength5) {
	  high = 0;
	} else {
	  high = Univcoordlist_last_value(path5->univdiagonals) - querylength5;
	}
	debug8(printf("high5 %u vs distal3 %u\n",high,distal3));

	if (0 && path5->transcriptome_method_p == true) {
	  /* Previously did not trim, but we could have added a splice to the end */
	} else if (high > distal3) { /* was path3->main_univdiagonal */
	  debug8(printf("Need to trim high part of 5' read, since %u > %u\n",high,distal3));
	  debug8(Path_print(path5));
	  debug8(Path_print(path3));
	  valid5p = Path_trim_qend_trimdiag(path5,query5_compress_fwd,query5_compress_rev,
					    intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,
					    querylength5,/*trimdiag*/distal3);
	  debug8(printf("Result of trimming\n"));
	  debug8(Path_print(path5));
	}
	
	low = Univcoordlist_head(path3->univdiagonals);
	debug8(printf("low3 %u vs distal5 %u\n",low,distal5));

	if (0 && path3->transcriptome_method_p == true) {
	  /* Previously did not trim, but we could have added a splice to the end */
	} else if (low < distal5) { /* was path5->main_univdiagonal */
	  debug8(printf("Need to trim low part of 3' read, since %u < %u\n",low,distal5));
	  debug8(Path_print(path5));
	  debug8(Path_print(path3));
	  valid3p = Path_trim_qstart_trimdiag(path3,query3_compress_fwd,query3_compress_rev,
					      intlistpool,univcoordlistpool,
					      listpool,pathpool,transcriptpool,
					      /*trimdiag*/distal5);
	  debug8(printf("Result of trimming\n"));
	  debug8(Path_print(path3));
	}
	
	if (valid5p == false || valid3p == false) {
	  /* Skip: delete below */
#if 0
	} else if (Path_resolved_qend_p(path5) == true && Path_resolved_qstart_p(path3) == true) {
	  debug8(printf("Have a resolved inner: ")); debug8(Path_print(path5)); debug8(Path_print(path3));
#endif
	} else {
	  /* Need to call even if resolved, because of overreach issues */
	  debug8(printf("Calling Pathpair_resolve\n"));
	  Pathpair_resolve(&(*found_score_5),&(*found_score_3),
			   pathpair,pathpair->plusp,/*genestrand*/0,
			   query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
			   /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,queryuc_ptr_5,queryuc_ptr_3,
			   /*querylengthL*/querylength5,/*querylengthH*/querylength3,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   stage1_5,stage1_3,knownsplicing,
			   nmismatches_allowed_5,nmismatches_allowed_3,
			   intlistpool,univcoordlistpool,listpool,pathpool,
			   transcriptpool,univdiagpool,vectorpool);
	  debug8(printf("Done with Pathpair_resolve\n"));
	}
	
      } else {
	/* minus */
	distal5 = Path_genomichigh(path5);
	distal3 = Path_genomiclow(path3);
	
	low = Univcoordlist_head(path5->univdiagonals);
	debug8(printf("low5 %u vs distal3 %u\n",low,distal3));

	if (0 && path3->transcriptome_method_p == true) {
	  /* Previously did not trim, but we could have added a splice to the end */
	} else if (low < distal3) { /* was path3->main_univdiagonal */
	  debug8(printf("Need to trim low part of 5' read, since %u < %u\n",low,distal3));
	  debug8(Path_print(path5));
	  debug8(Path_print(path3));
	  valid5p = Path_trim_qstart_trimdiag(path5,query5_compress_fwd,query5_compress_rev,
					      intlistpool,univcoordlistpool,
					      listpool,pathpool,transcriptpool,
					      /*trimdiag*/distal3);
	  debug8(printf("Result of trimming\n"));
	  debug8(Path_print(path5));
	}
	
	/* Need to shift 3' univdiagonals down by querylength3 */
	if (Univcoordlist_last_value(path3->univdiagonals) < (Univcoord_T) querylength3) {
	  high = 0;
	} else {
	  high = Univcoordlist_last_value(path3->univdiagonals) - querylength3;
	}
	debug8(printf("high3 %u vs distal5 %u\n",low,distal5));
	
	if (0 && path5->transcriptome_method_p == true) {
	  /* Previously did not trim, but we could have added a splice to the end */
	} else if (high > distal5) { /* was path5->main_univdiagonal */
	  debug8(printf("Need to trim high part of 3' read, since %u > %u\n",high,distal5));
	  debug8(Path_print(path5));
	  debug8(Path_print(path3));
	  valid3p = Path_trim_qend_trimdiag(path3,query3_compress_fwd,query3_compress_rev,
					    intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,
					    querylength3,/*trimdiag*/distal5);
	  debug8(printf("Result of trimming\n"));
	  debug8(Path_print(path3));
	}
	
	if (valid5p == false || valid3p == false) {
	  /* Skip: delete below */
#if 0
	} else if (Path_resolved_qend_p(path3) == true && Path_resolved_qstart_p(path5) == true) {
	  debug8(printf("Have a resolved inner: ")); debug8(Path_print(path5)); debug8(Path_print(path3));
#endif
	} else {
	  /* Need to call even if resolved, because of overreach issues */
	  debug8(printf("Calling Pathpair_resolve\n"));
	  Pathpair_resolve(&(*found_score_5),&(*found_score_3),
			   pathpair,pathpair->plusp,/*genestrand*/0,
			   query3_compress_rev,query3_compress_fwd,query3_compress_rev,
			   query5_compress_rev,query5_compress_fwd,query5_compress_rev,
			   /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,queryrc3,queryrc5,
			   /*querylengthL*/querylength3,/*querylengthH*/querylength5,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   stage1_3,stage1_5,knownsplicing,
			   nmismatches_allowed_3,nmismatches_allowed_5,
			   intlistpool,univcoordlistpool,listpool,pathpool,
			   transcriptpool,univdiagpool,vectorpool);
	  debug8(printf("Done with Pathpair_resolve\n"));
	}
      }
      
      if (valid5p == false || valid3p == false) {
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      } else {
	pathpairarray[k++] = pathpairarray[i];
      }
    }

    if ((npaths = k) == 0) {
      *npaths_primary = *npaths_altloc = 0;
      *first_absmq = 0;
      *second_absmq = 0;
      return (Pathpair_T *) NULL;
    }

    /* 2b.  Can get structurally identical results after Pathpair_resolve, so need to repeat step 1 */
    if (npaths > 1) {
      qsort(pathpairarray,npaths,sizeof(T),Pathpair_structure_cmp);
      
      k = 0;
      i = 0;
      while (i < npaths) {
	chop_qstart_L = Intlist_head(pathpairarray[i]->pathL->endpoints);
	chop_qend_L = Intlist_last_value(pathpairarray[i]->pathL->endpoints);
	chop_qstart_H = Intlist_head(pathpairarray[i]->pathH->endpoints);
	chop_qend_H = Intlist_last_value(pathpairarray[i]->pathH->endpoints);

	j = i + 1;
	while (j < npaths && Pathpair_structure_cmp(&(pathpairarray[j]),&(pathpairarray[i])) == 0) {
	  if ((qstart = Intlist_head(pathpairarray[j]->pathL->endpoints)) > chop_qstart_L) {
	    chop_qstart_L = qstart;
	  }
	  if ((qend = Intlist_last_value(pathpairarray[j]->pathL->endpoints)) < chop_qend_L) {
	    chop_qend_L = qend;
	  }
	  if ((qstart = Intlist_head(pathpairarray[j]->pathH->endpoints)) > chop_qstart_H) {
	    chop_qstart_H = qstart;
	  }
	  if ((qend = Intlist_last_value(pathpairarray[j]->pathH->endpoints)) < chop_qend_H) {
	    chop_qend_H = qend;
	  }
	  j++;
	}
	debug8(printf("Found an identical group by structure (except sensedir) of %d paths, chop region %d..%d and %d..%d => Re-sorting by method_cmp\n",
		      j - i,chop_qstart_L,chop_qend_L,chop_qstart_H,chop_qend_H));
	
	for (l = i; l < j; l++) {
	  pathpairarray[l]->pathL->chop_qstart = chop_qstart_L;
	  pathpairarray[l]->pathL->chop_qend = chop_qend_L;
	  pathpairarray[l]->pathH->chop_qstart = chop_qstart_H;
	  pathpairarray[l]->pathH->chop_qend = chop_qend_H;
	}

	qsort(&(pathpairarray[i]),j - i,sizeof(Pathpair_T),Pathpair_method_cmp);
	debug8(printf("(2) Keeping by method_cmp\n")); debug8(Pathpair_print(pathpairarray[i]));
	pathpairarray[k++] = pathpairarray[i];
	
	for (l = i + 1; l < j; l++) {
	  debug8(printf("(2) Eliminating by method_cmp\n")); debug8(Pathpair_print(pathpairarray[l]));
	  pathpair = pathpairarray[l];
	  Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	}
	
	i = j;
      }
      npaths = k;
    }
    debug8(printf("After removing duplicates (except for sensedir), have %d paths\n",npaths));
  }
  

  /* Note: Don't maximize on nmatches until we find fusions, because a fusion splice may need to truncate nmatches */

  /* 3.  Check all solutions for a possible outer fusion */
  if (splicingp == true) {
    /* 3a.  Look for possibilities */
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      if (pathpair->plusp == true) {
	if (Path_unextended_qstart_p(pathpair->path5,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
	  possible_fusion5_p = true;
	}
	if (Path_unextended_qend_p(pathpair->path3,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
	  possible_fusion3_p = true;
	}
      } else {
	if (Path_unextended_qend_p(pathpair->path5,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
	  possible_fusion5_p = true;
	}
	if (Path_unextended_qstart_p(pathpair->path3,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
	  possible_fusion3_p = true;
	}
      }
      i++;
    }
  }
  

  debug8(printf("possible_fusion5_p %d, possible_fusion3_p %d\n",
		possible_fusion5_p,possible_fusion3_p));
  
  if (possible_fusion5_p == true || possible_fusion3_p == true) {
    /* 3b.  Look for outer fusions */
    
    if (possible_fusion5_p == true) {
      debug8(Stage1_list_extension(stage1_5));
    }
    
    if (possible_fusion3_p == true) {
      debug8(Stage1_list_extension(stage1_3));
    }
    
    paths = (List_T) NULL;
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      
      /* Push original pathpair onto paths */
      paths = Hitlist_push(paths,hitlistpool, (void *) pathpair
			   hitlistpool_trace(__FILE__,__LINE__));
      
      /* Look for fusion paths */
      if (pathpair->plusp == true) {
	/* plus */
	if (Path_unextended_qstart_p(pathpair->path5,/*endtrim_allowed*/INNER_NEEDED,/*allow_ambig_p*/false) == false) {
	  fusion_paths_5 = (List_T) NULL;
	} else if (pathpair->path5->sensedir == SENSE_ANTI) {
	  /* Skip: requires sense */
	  fusion_paths_5 = (List_T) NULL;
	} else {
	  fusion_paths_5 = Path_fusion_outer_querystart_plus(&(*found_score_5),/*main*/pathpair->path5,stage1_5,

							     query5_compress_fwd,query5_compress_rev,
							     queryseq5,querylength5,knownsplicing,
							     pathpair->path5->genestrand,nmismatches_allowed_5,
							     intlistpool,uintlistpool,univcoordlistpool,
							     listpool,univdiagpool,pathpool,vectorpool,
							     transcriptpool,hitlistpool);
	}

	if (Path_unextended_qend_p(pathpair->path3,/*endtrim_allowed*/INNER_NEEDED,/*allow_ambig_p*/false) == false) {
	  fusion_paths_3 = (List_T) NULL;
	} else if (pathpair->path3->sensedir == SENSE_ANTI) {
	  /* Skip: requires sense */
	  fusion_paths_3 = (List_T) NULL;
	} else {
	  fusion_paths_3 = Path_fusion_outer_queryend_plus(&(*found_score_3),/*main*/pathpair->path3,stage1_3,

							   query3_compress_fwd,query3_compress_rev,
							   queryseq3,querylength3,knownsplicing,
							   pathpair->path3->genestrand,nmismatches_allowed_3,
							   intlistpool,uintlistpool,univcoordlistpool,
							   listpool,univdiagpool,pathpool,vectorpool,
							   transcriptpool,hitlistpool);
	}

	if (fusion_paths_5 == NULL && fusion_paths_3 == NULL) {
	  /* Skip */

	} else if (fusion_paths_5 != NULL && fusion_paths_3 == NULL) {
	  for (p = fusion_paths_5; p != NULL; p = List_next(p)) {
	    paths = Hitlist_push(paths,hitlistpool,
				 (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								  /*pathL*/(Path_T) List_head(p),/*pathH*/pathpair->path3,
								  /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
								  nmismatches_filter_5,nmismatches_filter_3,
								  mincoverage_filter_5,mincoverage_filter_3,
								  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								  transcriptpool,hitlistpool,/*check_inner_p*/false,
								  /*copyLp (was false)*/true,/*copyHp*/true)
				 hitlistpool_trace(__FILE__,__LINE__));
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_5,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

	} else if (fusion_paths_5 == NULL && fusion_paths_3 != NULL) {
	  for (p = fusion_paths_3; p != NULL; p = List_next(p)) {
	    paths = Hitlist_push(paths,hitlistpool,
				 (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								  /*pathL*/pathpair->path5,/*pathH*/(Path_T) List_head(p),
								  /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
								  nmismatches_filter_5,nmismatches_filter_3,
								  mincoverage_filter_5,mincoverage_filter_3,
								  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								  transcriptpool,hitlistpool,/*check_inner_p*/false,
								  /*copyLp*/true,/*copyHp (was false)*/true)
				 hitlistpool_trace(__FILE__,__LINE__));
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_3,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

	} else {
	  for (p = fusion_paths_5; p != NULL; p = List_next(p)) {
	    for (q = fusion_paths_3; q != NULL; q = List_next(q)) {
	      paths = Hitlist_push(paths,hitlistpool,
				   (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								    /*pathL*/(Path_T) List_head(p),/*pathH*/(Path_T) List_head(q),
								    /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
								    nmismatches_filter_5,nmismatches_filter_3,
								    mincoverage_filter_5,mincoverage_filter_3,
								    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								    transcriptpool,hitlistpool,/*check_inner_p*/false,
								    /*copyLp (was false)*/true,/*copyHp (was false)*/true)
				   hitlistpool_trace(__FILE__,__LINE__));
	    }
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_5,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
	  Path_gc(&fusion_paths_3,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
	}

      } else {
	/* minus */
	if (Path_unextended_qstart_p(pathpair->path3,/*endtrim_allowed*/INNER_NEEDED,/*allow_ambig_p*/false) == false) {
	  fusion_paths_3 = (List_T) NULL;
	} else if (pathpair->path3->sensedir == SENSE_FORWARD) {
	  /* Skip: requires antisense */
	  fusion_paths_3 = (List_T) NULL;
	} else {
	  fusion_paths_3 = Path_fusion_outer_querystart_minus(&(*found_score_3),/*main*/pathpair->path3,stage1_3,

							      query3_compress_fwd,query3_compress_rev,
							      queryseq3,querylength3,knownsplicing,
							      pathpair->path3->genestrand,nmismatches_allowed_3,
							      intlistpool,uintlistpool,univcoordlistpool,
							      listpool,univdiagpool,pathpool,vectorpool,
							      transcriptpool,hitlistpool);
	}

	if (Path_unextended_qend_p(pathpair->path5,/*endtrim_allowed*/INNER_NEEDED,/*allow_ambig_p*/false) == false) {
	  fusion_paths_5 = (List_T) NULL;
	} else if (pathpair->path5->sensedir == SENSE_FORWARD) {
	  /* Skip: requires antisense */
	  fusion_paths_5 = (List_T) NULL;
	} else {
	  fusion_paths_5 = Path_fusion_outer_queryend_minus(&(*found_score_5),/*main*/pathpair->path5,stage1_5,

							    query5_compress_fwd,query5_compress_rev,
							    queryseq5,querylength5,knownsplicing,
							    pathpair->path5->genestrand,nmismatches_allowed_5,
							    intlistpool,uintlistpool,univcoordlistpool,
							    listpool,univdiagpool,pathpool,vectorpool,
							    transcriptpool,hitlistpool);
	}

	if (fusion_paths_5 == NULL && fusion_paths_3 == NULL) {
	  /* Skip */

	} else if (fusion_paths_5 != NULL && fusion_paths_3 == NULL) {
	  for (p = fusion_paths_5; p != NULL; p = List_next(p)) {
	    paths = Hitlist_push(paths,hitlistpool,
				 (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								  /*pathL*/pathpair->path3,/*pathH*/(Path_T) List_head(p),
								  /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
								  nmismatches_filter_5,nmismatches_filter_3,
								  mincoverage_filter_5,mincoverage_filter_3,
								  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								  transcriptpool,hitlistpool,/*check_inner_p*/false,
								  /*copyLp*/true,/*copyHp (was false)*/true)
				 hitlistpool_trace(__FILE__,__LINE__));
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_5,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

	} else if (fusion_paths_5 == NULL && fusion_paths_3 != NULL) {
	  for (p = fusion_paths_3; p != NULL; p = List_next(p)) {
	    paths = Hitlist_push(paths,hitlistpool,
				 (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								  /*pathL*/(Path_T) List_head(p),/*pathH*/pathpair->path5,
								  /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
								  nmismatches_filter_5,nmismatches_filter_3,
								  mincoverage_filter_5,mincoverage_filter_3,
								  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								  transcriptpool,hitlistpool,/*check_inner_p*/false,
								  /*copyLp (was false)*/true,/*copyHp*/true)
				 hitlistpool_trace(__FILE__,__LINE__));
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_3,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

	} else {
	  for (p = fusion_paths_5; p != NULL; p = List_next(p)) {
	    for (q = fusion_paths_3; q != NULL; q = List_next(q)) {
	      paths = Hitlist_push(paths,hitlistpool,
				   (void *) Pathpair_new_concordant(&unresolved_pathpairs,
								    /*pathL*/(Path_T) List_head(q),/*pathH*/(Path_T) List_head(p),
								    /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
								    nmismatches_filter_5,nmismatches_filter_3,
								    mincoverage_filter_5,mincoverage_filter_3,
								    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
								    transcriptpool,hitlistpool,/*check_inner_p*/false,
								    /*copyLp (was false)*/true,/*copyHp (was false)*/true)
				   hitlistpool_trace(__FILE__,__LINE__));
	    }
	  }
	  /* Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
	     listpool,pathpool,transcriptpool,hitlistpool); */
	  Path_gc(&fusion_paths_5,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
	  Path_gc(&fusion_paths_3,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }

    /* Use new pathpairarray that integrates original paths and fusions */
    FREE_OUT(pathpairarray);
    pathpairarray = (Pathpair_T *) List_to_array_out_n(&npaths,paths);
    Hitlistpool_free_list(&paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    debug8(printf("HAVE %d PATHS\n",npaths));
  }


  /* 4.  Find best sensedir for paths in each local region */
  /* Do not reward unnecessary splices, so do not eliminate alignments with splice prob of 0.0 */
  if (npaths > 1) {
    qsort(pathpairarray,npaths,sizeof(Pathpair_T),Pathpair_interval_cmp);

    k = 0;
    i = 0;
    while (i < npaths) {
      j = i + 1;
      while (j < npaths && Pathpair_overlap_p(pathpairarray[j],pathpairarray[i]) == true) {
	j++;
      }
	
      debug8(printf("Found an overlapping group of %d => choosing sense vs antisense\n",j - i));
      best_avg_junction_sense_prob = best_avg_junction_antisense_prob = 0.0;
      best_end_sense_prob = best_end_antisense_prob = 0.0;
      for (l = i; l < j; l++) {
	pathpair = pathpairarray[l];
	if (Pathpair_sensedir(pathpair) == SENSE_FORWARD) {
	  if ((prob = Pathpair_avg_junction_splice_prob(pathpair)) > best_avg_junction_sense_prob) {
	    best_avg_junction_sense_prob = prob;
	  }
	  if ((prob = Pathpair_end_splice_prob(pathpair)) > best_end_sense_prob) {
	    best_end_sense_prob = prob;
	  }
	} else {
	  if ((prob = Pathpair_avg_junction_splice_prob(pathpair)) > best_avg_junction_antisense_prob) {
	    best_avg_junction_antisense_prob = prob;
	  }
	  if ((prob = Pathpair_end_splice_prob(pathpair)) > best_end_antisense_prob) {
	    best_end_antisense_prob = prob;
	  }
	}
      }
	
      debug8(printf("best avg junction sense prob %f, best avg junction antisense prob %f\n",
		    best_avg_junction_sense_prob,best_avg_junction_antisense_prob));
      debug8(printf("best end sense prob %f, best end antisense prob %f\n",
		    best_end_sense_prob,best_end_antisense_prob));
      if (best_avg_junction_sense_prob > best_avg_junction_antisense_prob ||
	  (best_avg_junction_sense_prob == best_avg_junction_antisense_prob && best_end_sense_prob > best_end_antisense_prob)) {
	/* Keep only sense */
	for (l = i; l < j; l++) {
	  pathpair = pathpairarray[l];
	  if (Pathpair_total_splice_prob(pathpair) == 0.0) {
	    debug8(printf("(4 sensedir) Keeping null ")); debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;
	  } else if (Pathpair_sensedir(pathpair) == SENSE_FORWARD) {
	    debug8(printf("(4 sensedir) Keeping sense ")); debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;
	  } else {
	    debug8(printf("(4 sensedir) Eliminating antisense ")); debug8(Pathpair_print(pathpair));
	    Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

      } else if (best_avg_junction_antisense_prob > best_avg_junction_sense_prob ||
		 (best_avg_junction_antisense_prob == best_avg_junction_sense_prob && best_end_antisense_prob > best_end_sense_prob)) {
	/* Keep only antisense */
	for (l = i; l < j; l++) {
	  pathpair = pathpairarray[l];
	  if (Pathpair_total_splice_prob(pathpair) == 0.0) {
	    debug8(printf("(4 sensedir) Keeping null ")); debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;
	  } else if (Pathpair_sensedir(pathpair) == SENSE_ANTI) {
	    debug8(printf("(4 sensedir) Keeping antisense ")); debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;
	  } else {
	    debug8(printf("(4 sensedir) Eliminating sense ")); debug8(Pathpair_print(pathpair));
	    Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

      } else {
	/* Keep both */
	for (l = i; l < j; l++) {
	  pathpair = pathpairarray[l];
	  pathpairarray[k++] = pathpair;
	}
      }
	
      i = j;
    }
    npaths = k;
  }


  /* 5.  Find best paths in each local region */
  if (npaths > 1) {
    /* Array should already be sorted */
    k = 0;
    i = 0;
    while (i < npaths) {
      chop_qstart_L = Intlist_head(pathpairarray[i]->pathL->endpoints);
      chop_qend_L = Intlist_last_value(pathpairarray[i]->pathL->endpoints);
      chop_qstart_H = Intlist_head(pathpairarray[i]->pathH->endpoints);
      chop_qend_H = Intlist_last_value(pathpairarray[i]->pathH->endpoints);

      j = i + 1;
      while (j < npaths && Pathpair_overlap_p(pathpairarray[j],pathpairarray[i]) == true) {
	if ((qstart = Intlist_head(pathpairarray[j]->pathL->endpoints)) > chop_qstart_L) {
	  chop_qstart_L = qstart;
	}
	if ((qend = Intlist_last_value(pathpairarray[j]->pathL->endpoints)) < chop_qend_L) {
	  chop_qend_L = qend;
	}
	if ((qstart = Intlist_head(pathpairarray[j]->pathH->endpoints)) > chop_qstart_H) {
	  chop_qstart_H = qstart;
	}
	if ((qend = Intlist_last_value(pathpairarray[j]->pathH->endpoints)) < chop_qend_H) {
	  chop_qend_H = qend;
	}
	j++;
      }
	
      debug8(printf("Found an overlapping group of %d, chop region %d..%d and %d..%d => re-sorting by Pathpair_local_cmp\n",
		    j - i,chop_qstart_L,chop_qend_L,chop_qstart_H,chop_qend_H));

      for (l = i; l < j; l++) {
	pathpairarray[l]->pathL->chop_qstart = chop_qstart_L;
	pathpairarray[l]->pathL->chop_qend = chop_qend_L;
	pathpairarray[l]->pathH->chop_qstart = chop_qstart_H;
	pathpairarray[l]->pathH->chop_qend = chop_qend_H;
      }
      
      /* Keep the best ones in the overlapping group */
      qsort(&(pathpairarray[i]),j - i,sizeof(Pathpair_T),Pathpair_local_cmp);
      debug8(printf("(5 best) Keeping by local_cmp ")); debug8(Pathpair_print(pathpairarray[i]));
      pathpairarray[k++] = pathpairarray[i];
	
      for (l = i + 1; l < j; l++) {
	pathpair = pathpairarray[l];
	if (Pathpair_structure_ignore_sense_cmp(&pathpair,&(pathpairarray[i])) == 0) {
	  debug8(printf("(5 identical) Eliminating by identical structure ")); debug8(Pathpair_print(pathpair));
	  Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	  
	} else if (Pathpair_local_cmp(&pathpair,&(pathpairarray[i])) == 0) {
	  debug8(printf("(5 tie) Keeping by local_cmp ")); debug8(Pathpair_print(pathpair));
	  pathpairarray[k++] = pathpair;
	    
	} else {
	  debug8(printf("(5 worse) Eliminating by local_cmp ")); debug8(Pathpair_print(pathpair));
	  Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
	
      i = j;
    }
    npaths = k;
  }


  /* 6.  Find best solution globally */
  /* TODO: Can make this O(n) rather than O(n*log n) */
  if (npaths > 1) {
    chop_qstart_L = Intlist_head(pathpairarray[0]->pathL->endpoints);
    chop_qend_L = Intlist_last_value(pathpairarray[0]->pathL->endpoints);
    chop_qstart_H = Intlist_head(pathpairarray[0]->pathH->endpoints);
    chop_qend_H = Intlist_last_value(pathpairarray[0]->pathH->endpoints);

    for (i = 1; i < npaths; i++) {
      if ((qstart = Intlist_head(pathpairarray[i]->pathL->endpoints)) > chop_qstart_L) {
	chop_qstart_L = qstart;
      }
      if ((qend = Intlist_last_value(pathpairarray[i]->pathL->endpoints)) < chop_qend_L) {
	chop_qend_L = qend;
      }
      if ((qstart = Intlist_head(pathpairarray[i]->pathH->endpoints)) > chop_qstart_H) {
	chop_qstart_H = qstart;
      }
      if ((qend = Intlist_last_value(pathpairarray[i]->pathH->endpoints)) < chop_qend_H) {
	chop_qend_H = qend;
      }
    }

    for (i = 0; i < npaths; i++) {
      pathpairarray[i]->pathL->chop_qstart = chop_qstart_L;
      pathpairarray[i]->pathL->chop_qend = chop_qend_L;
      pathpairarray[i]->pathH->chop_qstart = chop_qstart_H;
      pathpairarray[i]->pathH->chop_qend = chop_qend_H;
    }

    qsort(pathpairarray,npaths,sizeof(T),Pathpair_global_cmp);
  }
  debug8(printf("Found the global best solution across %d paths with %d + %d nmatches\n",
		npaths,pathpairarray[0]->path5->nmatches,pathpairarray[0]->path3->nmatches));


  /* 6.  Check if we should be filtering the result */
  pathpair = pathpairarray[0];
  if (filterp == true &&
      (pathpair->path5->score_within_trims > nmismatches_filter_5 || pathpair->path3->score_within_trims > nmismatches_filter_3)) {
    debug8(printf("(6 filter) Best solution has too many %d + %d nmismatches, so returning early\n",
		  pathpair->path5->score_within_trims,pathpair->path3->score_within_trims));
    for (k = 0; k < npaths; k++) {
      debug8(printf("(6 filter) Eliminating ")); debug8(Pathpair_print(pathpairarray[k]));
      pathpair = pathpairarray[k];
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    return (Pathpair_T *) NULL;

  } else if (filterp == true && (Path_coverage(pathpair->path5) < mincoverage_filter_5 &&
				 Path_coverage(pathpair->path3) < mincoverage_filter_3)) {
    /* Was || on coverage */
    debug8(printf("Best solution has too little %d and %d coverage, so eliminating all\n",
		  Path_coverage(pathpair->path5),Path_coverage(pathpair->path3)));
    for (k = 0; k < npaths; k++) {
      debug8(printf("(6 filter) Eliminating ")); debug8(Pathpair_print(pathpairarray[k]));
      pathpair = pathpairarray[k];
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    return (Pathpair_T *) NULL;

  } else if (filterp == true && allow_soft_clips_p == false &&
	     (Intlist_head(pathpair->path5->endpoints) != 0 || Intlist_last_value(pathpair->path5->endpoints) != querylength5 ||
	      Intlist_head(pathpair->path3->endpoints) != 0 || Intlist_last_value(pathpair->path3->endpoints) != querylength3)) {
    debug8(printf("Best solution has soft clips, so eliminating all\n"));
    for (k = 0; k < npaths; k++) {
      debug8(printf("(6 filter) Eliminating ")); debug8(Pathpair_print(pathpairarray[k]));
      pathpair = pathpairarray[k];
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    return (Pathpair_T *) NULL;

  } else {
    /* Otherwise, keep all paths equivalent to the best path */
    debug8(printf("(6 best) Keeping best by global_cmp\n")); debug8(Pathpair_print(pathpairarray[0]));
    i = 1;		/* Skip the best path */
      
    for (k = 1; k < npaths; k++) {
      pathpair = pathpairarray[k];
      if (Pathpair_global_cmp(&pathpair,&(pathpairarray[0])) == 0) {
	debug8(printf("(6 tie) Keeping by global_cmp\n")); debug8(Pathpair_print(pathpair));

	pathpairarray[i++] = pathpair;
      } else {
	debug8(printf("(6 worse) Eliminating by global_cmp\n")); debug8(Pathpair_print(pathpair));
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
    }

    npaths = i;			
  }


  if (want_random_p && npaths > 1) {
    /* Randomize among best alignments */
    /* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
    randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) npaths));
    /* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
    debug8(printf("Swapping random path (%d out of %d) %p with best path %p\n",
		  randomi,npaths,pathpairarray[randomi],pathpairarray[0]));
    temp = pathpairarray[0];
    pathpairarray[0] = pathpairarray[randomi];
    pathpairarray[randomi] = temp;
  }


  /* Trim alignments at chromosomal bounds */
  for (i = 0; i < npaths; i++) {
    pathpair = pathpairarray[i];
    if (circularp[pathpair->pathL->chrnum] == false) {
#if 0
      /* Now done at the start of the procedure */
      Path_trim_qstart_trimbounds(pathpair->pathL,intlistpool,univcoordlistpool,listpool,pathpool,
				  /*trimbounds*/pathpair->pathL->chroffset);
#endif
    } else if (output_type == STD_OUTPUT || output_type == M8_OUTPUT) {
      /* If output type is alignment or m8, then don't want to split up the parts */
    } else if (pathpair->pathL->plusp == true) {
      Path_trim_circular(pathpair->pathL,query5_compress_fwd,intlistpool,univcoordlistpool,listpool);
    } else {
      Path_trim_circular(pathpair->pathL,query3_compress_rev,intlistpool,univcoordlistpool,listpool);
    }

    if (circularp[pathpair->pathH->chrnum] == false) {
#if 0
      Path_trim_qend_trimbounds(pathpair->pathH,intlistpool,univcoordlistpool,listpool,pathpool,
				/*trimbounds*/pathpair->pathH->chrhigh);
#endif
    } else if (output_type == STD_OUTPUT || output_type == M8_OUTPUT) {
      /* If output type is alignment or m8, then don't want to split up the parts */
    } else if (pathpair->pathH->plusp == true) {
      Path_trim_circular(pathpair->pathH,query3_compress_fwd,intlistpool,univcoordlistpool,listpool);
    } else {
      Path_trim_circular(pathpair->pathH,query5_compress_rev,intlistpool,univcoordlistpool,listpool);
    }
  }


  /* Compute mapq_loglik */
  if (npaths == 1) {
    pathpair = pathpairarray[0];

    Path_mark_alignment(pathpair->path5,query5_compress_fwd,queryuc_ptr_5,query5_compress_rev,queryrc5,
			pathpool);
    Path_mark_alignment(pathpair->path3,query3_compress_fwd,queryuc_ptr_3,query3_compress_rev,queryrc3,
			pathpool);

    pathpair->mapq_loglik = MAPQ_MAXIMUM_SCORE;
    pathpair->mapq_score = MAPQ_max_quality_score(quality_string_5,querylength5);
    if ((mapq_score = MAPQ_max_quality_score(quality_string_3,querylength3)) > pathpairarray[0]->mapq_score) {
      pathpair->mapq_score = mapq_score;
    }
    pathpair->absmq_score = MAPQ_MAXIMUM_SCORE;

    *first_absmq = pathpair->absmq_score;
    *second_absmq = 0;

  } else {
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      path5 = pathpair->path5;
      path3 = pathpair->path3;
      Path_mark_alignment(path5,query5_compress_fwd,queryuc_ptr_5,query5_compress_rev,queryrc5,
			  pathpool);
      Path_mark_alignment(path3,query3_compress_fwd,queryuc_ptr_3,query3_compress_rev,queryrc3,
			  pathpool);
      pathpair->mapq_loglik =
	MAPQ_loglik_string(path5->genomic_diff,quality_string_5,querylength5,path5->plusp);
      pathpair->mapq_loglik +=
	MAPQ_loglik_string(path3->genomic_diff,quality_string_3,querylength3,path3->plusp);
    }


    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (pathpairarray[i-1]->mapq_loglik < pathpairarray[i]->mapq_loglik) {
	pathpairarray[i-1]->mapq_loglik = pathpairarray[i]->mapq_loglik;
      }
    }
    maxlik = pathpairarray[0]->mapq_loglik;

    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      pathpairarray[i]->mapq_loglik -= maxlik;
    }

    /* Compute absolute mapq */
    for (i = 0; i < npaths; i++) {
      loglik = pathpairarray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      pathpairarray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = pathpairarray[0]->absmq_score;
    if (npaths == 1) {
      *second_absmq = 0;
    } else {
      *second_absmq = pathpairarray[1]->absmq_score;
    }

    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (pathpairarray[i]->mapq_loglik = fasterexp(pathpairarray[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < npaths; i++) {
      pathpairarray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < npaths; i++) {
      if ((qual = 1.0 - pathpairarray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	pathpairarray[i]->mapq_score = 40;
      } else {
	pathpairarray[i]->mapq_score = rint(-10.0 * log10(qual));
      }
    }
  }


  /* Compute transcripts */
  if (transcriptome != NULL) {
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      path5 = pathpair->path5;
      path3 = pathpair->path3;

      Transcript_velocity_paired(path5,path3);
      if ((tr_insertlength = Transcript_fragment_length(path5,path3,queryseq5,queryseq3)) > 0) {
	pathpair->insertlength = tr_insertlength;
      }
    }
  }


  /* Filter for chrsubset */
  /* Want to allow other alignments to be found before filtering */
  *npaths_primary = *npaths_altloc = 0;
  if (chrsubsetp == NULL) {
    k = 0;
    for (pathi = 0; pathi < npaths; pathi++) {
      pathpair = pathpairarray[pathi];
      if (altlocp[pathpair->path5->chrnum] == true ||
	  altlocp[pathpair->path3->chrnum] == true) {
	(*npaths_altloc) += 1;
      } else {
	(*npaths_primary) += 1;
      }
    }

    return pathpairarray;

  } else {
    k = 0;
    for (pathi = 0; pathi < npaths; pathi++) {
      pathpair = pathpairarray[pathi];
      if (chrsubsetp[pathpair->path5->chrnum] == false ||
	  chrsubsetp[pathpair->path3->chrnum] == false) {
	/* Do not save this pathpair */
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	
      } else {
	/* Save this pathpair.  Re-use existing array */
	if (altlocp[pathpair->path5->chrnum] == true ||
	    altlocp[pathpair->path3->chrnum] == true) {
	  (*npaths_altloc) += 1;
	} else {
	  (*npaths_primary) += 1;
	}
	pathpairarray[k++] = pathpair;
      }
    }

    if ((*npaths_primary) + (*npaths_altloc) == 0) {
      FREE_OUT(pathpairarray);
      return (Pathpair_T *) NULL;
    } else {
      return pathpairarray;
    }
  }
}


void
Pathpair_eval_setup (int expected_pairlength, int pairlength_deviation,
		     Transcriptome_T transcriptome_in,
		     bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		     Outputtype_T output_type_in,
		     bool splicingp_in, bool resolve_inner_p_in, bool want_random_p_in,
		     bool allow_soft_clips_p_in) {

  max_insertlength = expected_pairlength + 2*pairlength_deviation;
  transcriptome = transcriptome_in;

  circularp = circularp_in;
  chrsubsetp = chrsubsetp_in;
  altlocp = altlocp_in;

  output_type = output_type_in;
  splicingp = splicingp_in;
  resolve_inner_p = resolve_inner_p_in;
  want_random_p = want_random_p_in;
  allow_soft_clips_p = allow_soft_clips_p_in;

  return;
}

