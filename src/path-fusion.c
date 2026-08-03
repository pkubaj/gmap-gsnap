#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"

#include "splice.h"
#include "junction.h"
#include "transcript.h"
#include "transcript-remap.h"
#include "spliceends.h"
#include "genomebits_trim.h"

#include "univdiagdef.h"
#include "extension-search.h"	/* For Elt_T */
#include "path-eval.h"


/* qstart and qend are represented in path->endpoints and path->fusion_endpoints.
   For anchor and main, qstart is the lowest genomic position, and qend is the highest genomic position.
   For fusion, qstart and qend represent the splice_qpos coordinate that matches main.

   Path endpoints shown below are 0, s for splice_querypos, and L for querylength

   Paired-end reads:

   Case 1:  fusion5+ [0,s)   (+D) .. (+A) main5+ (L-s,L) ; anchor3+  [fusion_querystart_plus]
   Case 2:  fusion5- (L-s,L] (+D) .. (+A) main5+ (L-s,L) ; anchor3+  [fusion_querystart_plus]
   Case 3:                                      anchor5+ ; (0,s) main3+ (+D) .. (+A) fusion3+ (L-s,L]  [fusion_queryend_plus]
   Case 4:                                      anchor5+ ; (0,s) main3+ (+D) .. (+A) fusion3- [0,s)    [fusion_queryend_plus]

   Case 4':                                     anchor3- ; (0,s) main5- (-A) .. (-D) fusion5+ [0,s)    [fusion_querystart_minus]
   Case 3':                                     anchor3- ; (0,s) main5- (-A) .. (-D) fusion5- (L-s,L]  [fusion_querystart_minus]
   Case 2': fusion3+ (L-s,L] (-A) .. (-D) main3- (L-s,L) ; anchor5-  [fusion_queryend_minus]
   Case 1': fusion3- [0,s)   (-A) .. (-D) main3- (L-s,L) ; anchor5-  [fusion_queryend_minus]


   Single-end reads:

   Case 1:  fusion+ (0,s)   (+D) .. (+A) main+ (L-s,L)    [fusion_querystart]
   Case 2:  fusion- (L-s,L) (+D) .. (+A) main+ (L-s,L)    [fusion_querystart]
   Case 3:  (0,s) main+     (+D) .. (+A) fusion+ (L-s,L)  [fusion_queryend]
   Case 4:  (0,s) main+     (+D) .. (+A) fusion- (0,s)    [fusion_queryend]

   Case 4': (0,s) main-     (-A) .. (-D) fusion+ (0,s)    [fusion_querystart]
   Case 3': (0,s) main-     (-A) .. (-D) fusion- (L-s,L)  [fusion_querystart]
   Case 2': fusion+ (L-s,L) (-A) .. (-D) main- (L-s,L)    [fusion_queryend]
   Case 1': fusion- (0,s)   (-A) .. (-D) main- (L-s,L)    [fusion_queryend]

*/


/* We are disallowing translocations involving circular chromosomes,
   but allowing splicing within circular chromosomes, which might be
   found through fusion methods */
static bool *circularp;

static EF64_T chromosome_ef64;
static Univcoord_T genomelength;
static Chrpos_T shortsplicedist;

static Transcriptome_T transcriptome;

static Genomebits_T genomebits;


#define MIN_OUTER_END 20

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Path creation */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Path_fusion_copy */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* For inner fusions, we do need to bound univdiagonals at chromosomal
   boundaries, or else newpath will not have the same chrnum as the
   anchor */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))

#define T Path_T


static Splicetype_T
flip_splicetype (Splicetype_T splicetype) {
  if (splicetype == DONOR) {
    return ACCEPTOR;
  } else if (splicetype == ACCEPTOR) {
    return DONOR;
  } else if (splicetype == ANTIDONOR) {
    return ANTIACCEPTOR;
  } else if (splicetype == ANTIACCEPTOR) {
    return ANTIDONOR;
  } else {
    return splicetype;
  }
}
    

static bool
querystart_fusion_p (Chrpos_T *splice_distance, T old, T fusion, int splice_querypos) {
  Univcoord_T univdiagonalH, univdiagonalL;

  if (old->chrnum != fusion->chrnum) {
    return true;

  } else if (old->plusp != fusion->plusp) {
    return true;

  } else if (circularp[old->chrnum] == true) {
    return false;

  } else if (old->plusp == true) {
    /* fusion is L; old is H */
    if ((univdiagonalH = Univcoordlist_head(old->univdiagonals)) <
	(univdiagonalL = Univcoordlist_last_value(fusion->univdiagonals))) {
      debug2(printf("(1) Checking univdiagonalL %u < univdiagonalH %u => true\n",univdiagonalL,univdiagonalH));
      return true;
    } else if ((*splice_distance = (Chrpos_T) (univdiagonalH - univdiagonalL)) <= shortsplicedist) {
      debug2(printf("(2) Checking univdiagonalH %u - univdiagonalL %u vs shortsplicedist => false\n",univdiagonalH,univdiagonalL));
      return false;
    } else if (Transcript_remap_matchp((Chrpos_T) (univdiagonalL - old->chroffset - old->querylength + splice_querypos),
				       (Chrpos_T) (univdiagonalH - old->chroffset - old->querylength + splice_querypos),
				       old->chrnum) == true) {
      return false;
    } else {
      return true;
    }

  } else {
    /* fusion is H; old is L */
    if ((univdiagonalH = Univcoordlist_head(fusion->univdiagonals)) <
	(univdiagonalL = Univcoordlist_last_value(old->univdiagonals))) {
      debug2(printf("(3) Checking univdiagonalH %u < univdiagonalL %u => true\n",univdiagonalH,univdiagonalL));
      return true;
    } else if ((*splice_distance = (Chrpos_T) (univdiagonalH - univdiagonalL)) <= shortsplicedist) {
      debug2(printf("(4) Checking univdiagonalH %u - univdiagonalL %u vs shortsplicedist => false\n",univdiagonalH,univdiagonalL));
      return false;
    } else if (Transcript_remap_matchp((Chrpos_T) (univdiagonalL - old->chroffset - splice_querypos),
				       (Chrpos_T) (univdiagonalH - old->chroffset - splice_querypos),
				       old->chrnum) == true) {
      return false;
    } else {
      return true;
    }
  }
}


static bool
queryend_fusion_p (Chrpos_T *splice_distance, T old, T fusion, int splice_querypos) {
  Univcoord_T univdiagonalH, univdiagonalL;

  if (old->chrnum != fusion->chrnum) {
    return true;

  } else if (old->plusp != fusion->plusp) {
    return true;

  } else if (circularp[old->chrnum] == true) {
    return false;

  } else if (old->plusp == true) {
    /* old is L; fusion is H */
    if ((univdiagonalH = Univcoordlist_head(fusion->univdiagonals)) <
	(univdiagonalL = Univcoordlist_last_value(old->univdiagonals))) {
      debug2(printf("(5) Checking univdiagonalH %u < univdiagonalL %u => true\n",univdiagonalH,univdiagonalL));
      return true;
    } else if ((*splice_distance = (Chrpos_T) (univdiagonalH - univdiagonalL)) <= shortsplicedist) {
      debug2(printf("(6) Checking univdiagonalH %u - univdiagonalL %u vs shortsplicedist => false\n",univdiagonalH,univdiagonalL));
      return false;
    } else if (Transcript_remap_matchp((Chrpos_T) (univdiagonalL - old->chroffset - old->querylength + splice_querypos),
				       (Chrpos_T) (univdiagonalH - old->chroffset - old->querylength + splice_querypos),
				       old->chrnum) == true) {
      return false;
    } else {
      return true;
    }

  } else {
    /* old is H; fusion is L */
    if ((univdiagonalH = Univcoordlist_head(old->univdiagonals)) < (univdiagonalL = Univcoordlist_last_value(fusion->univdiagonals))) {
      debug2(printf("(7) Checking univdiagonalH %u < univdiagonalL %u => true\n",univdiagonalH,univdiagonalL));
      return true;
    } else if ((*splice_distance = (Chrpos_T) (univdiagonalH - univdiagonalL)) <= shortsplicedist) {
      debug2(printf("(8) Checking univdiagonalH %u - univdiagonalL %u vs shortsplicedist => false\n",univdiagonalH,univdiagonalL));
      return false;
    } else if (Transcript_remap_matchp((Chrpos_T) (univdiagonalL - old->chroffset - splice_querypos),
				       (Chrpos_T) (univdiagonalH - old->chroffset - splice_querypos),
				       old->chrnum) == true) {
      return false;
    } else {
      return true;
    }
  }
}



static T
Path_fusion_copy_querystart (T old, T fusion, int splice_querypos, int querylength,
			     int nmismatches_main, int ref_nmismatches_main,
			     int nmismatches_fusion, int ref_nmismatches_fusion,
			     char donor1, char donor2, char acceptor1, char acceptor2,
			     double donor_prob, double acceptor_prob, Shortread_T queryseq,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			     bool fusion_created_p) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new;
  Chrpos_T splice_distance;
  Intlist_T qend_endpoints, qend_nmismatches, qend_ref_nmismatches;
  int ignore_int;
  int desired_genestrand;
  bool extend_qstart_p, extend_qend_p;

  debug2(printf("\nEntered Path_fusion_copy_querystart with nmismatches main %d and fusion %d, splice querypos %d, fusion created %d\n",
		nmismatches_main,nmismatches_fusion,splice_querypos,fusion_created_p));
  debug2(Path_print(old));
  debug2(Path_print(fusion));
  
  if (querystart_fusion_p(&splice_distance,old,fusion,splice_querypos) == true) {
    new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
    debug0(printf("%d: Creating path %p by Path_fusion_copy_querystart from mainpath %p and fusion %p\n",
		  ++call_i,new,old,fusion));

    new->query_compress_fwd = old->query_compress_fwd;
    new->query_compress_rev = old->query_compress_rev;

    new->nmatches = -1;
    new->ref_nmatches = -1;

    new->junction_splice_prob = old->junction_splice_prob;
    new->total_splice_prob = old->total_splice_prob;
    new->found_score = old->found_score;
    new->score_within_trims = old->score_within_trims;

    new->genomic_diff = (char *) NULL;

    new->plusp = old->plusp;
    new->genestrand = old->genestrand;
    new->sensedir = old->sensedir;
    new->querylength = old->querylength;

    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;
    new->chrhigh = old->chrhigh;

    if (new->plusp == true) {
      /* Revise the qstart of old */
      debug2(printf("Revising qstart of mainpath\n"));
      new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
      new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
      new->ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);

      Intlist_head_set(new->endpoints,splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);

      new->qstart_alts = (Altsplice_T) NULL;
      new->splice5p = false;
      new->splicetype5 = NO_SPLICE;
      new->ambig_prob_5 = 0.0;

      new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
      new->splice3p = old->splice3p;
      new->splicetype3 = old->splicetype3;
      new->ambig_prob_3 = old->ambig_prob_3;

    } else {
      /* Revise the qend of old */
      debug2(printf("Revising qend of mainpath\n"));
      new->endpoints = Intlist_reverse(Intlistpool_copy(old->endpoints,intlistpool));
      new->nmismatches = Intlist_reverse(Intlistpool_copy(old->nmismatches,intlistpool));
      new->ref_nmismatches = Intlist_reverse(Intlistpool_copy(old->ref_nmismatches,intlistpool));

      Intlist_head_set(new->endpoints,querylength - splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);

      new->endpoints = Intlist_reverse(new->endpoints);
      new->nmismatches = Intlist_reverse(new->nmismatches);
      new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

      new->qstart_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
      new->splice5p = old->splice5p;
      new->splicetype5 = old->splicetype5;
      new->ambig_prob_5 = old->ambig_prob_5;

      new->qend_alts = (Altsplice_T) NULL;
      new->splice3p = false;
      new->splicetype3 = NO_SPLICE;
      new->ambig_prob_3 = 0.0;
    }

    new->main_univdiagonal = old->main_univdiagonal;
    new->univdiagonals = Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool);
    new->junctions = Junction_copy_list(old->junctions,listpool,pathpool);

    if (old->circular_endpoints == NULL) {
      new->circular_endpoints = (Intlist_T) NULL;
    } else {
      new->circular_high_p = old->circular_high_p;
      new->circular_endpoints = Intlistpool_copy(old->circular_endpoints,intlistpool);
      new->circular_univdiagonals = Univcoordlistpool_copy(old->circular_univdiagonals,univcoordlistpool);
      new->circular_nmismatches = Intlistpool_copy(old->circular_nmismatches,intlistpool);
      new->circular_ref_nmismatches = Intlistpool_copy(old->circular_ref_nmismatches,intlistpool);
      new->circular_junctions = Junction_copy_list(old->circular_junctions,listpool,pathpool);
    }

    new->fusion_querystart_junction = Junction_new_chimera(donor1,donor2,acceptor1,acceptor2,
							   donor_prob,acceptor_prob,pathpool);
    new->fusion_queryend_junction = (Junction_T) NULL;

    new->fusion_chrnum = fusion->chrnum;
    new->fusion_chroffset = fusion->chroffset;
    new->fusion_chrhigh = fusion->chrhigh;
    new->fusion_plusp = fusion->plusp;

    if (fusion_created_p == true) {
      /* Already computed by Path_create */
      new->fusion_endpoints = Intlistpool_copy(fusion->endpoints,intlistpool);
      new->fusion_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
      new->fusion_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

    } else if (fusion->plusp == true) {
      /* Revise the qend of fusion */
      debug2(printf("Revising the qend of fusion\n"));
      new->fusion_endpoints = Intlist_reverse(Intlistpool_copy(fusion->endpoints,intlistpool));
      new->fusion_nmismatches = Intlist_reverse(Intlistpool_copy(fusion->nmismatches,intlistpool));
      new->fusion_ref_nmismatches = Intlist_reverse(Intlistpool_copy(fusion->ref_nmismatches,intlistpool));

      Intlist_head_set(new->fusion_endpoints,splice_querypos);
      Intlist_head_set(new->fusion_nmismatches,nmismatches_fusion);
      Intlist_head_set(new->fusion_ref_nmismatches,ref_nmismatches_fusion);

      new->fusion_endpoints = Intlist_reverse(new->fusion_endpoints);
      new->fusion_nmismatches = Intlist_reverse(new->fusion_nmismatches);
      new->fusion_ref_nmismatches = Intlist_reverse(new->fusion_ref_nmismatches);

    } else {
      /* Revise the qstart of fusion */
      debug2(printf("Revising the qstart of fusion\n"));
      new->fusion_endpoints = Intlistpool_copy(fusion->endpoints,intlistpool);
      new->fusion_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
      new->fusion_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

      Intlist_head_set(new->fusion_endpoints,querylength - splice_querypos);
      Intlist_head_set(new->fusion_nmismatches,nmismatches_fusion);
      Intlist_head_set(new->fusion_ref_nmismatches,ref_nmismatches_fusion);
    }

    if (fusion->plusp == true) {
      new->fusion_alts = Altsplice_copy(fusion->qstart_alts,pathpool,vectorpool);
      new->fusion_splicep = fusion->splice5p;
      new->fusion_splicetype = fusion->splicetype5;
      new->fusion_ambig_prob = fusion->ambig_prob_5;

    } else {
      new->fusion_alts = Altsplice_copy(fusion->qend_alts,pathpool,vectorpool);
      new->fusion_splicep = fusion->splice3p;
      new->fusion_splicetype = flip_splicetype(fusion->splicetype3);
      new->fusion_ambig_prob = fusion->ambig_prob_3;
    }

    new->fusion_univdiagonals = Univcoordlistpool_copy(fusion->univdiagonals,univcoordlistpool);
    new->fusion_junctions = Junction_copy_list(fusion->junctions,listpool,pathpool);

    /* Recompute transcripts based on new endpoints */
    if (transcriptome == NULL) {
      new->transcripts = (List_T) NULL;
      new->invalid_transcripts = (List_T) NULL;
      new->fusion_transcripts = (List_T) NULL;
      new->fusion_invalid_transcripts = (List_T) NULL;
    } else {
      if (new->plusp == true) {
	extend_qstart_p = false;
	extend_qend_p = true;
	desired_genestrand = (new->sensedir == SENSE_FORWARD) ? +1 : -1;
      } else {
	extend_qstart_p = true;
	extend_qend_p = false;
	desired_genestrand = (new->sensedir == SENSE_ANTI) ? +1 : -1;
      }
      Transcript_remap_all(&new->transcripts,&new->invalid_transcripts,
			   new->endpoints,new->univdiagonals,new->junctions,
			   queryseq,new->querylength,new->plusp,
			   new->chrnum,new->chroffset,new->chrhigh,
			   uintlistpool,listpool,transcriptpool,desired_genestrand,
			   extend_qstart_p,extend_qend_p,/*repairp*/false);

      if (fusion->transcripts != NULL || fusion->invalid_transcripts != NULL) {
	new->fusion_transcripts = Transcript_copy_list(fusion->transcripts,transcriptpool,listpool);
	new->fusion_invalid_transcripts = Transcript_copy_list(fusion->invalid_transcripts,transcriptpool,listpool);
      } else {
	if (new->fusion_plusp == true) {
	  extend_qstart_p = true;
	  extend_qend_p = false;
	  desired_genestrand = (new->sensedir == SENSE_FORWARD) ? +1 : -1;
	} else {
	  extend_qstart_p = false;
	  extend_qend_p = true;
	  desired_genestrand = (new->sensedir == SENSE_ANTI) ? +1 : -1;
	}
	Transcript_remap_all(&new->fusion_transcripts,&new->fusion_invalid_transcripts,
			     new->fusion_endpoints,new->fusion_univdiagonals,new->fusion_junctions,
			     queryseq,new->querylength,new->fusion_plusp,
			     new->fusion_chrnum,new->fusion_chroffset,new->fusion_chrhigh,
			     uintlistpool,listpool,transcriptpool,desired_genestrand,
			     extend_qstart_p,extend_qend_p,/*repairp*/false);
      }
    }

    new->childp = old->childp;
    new->extendedp = false;

    new->method = old->method;
    new->transcriptome_method_p = false;

    /* Performed by caller */
    /* Path_eval_nmatches(&(*found_score),new,query_compress_fwd,query_compress_rev); */

    debug2(printf("Result of Path_fusion_copy_querystart\n"));
    debug2(Path_print(new));

    return new;

  } else if (old->plusp == true) {
    /* plus */
    debug2(printf("Standard (non-fusion) case 1: fusion -> old, splice querypos %d, splice distance %u\n",
		  splice_querypos,splice_distance));
    new = Path_copy(fusion,intlistpool,univcoordlistpool,listpool,
		    pathpool,vectorpool,transcriptpool,hitlistpool);
    new->nmatches = -1;
    new->ref_nmatches = -1;
    
    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

    if (splice_distance == 0) {
      /* Merge pieces instead */
      new->endpoints = Intlistpool_pop(new->endpoints,intlistpool,&ignore_int
				       intlistpool_trace(__FILE__,__LINE__));
      new->nmismatches = Intlistpool_pop(new->nmismatches,intlistpool,&ignore_int
					 intlistpool_trace(__FILE__,__LINE__));
      new->ref_nmismatches = Intlistpool_pop(new->ref_nmismatches,intlistpool,&ignore_int
					     intlistpool_trace(__FILE__,__LINE__));
    } else {
      Intlist_head_set(new->endpoints,splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_fusion);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_fusion);
    }

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);


    qend_endpoints = Intlistpool_copy(Intlist_next(old->endpoints),intlistpool); /* Skip first elt in common with main */
    qend_nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
    qend_ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);

    if (splice_distance == 0) {
      Intlist_head_set(qend_nmismatches,-1);
      Intlist_head_set(qend_ref_nmismatches,-1);
    } else {
      /* Intlist_head_set(qend_endpoints,splice_querypos); -- first elt skipped */
      Intlist_head_set(qend_nmismatches,nmismatches_main);
      Intlist_head_set(qend_ref_nmismatches,ref_nmismatches_main);
    }

    new->endpoints = Intlist_append(new->endpoints,qend_endpoints);
    new->nmismatches = Intlist_append(new->nmismatches,qend_nmismatches);
    new->ref_nmismatches = Intlist_append(new->ref_nmismatches,qend_ref_nmismatches);

    /* new->qstart_alts should have be a copy of fusion->qstart_alts */
    if (new->qend_alts != NULL) {
      Altsplice_free(&new->qend_alts,pathpool);
    }
    new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);

    new->splice3p = old->splice3p;
    new->splicetype3 = old->splicetype3;
    new->ambig_prob_3 = old->ambig_prob_3;

    new->main_univdiagonal = old->main_univdiagonal;
    if (splice_distance == 0) {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(Univcoordlist_next(old->univdiagonals),univcoordlistpool));
      new->junctions = List_append(new->junctions,Junction_copy_list(old->junctions,listpool,pathpool));
    } else {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool));

      new->junctions = List_reverse(new->junctions);
      new->junctions = Listpool_push(new->junctions,listpool,
				     (void *) Junction_new_splice(splice_distance,new->sensedir,
								  donor_prob,acceptor_prob,pathpool)
				     listpool_trace(__FILE__,__LINE__));
      new->junctions = List_reverse(new->junctions);
      new->junctions = List_append(new->junctions,Junction_copy_list(old->junctions,listpool,pathpool));
    }

    debug2(printf("Standard case 1:\n"));
    debug2(Path_print(new));
    return new;

  } else {
    /* minus */
    debug2(printf("Standard (non-fusion) case 2: old -> fusion, splice querypos %d, splice distance %u\n",
		  splice_querypos,splice_distance));
    new = Path_copy(old,intlistpool,univcoordlistpool,listpool,
		    pathpool,vectorpool,transcriptpool,hitlistpool);
    new->nmatches = -1;
    new->ref_nmatches = -1;

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

    if (splice_distance == 0) {
      /* Merge pieces instead */
      new->endpoints = Intlistpool_pop(new->endpoints,intlistpool,&ignore_int
				       intlistpool_trace(__FILE__,__LINE__));
      new->nmismatches = Intlistpool_pop(new->nmismatches,intlistpool,&ignore_int
					 intlistpool_trace(__FILE__,__LINE__));
      new->ref_nmismatches = Intlistpool_pop(new->ref_nmismatches,intlistpool,&ignore_int
					     intlistpool_trace(__FILE__,__LINE__));
    } else {
      Intlist_head_set(new->endpoints,querylength - splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);
    }

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);


    qend_endpoints = Intlistpool_copy(Intlist_next(fusion->endpoints),intlistpool); /* Skip first elt in common with main */
    qend_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
    qend_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

    if (splice_distance == 0) {
      Intlist_head_set(qend_nmismatches,-1);
      Intlist_head_set(qend_ref_nmismatches,-1);
    } else {
      /* Intlist_head_set(qend_endpoints,querylength - splice_querypos); -- first elt skipped */
      Intlist_head_set(qend_nmismatches,nmismatches_fusion);
      Intlist_head_set(qend_ref_nmismatches,ref_nmismatches_fusion);
    }

    new->endpoints = Intlist_append(new->endpoints,qend_endpoints);
    new->nmismatches = Intlist_append(new->nmismatches,qend_nmismatches);
    new->ref_nmismatches = Intlist_append(new->ref_nmismatches,qend_ref_nmismatches);

    /* new->qstart_alts should be a copy of old->qstart_alts */
    if (new->qend_alts != NULL) {
      Altsplice_free(&new->qend_alts,pathpool);
    }
    new->qend_alts = Altsplice_copy(fusion->qend_alts,pathpool,vectorpool);

    new->splice3p = fusion->splice3p;
    new->splicetype3 = fusion->splicetype3;
    new->ambig_prob_3 = fusion->ambig_prob_3;

    new->main_univdiagonal = old->main_univdiagonal;
    if (splice_distance == 0) {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(Univcoordlist_next(fusion->univdiagonals),univcoordlistpool));
      new->junctions = List_append(new->junctions,Junction_copy_list(fusion->junctions,listpool,pathpool));
    } else {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(fusion->univdiagonals,univcoordlistpool));

      new->junctions = List_reverse(new->junctions);
      new->junctions = Listpool_push(new->junctions,listpool,
				     (void *) Junction_new_splice(splice_distance,new->sensedir,
								  donor_prob,acceptor_prob,pathpool)
				     listpool_trace(__FILE__,__LINE__));
      new->junctions = List_reverse(new->junctions);
      new->junctions = List_append(new->junctions,Junction_copy_list(fusion->junctions,listpool,pathpool));
    }

    debug2(printf("Standard case 2:\n"));
    debug2(Path_print(new));
    return new;
  }
}


static T
Path_fusion_copy_queryend (T old, T fusion, int splice_querypos, int querylength,
			   int nmismatches_main, int ref_nmismatches_main,
			   int nmismatches_fusion, int ref_nmismatches_fusion,
			   char donor1, char donor2, char acceptor1, char acceptor2,
			   double donor_prob, double acceptor_prob, Shortread_T queryseq,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			   bool fusion_created_p) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new;
  Chrpos_T splice_distance;
  Intlist_T qend_endpoints, qend_nmismatches, qend_ref_nmismatches;
  int ignore_int;
  int desired_genestrand;
  bool extend_qstart_p, extend_qend_p;


  debug2(printf("\nEntered Path_fusion_copy_queryend with nmismatches main %d and fusion %d, splice querypos %d, fusion created %d\n",
		nmismatches_main,nmismatches_fusion,splice_querypos,fusion_created_p));
  debug2(Path_print(old));
  debug2(Path_print(fusion));

  if (queryend_fusion_p(&splice_distance,old,fusion,splice_querypos) == true) {
    new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
    debug0(printf("%d: Creating path %p by Path_fusion_copy_queryend from mainpath %p and fusion %p\n",
		  ++call_i,new,old,fusion));

    new->query_compress_fwd = old->query_compress_fwd;
    new->query_compress_rev = old->query_compress_rev;

    new->nmatches = -1;
    new->ref_nmatches = -1;

    new->junction_splice_prob = old->junction_splice_prob;
    new->total_splice_prob = old->total_splice_prob;
    new->found_score = old->found_score;
    new->score_within_trims = old->score_within_trims;

    new->genomic_diff = (char *) NULL;

    new->plusp = old->plusp;
    new->genestrand = old->genestrand;
    new->sensedir = old->sensedir;
    new->querylength = old->querylength;

    new->chrnum = old->chrnum;
    new->chroffset = old->chroffset;
    new->chrhigh = old->chrhigh;

    if (new->plusp == true) {
      /* Revise the qend of old */
      debug2(printf("Revising qend of mainpath\n"));
      new->endpoints = Intlist_reverse(Intlistpool_copy(old->endpoints,intlistpool));
      new->nmismatches = Intlist_reverse(Intlistpool_copy(old->nmismatches,intlistpool));
      new->ref_nmismatches = Intlist_reverse(Intlistpool_copy(old->ref_nmismatches,intlistpool));

      Intlist_head_set(new->endpoints,splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);

      new->endpoints = Intlist_reverse(new->endpoints);
      new->nmismatches = Intlist_reverse(new->nmismatches);
      new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

      new->qstart_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
      new->splice5p = old->splice5p;
      new->splicetype5 = old->splicetype5;
      new->ambig_prob_5 = old->ambig_prob_5;

      new->qend_alts = (Altsplice_T) NULL;
      new->splice3p = false;
      new->splicetype3 = NO_SPLICE;
      new->ambig_prob_3 = 0.0;

    } else {
      /* Revise the qstart of old */
      debug2(printf("Revising qstart of mainpath\n"));
      new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
      new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
      new->ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);

      Intlist_head_set(new->endpoints,querylength - splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);

      new->qstart_alts = (Altsplice_T) NULL;
      new->splice5p = false;
      new->splicetype5 = NO_SPLICE;
      new->ambig_prob_5 = 0.0;

      new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
      new->splice3p = old->splice3p;
      new->splicetype3 = old->splicetype3;
      new->ambig_prob_3 = old->ambig_prob_3;
    }

    new->main_univdiagonal = old->main_univdiagonal;
    new->univdiagonals = Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool);
    new->junctions = Junction_copy_list(old->junctions,listpool,pathpool);

    if (old->circular_endpoints == NULL) {
      new->circular_endpoints = (Intlist_T) NULL;
    } else {
      new->circular_high_p = old->circular_high_p;
      new->circular_endpoints = Intlistpool_copy(old->circular_endpoints,intlistpool);
      new->circular_univdiagonals = Univcoordlistpool_copy(old->circular_univdiagonals,univcoordlistpool);
      new->circular_nmismatches = Intlistpool_copy(old->circular_nmismatches,intlistpool);
      new->circular_ref_nmismatches = Intlistpool_copy(old->circular_ref_nmismatches,intlistpool);
      new->circular_junctions = Junction_copy_list(old->circular_junctions,listpool,pathpool);
    }

    new->fusion_querystart_junction = (Junction_T) NULL;
    new->fusion_queryend_junction = Junction_new_chimera(donor1,donor2,acceptor1,acceptor2,
							 donor_prob,acceptor_prob,pathpool);

    new->fusion_chrnum = fusion->chrnum;
    new->fusion_chroffset = fusion->chroffset;
    new->fusion_chrhigh = fusion->chrhigh;
    new->fusion_plusp = fusion->plusp;

    if (fusion_created_p == true) {
      /* Already computed by Path_create */
      new->fusion_endpoints = Intlistpool_copy(fusion->endpoints,intlistpool);
      new->fusion_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
      new->fusion_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

    } else if (fusion->plusp == true) {
      /* Revise the qstart of fusion */
      debug2(printf("Revising the qstart of fusion\n"));
      new->fusion_endpoints = Intlistpool_copy(fusion->endpoints,intlistpool);
      new->fusion_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
      new->fusion_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

      Intlist_head_set(new->fusion_endpoints,splice_querypos);
      Intlist_head_set(new->fusion_nmismatches,nmismatches_fusion);
      Intlist_head_set(new->fusion_ref_nmismatches,ref_nmismatches_fusion);

    } else {
      /* Revise the qend of fusion */
      debug2(printf("Revising the qend of fusion\n"));
      new->fusion_endpoints = Intlist_reverse(Intlistpool_copy(fusion->endpoints,intlistpool));
      new->fusion_nmismatches = Intlist_reverse(Intlistpool_copy(fusion->nmismatches,intlistpool));
      new->fusion_ref_nmismatches = Intlist_reverse(Intlistpool_copy(fusion->ref_nmismatches,intlistpool));

      Intlist_head_set(new->fusion_endpoints,querylength - splice_querypos);
      Intlist_head_set(new->fusion_nmismatches,nmismatches_fusion);
      Intlist_head_set(new->fusion_ref_nmismatches,ref_nmismatches_fusion);

      new->fusion_endpoints = Intlist_reverse(new->fusion_endpoints);
      new->fusion_nmismatches = Intlist_reverse(new->fusion_nmismatches);
      new->fusion_ref_nmismatches = Intlist_reverse(new->fusion_ref_nmismatches);
    }

    if (fusion->plusp == true) {
      new->fusion_alts = Altsplice_copy(fusion->qend_alts,pathpool,vectorpool);
      new->fusion_splicep = fusion->splice3p;
      new->fusion_splicetype = fusion->splicetype3;
      new->fusion_ambig_prob = fusion->ambig_prob_3;
    } else {
      new->fusion_alts = Altsplice_copy(fusion->qstart_alts,pathpool,vectorpool);
      new->fusion_splicep = fusion->splice5p;
      new->fusion_splicetype = flip_splicetype(fusion->splicetype5);
      new->fusion_ambig_prob = fusion->ambig_prob_5;
    }

    new->fusion_univdiagonals = Univcoordlistpool_copy(fusion->univdiagonals,univcoordlistpool);
    new->fusion_junctions = Junction_copy_list(fusion->junctions,listpool,pathpool);

    /* Recompute transcripts based on new endpoints */
    if (transcriptome == NULL) {
      new->transcripts = (List_T) NULL;
      new->invalid_transcripts = (List_T) NULL;
      new->fusion_transcripts = (List_T) NULL;
      new->fusion_invalid_transcripts = (List_T) NULL;
    } else {
      if (new->plusp == true) {
	extend_qstart_p = true;
	extend_qend_p = false;
	desired_genestrand = (new->sensedir == SENSE_FORWARD) ? +1 : -1;
      } else {
	extend_qstart_p = false;
	extend_qend_p = true;
	desired_genestrand = (new->sensedir == SENSE_ANTI) ? +1 : -1;
      }
      Transcript_remap_all(&new->transcripts,&new->invalid_transcripts,
			   new->endpoints,new->univdiagonals,new->junctions,
			   queryseq,new->querylength,new->plusp,
			   new->chrnum,new->chroffset,new->chrhigh,
			   uintlistpool,listpool,transcriptpool,desired_genestrand,
			   extend_qstart_p,extend_qend_p,/*repairp*/false);
    
      if (fusion->transcripts != NULL || fusion->invalid_transcripts != NULL) {
	new->fusion_transcripts = Transcript_copy_list(fusion->transcripts,transcriptpool,listpool);
	new->fusion_invalid_transcripts = Transcript_copy_list(fusion->invalid_transcripts,transcriptpool,listpool);
      } else {
	if (new->fusion_plusp == true) {
	  extend_qstart_p = false;
	  extend_qend_p = true;
	  desired_genestrand = (new->sensedir == SENSE_FORWARD) ? +1 : -1;
	} else {
	  extend_qstart_p = true;
	  extend_qend_p = false;
	  desired_genestrand = (new->sensedir == SENSE_ANTI) ? +1 : -1;
	}
	Transcript_remap_all(&new->fusion_transcripts,&new->fusion_invalid_transcripts,
			     new->fusion_endpoints,new->fusion_univdiagonals,new->fusion_junctions,
			     queryseq,new->querylength,new->fusion_plusp,
			     new->fusion_chrnum,new->fusion_chroffset,new->fusion_chrhigh,
			     uintlistpool,listpool,transcriptpool,desired_genestrand,
			     extend_qstart_p,extend_qend_p,/*repairp*/false);
      }
    }

    new->childp = old->childp;
    new->extendedp = false;

    new->method = old->method;
    new->transcriptome_method_p = false;

    /* Performed by caller */
    /* Path_eval_nmatches(&(*found_score),new,query_compress_fwd,query_compress_rev); */

    debug2(printf("Result of Path_fusion_copy_queryend\n"));
    debug2(Path_print(new));

    return new;

  } else if (old->plusp == true) {
    /* plus */
    debug2(printf("Standard (non-fusion) case 3: old -> fusion, splice querypos %d, splice distance %u\n",
		  splice_querypos,splice_distance));
    new = Path_copy(old,intlistpool,univcoordlistpool,listpool,
		    pathpool,vectorpool,transcriptpool,hitlistpool);
    new->nmatches = -1;
    new->ref_nmatches = -1;

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

    if (splice_distance == 0) {
      /* Merge pieces instead */
      new->endpoints = Intlistpool_pop(new->endpoints,intlistpool,&ignore_int
				       intlistpool_trace(__FILE__,__LINE__));
      new->nmismatches = Intlistpool_pop(new->nmismatches,intlistpool,&ignore_int
					 intlistpool_trace(__FILE__,__LINE__));
      new->ref_nmismatches = Intlistpool_pop(new->ref_nmismatches,intlistpool,&ignore_int
					     intlistpool_trace(__FILE__,__LINE__));
    } else {
      Intlist_head_set(new->endpoints,splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_main);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_main);
    }

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);


    qend_endpoints = Intlistpool_copy(Intlist_next(fusion->endpoints),intlistpool); /* Skip first elt in common with main */
    qend_nmismatches = Intlistpool_copy(fusion->nmismatches,intlistpool);
    qend_ref_nmismatches = Intlistpool_copy(fusion->ref_nmismatches,intlistpool);

    if (splice_distance == 0) {
      Intlist_head_set(qend_nmismatches,-1);
      Intlist_head_set(qend_ref_nmismatches,-1);
    } else {
      /* Intlist_head_set(qend_endpoints,splice_querypos); -- first elt skipped */
      Intlist_head_set(qend_nmismatches,nmismatches_fusion);
      Intlist_head_set(qend_ref_nmismatches,ref_nmismatches_fusion);
    }

    new->endpoints = Intlist_append(new->endpoints,qend_endpoints);
    new->nmismatches = Intlist_append(new->nmismatches,qend_nmismatches);
    new->ref_nmismatches = Intlist_append(new->ref_nmismatches,qend_ref_nmismatches);

    /* new->qstart_alts should be a copy of old->qstart_alts */
    if (new->qend_alts != NULL) {
      Altsplice_free(&new->qend_alts,pathpool);
    }
    new->qend_alts = Altsplice_copy(fusion->qend_alts,pathpool,vectorpool);

    new->splice3p = fusion->splice3p;
    new->splicetype3 = fusion->splicetype3;
    new->ambig_prob_3 = fusion->ambig_prob_3;

    new->main_univdiagonal = old->main_univdiagonal;
    if (splice_distance == 0) {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(Univcoordlist_next(fusion->univdiagonals),univcoordlistpool));
      new->junctions = List_append(new->junctions,Junction_copy_list(fusion->junctions,listpool,pathpool));

    } else {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(fusion->univdiagonals,univcoordlistpool));

      new->junctions = List_reverse(new->junctions);
      new->junctions = Listpool_push(new->junctions,listpool,
				     (void *) Junction_new_splice(splice_distance,new->sensedir,
								  donor_prob,acceptor_prob,pathpool)
				     listpool_trace(__FILE__,__LINE__));
      new->junctions = List_reverse(new->junctions);
      new->junctions = List_append(new->junctions,Junction_copy_list(fusion->junctions,listpool,pathpool));
    }

    debug2(printf("Standard case 3:\n"));
    debug2(Path_print(new));
    return new;

  } else {
    /* minus */
    debug2(printf("Standard (non-fusion) case 4: fusion -> old, splice querypos %d, splice distance %u\n",
		  splice_querypos,splice_distance));
    new = Path_copy(fusion,intlistpool,univcoordlistpool,listpool,
		    pathpool,vectorpool,transcriptpool,hitlistpool);
    new->nmatches = -1;
    new->ref_nmatches = -1;

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);

    if (splice_distance == 0) {
      /* Merge pieces instead */
      new->endpoints = Intlistpool_pop(new->endpoints,intlistpool,&ignore_int
				       intlistpool_trace(__FILE__,__LINE__));
      new->nmismatches = Intlistpool_pop(new->nmismatches,intlistpool,&ignore_int
					 intlistpool_trace(__FILE__,__LINE__));
      new->ref_nmismatches = Intlistpool_pop(new->ref_nmismatches,intlistpool,&ignore_int
					     intlistpool_trace(__FILE__,__LINE__));
    } else {
      Intlist_head_set(new->endpoints,querylength - splice_querypos);
      Intlist_head_set(new->nmismatches,nmismatches_fusion);
      Intlist_head_set(new->ref_nmismatches,ref_nmismatches_fusion);
    }

    new->endpoints = Intlist_reverse(new->endpoints);
    new->nmismatches = Intlist_reverse(new->nmismatches);
    new->ref_nmismatches = Intlist_reverse(new->ref_nmismatches);


    qend_endpoints = Intlistpool_copy(Intlist_next(old->endpoints),intlistpool); /* Skip first elt in common with main */
    qend_nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
    qend_ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);

    if (splice_distance == 0) {
      Intlist_head_set(qend_nmismatches,-1);
      Intlist_head_set(qend_ref_nmismatches,-1);
    } else {
      /* Intlist_head_set(qend_endpoints,querylength - splice_querypos); -- first elt skipped */
      Intlist_head_set(qend_nmismatches,nmismatches_main);
      Intlist_head_set(qend_ref_nmismatches,ref_nmismatches_main);
    }

    new->endpoints = Intlist_append(new->endpoints,qend_endpoints);
    new->nmismatches = Intlist_append(new->nmismatches,qend_nmismatches);
    new->ref_nmismatches = Intlist_append(new->ref_nmismatches,qend_ref_nmismatches);

    /* new->qstart_alts should be a copy of fusion->qstart_alts */
    if (new->qend_alts != NULL) {
      Altsplice_free(&new->qend_alts,pathpool);
    }
    new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);

    new->splice3p = old->splice3p;
    new->splicetype3 = old->splicetype3;
    new->ambig_prob_3 = old->ambig_prob_3;

    new->main_univdiagonal = old->main_univdiagonal;
    if (splice_distance == 0) {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(Univcoordlist_next(old->univdiagonals),univcoordlistpool));
      new->junctions = List_append(new->junctions,Junction_copy_list(old->junctions,listpool,pathpool));
    } else {
      new->univdiagonals = Univcoordlist_append(new->univdiagonals,
						Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool));

      new->junctions = List_reverse(new->junctions);
      new->junctions = Listpool_push(new->junctions,listpool,
				     (void *) Junction_new_splice(splice_distance,new->sensedir,
								  donor_prob,acceptor_prob,pathpool)
				     listpool_trace(__FILE__,__LINE__));
      new->junctions = List_reverse(new->junctions);
      new->junctions = List_append(new->junctions,Junction_copy_list(old->junctions,listpool,pathpool));
    }

    debug2(printf("Standard case 4:\n"));
    debug2(Path_print(new));
    return new;
  }
}



#if 0
static void
collect_univdiags (List_T *gplus_univdiags, List_T *gminus_univdiags, Stage1_T stage1,
		   Univdiagpool_T univdiagpool) {
  List_T univdiags, p;
  Elt_T elt;
  Univdiag_T *univdiag_array;
  Univcoord_T univdiagonal;
  int qstart, qend;
  int ndiags, k, i, j;


  /* plus */
  univdiags = (List_T) NULL;
  for (p = stage1->queryfwd_plus_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    for (k = 0; k < elt->n_all_univdiags; k++) {
      univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,elt->all_univdiags[k]
					     univdiagpool_trace(__FILE__,__LINE__));
    }
  }
    
  for (p = stage1->queryrev_plus_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    for (k = 0; k < elt->n_all_univdiags; k++) {
      univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,elt->all_univdiags[k]
					     univdiagpool_trace(__FILE__,__LINE__));
    }
  }

  /* Consolidate duplicates, using code taken from extension-search.c */
  /* Needed to avoid multiple fusion paths, which would yield no results */
  /* univdiag_array = (Univdiag_T *) List_to_array_n(&ndiags,univdiags); */

  ndiags = List_length(univdiags);
  univdiag_array = Univdiag_new_array(ndiags);
  Univdiag_transfer_list_to_array(univdiag_array,univdiags);

  qsort(univdiag_array[0],ndiags,sizeof(struct Univdiag_T),Univdiag_struct_diagonal_fwd_cmp);
  Univdiagpool_free_list(&univdiags,univdiagpool
			 univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */

  *gplus_univdiags = (List_T) NULL;
  i = 0;
  while (i < ndiags) {
    univdiagonal = univdiag_array[i]->univdiagonal;
    qstart = univdiag_array[i]->qstart;
    qend = univdiag_array[i]->qend;
      
    j = i+1;
    while (j < ndiags && univdiag_array[j]->univdiagonal == univdiagonal) {
      debug(printf("At left diagonal %u, combining %d..%d with %d..%d\n",
		   univdiagonal,qstart,qend,univdiag_array[j]->qstart,univdiag_array[j]->qend));
      if (univdiag_array[j]->qstart < qstart) {
	qstart = univdiag_array[j]->qstart;
      }
      if (univdiag_array[j]->qend > qend) {
	qend = univdiag_array[j]->qend;
      }
      j++;
    }
      
    debug(printf("Pushing %u, %d..%d onto gplus_univdiags\n",univdiagonal,qstart,qend));
    if (qstart == univdiag_array[i]->qstart && qend == univdiag_array[i]->qend) {
      *gplus_univdiags = Univdiagpool_push_existing(*gplus_univdiags,univdiagpool,univdiag_array[i]
						    univdiagpool_trace(__FILE__,__LINE__));
    } else {
      *gplus_univdiags = Univdiagpool_push(*gplus_univdiags,univdiagpool,qstart,qend,/*nmismatches*/-1,univdiagonal
					   univdiagpool_trace(__FILE__,__LINE__));
    }	

    i = j;
  }
  Univdiag_free_array(&univdiag_array);


  /* minus */
  univdiags = (List_T) NULL;
  for (p = stage1->queryfwd_minus_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    for (k = 0; k < elt->n_all_univdiags; k++) {
      univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,elt->all_univdiags[k]
					     univdiagpool_trace(__FILE__,__LINE__));
    }
  }
    
  for (p = stage1->queryrev_minus_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    for (k = 0; k < elt->n_all_univdiags; k++) {
      univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,elt->all_univdiags[k]
					     univdiagpool_trace(__FILE__,__LINE__));
    }
  }


  /* Consolidate duplicates, using code taken from extension-search.c */
  /* Needed to avoid multiple fusion paths, which would yield no results */
  /* univdiag_array = (Univdiag_T *) List_to_array_n(&ndiags,univdiags); */

  ndiags = List_length(univdiags);
  univdiag_array = Univdiag_new_array(ndiags);
  Univdiag_transfer_list_to_array(univdiag_array,univdiags);

  qsort(univdiag_array[0],ndiags,sizeof(struct Univdiag_T),Univdiag_struct_diagonal_fwd_cmp);
  Univdiagpool_free_list(&univdiags,univdiagpool
			 univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */

  *gminus_univdiags = (List_T) NULL;
  i = 0;
  while (i < ndiags) {
    univdiagonal = univdiag_array[i]->univdiagonal;
    qstart = univdiag_array[i]->qstart;
    qend = univdiag_array[i]->qend;
      
    j = i+1;
    while (j < ndiags && univdiag_array[j]->univdiagonal == univdiagonal) {
      debug(printf("At left diagonal %u, combining %d..%d with %d..%d\n",
		   univdiagonal,qstart,qend,univdiag_array[j]->qstart,univdiag_array[j]->qend));
      if (univdiag_array[j]->qstart < qstart) {
	qstart = univdiag_array[j]->qstart;
      }
      if (univdiag_array[j]->qend > qend) {
	qend = univdiag_array[j]->qend;
      }
      j++;
    }
      
    debug(printf("Pushing %u, %d..%d onto gminus_univdiags\n",univdiagonal,qstart,qend));
    if (qstart == univdiag_array[i]->qstart && qend == univdiag_array[i]->qend) {
      *gminus_univdiags = Univdiagpool_push_existing(*gminus_univdiags,univdiagpool,univdiag_array[i]
						     univdiagpool_trace(__FILE__,__LINE__));
    } else {
      *gminus_univdiags = Univdiagpool_push(*gminus_univdiags,univdiagpool,qstart,qend,/*nmismatches*/-1,univdiagonal
					    univdiagpool_trace(__FILE__,__LINE__));
    }	

    i = j;
  }
  Univdiag_free_array(&univdiag_array);

  return;
}
#endif



/* We first try finding outer fusions against complete paths (which
   might have an ambig splice end) and unextended paths, which can
   contain indels and splicing.  These procedures are likely to work
   if the fusion fragment is long.  But extension search may have
   suboptimal elts that were not searched, especially if the fusion
   fragment is short, so we then try using univdiags from extension
   search */

#if 0
static List_T
compute_fusions_case1 (Path_T fusion, Path_T mainpath, List_T fusion_paths,
		       Univcoord_T univdiagonal3, int querystart3, int queryend3,
		       Shortread_T queryseq, int querylength,
		       Compress_T query_compress_fusion, Compress_T query_compress_main,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,
		       Stage1_T stage1, Knownsplicing_T knownsplicing,
		       int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
		       bool sense_forward_p, int genestrand, int endtrim_allowed,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal5;
  int querystart5, queryend5;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal5 = Univcoordlist_last_value(fusion->univdiagonals); /* Want highest */
  querystart5 = Intlist_penultimate_value(fusion->endpoints);
  if (fusion->junctions != NULL) {
    querystart5 += Junction_ninserts((Junction_T) List_last_value(fusion->junctions,NULL));
  }
  queryend5 = Intlist_last_value(fusion->endpoints);
  /* query_compress_fusion = query_compress_fwd; */
  /* queryptr_fusion = Shortread_queryuc_ptr(queryseq); */
  debug(printf("CASE 1: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
  
  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = Intlist_head(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_fusion,/*plusp_5*/true,
				      /*chroffset5*/fusion->chroffset,
				      /*query_compress_3*/query_compress_main,/*plusp_3:mainpath->plusp*/true,
				      /*chroffset3*/mainpath->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
	    
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 1: splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_querystart(mainpath,fusion,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					    /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case1p (Path_T fusion, Path_T mainpath, List_T fusion_paths,
			Univcoord_T univdiagonal5, int querystart5, int queryend5,
			Shortread_T queryseq, int querylength,
			Compress_T query_compress_fusion, Compress_T query_compress_main,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Stage1_T stage1, Knownsplicing_T knownsplicing,
			int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
			bool sense_forward_p, int genestrand, int endtrim_allowed,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal3;
  int querystart3, queryend3;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal3 = Univcoordlist_last_value(fusion->univdiagonals); /* Want highest */
  querystart3 = querylength - Intlist_last_value(fusion->endpoints);
  queryend3 = querylength - Intlist_penultimate_value(fusion->endpoints);
  if (fusion->junctions != NULL) {
    queryend3 -= Junction_ninserts((Junction_T) List_last_value(fusion->junctions,NULL));
  }
  /* query_compress_fusion = query_compress_rev; */
  /* queryptr_fusion = Shortread_queryrc(queryseq); */
  debug(printf("CASE 1': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));

  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = Intlist_head(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_main,/*plusp_5:mainpath->plusp*/false,
				      /*chroffset5*/mainpath->chroffset,
				      /*query_compress_3*/query_compress_fusion,/*plusp_3*/false,
				      /*chroffset3*/fusion->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
      
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 1': splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_queryend(mainpath,fusion,splice_querypos,querylength,
					  /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					  /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					  donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					  queryseq,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case2 (Path_T fusion, Path_T mainpath, List_T fusion_paths,
		       Univcoord_T univdiagonal3, int querystart3, int queryend3,
		       Shortread_T queryseq, int querylength,
		       Compress_T query_compress_fusion, Compress_T query_compress_main,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,
		       Stage1_T stage1, Knownsplicing_T knownsplicing,
		       int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
		       bool sense_forward_p, int genestrand, int endtrim_allowed,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal5;
  int querystart5, queryend5;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal5 = Univcoordlist_head(fusion->univdiagonals); /* Want lowest */
  querystart5 = querylength - Intlist_second_value(fusion->endpoints);
  /* No need to revise querystart5 for Junction_ninserts in reverse direction */
  queryend5 = querylength - Intlist_head(fusion->endpoints);
  /* query_compress_fusion = query_compress_rev; */
  /* queryptr_fusion = Shortread_queryrc(queryseq); */
  debug(printf("CASE 2: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
  
  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = querylength - Intlist_last_value(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_fusion,/*plusp_5*/false,
				      /*chroffset5*/fusion->chroffset,
				      /*query_compress_3*/query_compress_main,/*plusp_3:mainpath->plusp*/true,
				      /*chroffset3*/mainpath->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
      
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 2: splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_querystart(mainpath,fusion,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					    /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case2p (Path_T fusion, Path_T mainpath, List_T fusion_paths,
			Univcoord_T univdiagonal5, int querystart5, int queryend5,
			Shortread_T queryseq, int querylength,
			Compress_T query_compress_fusion, Compress_T query_compress_main,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Stage1_T stage1, Knownsplicing_T knownsplicing,
			int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
			bool sense_forward_p, int genestrand, int endtrim_allowed,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal3;
  int querystart3, queryend3;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal3 = Univcoordlist_head(fusion->univdiagonals); /* Want lowest */
  querystart3 = Intlist_head(fusion->endpoints);
  queryend3 = Intlist_second_value(fusion->endpoints);
  /* No need to revise queryend3 for Junction_ninserts in forward direction */
  /* query_compress_fusion = query_compress_fwd; */
  /* queryptr_fusion = Shortread_queryuc_ptr(queryseq); */
  debug(printf("CASE 2': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));

  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = querylength - Intlist_last_value(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_main,/*plusp_5:mainpath->plusp*/false,
				      /*chroffset5*/mainpath->chroffset,
				      /*query_compress_3*/query_compress_fusion,/*plusp_3*/true,
				      /*chroffset3*/fusion->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
	
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 2': splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_queryend(mainpath,fusion,splice_querypos,querylength,
					  /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					  /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					  donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					  queryseq,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case3 (Path_T fusion, Path_T mainpath, List_T fusion_paths,
		       Univcoord_T univdiagonal5, int querystart5, int queryend5,
		       Shortread_T queryseq, int querylength,
		       Compress_T query_compress_fusion, Compress_T query_compress_main,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,
		       Stage1_T stage1, Knownsplicing_T knownsplicing,
		       int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
		       bool sense_forward_p, int genestrand, int endtrim_allowed,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool) {
  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal3;
  int querystart3, queryend3;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal3 = Univcoordlist_head(fusion->univdiagonals); /* Want lowest */
  querystart3 = Intlist_head(fusion->endpoints);
  queryend3 = Intlist_second_value(fusion->endpoints);
  /* No need to revise queryend3 for Junction_ninserts in forward direction */
  /* query_compress_fusion = query_compress_fwd; */
  /* queryptr_fusion = Shortread_queryuc_ptr(queryseq); */
  debug(printf("CASE 3: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
  
  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = querylength - Intlist_last_value(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_main,/*plusp_5:mainpath->plusp*/true,
				      /*chroffset5*/mainpath->chroffset,
				      /*query_compress_3*/query_compress_fusion,/*plusp_3*/true,
				      /*chroffset3*/fusion->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
      
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 3: splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_queryend(mainpath,fusion,splice_querypos,querylength,
					  /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					  /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					  donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					  queryseq,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case3p (Path_T fusion, Path_T mainpath, List_T fusion_paths,
			Univcoord_T univdiagonal3, int querystart3, int queryend3,
			Shortread_T queryseq, int querylength,
			Compress_T query_compress_fusion, Compress_T query_compress_main,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Stage1_T stage1, Knownsplicing_T knownsplicing,
			int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
			bool sense_forward_p, int genestrand, int endtrim_allowed,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal5;
  int querystart5, queryend5;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal5 = Univcoordlist_head(fusion->univdiagonals); /* Want lowest */
  querystart5 = querylength - Intlist_second_value(fusion->endpoints);
  /* No need to revise querystart5 for Junction_ninserts in reverse direction */
  queryend5 = querylength - Intlist_head(fusion->endpoints);
  /* query_compress_fusion = query_compress_rev; */
  /* queryptr_fusion = Shortread_queryrc(queryseq); */
  debug(printf("CASE 3': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
  
  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = querylength - Intlist_last_value(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_fusion,/*plusp_5*/false,
				      /*chroffset5*/fusion->chroffset,
				      /*query_compress_3*/query_compress_main,/*plusp_3:mainpath->plusp*/false,
				      /*chroffset3*/mainpath->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
      
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 3': splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_querystart(mainpath,fusion,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					    /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case4 (Path_T fusion, Path_T mainpath, List_T fusion_paths,
		       Univcoord_T univdiagonal5, int querystart5, int queryend5,
		       Shortread_T queryseq, int querylength,
		       Compress_T query_compress_fusion, Compress_T query_compress_main,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,
		       Stage1_T stage1, Knownsplicing_T knownsplicing,
		       int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
		       bool sense_forward_p, int genestrand, int endtrim_allowed,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal3;
  int querystart3, queryend3;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */
    
  univdiagonal3 = Univcoordlist_last_value(fusion->univdiagonals); /* Want highest */
  querystart3 = querylength - Intlist_last_value(fusion->endpoints);
  queryend3 = querylength - Intlist_penultimate_value(fusion->endpoints);
  if (fusion->junctions != NULL) {
    queryend3 -= Junction_ninserts((Junction_T) List_last_value(fusion->junctions,NULL));
  }
  /* query_compress_fusion = query_compress_rev; */
  /* queryptr_fusion = Shortread_queryrc(queryseq); */
  debug(printf("CASE 4: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));

  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = Intlist_head(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_main,/*plusp_5:mainpath->plusp*/true,
				      /*chroffset5*/mainpath->chroffset,
				      /*query_compress_3*/query_compress_fusion,/*plusp_3*/false,
				      /*chroffset3*/fusion->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 4: splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_queryend(mainpath,fusion,splice_querypos,querylength,
					  /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					  /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					  donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					  queryseq,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
static List_T
compute_fusions_case4p (Path_T fusion, Path_T mainpath, List_T fusion_paths,
			Univcoord_T univdiagonal3, int querystart3, int queryend3,
			Shortread_T queryseq, int querylength,
			Compress_T query_compress_fusion, Compress_T query_compress_main,
			Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Stage1_T stage1, Knownsplicing_T knownsplicing,
			int nmismatches_allowed, int max_insertionlen, int max_deletionlen,
			bool sense_forward_p, int genestrand, int endtrim_allowed,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			Hitlistpool_T hitlistpool) {

  Path_T newpath;
  int splice_querypos, indel_pos, nindels, distal_trim;
  Univcoord_T univdiagonal5;
  int querystart5, queryend5;
  int nmismatches_5, nmismatches_3, ref_nmismatches_5, ref_nmismatches_3;
  double donor_prob, acceptor_prob;
  char donor1, donor2, acceptor1, acceptor2;

  int found_score_ignore = querylength;

  debug(Path_print(fusion));
  /* assert(old->extendedp == false); -- could be called again by Pathpair_eval_and_sort */

  univdiagonal5 = Univcoordlist_last_value(fusion->univdiagonals); /* Want highest */
  querystart5 = Intlist_penultimate_value(fusion->endpoints);
  if (fusion->junctions != NULL) {
    querystart5 += Junction_ninserts((Junction_T) List_last_value(fusion->junctions,NULL));
  }
  queryend5 = Intlist_last_value(fusion->endpoints);
  /* query_compress_fusion = query_compress_fwd; */
  /* queryptr_fusion = Shortread_queryuc_ptr(queryseq); */
  debug(printf("CASE 4' (unextended): %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));

  if (querystart5 >= querystart3 || queryend5 >= queryend3) {
    /* Skip */
  } else if (mainpath->chrnum != fusion->chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion->chrnum] == true)) {
    /* Skip */
  } else {
    distal_trim = Intlist_head(fusion->endpoints);
    if (distal_trim > endtrim_allowed) {
      /* Skip */
    } else if ((splice_querypos =
		Splice_resolve_fusion(&donor1,&donor2,&acceptor1,&acceptor2,&nindels,&indel_pos,
				      &nmismatches_5,&nmismatches_3,
				      &ref_nmismatches_5,&ref_nmismatches_3,
				      &donor_prob,&acceptor_prob,univdiagonal5,univdiagonal3,
				      /*query_compress_5*/query_compress_fusion,/*plusp_5*/true,
				      /*chroffset5*/fusion->chroffset,
				      /*query_compress_3*/query_compress_main,/*plusp_3:mainpath->plusp*/false,
				      /*chroffset3*/mainpath->chroffset,
				      querystart5,queryend3,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,sense_forward_p,genestrand,
				      nmismatches_allowed,max_insertionlen,max_deletionlen)) < 0) {
      /* Skip: no fusion found */
      
    } else if (nindels == 0) {
      /* Splice only */
      debug(printf("CASE 4': splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		   querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
      debug(Path_print(mainpath));
      
      newpath = Path_fusion_copy_querystart(mainpath,fusion,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					    /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/false);
      debug(Path_print(newpath));
      Path_eval_nmatches(&found_score_ignore,&end_score_ignore,newpath,query_compress_fwd,query_compress_rev);
      fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  return fusion_paths;
}
#endif


#if 0
/* Returns or replaces the original mainpath */
/* Handles cases 1 and 2 */
List_T
Path_fusion_outer_querystart_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				   Compress_T query_compress_fwd, Compress_T query_compress_rev,
				   Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				   int genestrand, int nmismatches_allowed,
				   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				   Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  /* Path_T *complete_paths_gplus, *complete_paths_gminus; */
  /* int n_complete_paths_gplus, n_complete_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;

  int type;
  int splice_querypos, trim_querypos;
  int nmismatches_5, nmismatches_3, nmismatches_to_trimpos;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_outer_querystart_plus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  /* Unextended_paths could have been run through Path_extend, but did not extend completely */
#if 0
  if (mainpath->sensedir != SENSE_FORWARD) {
    return (List_T) NULL;
  } else {
    sense_forward_p = true;
  }
#else
  assert(mainpath->sensedir == SENSE_FORWARD);
#endif


  /* path3 = mainpath; */
  univdiagonal3 = Univcoordlist_head(mainpath->univdiagonals);
  querystart3 = Intlist_head(mainpath->endpoints);
  queryend3 = Intlist_second_value(mainpath->endpoints);
  /* No need to revise queryend3 for Junction_ninserts in forward direction */
  query_compress_main = query_compress_fwd;
    

  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal5 = stage1->extension_gplus[i];
    querystart5 = stage1->extension_qstart_gplus[i];
    queryend5 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 1 (univdiag): %d..%d vs %s => %d..%d\n",
		   querystart5,queryend5,Intlist_to_string(mainpath->endpoints),querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */
      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_fusion,/*plusDp*/true,
				      /*chroffset_D*/fusion_chroffset,
				      /*query_compress_A*/query_compress_main,/*pluspAp:mainpath->plusp*/true,
				      /*chroffset_A*/mainpath->chroffset,
				      /*D*/querystart5,/*A*/queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_querystart_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
					       /*univdiagonal*/univdiagonal5,querylength,
					       /*pos5*/0,/*pos3*/splice_querypos,
					       /*plusp*/true,genestrand);
	debug(printf("TRIM_QUERYPOS %d, SPLICE_QUERYPOS %d\n",trim_querypos,splice_querypos));
	if (splice_querypos - trim_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  endpoints = Intlistpool_push(endpoints,intlistpool,querystart5
				       intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querystart5 + splice_querypos)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/true,/*first_read_p*/true,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);
			       
	  newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
						/*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
						/*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
						donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
						queryseq,intlistpool,uintlistpool,univcoordlistpool,
						listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
						/*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal5 = stage1->extension_gminus[i];
    querystart5 = querylength - stage1->extension_qend_gminus[i];
    queryend5 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 2 (univdiag): %d..%d vs %s => %d..%d\n",
		   querystart5,queryend5,
		   Intlist_to_string(mainpath->endpoints),querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */
      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_fusion,/*plusDp*/false,
				      /*chroffset_D*/fusion_chroffset,
				      /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/true,
				      /*chroffset_A*/mainpath->chroffset,
				      /*D*/querystart5,/*A*/queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_querystart_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = querylength - Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
							   /*univdiagonal*/univdiagonal5,querylength,
							   /*pos5*/querylength - splice_querypos,/*pos3*/querylength,
							   /*plusp*/false,genestrand);
	debug(printf("TRIM_QUERYPOS %d, SPLICE_QUERYPOS %d\n",trim_querypos,splice_querypos));
	if (splice_querypos - trim_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,querylength - querystart5
				       intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	  endpoints = Intlistpool_push(endpoints,intlistpool,querylength - splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - splice_querypos + querylength - querystart5)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/false,/*first_read_p*/true,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);
	
	  newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
						/*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
						/*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
						donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
						queryseq,intlistpool,uintlistpool,univcoordlistpool,
						listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
						/*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(1) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
/* Returns or replaces the original mainpath */
/* Handles cases 2' and 1' */
List_T
Path_fusion_outer_querystart_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				    Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				    int genestrand, int nmismatches_allowed,
				    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				    Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  /* Path_T *complete_paths_gplus, *complete_paths_gminus; */
  /* int n_complete_paths_gplus, n_complete_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;

  int type;
  int splice_querypos, trim_querypos;
  int nmismatches_5, nmismatches_3, nmismatches_to_trimpos;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_outer_querystart_minus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  /* Unextended_paths could have been run through Path_extend, but did not extend completely */
#if 0
  if (mainpath->sensedir == SENSE_FORWARD) {
    return (List_T) NULL;
  } else {
    sense_forward_p = false;
  }
#else
  assert(mainpath->sensedir == SENSE_ANTI);
#endif

  /* path5 = mainpath; */
  univdiagonal5 = Univcoordlist_head(mainpath->univdiagonals);
  querystart5 = querylength - Intlist_second_value(mainpath->endpoints);
  /* No need to revise querystart5 for Junction_ninserts in reverse direction */
  queryend5 = querylength - Intlist_head(mainpath->endpoints);
  query_compress_main = query_compress_rev;


  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal3 = stage1->extension_gplus[i];
    querystart3 = stage1->extension_qstart_gplus[i];
    queryend3 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 2' (univdiag): %s => %d..%d vs %d..%d\n",
		   Intlist_to_string(mainpath->endpoints),querystart5,queryend5,
		   querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/false,
					  /*chroffset_A*/mainpath->chroffset,
					  /*query_compress_D*/query_compress_fusion,/*plusDp*/true,
					  /*chroffset_D*/fusion_chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_querystart_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
					     /*univdiagonal*/univdiagonal3,querylength,
					     /*pos5*/splice_querypos,/*pos3*/querylength,
					     /*plusp*/true,genestrand);
	debug(printf("SPLICE_QUERYPOS %d, TRIM_QUERYPOS %d\n",splice_querypos,trim_querypos));
	if (trim_querypos - splice_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,queryend3
				       intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	  endpoints = Intlistpool_push(endpoints,intlistpool,splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (splice_querypos + queryend3)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/true,/*first_read_p*/false,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);
	
	  newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					      /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal3 = stage1->extension_gminus[i];
    querystart3 = querylength - stage1->extension_qend_gminus[i];
    queryend3 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 1' (univdiag): %s => %d..%d vs %d..%d\n",
		   Intlist_to_string(mainpath->endpoints),querystart5,queryend5,
		   querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/false,
					  /*chroffset_A*/mainpath->chroffset,
					  /*query_compress_D*/query_compress_fusion,/*plusDp*/false,
					  /*chroffset_D*/fusion_chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_querystart_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = querylength - Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
							     /*univdiagonal*/univdiagonal3,querylength,
							     /*pos5*/0,/*pos3*/querylength - splice_querypos,
							     /*plusp*/false,genestrand);
	debug(printf("SPLICE_QUERYPOS %d, TRIM_QUERYPOS %d\n",splice_querypos,trim_querypos));
	if (trim_querypos - splice_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,querylength - splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  endpoints = Intlistpool_push(endpoints,intlistpool,querylength - queryend3
				       intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - queryend3 + querylength - splice_querypos)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/false,/*first_read_p*/false,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	  newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					      /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(2) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
/* Returns or replaces the original mainpath */
/* Handles cases 3 and 4 */
List_T
Path_fusion_outer_queryend_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				 Compress_T query_compress_fwd, Compress_T query_compress_rev,
				 Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				 int genestrand, int nmismatches_allowed,
				 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				 Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  /* Path_T *complete_paths_gplus, *complete_paths_gminus; */
  /* int n_complete_paths_gplus, n_complete_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;

  int type;
  int splice_querypos, trim_querypos;
  int nmismatches_5, nmismatches_3, nmismatches_to_trimpos;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_outer_queryend_plus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  /* Unextended_paths could have been run through Path_extend, but did not extend completely */
#if 0
  if (mainpath->sensedir != SENSE_FORWARD) {
    return (List_T) NULL;
  } else {
    sense_forward_p = true;
  }
#else
  assert(mainpath->sensedir == SENSE_FORWARD);
#endif

  /* plus: path5 = mainpath */
  univdiagonal5 = Univcoordlist_last_value(mainpath->univdiagonals);
  querystart5 = Intlist_penultimate_value(mainpath->endpoints);
  if (mainpath->junctions != NULL) {
    querystart5 += Junction_ninserts((Junction_T) List_last_value(mainpath->junctions,NULL));
  }
  queryend5 = Intlist_last_value(mainpath->endpoints);
  query_compress_main = query_compress_fwd;
    

  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal3 = stage1->extension_gplus[i];
    querystart3 = stage1->extension_qstart_gplus[i];
    queryend3 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 3 (univdiag): %s => %d..%d vs %d..%d\n",
		   Intlist_to_string(mainpath->endpoints),querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/true,
				      /*chroffset_D*/mainpath->chroffset,
				      /*query_compress_A*/query_compress_fusion,/*plusAp*/true,
				      /*chroffset_A*/fusion_chroffset,
				      querystart5,queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_queryend_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
					     /*univdiagonal*/univdiagonal3,querylength,
					     /*pos5*/splice_querypos,/*pos3*/querylength,
					     /*plusp*/true,genestrand);
	debug(printf("SPLICE_QUERYPOS %d, TRIM_QUERYPOS %d\n",splice_querypos,trim_querypos));
	if (trim_querypos - splice_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,queryend3
				       intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	  endpoints = Intlistpool_push(endpoints,intlistpool,splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (splice_querypos + queryend3)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/true,/*first_read_p*/false,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	  newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					      /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal3 = stage1->extension_gminus[i];
    querystart3 = querylength - stage1->extension_qend_gminus[i];
    queryend3 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 4 (univdiag): %s => %d..%d vs %d..%d\n",
		   Intlist_to_string(mainpath->endpoints),querystart5,queryend5,
		   querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/true,
				      /*chroffset_D*/mainpath->chroffset,
				      /*query_compress_A*/query_compress_fusion,/*plusAp*/false,
				      /*chroffset_A*/fusion_chroffset,
				      querystart5,queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_queryend_plus (univdiags): splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = querylength - Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
							     /*univdiagonal*/univdiagonal3,querylength,
							     /*pos5*/0,/*pos3*/querylength - splice_querypos,
							     /*plusp*/false,genestrand);
	debug(printf("SPLICE_QUERYPOS %d, TRIM_QUERYPOS %d\n",splice_querypos,trim_querypos));
	if (trim_querypos - splice_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,querylength - splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  endpoints = Intlistpool_push(endpoints,intlistpool,querylength - queryend3
				       intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - queryend3 + querylength - splice_querypos)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/false,/*first_read_p*/false,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	  newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					      /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(3) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
List_T
Path_fusion_outer_queryend_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
				  int genestrand,  int nmismatches_allowed,
				  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				  Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  /* Path_T *complete_paths_gplus, *complete_paths_gminus; */
  /* int n_complete_paths_gplus, n_complete_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;

  int type;
  int splice_querypos, trim_querypos;
  int nmismatches_5, nmismatches_3, nmismatches_to_trimpos;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;

  debug(printf("Entered Path_fusion_outer_queryend_minus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  /* Unextended_paths could have been run through Path_extend, but did not extend completely */
#if 0
  if (mainpath->sensedir == SENSE_FORWARD) {
    sense_forward_p = true;
  } else {
    sense_forward_p = false;
  }
#else
  assert(mainpath->sensedir == SENSE_ANTI);
#endif


  /* minus: path3 = mainpath */
  univdiagonal3 = Univcoordlist_last_value(mainpath->univdiagonals);
  querystart3 = querylength - Intlist_last_value(mainpath->endpoints);
  queryend3 = querylength - Intlist_penultimate_value(mainpath->endpoints);
  if (mainpath->junctions != NULL) {
    queryend3 -= Junction_ninserts((Junction_T) List_last_value(mainpath->junctions,NULL));
  }
  query_compress_main = query_compress_rev;


  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal5 = stage1->extension_gplus[i];
    querystart5 = stage1->extension_qstart_gplus[i];
    queryend5 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 4' (univdiag): %d..%d vs %s => %d..%d\n",
		   querystart5,queryend5,Intlist_to_string(mainpath->endpoints),querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_fusion,/*plusAp*/true,
					  /*chroffset_A*/fusion_chroffset,
					  /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/false,
					  /*chroffset_D*/mainpath->chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_queryend_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
					       /*univdiagonal*/univdiagonal5,querylength,
					       /*pos5*/0,/*pos3*/splice_querypos,
					       /*plusp*/true,genestrand);
	debug(printf("TRIM_QUERYPOS %d, SPLICE_QUERYPOS %d\n",trim_querypos,splice_querypos));
	if (splice_querypos - trim_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  endpoints = Intlistpool_push(endpoints,intlistpool,querystart5
				       intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querystart5 + splice_querypos)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/true,/*first_read_p*/true,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	  newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
						/*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
						/*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
						donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
						queryseq,intlistpool,uintlistpool,univcoordlistpool,
						listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
						/*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }
  
  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal5 = stage1->extension_gminus[i];
    querystart5 = querylength - stage1->extension_qend_gminus[i];
    queryend5 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 3' (univdiag): %d..%d vs %s => %d..%d\n",
		   querystart5,queryend5,Intlist_to_string(mainpath->endpoints),querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_fusion,/*plusAp*/false,
					  /*chroffset_A*/fusion_chroffset,
					  /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/false,
					  /*chroffset_D*/mainpath->chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_outer_queryend_minus (univdiags): splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	trim_querypos = querylength - Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_fusion,genomebits,
							   /*univdiagonal*/univdiagonal5,querylength,
							   /*pos5*/querylength - splice_querypos,/*pos3*/querylength,
							   /*plusp*/false,genestrand);
	debug(printf("TRIM_QUERYPOS %d, SPLICE_QUERYPOS %d\n",trim_querypos,splice_querypos));
	if (splice_querypos - trim_querypos < MIN_OUTER_END) {
	  /* Skip: not a good alignment */
	} else {
	  endpoints = Intlistpool_push(NULL,intlistpool,querylength - querystart5
				       intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	  endpoints = Intlistpool_push(endpoints,intlistpool,querylength - splice_querypos
				       intlistpool_trace(__FILE__,__LINE__));
	  /* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - splice_querypos + querylength - querystart5)/2; */
	  univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
					 intlistpool_trace(__FILE__,__LINE__));
	  ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					     intlistpool_trace(__FILE__,__LINE__));

	  middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				    endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				    /*plusp*/false,/*first_read_p*/true,genestrand,
				    /*sensedir*/mainpath->sensedir,querylength,
				    /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				    /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				    query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	  newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
						/*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
						/*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
						donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
						queryseq,intlistpool,uintlistpool,univcoordlistpool,
						listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
						/*fusion_created_p*/true);
	  debug(Path_print(newpath));
	  debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	  Path_free(&middle_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				      hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(4) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


T
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
			Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Path_T newpath, mainpath;
  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoord_T goal_univdiagonal, low_univdiagonal, high_univdiagonal;
  Univcoordlist_T univdiagonals;

  Univcoord_T *new_univdiagonals = NULL;
  int ndiagonals = 0;

  Univcoord_T univdiagonal5, univdiagonal3;
  Compress_T query_compress_fusion5;
  int pos5, pos3, querystart5;

  int type;
  int splice_querypos;
  int local_nmismatches, nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


#ifdef DEBUG
  printf("Entered Path_fusion_inner_qend, sensedir %d\n",fusion5->sensedir);
  printf("Anchor: "); Path_print(anchor3);
  printf("Fusion: "); Path_print(fusion5);
#endif

  if (main_plusp == true) {
    if (fusion5->plusp == true) {
      debug(printf("CASE 1: %s\n",Intlist_to_string(fusion5->endpoints)));

      /* distal_trim = Intlist_head(fusion5->endpoints); */
      univdiagonal5 = Univcoordlist_last_value(fusion5->univdiagonals); /* Want highest */
      query_compress_fusion5 = query_compress_fwd;
      pos5 = Intlist_last_value(fusion5->endpoints);

      querystart5 = Intlist_penultimate_value(fusion5->endpoints);
      if (fusion5->junctions != NULL) {
	querystart5 += Junction_ninserts((Junction_T) List_last_value(fusion5->junctions,NULL));
      }
      
    } else {
      debug(printf("CASE 2: %s\n",Intlist_to_string(fusion5->endpoints)));

      /* distal_trim = querylength - Intlist_last_value(fusion5->endpoints); */
      univdiagonal5 = Univcoordlist_head(fusion5->univdiagonals); /* Want lowest */
      query_compress_fusion5 = query_compress_rev;
      pos5 = querylength - Intlist_head(fusion5->endpoints);

      querystart5 = querylength - Intlist_second_value(fusion5->endpoints);
      /* No need to revise querystart5 for Junction_ninserts in reverse direction */
    }

    goal_univdiagonal = subtract_bounded(Path_genomiclow(anchor3) + fusion5->querylength,
					 /*expected_pairlength*/30000,anchor3->chroffset + anchor3->querylength);
    debug(printf("Goal: Path_genomiclow3 %u + querylength5 %d - 30000, bounded by chroffset %u + querylength3 %d\n",
		 Path_genomiclow(anchor3),fusion5->querylength,anchor3->chroffset,anchor3->querylength));
    
    low_univdiagonal = subtract_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
					anchor3->chroffset + anchor3->querylength);
    high_univdiagonal = add_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
				    Univcoordlist_head(anchor3->univdiagonals));

    debug(printf("Looking for qend hit in %d..%d at %u..%u\n",
		 pos5,querylength,low_univdiagonal,high_univdiagonal));

    if (low_univdiagonal > high_univdiagonal) {
      debug(printf("Skipping because %u > %u\n",low_univdiagonal,high_univdiagonal));
      ndiagonals = 0;
      new_univdiagonals = (Univcoord_T *) NULL;
    } else {
      new_univdiagonals =
	Spliceends_qend_resolve(&ndiagonals,pos5,querylength,
				low_univdiagonal,high_univdiagonal,
				query_compress_main,queryptr_main,main_plusp,genestrand,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,stage1,
				nmismatches_allowed);
    }

  } else {
    if (fusion5->plusp == true) {
      debug(printf("CASE 4': %s\n",Intlist_to_string(fusion5->endpoints)));

      /* distal_trim = Intlist_head(fusion5->endpoints); */
      univdiagonal5 = Univcoordlist_last_value(fusion5->univdiagonals); /* Want highest */
      query_compress_fusion5 = query_compress_fwd;
      pos3 = querylength - Intlist_last_value(fusion5->endpoints);

      querystart5 = Intlist_penultimate_value(fusion5->endpoints);
      if (fusion5->junctions != NULL) {
	querystart5 += Junction_ninserts((Junction_T) List_last_value(fusion5->junctions,NULL));
      }

    } else {
      debug(printf("CASE 3': %s\n",Intlist_to_string(fusion5->endpoints)));

      /* distal_trim = querylength - Intlist_last_value(fusion5->endpoints); */
      univdiagonal5 = Univcoordlist_head(fusion5->univdiagonals); /* Want lowest */
      query_compress_fusion5 = query_compress_rev;
      pos3 = Intlist_head(fusion5->endpoints);

      querystart5 = querylength - Intlist_second_value(fusion5->endpoints);
      /* No need to revise querystart5 for Junction_ninserts in reverse direction */
    }

    goal_univdiagonal = add_bounded(Path_genomichigh(anchor3) - anchor3->querylength,
				    /*expected_pairlength*/30000,anchor3->chrhigh);
    low_univdiagonal = subtract_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
					Path_genomichigh(anchor3));
    high_univdiagonal = add_bounded(goal_univdiagonal,/*pairlength_delta*/30000,anchor3->chrhigh);

    debug(printf("Looking for qstart hit in 0..%d at %u..%u\n",
		 pos3,low_univdiagonal,high_univdiagonal));

    if (low_univdiagonal > high_univdiagonal) {
      ndiagonals = 0;
      new_univdiagonals = (Univcoord_T *) NULL;
    } else {
      new_univdiagonals =
	Spliceends_qstart_resolve(&ndiagonals,pos3,querylength,
				  low_univdiagonal,high_univdiagonal,
				  query_compress_main,queryptr_main,main_plusp,genestrand,
				  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				  localdb_tally,localdb_overflowp,max_nblocks,stage1,
				  nmismatches_allowed);
    }
  }

  if (ndiagonals == 0) {
    return (T) NULL;

  } else {
    debug(printf("FOUND: %u\n",new_univdiagonals[0]));
    univdiagonal3 = new_univdiagonals[0];
    debug(printf("(1) Looking for splice_querypos in %d..%d at %u..%u\n",
		 querystart5,querylength,univdiagonal5,univdiagonal3));
    if (main_chrnum != fusion5->chrnum && (circularp[main_chrnum] == true || circularp[fusion5->chrnum] == true)) {
      return (T) NULL;

    } else if (main_plusp == true) {
      /* Case 1, 2 */
      splice_querypos =
	Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
			    /*D*/&nmismatches_5,/*A*/&nmismatches_3,
			    /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
			    /*D*/univdiagonal5,/*A*/univdiagonal3,
			    /*query_compress_D*/query_compress_fusion5,/*plusDp*/fusion5->plusp,
			    /*chroffset_D*/fusion5->chroffset,
			    /*query_compress_A*/query_compress_main,/*plusAp*/main_plusp,
			    /*chroffset_A*/main_chroffset,
			    querystart5,/*queryend3*/querylength,querylength,
			    spliceinfo,knownsplicing,intlistpool,genestrand);
    } else {
      /* Case 3', 4' */
      splice_querypos =
	Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				/*A*/&nmismatches_5,/*D*/&nmismatches_3,
				/*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
				/*A*/univdiagonal5,/*D*/univdiagonal3,
				/*query_compress_A*/query_compress_fusion5,/*plusAp*/fusion5->plusp,
				/*chroffset_A*/fusion5->chroffset,
				/*query_compress_D*/query_compress_main,/*plusDp*/main_plusp,
				/*chroffset_D*/main_chroffset,
				querystart5,/*queryend3*/querylength,querylength,
				spliceinfo,knownsplicing,intlistpool,genestrand);
    }

    if (splice_querypos < 0) {
      /* Skip: no fusion found */
      return (T) NULL;

    } else {
      debug(printf("SPLICE_QUERYPOS %d\n",splice_querypos));
      if (main_plusp == true) {
	endpoints = Intlistpool_push(NULL,intlistpool,querylength
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (new_univdiagonals[0] - querylength) + (querylength + splice_querypos)/2; */
      } else {
	endpoints = Intlistpool_push(NULL,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,0
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (new_univdiagonals[0[ - querylength) + (querylength - splice_querypos)/2;*/
      }
      univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,new_univdiagonals[0]
					     univcoordlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
				     intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					 intlistpool_trace(__FILE__,__LINE__));

      debug(printf("Creating mainpath with univdiagonal %u\n",new_univdiagonals[0]));
      mainpath = Path_create(/*main_univdiagonal*/new_univdiagonals[0],
			     endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
			     main_plusp,/*first_read_p*/true,genestrand,
			     fusion5->sensedir,querylength,fusion5->method,
			     main_chrnum,main_chroffset,main_chrhigh,
			     /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			     /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			     /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
			     query_compress_fwd,query_compress_rev,pathpool,vectorpool);

      newpath = Path_fusion_copy_querystart(mainpath,fusion5,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					    /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/false);
      Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
      if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
	debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
	Path_free(&newpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	Path_free(&mainpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (T) NULL;
      } else {
	Path_free(&mainpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	return newpath;
      }

    }
  }
}


T
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
			  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Path_T newpath, mainpath;
  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoord_T goal_univdiagonal, low_univdiagonal, high_univdiagonal;
  Univcoordlist_T univdiagonals;

  Univcoord_T *new_univdiagonals = NULL;
  int ndiagonals = 0;

  Univcoord_T univdiagonal5, univdiagonal3;
  Compress_T query_compress_fusion3;
  int pos5, pos3, queryend3;

  int type;
  int splice_querypos;
  int local_nmismatches, nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


#ifdef DEBUG
  printf("Entered Path_fusion_inner_qstart, sensedir %d\n",fusion3->sensedir);
  printf("Anchor: "); Path_print(anchor5);
  printf("Fusion: "); Path_print(fusion3);
#endif

  if (main_plusp == true) {
    if (fusion3->plusp == true) {
      debug(printf("CASE 3: %s\n",Intlist_to_string(fusion3->endpoints)));

      /* distal_trim = querylength - Intlist_last_value(fusion3->endpoints); */
      univdiagonal3 = Univcoordlist_head(fusion3->univdiagonals); /* Want lowest */
      query_compress_fusion3 = query_compress_fwd;
      pos3 = Intlist_head(fusion3->endpoints);

      queryend3 = Intlist_second_value(fusion3->endpoints);
      /* No need to revise queryend3 for Junction_ninserts in forward direction */
      
    } else {
      debug(printf("CASE 4: %s\n",Intlist_to_string(fusion3->endpoints)));

      /* distal_trim = Intlist_head(fusion3->endpoints); */
      univdiagonal3 = Univcoordlist_last_value(fusion3->univdiagonals); /* Want highest */
      query_compress_fusion3 = query_compress_rev;
      pos3 = querylength - Intlist_last_value(fusion3->endpoints);

      queryend3 = querylength - Intlist_penultimate_value(fusion3->endpoints);
      if (fusion3->junctions != NULL) {
	queryend3 -= Junction_ninserts((Junction_T) List_last_value(fusion3->junctions,NULL));
      }
    }

    goal_univdiagonal = add_bounded(Path_genomichigh(anchor5) - anchor5->querylength,
				    /*expected_pairlength*/30000,anchor5->chrhigh);
    low_univdiagonal = subtract_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
					Path_genomichigh(anchor5));
    high_univdiagonal = add_bounded(goal_univdiagonal,/*pairlength_delta*/30000,anchor5->chrhigh);

    debug(printf("Looking for qstart hit in %d..%d at %u..%u\n",
		 0,pos3,low_univdiagonal,high_univdiagonal));

    if (low_univdiagonal > high_univdiagonal) {
      ndiagonals = 0;
      new_univdiagonals = (Univcoord_T *) NULL;
    } else {
      new_univdiagonals =
	Spliceends_qstart_resolve(&ndiagonals,pos3,querylength,
				  low_univdiagonal,high_univdiagonal,
				  query_compress_main,queryptr_main,main_plusp,genestrand,
				  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				  localdb_tally,localdb_overflowp,max_nblocks,stage1,
				  nmismatches_allowed);
    }

  } else {
    if (fusion3->plusp == true) {
      debug(printf("CASE 2': %s\n",Intlist_to_string(fusion3->endpoints)));

      /* distal_trim = querylength - Intlist_last_value(fusion3->endpoints); */
      univdiagonal3 = Univcoordlist_head(fusion3->univdiagonals); /* Want lowest */
      query_compress_fusion3 = query_compress_fwd;
      pos5 = querylength - Intlist_head(fusion3->endpoints);

      queryend3 = Intlist_second_value(fusion3->endpoints);
      /* No need to revise queryend3 for Junction_ninserts in forward direction */

    } else {
      debug(printf("CASE 1': %s\n",Intlist_to_string(fusion3->endpoints)));

      /* distal_trim = Intlist_head(fusion3->endpoints) */
      univdiagonal3 = Univcoordlist_last_value(fusion3->univdiagonals); /* Want highest */
      query_compress_fusion3 = query_compress_rev;
      pos5 = Intlist_last_value(fusion3->endpoints);

      queryend3 = querylength - Intlist_penultimate_value(fusion3->endpoints);
      if (fusion3->junctions != NULL) {
	queryend3 -= Junction_ninserts((Junction_T) List_last_value(fusion3->junctions,NULL));
      }
    }

    goal_univdiagonal = subtract_bounded(Path_genomiclow(anchor5) + anchor5->querylength + fusion3->querylength,
					 /*expected_pairlength*/30000,anchor5->chroffset + anchor5->querylength);
    low_univdiagonal = subtract_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
					anchor5->chroffset + anchor5->querylength);
    high_univdiagonal = add_bounded(goal_univdiagonal,/*pairlength_delta*/30000,
				    Univcoordlist_head(anchor5->univdiagonals));

    debug(printf("Looking for qend hit in %d..%d at %u..%u\n",
		 pos5,querylength,low_univdiagonal,high_univdiagonal));

    if (low_univdiagonal > high_univdiagonal) {
      ndiagonals = 0;
      new_univdiagonals = (Univcoord_T *) NULL;
    } else {
      new_univdiagonals =
	Spliceends_qend_resolve(&ndiagonals,pos5,querylength,
				low_univdiagonal,high_univdiagonal,
				query_compress_main,queryptr_main,main_plusp,genestrand,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,stage1,
				nmismatches_allowed);
    }
  }

  if (ndiagonals == 0) {
    return (T) NULL;

  } else {
    debug(printf("FOUND %u\n",new_univdiagonals[0]));
    univdiagonal5 = new_univdiagonals[0];
    debug(printf("(2) Looking for splice_querypos in %d..%d at %u..%u\n",
		 0,queryend3,univdiagonal5,univdiagonal3));
    if (main_chrnum != fusion3->chrnum && (circularp[main_chrnum] == true || circularp[fusion3->chrnum] == true)) {
      return (T) NULL;

    } else if (main_plusp == true) {
      /* Case 3, 4 */
      splice_querypos =
	Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
			    /*D*/&nmismatches_5,/*A*/&nmismatches_3,
			    /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
			    /*D*/univdiagonal5,/*A*/univdiagonal3,
			    /*query_compress_D*/query_compress_main,/*plusDp*/main_plusp,
			    /*chroffset_D*/main_chroffset,
			    /*query_compress_A*/query_compress_fusion3,/*plusAp*/fusion3->plusp,
			    /*chroffset_A*/fusion3->chroffset,
			    /*querystart5*/0,queryend3,querylength,
			    spliceinfo,knownsplicing,intlistpool,genestrand);

    } else {
      /* Case 1', 2' */
      splice_querypos =
	Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				/*A*/&nmismatches_5,/*D*/&nmismatches_3,
				/*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
				/*A*/univdiagonal5,/*D*/univdiagonal3,
				/*query_compress_A*/query_compress_main,/*plusAp*/main_plusp,
				/*chroffset_A*/main_chroffset,
				/*query_compress_D*/query_compress_fusion3,/*plusDp*/fusion3->plusp,
				/*chroffset_D*/fusion3->chroffset,
				/*querystart5*/0,queryend3,querylength,
				spliceinfo,knownsplicing,intlistpool,genestrand);
    }

    if (splice_querypos < 0) {
      /* Skip: no fusion found */
      return (T) NULL;

    } else {
      debug(printf("SPLICE_QUERYPOS %d\n",splice_querypos));
      if (main_plusp == true) {
	endpoints = Intlistpool_push(NULL,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,0
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (new_univdiagonals[0] - querylength) + splice_querypos/2; */
      } else {
	endpoints = Intlistpool_push(NULL,intlistpool,querylength
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (new_univdiagonals[0] - querylength) + (querylength + querylength - splice_querypos)/2; */
      }
      univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,new_univdiagonals[0]
					     univcoordlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
				     intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					 intlistpool_trace(__FILE__,__LINE__));

      debug(printf("Creating mainpath with univdiagonal %u\n",new_univdiagonals[0]));
      mainpath = Path_create(/*main_univdiagonal*/new_univdiagonals[0],
			     endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
			     main_plusp,/*first_read_p*/false,genestrand,
			     fusion3->sensedir,querylength,fusion3->method,
			     main_chrnum,main_chroffset,main_chrhigh,
			     /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			     /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			     /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
			     query_compress_fwd,query_compress_rev,pathpool,vectorpool);

      newpath = Path_fusion_copy_queryend(mainpath,fusion3,splice_querypos,querylength,
					  /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					  /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					  donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					  queryseq,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					  /*fusion_created_p*/false);
      Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
      if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
	debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
	Path_free(&newpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	Path_free(&mainpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (T) NULL;
      } else {
	Path_free(&mainpath,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	return newpath;
      }
    }
  }
}


/************************************************************************
 *  For single-end reads.  Procedures return only the fusion path, if found.
 ************************************************************************/

#if 0
List_T
Path_fusion_querystart_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			     int genestrand, int nmismatches_allowed,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  bool sense_forward_p;

  int type;
  int splice_querypos;
  int nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_querystart_plus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  if (mainpath->sensedir == SENSE_ANTI) {
    return (List_T) NULL;
  } else {
    sense_forward_p = true;
  }

  /* path3 = mainpath; */
  univdiagonal3 = Univcoordlist_head(mainpath->univdiagonals);
  querystart3 = Intlist_head(mainpath->endpoints);
  queryend3 = Intlist_second_value(mainpath->endpoints);
  /* No need to revise queryend3 for Junction_ninserts in forward direction */
  query_compress_main = query_compress_fwd;
    

  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal5 = stage1->extension_gplus[i];
    querystart5 = stage1->extension_qstart_gplus[i];
    queryend5 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 1: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_fusion,/*plusDp*/true,
				      /*chroffset_D*/fusion_chroffset,
				      /*query_compress_A*/query_compress_main,/*pluspAp:mainpath->plusp*/true,
				      /*chroffset_A*/mainpath->chroffset,
				      /*D*/querystart5,/*A*/queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_querystart_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,querystart5
				     intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querystart5 + splice_querypos)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/true,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);
			       
	newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					      /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal5 = stage1->extension_gminus[i];
    querystart5 = querylength - stage1->extension_qend_gminus[i];
    queryend5 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 2: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_fusion,/*plusDp*/false,
				      /*chroffset_D*/fusion_chroffset,
				      /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/true,
				      /*chroffset_A*/mainpath->chroffset,
				      /*D*/querystart5,/*A*/queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_querystart_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,querylength - querystart5
				     intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	endpoints = Intlistpool_push(endpoints,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - splice_querypos + querylength - querystart5)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/false,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);
	
	newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					      /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(5) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
/* Returns or replaces the original mainpath */
/* Handles cases 2' and 1' */
List_T
Path_fusion_querystart_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			      int genestrand, int nmismatches_allowed,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  bool sense_forward_p;

  int type;
  int splice_querypos;
  int nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_querystart_minus, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  if (mainpath->sensedir == SENSE_FORWARD) {
    return (List_T) NULL;
  } else {
    sense_forward_p = false;
  }

  /* path5 = mainpath; */
  univdiagonal5 = Univcoordlist_head(mainpath->univdiagonals);
  querystart5 = querylength - Intlist_second_value(mainpath->endpoints);
  /* No need to revise querystart5 for Junction_ninserts in reverse direction */
  queryend5 = querylength - Intlist_head(mainpath->endpoints);
  query_compress_main = query_compress_rev;


  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal3 = stage1->extension_gplus[i];
    querystart3 = stage1->extension_qstart_gplus[i];
    queryend3 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 2': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/false,
					  /*chroffset_A*/mainpath->chroffset,
					  /*query_compress_D*/query_compress_fusion,/*plusDp*/true,
					  /*chroffset_D*/fusion_chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_querystart_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,queryend3
				     intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	endpoints = Intlistpool_push(endpoints,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (splice_querypos + queryend3)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/true,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);
	
	newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					    /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal3 = stage1->extension_gminus[i];
    querystart3 = querylength - stage1->extension_qend_gminus[i];
    queryend3 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 1': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_main,/*plusAp:mainpath->plusp*/false,
					  /*chroffset_A*/mainpath->chroffset,
					  /*query_compress_D*/query_compress_fusion,/*plusDp*/false,
					  /*chroffset_D*/fusion_chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_querystart_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,querylength - queryend3
				     intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - queryend3 + querylength - splice_querypos)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/false,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					    /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(6) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
/* Returns or replaces the original mainpath */
/* Handles cases 3 and 4 */
List_T
Path_fusion_queryend_plus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			   int genestrand, int nmismatches_allowed,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  bool sense_forward_p;

  int type;
  int splice_querypos;
  int nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;


  debug(printf("Entered Path_fusion_queryend, sensedir %d, with path\n",mainpath->sensedir));
  debug(Path_print(mainpath));

  if (mainpath->sensedir == SENSE_ANTI) {
    return (List_T) NULL;
  } else {
    sense_forward_p = true;
  }

  /* plus: path5 = mainpath */
  univdiagonal5 = Univcoordlist_last_value(mainpath->univdiagonals);
  querystart5 = Intlist_penultimate_value(mainpath->endpoints);
  if (mainpath->junctions != NULL) {
    querystart5 += Junction_ninserts((Junction_T) List_last_value(mainpath->junctions,NULL));
  }
  queryend5 = Intlist_last_value(mainpath->endpoints);
  query_compress_main = query_compress_fwd;
    

  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal3 = stage1->extension_gplus[i];
    querystart3 = stage1->extension_qstart_gplus[i];
    queryend3 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 3: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/true,
				      /*chroffset_D*/mainpath->chroffset,
				      /*query_compress_A*/query_compress_fusion,/*plusAp*/true,
				      /*chroffset_A*/fusion_chroffset,
				      querystart5,queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_queryend_plus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,queryend3
				     intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	endpoints = Intlistpool_push(endpoints,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (splice_querypos + queryend3)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/true,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					    /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal3 = stage1->extension_gminus[i];
    querystart3 = querylength - stage1->extension_qend_gminus[i];
    queryend3 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal3 - querylength + querystart3,
				univdiagonal3 - querylength + queryend3);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querylength - queryend3 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 4: %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
				      /*D*/&nmismatches_5,/*A*/&nmismatches_3,
				      /*D*/&ref_nmismatches_5,/*A*/&ref_nmismatches_3,
				      /*D*/univdiagonal5,/*A*/univdiagonal3,
				      /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/true,
				      /*chroffset_D*/mainpath->chroffset,
				      /*query_compress_A*/query_compress_fusion,/*plusAp*/false,
				      /*chroffset_A*/fusion_chroffset,
				      querystart5,queryend3,querylength,
				      stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_queryend_plus (univdiags): splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,querylength - queryend3
				     intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - queryend3 + querylength - splice_querypos)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_3
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_3
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal3,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/false,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	newpath = Path_fusion_copy_queryend(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					    /*nmismatches_main*/nmismatches_5,/*ref_nmismatches_main*/ref_nmismatches_5,
					    /*nmismatches_fusion*/nmismatches_3,/*ref_nmismatches_fusion*/ref_nmismatches_3,
					    donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					    queryseq,intlistpool,uintlistpool,univcoordlistpool,
					    listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					    /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */

#ifdef DEBUG
  printf("(7) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


#if 0
List_T
Path_fusion_queryend_minus (int *found_score, int *end_score, T mainpath, Stage1_T stage1,

			    Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    Shortread_T queryseq, int querylength, Knownsplicing_T knownsplicing,
			    int genestrand, int nmismatches_allowed,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Univdiagpool_T univdiagpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  
  List_T result = NULL;
  List_T fusion_paths, gplus_univdiags = NULL, gminus_univdiags = NULL, p;
  /* List_T unextended_paths_gplus, unextended_paths_gminus; */
  int i;

  Intlist_T endpoints, nmismatches, ref_nmismatches;
  Univcoordlist_T univdiagonals;
  T newpath, middle_path;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset, fusion_chrhigh;
  /* char *queryptr_fusion; */

  Univcoord_T univdiagonal5, univdiagonal3;
  int querystart5, queryend5, querystart3, queryend3;
  Compress_T query_compress_main, query_compress_fusion;
  bool sense_forward_p;

  int type;
  int splice_querypos;
  int nmismatches_5, nmismatches_3;
  int ref_nmismatches_5, ref_nmismatches_3;
  char donor1, donor2, acceptor1, acceptor2;
  double donor_prob, acceptor_prob;

  if (mainpath->sensedir == SENSE_FORWARD) {
    return (List_T) NULL;
  } else {
    sense_forward_p = false;
  }

  /* minus: path3 = mainpath */
  univdiagonal3 = Univcoordlist_last_value(mainpath->univdiagonals);
  querystart3 = querylength - Intlist_last_value(mainpath->endpoints);
  queryend3 = querylength - Intlist_penultimate_value(mainpath->endpoints);
  if (mainpath->junctions != NULL) {
    queryend3 -= Junction_ninserts((Junction_T) List_last_value(mainpath->junctions,NULL));
  }
  query_compress_main = query_compress_rev;


  /* Previously looked at unextended paths, but we are no longer generating those */

  /* Try using univdiagonals from extension search */
  fusion_paths = (List_T) NULL;
  /* collect_univdiags(&gplus_univdiags,&gminus_univdiags,stage1,univdiagpool); */

  for (i = 0; i < stage1->nextension_gplus; i++) {
    univdiagonal5 = stage1->extension_gplus[i];
    querystart5 = stage1->extension_qstart_gplus[i];
    queryend5 = stage1->extension_qend_gplus[i];
    query_compress_fusion = query_compress_fwd;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 4': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_fusion,/*plusAp*/true,
					  /*chroffset_A*/fusion_chroffset,
					  /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/false,
					  /*chroffset_D*/mainpath->chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_queryend_minus (univdiags): splice_querypos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	endpoints = Intlistpool_push(endpoints,intlistpool,querystart5
				     intlistpool_trace(__FILE__,__LINE__)); /* 0 */
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querystart5 + splice_querypos)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/true,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					      /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }
  
  for (i = 0; i < stage1->nextension_gminus; i++) {
    univdiagonal5 = stage1->extension_gminus[i];
    querystart5 = querylength - stage1->extension_qend_gminus[i];
    queryend5 = querylength - stage1->extension_qstart_gminus[i];
    query_compress_fusion = query_compress_rev;

    fusion_chrnum = EF64_chrnum(&fusion_chroffset,&fusion_chrhigh,chromosome_ef64,
				univdiagonal5 - querylength + querystart5,
				univdiagonal5 - querylength + queryend5);

    /* Cannot rely on qstart or qend as being the endtrim because they are approximate from extension search */
    if (/* querystart5 <= endtrim_allowed && */
	querystart5 < querystart3 && queryend5 < queryend3) {
      debug(printf("CASE 3': %d..%d vs %d..%d\n",querystart5,queryend5,querystart3,queryend3));
      if (mainpath->chrnum != fusion_chrnum && (circularp[mainpath->chrnum] == true || circularp[fusion_chrnum] == true)) {
	/* Skip */

      } else if ((splice_querypos =
		  Splice_fusion_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&donor_prob,&acceptor_prob,
					  /*A*/&nmismatches_5,/*D*/&nmismatches_3,
					  /*A*/&ref_nmismatches_5,/*D*/&ref_nmismatches_3,
					  /*A*/univdiagonal5,/*D*/univdiagonal3,
					  /*query_compress_A*/query_compress_fusion,/*plusAp*/false,
					  /*chroffset_A*/fusion_chroffset,
					  /*query_compress_D*/query_compress_main,/*plusDp:mainpath->plusp*/false,
					  /*chroffset_D*/mainpath->chroffset,
					  querystart5,queryend3,querylength,
					  stage1->spliceinfo,knownsplicing,intlistpool,genestrand)) < 0) {
	/* Skip: no fusion found */

      } else {
	debug(printf("Path_fusion_queryend_minus (univdiags): splice_qpos in range %d..%d is %d with fusion, mismatches %d+%d, and probs %f and %f\n",
		     querystart5,queryend3,splice_querypos,nmismatches_5,nmismatches_3,donor_prob,acceptor_prob));
	debug(Path_print(mainpath));

	endpoints = Intlistpool_push(NULL,intlistpool,querylength - querystart5
				     intlistpool_trace(__FILE__,__LINE__)); /* querylength */
	endpoints = Intlistpool_push(endpoints,intlistpool,querylength - splice_querypos
				     intlistpool_trace(__FILE__,__LINE__));
	/* middle_univcoord = (univdiag->univdiagonal - querylength) + (querylength - splice_querypos + querylength - querystart5)/2; */
	univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal5
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches_5
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches_5
					   intlistpool_trace(__FILE__,__LINE__));

	middle_path = Path_create(/*main_univdiagonal*/univdiagonal5,
				  endpoints,univdiagonals,nmismatches,ref_nmismatches,/*junctions*/NULL,
				  /*plusp*/false,/*first_read_p*/true,genestrand,
				  /*sensedir*/mainpath->sensedir,querylength,
				  /*method*/FUSION,fusion_chrnum,fusion_chroffset,fusion_chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  /*qstart_alts*/(Altsplice_T) NULL,/*qend_alts*/(Altsplice_T) NULL,
				  query_compress_fwd,query_compress_rev,pathpool,vectorpool);

	newpath = Path_fusion_copy_querystart(mainpath,/*fusion*/middle_path,splice_querypos,querylength,
					      /*nmismatches_main*/nmismatches_3,/*ref_nmismatches_main*/ref_nmismatches_3,
					      /*nmismatches_fusion*/nmismatches_5,/*ref_nmismatches_fusion*/ref_nmismatches_5,
					      donor1,donor2,acceptor1,acceptor2,donor_prob,acceptor_prob,
					      queryseq,intlistpool,uintlistpool,univcoordlistpool,
					      listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
					      /*fusion_created_p*/true);

	Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
	debug(Path_print(newpath));
	debug(printf("Freeing middle path %p and pushing newpath %p\n",middle_path,newpath));
	Path_free(&middle_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	fusion_paths = Hitlist_push(fusion_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  Univdiagpool_gc(&gminus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  Univdiagpool_gc(&gplus_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_T */
  
#ifdef DEBUG
  printf("(8) From univdiags, got %d fusion_paths\n",List_length(fusion_paths));
  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_print(newpath);
  }
#endif

  for (p = fusion_paths; p != NULL; p = List_next(p)) {
    newpath = (T) List_head(p);
    Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev);
    if (0 && querylength - Path_eval_nmatches(&(*found_score),&(*end_score),newpath,query_compress_fwd,query_compress_rev) > nmismatches_allowed) {
      debug(printf("querylength %d - nmatches %d > nmismatches_allowed %d\n",querylength,newpath->nmatches,nmismatches_allowed));
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      result = Hitlist_push(result,hitlistpool,(void *) newpath
			    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&fusion_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by hitlistpool */

  return result;
}
#endif


void
Path_fusion_setup (bool *circularp_in, EF64_T chromosome_ef64_in, Univcoord_T genomelength_in,
		   Chrpos_T shortsplicedist_in, Transcriptome_T transcriptome_in,
		   Genomebits_T genomebits_in) {

  circularp = circularp_in;
  chromosome_ef64 = chromosome_ef64_in;
  genomelength = genomelength_in;
  shortsplicedist = shortsplicedist_in;
  transcriptome = transcriptome_in;
  genomebits = genomebits_in;
  
  return;
}

