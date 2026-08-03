#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "pathpair.h"
#include "debug_groups.h"
#include "path-eval.h"
#include "path-trim.h"
#include "path-solve.h"
#include "path-print-alignment.h"
#include "path-print-m8.h"
#include "spliceends.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "junction.h"


static Outputtype_T output_type;
static bool only_concordant_p = false; /* CONCORDANT_UNIQ, CONCORDANT_MULT, and CONCORDANT_CIRC */
static bool omit_concordant_uniq_p = false;
static bool omit_concordant_mult_p = false;

static Chrpos_T shortsplicedist;
static int max_insertlength;

static Univcoord_T genomelength;



#define CONCORDANT_TEXT "concordant"
#define PAIRED_TEXT "paired"
#define UNPAIRED_TEXT "unpaired"

#define INSERTLENGTH_FACTOR 1.5


#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
#define debug0(x) x
#else
#define debug0(x)
#endif

/* compute_insertlength */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Pathpair_resolve */
#if defined(DEBUG3) || defined(DEBUG_PAIRED_DEEP)
#define debug3(x) x
#else
#define debug3(x)
#endif


#define T Pathpair_T


bool
Pathpair_insertlength_knownp (T this) {

  if (Intlist_last_value(this->pathL->endpoints) < this->pathL->querylength) {
    return false;
  } else if (Intlist_head(this->pathH->endpoints) > 0) {
    return false;
  } else if (this->insertlength == 0) {
    return false;
  } else {
    return true;
  }
}


bool
Pathpair_outerlength_knownp (T this) {

  if (Intlist_head(this->pathL->endpoints) > 0) {
    return false;
  } else if (Intlist_last_value(this->pathH->endpoints) < this->pathH->querylength) {
    return false;
  } else {
    return true;
  }
}


void
Pathpair_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool) {
#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
  static int call_i = 0;
#endif

  debug0(printf("%d: Freeing pathpair %p\n",++call_i,*old));

  Path_free(&(*old)->pathL,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);
  Path_free(&(*old)->pathH,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);
  FREE(*old);

  return;
}


void
Pathpair_gc (List_T *list, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool) {
  List_T p;
  T old;
  
  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Pathpair_free(&old,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
  }
  Hitlistpool_free_list(&(*list),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* allocated by Hitlistpool_push */
  
  return;
}


bool
Pathpair_transcript_intersectp (Path_T pathL, Path_T pathH) {
  if (Transcript_intersectp(pathL->transcripts,pathH->transcripts) == true) {
    return true;
  } else if (Transcript_intersectp(pathL->transcripts,pathH->invalid_transcripts) == true) {
    return true;
  } else if (Transcript_intersectp(pathL->invalid_transcripts,pathH->transcripts) == true) {
    return true;
#if 0
  } else if (Transcript_intersectp(pathL->invalid_transcripts,pathH->invalid_transcripts) == true) {
    return true;
#endif
  } else {
    return false;
  }
}


static Chrpos_T
compute_insertlength (int *pair_relationship, Path_T pathL, Path_T pathH,
		      Shortread_T queryseqL, Shortread_T queryseqH, bool plusp) {
  Univcoord_T insert_start, insert_end;

  insert_end = Path_genomiclow(pathH) + pathH->querylength;
  if (plusp == true) {
    insert_end -= Shortread_left_choplength(queryseqH);
    insert_end += Shortread_right_choplength(queryseqH);
  } else {
    insert_end -= Shortread_right_choplength(queryseqH);
    insert_end += Shortread_left_choplength(queryseqH);
  }

  insert_start = Path_genomichigh(pathL) - pathL->querylength;
  if (plusp == true) {
    insert_start -= Shortread_left_choplength(queryseqL);
    insert_start += Shortread_right_choplength(queryseqL);
  } else {
    insert_start -= Shortread_right_choplength(queryseqL);
    insert_start += Shortread_left_choplength(queryseqL);
  }

#ifdef DEBUG1
  printf("Computing insertlength for\n");
  Path_print(pathL);
  printf("insert_start %u\n",insert_start);
  Path_print(pathH);
  printf("insert_end %u\n",insert_end);
#endif  

  if (insert_start > insert_end) {
    *pair_relationship = 0;
    debug1(printf("%p and %p: insert_start %u, insert_end %u => returning -1U\n",
		  pathL,pathH,insert_start,insert_end));
    return (Chrpos_T) -1;
  } else if (plusp == true) {
    *pair_relationship = +1;
    debug1(printf("%p and %p: insert_start %u, insert_end %u => returning %u\n",
		  pathL,pathH,insert_start,insert_end,insert_end - insert_start));
    return (Chrpos_T) (insert_end - insert_start);
  } else {
    *pair_relationship = -1;
    debug1(printf("%p and %p: insert_start %u, insert_end %u => returning %u\n",
		  pathL,pathH,insert_start,insert_end,insert_end - insert_start));
    return (Chrpos_T) (insert_end - insert_start);
  }
}


bool
Pathpair_overlappingp (T this) {
  Chrpos_T insertlength;
  int pair_relationship;

  insertlength = compute_insertlength(&pair_relationship,this->pathL,this->pathH,
				      this->queryseqL,this->queryseqH,this->plusp);
  if (insertlength == (Chrpos_T) -1) {
    return true;
  } else {
    return false;
  }
}


Chrpos_T
Pathpair_insertlength (T this) {
  int pair_relationship;

  return compute_insertlength(&pair_relationship,this->pathL,this->pathH,
			      this->queryseqL,this->queryseqH,this->plusp);
}


Chrpos_T
Pathpair_outerlength (T this) {
  return Path_genomichigh(this->pathH) - Path_genomiclow(this->pathL);
}


int
Pathpair_nbadsplices (T this) {
  return Path_nbadsplices(this->pathL) + Path_nbadsplices(this->pathH);
}


void
Pathpair_print (T this) {

  printf("chopped_nmatches %d+%d (region %d..%d and %d..%d), nmatches %d+%d, insertlength %u, outerlength %u\n",
	 Path_eval_nmatches_chopped(this->pathL),Path_eval_nmatches_chopped(this->pathH),
	 this->pathL->chop_qstart,this->pathL->chop_qend,this->pathH->chop_qstart,this->pathH->chop_qend,
	 this->pathL->nmatches,this->pathH->nmatches,Pathpair_insertlength(this),Pathpair_outerlength(this));
  Path_print(this->pathL);
  Path_print(this->pathH);
  return;
}


/* Previously returned some pathpairs in unresolved_pathpairs, but this does not work well when sequence quality is low */
T
Pathpair_new_concordant (List_T *unresolved_pathpairs,
			 Path_T pathL, Path_T pathH, Shortread_T queryseqL, Shortread_T queryseqH, bool plusp,

			 int nmismatches_filter_5, int nmismatches_filter_3,
			 int mincoverage_filter_5, int mincoverage_filter_3,

			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			 bool check_inner_p, bool copyLp, bool copyHp) {
#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
  static int call_i = 0;
#endif

  T new;
  Path_T newpathL, newpathH;


  debug0(printf("\nConsidering pathpair %d for L,H: plusp: %d, copyL %d, copy H %d\n",
		++call_i,plusp,copyLp,copyHp));
  debug0(Path_print(pathL));
  debug0(Path_print(pathH));

  if (pathL->chrnum != pathH->chrnum) {
    debug0(printf("chrnum %d != chrnum %d => returning NULL\n",
		  pathL->chrnum,pathH->chrnum));
    return (T) NULL;

  } else if (Path_genomiclow(pathL) > Path_genomichigh(pathH)) {
    debug0(printf("genomiclow %u > genomichigh %u => returning NULL\n",
		  Path_genomiclow(pathL),Path_genomichigh(pathH)));
    return (T) NULL;

  } else {
    if (copyLp == true) {
      newpathL = Path_copy(pathL,intlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool);
    } else {
      newpathL = pathL;
    }
    if (copyHp == true) {
      newpathH = Path_copy(pathH,intlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool);
    } else {
      newpathH = pathH;
    }
  }


  new = (T) MALLOC(sizeof(*new));
  new->pairtype = CONCORDANT;
  new->plusp = plusp;
  new->insertlength = compute_insertlength(&new->pair_relationship,newpathL,newpathH,
					   queryseqL,queryseqH,plusp);

  if (plusp == true) {
    new->path5 = new->pathL = newpathL;
    new->path3 = new->pathH = newpathH;
  } else {
    new->path3 = new->pathL = newpathL;
    new->path5 = new->pathH = newpathH;
  }

  new->queryseqL = queryseqL;
  new->queryseqH = queryseqH;

  /* Modifies path5 and path3 */
  new->transcript_concordant_p = 
    Transcript_intersection(new->path5,new->path3,listpool,transcriptpool);


  new->outerlength = Path_genomichigh(pathH) - Path_genomiclow(pathL);

  debug0(printf("Creating pathpair %p for 5' and 3' with insertlength %d:\n",
		new,new->insertlength));
  debug0(Pathpair_print(new));

#if 0
  if (check_inner_p == false) {
    return new;
    
  } else if (Path_unextended_qend_p(pathL,/*endtrim_allowed*/20,/*allow_ambig_p*/true) == true) {
    /* Pick 20 to allow for a kmer to not match */
    /* Set allow_ambig_p to be true to allow for Pathpair_resolve later, but may need to add a check in Pathpair_resolve */
    /* Allow outer parts to be unaligned to find outer fusions */
    debug0(printf("pathL is unextended at qend => returning NULL\n"));
    *unresolved_pathpairs = Hitlist_push(*unresolved_pathpairs,hitlistpool,(void *) new
					 hitlistpool_trace(__FILE__,__LINE__));
    return (T) NULL;

  } else if (Path_unextended_qstart_p(pathH,/*endtrim_allowed*/20,/*allow_ambig_p*/true) == true) {
    /* Pick 20 to allow for a kmer to not match */
    /* Set allow_ambig_p to be true to allow for Pathpair_resolve later, but may need to add a check in Pathpair_resolve */
    /* Allow outer parts to be unaligned to find outer fusions */
    debug0(printf("pathH is unextended at qstart => returning NULL\n"));
    *unresolved_pathpairs = Hitlist_push(*unresolved_pathpairs,hitlistpool,(void *) new
					 hitlistpool_trace(__FILE__,__LINE__));
    return (T) NULL;

  } else if (new->path5->score_within_trims > nmismatches_filter_5 ||
	     new->path3->score_within_trims > nmismatches_filter_3) {
    debug0(printf("Low scores within trims => returning NULL\n"));
    *unresolved_pathpairs = Hitlist_push(*unresolved_pathpairs,hitlistpool,(void *) new
					 hitlistpool_trace(__FILE__,__LINE__));
    return (T) NULL;

  } else if (Path_coverage(new->path5) < mincoverage_filter_5 ||
	     Path_coverage(new->path3) < mincoverage_filter_3) {
    debug0(printf("Low coverage => returning NULL\n"));
    *unresolved_pathpairs = Hitlist_push(*unresolved_pathpairs,hitlistpool,(void *) new
					 hitlistpool_trace(__FILE__,__LINE__));
    return (T) NULL;
    
  } else {
    return new;
  }
#else
  return new;
#endif
}


/* Does not need to resolve inner regions */
T
Pathpair_new_inner_fusion (Path_T pathL, Path_T pathH, Shortread_T queryseqL, Shortread_T queryseqH, bool plusp,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			   bool copyLp, bool copyHp) {
  T new;
  Path_T newpathL, newpathH;


  debug0(printf("Considering pathpair for\n"));
  debug0(Path_print(pathL));
  debug0(Path_print(pathH));

  assert(pathL->chrnum == pathH->chrnum);

  if (Path_genomichigh(pathL) > Path_genomichigh(pathH)) {
    debug0(printf("Returning NULL\n"));
    return (T) NULL;
  } else if (Path_genomiclow(pathH) < Path_genomiclow(pathL)) {
    debug0(printf("Returning NULL\n"));
    return (T) NULL;
  } else {
    if (copyLp == true) {
      newpathL = Path_copy(pathL,intlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool);
    } else {
      newpathL = pathL;
    }
    if (copyHp == true) {
      newpathH = Path_copy(pathH,intlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool);
    } else {
      newpathH = pathH;
    }
  }


  new = (T) MALLOC(sizeof(*new));
  new->pairtype = CONCORDANT;
  new->plusp = plusp;
  new->transcript_concordant_p = false;
  new->insertlength = compute_insertlength(&new->pair_relationship,newpathL,newpathH,
					   queryseqL,queryseqH,plusp);

  if (plusp == true) {
    new->path5 = new->pathL = newpathL;
    new->path3 = new->pathH = newpathH;
  } else {
    new->path3 = new->pathL = newpathL;
    new->path5 = new->pathH = newpathH;
  }

  new->queryseqL = queryseqL;
  new->queryseqH = queryseqH;

  /* outerlength is not defined for a fusion */
  new->outerlength = 0;

  debug0(printf("Creating pathpair %p\n",new));
  debug0(Path_print(new->path5));
  debug0(Path_print(new->path3));

  return new;
}



static Junction_T
find_overlap_for_pathL (bool *extendL_p, Univcoord_T *segmentH_high, Univcoord_T segmentL_low, Univcoord_T segmentL_high,
			Path_T pathH, int querylengthH) {
  Junction_T junctionH = NULL;
  int ninserts;
  Univcoord_T univdiagonal_H, segmentH_low;
  Intlist_T p;
  Univcoordlist_T q;
  List_T j;
  
  *extendL_p = false;

  ninserts = 0;		/* For first iteration */
  p = pathH->endpoints;
  q = pathH->univdiagonals;
  j = pathH->junctions;
  while (j != NULL) {
    univdiagonal_H = Univcoordlist_head(q);
    segmentH_low = univdiagonal_H - querylengthH + Intlist_head(p) + ninserts;
    *segmentH_high = univdiagonal_H - querylengthH + Intlist_second_value(p);
    
    junctionH = (Junction_T) List_head(j);
    ninserts = Junction_ninserts(junctionH); /* For next iteration */
    
    if (segmentL_high <= segmentH_low) {
      /* No overlap.  Far apart (as soon as we see this, we should be able to terminate) */
    } else if (*segmentH_high <= segmentL_low) {
      /* No overlap.  Overreach */
    } else {
      debug3(printf("Found overlap with segmentH %u..%u and %u..%u and a junction:  ",
		    segmentL_low,segmentL_high,segmentH_low,*segmentH_high));
      debug3(Junction_print(junctionH));
      return junctionH;
    }

    p = Intlist_next(p);
    q = Univcoordlist_next(q);
    j = List_next(j);
  }
    
  /* Final segment */
  univdiagonal_H = Univcoordlist_head(q);
  segmentH_low = univdiagonal_H - querylengthH + Intlist_head(p) + ninserts;
  *segmentH_high = univdiagonal_H - querylengthH + Intlist_second_value(p);
      
  if (segmentL_high <= segmentH_low) {
    /* No overlap.  Far apart (as soon as we see this, we should be able to terminate) */
  } else if (*segmentH_high <= segmentL_low) {
    /* No overlap.  Overreach */
  } else {
    debug3(printf("Found overlap with final segmentH %u..%u and %u..%u\n",
		  segmentL_low,segmentL_high,segmentH_low,*segmentH_high));
    *extendL_p = true;
  }

  return (Junction_T) NULL;
}


static Junction_T
find_overlap_for_pathH (bool *extendH_p,Univcoord_T *segmentL_low, Univcoord_T segmentH_low, Univcoord_T segmentH_high,
			Path_T pathL, int querylengthL) {
  Junction_T junctionL;
  int ninserts;
  Univcoord_T univdiagonal_L, segmentL_high;
  Intlist_T p;
  Univcoordlist_T q;
  List_T j;
  
  *extendH_p = false;

  pathL->endpoints = Intlist_reverse(pathL->endpoints);
  pathL->univdiagonals = Univcoordlist_reverse(pathL->univdiagonals);
  pathL->junctions = List_reverse(pathL->junctions);

  p = pathL->endpoints;
  q = pathL->univdiagonals;
  j = pathL->junctions;
  while (j != NULL) {
    junctionL = (Junction_T) List_head(j);
    ninserts = Junction_ninserts(junctionL); /* For this iteration */

    univdiagonal_L = Univcoordlist_head(q);
    segmentL_high = univdiagonal_L - querylengthL + Intlist_head(p);
    *segmentL_low = univdiagonal_L - querylengthL + Intlist_second_value(p) + ninserts;
    
    if (segmentL_high <= segmentH_low) {
      /* No overlap.  Far apart (as soon as we see this, we should be able to terminate) */
    } else if (segmentH_high <= *segmentL_low) {
      /* No overlap.  Overreach */
    } else {
      debug3(printf("Found overlap with segmentL %u..%u and %u..%u to a junction: ",
		    *segmentL_low,segmentL_high,segmentH_low,segmentH_high));
      debug3(Junction_print(junctionL));
      
      pathL->endpoints = Intlist_reverse(pathL->endpoints);
      pathL->univdiagonals = Univcoordlist_reverse(pathL->univdiagonals);
      pathL->junctions = List_reverse(pathL->junctions);
      
      return junctionL;
    }

    p = Intlist_next(p);
    q = Univcoordlist_next(q);
    j = List_next(j);
  }

  /* Final segment */
  univdiagonal_L = Univcoordlist_head(q);
  segmentL_high = univdiagonal_L - querylengthL + Intlist_head(p);
  *segmentL_low = univdiagonal_L - querylengthL + Intlist_second_value(p) /*+ ninserts (0)*/;
  
  if (segmentL_high <= segmentH_low) {
    /* No overlap.  Far apart (as soon as we see this, we should be able to terminate) */
  } else if (segmentH_high <= *segmentL_low) {
    /* No overlap.  Overreach */
  } else {
    debug3(printf("Found overlap with the final segmentL %u..%u and %u..%u\n",
		  *segmentL_low,segmentL_high,segmentH_low,segmentH_high));
    *extendH_p = true;
  }
      
  pathL->endpoints = Intlist_reverse(pathL->endpoints);
  pathL->univdiagonals = Univcoordlist_reverse(pathL->univdiagonals);
  pathL->junctions = List_reverse(pathL->junctions);
  
  return (Junction_T) NULL;
}


static bool
preserve_scored_local_mate_softclip_p (Path_T path) {
  if (path->method != LOCAL_MATE) {
    return false;
  } else if (path->junctions != NULL) {
    return false;
  } else if (Path_softclippedp(path) == false) {
    return false;
  } else if (path->nmismatches == NULL || Intlist_head(path->nmismatches) < 0) {
    return false;
  } else {
    return true;
  }
}



/* TODO: Handle the case where both pathL and pathH need to be
   resolved.  Could get the results of Localdb_get on each end, and if
   multiple, then add an altsplice */
void
Pathpair_resolve (int *found_score_5, int *found_score_3, int *end_score_5, int *end_score_3,
		  T this, bool plusp, int genestrand,
		  Compress_T query_compress_L, Compress_T queryL_compress_fwd, Compress_T queryL_compress_rev,
		  Compress_T query_compress_H, Compress_T queryH_compress_fwd, Compress_T queryH_compress_rev,
		  Shortread_T queryseqL, Shortread_T queryseqH, char *queryptrL, char *queryptrH,
		  int querylengthL, int querylengthH, 
		  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		  int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		  int localdb_nmismatches_allowed,
		  
		  Stage1_T stage1L, Stage1_T stage1H, Knownsplicing_T knownsplicing,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Univdiagpool_T univdiagpool, Vectorpool_T vectorpool) {

  Path_T pathL, pathH;
  int pos5_L, pos3_H;		/* Used for anchor_qpos */
#ifdef USE_EXPECTED
  Univcoord_T goal_univdiagonal;
#endif
  Univcoord_T low_univdiagonal, high_univdiagonal;
  int *found_score_ptr_L, *found_score_ptr_H, *end_score_ptr_L, *end_score_ptr_H;

  Univcoord_T univdiagonal_L, univdiagonal_H, distalL, distalH;
  Univcoord_T segmentL_low, segmentL_high, segmentH_low, segmentH_high, segment_length;
  int splice_qpos_L, splice_qpos_H, trimpos_L, trimpos_H;
  int medial_nmismatches_L, distal_nmismatches_L, medial_nmismatches_H, distal_nmismatches_H;
  double medial_prob_L, distal_prob_L, medial_prob_H, distal_prob_H;
  Chrpos_T splice_distance;
  Univcoord_T cutpointL, cutpointH;
  Junction_T junctionL, junctionH;
  int ninserts;
  bool sense_forward_p, extendL_p, extendH_p;
  bool preserve_localL_p, preserve_localH_p;


  pathL = this->pathL;
  pathH = this->pathH;
  preserve_localL_p = preserve_scored_local_mate_softclip_p(pathL);
  preserve_localH_p = preserve_scored_local_mate_softclip_p(pathH);
  distalL = Path_genomiclow(pathL);
  distalH = Path_genomichigh(pathH);
  if (plusp == true) {
    found_score_ptr_L = found_score_5;
    end_score_ptr_L = end_score_5;
    found_score_ptr_H = found_score_3;
    end_score_ptr_H = end_score_3;
  } else {
    found_score_ptr_L = found_score_3;
    end_score_ptr_L = end_score_3;
    found_score_ptr_H = found_score_5;
    end_score_ptr_H = end_score_5;
  }

  debug3(printf("\nEntering Pathpair_resolve with\n"));
  debug3(Path_print(pathL));
  debug3(Path_print(pathH));
  debug3(printf("insertlength %u\n",Pathpair_insertlength(this)));

  if (Intlist_last_value(pathL->endpoints) == querylengthL) {
    extendL_p = false;
    junctionH = (Junction_T) NULL;
  } else if (preserve_localL_p == true || preserve_localH_p == true) {
    extendL_p = false;
    junctionH = (Junction_T) NULL;
  } else {
    univdiagonal_L = Univcoordlist_last_value(pathL->univdiagonals);
    if (pathL->junctions == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_last_value(pathL->junctions,NULL));
    }
    segmentL_low = univdiagonal_L - querylengthL + Intlist_penultimate_value(pathL->endpoints) + ninserts;
    segmentL_high = univdiagonal_L - querylengthL + Intlist_last_value(pathL->endpoints);

    junctionH = find_overlap_for_pathL(&extendL_p,&cutpointH,segmentL_low,segmentL_high,pathH,querylengthH);
#if defined(DEBUG3) || defined(DEBUG_PAIRED_DEEP)
    printf("Last segment of pathL: %u..%u\n",segmentL_low,segmentL_high);
    if (junctionH != NULL) {
      printf("junctionH: ");
      Junction_print(junctionH);
      printf("\n");
    }
#endif
  }

  if (Intlist_head(pathH->endpoints) == 0) {
    extendH_p = false;
    junctionL = (Junction_T) NULL;
  } else if (preserve_localH_p == true || preserve_localL_p == true) {
    extendH_p = false;
    junctionL = (Junction_T) NULL;
  } else {
    univdiagonal_H = Univcoordlist_head(pathH->univdiagonals);
    segmentH_low = univdiagonal_H - querylengthH + Intlist_head(pathH->endpoints); 
    segmentH_high = univdiagonal_H - querylengthH + Intlist_second_value(pathH->endpoints); 

    junctionL = find_overlap_for_pathH(&extendH_p,&cutpointL,segmentH_low,segmentH_high,pathL,querylengthL);
#if defined(DEBUG3) || defined(DEBUG_PAIRED_DEEP)
    printf("First segment of pathH: %u..%u\n",segmentH_low,segmentH_high);
    if (junctionL != NULL) {
      printf("junctionL: ");
      Junction_print(junctionL);
      printf("\n");
    }
#endif    
  }


  if (extendL_p == false && junctionH == (Junction_T) NULL) {
    /* Skip */

  } else {
    if (junctionH != NULL && univdiagonal_L > cutpointH + (ninserts = Junction_ninserts(junctionH))) {
      /* Have junction and we can reach it */
      debug3(printf("Can extend pathL with junction from %u to %u\n",
		    cutpointH,cutpointH + Junction_adj(junctionH)));
      debug3(Junction_print(junctionH));

      segment_length = univdiagonal_L - cutpointH - ninserts;
      debug3(printf("segment_length on pathL: %u\n",segment_length));

      pathL->univdiagonals = Univcoordlist_reverse(pathL->univdiagonals);
      pathL->univdiagonals = Univcoordlistpool_push(pathL->univdiagonals,univcoordlistpool,
						    univdiagonal_L + Junction_adj(junctionH)
						    univcoordlistpool_trace(__FILE__,__LINE__));
      pathL->univdiagonals = Univcoordlist_reverse(pathL->univdiagonals);

      pathL->endpoints = Intlist_reverse(pathL->endpoints);
      pathL->nmismatches = Intlist_reverse(pathL->nmismatches);
      pathL->ref_nmismatches = Intlist_reverse(pathL->ref_nmismatches);

      /* Tests for endpoints and nmismatches are the same */
      if (Intlist_head(pathL->nmismatches) == 0 &&
	  Intlist_head(pathL->endpoints) >= querylengthL - (int) segment_length - ninserts) {
	/* Can keep nmismatches of 0 */
	Intlist_head_set(pathL->endpoints,querylengthL - (int) segment_length - ninserts); /* Need to adjust for ninserts */

      } else if (Intlist_head(pathL->endpoints) != querylengthL - (int) segment_length - ninserts) {
	Intlist_head_set(pathL->nmismatches,-1);
	Intlist_head_set(pathL->ref_nmismatches,-1);
	Intlist_head_set(pathL->endpoints,querylengthL - (int) segment_length - ninserts); /* Need to adjust for ninserts */
      }

      pathL->endpoints = Intlistpool_push(pathL->endpoints,intlistpool,querylengthL
					  intlistpool_trace(__FILE__,__LINE__));
      pathL->nmismatches = Intlistpool_push(pathL->nmismatches,intlistpool,-1
					    intlistpool_trace(__FILE__,__LINE__));
      pathL->ref_nmismatches = Intlistpool_push(pathL->ref_nmismatches,intlistpool,-1
						intlistpool_trace(__FILE__,__LINE__));
	
      pathL->ref_nmismatches = Intlist_reverse(pathL->ref_nmismatches);
      pathL->nmismatches = Intlist_reverse(pathL->nmismatches);
      pathL->endpoints = Intlist_reverse(pathL->endpoints);
      
      pathL->junctions = List_reverse(pathL->junctions);
      pathL->junctions = Listpool_push(pathL->junctions,listpool,(void *) Junction_copy(junctionH,pathpool)
				       listpool_trace(__FILE__,__LINE__));
      pathL->junctions = List_reverse(pathL->junctions);
      
    } else {
      pathL->endpoints = Intlist_reverse(pathL->endpoints);
      pathL->nmismatches = Intlist_reverse(pathL->nmismatches);
      pathL->ref_nmismatches = Intlist_reverse(pathL->ref_nmismatches);
      
      segmentL_high += querylengthL - Intlist_head(pathL->endpoints); /* If we go to endpoint querylengthL */
      debug3(printf("Can extend pathL with endpoint.  Extend from %u up to %u\n",segmentL_high,cutpointH));

      if (segmentL_high <= cutpointH) {
	Intlist_head_set(pathL->endpoints,querylengthL);
      } else {
	debug3(printf("Cannot extend all the way\n"));
	Intlist_head_set(pathL->endpoints,querylengthL - (segmentL_high - cutpointH));
      }
      Intlist_head_set(pathL->nmismatches,-1);
      Intlist_head_set(pathL->ref_nmismatches,-1);

      pathL->ref_nmismatches = Intlist_reverse(pathL->ref_nmismatches);
      pathL->nmismatches = Intlist_reverse(pathL->nmismatches);
      pathL->endpoints = Intlist_reverse(pathL->endpoints);
    }

    pathL->splice3p = false;
    pathL->splicetype3 = NO_SPLICE;
    pathL->ambig_prob_3 = 0.0;
    Altsplice_free(&pathL->qend_alts,pathpool);
    pathL->qend_alts = (Altsplice_T) NULL;
	
    debug3(printf("Freeing transcripts\n"));
    Transcript_list_gc(&pathL->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&pathL->invalid_transcripts,listpool,transcriptpool);

    debug3(Path_print(pathL));

#if 0
    pathL->nmatches = -1;
    Path_eval_nmatches(&(*found_score_ptr_L),&(*end_score_ptr_L),pathL,queryL_compress_fwd,queryL_compress_rev);
#else
    pathL->nmatches = -1;	/* Calls Path_eval_nmatches */
    Path_trim_chrbounds(&(*found_score_ptr_L),&(*end_score_ptr_L),
			pathL,queryL_compress_fwd,queryL_compress_rev,
			intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif    
  }


  if (extendH_p == false && junctionL == (Junction_T) NULL) {
    /* Skip */

  } else {
    if (junctionL != NULL && cutpointL + querylengthH > univdiagonal_H + (ninserts = Junction_ninserts(junctionL))) {
      /* Have junction and we can reach it */
      debug3(printf("Can extend pathH with junction from %u to %u\n",
		    cutpointL,cutpointL - Junction_adj(junctionL)));
      debug3(Junction_print(junctionL));
      
      segment_length = cutpointL - (univdiagonal_H - querylengthH) - ninserts;
      debug3(printf("segment_length on pathH: %u\n",segment_length));
      
      pathH->univdiagonals = Univcoordlistpool_push(pathH->univdiagonals,univcoordlistpool,
						    univdiagonal_H - Junction_adj(junctionL)
						    univcoordlistpool_trace(__FILE__,__LINE__));
      
      /* Tests for endpoints and nmismatches are different */
      if (Intlist_head(pathH->nmismatches) == 0 &&
	  Intlist_head(pathH->endpoints) <= (int) segment_length + ninserts) {
	/* Can keep nmismatches of 0 */
	
      } else if (Intlist_head(pathH->endpoints) != (int) segment_length + ninserts) {
	Intlist_head_set(pathH->nmismatches,-1);
	Intlist_head_set(pathH->ref_nmismatches,-1);
      }
      Intlist_head_set(pathH->endpoints,(int) segment_length); /* reduced by ninserts above */
      
      pathH->endpoints = Intlistpool_push(pathH->endpoints,intlistpool,0
					  intlistpool_trace(__FILE__,__LINE__));
      pathH->nmismatches = Intlistpool_push(pathH->nmismatches,intlistpool,-1
					    intlistpool_trace(__FILE__,__LINE__));
      pathH->ref_nmismatches = Intlistpool_push(pathH->ref_nmismatches,intlistpool,-1
						intlistpool_trace(__FILE__,__LINE__));
      
      pathH->junctions = Listpool_push(pathH->junctions,listpool,(void *) Junction_copy(junctionL,pathpool)
				       listpool_trace(__FILE__,__LINE__));

    } else {
      segmentH_low -= Intlist_head(pathH->endpoints); /* If we go to endpoint 0 */
      debug3(printf("Can extend pathH with endpoint.  Extend from %u down to %u\n",segmentH_low,cutpointL));

      if (segmentH_low >= cutpointL) {
	Intlist_head_set(pathH->endpoints,0);
      } else {
	debug3(printf("Cannot extend all the way\n"));
	Intlist_head_set(pathH->endpoints,cutpointL - segmentH_low);
      }
      Intlist_head_set(pathH->nmismatches,-1);
      Intlist_head_set(pathH->ref_nmismatches,-1);
    }

    pathH->splice5p = false;
    pathH->splicetype5 = NO_SPLICE;
    pathH->ambig_prob_5 = 0.0;
    Altsplice_free(&pathH->qstart_alts,pathpool);
    pathH->qstart_alts = (Altsplice_T) NULL;
	
    debug3(printf("Freeing transcripts\n"));
    Transcript_list_gc(&pathH->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&pathH->invalid_transcripts,listpool,transcriptpool);

    debug3(Path_print(pathH));

#if 0
    pathH->nmatches = -1;
    Path_eval_nmatches(&(*found_score_ptr_H),&(*end_score_ptr_H),pathH,queryH_compress_fwd,queryH_compress_rev);
#else
    pathH->nmatches = -1;	/* Calls Path_eval_nmatches */
    Path_trim_chrbounds(&(*found_score_ptr_H),&(*end_score_ptr_H),
			pathH,queryH_compress_fwd,queryH_compress_rev,
			intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif
  }


  /* Remove qend_alts if it did not extend far enough (boundedp) or overextends */
  if (pathL->qend_alts != NULL &&
      ((Altsplice_boundedp(pathL->qend_alts) == true &&
	Altsplice_last_univdiagonal(pathL->qend_alts) + pathL->qend_alts->splice_qpos < distalH + querylengthL) ||
       Altsplice_first_univdiagonal(pathL->qend_alts) > distalH + querylengthL)) {
    /* Remove altsplice */
    debug3(printf("Removing altsplice of pathL: boundedp %d, coords %u..%u vs distal %u\n",
		  Altsplice_boundedp(pathL->qend_alts),Altsplice_first_univdiagonal(pathL->qend_alts),
		  Altsplice_last_univdiagonal(pathL->qend_alts),distalH));
    pathL->splice3p = true;
    if (pathL->sensedir == SENSE_FORWARD) {
      if (pathL->plusp == true) {
	pathL->splicetype3 = DONOR;
      } else {
	pathL->splicetype3 = ANTIACCEPTOR;
      }
    } else {
      /* SENSE_ANTI */
      if (pathL->plusp == true) {
	pathL->splicetype3 = ANTIACCEPTOR;
      } else {
	pathL->splicetype3 = DONOR;
      }
    }

    /* assert(Intlist_last_value(pathL->endpoints) == pathL->querylength - pathL->qend_alts->best_distal_length); */
    pathL->ambig_prob_3 = pathL->qend_alts->medial_prob;

    Altsplice_free(&pathL->qend_alts,pathpool);
    pathL->qend_alts = (Altsplice_T) NULL;

    debug3(printf("Freeing transcripts\n"));
    Transcript_list_gc(&pathL->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&pathL->invalid_transcripts,listpool,transcriptpool);
  }


  /* Remove qstart_alts if it did not extend far enough (boundedp) or overextends */
  if (pathH->qstart_alts != NULL &&
      ((Altsplice_boundedp(pathH->qstart_alts) == true &&
	Altsplice_last_univdiagonal(pathH->qstart_alts) + pathH->qstart_alts->splice_qpos > distalL + querylengthH) ||
       Altsplice_first_univdiagonal(pathH->qstart_alts) < distalL + querylengthH)) {
    /* Remove altsplice */
    debug3(printf("Removing altsplice of pathH: boundedp %d, coords %u..%u vs distal %u\n",
		  Altsplice_boundedp(pathH->qstart_alts),Altsplice_first_univdiagonal(pathH->qstart_alts),
		  Altsplice_last_univdiagonal(pathH->qstart_alts),distalL));
    pathH->splice5p = true;
    if (pathH->sensedir == SENSE_FORWARD) {
      if (pathH->plusp == true) {
	pathH->splicetype5 = ACCEPTOR;
      } else {
	pathH->splicetype5 = ANTIDONOR;
      }

    } else {
      /* SENSE_ANTI */
      if (pathH->plusp == true) {
	pathH->splicetype5 = ANTIDONOR;
      } else {
	pathH->splicetype5 = ACCEPTOR;
      }
    }

    /* assert(Intlist_head(pathH->endpoints) == pathH->qstart_alts->best_distal_length); */
    pathH->ambig_prob_5 = pathH->qstart_alts->medial_prob;

    Altsplice_free(&pathH->qstart_alts,pathpool);
    pathH->qstart_alts = (Altsplice_T) NULL;

    debug3(printf("Freeing transcripts\n"));
    Transcript_list_gc(&pathH->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&pathH->invalid_transcripts,listpool,transcriptpool);
  }


  /* Resolve pathL */
  if (pathL->qend_alts != NULL) {
    /* Already have alts to consider */
    debug3(printf("Already have alts at qend of pathL\n"));
    
  } else if (/*qend*/Intlist_last_value(pathL->endpoints) == pathL->querylength ||
	     preserve_localL_p == true || preserve_localH_p == true) {
    /* No need to resolve */
    debug3(printf("No need to resolve qend of pathL\n"));

#if 0
  } else if (0 && Path_genomiclow(pathH) + pathH->querylength < (Univcoord_T) (expected_pairlength + pairlength_delta)) {
    /* Not enough room to perform resolve */
#endif

#if 0
  } else if (pathL->spliceends_qend >= 0 && pathL->spliceends_qend <= Intlist_last_value(pathL->endpoints)) {
    debug3(printf("Already called spliceends at qend\n"));
#endif
    
  } else {
#ifdef USE_EXPECTED
    goal_univdiagonal = Path_genomiclow(pathH) + pathH->querylength + pathL->querylength - expected_pairlength;
    low_univdiagonal = goal_univdiagonal - pairlength_delta;
    high_univdiagonal = goal_univdiagonal + pairlength_delta;
#elif 0
    /* Searches just near the other end */
    high_univdiagonal = Univcoordlist_head(pathH->univdiagonals);
    low_univdiagonal = subtract_bounded(high_univdiagonal,max_insertlength,Univcoordlist_last_value(pathL->univdiagonals));    
#else
    /* Searches the entire space.  Subtract querylengthL so we don't overlap */
    low_univdiagonal = Univcoordlist_last_value(pathL->univdiagonals);
    high_univdiagonal = add_bounded(low_univdiagonal,shortsplicedist,Univcoordlist_head(pathH->univdiagonals) - querylengthL);
#endif
    debug3(printf("Resolving qend of pathL: resolve region is %u..%u (%u..%u)\n",
		  low_univdiagonal,high_univdiagonal,
		  low_univdiagonal - pathH->chroffset,high_univdiagonal - pathH->chroffset));

    Path_qend_resolve(&(*found_score_ptr_L),&(*end_score_ptr_L),
		      pathL,low_univdiagonal,high_univdiagonal,
		      queryptrL,pathL->querylength,
		      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		      stage1L,knownsplicing,stage1L->indelinfo,
		      query_compress_L,queryL_compress_fwd,queryL_compress_rev,genestrand,
		      intlistpool,univcoordlistpool,listpool,pathpool,
		      transcriptpool,univdiagpool,vectorpool);
    debug3(printf("Result of Path_qend_resolve: ")); debug3(Path_print(pathL));
  }


  /* Resolve pathH */
  if (pathH->qstart_alts != NULL) {
    /* Already have alts to consider */
    debug3(printf("Already have alts at qstart of pathH\n"));

  } else if (/*qstart*/Intlist_head(pathH->endpoints) == 0 ||
	     preserve_localH_p == true || preserve_localL_p == true) {
    /* No need to resolve */
    debug3(printf("No need to resolve qstart of pathH\n"));

#if 0
  } else if (0 && Path_genomichigh(pathL) + pathH->querylength + expected_pairlength + pairlength_delta >= genomelength + pathL->querylength) {
    /* Not enough room to perform resolve */
#endif

#if 0
  } else if (pathH->spliceends_qstart >= 0 && pathH->spliceends_qstart >= Intlist_head(pathH->endpoints)) {
    debug3(printf("Already called spliceends at qstart\n"));
#endif

  } else {
#ifdef USE_EXPECTED
    goal_univdiagonal = Path_genomichigh(pathL) - pathL->querylength + expected_pairlength;
    low_univdiagonal = goal_univdiagonal - pairlength_delta;
    high_univdiagonal = goal_univdiagonal + pairlength_delta;
#elif 0
    /* Searches just near the other end */
    low_univdiagonal = Univcoordlist_last_value(pathL->univdiagonals);
    high_univdiagonal = add_bounded(low_univdiagonal,max_insertlength,Univcoordlist_head(pathH->univdiagonals));
#else
    /* Searches the entire space.  Add querylengthH so we don't overlap */
    high_univdiagonal = Univcoordlist_head(pathH->univdiagonals);
    low_univdiagonal = subtract_bounded(high_univdiagonal,shortsplicedist,Univcoordlist_last_value(pathL->univdiagonals) + querylengthH);
#endif
    debug3(printf("Resolving qstart of pathH: resolve region is %u..%u (%u..%u)\n",
		  low_univdiagonal,high_univdiagonal,
		  low_univdiagonal - pathL->chroffset,high_univdiagonal - pathL->chroffset));
      
    Path_qstart_resolve(&(*found_score_ptr_H),&(*end_score_ptr_H),
			pathH,low_univdiagonal,high_univdiagonal,
			queryptrH,pathH->querylength,
			novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			stage1H,knownsplicing,stage1H->indelinfo,
			query_compress_H,queryH_compress_fwd,queryH_compress_rev,genestrand,
			intlistpool,univcoordlistpool,listpool,pathpool,
			transcriptpool,univdiagpool,vectorpool);

    debug3(printf("Result of Path_qstart_resolve: ")); debug3(Path_print(pathH));
  }


  /* Could have original altsplices, or those resulting from resolve procedures above */
  if (pathL->sensedir == SENSE_FORWARD) {
    sense_forward_p = true;
  } else {
    sense_forward_p = false;
  }

  if (pathL->qend_alts != NULL && pathH->qstart_alts != NULL) {
    debug3(printf("Resolving altsplices, case A: resolve both\n"));
    debug3(Path_print(pathL));
    debug3(Path_print(pathH));

    pos5_L = Intlist_penultimate_value(pathL->endpoints); /* Altsplice_qend */
    if (pathL->junctions != NULL) {
      pos5_L += Junction_ninserts((Junction_T) List_last_value(pathL->junctions,NULL));
    }
    pos3_H = Intlist_second_value(pathH->endpoints); /* Altsplice_qstart */

    if (Altsplice_resolve_both(&univdiagonal_L,&splice_qpos_L,&trimpos_L,
			       &medial_nmismatches_L,&distal_nmismatches_L,&medial_prob_L,&distal_prob_L,
			       &univdiagonal_H,&splice_qpos_H,&trimpos_H,
			       &medial_nmismatches_H,&distal_nmismatches_H,&medial_prob_H,&distal_prob_H,
			       pathL->qend_alts,/*anchor_qpos_L*/pos5_L,
			       pathH->qstart_alts,/*anchor_qpos_H*/pos3_H,
			       pathL->querylength,pathH->querylength,
			       query_compress_L,query_compress_H,plusp,genestrand,sense_forward_p) == true) {
      debug3(printf("Case A got univdiagonals %u and %u with trimpos %d and %d\n",
		    univdiagonal_L,univdiagonal_H,trimpos_L,trimpos_H));

      /* Modify pathL */
      pathL = Path_reverse(pathL,/*expect_fwd_p*/false);
      Intlist_head_set(pathL->endpoints,splice_qpos_L);
      pathL->endpoints = Intlistpool_push(pathL->endpoints,intlistpool,trimpos_L  /* was pathL->querylength */
					  intlistpool_trace(__FILE__,__LINE__));
      Intlist_head_set(pathL->nmismatches,medial_nmismatches_L);
      Intlist_head_set(pathL->ref_nmismatches,medial_nmismatches_L);
      pathL->nmismatches = Intlistpool_push(pathL->nmismatches,intlistpool,distal_nmismatches_L
					    intlistpool_trace(__FILE__,__LINE__));
      pathL->ref_nmismatches = Intlistpool_push(pathL->ref_nmismatches,intlistpool,distal_nmismatches_L
						intlistpool_trace(__FILE__,__LINE__));

      splice_distance = univdiagonal_L - Univcoordlist_head(pathL->univdiagonals);
      if (pathL->sensedir == SENSE_FORWARD) {
	sense_forward_p = true;
      } else {
	sense_forward_p = false;
      }
      if (pathL->plusp == sense_forward_p) {
	pathL->junctions =
	  Listpool_push(pathL->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathL->sensedir,
						     /*donor_prob*/medial_prob_L,
						     /*acceptor_prob*/distal_prob_L,pathpool)
			listpool_trace(__FILE__,__LINE__));
      } else {
	pathL->junctions =
	  Listpool_push(pathL->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathL->sensedir,
						     /*donor_prob*/distal_prob_L,
						     /*acceptor_prob*/medial_prob_L,
						     pathpool)
			listpool_trace(__FILE__,__LINE__));
      }

      pathL->univdiagonals = Univcoordlistpool_push(pathL->univdiagonals,univcoordlistpool,univdiagonal_L
						    univcoordlistpool_trace(__FILE__,__LINE__));

      Altsplice_free(&pathL->qend_alts,pathpool);
      pathL->qend_alts = (Altsplice_T) NULL; /* Important if Altsplice_free is an empty procedure */
      pathL->splice3p = false;
      pathL->splicetype3 = NO_SPLICE;
      pathL->ambig_prob_3 = 0.0;

      debug3(printf("Freeing transcripts\n"));
      Transcript_list_gc(&pathL->transcripts,listpool,transcriptpool);
      Transcript_list_gc(&pathL->invalid_transcripts,listpool,transcriptpool);

      pathL = Path_reverse(pathL,/*expect_fwd_p*/true);


      /* Modify pathH */
      Intlist_head_set(pathH->endpoints,splice_qpos_H);
      pathH->endpoints = Intlistpool_push(pathH->endpoints,intlistpool,trimpos_H /* was 0 */
					  intlistpool_trace(__FILE__,__LINE__));
      Intlist_head_set(pathH->nmismatches,medial_nmismatches_H);
      Intlist_head_set(pathH->ref_nmismatches,medial_nmismatches_H);
      pathH->nmismatches = Intlistpool_push(pathH->nmismatches,intlistpool,distal_nmismatches_H
					    intlistpool_trace(__FILE__,__LINE__));
      pathH->ref_nmismatches = Intlistpool_push(pathH->ref_nmismatches,intlistpool,distal_nmismatches_H
						intlistpool_trace(__FILE__,__LINE__));

      splice_distance = Univcoordlist_head(pathH->univdiagonals) - univdiagonal_H;
      if (pathH->sensedir == SENSE_FORWARD) {
	sense_forward_p = true;
      } else {
	sense_forward_p = false;
      }
      if (pathH->plusp == sense_forward_p) {
	pathH->junctions =
	  Listpool_push(pathH->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathH->sensedir,
						     /*donor_prob*/distal_prob_H,
						     /*acceptor_prob*/medial_prob_H,
						     pathpool)
			listpool_trace(__FILE__,__LINE__));
      } else {
	pathH->junctions =
	  Listpool_push(pathH->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathH->sensedir,
						     /*donor_prob*/medial_prob_H,
						     /*acceptor_prob*/distal_prob_H,pathpool)
			listpool_trace(__FILE__,__LINE__));
      }
	
      pathH->univdiagonals = Univcoordlistpool_push(pathH->univdiagonals,univcoordlistpool,univdiagonal_H
						    univcoordlistpool_trace(__FILE__,__LINE__));

      Altsplice_free(&pathH->qstart_alts,pathpool);
      pathH->qstart_alts = (Altsplice_T) NULL; /* Important if Altsplice_free is an empty procedure */
      pathH->splice5p = false;
      pathH->splicetype5 = NO_SPLICE;
      pathH->ambig_prob_5 = 0.0;

      debug3(printf("Freeing transcripts\n"));
      Transcript_list_gc(&pathH->transcripts,listpool,transcriptpool);
      Transcript_list_gc(&pathH->invalid_transcripts,listpool,transcriptpool);

      if (medial_nmismatches_L < 0 || distal_nmismatches_L < 0) {
	found_score_ptr_L = (plusp == true) ? found_score_5 : found_score_3;
	end_score_ptr_L = (plusp == true) ? end_score_5 : end_score_3;
#if 0
	pathL->nmatches = -1;
	Path_eval_nmatches(&(*found_score_ptr_L),&(*end_score_ptr_L),pathL,queryL_compress_fwd,queryL_compress_rev);
#else
	pathL->nmatches = -1;	/* Calls Path_eval_nmatches */
	Path_trim_chrbounds(&(*found_score_ptr_L),&(*end_score_ptr_L),
			    pathL,queryL_compress_fwd,queryL_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif
      }
      if (medial_nmismatches_H < 0 || distal_nmismatches_H < 0) {
	found_score_ptr_H = (plusp == true) ? found_score_3 : found_score_5;
	end_score_ptr_H = (plusp == true) ? end_score_3 : end_score_5;
#if 0
	pathH->nmatches = -1;
	Path_eval_nmatches(&(*found_score_ptr_H),&(*end_score_ptr_H),pathH,queryH_compress_fwd,queryH_compress_rev);
#else
	pathH->nmatches = -1;	/* Calls Path_eval_nmatches */
	Path_trim_chrbounds(&(*found_score_ptr_H),&(*end_score_ptr_H),
			    pathH,queryH_compress_fwd,queryH_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif
      }
    }


  } else if (pathL->qend_alts != NULL) {
    debug3(printf("Resolving altsplices, case B: resolve qend of pathL\n"));
    debug3(Path_print(pathL));

    pos5_L = Intlist_penultimate_value(pathL->endpoints); /* Altsplice_qend */
    if (pathL->junctions != NULL) {
      pos5_L += Junction_ninserts((Junction_T) List_last_value(pathL->junctions,NULL));
    }

    if (Altsplice_resolve_qend(&univdiagonal_L,&splice_qpos_L,&trimpos_L,
			       &medial_nmismatches_L,&distal_nmismatches_L,&medial_prob_L,&distal_prob_L,
			       pathL->qend_alts,/*anchor_qpos_L*/pos5_L,
			       pathL->querylength,pathH->querylength,
			       /*genomiclowH*/Path_genomiclow(pathH),
			       query_compress_L,plusp,genestrand,sense_forward_p) == true) {
      debug3(printf("Case B got univdiagonal %u with trimpos %d\n",univdiagonal_L,trimpos_L));

      /* Modify pathL */
      pathL = Path_reverse(pathL,/*expect_fwd_p*/false);
      Intlist_head_set(pathL->endpoints,splice_qpos_L);
      pathL->endpoints = Intlistpool_push(pathL->endpoints,intlistpool,trimpos_L   /* was pathL->querylength */
					  intlistpool_trace(__FILE__,__LINE__));
      Intlist_head_set(pathL->nmismatches,medial_nmismatches_L);
      Intlist_head_set(pathL->ref_nmismatches,medial_nmismatches_L);
      pathL->nmismatches = Intlistpool_push(pathL->nmismatches,intlistpool,distal_nmismatches_L
					    intlistpool_trace(__FILE__,__LINE__));
      pathL->ref_nmismatches = Intlistpool_push(pathL->ref_nmismatches,intlistpool,distal_nmismatches_L
						intlistpool_trace(__FILE__,__LINE__));

      splice_distance = univdiagonal_L - Univcoordlist_head(pathL->univdiagonals);
      if (pathL->sensedir == SENSE_FORWARD) {
	sense_forward_p = true;
      } else {
	sense_forward_p = false;
      }
      if (pathL->plusp == sense_forward_p) {
	pathL->junctions =
	  Listpool_push(pathL->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathL->sensedir,
						     /*donor_prob*/medial_prob_L,
						     /*acceptor_prob*/distal_prob_L,pathpool)
			listpool_trace(__FILE__,__LINE__));
      } else {
	pathL->junctions =
	  Listpool_push(pathL->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathL->sensedir,
						     /*donor_prob*/distal_prob_L,
						     /*acceptor_prob*/medial_prob_L,
						     pathpool)
			listpool_trace(__FILE__,__LINE__));
      }

      pathL->univdiagonals = Univcoordlistpool_push(pathL->univdiagonals,univcoordlistpool,univdiagonal_L
						    univcoordlistpool_trace(__FILE__,__LINE__));

      Altsplice_free(&pathL->qend_alts,pathpool);
      pathL->qend_alts = (Altsplice_T) NULL; /* Important if Altsplice_free is an empty procedure */
      pathL->splice3p = false;
      pathL->splicetype3 = NO_SPLICE;
      pathL->ambig_prob_3 = 0.0;

      debug3(printf("Freeing transcripts\n"));
      Transcript_list_gc(&pathL->transcripts,listpool,transcriptpool);
      Transcript_list_gc(&pathL->invalid_transcripts,listpool,transcriptpool);

      pathL = Path_reverse(pathL,/*expect_fwd_p*/true);

      if (medial_nmismatches_L < 0 || distal_nmismatches_L < 0) {
	found_score_ptr_L = (plusp == true) ? found_score_5 : found_score_3;
	end_score_ptr_L = (plusp == true) ? end_score_5 : end_score_3;
#if 0
	pathL->nmatches = -1;
	Path_eval_nmatches(&(*found_score_ptr_L),&(*end_score_ptr_L),pathL,queryL_compress_fwd,queryL_compress_rev);
#else
	pathL->nmatches = -1;	/* Calls Path_eval_nmatches */
	Path_trim_chrbounds(&(*found_score_ptr_L),&(*end_score_ptr_L),
			    pathL,queryL_compress_fwd,queryL_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif
      }
    }
    
  } else if (pathH->qstart_alts != NULL) {
    debug3(printf("Resolving altsplices, case C: resolve qstart of pathH\n"));
    debug3(Path_print(pathH));

    pos3_H = Intlist_second_value(pathH->endpoints); /* Altsplice_qstart */

    if (Altsplice_resolve_qstart(&univdiagonal_H,&splice_qpos_H,&trimpos_H,
				 &medial_nmismatches_H,&distal_nmismatches_H,
				 &medial_prob_H,&distal_prob_H,
				 pathH->qstart_alts,/*anchor_qpos_H*/pos3_H,
				 pathL->querylength,pathH->querylength,
				 /*genomichighL*/Path_genomichigh(pathL),
				 query_compress_H,plusp,genestrand,sense_forward_p) == true) {
      debug3(printf("Case C got univdiagonal %u with trimpos %d\n",univdiagonal_H,trimpos_H));

      /* Modify pathH */
      Intlist_head_set(pathH->endpoints,splice_qpos_H);
      pathH->endpoints = Intlistpool_push(pathH->endpoints,intlistpool,trimpos_H /* was 0 */
					  intlistpool_trace(__FILE__,__LINE__));
      Intlist_head_set(pathH->nmismatches,medial_nmismatches_H);
      Intlist_head_set(pathH->ref_nmismatches,medial_nmismatches_H);
      pathH->nmismatches = Intlistpool_push(pathH->nmismatches,intlistpool,distal_nmismatches_H
					    intlistpool_trace(__FILE__,__LINE__));
      pathH->ref_nmismatches = Intlistpool_push(pathH->ref_nmismatches,intlistpool,distal_nmismatches_H
						intlistpool_trace(__FILE__,__LINE__));

      splice_distance = Univcoordlist_head(pathH->univdiagonals) - univdiagonal_H;
      if (pathH->sensedir == SENSE_FORWARD) {
	sense_forward_p = true;
      } else {
	sense_forward_p = false;
      }
      if (pathH->plusp == sense_forward_p) {
	pathH->junctions =
	  Listpool_push(pathH->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathH->sensedir,
						     /*donor_prob*/distal_prob_H,
						     /*acceptor_prob*/medial_prob_H,
						     pathpool)
			listpool_trace(__FILE__,__LINE__));
      } else {
	pathH->junctions =
	  Listpool_push(pathH->junctions,listpool,
			(void *) Junction_new_splice(splice_distance,pathH->sensedir,
						     /*donor_prob*/medial_prob_H,
						     /*acceptor_prob*/distal_prob_H,pathpool)
			listpool_trace(__FILE__,__LINE__));
      }
	
      pathH->univdiagonals = Univcoordlistpool_push(pathH->univdiagonals,univcoordlistpool,univdiagonal_H
						    univcoordlistpool_trace(__FILE__,__LINE__));

      Altsplice_free(&pathH->qstart_alts,pathpool);
      pathH->qstart_alts = (Altsplice_T) NULL; /* Important if Altsplice_free is an empty procedure */
      pathH->splice5p = false;
      pathH->splicetype5 = NO_SPLICE;
      pathH->ambig_prob_5 = 0.0;

      debug3(printf("Freeing transcripts\n"));
      Transcript_list_gc(&pathH->transcripts,listpool,transcriptpool);
      Transcript_list_gc(&pathH->invalid_transcripts,listpool,transcriptpool);

      if (medial_nmismatches_H < 0 || distal_nmismatches_H < 0) {
	found_score_ptr_H = (plusp == true) ? found_score_3 : found_score_5;
	end_score_ptr_H = (plusp == true) ? end_score_3 : end_score_5;
#if 0
	pathH->nmatches = -1;
	Path_eval_nmatches(&(*found_score_ptr_H),&(*end_score_ptr_H),pathH,queryH_compress_fwd,queryH_compress_rev);
#else
	pathH->nmatches = -1;	/* Calls Path_eval_nmatches */
	Path_trim_chrbounds(&(*found_score_ptr_H),&(*end_score_ptr_H),
			    pathH,queryH_compress_fwd,queryH_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif
      }
    }

  } else {
    debug3(printf("No altsplices\n"));
  }

  this->insertlength = compute_insertlength(&this->pair_relationship,pathL,pathH,
					    queryseqL,queryseqH,this->plusp);

  debug3(printf("Exiting Pathpair_resolve with\n"));
  debug3(Path_print(pathL));
  debug3(Path_print(pathH));
  debug3(printf("insertlength %u\n\n",Pathpair_insertlength(this)));

  /* outerlength should remain the same */
  /* this->outerlength = Path_genomichigh(pathH) - Path_genomiclow(pathL); */

  return;
}



int
Pathpair_interval_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  Univcoord_T genomiclow_x, genomiclow_y,
    genomichigh_x, genomichigh_y;

  genomiclow_x = Path_genomiclow(a->pathL);
  genomiclow_y = Path_genomiclow(b->pathL);

  if (genomiclow_x < genomiclow_y) {
    return -1;
  } else if (genomiclow_y < genomiclow_x) {
    return +1;
  } else {
    genomichigh_x = Path_genomichigh(a->pathH);
    genomichigh_y = Path_genomichigh(b->pathH);
    if (genomichigh_x > genomichigh_y) {
      return -1;
    } else if (genomichigh_y > genomichigh_x) {
      return +1;
    } else {
      return 0;
    }
  }
}


int
Pathpair_structure_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  /* Does not consider main_univdiagonals, which means that the choice
     of main_univdiagonal is arbitrary */

  if ((cmp = Path_structure_cmp(&(a->pathL),&(b->pathL))) != 0) {
    return cmp;
  } else if ((cmp = Path_structure_cmp(&(a->pathH),&(b->pathH))) != 0) {
    return cmp;
  } else {
    return 0;
  }
}


int
Pathpair_structure_ignore_sense_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  /* Does not consider main_univdiagonals, which means that the choice
     of main_univdiagonal is arbitrary */

  if ((cmp = Path_structure_ignore_sense_cmp(&(a->pathL),&(b->pathL))) != 0) {
    return cmp;
  } else if ((cmp = Path_structure_ignore_sense_cmp(&(a->pathH),&(b->pathH))) != 0) {
    return cmp;
  } else {
    return 0;
  }
}


bool
Pathpair_overlap_p (T x, T y) {
  Univcoord_T genomiclow_x, genomiclow_y,
    genomichigh_x, genomichigh_y;

  genomiclow_x = Path_genomiclow(x->pathL);
  genomiclow_y = Path_genomiclow(y->pathL);
  genomichigh_x = Path_genomichigh(x->pathH);
  genomichigh_y = Path_genomichigh(y->pathH);
  
  if (genomichigh_x < genomiclow_y) {
    return false;
  } else if (genomichigh_y < genomiclow_x) {
    return false;
  } else {
    return true;
  }
}


#if 0
List_T
Pathpair_filter (List_T pathpairs, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		 Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		 Hitlistpool_T hitlistpool) {
  List_T filtered, p;
  /* Chrpos_T min_insertlength, insertlength; */
  double max_splice_prob = 0.0, splice_prob;
  int max_nmatches = 0, nmatches;
  T pathpair;


  debug0(printf("Entered Pathpair_filter with %d pathpairs\n",List_length(pathpairs)));

  /* Causes problems because it can lead to bad splices */
  /* Need to have thresholds on splice probabilities */
  /* 1.  Filter by splice prob */
  filtered = (List_T) NULL;
  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (T) List_head(p);
    if ((splice_prob = pathpair->pathL->total_splice_prob + pathpair->pathH->total_splice_prob) > max_splice_prob) {
      max_splice_prob = splice_prob;
    }
    if ((nmatches = pathpair->pathL->nmatches + pathpair->pathH->nmatches) > max_nmatches) {
      max_nmatches = nmatches;
    }
  }
  
  debug0(printf("Have max_splice_prob %f and max_nmatches %d\n",max_splice_prob,max_nmatches));
  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (T) List_head(p);
#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
    printf("This path has %f+%f splice prob\n",
	   pathpair->pathL->total_splice_prob,pathpair->pathH->total_splice_prob);
    printf("This path has %d+%d nmatches\n",
	   pathpair->pathL->nmatches,pathpair->pathH->nmatches);
    Path_print(pathpair->pathL);
    Path_print(pathpair->pathH);
#endif

    if (pathpair->pathL->total_splice_prob + pathpair->pathH->total_splice_prob == max_splice_prob ||
	pathpair->pathL->nmatches + pathpair->pathH->nmatches == max_nmatches) {
      filtered = Hitlist_push(filtered,hitlistpool,(void *) pathpair
			      hitlistpool_trace(__FILE__,__LINE__));
    } else {
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
  }

  Hitlistpool_free_list(&pathpairs,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));

#ifdef FILTER_BY_INSERTLENGTH
  /* Causes problems because it does not account for nmatches */
  /* Speed up is achieved instead by limiting hits per rare in Intersect_concordance */
  /* 2.  Filter by insertlength */
  if (List_length(filtered) > 1) {
    pathpairs = filtered;
    filtered = (List_T) NULL;

    pathpair = (T) List_head(pathpairs);
    min_insertlength = pathpair->pathH->main_univdiagonal - pathpair->pathL->main_univdiagonal;
    for (p = List_next(pathpairs); p != NULL; p = List_next(p)) {
      pathpair = (T) List_head(p);
      if ((insertlength = pathpair->pathH->main_univdiagonal - pathpair->pathL->main_univdiagonal) < min_insertlength) {
	min_insertlength = insertlength;
      }
    }
  
    debug0(printf("Have min_insertlength %u\n",min_insertlength));
    for (p = pathpairs; p != NULL; p = List_next(p)) {
      pathpair = (T) List_head(p);
#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
      printf("This path has %u insertlength\n",
	     pathpair->pathH->main_univdiagonal - pathpair->pathL->main_univdiagonal);
      Path_print(pathpair->pathL);
      Path_print(pathpair->pathH);
#endif

      insertlength = pathpair->pathH->main_univdiagonal - pathpair->pathL->main_univdiagonal;
      if (insertlength <= min_insertlength * INSERTLENGTH_FACTOR) {
	filtered = Hitlist_push(filtered,hitlistpool,(void *) pathpair
				hitlistpool_trace(__FILE__,__LINE__));
      } else {
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
    }

    Hitlistpool_free_list(&pathpairs,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }
#endif

  debug0(printf("Exiting Pathpair_filter with %d pathpairs\n",List_length(filtered)));

  return filtered;
}
#endif


#ifdef TO_FIX
/* Returns true if ilengths are valid */
static bool
find_ilengths (int *ilength_low, int *ilength_high, T path, Univcoord_T common_genomicpos) {
  List_T p, q;
  Substring_T substring;
  Junction_T junction;


  debug15(printf("Finding ilengths for common_genomicpos %u\n",(Chrpos_T) (common_genomicpos - chroffset)));
  if (path->plusp == true) {
#ifdef DEBUG15
    printf("plus.  Checking common genomicpos %llu against\n",common_genomicpos - path->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u, trimmed %d..%d\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset,
	     Substring_querystart_trimmed(substring),Substring_queryend_trimmed(substring));
    }
    printf("\n");
#endif
    /* Plus: Subtract 1 from alignend */
    *ilength_low = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring %p: %u..%u, trimmed %d..%d\n",substring,
		     Substring_alignstart_trim(substring) - hit->chroffset,
		     Substring_alignend_trim(substring) - 1U - hit->chroffset,
		     Substring_querystart_trimmed(substring),Substring_queryend_trimmed(substring)));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_low += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_low += (common_genomicpos - Substring_alignstart_trim(substring) + 1);
	*ilength_high = ((Substring_alignend_trim(substring) - 1) - common_genomicpos + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_high += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Plus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  } else {
#ifdef DEBUG15
    printf("minus.  Checking common genomicpos %llu against\n",common_genomicpos - hit->chroffset);
    for (p = hit->substrings_1toN; p != NULL; p = List_next(p)) {
      substring = (Substring_T) List_head(p);
      printf("substring %p: %u..%u, trimmed %d..%d\n",
	     substring,Substring_alignstart_trim(substring) - hit->chroffset,
	     Substring_alignend_trim(substring) - 1U - hit->chroffset,
	     Substring_querystart_trimmed(substring),Substring_queryend_trimmed(substring));
    }
    printf("\n");
#endif
    /* Minus: Subtract 1 from alignstart */
    *ilength_high = 0;
    for (p = hit->substrings_1toN, q = hit->junctions_1toN; p != NULL; p = List_next(p), q = List_next(q)) {
      substring = (Substring_T) List_head(p);
      debug15(printf("substring: %u..%u\n",
		     Substring_alignstart_trim(substring) - 1U - hit->chroffset,
		     Substring_alignend_trim(substring) - hit->chroffset));
      if (Substring_overlap_point_trimmed_p(substring,common_genomicpos) == false) {
	*ilength_high += Substring_genomic_alignment_length(substring);
	if (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_high += Junction_nindels(junction);
	  }
	}

      } else {
	*ilength_high += ((Substring_alignstart_trim(substring) - 1) - common_genomicpos + 1);
	*ilength_low = (common_genomicpos - (Substring_alignend_trim(substring) /*+ 1*/) + 1);
	p = List_next(p);
	while (p != NULL) {
	  substring = (Substring_T) List_head(p);
	  *ilength_low += Substring_genomic_alignment_length(substring);
	  p = List_next(p);
	}
	while (q != NULL) {
	  junction = (Junction_T) List_head(q);
	  if (Junction_type(junction) == INS_JUNCTION) {
	    *ilength_low += Junction_nindels(junction);
	  }
	  q = List_next(q);
	}
	debug15(printf("Minus: Have ilength_low %d and ilength_high %d\n",*ilength_low,*ilength_high));
	return true;
      }
    }
  }

  return false;
}
#endif


#ifdef TO_FIX
/* Needed to compute overlap properly.  Based on pair_insert_length below, plus code for handling GMAP. */
static Univcoord_T
pair_common_genomicpos (Stage3end_T hit5, Stage3end_T hit3) {
  Univcoord_T common_genomicpos;
  Univcoord_T start5, end5, start3, end3;
  List_T p, q;
  Substring_T substring, substring5, substring3;

  if (hit5->plusp == true && hit3->plusp == true) {
    /* plus/plus */
    debug15(printf("Computing overlap using substrings plus/plus\n"));

    start5 = hit5->genomiclow + hit5->querystart_trimmed + start_amb_length(hit5);
    end5 = (hit5->genomichigh - 1) - (hit5->querylength - hit5->queryend_trimmed) - end_amb_length(hit5);
    start3 = hit3->genomiclow + hit3->querystart_trimmed + start_amb_length(hit3);
    end3 = (hit3->genomichigh - 1) - (hit3->querylength - hit3->queryend_trimmed) - end_amb_length(hit3);
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 < start5) {
      /* Case 1 */
      return false;
    } else if (end5 < start3) {
      /* Case 6 */
      return false;
    } else if (start3 < start5) {
      if (end3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus/plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus/plus case 3\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus/plus case 4\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus/plus case 5a\n"));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus/plus case 5b\n"));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus/plus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else if (hit5->plusp == true && hit3->plusp == false) {
    /* plus/minus */
    debug15(printf("Computing overlap using substrings plus/minus\n"));
    return 0;

#if 0
    start5 = hit5->genomiclow + hit5->querystart_trimmed + start_amb_length(hit5);
    end5 = hit5->genomichigh - (hit5->querylength - hit5->queryend_trimmed) - end_amb_length(hit5);
    start3 = hit3->genomiclow - hit3->querystart_trimmed - start_amb_length(hit3);
    end3 = hit3->genomichigh + (hit3->querylength - hit3->queryend_trimmed) + end_amb_length(hit3);

    if (start3 < start5) {
      /* Case 1 */
      return 0;
    } else if (end5 < end3) {
      /* Case 6 */
      return 0;
    } else if (end3 < start5) {
      if (start3 < end5) {
	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2a: start5 %u\n",start5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	}

	/* Case 2: Tails overlap.  Go from start5 to start3 */
	debug15(printf("plus case 2b: start3 %u\n",start3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (start3 < end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0U;
#endif

  } else if (hit5->plusp == false && hit3->plusp == true) {
    /* minus/plus */
    debug15(printf("Computing overlap using substrings minus/plus\n"));
    return 0;

#if 0
    start5 = hit5->genomiclow - hit5->querystart_trimmed - start_amb_length(hit5);
    end5 = hit5->genomichigh + (hit5->querylength - hit5->queryend_trimmed) + end_amb_length(hit5);
    start3 = hit3->genomiclow + hit3->querystart_trimmed + start_amb_length(hit3);
    end3 = hit3->genomichigh - (hit3->querylength - hit3->queryend_trimmed) - end_amb_length(hit3);

    if (end3 < end5) {
      /* Case 1 */
      return 0;
    } else if (start5 < start3) {
      /* Case 6 */
      return 0;
    } else if (start3 < end5) {
      if (end3 < start5) {
	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2a: end5 %u\n",end5 - hit5->chroffset));
	if (Substring_overlap_point_trimmed_p(hit3->substring0,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,end5)) {
	  return end5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring2,end5)) {
	  return end5;
	}

	/* Case 2: Tails overlap.  Go from end5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	if (Substring_overlap_point_trimmed_p(hit5->substring2,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,end3)) {
	  return end3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring0,end3)) {
	  return end3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("plus case 3\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }

    } else {
      if (end3 < start5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("plus case 4\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("plus case 5a\n"));
	if (Substring_overlap_point_trimmed_p(hit5->substring0,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring1,start3)) {
	  return start3;
	} else if (Substring_overlap_point_trimmed_p(hit5->substring2,start3)) {
	  return start3;
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("plus case 5b\n"));
	if (Substring_overlap_point_trimmed_p(hit3->substring2,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring1,start5)) {
	  return start5;
	} else if (Substring_overlap_point_trimmed_p(hit3->substring0,start5)) {
	  return start5;
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("plus general: hit3->substring1\n"));
    if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring2 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring1)) != 0) {
      return common_genomicpos;
    } else if (hit5->substring0 != NULL &&
	       (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring1)) != 0) {
      return common_genomicpos;
    }

    if (hit3->substring2 != NULL) {
      debug15(printf("plus general: hit3->substring2\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring2)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring2)) != 0) {
	return common_genomicpos;
      }
    }

    if (hit3->substring0 != NULL) {
      debug15(printf("plus general: hit3->substring0\n"));
      if ((common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring1,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring2 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring2,hit3->substring0)) != 0) {
	return common_genomicpos;
      } else if (hit5->substring0 != NULL &&
		 (common_genomicpos = Substring_overlap_segment_trimmed(hit5->substring0,hit3->substring0)) != 0) {
	return common_genomicpos;
      }
    }

    return 0;
#endif

  } else if (hit5->plusp == false && hit3->plusp == false) {
    /* minus/minus */
    debug15(printf("Computing overlap using substrings minus/minus\n"));

    start5 = (hit5->genomiclow - 1) - hit5->querystart_trimmed /*- start_amb_length(hit5)*/;
    end5 = hit5->genomichigh + (hit5->querylength - hit5->queryend_trimmed) /*+ end_amb_length(hit5)*/;
    start3 = (hit3->genomiclow - 1) - hit3->querystart_trimmed /*- start_amb_length(hit3)*/;
    end3 = hit3->genomichigh + (hit3->querylength - hit3->queryend_trimmed) /*+ end_amb_length(hit3)*/;
    debug15(printf("hit5 endpoints are %u..%u.  hit3 endpoints are %u..%u\n",
		   start5-hit5->chroffset,end5-hit5->chroffset,start3-hit3->chroffset,end3-hit3->chroffset));

    if (end3 > start5) {
      /* Case 1 */
      return 0;
    } else if (end5 > start3) {
      /* Case 6 */
      return 0;
    } else if (start3 > start5) {
      if (end3 > end5) {
	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("minus/minus case 2a: start5 %llu (%u)\n",start5,start5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start5)) {
	    return start5;
	  }
	}

	/* Case 2: Tails overlap.  Go from start5 to end3 */
	debug15(printf("plus case 2b: end3 %u\n",end3 - hit3->chroffset));
	for (p = hit5->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end3)) {
	    return end3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 3: hit3 subsumes hit5 */
	debug15(printf("minus/minus case 3: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}

	/* Fall through to general algorithm */
      }

    } else {
      if (end3 > end5) {
	/* Case 4: hit5 subsumes hit3 */
	debug15(printf("minus/minus case 4: start3 %u\n",(Chrpos_T) (start3 - hit3->chroffset)));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}
	/* Fall through to general algorithm */

      } else {
	/* Case 5: Based on hit3_trimmed_length */
	debug15(printf("minus case 5a: start3 %u\n",start3 - hit3->chroffset));
	for (p = hit5->substrings_1toN; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,start3)) {
	    return start3;
	  }
	}

	/* Case 5: Based on hit5_trimmed_length */
	debug15(printf("minus case 5b: end5 %u\n",end5 - hit5->chroffset));
	for (p = hit3->substrings_Nto1; p != NULL; p = List_next(p)) {
	  substring = (Substring_T) List_head(p);
	  if (Substring_overlap_point_trimmed_p(substring,end5)) {
	    return end5;
	  }
	}
	/* Fall through to general algorithm */
      }
    }

    /* General algorithm */
    debug15(printf("minus/minus general\n"));
    for (p = hit3->substrings_1toN; p != NULL; p = List_next(p)) {
      substring3 = (Substring_T) List_head(p);
      for (q = hit5->substrings_1toN; q != NULL; q = List_next(q)) {
	substring5 = (Substring_T) List_head(q);
	if ((common_genomicpos = Substring_overlap_segment_trimmed(substring5,substring3)) != 0) {
	  return common_genomicpos;
	}
      }
    }

    return 0;

  } else {
    abort();
    return 0;
  }
}
#endif


#ifdef TO_FIX
/* Note: Do not alter this->insertlength, which is used for SAM
   output.  The insertlength computed here is used only for performing
   --clip-overlap or --merge-overlap */
int
Pathpair_overlap (int *hardclip5_low, int *hardclip5_high, int *hardclip3_low, int *hardclip3_high, T this) {
  Pair_T path5, path3;
  int clipdir;
  int ilength53, ilength35, ilength5_low, ilength5_high, ilength3_low, ilength3_high;
  int common_shift, common_left, common_right;
  Univcoord_T common_genomicpos, common_genomicpos_right, common_genomicpos_left;
  int shift_right, shift_left;
#ifdef DEBUG15
  int overlap;
#endif


  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;

  path5 = this->path5;
  path3 = this->path3;

  debug15(printf("Entered Stage3pair_overlap with hittype %s and %s\n",
		 hittype_string(path5->hittype),hittype_string(path3->hittype)));
  if (path5->hittype == SAMECHR_SPLICE || path5->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (path3->hittype == SAMECHR_SPLICE || path3->hittype == TRANSLOC_SPLICE) {
    return 0;
  } else if (path5->plusp != path3->plusp) {
    debug15(printf("The two ends are not on the same strand, so returning 0\n"));
    return 0;
  } else {
    debug15(printf("path5 querystart_trimmed %d + amb_start %d, queryend_trimmed %d + amb_end %d, path3 querystart_trimmed %d + amb_start %d, queryend_trimmed %d + amb_end %d\n",
		   path5->querystart_trimmed,start_amb_length(path5),path5->queryend_trimmed,end_amb_length(path5),
		   path3->querystart_trimmed,start_amb_length(path3),path3->queryend_trimmed,end_amb_length(path3)));
    if (path5->plusp == true) {
      /* plus */
#if 0
      path5_trimmed_length = path5->querylength - path5->querystart_trimmed - path5->queryend_trimmed - start_amb_length(path5) - end_amb_length(path5);
      path3_trimmed_length = path3->querylength - path3->querystart_trimmed - path3->queryend_trimmed - start_amb_length(path3) - end_amb_length(path3);
      totallength = path5_trimmed_length + path3_trimmed_length;
      debug15(printf("totallength = %d, path5 trimmed length = %d, path3 trimmed length = %d\n",
		     totallength,path5_trimmed_length,path3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,path5->querystart_trimmed + start_amb_length(path5),
		     path5->queryend_trimmed + end_amb_length(path5),path3->querystart_trimmed + start_amb_length(path3),
		     path3->queryend_trimmed + end_amb_length(path3)));
#endif

      if ((common_genomicpos = pair_common_genomicpos(path5,path3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53 is %d, ilength 35 is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_low > 0);
	  assert(ilength3_low > 0);
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_high > 0);
	  assert(ilength3_high > 0);
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));


	if ((ilength53 = ilength5_low + ilength3_high) >= (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >=, not >, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("plus, ilength53 is longer.  Clipping heads.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = +1;

	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
	  *hardclip3_low += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
	  debug15(printf("Ambig clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
		*hardclip3_low += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - path3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_low > 0);
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_high > 0);
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    debug15(printf("Initial computation of clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_high += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
	    *hardclip3_low += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
	    debug15(printf("Recomputed clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("plus, ilength35 is longer.  Clipping tails.\n"));
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 common_left+common_right-1,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
	  *hardclip3_high += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
	  debug15(printf("Ambig clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
		*hardclip3_high += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - path3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_low > 0);
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_high > 0);
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    debug15(printf("Initial computation of clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_low += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
	    *hardclip3_high += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
	    debug15(printf("Recomputed clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 plus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}

	debug15(printf("returning clipdir %d\n",clipdir));
	return clipdir;
      }

    } else {
      /* minus */
#if 0
      path5_trimmed_length = path5->querylength - path5->querystart_trimmed - (path5->querylength - path5->queryend_trimmed) - start_amb_length(path5) - end_amb_length(path5);
      path3_trimmed_length = path3->querylength - path3->querystart_trimmed - (path3->querylength - path3->queryend_trimmed) - start_amb_length(path3) - end_amb_length(path3);
      totallength = path5_trimmed_length + path3_trimmed_length;
      debug15(printf("totallength = %d, path5 trimmed length = %d, path3 trimmed length = %d\n",
		     totallength,path5_trimmed_length,path3_trimmed_length));
      debug15(printf("original insertlength: %d, trim+amb5: %d..%d, trim+amb3: %d..%d\n",
		     this->insertlength,path5->querystart_trimmed + start_amb_length(path5),
		     path5->queryend_trimmed + path5->end_amb_length,path3->querystart_trimmed + start_amb_length(path3),
		     path3->queryend_trimmed + path3->end_amb_length));
#endif

      if ((common_genomicpos = pair_common_genomicpos(path5,path3)) == 0) {
	debug15(printf("Cannot determine a common point, so returning 0\n"));
	return 0;

      } else if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		 find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	debug15(printf("Cannot determine ilengths, so returning 0\n"));
	return 0;

      } else {
	debug15(printf("Inclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	debug15(printf("ilength53lh is %d, ilength35lh is %d\n",ilength5_low + ilength3_high - 1,ilength3_low + ilength5_high - 1));

	common_left = (ilength5_low < ilength3_low) ? ilength5_low : ilength3_low;
	common_right = (ilength5_high < ilength3_high) ? ilength5_high : ilength3_high;
	if (common_right > common_left) {
	  common_shift = common_right/2 - (common_left - 1)/2;
	  debug15(printf("Common shift is %d = common_right %d/2 - (common_left %d - 1)/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_low > 0);
	  assert(ilength3_low > 0);
	  ilength5_low -= 1;
	  ilength3_low -= 1;
	} else {
	  common_shift = (common_right - 1)/2 - common_left/2;
	  debug15(printf("Common shift is %d = (common_right %d - 1)/2 - common_left %d/2\n",
			 common_shift,common_right,common_left));
	  assert(ilength5_high > 0);
	  assert(ilength3_high > 0);
	  ilength5_high -= 1;
	  ilength3_high -= 1;
	}
	debug15(printf("Exclusive: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	if ((ilength53 = ilength5_low + ilength3_high) > (ilength35 = ilength3_low + ilength5_high)) {
	  /* Use >, not >=, so we favor clipping heads over clipping tails in case of a tie */
	  debug15(printf("minus, ilength53 is longer.  Clipping tails.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = +1;


	  /* Want to clip 5 high and 3 low */
	  *hardclip5_high = ilength5_high - common_shift;
	  *hardclip3_low = ilength3_low + common_shift;
	  debug15(printf("Overlap clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_high += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
	  *hardclip3_low += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
	  debug15(printf("Ambig clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_high = ilength5_high /*- common_shift*/;
		*hardclip3_low = ilength3_low /*+ common_shift*/;
		*hardclip5_high += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
		*hardclip3_low += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - path3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength3_low > ilength5_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_low > 0);
	      ilength3_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_high > 0);
	      ilength5_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_high = ilength5_high /*- common_shift*/;
	    *hardclip3_low = ilength3_low /*+ common_shift*/;
	    debug15(printf("Initial computation of clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_high += path5->querystart_trimmed /*+ start_amb_length(path5)*/;
	    *hardclip3_low += (path3->querylength - path3->queryend_trimmed) /*+ end_amb_length(path3)*/;
	    debug15(printf("Recomputed clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_high < 0) {
	    *hardclip5_high = 0;
	  }
	  if (*hardclip3_low < 0) {
	    *hardclip3_low = 0;
	  }
	  debug15(printf("Positive clip for ilength53 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif

	} else {
	  debug15(printf("minus, ilength35 is longer.  Clipping heads.\n"));
	  debug15(overlap = common_left + common_right - 1);
	  debug15(printf("Overlap is %d = common_left %d + common_right %d - 1\n",
			 overlap,common_left,common_right));
	  clipdir = -1;

	  /* Want to clip 5 low and 3 high */
	  *hardclip5_low = ilength5_low + common_shift;
	  *hardclip3_high = ilength3_high - common_shift;
	  debug15(printf("Overlap clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  *hardclip5_low += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
	  *hardclip3_high += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
	  debug15(printf("Ambig clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	  if (common_shift != 0) {
	    if (test_hardclips(&common_genomicpos,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset) == true) {
	      /* No adjustment needed, but need to recompute ilengths for shifted common_genomicpos */
	    } else {
	      common_genomicpos_right = adjust_hardclips_right(&shift_right,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset);
	      common_genomicpos_left = adjust_hardclips_left(&shift_left,*hardclip5_low,path5,*hardclip3_high,path3,path3->chroffset);
	      debug15(printf("shift_right %d, shift_left %d\n",shift_right,shift_left));
	      if (shift_right == 0 && shift_left == 0) {
		/* Try original position without a shift */
		*hardclip5_low = ilength5_low /*+ common_shift*/;
		*hardclip3_high = ilength3_high /*- common_shift*/;
		*hardclip5_low += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
		*hardclip3_high += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
		if (test_hardclips(&common_genomicpos,*hardclip3_low,path3,*hardclip5_high,path5,path3->chroffset) == false) {
		  *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
		  return 0;
		}
	      } else if (shift_left == 0) {
		common_genomicpos = common_genomicpos_right;
	      } else if (shift_right == 0) {
		common_genomicpos = common_genomicpos_left;
	      } else if (shift_right <= shift_left) {
		common_genomicpos = common_genomicpos_right;
	      } else {
		common_genomicpos = common_genomicpos_left;
	      }
	    }

	    debug15(printf("New common point is %u\n",common_genomicpos - path3->chroffset));
	    /* Recompute hardclips */
	    if (find_ilengths(&ilength5_low,&ilength5_high,path5,common_genomicpos) == false ||
		find_ilengths(&ilength3_low,&ilength3_high,path3,common_genomicpos) == false) {
	      *hardclip5_low = *hardclip5_high = *hardclip3_low = *hardclip3_high = 0;
	      return 0;
	    } else if (ilength5_low > ilength3_high) {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength5_low > 0);
	      ilength5_low -= 1;
	    } else {
	      debug15(printf("Uneven: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));
	      assert(ilength3_high > 0);
	      ilength3_high -= 1;
	    }
	    debug15(printf("Even: ilengths5: %d|%d.  ilengths3: %d|%d\n",ilength5_low,ilength5_high,ilength3_low,ilength3_high));

	    *hardclip5_low = ilength5_low /*+ common_shift*/;
	    *hardclip3_high = ilength3_high /*- common_shift*/;
	    debug15(printf("Initial computation of clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));

	    *hardclip5_low += (path5->querylength - path5->queryend_trimmed) /*+ end_amb_length(path5)*/;
	    *hardclip3_high += path3->querystart_trimmed /*+ start_amb_length(path3)*/;
	    debug15(printf("Recomputed clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			   *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
	  }

#if 0
	  if (*hardclip5_low < 0) {
	    *hardclip5_low = 0;
	  }
	  if (*hardclip3_high < 0) {
	    *hardclip3_high = 0;
	  }
	  debug15(printf("Positive clip for ilength35 minus is hardclip5 %d..%d and hardclip3 %d..%d\n",
			 *hardclip5_low,*hardclip5_high,*hardclip3_low,*hardclip3_high));
#endif
	}
      }

      debug15(printf("returning clipdir %d\n",clipdir));
      return clipdir;
    }
  }
}
#endif



static void
print_query_header (Filestring_T fp, char initchar, Shortread_T queryseq, bool invertp) {
  FPRINTF(fp,"%c",initchar);
  if (invertp == false) {
    Shortread_print_oneline(fp,queryseq);
  } else {
    Shortread_print_oneline_revcomp(fp,queryseq);
  }

  return;
}



static void
print_barcode_and_quality (Filestring_T fp, Shortread_T queryseq, bool invertp, int quality_shift) {
  char *barcode;

  if ((barcode = Shortread_barcode(queryseq)) != NULL) {
    FPRINTF(fp,"\tbarcode:%s",barcode);
  }

  if (Shortread_quality_string(queryseq) != NULL) {
    FPRINTF(fp,"\t");
    if (invertp == false) {
      Shortread_print_quality(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
			      quality_shift,/*show_chopped_p*/true);
    } else {
      Shortread_print_quality_revcomp(fp,queryseq,/*hardclip_low*/0,/*hardclip_high*/0,
				      quality_shift,/*show_chopped_p*/true);
    }
  }

  return;
}

bool
Pathpair_print_end_alignment (Filestring_T fp, Result_T result, Resulttype_T resulttype,
			      char initchar, bool firstp, 
			      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
			      int maxpaths, bool quiet_if_excessive_p, bool invertp, int quality_shift,
			      Listpool_T listpool) {
  bool printp = false, excessivep, translocationp, concordant_softclipped_p;
  Pathpair_T *pathpairarray, pathpair;
  Path_T *patharray, path5, path3, path;
  int npaths_primary, npaths_altloc, pathnum;
  int first_absmq, second_absmq;

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (only_concordant_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NM);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t0 %s",UNPAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
    
      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == CONCORDANT_UNIQ) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    pathpair = pathpairarray[0];
    path5 = pathpair->path5;
    path3 = pathpair->path3;

    if (Path_softclippedp(path5) == true || Path_softclippedp(path3) == true) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/true,OUTPUT_CU);
    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CU);
    }

    if (omit_concordant_uniq_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

    } else {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t1 %s",CONCORDANT_TEXT);
    
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
    
      printp = true;
      if (firstp == true) {
	Path_print_alignment(fp,path5,pathpair,queryseq,invertp,listpool);
      } else {
	Path_print_alignment(fp,path3,pathpair,queryseq,invertp,listpool);
      }

      FPRINTF(fp,"\n");
    }

  } else if (resulttype == CONCORDANT_MULT) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (omit_concordant_mult_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CX);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);

      /* No further output */
      FPRINTF(fp,"\n");
      printp = false;

    } else {
      concordant_softclipped_p = false;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	if (Path_softclippedp(path5) == true || Path_softclippedp(path3) == true) {
	  concordant_softclipped_p = true;
	}
      }
      Filestring_set_split_output(fp,concordant_softclipped_p,OUTPUT_CM);
      
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
	
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
	
      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);
      
      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	
	if (firstp == true) {
	  Path_print_alignment(fp,path5,pathpair,queryseq,invertp,listpool);
	} else {
	  Path_print_alignment(fp,path3,pathpair,queryseq,invertp,listpool);
	}
      }

      FPRINTF(fp,"\n");
    }

  } else if (only_concordant_p == true) {
    /* Skip printing */
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

  } else if (resulttype == CONCORDANT_TRANSLOC) {
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CT);
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      /* No xs category for transloc, so ignore quiet-if-excessive_p */
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
      FPRINTF(fp," (transloc)");
	
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);

      /* No further output */
      FPRINTF(fp,"\n");

    } else {
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,CONCORDANT_TEXT);
      FPRINTF(fp," (transloc)");

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);
	
      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);

      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	
	if (firstp == true) {
	  Path_print_alignment(fp,path5,pathpair,queryseq,invertp,listpool);
	} else {
	  Path_print_alignment(fp,path3,pathpair,queryseq,invertp,listpool);
	}
      }
      
      FPRINTF(fp,"\n");
    }

  } else if (resulttype == PAIRED_UNIQ_INV || resulttype == PAIRED_UNIQ_SCR || resulttype == PAIRED_UNIQ_TOOLONG) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    pathpair = pathpairarray[0];
    path5 = pathpair->path5;
    path3 = pathpair->path3;

    if (resulttype == PAIRED_UNIQ_INV) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PI);
    } else if (resulttype == PAIRED_UNIQ_SCR) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PS);
    } else if (resulttype == PAIRED_UNIQ_TOOLONG) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PL);
    } else {
      fprintf(stderr,"Unexpected resulttype %d\n",resulttype);
      abort();
    }
    
    print_query_header(fp,initchar,queryseq,invertp);
    FPRINTF(fp,"\t1 %s",PAIRED_TEXT);

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    FPRINTF(fp,"\t");
    Shortread_print_header(fp,headerseq1,headerseq2);

    printp = true;
    if (firstp == true) {
      Path_print_alignment(fp,path5,pathpair,queryseq,invertp,listpool);
    } else {
      Path_print_alignment(fp,path3,pathpair,queryseq,invertp,listpool);
    }

    FPRINTF(fp,"\n");

  } else if (resulttype == PAIRED_MULT) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PX);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,PAIRED_TEXT);
	
      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);

      /* No further output */
      FPRINTF(fp,"\n");

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PM);
      print_query_header(fp,initchar,queryseq,invertp);
      FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,PAIRED_TEXT);

      print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

      FPRINTF(fp,"\t");
      Shortread_print_header(fp,headerseq1,headerseq2);

      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;

	if (firstp == true) {
	  Path_print_alignment(fp,path5,pathpair,queryseq,invertp,listpool);
	} else {
	  Path_print_alignment(fp,path3,pathpair,queryseq,invertp,listpool);
	}
      }

      FPRINTF(fp,"\n");
    }
    
  } else {
    /* Print as singles */
    if (firstp == true) {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      /* patharray_mate = (Path_T *) Result_array2(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result); */
      patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    } else {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      /* patharray_mate = (Path_T *) Result_array(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result); */
      patharray = (Path_T *) Result_array2(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    }

    excessivep = false;
    translocationp = false;
    if (resulttype == HALFMAPPING_UNIQ) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HU);

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HT);
      translocationp = true;

    } else if (resulttype == HALFMAPPING_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HM);
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);

    } else if (resulttype == UNPAIRED_TRANSLOC) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);
      translocationp = true;

    } else if (resulttype == UNPAIRED_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      }

    } else {
      fprintf(stderr,"Resulttype is %s\n",Resulttype_string(resulttype));
      abort();
    }

    print_query_header(fp,initchar,queryseq,invertp);
    FPRINTF(fp,"\t%d %s",npaths_primary + npaths_altloc,UNPAIRED_TEXT);
    if (translocationp == true) {
      FPRINTF(fp," (transloc)");
    }

    print_barcode_and_quality(fp,queryseq,invertp,quality_shift);

    FPRINTF(fp,"\t");
    Shortread_print_header(fp,headerseq1,headerseq2);

    if (excessivep == true) {
      /* No output */
					      
    } else {
      printp = true;
      if (firstp == true) {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  path = patharray[pathnum-1];
	  Path_print_alignment(fp,path,/*pathpair*/NULL,queryseq,invertp,listpool);
	}
      } else {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  path = patharray[pathnum-1];
	  Path_print_alignment(fp,path,/*pathpair*/NULL,queryseq,invertp,listpool);
	}
      }
    }

    FPRINTF(fp,"\n");
  }

  return printp;
}


bool
Pathpair_print_end_m8 (Filestring_T fp, Result_T result, Resulttype_T resulttype,
		       bool firstp, char *accession,
		       int maxpaths, bool quiet_if_excessive_p, bool invertp,
		       Listpool_T listpool) {
  bool printp = false, excessivep, concordant_softclipped_p;
  Pathpair_T *pathpairarray, pathpair;
  Path_T *patharray, path5, path3, path;
  int npaths_primary, npaths_altloc, pathnum;
  int first_absmq, second_absmq;

  if (resulttype == PAIREDEND_NOMAPPING) {
    if (only_concordant_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);
    }

  } else if (resulttype == CONCORDANT_UNIQ) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    pathpair = pathpairarray[0];
    path5 = pathpair->path5;
    path3 = pathpair->path3;

    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CU);

    if (omit_concordant_uniq_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

    } else {
      printp = true;
      if (firstp == true) {
	Path_print_m8(fp,path5,accession,/*acc_suffix*/"/1",invertp,listpool);
      } else {
	Path_print_m8(fp,path3,accession,/*acc_suffix*/"/2",invertp,listpool);
      }
    }

  } else if (resulttype == CONCORDANT_MULT) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (omit_concordant_mult_p == true) {
      /* Skip printing */
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

    } else if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CX);

    } else {
      concordant_softclipped_p = false;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	if (Path_softclippedp(path5) == true || Path_softclippedp(path3) == true) {
	  concordant_softclipped_p = true;
	}
      }
      Filestring_set_split_output(fp,concordant_softclipped_p,OUTPUT_CM);
      
      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	
	if (firstp == true) {
	  Path_print_m8(fp,path5,accession,/*acc_suffix*/"/1",invertp,listpool);
	} else {
	  Path_print_m8(fp,path3,accession,/*acc_suffix*/"/2",invertp,listpool);
	}
      }
    }

  } else if (only_concordant_p == true) {
    /* Skip printing */
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_NONE);

  } else if (resulttype == CONCORDANT_TRANSLOC) {
    Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_CT);
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      /* Skip */
    } else {
      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;
	
	if (firstp == true) {
	  Path_print_m8(fp,path5,accession,/*acc_suffix*/"/1",invertp,listpool);
	} else {
	  Path_print_m8(fp,path3,accession,/*acc_suffix*/"/2",invertp,listpool);
	}
      }
    }

  } else if (resulttype == PAIRED_UNIQ_INV || resulttype == PAIRED_UNIQ_SCR || resulttype == PAIRED_UNIQ_TOOLONG) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    pathpair = pathpairarray[0];
    path5 = pathpair->path5;
    path3 = pathpair->path3;

    if (resulttype == PAIRED_UNIQ_INV) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PI);
    } else if (resulttype == PAIRED_UNIQ_SCR) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PS);
    } else if (resulttype == PAIRED_UNIQ_TOOLONG) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PL);
    } else {
      fprintf(stderr,"Unexpected resulttype %d\n",resulttype);
      abort();
    }
    
    printp = true;
    if (firstp == true) {
      Path_print_m8(fp,path5,accession,/*acc_suffix*/"/1",invertp,listpool);
    } else {
      Path_print_m8(fp,path3,accession,/*acc_suffix*/"/2",invertp,listpool);
    }

  } else if (resulttype == PAIRED_MULT) {
    pathpairarray = (Pathpair_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);

    if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PX);

    } else {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_PM);

      printp = true;
      for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	pathpair = pathpairarray[pathnum-1];
	path5 = pathpair->path5;
	path3 = pathpair->path3;

	if (firstp == true) {
	  Path_print_m8(fp,path5,accession,/*acc_suffix*/"/1",invertp,listpool);
	} else {
	  Path_print_m8(fp,path3,accession,/*acc_suffix*/"/2",invertp,listpool);
	}
      }
    }
    
  } else {
    /* Print as singles */
    if (firstp == true) {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      /* patharray_mate = (Path_T *) Result_array2(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result); */
      patharray = (Path_T *) Result_array(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    } else {
      /* Get stage3array_mate first to avoid incorrect values for npaths */
      /* patharray_mate = (Path_T *) Result_array(&npaths_mate_primary,&npaths_mate_altloc,&first_absmq,&second_absmq,result); */
      patharray = (Path_T *) Result_array2(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,result);
    }

    excessivep = false;
    /* translocationp = false; */
    if (resulttype == HALFMAPPING_UNIQ) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HU);

    } else if (resulttype == HALFMAPPING_TRANSLOC) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HT);
      /* translocationp = true; */

    } else if (resulttype == HALFMAPPING_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_HM);
      }

    } else if (resulttype == UNPAIRED_UNIQ) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UU);

    } else if (resulttype == UNPAIRED_TRANSLOC) {
      Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UT);
      /* translocationp = true; */

    } else if (resulttype == UNPAIRED_MULT) {
      if (quiet_if_excessive_p && npaths_primary + npaths_altloc > maxpaths) {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UX);
	excessivep = true;
      } else {
	Filestring_set_split_output(fp,/*concordant_softclipped_p*/false,OUTPUT_UM);
      }

    } else {
      fprintf(stderr,"Resulttype is %s\n",Resulttype_string(resulttype));
      abort();
    }

    if (excessivep == true) {
      /* No output */
					      
    } else {
      printp = true;
      if (firstp == true) {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  path = patharray[pathnum-1];
	  Path_print_m8(fp,path,accession,/*acc_suffix*/"/1",invertp,listpool);
	}
      } else {
	for (pathnum = 1; pathnum <= npaths_primary + npaths_altloc && pathnum <= maxpaths; pathnum++) {
	  path = patharray[pathnum-1];
	  Path_print_m8(fp,path,accession,/*acc_suffix*/"/2",invertp,listpool);
	}
      }
    }
  }

  return printp;
}



void
Pathpair_setup (Outputtype_T output_type_in, Univcoord_T genomelength_in,
		Chrpos_T shortsplicedist_in) {

  output_type = output_type_in;
  genomelength = genomelength_in;
  shortsplicedist = shortsplicedist_in;

  return;
}

#if 0
void
Pathpair_pass2_setup (int expected_pairlength, int pairlength_deviation) {
  int max_insertlength;

  /* expected_pairlength = expected_pairlength_in; */
  /* pairlength_delta = 3*pairlength_deviation_in; */

  max_insertlength = expected_pairlength + 2*pairlength_deviation;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  return;
}
#endif
