static char rcsid[] = "$Id: splicetrie.c 223511 2020-11-14 15:50:08Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "splicetrie.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"
#include "univcoord.h"

#include "iitdef.h"
#include "interval.h"
#include "splicetrie_build.h"	/* For single_leaf_p and multiple_leaf_p macros */
#include "dynprog_end.h"



/* Finding short-overlap splicing.  Also may want to turn on DEBUG4H in stage1hr.c.  No longer working; try debug2 */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Finding short exons and short overlap.  Also may want to turn on DEBUG4K in stage1hr.c. */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Dump coords.  May want to turn on DEBUG7 in dynprog.c. */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Solve ends */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif


#define clear_start(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))


static Univcoord_T *splicesites;
static Genomecomp_T *splicefrags_ref;
static Genomecomp_T *splicefrags_alt;

static Trieoffset_T *trieoffsets_obs;
static Triecontent_T *triecontents_obs;
static Trieoffset_T *trieoffsets_max;
static Triecontent_T *triecontents_max;

static bool snpp;
static bool amb_closest_p;
static int min_shortend;
static bool amb_clip_p;


void
Splicetrie_setup (Univcoord_T *splicesites_in, Genomecomp_T *splicefrags_ref_in, Genomecomp_T *splicefrags_alt_in,
		  Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		  Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
		  bool snpp_in, bool amb_closest_p_in, bool amb_clip_p_in, int min_shortend_in) {

  splicesites = splicesites_in;
  splicefrags_ref = splicefrags_ref_in;
  splicefrags_alt = splicefrags_alt_in;

  trieoffsets_obs = trieoffsets_obs_in;
  triecontents_obs = triecontents_obs_in;
  trieoffsets_max = trieoffsets_max_in;
  triecontents_max = triecontents_max_in;

  snpp = snpp_in;
  amb_closest_p = amb_closest_p_in;
  amb_clip_p = amb_clip_p_in;
  min_shortend = min_shortend_in;

  return;
}



/************************************************************************
 *   Using splicetries
 ************************************************************************/

#ifdef USE_2BYTE_RELOFFSETS
static void
get_offsets (int *offseta, int *offsetc, int *offsetg, int *offsett,
	     Trieoffset_T offsets1, Trieoffset_T offsets2) {

  *offsetc = (int) (offsets1 & 0xffff);
  *offseta = (int) ((offsets1 >>= 16) & 0xffff);

  *offsett = (int) (offsets2 & 0xffff);
  *offsetg = (int) ((offsets2 >>= 16) & 0xffff);

  return;
}
#endif


/* Modified from Splicetrie_dump in splicetrie_build.c */
static int
splicetrie_size (Triecontent_T *triestart) {
  int size;
  Triecontent_T leaf;
  int nleaves;
  int offseta, offsetc, offsetg, offsett;

  if (single_leaf_p(leaf = triestart[0])) {
    return 1;

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    return nleaves;

  } else {
#ifdef USE_2BYTE_RELOFFSETS
    get_offsets(&offseta,&offsetc,&offsetg,&offsett,triestart[1],triestart[2]);
#else
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];
#endif

    size = 0;
    if (offseta > 0) {
      size += splicetrie_size(&(triestart[-offseta]));
    }
    if (offsetc > 0) {
      size += splicetrie_size(&(triestart[-offsetc]));
    }
    if (offsetg > 0) {
      size += splicetrie_size(&(triestart[-offsetg]));
    }
    if (offsett > 0) {
      size += splicetrie_size(&(triestart[-offsett]));
    }

    return size;
  }
}



static List_T
solve_end5_aux (Univcoord_T **coordsptr, Univcoord_T *coords,
		List_T best_pairs, Triecontent_T *triestart,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,
		       
		int *finalscore, int *nmatches, int *nmismatches,
		int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		int splicelength, int contlength, Splicetype_T far_splicetype,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int *dynprogindex, Dynprog_T dynprog, 
		char *revsequence1, char *revsequenceuc1,
		int length1, int length2, int revoffset1, int revoffset2,
		int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		int extraband_end, double defect_rate) {
  Triecontent_T leaf;
  int nleaves, i;
  Univcoord_T splicecoord;
  Chrpos_T shortest_intron_length = -1U, intron_length;
  int offseta, offsetc, offsetg, offsett;
  int spliceoffset2_anchor, spliceoffset2_far;

  int score, miss_score, nmatches0, nmismatches0, nopens0, nindels0;
  int length_distal;
  List_T pairs;
    
  if (single_leaf_p(leaf = triestart[0])) {
    splicecoord = splicesites[leaf];
    debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
    if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	Dynprog_make_splicejunction_5(splicejunction,splicejunction_alt,splicecoord,
				      splicelength,contlength,far_splicetype,genome,genomealt,watsonp) == true) {
      debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
      debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		    length1,length2,chroffset,splicecoord));
      if (watsonp) {
	spliceoffset2_anchor = revoffset2;
	spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
      } else {
	spliceoffset2_anchor = revoffset2;
	spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
      }
      if ((pairs = Dynprog_end5_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
					       &nopens0,&nindels0,dynprog,revsequence1,revsequenceuc1,
					       /*revsequence2*/&(splicejunction[length2-1]),
					       /*revsequencealt2*/&(splicejunction_alt[length2-1]),
					       length1,length2,revoffset1,spliceoffset2_anchor,spliceoffset2_far,
					       chroffset,chrhigh,watsonp,genestrand,jump_late_p,
					       genome,genomealt,pairpool,
					       extraband_end,defect_rate,contlength)) != NULL) {

	/* miss_score = perfect_score - score; */
	assert(miss_score <= 0);
	debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
		      score,perfect_score,score-perfect_score,miss_score));
	debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	  debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	  /* Just use results from Dynprog_end5_splicejunction */
	  pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,
					     watsonp,pairpool);
#else
	  length_distal = length1 - contlength;
#endif
	  best_pairs = pairs;
	  *finalscore = score;
	  *nmatches = nmatches0;
	  *nmismatches = nmismatches0;
	  *nopens = nopens0;
	  *nindels = nindels0;
	  *knownsplicep = true;
	  *ambig_end_length = length_distal;
	  *threshold_miss_score = miss_score - obsmax_penalty;
	  shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	  *coordsptr = coords;
	  *(*coordsptr)++ = splicecoord;

	} else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		   && Univcoordlist_find(*coords,splicecoord) == false
#endif
		   ) {
	  if (amb_closest_p == false) {
	    debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	    /* best_pairs = (List_T) NULL; */
	    *(*coordsptr)++ = splicecoord;
	  } else {
	    intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    if (intron_length > shortest_intron_length) {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	    } else {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
	      /* Just use results from Dynprog_end5_splicejunction */
	      pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,
						 watsonp,pairpool);
#else
	      length_distal = length1 - contlength;
#endif
	      best_pairs = pairs;
	      *finalscore = score;
	      *nmatches = nmatches0;
	      *nmismatches = nmismatches0;
	      *nopens = nopens0;
	      *nindels = nindels0;
	      *knownsplicep = true;
	      *ambig_end_length = length_distal;
	      *threshold_miss_score = miss_score - obsmax_penalty;
	      shortest_intron_length = intron_length;
	      *coordsptr = coords;
	      *(*coordsptr)++ = splicecoord;
	    }
	  }
	}
      }
    }
    debug7(printf("\n"));

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      splicecoord = splicesites[leaf];
      debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
      if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	  Dynprog_make_splicejunction_5(splicejunction,splicejunction_alt,splicecoord,
					splicelength,contlength,far_splicetype,genome,genomealt,watsonp) == true) {
	debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
	debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		      length1,length2,chroffset,splicecoord));
	if (watsonp) {
	  spliceoffset2_anchor = revoffset2;
	  spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
	} else {
	  spliceoffset2_anchor = revoffset2;
	  spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
	}
	if ((pairs = Dynprog_end5_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
						 &nopens0,&nindels0,dynprog,revsequence1,revsequenceuc1,
						 /*revsequence2*/&(splicejunction[length2-1]),
						 /*revsequencealt2*/&(splicejunction_alt[length2-1]),
						 length1,length2,revoffset1,spliceoffset2_anchor,spliceoffset2_far,
						 chroffset,chrhigh,watsonp,genestrand,jump_late_p,
						 genome,genomealt,pairpool,
						 extraband_end,defect_rate,contlength)) != NULL) {

	  /* miss_score = perfect_score - score; */
	  assert(miss_score <= 0);
	  debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
			score,perfect_score,score-perfect_score,miss_score));
	  debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	  if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	    debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	    /* Just use results from Dynprog_end5_splicejunction */
	    pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	    length_distal = length1 - contlength;
#endif
	    best_pairs = pairs;
	    *finalscore = score;
	    *nmatches = nmatches0;
	    *nmismatches = nmismatches0;
	    *nopens = nopens0;
	    *nindels = nindels0;
	    *knownsplicep = true;
	    *ambig_end_length = length_distal;
	    *threshold_miss_score = miss_score - obsmax_penalty;
	    shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    *coordsptr = coords;
	    *(*coordsptr)++ = splicecoord;

	  } else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		     /* Filter for duplicates later */
		     && Univcoordlist_find(*coords,splicecoord) == false
#endif
		     ) {
	    if (amb_closest_p == false) {
	      debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	      /* best_pairs = (List_T) NULL; */
	      *(*coordsptr)++ = splicecoord;
	    } else {
	      intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	      if (intron_length > shortest_intron_length) {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	      } else {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
		/* Just use results from Dynprog_end5_splicejunction */
		pairs = Dynprog_add_known_splice_5(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
		length_distal = length1 - contlength;
#endif
		best_pairs = pairs;
		*finalscore = score;
		*nmatches = nmatches0;
		*nmismatches = nmismatches0;
		*nopens = nopens0;
		*nindels = nindels0;
		*knownsplicep = true;
		*ambig_end_length = length_distal;
		*threshold_miss_score = miss_score - obsmax_penalty;
		shortest_intron_length = intron_length;
		*coordsptr = coords;
		*(*coordsptr)++ = splicecoord;
	      }
	    }
	  }
	}
      }
      debug7(printf("\n"));

    }

  } else {
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];

    if (offseta > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offseta]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsetc > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetc]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsetg > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetg]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsett > 0) {
      best_pairs = solve_end5_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsett]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
  }

  return best_pairs;
}


static List_T
solve_end3_aux (Univcoord_T **coordsptr, Univcoord_T *coords,
		List_T best_pairs, Triecontent_T *triestart,
		Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		int *finalscore, int *nmatches, int *nmismatches,
		int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		int splicelength, int contlength, Splicetype_T far_splicetype,
		Univcoord_T chroffset, Univcoord_T chrhigh,
		int *dynprogindex, Dynprog_T dynprog, 
		char *sequence1, char *sequenceuc1,
		int length1, int length2, int offset1, int offset2,
		int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		int extraband_end, double defect_rate) {
  Triecontent_T leaf;
  int nleaves, i;
  Univcoord_T splicecoord;
  Chrpos_T shortest_intron_length = -1U, intron_length;
  int offseta, offsetc, offsetg, offsett;
  int spliceoffset2_anchor, spliceoffset2_far;

  int score, miss_score, nmatches0, nmismatches0, nopens0, nindels0;
  int length_distal;
  List_T pairs;
    
  if (single_leaf_p(leaf = triestart[0])) {
    splicecoord = splicesites[leaf];
    debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
    if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	Dynprog_make_splicejunction_3(splicejunction,splicejunction_alt,splicecoord,
				      splicelength,contlength,far_splicetype,genome,genomealt,watsonp) == true) {
      debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
      debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		    length1,length2,chroffset,splicecoord));
      if (watsonp) {
	spliceoffset2_anchor = offset2;
	spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
      } else {
	spliceoffset2_anchor = offset2;
	spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
      }
      if ((pairs = Dynprog_end3_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
					       &nopens0,&nindels0,dynprog,sequence1,sequenceuc1,
					       /*sequence2*/splicejunction,/*sequencealt2*/splicejunction_alt,
					       length1,length2,offset1,spliceoffset2_anchor,spliceoffset2_far,
					       chroffset,chrhigh,watsonp,genestrand,jump_late_p,
					       genome,genomealt,pairpool,
					       extraband_end,defect_rate,contlength)) != NULL) {

	/* miss_score = perfect_score - score; */
	assert(miss_score <= 0);
	debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
		      score,perfect_score,score-perfect_score,miss_score));
	debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	  debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	  /* Just results of Dynprog_end3_splicejunction */
	  pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	  length_distal = length1 - contlength;
#endif
	  best_pairs = pairs;
	  *finalscore = score;
	  *nmatches = nmatches0;
	  *nmismatches = nmismatches0;
	  *nopens = nopens0;
	  *nindels = nindels0;
	  *knownsplicep = true;
	  *ambig_end_length = length_distal;
	  *threshold_miss_score = miss_score - obsmax_penalty;
	  shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	  *coordsptr = coords;
	  *(*coordsptr)++ = splicecoord;

	} else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		   && Univcoordlist_find(*coords,splicecoord) == false
#endif
		   ) {
	  if (amb_closest_p == false) {
	    debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	    /* best_pairs = (List_T) NULL; */
	    *(*coordsptr)++ = splicecoord;
	  } else {
	    intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    if (intron_length > shortest_intron_length) {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	    } else {
	      debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			    miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
	      /* Just use results of Dynprog_end3_splicejunction */
	      pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	      length_distal = length1 - contlength;
#endif
	      best_pairs = pairs;
	      *finalscore = score;
	      *nmatches = nmatches0;
	      *nmismatches = nmismatches0;
	      *nopens = nopens0;
	      *nindels = nindels0;
	      *knownsplicep = true;
	      *ambig_end_length = length_distal;
	      *threshold_miss_score = miss_score - obsmax_penalty;
	      shortest_intron_length = intron_length;
	      *coordsptr = coords;
	      *(*coordsptr)++ = splicecoord;
	    }
	  }
	}
      }
    }
    debug7(printf("\n"));

  } else if (multiple_leaf_p(leaf)) {
    nleaves = (int) (-leaf);
    for (i = 1; i <= nleaves; i++) {
      leaf = triestart[i];
      splicecoord = splicesites[leaf];
      debug7(printf("Checking leaf %d at %u: ",(int) leaf,splicecoord));
      if (splicecoord >= knownsplice_limit_low && splicecoord <= knownsplice_limit_high &&
	  Dynprog_make_splicejunction_3(splicejunction,splicejunction_alt,splicecoord,
					splicelength,contlength,far_splicetype,genome,genomealt,watsonp) == true) {
	debug7(printf("intron length %d, ",splicecoord - anchor_splicesite));
	debug7(printf("length1 = %d, length2 = %d, chroffset = %u, splicecoord = %u\n",
		      length1,length2,chroffset,splicecoord));
	if (watsonp) {
	  spliceoffset2_anchor = offset2;
	  spliceoffset2_far = spliceoffset2_anchor - anchor_splicesite + splicecoord;
	} else {
	  spliceoffset2_anchor = offset2;
	  spliceoffset2_far = spliceoffset2_anchor + anchor_splicesite - splicecoord;
	}
	if ((pairs = Dynprog_end3_splicejunction(&(*dynprogindex),&score,&miss_score,&nmatches0,&nmismatches0,
						 &nopens0,&nindels0,dynprog,sequence1,sequenceuc1,
						 /*sequence2*/splicejunction,/*sequencealt2*/splicejunction_alt,
						 length1,length2,offset1,spliceoffset2_anchor,spliceoffset2_far,
						 chroffset,chrhigh,watsonp,genestrand,jump_late_p,
						 genome,genomealt,pairpool,
						 extraband_end,defect_rate,contlength)) != NULL) {

	  /* miss_score = perfect_score - score; */
	  assert(miss_score <= 0);
	  debug7(printf("score %d - perfect score %d = miss %d expected vs %d returned.  ",
			score,perfect_score,score-perfect_score,miss_score));
	  debug7(printf("  comparing against threshold_miss %d + obsmax_penalty %d\n",*threshold_miss_score,obsmax_penalty));
	  if (score > 0 && miss_score > *threshold_miss_score + obsmax_penalty) {
	    debug7(printf("miss %d > threshold %d + %d",miss_score,*threshold_miss_score,obsmax_penalty));
#if 0
	    /* Just use results of Dynprog_end3_splicejunction */
	    pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
	    length_distal = length1 - contlength;
#endif
	    best_pairs = pairs;
	    *finalscore = score;
	    *nmatches = nmatches0;
	    *nmismatches = nmismatches0;
	    *nopens = nopens0;
	    *nindels = nindels0;
	    *knownsplicep = true;
	    *ambig_end_length = length_distal;
	    *threshold_miss_score = miss_score - obsmax_penalty;
	    shortest_intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	    *coordsptr = coords;
	    *(*coordsptr)++ = splicecoord;
	  } else if (miss_score == *threshold_miss_score + obsmax_penalty
#if 0
		     && Univcoordlist_find(*coords,splicecoord) == false
#endif
		     ) {
	    if (amb_closest_p == false) {
	      debug7(printf("miss %d == threshold %d + %d, so ambiguous",miss_score,*threshold_miss_score,obsmax_penalty));
	      /* best_pairs = (List_T) NULL; */
	      *(*coordsptr)++ = splicecoord;
	    } else {
	      intron_length = (splicecoord > anchor_splicesite) ? splicecoord - anchor_splicesite : anchor_splicesite - splicecoord;
	      if (intron_length > shortest_intron_length) {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d > shortest %d, so ignore",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
	      } else {
		debug7(printf("miss %d == threshold %d + %d, but intron_length %d < shortest %d, so new best",
			      miss_score,*threshold_miss_score,obsmax_penalty,intron_length,shortest_intron_length));
#if 0
		/* Just use results of Dynprog_end3_splicejunction */
		pairs = Dynprog_add_known_splice_3(&length_distal,pairs,anchor_splicesite,splicecoord,chroffset,watsonp,pairpool);
#else
		length_distal = length1 - contlength;
#endif
		best_pairs = pairs;
		*finalscore = score;
		*nmatches = nmatches0;
		*nmismatches = nmismatches0;
		*nopens = nopens0;
		*nindels = nindels0;
		*knownsplicep = true;
		*ambig_end_length = length_distal;
		*threshold_miss_score = miss_score - obsmax_penalty;
		shortest_intron_length = intron_length;
		*coordsptr = coords;
		*(*coordsptr)++ = splicecoord;
	      }
	    }
	  }
	}
      }
      debug7(printf("\n"));

    }

  } else {
    offseta = (int) triestart[1];
    offsetc = (int) triestart[2];
    offsetg = (int) triestart[3];
    offsett = (int) triestart[4];

    if (offseta > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offseta]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsetc > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetc]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsetg > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsetg]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
    if (offsett > 0) {
      best_pairs = solve_end3_aux(&(*coordsptr),coords,best_pairs,&(triestart[-offsett]),
				  knownsplice_limit_low,knownsplice_limit_high,
				  &(*finalscore),&(*nmatches),&(*nmismatches),
				  &(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				  &(*threshold_miss_score),obsmax_penalty,perfect_score,
				  anchor_splicesite,splicejunction,splicejunction_alt,
				  splicelength,contlength,far_splicetype,
				  chroffset,chrhigh,&(*dynprogindex),dynprog,
				  sequence1,sequenceuc1,length1,length2,offset1,offset2,
				  cdna_direction,watsonp,genestrand,jump_late_p,
				  genome,genomealt,pairpool,extraband_end,defect_rate);
    }
  }

  return best_pairs;
}


List_T
Splicetrie_solve_end5 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *revsequence1, char *revsequenceuc1,
		       int length1, int length2, int revoffset1, int revoffset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		       int extraband_end, double defect_rate) {
  Univcoord_T *coordsptr, *coords, splicecoord0;
  unsigned int *triestart;
  int size, i;

  debug7(printf("Entering Splicetrie_solve_end5 with limits %u..%u, anchor splicesite %u (%u)\n",
		knownsplice_limit_low,knownsplice_limit_high,anchor_splicesite,anchor_splicesite - chroffset));
  if (trieoffsets[j] == NULL_POINTER) {
    return best_pairs;
  } else {
    triestart = &(triecontents[trieoffsets[j]]);
  }

  if ((size = splicetrie_size(triestart)) == 0) {
    return best_pairs;
  } else {
    coordsptr = coords = (Univcoord_T *) MALLOCA(size * sizeof(Univcoord_T));

    best_pairs = solve_end5_aux(&coordsptr,coords,best_pairs,triestart,
				knownsplice_limit_low,knownsplice_limit_high,
				&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				&(*threshold_miss_score),obsmax_penalty,perfect_score,
				anchor_splicesite,splicejunction,splicejunction_alt,
				splicelength,contlength,far_splicetype,
				chroffset,chrhigh,&(*dynprogindex),dynprog,
				revsequence1,revsequenceuc1,length1,length2,revoffset1,revoffset2,
				cdna_direction,watsonp,genestrand,jump_late_p,
				genome,genomealt,pairpool,extraband_end,defect_rate);
    debug7(printf("\n"));

    /* Check for unique splicecoord */
    size = (int) (coordsptr - coords);
    if (size > 1) {
      splicecoord0 = coords[0];
      i = 1;
      while (i < size && coords[i] == splicecoord0) {
	i++;
      }
      if (i < size) {
	/* Signal non-uniqueness or ambiguity */
	best_pairs = (List_T) NULL;
      }
    }
    
    FREEA(coords);
    return best_pairs;
  }
}


List_T
Splicetrie_solve_end3 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2, int offset1, int offset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		       int extraband_end, double defect_rate) {
  Univcoord_T *coordsptr, *coords, splicecoord0;
  unsigned int *triestart;
  int size, i;

  debug7(printf("Entering Splicetrie_solve_end3 with limits %u..%u, anchor splicesite %u (%u)\n",
		knownsplice_limit_low,knownsplice_limit_high,anchor_splicesite,anchor_splicesite - chroffset));
  if (trieoffsets[j] == NULL_POINTER) {
    return best_pairs;
  } else {
    triestart = &(triecontents[trieoffsets[j]]);
  }

  if ((size = splicetrie_size(triestart)) == 0) {
    return best_pairs;
  } else {
    coordsptr = coords = (Univcoord_T *) MALLOCA(size * sizeof(Univcoord_T));

    best_pairs = solve_end3_aux(&coordsptr,coords,best_pairs,triestart,
				knownsplice_limit_low,knownsplice_limit_high,
				&(*finalscore),&(*nmatches),&(*nmismatches),
				&(*nopens),&(*nindels),&(*knownsplicep),&(*ambig_end_length),
				&(*threshold_miss_score),obsmax_penalty,perfect_score,
				anchor_splicesite,splicejunction,splicejunction_alt,
				splicelength,contlength,far_splicetype,
				chroffset,chrhigh,&(*dynprogindex),dynprog,
				sequence1,sequenceuc1,length1,length2,offset1,offset2,
				cdna_direction,watsonp,genestrand,jump_late_p,
				genome,genomealt,pairpool,extraband_end,defect_rate);
    debug7(printf("\n"));

    /* Check for unique splicecoord */
    size = (int) (coordsptr - coords);
    if (size > 1) {
      splicecoord0 = coords[0];
      i = 1;
      while (i < size && coords[i] == splicecoord0) {
	i++;
      }
      if (i < size) {
	/* Signal non-uniqueness or ambiguity */
	best_pairs = (List_T) NULL;
      }
    }

    FREEA(coords);
    return best_pairs;
  }
}


