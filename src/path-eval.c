#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-eval.h"
#include "debug_groups.h"
#include "path-solve.h"
#include "path-trim.h"
#include "transcript-remap.h"

#include <stdio.h>
#include <string.h>		/* For strcpy */
#include <math.h>		/* For rint */
#include <ctype.h>		/* For islower */
#include <limits.h>
#include "fastlog.h"		/* For fasterexp */

#include "assert.h"
#include "list.h"
#include "genomebits_count.h"
#include "junction.h"
#include "mapq.h"

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Transcriptome_T transcriptome;

static bool *circularp;
static bool *chrsubsetp;
static bool *altlocp;

static int index1part;
static int index1interval;

static Outputtype_T output_type;
static bool md_report_snps_p;
static bool want_random_p;

static bool allow_soft_clips_start_p;
static bool allow_soft_clips_end_p;

/* Set to index1part, because segment search cannot do any better than this */
static int end_ambig_allowed;


/* Splice goodness now determined by Splice_accept_p */
/* #define OMIT_BAD_SPLICES 1 */

#define NMATCHES_FACTOR 8
#define CHOPPED_SLOP 1

static bool quality_ranker_p = false;
static double quality_match_reward = 1.0;
static double quality_mismatch_penalty = 2.0;
static double quality_indel_penalty = 2.0;
static double quality_softclip_penalty = 0.0;
static double quality_error_rate_penalty = 0.0;
static double quality_rank_margin = 0.20;
static bool quality_split_ranker_p = false;
static double quality_split_gap_penalty = 0.0;
static double quality_split_segment_penalty = 0.0;


#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Path_common_structure_trim_p */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Path_eval_nmatches_chopped */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Path_eval_nmatches */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Path_eval_and_sort */
#if defined(DEBUG8) || defined(DEBUG_TRANSCRIPTOME)
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Path_consolidate */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define T Path_T


/* Previously, some procedures removed a read if path->nmatches <
   nmismatches_allowed.  However, since this procedure sets
   found_score in such cases, found_score is no longer accurate.  One
   solution would be for caller to revise found_score.  Another
   solution is not to compare path->nmatches against
   nmismatches_allowed.  Also, could consider using found_score to
   constrain subsequent procedures. */

/* Sets found_score, score_within_trims, nmatches, ref_nmatches, junction_splice_prob, total_splice_prob */
int
Path_eval_nmatches (int *found_score, int *end_score,
		    T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int qstart, qend, ninserts;
  Univcoord_T univdiagonal;
  Intlist_T r, x, y;
  Univcoordlist_T q;
  Junction_T junction;
  List_T j;
  Altsplice_T altsplice;
  /* bool insertionp = false; */
  /* int adj0; deletions - insertions */
  /* int total_ninserts = 0; */
  int nmismatches, ref_nmismatches;
#ifdef CHECK_ASSERTIONS
  int orig_nmatches;
#endif


  debug7(printf("\nEntering Path_eval_nmatches on path %p, and comparing against %d\n",this,*found_score));
  debug7(Path_print(this));

  Path_expect_fwd(this);

#ifdef CHECK_ASSERTIONS
  orig_nmatches = this->nmatches;
#endif

  this->found_score = 0;
  this->score_within_trims = 0;
  this->nmatches = 0;
  this->ref_nmatches = 0;
  this->junction_splice_prob = 0.0;
  this->total_splice_prob = 0.0;


  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);


  qstart = Intlist_head(this->endpoints);
  qend = Intlist_second_value(this->endpoints);
  nmismatches = Intlist_head(this->nmismatches);
  ref_nmismatches = Intlist_head(this->ref_nmismatches);
  ninserts = 0;
  this->end_score = 0;
  
  if (this->plusp == true) {
    /* plus */
    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.9) {
      /* Count ambig splice site as matches, but counts against found_score, unless we are close to 5' end */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      if (qstart >= end_ambig_allowed) {
	this->found_score += qstart;
	this->end_score = qstart;
      }
      debug7(printf("(1) found score now %d\n",this->found_score));
      /* Does not affect score_within_trims */

    } else if ((altsplice = this->qstart_alts) != NULL && altsplice->nmismatches_known_p == true &&
	       (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      this->nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("(1) Adding %d matches for qstart_alts => total %d\n",
		    Altsplice_best_nmatches(this->qstart_alts),this->nmatches));
      this->ref_nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->total_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->found_score += Altsplice_distal_nmismatches(this->qstart_alts);
      this->end_score = Altsplice_distal_trimpos(this->qstart_alts);
      /* this->score_within_trims += Altsplice_distal_nmismatches(this->qstart_alts); */
      debug7(printf("(2) found score now %d\n",this->found_score));

    } else if (this->fusion_querystart_junction != NULL) {
      /* Rest of fusions handled below */

    } else {
      this->found_score += qstart;
      this->end_score = qstart;
      debug7(printf("(3) found score now %d\n",this->found_score));
      /* Does not affect score_within_trims */
    }

    /* Add qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
#if 0
      if (insertionp == true) {
	nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
      } else {
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
      }
#else
      nmismatches = Intlist_head(x);
      ref_nmismatches = Intlist_head(y);
#endif

      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */
      debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		    ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

      if (nmismatches >= 0 && ref_nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",univdiagonal - this->chroffset,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,
		      Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
								    univdiagonal,this->querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
		      nmismatches,ref_nmismatches,qstart,qend));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      this->nmatches += (qend - qstart) - nmismatches;
      debug7(printf("(1) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		    (qend-qstart)-nmismatches,qend,qstart,nmismatches,this->nmatches));
      this->ref_nmatches += (qend - qstart) - ref_nmismatches;
      this->found_score += nmismatches;
      this->score_within_trims += nmismatches;
      debug7(printf("(4) found score now %d\n",this->found_score));

      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	ninserts = Junction_ninserts(junction);
	this->junction_splice_prob += Junction_prob(junction);
	this->total_splice_prob += Junction_prob(junction);
      }
    }

  } else {
    /* minus */

    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.9) {
      /* Count ambig splice site as matches, but counts against found_score, unless we are close to 5' end */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      if (qstart >= end_ambig_allowed) {
	this->found_score += qstart;
	this->end_score = qstart;
      }
      debug7(printf("(5) found score now %d\n",this->found_score));
      /* Does not affect score_within_trims */

    } else if ((altsplice = this->qstart_alts )!= NULL && altsplice->nmismatches_known_p == true &&
	       (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      this->nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("(2) Adding %d matches for qstart_alts => total %d\n",
		    Altsplice_best_nmatches(this->qstart_alts),this->nmatches));
      this->ref_nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->total_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->found_score += Altsplice_distal_nmismatches(this->qstart_alts);
      this->end_score = Altsplice_distal_trimpos(this->qstart_alts);
      /* this->score_within_trims += Altsplice_distal_nmismatches(this->qstart_alts); */
      debug7(printf("(6) found score now %d\n",this->found_score));

    } else if (this->fusion_querystart_junction != NULL) {
      /* Fusions handled below */

    } else {
      this->found_score += qstart;
      this->end_score = qstart;
      debug7(printf("(7) found score now %d\n",this->found_score));
      /* Does not affect score_within_trims */
    }

    /* Subtract qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
#if 0
      if (insertionp == true) {
	nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
      } else {
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
      }
#else
      nmismatches = Intlist_head(x);
      ref_nmismatches = Intlist_head(y);
#endif
      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */
      debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		    ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

      if (nmismatches >= 0 && ref_nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",univdiagonal - this->chroffset,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
									univdiagonal,this->querylength,
									/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
								    univdiagonal,this->querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
							    univdiagonal,this->querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
		      nmismatches,ref_nmismatches,this->querylength - qend,this->querylength - qstart));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      this->nmatches += (qend - qstart) - nmismatches;
      debug7(printf("(2) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		    (qend-qstart)-nmismatches,this->querylength - qstart,this->querylength - qend,nmismatches,this->nmatches));
      this->ref_nmatches += (qend - qstart) - ref_nmismatches;
      this->found_score += nmismatches;
      this->score_within_trims += nmismatches;
      debug7(printf("(8) found score now %d\n",this->found_score));

      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	ninserts = Junction_ninserts(junction);
	this->junction_splice_prob += Junction_prob(junction);
	this->total_splice_prob += Junction_prob(junction);
      }
    }
  }

  if (this->ambig_prob_3 > 0.9) {
    /* Count ambig splice site as matches, but counts against found_score, unless we are close to 3' end */
    this->total_splice_prob += this->ambig_prob_3;
    /* this->nmatches += this->querylength - qend; (inferred) */
    /* this->ref_nmatches += this->querylength - qend; (inferred) */
    if (this->querylength - qend >= end_ambig_allowed) {
      this->found_score += this->querylength - qend;
      if (this->querylength - qend > this->end_score) {
	this->end_score = this->querylength - qend;
      }
    }
    debug7(printf("(9) found score now %d\n",this->found_score));
    /* Does not affect score_within_trims */

  } else if ((altsplice = this->qend_alts) != NULL && altsplice->nmismatches_known_p == true &&
	     (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
    this->nmatches += Altsplice_best_nmatches(this->qend_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
    debug7(printf("(3) Adding %d matches for qend_alts => total %d\n",
		  Altsplice_best_nmatches(this->qend_alts),this->nmatches));
    this->ref_nmatches += Altsplice_best_nmatches(this->qend_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
    this->junction_splice_prob += this->qend_alts->medial_prob + Altsplice_best_distal_prob(this->qend_alts);
    this->total_splice_prob += this->qend_alts->medial_prob + Altsplice_best_distal_prob(this->qend_alts);
    this->found_score += Altsplice_distal_nmismatches(this->qend_alts);
    if (Altsplice_distal_trimpos(this->qend_alts) > this->end_score) {
      this->end_score = Altsplice_distal_trimpos(this->qend_alts);
    }

    /* this->score_within_trims += Altsplice_distal_nmismatches(this->qend_alts); */
    debug7(printf("(10) found score now %d\n",this->found_score));

  } else if (this->fusion_queryend_junction != NULL) {
    /* Rest of fusions handled below */

  } else {
    this->found_score += this->querylength - qend;
    if (this->querylength - qend > this->end_score) {
      this->end_score = this->querylength - qend;
    }
    debug7(printf("(11) found score now %d\n",this->found_score));
    /* Does not affect score_within_trims */
  }



  /* Fusion */
  if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
    this->junction_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->junction_splice_prob += Junction_prob(this->fusion_queryend_junction);
    this->total_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->total_splice_prob += Junction_prob(this->fusion_queryend_junction);

    qstart = Intlist_head(this->fusion_endpoints);
    nmismatches = Intlist_head(this->fusion_nmismatches);
    ref_nmismatches = Intlist_head(this->fusion_ref_nmismatches);
    ninserts = 0;

    if (this->fusion_querystart_junction != NULL && this->fusion_plusp == this->plusp) {
      this->found_score += qstart;
      this->end_score = qstart;
      debug7(printf("(12) found score now %d\n",this->found_score));
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp != this->plusp) {
      this->found_score += qstart;
      this->end_score = qstart;
      debug7(printf("(13) found score now %d\n",this->found_score));
    }

    if (this->fusion_plusp == true) {
      j = this->fusion_junctions;
      /* insertionp = false; */

      /* Add qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);
#if 0
	if (insertionp == true) {
	  nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
	} else {
	  nmismatches = Intlist_head(x);
	  ref_nmismatches = Intlist_head(y);
	}
#else
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
#endif
	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */
	debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		      ninserts,qstart,qend,univdiagonal,univdiagonal - this->chroffset));

	if (nmismatches >= 0 && ref_nmismatches >= 0) {
	  debug7(printf("Checking fusion mismatches, plus, at %u from querystart %d to queryend %d\n",
			univdiagonal - this->fusion_chroffset,qstart,qend));
	  debug7(printf("%d mismatches expected vs %d measured\n",
			nmismatches,
			Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	  assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
								      univdiagonal,this->querylength,
								      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand));
#endif
	} else {
	  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,this->genestrand);
	  Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	  Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	  debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
			nmismatches,ref_nmismatches,qstart,qend));
	}

	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	this->nmatches += (qend - qstart) - nmismatches;
	this->ref_nmatches += (qend - qstart) - ref_nmismatches;
	debug7(printf("(3) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		      (qend-qstart)-nmismatches,qend,qstart,nmismatches,this->nmatches));
	this->found_score += nmismatches;
	this->score_within_trims += nmismatches;
	debug7(printf("(14) found score now %d\n",this->found_score));

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	  ninserts = Junction_ninserts(junction);
	  this->junction_splice_prob += Junction_prob(junction);
	  this->total_splice_prob += Junction_prob(junction);
	}
      }

    } else {
      j = this->fusion_junctions; /* Put here before we handle querystart_alts */

      /* Subtract qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);
#if 0
	if (insertionp == true) {
	  nmismatches = ref_nmismatches = -1; /* Recompute nmismatches */
	} else {
	  nmismatches = Intlist_head(x);
	  ref_nmismatches = Intlist_head(y);
	}
#else
	nmismatches = Intlist_head(x);
	ref_nmismatches = Intlist_head(y);
#endif

	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */
	debug7(printf("Path_eval_nmatches: ninserts %d, qstart %d..qend %d at univdiagonal %u [%u]\n",
		      ninserts, qstart,qend,univdiagonal,univdiagonal - this->chroffset));

	if (nmismatches >= 0 && ref_nmismatches >= 0) {
	  debug7(printf("Checking fusion mismatches at %u, minus, from querystart %d to queryend %d\n",
			univdiagonal - this->fusion_chroffset,qstart,qend));
	  debug7(printf("%d mismatches expected vs %d measured\n",
			nmismatches,Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
									  univdiagonal,this->querylength,
									  /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand)));
#ifdef CHECK_NMISMATCHES
	  assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
								      univdiagonal,this->querylength,
								      /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand));
#endif
	} else {
	  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
							      univdiagonal,this->querylength,
							      /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,this->genestrand);
	  Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	  Intlist_head_set(y,ref_nmismatches);		/* Save for Stage3end_new_substrings */
	  debug7(printf("%d (%d ref) mismatches from genome over querypos %d..%d\n",
			nmismatches,ref_nmismatches,this->querylength - qend,this->querylength - qstart));
	}

	/* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
	this->nmatches += (qend - qstart) - nmismatches;
	this->ref_nmatches += (qend - qstart) - ref_nmismatches;
	debug7(printf("(4) Shortcut adds matches of %d = (%d - %d) - nmismatches %d => total %d\n",
		      (qend-qstart)-nmismatches,this->querylength - qstart,this->querylength - qend,nmismatches,this->nmatches));
	this->found_score += nmismatches;
	this->score_within_trims += nmismatches;
	debug7(printf("(15) found score now %d\n",this->found_score));

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  debug7(printf("Junction: ")); debug7(Junction_print(junction)); debug7(printf("\n"));
	  ninserts = Junction_ninserts(junction);
	  this->junction_splice_prob += Junction_prob(junction);
	  this->total_splice_prob += Junction_prob(junction);
	}
      }
    }

    if (this->fusion_querystart_junction != NULL && this->fusion_plusp != this->plusp) {
      this->found_score += this->querylength - qend;
      debug7(printf("(16) found score now %d\n",this->found_score));
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp == this->plusp) {
      this->found_score += this->querylength - qend;
      debug7(printf("(17) found score now %d\n",this->found_score));
    }
  }


  /* this->nmatches += total_ninserts; */

  if (this->found_score < *found_score) {
    *found_score = this->found_score;
  }
  if (this->end_score < *end_score) {
    *end_score = this->end_score;
  }

  debug7(printf("Path_eval_nmatches returning %d matches (orig %d) in %s for score of %d within trims\n",
		this->nmatches,orig_nmatches,Intlist_to_string(this->endpoints),this->score_within_trims));

  assert(orig_nmatches < 0 || this->nmatches == orig_nmatches);

  assert(this->nmatches <= this->querylength);
  return this->nmatches;
}



int
Path_eval_nmatches_chopped (T this) {
  Compress_T query_compress_fwd = this->query_compress_fwd;
  Compress_T query_compress_rev = this->query_compress_rev;
  int chop_qstart = this->chop_qstart;
  int chop_qend = this->chop_qend;

  int nmatches = 0, ref_nmatches = 0;
  int overlap_qstart, overlap_qend, qstart, qend, ninserts;
  Univcoord_T univdiagonal;
  Intlist_T r, x, y;
  Univcoordlist_T q;
  Junction_T junction;
  List_T j;
  /* bool insertionp = false; */
  /* int adj0; deletions - insertions */
  /* int total_ninserts = 0; */
  int ref_nmismatches;


  debug6(printf("\nEntering Path_eval_nmatches_chopped on path %p\n",this));
  debug6(Path_print(this));

  Path_expect_fwd(this);

  if (chop_qend <= chop_qstart) {
    return 0;
  }

  qstart = Intlist_head(this->endpoints);
  ninserts = 0;
  
  if (this->plusp == true) {
    /* plus */
    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.9) {
      /* Count ambig splice site as matches, but counts against found_score */
      /* nmatches += qstart; (inferred) */
      /* ref_nmatches += qstart; (inferred) */
      /* Does not affect score_within_trims */

#if 0
    } else if ((altsplice = this->qstart_alts) != NULL && altsplice->nmismatches_known_p == true &&
	       (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      /* Not considering altsplices */
      nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug6(printf("Adding %d matches for qstart_alts => total %d\n",this->qstart_alts->best_distal_nmatches,nmatches));
      ref_nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
#endif

    } else if (this->fusion_querystart_junction != NULL) {
      /* Rest of fusions handled below */

    } else {
      /* Does not affect score_within_trims */
    }

    /* Add qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);

      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */

      if (qend <= chop_qstart) {
	/* Segment is to the left of the chop boundary */

      } else if (qstart > chop_qend) {
	/* Segment is to the right of the chop boundary */

      } else if (qstart > chop_qstart && qend <= chop_qend) {
	/* Segment is within chop boundaries */
	assert(Intlist_head(x) >= 0);
	nmatches += (qend - qstart) - Intlist_head(x);
	ref_nmatches += (qend - qstart) - Intlist_head(y);

      } else {
	overlap_qstart = (qstart > chop_qstart) ? qstart : chop_qstart;
	overlap_qend = (qend <= chop_qend) ? qend : chop_qend;

	/* Count only in the overlapping region */
	nmatches += (overlap_qend - overlap_qstart) - 
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
						univdiagonal,this->querylength,
						/*pos5*/overlap_qstart,/*pos3*/overlap_qend,/*plusp*/true,this->genestrand);
	ref_nmatches += (overlap_qend - overlap_qstart) - ref_nmismatches;
      }

      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	ninserts = Junction_ninserts(junction);
      }
    }

  } else {
    /* minus */

    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.9) {
      /* Count ambig splice site as matches, but counts against found_score */
      /* nmatches += qstart; (inferred) */
      /* ref_nmatches += qstart; (inferred) */
      /* Does not affect score_within_trims */

#if 0
    } else if ((altsplice = this->qstart_alts )!= NULL && altsplice->nmismatches_known_p == true &&
	       (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      /* Not considering altsplices */
      nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
      ref_nmatches += this->qstart_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
#endif

    } else if (this->fusion_querystart_junction != NULL) {
      /* Fusions handled below */

    } else {
      /* Does not affect score_within_trims */
    }

    /* Subtract qpos to get alignstart/alignend */
    for (q = this->univdiagonals, x = this->nmismatches, y = this->ref_nmismatches, r = Intlist_next(this->endpoints); q != NULL;
	 q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);

      univdiagonal = Univcoordlist_head(q);
      /* left = univdiagonal - (Univcoord_T) this->querylength; */

      if (qend <= chop_qstart) {
	/* Segment is to the left of the chop boundary */

      } else if (qstart > chop_qend) {
	/* Segment is to the right of the chop boundary */

      } else if (qstart > chop_qstart && qend <= chop_qend) {
	/* Segment is within chop boundaries */
	assert(Intlist_head(x) >= 0);
	nmatches += (qend - qstart) - Intlist_head(x);
	ref_nmatches += (qend - qstart) - Intlist_head(y);

      } else {
	overlap_qstart = (qstart > chop_qstart) ? qstart : chop_qstart;
	overlap_qend = (qend <= chop_qend) ? qend : chop_qend;

	/* Count only in the overlapping region */
	nmatches += (overlap_qend - overlap_qstart) - 
	  Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
							    univdiagonal,this->querylength,
							    /*pos5*/overlap_qstart,/*pos3*/overlap_qend,/*plusp*/false,this->genestrand);
	ref_nmatches += (overlap_qend - overlap_qstart) - ref_nmismatches;
      }


      /* Prepare for next iteration */
      qstart = qend;

      if (j == NULL) {
	ninserts = 0;
      } else if ((junction = (Junction_T) List_head(j)) == NULL) {
	/* qstart_junction */
	ninserts = 0;
      } else {
	ninserts = Junction_ninserts(junction);
      }
    }
  }

  if (this->ambig_prob_3 > 0.9) {
    /* Count ambig splice site as matches, but counts against found_score */
    /* nmatches += this->querylength - qend; (inferred) */
    /* ref_nmatches += this->querylength - qend; (inferred) */
    /* Does not affect score_within_trims */

#if 0
  } else if ((altsplice = this->qend_alts) != NULL && altsplice->nmismatches_known_p == true &&
	     (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
    /* Not considering altsplices */
    nmatches += this->qend_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
    ref_nmatches += this->qend_alts->best_distal_nmatches; /* Not nmatches_to_trims, which is 0 for alts_substring */
#endif

  } else if (this->fusion_queryend_junction != NULL) {
    /* Rest of fusions handled below */

  } else {
    /* Does not affect score_within_trims */
  }



  /* Fusion */
  if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
    this->junction_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->junction_splice_prob += Junction_prob(this->fusion_queryend_junction);
    this->total_splice_prob += Junction_prob(this->fusion_querystart_junction);
    this->total_splice_prob += Junction_prob(this->fusion_queryend_junction);

    qstart = Intlist_head(this->fusion_endpoints);
    ninserts = 0;

    if (this->fusion_plusp == true) {
      j = this->fusion_junctions;
      /* insertionp = false; */

      /* Add qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);

	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */

	if (qend <= chop_qstart) {
	  /* Segment is to the left of the chop boundary */
	  
	} else if (qstart > chop_qend) {
	  /* Segment is to the right of the chop boundary */
	  
	} else if (qstart > chop_qstart && qend <= chop_qend) {
	  /* Segment is within chop boundaries */
	  assert(Intlist_head(x) >= 0);
	  nmatches += (qend - qstart) - Intlist_head(x);
	  ref_nmatches += (qend - qstart) - Intlist_head(y);

	} else {
	  overlap_qstart = (qstart > chop_qstart) ? qstart : chop_qstart;
	  overlap_qend = (qend <= chop_qend) ? qend : chop_qend;

	  /* Count only in the overlapping region */
	  nmatches += (overlap_qend - overlap_qstart) - 
	    Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_fwd,
						  univdiagonal,this->querylength,
						  /*pos5*/overlap_qstart,/*pos3*/overlap_qend,/*plusp*/true,this->genestrand);
	  ref_nmatches += (overlap_qend - overlap_qstart) - ref_nmismatches;
	}

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  ninserts = Junction_ninserts(junction);
	}
      }

    } else {
      j = this->fusion_junctions; /* Put here before we handle querystart_alts */

      /* Subtract qpos to get alignstart/alignend */
      for (q = this->fusion_univdiagonals, x = this->fusion_nmismatches,
	     y = this->fusion_ref_nmismatches, r = Intlist_next(this->fusion_endpoints); q != NULL;
	   q = Univcoordlist_next(q), x = Intlist_next(x), y = Intlist_next(y), r = Intlist_next(r), j = List_next(j)) {
	qstart += ninserts;
	qend = Intlist_head(r);

	univdiagonal = Univcoordlist_head(q);
	/* left = univdiagonal - (Univcoord_T) this->querylength; */

	if (qend <= chop_qstart) {
	  /* Segment is to the left of the chop boundary */
	  
	} else if (qstart > chop_qend) {
	  /* Segment is to the right of the chop boundary */
	  
	} else if (qstart > chop_qstart && qend <= chop_qend) {
	  /* Segment is within chop boundaries */
	  assert(Intlist_head(x) >= 0);
	  nmatches += (qend - qstart) - Intlist_head(x);
	  ref_nmatches += (qend - qstart) - Intlist_head(y);
	  
	} else {
	  overlap_qstart = (qstart > chop_qstart) ? qstart : chop_qstart;
	  overlap_qend = (qend <= chop_qend) ? qend : chop_qend;

	  /* Count only in the overlapping region */
	  nmatches += (overlap_qend - overlap_qstart) - 
	    Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress_rev,
						  univdiagonal,this->querylength,
						  /*pos5*/overlap_qstart,/*pos3*/overlap_qend,/*plusp*/false,this->genestrand);
	  ref_nmatches += (overlap_qend - overlap_qstart) - ref_nmismatches;
	}

	/* Prepare for next iteration */
	qstart = qend;

	if (j == NULL) {
	  ninserts = 0;
	} else if ((junction = (Junction_T) List_head(j)) == NULL) {
	  /* qstart_junction */
	  ninserts = 0;
	} else {
	  ninserts = Junction_ninserts(junction);
	}
      }
    }
  }

  /* nmatches += total_ninserts; */

  debug6(printf("Path_eval_nmatches_chopped returning %d matches in %s, chopped region %d..%d\n",
		nmatches,Intlist_to_string(this->endpoints),chop_qstart,chop_qend));

  return nmatches;
}


int
Path_nmatches_cmp (T a, T b) {

  int chopped_nmatches_a, chopped_nmatches_b;
  int qstart_a = Intlist_head(a->endpoints);
  int qstart_b = Intlist_head(b->endpoints);
  int qend_a = Intlist_last_value(a->endpoints);
  int qend_b = Intlist_last_value(b->endpoints);

  if (a->chop_qstart == 0 && a->chop_qend == a->querylength) {
    /* No need to chop, since we can use the entire region */
  } else if (qstart_a == qstart_b && qend_a == qend_b) {
    /* No need to chop, since the two chopped regions are the same */
  } else if (a->nmatches > b->nmatches + NMATCHES_FACTOR) {
    /* Skip because one is clearly better */
  } else if (b->nmatches > a->nmatches + NMATCHES_FACTOR) {
    /* Skip because one is clearly better */
  } else {
    chopped_nmatches_a = Path_eval_nmatches_chopped(a);
    chopped_nmatches_b = Path_eval_nmatches_chopped(b);

    if (chopped_nmatches_a > chopped_nmatches_b + CHOPPED_SLOP) {
      return -1;
    } else if (chopped_nmatches_b > chopped_nmatches_a + CHOPPED_SLOP) {
      return +1;
    } else {
      /* Fall through to the comparison below */
    }
  }

  if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_quality_softclip (T path) {
  int softclip;

  softclip = Intlist_head(path->endpoints) + path->querylength - Intlist_last_value(path->endpoints);
  if (path->found_score - path->score_within_trims > softclip) {
    softclip = path->found_score - path->score_within_trims;
  }
  if (path->querylength - (path->nmatches + path->score_within_trims) > softclip) {
    softclip = path->querylength - (path->nmatches + path->score_within_trims);
  }

  return softclip;
}


static double
Path_quality_score_components (int *mismatches, int *nindels, int *softclip,
			       double *error_rate, T path) {
  double aligned_bases, errors;

  *mismatches = path->score_within_trims;
  *nindels = Path_nindels(path);
  *softclip = Path_quality_softclip(path);

  aligned_bases = (double) Path_coverage(path);
  if (aligned_bases < 1.0) {
    aligned_bases = 1.0;
  }
  errors = (double) (*mismatches + *nindels);
  *error_rate = errors/aligned_bases;

  return quality_match_reward * (double) path->nmatches
    - quality_mismatch_penalty * (double) (*mismatches)
    - quality_indel_penalty * (double) (*nindels)
    - quality_softclip_penalty * (double) (*softclip)
    - quality_error_rate_penalty * errors * (*error_rate);
}


static double
Path_quality_score (T path) {
  int mismatches, nindels, softclip;
  double error_rate;

  return Path_quality_score_components(&mismatches,&nindels,&softclip,&error_rate,path);
}


static void
Path_quality_debug (const char *label, const char *action, double score, T path) {
  int mismatches, nindels, softclip;
  double error_rate;

  (void) Path_quality_score_components(&mismatches,&nindels,&softclip,&error_rate,path);
  debug8(printf("(%s quality ranker) %s score %.3f [matches=%d mismatches=%d indels=%d softclip=%d error_rate=%.4f] ",
		label,action,score,path->nmatches,mismatches,nindels,softclip,error_rate));
  debug8(Path_print(path));
  return;
}


static int
Path_quality_cmp (T a, T b) {
  double quality_a = Path_quality_score(a);
  double quality_b = Path_quality_score(b);

  if (quality_a > quality_b) {
    return -1;
  } else if (quality_b > quality_a) {
    return +1;
  } else {
    return 0;
  }
}


static double
Path_quality_threshold (double best_score) {
  return best_score - fabs(best_score) * quality_rank_margin;
}


static int
Path_filter_quality_score (T *patharray, int npaths, const char *label,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  int i, k = 0;
  double score, best_score, threshold;
  T path;

  if (quality_ranker_p == false || npaths <= 1) {
    return npaths;
  }

  best_score = Path_quality_score(patharray[0]);
  for (i = 1; i < npaths; i++) {
    score = Path_quality_score(patharray[i]);
    if (score > best_score) {
      best_score = score;
    }
  }
  threshold = Path_quality_threshold(best_score);

  debug8(printf("(%s quality ranker) best score %.3f, threshold %.3f over %d paths\n",
		label,best_score,threshold,npaths));

  for (i = 0; i < npaths; i++) {
    path = patharray[i];
    score = Path_quality_score(path);
    if (score >= threshold) {
      Path_quality_debug(label,/*action*/"Keeping",score,path);
      patharray[k++] = path;
    } else {
      Path_quality_debug(label,/*action*/"Eliminating",score,path);
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
  }

  return k;
}


#if 0
void
Path_eval_extend (T this, int querystart, int queryend,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int qstart, qend;
  bool changedp = false;
  int ignore_score = 0;

  if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
    /* Cannot extend */
    return;
  } else if (this->plusp == true) {
    qstart = querystart;
    qend = queryend;
  } else {
    qstart = this->querylength - queryend;
    qend = this->querylength - querystart;
  }

  if (this->qstart_alts != NULL) {
    /* Cannot extend */
  } else if (Intlist_head(this->endpoints) == qstart) {
    /* No need to change */
  } else if (Univcoordlist_head(this->univdiagonals) + qstart < this->chroffset + this->querylength) {
    /* Cannot extend past chromosomal bounds */
  } else {
    Intlist_head_set(this->endpoints,qstart);
    Intlist_head_set(this->nmismatches,-1);
    Intlist_head_set(this->ref_nmismatches,-1);
    changedp = true;
  }

  if (this->qend_alts != NULL) {
    /* Cannot extend */
  } else if (Intlist_last_value(this->endpoints) == qend) {
    /* No need to change */
  } else if (Univcoordlist_last_value(this->univdiagonals) + qend > this->chrhigh + this->querylength) {
    /* Cannot extend past chromosomal bounds */
  } else {
    Intlist_last_value_set(this->endpoints,qend);
    Intlist_last_value_set(this->nmismatches,-1);
    Intlist_last_value_set(this->ref_nmismatches,-1);
    changedp = true;
  }

  if (changedp == true) {
    Path_eval_nmatches(&ignore_score,&ignore_score,this,query_compress_fwd,query_compress_rev);
  }

  return;
}
#endif



void
Path_mark_alignment (T path, Compress_T query_compress_fwd, char *queryuc_ptr,
		     Compress_T query_compress_rev, char *queryrc, Pathpool_T pathpool) {

  Compress_T query_compress;
  Univcoord_T univdiagonal;
  Intlist_T q, n;
  Univcoordlist_T p;
  List_T j;
  int nmatches_exonic;
  int pos5, pos3, ninserts;

  if (path->genomic_diff != NULL) {
    /* Already marked */
    return;

  } else {
    path->genomic_diff = Pathpool_new_string(pathpool,path->querylength+1);
    /* genomic_bothdiff = (char *) MALLOC((path->querylength+1) * sizeof(char)); */
    /* Genome_fill_buffer(left,querylength,genomic_diff); */
  }

  if (path->plusp == true) {
    query_compress = query_compress_fwd;
    strcpy(path->genomic_diff,queryuc_ptr); /* Start with query sequence on genomic plus strand */
  } else {
    query_compress = query_compress_rev;
    strcpy(path->genomic_diff,queryrc); /* Start with query sequence on genomic plus strand */
  }

  ninserts = 0;
  for (p = path->univdiagonals, q = path->endpoints, n = path->nmismatches, j = path->junctions; p != NULL;
       p = Univcoordlist_next(p), q = Intlist_next(q), n = Intlist_next(n), j = List_next(j)) {
    assert(Intlist_head(n) >= 0);
    if (Intlist_head(n) > 0 || Compress_non_acgt(query_compress) == true || md_report_snps_p == true) {
      univdiagonal = Univcoordlist_head(p);
      pos5 = Intlist_head(q) + ninserts;
      pos3 = Intlist_head(Intlist_next(q));
      Genomebits_mark_mismatches(&nmatches_exonic,path->genomic_diff,
				 query_compress,univdiagonal,path->querylength,
				 pos5,pos3,/*segment_plusp*/path->plusp,
				 /*query_plusp*/path->plusp,path->genestrand);
    }

    if (j != NULL) {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }


  /* Fusion */
  if (path->fusion_querystart_junction != NULL || path->fusion_queryend_junction != NULL) {
    if (path->fusion_plusp == true) {
      query_compress = query_compress_fwd;
    } else {
      query_compress = query_compress_rev;
    }

    ninserts = 0;
    for (p = path->fusion_univdiagonals, q = path->fusion_endpoints,
	   n = path->fusion_nmismatches, j = path->fusion_junctions; p != NULL;
	 p = Univcoordlist_next(p), q = Intlist_next(q), n = Intlist_next(n), j = List_next(j)) {
      if (Intlist_head(n) > 0 || Compress_non_acgt(query_compress) == true || md_report_snps_p == true) {
	univdiagonal = Univcoordlist_head(p);
	pos5 = Intlist_head(q) + ninserts;
	pos3 = Intlist_head(Intlist_next(q));
	Genomebits_mark_mismatches(&nmatches_exonic,path->genomic_diff,
				   query_compress,univdiagonal,path->querylength,
				   pos5,pos3,/*segment_plusp*/path->fusion_plusp,
				   /*query_plusp*/path->plusp,path->genestrand);
      }
    }

    if (j != NULL) {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }

  return;
}


bool
Path_eval_perfect_ends_p (T this, Compress_T query_compress_fwd, char *queryuc_ptr,
			  Compress_T query_compress_rev, char *queryrc,
			  int querystart, int queryend, Pathpool_T pathpool) {
  int mod;
  int querypos5, querypos3;

  if (this->found_score == 0) {
    return true;
  } else {
    Path_mark_alignment(this,query_compress_fwd,queryuc_ptr,
			query_compress_rev,queryrc,pathpool);
    for (mod = 0; mod < index1interval; mod++) {
      querypos5 = querystart + mod;
      querypos3 = queryend - mod; /* Typically query_lastpos - mod */
      if (querypos3 < querypos5) {
	/* For short querystart..queryend, querypos3 could be negative */
	return false;
      } else if (islower(this->genomic_diff[querypos5])) {
	return false;
      } else if (islower(this->genomic_diff[querypos5 + index1part - 1])) {
	return false;
      } else if (islower(this->genomic_diff[querypos3])) {
	return false;
      } else if (islower(this->genomic_diff[querypos3 + index1part - 1])) {
	return false;
      }
    }

    return true;
  }
}


static int
Path_method_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp, coverage_a, coverage_b, nindels_a, nindels_b;
  int qstart_a, qstart_b, qend_a, qend_b;

  if (a->method > b->method) {
    return -1;
  } else if (b->method > a->method) {
    return +1;
  } else if ((cmp = Path_quality_cmp(a,b)) != 0) {
    return cmp;
  } else if ((coverage_a = Path_coverage(a)) > (coverage_b = Path_coverage(b))) {
    return -1;
  } else if (coverage_b > coverage_a) {
    return +1;
  } else if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;
  } else if (a->score_within_trims < b->score_within_trims) {
    return -1;
  } else if (b->score_within_trims < a->score_within_trims) {
    return +1;
  } else if ((nindels_a = Path_nindels(a)) < (nindels_b = Path_nindels(b))) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  } else {
    /* Path_structure_cmp ignores terminal trims, so make this deterministic. */
    qstart_a = (a->qstart_alts != NULL) ? 0 : Intlist_min(a->endpoints);
    qstart_b = (b->qstart_alts != NULL) ? 0 : Intlist_min(b->endpoints);
    qend_a = (a->qend_alts != NULL) ? a->querylength : Intlist_max(a->endpoints);
    qend_b = (b->qend_alts != NULL) ? b->querylength : Intlist_max(b->endpoints);

    if (qstart_a < qstart_b) {
      return -1;
    } else if (qstart_b < qstart_a) {
      return +1;
    } else if (qend_a > qend_b) {
      return -1;
    } else if (qend_b > qend_a) {
      return +1;
    } else {
      return 0;
    }
  }
}


/* Used to identify best paths for pathstore */
/* Assumes that chop_qstart and chop_qend have been set for the two Pair_T objects */
int
Path_best_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int cmp;
  int coverage_a, coverage_b;
  int nindels_a, nindels_b;
#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

#ifdef OMIT_BAD_SPLICES
  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);
#endif

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  }

#ifdef OMIT_BAD_SPLICES
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif

  if ((cmp = Path_nmatches_cmp(a,b)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  }

  nindels_a = Path_nindels(a);
  nindels_b = Path_nindels(b);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  } else {
    return 0;
  }
}


/* Used to identify best paths in comparing partners in stage1hr-paired.c */
int
Path_best_pair_cmp (T a, T b) {

  int cmp;
  int chop_qstart, chop_qend, qstart_a, qstart_b, qend_a, qend_b;

  int coverage_a, coverage_b;
  int nindels_a, nindels_b;
#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

#ifdef OMIT_BAD_SPLICES
  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);
#endif

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  }

#ifdef OMIT_BAD_SPLICES
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif

  /* Need to set chop_qstart and chop_qend for the two paths */
  qstart_a = Intlist_head(a->endpoints);
  qstart_b = Intlist_head(b->endpoints);
  qend_a = Intlist_last_value(a->endpoints);
  qend_b = Intlist_last_value(b->endpoints);

  chop_qstart = /*max*/ (qstart_a > qstart_b) ? qstart_a : qstart_b;
  chop_qend = /*min*/ (qend_a < qend_b) ? qend_a : qend_b;

  a->chop_qstart = chop_qstart;
  a->chop_qend = chop_qend;
  b->chop_qstart = chop_qstart;
  b->chop_qend = chop_qend;

  if ((cmp = Path_nmatches_cmp(a,b)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  }

  nindels_a = Path_nindels(a);
  nindels_b = Path_nindels(b);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  } else {
    return 0;
  }
}


/* Duplicates with respect to method have already been taken care of */
/* Ignore sensedir, so we keep both sensedirs if they have equivalent matches */
static int
Path_local_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int cmp;
  int coverage_a, coverage_b;
  /* int nbadsplices_a, nbadsplices_b; */
  int nindels_a, nindels_b;
  int nalts_a = 0, nalts_b = 0;

#if 0
  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);
#endif

  bool fusionp_x = false, fusionp_y = false;

  if (a->fusion_querystart_junction != NULL ||
      a->fusion_queryend_junction != NULL) {
    fusionp_x = true;
  }
  if (b->fusion_querystart_junction != NULL ||
      b->fusion_queryend_junction != NULL) {
    fusionp_y = true;
  }

  if (fusionp_x == false && fusionp_y == true) {
    return -1;
  } else if (fusionp_y == false && fusionp_x == true) {
    return +1;
  }

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  }
    
#if 0
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif

  if ((cmp = Path_nmatches_cmp(a,b)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;

  } else if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;
  }

  nindels_a = Path_nindels(a);
  nindels_b = Path_nindels(b);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  }

  if (a->qstart_alts != NULL) {
    nalts_a++;
  }
  if (a->qend_alts != NULL) {
    nalts_a++;
  }
  if (b->qstart_alts != NULL) {
    nalts_b++;
  }
  if (b->qend_alts != NULL) {
    nalts_b++;
  }

  if (nalts_a < nalts_b) {
    return -1;
  } else if (nalts_b < nalts_a) {
    return +1;

#if 0
    /* Already compared sense/antisense */
    /* We want to keep all possible transcript results at a locus */
  } else if (a->junction_splice_prob > b->junction_splice_prob) {
    return -1;
  } else if (b->junction_splice_prob > a->junction_splice_prob) {
    return +1;
  } else if (a->total_splice_prob > b->total_splice_prob) {
    return -1;
  } else if (b->total_splice_prob > a->total_splice_prob) {
    return +1;
#endif

  } else {
    return 0;
  }
}


static int
Path_global_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  int cmp;
  int coverage_a, coverage_b;
#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif

  coverage_a = Path_coverage(a);
  coverage_b = Path_coverage(b);

#ifdef OMIT_BAD_SPLICES
  nbadsplices_a = Path_nbadsplices(a);
  nbadsplices_b = Path_nbadsplices(b);
#endif

  bool fusionp_x = false, fusionp_y = false;

  if (a->fusion_querystart_junction != NULL ||
      a->fusion_queryend_junction != NULL) {
    fusionp_x = true;
  }
  if (b->fusion_querystart_junction != NULL ||
      b->fusion_queryend_junction != NULL) {
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
    /* Do not consider for global cmp */
  if (nbadsplices_a < nbadsplices_b) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
  }
#endif

  if ((cmp = Path_nmatches_cmp(a,b)) < 0) {
    return -1;
  } else if (cmp > 0) {
    return +1;
  } else if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_fusion_cmp (T a, T b) {
  bool fusionp_x = false, fusionp_y = false;

  if (a->fusion_querystart_junction != NULL ||
      a->fusion_queryend_junction != NULL) {
    fusionp_x = true;
  }
  if (b->fusion_querystart_junction != NULL ||
      b->fusion_queryend_junction != NULL) {
    fusionp_y = true;
  }

  if (fusionp_x == false && fusionp_y == true) {
    return -1;
  } else if (fusionp_y == false && fusionp_x == true) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_transcript_cmp (T a, T b) {
  if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_nalts_cmp (T a, T b) {
  int nalts_a = 0, nalts_b = 0;

  if (a->qstart_alts != NULL) {
    nalts_a++;
  }
  if (a->qend_alts != NULL) {
    nalts_a++;
  }
  if (b->qstart_alts != NULL) {
    nalts_b++;
  }
  if (b->qend_alts != NULL) {
    nalts_b++;
  }

  if (nalts_a < nalts_b) {
    return -1;
  } else if (nalts_b < nalts_a) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_nindels_cmp (T a, T b) {
  int nindels_a = Path_nindels(a);
  int nindels_b = Path_nindels(b);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_coverage_cmp (T a, T b) {
  int coverage_a = Path_coverage(a);
  int coverage_b = Path_coverage(b);

  if (coverage_a > coverage_b + 20) {
    return -1;
  } else if (coverage_b > coverage_a + 20) {
    return +1;
  } else {
    return 0;
  }
}


static int
Path_local_quality_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;

  if ((cmp = Path_fusion_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_quality_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_nmatches_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_transcript_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_nindels_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_nalts_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_coverage_cmp(a,b)) != 0) {
    return cmp;
  } else {
    return 0;
  }
}


static int
Path_global_quality_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  int cmp;
#ifdef OMIT_BAD_SPLICES
  int nbadsplices_a, nbadsplices_b;
#endif

  if ((cmp = Path_fusion_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_quality_cmp(a,b)) != 0) {
    return cmp;

#ifdef OMIT_BAD_SPLICES
  } else if ((nbadsplices_a = Path_nbadsplices(a)) <
	     (nbadsplices_b = Path_nbadsplices(b))) {
    return -1;
  } else if (nbadsplices_b < nbadsplices_a) {
    return +1;
#endif

  } else if ((cmp = Path_nmatches_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_transcript_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_nindels_cmp(a,b)) != 0) {
    return cmp;
  } else if ((cmp = Path_coverage_cmp(a,b)) != 0) {
    return cmp;
  } else {
    return 0;
  }
}


static int
Path_local_cmp_active (const void *x, const void *y) {
  if (quality_ranker_p == true) {
    return Path_local_quality_cmp(x,y);
  } else {
    return Path_local_cmp(x,y);
  }
}


static int
Path_global_cmp_active (const void *x, const void *y) {
  if (quality_ranker_p == true) {
    return Path_global_quality_cmp(x,y);
  } else {
    return Path_global_cmp(x,y);
  }
}


static void
Path_query_span (int *qstart, int *qend, T path) {
  int path_qstart, path_qend;

  if (path->qstart_alts != NULL) {
    path_qstart = 0;
  } else {
    path_qstart = Intlist_min(path->endpoints);
  }

  if (path->qend_alts != NULL) {
    path_qend = path->querylength;
  } else {
    path_qend = Intlist_max(path->endpoints);
  }

  if (path->plusp == true) {
    *qstart = path_qstart;
    *qend = path_qend;
  } else {
    *qstart = path->querylength - path_qend;
    *qend = path->querylength - path_qstart;
  }

  return;
}


static int
Path_query_span_overlap (T a, T b) {
  int a_qstart, a_qend, b_qstart, b_qend, qstart, qend;

  Path_query_span(&a_qstart,&a_qend,a);
  Path_query_span(&b_qstart,&b_qend,b);

  qstart = (a_qstart > b_qstart) ? a_qstart : b_qstart;
  qend = (a_qend < b_qend) ? a_qend : b_qend;

  if (qend > qstart) {
    return qend - qstart;
  } else {
    return 0;
  }
}


static int
Path_query_span_union_coverage (T *patharray, int npaths) {
  bool foundp;
  int i, qstart, qend, trimmed_qstart, best_qstart, best_qend;
  int coverage = 0, covered_qend = 0;

  while (1) {
    foundp = false;
    best_qstart = best_qend = 0;

    for (i = 0; i < npaths; i++) {
      Path_query_span(&qstart,&qend,patharray[i]);
      if (qend <= covered_qend) {
	/* This path is already covered by the union built so far */
      } else {
	trimmed_qstart = (qstart < covered_qend) ? covered_qend : qstart;
	if (foundp == false || trimmed_qstart < best_qstart ||
	    (trimmed_qstart == best_qstart && qend > best_qend)) {
	  foundp = true;
	  best_qstart = trimmed_qstart;
	  best_qend = qend;
	}
      }
    }

    if (foundp == false) {
      return coverage;
    } else {
      coverage += best_qend - best_qstart;
      covered_qend = best_qend;
    }
  }
}


static bool
Path_split_complementary_p (T *patharray, int npaths, int split_path_max_overlap) {
  int i, j;

  for (i = 0; i < npaths; i++) {
    for (j = i + 1; j < npaths; j++) {
      if (Path_query_span_overlap(patharray[i],patharray[j]) <= split_path_max_overlap) {
	return true;
      }
    }
  }

  return false;
}


static bool
Path_split_compatible_with_kept_p (T path, T *patharray, int kept_start, int kept_end,
				   int split_path_max_overlap) {
  int i;

  for (i = kept_start; i < kept_end; i++) {
    if (Path_query_span_overlap(path,patharray[i]) > split_path_max_overlap) {
      return false;
    }
  }

  return true;
}


static T
Path_repair_terminal_suffix (T full_path, T suffix_path,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool) {
  Intlist_T endpoints = NULL, nmismatches = NULL, ref_nmismatches = NULL;
  Univcoordlist_T univdiagonals = NULL;
  List_T junctions = NULL;
  Intlist_T fq, fs, fr, sq, ss, sr;
  Univcoordlist_T fu, su;
  List_T fj, sj;
  T repaired = NULL;
  int suffix_qstart, suffix_qend, full_qend, cut_index = -1, segmenti, bridge_ninserts = 0;
  int endpoint, score;
  Univcoord_T suffix_univdiagonal, full_univdiagonal;

  if (full_path->plusp == false || suffix_path->plusp == false) {
    return (T) NULL;
  } else if (full_path->chrnum != suffix_path->chrnum ||
	     full_path->chroffset != suffix_path->chroffset ||
	     full_path->chrhigh != suffix_path->chrhigh ||
	     full_path->plusp != suffix_path->plusp ||
	     full_path->genestrand != suffix_path->genestrand) {
    return (T) NULL;
  } else if (full_path->fusion_querystart_junction != NULL ||
	     full_path->fusion_queryend_junction != NULL ||
	     suffix_path->fusion_querystart_junction != NULL ||
	     suffix_path->fusion_queryend_junction != NULL) {
    return (T) NULL;
  } else if (full_path->qstart_alts != NULL || full_path->qend_alts != NULL ||
	     suffix_path->qstart_alts != NULL || suffix_path->qend_alts != NULL) {
    return (T) NULL;
  } else if (Intlist_length(suffix_path->endpoints) < 3) {
    return (T) NULL;
  }

  suffix_qstart = Intlist_head(suffix_path->endpoints);
  suffix_qend = Intlist_last_value(suffix_path->endpoints);
  full_qend = Intlist_last_value(full_path->endpoints);
  if (suffix_qstart <= Intlist_head(full_path->endpoints) || suffix_qend != full_qend) {
    return (T) NULL;
  }

  suffix_univdiagonal = Univcoordlist_head(suffix_path->univdiagonals);
  segmenti = 0;
  for (fq = full_path->endpoints, fu = full_path->univdiagonals;
       fu != NULL && Intlist_next(fq) != NULL;
       fq = Intlist_next(fq), fu = Univcoordlist_next(fu), segmenti++) {
    endpoint = Intlist_head(Intlist_next(fq));
    full_univdiagonal = Univcoordlist_head(fu);

    if (endpoint < suffix_qstart &&
	full_univdiagonal > suffix_univdiagonal &&
	full_univdiagonal - suffix_univdiagonal <= 10 &&
	endpoint + (int) (full_univdiagonal - suffix_univdiagonal) == suffix_qstart) {
      cut_index = segmenti;
      bridge_ninserts = (int) (full_univdiagonal - suffix_univdiagonal);
    }
  }

  if (cut_index < 0) {
    return (T) NULL;
  }

  fq = full_path->endpoints;
  fu = full_path->univdiagonals;
  fs = full_path->nmismatches;
  fr = full_path->ref_nmismatches;
  fj = full_path->junctions;

  endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(fq)
			       intlistpool_trace(__FILE__,__LINE__));
  for (segmenti = 0; segmenti <= cut_index; segmenti++) {
    fq = Intlist_next(fq);
    endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(fq)
				 intlistpool_trace(__FILE__,__LINE__));
    univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(fu)
					   univcoordlistpool_trace(__FILE__,__LINE__));
    nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(fs)
				   intlistpool_trace(__FILE__,__LINE__));
    ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(fr)
				       intlistpool_trace(__FILE__,__LINE__));
    if (segmenti < cut_index) {
      junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(fj),pathpool)
				listpool_trace(__FILE__,__LINE__));
      fj = List_next(fj);
    }
    fu = Univcoordlist_next(fu);
    fs = Intlist_next(fs);
    fr = Intlist_next(fr);
  }

  junctions = Listpool_push(junctions,listpool,(void *) Junction_new_insertion(bridge_ninserts,pathpool)
			    listpool_trace(__FILE__,__LINE__));

  sq = Intlist_next(suffix_path->endpoints);
  su = suffix_path->univdiagonals;
  ss = suffix_path->nmismatches;
  sr = suffix_path->ref_nmismatches;
  sj = suffix_path->junctions;
  while (su != NULL) {
    endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(sq)
				 intlistpool_trace(__FILE__,__LINE__));
    univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(su)
					   univcoordlistpool_trace(__FILE__,__LINE__));
    nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(ss)
				   intlistpool_trace(__FILE__,__LINE__));
    ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(sr)
				       intlistpool_trace(__FILE__,__LINE__));

    su = Univcoordlist_next(su);
    ss = Intlist_next(ss);
    sr = Intlist_next(sr);
    sq = Intlist_next(sq);
    if (sj != NULL) {
      junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(sj),pathpool)
				listpool_trace(__FILE__,__LINE__));
      sj = List_next(sj);
    }
  }

  endpoints = Intlist_reverse(endpoints);
  univdiagonals = Univcoordlist_reverse(univdiagonals);
  nmismatches = Intlist_reverse(nmismatches);
  ref_nmismatches = Intlist_reverse(ref_nmismatches);
  junctions = List_reverse(junctions);

  if (Path_endpoints_acceptable_p(endpoints,junctions) == false) {
    Intlistpool_free_list(&endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Junction_list_gc(&junctions,listpool,pathpool);
    return (T) NULL;
  }

  repaired = Path_create(full_path->main_univdiagonal,
			 endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,
			 full_path->plusp,/*first_read_p*/true,full_path->genestrand,
			 full_path->sensedir,full_path->querylength,full_path->method,
			 full_path->chrnum,full_path->chroffset,full_path->chrhigh,
			 full_path->splice5p,full_path->splicetype5,full_path->ambig_prob_5,
			 full_path->splice3p,full_path->splicetype3,full_path->ambig_prob_3,
			 /*qstart_alts*/NULL,/*qend_alts*/NULL,
			 full_path->query_compress_fwd,full_path->query_compress_rev,
			 pathpool,/*vectorpool*/NULL);
  repaired->main_univdiagonal = full_path->main_univdiagonal;
  repaired->junction_splice_prob = full_path->junction_splice_prob;
  repaired->total_splice_prob = full_path->total_splice_prob;
  repaired->completep = full_path->completep;
  repaired->childp = full_path->childp;
  repaired->extendedp = full_path->extendedp;

  score = repaired->score_within_trims;
  repaired->nmatches = repaired->ref_nmatches =
    Path_coverage(repaired) - score - Junction_total_ninserts(repaired->junctions);

  return repaired;
}


static T
Path_merge_split_set (T *patharray, int npaths,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, Pathpool_T pathpool) {
  int *qstarts = NULL, *qends = NULL, *order = NULL, *bridge_kind = NULL, *bridge_length = NULL;
  Intlist_T endpoints = NULL, nmismatches = NULL, ref_nmismatches = NULL;
  Univcoordlist_T univdiagonals = NULL;
  List_T junctions = NULL;
  Intlist_T q, s, r;
  Univcoordlist_T u;
  List_T j;
  T first, prev, curr, merged = NULL;
  int i, k, temp_idx, idx, prev_idx, last_idx, qgap, score;
  long long prev_ref_end, curr_ref_start, ref_gap;

  if (npaths <= 1) {
    return (T) NULL;
  }

  first = patharray[0];
  if (first->fusion_querystart_junction != NULL ||
      first->fusion_queryend_junction != NULL ||
      first->circular_endpoints != NULL ||
      first->qstart_alts != NULL || first->qend_alts != NULL ||
      first->transcripts != NULL || first->invalid_transcripts != NULL ||
      first->fusion_transcripts != NULL || first->fusion_invalid_transcripts != NULL) {
    return (T) NULL;
  }

  qstarts = (int *) MALLOC(npaths * sizeof(int));
  qends = (int *) MALLOC(npaths * sizeof(int));
  order = (int *) MALLOC(npaths * sizeof(int));
  bridge_kind = (int *) MALLOC((npaths - 1) * sizeof(int));
  bridge_length = (int *) MALLOC((npaths - 1) * sizeof(int));

  for (i = 0; i < npaths; i++) {
    curr = patharray[i];
    if (curr->chrnum != first->chrnum ||
	curr->chroffset != first->chroffset ||
	curr->chrhigh != first->chrhigh ||
	curr->plusp != first->plusp ||
	curr->genestrand != first->genestrand ||
	curr->sensedir != first->sensedir ||
	curr->querylength != first->querylength ||
	curr->fusion_querystart_junction != NULL ||
	curr->fusion_queryend_junction != NULL ||
	curr->circular_endpoints != NULL ||
	curr->qstart_alts != NULL || curr->qend_alts != NULL ||
	curr->transcripts != NULL || curr->invalid_transcripts != NULL ||
	curr->fusion_transcripts != NULL || curr->fusion_invalid_transcripts != NULL) {
      FREE(bridge_length);
      FREE(bridge_kind);
      FREE(order);
      FREE(qends);
      FREE(qstarts);
      return (T) NULL;
    }

    qstarts[i] = Intlist_head(curr->endpoints);
    qends[i] = Intlist_last_value(curr->endpoints);
    order[i] = i;
  }

  for (i = 1; i < npaths; i++) {
    temp_idx = order[i];
    k = i - 1;
    while (k >= 0 &&
	   (qstarts[order[k]] > qstarts[temp_idx] ||
	    (qstarts[order[k]] == qstarts[temp_idx] && qends[order[k]] > qends[temp_idx]))) {
      order[k + 1] = order[k];
      k--;
    }
    order[k + 1] = temp_idx;
  }

  for (i = 1; i < npaths; i++) {
    prev_idx = order[i - 1];
    idx = order[i];
    prev = patharray[prev_idx];
    curr = patharray[idx];

    qgap = qstarts[idx] - qends[prev_idx];
    if (qgap < 0) {
      FREE(bridge_length);
      FREE(bridge_kind);
      FREE(order);
      FREE(qends);
      FREE(qstarts);
      return (T) NULL;
    }

    prev_ref_end = (long long) Univcoordlist_last_value(prev->univdiagonals) -
      (long long) prev->querylength + (long long) qends[prev_idx];
    curr_ref_start = (long long) Univcoordlist_head(curr->univdiagonals) -
      (long long) curr->querylength + (long long) qstarts[idx];
    ref_gap = curr_ref_start - prev_ref_end;

    if (qgap == 0 && ref_gap > 0 && ref_gap <= (long long) INT_MAX) {
      bridge_kind[i - 1] = SPLICE_JUNCTION;
      bridge_length[i - 1] = (int) ref_gap;
    } else if (qgap > 0 && ref_gap == 0) {
      bridge_kind[i - 1] = INS_JUNCTION;
      bridge_length[i - 1] = qgap;
    } else {
      FREE(bridge_length);
      FREE(bridge_kind);
      FREE(order);
      FREE(qends);
      FREE(qstarts);
      return (T) NULL;
    }
  }

  first = patharray[order[0]];
  idx = order[0];
  q = patharray[idx]->endpoints;
  endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(q)
			       intlistpool_trace(__FILE__,__LINE__));

  for (i = 0; i < npaths; i++) {
    idx = order[i];
    curr = patharray[idx];

    if (i > 0) {
      if (bridge_kind[i - 1] == INS_JUNCTION) {
	junctions = Listpool_push(junctions,listpool,
				  (void *) Junction_new_insertion(bridge_length[i - 1],pathpool)
				  listpool_trace(__FILE__,__LINE__));
      } else {
	junctions = Listpool_push(junctions,listpool,
				  (void *) Junction_new_splice((Chrpos_T) bridge_length[i - 1],
							       first->sensedir,0.0,0.0,pathpool)
				  listpool_trace(__FILE__,__LINE__));
      }
    }

    q = Intlist_next(curr->endpoints);
    u = curr->univdiagonals;
    s = curr->nmismatches;
    r = curr->ref_nmismatches;
    j = curr->junctions;
    while (u != NULL) {
      endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(q)
				   intlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,
					     Univcoordlist_head(u)
					     univcoordlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(s)
				     intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(r)
					 intlistpool_trace(__FILE__,__LINE__));

      u = Univcoordlist_next(u);
      s = Intlist_next(s);
      r = Intlist_next(r);
      q = Intlist_next(q);
      if (j != NULL) {
	junctions = Listpool_push(junctions,listpool,
				  (void *) Junction_copy((Junction_T) List_head(j),pathpool)
				  listpool_trace(__FILE__,__LINE__));
	j = List_next(j);
      }
    }
  }

  endpoints = Intlist_reverse(endpoints);
  univdiagonals = Univcoordlist_reverse(univdiagonals);
  nmismatches = Intlist_reverse(nmismatches);
  ref_nmismatches = Intlist_reverse(ref_nmismatches);
  junctions = List_reverse(junctions);
  last_idx = order[npaths - 1];

  FREE(bridge_length);
  FREE(bridge_kind);
  FREE(order);
  FREE(qends);
  FREE(qstarts);

  if (Path_endpoints_acceptable_p(endpoints,junctions) == false) {
    Intlistpool_free_list(&endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Junction_list_gc(&junctions,listpool,pathpool);
    return (T) NULL;
  }

  merged = Path_create(first->main_univdiagonal,
		       endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,
		       first->plusp,/*first_read_p*/true,first->genestrand,
		       first->sensedir,first->querylength,first->method,
		       first->chrnum,first->chroffset,first->chrhigh,
		       first->splice5p,first->splicetype5,first->ambig_prob_5,
		       patharray[last_idx]->splice3p,
		       patharray[last_idx]->splicetype3,
		       patharray[last_idx]->ambig_prob_3,
		       /*qstart_alts*/NULL,/*qend_alts*/NULL,
		       first->query_compress_fwd,first->query_compress_rev,
		       pathpool,/*vectorpool*/NULL);
  merged->main_univdiagonal = first->main_univdiagonal;
  merged->junction_splice_prob = first->junction_splice_prob;
  merged->total_splice_prob = first->total_splice_prob;
  merged->completep = true;
  merged->childp = first->childp;
  merged->extendedp = first->extendedp;

  score = merged->score_within_trims;
  merged->nmatches = merged->ref_nmatches =
    Path_coverage(merged) - score - Junction_total_ninserts(merged->junctions);

  return merged;
}


static int
Path_filter_standard_global (T *patharray, int npaths, int mincoverage_filter,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool,
			     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  int i, k;
  T path;

  k = 0;
  for (i = 0; i < npaths; i++) {
    path = patharray[i];
    if (Path_coverage(path) >= mincoverage_filter) {
      patharray[k++] = path;
    } else {
      debug8(printf("(6 split fallback) Eliminating sub-mincoverage ")); debug8(Path_print(path));
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
  }
  npaths = k;

  if (npaths == 0) {
    return 0;
  } else if (npaths > 1) {
    qsort(patharray,npaths,sizeof(T),Path_global_cmp_active);

    debug8(printf("(6 split fallback) Keeping best ")); debug8(Path_print(patharray[0]));
    i = 1;
    for (k = 1; k < npaths; k++) {
      path = patharray[k];
      if (Path_global_cmp_active(&path,&(patharray[0])) == 0) {
	debug8(printf("(6 split fallback tie) Keeping ")); debug8(Path_print(path));
	patharray[i++] = path;
      } else {
	debug8(printf("(6 split fallback worse) Eliminating ")); debug8(Path_print(path));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
    }

    return i;
  } else {
    return 1;
  }
}


static int
Path_filter_quality_split_global (T *patharray, int npaths, int split_path_max_overlap,
				  int split_path_min_segment, int split_path_min_combined_coverage,
				  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				  Listpool_T listpool, Pathpool_T pathpool,
				  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  bool *selected, *eligible;
  int *qstarts, *qends, *order, *prev;
  double *scores, *best_scores;
  T *selected_paths;
  bool movedp;
  int i, j, k, oi, oj, idx, jdx, temp_idx, selected_count, union_coverage;
  int best_oi = -1, best_single_i, gap;
  double score, repaired_score, candidate_score, best_set_score, best_single_score, best_single_margin_score;
  T path, repaired_path, merged_path;

  if (quality_split_ranker_p == false || npaths <= 1) {
    return 0;
  }

  debug8(printf("(split quality) Evaluating %d paths for a compatible split set\n",npaths));

  selected = (bool *) MALLOC(npaths * sizeof(bool));
  eligible = (bool *) MALLOC(npaths * sizeof(bool));
  qstarts = (int *) MALLOC(npaths * sizeof(int));
  qends = (int *) MALLOC(npaths * sizeof(int));
  order = (int *) MALLOC(npaths * sizeof(int));
  prev = (int *) MALLOC(npaths * sizeof(int));
  scores = (double *) MALLOC(npaths * sizeof(double));
  best_scores = (double *) MALLOC(npaths * sizeof(double));

  for (i = 0; i < npaths; i++) {
    Path_query_span(&qstarts[i],&qends[i],patharray[i]);
    scores[i] = Path_quality_score(patharray[i]);
    selected[i] = false;
    eligible[i] = (Path_coverage(patharray[i]) >= split_path_min_segment) ? true : false;
    order[i] = i;
  }

  /* Insertion sort by query end, then query start.  npaths is already small here. */
  for (i = 1; i < npaths; i++) {
    temp_idx = order[i];
    j = i - 1;
    while (j >= 0 &&
	   (qends[order[j]] > qends[temp_idx] ||
	    (qends[order[j]] == qends[temp_idx] && qstarts[order[j]] > qstarts[temp_idx]))) {
      order[j + 1] = order[j];
      j--;
    }
    order[j + 1] = temp_idx;
  }

  best_set_score = -HUGE_VAL;
  for (oi = 0; oi < npaths; oi++) {
    idx = order[oi];
    prev[oi] = -1;

    if (eligible[idx] == false) {
      best_scores[oi] = -HUGE_VAL;
      continue;
    } else {
      best_scores[oi] = scores[idx];
    }

    for (oj = 0; oj < oi; oj++) {
      jdx = order[oj];
      if (best_scores[oj] == -HUGE_VAL) {
	/* Previous path was too short for split-set ranking */
      } else if (qends[jdx] <= qstarts[idx] + split_path_max_overlap) {
	gap = qstarts[idx] - qends[jdx];
	if (gap < 0) {
	  gap = 0;
	}
	candidate_score = best_scores[oj] + scores[idx] -
	  quality_split_gap_penalty * (double) gap - quality_split_segment_penalty;
	if (candidate_score > best_scores[oi]) {
	  best_scores[oi] = candidate_score;
	  prev[oi] = oj;
	}
      }
    }

    if (best_oi < 0 || best_scores[oi] > best_set_score) {
      best_set_score = best_scores[oi];
      best_oi = oi;
    }
  }

  if (best_oi < 0) {
    FREE(best_scores);
    FREE(scores);
    FREE(prev);
    FREE(order);
    FREE(qends);
    FREE(qstarts);
    FREE(eligible);
    FREE(selected);
    return 0;
  }

  selected_count = 0;
  for (oi = best_oi; oi >= 0; oi = prev[oi]) {
    idx = order[oi];
    if (selected[idx] == false) {
      selected[idx] = true;
      selected_count++;
    }
    if (prev[oi] < 0) {
      break;
    }
  }

  if (selected_count <= 1) {
    FREE(best_scores);
    FREE(scores);
    FREE(prev);
    FREE(order);
    FREE(qends);
    FREE(qstarts);
    FREE(eligible);
    FREE(selected);
    return 0;
  }

  selected_paths = (T *) MALLOC(selected_count * sizeof(T));
  k = 0;
  for (i = 0; i < npaths; i++) {
    if (selected[i] == true) {
      selected_paths[k++] = patharray[i];
    }
  }
  union_coverage = Path_query_span_union_coverage(selected_paths,selected_count);
  FREE(selected_paths);

  best_single_i = 0;
  best_single_score = scores[0];
  for (i = 1; i < npaths; i++) {
    if (scores[i] > best_single_score) {
      best_single_i = i;
      best_single_score = scores[i];
    }
  }

  for (i = 0; i < npaths; i++) {
    if (selected[i] == true && i != best_single_i &&
	(repaired_path = Path_repair_terminal_suffix(patharray[best_single_i],patharray[i],
						     intlistpool,univcoordlistpool,
						     listpool,pathpool)) != NULL) {
      repaired_score = Path_quality_score(repaired_path);
      if (repaired_score > best_single_score) {
	path = patharray[best_single_i];
	debug8(printf("(split quality) Repaired best single with split suffix; old score %.3f, repaired score %.3f\n",
		      best_single_score,repaired_score));
	debug8(Path_print(repaired_path));
	patharray[best_single_i] = repaired_path;
	scores[best_single_i] = best_single_score = repaired_score;
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      } else {
	Path_free(&repaired_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
  }

  best_single_margin_score = best_single_score + fabs(best_single_score) * quality_rank_margin;
  if (union_coverage < split_path_min_combined_coverage ||
      best_set_score <= best_single_margin_score) {
    debug8(printf("(split quality) Candidate set score %.3f with %d paths and union coverage %d rejected against best single %.3f, margin threshold %.3f, and min coverage %d\n",
		  best_set_score,selected_count,union_coverage,best_single_score,
		  best_single_margin_score,split_path_min_combined_coverage));
    FREE(best_scores);
    FREE(scores);
    FREE(prev);
    FREE(order);
    FREE(qends);
    FREE(qstarts);
    FREE(eligible);
    FREE(selected);
    return 0;
  }

  debug8(printf("(split quality) Keeping set score %.3f with %d paths and union coverage %d >= %d over best single %.3f, margin threshold %.3f\n",
		best_set_score,selected_count,union_coverage,split_path_min_combined_coverage,
		best_single_score,best_single_margin_score));

  k = 0;
  for (i = 0; i < npaths; i++) {
    path = patharray[i];
    if (selected[i] == true) {
      score = Path_quality_score(path);
      Path_quality_debug(/*label*/"split",/*action*/"Keeping set member",score,path);
      patharray[k++] = path;
    } else {
      score = Path_quality_score(path);
      Path_quality_debug(/*label*/"split",/*action*/"Eliminating outside set",score,path);
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
  }

  if (k > 1 &&
      (merged_path = Path_merge_split_set(patharray,k,
					  intlistpool,univcoordlistpool,
					  listpool,pathpool)) != NULL) {
    debug8(printf("(split quality) Merged %d split paths into ",k)); debug8(Path_print(merged_path));
    for (i = 0; i < k; i++) {
      path = patharray[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    patharray[0] = merged_path;
    k = 1;
  }

  if (k > 1) {
    do {
      movedp = false;
      for (i = 1; i < k; i++) {
	if (Path_global_cmp_active(&(patharray[i]),&(patharray[i - 1])) < 0) {
	  path = patharray[i - 1];
	  patharray[i - 1] = patharray[i];
	  patharray[i] = path;
	  movedp = true;
	}
      }
    } while (movedp == true);
  }

  FREE(best_scores);
  FREE(scores);
  FREE(prev);
  FREE(order);
  FREE(qends);
  FREE(qstarts);
  FREE(eligible);
  FREE(selected);

  return k;
}


static int
Path_filter_split_global (T *patharray, int npaths, int mincoverage_filter,
			  int split_path_max_overlap, int split_path_min_segment,
			  int split_path_min_combined_coverage,
			  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool,
			  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  bool overlapping_kept_p, equivalent_overlap_p;
  int i, j, k, union_coverage;
  T path, merged_path;

  if ((k = Path_filter_quality_split_global(patharray,npaths,split_path_max_overlap,
					    split_path_min_segment,split_path_min_combined_coverage,
					    intlistpool,univcoordlistpool,listpool,pathpool,
					    transcriptpool,hitlistpool)) > 0) {
    return k;
  }

  debug8(printf("(6 split best) Keeping best ")); debug8(Path_print(patharray[0]));
  i = 1;

  for (k = 1; k < npaths; k++) {
    path = patharray[k];
    overlapping_kept_p = equivalent_overlap_p = false;

    for (j = 0; j < i; j++) {
      if (Path_query_span_overlap(path,patharray[j]) > split_path_max_overlap) {
	overlapping_kept_p = true;
	if (Path_global_cmp_active(&path,&(patharray[j])) == 0) {
	  equivalent_overlap_p = true;
	}
      }
    }

    if (overlapping_kept_p == false) {
      debug8(printf("(6 split complementary) Keeping ")); debug8(Path_print(path));
      patharray[i++] = path;
    } else if (equivalent_overlap_p == true) {
      debug8(printf("(6 split tie) Keeping overlapping tie ")); debug8(Path_print(path));
      patharray[i++] = path;
    } else {
      debug8(printf("(6 split worse) Eliminating overlapping worse ")); debug8(Path_print(path));
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
  }
  npaths = i;

  if (Path_split_complementary_p(patharray,npaths,split_path_max_overlap) == false) {
    debug8(printf("(6 split) No complementary query spans retained; falling back to standard coverage filter\n"));
    return Path_filter_standard_global(patharray,npaths,mincoverage_filter,
				       intlistpool,univcoordlistpool,listpool,pathpool,
				       transcriptpool,hitlistpool);
  }

  union_coverage = Path_query_span_union_coverage(patharray,npaths);
  if (union_coverage < split_path_min_combined_coverage) {
    debug8(printf("(6 split) Combined query span coverage %d < %d; falling back to standard coverage filter\n",
		  union_coverage,split_path_min_combined_coverage));
    return Path_filter_standard_global(patharray,npaths,mincoverage_filter,
				       intlistpool,univcoordlistpool,listpool,pathpool,
				       transcriptpool,hitlistpool);
  } else if ((merged_path = Path_merge_split_set(patharray,npaths,
						 intlistpool,univcoordlistpool,
						 listpool,pathpool)) != NULL) {
    debug8(printf("(6 split) Merged %d paths with combined query span coverage %d into ",
		  npaths,union_coverage)); debug8(Path_print(merged_path));
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    patharray[0] = merged_path;
    return 1;
  } else {
    debug8(printf("(6 split) Keeping %d paths with combined query span coverage %d >= %d\n",
		  npaths,union_coverage,split_path_min_combined_coverage));
    return npaths;
  }
}
  

static bool
Path_main_fusion_equal (T a, T b) {

  if (a->fusion_endpoints == NULL) {
    return false;
  } else if (b->fusion_endpoints == NULL) {
    return false;
  } else if (Intlist_equal(a->endpoints,b->fusion_endpoints) == false) {
    return false;
  } else if (Intlist_equal(a->fusion_endpoints,b->endpoints) == false) {
    return false;
  } else if (Univcoordlist_equal(a->univdiagonals,b->fusion_univdiagonals) == false) {
    return false;
  } else if (Univcoordlist_equal(a->fusion_univdiagonals,b->univdiagonals) == false) {
    return false;
  } else {
    return true;
  }
}


static bool
Path_fullspan_coverage_dominant_p (T path, int best_opposite_coverage) {

  if (path->qstart_alts != NULL || path->qend_alts != NULL) {
    return false;
  } else if (path->endpoints == NULL) {
    return false;
  } else if (Intlist_head(path->endpoints) != 0) {
    return false;
  } else if (Intlist_last_value(path->endpoints) != path->querylength) {
    return false;
  } else if (Path_coverage(path) > best_opposite_coverage + 20) {
    return true;
  } else {
    return false;
  }
}


T *
Path_eval_and_sort (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T *patharray, int npaths, 
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    char *queryuc_ptr, char *queryrc, char *quality_string,
		    int nmismatches_filter, int mincoverage_filter, int max_softclip_filter,
		    bool split_paths_p, int split_path_max_overlap,
		    int split_path_min_segment, int split_path_min_combined_coverage,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Hitlistpool_T hitlistpool, bool filterp) {
  int pathi;
  float maxlik, loglik;
  float total, qual;		/* For Bayesian mapq calculation */

  int randomi, i, j, k, l, group_end, group_kept_start, local_npaths, softclip;
  T temp, path;
  int coverage, best_coverage_sense, best_coverage_antisense;
  double best_avg_junction_sense_prob, best_avg_junction_antisense_prob,
    best_end_sense_prob, best_end_antisense_prob, prob;
  /* int best_querystart, best_queryend, querystart, queryend; */
  int chop_qstart, chop_qend, qstart, qend;
  int ignore_found_score = 0;
  int ignore_end_score = 0;

  bool allow_soft_clips_qstart_p, allow_soft_clips_qend_p;


  debug8(printf("Entered Path_eval_and_sort with patharray %p and npaths %d\n",patharray,npaths));
  if (npaths == 0) {
    /* Skip */
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;

  } else {
    /* Have already called Path_extend */
    
    /* 0.  Unalias circular alignments and trim chrbounds */
    k = 0;
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      if (circularp[path->chrnum] == true) {
	if (path->plusp == true) {
	  Path_trim_circular_unalias(path,query_compress_fwd,
				     query_compress_fwd,query_compress_rev);
	} else {
	  Path_trim_circular_unalias(path,query_compress_rev,
				     query_compress_fwd,query_compress_rev);
	}
      }
      if (Path_trim_chrbounds(&ignore_found_score,&ignore_end_score,
			      path,query_compress_fwd,query_compress_rev,
			      intlistpool,univcoordlistpool,listpool,pathpool,
			      transcriptpool) == true) {
	patharray[k++] = path;
      } else {
	debug8(printf("Eliminating due to chrbounds: ")); debug8(Path_print(path));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
    npaths = k;
    debug8(printf("After trimming, have %d paths\n",npaths));
    if (npaths == 0) {
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;
    }


    /* 1.  Sort by structure to remove duplicates */
    if (npaths > 1) {
      qsort(patharray,npaths,sizeof(T),Path_structure_cmp);

      k = 0;
      i = 0;
      while (i < npaths) {
	j = i + 1;
	while (j < npaths && Path_structure_cmp(&(patharray[j]),&(patharray[i])) == 0) {
	  j++;
	}
	debug8(printf("Found an identical group by structure (except sensedir) of %d paths => Re-sorting by method_cmp\n",j - i));
	
	qsort(&(patharray[i]),j - i,sizeof(T),Path_method_cmp);
	debug8(printf("(1) Keeping by method_cmp ")); debug8(Path_print(patharray[i]));
	patharray[k++] = patharray[i];
	
	for (l = i + 1; l < j; l++) {
	  debug8(printf("(1) Eliminating by method_cmp ")); debug8(Path_print(patharray[l]));
	  path = patharray[l];
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
	
	i = j;
      }
      npaths = k;
    }
    

    /* 2.  Sort by intervals to find best sensedir at a given locus */
    /* Do not reward unnecessary splices, so do not eliminate alignments with splice prob of 0.0 */
    if (npaths > 1) {
      qsort(patharray,npaths,sizeof(T),Path_interval_cmp);
    
      k = 0;
      i = 0;
      while (i < npaths) {
	j = i + 1;
	while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	  j++;
	}
	debug8(printf("Found an overlapping group of %d => choosing sense vs antisense\n",j - i));

	best_avg_junction_sense_prob = best_avg_junction_antisense_prob = 0.0;
	best_end_sense_prob = best_end_antisense_prob = 0.0;
	best_coverage_sense = best_coverage_antisense = 0;
	for (l = i; l < j; l++) {
	  path = patharray[l];
	  coverage = Path_coverage(path);
	  if (path->sensedir == SENSE_FORWARD) {
	    if (coverage > best_coverage_sense) {
	      best_coverage_sense = coverage;
	    }
	    if ((prob = Path_avg_junction_splice_prob(path)) > best_avg_junction_sense_prob) {
	      best_avg_junction_sense_prob = prob;
	    }
	    if ((prob = Path_end_splice_prob(path)) > best_end_sense_prob) {
	      best_end_sense_prob = prob;
	    }

	  } else {
	    if (coverage > best_coverage_antisense) {
	      best_coverage_antisense = coverage;
	    }
	    if ((prob = Path_avg_junction_splice_prob(path)) > best_avg_junction_antisense_prob) {
	      best_avg_junction_antisense_prob = prob;
	    }
	    if ((prob = Path_end_splice_prob(path)) > best_end_antisense_prob) {
	      best_end_antisense_prob = prob;
	    }
	  }
	}

	if (best_avg_junction_sense_prob > best_avg_junction_antisense_prob ||
	    (best_avg_junction_sense_prob == best_avg_junction_antisense_prob && best_end_sense_prob > best_end_antisense_prob)) {
	  /* Keep only sense */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    if (path->total_splice_prob == 0.0) {
	      debug8(printf("(2 sensedir) Keeping null ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else if (path->sensedir == SENSE_FORWARD) {
	      debug8(printf("(2 sensedir) Keeping sense ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else if (Path_fullspan_coverage_dominant_p(path,best_coverage_sense) == true) {
	      debug8(printf("(2 sensedir) Keeping coverage-dominant antisense ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else {
	      debug8(printf("(2 sensedir) Eliminating antisense ")); debug8(Path_print(path));
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }

	} else if (best_avg_junction_antisense_prob > best_avg_junction_sense_prob ||
		   (best_avg_junction_antisense_prob == best_avg_junction_sense_prob && best_end_antisense_prob > best_end_sense_prob)) {
	  /* Keep only antisense */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    if (path->total_splice_prob == 0.0) {
	      debug8(printf("(2 sensedir) Keeping null ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else if (path->sensedir == SENSE_ANTI) {
	      debug8(printf("(2 sensedir) Keeping anti ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else if (Path_fullspan_coverage_dominant_p(path,best_coverage_antisense) == true) {
	      debug8(printf("(2 sensedir) Keeping coverage-dominant sense ")); debug8(Path_print(path));
	      patharray[k++] = path;
	    } else {
	      debug8(printf("(2 sensedir) Eliminating sense ")); debug8(Path_print(path));
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }
	} else {
	  /* Keep both */
	  for (l = i; l < j; l++) {
	    path = patharray[l];
	    patharray[k++] = path;
	  }
	}

	i = j;
      }
      npaths = k;
    }


    /* 3.  Sort by intervals to best paths at a given locus */
    if (npaths > 1) {
      /* Should already be sorted */
      /* qsort(patharray,npaths,sizeof(T),Path_interval_cmp); */
    
      k = 0;
      i = 0;
      while (i < npaths) {
	chop_qstart = Intlist_head(patharray[i]->endpoints);
	chop_qend = Intlist_last_value(patharray[i]->endpoints);

	j = i + 1;
	while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	  if ((qstart = Intlist_head(patharray[j]->endpoints)) > chop_qstart) {
	    chop_qstart = qstart;
	  }
	  if ((qend = Intlist_last_value(patharray[j]->endpoints)) < chop_qend) {
	    chop_qend = qend;
	  }
	  j++;
	}
	group_end = j;
	debug8(printf("Found an overlapping group of %d => re-sorting by Path_local_cmp\n",j - i));

	for (l = i; l < j; l++) {
	  patharray[l]->chop_qstart = chop_qstart;
	  patharray[l]->chop_qend = chop_qend;
	}

	if (split_paths_p == true &&
	    (local_npaths = Path_filter_quality_split_global(&(patharray[i]),j - i,
							     split_path_max_overlap,split_path_min_segment,
							     split_path_min_combined_coverage,
							     intlistpool,univcoordlistpool,listpool,pathpool,
							     transcriptpool,hitlistpool)) > 0) {
	  j = i + local_npaths;
	} else {
	  j = i + Path_filter_quality_score(&(patharray[i]),j - i,/*label*/"3 local",
					    intlistpool,univcoordlistpool,listpool,pathpool,
					    transcriptpool,hitlistpool);
	}

	/* Keep the best ones in the overlapping group */
	qsort(&(patharray[i]),j - i,sizeof(T),Path_local_cmp_active);
	debug8(printf("(3 best) Keeping by local_cmp ")); debug8(Path_print(patharray[i]));
	group_kept_start = k;
	patharray[k++] = patharray[i];
      
	for (l = i + 1; l < j; l++) {
	  if (Path_structure_ignore_sense_cmp(&(patharray[l]),&(patharray[i])) == 0) {
	    debug8(printf("(3 identical) Eliminating by identical structure ")); debug8(Path_print(patharray[l]));
	    path = patharray[l];
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

	  } else if (Path_local_cmp_active(&(patharray[l]),&(patharray[i])) == 0) {
	    debug8(printf("(3 tie) Keeping by local_cmp ")); debug8(Path_print(patharray[l]));
	    patharray[k++] = patharray[l];

	  } else if (quality_split_ranker_p == true && split_paths_p == true &&
		     Path_coverage(patharray[l]) >= split_path_min_segment &&
		     Path_split_compatible_with_kept_p(patharray[l],patharray,group_kept_start,k,
						       split_path_max_overlap) == true) {
	    debug8(printf("(3 split complementary) Keeping for split-set ranking ")); debug8(Path_print(patharray[l]));
	    patharray[k++] = patharray[l];

	  } else {
	    debug8(printf("(3 worse) Eliminating by local_cmp ")); debug8(Path_print(patharray[l]));
	    path = patharray[l];
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

	i = group_end;
      }
      npaths = k;
    }


    /* 4.  Find best solution globally */
    /* TODO: Can make this O(n) rather than O(n*log n) */
    if (npaths > 1) {
      chop_qstart = Intlist_head(patharray[0]->endpoints);
      chop_qend = Intlist_last_value(patharray[0]->endpoints);

      for (i = 1; i < npaths; i++) {
	if ((qstart = Intlist_head(patharray[i]->endpoints)) > chop_qstart) {
	  chop_qstart = qstart;
	}
	if ((qend = Intlist_last_value(patharray[i]->endpoints)) < chop_qend) {
	  chop_qend = qend;
	}
      }

      for (i = 0; i < npaths; i++) {
	patharray[i]->chop_qstart = chop_qstart;
	patharray[i]->chop_qend = chop_qend;
      }

      if (split_paths_p == true &&
	  (k = Path_filter_quality_split_global(patharray,npaths,split_path_max_overlap,
						split_path_min_segment,split_path_min_combined_coverage,
						intlistpool,univcoordlistpool,listpool,pathpool,
						transcriptpool,hitlistpool)) > 0) {
	npaths = k;
      } else {
	npaths = Path_filter_quality_score(patharray,npaths,/*label*/"4 global",
					   intlistpool,univcoordlistpool,listpool,pathpool,
					   transcriptpool,hitlistpool);
	if (npaths > 1) {
	  qsort(patharray,npaths,sizeof(T),Path_global_cmp_active);
	}
      }
    }
    debug8(printf("Found the global best solution with %d nmatches\n",patharray[0]->nmatches));


    /* 5.  Filter the results */
    if (filterp == true) {
      k = 0;
      for (i = 0; i < npaths; i++) {
	path = patharray[i];
	softclip = Intlist_head(path->endpoints) + path->querylength - Intlist_last_value(path->endpoints);
	if (path->found_score - path->score_within_trims > softclip) {
	  softclip = path->found_score - path->score_within_trims;
	}
	if (path->querylength - (path->nmatches + path->score_within_trims) > softclip) {
	  softclip = path->querylength - (path->nmatches + path->score_within_trims);
	}
	
	if (path->plusp == true) {
	  allow_soft_clips_qstart_p = allow_soft_clips_start_p;
	  allow_soft_clips_qend_p = allow_soft_clips_end_p;
	} else {
	  allow_soft_clips_qstart_p = allow_soft_clips_end_p;
	  allow_soft_clips_qend_p = allow_soft_clips_start_p;
	}

	if (path->score_within_trims > nmismatches_filter) {
	  debug8(printf("(5 filter) Eliminating: solution has too many %d nmismatches\n",
			path->score_within_trims));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else if (max_softclip_filter >= 0 && softclip > max_softclip_filter) {
	  debug8(printf("(5 filter) Eliminating: solution has too many %d soft-clipped bases\n",
			softclip));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	  
	} else if (allow_soft_clips_qstart_p == false && Intlist_head(path->endpoints) != 0) {
	  debug8(printf("(5 filter) Eliminating: soft clips %d at qstart; %d at qend\n",
			Intlist_head(path->endpoints),path->querylength - Intlist_last_value(path->endpoints)));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else if (allow_soft_clips_qend_p == false && Intlist_last_value(path->endpoints) != path->querylength) {
	  debug8(printf("(5 filter) Eliminating: soft clips %d at qstart; %d at qend\n",
			Intlist_head(path->endpoints),path->querylength - Intlist_last_value(path->endpoints)));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else if (Path_coverage(path) < mincoverage_filter &&
		   (split_paths_p == false || Path_coverage(path) < split_path_min_segment)) {
	  debug8(printf("(5 filter) solution has too little %d coverage\n",
		    Path_coverage(path)));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else if (Path_coverage(path) < mincoverage_filter) {
	  debug8(printf("(5 split filter) Keeping partial path with coverage %d for combined evaluation\n",
			Path_coverage(path)));
	  patharray[k++] = path;

	} else {
	  debug8(printf("(5 filter) Keeping\n"));
	  patharray[k++] = path;
	}
      }

      npaths = k;
    }

    if (npaths == 0) {
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;
    }


    /* 6.  Find best solution globally */
    /* TODO: Can make this O(n) rather than O(n*log n) */

    if (npaths > 1) {
      chop_qstart = Intlist_head(patharray[0]->endpoints);
      chop_qend = Intlist_last_value(patharray[0]->endpoints);

      for (i = 1; i < npaths; i++) {
	if ((qstart = Intlist_head(patharray[i]->endpoints)) > chop_qstart) {
	  chop_qstart = qstart;
	}
	if ((qend = Intlist_last_value(patharray[i]->endpoints)) < chop_qend) {
	  chop_qend = qend;
	}
      }

      for (i = 0; i < npaths; i++) {
	patharray[i]->chop_qstart = chop_qstart;
	patharray[i]->chop_qend = chop_qend;
      }
      
      if (split_paths_p == true &&
	  (k = Path_filter_quality_split_global(patharray,npaths,split_path_max_overlap,
						split_path_min_segment,split_path_min_combined_coverage,
						intlistpool,univcoordlistpool,listpool,pathpool,
						transcriptpool,hitlistpool)) > 0) {
	npaths = k;
      } else {
	npaths = Path_filter_quality_score(patharray,npaths,/*label*/"6 global",
					   intlistpool,univcoordlistpool,listpool,pathpool,
					   transcriptpool,hitlistpool);
	if (npaths > 1) {
	  qsort(patharray,npaths,sizeof(T),Path_global_cmp_active);
	}

	debug8(printf("Found the global best solution across %d paths with %d nmatches\n",
		      npaths,patharray[0]->nmatches));

	if (split_paths_p == true) {
	  npaths = Path_filter_split_global(patharray,npaths,mincoverage_filter,
					    split_path_max_overlap,split_path_min_segment,
					    split_path_min_combined_coverage,
					    intlistpool,univcoordlistpool,listpool,pathpool,
					    transcriptpool,hitlistpool);
	} else {
	  /* Keep all paths equivalent to the best path */
	  debug8(printf("(6 best) Keeping best ")); debug8(Path_print(patharray[0]));
	  i = 1;			/* Skip the best path */

	  for (k = 1; k < npaths; k++) {
	    path = patharray[k];
	    if (Path_global_cmp_active(&path,&(patharray[0])) == 0) {
	      debug8(printf("(6 tie) Keeping ")); debug8(Path_print(path));
	      patharray[i++] = path;
	    } else {
	      debug8(printf("(6 worse) Eliminating ")); debug8(Path_print(path));
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }

	  npaths = i;
	}
      }
    }

    if (split_paths_p == true && filterp == true && npaths == 1 &&
	Path_coverage(patharray[0]) < mincoverage_filter) {
      debug8(printf("(6 split singleton) Eliminating partial path with coverage %d < %d\n",
		    Path_coverage(patharray[0]),mincoverage_filter));
      path = patharray[0];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      npaths = 0;
    }

    if (npaths == 0) {
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;
    }


    /* 7.  Check for main <=> fusion identity (can occur in single-end
       reads without a paired-end as anchor) */
    if (npaths > 1) {
      for (i = 0; i < npaths; i++) {
	if (patharray[i] != NULL) {
	  for (j = i + 1; j < npaths; j++) {
	    if (patharray[j] != NULL) {
	      if (Path_main_fusion_equal(patharray[j],patharray[i]) == true) {
		debug8(printf("(6) Eliminating ")); debug8(Path_print(patharray[j]));
		path = patharray[j];
		Path_free(&path,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
		patharray[j] = (Path_T) NULL;
	      }
	    }
	  }
	}
      }

      k = 0;
      for (i = 0; i < npaths; i++) {
	if (patharray[i] != NULL) {
	  patharray[k++] = patharray[i];
	}
      }
      npaths = k;
    }

    if (want_random_p && npaths > 1) {
      /* Randomize among best alignments */
      /* randomi = (int) ((double) i * rand()/((double) RAND_MAX + 1.0)); */
      randomi = (int) (rand() / (((double) RAND_MAX + 1.0) / (double) npaths));
      /* fprintf(stderr,"%d dups => random %d\n",i,randomi); */
      temp = patharray[0];
      patharray[0] = patharray[randomi];
      patharray[randomi] = temp;
    }

    /* Trim alignments at chromosomal bounds */
    for (i = 0; i < npaths; i++) {
      debug8(printf("(8) Trimming ")); debug8(Path_print(patharray[i]));
      path = patharray[i];
      if (circularp[path->chrnum] == false) {
#if 0
	/* Now done at the start of this procedure */
	debug8(printf("Chrnum %d is not circular\n",path->chrnum));
	Path_trim_qstart_trimbounds(path,intlistpool,univcoordlistpool,listpool,pathpool,
				    /*trimbounds*/path->chroffset);
	Path_trim_qend_trimbounds(path,intlistpool,univcoordlistpool,listpool,pathpool,
				  /*trimbounds*/path->chrhigh);
#endif

      } else if (output_type == STD_OUTPUT || output_type == M8_OUTPUT) {
	/* If output type is alignment or m8, then don't want to split up the parts */

      } else if (path->plusp == true) {
	Path_trim_circular(path,/*query_compress*/query_compress_fwd,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool);

      } else {
	Path_trim_circular(path,/*query_compress*/query_compress_rev,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool);
      }
    }

    /* Compute mapq_loglik */
    if (npaths == 1) {
      path = patharray[0];
      Path_mark_alignment(path,query_compress_fwd,queryuc_ptr,query_compress_rev,queryrc,
			  pathpool);
      path->mapq_loglik = MAPQ_MAXIMUM_SCORE;
      path->mapq_score = MAPQ_max_quality_score(quality_string,path->querylength);
      path->absmq_score = MAPQ_MAXIMUM_SCORE;

      *first_absmq = path->absmq_score;
      *second_absmq = 0;

    } else {
      for (i = 0; i < npaths; i++) {
	path = patharray[i];
	Path_mark_alignment(path,query_compress_fwd,queryuc_ptr,query_compress_rev,queryrc,
			    pathpool);
	path->mapq_loglik =
	  MAPQ_loglik_string(path->genomic_diff,quality_string,path->querylength,path->plusp);
      }
    }
    

    /* Enforce monotonicity */
    for (i = npaths - 1; i > 0; i--) {
      if (patharray[i-1]->mapq_loglik < patharray[i]->mapq_loglik) {
	patharray[i-1]->mapq_loglik = patharray[i]->mapq_loglik;
      }
    }
    maxlik = patharray[0]->mapq_loglik;
    
    /* Subtract maxlik to avoid underflow */
    for (i = 0; i < npaths; i++) {
      patharray[i]->mapq_loglik -= maxlik;
    }

    /* Compute absolute mapq */
    for (i = 0; i < npaths; i++) {
      loglik = patharray[i]->mapq_loglik + MAPQ_MAXIMUM_SCORE;
      if (loglik < 0.0) {
	loglik = 0.0;
      }
      patharray[i]->absmq_score = rint(loglik);
    }
    *first_absmq = patharray[0]->absmq_score;
    if (npaths == 1) {
      *second_absmq = 0;
    } else {
      *second_absmq = patharray[1]->absmq_score;
    }

    /* Compute Bayesian mapq */
    total = 0.0;
    for (i = 0; i < npaths; i++) {
      total += (patharray[i]->mapq_loglik = fasterexp(patharray[i]->mapq_loglik));
    }

    /* Obtain posterior probabilities of being true */
    for (i = 0; i < npaths; i++) {
      patharray[i]->mapq_loglik /= total;
    }

    /* Convert to Phred scores */
    for (i = 0; i < npaths; i++) {
      if ((qual = 1.0 - patharray[i]->mapq_loglik) < 2.5e-10 /* 10^-9.6 */) {
	patharray[i]->mapq_score = 40;
      } else {
	patharray[i]->mapq_score = rint(-10.0 * log10(qual));
      }
    }
  }

  debug8(printf("Exiting Path_eval_and_sort with %d paths\n",npaths));
  
  if (transcriptome != NULL) {
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      Transcript_velocity_single(path);
    }
  }

  /* Filter for chrsubset */
  /* Want to allow other alignments to be found before filtering */
  *npaths_primary = *npaths_altloc = 0;
  if (chrsubsetp == NULL) {
    for (pathi = 0; pathi < npaths; pathi++) {
      path = patharray[pathi];
      if (altlocp[path->chrnum] == true) {
	(*npaths_altloc) += 1;
      } else {
	(*npaths_primary) += 1;
      }
    }
    return patharray;

  } else {
    k = 0;
    for (pathi = 0; pathi < npaths; pathi++) {
      path = patharray[pathi];
      if (chrsubsetp[path->chrnum] == false) {
	/* Do not save this hit */
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	
      } else {
	/* Save this hit.  Re-use existing array */
	if (altlocp[path->chrnum] == true) {
	  (*npaths_altloc) += 1;
	} else {
	  (*npaths_primary) += 1;
	}
	patharray[k++] = path;
      }
    }

    if ((*npaths_primary) + (*npaths_altloc) == 0) {
      FREE_OUT(patharray);
      return (T *) NULL;
    } else {
      return patharray;
    }
  }
}


#if 0
static bool
trims_compatible_p (int common_trim_qstart, int common_trim_qend,
		    int qstart_trim, int qend_trim) {
  if (common_trim_qstart == 0 && common_trim_qend == 0) {
    return true;
  } else if (common_trim_qstart > 0 && common_trim_qend > 0) {
    /* Trimming on both ends by a single path */
    return true;
  } else if (qend_trim > 0 && common_trim_qstart > 0) {
    /* Trimming on a different end from the equivalence class */
    return false;
  } else if (qstart_trim > 0 && common_trim_qend > 0) {
    /* Trimming on a different end from the equivalence class */
    return false;
  } else {
    return true;
  }
}
#endif


#if 0
/* Now using Path_structure_cmp instead */

/* Algorithm: Find the lowest common univdiagonal, keeping track of
   qstart trims.  Continue through common univdiagonals, returning
   false if different endpoints are found, and then initializing qend
   trims.  If a common univdiagonal is found after that, return false.
   Otherwise, return true and report the trims. */

/* If we allow unequal ends, then we could trim one path
   unnecessarily, which could penalize it in the global_cmp
   comparison */

/* #define ALLOW_UNEQUAL_ENDS 1 */

static bool
Path_common_structure_trim_p (int *common_trim_qstart, int *common_trim_qend, T a, T b) {
  bool common_univdiagonal_p = false, same_univdiagonal_p;
  int a_leading_univdiagonals = 0, a_trailing_univdiagonals = 0,
    b_leading_univdiagonals = 0, b_trailing_univdiagonals = 0;
  bool a_common_longest_p = false, b_common_longest_p = false;
  Univcoordlist_T p0, p1;
  Intlist_T q0, q1;
  Univcoord_T a_longest_univdiagonal, b_longest_univdiagonal;
  Univcoord_T a_univdiagonal, b_univdiagonal;
  int a_qstart, a_qend, b_qstart, b_qend;
  int maxlength, length;

  int a_trim_qstart = 0, a_trim_qend = 0;
  int b_trim_qstart = 0, b_trim_qend = 0;

#ifdef DEBUG5
  printf("Entered Path_common_structure_p\n");
  Path_print(a);
  Path_print(b);
#endif

  /* Find longest univdiagonals on each path */
  maxlength = 0;
  p0 = a->univdiagonals;
  q0 = a->endpoints;
  while (p0 != NULL) {
    length = Intlist_second_value(q0) - Intlist_head(q0);
    if (length > maxlength) {
      a_longest_univdiagonal = Univcoordlist_head(p0);
      maxlength = length;
    }
    p0 = Univcoordlist_next(p0);
    q0 = Intlist_next(q0);
  }

  maxlength = 0;
  p1 = b->univdiagonals;
  q1 = b->endpoints;
  while (p1 != NULL) {
    length = Intlist_second_value(q1) - Intlist_head(q1);
    if (length > maxlength) {
      b_longest_univdiagonal = Univcoordlist_head(p1);
      maxlength = length;
    }
    p1 = Univcoordlist_next(p1);
    q1 = Intlist_next(q1);
  }


  /* Find common structure */
  p0 = a->univdiagonals;
  p1 = b->univdiagonals;
  q0 = a->endpoints;
  q1 = b->endpoints;

  a_qend = Intlist_second_value(q0);
  b_qend = Intlist_second_value(q1);

  /* Find qstart trims until the first common univdiagonal */
  while (p0 != NULL && p1 != NULL && common_univdiagonal_p == false) {
    a_univdiagonal = Univcoordlist_head(p0);
    b_univdiagonal = Univcoordlist_head(p1);

    if (a_univdiagonal < b_univdiagonal) {
      debug5(printf("pre: a segment is entirely before b segment\n"));
      a_trim_qstart += Intlist_second_value(q0) - Intlist_head(q0);
      p0 = Univcoordlist_next(p0);
      q0 = Intlist_next(q0);
      a_leading_univdiagonals += 1;

    } else if (b_univdiagonal < a_univdiagonal) {
      debug5(printf("pre: b segment is entirely before a segment\n"));
      b_trim_qstart += Intlist_second_value(q1) - Intlist_head(q1);
      p1 = Univcoordlist_next(p1);
      q1 = Intlist_next(q1);
      b_leading_univdiagonals += 1;

    } else {
      /* Found first common univdiagonal */
      a_qstart = Intlist_head(q0);
      a_qend = Intlist_second_value(q0);
      b_qstart = Intlist_head(q1);
      b_qend = Intlist_second_value(q1);
      debug5(printf("pre: %u %d..%d vs %u %d..%d, with the same univdiagonal\n",
		    a_univdiagonal,a_qstart,a_qend,b_univdiagonal,b_qstart,b_qend));
      
      if (a_qstart < b_qstart) {
	debug5(printf("pre: trimming a only\n"));
#ifndef ALLOW_UNEQUAL_ENDS
	debug5(printf("Requiring equal ends, so returning false\n"));
	return false;
#else
	a_trim_qstart += (b_qstart - a_qstart);
#endif	  
      } else if (b_qstart < a_qstart) {
	debug5(printf("pre: trimming b only\n"));
#ifndef ALLOW_UNEQUAL_ENDS
	debug5(printf("Requiring equal ends, so returning false\n"));
	return false;
#else
	b_trim_qstart += (a_qstart - b_qstart);
#endif
      } else {
	debug5(printf("pre: no trimming\n"));
      }
      common_univdiagonal_p = true;
      if (a_univdiagonal == a_longest_univdiagonal) {
	a_common_longest_p = true;
      }
      if (b_univdiagonal == b_longest_univdiagonal) {
	b_common_longest_p = true;
      }

      p0 = Univcoordlist_next(p0);
      q0 = Intlist_next(q0);
      p1 = Univcoordlist_next(p1);
      q1 = Intlist_next(q1);
    }
  }

  debug5(printf("done with pre: a univdiagonals %p, b univdiagonals %p, common_univdiagonal_p %d, a_trim_qstart %d, b_trim_qstart %d\n",
		p0,p1,common_univdiagonal_p,a_trim_qstart,b_trim_qstart));

  /* Continue through common internal structure */
  same_univdiagonal_p = true;
  while (p0 != NULL && p1 != NULL && same_univdiagonal_p == true) {
    a_univdiagonal = Univcoordlist_head(p0);
    b_univdiagonal = Univcoordlist_head(p1);

    if (a_univdiagonal < b_univdiagonal) {
      debug5(printf("mid: a segment is entirely before b segment\n"));
      a_trim_qend += Intlist_second_value(q0) - Intlist_head(q0);
      p0 = Univcoordlist_next(p0);
      q0 = Intlist_next(q0);
      same_univdiagonal_p = false;

    } else if (b_univdiagonal < a_univdiagonal) {
      debug5(printf("mid: b segment is entirely before a segment\n"));
      b_trim_qend += Intlist_second_value(q1) - Intlist_head(q1);
      p1 = Univcoordlist_next(p1);
      q1 = Intlist_next(q1);
      same_univdiagonal_p = false;

    } else if (a_qend != b_qend) {
      /* Another common univdiagonal but previous one has different ends */
      debug5(printf("mid: Previous univdiagonal had different ends, so returning false\n"));
      return false;

    } else if ((a_qstart = Intlist_head(q0)) != (b_qstart = Intlist_head(q1))) {
      /* Another common univdiagonal but this one has different starts */
      debug5(printf("mid: This univdiagonal had different starts, so returning false\n"));
      return false;
      
    } else {
      a_qend = Intlist_second_value(q0);
      b_qend = Intlist_second_value(q1);
      debug5(printf("mid: %u %d..%d vs %d..%d\n",a_univdiagonal,a_qstart,a_qend,b_qstart,b_qend));

      if (a_univdiagonal == a_longest_univdiagonal) {
	a_common_longest_p = true;
      }
      if (b_univdiagonal == b_longest_univdiagonal) {
	b_common_longest_p = true;
      }

      p0 = Univcoordlist_next(p0);
      q0 = Intlist_next(q0);
      p1 = Univcoordlist_next(p1);
      q1 = Intlist_next(q1);
    }
  }

  /* Handle the final common univdiagonal */
  if (a_qend > b_qend) {
    debug5(printf("mid: trimming a only\n"));
#ifndef ALLOW_UNEQUAL_ENDS
    debug5(printf("Requiring equal ends, so returning false\n"));
    return false;
#else
    a_trim_qend += (a_qend - b_qend);
#endif
	
  } else if (b_qend > a_qend) {
    debug5(printf("mid: trimming b only\n"));
#ifndef ALLOW_UNEQUAL_ENDS
    debug5(printf("Requiring equal ends, so returning false\n"));
    return false;
#else
    b_trim_qend += (b_qend - a_qend);
#endif

  } else {
    debug5(printf("mid: no trimming\n"));
  }
    
  debug5(printf("done with mid: a univdiagonals %p, b univdiagonals %p, same_univdiagonal_p %d\n",
		p0,p1,same_univdiagonal_p));


  /* Continue after the common internal structure */
  while (p0 != NULL && p1 != NULL) {
    a_univdiagonal = Univcoordlist_head(p0);
    b_univdiagonal = Univcoordlist_head(p1);

    if (a_univdiagonal < b_univdiagonal) {
      debug5(printf("post: a segment is entirely before b segment\n"));
      a_trim_qend += Intlist_second_value(q0) - Intlist_head(q0);
      p0 = Univcoordlist_next(p0);
      q0 = Intlist_next(q0);
      a_trailing_univdiagonals += 1;

    } else if (b_univdiagonal < a_univdiagonal) {
      debug5(printf("post: b segment is entirely before a segment\n"));
      b_trim_qend += Intlist_second_value(q1) - Intlist_head(q1);
      p1 = Univcoordlist_next(p1);
      q1 = Intlist_next(q1);
      b_trailing_univdiagonals += 1;

    } else {
      /* Found a common univdiagonal after different univdiagonals */
      debug5(printf("post: found a common univdiagonal after different ones, so returning false\n"));
      return false;
    }
  }

  while (p0 != NULL) {
    debug5(printf("post: have remaining segments for a\n"));
    a_trim_qend += Intlist_second_value(q0) - Intlist_head(q0);
    p0 = Univcoordlist_next(p0);
    q0 = Intlist_next(q0);
    a_trailing_univdiagonals += 1;
  }

  while (p1 != NULL) {
    debug5(printf("post: have remaining segments for b\n"));
    b_trim_qend += Intlist_second_value(q1) - Intlist_head(q1);
    p1 = Univcoordlist_next(p1);
    q1 = Intlist_next(q1);
    b_trailing_univdiagonals += 1;
  }
      
  debug5(printf("end: Have trims %d vs %d and %d vs %d\n",
		a_trim_qstart,b_trim_qstart,a_trim_qend,b_trim_qend));
  debug5(printf("leading %d and %d.  trailing %d and %d\n",
		a_leading_univdiagonals,b_leading_univdiagonals,
		a_trailing_univdiagonals,b_trailing_univdiagonals));

  if (common_univdiagonal_p == false) {
    debug5(printf("No common univdiagonal, so returning false\n"));
    return false;
  } else if (a_common_longest_p == false) {
    debug5(printf("Common univdiagonal for path a is not the longest, so returning false\n"));
    return false;
  } else if (b_common_longest_p == false) {
    debug5(printf("Common univdiagonal for path b is not the longest, so returning false\n"));
    return false;
  } else if (a_trim_qstart != b_trim_qstart) {
    debug5(printf("qstarts unequal, so returning false\n"));
    return false;
  } else if (a_trim_qend != b_trim_qend) {
    debug5(printf("qends unequal, so returning false\n"));
    return false;
  } else if (a_leading_univdiagonals > 1 || b_leading_univdiagonals > 1) {
    debug5(printf("too many leading univdiagonals\n"));
    return false;
  } else if (a_trailing_univdiagonals > 1 || b_trailing_univdiagonals > 1) {
    debug5(printf("too many trailing univdiagonals\n"));
    return false;
  } else if (a_leading_univdiagonals == 1 && b_trailing_univdiagonals == 1) {
    debug5(printf("mismatch leading/trailing univdiagonals\n"));
    return false;
  } else if (b_leading_univdiagonals == 1 && a_trailing_univdiagonals == 1) {
    debug5(printf("mismatch leading/trailing univdiagonals\n"));
    return false;
  } else if (a_trim_qstart + a_trim_qend >= Path_coverage(a)) {
    debug5(printf("Trims are too long, so returning false\n"));
    return false;
  } else {
    *common_trim_qstart = a_trim_qstart;
    *common_trim_qend = a_trim_qend;
    debug5(printf("Returning true with trims %d and %d\n",*common_trim_qstart,*common_trim_qend));
    return true;
  }
}
#endif



/* Allowing invalid transcripts but remapping them to the trimmed
   transcript and moving them to valid if appropriate */
/* We have to allow invalid transcripts because Trpath_convert
   procedures call Path_solve_junctions on unsolved paths, resulting
   in invalid transcripts */
#define ALLOW_INVALID_TRANSCRIPTS 1


static List_T
find_equivalence_class (T path, List_T equiv_classes) {
  List_T p = equiv_classes;
  T representative_path;
  
  while (p != NULL) {
    representative_path = (T) List_head((List_T) List_head(p));
    if (path->transcriptome_method_p == representative_path->transcriptome_method_p &&
	Path_structure_cmp(&path,&representative_path) == 0) {
      return p;
    } else {
      p = List_next(p);
    }
  }

  return (List_T) NULL;
}

static bool
path_group_has_transcriptome_method_p (T *patharray, int lowi, int highi) {
  int i;

  for (i = lowi; i < highi; i++) {
    if (patharray[i]->transcriptome_method_p == true) {
      return true;
    }
  }

  return false;
}


List_T
Path_consolidate (List_T paths, Shortread_T queryseq,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		  Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool,
		  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  T *patharray;
  List_T class, c;
  List_T equiv_classes, p;
  /* Intlist_T qstart_trims, qend_trims, q, r; */
  T bestpath, path;
  /* int qstart_trim, qend_trim, common_trim_qstart, common_trim_qend; */
  int npaths, i, j, k;
  int best_nmatches;
  int nbest;
  bool transcriptome_group_p;

  List_T invalid_transcripts, t;
  Transcript_T transcript;
  bool assignedp;


  debug9(printf("Entered Path_consolidate with %d paths\n",List_length(paths)));
  if ((npaths = List_length(paths)) > 1) {
    patharray = (T *) List_to_array(paths,NULL);
    paths = (List_T) NULL;

    qsort(patharray,npaths,sizeof(T),Path_interval_cmp);
    
    i = 0;
    while (i < npaths) {
      /* Consider only paths with the most nmatches.  Use a linear pass instead of sorting */
      best_nmatches = patharray[i]->nmatches;
      bestpath = patharray[i];
      nbest = 1;

      j = i + 1;
      while (j < npaths && Path_overlap_p(patharray[j],patharray[i]) == true) {
	if (patharray[j]->nmatches > best_nmatches) {
	  bestpath = patharray[j];
	  nbest = 1;
	} else if (patharray[j]->nmatches == best_nmatches) {
	  nbest++;
	}
	j++;
      }
      debug9(printf("Found an overlapping group of %d starting at %d, best nmatches %d with %d paths\n",
		    j - i,i,best_nmatches,nbest));

      transcriptome_group_p = path_group_has_transcriptome_method_p(patharray,i,j);
      
      /* Handle the best paths in this overlapping group */
      if (nbest == 1 && transcriptome_group_p == false) {
	/* Only a single best path */
	debug9(printf("Keeping by best_nmatches ")); debug9(Path_print(bestpath));
	paths = Hitlist_push(paths,hitlistpool,(void *) bestpath
			     hitlistpool_trace(__FILE__,__LINE__));

	/* Remap path here if transcriptome_method_p is false? */

	for (k = i; k < j; k++) {
	  if (patharray[k]->nmatches < best_nmatches) {
	    path = patharray[k];
	    debug9(printf("(1) Eliminating by best_nmatches ")); debug9(Path_print(path));
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

      } else {
	/* Among the nbest paths, find equivalence classes based on
	   Path_structure_cmp (previously used Path_common_structure_trim_p,
	   which disregads trims at the end) */
	equiv_classes = (List_T) NULL;
	/* qstart_trims = (Intlist_T) NULL; */
	/* qend_trims = (Intlist_T) NULL; */

	for (k = i; k < j; k++) {
	  path = patharray[k];
	  if (path->nmatches < best_nmatches &&
	      (transcriptome_group_p == false || path->transcriptome_method_p == false)) {
	    debug9(printf("(2) Eliminating by best_nmatches ")); debug9(Path_print(path));
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

	  } else if ((p = find_equivalence_class(path,equiv_classes)) == NULL) {
	    /* Create new equivalence class */
	    debug9(printf("Creating new equivalence class\n"));
	    equiv_classes = Listpool_push(equiv_classes,listpool,
					  (void *) Listpool_push(NULL,listpool,(void *) path
								 listpool_trace(__FILE__,__LINE__))
					  listpool_trace(__FILE__,__LINE__));

	    /* qstart_trims = Intlistpool_push(qstart_trims,intlistpool,0
	       intlistpool_trace(__FILE__,__LINE__)); */
	    /* qend_trims = Intlistpool_push(qend_trims,intlistpool,0
	       intlistpool_trace(__FILE__,__LINE__)); */

	  } else {
	    /* Insert path into equivalence class and revise trims */
	    debug9(printf("Putting path existing equivalence class\n"));
	    
	    class = (List_T) List_head(p);
	    class = Listpool_push(class,listpool,(void *) path
				  listpool_trace(__FILE__,__LINE__));
	    List_head_set(p,(void *) class);
#if 0
	    debug9(printf("Revising trims from %d and %d",Intlist_head(q),Intlist_head(r)));
	    if (common_trim_qstart > Intlist_head(q)) {
	      Intlist_head_set(q,common_trim_qstart);
	    }
	    if (common_trim_qend > Intlist_head(r)) {
	      Intlist_head_set(r,common_trim_qend);
	    }
	    debug9(printf(" to %d and %d\n",Intlist_head(q),Intlist_head(r)));
#endif
	  }
	}
	      
	debug9(printf("Found %d equivalence classes\n",List_length(equiv_classes)));
	for (p = equiv_classes /*, q = qstart_trims, r = qend_trims*/;
	     p != NULL; p = List_next(p) /*, q = Intlist_next(q), r = Intlist_next(r)*/) {
	  /* Process each equivalence class */
	  class = (List_T) List_head(p);

#if 0
	  /* Not needed, since we have separate equivalence classes for transcript and non-transcript methods */
	  /* Find a representative from a transcriptome method */
	  /* non_transcriptomep = false; */
	  bestpath = (T) NULL;
	  for (c = class; c != NULL; c = List_next(c)) {
	    path = (T) List_head(c);
	    if (path->transcriptome_method_p == true) {
	      bestpath = path;
	    } else {
	      /* non_transcriptomep = true; */
	    }
	  }
	  if (bestpath == NULL) {
	    bestpath = (T) List_head(class);
	  }
#else
	  bestpath = (T) List_head(class);
	  for (c = List_next(class); c != NULL; c = List_next(c)) {
	    path = (T) List_head(c);
	    if (Path_method_cmp(&path,&bestpath) < 0) {
	      bestpath = path;
	    }
	  }
	  debug9(printf("Best path: ")); debug9(Path_print(bestpath));
#endif
	  for (c = class; c != NULL; c = List_next(c)) {
	    path = (T) List_head(c);
	    if (path != bestpath) {
	      debug9(printf("Eliminating ")); debug9(Path_print(path));
	      bestpath->transcripts = List_append(path->transcripts,bestpath->transcripts);
	      path->transcripts = (List_T) NULL;
#ifdef ALLOW_INVALID_TRANSCRIPTS
	      /* If we allow invalid transcripts from transcriptome methods */
	      bestpath->invalid_transcripts = List_append(path->invalid_transcripts,bestpath->invalid_transcripts);
	      path->invalid_transcripts = (List_T) NULL;
#endif
	      Path_free(&path,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	    }
	  }
	  
#if 0
	  qstart_trim = Intlist_head(q);
	  qend_trim = Intlist_head(r);
	  debug9(printf("qstart_trim %d, qend_trim %d\n",qstart_trim,qend_trim));
	  if (qstart_trim > 0) {
	    Path_trim_qstart_n(qstart_trim,bestpath,
			       query_compress_fwd,query_compress_rev,
			       intlistpool,univcoordlistpool,listpool,pathpool);
	  }
	  if (qend_trim > 0) {
	    Path_trim_qend_n(qend_trim,bestpath,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,univcoordlistpool,listpool,pathpool);
	  }
	  if (bestpath->genestrand > 0) {
	    Transcript_list_trim(bestpath->transcripts,/*trim5*/qstart_trim,/*trim3*/qend_trim,
				 transcriptome,listpool,transcriptpool);
	  } else {
	    Transcript_list_trim(bestpath->transcripts,/*trim5*/qend_trim,/*trim3*/qstart_trim,
				 transcriptome,listpool,transcriptpool);
	  }
#endif

#ifdef ALLOW_INVALID_TRANSCRIPTS
	  /* Hard to trim invalid transcripts, so remap them */
	  invalid_transcripts = bestpath->invalid_transcripts;
	  bestpath->invalid_transcripts = (List_T) NULL;
	  for (t = invalid_transcripts; t != NULL; t = List_next(t)) {
	    transcript = (Transcript_T) List_head(t);
	    if ((assignedp = Transcript_remap_invalid(transcript,bestpath,transcriptome,queryseq,
						      uintlistpool,listpool,transcriptpool)) == true) {
	      Transcript_free(&transcript,listpool,transcriptpool);
	    } else {
	      bestpath->invalid_transcripts = Listpool_push(bestpath->invalid_transcripts,listpool,(void *) transcript
							    listpool_trace(__FILE__,__LINE__));
	    }
	  }
	  Listpool_free_list(&invalid_transcripts,listpool
			     listpool_trace(__FILE__,__LINE__));
#endif

	  bestpath->transcripts = Transcript_list_sort(bestpath->transcripts,listpool,transcriptpool);
	  bestpath->invalid_transcripts = Transcript_list_sort(bestpath->invalid_transcripts,listpool,transcriptpool);

	  debug9(printf("Keeping representative of class ")); debug9(Path_print(bestpath));
	  paths = Hitlist_push(paths,hitlistpool,(void *) bestpath
			       hitlistpool_trace(__FILE__,__LINE__));
	  Listpool_free_list(&class,listpool
			     listpool_trace(__FILE__,__LINE__));
	}

#if 0
	Intlistpool_free_list(&qend_trims,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&qstart_trims,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));			      
#endif
	Listpool_free_list(&equiv_classes,listpool
			   listpool_trace(__FILE__,__LINE__));
      }

      i = j;
    }
	  
    FREE(patharray);
  }

#ifdef DEBUG9
  printf("Exiting Path_consolidate with %d paths:\n",List_length(paths));
  for (p = paths; p != NULL; p = List_next(p)) {
    Path_print((T) List_head(p));
  }
#endif

  return paths;
}



void
Path_eval_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		 Transcriptome_T transcriptome_in,
		 bool *circularp_in, bool *chrsubsetp_in, bool *altlocp_in,
		 int index1part_in, int index1interval_in,
		 Outputtype_T output_type_in, bool md_report_snps_p_in,
		 bool want_random_p_in, bool allow_soft_clips_start_p_in,
		 bool allow_soft_clips_end_p_in,
		 bool quality_ranker_p_in, double quality_match_reward_in,
		 double quality_mismatch_penalty_in,
		 double quality_indel_penalty_in,
		 double quality_softclip_penalty_in,
		 double quality_error_rate_penalty_in,
		 double quality_rank_margin_in,
		 bool quality_split_ranker_p_in,
		 double quality_split_gap_penalty_in,
		 double quality_split_segment_penalty_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  transcriptome = transcriptome_in;

  circularp = circularp_in;
  chrsubsetp = chrsubsetp_in;
  altlocp = altlocp_in;

  index1part = index1part_in;
  index1interval = index1interval_in;

  output_type = output_type_in;
  md_report_snps_p = md_report_snps_p_in;
  want_random_p = want_random_p_in;

  allow_soft_clips_start_p = allow_soft_clips_start_p_in;
  allow_soft_clips_end_p = allow_soft_clips_end_p_in;

  quality_ranker_p = quality_ranker_p_in;
  quality_match_reward = quality_match_reward_in;
  quality_mismatch_penalty = quality_mismatch_penalty_in;
  quality_indel_penalty = quality_indel_penalty_in;
  quality_softclip_penalty = quality_softclip_penalty_in;
  quality_error_rate_penalty = quality_error_rate_penalty_in;
  quality_rank_margin = quality_rank_margin_in;
  quality_split_ranker_p = quality_split_ranker_p_in;
  quality_split_gap_penalty = quality_split_gap_penalty_in;
  quality_split_segment_penalty = quality_split_segment_penalty_in;

  end_ambig_allowed = index1part_in;

  return;
}
