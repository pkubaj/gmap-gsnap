static char rcsid[] = "$Id: e7f61e0ccbcbcbee8c6dfc17d5225a1c092f4353 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-eval.h"
#include "path-solve.h"
#include "path-trim.h"
#include "transcript-remap.h"

#include <stdio.h>
#include <string.h>		/* For strcpy */
#include <math.h>		/* For rint */
#include <ctype.h>		/* For islower */
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
static bool allow_soft_clips_p;


/* Splice goodness now determined by Splice_accept_p */
/* #define OMIT_BAD_SPLICES 1 */

#define NMATCHES_FACTOR 8
#define CHOPPED_SLOP 1


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
#ifdef DEBUG8
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
Path_eval_nmatches (int *found_score, T this, Compress_T query_compress_fwd, Compress_T query_compress_rev) {
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


  debug7(printf("\nEntering Path_eval_nmatches on path %p\n",this));
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
  nmismatches = Intlist_head(this->nmismatches);
  ref_nmismatches = Intlist_head(this->ref_nmismatches);
  ninserts = 0;
  
  if (this->plusp == true) {
    /* plus */
    j = this->junctions;		/* Put here before we handle querystart_alts */
    if (this->ambig_prob_5 > 0.9) {
      /* Count ambig splice site as matches, but counts against found_score */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      this->found_score += qstart;
      /* Does not affect score_within_trims */

    } else if ((altsplice = this->qstart_alts) != NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      this->nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("Adding %d matches for qstart_alts => total %d\n",
		    Altsplice_best_nmatches(this->qstart_alts),this->nmatches));
      this->ref_nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->total_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->found_score += Altsplice_found_score(this->qstart_alts);
      this->score_within_trims += Altsplice_found_score(this->qstart_alts);

    } else if (this->fusion_querystart_junction != NULL) {
      /* Rest of fusions handled below */

    } else {
      this->found_score += qstart;
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
      /* Count ambig splice site as matches, but counts against found_score */
      this->total_splice_prob += this->ambig_prob_5;
      /* this->nmatches += qstart; (inferred) */
      /* this->ref_nmatches += qstart; (inferred) */
      this->found_score += qstart;
      /* Does not affect score_within_trims */

    } else if ((altsplice = this->qstart_alts )!= NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
      this->nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      debug7(printf("Adding %d matches for qstart_alts => total %d\n",
		    Altsplice_best_nmatches(this->qstart_alts),this->nmatches));
      this->ref_nmatches += Altsplice_best_nmatches(this->qstart_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
      this->junction_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->total_splice_prob += this->qstart_alts->medial_prob + Altsplice_best_distal_prob(this->qstart_alts);
      this->found_score += Altsplice_found_score(this->qstart_alts);
      this->score_within_trims += Altsplice_found_score(this->qstart_alts);

    } else if (this->fusion_querystart_junction != NULL) {
      /* Fusions handled below */

    } else {
      this->found_score += qstart;
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
    /* Count ambig splice site as matches, but counts against found_score */
    this->total_splice_prob += this->ambig_prob_3;
    /* this->nmatches += this->querylength - qend; (inferred) */
    /* this->ref_nmatches += this->querylength - qend; (inferred) */
    this->found_score += this->querylength - qend;
    /* Does not affect score_within_trims */

  } else if ((altsplice = this->qend_alts) != NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
    this->nmatches += Altsplice_best_nmatches(this->qend_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
    debug7(printf("Adding %d matches for qend_alts => total %d\n",
		  Altsplice_best_nmatches(this->qend_alts),this->nmatches));
    this->ref_nmatches += Altsplice_best_nmatches(this->qend_alts); /* Not nmatches_to_trims, which is 0 for alts_substring */
    this->junction_splice_prob += this->qend_alts->medial_prob + Altsplice_best_distal_prob(this->qend_alts);
    this->total_splice_prob += this->qend_alts->medial_prob + Altsplice_best_distal_prob(this->qend_alts);
    this->found_score += Altsplice_found_score(this->qend_alts);
    this->score_within_trims += Altsplice_found_score(this->qend_alts);

  } else if (this->fusion_queryend_junction != NULL) {
    /* Rest of fusions handled below */

  } else {
    this->found_score += this->querylength - qend;
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
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp != this->plusp) {
      this->found_score += qstart;
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
    } else if (this->fusion_queryend_junction != NULL && this->fusion_plusp == this->plusp) {
      this->found_score += this->querylength - qend;
    }
  }


  /* this->nmatches += total_ninserts; */

  if (this->found_score < *found_score) {
    *found_score = this->found_score;
  }

  debug7(printf("Path_eval_nmatches returning %d matches in %s for score of %d within trims\n",
		this->nmatches,Intlist_to_string(this->endpoints),this->score_within_trims));

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
    } else if ((altsplice = this->qstart_alts) != NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
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
    } else if ((altsplice = this->qstart_alts )!= NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
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
  } else if ((altsplice = this->qend_alts) != NULL && (altsplice->medial_prob > 0.9 || Altsplice_best_distal_prob(altsplice) > 0.9)) {
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
    Path_eval_nmatches(&ignore_score,this,query_compress_fwd,query_compress_rev);
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

  if (a->method > b->method) {
    return -1;
  } else if (b->method > a->method) {
    return +1;
  } else {
    return 0;
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

  if (a->transcripts != NULL && b->transcripts == NULL) {
    return -1;
  } else if (b->transcripts != NULL && a->transcripts == NULL) {
    return +1;
  }

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
  } else {
    return 0;
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


T *
Path_eval_and_sort (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T *patharray, int npaths, 
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    char *queryuc_ptr, char *queryrc, char *quality_string,
		    int nmismatches_filter, int mincoverage_filter,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Hitlistpool_T hitlistpool, bool filterp) {
  int pathi;
  float maxlik, loglik;
  float total, qual;		/* For Bayesian mapq calculation */

  int randomi, i, j, k, l;
  T temp, path;
  double best_avg_junction_sense_prob, best_avg_junction_antisense_prob,
    best_end_sense_prob, best_end_antisense_prob, prob;
  /* int best_querystart, best_queryend, querystart, queryend; */
  int chop_qstart, chop_qend, qstart, qend;


  debug8(printf("Entered Path_eval_and_sort with patharray %p and npaths %d\n",patharray,npaths));
  if (npaths == 0) {
    /* Skip */
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;

  } else {
    /* Have already called Path_extend */
    
    /* 0.  Unalias circular alignments and trim chrbounds */
    for (i = 0; i < npaths; i++) {
      path = patharray[i];
      if (circularp[path->chrnum] == true) {
	Path_trim_circular_unalias(path);
      }
      Path_trim_chrbounds(path,query_compress_fwd,query_compress_rev,
			  intlistpool,univcoordlistpool,listpool,pathpool,
			  transcriptpool);
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
	for (l = i; l < j; l++) {
	  path = patharray[l];
	  if (path->sensedir == SENSE_FORWARD) {
	    if ((prob = Path_avg_junction_splice_prob(path)) > best_avg_junction_sense_prob) {
	      best_avg_junction_sense_prob = prob;
	    }
	    if ((prob = Path_end_splice_prob(path)) > best_end_sense_prob) {
	      best_end_sense_prob = prob;
	    }

	  } else {
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
	debug8(printf("Found an overlapping group of %d => re-sorting by Path_local_cmp\n",j - i));

	for (l = i; l < j; l++) {
	  patharray[l]->chop_qstart = chop_qstart;
	  patharray[l]->chop_qend = chop_qend;
	}

	/* Keep the best ones in the overlapping group */
	qsort(&(patharray[i]),j - i,sizeof(T),Path_local_cmp);
	debug8(printf("(3 best) Keeping by local_cmp ")); debug8(Path_print(patharray[i]));
	patharray[k++] = patharray[i];
      
	for (l = i + 1; l < j; l++) {
	  if (Path_structure_ignore_sense_cmp(&(patharray[l]),&(patharray[i])) == 0) {
	    debug8(printf("(3 identical) Eliminating by identical structure ")); debug8(Path_print(patharray[l]));
	    path = patharray[l];
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

	  } else if (Path_local_cmp(&(patharray[l]),&(patharray[i])) == 0) {
	    debug8(printf("(3 tie) Keeping by local_cmp ")); debug8(Path_print(patharray[l]));
	    patharray[k++] = patharray[l];

	  } else {
	    debug8(printf("(3 worse) Eliminating by local_cmp ")); debug8(Path_print(patharray[l]));
	    path = patharray[l];
	    Path_free(&path,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}

	i = j;
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

      qsort(patharray,npaths,sizeof(T),Path_global_cmp);
    }
    debug8(printf("Found the global best solution with %d nmatches\n",patharray[0]->nmatches));

    /* 5.  Check if we should be filtering the result */
    path = patharray[0];
    if (filterp == true && path->score_within_trims > nmismatches_filter) {
      debug8(printf("(4 filter) Best solution has too many %d nmismatches, so eliminating all\n",
		    path->score_within_trims));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else if (filterp == true && Path_coverage(path) < mincoverage_filter) {
      debug8(printf("(4 filter) Best solution has too little %d coverage, so eliminating all\n",
		    Path_coverage(path)));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else if (filterp == true && allow_soft_clips_p == false &&
	       (Intlist_head(path->endpoints) != 0 || Intlist_last_value(path->endpoints) != path->querylength)) {
      debug8(printf("(4 filter) Best solution has soft clips, so eliminating all\n",
		    Path_coverage(path)));
      for (k = 0; k < npaths; k++) {
	debug8(printf("(4 filter) Eliminating ")); debug8(Path_print(patharray[k]));
	path = patharray[k];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray);
      *npaths_primary = *npaths_altloc = 0;
      return (T *) NULL;

    } else {
      /* Otherwise, keep all paths equivalent to the best path */
      debug8(printf("(4 best) Keeping best ")); debug8(Path_print(patharray[0]));
      i = 1;			/* Skip the best path */

      for (k = 1; k < npaths; k++) {
	path = patharray[k];
	if (Path_global_cmp(&path,&(patharray[0])) == 0) {
	  debug8(printf("(4 tie) Keeping ")); debug8(Path_print(path));
	  patharray[i++] = path;
	} else {
	  debug8(printf("(4 worse) Eliminating ")); debug8(Path_print(path));
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      npaths = i;			
    }


    /* 4.  Check for main <=> fusion identity (can occur in single-end
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
	Path_trim_circular(path,query_compress_fwd,intlistpool,univcoordlistpool,listpool);

      } else {
	Path_trim_circular(path,query_compress_rev,intlistpool,univcoordlistpool,listpool);
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
      
      /* Handle the best paths in this overlapping group */
      if (nbest == 1) {
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
	  if (path->nmatches < best_nmatches) {
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
		 bool want_random_p_in, bool allow_soft_clips_p_in) {

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
  allow_soft_clips_p = allow_soft_clips_p_in;

  return;
}

