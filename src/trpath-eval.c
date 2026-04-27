static char rcsid[] = "$Id: f97f3d20310ada1ff47dfa2de217a847f9d3e7ed $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath-eval.h"

#include <stdio.h>
#include <string.h>		/* For strcpy */
#include <math.h>		/* For rint */
#include <ctype.h>		/* For islower */

#include "assert.h"
#include "genomebits_count.h"
#include "junction.h"

static Genomebits_T transcriptomebits;

#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif


/* Trpath_eval_nmatches */
#ifdef DEBUG7
#define debug7(x) x
#else
#define debug7(x)
#endif

#define T Trpath_T


/* Simplified version of Trpath_eval_nmatches */
/* Returns score */
int
Trpath_eval_score (Uintlist_T trdiagonals, Intlist_T tr_endpoints, Intlist_T nmismatches_list, List_T junctions,
		   Compress_T query_compress_tr, int querylength, bool tplusp) {
  int score = 0;
  int qstart, qend, ninserts;
  Trcoord_T trdiagonal;
  Intlist_T r, x;
  Uintlist_T q;
  Junction_T junction;
  List_T j;
  /* bool insertionp = false; */
  /* int adj0; deletions - insertions */
  int nmismatches, ref_nmismatches;


  debug7(printf("\nEntering Trpath_eval_nmatches on path %p\n",this));
  debug7(Trpath_print(this));

  assert(Uintlist_length(trdiagonals) == Intlist_length(tr_endpoints) - 1);
  assert(Intlist_length(nmismatches_list) == Intlist_length(tr_endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(tr_endpoints) - 2);

  qstart = Intlist_head(tr_endpoints);
  nmismatches = Intlist_head(nmismatches_list);
  ninserts = 0;
  
  if (tplusp == true) {
    /* tplus */
    j = junctions;		/* Put here before we handle querystart_alts */
    score += qstart;

    /* Add qpos to get alignstart/alignend */
    for (q = trdiagonals, x = nmismatches_list, r = Intlist_next(tr_endpoints); q != NULL;
	 q = Uintlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);

      trdiagonal = Uintlist_head(q);
      /* left = univdiagonal - (Univcoord_T) querylength; */
      debug7(printf("Trpath_eval_nmatches: ninserts %d, qstart %d..qend %d at trdiagonal %u\n",
		    ninserts,qstart,qend,trdiagonal));

      if (nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",trdiagonal,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,
		      Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
							    trdiagonal,querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,/*genestrand*/0)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
								    trdiagonal,querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,/*genestrand*/0));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
							    trdiagonal,querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/true,/*genestrand*/0);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d mismatches from transcriptome over querypos %d..%d\n",
		      nmismatches,qstart,qend));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      /* this->nmatches += (qend - qstart) - nmismatches; */
      score += nmismatches;

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
      }
    }

  } else {
    /* minus */

    j = junctions;		/* Put here before we handle querystart_alts */
    score += qstart;

    /* Subtract qpos to get alignstart/alignend */
    for (q = trdiagonals, x = nmismatches_list, r = Intlist_next(tr_endpoints); q != NULL;
	 q = Uintlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
      qstart += ninserts;
      qend = Intlist_head(r);
      nmismatches = Intlist_head(x);

      trdiagonal = Uintlist_head(q);
      /* left = univdiagonal - (Univcoord_T) querylength; */
      debug7(printf("Trpath_eval_nmatches: ninserts %d, qstart %d..qend %d at trdiagonal %u\n",
		    ninserts,qstart,qend,trdiagonal));

      if (nmismatches >= 0) {
	debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",trdiagonal,qstart,qend));
	debug7(printf("%d mismatches expected vs %d measured\n",
		      nmismatches,Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
									trdiagonal,querylength,
									/*pos5*/qstart,/*pos3*/qend,/*plusp*/false,/*genestrand*/0)));
#ifdef CHECK_NMISMATCHES
	assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
								    trdiagonal,querylength,
								    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,/*genestrand*/0));
#endif
      } else {
	nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
							    trdiagonal,querylength,
							    /*pos5*/qstart,/*pos3*/qend,/*plusp*/false,/*genestrand*/0);
	Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
	debug7(printf("%d (%d ref) mismatches from transcriptome over querypos %d..%d\n",
		      nmismatches,ref_nmismatches,querylength - qend,querylength - qstart));
      }

      /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
      score += nmismatches;

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
      }
    }
  }

  score += querylength - qend;

  return score;
}


void
Trpath_eval_nmatches (int *found_score, T this, int querylength, Compress_T query_compress_tr) {
  int qstart, qend, ninserts;
  Trcoord_T trdiagonal;
  Intlist_T r, x;
  Uintlist_T q;
  Junction_T junction;
  List_T j;
  int nmismatches, ref_nmismatches;


  debug7(printf("\nEntering Trpath_eval_nmatches on trpath %p\n",this));
  debug7(Trpath_print(this));

  this->found_score = 0;

  assert(Uintlist_length(this->trdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);


  qstart = Intlist_head(this->endpoints);
  nmismatches = Intlist_head(this->nmismatches);
  ninserts = 0;
  
  j = this->junctions;		/* Put here before we handle querystart_alts */
  this->found_score += qstart;

  /* Add qpos to get alignstart/alignend */
  for (q = this->trdiagonals, x = this->nmismatches, r = Intlist_next(this->endpoints); q != NULL;
       q = Uintlist_next(q), x = Intlist_next(x), r = Intlist_next(r), j = List_next(j)) {
    qstart += ninserts;
    qend = Intlist_head(r);
    nmismatches = Intlist_head(x);

    trdiagonal = Uintlist_head(q);
    /* left = trdiagonal - (Univcoord_T) querylength; */
    debug7(printf("Trpath_eval_nmatches: ninserts %d, qstart %d..qend %d at trdiagonal %u [%u]\n",
		  ninserts,qstart,qend,trdiagonal,trdiagonal - this->troffset));

    if (nmismatches >= 0) {
      debug7(printf("Checking mismatches at %u from querystart %d to queryend %d\n",trdiagonal - this->troffset,qstart,qend));
      debug7(printf("%d mismatches expected vs %d measured\n",
		    nmismatches,
		    Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
							  trdiagonal,querylength,
							  /*pos5*/qstart,/*pos3*/qend,/*plusp*/this->tplusp,/*genestrand*/0)));
#ifdef CHECK_NMISMATCHES
      assert(nmismatches == Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
								  trdiagonal,querylength,
								  /*pos5*/qstart,/*pos3*/qend,/*plusp*/this->tplusp,/*genestrand*/0));
#endif
    } else {
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,/*alt*/NULL,query_compress_tr,
							  trdiagonal,querylength,
							  /*pos5*/qstart,/*pos3*/qend,/*plusp*/this->tplusp,/*genestrand*/0);
      Intlist_head_set(x,nmismatches);		/* Save for Stage3end_new_substrings */
      debug7(printf("%d mismatches from genome over querypos %d..%d\n",
		    nmismatches,qstart,qend));
    }
    
    /* Could potentially check here if qstart < qend, but relying upon caller to use endpoints_acceptable_p */
    /* this->nmatches += (qend - qstart) - nmismatches; */
    this->found_score += nmismatches;

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
    }
  }

  this->found_score += querylength - qend;

  if (this->found_score < *found_score) {
    *found_score = this->found_score;
  }

  debug7(printf("Trpath_eval_nmatches returning %s for found_score %d\n",
		Intlist_to_string(this->endpoints),this->found_score));

  return;
}


void
Trpath_eval_setup (Genomebits_T transcriptomebits_in) {

  transcriptomebits = transcriptomebits_in;

  return;
}

