
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath-eval.h"
#include "debug_groups.h"

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
#if defined(DEBUG7) || defined(DEBUG_TRPATH)
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


  debug7(printf("\nEntering Trpath_eval_score (tplusp=%d, querylength=%d, ndiagonals=%d)\n",
		tplusp,querylength,Uintlist_length(trdiagonals)));

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
Trpath_score_read (Trpath_score_T *components, T this, int querylength, Compress_T query_compress_tr,
		   double match_reward, double mismatch_penalty,
		   double indel_open_penalty, double indel_extend_penalty,
		   double softclip_penalty) {
  int qstart, qend, ninserts;
  Trcoord_T trdiagonal;
  Intlist_T r, x;
  Uintlist_T q;
  Junction_T junction;
  Junctiontype_T junction_type;
  List_T j;
  int nmismatches, ref_nmismatches, segment_matches, nindels;


  debug7(printf("\nEntering Trpath_score_read on trpath %p\n",this));
  debug7(Trpath_print(this));

  components->matches = 0;
  components->mismatches = 0;
  components->indels = 0;
  components->indel_opens = 0;
  components->indel_extends = 0;
  components->softclip = 0;
  components->coverage = 0;
  components->transcript_discordant_p = false;
  components->found_score = 0;
  components->end_score = 0;
  components->score = 0.0;

  this->found_score = 0;

  assert(Uintlist_length(this->trdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);


  qstart = Intlist_head(this->endpoints);
  nmismatches = Intlist_head(this->nmismatches);
  ninserts = 0;
  
  j = this->junctions;		/* Put here before we handle querystart_alts */
  this->found_score += qstart;
  this->end_score = qstart;
  components->softclip += qstart;

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
    segment_matches = (qend - qstart) - nmismatches;
    components->matches += segment_matches;
    components->mismatches += nmismatches;
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
      junction_type = Junction_type(junction);
      if (junction_type == INS_JUNCTION || junction_type == DEL_JUNCTION) {
	nindels = Junction_nindels(junction);
	components->indels += nindels;
	components->indel_opens += 1;
	components->indel_extends += nindels - 1;
      }
    }
  }

  this->found_score += querylength - qend;
  components->softclip += querylength - qend;
  components->coverage = Intlist_last_value(this->endpoints) - Intlist_head(this->endpoints);
  components->found_score = this->found_score;

  components->score = match_reward * (double) components->matches
    - mismatch_penalty * (double) components->mismatches
    - indel_open_penalty * (double) components->indel_opens
    - indel_extend_penalty * (double) components->indel_extends
    - softclip_penalty * (double) components->softclip;

  if (querylength - qend > this->end_score) {
    this->end_score = querylength - qend;
  }
  components->end_score = this->end_score;

  debug7(printf("Trpath_score_read returning %s with score %.3f [matches=%d mismatches=%d indels=%d indel_opens=%d indel_extends=%d softclip=%d coverage=%d found_score=%d end_score=%d]\n",
		Intlist_to_string(this->endpoints),components->score,components->matches,
		components->mismatches,components->indels,components->indel_opens,
		components->indel_extends,components->softclip,components->coverage,
		components->found_score,components->end_score));

  return;
}


double
Trpath_score_pair (Trpath_score_T *pair_components, Trpath_score_T *components5,
		   Trpath_score_T *components3, T trpath5, T trpath3,
		   int querylength5, int querylength3,
		   Compress_T query5_compress_tr, Compress_T query3_compress_tr,
		   double match_reward, double mismatch_penalty,
		   double indel_open_penalty, double indel_extend_penalty,
		   double softclip_penalty, double transcript_discordance_penalty) {
  Trpath_score_T local5, local3;
  bool transcript_discordant_p;

  if (components5 == NULL) {
    components5 = &local5;
  }
  if (components3 == NULL) {
    components3 = &local3;
  }

  Trpath_score_read(components5,trpath5,querylength5,query5_compress_tr,
		    match_reward,mismatch_penalty,indel_open_penalty,
		    indel_extend_penalty,softclip_penalty);
  Trpath_score_read(components3,trpath3,querylength3,query3_compress_tr,
		    match_reward,mismatch_penalty,indel_open_penalty,
		    indel_extend_penalty,softclip_penalty);

  transcript_discordant_p = (trpath5->trnum == trpath3->trnum) ? false : true;

  if (pair_components != NULL) {
    pair_components->matches = components5->matches + components3->matches;
    pair_components->mismatches = components5->mismatches + components3->mismatches;
    pair_components->indels = components5->indels + components3->indels;
    pair_components->indel_opens = components5->indel_opens + components3->indel_opens;
    pair_components->indel_extends = components5->indel_extends + components3->indel_extends;
    pair_components->softclip = components5->softclip + components3->softclip;
    pair_components->coverage = components5->coverage + components3->coverage;
    pair_components->transcript_discordant_p = transcript_discordant_p;
    pair_components->found_score = components5->found_score + components3->found_score;
    pair_components->end_score = components5->end_score + components3->end_score;
    pair_components->score = components5->score + components3->score;
    if (transcript_discordant_p == true) {
      pair_components->score -= transcript_discordance_penalty;
    }
    return pair_components->score;
  } else {
    if (transcript_discordant_p == true) {
      return components5->score + components3->score - transcript_discordance_penalty;
    } else {
      return components5->score + components3->score;
    }
  }
}


void
Trpath_eval_nmatches (int *found_score, int *end_score, T this, int querylength, Compress_T query_compress_tr) {
  Trpath_score_T components;

  debug7(printf("\nEntering Trpath_eval_nmatches on trpath %p\n",this));
  debug7(Trpath_print(this));

  Trpath_score_read(&components,this,querylength,query_compress_tr,
		    /*match_reward*/1.0,/*mismatch_penalty*/1.0,
		    /*indel_open_penalty*/0.0,/*indel_extend_penalty*/0.0,
		    /*softclip_penalty*/1.0);

  if (this->found_score < *found_score) {
    *found_score = this->found_score;
  }
  if (this->end_score < *end_score) {
    *end_score = this->end_score;
  }

  debug7(printf("Trpath_eval_nmatches returning %s for found_score %d, end_score %d\n",
		Intlist_to_string(this->endpoints),this->found_score,this->end_score));

  return;
}


void
Trpath_eval_setup (Genomebits_T transcriptomebits_in) {

  transcriptomebits = transcriptomebits_in;

  return;
}
