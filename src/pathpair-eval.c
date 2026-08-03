#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "pathpair-eval.h"
#include "debug_groups.h"
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
#include "transcript.h"


static bool *circularp;
static bool *chrsubsetp;
static bool *altlocp;

static Outputtype_T output_type;
static bool splicingp;
static bool resolve_inner_p;
static bool want_random_p;

static bool allow_soft_clips_start_p;
static bool allow_soft_clips_end_p;
static int paired_terminal_anchor_min_score;
static int paired_terminal_anchor_slop;
static bool mincoverage_filter_per_read_p;
static bool pathpair_equal_coverage_diff_tiebreak_p;
static double score_match_reward;
static double score_mismatch_penalty;
static double score_indel_open_penalty;
static double score_indel_extend_penalty;
static double score_splice_unannotated_penalty;
static double score_pair_transcript_discordance_penalty;
static double score_pair_margin;
static double score_min_read;
static double score_min_pair;

static Chrpos_T max_insertlength;

static Transcriptome_T transcriptome;


/* Splice goodness now determined by Splice_accept_p */
/* #define OMIT_BAD_SPLICES 1 */

#define NMATCHES_FACTOR 8
#define INSERTLENGTH_FACTOR 1.5
#define OUTERLENGTH_FACTOR 1.5

#define INNER_NEEDED 15
#define POSTFALLBACK_TERMINAL_ANCHOR_MAX_LOW_SCORE 0
#define POSTFALLBACK_TERMINAL_ANCHOR_MAX_HIGH_SCORE 3
#define POSTFALLBACK_TERMINAL_ANCHOR_MAX_HIGH_INDELS 1
#define POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_LOW_SCORE 1
#define POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_LOW_INDELS 1
#define POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_HIGH_SCORE 1


/* Pathpair_preeval_and_uniq */
#if defined(DEBUG7) || defined(DEBUG_PATH_CORE)
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Pathpair_eval_and_sort */
#if defined(DEBUG8) || defined(DEBUG_PATH_CORE)
#define debug8(x) x
#else
#define debug8(x)
#endif

#define T Path_T

static bool
path_single_filter_pass_p (Path_T path, int coverage, int nmismatches_filter,
			   int mincoverage_filter, int querylength,
			   bool allow_soft_clips_qstart_p, bool allow_soft_clips_qend_p) {
  if (path->score_within_trims > nmismatches_filter) {
    return false;
  } else if (coverage < mincoverage_filter) {
    return false;
  } else if (allow_soft_clips_qstart_p == false &&
	     Intlist_head(path->endpoints) != 0) {
    return false;
  } else if (allow_soft_clips_qend_p == false &&
	     Intlist_last_value(path->endpoints) != querylength) {
    return false;
  } else {
    return true;
  }
}

static void
save_rejected_single_path (List_T *paths, Path_T path,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  Path_T copy;

  if (paths != NULL) {
    copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		     pathpool,vectorpool,transcriptpool,hitlistpool);
    *paths = Hitlist_push(*paths,hitlistpool,(void *) copy
			  hitlistpool_trace(__FILE__,__LINE__));
  }

  return;
}

typedef struct Pathpair_one_low_side_T {
  bool low5p;
  Path_T low_path;
  Path_T high_path;
  int low_coverage;
  int high_coverage;
  int low_mincoverage_filter;
  int high_mincoverage_filter;
  int low_querylength;
  int high_querylength;
} Pathpair_one_low_side_T;


static int
Pathpair_global_cmp (const void *x, const void *y);


static void
Path_final_indel_counts (int *indel_opens, int *indel_extends, Path_T path) {
  List_T p;
  Junction_T junction;
  int nindels;

  *indel_opens = 0;
  *indel_extends = 0;

  for (p = path->junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == INS_JUNCTION ||
	Junction_type(junction) == DEL_JUNCTION) {
      nindels = Junction_nindels(junction);
      if (nindels > 0) {
	*indel_opens += 1;
	*indel_extends += nindels - 1;
      }
    }
  }

  return;
}


static int
Path_unannotated_splice_count (Path_T path) {
  List_T p;
  Junction_T junction;
  int count = 0;

  for (p = path->junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION &&
	Junction_annotatedp(junction) == false) {
      count++;
    }
  }

  return count;
}


static bool
Path_transcript_evidence_p (Path_T path) {
  return path->transcripts != NULL || path->invalid_transcripts != NULL;
}


double
Path_final_score (Path_T path) {
  int indel_opens, indel_extends, unannotated_splices;

  Path_final_indel_counts(&indel_opens,&indel_extends,path);
  unannotated_splices = Path_unannotated_splice_count(path);

  return score_match_reward * (double) path->nmatches
    - score_mismatch_penalty * (double) path->score_within_trims
    - score_indel_open_penalty * (double) indel_opens
    - score_indel_extend_penalty * (double) indel_extends
    - score_splice_unannotated_penalty * (double) unannotated_splices;
}


bool
Path_final_score_reportable_p (Path_T path) {
  return Path_final_score(path) >= score_min_read ? true : false;
}


void
Path_set_final_score (Path_T path, int pair_score) {
  path->final_score_p = true;
  path->final_read_score = (int) rint(Path_final_score(path));
  path->final_pair_score = pair_score;
  return;
}


static double
Pathpair_transcript_discordance_penalty (Pathpair_T pathpair) {
  if (Path_transcript_evidence_p(pathpair->path5) == true &&
      Path_transcript_evidence_p(pathpair->path3) == true &&
      Transcript_intersectp(pathpair->path5->transcripts,pathpair->path3->transcripts) == false) {
    return score_pair_transcript_discordance_penalty;
  } else {
    return 0.0;
  }
}


static double
Pathpair_final_score (Pathpair_T pathpair) {
  return Path_final_score(pathpair->path5) + Path_final_score(pathpair->path3)
    - Pathpair_transcript_discordance_penalty(pathpair);
}


static bool
Pathpair_final_score_reportable_p (Pathpair_T pathpair) {
  double read_score_5 = Path_final_score(pathpair->path5);
  double read_score_3 = Path_final_score(pathpair->path3);
  double pair_score = read_score_5 + read_score_3
    - Pathpair_transcript_discordance_penalty(pathpair);

  if (read_score_5 < score_min_read || read_score_3 < score_min_read) {
    return false;
  } else if (pair_score < score_min_pair) {
    return false;
  } else {
    return true;
  }
}


static void
Pathpair_set_final_scores (Pathpair_T pathpair) {
  int pair_score = (int) rint(Pathpair_final_score(pathpair));

  Path_set_final_score(pathpair->path5,pair_score);
  Path_set_final_score(pathpair->path3,pair_score);

  return;
}


static int
Pathpair_final_score_cmp (const void *x, const void *y) {
  Pathpair_T a = * (Pathpair_T *) x;
  Pathpair_T b = * (Pathpair_T *) y;
  double score_a = Pathpair_final_score(a);
  double score_b = Pathpair_final_score(b);

  if (score_a > score_b) {
    return -1;
  } else if (score_b > score_a) {
    return +1;
  } else {
    return Pathpair_global_cmp(x,y);
  }
}

static bool
pathpair_one_low_side_setup (Pathpair_one_low_side_T *side, Pathpair_T pathpair,
			     int coverage5, int coverage3,
			     int mincoverage_filter_5, int mincoverage_filter_3,
			     int querylength5, int querylength3) {
  bool low5p = coverage5 < mincoverage_filter_5;
  bool low3p = coverage3 < mincoverage_filter_3;

  if (low5p == low3p) {
    return false;
  } else if (low5p == true) {
    side->low5p = true;
    side->low_path = pathpair->path5;
    side->high_path = pathpair->path3;
    side->low_coverage = coverage5;
    side->high_coverage = coverage3;
    side->low_mincoverage_filter = mincoverage_filter_5;
    side->high_mincoverage_filter = mincoverage_filter_3;
    side->low_querylength = querylength5;
    side->high_querylength = querylength3;
  } else {
    side->low5p = false;
    side->low_path = pathpair->path3;
    side->high_path = pathpair->path5;
    side->low_coverage = coverage3;
    side->high_coverage = coverage5;
    side->low_mincoverage_filter = mincoverage_filter_3;
    side->high_mincoverage_filter = mincoverage_filter_5;
    side->low_querylength = querylength3;
    side->high_querylength = querylength5;
  }

  return true;
}

static bool
pathpair_mincoverage_rescue_side_p (const Pathpair_one_low_side_T *side) {
  if (side->low_path->score_within_trims != 0 ||
      side->high_path->score_within_trims != 0) {
    return false;
  } else if (Path_nindels(side->low_path) + Path_nindels(side->high_path) != 0) {
    return false;
  } else if (side->low_coverage + side->high_coverage <
	     side->low_mincoverage_filter + side->high_mincoverage_filter) {
    return false;
  } else if (side->low_coverage < (side->low_mincoverage_filter * 3)/4) {
    return false;
  } else {
    return true;
  }
}

static bool
path_terminal_anchor_p (Path_T path, int coverage, int querylength) {
  int qstart, qend;

  if (coverage < paired_terminal_anchor_min_score) {
    return false;
  } else {
    qstart = Intlist_head(path->endpoints);
    qend = Intlist_last_value(path->endpoints);
    return (qstart <= paired_terminal_anchor_slop ||
	    qend >= querylength - paired_terminal_anchor_slop) ? true : false;
  }
}

static bool
path_clean_full_partner_p (Path_T path, int coverage, int mincoverage_filter, int querylength) {
  if (path->score_within_trims != 0) {
    return false;
  } else if (Path_nindels(path) != 0) {
    return false;
  } else if (coverage < mincoverage_filter) {
    return false;
  } else if (Intlist_head(path->endpoints) != 0 ||
	     Intlist_last_value(path->endpoints) != querylength) {
    return false;
  } else {
    return true;
  }
}

bool
Pathpair_near_exact_terminal_anchor_partner_p (Path_T path) {
  return (path->score_within_trims <= POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_HIGH_SCORE &&
	  Path_nindels(path) == 0) ? true : false;
}

static bool
path_near_exact_terminal_anchor_core_p (Path_T path, int coverage, int mincoverage_filter,
					int querylength, bool clean_full_partner_p) {
  int mincoverage, nindels;
  bool one_error_no_indel_p, transcript_backed_one_error_p, paired_supported_one_error_p;

  if (path->score_within_trims < 0 ||
      path->score_within_trims > POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_LOW_SCORE) {
    return false;
  } else if ((nindels = Path_nindels(path)) > POSTFALLBACK_NEAR_EXACT_TERMINAL_ANCHOR_MAX_LOW_INDELS) {
    return false;
  } else {
    one_error_no_indel_p =
      (path->score_within_trims > POSTFALLBACK_TERMINAL_ANCHOR_MAX_LOW_SCORE && nindels == 0) ?
      true : false;
    transcript_backed_one_error_p =
      (one_error_no_indel_p == true && path->transcriptome_method_p == true &&
       path->transcripts != NULL) ? true : false;
    paired_supported_one_error_p =
      (one_error_no_indel_p == true && clean_full_partner_p == true) ? true : false;

    mincoverage = paired_terminal_anchor_min_score;
    /* A one-error terminal anchor needs more support than an exact one. */
    if ((path->score_within_trims > POSTFALLBACK_TERMINAL_ANCHOR_MAX_LOW_SCORE || nindels > 0) &&
	transcript_backed_one_error_p == false &&
	paired_supported_one_error_p == false &&
	mincoverage_filter/2 > mincoverage) {
      mincoverage = mincoverage_filter/2;
    }
    /* Paired/deferred callers may keep a one-mismatch genomic terminal anchor
       with one-third coverage if the opposite mate is a clean full-length path. */
    if (paired_supported_one_error_p == true &&
	transcript_backed_one_error_p == false &&
	mincoverage_filter/3 > mincoverage) {
      mincoverage = mincoverage_filter/3;
    }

    if (coverage < mincoverage) {
      return false;
    } else {
      return path_terminal_anchor_p(path,coverage,querylength);
    }
  }
}

static bool
path_near_exact_terminal_anchor_with_partner_p (Path_T path, int coverage, int mincoverage_filter,
						int querylength, Path_T partner, int partner_coverage,
						int partner_mincoverage_filter, int partner_querylength) {
  bool clean_full_partner_p =
    path_clean_full_partner_p(partner,partner_coverage,partner_mincoverage_filter,partner_querylength);

  return path_near_exact_terminal_anchor_core_p(path,coverage,mincoverage_filter,
						querylength,clean_full_partner_p);
}

static bool
pathpair_terminal_anchor_rescue_side_p (const Pathpair_one_low_side_T *side) {
  if (side->low_path->score_within_trims != 0 ||
      side->high_path->score_within_trims != 0) {
    return false;
  } else if (Path_nindels(side->low_path) + Path_nindels(side->high_path) != 0) {
    return false;
  } else {
    return (side->high_coverage >= side->high_mincoverage_filter &&
	    path_terminal_anchor_p(side->low_path,side->low_coverage,side->low_querylength)) ?
      true : false;
  }
}

static bool
pathpair_terminal_anchor_rescue_p (Pathpair_T pathpair, int coverage5, int coverage3,
				   int mincoverage_filter_5, int mincoverage_filter_3,
				   int querylength5, int querylength3) {
  Pathpair_one_low_side_T side;

  return (pathpair_one_low_side_setup(&side,pathpair,coverage5,coverage3,
				      mincoverage_filter_5,mincoverage_filter_3,
				      querylength5,querylength3) == true &&
	  pathpair_terminal_anchor_rescue_side_p(&side) == true) ? true : false;
}

static bool
pathpair_low_coverage_side_nonexact_p (Pathpair_T pathpair, int coverage5, int coverage3,
				       int mincoverage_filter_5, int mincoverage_filter_3) {
  if (coverage5 < mincoverage_filter_5 &&
      (pathpair->path5->score_within_trims > 0 || Path_nindels(pathpair->path5) > 0)) {
    return true;
  } else if (coverage3 < mincoverage_filter_3 &&
	     (pathpair->path3->score_within_trims > 0 || Path_nindels(pathpair->path3) > 0)) {
    return true;
  } else {
    return false;
  }
}

static bool
pathpair_low_coverage_side_substantial_p (int coverage5, int coverage3,
					  int mincoverage_filter_5, int mincoverage_filter_3,
					  int querylength5, int querylength3) {
  if (coverage5 < mincoverage_filter_5) {
    return (coverage5 * 3 >= querylength5) ? true : false;
  } else if (coverage3 < mincoverage_filter_3) {
    return (coverage3 * 3 >= querylength3) ? true : false;
  } else {
    return false;
  }
}

static bool
pathpair_postfallback_terminal_anchor_rescue_side_p (const Pathpair_one_low_side_T *side) {
  int nindels_high;

  nindels_high = Path_nindels(side->high_path);
  return (side->low_path->score_within_trims <= POSTFALLBACK_TERMINAL_ANCHOR_MAX_LOW_SCORE &&
	  side->high_path->score_within_trims <= POSTFALLBACK_TERMINAL_ANCHOR_MAX_HIGH_SCORE &&
	  Path_nindels(side->low_path) == 0 &&
	  nindels_high <= POSTFALLBACK_TERMINAL_ANCHOR_MAX_HIGH_INDELS &&
	  (nindels_high == 0 || side->low_coverage >= side->low_mincoverage_filter/2) &&
	  side->high_coverage >= side->high_mincoverage_filter &&
	  path_terminal_anchor_p(side->low_path,side->low_coverage,side->low_querylength)) ?
    true : false;
}

static bool
pathpair_postfallback_terminal_anchor_rescue_p (Pathpair_T pathpair, int coverage5, int coverage3,
						int mincoverage_filter_5, int mincoverage_filter_3,
						int querylength5, int querylength3) {
  Pathpair_one_low_side_T side;

  return (pathpair_one_low_side_setup(&side,pathpair,coverage5,coverage3,
				      mincoverage_filter_5,mincoverage_filter_3,
				      querylength5,querylength3) == true &&
	  pathpair_postfallback_terminal_anchor_rescue_side_p(&side) == true) ? true : false;
}

typedef enum Pathpair_local_anchor_rescue_T {
  PATHPAIR_LOCAL_ANCHOR_RESCUE_NONE,
  PATHPAIR_LOCAL_ANCHOR_RESCUE_TERMINAL,
  PATHPAIR_LOCAL_ANCHOR_RESCUE_POSTFALLBACK
} Pathpair_local_anchor_rescue_T;

static Pathpair_local_anchor_rescue_T
pathpair_local_terminal_anchor_rescue (int *coverage5_out, int *coverage3_out,
				       Pathpair_T pathpair,
				       bool best_low_coverage_p,
				       bool best_low_coverage_side_nonexact_p,
				       bool best_terminal_anchor_rescue_p,
				       bool best_postfallback_terminal_anchor_rescue_p,
				       bool require_best_low_coverage_side_nonexact_p,
				       bool require_low_coverage_side_substantial_p,
				       int mincoverage_filter_5, int mincoverage_filter_3,
				       int querylength5, int querylength3,
				       List_T *terminal_anchor_pathpairs_deferred,
				       bool allow_postfallback_terminal_anchor_rescue_p) {
  int coverage5, coverage3;

  if (best_low_coverage_p == false ||
      best_terminal_anchor_rescue_p == true) {
    return PATHPAIR_LOCAL_ANCHOR_RESCUE_NONE;
  } else if (require_best_low_coverage_side_nonexact_p == true &&
	     best_low_coverage_side_nonexact_p == false) {
    return PATHPAIR_LOCAL_ANCHOR_RESCUE_NONE;
  } else {
    coverage5 = Path_coverage(pathpair->path5);
    coverage3 = Path_coverage(pathpair->path3);
  }

  if (pathpair_terminal_anchor_rescue_p(pathpair,coverage5,coverage3,
					mincoverage_filter_5,mincoverage_filter_3,
					querylength5,querylength3) == true &&
      (require_low_coverage_side_substantial_p == false ||
       pathpair_low_coverage_side_substantial_p(coverage5,coverage3,
						mincoverage_filter_5,mincoverage_filter_3,
						querylength5,querylength3) == true)) {
    *coverage5_out = coverage5;
    *coverage3_out = coverage3;
    return PATHPAIR_LOCAL_ANCHOR_RESCUE_TERMINAL;

  } else if (best_postfallback_terminal_anchor_rescue_p == false &&
	     (terminal_anchor_pathpairs_deferred != NULL ||
	      allow_postfallback_terminal_anchor_rescue_p == true) &&
	     pathpair_postfallback_terminal_anchor_rescue_p(pathpair,coverage5,coverage3,
							   mincoverage_filter_5,mincoverage_filter_3,
							   querylength5,querylength3) == true &&
	     (require_low_coverage_side_substantial_p == false ||
	      pathpair_low_coverage_side_substantial_p(coverage5,coverage3,
						       mincoverage_filter_5,mincoverage_filter_3,
						       querylength5,querylength3) == true)) {
    *coverage5_out = coverage5;
    *coverage3_out = coverage3;
    return PATHPAIR_LOCAL_ANCHOR_RESCUE_POSTFALLBACK;

  } else {
    return PATHPAIR_LOCAL_ANCHOR_RESCUE_NONE;
  }
}

static bool
pathpair_near_exact_terminal_anchor_rescue_side_p (const Pathpair_one_low_side_T *side) {
  return (Pathpair_near_exact_terminal_anchor_partner_p(side->high_path) == true &&
	  side->high_coverage >= side->high_mincoverage_filter &&
	  path_near_exact_terminal_anchor_with_partner_p(side->low_path,side->low_coverage,
							 side->low_mincoverage_filter,side->low_querylength,
							 side->high_path,side->high_coverage,
							 side->high_mincoverage_filter,
							 side->high_querylength)) ? true : false;
}

static bool
pathpair_near_exact_low_side_path_p (const Pathpair_one_low_side_T *side) {
  return (side->high_coverage >= side->high_mincoverage_filter &&
	  path_near_exact_terminal_anchor_with_partner_p(side->low_path,side->low_coverage,
							 side->low_mincoverage_filter,side->low_querylength,
							 side->high_path,side->high_coverage,
							 side->high_mincoverage_filter,
							 side->high_querylength)) ? true : false;
}

static bool
pathpair_clean_transcript_overlap_p (Pathpair_T pathpair, int coverage5, int coverage3,
				     int mincoverage_filter_5, int mincoverage_filter_3) {
  if (pathpair->outerlength > INNER_NEEDED) {
    return false;
  } else if (pathpair->path5->transcripts == NULL || pathpair->path3->transcripts == NULL) {
    return false;
  } else if (pathpair->transcript_concordant_p == false) {
    return false;
  } else if (pathpair->path5->score_within_trims != 0 ||
	     pathpair->path3->score_within_trims != 0) {
    return false;
  } else if (Path_nindels(pathpair->path5) + Path_nindels(pathpair->path3) != 0) {
    return false;
  } else if (coverage5 < mincoverage_filter_5 || coverage3 < mincoverage_filter_3) {
    return false;
  } else {
    return true;
  }
}

static int
Pathpair_equal_coverage_diff_cmp (Pathpair_T a, Pathpair_T b) {
  int ndiffs_a, ndiffs_b, nindels_a, nindels_b;

  ndiffs_a = Path_ndiffs(a->pathL) + Path_ndiffs(a->pathH);
  ndiffs_b = Path_ndiffs(b->pathL) + Path_ndiffs(b->pathH);

  if (ndiffs_a < ndiffs_b) {
    return -1;
  } else if (ndiffs_b < ndiffs_a) {
    return +1;
  }

  nindels_a = Path_nindels(a->pathL) + Path_nindels(a->pathH);
  nindels_b = Path_nindels(b->pathL) + Path_nindels(b->pathH);

  if (nindels_a < nindels_b) {
    return -1;
  } else if (nindels_b < nindels_a) {
    return +1;
  } else {
    return 0;
  }
}

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
  int diff_cmp;
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

#if defined(DEBUG8) || defined(DEBUG_PATH_CORE)
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

  if (pathpair_equal_coverage_diff_tiebreak_p == true && coverage_a == coverage_b &&
      (diff_cmp = Pathpair_equal_coverage_diff_cmp(a,b)) != 0) {
    return diff_cmp;
  }
    
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
  int diff_cmp;
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

  if (pathpair_equal_coverage_diff_tiebreak_p == true && coverage_a == coverage_b &&
      (diff_cmp = Pathpair_equal_coverage_diff_cmp(a,b)) != 0) {
    return diff_cmp;
  }
    
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


#if 0
/* Accumulates pathpairs into the given list called pathpairs */
List_T
Pathpair_preeval_and_uniq (List_T pathpairs, Pathpair_T *pathpairarray, int npaths,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Hitlistpool_T hitlistpool) {
  int best_nmatches, nmatches;

  int i, j, k, l, ll;
  Pathpair_T pathpair;

  bool *duplicatep;

  int chop_qstart_L, chop_qend_L, chop_qstart_H, chop_qend_H, qstart, qend;

  
#if defined(DEBUG7) || defined(DEBUG_PATH_CORE)
  printf("Entered Pathpair_preeval_and_uniq with %d pathpairs\n",npaths);
  for (int pathi = 0; pathi < npaths; pathi++) {
    printf("Pathpair (sense %d and %d): %p and %p\n",
	   pathpairarray[pathi]->pathL->sensedir,pathpairarray[pathi]->pathH->sensedir,
	   pathpairarray[pathi]->pathL,pathpairarray[pathi]->pathH);
  }
  printf("\n");
#endif

  /* 0.  Sort by structure to remove duplicates */
  if (npaths == 0) {
    return pathpairs;

  } else {
    duplicatep = (bool *) CALLOC(npaths,sizeof(bool));
    qsort(pathpairarray,npaths,sizeof(Pathpair_T),Pathpair_interval_cmp);

    k = 0;
    i = 0;
    while (i < npaths) {
      j = i + 1;
      while (j < npaths && Pathpair_overlap_p(pathpairarray[j],pathpairarray[i]) == true) {
	j++;
      }
      
      if (j - i == 1) {
	pathpairarray[k++] = pathpairarray[i];

      } else {
	debug7(printf("Found an overlapping group of %d => eliminating duplicates\n",j - i));
	for (l = i; l < j; l++) {
	  if (duplicatep[l] == true) {
	    /* Already found to be a duplicate */
	  } else {
	    for (ll = l + 1; ll < j; ll++) {
	      /* Need to keep sensedirs separate so we can extend them separately */
	      if (Pathpair_structure_cmp(&(pathpairarray[l]),&(pathpairarray[ll])) == 0) {
		debug7(printf("Path %d is a duplicate of %d\n",ll,l));
		duplicatep[ll] = true;
	      }
	    }
	  }
	}

	for (l = i; l < j; l++) {
	  pathpair = pathpairarray[l];
	  if (duplicatep[l] == false) {
	    pathpairarray[k++] = pathpair;
	  } else {
	    Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}
      }
      
      i = j;
    }

    FREE(duplicatep);
    npaths = k;
  }

  debug7(printf("(1) After removing duplicates (except for sensedir), have %d paths\n",npaths));


  /* Find best chopped solution globally (among given sensedirs) */
  /* This code is designed to help speed up the algorithm, but it might be too greedy */
  if (npaths == 1) {
    pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpairarray[0]
			     hitlistpool_trace(__FILE__,__LINE__));
    return pathpairs;

  } else {
    /* Find chopped regions */
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

    best_nmatches = 0;
    for (i = 0; i < npaths; i++) {
      /* Set chopped regions and evaluate for nmatches */
      pathpair = pathpairarray[i];
      pathpair->pathL->chop_qstart = chop_qstart_L;
      pathpair->pathL->chop_qend = chop_qend_L;
      pathpair->pathH->chop_qstart = chop_qstart_H;
      pathpair->pathH->chop_qend = chop_qend_H;

      if ((nmatches = Path_eval_nmatches_chopped(pathpair->pathL) + 
	   Path_eval_nmatches_chopped(pathpair->pathH)) > best_nmatches) {
	best_nmatches = nmatches;
      }
    }

    debug7(printf("best_nmatches is %d\n",best_nmatches));

    if (best_nmatches == 0) {
      /* Unable to find global best solutions */
      for (i = 0; i < npaths; i++) {
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpairarray[i]
				 hitlistpool_trace(__FILE__,__LINE__));
      }

    } else {
      for (i = 0; i < npaths; i++) {
	/* Set chopped regions and evaluate for nmatches */
	pathpair = pathpairarray[i];
	if (Path_eval_nmatches_chopped(pathpair->pathL) + 
	    Path_eval_nmatches_chopped(pathpair->pathH) == best_nmatches) {
	  debug7(printf("Keeping by nmatches ")); debug7(Pathpair_print(pathpair));
	  pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				   hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  debug7(printf("Eliminating by nmatches ")); debug7(Pathpair_print(pathpair));
	  Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  return pathpairs;
}
#endif



/* Accumulates pathpairs into the given list called pathpairs */
List_T
Pathpair_uniq (List_T pathpairs, Pathpair_T *pathpairarray, int npaths,
	       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool) {
  int best_nmatches, nmatches;

  int i, j, k, l, ll;
  Pathpair_T pathpair;

  bool *duplicatep;

  int chop_qstart_L, chop_qend_L, chop_qstart_H, chop_qend_H, qstart, qend;

  
#if defined(DEBUG7) || defined(DEBUG_PATH_CORE)
  printf("Entered Pathpair_preeval_and_uniq with %d pathpairs\n",npaths);
  for (int pathi = 0; pathi < npaths; pathi++) {
    printf("Pathpair (sense %d and %d): %p and %p\n",
	   pathpairarray[pathi]->pathL->sensedir,pathpairarray[pathi]->pathH->sensedir,
	   pathpairarray[pathi]->pathL,pathpairarray[pathi]->pathH);
    Pathpair_print(pathpairarray[pathi]);
  }
  printf("\n");
#endif

  /* 0.  Sort by structure to remove duplicates */
  if (npaths == 0) {
    return pathpairs;

  } else {
    duplicatep = (bool *) CALLOC(npaths,sizeof(bool));
    qsort(pathpairarray,npaths,sizeof(Pathpair_T),Pathpair_interval_cmp);

    k = 0;
    i = 0;
    while (i < npaths) {
      j = i + 1;
      while (j < npaths && Pathpair_overlap_p(pathpairarray[j],pathpairarray[i]) == true) {
	j++;
      }
      
      if (j - i == 1) {
	pathpairarray[k++] = pathpairarray[i];

      } else {
	debug7(printf("Found an overlapping group of %d => eliminating duplicates\n",j - i));
	for (l = i; l < j; l++) {
	  if (duplicatep[l] == true) {
	    /* Already found to be a duplicate */
	  } else {
	    for (ll = l + 1; ll < j; ll++) {
	      /* Need to keep sensedirs separate so we can extend them separately */
	      if (Pathpair_structure_cmp(&(pathpairarray[l]),&(pathpairarray[ll])) == 0) {
		debug7(printf("Path %d is a duplicate of %d\n",ll,l));
		duplicatep[ll] = true;
	      }
	    }
	  }
	}

	for (l = i; l < j; l++) {
	  pathpair = pathpairarray[l];
	  if (duplicatep[l] == false) {
	    pathpairarray[k++] = pathpair;
	  } else {
	    Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}
      }
      
      i = j;
    }

    FREE(duplicatep);
    npaths = k;
  }

  debug7(printf("(1) After removing duplicates (except for sensedir), have %d paths\n",npaths));


  for (i = 0; i < npaths; i++) {
    /* Set chopped regions and evaluate for nmatches */
    pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpairarray[i]
			     hitlistpool_trace(__FILE__,__LINE__));
  }

  return pathpairs;
}




/* Taken from step 4 in Pathpair_eval_and_sort, which finds the best
   sensedir for paths in each local region */

#if 0
/* Can be led astray by end splice probs */
List_T
Pathpair_pick_dir (List_T pathpairs, Pathpair_T *pathpairarray, int npaths,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		   Hitlistpool_T hitlistpool) {

  Pathpair_T pathpair;
  int i, j, k, l;
  double best_avg_junction_sense_prob, best_avg_junction_antisense_prob,
    best_end_sense_prob, best_end_antisense_prob, prob;


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

  for (i = 0; i < npaths; i++) {
    pathpair = pathpairarray[i];
    pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
			     hitlistpool_trace(__FILE__,__LINE__));
  }

  return pathpairs;
}
#endif


/* filtering is creating an issue between pathpairs != NULL, stopping
   path creation, and then an empty pathpairarray */

Pathpair_T *
Pathpair_eval_and_sort (int *found_score_5, int *found_score_3, int *end_score_5, int *end_score_3,
			int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
			List_T *singlepaths5_from_rejected, List_T *singlepaths3_from_rejected,
			Pathpair_T *pathpairarray, int npaths, Stage1_T stage1_5, Stage1_T stage1_3,
			Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			Shortread_T queryseq5, Shortread_T queryseq3,
			char *queryuc_ptr_5, char *queryrc5, char *queryuc_ptr_3, char *queryrc3, 
			char *quality_string_5, char *quality_string_3,

			int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			int max_nblocks, int localdb_nmismatches_allowed,

			Knownsplicing_T knownsplicing, Knownindels_T knownindels,

			int nmismatches_filter_5, int nmismatches_filter_3,
			int mincoverage_filter_5, int mincoverage_filter_3,
			int querylength5, int querylength3,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
			Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			List_T *terminal_anchor_pathpairs_deferred,
			List_T *terminal_anchor_paths5_deferred,
			List_T *terminal_anchor_paths3_deferred,
			bool filterp, bool allow_postfallback_terminal_anchor_rescue_p,
			bool allow_near_exact_terminal_anchor_rescue_p) {
  int pathi;
  float maxlik, loglik;

  float total, qual;
  int mapq_score;

  int randomi, i, j, k, l, structure_cmp;
  int rescue_npaths;
  Pathpair_T temp, pathpair, *rescue_pathpairarray;
  Path_T path5, path3;
  Univcoord_T high, low, distal5, distal3;
  int tr_insertlength;

  bool valid5p, valid3p;
  bool allow_soft_clips_qstart_p, allow_soft_clips_qend_p;
  bool possible_fusion5_p = false, possible_fusion3_p = false;
  bool best_low_coverage_p, best_low_coverage_side_nonexact_p;
  bool best_terminal_anchor_rescue_p, best_postfallback_terminal_anchor_rescue_p;
  Pathpair_local_anchor_rescue_T local_anchor_rescue;
  bool one_low_side_p;
  Pathpair_one_low_side_T low_side;
  int coverage5, coverage3, best_coverage5, best_coverage3;

  List_T paths, fusion_paths_5, fusion_paths_3, p, q;
  List_T unresolved_pathpairs;
  double best_pair_score, accepted_pair_score, pair_score;
  double best_avg_junction_sense_prob, best_avg_junction_antisense_prob,
    best_end_sense_prob, best_end_antisense_prob, prob;
  /* int best_querystart5, best_queryend5, best_querystart3, best_queryend3, querystart, queryend; */
  int chop_qstart_L, chop_qend_L, chop_qstart_H, chop_qend_H, qstart, qend;
  int ignore_found_score = 0, ignore_end_score = 0;

  
#if defined(DEBUG8) || defined(DEBUG_PATH_CORE)
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
  k = 0;
  for (i = 0; i < npaths; i++) {
    pathpair = pathpairarray[i];
    path5 = pathpair->path5;
    path3 = pathpair->path3;
    if (circularp[pathpair->pathL->chrnum] == true) {
      Path_trim_circular_unalias_pair(pathpair->path5,pathpair->path3,
				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev);
    }
    /* true means that endpoints were acceptable */
    if (Path_trim_chrbounds(&ignore_found_score,&ignore_end_score,
			    path5,query5_compress_fwd,query5_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,
			    transcriptpool) == true &&
	Path_trim_chrbounds(&ignore_found_score,&ignore_end_score,
			    path3,query3_compress_fwd,query3_compress_rev,
			    intlistpool,univcoordlistpool,listpool,pathpool,
			    transcriptpool) == true) {
      pathpairarray[k++] = pathpair;
    } else {
      debug8(printf("Eliminating due to chrbounds: ")); debug8(Pathpair_print(pathpair));
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
  }
  
  npaths = k;
  debug8(printf("(1) After trimming chrbounds, have %d paths\n",npaths));
  if (npaths == 0) {
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;
    return (Pathpair_T *) NULL;
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
      coverage5 = Path_coverage(path5);
      coverage3 = Path_coverage(path3);
      if (pathpair_clean_transcript_overlap_p(pathpair,coverage5,coverage3,
					      mincoverage_filter_5,mincoverage_filter_3) == true) {
	/* This short-overlap transcript pair is already scored; trimming can erase valid exonic coverage. */
	debug8(printf("Skipping inner overlap resolve for clean transcript pathpair: %d and %d coverage, outerlength %u\n",
		      coverage5,coverage3,pathpair->outerlength));

      } else if (pathpair->plusp == true) {
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
	  Pathpair_resolve(&(*found_score_5),&(*found_score_3),&(*end_score_5),&(*end_score_3),
			   pathpair,pathpair->plusp,/*genestrand*/0,
			   query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
			   /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,queryuc_ptr_5,queryuc_ptr_3,
			   /*querylengthL*/querylength5,/*querylengthH*/querylength3,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   localdb_nmismatches_allowed,
			   
			   stage1_5,stage1_3,knownsplicing,
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
	  Pathpair_resolve(&(*found_score_5),&(*found_score_3),&(*end_score_5),&(*end_score_3),
			   pathpair,pathpair->plusp,/*genestrand*/0,
			   query3_compress_rev,query3_compress_fwd,query3_compress_rev,
			   query5_compress_rev,query5_compress_fwd,query5_compress_rev,
			   /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,queryrc3,queryrc5,
			   /*querylengthL*/querylength3,/*querylengthH*/querylength5,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   localdb_nmismatches_allowed,

			   stage1_3,stage1_5,knownsplicing,
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
      FREE_OUT(pathpairarray);
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
      
#ifdef FUSION
      /* Look for fusion paths */
      if (pathpair->plusp == true) {
	/* plus */
	if (Path_unextended_qstart_p(pathpair->path5,/*endtrim_allowed*/INNER_NEEDED,/*allow_ambig_p*/false) == false) {
	  fusion_paths_5 = (List_T) NULL;
	} else if (pathpair->path5->sensedir == SENSE_ANTI) {
	  /* Skip: requires sense */
	  fusion_paths_5 = (List_T) NULL;
	} else {
	  fusion_paths_5 = Path_fusion_outer_querystart_plus(&(*found_score_5),&(*end_score_5),
							     /*main*/pathpair->path5,stage1_5,

							     query5_compress_fwd,query5_compress_rev,
							     queryseq5,querylength5,knownsplicing,
							     pathpair->path5->genestrand,/*nmismatches_allowed:not used*/-1,
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
	  fusion_paths_3 = Path_fusion_outer_queryend_plus(&(*found_score_3),&(*end_score_3),
							   /*main*/pathpair->path3,stage1_3,

							   query3_compress_fwd,query3_compress_rev,
							   queryseq3,querylength3,knownsplicing,
							   pathpair->path3->genestrand,/*nmismatches_allowed:not used*/-1,
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
	  fusion_paths_3 = Path_fusion_outer_querystart_minus(&(*found_score_3),&(*end_score_3),
							      /*main*/pathpair->path3,stage1_3,

							      query3_compress_fwd,query3_compress_rev,
							      queryseq3,querylength3,knownsplicing,
							      pathpair->path3->genestrand,/*nmismatches_allowed:not used*/-1,
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
	  fusion_paths_5 = Path_fusion_outer_queryend_minus(&(*found_score_5),&(*end_score_5),
							    /*main*/pathpair->path5,stage1_5,

							    query5_compress_fwd,query5_compress_rev,
							    queryseq5,querylength5,knownsplicing,
							    pathpair->path5->genestrand,/*nmismatches_allowed:not used*/-1,
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
#endif	/* FUSION */
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
      temp = pathpairarray[i];
      debug8(printf("(5 best) Keeping by local_cmp ")); debug8(Pathpair_print(temp));
      pathpairarray[k++] = temp;

      if (mincoverage_filter_per_read_p == true) {
	best_coverage5 = Path_coverage(temp->path5);
	best_coverage3 = Path_coverage(temp->path3);
	best_low_coverage_p =
	  (best_coverage5 < mincoverage_filter_5 ||
	   best_coverage3 < mincoverage_filter_3) ? true : false;
	best_terminal_anchor_rescue_p =
	  pathpair_terminal_anchor_rescue_p(temp,best_coverage5,best_coverage3,
					    mincoverage_filter_5,mincoverage_filter_3,
					    querylength5,querylength3);
	best_postfallback_terminal_anchor_rescue_p =
	  pathpair_postfallback_terminal_anchor_rescue_p(temp,best_coverage5,best_coverage3,
							 mincoverage_filter_5,mincoverage_filter_3,
							 querylength5,querylength3);
	best_low_coverage_side_nonexact_p =
	  pathpair_low_coverage_side_nonexact_p(temp,best_coverage5,best_coverage3,
					       mincoverage_filter_5,mincoverage_filter_3);
      } else {
	best_low_coverage_p = best_low_coverage_side_nonexact_p =
	  best_terminal_anchor_rescue_p =
	  best_postfallback_terminal_anchor_rescue_p = false;
      }
	
      for (l = i + 1; l < j; l++) {
	pathpair = pathpairarray[l];
	structure_cmp = Pathpair_structure_ignore_sense_cmp(&pathpair,&temp);
	if (structure_cmp == 0) {
	  local_anchor_rescue =
	    pathpair_local_terminal_anchor_rescue(&coverage5,&coverage3,pathpair,
						  best_low_coverage_p,
						  best_low_coverage_side_nonexact_p,
						  best_terminal_anchor_rescue_p,
						  best_postfallback_terminal_anchor_rescue_p,
						  /*require_best_low_coverage_side_nonexact_p*/true,
						  /*require_low_coverage_side_substantial_p*/true,
						  mincoverage_filter_5,mincoverage_filter_3,
						  querylength5,querylength3,
						  terminal_anchor_pathpairs_deferred,
						  allow_postfallback_terminal_anchor_rescue_p);
	  if (local_anchor_rescue == PATHPAIR_LOCAL_ANCHOR_RESCUE_TERMINAL) {
	    debug8(printf("(5 rescue) Keeping paired terminal-anchor rescue before identical-structure deletion: %d and %d coverage ",
			  coverage5,coverage3));
	    debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;

	  } else if (local_anchor_rescue == PATHPAIR_LOCAL_ANCHOR_RESCUE_POSTFALLBACK) {
	    debug8(printf("(5 rescue) Keeping post-fallback terminal-anchor rescue before identical-structure deletion: %d and %d coverage ",
			  coverage5,coverage3));
	    debug8(Pathpair_print(pathpair));
	    pathpairarray[k++] = pathpair;

	  } else {
	    debug8(printf("(5 identical) Eliminating by identical structure ")); debug8(Pathpair_print(pathpair));
	    Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
	  }
	  
	} else if (Pathpair_local_cmp(&pathpair,&temp) == 0) {
	  debug8(printf("(5 tie) Keeping by local_cmp ")); debug8(Pathpair_print(pathpair));
	  pathpairarray[k++] = pathpair;

	} else if ((local_anchor_rescue =
		    pathpair_local_terminal_anchor_rescue(&coverage5,&coverage3,pathpair,
							  best_low_coverage_p,
							  best_low_coverage_side_nonexact_p,
							  best_terminal_anchor_rescue_p,
							  best_postfallback_terminal_anchor_rescue_p,
							  /*require_best_low_coverage_side_nonexact_p*/false,
							  /*require_low_coverage_side_substantial_p*/false,
							  mincoverage_filter_5,mincoverage_filter_3,
							  querylength5,querylength3,
							  terminal_anchor_pathpairs_deferred,
							  allow_postfallback_terminal_anchor_rescue_p)) ==
		   PATHPAIR_LOCAL_ANCHOR_RESCUE_TERMINAL) {
	  debug8(printf("(5 rescue) Keeping paired terminal-anchor rescue over low-coverage local_cmp best: %d and %d coverage ",
			coverage5,coverage3));
	  debug8(Pathpair_print(pathpair));
	  pathpairarray[k++] = pathpair;

	} else if (local_anchor_rescue == PATHPAIR_LOCAL_ANCHOR_RESCUE_POSTFALLBACK) {
	  debug8(printf("(5 rescue) Keeping post-fallback terminal-anchor rescue over low-coverage local_cmp best: %d and %d coverage ",
			coverage5,coverage3));
	  debug8(Pathpair_print(pathpair));
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


  /* 6.  Filter the results */
  if (filterp == true) {
    rescue_pathpairarray = (Pathpair_T *) MALLOC_OUT(npaths * sizeof(Pathpair_T));
    rescue_npaths = 0;
    k = 0;
    for (i = 0; i < npaths; i++) {
      pathpair = pathpairarray[i];
      
      if (pathpair->plusp == true) {
	allow_soft_clips_qstart_p = allow_soft_clips_start_p;
	allow_soft_clips_qend_p = allow_soft_clips_end_p;
      } else {
	allow_soft_clips_qstart_p = allow_soft_clips_end_p;
	allow_soft_clips_qend_p = allow_soft_clips_start_p;
      }
      coverage5 = Path_coverage(pathpair->path5);
      coverage3 = Path_coverage(pathpair->path3);
      one_low_side_p =
	(mincoverage_filter_per_read_p == true &&
	 pathpair_one_low_side_setup(&low_side,pathpair,coverage5,coverage3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     querylength5,querylength3) == true) ? true : false;

      if (pathpair->path5->score_within_trims > nmismatches_filter_5 || pathpair->path3->score_within_trims > nmismatches_filter_3) {
	debug8(printf("(6 filter) Eliminating: solution has too many %d + %d nmismatches\n",
		      pathpair->path5->score_within_trims,pathpair->path3->score_within_trims));
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

      } else if ((mincoverage_filter_per_read_p == true &&
		  (coverage5 < mincoverage_filter_5 || coverage3 < mincoverage_filter_3)) ||
		 (mincoverage_filter_per_read_p == false &&
		  coverage5 < mincoverage_filter_5 && coverage3 < mincoverage_filter_3)) {
	if (one_low_side_p == true &&
	    pathpair_mincoverage_rescue_side_p(&low_side) == true) {
	  debug8(printf("(6 filter) Saving clean paired coverage rescue: %d and %d coverage\n",
			coverage5,coverage3));
	  rescue_pathpairarray[rescue_npaths++] = pathpair;

	} else if (one_low_side_p == true &&
		   pathpair_terminal_anchor_rescue_side_p(&low_side) == true) {
	  debug8(printf("(6 filter) Saving paired terminal-anchor coverage rescue: %d and %d coverage\n",
			coverage5,coverage3));
	  rescue_pathpairarray[rescue_npaths++] = pathpair;

	} else if (allow_postfallback_terminal_anchor_rescue_p == true &&
		   one_low_side_p == true &&
		   pathpair_postfallback_terminal_anchor_rescue_side_p(&low_side) == true) {
	  debug8(printf("(6 filter) Saving post-fallback terminal-anchor coverage rescue: %d and %d coverage, scores %d and %d\n",
			coverage5,coverage3,pathpair->path5->score_within_trims,pathpair->path3->score_within_trims));
	  rescue_pathpairarray[rescue_npaths++] = pathpair;

	} else if (allow_near_exact_terminal_anchor_rescue_p == true &&
		   one_low_side_p == true &&
		   pathpair_near_exact_terminal_anchor_rescue_side_p(&low_side) == true) {
	  debug8(printf("(6 filter) Saving near-exact terminal-anchor coverage rescue: %d and %d coverage, scores %d and %d\n",
			coverage5,coverage3,pathpair->path5->score_within_trims,pathpair->path3->score_within_trims));
	  rescue_pathpairarray[rescue_npaths++] = pathpair;

	} else if (terminal_anchor_pathpairs_deferred != NULL &&
		   one_low_side_p == true &&
		   pathpair_postfallback_terminal_anchor_rescue_side_p(&low_side) == true) {
	  /* Keep this pair out of ranking until final singles show that one mate is absent. */
	  debug8(printf("(6 filter) Deferring post-fallback terminal-anchor coverage rescue: %d and %d coverage, scores %d and %d\n",
			coverage5,coverage3,pathpair->path5->score_within_trims,pathpair->path3->score_within_trims));
	  if (path_single_filter_pass_p(pathpair->path5,coverage5,nmismatches_filter_5,
					mincoverage_filter_5,querylength5,
					allow_soft_clips_qstart_p,allow_soft_clips_qend_p)) {
	    save_rejected_single_path(singlepaths5_from_rejected,pathpair->path5,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  if (path_single_filter_pass_p(pathpair->path3,coverage3,nmismatches_filter_3,
					mincoverage_filter_3,querylength3,
					allow_soft_clips_qstart_p,allow_soft_clips_qend_p)) {
	    save_rejected_single_path(singlepaths3_from_rejected,pathpair->path3,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  *terminal_anchor_pathpairs_deferred =
	    Hitlist_push(*terminal_anchor_pathpairs_deferred,hitlistpool,(void *) pathpair
			 hitlistpool_trace(__FILE__,__LINE__));

	} else {
	  debug8(printf("(6 filter) Eliminating: too little %d and %d coverage\n",
			coverage5,coverage3));
	  /* Save near-exact terminal paths before freeing the rejected pair. */
	  if (one_low_side_p == true &&
	      low_side.low5p == true && terminal_anchor_paths5_deferred != NULL &&
	      pathpair_near_exact_low_side_path_p(&low_side) == true) {
	    debug8(printf("(6 filter) Deferring near-exact 5' terminal-anchor path: coverage %d, score %d\n",
			  coverage5,pathpair->path5->score_within_trims));
	    save_rejected_single_path(terminal_anchor_paths5_deferred,pathpair->path5,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  if (one_low_side_p == true &&
	      low_side.low5p == false && terminal_anchor_paths3_deferred != NULL &&
	      pathpair_near_exact_low_side_path_p(&low_side) == true) {
	    debug8(printf("(6 filter) Deferring near-exact 3' terminal-anchor path: coverage %d, score %d\n",
			  coverage3,pathpair->path3->score_within_trims));
	    save_rejected_single_path(terminal_anchor_paths3_deferred,pathpair->path3,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  if (path_single_filter_pass_p(pathpair->path5,coverage5,nmismatches_filter_5,
					mincoverage_filter_5,querylength5,
					allow_soft_clips_qstart_p,allow_soft_clips_qend_p)) {
	    save_rejected_single_path(singlepaths5_from_rejected,pathpair->path5,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  if (path_single_filter_pass_p(pathpair->path3,coverage3,nmismatches_filter_3,
					mincoverage_filter_3,querylength3,
					allow_soft_clips_qstart_p,allow_soft_clips_qend_p)) {
	    save_rejected_single_path(singlepaths3_from_rejected,pathpair->path3,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
	  }
	  Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
			listpool,pathpool,transcriptpool,hitlistpool);
	}

      } else if (allow_soft_clips_qstart_p == false &&
		 (Intlist_head(pathpair->path5->endpoints) != 0 || Intlist_head(pathpair->path3->endpoints) != 0)) {
	debug8(printf("(6 filter) Eliminating: soft clips %d and %d at qstart; %d and %d at qend\n",
		      Intlist_head(pathpair->path5->endpoints),Intlist_head(pathpair->path5->endpoints),
		      querylength5 - Intlist_last_value(pathpair->path5->endpoints),
		      querylength3 - Intlist_last_value(pathpair->path3->endpoints)));
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

      } else if (allow_soft_clips_qend_p == false &&
		 (Intlist_last_value(pathpair->path5->endpoints) != querylength5 ||
		  Intlist_last_value(pathpair->path3->endpoints) != querylength3)) {
	debug8(printf("(6 filter) Eliminating: soft clips %d and %d at qstart; %d and %d at qend\n",
		      Intlist_head(pathpair->path5->endpoints),Intlist_head(pathpair->path5->endpoints),
		      querylength5 - Intlist_last_value(pathpair->path5->endpoints),
		      querylength3 - Intlist_last_value(pathpair->path3->endpoints)));
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);

      } else {
	debug8(printf("(6 filter) Keeping\n"));
	pathpairarray[k++] = pathpair;
      }
    }

    if (k == 0 && rescue_npaths > 0) {
      debug8(printf("(6 filter) Keeping %d clean paired coverage rescues because no full-coverage pair survived\n",
		    rescue_npaths));
      for (i = 0; i < rescue_npaths; i++) {
	pathpairarray[k++] = rescue_pathpairarray[i];
      }

    } else {
      for (i = 0; i < rescue_npaths; i++) {
	pathpair = rescue_pathpairarray[i];
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
    FREE_OUT(rescue_pathpairarray);

    npaths = k;
  }

  if (npaths == 0) {
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;
    return (Pathpair_T *) NULL;
  }

  /* 7.  Find best solution globally */
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

    qsort(pathpairarray,npaths,sizeof(T),Pathpair_final_score_cmp);

    best_pair_score = Pathpair_final_score(pathpairarray[0]);
    accepted_pair_score = best_pair_score * (1.0 - score_pair_margin);

    debug8(printf("Found the global best solution across %d paths with score %.3f, accepted score %.3f, and %d + %d nmatches\n",
		  npaths,best_pair_score,accepted_pair_score,
		  pathpairarray[0]->path5->nmatches,pathpairarray[0]->path3->nmatches));

    /* Keep all paths inside the best-score margin */
    debug8(printf("(7 best) Keeping best by final score %.3f\n",best_pair_score));
    debug8(Pathpair_print(pathpairarray[0]));
    i = 1;		/* Skip the best path */
      
    for (k = 1; k < npaths; k++) {
      pathpair = pathpairarray[k];
      pair_score = Pathpair_final_score(pathpair);
      if (pair_score >= accepted_pair_score) {
	debug8(printf("(7 margin) Keeping by final score %.3f >= %.3f\n",
		      pair_score,accepted_pair_score));
	debug8(Pathpair_print(pathpair));

	pathpairarray[i++] = pathpair;
      } else {
	debug8(printf("(7 worse) Eliminating by final score %.3f < %.3f\n",
		      pair_score,accepted_pair_score));
	debug8(Pathpair_print(pathpair));
	Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
    }

    npaths = i;			
  }

  for (i = k = 0; i < npaths; i++) {
    pathpair = pathpairarray[i];
    if (Pathpair_final_score_reportable_p(pathpair) == true) {
      pathpairarray[k++] = pathpair;
    } else {
      debug8(printf("(7 nonreportable) Eliminating by final scores read5 %.3f, read3 %.3f, pair %.3f below minima %.3f/%.3f\n",
		    Path_final_score(pathpair->path5),Path_final_score(pathpair->path3),
		    Pathpair_final_score(pathpair),score_min_read,score_min_pair));
      debug8(Pathpair_print(pathpair));
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
  }
  npaths = k;

  if (npaths == 0) {
    FREE_OUT(pathpairarray);
    *npaths_primary = *npaths_altloc = 0;
    *first_absmq = 0;
    *second_absmq = 0;
    return (Pathpair_T *) NULL;
  }

  for (i = 0; i < npaths; i++) {
    Pathpair_set_final_scores(pathpairarray[i]);
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
      Path_trim_circular(pathpair->pathL,/*query_compress*/query5_compress_fwd,
			 query5_compress_fwd,query5_compress_rev,
			 intlistpool,univcoordlistpool,listpool);
    } else {
      Path_trim_circular(pathpair->pathL,/*query_compress*/query3_compress_rev,
			 query3_compress_fwd,query3_compress_rev,
			 intlistpool,univcoordlistpool,listpool);
    }

    if (circularp[pathpair->pathH->chrnum] == false) {
#if 0
      Path_trim_qend_trimbounds(pathpair->pathH,intlistpool,univcoordlistpool,listpool,pathpool,
				/*trimbounds*/pathpair->pathH->chrhigh);
#endif
    } else if (output_type == STD_OUTPUT || output_type == M8_OUTPUT) {
      /* If output type is alignment or m8, then don't want to split up the parts */
    } else if (pathpair->pathH->plusp == true) {
      Path_trim_circular(pathpair->pathH,/*query_compress*/query3_compress_fwd,
			 query3_compress_fwd,query3_compress_rev,
			 intlistpool,univcoordlistpool,listpool);
    } else {
      Path_trim_circular(pathpair->pathH,/*query_compress*/query5_compress_rev,
			 query5_compress_fwd,query5_compress_rev,
			 intlistpool,univcoordlistpool,listpool);
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
		     bool allow_soft_clips_start_p_in, bool allow_soft_clips_end_p_in,
		     int paired_terminal_anchor_min_score_in,
		     int paired_terminal_anchor_slop_in,
		     bool mincoverage_filter_per_read_p_in,
		     bool pathpair_equal_coverage_diff_tiebreak_p_in,
		     double score_match_reward_in,
		     double score_mismatch_penalty_in,
		     double score_indel_open_penalty_in,
		     double score_indel_extend_penalty_in,
		     double score_splice_unannotated_penalty_in,
		     double score_pair_transcript_discordance_penalty_in,
		     double score_pair_margin_in,
		     double score_min_read_in,
		     double score_min_pair_in) {

  max_insertlength = expected_pairlength + 2*pairlength_deviation;
  transcriptome = transcriptome_in;

  circularp = circularp_in;
  chrsubsetp = chrsubsetp_in;
  altlocp = altlocp_in;

  output_type = output_type_in;
  splicingp = splicingp_in;
  resolve_inner_p = resolve_inner_p_in;
  want_random_p = want_random_p_in;

  allow_soft_clips_start_p = allow_soft_clips_start_p_in;
  allow_soft_clips_end_p = allow_soft_clips_end_p_in;
  paired_terminal_anchor_min_score = paired_terminal_anchor_min_score_in;
  paired_terminal_anchor_slop = paired_terminal_anchor_slop_in;
  mincoverage_filter_per_read_p = mincoverage_filter_per_read_p_in;
  pathpair_equal_coverage_diff_tiebreak_p = pathpair_equal_coverage_diff_tiebreak_p_in;
  score_match_reward = score_match_reward_in;
  score_mismatch_penalty = score_mismatch_penalty_in;
  score_indel_open_penalty = score_indel_open_penalty_in;
  score_indel_extend_penalty = score_indel_extend_penalty_in;
  score_splice_unannotated_penalty = score_splice_unannotated_penalty_in;
  score_pair_transcript_discordance_penalty = score_pair_transcript_discordance_penalty_in;
  score_pair_margin = score_pair_margin_in;
  score_min_read = score_min_read_in;
  score_min_pair = score_min_pair_in;

  return;
}
