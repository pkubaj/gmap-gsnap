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
#include "align-cascade.h"
#include "cascade-timing.h"
#include "cascade-stage-probe.h"
#include "accel-profile-stats.h"
#include "debug_groups.h"

/* Exact EXT-anchor root reuse is shelved pending a proven state-level
   equivalence relation (value purity + ordering/scoring parity).  It is enabled
   only with an explicit build flag; default and stats builds neither allocate
   nor consult the anchor cache.  It targets only same-head EXT repeats, a
   minority of overall LocalDB request repetition, so the general LocalDB
   exact-request memoization is the primary redundancy lever. */
#if defined(GSNAP_EXT_ANCHOR_REUSE_ENABLE) && !defined(GSNAP_STATS_LEGACY_REQUEST_GRAPH)
#define GSNAP_EXT_ANCHOR_REUSE 1
#endif
#if defined(GSNAP_EXT_ANCHOR_REUSE) || defined(GSNAP_EXT_ANCHOR_SUBSTITUTION_PROOF)
#define GSNAP_EXT_ANCHOR_CACHE_ACTIVE 1
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>		/* For rint */
#include <string.h>		/* For memset */
#ifdef GSNAP_PAIRED_SOLVE_STATS
#include <time.h>
#endif

#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"
#include "univdiag.h"
#if defined(GSNAP_PAIRED_SOLVE_STATS) || defined(GSNAP_ACCEL_PROFILE_STATS) || \
  defined(GSNAP_EXT_ANCHOR_SUBSTITUTION_PROOF) || defined(GSNAP_EXT_ANCHOR_REUSE) || \
  defined(GSNAP_EXT_CHAIN_PREFIX_STATS)
#include "univdiagdef.h"
#endif

#include "list.h"
#include "compress.h"

#include "genomebits_mismatches.h" /* For MISMATCH_EXTRA */
#include "genomebits_kmer.h" /* For MISMATCH_EXTRA */

#include "stage1hr.h"
#include "kmer-search.h"
#include "extension-search.h"
#include "stage1hr-single.h"
#include "transcriptome-search.h"
#include "tr-extension-search.h"
#include "tr-chain-search.h"

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



/* Also defined in stage1hr-single.c */
#define EXON_LENGTH 50		/* Only 5% of exons are less than 50 bp */

/* Mutual concordance is greedy and can miss optimal solutions, but
   needed when there are too many pairs */
#define MAX_OPTIMAL 100

#define LOCALDB_REGION_SIZE 65536
#define QUERYLENGTH_FOR_LOCALDB_MATE 50
#define LOCAL_MATE_MIN_SOFTCLIP 8
#define LOCAL_MATE_MAX_SOFTCLIP 20
#define LOCAL_MATE_MIN_ALIGNED 50
#define LOCAL_MATE_MAX_INDEL 3
#define LOCAL_MATE_MAX_TOTAL_INDEL 3
#define LOCAL_MATE_MIN_INDEL_SUPPORT 12


void
Stage1hr_paired_ext_isolated_flush (Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				     Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
				     Pathpool_T pathpool, Transcriptpool_T transcriptpool,
				     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
				     Hitlistpool_T hitlistpool) {
  (void) intlistpool; (void) uintlistpool; (void) univcoordlistpool;
  (void) listpool; (void) pathpool; (void) transcriptpool;
  (void) univdiagpool; (void) vectorpool; (void) hitlistpool;
}

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


#if defined(DEBUG0) || defined(DEBUG_PAIRED_FLOW)
#define debug0(x) x
#else
#define debug0(x)
#endif

#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
#define debug(x) x
#else
#define debug(x)
#endif

/* Trdiagonals */
#if defined(DEBUG1) || defined(DEBUG_PAIRED_DEEP)
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
static int paired_max_insertlength;
static int paired_segment_exact_max_mismatch_divisor =
  STAGE1HR_PAIRED_DEFAULT_SEGMENT_EXACT_MAX_MISMATCH_DIVISOR;
static bool paired_mate_rescue_gapped_p = false;
static bool paired_preeval_ungapped_mate_p = false;
static bool paired_preeval_rescue_skip_p = false;
static bool paired_segment_stream_chains_p = true;

#ifdef GSNAP_EXT_CHAIN_PREFIX_STATS
#define EXT_PREFIX_SEEN_SIZE (1U << 22)
#define EXT_PREFIX_SEEN_MASK (EXT_PREFIX_SEEN_SIZE - 1U)
static unsigned long long *ext_prefix_seen = NULL;
static unsigned long long ext_prefix_nodes = 0ULL;
static unsigned long long ext_prefix_unique = 0ULL;
static unsigned long long ext_prefix_repeats = 0ULL;
static unsigned long long ext_prefix_table_full = 0ULL;
static bool ext_prefix_report_registered_p = false;

static void
ext_prefix_report (void) {
  fprintf(stderr,"ext_chain_prefix_summary\tnodes=%llu\tunique=%llu\trepeats=%llu\t"
	  "table_full=%llu\n",ext_prefix_nodes,ext_prefix_unique,ext_prefix_repeats,
	  ext_prefix_table_full);
}
#endif
static bool paired_solve_all_anchor_first_p = false;
static bool paired_postsolve_reliable_anchors_first_p = false;
static bool paired_seg1_exhaustive_concordance_p = false;
static bool paired_duplicate_auxinfo_frontier_p = false;
static unsigned int paired_mate_minpure_floor = 20U;
static bool paired_transcript_completion_p = false;
static bool paired_tr_adjacent_indel_repair_p = false;
static bool tr_chain_search_p = false;
static double se_segment_min_kmer_avg_quality = -1.0;
static int se_segment_max_kmer_bad_bases = -1;
static int se_segment_bad_base_quality = 10;
static int segment_max_avg_kmer_positions = 0;
static Align_cascade_T align_cascade = ALIGN_CASCADE_DEFAULT;


#define T Stage1_T

static bool
segment_quality_filter_enabled_p (void) {
  return se_segment_min_kmer_avg_quality >= 0.0 || se_segment_max_kmer_bad_bases >= 0;
}


static bool
segment_kmer_quality_bad_p (char *quality_string, int querypos) {
  int i, q, nbad = 0, total = 0;
  double avg_quality;

  for (i = 0; i < index1part; i++) {
    q = (int) quality_string[querypos + i] - 33;
    total += q;
    if (se_segment_max_kmer_bad_bases >= 0 && q < se_segment_bad_base_quality) {
      nbad++;
    }
  }

  if (se_segment_max_kmer_bad_bases >= 0 && nbad > se_segment_max_kmer_bad_bases) {
    return true;
  }

  if (se_segment_min_kmer_avg_quality >= 0.0) {
    avg_quality = (double) total/(double) index1part;
    if (avg_quality < se_segment_min_kmer_avg_quality) {
      return true;
    }
  }

  return false;
}


static void
segment_apply_quality_filter (T this, char *quality_string, int querylength) {
  int query_lastpos, querypos, nfiltered = 0;

  if (segment_quality_filter_enabled_p() == false) {
    return;
  } else if (quality_string == NULL) {
    return;
  }

  query_lastpos = querylength - index1part;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->validp[querypos] == true &&
	segment_kmer_quality_bad_p(quality_string,querypos) == true) {
      this->validp[querypos] = false;
      this->plus_npositions[querypos] = 0;
      this->minus_npositions[querypos] = 0;
      nfiltered++;
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,
	  "[SEGMENT_QUALITY_FILTER] filtered_kmers=%d query_kmers=%d min_avg=%g max_bad=%d bad_base_q=%d\n",
	  nfiltered,query_lastpos + 1,se_segment_min_kmer_avg_quality,
	  se_segment_max_kmer_bad_bases,se_segment_bad_base_quality);
#endif

  return;
}

static int
segment_quality_seedable_querylength (char *quality_string, int querylength) {
  int query_lastpos, querypos;
  int interval_start = -1, interval_end = -1;
  int kmer_start, kmer_end, covered = 0;

  if (segment_quality_filter_enabled_p() == false) {
    return querylength;
  } else if (quality_string == NULL) {
    return querylength;
  } else if (querylength < index1part) {
    return querylength;
  }

  query_lastpos = querylength - index1part;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (segment_kmer_quality_bad_p(quality_string,querypos) == false) {
      kmer_start = querypos;
      kmer_end = querypos + index1part;

      if (interval_start < 0) {
	interval_start = kmer_start;
	interval_end = kmer_end;
      } else if (kmer_start <= interval_end) {
	if (kmer_end > interval_end) {
	  interval_end = kmer_end;
	}
      } else {
	covered += interval_end - interval_start;
	interval_start = kmer_start;
	interval_end = kmer_end;
      }
    }
  }

  if (interval_start >= 0) {
    covered += interval_end - interval_start;
  }

  return covered;
}

#define CASCADE_TIMING_START(timing) ((timing) != NULL ? Cascade_timing_now_us() : 0ULL)
#define CASCADE_TIMING_ADD(timing,field,start_us) \
  do { \
    if ((timing) != NULL) { \
      (timing)->field += Cascade_timing_now_us() - (start_us); \
    } \
  } while (0)

#define CASCADE_STAGE_PROBE_SET(stage,queryseq,role) \
  GSNAP_univdiagpool_stage_probe_set((stage),(queryseq) == NULL ? "unset" : Shortread_accession(queryseq),(role))
#define CASCADE_STAGE_PROBE_CLEAR() \
  GSNAP_univdiagpool_stage_probe_clear()

static unsigned long long
paired_cascade_stage_start (const char *stage, Shortread_T queryseq,
			    const char *role, Univdiagpool_T univdiagpool) {
  const char *qname = queryseq == NULL ? "unset" : Shortread_accession(queryseq);

  CASCADE_STAGE_PROBE_SET(stage,queryseq,role);
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
  int univdiag_chunks = 0, univdiag_list_chunks = 0;
  if (univdiagpool != NULL) {
    Univdiagpool_nchunks(&univdiag_chunks,&univdiag_list_chunks,univdiagpool);
  }
  fprintf(stderr,
	  "cascade_stage\tbegin\tstage=%s\trole=%s\tqname=%s\tunivdiag_chunks=%d\tunivdiag_list_chunks=%d\n",
	  stage,role,qname == NULL ? "unset" : qname,univdiag_chunks,univdiag_list_chunks);
#endif
  return Cascade_timing_now_us();
}

static unsigned long long
paired_cascade_stage_end (const char *stage, Shortread_T queryseq,
			  const char *role, Univdiagpool_T univdiagpool,
			  unsigned long long start_us) {
  const char *qname = queryseq == NULL ? "unset" : Shortread_accession(queryseq);
  unsigned long long elapsed_us = Cascade_timing_now_us() - start_us;

#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
  int univdiag_chunks = 0, univdiag_list_chunks = 0;
  if (univdiagpool != NULL) {
    Univdiagpool_nchunks(&univdiag_chunks,&univdiag_list_chunks,univdiagpool);
  }
  fprintf(stderr,
	  "cascade_stage\tend\tstage=%s\trole=%s\tqname=%s\telapsed_us=%llu\tunivdiag_chunks=%d\tunivdiag_list_chunks=%d\n",
	  stage,role,qname == NULL ? "unset" : qname,elapsed_us,
	  univdiag_chunks,univdiag_list_chunks);
#endif
  CASCADE_STAGE_PROBE_CLEAR();
  return elapsed_us;
}

#define CASCADE_STAGE_START(stage,queryseq,univdiagpool) \
  paired_cascade_stage_start((stage),(queryseq),"pair",(univdiagpool))
#define CASCADE_STAGE_END(stage,queryseq,univdiagpool,start_us) \
  paired_cascade_stage_end((stage),(queryseq),"pair",(univdiagpool),(start_us))
#define CASCADE_STAGE_START_ROLE(stage,role,queryseq,univdiagpool) \
  paired_cascade_stage_start((stage),(queryseq),(role),(univdiagpool))
#define CASCADE_STAGE_END_ROLE(stage,role,queryseq,univdiagpool,start_us) \
  paired_cascade_stage_end((stage),(queryseq),(role),(univdiagpool),(start_us))
#define CASCADE_STAGE_ADD(timing,field,stage,queryseq,univdiagpool,start_us) \
  do { \
    unsigned long long _elapsed_us = \
      CASCADE_STAGE_END((stage),(queryseq),(univdiagpool),(start_us)); \
    if ((timing) != NULL) { \
      (timing)->field += _elapsed_us; \
    } \
  } while (0)
#define CASCADE_STAGE_ADD_ROLE(timing,field,stage,role,queryseq,univdiagpool,start_us) \
  do { \
    unsigned long long _elapsed_us = \
      CASCADE_STAGE_END_ROLE((stage),(role),(queryseq),(univdiagpool),(start_us)); \
    if ((timing) != NULL) { \
      (timing)->field += _elapsed_us; \
    } \
  } while (0)

static void
paired_cascade_add_gen_method (Cascade_Timing_T *timing, Method_T method,
			       unsigned long long elapsed_us) {
  if (timing == NULL) {
    return;
  } else {
    switch (method) {
    case KMER_EXACT1: case KMER_ANYPAIR: case KMER_EXACT2:
      timing->gen_exact_us += elapsed_us;
      break;
    case EXT: case EXT_NOSPLICE: case EXT_SPLICE:
      timing->gen_extension_us += elapsed_us;
      break;
    case SEGMENT1:
      timing->gen_segment_us += elapsed_us;
      break;
    default:
      timing->gen_other_us += elapsed_us;
      break;
    }
  }
  return;
}

static const char *
paired_cascade_next_gen_stage_label (Method_T last_method, Method_T method_goal) {
  if (last_method < KMER_EXACT1 && method_goal >= KMER_EXACT1 &&
      Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXACT)) {
    return "gex";
  } else if (last_method < EXT && method_goal >= EXT &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT)) {
    return "gext";
  } else if (last_method < SEGMENT1 && method_goal >= SEGMENT1 &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_SEGMENT)) {
    return "gseg";
  } else {
    return "goth";
  }
}

#ifdef GSNAP_PAIRED_SOLVE_STATS
#ifndef GSNAP_PAIRED_SOLVE_STATS_MAX_PATHS
#define GSNAP_PAIRED_SOLVE_STATS_MAX_PATHS 80
#endif

static double
paired_solve_stats_now (void) {
  struct timespec ts;

  clock_gettime(CLOCK_MONOTONIC,&ts);
  return (double) ts.tv_sec + (double) ts.tv_nsec * 1.0e-9;
}

static void
paired_solve_stats_print_pools (const char *label, const char *qname,
				Listpool_T listpool, Univdiagpool_T univdiagpool,
				Pathpool_T pathpool, Hitlistpool_T hitlistpool) {
  int list_chunks, univdiag_chunks, univdiag_list_chunks;
  int path_chunks, junction_chunks, altsplice_chunks, string_chunks;
  int hitlist_chunks;

  list_chunks = Listpool_nchunks(listpool);
  Univdiagpool_nchunks(&univdiag_chunks,&univdiag_list_chunks,univdiagpool);
  Pathpool_nchunks(&path_chunks,&junction_chunks,&altsplice_chunks,&string_chunks,pathpool);
  hitlist_chunks = Hitlistpool_nchunks(hitlistpool);

  fprintf(stderr,
	  "paired_solve_stats\tpools\tstage=%s\tqname=%s\tlist=%d\tunivdiag=%d\t"
	  "univdiag_list=%d\tpath=%d\tjunction=%d\taltsplice=%d\t"
	  "string=%d\thitlist=%d\n",
	  label,qname,list_chunks,univdiag_chunks,univdiag_list_chunks,
	  path_chunks,junction_chunks,altsplice_chunks,string_chunks,hitlist_chunks);
  return;
}

static int
paired_solve_stats_count_pathstore_paths (Pathstore_T *pathstores, int n) {
  int i, count = 0;

  for (i = 0; i < n; i++) {
    count += List_length(pathstores[i]->best_sense_paths);
    count += List_length(pathstores[i]->best_antisense_paths);
  }

  return count;
}

static void
paired_solve_stats_print_path (const char *label, Path_T path, int *printed) {
  Univcoord_T low, high, chroffset;
  unsigned long long rel_low, rel_high;

  if (*printed >= GSNAP_PAIRED_SOLVE_STATS_MAX_PATHS) {
    return;
  } else {
    low = Path_genomiclow(path);
    high = Path_genomichigh(path);
    chroffset = Path_chroffset(path);
    rel_low = (low >= chroffset) ? (unsigned long long) (low - chroffset + 1) : 0ULL;
    rel_high = (high >= chroffset) ? (unsigned long long) (high - chroffset) : 0ULL;

    fprintf(stderr,
	    "paired_solve_stats\tpath\t%s\tchrnum=%d\trel=%llu..%llu\tabs=%llu..%llu\tplus=%d\tq=%d..%d\tfound=%d\tend=%d\tnmatches=%d\tmethod=%s\tcomplete=%d\ttr=%d\n",
	    label,path->chrnum,rel_low,rel_high,
	    (unsigned long long) low,(unsigned long long) high,
	    path->plusp,Path_qstart(path),Path_qend(path),
	    path->found_score,path->end_score,path->nmatches,
	    Method_string(path->method),path->completep,path->transcriptome_method_p);
    *printed += 1;
  }

  return;
}

static void
paired_solve_stats_print_pathstores (const char *label, Pathstore_T *pathstores, int n) {
  List_T p;
  int i, printed = 0, total = 0;

  for (i = 0; i < n; i++) {
    for (p = pathstores[i]->best_sense_paths; p != NULL; p = List_next(p)) {
      total++;
      paired_solve_stats_print_path(label,(Path_T) List_head(p),&printed);
    }
    for (p = pathstores[i]->best_antisense_paths; p != NULL; p = List_next(p)) {
      total++;
      paired_solve_stats_print_path(label,(Path_T) List_head(p),&printed);
    }
  }

  if (total > printed) {
    fprintf(stderr,"paired_solve_stats\tpath\t%s\tomitted=%d\tprinted=%d\tlimit=%d\n",
	    label,total - printed,printed,GSNAP_PAIRED_SOLVE_STATS_MAX_PATHS);
  }

  return;
}

static void
paired_solve_stats_print_path_list (const char *label, List_T paths) {
  List_T p;
  int printed = 0, total = 0;

  for (p = paths; p != NULL; p = List_next(p)) {
    total++;
    paired_solve_stats_print_path(label,(Path_T) List_head(p),&printed);
  }

  if (total > printed) {
    fprintf(stderr,"paired_solve_stats\tpath\t%s\tomitted=%d\tprinted=%d\tlimit=%d\n",
	    label,total - printed,printed,GSNAP_PAIRED_SOLVE_STATS_MAX_PATHS);
  }

  return;
}

static void
paired_solve_stats_print_pathpair_list (const char *label, List_T pathpairs) {
  List_T p;
  Pathpair_T pathpair;
  Path_T path5, path3;
  Univcoord_T low5, high5, offset5, low3, high3, offset3;
  unsigned long long rel5_low, rel5_high, rel3_low, rel3_high;
  int i = 0;

  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (Pathpair_T) List_head(p);
    path5 = pathpair->path5;
    path3 = pathpair->path3;
    low5 = Path_genomiclow(path5);
    high5 = Path_genomichigh(path5);
    offset5 = Path_chroffset(path5);
    low3 = Path_genomiclow(path3);
    high3 = Path_genomichigh(path3);
    offset3 = Path_chroffset(path3);
    rel5_low = (low5 >= offset5) ? (unsigned long long) (low5 - offset5 + 1) : 0ULL;
    rel5_high = (high5 >= offset5) ? (unsigned long long) (high5 - offset5) : 0ULL;
    rel3_low = (low3 >= offset3) ? (unsigned long long) (low3 - offset3 + 1) : 0ULL;
    rel3_high = (high3 >= offset3) ? (unsigned long long) (high3 - offset3) : 0ULL;

    fprintf(stderr,
	    "paired_solve_stats\tpathpair\t%s\tidx=%d\tpairtype=%d\trelationship=%d\tplus=%d\tinsert=%u\touter=%u\tfound=%d\tend=%d\tchr5=%d\trel5=%llu..%llu\tplus5=%d\tq5=%d..%d\tfound5=%d\tend5=%d\tnm5=%d\tmethod5=%s\tchr3=%d\trel3=%llu..%llu\tplus3=%d\tq3=%d..%d\tfound3=%d\tend3=%d\tnm3=%d\tmethod3=%s\n",
	    label,i,pathpair->pairtype,pathpair->pair_relationship,pathpair->plusp,
	    pathpair->insertlength,pathpair->outerlength,
	    Pathpair_found_score(pathpair),Pathpair_end_score(pathpair),
	    path5->chrnum,rel5_low,rel5_high,path5->plusp,Path_qstart(path5),Path_qend(path5),
	    path5->found_score,path5->end_score,path5->nmatches,Method_string(path5->method),
	    path3->chrnum,rel3_low,rel3_high,path3->plusp,Path_qstart(path3),Path_qend(path3),
	    path3->found_score,path3->end_score,path3->nmatches,Method_string(path3->method));
    i++;
  }

  return;
}

static void
paired_solve_stats_print_stage (const char *label, T this) {
  fprintf(stderr,
	  "paired_solve_stats\tstage\t%s\tunivdiag_plus=%d\tunivdiag_minus=%d\texhaustive_plus=%d\texhaustive_minus=%d\tbest_plus=%d\tbest_minus=%d\n",
	  label,this->all_nunivdiagonals_gplus,this->all_nunivdiagonals_gminus,
	  this->nexhaustive_gplus,this->nexhaustive_gminus,
	  paired_solve_stats_count_pathstore_paths(this->all_pathstores_gplus,this->all_nunivdiagonals_gplus),
	  paired_solve_stats_count_pathstore_paths(this->all_pathstores_gminus,this->all_nunivdiagonals_gminus));

  return;
}
#endif

#ifdef GSNAP_PAIRED_SOLVE_STATS
typedef struct Mate_rescue_stats_T {
  const char *stage;

  unsigned long long threshold_anchor_scanned;
  unsigned long long threshold_anchor_complete_skip;
  unsigned long long threshold_anchor_rule_reject;
  unsigned long long threshold_anchor_best;
  unsigned long long threshold_anchor_with_mate_count;

  unsigned long long rescue_anchor_scanned;
  unsigned long long rescue_anchor_complete_skip;
  unsigned long long rescue_anchor_rule_reject;
  unsigned long long rescue_anchor_best;
  unsigned long long rescue_anchor_threshold_pass;

  unsigned long long auxinfo_null;
  unsigned long long mate_bounds_skip;
  unsigned long long pathstore_mate_created;
  unsigned long long segment_solve_calls;
  unsigned long long chain_solve_calls;
  unsigned long long bounded_calls;
  unsigned long long bounded_created;
  unsigned long long bounded_replaced;
  unsigned long long mate_solved;
  unsigned long long make_pathpairs_calls;
  unsigned long long local_mate_singlepath_copy_calls;
  unsigned long long local_mate_singlepath_copy_skipped_unusable;
  unsigned long long local_mate_singlepaths_copied;
  unsigned long long segment_control_enabled;
  unsigned long long segment_control_dedup_activated;
  unsigned long long segment_control_cap_hit;
  unsigned long long segment_control_chains_seen;
  unsigned long long segment_control_chains_solved;
  unsigned long long segment_control_chains_skipped;

  long long pathpairs_delta;
  long long unresolved_delta;
  unsigned long long next_checkpoint_work;
} Mate_rescue_stats_T;

#define MATE_RESCUE_STATS_CHECKPOINT_STEP 4096ULL

#define MATE_RESCUE_STATS_INC(stats,field) \
  do { if ((stats) != NULL) { (stats)->field++; } } while (0)

#define MATE_RESCUE_STATS_ADD(stats,field,value) \
  do { if ((stats) != NULL) { (stats)->field += (value); } } while (0)

static void
mate_rescue_stats_init (Mate_rescue_stats_T *stats, const char *stage) {
  memset(stats,0,sizeof(*stats));
  stats->stage = stage;
  stats->next_checkpoint_work = MATE_RESCUE_STATS_CHECKPOINT_STEP;
  return;
}

static unsigned long long
mate_rescue_stats_work (Mate_rescue_stats_T *stats) {
  return stats->rescue_anchor_scanned + stats->pathstore_mate_created +
    stats->segment_solve_calls + stats->chain_solve_calls +
    stats->bounded_calls + stats->make_pathpairs_calls +
    stats->local_mate_singlepath_copy_calls +
    stats->local_mate_singlepaths_copied +
    stats->segment_control_chains_seen;
}

static void
mate_rescue_stats_checkpoint (Mate_rescue_stats_T *stats, const char *side,
			      const char *direction, bool forcep) {
  unsigned long long work;

  if (stats == NULL) {
    return;
  } else {
    work = mate_rescue_stats_work(stats);
    if (forcep == false && work < stats->next_checkpoint_work) {
      return;
    } else if (forcep == false) {
      while (stats->next_checkpoint_work <= work) {
	stats->next_checkpoint_work += MATE_RESCUE_STATS_CHECKPOINT_STEP;
      }
    }

    fprintf(stderr,
	    "paired_solve_stats\tmate_rescue_progress\tstage=%s\tside=%s\tdirection=%s\twork=%llu\t"
	    "threshold_scanned=%llu\tthreshold_complete_skip=%llu\tthreshold_rule_reject=%llu\t"
	    "threshold_best=%llu\tthreshold_with_mate_count=%llu\t"
	    "rescue_scanned=%llu\trescue_complete_skip=%llu\trescue_rule_reject=%llu\t"
	    "rescue_best=%llu\trescue_threshold_pass=%llu\t"
	    "auxinfo_null=%llu\tmate_bounds_skip=%llu\tpathstore_mate_created=%llu\t"
	    "segment_solve_calls=%llu\tchain_solve_calls=%llu\tbounded_calls=%llu\t"
	    "bounded_created=%llu\tbounded_replaced=%llu\tmate_solved=%llu\t"
	    "make_pathpairs_calls=%llu\tpathpairs_delta=%lld\tunresolved_delta=%lld\t"
	    "local_mate_singlepath_copy_calls=%llu\t"
	    "local_mate_singlepath_copy_skipped_unusable=%llu\t"
	    "local_mate_singlepaths_copied=%llu\t"
	    "segment_control_enabled=%llu\tsegment_control_dedup_activated=%llu\t"
	    "segment_control_cap_hit=%llu\tsegment_control_chains_seen=%llu\t"
	    "segment_control_chains_solved=%llu\tsegment_control_chains_skipped=%llu\n",
	    stats->stage,side,direction,work,
	    stats->threshold_anchor_scanned,stats->threshold_anchor_complete_skip,
	    stats->threshold_anchor_rule_reject,stats->threshold_anchor_best,
	    stats->threshold_anchor_with_mate_count,
	    stats->rescue_anchor_scanned,stats->rescue_anchor_complete_skip,
	    stats->rescue_anchor_rule_reject,stats->rescue_anchor_best,
	    stats->rescue_anchor_threshold_pass,
	    stats->auxinfo_null,stats->mate_bounds_skip,stats->pathstore_mate_created,
	    stats->segment_solve_calls,stats->chain_solve_calls,stats->bounded_calls,
	    stats->bounded_created,stats->bounded_replaced,stats->mate_solved,
	    stats->make_pathpairs_calls,stats->pathpairs_delta,stats->unresolved_delta,
	    stats->local_mate_singlepath_copy_calls,
	    stats->local_mate_singlepath_copy_skipped_unusable,
	    stats->local_mate_singlepaths_copied,
	    stats->segment_control_enabled,stats->segment_control_dedup_activated,
	    stats->segment_control_cap_hit,stats->segment_control_chains_seen,
	    stats->segment_control_chains_solved,stats->segment_control_chains_skipped);
    return;
  }
}

static void
mate_rescue_stats_add_pathpair_delta (Mate_rescue_stats_T *stats,
				      int before_pathpairs, int after_pathpairs,
				      int before_unresolved, int after_unresolved) {
  if (stats != NULL) {
    stats->make_pathpairs_calls++;
    stats->pathpairs_delta += (long long) after_pathpairs - (long long) before_pathpairs;
    stats->unresolved_delta += (long long) after_unresolved - (long long) before_unresolved;
  }
  return;
}

static void
mate_rescue_stats_note_bounded_result (Mate_rescue_stats_T *stats,
				       Pathstore_T before, Pathstore_T after) {
  if (stats != NULL) {
    stats->bounded_calls++;
    if (before == NULL && after != NULL) {
      stats->bounded_created++;
    } else if (before != NULL && after != NULL && before != after) {
      stats->bounded_replaced++;
    }
  }
  return;
}
#else
typedef void Mate_rescue_stats_T;
#define MATE_RESCUE_STATS_INC(stats,field) do { } while (0)
#define MATE_RESCUE_STATS_ADD(stats,field,value) do { } while (0)
#define mate_rescue_stats_checkpoint(stats,side,direction,forcep) do { } while (0)
#define mate_rescue_stats_add_pathpair_delta(stats,bp,ap,bu,au) do { } while (0)
#define mate_rescue_stats_note_bounded_result(stats,before,after) do { } while (0)
#endif

#define MATE_RESCUE_SEGMENT_DEDUP_CHAIN_TRIGGER 4096ULL
#define MATE_RESCUE_SEGMENT_DEDUP_COMPLETE_TRIGGER 8192
#define MATE_RESCUE_SEGMENT_SOFT_CHAIN_LIMIT 16384ULL
#define MATE_RESCUE_SEGMENT_HARD_CHAIN_LIMIT 65536ULL
#define MATE_RESCUE_SEGMENT_SOFT_COMPLETE_LIMIT 8192
#define MATE_RESCUE_SEGMENT_HARD_COMPLETE_LIMIT 32768
#define MATE_RESCUE_SEGMENT_CHECK_INTERVAL 256ULL

typedef struct Mate_rescue_segment_policy_T {
  bool activep;
  bool dedup_activated_p;
  bool cap_hit_p;
  bool stop_solving_p;
  const char *stage;
  const char *side;
  const char *direction;
  unsigned long long chains_seen;
  unsigned long long chains_solved;
  unsigned long long chains_skipped;
  unsigned long long dedup_chain_trigger;
  int dedup_complete_trigger;
  unsigned long long soft_chain_limit;
  unsigned long long hard_chain_limit;
  int soft_complete_path_limit;
  int hard_complete_path_limit;
  int complete_paths_at_cap;
  Path_segment_cache_T segment_cache;
  Mate_rescue_stats_T *rescue_stats;
} Mate_rescue_segment_policy_T;

static void
mate_rescue_segment_policy_init (Mate_rescue_segment_policy_T *policy,
				 Mate_rescue_stats_T *rescue_stats, const char *stage,
				 const char *side, const char *direction,
				 bool retain_local_mate_singlepaths_p) {
  memset(policy,0,sizeof(*policy));
  policy->activep = retain_local_mate_singlepaths_p == false ? true : false;
  policy->stage = stage;
  policy->side = side;
  policy->direction = direction;
  policy->dedup_chain_trigger = MATE_RESCUE_SEGMENT_DEDUP_CHAIN_TRIGGER;
  policy->dedup_complete_trigger = MATE_RESCUE_SEGMENT_DEDUP_COMPLETE_TRIGGER;
  policy->soft_chain_limit = MATE_RESCUE_SEGMENT_SOFT_CHAIN_LIMIT;
  policy->hard_chain_limit = MATE_RESCUE_SEGMENT_HARD_CHAIN_LIMIT;
  policy->soft_complete_path_limit = MATE_RESCUE_SEGMENT_SOFT_COMPLETE_LIMIT;
  policy->hard_complete_path_limit = MATE_RESCUE_SEGMENT_HARD_COMPLETE_LIMIT;
  policy->rescue_stats = rescue_stats;
  if (policy->activep == true) {
    MATE_RESCUE_STATS_INC(rescue_stats,segment_control_enabled);
  }
  return;
}

static int
mate_rescue_segment_complete_paths (Pathstore_T pathstore) {
  return List_length(pathstore->complete_sense_paths) +
    List_length(pathstore->complete_antisense_paths);
}

static void
mate_rescue_segment_policy_activate_dedup (Mate_rescue_segment_policy_T *policy) {
  if (policy->segment_cache == (Path_segment_cache_T) NULL) {
    policy->segment_cache = Path_segment_cache_new(/*emitted_path_dedup_p*/true);
  }
  if (policy->dedup_activated_p == false) {
    policy->dedup_activated_p = true;
    MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_dedup_activated);
  }
  return;
}

static void
mate_rescue_segment_policy_cap (Mate_rescue_segment_policy_T *policy,
				int complete_paths) {
  if (policy->cap_hit_p == false) {
    policy->cap_hit_p = true;
    policy->stop_solving_p = true;
    policy->complete_paths_at_cap = complete_paths;
    MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_cap_hit);
  }
  return;
}

static bool
mate_rescue_segment_policy_stop_p (Mate_rescue_segment_policy_T *policy,
				   Pathstore_T pathstore) {
  int complete_paths;

  if (policy->stop_solving_p == true) {
    return true;
  } else if (policy->chains_seen <= policy->soft_chain_limit) {
    return false;
  } else {
    complete_paths = mate_rescue_segment_complete_paths(pathstore);
    if (complete_paths >= policy->hard_complete_path_limit ||
	policy->chains_seen > policy->hard_chain_limit ||
	complete_paths >= policy->soft_complete_path_limit) {
      mate_rescue_segment_policy_cap(policy,complete_paths);
      return true;
    } else {
      return false;
    }
  }
}

static void
mate_rescue_segment_policy_free (Mate_rescue_segment_policy_T *policy,
				 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				 Listpool_T listpool, Pathpool_T pathpool,
				 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  if (policy != NULL && policy->segment_cache != (Path_segment_cache_T) NULL) {
    Path_segment_cache_free(&policy->segment_cache,intlistpool,univcoordlistpool,
			    listpool,pathpool,transcriptpool,hitlistpool);
  }
  return;
}

static Method_T
paired_gen_method_goal (void) {
  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_SEGMENT)) {
    return SEGMENT1;
  } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT)) {
    return EXT;
  } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXACT)) {
    return KMER_EXACT1;
  } else {
    return METHOD_INIT;
  }
}


static bool
paired_gen_method_available_p (Method_T last_method, Method_T method_goal) {
  return
    (last_method < KMER_EXACT1 && method_goal >= KMER_EXACT1 &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXACT)) ||
    (last_method < EXT && method_goal >= EXT &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT)) ||
    (last_method < SEGMENT1 && method_goal >= SEGMENT1 &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_SEGMENT));
}


static bool
tr_only_extension_still_available_p (Method_T last_method_5, Method_T last_method_3) {
  Method_T method_goal;

  if (paired_tr_adjacent_indel_repair_p == false) {
    return false;
  } else if (paired_gen_method_goal() != METHOD_INIT) {
    return false;
  } else if ((method_goal = single_read_tr_method_goal()) == METHOD_INIT) {
    return false;
  } else {
    return single_read_tr_method_available_p(last_method_5,method_goal) ||
      single_read_tr_method_available_p(last_method_3,method_goal);
  }
}

static bool
paired_transcript_completion_ready_p (int found_score_5, int found_score_3,
				      int end_score_5, int end_score_3,
				      Method_T last_method_5, Method_T last_method_3) {
  bool complete_mate3_p, complete_mate5_p;

  (void) last_method_5;
  (void) last_method_3;

  if (paired_transcript_completion_p == false) {
    return false;
  } else {
    complete_mate3_p = found_score_3 > found_score_5 ||
      (found_score_3 == found_score_5 && end_score_3 > end_score_5);
    complete_mate5_p = found_score_5 > found_score_3 ||
      (found_score_5 == found_score_3 && end_score_5 > end_score_3);

    return complete_mate3_p == true || complete_mate5_p == true;
  }
}


static int
paired_segment_exact_max_nmismatches (int querylength) {
  return querylength / paired_segment_exact_max_mismatch_divisor;
}


static int
positions_lower_bound (UINT4 *positions, int npositions, Trcoord_T goal) {
  int low = 0, high = npositions, mid;

  while (low < high) {
    mid = low + (high - low)/2;
    if ((Trcoord_T) positions[mid] < goal) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  return low;
}


static bool
trpath_list_has_trnum_p (List_T trpaths, Trnum_T trnum) {
  List_T p;
  Trpath_T trpath;

  for (p = trpaths; p != NULL; p = List_next(p)) {
    trpath = (Trpath_T) List_head(p);
    if (trpath->trnum == trnum) {
      return true;
    }
  }

  return false;
}

static bool
trpath_list_has_complete_trnum_p (List_T trpaths, Trnum_T trnum, int querylength) {
  List_T p;
  Trpath_T trpath;

  for (p = trpaths; p != NULL; p = List_next(p)) {
    trpath = (Trpath_T) List_head(p);
    if (trpath->trnum == trnum &&
	Intlist_head(trpath->endpoints) <= index1part_tr &&
	Intlist_last_value(trpath->endpoints) >= querylength - index1part_tr) {
      return true;
    }
  }

  return false;
}


static bool
trdiagonal_seen_p (Trcoord_T *tried_trdiagonals, int ntried, Trcoord_T trdiagonal) {
  int i;

  for (i = 0; i < ntried; i++) {
    if (tried_trdiagonals[i] == trdiagonal) {
      return true;
    }
  }

  return false;
}


static List_T
seeded_transcript_completion_for_mate (int *found_score, int *end_score, int sufficient_score,
				       List_T existing_trpaths, List_T anchor_trpaths,

				       T mate, int querylength, int *mismatch_positions_alloc,
				       Compress_T query_compress_fwd, Compress_T query_compress_rev,
				       bool tplusp, bool first_read_p,

				       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				       Listpool_T listpool, Trpathpool_T trpathpool,
				       Pathpool_T pathpool, Hitlistpool_T hitlistpool) {
  List_T new_trpaths = NULL, p;
  Trpath_T anchor, trpath;
  bool partialp;

  UINT4 *positions;
  int npositions, query_lastpos, querypos, diagterm, lowi, i;
  Trcoord_T position, trdiagonal, *tried_trdiagonals;
  int ntried = 0;

  query_lastpos = querylength - index1part_tr;
  if (query_lastpos < 0 || maxpaths_search <= 0) {
    return (List_T) NULL;
  } else {
    tried_trdiagonals = (Trcoord_T *) MALLOC(maxpaths_search * sizeof(Trcoord_T));
  }

  for (p = anchor_trpaths; p != NULL; p = List_next(p)) {
    anchor = (Trpath_T) List_head(p);

    if (trpath_list_has_complete_trnum_p(existing_trpaths,anchor->trnum,querylength) == true ||
	trpath_list_has_trnum_p(new_trpaths,anchor->trnum) == true) {
      /* Already have a complete mate candidate on this anchored transcript. */

    } else {
      for (querypos = 0; querypos <= query_lastpos && ntried < maxpaths_search; querypos++) {
	if (mate->tr_validp[querypos] == false) {
	  /* Skip */

	} else if (tplusp == true) {
	  positions = mate->tr_plus_positions[querypos];
	  npositions = mate->tr_plus_npositions[querypos];
	  diagterm = querylength - querypos;

	  lowi = positions_lower_bound(positions,npositions,anchor->troffset);
	  for (i = lowi; i < npositions && (position = (Trcoord_T) positions[i]) < anchor->trhigh &&
		 ntried < maxpaths_search; i++) {
	    trdiagonal = position + (Trcoord_T) diagterm;
	    if (trdiagonal_seen_p(tried_trdiagonals,ntried,trdiagonal) == false) {
	      tried_trdiagonals[ntried++] = trdiagonal;
	      if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&(*end_score),&partialp,
							 sufficient_score,trdiagonal,querypos,
							 querypos + index1part_tr,anchor->trnum,
							 anchor->troffset,anchor->trhigh,
							 query_compress_fwd,/*tplusp*/true,querylength,
							 mismatch_positions_alloc,anchor->chrnum,
							 anchor->transcript_genestrand,
							 /*paired_end_p*/true,first_read_p,mate->indelinfo,
							 intlistpool,uintlistpool,listpool,
							 trpathpool,pathpool,/*method*/TR_COMPLETE)) != NULL) {
		if (partialp == false) {
		  new_trpaths = Hitlist_push(new_trpaths,hitlistpool,(void *) trpath
					     hitlistpool_trace(__FILE__,__LINE__));
		  FREE(tried_trdiagonals);
		  return new_trpaths;
		} else {
		  Trpath_free(&trpath,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
		}
	      }
	    }
	  }

	} else {
	  positions = mate->tr_minus_positions[querypos];
	  npositions = mate->tr_minus_npositions[querypos];
	  diagterm = querypos + index1part_tr;

	  lowi = positions_lower_bound(positions,npositions,anchor->troffset);
	  for (i = lowi; i < npositions && (position = (Trcoord_T) positions[i]) < anchor->trhigh &&
		 ntried < maxpaths_search; i++) {
	    trdiagonal = position + (Trcoord_T) diagterm;
	    if (trdiagonal_seen_p(tried_trdiagonals,ntried,trdiagonal) == false) {
	      tried_trdiagonals[ntried++] = trdiagonal;
	      if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&(*end_score),&partialp,
							 sufficient_score,trdiagonal,querypos,
							 querypos + index1part_tr,anchor->trnum,
							 anchor->troffset,anchor->trhigh,
							 query_compress_rev,/*tplusp*/false,querylength,
							 mismatch_positions_alloc,anchor->chrnum,
							 anchor->transcript_genestrand,
							 /*paired_end_p*/true,first_read_p,mate->indelinfo,
							 intlistpool,uintlistpool,listpool,
							 trpathpool,pathpool,/*method*/TR_COMPLETE)) != NULL) {
		if (partialp == false) {
		  new_trpaths = Hitlist_push(new_trpaths,hitlistpool,(void *) trpath
					     hitlistpool_trace(__FILE__,__LINE__));
		  FREE(tried_trdiagonals);
		  return new_trpaths;
		} else {
		  Trpath_free(&trpath,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  FREE(tried_trdiagonals);

  return new_trpaths;
}


#if 0
/* Do not need to extend all paths, since fusion procedures will extend fusion candidates */
static void
eval_unextended_paths (int *found_score, int *end_score, T this,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  List_T p;

  for (p = this->unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),&(*end_score),(Path_T) List_head(p),
		       query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),&(*end_score),(Path_T) List_head(p),
		       query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),&(*end_score),(Path_T) List_head(p),
		       query_compress_fwd,query_compress_rev);
  }
  for (p = this->unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    Path_eval_nmatches(&(*found_score),&(*end_score),(Path_T) List_head(p),
		       query_compress_fwd,query_compress_rev);
  }

  return;
}
#endif


static List_T
paired_search_trdiagonals (int *found_score_paired, int *found_score_5, int *found_score_3,
			   int *end_score_paired, int *end_score_5, int *end_score_3,
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

			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			   unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			   int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

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

#if defined(DEBUG1) || defined(DEBUG_PAIRED_DEEP)
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
    if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
    if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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
    if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
    if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),
				 
			     pathpairs,/*new*/this5->sense_trpaths,/*new*/this3->sense_trpaths,
			     /*trpaths5*/NULL,/*trpaths3*/NULL,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
			     
			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD,paired_tr_adjacent_indel_repair_p);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),
				 
			     pathpairs,/*new*/this5->antisense_trpaths,/*new*/this3->antisense_trpaths,
			     /*trpaths5*/NULL,/*trpaths3*/NULL,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI,paired_tr_adjacent_indel_repair_p);

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
	if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
      if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
	if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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
      if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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
	if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
      if ((trpath5 = Trpath_solve_from_trdiagonal(&(*found_score_5),&(*end_score_5),&partialp,sufficient_score_5,
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
	if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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
      if ((trpath3 = Trpath_solve_from_trdiagonal(&(*found_score_3),&(*end_score_3),&partialp,sufficient_score_3,
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
			      int *end_score_paired, int *end_score_5, int *end_score_3,
			      int sufficient_score_5, Method_T *last_method_5,

			      List_T pathpairs, T this5, T this3, Knownsplicing_T knownsplicing,
			      
			      Shortread_T queryseq5, Shortread_T queryseq3,
			      int querylength5, int querylength3,

			      int *mismatch_positions_alloc_5, int nmismatches_allowed_5,

			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

			      int nmismatches_filter_5, int nmismatches_filter_3,
			      int mincoverage_filter_5, int mincoverage_filter_3,

			      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			      int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

			      int genestrand, Trdiagpool_T trdiagpool, Intlistpool_T intlistpool,
			      Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Transcriptpool_T transcriptpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T partial_sense_trpaths5, partial_antisense_trpaths5, sense_trpaths5, antisense_trpaths5;

  *last_method_5 = single_read_next_method_tr(&(*found_score_5),&(*end_score_5),sufficient_score_5,*last_method_5,
					     
					      &partial_sense_trpaths5,&partial_antisense_trpaths5,
					      &sense_trpaths5,&antisense_trpaths5,
					     
					      this5,genestrand,querylength5,
					      mismatch_positions_alloc_5,
					      Shortread_queryuc_ptr(queryseq5),Shortread_queryrc(queryseq5),
					      query5_compress_fwd,query5_compress_rev,
					      nmismatches_allowed_5,
					      trdiagpool,intlistpool,uintlistpool,
					      listpool,trpathpool,pathpool,hitlistpool,
					      /*paired_end_p*/true,/*first_read_p*/true,/*appendp*/false);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),
				 
			     pathpairs,/*new*/sense_trpaths5,/*sense_trpaths3*/NULL,
			     /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD,paired_tr_adjacent_indel_repair_p);
  this5->partial_sense_trpaths = List_append(partial_sense_trpaths5,this5->partial_sense_trpaths);
  this5->sense_trpaths = List_append(sense_trpaths5,this5->sense_trpaths);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),
				 
			     pathpairs,/*new*/antisense_trpaths5,/*antisense_trpaths3*/NULL,
			     /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI,paired_tr_adjacent_indel_repair_p);
  this5->partial_antisense_trpaths = List_append(partial_antisense_trpaths5,this5->partial_antisense_trpaths);
  this5->antisense_trpaths = List_append(antisense_trpaths5,this5->antisense_trpaths);

  return pathpairs;
}


static List_T
paired_read_next_method_tr_3 (int *found_score_paired, int *found_score_5, int *found_score_3,
			      int *end_score_paired, int *end_score_5, int *end_score_3,
			      int sufficient_score_3, Method_T *last_method_3,

			      List_T pathpairs, T this5, T this3, Knownsplicing_T knownsplicing,
			      
			      Shortread_T queryseq5, Shortread_T queryseq3,
			      int querylength5, int querylength3,

			      int *mismatch_positions_alloc_3, int nmismatches_allowed_3,
			      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			    
			      int nmismatches_filter_5, int nmismatches_filter_3,
			      int mincoverage_filter_5, int mincoverage_filter_3,

			      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			      int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

			      int genestrand, Trdiagpool_T trdiagpool, Intlistpool_T intlistpool,
			      Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Transcriptpool_T transcriptpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T partial_sense_trpaths3, partial_antisense_trpaths3, sense_trpaths3, antisense_trpaths3;

  *last_method_3 = single_read_next_method_tr(&(*found_score_3),&(*end_score_3),sufficient_score_3,*last_method_3,
					     
					      &partial_sense_trpaths3,&partial_antisense_trpaths3,
					      &sense_trpaths3,&antisense_trpaths3,
					     
					      this3,genestrand,querylength3,
					      mismatch_positions_alloc_3,
					      Shortread_queryuc_ptr(queryseq3),Shortread_queryrc(queryseq3),
					      query3_compress_fwd,query3_compress_rev,
					      nmismatches_allowed_3,
					      trdiagpool,intlistpool,uintlistpool,
					      listpool,trpathpool,pathpool,hitlistpool,
					      /*paired_end_p*/true,/*first_read_p*/false,/*appendp*/false);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),				 

			     pathpairs,/*sense_trpaths5*/NULL,/*new*/sense_trpaths3,
			     /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_FORWARD,paired_tr_adjacent_indel_repair_p);
  this3->partial_sense_trpaths = List_append(partial_sense_trpaths3,this3->partial_sense_trpaths);
  this3->sense_trpaths = List_append(sense_trpaths3,this3->sense_trpaths);

  pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
			     &(*end_score_paired),&(*end_score_5),&(*end_score_3),				 
				 
			     pathpairs,/*antisense_trpaths5*/NULL,/*new*/antisense_trpaths3,
			     /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,
			     
			     query5_compress_fwd,query5_compress_rev,
			     query3_compress_fwd,query3_compress_rev,
			     queryseq5,queryseq3,
			     querylength5,querylength3,this5,this3,knownsplicing,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*sensedir*/SENSE_ANTI,paired_tr_adjacent_indel_repair_p);
  this3->partial_antisense_trpaths = List_append(partial_antisense_trpaths3,this3->partial_antisense_trpaths);
  this3->antisense_trpaths = List_append(antisense_trpaths3,this3->antisense_trpaths);

  return pathpairs;
}


/* set_best_paths_p is set to false when computing from univdiagonals,
   and scanning over entire genome, but true when computing for mates */
/* Puts Path_T objects into pathstore */
static bool
path_satisfies_coverage_p (Path_T path, int coverage_threshold) {
  return Path_coverage(path) >= coverage_threshold ? true : false;
}

static bool
pathlist_all_satisfy_coverage_p (List_T paths, int coverage_threshold) {
  List_T p;

  for (p = paths; p != NULL; p = List_next(p)) {
    if (path_satisfies_coverage_p((Path_T) List_head(p),coverage_threshold) == false) {
      return false;
    }
  }

  return true;
}

static void
mark_pathstore_solved_after_auxinfo (bool *complete_sense_p, bool *complete_antisense_p,
				     Pathstore_T pathstore, int ncomplete_sense_before,
				     int ncomplete_antisense_before, int coverage_threshold) {
  bool new_complete_p = false;

  if (pathstore->complete_sense_paths != NULL) {
    *complete_sense_p = true;
    if (List_length(pathstore->complete_sense_paths) > ncomplete_sense_before) {
      new_complete_p = true;
    }
  }
  if (pathstore->complete_antisense_paths != NULL) {
    *complete_antisense_p = true;
    if (List_length(pathstore->complete_antisense_paths) > ncomplete_antisense_before) {
      new_complete_p = true;
    }
  }

  if (new_complete_p == true ||
      ((pathstore->complete_sense_paths != NULL || pathstore->complete_antisense_paths != NULL) &&
       pathlist_all_satisfy_coverage_p(pathstore->complete_sense_paths,coverage_threshold) == true &&
       pathlist_all_satisfy_coverage_p(pathstore->complete_antisense_paths,coverage_threshold) == true)) {
    pathstore->solvedp = true;
  }

  return;
}

struct Paired_segment_chain_apply_T {
  int *found_score;
  int *end_score;
  Pathstore_T pathstore;

  Shortread_T queryseq;
  char *queryptr;
  int querylength;
  int *mismatch_positions_alloc;
  Univdiag_T *novel_univdiags_alloc;
  Univcoord_T *novel_diagonals_alloc;
  unsigned short *localdb_alloc;
  int16_t *localdb_tally;
  int16_t *localdb_overflowp;
  int max_nblocks;
  int localdb_nmismatches_allowed;

  T stage1;
  Knownsplicing_T knownsplicing;
  Knownindels_T knownindels;
  Compress_T query_compress;
  Compress_T query_compress_fwd;
  Compress_T query_compress_rev;
  bool plusp;
  int genestrand;
  bool first_read_p;
  bool lowp;
  int selected_element_id;
  Mate_rescue_segment_policy_T *mate_rescue_segment_policy;
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
  Path_ext_anchor_cache_T anchor_cache;
  bool anchor_head_valid_p;
  Univcoord_T anchor_head_univdiagonal;
  int anchor_head_qstart;
  int anchor_head_qend;
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long profile_chain_count;
#endif

  Intlistpool_T intlistpool;
  Uintlistpool_T uintlistpool;
  Univcoordlistpool_T univcoordlistpool;
  Listpool_T listpool;
  Pathpool_T pathpool;
  Transcriptpool_T transcriptpool;
  Univdiagpool_T univdiagpool;
  Vectorpool_T vectorpool;
  Hitlistpool_T hitlistpool;
  Spliceendsgen_T spliceendsgen;

#ifdef GSNAP_PAIRED_SOLVE_STATS
#define PAIRED_CHAIN_HASH_SEEN_SIZE (1U << 20)
#define PAIRED_CHAIN_HASH_SEEN_MASK (PAIRED_CHAIN_HASH_SEEN_SIZE - 1U)
  unsigned long long chain_count;
  unsigned long long chain_hash;
  unsigned long long *chain_hash_seen;
  unsigned long long chain_hash_unique;
  unsigned long long chain_hash_repeats;
  unsigned long long chain_hash_table_full;
#endif
};

#if defined(GSNAP_PAIRED_SOLVE_STATS) || defined(GSNAP_ACCEL_PROFILE_STATS) || \
  defined(GSNAP_EXT_CHAIN_PREFIX_STATS) || defined(GSNAP_EXT_EXACT_CHAIN_DEDUP)
static unsigned long long
paired_segment_chain_hash_mix (unsigned long long hash, unsigned long long value) {
  int i;

  for (i = 0; i < 8; i++) {
    hash ^= value & 0xFFU;
    hash *= 1099511628211ULL;
    value >>= 8;
  }

  return hash;
}

static unsigned long long
paired_segment_chain_signature (List_T chain) {
  Univdiag_T univdiag;
  List_T p;
  unsigned long long hash = 1469598103934665603ULL;

  hash = paired_segment_chain_hash_mix(hash,0x434841494EULL);
  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) univdiag->univdiagonal);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) (unsigned int) univdiag->qstart);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) (unsigned int) univdiag->qend);
  }
  return paired_segment_chain_hash_mix(hash,0x454E44ULL);
}

#ifdef GSNAP_EXT_EXACT_CHAIN_DEDUP
static bool
ext_chain_equal_p (List_T x, List_T y) {
  Univdiag_T a, b;

  while (x != NULL && y != NULL) {
    a = (Univdiag_T) List_head(x);
    b = (Univdiag_T) List_head(y);
    /* solve_chain_qstart/qend consume only these three Univdiag fields.
       nmismatches is deliberately excluded: both solvers create their
       initial Path_T with nmismatches=-1 and never read the chain value. */
    if (a->univdiagonal != b->univdiagonal || a->qstart != b->qstart ||
	a->qend != b->qend) {
      return false;
    }
    x = List_next(x);
    y = List_next(y);
  }
  return x == NULL && y == NULL;
}
#endif
#endif

#ifdef GSNAP_EXT_CHAIN_PREFIX_STATS
static void
ext_chain_prefix_account (List_T chain, char *queryptr, bool plusp,
			  bool first_read_p, bool lowp) {
  Univdiag_T univdiag;
  List_T p;
  unsigned long long hash = 1469598103934665603ULL;
  unsigned int i, index;

  if (ext_prefix_seen == NULL) {
    ext_prefix_seen = (unsigned long long *) CALLOC(EXT_PREFIX_SEEN_SIZE,
						    sizeof(unsigned long long));
    if (ext_prefix_report_registered_p == false) {
      atexit(ext_prefix_report);
      ext_prefix_report_registered_p = true;
    }
  }
  hash = paired_segment_chain_hash_mix(hash,(unsigned long long) (uintptr_t) queryptr);
  hash = paired_segment_chain_hash_mix(hash,(unsigned long long) plusp);
  hash = paired_segment_chain_hash_mix(hash,(unsigned long long) first_read_p);
  hash = paired_segment_chain_hash_mix(hash,(unsigned long long) lowp);
  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) univdiag->univdiagonal);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) (unsigned int) univdiag->qstart);
    hash = paired_segment_chain_hash_mix(hash,(unsigned long long) (unsigned int) univdiag->qend);
    if (hash == 0ULL) hash = 1ULL;
    ext_prefix_nodes++;
    index = (unsigned int) hash & EXT_PREFIX_SEEN_MASK;
    for (i = 0; i < EXT_PREFIX_SEEN_SIZE; i++) {
      if (ext_prefix_seen[index] == 0ULL) {
	ext_prefix_seen[index] = hash;
	ext_prefix_unique++;
	break;
      } else if (ext_prefix_seen[index] == hash) {
	ext_prefix_repeats++;
	break;
      }
      index = (index + 1U) & EXT_PREFIX_SEEN_MASK;
    }
    if (i == EXT_PREFIX_SEEN_SIZE) ext_prefix_table_full++;
  }
}
#endif

#ifdef GSNAP_PAIRED_SOLVE_STATS
static void
paired_segment_chain_account (struct Paired_segment_chain_apply_T *data, List_T chain) {
  Univdiag_T univdiag;
  List_T p;
  unsigned long long chain_hash = 1469598103934665603ULL;
  unsigned int i, index;

  data->chain_hash = paired_segment_chain_hash_mix(data->chain_hash,0x434841494EULL);
  chain_hash = paired_segment_chain_hash_mix(chain_hash,0x434841494EULL);
  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    data->chain_hash = paired_segment_chain_hash_mix(data->chain_hash,
						     (unsigned long long) univdiag->univdiagonal);
    chain_hash = paired_segment_chain_hash_mix(chain_hash,
					       (unsigned long long) univdiag->univdiagonal);
    data->chain_hash = paired_segment_chain_hash_mix(data->chain_hash,
						     (unsigned long long) (unsigned int) univdiag->qstart);
    chain_hash = paired_segment_chain_hash_mix(chain_hash,
					       (unsigned long long) (unsigned int) univdiag->qstart);
    data->chain_hash = paired_segment_chain_hash_mix(data->chain_hash,
						     (unsigned long long) (unsigned int) univdiag->qend);
    chain_hash = paired_segment_chain_hash_mix(chain_hash,
					       (unsigned long long) (unsigned int) univdiag->qend);
  }
  data->chain_hash = paired_segment_chain_hash_mix(data->chain_hash,0x454E44ULL);
  chain_hash = paired_segment_chain_hash_mix(chain_hash,0x454E44ULL);
  data->chain_count++;

  if (chain_hash == 0) {
    chain_hash = 1;
  }
  index = (unsigned int) chain_hash & PAIRED_CHAIN_HASH_SEEN_MASK;
  for (i = 0; i < PAIRED_CHAIN_HASH_SEEN_SIZE; i++) {
    if (data->chain_hash_seen[index] == 0) {
      data->chain_hash_seen[index] = chain_hash;
      data->chain_hash_unique++;
      return;
    } else if (data->chain_hash_seen[index] == chain_hash) {
      data->chain_hash_repeats++;
      return;
    } else {
      index = (index + 1U) & PAIRED_CHAIN_HASH_SEEN_MASK;
    }
  }
  data->chain_hash_table_full++;
  return;
}

static void
paired_segment_pool_checkpoint (struct Paired_segment_chain_apply_T *data, const char *phase) {
  int univdiag_chunks, univdiag_list_chunks;
  int path_chunks, junction_chunks, altsplice_chunks, string_chunks;

  Univdiagpool_nchunks(&univdiag_chunks,&univdiag_list_chunks,data->univdiagpool);
  Pathpool_nchunks(&path_chunks,&junction_chunks,&altsplice_chunks,&string_chunks,data->pathpool);
  fprintf(stderr,
	  "paired_solve_stats\tsegment_progress\tphase=%s\tchain=%llu\tchain_hash_unique=%llu\t"
	  "chain_hash_repeats=%llu\tchain_hash_table_full=%llu\t"
	  "unextended_sense=%d\tunextended_antisense=%d\tcomplete_sense=%d\tcomplete_antisense=%d\t"
	  "list=%d\tunivdiag=%d\tunivdiag_list=%d\tpath=%d\tjunction=%d\taltsplice=%d\t"
	  "string=%d\thitlist=%d\n",
	  phase,data->chain_count,data->chain_hash_unique,data->chain_hash_repeats,
	  data->chain_hash_table_full,
	  List_length(data->pathstore->unextended_sense_paths),
	  List_length(data->pathstore->unextended_antisense_paths),
	  List_length(data->pathstore->complete_sense_paths),
	  List_length(data->pathstore->complete_antisense_paths),
	  Listpool_nchunks(data->listpool),univdiag_chunks,univdiag_list_chunks,
	  path_chunks,junction_chunks,altsplice_chunks,string_chunks,
	  Hitlistpool_nchunks(data->hitlistpool));
  return;
}
#endif

static void
paired_segment_chain_apply (List_T chain, void *x) {
  struct Paired_segment_chain_apply_T *data = (struct Paired_segment_chain_apply_T *) x;
  Pathstore_T pathstore = data->pathstore;
  Mate_rescue_segment_policy_T *policy = data->mate_rescue_segment_policy;
  Path_segment_cache_T segment_cache = (Path_segment_cache_T) NULL;
  bool emitted_path_dedup_p = false;
  int complete_paths;

  if (policy != (Mate_rescue_segment_policy_T *) NULL && policy->activep == true) {
    policy->chains_seen++;
    MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_chains_seen);

    if (mate_rescue_segment_policy_stop_p(policy,pathstore) == true) {
      policy->chains_skipped++;
      MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_chains_skipped);
      return;
    }
  }

#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
  Univdiag_T anchor_head = (Univdiag_T) List_head(chain);

  if (data->anchor_head_valid_p == true &&
      (anchor_head->univdiagonal != data->anchor_head_univdiagonal ||
       anchor_head->qstart != data->anchor_head_qstart ||
       anchor_head->qend != data->anchor_head_qend)) {
    Path_ext_anchor_cache_free(&data->anchor_cache,data->intlistpool,data->univcoordlistpool,
			       data->listpool,data->pathpool,data->transcriptpool,
			       data->hitlistpool);
    data->anchor_head_valid_p = false;
  }
  if (data->anchor_cache == NULL) {
    data->anchor_cache = Path_ext_anchor_cache_new();
    data->anchor_head_univdiagonal = anchor_head->univdiagonal;
    data->anchor_head_qstart = anchor_head->qstart;
    data->anchor_head_qend = anchor_head->qend;
    data->anchor_head_valid_p = true;
  }
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
  Univdiag_T profile_head = (Univdiag_T) List_head(chain);
#endif

#ifdef GSNAP_PAIRED_SOLVE_STATS
  paired_segment_chain_account(data,chain);
  if (data->chain_count == 1 || data->chain_count % 256 == 0) {
    paired_segment_pool_checkpoint(data,"before");
  }
#endif

  if (policy != (Mate_rescue_segment_policy_T *) NULL && policy->activep == true) {
    if (policy->segment_cache != (Path_segment_cache_T) NULL) {
      segment_cache = policy->segment_cache;
      emitted_path_dedup_p = true;
    } else if (policy->chains_seen >= policy->dedup_chain_trigger) {
      mate_rescue_segment_policy_activate_dedup(policy);
      segment_cache = policy->segment_cache;
      emitted_path_dedup_p = true;
    } else if (policy->chains_seen % MATE_RESCUE_SEGMENT_CHECK_INTERVAL == 0ULL) {
      complete_paths = mate_rescue_segment_complete_paths(pathstore);
      if (complete_paths >= policy->dedup_complete_trigger) {
	mate_rescue_segment_policy_activate_dedup(policy);
	segment_cache = policy->segment_cache;
	emitted_path_dedup_p = true;
      }
      if (complete_paths >= policy->hard_complete_path_limit) {
	mate_rescue_segment_policy_cap(policy,complete_paths);
      }
    }

    if (mate_rescue_segment_policy_stop_p(policy,pathstore) == true) {
      policy->chains_skipped++;
      MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_chains_skipped);
      return;
    }
  }

#ifdef GSNAP_ACCEL_PROFILE_STATS
  data->profile_chain_count++;
  Accel_profile_chain_begin(data->selected_element_id,data->profile_chain_count,
			    paired_segment_chain_signature(chain),List_length(chain),
			    profile_head->univdiagonal,profile_head->qstart,profile_head->qend,
			    data->first_read_p,data->plusp,data->lowp,"SEGMENT1");
#endif

  if (data->lowp == true) {
    if (policy != (Mate_rescue_segment_policy_T *) NULL && policy->activep == true) {
      segment_cache = policy->segment_cache;
      emitted_path_dedup_p = segment_cache != (Path_segment_cache_T) NULL ? true : false;
    }
    Path_solve_from_left_chain_cached(data->found_score,data->end_score,

			       &pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
			       &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,

			       chain,data->queryseq,data->queryptr,data->querylength,
			       data->mismatch_positions_alloc,
			       data->novel_univdiags_alloc,data->novel_diagonals_alloc,data->localdb_alloc,
			       data->localdb_tally,data->localdb_overflowp,
			       data->max_nblocks,data->localdb_nmismatches_allowed,

			       data->stage1,data->knownsplicing,data->knownindels,
			       data->query_compress,data->query_compress_fwd,data->query_compress_rev,
			       pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
			       data->plusp,data->genestrand,
			       /*paired_end_p*/true,data->first_read_p,
			       data->intlistpool,data->uintlistpool,data->univcoordlistpool,data->listpool,
			       data->pathpool,data->transcriptpool,data->univdiagpool,data->vectorpool,
			       data->hitlistpool,data->spliceendsgen,
			       /*method*/SEGMENT1,/*find_splices_p*/true,
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
			       data->anchor_cache,
#else
			       (Path_ext_anchor_cache_T) NULL,
#endif
			       segment_cache,
			       /*segment_prefix_skip_p*/false,
			       /*segment_qstart_anchor_cache_p*/false,
			       /*path_qend_local_cache_p*/false,
			       emitted_path_dedup_p);

  } else {
    if (policy != (Mate_rescue_segment_policy_T *) NULL && policy->activep == true) {
      segment_cache = policy->segment_cache;
      emitted_path_dedup_p = segment_cache != (Path_segment_cache_T) NULL ? true : false;
    }
    Path_solve_from_right_chain_cached(data->found_score,data->end_score,

				&pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
				&pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,

				chain,data->queryseq,data->queryptr,data->querylength,
				data->mismatch_positions_alloc,
				data->novel_univdiags_alloc,data->novel_diagonals_alloc,data->localdb_alloc,
				data->localdb_tally,data->localdb_overflowp,
				data->max_nblocks,data->localdb_nmismatches_allowed,

				data->stage1,data->knownsplicing,data->knownindels,
				data->query_compress,data->query_compress_fwd,data->query_compress_rev,
				pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
				data->plusp,data->genestrand,
				/*paired_end_p*/true,data->first_read_p,
				data->intlistpool,data->uintlistpool,data->univcoordlistpool,data->listpool,
				data->pathpool,data->transcriptpool,data->univdiagpool,data->vectorpool,
				data->hitlistpool,data->spliceendsgen,
				/*method*/SEGMENT1,/*find_splices_p*/true,
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
				data->anchor_cache,
#else
				(Path_ext_anchor_cache_T) NULL,
#endif
				segment_cache,
				/*segment_prefix_skip_p*/false,
				/*segment_qstart_anchor_cache_p*/false,
				/*path_qend_local_cache_p*/false,
				emitted_path_dedup_p);
  }

  if (policy != (Mate_rescue_segment_policy_T *) NULL && policy->activep == true) {
    policy->chains_solved++;
    MATE_RESCUE_STATS_INC(policy->rescue_stats,segment_control_chains_solved);
    if (policy->chains_solved % MATE_RESCUE_SEGMENT_CHECK_INTERVAL == 0ULL) {
      complete_paths = mate_rescue_segment_complete_paths(pathstore);
      if (complete_paths >= policy->hard_complete_path_limit ||
	  (policy->chains_seen >= policy->soft_chain_limit &&
	   complete_paths >= policy->soft_complete_path_limit)) {
	mate_rescue_segment_policy_cap(policy,complete_paths);
      }
    }
  }

#ifdef GSNAP_ACCEL_PROFILE_STATS
  Accel_profile_chain_end();
#endif

#ifdef GSNAP_PAIRED_SOLVE_STATS
  if (data->chain_count == 1 || data->chain_count % 256 == 0) {
    paired_segment_pool_checkpoint(data,"after");
  }
#endif

  return;
}

static void
solve_univdiagonal_auxinfo (bool *complete_sense_p, bool *complete_antisense_p,
			    int *found_score, int *end_score,
			    Univcoord_T univdiagonal, Auxinfo_T *auxinfo, Pathstore_T pathstore,
		  
			    Shortread_T queryseq,
			    char *queryptr, char *queryuc_ptr, char *queryrc, int querylength, int max_nsplices,
			    T this, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			    int *mismatch_positions_alloc, 
			    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			    unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			    int max_nblocks, int localdb_nmismatches_allowed,

			    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    int genestrand,
			     
			    Auxinfopool_T auxinfopool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			    Spliceendsgen_T spliceendsgen, bool plusp, bool first_read_p, bool lowp,
			    bool set_best_paths_p, bool force_pending_auxinfo_p,
			    Mate_rescue_segment_policy_T *mate_rescue_segment_policy) {
  List_T chains, chain, p;
  Singleelt_T *eltarray;
  struct Paired_segment_chain_apply_T segment_apply;
  int nsplices, narray, ncomplete_sense_before, ncomplete_antisense_before;
  Path_ext_anchor_cache_T ext_anchor_cache = (Path_ext_anchor_cache_T) NULL;
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
  Univdiag_T ext_head;
  Univcoord_T ext_cache_head_univdiagonal = 0U;
  int ext_cache_head_qstart = 0, ext_cache_head_qend = 0;
  bool ext_cache_head_valid_p = false;
#endif
#ifdef GSNAP_EXT_EXACT_CHAIN_DEDUP
  List_T ext_previous_chain = (List_T) NULL;
  int ext_previous_chaindir = CHAIN_NONE;
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long profile_ext_chain_count = 0ULL;
#endif
  
  if (univdiagonal < (Univcoord_T) querylength) {
    /* Skip */

  } else {
    while ((pathstore->solvedp == false || force_pending_auxinfo_p == true) && *auxinfo != NULL) {
      debug(printf("Method is %s\n",Method_string((*auxinfo)->method)));
      ncomplete_sense_before = List_length(pathstore->complete_sense_paths);
      ncomplete_antisense_before = List_length(pathstore->complete_antisense_paths);

      if ((*auxinfo)->method == KMER_EXACT1) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif
	debug(printf("Solving auxinfo %p using method EXACT1\n",*auxinfo));
	Path_solve_exact(&(*found_score),&(*end_score),
			 
			 &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,
			 
			 univdiagonal,*auxinfo,querylength,
			 /*max_nmismatches*/-1,
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
	mark_pathstore_solved_after_auxinfo(&(*complete_sense_p),&(*complete_antisense_p),
					    pathstore,ncomplete_sense_before,ncomplete_antisense_before,
					    querylength);
	
      } else if ((*auxinfo)->method == EXT) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif

	if (lowp == true) {
	  assert((*auxinfo)->chaindir == CHAIN_LEFT);
	  for (nsplices = 0; nsplices <= (*auxinfo)->max_nsplices; nsplices++) {
	    for (p = (*auxinfo)->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
	      chain = (List_T) List_head(p);
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
	      ext_head = (Univdiag_T) List_head(chain);
	      if (ext_cache_head_valid_p == true &&
		  (ext_head->univdiagonal != ext_cache_head_univdiagonal ||
		   ext_head->qstart != ext_cache_head_qstart ||
		   ext_head->qend != ext_cache_head_qend)) {
		Path_ext_anchor_cache_free(&ext_anchor_cache,
		  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
		ext_cache_head_valid_p = false;
	      }
	      if (ext_anchor_cache == NULL &&
		  List_next(p) != NULL) {
		List_T next_chain = (List_T) List_head(List_next(p));
		Univdiag_T next_head = (Univdiag_T) List_head(next_chain);
		if (next_head->univdiagonal == ext_head->univdiagonal &&
		    next_head->qstart == ext_head->qstart && next_head->qend == ext_head->qend) {
		  ext_anchor_cache = Path_ext_anchor_cache_new();
		  ext_cache_head_univdiagonal = ext_head->univdiagonal;
		  ext_cache_head_qstart = ext_head->qstart;
		  ext_cache_head_qend = ext_head->qend;
		  ext_cache_head_valid_p = true;
		}
	      }
#endif
#ifdef GSNAP_EXT_EXACT_CHAIN_DEDUP
	      if (ext_previous_chaindir == CHAIN_LEFT &&
		  ext_chain_equal_p(chain,ext_previous_chain) == true) {
		continue;
	      }
	      ext_previous_chain = chain;
	      ext_previous_chaindir = CHAIN_LEFT;
#endif
#ifdef GSNAP_EXT_CHAIN_PREFIX_STATS
	      ext_chain_prefix_account(chain,queryptr,plusp,first_read_p,lowp);
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
	      profile_ext_chain_count++;
	      Accel_profile_chain_begin((*auxinfo)->elti,profile_ext_chain_count,
				  paired_segment_chain_signature(chain),List_length(chain),
				  ((Univdiag_T) List_head(chain))->univdiagonal,
				  ((Univdiag_T) List_head(chain))->qstart,
				  ((Univdiag_T) List_head(chain))->qend,
				  first_read_p,plusp,lowp,"EXT");
#endif

	      Path_solve_from_left_chain_cached(&(*found_score),&(*end_score),
				       
					 &pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
					 &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,

					 chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					 this,knownsplicing,knownindels,
					 query_compress,query_compress_fwd,query_compress_rev,
					 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
					 /*chrhigh*/pathstore->chrhigh,plusp,genestrand,
					 /*paired_end_p*/true,first_read_p,
					 intlistpool,uintlistpool,univcoordlistpool,listpool,
					 pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
					 (*auxinfo)->method,/*find_splices_p*/true,
					 ext_anchor_cache,
					 /*segment_cache*/(Path_segment_cache_T) NULL,
					 /*segment_prefix_skip_p*/false,
					 /*segment_qstart_anchor_cache_p*/false,
					 /*path_qend_local_cache_p*/false,
					 /*emitted_path_dedup_p*/false);
#ifdef GSNAP_ACCEL_PROFILE_STATS
	      Accel_profile_chain_end();
#endif
	    }
	  }

	} else {
	  assert((*auxinfo)->chaindir == CHAIN_RIGHT);
	  for (nsplices = 0; nsplices <= (*auxinfo)->max_nsplices; nsplices++) {
	    for (p = (*auxinfo)->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
	      chain = (List_T) List_head(p);
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
	      ext_head = (Univdiag_T) List_head(chain);
	      if (ext_cache_head_valid_p == true &&
		  (ext_head->univdiagonal != ext_cache_head_univdiagonal ||
		   ext_head->qstart != ext_cache_head_qstart ||
		   ext_head->qend != ext_cache_head_qend)) {
		Path_ext_anchor_cache_free(&ext_anchor_cache,
		  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
		ext_cache_head_valid_p = false;
	      }
	      if (ext_anchor_cache == NULL &&
		  List_next(p) != NULL) {
		List_T next_chain = (List_T) List_head(List_next(p));
		Univdiag_T next_head = (Univdiag_T) List_head(next_chain);
		if (next_head->univdiagonal == ext_head->univdiagonal &&
		    next_head->qstart == ext_head->qstart && next_head->qend == ext_head->qend) {
		  ext_anchor_cache = Path_ext_anchor_cache_new();
		  ext_cache_head_univdiagonal = ext_head->univdiagonal;
		  ext_cache_head_qstart = ext_head->qstart;
		  ext_cache_head_qend = ext_head->qend;
		  ext_cache_head_valid_p = true;
		}
	      }
#endif
#ifdef GSNAP_EXT_EXACT_CHAIN_DEDUP
	      if (ext_previous_chaindir == CHAIN_RIGHT &&
		  ext_chain_equal_p(chain,ext_previous_chain) == true) {
		continue;
	      }
	      ext_previous_chain = chain;
	      ext_previous_chaindir = CHAIN_RIGHT;
#endif
#ifdef GSNAP_EXT_CHAIN_PREFIX_STATS
	      ext_chain_prefix_account(chain,queryptr,plusp,first_read_p,lowp);
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
	      profile_ext_chain_count++;
	      Accel_profile_chain_begin((*auxinfo)->elti,profile_ext_chain_count,
				  paired_segment_chain_signature(chain),List_length(chain),
				  ((Univdiag_T) List_head(chain))->univdiagonal,
				  ((Univdiag_T) List_head(chain))->qstart,
				  ((Univdiag_T) List_head(chain))->qend,
				  first_read_p,plusp,lowp,"EXT");
#endif
	      
	      Path_solve_from_right_chain_cached(&(*found_score),&(*end_score),
				       
					  &pathstore->unextended_sense_paths,&pathstore->unextended_antisense_paths,
					  &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,

					  chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					  this,knownsplicing,knownindels,
					  query_compress,query_compress_fwd,query_compress_rev,
					  /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
					  /*chrhigh*/pathstore->chrhigh,plusp,genestrand,
					  /*paired_end_p*/true,first_read_p,
					  intlistpool,uintlistpool,univcoordlistpool,listpool,
					  pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
					  (*auxinfo)->method,/*find_splices_p*/true,
					  ext_anchor_cache,
					  /*segment_cache*/
					  (Path_segment_cache_T) NULL,
					  /*segment_prefix_skip_p*/false,
					  /*segment_qstart_anchor_cache_p*/false,
					  /*path_qend_local_cache_p*/
					  false,
					  /*emitted_path_dedup_p*/false);
#ifdef GSNAP_ACCEL_PROFILE_STATS
	      Accel_profile_chain_end();
#endif
	    }
	  }
	}

	if (set_best_paths_p == true) {
	  Pathstore_set_best_paths(pathstore,hitlistpool);
	}
	mark_pathstore_solved_after_auxinfo(&(*complete_sense_p),&(*complete_antisense_p),
					    pathstore,ncomplete_sense_before,ncomplete_antisense_before,
					    querylength);
	
      } else if ((*auxinfo)->method == SEGMENT1) {
#ifdef INDIVIDUAL_CHRINFO
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     univdiagonal - querylength,univdiagonal);
#endif

	/* Preserve low-mismatch full-read segment hits before chain solving
	   can prefer off-diagonal paths. */
	Path_solve_exact(&(*found_score),&(*end_score),

			 &pathstore->complete_sense_paths,&pathstore->complete_antisense_paths,

			 univdiagonal,*auxinfo,querylength,
			 /*max_nmismatches*/paired_segment_exact_max_nmismatches(querylength),
			 plusp,first_read_p,genestrand,
			 query_compress,query_compress_fwd,query_compress_rev,
			 queryseq,queryuc_ptr,queryrc,
			 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			 /*chrhigh*/pathstore->chrhigh,
			 intlistpool,uintlistpool,univcoordlistpool,
			 listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			 /*method*/SEGMENT1);

	if (plusp == true) {
	  eltarray = this->exhaustive_eltarray_gplus;
	  narray = this->nexhaustive_gplus;
	} else {
	  eltarray = this->exhaustive_eltarray_gminus;
	  narray = this->nexhaustive_gminus;
	}

	memset(&segment_apply,0,sizeof(segment_apply));
	segment_apply.found_score = &(*found_score);
	segment_apply.end_score = &(*end_score);
	segment_apply.pathstore = pathstore;
	segment_apply.queryseq = queryseq;
	segment_apply.queryptr = queryptr;
	segment_apply.querylength = querylength;
	segment_apply.mismatch_positions_alloc = mismatch_positions_alloc;
	segment_apply.novel_univdiags_alloc = novel_univdiags_alloc;
	segment_apply.novel_diagonals_alloc = novel_diagonals_alloc;
	segment_apply.localdb_alloc = localdb_alloc;
	segment_apply.localdb_tally = localdb_tally;
	segment_apply.localdb_overflowp = localdb_overflowp;
	segment_apply.max_nblocks = max_nblocks;
	segment_apply.localdb_nmismatches_allowed = localdb_nmismatches_allowed;
	segment_apply.stage1 = this;
	segment_apply.knownsplicing = knownsplicing;
	segment_apply.knownindels = knownindels;
	segment_apply.query_compress = query_compress;
	segment_apply.query_compress_fwd = query_compress_fwd;
	segment_apply.query_compress_rev = query_compress_rev;
	segment_apply.plusp = plusp;
	segment_apply.genestrand = genestrand;
	segment_apply.first_read_p = first_read_p;
	segment_apply.lowp = lowp;
	segment_apply.selected_element_id = (*auxinfo)->elti;
	segment_apply.mate_rescue_segment_policy = mate_rescue_segment_policy;
	segment_apply.intlistpool = intlistpool;
	segment_apply.uintlistpool = uintlistpool;
	segment_apply.univcoordlistpool = univcoordlistpool;
	segment_apply.listpool = listpool;
	segment_apply.pathpool = pathpool;
	segment_apply.transcriptpool = transcriptpool;
	segment_apply.univdiagpool = univdiagpool;
	segment_apply.vectorpool = vectorpool;
	segment_apply.hitlistpool = hitlistpool;
	segment_apply.spliceendsgen = spliceendsgen;
#ifdef GSNAP_PAIRED_SOLVE_STATS
	segment_apply.chain_hash = 1469598103934665603ULL;
	segment_apply.chain_hash_seen = (unsigned long long *)
	  CALLOC(PAIRED_CHAIN_HASH_SEEN_SIZE,sizeof(unsigned long long));
#endif

	Pathpool_set_path_reuse(pathpool,true);
	if (paired_segment_stream_chains_p == true) {
	  if (lowp == true) {
	    Kmer_apply_eltarray_left(eltarray,/*elti*/(*auxinfo)->elti,narray,max_nsplices,
				     listpool,/*skip_chain_p*/NULL,
				     /*equiv_prev_prune_p*/false,/*max_chain_window*/0,
				     /*materialized_order_p*/true,
				     paired_segment_chain_apply,(void *) &segment_apply);
	  } else {
	    Kmer_apply_eltarray_right(eltarray,/*elti*/(*auxinfo)->elti,narray,max_nsplices,
				      listpool,/*skip_chain_p*/NULL,
				      /*equiv_prev_prune_p*/false,/*max_chain_window*/0,
				      /*materialized_order_p*/true,
				      paired_segment_chain_apply,(void *) &segment_apply);
	  }

	} else if (lowp == true) {
	  chains = Kmer_solve_eltarray_left(eltarray,/*elti*/(*auxinfo)->elti,narray,
					    max_nsplices,univdiagpool,listpool);
	  for (p = chains; p != NULL; p = List_next(p)) {
	    paired_segment_chain_apply((List_T) List_head(p),(void *) &segment_apply);
	  }

	} else {
	  chains = Kmer_solve_eltarray_right(eltarray,/*elti*/(*auxinfo)->elti,narray,
					     max_nsplices,univdiagpool,listpool);
	  for (p = chains; p != NULL; p = List_next(p)) {
	    paired_segment_chain_apply((List_T) List_head(p),(void *) &segment_apply);
	  }
	}
#ifdef GSNAP_EXT_ANCHOR_CACHE_ACTIVE
	if (segment_apply.anchor_cache != (Path_ext_anchor_cache_T) NULL) {
	  Path_ext_anchor_cache_free(&segment_apply.anchor_cache,intlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,hitlistpool);
	}
#endif
	mate_rescue_segment_policy_free(mate_rescue_segment_policy,intlistpool,univcoordlistpool,
					listpool,pathpool,transcriptpool,hitlistpool);
	Pathpool_set_path_reuse(pathpool,false);

#ifdef GSNAP_PAIRED_SOLVE_STATS
	fprintf(stderr,
		"paired_solve_stats\tsegment_chains\tstream=%d\tplus=%d\tlow=%d\t"
		"elti=%d\tnarray=%d\tchains=%llu\thash=%016llx\t"
		"mate_rescue_control=%d\tcontrol_dedup=%d\tcontrol_cap=%d\t"
		"control_chains_seen=%llu\tcontrol_chains_solved=%llu\t"
		"control_chains_skipped=%llu\tcontrol_complete_paths_at_cap=%d\n",
		paired_segment_stream_chains_p,plusp,lowp,(*auxinfo)->elti,narray,
		segment_apply.chain_count,segment_apply.chain_hash,
		mate_rescue_segment_policy != NULL && mate_rescue_segment_policy->activep == true ? 1 : 0,
		mate_rescue_segment_policy != NULL && mate_rescue_segment_policy->dedup_activated_p == true ? 1 : 0,
		mate_rescue_segment_policy != NULL && mate_rescue_segment_policy->cap_hit_p == true ? 1 : 0,
		mate_rescue_segment_policy != NULL ? mate_rescue_segment_policy->chains_seen : 0ULL,
		mate_rescue_segment_policy != NULL ? mate_rescue_segment_policy->chains_solved : 0ULL,
		mate_rescue_segment_policy != NULL ? mate_rescue_segment_policy->chains_skipped : 0ULL,
		mate_rescue_segment_policy != NULL ? mate_rescue_segment_policy->complete_paths_at_cap : 0);
	FREE(segment_apply.chain_hash_seen);
#endif

	if (set_best_paths_p == true) {
	  Pathstore_set_best_paths(pathstore,hitlistpool);
	}
	mark_pathstore_solved_after_auxinfo(&(*complete_sense_p),&(*complete_antisense_p),
					    pathstore,ncomplete_sense_before,ncomplete_antisense_before,
					    querylength);

      } else {
	fprintf(stderr,"Unexpected method %s\n",Method_string((*auxinfo)->method));
	abort();
      }

      *auxinfo = Auxinfo_pop(*auxinfo,univdiagpool,auxinfopool,listpool);
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

  if (ext_anchor_cache != (Path_ext_anchor_cache_T) NULL) {
    Path_ext_anchor_cache_free(&ext_anchor_cache,
      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  }

  return;
}


#define GENERATE_ALL_PAIRS 1

static Pairtype_T
determine_pairtype (Path_T path5, Path_T path3);

static bool
make_pathpair_concordant_p (Path_T pathL, Path_T pathH, bool plusp) {
  if (plusp == true) {
    return determine_pairtype(pathL,pathH) == CONCORDANT ? true : false;
  } else {
    return determine_pairtype(pathH,pathL) == CONCORDANT ? true : false;
  }
}

static List_T
make_pathpairs (int *found_score_paired, int *end_score_paired,
		List_T *unresolved_pathpairs, List_T pathpairs,
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

      if (make_pathpair_concordant_p(pathL,pathH,plusp) == false) {
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
	if (Pathpair_end_score(pathpair) < *end_score_paired) {
	  *end_score_paired = Pathpair_end_score(pathpair);
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

      if (make_pathpair_concordant_p(pathL,pathH,plusp) == false) {
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
	if (Pathpair_end_score(pathpair) < *end_score_paired) {
	  *end_score_paired = Pathpair_end_score(pathpair);
	}
	pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  return pathpairs;
}

static List_T
replay_near_exact_terminal_anchor_pathpairs (List_T *unresolved_pathpairs, bool replay5p,
					     Path_T *partner_patharray, int partner_npaths,
					     List_T terminal_anchor_paths_deferred,
					     Shortread_T queryseq5, Shortread_T queryseq3,
					     int nmismatches_filter_5, int nmismatches_filter_3,
					     int mincoverage_filter_5, int mincoverage_filter_3,
					     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
					     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
					     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  Pathpair_T pathpair;
  Path_T path5, path3, partner;
  List_T p, pathpairs = NULL;
  int i;

  for (i = 0; i < partner_npaths; i++) {
    partner = partner_patharray[i];
    if (Pathpair_near_exact_terminal_anchor_partner_p(partner) == true) {
      for (p = terminal_anchor_paths_deferred; p != NULL; p = List_next(p)) {
	if (replay5p == true) {
	  path5 = (Path_T) List_head(p);
	  path3 = partner;
	} else {
	  path5 = partner;
	  path3 = (Path_T) List_head(p);
	}

	if (determine_pairtype(path5,path3) != CONCORDANT) {
	  /* Skip */
	} else if (path5->plusp == true) {
	  if ((pathpair = Pathpair_new_concordant(&(*unresolved_pathpairs),
						  /*pathL*/path5,/*pathH*/path3,
						  /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,
						  /*plusp*/true,
						  nmismatches_filter_5,nmismatches_filter_3,
						  mincoverage_filter_5,mincoverage_filter_3,
						  intlistpool,univcoordlistpool,listpool,
						  pathpool,vectorpool,transcriptpool,hitlistpool,
						  /*check_inner_p*/false,
						  /*copyLp*/true,/*copyHp*/true)) != NULL) {
	    pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				     hitlistpool_trace(__FILE__,__LINE__));
	  }
	} else {
	  if ((pathpair = Pathpair_new_concordant(&(*unresolved_pathpairs),
						  /*pathL*/path3,/*pathH*/path5,
						  /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,
						  /*plusp*/false,
						  nmismatches_filter_5,nmismatches_filter_3,
						  mincoverage_filter_5,mincoverage_filter_3,
						  intlistpool,univcoordlistpool,listpool,
						  pathpool,vectorpool,transcriptpool,hitlistpool,
						  /*check_inner_p*/false,
						  /*copyLp*/true,/*copyHp*/true)) != NULL) {
	    pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				     hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
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
gather_pathpairs_mutual (bool *completeL_p, bool *completeH_p,
			 int *found_score_paired, int *end_score_paired,
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
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
  
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
  
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
      pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
      pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
gather_pathpairs_optimal (bool *completeL_p, bool *completeH_p,
			  int *found_score_paired, int *end_score_paired,
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
  
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
  
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
  last_pathstoreL = NULL;	/* Not necessary because we check for i == 0 */
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
  last_pathstoreH = NULL;	/* Not necessary because we check for i == 0 */
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
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
  
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
  printf("returning %d pathpairs\n",List_length(pathpairs));
  for (List_T p = pathpairs; p != NULL; p = List_next(p)) {
    Pathpair_print((Pathpair_T) List_head(p));
  }
#endif

  return pathpairs;
}


static List_T
concordance_univdiagonals (bool *completeL_p, bool *completeH_p,
			   int *found_score_paired, int *found_score_L, int *found_score_H,
			   int *end_score_paired, int *end_score_L, int *end_score_H,
		  
			   List_T *unresolved_pathpairs, List_T pathpairs, T thisL, T thisH,
	      
			   Univcoord_T *univdiagonalsL_array, Auxinfo_T *auxinfoL_array,
			   Auxinfo_T **auxinfoL_ptr_array,
			   Pathstore_T *pathstoreL_array, int nunivdiagonalsL_array,
			   Univcoord_T *univdiagonalsH_array, Auxinfo_T *auxinfoH_array,
			   Auxinfo_T **auxinfoH_ptr_array,
			   Pathstore_T *pathstoreH_array, int nunivdiagonalsH_array,

			   Shortread_T queryseqL, Shortread_T queryseqH,
			   char *queryptrL, char *queryuc_ptr_L, char *queryrcL, int querylengthL, int max_nsplices_L,
			   char *queryptrH, char *queryuc_ptr_H, char *queryrcH, int querylengthH, int max_nsplices_H,
			   Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			   int *mismatch_positions_alloc_L, int *mismatch_positions_alloc_H,
			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			   int localdb_nmismatches_allowed,

			   Compress_T queryL_compress, Compress_T queryL_compress_fwd, Compress_T queryL_compress_rev,
			   Compress_T queryH_compress, Compress_T queryH_compress_fwd, Compress_T queryH_compress_rev,
			   int sufficient_score_L, int sufficient_score_H, int genestrand,
			     
			   int nmismatches_filter_5, int nmismatches_filter_3,
			   int mincoverage_filter_5, int mincoverage_filter_3,

			   Auxinfopool_T auxinfopool, Intlistpool_T intlistpool,
			   Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			   Spliceendsgen_T spliceendsgenL, Spliceendsgen_T spliceendsgenH,
			   bool plusp, bool L_first_read_p, bool H_first_read_p,
			   bool force_pending_auxinfo_L_p, bool force_pending_auxinfo_H_p) {

  int *indices, index1, index2;
  int nindices, i, k;
  Univcoord_T univdiagonalL, univdiagonalH;
  Pathstore_T pathstoreL, pathstoreH, last_pathstoreL, last_pathstoreH;
  Auxinfo_T *auxinfo_ptr;
  bool complete_sense_p, complete_antisense_p;
  /* List_T paths; */
  /* int cmp; */
#ifdef GSNAP_PAIRED_SOLVE_STATS
  double concordance_t0 = paired_solve_stats_now();
#endif


  debug(printf("Entering concordance_univdiagonals with found_scores %d and %d, array sizes %d and %d\n",
	       *found_score_L,*found_score_H,nunivdiagonalsL_array,nunivdiagonalsH_array));

#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
	auxinfo_ptr = (auxinfoL_ptr_array != NULL) ? auxinfoL_ptr_array[index1] : &(auxinfoL_array[index1]);
	solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				   &(*found_score_L),&(*end_score_L),
				   univdiagonalL,auxinfo_ptr,pathstoreL,
				   
				   queryseqL,queryptrL,queryuc_ptr_L,queryrcL,querylengthL,max_nsplices_L,
				   thisL,knownsplicing,knownindels,
				   
				   mismatch_positions_alloc_L,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
				   
				   queryL_compress,queryL_compress_fwd,queryL_compress_rev,genestrand,
				   auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
				   hitlistpool,spliceendsgenL,plusp,L_first_read_p,/*lowp*/true,
				   /*set_best_paths_p*/false,force_pending_auxinfo_L_p,
				   /*mate_rescue_segment_policy*/NULL);
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
	auxinfo_ptr = (auxinfoH_ptr_array != NULL) ? auxinfoH_ptr_array[index2] : &(auxinfoH_array[index2]);
	solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
				   &(*found_score_H),&(*end_score_H),
				   univdiagonalH,auxinfo_ptr,pathstoreH,
				   
				   queryseqH,queryptrH,queryuc_ptr_H,queryrcH,querylengthH,max_nsplices_H,
				   thisH,knownsplicing,knownindels,
				   
				   mismatch_positions_alloc_H,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
				   
				   queryH_compress,queryH_compress_fwd,queryH_compress_rev,genestrand,
				   auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
				   hitlistpool,spliceendsgenH,plusp,H_first_read_p,/*lowp*/false,
				   /*set_best_paths_p*/false,force_pending_auxinfo_H_p,
				   /*mate_rescue_segment_policy*/NULL);
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
    
#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
	gather_pathpairs_mutual(&(*completeL_p),&(*completeH_p),&(*found_score_paired),&(*end_score_paired),
				 &(*unresolved_pathpairs),pathpairs,indices,nindices,
					   
				 univdiagonalsL_array,pathstoreL_array,univdiagonalsH_array,pathstoreH_array,
				 
				 queryseqL,queryseqH,nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 
				 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,plusp);
      
    } else {
      pathpairs =
	gather_pathpairs_optimal(&(*completeL_p),&(*completeH_p),&(*found_score_paired),&(*end_score_paired),
				 &(*unresolved_pathpairs),pathpairs,indices,nindices,
			       
				 univdiagonalsL_array,pathstoreL_array,univdiagonalsH_array,pathstoreH_array,
				 
				 queryseqL,queryseqH,nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 
				 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,plusp);
    }

#ifdef GSNAP_PAIRED_SOLVE_STATS
    fprintf(stderr,
	    "paired_solve_stats\tconcordance\tplus=%d\tL_first=%d\tH_first=%d\tL_in=%d\tH_in=%d\tnindices=%d\tpathpairs=%d\tunresolved=%d\tfound_pair=%d\tfoundL=%d\tfoundH=%d\telapsed=%.6f\n",
	    plusp,L_first_read_p,H_first_read_p,nunivdiagonalsL_array,nunivdiagonalsH_array,
	    nindices,List_length(pathpairs),List_length(*unresolved_pathpairs),
	    *found_score_paired,*found_score_L,*found_score_H,paired_solve_stats_now() - concordance_t0);
#endif

    FREE(indices);
  }


  return pathpairs;
}


/* _univdiagonals_gplus and _univdiagonals_gminus are aligned */
static void
single_read_gen_univdiagonals (Method_T *last_method, int *querystart, int *queryend, int end_score,

			       Univcoord_T **_univdiagonals_gplus, Univcoord_T **_univdiagonals_gminus,
			       Auxinfo_T **auxinfo_gplus, Auxinfo_T **auxinfo_gminus,
			       Pathstore_T **pathstores_gplus, Pathstore_T **pathstores_gminus,
			       int *nunivdiagonals_gplus, int *nunivdiagonals_gminus,

			       T this, int genestrand,			       
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       int querylength, char *quality_string, const char *read_name,
			       int max_nsplices, EF64_T repetitive_ef64,
			       
			       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			       Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
			       Univcoordlistpool_T univcoordlistpool, Intlistpool_T intlistpool, Listpool_T listpool,
			       bool first_read_p, Method_T method_goal) {

  int total_npositions_plus, total_npositions_minus;
  int **indices_by_nsplices_gplus, **indices_by_nsplices_gminus;
  int *nindices_by_nsplices_gplus, *nindices_by_nsplices_gminus;
  int nsplices;


  debug(printf("Entering single_read_gen_univdiagonals with last method %s\n",Method_string(*last_method)));

  if (*last_method < KMER_EXACT1 && method_goal >= KMER_EXACT1 &&
      Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXACT)) {
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

  } else if (*last_method < EXT && method_goal >= EXT &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT)) {
    /* 2. Extension search */
    debug(printf("%s Read: 2.  Running Extension search\n",first_read_p ? "5'" : "3'"));
    Stage1_fill_all_oligos_gen(this,querylength,genestrand);
    Extension_search(&indices_by_nsplices_gplus,&nindices_by_nsplices_gplus,
		     &indices_by_nsplices_gminus,&nindices_by_nsplices_gminus,

		     &(*_univdiagonals_gplus),&(*auxinfo_gplus),
		     &(*pathstores_gplus),&(*nunivdiagonals_gplus),
		     &(*_univdiagonals_gminus),&(*auxinfo_gminus),
		     &(*pathstores_gminus),&(*nunivdiagonals_gminus),

			     this,query_compress_fwd,query_compress_rev,
			     querylength,max_nsplices,first_read_p,/*retain_terminal_p*/true,
			     /*prune_equal_score_p*/false,
			     segmenteltpool,auxinfopool,pathstorepool,
		     univcoordlistpool,univdiagpool,intlistpool,listpool);

    Pathstore_assign_chrinfo(*_univdiagonals_gplus,*pathstores_gplus,*nunivdiagonals_gplus,querylength);
    Pathstore_assign_chrinfo(*_univdiagonals_gminus,*pathstores_gminus,*nunivdiagonals_gminus,querylength);

    /* For paired-end, using concordance rather than trying all values of nsplices */
    for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
      FREE(indices_by_nsplices_gminus[nsplices]);
      FREE(indices_by_nsplices_gplus[nsplices]);
    }
    FREE(indices_by_nsplices_gminus);
    FREE(indices_by_nsplices_gplus);
    FREE(nindices_by_nsplices_gminus);
    FREE(nindices_by_nsplices_gplus);

    debug(printf("Extension search returning %d plus and %d minus univdiagonals\n",
		 *nunivdiagonals_gplus,*nunivdiagonals_gminus));
    *last_method = EXT;
    return;

  } else if (*last_method < SEGMENT1 && method_goal >= SEGMENT1 &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_SEGMENT)) {
    /* 3. Segment search */
    debug(printf("%s Read: 3.  Running Segment search\n",first_read_p ? "5'" : "3'"));
    if (this->all_oligos_gen_filledp == false) {
      Stage1_fill_all_oligos_gen(this,querylength,genestrand);
    }
    segment_apply_quality_filter(this,quality_string,querylength);
    Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				  this,querylength,genestrand);
    Stage1_filter_segment_kmer_positions_gen(this,querylength,segment_max_avg_kmer_positions,
					     first_read_p ? "paired-5" : "paired-3");
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    Stage1_print_segment_kmer_positions_gen(this,querylength,first_read_p ? "paired-5" : "paired-3");
#endif

    Kmer_segment(&(*_univdiagonals_gplus),&(*auxinfo_gplus),
		 &(*pathstores_gplus),&(*nunivdiagonals_gplus),
		 &(*_univdiagonals_gminus),&(*auxinfo_gminus),
		 &(*pathstores_gminus),&(*nunivdiagonals_gminus),

		 this,querylength,max_nsplices,first_read_p,read_name,
		 first_read_p ? "paired-5" : "paired-3",
		 repetitive_ef64,univdiagpool,listpool,
		 univcoordlistpool,auxinfopool,pathstorepool);

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
    fprintf(stderr,"No enabled genomic method after %s\n",Method_string(*last_method));
    abort();
  }
}


/* A later method can add auxinfo to a diagonal already seen by an
   earlier method.  Path_setdiff merges that auxinfo into the retained
   all_* arrays and returns only truly novel diagonals, so we build a
   transient concordance frontier that includes both novel diagonals and
   changed old diagonals. */
static bool
make_univdiagonal_frontier (Univcoord_T **frontier_univdiagonals, Auxinfo_T **frontier_auxinfo,
			    Auxinfo_T ***frontier_auxinfo_ptrs, Pathstore_T **frontier_pathstores,
			    int *nfrontier,
			    Univcoord_T *novel_univdiagonals, Auxinfo_T *novel_auxinfo,
			    Pathstore_T *novel_pathstores, int nnovel,
			    Univcoord_T *old_univdiagonals, Auxinfo_T *old_auxinfo,
			    Pathstore_T *old_pathstores, int *duplicate_old_indices,
			    int nduplicates) {
  int i, j, k, old_i;

  if (nduplicates == 0) {
    *frontier_univdiagonals = novel_univdiagonals;
    *frontier_auxinfo = novel_auxinfo;
    *frontier_auxinfo_ptrs = (Auxinfo_T **) NULL;
    *frontier_pathstores = novel_pathstores;
    *nfrontier = nnovel;
    return false;

  } else {
    *nfrontier = nnovel + nduplicates;
    *frontier_univdiagonals = (Univcoord_T *) MALLOC_ALIGN((*nfrontier)*sizeof(Univcoord_T));
    *frontier_auxinfo = (Auxinfo_T *) MALLOC((*nfrontier)*sizeof(Auxinfo_T));
    *frontier_auxinfo_ptrs = (Auxinfo_T **) MALLOC((*nfrontier)*sizeof(Auxinfo_T *));
    *frontier_pathstores = (Pathstore_T *) MALLOC((*nfrontier)*sizeof(Pathstore_T));

    i = j = k = 0;
    while (i < nnovel && j < nduplicates) {
      old_i = duplicate_old_indices[j];
      if (novel_univdiagonals[i] < old_univdiagonals[old_i]) {
	(*frontier_univdiagonals)[k] = novel_univdiagonals[i];
	(*frontier_auxinfo)[k] = novel_auxinfo[i];
	(*frontier_auxinfo_ptrs)[k] = &(novel_auxinfo[i]);
	(*frontier_pathstores)[k] = novel_pathstores[i];
	i++;
      } else {
	(*frontier_univdiagonals)[k] = old_univdiagonals[old_i];
	(*frontier_auxinfo)[k] = old_auxinfo[old_i];
	(*frontier_auxinfo_ptrs)[k] = &(old_auxinfo[old_i]);
	(*frontier_pathstores)[k] = old_pathstores[old_i];
	j++;
      }
      k++;
    }
    while (i < nnovel) {
      (*frontier_univdiagonals)[k] = novel_univdiagonals[i];
      (*frontier_auxinfo)[k] = novel_auxinfo[i];
      (*frontier_auxinfo_ptrs)[k] = &(novel_auxinfo[i]);
      (*frontier_pathstores)[k] = novel_pathstores[i];
      i++; k++;
    }
    while (j < nduplicates) {
      old_i = duplicate_old_indices[j];
      (*frontier_univdiagonals)[k] = old_univdiagonals[old_i];
      (*frontier_auxinfo)[k] = old_auxinfo[old_i];
      (*frontier_auxinfo_ptrs)[k] = &(old_auxinfo[old_i]);
      (*frontier_pathstores)[k] = old_pathstores[old_i];
      j++; k++;
    }

    assert(k == *nfrontier);
    return true;
  }
}


static void
free_univdiagonal_frontier (Univcoord_T **frontier_univdiagonals, Auxinfo_T **frontier_auxinfo,
			    Auxinfo_T ***frontier_auxinfo_ptrs,
			    Pathstore_T **frontier_pathstores, bool allocatedp) {
  if (allocatedp == true) {
    FREE_ALIGN(*frontier_univdiagonals);
    FREE(*frontier_auxinfo);
    FREE(*frontier_auxinfo_ptrs);
    FREE(*frontier_pathstores);
  }

  *frontier_univdiagonals = (Univcoord_T *) NULL;
  *frontier_auxinfo = (Auxinfo_T *) NULL;
  *frontier_auxinfo_ptrs = (Auxinfo_T **) NULL;
  *frontier_pathstores = (Pathstore_T *) NULL;
  return;
}


static bool
univdiagonal_member_p (Univcoord_T univdiagonal, Univcoord_T *univdiagonals, int n) {
  int low = 0, high = n, mid;

  while (low < high) {
    mid = low + (high - low)/2;
    if (univdiagonals[mid] < univdiagonal) {
      low = mid + 1;
    } else {
      high = mid;
    }
  }

  return (low < n && univdiagonals[low] == univdiagonal) ? true : false;
}

static Auxinfo_T
auxinfo_from_exhaustive_elt (Singleelt_T *eltarray, int narray, int k,
			     Auxinfopool_T auxinfopool);

/* SEGMENT1 keeps only the best global segment candidates.  In paired
   mode, a mate can supply the missing context, so admit exhaustive
   SEGMENT1 candidates that are already concordant with retained mate
   candidates. */
static int
make_segment_exhaustive_concordance_frontier (Univcoord_T **target_univdiagonals,
					     Auxinfo_T **target_auxinfo,
					     Pathstore_T **target_pathstores,
					     bool target_is_L_p,
					     Univcoord_T *L_univdiagonals, int nL, int querylengthL,
					     Univcoord_T *H_univdiagonals, int nH,
					     Univcoord_T *target_all_univdiagonals, int ntarget_all,
					     Univcoord_T *target_selected_univdiagonals, int ntarget_selected,
					     Singleelt_T *target_eltarray, int ntarget_exhaustive,
					     int *target_counts,
					     int target_querylength,
					     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {
  int *indices;
  int nindices, i, k, anchor_i, target_i, best_target_i, best_anchor_i, best_count, best_width;
  int count, width, neligible;
  Chrpos_T promotion_distance;
  Univcoord_T univdiagonal;
  Singleelt_T elt;

  *target_univdiagonals = (Univcoord_T *) NULL;
  *target_auxinfo = (Auxinfo_T *) NULL;
  *target_pathstores = (Pathstore_T *) NULL;

  if (nL == 0 || nH == 0 || ntarget_exhaustive == 0 || target_eltarray == NULL) {
    return 0;
  } else {
    /* This repair is for local paired fragments that SEGMENT1 ranked out.
       Keep splice-aware concordance_distance for the main paired search, but
       use the configured fragment insert bound here so a distant repeat pair
       does not preempt the later splice/solve-all machinery. */
    promotion_distance = paired_max_insertlength > 0 ?
      (Chrpos_T) paired_max_insertlength : concordance_distance;
#ifdef LARGE_GENOMES
    indices = Intersect_approx_indices_uint8(&nindices,
					     L_univdiagonals,nL,/*diagterm1*/-querylengthL,
					     H_univdiagonals,nH,/*diagterm2*/0,
					     /*below_slop*/0,/*above_slop*/promotion_distance);
#else
    indices = Intersect_approx_indices_uint4(&nindices,
					     L_univdiagonals,nL,/*diagterm1*/-querylengthL,
					     H_univdiagonals,nH,/*diagterm2*/0,
					     /*below_slop*/0,/*above_slop*/promotion_distance);
#endif
  }

  if (nindices == 0) {
    FREE(indices);
    return 0;
  }

  best_target_i = best_anchor_i = -1;
  best_count = best_width = -1;
  neligible = 0;
  for (i = 0, k = 0; i < nindices; i++, k += 2) {
    target_i = target_is_L_p == true ? indices[k] : indices[k+1];
    anchor_i = target_is_L_p == true ? indices[k+1] : indices[k];
    if (target_i < 0 || target_i >= ntarget_exhaustive || target_eltarray[target_i] == NULL) {
      /* Skip */
    } else {
      univdiagonal = target_is_L_p == true ? L_univdiagonals[target_i] : H_univdiagonals[target_i];
      if (univdiagonal_member_p(univdiagonal,target_all_univdiagonals,ntarget_all) == true) {
	/* Already retained */
      } else if (univdiagonal_member_p(univdiagonal,target_selected_univdiagonals,ntarget_selected) == true) {
	/* Already selected by SEGMENT1 */
      } else {
	elt = target_eltarray[target_i];
	count = target_counts == NULL ? 0 : target_counts[target_i];
	width = elt->qend - elt->qstart;
	neligible++;
	if (count > best_count ||
	    (count == best_count && width > best_width) ||
	    (count == best_count && width == best_width && target_i < best_target_i)) {
	  best_target_i = target_i;
	  best_anchor_i = anchor_i;
	  best_count = count;
	  best_width = width;
	}
      }
    }
  }
  FREE(indices);

  if (best_target_i < 0) {
    return 0;
  }

  *target_univdiagonals = (Univcoord_T *) MALLOC_ALIGN(sizeof(Univcoord_T));
  *target_auxinfo = (Auxinfo_T *) MALLOC(sizeof(Auxinfo_T));
  *target_pathstores = (Pathstore_T *) MALLOC(sizeof(Pathstore_T));

  (*target_univdiagonals)[0] =
    target_is_L_p == true ? L_univdiagonals[best_target_i] : H_univdiagonals[best_target_i];
  (*target_auxinfo)[0] =
    auxinfo_from_exhaustive_elt(target_eltarray,ntarget_exhaustive,best_target_i,auxinfopool);
  assert((*target_auxinfo)[0] != NULL);
  (*target_pathstores)[0] = Pathstore_new(pathstorepool);

  Pathstore_assign_chrinfo(*target_univdiagonals,*target_pathstores,1,target_querylength);

#ifdef GSNAP_PAIRED_SOLVE_STATS
  fprintf(stderr,
	  "paired_solve_stats\tsegment_exhaustive_frontier\ttarget_is_L=%d\tpairs=%d\teligible=%d\tselected=1\ttarget_i=%d\tanchor_i=%d\tcount=%d\twidth=%d\tslop=%u\n",
	  target_is_L_p,nindices,neligible,best_target_i,best_anchor_i,best_count,best_width,
	  promotion_distance);
#endif

  return 1;
}

static bool
augment_univdiagonal_frontier (Univcoord_T **frontier_univdiagonals,
			       Auxinfo_T **frontier_auxinfo,
			       Auxinfo_T ***frontier_auxinfo_ptrs,
			       Pathstore_T **frontier_pathstores,
			       int *nfrontier, bool frontier_allocated_p,
			       Univcoord_T *extra_univdiagonals, Auxinfo_T *extra_auxinfo,
			       Pathstore_T *extra_pathstores, int nextra) {
  Univcoord_T *old_univdiagonals, *new_univdiagonals;
  Auxinfo_T *old_auxinfo, *new_auxinfo;
  Auxinfo_T **old_auxinfo_ptrs, **new_auxinfo_ptrs;
  Pathstore_T *old_pathstores, *new_pathstores;
  int nold, i, j, k;

  if (nextra == 0) {
    return frontier_allocated_p;
  }

  old_univdiagonals = *frontier_univdiagonals;
  old_auxinfo = *frontier_auxinfo;
  old_auxinfo_ptrs = *frontier_auxinfo_ptrs;
  old_pathstores = *frontier_pathstores;
  nold = *nfrontier;

  new_univdiagonals = (Univcoord_T *) MALLOC_ALIGN((nold + nextra) * sizeof(Univcoord_T));
  new_auxinfo = (Auxinfo_T *) MALLOC((nold + nextra) * sizeof(Auxinfo_T));
  new_auxinfo_ptrs = (Auxinfo_T **) MALLOC((nold + nextra) * sizeof(Auxinfo_T *));
  new_pathstores = (Pathstore_T *) MALLOC((nold + nextra) * sizeof(Pathstore_T));

  i = j = k = 0;
  while (i < nold && j < nextra) {
    if (old_univdiagonals[i] < extra_univdiagonals[j]) {
      new_univdiagonals[k] = old_univdiagonals[i];
      new_auxinfo[k] = old_auxinfo[i];
      new_auxinfo_ptrs[k] = old_auxinfo_ptrs == NULL ? &(old_auxinfo[i]) : old_auxinfo_ptrs[i];
      new_pathstores[k] = old_pathstores[i];
      i++;
    } else {
      new_univdiagonals[k] = extra_univdiagonals[j];
      new_auxinfo[k] = extra_auxinfo[j];
      new_auxinfo_ptrs[k] = &(extra_auxinfo[j]);
      new_pathstores[k] = extra_pathstores[j];
      j++;
    }
    k++;
  }
  while (i < nold) {
    new_univdiagonals[k] = old_univdiagonals[i];
    new_auxinfo[k] = old_auxinfo[i];
    new_auxinfo_ptrs[k] = old_auxinfo_ptrs == NULL ? &(old_auxinfo[i]) : old_auxinfo_ptrs[i];
    new_pathstores[k] = old_pathstores[i];
    i++; k++;
  }
  while (j < nextra) {
    new_univdiagonals[k] = extra_univdiagonals[j];
    new_auxinfo[k] = extra_auxinfo[j];
    new_auxinfo_ptrs[k] = &(extra_auxinfo[j]);
    new_pathstores[k] = extra_pathstores[j];
    j++; k++;
  }
  assert(k == nold + nextra);

  if (frontier_allocated_p == true) {
    FREE_ALIGN(old_univdiagonals);
    FREE(old_auxinfo);
    FREE(old_auxinfo_ptrs);
    FREE(old_pathstores);
  }

  *frontier_univdiagonals = new_univdiagonals;
  *frontier_auxinfo = new_auxinfo;
  *frontier_auxinfo_ptrs = new_auxinfo_ptrs;
  *frontier_pathstores = new_pathstores;
  *nfrontier = nold + nextra;

  return true;
}


/* Follows paired_search_trdiagonals */
static List_T
paired_search_univdiagonals (bool *complete5_p, bool *complete3_p,
			     int *found_score_paired, int *found_score_5, int *found_score_3,
			     int *end_score_paired, int *end_score_5, int *end_score_3,

			     Method_T *last_method_5, Method_T *last_method_3,
			     bool *any_imperfect_ends5_p, bool *any_imperfect_ends3_p,
		  
			     List_T *unresolved_pathpairs, List_T pathpairs, T this5, T this3,
			     
			     Shortread_T queryseq5, Shortread_T queryseq3,
			     char *queryuc_ptr_5, char *queryrc5, int querylength5,
			     char *queryuc_ptr_3, char *queryrc3, int querylength3,
			     Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      
			     int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,

			     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			     int localdb_nmismatches_allowed,

			     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			     int sufficient_score_5, int sufficient_score_3,
			     int max_nsplices_5, int max_nsplices_3, int genestrand,
			     
			     int nmismatches_filter_5, int nmismatches_filter_3,
			     int mincoverage_filter_5, int mincoverage_filter_3,

			     EF64_T repetitive_ef64,
			     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			     Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool, Intlistpool_T intlistpool,
			     Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			     Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3,
			     Cascade_Timing_T *cascade_timing, Method_T method_goal) {
  
  /* Method_T method_goal = SEGMENT1; */
  int kmer5_querystart, kmer5_queryend, kmer3_querystart, kmer3_queryend;

  Univcoord_T *_univdiagonals5_gplus = NULL, *_univdiagonals5_gminus = NULL, *_univdiagonals3_gplus = NULL, *_univdiagonals3_gminus = NULL;
  Auxinfo_T *auxinfo5_gplus = NULL, *auxinfo5_gminus = NULL, *auxinfo3_gplus = NULL, *auxinfo3_gminus = NULL;
  Pathstore_T *pathstores5_gplus = NULL, *pathstores5_gminus = NULL, *pathstores3_gplus = NULL, *pathstores3_gminus = NULL;
  int nunivdiagonals5_gplus = 0, nunivdiagonals5_gminus = 0, nunivdiagonals3_gplus = 0, nunivdiagonals3_gminus = 0;
  int nnovel5_gplus, nnovel5_gminus, nnovel3_gplus, nnovel3_gminus;
  int *duplicate5_gplus_indices = NULL, *duplicate5_gminus_indices = NULL;
  int *duplicate3_gplus_indices = NULL, *duplicate3_gminus_indices = NULL;
  int nduplicate5_gplus = 0, nduplicate5_gminus = 0, nduplicate3_gplus = 0, nduplicate3_gminus = 0;
  Univcoord_T *frontier5_gplus = NULL, *frontier5_gminus = NULL, *frontier3_gplus = NULL, *frontier3_gminus = NULL;
  Auxinfo_T *frontier_auxinfo5_gplus = NULL, *frontier_auxinfo5_gminus = NULL;
  Auxinfo_T *frontier_auxinfo3_gplus = NULL, *frontier_auxinfo3_gminus = NULL;
  Auxinfo_T **frontier_auxinfo5_gplus_ptrs = NULL, **frontier_auxinfo5_gminus_ptrs = NULL;
  Auxinfo_T **frontier_auxinfo3_gplus_ptrs = NULL, **frontier_auxinfo3_gminus_ptrs = NULL;
  Pathstore_T *frontier_pathstores5_gplus = NULL, *frontier_pathstores5_gminus = NULL;
  Pathstore_T *frontier_pathstores3_gplus = NULL, *frontier_pathstores3_gminus = NULL;
  int nfrontier5_gplus = 0, nfrontier5_gminus = 0, nfrontier3_gplus = 0, nfrontier3_gminus = 0;
  Univcoord_T *extra5_gplus = NULL, *extra5_gminus = NULL, *extra3_gplus = NULL, *extra3_gminus = NULL;
  Auxinfo_T *extra_auxinfo5_gplus = NULL, *extra_auxinfo5_gminus = NULL;
  Auxinfo_T *extra_auxinfo3_gplus = NULL, *extra_auxinfo3_gminus = NULL;
  Pathstore_T *extra_pathstores5_gplus = NULL, *extra_pathstores5_gminus = NULL;
  Pathstore_T *extra_pathstores3_gplus = NULL, *extra_pathstores3_gminus = NULL;
  int nextra5_gplus = 0, nextra5_gminus = 0, nextra3_gplus = 0, nextra3_gminus = 0;
  bool frontier5_gplus_allocated_p = false, frontier5_gminus_allocated_p = false;
  bool frontier3_gplus_allocated_p = false, frontier3_gminus_allocated_p = false;
  bool force_duplicate5_p = false, force_duplicate3_p = false;
  bool force_duplicate5_gplus_p = false, force_duplicate5_gminus_p = false;
  bool force_duplicate3_gplus_p = false, force_duplicate3_gminus_p = false;
  bool run5p = false, run3p = false;
  const char *stage_label;
  const char *paired_read_name;
  unsigned long long search_start_us, stage_start_us, elapsed_us;
  unsigned long long generation_us = 0ULL, other_us = 0ULL, concordance_us = 0ULL, merge_us = 0ULL;
#ifdef GSNAP_PAIRED_SOLVE_STATS
  double search_t0 = paired_solve_stats_now();
#endif

  /* int ndense5_gplus, ndense3_gplus, ndense5_gminus, ndense3_gminus; */

  search_start_us = CASCADE_TIMING_START(cascade_timing);
  paired_read_name = Shortread_accession(queryseq5);
  if (paired_read_name == NULL || paired_read_name[0] == '\0') {
    paired_read_name = Shortread_accession(queryseq3);
  }

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
    stage_label = paired_cascade_next_gen_stage_label(*last_method_5,method_goal);
    stage_start_us = CASCADE_STAGE_START_ROLE(stage_label,"mate5",queryseq5,univdiagpool);
    single_read_gen_univdiagonals(&(*last_method_5),&kmer5_querystart,&kmer5_queryend,*end_score_5,
				  &_univdiagonals5_gplus,&_univdiagonals5_gminus,
				  &auxinfo5_gplus,&auxinfo5_gminus,
				  &pathstores5_gplus,&pathstores5_gminus,
				  &nunivdiagonals5_gplus,&nunivdiagonals5_gminus,
				  this5,genestrand,query5_compress_fwd,query5_compress_rev,
				  querylength5,Shortread_quality_string(queryseq5),
				  paired_read_name,
				  max_nsplices_5,repetitive_ef64,
				  auxinfopool,pathstorepool,
				  univdiagpool,segmenteltpool,univcoordlistpool,intlistpool,listpool,
				  /*first_read_p*/true,method_goal);
    elapsed_us = CASCADE_STAGE_END_ROLE(stage_label,"mate5",queryseq5,univdiagpool,stage_start_us);
    generation_us += elapsed_us;
    if (cascade_timing != NULL) {
      paired_cascade_add_gen_method(cascade_timing,*last_method_5,elapsed_us);
    }
  }
  if (run3p == true) {
    stage_label = paired_cascade_next_gen_stage_label(*last_method_3,method_goal);
    stage_start_us = CASCADE_STAGE_START_ROLE(stage_label,"mate3",queryseq5,univdiagpool);
    single_read_gen_univdiagonals(&(*last_method_3),&kmer3_querystart,&kmer3_queryend,*end_score_3,
				  &_univdiagonals3_gplus,&_univdiagonals3_gminus,
				  &auxinfo3_gplus,&auxinfo3_gminus,
				  &pathstores3_gplus,&pathstores3_gminus,
				  &nunivdiagonals3_gplus,&nunivdiagonals3_gminus,
				  this3,genestrand,query3_compress_fwd,query3_compress_rev,
				  querylength3,Shortread_quality_string(queryseq3),
				  paired_read_name,
				  max_nsplices_3,repetitive_ef64,
				  auxinfopool,pathstorepool,
				  univdiagpool,segmenteltpool,univcoordlistpool,intlistpool,listpool,
				  /*first_read_p*/false,method_goal);
    elapsed_us = CASCADE_STAGE_END_ROLE(stage_label,"mate3",queryseq5,univdiagpool,stage_start_us);
    generation_us += elapsed_us;
    if (cascade_timing != NULL) {
      paired_cascade_add_gen_method(cascade_timing,*last_method_3,elapsed_us);
    }
  }


#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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


  stage_start_us = CASCADE_STAGE_START("goth",queryseq5,univdiagpool);

  /* Identify novel univdiagonals, so we can reuse results from the old Auxinfo_T objects */
  if (run5p == false) {
    nnovel5_gplus = nnovel5_gminus = 0;
  } else {
    nnovel5_gplus = Path_setdiff_univdiagonals_auxinfo(paired_duplicate_auxinfo_frontier_p == true ? &duplicate5_gplus_indices : NULL,
						       paired_duplicate_auxinfo_frontier_p == true ? &nduplicate5_gplus : NULL,
						       this5->all_univdiagonals_gplus,this5->all_auxinfo_gplus,
						       this5->all_pathstores_gplus,this5->all_nunivdiagonals_gplus,
						       _univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,
						       nunivdiagonals5_gplus,pathstorepool,intlistpool);
    nnovel5_gminus = Path_setdiff_univdiagonals_auxinfo(paired_duplicate_auxinfo_frontier_p == true ? &duplicate5_gminus_indices : NULL,
							paired_duplicate_auxinfo_frontier_p == true ? &nduplicate5_gminus : NULL,
							this5->all_univdiagonals_gminus,this5->all_auxinfo_gminus,
							this5->all_pathstores_gminus,this5->all_nunivdiagonals_gminus,
							_univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,
							nunivdiagonals5_gminus,pathstorepool,intlistpool);
  }

  if (run3p == false) {
    nnovel3_gplus = nnovel3_gminus = 0;
  } else {
    nnovel3_gplus = Path_setdiff_univdiagonals_auxinfo(paired_duplicate_auxinfo_frontier_p == true ? &duplicate3_gplus_indices : NULL,
						       paired_duplicate_auxinfo_frontier_p == true ? &nduplicate3_gplus : NULL,
						       this3->all_univdiagonals_gplus,this3->all_auxinfo_gplus,
						       this3->all_pathstores_gplus,this3->all_nunivdiagonals_gplus,
						       _univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,
						       nunivdiagonals3_gplus,pathstorepool,intlistpool);
    nnovel3_gminus = Path_setdiff_univdiagonals_auxinfo(paired_duplicate_auxinfo_frontier_p == true ? &duplicate3_gminus_indices : NULL,
							paired_duplicate_auxinfo_frontier_p == true ? &nduplicate3_gminus : NULL,
							this3->all_univdiagonals_gminus,this3->all_auxinfo_gminus,
							this3->all_pathstores_gminus,this3->all_nunivdiagonals_gminus,
							_univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,
							nunivdiagonals3_gminus,pathstorepool,intlistpool);
  }

  force_duplicate5_p = (*found_score_paired > *found_score_5 + *found_score_3 ||
			*found_score_5 > sufficient_score_5) ? true : false;
  force_duplicate3_p = (*found_score_paired > *found_score_5 + *found_score_3 ||
			*found_score_3 > sufficient_score_3) ? true : false;

  force_duplicate5_gplus_p = (force_duplicate5_p == true && nduplicate5_gplus > 0) ? true : false;
  force_duplicate5_gminus_p = (force_duplicate5_p == true && nduplicate5_gminus > 0) ? true : false;
  force_duplicate3_gplus_p = (force_duplicate3_p == true && nduplicate3_gplus > 0) ? true : false;
  force_duplicate3_gminus_p = (force_duplicate3_p == true && nduplicate3_gminus > 0) ? true : false;

  frontier5_gplus_allocated_p =
    make_univdiagonal_frontier(&frontier5_gplus,&frontier_auxinfo5_gplus,
			       &frontier_auxinfo5_gplus_ptrs,&frontier_pathstores5_gplus,
			       &nfrontier5_gplus,
			       _univdiagonals5_gplus,auxinfo5_gplus,pathstores5_gplus,nnovel5_gplus,
			       this5->all_univdiagonals_gplus,this5->all_auxinfo_gplus,
			       this5->all_pathstores_gplus,duplicate5_gplus_indices,
			       force_duplicate5_gplus_p == true ? nduplicate5_gplus : 0);
  frontier5_gminus_allocated_p =
    make_univdiagonal_frontier(&frontier5_gminus,&frontier_auxinfo5_gminus,
			       &frontier_auxinfo5_gminus_ptrs,&frontier_pathstores5_gminus,
			       &nfrontier5_gminus,
			       _univdiagonals5_gminus,auxinfo5_gminus,pathstores5_gminus,nnovel5_gminus,
			       this5->all_univdiagonals_gminus,this5->all_auxinfo_gminus,
			       this5->all_pathstores_gminus,duplicate5_gminus_indices,
			       force_duplicate5_gminus_p == true ? nduplicate5_gminus : 0);
  frontier3_gplus_allocated_p =
    make_univdiagonal_frontier(&frontier3_gplus,&frontier_auxinfo3_gplus,
			       &frontier_auxinfo3_gplus_ptrs,&frontier_pathstores3_gplus,
			       &nfrontier3_gplus,
			       _univdiagonals3_gplus,auxinfo3_gplus,pathstores3_gplus,nnovel3_gplus,
			       this3->all_univdiagonals_gplus,this3->all_auxinfo_gplus,
			       this3->all_pathstores_gplus,duplicate3_gplus_indices,
			       force_duplicate3_gplus_p == true ? nduplicate3_gplus : 0);
  frontier3_gminus_allocated_p =
    make_univdiagonal_frontier(&frontier3_gminus,&frontier_auxinfo3_gminus,
			       &frontier_auxinfo3_gminus_ptrs,&frontier_pathstores3_gminus,
			       &nfrontier3_gminus,
			       _univdiagonals3_gminus,auxinfo3_gminus,pathstores3_gminus,nnovel3_gminus,
			       this3->all_univdiagonals_gminus,this3->all_auxinfo_gminus,
			       this3->all_pathstores_gminus,duplicate3_gminus_indices,
			       force_duplicate3_gminus_p == true ? nduplicate3_gminus : 0);

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

  if (paired_seg1_exhaustive_concordance_p == true &&
      *found_score_paired > *found_score_5 + *found_score_3) {
    if (run5p == true && *last_method_5 == SEGMENT1 && *last_method_3 >= EXT) {
      nextra5_gplus =
	make_segment_exhaustive_concordance_frontier(&extra5_gplus,&extra_auxinfo5_gplus,
						    &extra_pathstores5_gplus,
						    /*target_is_L_p*/true,
						    /*L*/this5->exhaustive_gplus,this5->nexhaustive_gplus,
						    querylength5,
						    /*H*/this3->all_univdiagonals_gplus,
						    this3->all_nunivdiagonals_gplus,
						    this5->all_univdiagonals_gplus,
						    this5->all_nunivdiagonals_gplus,
						    _univdiagonals5_gplus,nnovel5_gplus,
						    this5->exhaustive_eltarray_gplus,
						    this5->nexhaustive_gplus,
						    this5->exhaustive_counts_gplus,querylength5,
						    auxinfopool,pathstorepool);
      frontier5_gplus_allocated_p =
	augment_univdiagonal_frontier(&frontier5_gplus,&frontier_auxinfo5_gplus,
				      &frontier_auxinfo5_gplus_ptrs,&frontier_pathstores5_gplus,
				      &nfrontier5_gplus,frontier5_gplus_allocated_p,
				      extra5_gplus,extra_auxinfo5_gplus,extra_pathstores5_gplus,
				      nextra5_gplus);

      nextra5_gminus =
	make_segment_exhaustive_concordance_frontier(&extra5_gminus,&extra_auxinfo5_gminus,
						    &extra_pathstores5_gminus,
						    /*target_is_L_p*/false,
						    /*L*/this3->all_univdiagonals_gminus,
						    this3->all_nunivdiagonals_gminus,
						    querylength3,
						    /*H*/this5->exhaustive_gminus,
						    this5->nexhaustive_gminus,
						    this5->all_univdiagonals_gminus,
						    this5->all_nunivdiagonals_gminus,
						    _univdiagonals5_gminus,nnovel5_gminus,
						    this5->exhaustive_eltarray_gminus,
						    this5->nexhaustive_gminus,
						    this5->exhaustive_counts_gminus,querylength5,
						    auxinfopool,pathstorepool);
      frontier5_gminus_allocated_p =
	augment_univdiagonal_frontier(&frontier5_gminus,&frontier_auxinfo5_gminus,
				      &frontier_auxinfo5_gminus_ptrs,&frontier_pathstores5_gminus,
				      &nfrontier5_gminus,frontier5_gminus_allocated_p,
				      extra5_gminus,extra_auxinfo5_gminus,extra_pathstores5_gminus,
				      nextra5_gminus);
    }

    if (run3p == true && *last_method_3 == SEGMENT1 && *last_method_5 >= EXT) {
      nextra3_gplus =
	make_segment_exhaustive_concordance_frontier(&extra3_gplus,&extra_auxinfo3_gplus,
						    &extra_pathstores3_gplus,
						    /*target_is_L_p*/false,
						    /*L*/this5->all_univdiagonals_gplus,
						    this5->all_nunivdiagonals_gplus,
						    querylength5,
						    /*H*/this3->exhaustive_gplus,
						    this3->nexhaustive_gplus,
						    this3->all_univdiagonals_gplus,
						    this3->all_nunivdiagonals_gplus,
						    _univdiagonals3_gplus,nnovel3_gplus,
						    this3->exhaustive_eltarray_gplus,
						    this3->nexhaustive_gplus,
						    this3->exhaustive_counts_gplus,querylength3,
						    auxinfopool,pathstorepool);
      frontier3_gplus_allocated_p =
	augment_univdiagonal_frontier(&frontier3_gplus,&frontier_auxinfo3_gplus,
				      &frontier_auxinfo3_gplus_ptrs,&frontier_pathstores3_gplus,
				      &nfrontier3_gplus,frontier3_gplus_allocated_p,
				      extra3_gplus,extra_auxinfo3_gplus,extra_pathstores3_gplus,
				      nextra3_gplus);

      nextra3_gminus =
	make_segment_exhaustive_concordance_frontier(&extra3_gminus,&extra_auxinfo3_gminus,
						    &extra_pathstores3_gminus,
						    /*target_is_L_p*/true,
						    /*L*/this3->exhaustive_gminus,
						    this3->nexhaustive_gminus,
						    querylength3,
						    /*H*/this5->all_univdiagonals_gminus,
						    this5->all_nunivdiagonals_gminus,
						    this3->all_univdiagonals_gminus,
						    this3->all_nunivdiagonals_gminus,
						    _univdiagonals3_gminus,nnovel3_gminus,
						    this3->exhaustive_eltarray_gminus,
						    this3->nexhaustive_gminus,
						    this3->exhaustive_counts_gminus,querylength3,
						    auxinfopool,pathstorepool);
      frontier3_gminus_allocated_p =
	augment_univdiagonal_frontier(&frontier3_gminus,&frontier_auxinfo3_gminus,
				      &frontier_auxinfo3_gminus_ptrs,&frontier_pathstores3_gminus,
				      &nfrontier3_gminus,frontier3_gminus_allocated_p,
				      extra3_gminus,extra_auxinfo3_gminus,extra_pathstores3_gminus,
				      nextra3_gminus);
    }
  }

  elapsed_us = CASCADE_STAGE_END("goth",queryseq5,univdiagpool,stage_start_us);
  other_us += elapsed_us;
  if (cascade_timing != NULL) {
    cascade_timing->gen_other_us += elapsed_us;
  }

  stage_start_us = CASCADE_STAGE_START("gconc",queryseq5,univdiagpool);
  if (run5p == true && run3p == true) {
    /* Novel5 vs Novel3 plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),
					  &(*found_score_paired),&(*found_score_5),&(*found_score_3),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/frontier5_gplus,frontier_auxinfo5_gplus,
					  frontier_auxinfo5_gplus_ptrs,frontier_pathstores5_gplus,
					  nfrontier5_gplus,
					  /*H*/frontier3_gplus,frontier_auxinfo3_gplus,
					  frontier_auxinfo3_gplus_ptrs,frontier_pathstores3_gplus,
					  nfrontier3_gplus,

					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  knownsplicing,knownindels,
			      
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					  
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					  sufficient_score_5,sufficient_score_3,genestrand,
			     
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false,
					  /*force_pending_auxinfo_L_p*/force_duplicate5_gplus_p,
					  /*force_pending_auxinfo_H_p*/force_duplicate3_gplus_p);

    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/frontier3_gminus,frontier_auxinfo3_gminus,
					  frontier_auxinfo3_gminus_ptrs,frontier_pathstores3_gminus,
					  nfrontier3_gminus,
					  /*H*/frontier5_gminus,frontier_auxinfo5_gminus,
					  frontier_auxinfo5_gminus_ptrs,frontier_pathstores5_gminus,
					  nfrontier5_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					  
					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  sufficient_score_3,sufficient_score_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true,
					  /*force_pending_auxinfo_L_p*/force_duplicate3_gminus_p,
					  /*force_pending_auxinfo_H_p*/force_duplicate5_gminus_p);

  } else if (run5p == true) {
    /* Novel5 vs Old3, plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),&(*found_score_paired),
					  &(*found_score_5),&(*found_score_3),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/frontier5_gplus,frontier_auxinfo5_gplus,
					  frontier_auxinfo5_gplus_ptrs,frontier_pathstores5_gplus,
					  nfrontier5_gplus,
					  /*H*/this3->all_univdiagonals_gplus,this3->all_auxinfo_gplus,
					  (Auxinfo_T **) NULL,this3->all_pathstores_gplus,
					  this3->all_nunivdiagonals_gplus,

					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					  
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					  sufficient_score_5,sufficient_score_3,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false,
					  /*force_pending_auxinfo_L_p*/force_duplicate5_gplus_p,
					  /*force_pending_auxinfo_H_p*/false);

    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/this3->all_univdiagonals_gminus,this3->all_auxinfo_gminus,
					  (Auxinfo_T **) NULL,this3->all_pathstores_gminus,
					  this3->all_nunivdiagonals_gminus,
					  /*H*/frontier5_gminus,frontier_auxinfo5_gminus,
					  frontier_auxinfo5_gminus_ptrs,frontier_pathstores5_gminus,
					  nfrontier5_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  sufficient_score_3,sufficient_score_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true,
					  /*force_pending_auxinfo_L_p*/false,
					  /*force_pending_auxinfo_H_p*/force_duplicate5_gminus_p);

  } else if (run3p == true) {
    /* Old5 vs Novel3, plus */
    pathpairs = concordance_univdiagonals(&(*complete5_p),&(*complete3_p),
					  &(*found_score_paired),&(*found_score_5),&(*found_score_3),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this5,/*stage1H*/this3,
					  /*L*/this5->all_univdiagonals_gplus,this5->all_auxinfo_gplus,
					  (Auxinfo_T **) NULL,this5->all_pathstores_gplus,
					  this5->all_nunivdiagonals_gplus,
					  /*H*/frontier3_gplus,frontier_auxinfo3_gplus,
					  frontier_auxinfo3_gplus_ptrs,frontier_pathstores3_gplus,
					  nfrontier3_gplus,
					  
					  queryseq5,queryseq3,
					  /*queryptrL*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  /*queryptrH*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_5,mismatch_positions_alloc_3,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					  
					  /*queryL_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					  /*queryH_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					  sufficient_score_5,sufficient_score_3,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen5,/*H*/spliceendsgen3,/*plusp*/true,
					  /*L_first_read_p*/true,/*H_first_read_p*/false,
					  /*force_pending_auxinfo_L_p*/false,
					  /*force_pending_auxinfo_H_p*/force_duplicate3_gplus_p);
    
    /* minus */
    pathpairs = concordance_univdiagonals(&(*complete3_p),&(*complete5_p),
					  &(*found_score_paired),&(*found_score_3),&(*found_score_5),
					  &(*end_score_paired),&(*end_score_5),&(*end_score_3),
					  &(*unresolved_pathpairs),pathpairs,/*stage1L*/this3,/*stage1H*/this5,
					  /*L*/frontier3_gminus,frontier_auxinfo3_gminus,
					  frontier_auxinfo3_gminus_ptrs,frontier_pathstores3_gminus,
					  nfrontier3_gminus,
					  /*H*/this5->all_univdiagonals_gminus,this5->all_auxinfo_gminus,
					  (Auxinfo_T **) NULL,this5->all_pathstores_gminus,
					  this5->all_nunivdiagonals_gminus,
					  
					  queryseq3,queryseq5,
					  /*queryptrL*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,max_nsplices_3,
					  /*queryptrH*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,max_nsplices_5,
					  knownsplicing,knownindels,
					  
					  mismatch_positions_alloc_3,mismatch_positions_alloc_5,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					  /*queryL_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					  /*queryH_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					  sufficient_score_3,sufficient_score_5,genestrand,
					  
					  nmismatches_filter_5,nmismatches_filter_3,
					  mincoverage_filter_5,mincoverage_filter_3,
					  
					  auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					  /*L*/spliceendsgen3,/*H*/spliceendsgen5,/*plusp*/false,
					  /*L_first_read_p*/false,/*H_first_read_p*/true,
					  /*force_pending_auxinfo_L_p*/force_duplicate3_gminus_p,
					  /*force_pending_auxinfo_H_p*/false);

  }
  concordance_us = CASCADE_STAGE_END("gconc",queryseq5,univdiagpool,stage_start_us);
  if (cascade_timing != NULL) {
    cascade_timing->gen_concordance_us += concordance_us;
  }
  debug(printf("end_scores are %d and %d\n",*end_score_5,*end_score_3));

  free_univdiagonal_frontier(&frontier5_gplus,&frontier_auxinfo5_gplus,
			     &frontier_auxinfo5_gplus_ptrs,&frontier_pathstores5_gplus,
			     frontier5_gplus_allocated_p);
  free_univdiagonal_frontier(&frontier5_gminus,&frontier_auxinfo5_gminus,
			     &frontier_auxinfo5_gminus_ptrs,&frontier_pathstores5_gminus,
			     frontier5_gminus_allocated_p);
  free_univdiagonal_frontier(&frontier3_gplus,&frontier_auxinfo3_gplus,
			     &frontier_auxinfo3_gplus_ptrs,&frontier_pathstores3_gplus,
			     frontier3_gplus_allocated_p);
  free_univdiagonal_frontier(&frontier3_gminus,&frontier_auxinfo3_gminus,
			     &frontier_auxinfo3_gminus_ptrs,&frontier_pathstores3_gminus,
			     frontier3_gminus_allocated_p);

  FREE(duplicate5_gplus_indices);
  FREE(duplicate5_gminus_indices);
  FREE(duplicate3_gplus_indices);
  FREE(duplicate3_gminus_indices);


  stage_start_us = CASCADE_STAGE_START("gmerge",queryseq5,univdiagpool);
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

  Path_merge_univdiagonals_auxinfo(&this5->all_univdiagonals_gplus,&this5->all_auxinfo_gplus,
				   &this5->all_pathstores_gplus,&this5->all_nunivdiagonals_gplus,
				   extra5_gplus,extra_auxinfo5_gplus,extra_pathstores5_gplus,nextra5_gplus,
				   pathstorepool,intlistpool);
  Path_merge_univdiagonals_auxinfo(&this5->all_univdiagonals_gminus,&this5->all_auxinfo_gminus,
				   &this5->all_pathstores_gminus,&this5->all_nunivdiagonals_gminus,
				   extra5_gminus,extra_auxinfo5_gminus,extra_pathstores5_gminus,nextra5_gminus,
				   pathstorepool,intlistpool);
  Path_merge_univdiagonals_auxinfo(&this3->all_univdiagonals_gplus,&this3->all_auxinfo_gplus,
				   &this3->all_pathstores_gplus,&this3->all_nunivdiagonals_gplus,
				   extra3_gplus,extra_auxinfo3_gplus,extra_pathstores3_gplus,nextra3_gplus,
				   pathstorepool,intlistpool);
  Path_merge_univdiagonals_auxinfo(&this3->all_univdiagonals_gminus,&this3->all_auxinfo_gminus,
				   &this3->all_pathstores_gminus,&this3->all_nunivdiagonals_gminus,
				   extra3_gminus,extra_auxinfo3_gminus,extra_pathstores3_gminus,nextra3_gminus,
				   pathstorepool,intlistpool);
  merge_us = CASCADE_STAGE_END("gmerge",queryseq5,univdiagpool,stage_start_us);
  if (cascade_timing != NULL) {
    cascade_timing->gen_merge_us += merge_us;
    elapsed_us = Cascade_timing_now_us() - search_start_us;
    if (elapsed_us > generation_us + other_us + concordance_us + merge_us) {
      cascade_timing->gen_other_us += elapsed_us - generation_us - other_us - concordance_us - merge_us;
    }
  }

#ifdef GSNAP_PAIRED_SOLVE_STATS
  fprintf(stderr,
	  "paired_solve_stats\tpaired_search_univdiagonals\telapsed=%.6f\trun5=%d\trun3=%d\tlast5=%s\tlast3=%s\tnovel5_plus=%d\tnovel5_minus=%d\tnovel3_plus=%d\tnovel3_minus=%d\tdup5_plus=%d\tdup5_minus=%d\tdup3_plus=%d\tdup3_minus=%d\textra5_plus=%d\textra5_minus=%d\textra3_plus=%d\textra3_minus=%d\tpathpairs=%d\tunresolved=%d\tfound_pair=%d\tfound5=%d\tfound3=%d\n",
	  paired_solve_stats_now() - search_t0,run5p,run3p,
	  Method_string(*last_method_5),Method_string(*last_method_3),
	  nnovel5_gplus,nnovel5_gminus,nnovel3_gplus,nnovel3_gminus,
	  nduplicate5_gplus,nduplicate5_gminus,nduplicate3_gplus,nduplicate3_gminus,
	  nextra5_gplus,nextra5_gminus,nextra3_gplus,nextra3_gminus,
	  List_length(pathpairs),List_length(*unresolved_pathpairs),
	  *found_score_paired,*found_score_5,*found_score_3);
  paired_solve_stats_print_stage("paired_search_5_after_merge",this5);
  paired_solve_stats_print_stage("paired_search_3_after_merge",this3);
#endif

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

    Auxinfo_gc(auxinfo5_gplus,nnovel5_gplus,univdiagpool,auxinfopool,listpool);
    Auxinfo_gc(auxinfo5_gminus,nnovel5_gminus,univdiagpool,auxinfopool,listpool);
    Auxinfo_gc(auxinfo3_gplus,nnovel3_gplus,univdiagpool,auxinfopool,listpool);
    Auxinfo_gc(auxinfo3_gminus,nnovel3_gminus,univdiagpool,auxinfopool,listpool);

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
find_inner_fusions_path5 (int *found_score_5, int *end_score_5, List_T pathpairs, Path_T path5, List_T singlepaths3,
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
      if ((newpath = Path_fusion_inner_qend(&(*found_score_5),&(*end_score_5),/*fusion5*/path5,/*anchor3*/path3,
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
      if ((newpath = Path_fusion_inner_qend(&(*found_score_5),&(*end_score_5),/*fusion5*/path5,/*anchor3*/path3,
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
find_inner_fusions_path3 (int *found_score_3, int *end_score_3, List_T pathpairs, Path_T path3, List_T singlepaths5,
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
      if ((newpath = Path_fusion_inner_qstart(&(*found_score_3),&(*end_score_3),/*fusion3*/path3,/*anchor5*/path5,
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
      if ((newpath = Path_fusion_inner_qstart(&(*found_score_3),&(*end_score_3),/*fusion3*/path3,/*anchor5*/path5,
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


#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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
push_converted_tr_path (Pathstore_T pathstore, Path_T path, bool sensep,
			Hitlistpool_T hitlistpool) {

  if (sensep == true) {
    pathstore->complete_sense_paths =
      Hitlist_push(pathstore->complete_sense_paths,hitlistpool,(void *) path
		   hitlistpool_trace(__FILE__,__LINE__));
  } else {
    pathstore->complete_antisense_paths =
      Hitlist_push(pathstore->complete_antisense_paths,hitlistpool,(void *) path
		   hitlistpool_trace(__FILE__,__LINE__));
  }

  return;
}

static void
finish_converted_tr_pathstore (Pathstore_T pathstore, int mincoverage_filter, Hitlistpool_T hitlistpool) {

  if (pathstore->complete_sense_paths != NULL) {
    pathstore->complete_sense_p = true;
    pathstore->best_sense_paths = Hitlist_copy(pathstore->complete_sense_paths,hitlistpool);
  } else {
    pathstore->best_sense_paths = Hitlist_copy(pathstore->unextended_sense_paths,hitlistpool);
  }

  if (pathstore->complete_antisense_paths != NULL) {
    pathstore->complete_antisense_p = true;
    pathstore->best_antisense_paths = Hitlist_copy(pathstore->complete_antisense_paths,hitlistpool);
  } else {
    pathstore->best_antisense_paths = Hitlist_copy(pathstore->unextended_antisense_paths,hitlistpool);
  }

  /* Below-threshold transcript paths remain output candidates, but should
     not prevent genomic auxinfo rescue. */
  pathstore->solvedp =
    (pathstore->complete_sense_paths != NULL || pathstore->complete_antisense_paths != NULL) &&
    pathlist_all_satisfy_coverage_p(pathstore->complete_sense_paths,mincoverage_filter) == true &&
    pathlist_all_satisfy_coverage_p(pathstore->complete_antisense_paths,mincoverage_filter) == true;

  return;
}


static void
create_univdiagonals_auxinfo_from_tr (int *found_score, int *end_score, Stage1_T this, Knownsplicing_T knownsplicing,
				      Shortread_T queryseq, int querylength,
				      List_T partial_sense_trpaths, List_T partial_antisense_trpaths,
				      List_T sense_trpaths, List_T antisense_trpaths,
				      Compress_T query_compress_fwd, Compress_T query_compress_rev,
				      int mincoverage_filter,

				      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
				      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
				      int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

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


  /* Convert partial sense_trpaths for genomic completion candidates */
  Trpath_convert_sense(&(*found_score),&(*end_score),
		       &unsolved_sense_pathlist_gplus,&unsolved_sense_pathlist_gminus,
		       &unsolved_antisense_pathlist_gplus,&unsolved_antisense_pathlist_gminus,

		       partial_sense_trpaths,first_read_p,
		       queryseq,querylength,this,knownsplicing,
		       query_compress_fwd,query_compress_rev,

		       novel_univdiags_alloc,novel_diagonals_alloc,
		       localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
		       localdb,localdb_nmismatches_allowed,

		       intlistpool,uintlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool,
		       /*crossover_sense_p*/true);

  /* Convert partial antisense_trpaths for genomic completion candidates */
  Trpath_convert_antisense(&(*found_score),&(*end_score),
			   &unsolved_antisense_pathlist_gplus,&unsolved_antisense_pathlist_gminus,
			   &unsolved_sense_pathlist_gplus,&unsolved_sense_pathlist_gminus,

			   partial_antisense_trpaths,first_read_p,
			   queryseq,querylength,this,knownsplicing,
			   query_compress_fwd,query_compress_rev,

			   novel_univdiags_alloc,novel_diagonals_alloc,
			   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			   localdb,localdb_nmismatches_allowed,

			   intlistpool,uintlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool,
			   /*crossover_sense_p*/true);

  /* Convert sense_trpaths */
  Trpath_convert_sense(&(*found_score),&(*end_score),
		       &sense_pathlist_gplus,&sense_pathlist_gminus,
		       &antisense_pathlist_gplus,&antisense_pathlist_gminus,

		       sense_trpaths,first_read_p,
		       queryseq,querylength,this,knownsplicing,
		       query_compress_fwd,query_compress_rev,

		       novel_univdiags_alloc,novel_diagonals_alloc,
		       localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
		       localdb,localdb_nmismatches_allowed,

		       intlistpool,uintlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool,
		       /*crossover_sense_p*/true);

  /* Convert antisense_trpaths */
  Trpath_convert_antisense(&(*found_score),&(*end_score),
			   &antisense_pathlist_gplus,&antisense_pathlist_gminus,
			   &sense_pathlist_gplus,&sense_pathlist_gminus,

			   antisense_trpaths,first_read_p,
			   queryseq,querylength,this,knownsplicing,
			   query_compress_fwd,query_compress_rev,

			   novel_univdiags_alloc,novel_diagonals_alloc,
			   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			   localdb,localdb_nmismatches_allowed,

			   intlistpool,uintlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool,
			   /*crossover_sense_p*/true);

  sense_pathlist_gplus = List_append(unsolved_sense_pathlist_gplus,sense_pathlist_gplus);
  sense_pathlist_gminus = List_append(unsolved_sense_pathlist_gminus,sense_pathlist_gminus);
  
  sense_paths_gplus = (Path_T *) List_to_array_n(&sense_npaths_gplus,sense_pathlist_gplus);
  sense_paths_gminus = (Path_T *) List_to_array_n(&sense_npaths_gminus,sense_pathlist_gminus);
  if (sense_npaths_gplus > 0) {
    qsort(sense_paths_gplus,sense_npaths_gplus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  }
  if (sense_npaths_gminus > 0) {
    qsort(sense_paths_gminus,sense_npaths_gminus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  }
  Hitlistpool_free_list(&sense_pathlist_gplus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&sense_pathlist_gminus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  
  antisense_pathlist_gplus = List_append(unsolved_antisense_pathlist_gplus,antisense_pathlist_gplus);
  antisense_pathlist_gminus = List_append(unsolved_antisense_pathlist_gminus,antisense_pathlist_gminus);

  antisense_paths_gplus = (Path_T *) List_to_array_n(&antisense_npaths_gplus,antisense_pathlist_gplus);
  antisense_paths_gminus = (Path_T *) List_to_array_n(&antisense_npaths_gminus,antisense_pathlist_gminus);
  if (antisense_npaths_gplus > 0) {
    qsort(antisense_paths_gplus,antisense_npaths_gplus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  }
  if (antisense_npaths_gminus > 0) {
    qsort(antisense_paths_gminus,antisense_npaths_gminus,sizeof(Path_T),Path_main_univdiagonal_cmp);
  }
  Hitlistpool_free_list(&antisense_pathlist_gplus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&antisense_pathlist_gminus,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));


  /* Build plus */
  if ((npaths = sense_npaths_gplus + antisense_npaths_gplus) > 0) {
    this->all_univdiagonals_gplus = (Univcoord_T *) MALLOC_ALIGN(npaths*sizeof(Univcoord_T));
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
	  push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
	}
	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	  push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
	}
	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	  push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
	}

	jj = j + 1;
	while (jj < antisense_npaths_gplus && antisense_paths_gplus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gplus[l];
	  push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
	}

	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
      }
      finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
      }
      finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

      j = jj; k++;
    }

    this->all_nunivdiagonals_gplus = k;
  }

  /* Build minus */
  if ((npaths = sense_npaths_gminus + antisense_npaths_gminus) > 0) {
    this->all_univdiagonals_gminus = (Univcoord_T *) MALLOC_ALIGN(npaths*sizeof(Univcoord_T));
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
	  push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
	}
	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	  push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
	}
	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	  push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
	}

	jj = j + 1;
	while (jj < antisense_npaths_gminus && antisense_paths_gminus[jj]->main_univdiagonal == univdiagonal_j) {
	  jj++;
	}
	for (l = j; l < jj; l++) {
	  path = antisense_paths_gminus[l];
	  push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
	}

	finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	push_converted_tr_path(pathstore,path,/*sensep*/true,hitlistpool);
      }
      finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
	push_converted_tr_path(pathstore,path,/*sensep*/false,hitlistpool);
      }
      finish_converted_tr_pathstore(pathstore,mincoverage_filter,hitlistpool);

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
		  int *end_score_paired, int *end_score_5, int *end_score_3,
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

		  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		  int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

		  Trdiagpool_T trdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		  Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {
  
  List_T pathpairs = NULL;
  Method_T last_method_5 = METHOD_INIT, last_method_3 = METHOD_INIT;
  Method_T method_goal;

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


  if (tr_chain_search_p == false) {
    /* A.  Initial search.  Try each side until we get a result */
    /* Disjunction: Stops when we have found pathpairs or both ends have gone to TR_PREVALENT (conjunction would mean either end) */
    method_goal = single_read_trdiagonal_method_goal();
    while (single_read_trdiagonal_method_available_p(last_method_5,method_goal) ||
	   single_read_trdiagonal_method_available_p(last_method_3,method_goal)) {
      if ((pathpairs = paired_search_trdiagonals(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						 &(*end_score_paired),&(*end_score_5),&(*end_score_3),
						 sufficient_score_5,sufficient_score_3,
						 &last_method_5,&last_method_3,

						 pathpairs,this5,this3,knownsplicing,
						 queryseq5,queryseq3,querylength5,querylength3,
						 query5_compress_fwd,query5_compress_rev,
						 query3_compress_fwd,query3_compress_rev,
						 mismatch_positions_alloc_5,mismatch_positions_alloc_3,

						 nmismatches_filter_5,nmismatches_filter_3,
						 mincoverage_filter_5,mincoverage_filter_3,

						 novel_univdiags_alloc,novel_diagonals_alloc,
						 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
						 localdb,localdb_nmismatches_allowed,

						 intlistpool,uintlistpool,univcoordlistpool,
						 listpool,trpathpool,pathpool,transcriptpool,vectorpool,
						 hitlistpool,method_goal)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	/* Path_print(((Pathpair_T) pathpairs->first)->path5); */
	/* Path_print(((Pathpair_T) pathpairs->first)->path3); */

	if (tr_only_extension_still_available_p(last_method_5,last_method_3) == true) {
	  debug(printf("Continuing TR-only search to enabled TR_EXT before accepting pathpairs\n"));
	  break;
	} else if (paired_transcript_completion_ready_p(*found_score_5,*found_score_3,
							*end_score_5,*end_score_3,
							last_method_5,last_method_3) == true) {
	  debug(printf("Deferring pathpair acceptance for seeded paired transcript completion\n"));
	  break;
	} else {
	  return pathpairs;
	}

      } else {
	debug1(printf("Continuing with found scores %d and %d\n",*found_score_5,*found_score_3));
      }
    }
  }
  
  
  /* B.  Advance by side separately to the selected TR path method */
  /* Helps to avoid genomic search */
  /* Disjunction: Stops when we have found pathpairs or both ends have reached the selected TR path method */
  method_goal = single_read_tr_method_goal();
  while (single_read_tr_method_available_p(last_method_5,method_goal) ||
	 single_read_tr_method_available_p(last_method_3,method_goal)) {
    if (last_method_3 >= method_goal) {
      if ((pathpairs = paired_read_next_method_tr_5(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    &(*end_score_paired),&(*end_score_5),&(*end_score_3),
						    sufficient_score_5,&last_method_5,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_5,nmismatches_allowed_5,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,

						    novel_univdiags_alloc,novel_diagonals_alloc,
						    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
						    localdb,localdb_nmismatches_allowed,

						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	if (tr_only_extension_still_available_p(last_method_5,last_method_3) == true) {
	  debug(printf("Continuing TR-only search to enabled %s before accepting pathpairs\n",
		       Method_string(method_goal)));
	} else if (paired_transcript_completion_ready_p(*found_score_5,*found_score_3,
							*end_score_5,*end_score_3,
							last_method_5,last_method_3) == true) {
	  debug(printf("Deferring pathpair acceptance for seeded paired transcript completion\n"));
	  break;
	} else if (tr_chain_search_p == true) {
	  debug(printf("Deferring TR-chain pathpair acceptance to convert TR-chain paths for genomic pairing\n"));
	  break;
	} else {
	  return pathpairs;
	}
      }
      
    } else if (last_method_5 >= method_goal) {
      if ((pathpairs = paired_read_next_method_tr_3(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    &(*end_score_paired),&(*end_score_5),&(*end_score_3),
						    sufficient_score_3,&last_method_3,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_3,nmismatches_allowed_3,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,

						    novel_univdiags_alloc,novel_diagonals_alloc,
						    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
						    localdb,localdb_nmismatches_allowed,

						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	if (tr_only_extension_still_available_p(last_method_5,last_method_3) == true) {
	  debug(printf("Continuing TR-only search to enabled %s before accepting pathpairs\n",
		       Method_string(method_goal)));
	} else if (paired_transcript_completion_ready_p(*found_score_5,*found_score_3,
							*end_score_5,*end_score_3,
							last_method_5,last_method_3) == true) {
	  debug(printf("Deferring pathpair acceptance for seeded paired transcript completion\n"));
	  break;
	} else if (tr_chain_search_p == true) {
	  debug(printf("Deferring TR-chain pathpair acceptance to convert TR-chain paths for genomic pairing\n"));
	  break;
	} else {
	  return pathpairs;
	}
      }
      
    } else if ((*found_score_5) >= (*found_score_3)) {
      if ((pathpairs = paired_read_next_method_tr_5(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    &(*end_score_paired),&(*end_score_5),&(*end_score_3),
						    sufficient_score_5,&last_method_5,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_5,nmismatches_allowed_5,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,

						    novel_univdiags_alloc,novel_diagonals_alloc,
						    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
						    localdb,localdb_nmismatches_allowed,

						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	if (tr_only_extension_still_available_p(last_method_5,last_method_3) == true) {
	  debug(printf("Continuing TR-only search to enabled %s before accepting pathpairs\n",
		       Method_string(method_goal)));
	} else if (paired_transcript_completion_ready_p(*found_score_5,*found_score_3,
							*end_score_5,*end_score_3,
							last_method_5,last_method_3) == true) {
	  debug(printf("Deferring pathpair acceptance for seeded paired transcript completion\n"));
	  break;
	} else if (tr_chain_search_p == true) {
	  debug(printf("Deferring TR-chain pathpair acceptance to convert TR-chain paths for genomic pairing\n"));
	  break;
	} else {
	  return pathpairs;
	}
      }
      
    } else {
      if ((pathpairs = paired_read_next_method_tr_3(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
						    &(*end_score_paired),&(*end_score_5),&(*end_score_3),
						    sufficient_score_3,&last_method_3,
						    
						    pathpairs,this5,this3,knownsplicing,queryseq5,queryseq3,
						    querylength5,querylength3,mismatch_positions_alloc_3,nmismatches_allowed_3,
						    query5_compress_fwd,query5_compress_rev,
						    query3_compress_fwd,query3_compress_rev,
						    nmismatches_filter_5,nmismatches_filter_3,
						    mincoverage_filter_5,mincoverage_filter_3,

						    novel_univdiags_alloc,novel_diagonals_alloc,
						    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
						    localdb,localdb_nmismatches_allowed,

						    genestrand,trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
						    listpool,trpathpool,pathpool,transcriptpool,
						    vectorpool,hitlistpool)) != NULL) {
	debug(printf("Exiting after methods %s and %s with %d pathpairs.  found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     Method_string(last_method_5),Method_string(last_method_3),List_length(pathpairs),
		     *found_score_paired,*found_score_5,*found_score_3));
	if (tr_only_extension_still_available_p(last_method_5,last_method_3) == true) {
	  debug(printf("Continuing TR-only search to enabled %s before accepting pathpairs\n",
		       Method_string(method_goal)));
	} else if (paired_transcript_completion_ready_p(*found_score_5,*found_score_3,
							*end_score_5,*end_score_3,
							last_method_5,last_method_3) == true) {
	  debug(printf("Deferring pathpair acceptance for seeded paired transcript completion\n"));
	  break;
	} else if (tr_chain_search_p == true) {
	  debug(printf("Deferring TR-chain pathpair acceptance to convert TR-chain paths for genomic pairing\n"));
	  break;
	} else {
	  return pathpairs;
	}
      }
    }
  }
  

  /* D.  Prep for genomic search.  Previously converted all trpaths to
     paths.  Now create univdiagonals and auxinfo in Stage1_T
     objects */
  if (paired_transcript_completion_p == true) {
    List_T new_sense_trpaths5, new_sense_trpaths3, new_antisense_trpaths5, new_antisense_trpaths3;
    bool complete_sense3_p, complete_sense5_p, complete_antisense3_p, complete_antisense5_p;
    bool complete_mate3_p, complete_mate5_p;

    debug(printf("Trying seeded paired transcript completion\n"));

    complete_mate3_p = (*found_score_3) > (*found_score_5) ||
      ((*found_score_3) == (*found_score_5) && (*end_score_3) > (*end_score_5));
    complete_mate5_p = (*found_score_5) > (*found_score_3) ||
      ((*found_score_5) == (*found_score_3) && (*end_score_5) > (*end_score_3));

    complete_sense3_p = complete_mate3_p == true && this5->sense_trpaths != NULL;
    complete_sense5_p = complete_mate5_p == true && this3->sense_trpaths != NULL;
    complete_antisense3_p = complete_mate3_p == true && this5->antisense_trpaths != NULL;
    complete_antisense5_p = complete_mate5_p == true && this3->antisense_trpaths != NULL;

    if (complete_sense3_p == true || complete_antisense3_p == true) {
      Stage1_fill_all_oligos_tr(this3,querylength3);
      Stage1_fill_all_positions_tr(this3,querylength3);
    }
    if (complete_sense5_p == true || complete_antisense5_p == true) {
      Stage1_fill_all_oligos_tr(this5,querylength5);
      Stage1_fill_all_positions_tr(this5,querylength5);
    }

    if (complete_sense3_p == true &&
	(new_sense_trpaths3 =
	 seeded_transcript_completion_for_mate(&(*found_score_3),&(*end_score_3),sufficient_score_3,
					       this3->sense_trpaths,this5->sense_trpaths,
					       this3,querylength3,mismatch_positions_alloc_3,
					       query3_compress_fwd,query3_compress_rev,
					       /*tplusp*/true,/*first_read_p*/false,
					       intlistpool,uintlistpool,listpool,
					       trpathpool,pathpool,hitlistpool)) != NULL) {
      pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
				 &(*end_score_paired),&(*end_score_5),&(*end_score_3),

				 pathpairs,/*sense_trpaths5*/NULL,/*new*/new_sense_trpaths3,
				 /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryseq3,
				 querylength5,querylength3,this5,this3,knownsplicing,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,

				 novel_univdiags_alloc,novel_diagonals_alloc,
				 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				 localdb,localdb_nmismatches_allowed,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				 /*sensedir*/SENSE_FORWARD,paired_tr_adjacent_indel_repair_p);
      this3->sense_trpaths = List_append(new_sense_trpaths3,this3->sense_trpaths);
    }

    if (complete_sense5_p == true &&
	(new_sense_trpaths5 =
	 seeded_transcript_completion_for_mate(&(*found_score_5),&(*end_score_5),sufficient_score_5,
					       this5->sense_trpaths,this3->sense_trpaths,
					       this5,querylength5,mismatch_positions_alloc_5,
					       query5_compress_fwd,query5_compress_rev,
					       /*tplusp*/true,/*first_read_p*/true,
					       intlistpool,uintlistpool,listpool,
					       trpathpool,pathpool,hitlistpool)) != NULL) {
      pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
				 &(*end_score_paired),&(*end_score_5),&(*end_score_3),

				 pathpairs,/*new*/new_sense_trpaths5,/*sense_trpaths3*/NULL,
				 /*trpaths5*/this5->sense_trpaths,/*trpaths3*/this3->sense_trpaths,

				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryseq3,
				 querylength5,querylength3,this5,this3,knownsplicing,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,

				 novel_univdiags_alloc,novel_diagonals_alloc,
				 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				 localdb,localdb_nmismatches_allowed,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				 /*sensedir*/SENSE_FORWARD,paired_tr_adjacent_indel_repair_p);
      this5->sense_trpaths = List_append(new_sense_trpaths5,this5->sense_trpaths);
    }

    if (complete_antisense3_p == true &&
	(new_antisense_trpaths3 =
	 seeded_transcript_completion_for_mate(&(*found_score_3),&(*end_score_3),sufficient_score_3,
					       this3->antisense_trpaths,this5->antisense_trpaths,
					       this3,querylength3,mismatch_positions_alloc_3,
					       query3_compress_fwd,query3_compress_rev,
					       /*tplusp*/false,/*first_read_p*/false,
					       intlistpool,uintlistpool,listpool,
					       trpathpool,pathpool,hitlistpool)) != NULL) {
      pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
				 &(*end_score_paired),&(*end_score_5),&(*end_score_3),

				 pathpairs,/*antisense_trpaths5*/NULL,/*new*/new_antisense_trpaths3,
				 /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,

				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryseq3,
				 querylength5,querylength3,this5,this3,knownsplicing,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,

				 novel_univdiags_alloc,novel_diagonals_alloc,
				 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				 localdb,localdb_nmismatches_allowed,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				 /*sensedir*/SENSE_ANTI,paired_tr_adjacent_indel_repair_p);
      this3->antisense_trpaths = List_append(new_antisense_trpaths3,this3->antisense_trpaths);
    }

    if (complete_antisense5_p == true &&
	(new_antisense_trpaths5 =
	 seeded_transcript_completion_for_mate(&(*found_score_5),&(*end_score_5),sufficient_score_5,
					       this5->antisense_trpaths,this3->antisense_trpaths,
					       this5,querylength5,mismatch_positions_alloc_5,
					       query5_compress_fwd,query5_compress_rev,
					       /*tplusp*/false,/*first_read_p*/true,
					       intlistpool,uintlistpool,listpool,
					       trpathpool,pathpool,hitlistpool)) != NULL) {
      pathpairs = Concordance_tr(&(*found_score_paired),&(*found_score_5),&(*found_score_3),
				 &(*end_score_paired),&(*end_score_5),&(*end_score_3),

				 pathpairs,/*new*/new_antisense_trpaths5,/*antisense_trpaths3*/NULL,
				 /*trpaths5*/this5->antisense_trpaths,/*trpaths3*/this3->antisense_trpaths,

				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 queryseq5,queryseq3,
				 querylength5,querylength3,this5,this3,knownsplicing,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,

				 novel_univdiags_alloc,novel_diagonals_alloc,
				 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				 localdb,localdb_nmismatches_allowed,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				 /*sensedir*/SENSE_ANTI,paired_tr_adjacent_indel_repair_p);
      this5->antisense_trpaths = List_append(new_antisense_trpaths5,this5->antisense_trpaths);
    }
  }

  if (pathpairs != NULL) {
    if (tr_chain_search_p == false) {
      debug(printf("Returning from paired_search_tr after seeded transcript completion with %d pathpairs\n",
		   List_length(pathpairs)));
      return pathpairs;
    } else if (genome_align_p == false) {
      debug(printf("Returning from paired_search_tr with %d TR-chain pathpairs because genome stages are disabled\n",
		   List_length(pathpairs)));
      return pathpairs;
    } else {
      debug(printf("Keeping %d TR-chain pathpairs and preparing converted TR-chain paths for genomic pairing\n",
		   List_length(pathpairs)));
    }
  }

  debug(printf("Preparing for genomic search\n"));

  create_univdiagonals_auxinfo_from_tr(&(*found_score_5),&(*end_score_5),this5,knownsplicing,
				       queryseq5,querylength5,
				       this5->partial_sense_trpaths,this5->partial_antisense_trpaths,
				       this5->sense_trpaths,this5->antisense_trpaths,
				       query5_compress_fwd,query5_compress_rev,
				       mincoverage_filter_5,

				       novel_univdiags_alloc,novel_diagonals_alloc,
				       localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				       localdb,localdb_nmismatches_allowed,

				       auxinfopool,pathstorepool,intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
				       /*first_read_p*/true);

  create_univdiagonals_auxinfo_from_tr(&(*found_score_3),&(*end_score_3),this3,knownsplicing,
				       queryseq3,querylength3,
				       this3->partial_sense_trpaths,this3->partial_antisense_trpaths,
				       this3->sense_trpaths,this3->antisense_trpaths,
				       query3_compress_fwd,query3_compress_rev,
				       mincoverage_filter_3,

				       novel_univdiags_alloc,novel_diagonals_alloc,
				       localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				       localdb,localdb_nmismatches_allowed,

				       auxinfopool,pathstorepool,intlistpool,uintlistpool,univcoordlistpool,
				       listpool,pathpool,vectorpool,transcriptpool,hitlistpool,
				       /*first_read_p*/false);

  debug(printf("Returning from paired_search_tr with %d pathpairs\n",List_length(pathpairs)));
  return pathpairs;
}


#if 0
static List_T
find_inner_fusions (int *found_score_5, int *found_score_3, int *end_score_5, int *end_score_3,

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
      pathpairs = find_inner_fusions_path5(&(*found_score_5),&(*end_score_5),
					   pathpairs,path5,local_sense_paths3,
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
      pathpairs = find_inner_fusions_path5(&(*found_score_5),&(*end_score_5),
					   pathpairs,path5,local_antisense_paths3,
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
      pathpairs = find_inner_fusions_path3(&(*found_score_3),&(*end_score_3),
					   pathpairs,path3,local_sense_paths5,
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
      pathpairs = find_inner_fusions_path3(&(*found_score_3),&(*end_score_3),
					   pathpairs,path3,local_antisense_paths5,
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
    genomicstart5 = Path_genomicstart(path5);
    genomicend5 = Path_genomicend(path5);
    genomicstart3 = Path_genomicstart(path3);
    genomicend3 = Path_genomicend(path3);

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
    genomicstart5 = Path_genomicstart(path5);
    genomicend5 = Path_genomicend(path5);
    genomicstart3 = Path_genomicstart(path3);
    genomicend3 = Path_genomicend(path3);

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

static bool
paired_reads_are_exact_full_overlap_p (Shortread_T queryseq5, Shortread_T queryseq3,
				       char *queryuc_ptr_5, char *queryrc5,
				       char *queryuc_ptr_3, char *queryrc3,
				       int querylength5, int querylength3);

static List_T
append_overlap_mate_mirror_paths (int *ncopied, List_T dest, List_T src,
				  int mincoverage_filter, int querylength, bool first_read_p,
				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

static Path_T *
filter_reportable_single_paths (Path_T *patharray, int *npaths_primary, int *npaths_altloc,
				Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				Listpool_T listpool, Pathpool_T pathpool,
				Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  Path_T path;
  int old_npaths_primary = *npaths_primary, old_npaths_altloc = *npaths_altloc;
  int i, k = 0, kept_primary = 0, kept_altloc = 0, total_npaths;

  total_npaths = old_npaths_primary + old_npaths_altloc;
  for (i = 0; i < total_npaths; i++) {
    path = patharray[i];
    if (Path_final_score_reportable_p(path) == true) {
      Path_set_final_score(path,/*pair_score*/(int) rint(Path_final_score(path)));
      patharray[k++] = path;
      if (i < old_npaths_primary) {
	kept_primary++;
      } else {
	kept_altloc++;
      }
    } else {
      debug(printf("Eliminating fallback single path by final score %.3f < score-min-read\n",
		   Path_final_score(path)));
      debug(Path_print(path));
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
  }

  *npaths_primary = kept_primary;
  *npaths_altloc = kept_altloc;
  if (k == 0) {
    FREE_OUT(patharray);
    return (Path_T *) NULL;
  } else {
    return patharray;
  }
}

static Pathpair_T *
consolidate_results (int *found_score_5, int *found_score_3, int *end_score_5, int *end_score_3, 
		     Pairtype_T *final_pairtype,

		     int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		     Path_T **patharray5, int *npaths5_primary, int *npaths5_altloc, int *first_absmq5, int *second_absmq5,
		     Path_T **patharray3, int *npaths3_primary, int *npaths3_altloc, int *first_absmq3, int *second_absmq3,

		     List_T sense_paths5, List_T antisense_paths5,
		     List_T sense_paths3, List_T antisense_paths3,
		     List_T rejected_singlepaths5, List_T rejected_singlepaths3,
		     List_T local_mate_singlepaths5, List_T local_mate_singlepaths3,

		     Shortread_T queryseq5, Shortread_T queryseq3,
		     char *queryuc_ptr_5, char *queryrc5, int querylength5,
		     char *queryuc_ptr_3, char *queryrc3, int querylength3,

		     Knownsplicing_T knownsplicing, Knownindels_T knownindels,

		     int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
		     Univdiag_T *novel_univdiags_alloc,Univcoord_T *novel_diagonals_alloc,
		     unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		     int max_nblocks, int localdb_nmismatches_allowed,

		     T this5, T this3,
		     Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		     Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

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
  int i, j, ncopied;
  bool use_local_mate_singlepaths5_p = false, use_local_mate_singlepaths3_p = false;


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

    } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT) &&
	       (complete_paths = Path_extend(&(*found_score_5),&(*end_score_5),&unextended_paths5,
					     /*original_path*/path,queryseq5,queryptr,querylength5,
					     this5->mismatch_positions_alloc,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,

					     this5,knownsplicing,knownindels,
					     query_compress,query5_compress_fwd,query5_compress_rev,/*genestrand*/0,
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

    } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT) &&
	       (complete_paths = Path_extend(&(*found_score_3),&(*end_score_3),&unextended_paths3,
					     /*original_path*/path,queryseq3,queryptr,querylength3,
					     this3->mismatch_positions_alloc,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,

					     this3,knownsplicing,knownindels,
					     query_compress,query3_compress_fwd,query3_compress_rev,/*genestrand*/0,
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

  if (paired_reads_are_exact_full_overlap_p(queryseq5,queryseq3,
					    queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					    querylength5,querylength3) == true) {
    if (singlepaths5 == NULL && rejected_singlepaths5 == NULL && local_mate_singlepaths5 == NULL &&
	local_mate_singlepaths3 != NULL) {
      ncopied = 0;
      singlepaths5 =
	append_overlap_mate_mirror_paths(&ncopied,singlepaths5,local_mate_singlepaths3,
					 mincoverage_filter_5,querylength5,/*first_read_p*/true,
					 query5_compress_fwd,query5_compress_rev,
					 intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
					 transcriptpool,hitlistpool);
      use_local_mate_singlepaths3_p = (ncopied > 0) ? true : false;
      debug16(printf("Mirrored %d complete overlap paths from mate3 to mate5\n",ncopied));
    } else if (singlepaths3 == NULL && rejected_singlepaths3 == NULL && local_mate_singlepaths3 == NULL &&
	       local_mate_singlepaths5 != NULL) {
      ncopied = 0;
      singlepaths3 =
	append_overlap_mate_mirror_paths(&ncopied,singlepaths3,local_mate_singlepaths5,
					 mincoverage_filter_3,querylength3,/*first_read_p*/false,
					 query3_compress_fwd,query3_compress_rev,
					 intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
					 transcriptpool,hitlistpool);
      use_local_mate_singlepaths5_p = (ncopied > 0) ? true : false;
      debug16(printf("Mirrored %d complete overlap paths from mate5 to mate3\n",ncopied));
    }
  }

  if (use_local_mate_singlepaths5_p == true) {
    singlepaths5 = List_append(local_mate_singlepaths5,singlepaths5);
  } else {
    Path_gc(&local_mate_singlepaths5,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
  }
  singlepaths5 = List_append(rejected_singlepaths5,singlepaths5);
  if (use_local_mate_singlepaths3_p == true) {
    singlepaths3 = List_append(local_mate_singlepaths3,singlepaths3);
  } else {
    Path_gc(&local_mate_singlepaths3,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
  }
  singlepaths3 = List_append(rejected_singlepaths3,singlepaths3);

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
				     /*max_softclip_filter*/-1,
				     /*split_paths_p*/false,/*split_path_max_overlap*/0,
				     /*split_path_min_segment*/0,/*split_path_min_combined_coverage*/0,
				     intlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,hitlistpool,/*filterp*/true);
    if (*patharray5 != NULL) {
      *patharray5 = filter_reportable_single_paths(*patharray5,&(*npaths5_primary),&(*npaths5_altloc),
						   intlistpool,univcoordlistpool,listpool,pathpool,
						   transcriptpool,hitlistpool);
    }

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
				     /*max_softclip_filter*/-1,
				     /*split_paths_p*/false,/*split_path_max_overlap*/0,
				     /*split_path_min_segment*/0,/*split_path_min_combined_coverage*/0,
				     intlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,hitlistpool,/*filterp*/true);
    if (*patharray3 != NULL) {
      *patharray3 = filter_reportable_single_paths(*patharray3,&(*npaths3_primary),&(*npaths3_altloc),
						   intlistpool,univcoordlistpool,listpool,pathpool,
						   transcriptpool,hitlistpool);
    }

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
    pathpairarray = Pathpair_eval_and_sort(&(*found_score_5),&(*found_score_3),&(*end_score_5),&(*end_score_3),
					   &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					   /*singlepaths5_from_rejected*/NULL,/*singlepaths3_from_rejected*/NULL,
					   pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					   /*quality_string_5*/Shortread_quality_string(queryseq5),
					   /*quality_string_3*/Shortread_quality_string(queryseq3),
					   this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					   localdb_tally,localdb_overflowp,max_nblocks,
					   localdb_nmismatches_allowed,

					   knownsplicing,knownindels,
					   nmismatches_filter_5,nmismatches_filter_3,
					   mincoverage_filter_5,mincoverage_filter_3,
					   querylength5,querylength3,
					   univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					   /*terminal_anchor_pathpairs_deferred*/NULL,
					   /*terminal_anchor_paths5_deferred*/NULL,
					   /*terminal_anchor_paths3_deferred*/NULL,
					   /*filterp*/true,
					   /*allow_postfallback_terminal_anchor_rescue_p*/false,
					   /*allow_near_exact_terminal_anchor_rescue_p*/false);

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


#if defined(DEBUG) || defined(DEBUG_PAIRED_FLOW)
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


#if 0
static List_T
paired_search_exhaustive (int *found_score_paired, int *found_score_5, int *found_score_3,
			  int *end_score_paired, int *end_score_5, int *end_score_3,
			  List_T *unresolved_pathpairs,

			  Shortread_T queryseq5, char *queryuc_ptr_5, char *queryrc5,
			  Shortread_T queryseq3, char *queryuc_ptr_3, char *queryrc3,
			  int querylength5, int querylength3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
			  
			  int *mismatch_positions_alloc_5, int *mismatch_positions_alloc_3,
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			  int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			  int localdb_nmismatches_allowed,

			  T this5, T this3, Knownsplicing_T knownsplicing, Knownindels_T knownindels,

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
  List_T chain, p;

  int *indices_gplus, *indices_gminus;
  int nindices_gplus, nindices_gminus, i, k;
  int index1, index2;
  int max_count, max_count_gplus, max_count_gminus, count;
  
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
	auxinfo5 = this5->exhaustive_auxinfo_gplus[index1];
	pathstore5 = Pathstore_new(pathstorepool);

	if (univdiagonal5 + auxinfo5->qstart < (Univcoord_T) querylength5) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/0,/*high_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qend);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qstart,
			       /*high_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qend);
	}

	assert(auxinfo5->chaindir == CHAIN_LEFT);
	for (p = auxinfo5->chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_left_chain(&(*found_score_5),&(*end_score_5),
				  
				     &pathstore5->unextended_sense_paths,&pathstore5->unextended_antisense_paths,
				     &pathstore5->complete_sense_paths,&pathstore5->complete_antisense_paths,

				     chain,queryseq5,/*queryptr*/queryuc_ptr_5,querylength5,mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,
				     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				     localdb_nmismatches_allowed,
				      
				     this5,knownsplicing,knownindels,
				     /*query_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
				     chrnum,chroffset,chrhigh,/*plusp*/true,/*genestrand*/0,
				     /*paired_end_p*/true,/*first_read_p*/true,intlistpool,
				     uintlistpool,univcoordlistpool,listpool,pathpool,
				     transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen5,
				     /*method*/EXHAUSTIVE,/*find_splices_p*/true);
	}

	pathstore5->solvedp = true;
	Pathstore_set_best_paths(pathstore5,hitlistpool); /* Needed for make_pathpairs to work */
      

	univdiagonal3 = this3->exhaustive_gplus[index2];
	auxinfo3 = this3->exhaustive_auxinfo_gplus[index2];
	pathstore3 = Pathstore_new(pathstorepool);

	if (univdiagonal3 + auxinfo3->qstart < (Univcoord_T) querylength3) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/0,/*high_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qend);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qstart,
			       /*high_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qend);
	}

	assert(auxinfo3->chaindir == CHAIN_RIGHT);
	for (p = auxinfo3->chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_right_chain(&(*found_score_3),&(*end_score_3),
				  
				     &pathstore3->unextended_sense_paths,&pathstore3->unextended_antisense_paths,
				     &pathstore3->complete_sense_paths,&pathstore3->complete_antisense_paths,

				     chain,queryseq3,/*queryptr*/queryuc_ptr_3,querylength3,mismatch_positions_alloc_3,
				     novel_univdiags_alloc,novel_diagonals_alloc,
				     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				     localdb_nmismatches_allowed,
				      
				     this3,knownsplicing,knownindels,
				     /*query_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
				     chrnum,chroffset,chrhigh,/*plusp*/true,/*genestrand*/0,
				     /*paired_end_p*/true,/*first_read_p*/false,intlistpool,
				     uintlistpool,univcoordlistpool,listpool,pathpool,
				     transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen3,
				     /*method*/EXHAUSTIVE,/*find_splices_p*/true);
	}

	pathstore3->solvedp = true;
	Pathstore_set_best_paths(pathstore3,hitlistpool); /* Needed for make_pathpairs to work */
      
	debug(printf("plus max count %d, %u %d..%d to %u %d..%d\n",
		     max_count,univdiagonal5,auxinfo5->qstart,auxinfo5->qend,
		     univdiagonal3,auxinfo3->qstart,auxinfo3->qend));
	debug(printf("(2) Calling make_pathpairs\n"));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstore5,pathstore3,queryseq5,queryseq3,/*plusp*/true,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/false,/*only_completeH_p*/false);

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
	auxinfo3 = this3->exhaustive_auxinfo_gminus[index1];
	pathstore3 = Pathstore_new(pathstorepool);

	if (univdiagonal3 + auxinfo3->qstart < (Univcoord_T) querylength3) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/0,/*high_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qend);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qstart,
			       /*high_univcoord*/univdiagonal3 - querylength3 + auxinfo3->qend);
	}

	assert(auxinfo5->chaindir == CHAIN_RIGHT);
	for (p = auxinfo5->chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_right_chain(&(*found_score_5),&(*end_score_5),
				  
				      &pathstore5->unextended_sense_paths,&pathstore5->unextended_antisense_paths,
				      &pathstore5->complete_sense_paths,&pathstore5->complete_antisense_paths,
				      
				      chain,queryseq3,/*queryptr*/queryrc3,querylength3,mismatch_positions_alloc_3,

				      novel_univdiags_alloc,novel_diagonals_alloc,
				      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				      localdb_nmismatches_allowed,
				      
				      this3,knownsplicing,knownindels,
				      /*query_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
				      chrnum,chroffset,chrhigh,/*plusp*/false,/*genestrand*/0,
				      /*paired_end_p*/true,/*first_read_p*/false,intlistpool,
				      uintlistpool,univcoordlistpool,listpool,pathpool,
				      transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen3,
				      /*method*/EXHAUSTIVE,/*find_splices_p*/true);
	}

	pathstore3->solvedp = true;
	Pathstore_set_best_paths(pathstore3,hitlistpool); /* Needed for make_pathpairs to work */


	univdiagonal5 = this5->exhaustive_gminus[index2];
	auxinfo5 = this5->exhaustive_auxinfo_gminus[index2];
	pathstore5 = Pathstore_new(pathstorepool);

	if (univdiagonal5 + auxinfo5->qstart < (Univcoord_T) querylength5) {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/0,/*high_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qend);
	} else {
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*low_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qstart,
			       /*high_univcoord*/univdiagonal5 - querylength5 + auxinfo5->qend);
	}

	assert(auxinfo3->chaindir == CHAIN_LEFT);
	for (p = auxinfo3->chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_left_chain(&(*found_score_3),&(*end_score_3),
				  
				     &pathstore3->unextended_sense_paths,&pathstore3->unextended_antisense_paths,
				     &pathstore3->complete_sense_paths,&pathstore3->complete_antisense_paths,

				     chain,queryseq5,/*queryptr*/queryrc5,querylength5,mismatch_positions_alloc_5,
				     novel_univdiags_alloc,novel_diagonals_alloc,
				     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				     localdb_nmismatches_allowed,
				      
				     this5,knownsplicing,knownindels,
				     /*query_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
				     chrnum,chroffset,chrhigh,/*plusp*/false,/*genestrand*/0,
				     /*paired_end_p*/true,/*first_read_p*/true,intlistpool,
				     uintlistpool,univcoordlistpool,listpool,pathpool,
				     transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen5,
				     /*method*/EXHAUSTIVE,/*find_splices_p*/true);
	}

	pathstore5->solvedp = true;
	Pathstore_set_best_paths(pathstore5,hitlistpool); /* Needed for make_pathpairs to work */


	debug(printf("minus max count %d, %u %d..%d to %u %d..%d\n",
		     max_count,univdiagonal3,auxinfo3->qstart,auxinfo3->qend,
		     univdiagonal5,auxinfo5->qstart,auxinfo5->qend));
	debug(printf("(3) Calling make_pathpairs\n"));
	pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),
				   pathpairs,pathstore3,pathstore5,queryseq3,queryseq5,/*plusp*/false,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,
				   intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				   transcriptpool,hitlistpool,
				   /*only_completeL_p*/false,/*only_completeH_p*/false);

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
#endif


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

static unsigned int
mate_minpure_specificity (bool *observedp, unsigned int *specificity_array, int nexhaustive, int k) {
  if (specificity_array == NULL || k < 0 || k >= nexhaustive || specificity_array[k] == 0U) {
    *observedp = false;
    return 0U;
  } else {
    *observedp = true;
    return specificity_array[k];
  }
}

static int
find_minpure_ranked_upstream (int *bestk, unsigned int *best_specificity,
			      bool *best_observedp,
			      Univcoord_T *exhaustive, int *counts, unsigned int *specificity_array,
			      int nexhaustive, int k, Univcoord_T high_univdiagonal,
			      int min_count, unsigned int min_specificity) {
  unsigned long long score, best_score = 0ULL;
  unsigned int specificity;
  bool observedp;
  int best_count = 0;

  *bestk = -1;
  *best_specificity = 0U;
  *best_observedp = false;
  while (k < nexhaustive && exhaustive[k] < high_univdiagonal) {
    if (counts[k] >= min_count) {
      specificity = mate_minpure_specificity(&observedp,specificity_array,nexhaustive,k);
      if (observedp == true && specificity >= min_specificity) {
	score = (unsigned long long) counts[k] * (unsigned long long) specificity;
	if (*bestk < 0 || score > best_score ||
	    (score == best_score && counts[k] > best_count)) {
	  *bestk = k;
	  best_count = counts[k];
	  best_score = score;
	  *best_specificity = specificity;
	  *best_observedp = observedp;
	}
      }
    }
    k++;
  }

  return best_count;
}

static int
find_minpure_ranked_downstream (int *bestk, unsigned int *best_specificity,
				bool *best_observedp,
				Univcoord_T *exhaustive, int *counts, unsigned int *specificity_array,
				int nexhaustive, int k, Univcoord_T low_univdiagonal,
				int min_count, unsigned int min_specificity) {
  unsigned long long score, best_score = 0ULL;
  unsigned int specificity;
  bool observedp;
  int best_count = 0;

  *bestk = -1;
  *best_specificity = 0U;
  *best_observedp = false;
  while (k >= 0 && exhaustive[k] > low_univdiagonal) {
    if (counts[k] >= min_count) {
      specificity = mate_minpure_specificity(&observedp,specificity_array,nexhaustive,k);
      if (observedp == true && specificity >= min_specificity) {
	score = (unsigned long long) counts[k] * (unsigned long long) specificity;
	if (*bestk < 0 || score > best_score ||
	    (score == best_score && counts[k] > best_count)) {
	  *bestk = k;
	  best_count = counts[k];
	  best_score = score;
	  *best_specificity = specificity;
	  *best_observedp = observedp;
	}
      }
    }
    k--;
  }

  return best_count;
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


static bool
pathstore_has_best_anchor_p (Pathstore_T pathstore) {
  return pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL;
}

static bool
pathstore_reliable_local_mate_anchor_p (Pathstore_T pathstore);

static bool
pathstore_mate_rescue_anchor_p (Pathstore_T pathstore, bool reliable_anchors_only_p) {
  /* The mate-rescue algorithm intentionally uses the broad 2025 anchor rule:
     any retained best path may seed a bounded local search for the mate. */
  if (reliable_anchors_only_p == true) {
    return pathstore_reliable_local_mate_anchor_p(pathstore);
  } else {
    return pathstore_has_best_anchor_p(pathstore);
  }
}

static int
local_mate_max_softclip (int querylength);

static bool
pathlist_has_reliable_local_mate_anchor_p (List_T paths) {
  Path_T path;
  List_T p;
  int min_anchor_coverage, max_anchor_score;

  for (p = paths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    min_anchor_coverage = path->querylength - local_mate_max_softclip(path->querylength);
    if (min_anchor_coverage < LOCAL_MATE_MIN_ALIGNED) {
      min_anchor_coverage = LOCAL_MATE_MIN_ALIGNED;
    }
    max_anchor_score = path->querylength/5;
    if (max_anchor_score < 12) {
      max_anchor_score = 12;
    }

    if (Path_coverage(path) >= min_anchor_coverage && path->found_score <= max_anchor_score) {
      return true;
    }
  }

  return false;
}

static bool
pathstore_reliable_local_mate_anchor_p (Pathstore_T pathstore) {
  return (pathlist_has_reliable_local_mate_anchor_p(pathstore->best_sense_paths) == true ||
	  pathlist_has_reliable_local_mate_anchor_p(pathstore->best_antisense_paths) == true) ?
    true : false;
}

static Auxinfo_T
auxinfo_from_exhaustive_elt (Singleelt_T *eltarray, int narray, int k,
			     Auxinfopool_T auxinfopool) {
  Singleelt_T elt;

  if (k < 0 || k >= narray || eltarray == NULL || eltarray[k] == NULL) {
    return (Auxinfo_T) NULL;
  } else {
    elt = eltarray[k];
    return Auxinfo_new_elt(SEGMENT1,elt->qstart,elt->qend,/*elti*/k,auxinfopool);
  }
}

static int
local_mate_max_softclip (int querylength) {
  int max_softclip = querylength/5;

  if (max_softclip < LOCAL_MATE_MIN_SOFTCLIP) {
    return LOCAL_MATE_MIN_SOFTCLIP;
  } else if (max_softclip > LOCAL_MATE_MAX_SOFTCLIP) {
    return LOCAL_MATE_MAX_SOFTCLIP;
  } else {
    return max_softclip;
  }
}

static bool
local_mate_segment_in_bounds_p (Univcoord_T univdiagonal, int qstart, int qend,
				int querylength, Pathstore_T anchor_pathstore) {
  return (univdiagonal >= anchor_pathstore->chroffset + querylength - qstart &&
	  univdiagonal <= anchor_pathstore->chrhigh + querylength - qend) ? true : false;
}

static bool
local_mate_better_p (int weighted_score, int nmismatches, int nindels,
		     int found_score_mate, int qstart, int qend,
		     int best_weighted_score, int best_nmismatches, int best_nindels,
		     int best_found_score, int best_qstart, int best_qend) {
  int ndiffs = nmismatches + nindels;
  int best_ndiffs = best_nmismatches + best_nindels;
  int coverage = qend - qstart;
  int best_coverage = best_qend - best_qstart;

  if (weighted_score < best_weighted_score) {
    return true;
  } else if (weighted_score > best_weighted_score) {
    return false;
  } else if (ndiffs < best_ndiffs) {
    return true;
  } else if (ndiffs > best_ndiffs) {
    return false;
  } else if (found_score_mate < best_found_score) {
    return true;
  } else if (found_score_mate > best_found_score) {
    return false;
  } else if (coverage > best_coverage) {
    return true;
  } else {
    return false;
  }
}

static Pathstore_T
bounded_mate_pathstore (Univcoord_T *univdiagonal_mate_out,
			int *found_score, int *end_score,
			int *mate_qstart_out, int *mate_qend_out,
			int *mate_nmismatches_out, int *mate_found_score_out,
			Shortread_T queryseq, int querylength,
			Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			Pathstore_T anchor_pathstore,
			Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			bool scan_from_low_p, bool plusp, bool first_read_p,
			bool allow_gapped_p,
			Pathstorepool_T pathstorepool,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {
  Path_T path, antipath;
  Pathstore_T pathstore_mate;
  Intlist_T endpoints, nmismatches_list, ref_nmismatches_list;
  Univcoordlist_T univdiagonals;
  List_T junctions;
  Junction_T junction;
  Univcoord_T scan_low, scan_high, refine_low, refine_high;
  Univcoord_T univdiagonal, univdiagonal2, univdiagonal3;
  Chrpos_T scan_span;
  int max_softclip, min_aligned, max_ungapped_weighted_score, max_gapped_weighted_score, max_gapped_ndiffs;
  int qstart, qend, softclip, candidate_softclip, clip_side;
  int nmismatches, ref_nmismatches, weighted_score, found_score_mate;
  int nindels, nindels1, nindels2, total_nindels;
  int indel_pos, indel_pos1, indel_pos2, nmismatches1, nmismatches2, nmismatches3;
  int ref_nmismatches1, ref_nmismatches2, ref_nmismatches3, qstart2;
  int junction_type1, junction_type2, qstart2_trial, qstart3_trial;
  int ignore_found_score, ignore_end_score;
  int best_qstart = -1, best_qend = -1, best_indel_pos = -1, best_indel_pos2 = -1;
  int best_nmismatches = querylength + 1, best_ref_nmismatches = querylength + 1;
  int best_nmismatches1 = querylength + 1, best_nmismatches2 = querylength + 1, best_nmismatches3 = querylength + 1;
  int best_ref_nmismatches1 = querylength + 1, best_ref_nmismatches2 = querylength + 1, best_ref_nmismatches3 = querylength + 1;
  int best_nindels = 0, best_nindels2 = 0, best_weighted_score = 2*querylength + 1;
  int best_found_score = querylength + 1, best_junction_type = 0, best_junction_type2 = 0;
  Univcoord_T best_univdiagonal = 0, best_univdiagonal2 = 0, best_univdiagonal3 = 0;

  if (high_univdiagonal < low_univdiagonal || paired_max_insertlength <= 0) {
    return (Pathstore_T) NULL;
  } else if (pathstore_reliable_local_mate_anchor_p(anchor_pathstore) == false) {
    debug11(printf("bounded mate skip: weak anchor\n"));
    return (Pathstore_T) NULL;
  } else {
    max_softclip = local_mate_max_softclip(querylength);
    min_aligned = querylength - max_softclip;
    if (min_aligned < LOCAL_MATE_MIN_ALIGNED) {
      min_aligned = LOCAL_MATE_MIN_ALIGNED;
    }
    max_ungapped_weighted_score = querylength/5;
    if (max_ungapped_weighted_score < 12) {
      max_ungapped_weighted_score = 12;
    }
    max_gapped_weighted_score = 2*querylength/5;
    if (max_gapped_weighted_score < 24) {
      max_gapped_weighted_score = 24;
    }
    max_gapped_ndiffs = querylength/8;
    if (max_gapped_ndiffs < 10) {
      max_gapped_ndiffs = 10;
    }

    scan_span = (Chrpos_T) paired_max_insertlength + 2*querylength +
      max_softclip + LOCAL_MATE_MAX_INDEL + LOCAL_MATE_MIN_SOFTCLIP;
    if (high_univdiagonal - low_univdiagonal > scan_span) {
      if (scan_from_low_p == true) {
	scan_low = low_univdiagonal;
	scan_high = add_bounded(scan_low,scan_span,high_univdiagonal + 1);
      } else {
	scan_high = high_univdiagonal;
	scan_low = subtract_bounded(scan_high,scan_span,low_univdiagonal);
      }
    } else {
      scan_low = low_univdiagonal;
      scan_high = high_univdiagonal;
    }
  }

  debug11(printf("bounded mate scan: %u..%u from %u..%u, plusp %d, first_read_p %d\n",
		 scan_low,scan_high,low_univdiagonal,high_univdiagonal,plusp,first_read_p));

  for (univdiagonal = scan_low; univdiagonal <= scan_high; univdiagonal++) {
    for (softclip = 0; softclip <= max_softclip; softclip++) {
      /* Left softclip, including softclip 0 for full-length ungapped alignment. */
      qstart = softclip;
      qend = querylength;
      if (qend - qstart >= min_aligned &&
	  local_mate_segment_in_bounds_p(univdiagonal,qstart,qend,querylength,anchor_pathstore) == true) {
	nmismatches = Path_solve_count_mismatches_substring(&ref_nmismatches,query_compress,
							    univdiagonal,querylength,
							    qstart,qend,plusp,/*genestrand*/0);
	weighted_score = 2*nmismatches + softclip;
	found_score_mate = nmismatches + softclip;
	if (local_mate_better_p(weighted_score,nmismatches,/*nindels*/0,
				found_score_mate,qstart,qend,
				best_weighted_score,best_nmismatches,best_nindels,
				best_found_score,best_qstart,best_qend) == true) {
	  best_weighted_score = weighted_score;
	  best_found_score = found_score_mate;
	  best_nmismatches = nmismatches;
	  best_ref_nmismatches = ref_nmismatches;
	  best_qstart = qstart;
	  best_qend = qend;
	  best_univdiagonal = univdiagonal;
	  best_junction_type = 0;
	  best_junction_type2 = 0;
	  best_nindels = best_nindels2 = 0;
	  best_indel_pos = best_indel_pos2 = -1;
	  best_univdiagonal2 = best_univdiagonal3 = 0;
	}

	if (allow_gapped_p == true) {
	  for (nindels = 1; nindels <= LOCAL_MATE_MAX_INDEL; nindels++) {
	    /* Insertion in the read. */
	    if (univdiagonal >= (Univcoord_T) nindels) {
	      univdiagonal2 = univdiagonal - nindels;
	      for (indel_pos = qstart + LOCAL_MATE_MIN_INDEL_SUPPORT;
		   indel_pos + nindels + LOCAL_MATE_MIN_INDEL_SUPPORT <= qend;
		   indel_pos++) {
		qstart2 = indel_pos + nindels;
		if (local_mate_segment_in_bounds_p(univdiagonal2,qstart2,qend,querylength,anchor_pathstore) == true) {
		  nmismatches1 = Path_solve_count_mismatches_substring(&ref_nmismatches1,query_compress,
								       univdiagonal,querylength,
								       qstart,indel_pos,plusp,/*genestrand*/0);
		  nmismatches2 = Path_solve_count_mismatches_substring(&ref_nmismatches2,query_compress,
								       univdiagonal2,querylength,
								       qstart2,qend,plusp,/*genestrand*/0);
		  nmismatches = nmismatches1 + nmismatches2;
		  if (nmismatches + nindels <= max_gapped_ndiffs) {
		    weighted_score = 2*(nmismatches + nindels) + softclip;
		    found_score_mate = nmismatches + softclip;
		    if (local_mate_better_p(weighted_score,nmismatches,nindels,
					    found_score_mate,qstart,qend,
					    best_weighted_score,best_nmismatches,best_nindels,
					    best_found_score,best_qstart,best_qend) == true) {
		      best_weighted_score = weighted_score;
		      best_found_score = found_score_mate;
		      best_nmismatches = nmismatches;
		      best_ref_nmismatches = ref_nmismatches1 + ref_nmismatches2;
		      best_nmismatches1 = nmismatches1;
		      best_nmismatches2 = nmismatches2;
		      best_ref_nmismatches1 = ref_nmismatches1;
		      best_ref_nmismatches2 = ref_nmismatches2;
		      best_qstart = qstart;
		      best_qend = qend;
		      best_indel_pos = indel_pos;
		      best_nindels = nindels;
		      best_univdiagonal = univdiagonal;
		      best_univdiagonal2 = univdiagonal2;
		      best_junction_type = INS_JUNCTION;
		      best_junction_type2 = 0;
		      best_nindels2 = 0;
		      best_indel_pos2 = -1;
		      best_univdiagonal3 = 0;
		    }
		  }
		}
	      }
	    }

	    /* Deletion in the reference. */
	    univdiagonal2 = univdiagonal + nindels;
	    for (indel_pos = qstart + LOCAL_MATE_MIN_INDEL_SUPPORT;
		 indel_pos + LOCAL_MATE_MIN_INDEL_SUPPORT <= qend;
		 indel_pos++) {
	      if (local_mate_segment_in_bounds_p(univdiagonal2,indel_pos,qend,querylength,anchor_pathstore) == true) {
		nmismatches1 = Path_solve_count_mismatches_substring(&ref_nmismatches1,query_compress,
								     univdiagonal,querylength,
								     qstart,indel_pos,plusp,/*genestrand*/0);
		nmismatches2 = Path_solve_count_mismatches_substring(&ref_nmismatches2,query_compress,
								     univdiagonal2,querylength,
								     indel_pos,qend,plusp,/*genestrand*/0);
		nmismatches = nmismatches1 + nmismatches2;
		if (nmismatches + nindels <= max_gapped_ndiffs) {
		  weighted_score = 2*(nmismatches + nindels) + softclip;
		  found_score_mate = nmismatches + softclip;
		  if (local_mate_better_p(weighted_score,nmismatches,nindels,
					  found_score_mate,qstart,qend,
					  best_weighted_score,best_nmismatches,best_nindels,
					  best_found_score,best_qstart,best_qend) == true) {
		    best_weighted_score = weighted_score;
		    best_found_score = found_score_mate;
		    best_nmismatches = nmismatches;
		    best_ref_nmismatches = ref_nmismatches1 + ref_nmismatches2;
		    best_nmismatches1 = nmismatches1;
		    best_nmismatches2 = nmismatches2;
		    best_ref_nmismatches1 = ref_nmismatches1;
		    best_ref_nmismatches2 = ref_nmismatches2;
		    best_qstart = qstart;
		    best_qend = qend;
		    best_indel_pos = indel_pos;
		    best_nindels = nindels;
		    best_univdiagonal = univdiagonal;
		    best_univdiagonal2 = univdiagonal2;
		    best_junction_type = DEL_JUNCTION;
		    best_junction_type2 = 0;
		    best_nindels2 = 0;
		    best_indel_pos2 = -1;
		    best_univdiagonal3 = 0;
		  }
		}
	      }
	    }
	  }
	}
      }

      if (softclip > 0) {
	/* Right softclip. */
	qstart = 0;
	qend = querylength - softclip;
	if (qend - qstart >= min_aligned &&
	    local_mate_segment_in_bounds_p(univdiagonal,qstart,qend,querylength,anchor_pathstore) == true) {
	  nmismatches = Path_solve_count_mismatches_substring(&ref_nmismatches,query_compress,
							      univdiagonal,querylength,
							      qstart,qend,plusp,/*genestrand*/0);
	  weighted_score = 2*nmismatches + softclip;
	  found_score_mate = nmismatches + softclip;
	  if (local_mate_better_p(weighted_score,nmismatches,/*nindels*/0,
				  found_score_mate,qstart,qend,
				  best_weighted_score,best_nmismatches,best_nindels,
				  best_found_score,best_qstart,best_qend) == true) {
	    best_weighted_score = weighted_score;
	    best_found_score = found_score_mate;
	    best_nmismatches = nmismatches;
	    best_ref_nmismatches = ref_nmismatches;
	    best_qstart = qstart;
	    best_qend = qend;
	    best_univdiagonal = univdiagonal;
	    best_junction_type = 0;
	    best_junction_type2 = 0;
	    best_nindels = best_nindels2 = 0;
	    best_indel_pos = best_indel_pos2 = -1;
	    best_univdiagonal2 = best_univdiagonal3 = 0;
	  }

	  if (allow_gapped_p == true) {
	    for (nindels = 1; nindels <= LOCAL_MATE_MAX_INDEL; nindels++) {
	      /* Insertion in the read. */
	      if (univdiagonal >= (Univcoord_T) nindels) {
		univdiagonal2 = univdiagonal - nindels;
		for (indel_pos = qstart + LOCAL_MATE_MIN_INDEL_SUPPORT;
		     indel_pos + nindels + LOCAL_MATE_MIN_INDEL_SUPPORT <= qend;
		     indel_pos++) {
		  qstart2 = indel_pos + nindels;
		  if (local_mate_segment_in_bounds_p(univdiagonal2,qstart2,qend,querylength,anchor_pathstore) == true) {
		    nmismatches1 = Path_solve_count_mismatches_substring(&ref_nmismatches1,query_compress,
									 univdiagonal,querylength,
									 qstart,indel_pos,plusp,/*genestrand*/0);
		    nmismatches2 = Path_solve_count_mismatches_substring(&ref_nmismatches2,query_compress,
									 univdiagonal2,querylength,
									 qstart2,qend,plusp,/*genestrand*/0);
		    nmismatches = nmismatches1 + nmismatches2;
		    if (nmismatches + nindels <= max_gapped_ndiffs) {
		      weighted_score = 2*(nmismatches + nindels) + softclip;
		      found_score_mate = nmismatches + softclip;
		      if (local_mate_better_p(weighted_score,nmismatches,nindels,
					      found_score_mate,qstart,qend,
					      best_weighted_score,best_nmismatches,best_nindels,
					      best_found_score,best_qstart,best_qend) == true) {
			best_weighted_score = weighted_score;
			best_found_score = found_score_mate;
			best_nmismatches = nmismatches;
			best_ref_nmismatches = ref_nmismatches1 + ref_nmismatches2;
			best_nmismatches1 = nmismatches1;
			best_nmismatches2 = nmismatches2;
			best_ref_nmismatches1 = ref_nmismatches1;
			best_ref_nmismatches2 = ref_nmismatches2;
			best_qstart = qstart;
			best_qend = qend;
			best_indel_pos = indel_pos;
			best_nindels = nindels;
			best_univdiagonal = univdiagonal;
			best_univdiagonal2 = univdiagonal2;
			best_junction_type = INS_JUNCTION;
			best_junction_type2 = 0;
			best_nindels2 = 0;
			best_indel_pos2 = -1;
			best_univdiagonal3 = 0;
		      }
		    }
		  }
		}
	      }

	      /* Deletion in the reference. */
	      univdiagonal2 = univdiagonal + nindels;
	      for (indel_pos = qstart + LOCAL_MATE_MIN_INDEL_SUPPORT;
		   indel_pos + LOCAL_MATE_MIN_INDEL_SUPPORT <= qend;
		   indel_pos++) {
		if (local_mate_segment_in_bounds_p(univdiagonal2,indel_pos,qend,querylength,anchor_pathstore) == true) {
		  nmismatches1 = Path_solve_count_mismatches_substring(&ref_nmismatches1,query_compress,
								       univdiagonal,querylength,
								       qstart,indel_pos,plusp,/*genestrand*/0);
		  nmismatches2 = Path_solve_count_mismatches_substring(&ref_nmismatches2,query_compress,
								       univdiagonal2,querylength,
								       indel_pos,qend,plusp,/*genestrand*/0);
		  nmismatches = nmismatches1 + nmismatches2;
		  if (nmismatches + nindels <= max_gapped_ndiffs) {
		    weighted_score = 2*(nmismatches + nindels) + softclip;
		    found_score_mate = nmismatches + softclip;
		    if (local_mate_better_p(weighted_score,nmismatches,nindels,
					    found_score_mate,qstart,qend,
					    best_weighted_score,best_nmismatches,best_nindels,
					    best_found_score,best_qstart,best_qend) == true) {
		      best_weighted_score = weighted_score;
		      best_found_score = found_score_mate;
		      best_nmismatches = nmismatches;
		      best_ref_nmismatches = ref_nmismatches1 + ref_nmismatches2;
		      best_nmismatches1 = nmismatches1;
		      best_nmismatches2 = nmismatches2;
		      best_ref_nmismatches1 = ref_nmismatches1;
		      best_ref_nmismatches2 = ref_nmismatches2;
		      best_qstart = qstart;
		      best_qend = qend;
		      best_indel_pos = indel_pos;
		      best_nindels = nindels;
		      best_univdiagonal = univdiagonal;
		      best_univdiagonal2 = univdiagonal2;
		      best_junction_type = DEL_JUNCTION;
		      best_junction_type2 = 0;
		      best_nindels2 = 0;
		      best_indel_pos2 = -1;
		      best_univdiagonal3 = 0;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if (allow_gapped_p == true &&
      best_univdiagonal != 0 &&
      best_weighted_score <= 3*max_gapped_weighted_score &&
      ((best_junction_type == 0 && best_weighted_score >= max_ungapped_weighted_score) ||
       (best_junction_type != 0 && best_weighted_score > max_gapped_weighted_score))) {
    refine_low = subtract_bounded(best_univdiagonal,LOCAL_MATE_MAX_TOTAL_INDEL,scan_low);
    refine_high = add_bounded(best_univdiagonal,LOCAL_MATE_MAX_TOTAL_INDEL,scan_high + 1);

    for (candidate_softclip = 0; candidate_softclip <= max_softclip; candidate_softclip++) {
      for (clip_side = 0; clip_side < (candidate_softclip == 0 ? 1 : 2); clip_side++) {
	if (clip_side == 0) {
	  qstart = candidate_softclip;
	  qend = querylength;
	} else {
	  qstart = 0;
	  qend = querylength - candidate_softclip;
	}
	if (qend - qstart < min_aligned) {
	  continue;
	}

	for (univdiagonal = refine_low; univdiagonal <= refine_high; univdiagonal++) {
	  for (nindels1 = 1; nindels1 <= LOCAL_MATE_MAX_TOTAL_INDEL; nindels1++) {
	    for (nindels2 = 1; nindels1 + nindels2 <= LOCAL_MATE_MAX_TOTAL_INDEL; nindels2++) {
	      total_nindels = nindels1 + nindels2;

	      for (junction_type1 = INS_JUNCTION; junction_type1 <= DEL_JUNCTION; junction_type1++) {
		if (junction_type1 == INS_JUNCTION) {
		  if (univdiagonal < (Univcoord_T) nindels1) {
		    continue;
		  } else {
		    univdiagonal2 = univdiagonal - nindels1;
		  }
		} else {
		  univdiagonal2 = univdiagonal + nindels1;
		}

		for (junction_type2 = INS_JUNCTION; junction_type2 <= DEL_JUNCTION; junction_type2++) {
		  for (indel_pos1 = qstart + LOCAL_MATE_MIN_INDEL_SUPPORT;
		       indel_pos1 + LOCAL_MATE_MIN_INDEL_SUPPORT < qend;
		       indel_pos1++) {
		    if (junction_type1 == INS_JUNCTION) {
		      qstart2_trial = indel_pos1 + nindels1;
		    } else {
		      qstart2_trial = indel_pos1;
		    }

		    if (junction_type2 == INS_JUNCTION) {
		      if (univdiagonal2 < (Univcoord_T) nindels2) {
			continue;
		      } else {
			univdiagonal3 = univdiagonal2 - nindels2;
		      }
		    } else {
		      univdiagonal3 = univdiagonal2 + nindels2;
		    }

		    for (indel_pos2 = qstart2_trial + 1;
			 indel_pos2 + LOCAL_MATE_MIN_INDEL_SUPPORT < qend;
			 indel_pos2++) {
		      if (junction_type2 == INS_JUNCTION) {
			qstart3_trial = indel_pos2 + nindels2;
		      } else {
			qstart3_trial = indel_pos2;
		      }

		      if (qstart2_trial >= indel_pos2 ||
			  qstart3_trial + LOCAL_MATE_MIN_INDEL_SUPPORT > qend) {
			continue;
		      } else if (local_mate_segment_in_bounds_p(univdiagonal,qstart,indel_pos1,
								querylength,anchor_pathstore) == false ||
				 local_mate_segment_in_bounds_p(univdiagonal2,qstart2_trial,indel_pos2,
								querylength,anchor_pathstore) == false ||
				 local_mate_segment_in_bounds_p(univdiagonal3,qstart3_trial,qend,
								querylength,anchor_pathstore) == false) {
			continue;
		      }

		      nmismatches1 =
			Path_solve_count_mismatches_substring(&ref_nmismatches1,query_compress,
							      univdiagonal,querylength,
							      qstart,indel_pos1,plusp,/*genestrand*/0);
		      nmismatches2 =
			Path_solve_count_mismatches_substring(&ref_nmismatches2,query_compress,
							      univdiagonal2,querylength,
							      qstart2_trial,indel_pos2,plusp,/*genestrand*/0);
		      nmismatches3 =
			Path_solve_count_mismatches_substring(&ref_nmismatches3,query_compress,
							      univdiagonal3,querylength,
							      qstart3_trial,qend,plusp,/*genestrand*/0);
		      nmismatches = nmismatches1 + nmismatches2 + nmismatches3;

		      if (nmismatches + total_nindels <= max_gapped_ndiffs) {
			softclip = qstart + querylength - qend;
			weighted_score = 2*(nmismatches + total_nindels) + softclip;
			found_score_mate = nmismatches + softclip;
			if (local_mate_better_p(weighted_score,nmismatches,total_nindels,
						found_score_mate,qstart,qend,
						best_weighted_score,best_nmismatches,
						best_nindels + best_nindels2,
						best_found_score,best_qstart,best_qend) == true) {
			  best_weighted_score = weighted_score;
			  best_found_score = found_score_mate;
			  best_nmismatches = nmismatches;
			  best_ref_nmismatches = ref_nmismatches1 + ref_nmismatches2 + ref_nmismatches3;
			  best_nmismatches1 = nmismatches1;
			  best_nmismatches2 = nmismatches2;
			  best_nmismatches3 = nmismatches3;
			  best_ref_nmismatches1 = ref_nmismatches1;
			  best_ref_nmismatches2 = ref_nmismatches2;
			  best_ref_nmismatches3 = ref_nmismatches3;
			  best_qstart = qstart;
			  best_qend = qend;
			  best_indel_pos = indel_pos1;
			  best_indel_pos2 = indel_pos2;
			  best_nindels = nindels1;
			  best_nindels2 = nindels2;
			  best_univdiagonal = univdiagonal;
			  best_univdiagonal2 = univdiagonal2;
			  best_univdiagonal3 = univdiagonal3;
			  best_junction_type = junction_type1;
			  best_junction_type2 = junction_type2;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if (best_univdiagonal == 0) {
    debug11(printf("bounded mate rejected: no candidate\n"));
    return (Pathstore_T) NULL;
  } else if (best_junction_type == 0 && best_weighted_score >= max_ungapped_weighted_score) {
    debug11(printf("bounded mate rejected: ungapped weighted %d >= max %d\n",
		   best_weighted_score,max_ungapped_weighted_score));
    return (Pathstore_T) NULL;
  } else if (best_junction_type != 0 && best_weighted_score > max_gapped_weighted_score) {
    debug11(printf("bounded mate rejected: gapped weighted %d > max %d\n",
		   best_weighted_score,max_gapped_weighted_score));
    return (Pathstore_T) NULL;
  } else if (*found_score < querylength && best_found_score > *found_score + 2) {
    debug11(printf("bounded mate rejected: found %d > current %d + 2\n",
		   best_found_score,*found_score));
    return (Pathstore_T) NULL;
  }

  debug11(printf("bounded mate selected: %u/%u/%u q%d..q%d indels %d:%d at %d, %d:%d at %d nm %d weighted %d found %d\n",
		 best_univdiagonal,best_univdiagonal2,best_univdiagonal3,best_qstart,best_qend,
		 best_junction_type,best_nindels,best_indel_pos,
		 best_junction_type2,best_nindels2,best_indel_pos2,best_nmismatches,
		 best_weighted_score,best_found_score));

  pathstore_mate = Pathstore_new(pathstorepool);
  pathstore_mate->chrnum = anchor_pathstore->chrnum;
  pathstore_mate->chroffset = anchor_pathstore->chroffset;
  pathstore_mate->chrhigh = anchor_pathstore->chrhigh;

  if (best_junction_type == 0) {
    path = Path_new_exact(best_univdiagonal,best_qstart,best_qend,
			  best_nmismatches,best_ref_nmismatches,
			  plusp,first_read_p,/*genestrand*/0,querylength,best_found_score,
			  anchor_pathstore->chrnum,anchor_pathstore->chroffset,anchor_pathstore->chrhigh,
			  query_compress_fwd,query_compress_rev,
			  intlistpool,univcoordlistpool,pathpool,/*method*/LOCAL_MATE);
  } else {
    endpoints = Intlistpool_push(NULL,intlistpool,best_qend
				 intlistpool_trace(__FILE__,__LINE__));
    if (best_junction_type2 != 0) {
      endpoints = Intlistpool_push(endpoints,intlistpool,best_indel_pos2
				   intlistpool_trace(__FILE__,__LINE__));
    }
    endpoints = Intlistpool_push(endpoints,intlistpool,best_indel_pos
				 intlistpool_trace(__FILE__,__LINE__));
    endpoints = Intlistpool_push(endpoints,intlistpool,best_qstart
				 intlistpool_trace(__FILE__,__LINE__));

    if (best_junction_type2 != 0) {
      univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,best_univdiagonal3
					     univcoordlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,best_univdiagonal2
					     univcoordlistpool_trace(__FILE__,__LINE__));
    } else {
      univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,best_univdiagonal2
					     univcoordlistpool_trace(__FILE__,__LINE__));
    }
    univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,best_univdiagonal
					   univcoordlistpool_trace(__FILE__,__LINE__));

    if (best_junction_type2 != 0) {
      nmismatches_list = Intlistpool_push(NULL,intlistpool,best_nmismatches3
					  intlistpool_trace(__FILE__,__LINE__));
      nmismatches_list = Intlistpool_push(nmismatches_list,intlistpool,best_nmismatches2
					  intlistpool_trace(__FILE__,__LINE__));
    } else {
      nmismatches_list = Intlistpool_push(NULL,intlistpool,best_nmismatches2
					  intlistpool_trace(__FILE__,__LINE__));
    }
    nmismatches_list = Intlistpool_push(nmismatches_list,intlistpool,best_nmismatches1
					intlistpool_trace(__FILE__,__LINE__));
    if (best_junction_type2 != 0) {
      ref_nmismatches_list = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches3
					      intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches_list = Intlistpool_push(ref_nmismatches_list,intlistpool,best_ref_nmismatches2
					      intlistpool_trace(__FILE__,__LINE__));
    } else {
      ref_nmismatches_list = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches2
					      intlistpool_trace(__FILE__,__LINE__));
    }
    ref_nmismatches_list = Intlistpool_push(ref_nmismatches_list,intlistpool,best_ref_nmismatches1
					    intlistpool_trace(__FILE__,__LINE__));

    junctions = (List_T) NULL;
    if (best_junction_type2 == INS_JUNCTION) {
      junction = Junction_new_insertion(best_nindels2,pathpool);
      junctions = Listpool_push(junctions,listpool,(void *) junction
				listpool_trace(__FILE__,__LINE__));
    } else if (best_junction_type2 == DEL_JUNCTION) {
      junction = Junction_new_deletion(best_nindels2,pathpool);
      junctions = Listpool_push(junctions,listpool,(void *) junction
				listpool_trace(__FILE__,__LINE__));
    }

    if (best_junction_type == INS_JUNCTION) {
      junction = Junction_new_insertion(best_nindels,pathpool);
    } else {
      junction = Junction_new_deletion(best_nindels,pathpool);
    }
    junctions = Listpool_push(junctions,listpool,(void *) junction
			      listpool_trace(__FILE__,__LINE__));

    path = Path_convert_simple(anchor_pathstore->chrnum,anchor_pathstore->chroffset,
			       anchor_pathstore->chrhigh,endpoints,univdiagonals,
			       nmismatches_list,ref_nmismatches_list,junctions,/*evalp*/true,
			       plusp,first_read_p,/*genestrand*/0,SENSE_FORWARD,querylength,
			       query_compress_fwd,query_compress_rev,listpool,pathpool,/*method*/LOCAL_MATE);
    ignore_found_score = ignore_end_score = querylength;
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,path,query_compress_fwd,query_compress_rev);
  }
  path->completep = true;
  pathstore_mate->complete_sense_paths =
    Hitlist_push(pathstore_mate->complete_sense_paths,hitlistpool,(void *) path
		 hitlistpool_trace(__FILE__,__LINE__));

  if (splicingp == true) {
    antipath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			 pathpool,vectorpool,transcriptpool,hitlistpool);
    antipath->sensedir = SENSE_ANTI;
    antipath->completep = true;
    pathstore_mate->complete_antisense_paths =
      Hitlist_push(pathstore_mate->complete_antisense_paths,hitlistpool,(void *) antipath
		   hitlistpool_trace(__FILE__,__LINE__));
  }

  pathstore_mate->solvedp = true;
  Pathstore_set_best_paths(pathstore_mate,hitlistpool);

  *found_score = best_found_score < *found_score ? best_found_score : *found_score;
  *end_score = 0;
  *univdiagonal_mate_out = best_univdiagonal;
  *mate_qstart_out = best_qstart;
  *mate_qend_out = best_qend;
  *mate_nmismatches_out = best_nmismatches;
  *mate_found_score_out = best_found_score;

  (void) queryseq;
  (void) uintlistpool;
  return pathstore_mate;
}

static bool
pathlist_has_local_mate_coverage_p (List_T paths, int querylength) {
  Path_T path;
  List_T p;
  int min_aligned = querylength - local_mate_max_softclip(querylength);
  int max_found_score = querylength/5;

  if (min_aligned < LOCAL_MATE_MIN_ALIGNED) {
    min_aligned = LOCAL_MATE_MIN_ALIGNED;
  }
  if (max_found_score < 12) {
    max_found_score = 12;
  }

  for (p = paths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (Path_coverage(path) >= min_aligned &&
	path->score_within_trims >= 0 && path->found_score <= max_found_score) {
      return true;
    }
  }

  return false;
}

static bool
pathstore_has_local_mate_coverage_p (Pathstore_T pathstore, int querylength) {
  return (pathlist_has_local_mate_coverage_p(pathstore->best_sense_paths,querylength) == true ||
	  pathlist_has_local_mate_coverage_p(pathstore->best_antisense_paths,querylength) == true) ?
    true : false;
}

static List_T
append_local_mate_singlepath_list (int *ncopied, List_T singlepaths, List_T paths, int mincoverage_filter,
				   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  List_T p;
  Path_T path, copy;

  for (p = paths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (Path_coverage(path) >= mincoverage_filter && path->score_within_trims >= 0) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      singlepaths = Hitlist_push(singlepaths,hitlistpool,(void *) copy
				 hitlistpool_trace(__FILE__,__LINE__));
      *ncopied += 1;
    }
  }

  return singlepaths;
}

static List_T
append_local_mate_singlepaths (int *ncopied, List_T singlepaths, Pathstore_T pathstore, int mincoverage_filter,
			       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			       Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  singlepaths =
    append_local_mate_singlepath_list(&(*ncopied),singlepaths,pathstore->complete_sense_paths,mincoverage_filter,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);
  singlepaths =
    append_local_mate_singlepath_list(&(*ncopied),singlepaths,pathstore->complete_antisense_paths,mincoverage_filter,
				      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				      transcriptpool,hitlistpool);

  return singlepaths;
}

static List_T
append_local_mate_singlepaths_policy (int *ncopied, List_T singlepaths, Pathstore_T pathstore,
				      int mincoverage_filter, bool retain_local_mate_singlepaths_p,
				      Mate_rescue_stats_T *rescue_stats,
				      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				      Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  *ncopied = 0;
  MATE_RESCUE_STATS_INC(rescue_stats,local_mate_singlepath_copy_calls);

  if (retain_local_mate_singlepaths_p == false) {
    MATE_RESCUE_STATS_INC(rescue_stats,local_mate_singlepath_copy_skipped_unusable);
    return singlepaths;
  } else {
    singlepaths =
      append_local_mate_singlepaths(&(*ncopied),singlepaths,pathstore,mincoverage_filter,
				    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				    transcriptpool,hitlistpool);
    MATE_RESCUE_STATS_ADD(rescue_stats,local_mate_singlepaths_copied,(unsigned long long) *ncopied);
    return singlepaths;
  }
}

static bool
paired_reads_are_exact_full_overlap_p (Shortread_T queryseq5, Shortread_T queryseq3,
				       char *queryuc_ptr_5, char *queryrc5,
				       char *queryuc_ptr_3, char *queryrc3,
				       int querylength5, int querylength3) {
  if (querylength5 != querylength3) {
    return false;
  } else if (Shortread_invertedp(queryseq5) != Shortread_invertedp(queryseq3)) {
    return (strncmp(queryuc_ptr_5,queryuc_ptr_3,querylength5) == 0 &&
	    strncmp(queryrc5,queryrc3,querylength5) == 0) ? true : false;
  } else {
    return (strncmp(queryuc_ptr_5,queryrc3,querylength5) == 0 &&
	    strncmp(queryrc5,queryuc_ptr_3,querylength5) == 0) ? true : false;
  }
}

static Univcoord_T
mate_main_univdiagonal (Path_T path, int querylength, bool first_read_p) {
  if (path->plusp == true) {
    if (first_read_p == false) {
      return Univcoordlist_last_value(path->univdiagonals);
    } else if (Univcoordlist_head(path->univdiagonals) < (Univcoord_T) querylength) {
      return 0;
    } else {
      return Univcoordlist_head(path->univdiagonals) - querylength;
    }
  } else {
    if (first_read_p == true) {
      return Univcoordlist_last_value(path->univdiagonals);
    } else if (Univcoordlist_head(path->univdiagonals) < (Univcoord_T) querylength) {
      return 0;
    } else {
      return Univcoordlist_head(path->univdiagonals) - querylength;
    }
  }
}

static List_T
append_overlap_mate_mirror_paths (int *ncopied, List_T dest, List_T src,
				  int mincoverage_filter, int querylength, bool first_read_p,
				  Compress_T query_compress_fwd, Compress_T query_compress_rev,
				  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
				  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
				  Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  List_T p;
  Path_T path, copy;

  for (p = src; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->completep == true && Path_coverage(path) >= mincoverage_filter &&
	path->score_within_trims >= 0 && path->querylength == querylength) {
      copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool);
      copy->query_compress_fwd = query_compress_fwd;
      copy->query_compress_rev = query_compress_rev;
      copy->main_univdiagonal = mate_main_univdiagonal(copy,querylength,first_read_p);
      if (copy->genomic_diff != NULL) {
	copy->genomic_diff = (char *) NULL;
      }
      dest = Hitlist_push(dest,hitlistpool,(void *) copy
			  hitlistpool_trace(__FILE__,__LINE__));
      *ncopied += 1;
    }
  }

  return dest;
}

static Pathstore_T
improve_bounded_mate_pathstore (Pathstore_T pathstore_mate,
				Univcoord_T *univdiagonal_mate_out,
				int *found_score, int *end_score,
				int *mate_qstart_out, int *mate_qend_out,
				int *mate_nmismatches_out, int *mate_found_score_out,
				Shortread_T queryseq, int querylength,
				Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				Pathstore_T anchor_pathstore,
				Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
				bool scan_from_low_p, bool plusp, bool first_read_p,
				Pathstorepool_T pathstorepool,
				Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
				Pathpool_T pathpool, Transcriptpool_T transcriptpool,
				Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {
  Pathstore_T replacement;

  if (pathstore_mate != NULL &&
      pathstore_has_local_mate_coverage_p(pathstore_mate,querylength) == true) {
    return pathstore_mate;
  } else if ((replacement =
	      bounded_mate_pathstore(univdiagonal_mate_out,found_score,end_score,
				     mate_qstart_out,mate_qend_out,
				     mate_nmismatches_out,mate_found_score_out,
				     queryseq,querylength,query_compress,
				     query_compress_fwd,query_compress_rev,anchor_pathstore,
				     low_univdiagonal,high_univdiagonal,scan_from_low_p,plusp,first_read_p,
				     paired_mate_rescue_gapped_p,
				     pathstorepool,intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,vectorpool,hitlistpool)) == NULL) {
    return pathstore_mate;
  } else {
    if (pathstore_mate != NULL) {
      Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
		     listpool,pathpool,transcriptpool,hitlistpool);
    }
    return replacement;
  }
}

static List_T
paths5_mates (int *found_score_paired, int *found_score_3,
	      int *end_score_paired, int *end_score_3,
	      List_T *unresolved_pathpairs, List_T pathpairs,
	      List_T *local_mate_singlepaths3,
		     
	      Shortread_T queryseq5, Shortread_T queryseq3,
	      char *queryuc_ptr_3, char *queryrc3, int querylength3,
	      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	      int *mismatch_positions_alloc_3,
	      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
	      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
	      int max_nblocks, int localdb_nmismatches_allowed,

	      T this5, T this3,
	      Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		  
	      int nmismatches_filter_5, int nmismatches_filter_3,
	      int mincoverage_filter_5, int mincoverage_filter_3,

	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
	      Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
	      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
	      Spliceendsgen_T spliceendsgen3, bool only_complete5_p,
	      bool retain_local_mate_singlepaths_p,
	      bool allow_ungapped_mate_p, bool reliable_anchors_only_p,
	      Mate_rescue_stats_T *rescue_stats) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo_mate, auxinfo_copy;
  Pathstore_T pathstore, pathstore_mate;
  bool complete_sense_p, complete_antisense_p;
  int qstart, qend;

  List_T chain, p;
  Mate_rescue_segment_policy_T segment_policy;

  int i, nsplices;
  int k_gplus3, k_gminus3, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  int ignore;
  int mate_qstart, mate_qend, mate_nmismatches, mate_found_score, ncopied;
  int mate_selected_count;
  unsigned int mate_minpure_specificity_score;
  bool mate_minpure_observed_p;
  bool segment_control_stop_p = false;
#ifdef GSNAP_PAIRED_SOLVE_STATS
  int before_pathpairs, after_pathpairs, before_unresolved, after_unresolved;
  Pathstore_T pathstore_mate_before_bounded;
#endif


  /* From solved paths, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered paths5_mates with 5' all plus %d minus %d; 3' exhaustive plus %d minus %d\n",
		 this5->all_nunivdiagonals_gplus,this5->all_nunivdiagonals_gminus,
		 this3->nexhaustive_gplus,this3->nexhaustive_gminus));

  /* First, find threshold for mates */
  k_gplus3 = 0;
  k_gminus3 = this3->nexhaustive_gminus - 1;

  /* Go upstream */
  for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
    MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_scanned);
    univdiagonal = this5->all_univdiagonals_gplus[i];
    while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
      k_gplus3++;
    }

    /* auxinfo = this5->all_auxinfo_gplus[i]; */
    pathstore = this5->all_pathstores_gplus[i];

    if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_complete_skip);
    } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_rule_reject);
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_best);
      /* (1,2) Look upstream for 3' end matches */
      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
      if ((pathstore->mate_count =
	   find_max_count_upstream(&pathstore->mate_bestk,this3->exhaustive_gplus,
				   this3->exhaustive_counts_gplus,this3->nexhaustive_gplus,
				   k_gplus3,high_univdiagonal)) > max_count_gplus) {
	max_count_gplus = pathstore->mate_count;
      }
      if (pathstore->mate_count > 0) {
	MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_with_mate_count);
      }
    }
  }


  debug11(printf("paths5_mates upstream threshold max_count_gplus %d\n",max_count_gplus));

  /* Go downstream */
  for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
    MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_scanned);
    univdiagonal = this5->all_univdiagonals_gminus[i];
    while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
      k_gminus3--;
    }

    /* auxinfo = this5->all_auxinfo_gminus[i]; */
    pathstore = this5->all_pathstores_gminus[i];

    if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_complete_skip);
    } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_rule_reject);
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_best);
      /* (3,4) Look downstream for 3' end matches */
      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
      if ((pathstore->mate_count =
	   find_max_count_downstream(&pathstore->mate_bestk,
				     this3->exhaustive_gminus,this3->exhaustive_counts_gminus,
				     k_gminus3,low_univdiagonal)) > max_count_gminus) {
	max_count_gminus = pathstore->mate_count;
      }
      if (pathstore->mate_count > 0) {
	MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_with_mate_count);
      }
    }
  }
  

  /* Find threshold */
  debug11(printf("paths5_mates got max_count_gplus %d, max_count_gminus %d\n",
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

  if (max_count_gplus >= max_count || allow_ungapped_mate_p == true) {
    /* Go upstream */
    for (i = 0; i < this5->all_nunivdiagonals_gplus; i++) {
      if (segment_control_stop_p == true) {
	break;
      }
      MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_scanned);
      mate_rescue_stats_checkpoint(rescue_stats,"paths5_mates","upstream",/*forcep*/false);
      univdiagonal = this5->all_univdiagonals_gplus[i];
      while (k_gplus3 < this3->nexhaustive_gplus && this3->exhaustive_gplus[k_gplus3] < univdiagonal) {
	k_gplus3++;
      }

      debug11(printf("Finding upstream mates for 5' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this5->all_auxinfo_gplus[i]; */
      pathstore = this5->all_pathstores_gplus[i];

      if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_complete_skip);
      } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_rule_reject);

      } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_best);
	/* (1-2) Look upstream for 3' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
	pathstore_mate = (Pathstore_T) NULL;

	if (pathstore->mate_count >= max_count) {
	  MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_threshold_pass);
	  mate_selected_count =
	    find_minpure_ranked_upstream(&k,&mate_minpure_specificity_score,&mate_minpure_observed_p,
					 this3->exhaustive_gplus,this3->exhaustive_counts_gplus,
					 this3->exhaustive_minpure_specificity_gplus,this3->nexhaustive_gplus,
					 k_gplus3,high_univdiagonal,max_count,paired_mate_minpure_floor);
	  pathstore->mate_count = mate_selected_count;
	  debug11(printf("paths5_mates upstream minpure selected k %d count %d specificity %u observed %d floor %u\n",
			 k,mate_selected_count,mate_minpure_specificity_score,mate_minpure_observed_p,
			 paired_mate_minpure_floor));

	  if (k < 0) {
	    /* Skip */
	    debug11(printf("paths5_mates upstream minpure found no support/specificity-passing mate\n"));
	  } else {
	    assert(k >= 0);
	    univdiagonal_mate = this3->exhaustive_gplus[k];
	    auxinfo_mate = auxinfo_from_exhaustive_elt(this3->exhaustive_eltarray_gplus,
						       this3->nexhaustive_gplus,k,auxinfopool);
	    if (auxinfo_mate == NULL) {
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,auxinfo_null);
	    } else
	    if (univdiagonal_mate + auxinfo_mate->qstart >= pathstore->chrhigh + querylength3) { /* Upstream */
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,mate_bounds_skip);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);
	    } else {
	    pathstore_mate = Pathstore_new(pathstorepool);
	    MATE_RESCUE_STATS_INC(rescue_stats,pathstore_mate_created);
	    pathstore_mate->chrnum = pathstore->chrnum;
	    pathstore_mate->chroffset = pathstore->chroffset;
	    pathstore_mate->chrhigh = pathstore->chrhigh;
	    debug11(printf("(1-2) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,
			   auxinfo_mate->qstart,auxinfo_mate->qend,
			   this3->exhaustive_counts_gplus[k]));

	    if (auxinfo_mate->method == SEGMENT1) {
	      MATE_RESCUE_STATS_INC(rescue_stats,segment_solve_calls);
	      complete_sense_p = complete_antisense_p = false;
	      auxinfo_copy = auxinfo_mate;
	      mate_rescue_segment_policy_init(&segment_policy,rescue_stats,
					      "mate-rescue",
					      "paths5_mates","upstream",
					      retain_local_mate_singlepaths_p);
	      solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
					 &(*found_score_3),&(*end_score_3),
					 univdiagonal_mate,&auxinfo_copy,pathstore_mate,
					 queryseq3,/*queryptr*/queryuc_ptr_3,queryuc_ptr_3,queryrc3,querylength3,
					 querylength3/EXON_LENGTH + 1,this3,knownsplicing,knownindels,
					 mismatch_positions_alloc_3,
					 novel_univdiags_alloc,novel_diagonals_alloc,
					 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
					 /*query_compress*/query3_compress_fwd,
					 query3_compress_fwd,query3_compress_rev,/*genestrand*/0,
					 auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen3,/*plusp*/true,/*first_read_p*/false,
					 /*lowp*/false,/*set_best_paths_p*/true,
					 /*force_pending_auxinfo_p*/false,&segment_policy);
	      if (segment_policy.cap_hit_p == true) {
		segment_control_stop_p = true;
	      }

	    } else {
	      assert(auxinfo_mate->chaindir == CHAIN_RIGHT);
	      for (nsplices = 0; nsplices <= auxinfo_mate->max_nsplices; nsplices++) {
		for (p = auxinfo_mate->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
		  chain = (List_T) List_head(p);
		  MATE_RESCUE_STATS_INC(rescue_stats,chain_solve_calls);
		  Path_solve_from_right_chain(&(*found_score_3),&(*end_score_3),
				  
					      &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					      &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,

					      chain,queryseq3,/*queryptr*/queryuc_ptr_3,querylength3,mismatch_positions_alloc_3,
					  
					      novel_univdiags_alloc,novel_diagonals_alloc,
					      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					      localdb_nmismatches_allowed,
					  
					      this3,knownsplicing,knownindels,
					      /*query_compress*/query3_compress_fwd,query3_compress_fwd,query3_compress_rev,
					      pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
					      /*plusp*/true,/*genestrand*/0,
					      /*paired_end_p*/true,/*first_read_p*/false,intlistpool,
					      uintlistpool,univcoordlistpool,listpool,pathpool,
					      transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen3,
					      /*method*/LOCAL_MATE,/*find_splices_p*/true);
		}
	      }
	    }

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	    }
	  }
	    
#ifdef USE_SINGLE_UNIVDIAGONAL
	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,
							    /*queryptr*/queryuc_ptr_3,querylength3,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query3_compress_fwd,
							    /*plusp*/true,/*genestrand*/0,
							    genomebits)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_fwd,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  if (qstart < qend) {
	    debug11(printf("(1) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	    auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	    pathstore_mate = Pathstore_new(pathstorepool);

	    Path_solve_from_univdiagonal(&(*found_score_3),&(*end_score_3),
					 
					 &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					 &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
					 
					 univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryuc_ptr_3,
					 /*query_compress*/query3_compress_fwd,
					 query3_compress_fwd,query3_compress_rev,
					 /*plusp*/true,querylength3,mismatch_positions_alloc_3,

					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,

					 this3,knownsplicing,knownindels,/*paired_end_p*/true,/*first_read_p*/false,
					 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
					 intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  }
#endif
	}

	if (allow_ungapped_mate_p == true) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  pathstore_mate_before_bounded = pathstore_mate;
#endif
	  pathstore_mate =
	    improve_bounded_mate_pathstore(pathstore_mate,
					   &univdiagonal_mate,&(*found_score_3),&(*end_score_3),
					   &mate_qstart,&mate_qend,&mate_nmismatches,&mate_found_score,
					   queryseq3,querylength3,
					   /*query_compress*/query3_compress_fwd,
					   query3_compress_fwd,query3_compress_rev,pathstore,
					   low_univdiagonal,high_univdiagonal,/*scan_from_low_p*/true,
					   /*plusp*/true,/*first_read_p*/false,pathstorepool,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
	  mate_rescue_stats_note_bounded_result(rescue_stats,pathstore_mate_before_bounded,pathstore_mate);
	}

	if (pathstore_mate != NULL) {
	  MATE_RESCUE_STATS_INC(rescue_stats,mate_solved);
	  debug(printf("(4) Calling make_pathpairs\n"));
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  before_pathpairs = List_length(pathpairs);
	  before_unresolved = List_length(*unresolved_pathpairs);
#endif
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete5_p,/*only_completeH_p*/false);
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  after_pathpairs = List_length(pathpairs);
	  after_unresolved = List_length(*unresolved_pathpairs);
	  mate_rescue_stats_add_pathpair_delta(rescue_stats,before_pathpairs,after_pathpairs,
					       before_unresolved,after_unresolved);
#endif
	  ncopied = 0;
	  *local_mate_singlepaths3 =
	    append_local_mate_singlepaths_policy(&ncopied,*local_mate_singlepaths3,pathstore_mate,
						 mincoverage_filter_3,retain_local_mate_singlepaths_p,
						 rescue_stats,intlistpool,univcoordlistpool,listpool,
						 pathpool,vectorpool,transcriptpool,hitlistpool);
	  mate_rescue_stats_checkpoint(rescue_stats,"paths5_mates","upstream",/*forcep*/false);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }


  if (max_count_gminus >= max_count || allow_ungapped_mate_p == true) {
    /* Go downstream */
    for (i = this5->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
      if (segment_control_stop_p == true) {
	break;
      }
      MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_scanned);
      mate_rescue_stats_checkpoint(rescue_stats,"paths5_mates","downstream",/*forcep*/false);
      univdiagonal = this5->all_univdiagonals_gminus[i];
      while (k_gminus3 >= 0 && this3->exhaustive_gminus[k_gminus3] > univdiagonal) {
	k_gminus3--;
      }

      debug11(printf("Finding downstream mates for 5' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this5->all_auxinfo_gminus[i]; */
      pathstore = this5->all_pathstores_gminus[i];

      if (only_complete5_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_complete_skip);
      } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_rule_reject);
      } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_best);
	/* (3-4) Look downstream for 3' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
	pathstore_mate = (Pathstore_T) NULL;

	if (pathstore->mate_count >= max_count) {
	  MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_threshold_pass);
	  mate_selected_count =
	    find_minpure_ranked_downstream(&k,&mate_minpure_specificity_score,&mate_minpure_observed_p,
					   this3->exhaustive_gminus,this3->exhaustive_counts_gminus,
					   this3->exhaustive_minpure_specificity_gminus,this3->nexhaustive_gminus,
					   k_gminus3,low_univdiagonal,max_count,paired_mate_minpure_floor);
	  pathstore->mate_count = mate_selected_count;
	  debug11(printf("paths5_mates downstream minpure selected k %d count %d specificity %u observed %d floor %u\n",
			 k,mate_selected_count,mate_minpure_specificity_score,mate_minpure_observed_p,
			 paired_mate_minpure_floor));

	  if (k < 0) {
	    /* Skip */
	    debug11(printf("paths5_mates downstream minpure found no support/specificity-passing mate\n"));
	  } else {
	    assert(k >= 0);
	    univdiagonal_mate = this3->exhaustive_gminus[k];
	    auxinfo_mate = auxinfo_from_exhaustive_elt(this3->exhaustive_eltarray_gminus,
						       this3->nexhaustive_gminus,k,auxinfopool);
	    if (auxinfo_mate == NULL) {
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,auxinfo_null);
	    } else
	    if (univdiagonal_mate + auxinfo_mate->qend < pathstore->chroffset + querylength3) { /* Downstream */
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,mate_bounds_skip);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);

	    } else {
	    pathstore_mate = Pathstore_new(pathstorepool);
	    MATE_RESCUE_STATS_INC(rescue_stats,pathstore_mate_created);
	    pathstore_mate->chrnum = pathstore->chrnum;
	    pathstore_mate->chroffset = pathstore->chroffset;
	    pathstore_mate->chrhigh = pathstore->chrhigh;
	    debug11(printf("(3-4) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,
			   auxinfo_mate->qstart,auxinfo_mate->qend,
			   this3->exhaustive_counts_gminus[k]));

	    if (auxinfo_mate->method == SEGMENT1) {
	      MATE_RESCUE_STATS_INC(rescue_stats,segment_solve_calls);
	      complete_sense_p = complete_antisense_p = false;
	      auxinfo_copy = auxinfo_mate;
	      mate_rescue_segment_policy_init(&segment_policy,rescue_stats,
					      "mate-rescue",
					      "paths5_mates","downstream",
					      retain_local_mate_singlepaths_p);
	      solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
					 &(*found_score_3),&(*end_score_3),
					 univdiagonal_mate,&auxinfo_copy,pathstore_mate,
					 queryseq3,/*queryptr*/queryrc3,queryuc_ptr_3,queryrc3,querylength3,
					 querylength3/EXON_LENGTH + 1,this3,knownsplicing,knownindels,
					 mismatch_positions_alloc_3,
					 novel_univdiags_alloc,novel_diagonals_alloc,
					 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
					 /*query_compress*/query3_compress_rev,
					 query3_compress_fwd,query3_compress_rev,/*genestrand*/0,
					 auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen3,/*plusp*/false,/*first_read_p*/false,
					 /*lowp*/true,/*set_best_paths_p*/true,
					 /*force_pending_auxinfo_p*/false,&segment_policy);
	      if (segment_policy.cap_hit_p == true) {
		segment_control_stop_p = true;
	      }

	    } else {
	      assert(auxinfo_mate->chaindir == CHAIN_LEFT);
	      for (nsplices = 0; nsplices <= auxinfo_mate->max_nsplices; nsplices++) {
		for (p = auxinfo_mate->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
		  chain = (List_T) List_head(p);
		  MATE_RESCUE_STATS_INC(rescue_stats,chain_solve_calls);
		  Path_solve_from_left_chain(&(*found_score_3),&(*end_score_3),

					     &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					     &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,

					     chain,queryseq3,/*queryptr*/queryrc3,querylength3,mismatch_positions_alloc_3,

					     novel_univdiags_alloc,novel_diagonals_alloc,
					     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,
					  
					     this3,knownsplicing,knownindels,
					     /*query_compress*/query3_compress_rev,query3_compress_fwd,query3_compress_rev,
					     pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
					     /*plusp*/false,/*genestrand*/0,
					     /*paired_end_p*/true,/*first_read_p*/false,intlistpool,
					     uintlistpool,univcoordlistpool,listpool,pathpool,
					     transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen3,
					     /*method*/LOCAL_MATE,/*find_splices_p*/true);
		}
	      }
	    }

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	    }
	  }

#ifdef USE_SINGLE_UNIVDIAGONAL
	} else if (localdb != NULL &&
		   querylength3 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,
							     /*queryptr*/queryrc3,querylength3,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query3_compress_rev,
							     /*plusp*/false,/*genestrand*/0,
							     genomebits)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					      univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query3_compress_rev,
					     univdiagonal_mate,querylength3,/*pos5*/0,/*pos3*/querylength3,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  if (qstart < qend) {
	    debug11(printf("(3) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	    auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	    pathstore_mate = Pathstore_new(pathstorepool);
	
	    Path_solve_from_univdiagonal(&(*found_score_3),&(*end_score_3),
				       
					 &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					 &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
					 univdiagonal_mate,auxinfo_mate,queryseq3,/*queryptr*/queryrc3,
					 /*query_compress*/query3_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 /*plusp*/false,querylength3,mismatch_positions_alloc_3,

					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
				       
					 this3,knownsplicing,knownindels,/*paired_end_p*/true,/*first_read_p*/false,
					 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
					 intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen3,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  }
#endif
	}

	if (allow_ungapped_mate_p == true) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  pathstore_mate_before_bounded = pathstore_mate;
#endif
	  pathstore_mate =
	    improve_bounded_mate_pathstore(pathstore_mate,
					   &univdiagonal_mate,&(*found_score_3),&(*end_score_3),
					   &mate_qstart,&mate_qend,&mate_nmismatches,&mate_found_score,
					   queryseq3,querylength3,
					   /*query_compress*/query3_compress_rev,
					   query3_compress_fwd,query3_compress_rev,pathstore,
					   low_univdiagonal,high_univdiagonal,/*scan_from_low_p*/false,
					   /*plusp*/false,/*first_read_p*/false,pathstorepool,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
	  mate_rescue_stats_note_bounded_result(rescue_stats,pathstore_mate_before_bounded,pathstore_mate);
	}
	  
	if (pathstore_mate != NULL) {
	  MATE_RESCUE_STATS_INC(rescue_stats,mate_solved);
	  debug(printf("(6) Calling make_pathpairs\n"));
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  before_pathpairs = List_length(pathpairs);
	  before_unresolved = List_length(*unresolved_pathpairs);
#endif
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete5_p);
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  after_pathpairs = List_length(pathpairs);
	  after_unresolved = List_length(*unresolved_pathpairs);
	  mate_rescue_stats_add_pathpair_delta(rescue_stats,before_pathpairs,after_pathpairs,
					       before_unresolved,after_unresolved);
#endif
	  ncopied = 0;
	  *local_mate_singlepaths3 =
	    append_local_mate_singlepaths_policy(&ncopied,*local_mate_singlepaths3,pathstore_mate,
						 mincoverage_filter_3,retain_local_mate_singlepaths_p,
						 rescue_stats,intlistpool,univcoordlistpool,listpool,
						 pathpool,vectorpool,transcriptpool,hitlistpool);
	  mate_rescue_stats_checkpoint(rescue_stats,"paths5_mates","downstream",/*forcep*/false);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  mate_rescue_stats_checkpoint(rescue_stats,"paths5_mates","done",/*forcep*/true);
  return pathpairs;
}


static List_T
paths3_mates (int *found_score_paired, int *found_score_5,
	      int *end_score_paired, int *end_score_5,
	      List_T *unresolved_pathpairs, List_T pathpairs,
	      List_T *local_mate_singlepaths5,

	      Shortread_T queryseq5, Shortread_T queryseq3,
	      char *queryuc_ptr_5, char *queryrc5, int querylength5,
	      Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	      int *mismatch_positions_alloc_5,
	      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
	      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
	      int max_nblocks, int localdb_nmismatches_allowed,

	      T this5, T this3,
	      Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		  
	      int nmismatches_filter_5, int nmismatches_filter_3,
	      int mincoverage_filter_5, int mincoverage_filter_3,

	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
	      Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
	      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	      Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
	      Spliceendsgen_T spliceendsgen5, bool only_complete3_p,
	      bool retain_local_mate_singlepaths_p,
	      bool allow_ungapped_mate_p, bool reliable_anchors_only_p,
	      Mate_rescue_stats_T *rescue_stats) {
  
  Univcoord_T low_univdiagonal, high_univdiagonal, univdiagonal_mate, univdiagonal;
  Auxinfo_T auxinfo_mate, auxinfo_copy;
  Pathstore_T pathstore, pathstore_mate;
  bool complete_sense_p, complete_antisense_p;
  int qstart, qend;

  List_T chain, p;
  Mate_rescue_segment_policy_T segment_policy;

  int i, nsplices;
  int k_gplus5, k_gminus5, k;
  int max_count, max_count_gplus = 0, max_count_gminus = 0;
  int ignore;
  int mate_qstart, mate_qend, mate_nmismatches, mate_found_score, ncopied;
  int mate_selected_count;
  unsigned int mate_minpure_specificity_score;
  bool mate_minpure_observed_p;
  bool segment_control_stop_p = false;
#ifdef GSNAP_PAIRED_SOLVE_STATS
  int before_pathpairs, after_pathpairs, before_unresolved, after_unresolved;
  Pathstore_T pathstore_mate_before_bounded;
#endif


  /* From solved paths, find local mates */

  /* Ideally would like to search from the distal end to allow for
     cases where only tails overlap, but then we would have no
     guarantee that high_univdiagonal - low_univdiagonal <
     concordance_distance, and could overflow the stream or mergeinfo
     buffers */

  debug11(printf("Entered paths3_mates with 3' all plus %d minus %d; 5' exhaustive plus %d minus %d\n",
		 this3->all_nunivdiagonals_gplus,this3->all_nunivdiagonals_gminus,
		 this5->nexhaustive_gplus,this5->nexhaustive_gminus));

  /* First, find threshold for mates */
  k_gplus5 = this5->nexhaustive_gplus - 1;
  k_gminus5 = 0;

  /* Go downstream */
  for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
    MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_scanned);
    univdiagonal = this3->all_univdiagonals_gplus[i];
    while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
      k_gplus5--;
    }

    /* auxinfo = this3->all_auxinfo_gplus[i]; */
    pathstore = this3->all_pathstores_gplus[i];

    if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_complete_skip);
    } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_rule_reject);
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_best);
      /* (5,6) Look downstream for 5' end matches */
      high_univdiagonal = univdiagonal;
      low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
      if ((pathstore->mate_count =
	   find_max_count_downstream(&pathstore->mate_bestk,
				     this5->exhaustive_gplus,this5->exhaustive_counts_gplus,
				     k_gplus5,low_univdiagonal)) > max_count_gplus) {
	max_count_gplus = pathstore->mate_count;
      }
      if (pathstore->mate_count > 0) {
	MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_with_mate_count);
      }
    }
  }

  debug11(printf("paths3_mates downstream threshold max_count_gplus %d\n",max_count_gplus));

  /* Go upstream */
  for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
    MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_scanned);
    univdiagonal = this3->all_univdiagonals_gminus[i];
    while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
      k_gminus5++;
    }

    /* auxinfo = this3->all_auxinfo_gminus[i]; */
    pathstore = this3->all_pathstores_gminus[i];

    if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_complete_skip);
    } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
      /* Skip */
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_rule_reject);
    } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
      MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_best);
      /* (7,8) Look upstream for 5' end matches */
      low_univdiagonal = univdiagonal;
      high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
      if ((pathstore->mate_count =
	   find_max_count_upstream(&pathstore->mate_bestk,
				   this5->exhaustive_gminus,this5->exhaustive_counts_gminus,this5->nexhaustive_gminus,
				   k_gminus5,high_univdiagonal)) > max_count_gminus) {
	max_count_gminus = pathstore->mate_count;
      }
      if (pathstore->mate_count > 0) {
	MATE_RESCUE_STATS_INC(rescue_stats,threshold_anchor_with_mate_count);
      }
    }
  }

  /* Find threshold */
  debug11(printf("paths3_mates got max_count_gplus %d, max_count_gminus %d\n",
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

  if (max_count_gplus >= max_count || allow_ungapped_mate_p == true) {
    /* Go downstream */
    for (i = this3->all_nunivdiagonals_gplus - 1; i >= 0; i--) {
      if (segment_control_stop_p == true) {
	break;
      }
      MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_scanned);
      mate_rescue_stats_checkpoint(rescue_stats,"paths3_mates","downstream",/*forcep*/false);
      univdiagonal = this3->all_univdiagonals_gplus[i];
      while (k_gplus5 >= 0 && this5->exhaustive_gplus[k_gplus5] > univdiagonal) {
	k_gplus5--;
      }

      debug11(printf("Finding downstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gplus[i]; */
      pathstore = this3->all_pathstores_gplus[i];
      if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_complete_skip);
      } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_rule_reject);
      } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_best);
	/* (5-6) Look downstream for 5' end matches */
	high_univdiagonal = univdiagonal;
	low_univdiagonal = subtract_bounded(high_univdiagonal,concordance_distance,pathstore->chroffset);
	pathstore_mate = (Pathstore_T) NULL;

	if (pathstore->mate_count >= max_count) {
	  MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_threshold_pass);
	  mate_selected_count =
	    find_minpure_ranked_downstream(&k,&mate_minpure_specificity_score,&mate_minpure_observed_p,
					   this5->exhaustive_gplus,this5->exhaustive_counts_gplus,
					   this5->exhaustive_minpure_specificity_gplus,this5->nexhaustive_gplus,
					   k_gplus5,low_univdiagonal,max_count,paired_mate_minpure_floor);
	  pathstore->mate_count = mate_selected_count;
	  debug11(printf("paths3_mates downstream minpure selected k %d count %d specificity %u observed %d floor %u\n",
			 k,mate_selected_count,mate_minpure_specificity_score,mate_minpure_observed_p,
			 paired_mate_minpure_floor));

	  if (k < 0) {
	    /* Skip */
	    debug11(printf("paths3_mates downstream minpure found no support/specificity-passing mate\n"));
	  } else {
	    assert(k >= 0);
	    univdiagonal_mate = this5->exhaustive_gplus[k];
	    auxinfo_mate = auxinfo_from_exhaustive_elt(this5->exhaustive_eltarray_gplus,
						       this5->nexhaustive_gplus,k,auxinfopool);
	    debug11(printf("paths3_mates plus rescue anchor %u mate index %d mate %u auxinfo %p\n",
			   univdiagonal,k,univdiagonal_mate,auxinfo_mate));
	    if (auxinfo_mate == NULL) {
	    pathstore_mate = (Pathstore_T) NULL;
	    MATE_RESCUE_STATS_INC(rescue_stats,auxinfo_null);
	    } else
	    if (univdiagonal_mate + auxinfo_mate->qend < pathstore->chroffset + querylength5) { /* Downstream */
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,mate_bounds_skip);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);
	    } else {
	    pathstore_mate = Pathstore_new(pathstorepool);
	    MATE_RESCUE_STATS_INC(rescue_stats,pathstore_mate_created);
	    pathstore_mate->chrnum = pathstore->chrnum;
	    pathstore_mate->chroffset = pathstore->chroffset;
	    pathstore_mate->chrhigh = pathstore->chrhigh;
	    debug11(printf("(5-6) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,
			   auxinfo_mate->qstart,auxinfo_mate->qend,
			   this5->exhaustive_counts_gplus[k]));

	    if (auxinfo_mate->method == SEGMENT1) {
	      MATE_RESCUE_STATS_INC(rescue_stats,segment_solve_calls);
	      complete_sense_p = complete_antisense_p = false;
	      auxinfo_copy = auxinfo_mate;
	      mate_rescue_segment_policy_init(&segment_policy,rescue_stats,
					      "mate-rescue",
					      "paths3_mates","downstream",
					      retain_local_mate_singlepaths_p);
	      solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
					 &(*found_score_5),&(*end_score_5),
					 univdiagonal_mate,&auxinfo_copy,pathstore_mate,
					 queryseq5,/*queryptr*/queryuc_ptr_5,queryuc_ptr_5,queryrc5,querylength5,
					 querylength5/EXON_LENGTH + 1,this5,knownsplicing,knownindels,
					 mismatch_positions_alloc_5,
					 novel_univdiags_alloc,novel_diagonals_alloc,
					 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
					 /*query_compress*/query5_compress_fwd,
					 query5_compress_fwd,query5_compress_rev,/*genestrand*/0,
					 auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen5,/*plusp*/true,/*first_read_p*/true,
					 /*lowp*/true,/*set_best_paths_p*/true,
					 /*force_pending_auxinfo_p*/false,&segment_policy);
	      if (segment_policy.cap_hit_p == true) {
		segment_control_stop_p = true;
	      }

	    } else {
	      assert(auxinfo_mate->chaindir == CHAIN_LEFT);
	      for (nsplices = 0; nsplices <= auxinfo_mate->max_nsplices; nsplices++) {
		for (p = auxinfo_mate->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
		  chain = (List_T) List_head(p);
		  MATE_RESCUE_STATS_INC(rescue_stats,chain_solve_calls);
		  Path_solve_from_left_chain(&(*found_score_5),&(*end_score_5),

					     &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					     &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,

					     chain,queryseq5,/*queryptr*/queryuc_ptr_5,querylength5,mismatch_positions_alloc_5,

					     novel_univdiags_alloc,novel_diagonals_alloc,
					     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,
					  
					     this5,knownsplicing,knownindels,
					     /*query_compress*/query5_compress_fwd,query5_compress_fwd,query5_compress_rev,
					     pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
					     /*plusp*/true,/*genestrand*/0,
					     /*paired_end_p*/true,/*first_read_p*/true,intlistpool,
					     uintlistpool,univcoordlistpool,listpool,pathpool,
					     transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen5,
					     /*method*/LOCAL_MATE,/*find_splices_p*/true);
		}
	      }
	    }

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	    }
	  }

#ifdef USE_SINGLE_UNIVDIAGONAL
	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE && 
		   (univdiagonal_mate = Localdb_get_one_high(novel_univdiags_alloc,localdb,
							     novel_diagonals_alloc,localdb_alloc,
							     /*queryptr*/queryuc_ptr_5,querylength5,
							     low_univdiagonal,high_univdiagonal,
							     /*query_compress*/query5_compress_fwd,
							     /*plusp*/true,/*genestrand*/0,
							     genomebits)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_fwd,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/true,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  if (qstart < qend) {
	    debug11(printf("(5) highest mate: %u %u %d..%d\n",univdiagonal_mate,high_univdiagonal - univdiagonal_mate,qstart,qend));
	    auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	    pathstore_mate = Pathstore_new(pathstorepool);
	
	    Path_solve_from_univdiagonal(&(*found_score_5),&(*end_score_5),
				       
					 &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					 &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
					 univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryuc_ptr_5,
					 /*query_compress*/query5_compress_fwd,
					 query5_compress_fwd,query5_compress_rev,
					 /*plusp*/true,querylength5,mismatch_positions_alloc_5,

					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
				       
					 this5,knownsplicing,knownindels,/*paired_end_p*/true,/*first_read_p*/true,
					 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
					 intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  }
#endif
	}

	if (allow_ungapped_mate_p == true) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  pathstore_mate_before_bounded = pathstore_mate;
#endif
	  pathstore_mate =
	    improve_bounded_mate_pathstore(pathstore_mate,
					   &univdiagonal_mate,&(*found_score_5),&(*end_score_5),
					   &mate_qstart,&mate_qend,&mate_nmismatches,&mate_found_score,
					   queryseq5,querylength5,
					   /*query_compress*/query5_compress_fwd,
					   query5_compress_fwd,query5_compress_rev,pathstore,
					   low_univdiagonal,high_univdiagonal,/*scan_from_low_p*/false,
					   /*plusp*/true,/*first_read_p*/true,pathstorepool,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
	  mate_rescue_stats_note_bounded_result(rescue_stats,pathstore_mate_before_bounded,pathstore_mate);
	}
	  
	if (pathstore_mate != NULL) {
	  MATE_RESCUE_STATS_INC(rescue_stats,mate_solved);
	  debug(printf("(8) Calling make_pathpairs\n"));
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  before_pathpairs = List_length(pathpairs);
	  before_unresolved = List_length(*unresolved_pathpairs);
#endif
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore_mate,/*pathstoreH*/pathstore,
				     /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,/*plusp*/true,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/false,/*only_completeH_p*/only_complete3_p);
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  after_pathpairs = List_length(pathpairs);
	  after_unresolved = List_length(*unresolved_pathpairs);
	  mate_rescue_stats_add_pathpair_delta(rescue_stats,before_pathpairs,after_pathpairs,
					       before_unresolved,after_unresolved);
#endif
	  ncopied = 0;
	  *local_mate_singlepaths5 =
	    append_local_mate_singlepaths_policy(&ncopied,*local_mate_singlepaths5,pathstore_mate,
						 mincoverage_filter_5,retain_local_mate_singlepaths_p,
						 rescue_stats,intlistpool,univcoordlistpool,listpool,
						 pathpool,vectorpool,transcriptpool,hitlistpool);
	  mate_rescue_stats_checkpoint(rescue_stats,"paths3_mates","downstream",/*forcep*/false);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }


  if (max_count_gminus >= max_count || allow_ungapped_mate_p == true) {
    /* Go upstream */
    for (i = 0; i < this3->all_nunivdiagonals_gminus; i++) {
      if (segment_control_stop_p == true) {
	break;
      }
      MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_scanned);
      mate_rescue_stats_checkpoint(rescue_stats,"paths3_mates","upstream",/*forcep*/false);
      univdiagonal = this3->all_univdiagonals_gminus[i];
      while (k_gminus5 < this5->nexhaustive_gminus && this5->exhaustive_gminus[k_gminus5] < univdiagonal) {
	k_gminus5++;
      }
      
      debug11(printf("Finding upstream mates for 3' univdiagonal %u\n",univdiagonal));
      /* auxinfo = this3->all_auxinfo_gminus[i]; */
      pathstore = this3->all_pathstores_gminus[i];
      if (only_complete3_p == true && pathstore->complete_sense_p == false && pathstore->complete_antisense_p == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_complete_skip);
      } else if (pathstore_mate_rescue_anchor_p(pathstore,reliable_anchors_only_p) == false) {
	/* Skip */
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_rule_reject);
      } else if (pathstore->best_sense_paths != NULL || pathstore->best_antisense_paths != NULL) {
	MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_best);
	/* (7-8) Look upstream for 5' end matches */
	low_univdiagonal = univdiagonal;
	high_univdiagonal = add_bounded(low_univdiagonal,concordance_distance,pathstore->chrhigh);
	pathstore_mate = (Pathstore_T) NULL;

	if (pathstore->mate_count >= max_count) {
	  MATE_RESCUE_STATS_INC(rescue_stats,rescue_anchor_threshold_pass);
	  mate_selected_count =
	    find_minpure_ranked_upstream(&k,&mate_minpure_specificity_score,&mate_minpure_observed_p,
					 this5->exhaustive_gminus,this5->exhaustive_counts_gminus,
					 this5->exhaustive_minpure_specificity_gminus,this5->nexhaustive_gminus,
					 k_gminus5,high_univdiagonal,max_count,paired_mate_minpure_floor);
	  pathstore->mate_count = mate_selected_count;
	  debug11(printf("paths3_mates upstream minpure selected k %d count %d specificity %u observed %d floor %u\n",
			 k,mate_selected_count,mate_minpure_specificity_score,mate_minpure_observed_p,
			 paired_mate_minpure_floor));

	  if (k < 0) {
	    /* Skip */
	    debug11(printf("paths3_mates upstream minpure found no support/specificity-passing mate\n"));
	  } else {
	    assert(k >= 0);
	    univdiagonal_mate = this5->exhaustive_gminus[k];
	    auxinfo_mate = auxinfo_from_exhaustive_elt(this5->exhaustive_eltarray_gminus,
						       this5->nexhaustive_gminus,k,auxinfopool);
	    debug11(printf("paths3_mates minus rescue anchor %u mate index %d mate %u auxinfo %p\n",
			   univdiagonal,k,univdiagonal_mate,auxinfo_mate));
	    if (auxinfo_mate == NULL) {
	    pathstore_mate = (Pathstore_T) NULL;
	    MATE_RESCUE_STATS_INC(rescue_stats,auxinfo_null);
	    } else
	    if (univdiagonal_mate + auxinfo_mate->qstart >= pathstore->chrhigh + querylength5) { /* Upstream */
	    /* Skip */
	    MATE_RESCUE_STATS_INC(rescue_stats,mate_bounds_skip);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);
	    } else {
	    pathstore_mate = Pathstore_new(pathstorepool);
	    MATE_RESCUE_STATS_INC(rescue_stats,pathstore_mate_created);
	    pathstore_mate->chrnum = pathstore->chrnum;
	    pathstore_mate->chroffset = pathstore->chroffset;
	    pathstore_mate->chrhigh = pathstore->chrhigh;
	    debug11(printf("(7-8) prevalent mate: %u %d..%d %d\n",univdiagonal_mate,
			   auxinfo_mate->qstart,auxinfo_mate->qend,
			   this5->exhaustive_counts_gminus[k]));

	    if (auxinfo_mate->method == SEGMENT1) {
	      MATE_RESCUE_STATS_INC(rescue_stats,segment_solve_calls);
	      complete_sense_p = complete_antisense_p = false;
	      auxinfo_copy = auxinfo_mate;
	      mate_rescue_segment_policy_init(&segment_policy,rescue_stats,
					      "mate-rescue",
					      "paths3_mates","upstream",
					      retain_local_mate_singlepaths_p);
	      solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
					 &(*found_score_5),&(*end_score_5),
					 univdiagonal_mate,&auxinfo_copy,pathstore_mate,
					 queryseq5,/*queryptr*/queryrc5,queryuc_ptr_5,queryrc5,querylength5,
					 querylength5/EXON_LENGTH + 1,this5,knownsplicing,knownindels,
					 mismatch_positions_alloc_5,
					 novel_univdiags_alloc,novel_diagonals_alloc,
					 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
					 /*query_compress*/query5_compress_rev,
					 query5_compress_fwd,query5_compress_rev,/*genestrand*/0,
					 auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen5,/*plusp*/false,/*first_read_p*/true,
					 /*lowp*/false,/*set_best_paths_p*/true,
					 /*force_pending_auxinfo_p*/false,&segment_policy);
	      if (segment_policy.cap_hit_p == true) {
		segment_control_stop_p = true;
	      }

	    } else {
	      assert(auxinfo_mate->chaindir == CHAIN_RIGHT);
	      for (nsplices = 0; nsplices <= auxinfo_mate->max_nsplices; nsplices++) {
		for (p = auxinfo_mate->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
		  chain = (List_T) List_head(p);
		  MATE_RESCUE_STATS_INC(rescue_stats,chain_solve_calls);
		  Path_solve_from_right_chain(&(*found_score_5),&(*end_score_5),

					      &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					      &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
					  
					      chain,queryseq5,/*queryptr*/queryrc5,querylength5,mismatch_positions_alloc_5,
					  
					      novel_univdiags_alloc,novel_diagonals_alloc,
					      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					      localdb_nmismatches_allowed,
					  
					      this5,knownsplicing,knownindels,
					      /*query_compress*/query5_compress_rev,query5_compress_fwd,query5_compress_rev,
					      pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,
					      /*plusp*/false,/*genestrand*/0,
					      /*paired_end_p*/true,/*first_read_p*/true,intlistpool,
					      uintlistpool,univcoordlistpool,listpool,pathpool,
					      transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen5,
					      /*method*/LOCAL_MATE,/*find_splices_p*/true);
		}
	      }
	    }

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	    }
	  }

#ifdef USE_SINGLE_UNIVDIAGONAL
	} else if (localdb != NULL &&
		   querylength5 < QUERYLENGTH_FOR_LOCALDB_MATE &&
		   (univdiagonal_mate = Localdb_get_one_low(novel_univdiags_alloc,localdb,
							    novel_diagonals_alloc,localdb_alloc,
							    /*queryptr*/queryrc5,querylength5,
							    low_univdiagonal,high_univdiagonal,
							    /*query_compress*/query5_compress_rev,
							    /*plusp*/false,/*genestrand*/0,
							    genomebits)) != 0) {
	  qstart = Genomebits_first_kmer_left(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					      univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					      /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  qend = Genomebits_first_kmer_right(&ignore,genomebits,/*query_compress*/query5_compress_rev,
					     univdiagonal_mate,querylength5,/*pos5*/0,/*pos3*/querylength5,
					     /*plusp*/false,/*genestrand*/0,/*query_unk_mismatch_p*/true,/*kmer*/4);
	  if (qstart < qend) {
	    debug11(printf("(7) lowest mate: %u %u %d..%d\n",univdiagonal_mate,univdiagonal_mate - low_univdiagonal,qstart,qend));
	    auxinfo_mate = Auxinfo_new(LOCAL_MATE,qstart,qend,auxinfopool);
	    pathstore_mate = Pathstore_new(pathstorepool);

	    Path_solve_from_univdiagonal(&(*found_score_5),&(*end_score_5),
				       
					 &pathstore_mate->unextended_sense_paths,&pathstore_mate->unextended_antisense_paths,
					 &pathstore_mate->complete_sense_paths,&pathstore_mate->complete_antisense_paths,
				       
					 univdiagonal_mate,auxinfo_mate,queryseq5,/*queryptr*/queryrc5,
					 /*query_compress*/query5_compress_rev,
					 query5_compress_fwd,query5_compress_rev,
					 /*plusp*/false,querylength5,mismatch_positions_alloc_5,

					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,
					 localdb_nmismatches_allowed,
				       
					 this5,knownsplicing,knownindels,/*paired_end_p*/true,/*first_read_p*/true,
					 /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,/*chrhigh*/pathstore->chrhigh,
					 intlistpool,uintlistpool,univcoordlistpool,
					 listpool,pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
					 spliceendsgen5,/*method*/LOCAL_MATE,/*find_splices_p*/true);
	    Auxinfo_free(&auxinfo_mate,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);

	    pathstore_mate->solvedp = true;
	    Pathstore_set_best_paths(pathstore_mate,hitlistpool); /* Needed for make_pathpairs to work */
	  }
#endif
	}

	if (allow_ungapped_mate_p == true) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  pathstore_mate_before_bounded = pathstore_mate;
#endif
	  pathstore_mate =
	    improve_bounded_mate_pathstore(pathstore_mate,
					   &univdiagonal_mate,&(*found_score_5),&(*end_score_5),
					   &mate_qstart,&mate_qend,&mate_nmismatches,&mate_found_score,
					   queryseq5,querylength5,
					   /*query_compress*/query5_compress_rev,
					   query5_compress_fwd,query5_compress_rev,pathstore,
					   low_univdiagonal,high_univdiagonal,/*scan_from_low_p*/true,
					   /*plusp*/false,/*first_read_p*/true,pathstorepool,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
	  mate_rescue_stats_note_bounded_result(rescue_stats,pathstore_mate_before_bounded,pathstore_mate);
	}

	if (pathstore_mate != NULL) {
	  MATE_RESCUE_STATS_INC(rescue_stats,mate_solved);
	  debug(printf("(10) Calling make_pathpairs\n"));
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  before_pathpairs = List_length(pathpairs);
	  before_unresolved = List_length(*unresolved_pathpairs);
#endif
	  pathpairs = make_pathpairs(&(*found_score_paired),&(*end_score_paired),&(*unresolved_pathpairs),pathpairs,
				     /*pathstoreL*/pathstore,/*pathstoreH*/pathstore_mate,
				     /*queryseqL*/queryseq3,/*queryseqH*/queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
				     transcriptpool,hitlistpool,
				     /*only_completeL_p*/only_complete3_p,/*only_completeH_p*/false);
#ifdef GSNAP_PAIRED_SOLVE_STATS
	  after_pathpairs = List_length(pathpairs);
	  after_unresolved = List_length(*unresolved_pathpairs);
	  mate_rescue_stats_add_pathpair_delta(rescue_stats,before_pathpairs,after_pathpairs,
					       before_unresolved,after_unresolved);
#endif
	  ncopied = 0;
	  *local_mate_singlepaths5 =
	    append_local_mate_singlepaths_policy(&ncopied,*local_mate_singlepaths5,pathstore_mate,
						 mincoverage_filter_5,retain_local_mate_singlepaths_p,
						 rescue_stats,intlistpool,univcoordlistpool,listpool,
						 pathpool,vectorpool,transcriptpool,hitlistpool);
	  mate_rescue_stats_checkpoint(rescue_stats,"paths3_mates","upstream",/*forcep*/false);
	  Pathstore_free(&pathstore_mate,pathstorepool,intlistpool,univcoordlistpool,
			 listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  mate_rescue_stats_checkpoint(rescue_stats,"paths3_mates","done",/*forcep*/true);
  return pathpairs;
}


static List_T
paired_local_mate_rescue (const char *stage,
			  int *found_score_paired, int *found_score_5, int *found_score_3,
			  int *end_score_paired, int *end_score_5, int *end_score_3,
			  List_T *unresolved_pathpairs, List_T pathpairs,
			  List_T *local_mate_singlepaths5, List_T *local_mate_singlepaths3,

			  Shortread_T queryseq5, Shortread_T queryseq3,
			  char *queryuc_ptr_5, char *queryrc5, int querylength5,
			  char *queryuc_ptr_3, char *queryrc3, int querylength3,
			  Knownsplicing_T knownsplicing, Knownindels_T knownindels,

			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			  int max_nblocks, int localdb_nmismatches_allowed,

			  T this5, T this3,
			  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
			  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

			  int nmismatches_filter_5, int nmismatches_filter_3,
			  int mincoverage_filter_5, int mincoverage_filter_3,

			  Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			  Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			  Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			  Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3,
			  bool complete5_p, bool complete3_p,
			  bool retain_local_mate_singlepaths_p,
			  bool allow_ungapped_mate_p, bool reliable_anchors_only_p) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
  double rescue_t0;
  int before_pathpairs, before_unresolved;
  Mate_rescue_stats_T rescue_stats_body, *rescue_stats = &rescue_stats_body;
#else
  Mate_rescue_stats_T *rescue_stats = (Mate_rescue_stats_T *) NULL;
#endif

  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE) == false) {
    return pathpairs;
  }

#ifdef GSNAP_PAIRED_SOLVE_STATS
  mate_rescue_stats_init(rescue_stats,stage);
  rescue_t0 = paired_solve_stats_now();
  before_pathpairs = List_length(pathpairs);
  before_unresolved = List_length(*unresolved_pathpairs);
  fprintf(stderr,
	  "paired_solve_stats\tmate_rescue\tstage=%s\tbegin_pathpairs=%d\tbegin_unresolved=%d\tcomplete5=%d\tcomplete3=%d\tretain_local_mate_singlepaths=%d\tallow_ungapped=%d\treliable_anchors_only=%d\n",
	  stage,before_pathpairs,before_unresolved,complete5_p,complete3_p,
	  retain_local_mate_singlepaths_p,allow_ungapped_mate_p,reliable_anchors_only_p);
  paired_solve_stats_print_stage("mate_rescue_5_begin",this5);
  paired_solve_stats_print_stage("mate_rescue_3_begin",this3);
#endif

  debug11(printf("%s local mate rescue: 5' exhaustive plus %d, minus %d; 3' exhaustive plus %d, minus %d\n",
		 stage,this5->nexhaustive_gplus,this5->nexhaustive_gminus,
		 this3->nexhaustive_gplus,this3->nexhaustive_gminus));

  pathpairs = paths5_mates(found_score_paired,found_score_3,end_score_paired,end_score_3,
			   unresolved_pathpairs,pathpairs,local_mate_singlepaths3,
			   queryseq5,queryseq3,
			   queryuc_ptr_3,queryrc3,querylength3,
			   knownsplicing,knownindels,this3->mismatch_positions_alloc,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
			   this5,this3,query3_compress_fwd,query3_compress_rev,
			   nmismatches_filter_5,nmismatches_filter_3,
			   mincoverage_filter_5,mincoverage_filter_3,
			   auxinfopool,pathstorepool,univdiagpool,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			   spliceendsgen3,/*only_complete5_p*/complete5_p,
			   retain_local_mate_singlepaths_p,allow_ungapped_mate_p,
			   reliable_anchors_only_p,rescue_stats);

  pathpairs = paths3_mates(found_score_paired,found_score_5,end_score_paired,end_score_5,
			   unresolved_pathpairs,pathpairs,local_mate_singlepaths5,
			   queryseq5,queryseq3,
			   queryuc_ptr_5,queryrc5,querylength5,
			   knownsplicing,knownindels,this5->mismatch_positions_alloc,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
			   this5,this3,query5_compress_fwd,query5_compress_rev,
			   nmismatches_filter_5,nmismatches_filter_3,
			   mincoverage_filter_5,mincoverage_filter_3,
			   auxinfopool,pathstorepool,univdiagpool,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			   spliceendsgen5,/*only_complete3_p*/complete3_p,
			   retain_local_mate_singlepaths_p,allow_ungapped_mate_p,
			   reliable_anchors_only_p,rescue_stats);

#ifdef GSNAP_PAIRED_SOLVE_STATS
  mate_rescue_stats_checkpoint(rescue_stats,"paired_local_mate_rescue","end",/*forcep*/true);
  fprintf(stderr,
	  "paired_solve_stats\tmate_rescue\tstage=%s\telapsed=%.6f\tend_pathpairs=%d\tend_unresolved=%d\tdelta_pathpairs=%d\tdelta_unresolved=%d\n",
	  stage,paired_solve_stats_now() - rescue_t0,List_length(pathpairs),List_length(*unresolved_pathpairs),
	  List_length(pathpairs) - before_pathpairs,List_length(*unresolved_pathpairs) - before_unresolved);
  paired_solve_stats_print_pathpair_list(stage,pathpairs);
#endif

  return pathpairs;
}

static bool
pathpairs_have_resolved_best_pair_p (List_T pathpairs, int found_score_paired) {
  Pathpair_T pathpair;
  List_T p;

  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (Pathpair_T) List_head(p);
    if (Pathpair_found_score(pathpair) == found_score_paired &&
	Pathpair_end_score(pathpair) == 0) {
      return true;
    }
  }

  return false;
}

static void
solve_all (bool *completep, int *found_score, int *end_score, Shortread_T queryseq,
	   char *queryuc_ptr, char *queryrc, int querylength, int max_nsplices,
	   T this, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  
	   int *mismatch_positions_alloc,
	   Univdiag_T *novel_univdiags_alloc,
	   Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
	   int localdb_nmismatches_allowed, 
	   
	   Compress_T query_compress_fwd, Compress_T query_compress_rev,
	   bool first_read_p,

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
#ifdef GSNAP_PAIRED_SOLVE_STATS
  double solve_t0, plus_t1, minus_t1;
  const char *mate_label = first_read_p == true ? "5" : "3";

  solve_t0 = paired_solve_stats_now();
  paired_solve_stats_print_stage(first_read_p == true ? "solve_all_5_begin" : "solve_all_3_begin",this);
#endif

  /* plus */
  lowp = (first_read_p == true) ? true : false;
  for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
    univdiagonal = this->all_univdiagonals_gplus[i];
    pathstore = this->all_pathstores_gplus[i];
    solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
			       &(*found_score),&(*end_score),
			       univdiagonal,&(this->all_auxinfo_gplus[i]),pathstore,

			       queryseq,/*queryptr*/queryuc_ptr,queryuc_ptr,queryrc,querylength,
			       max_nsplices,this,knownsplicing,knownindels,
				 
			       this->mismatch_positions_alloc,
			       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			       localdb_tally,localdb_overflowp,max_nblocks,
			       localdb_nmismatches_allowed,

			       /*query_compress*/query_compress_fwd,
			       query_compress_fwd,query_compress_rev,/*genestrand*/0,
			       auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
			       hitlistpool,spliceendsgen,/*plusp*/true,first_read_p,lowp,
			       /*set_best_paths_p*/false,/*force_pending_auxinfo_p*/false,
			       /*mate_rescue_segment_policy*/NULL);
    if (complete_sense_p == true || complete_antisense_p == true) {
      *completep = true;
    }
  }
#ifdef GSNAP_PAIRED_SOLVE_STATS
  plus_t1 = paired_solve_stats_now();
  fprintf(stderr,"paired_solve_stats\tsolve_all\tmate=%s\tstrand=plus\telapsed=%.6f\tcomplete_sense=%d\tcomplete_antisense=%d\tfound=%d\tend=%d\n",
	  mate_label,plus_t1 - solve_t0,complete_sense_p,complete_antisense_p,*found_score,*end_score);
#endif

  /* minus */
  lowp = (first_read_p == true) ? false : true;
  for (i = this->all_nunivdiagonals_gminus - 1; i >= 0; i--) {
    univdiagonal = this->all_univdiagonals_gminus[i];
    pathstore = this->all_pathstores_gminus[i];
    solve_univdiagonal_auxinfo(&complete_sense_p,&complete_antisense_p,
			       &(*found_score),&(*end_score),
			       univdiagonal,&(this->all_auxinfo_gminus[i]),pathstore,

			       queryseq,/*queryptr*/queryrc,queryuc_ptr,queryrc,querylength,
			       max_nsplices,this,knownsplicing,knownindels,

			       this->mismatch_positions_alloc,
			       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			       localdb_tally,localdb_overflowp,max_nblocks,
			       localdb_nmismatches_allowed,

			       /*query_compress*/query_compress_rev,
			       query_compress_fwd,query_compress_rev,/*genestrand*/0,
			       auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
			       hitlistpool,spliceendsgen,/*plusp*/false,first_read_p,lowp,
			       /*set_best_paths_p*/false,/*force_pending_auxinfo_p*/false,
			       /*mate_rescue_segment_policy*/NULL);
  }
#ifdef GSNAP_PAIRED_SOLVE_STATS
  minus_t1 = paired_solve_stats_now();
  fprintf(stderr,"paired_solve_stats\tsolve_all\tmate=%s\tstrand=minus\telapsed=%.6f\tcomplete_sense=%d\tcomplete_antisense=%d\tfound=%d\tend=%d\n",
	  mate_label,minus_t1 - plus_t1,complete_sense_p,complete_antisense_p,*found_score,*end_score);
#endif

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

#ifdef GSNAP_PAIRED_SOLVE_STATS
  paired_solve_stats_print_stage(first_read_p == true ? "solve_all_5_after_best" : "solve_all_3_after_best",this);
  paired_solve_stats_print_pathstores(first_read_p == true ? "solve_all_5_plus" : "solve_all_3_plus",
				      this->all_pathstores_gplus,this->all_nunivdiagonals_gplus);
  paired_solve_stats_print_pathstores(first_read_p == true ? "solve_all_5_minus" : "solve_all_3_minus",
				      this->all_pathstores_gminus,this->all_nunivdiagonals_gminus);
  fprintf(stderr,"paired_solve_stats\tsolve_all\tmate=%s\ttotal_elapsed=%.6f\tcomplete=%d\tfound=%d\tend=%d\n",
	  mate_label,paired_solve_stats_now() - solve_t0,*completep,*found_score,*end_score);
#endif

  return;
}

static void
collect_single_end_paths (List_T *sense_paths5, List_T *antisense_paths5,
			  List_T *sense_paths3, List_T *antisense_paths3,
			  T this5, T this3, Hitlistpool_T hitlistpool) {
  List_T sense_paths5_gplus = NULL, sense_paths5_gminus = NULL,
    antisense_paths5_gplus = NULL, antisense_paths5_gminus = NULL,
    sense_paths3_gplus = NULL, sense_paths3_gminus = NULL,
    antisense_paths3_gplus = NULL, antisense_paths3_gminus = NULL;

  Stage1_collect_paths(&sense_paths5_gplus,&sense_paths5_gminus,
		       &antisense_paths5_gplus,&antisense_paths5_gminus,
		       this5,hitlistpool);
  Stage1_collect_paths(&sense_paths3_gplus,&sense_paths3_gminus,
		       &antisense_paths3_gplus,&antisense_paths3_gminus,
		       this3,hitlistpool);

  *sense_paths5 = List_append(sense_paths5_gplus,sense_paths5_gminus);
  *antisense_paths5 = List_append(antisense_paths5_gplus,antisense_paths5_gminus);
  *sense_paths3 = List_append(sense_paths3_gplus,sense_paths3_gminus);
  *antisense_paths3 = List_append(antisense_paths3_gplus,antisense_paths3_gminus);

#ifdef GSNAP_PAIRED_SOLVE_STATS
  fprintf(stderr,
	  "paired_solve_stats\tcollect_single_end_paths\tsense5=%d\tantisense5=%d\tsense3=%d\tantisense3=%d\n",
	  List_length(*sense_paths5),List_length(*antisense_paths5),
	  List_length(*sense_paths3),List_length(*antisense_paths3));
  paired_solve_stats_print_path_list("fallback_sense5",*sense_paths5);
  paired_solve_stats_print_path_list("fallback_antisense5",*antisense_paths5);
  paired_solve_stats_print_path_list("fallback_sense3",*sense_paths3);
  paired_solve_stats_print_path_list("fallback_antisense3",*antisense_paths3);
#endif

  return;
}

static int
stage1_all_univdiagonal_count (T this) {
  return this->all_nunivdiagonals_gplus + this->all_nunivdiagonals_gminus;
}

static bool
pathpair_uses_local_mate_p (Pathpair_T pathpair) {
  return pathpair->path5->method == LOCAL_MATE || pathpair->path3->method == LOCAL_MATE;
}

static List_T
keep_local_mate_pathpairs (List_T pathpairs,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Hitlistpool_T hitlistpool) {
  List_T p, kept = NULL, rejected = NULL;
  Pathpair_T pathpair;
#ifdef GSNAP_PAIRED_SOLVE_STATS
  int nkept = 0, nrejected = 0;
#endif

  for (p = pathpairs; p != NULL; p = List_next(p)) {
    pathpair = (Pathpair_T) List_head(p);
    if (pathpair_uses_local_mate_p(pathpair) == true) {
      kept = Hitlist_push(kept,hitlistpool,(void *) pathpair
			  hitlistpool_trace(__FILE__,__LINE__));
#ifdef GSNAP_PAIRED_SOLVE_STATS
      nkept++;
#endif
    } else {
      rejected = Hitlist_push(rejected,hitlistpool,(void *) pathpair
			      hitlistpool_trace(__FILE__,__LINE__));
#ifdef GSNAP_PAIRED_SOLVE_STATS
      nrejected++;
#endif
    }
  }

  Hitlistpool_free_list(&pathpairs,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Pathpair_gc(&rejected,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

#ifdef GSNAP_PAIRED_SOLVE_STATS
  fprintf(stderr,"paired_solve_stats\tanchor_first_filter\tkept_local_mate=%d\trejected_non_mate=%d\n",
	  nkept,nrejected);
#endif

  return List_reverse(kept);
}


/* final_pairtype can be CONCORDANT_TRANSLOCATIONS, CONCORDANT, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG, UNPAIRED */
Pathpair_T *
Stage1_paired_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, Pairtype_T *final_pairtype,
		    Path_T **patharray5, int *npaths5_primary, int *npaths5_altloc, int *first_absmq5, int *second_absmq5,
		    Path_T **patharray3, int *npaths3_primary, int *npaths3_altloc, int *first_absmq3, int *second_absmq3,
		    T this5, T this3, Shortread_T queryseq5, Shortread_T queryseq3, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Chrpos_T pairmax_linear,

		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		    int localdb_nmismatches_allowed,

		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		    Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3,
		    Cascade_Timing_T *cascade_timing, Pass_T pass) {
  Pathpair_T *pathpairarray = NULL;
  List_T pathpairs = NULL, unresolved_pathpairs = NULL, unresolved_pathpairs_from_exhaustive = NULL;
  List_T terminal_anchor_pathpairs_deferred = NULL;
  List_T terminal_anchor_paths5_deferred = NULL;
  List_T terminal_anchor_paths3_deferred = NULL;
  Method_T last_method_5, last_method_3, method_goal;

  List_T sense_paths5 = NULL, antisense_paths5 = NULL, sense_paths3 = NULL, antisense_paths3 = NULL;
  List_T rejected_singlepaths5 = NULL, rejected_singlepaths3 = NULL;
  List_T local_mate_singlepaths5 = NULL, local_mate_singlepaths3 = NULL;

  /* List_T local_sense_paths5 = NULL, local_antisense_paths5 = NULL,
     local_sense_paths3 = NULL, local_antisense_paths3 = NULL; */

  /* List_T geneplus_pathpairs, geneminus_pathpairs; */
  /* List_T geneplus_paths5, geneplus_paths3, geneminus_paths5, geneminus_paths3; */

  /* Maximum number of localdb regions possible, based on region size of 65536 bp */
  int max_localdb_nregions = (positive_gap_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  /* int max_nintersections = max_localdb_nregions * LOCALDB_REGION_SIZE; */
  int sufficient_score_5, sufficient_score_3, max_nsplices_5, max_nsplices_3, querylength5, querylength3;

  bool complete5_p, complete3_p;
  int found_score_paired, found_score_5, found_score_3;
  int end_score_paired, end_score_5, end_score_3;
  bool any_imperfect_ends5_p = false, any_imperfect_ends3_p = false;
  bool solved_all_p = false;
  bool solved_all_5_p = false, solved_all_3_p = false;
  bool anchor_first_solve_attempted_p = false;
  bool post_solve_fallback_entered_p = false;
  bool exact_full_overlap_p = false;

  int nmismatches_filter_5, nmismatches_filter_3;
  int mincoverage_filter_5, mincoverage_filter_3;
  int mincoverage_querylength_5, mincoverage_querylength_3;
  int nmismatches_allowed_5, nmismatches_allowed_3;
  int i;
  Path_T path5, path3;
  char *queryuc_ptr_5, *queryuc_ptr_3, *queryrc5, *queryrc3;
  Compress_T query5_compress_fwd, query5_compress_rev, query3_compress_fwd, query3_compress_rev;
  unsigned long long stage_start_us = CASCADE_STAGE_START("setup",queryseq5,univdiagpool);


  *npaths_primary = *npaths_altloc = 0;
  *first_absmq = *second_absmq = 0;
  *final_pairtype = UNPAIRED;

  *patharray5 = *patharray3 = (Path_T *) NULL;
  *npaths5_primary = *npaths5_altloc = 0;
  *first_absmq5 = *second_absmq5 = 0;
  *npaths3_primary = *npaths3_altloc = 0;
  *first_absmq3 = *second_absmq3 = 0;

  if ((querylength5 = Shortread_fulllength(queryseq5)) < index1part + index1interval - 1 ||
      (querylength3 = Shortread_fulllength(queryseq3)) < index1part + index1interval - 1) {
    CASCADE_STAGE_ADD(cascade_timing,setup_us,"setup",queryseq5,univdiagpool,stage_start_us);
    return (Pathpair_T *) NULL;
  }

  queryuc_ptr_5 = Shortread_queryuc_ptr(queryseq5);
  queryuc_ptr_3 = Shortread_queryuc_ptr(queryseq3);
  queryrc5 = Shortread_queryrc(queryseq5);
  queryrc3 = Shortread_queryrc(queryseq3);
  exact_full_overlap_p =
    paired_reads_are_exact_full_overlap_p(queryseq5,queryseq3,
					  queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					  querylength5,querylength3);

#if 0
  this5 = Stage1_new(queryuc_ptr_5,querylength5,/*first_read_p*/true);
  this3 = Stage1_new(queryuc_ptr_3,querylength3,/*first_read_p*/false);
#endif

  found_score_5 = end_score_5 = querylength5;
  found_score_3 = end_score_3 = querylength3;

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
    mincoverage_querylength_5 = querylength5;
    mincoverage_querylength_3 = querylength3;
    mincoverage_filter_5 = 0;
    mincoverage_filter_3 = 0;
  } else if (user_mincoverage_filter_float <= 1.0) {
    /* Assuming that --min-coverage=1 must mean 1.0 and not a coverage of 1 bp */
    mincoverage_querylength_5 =
      segment_quality_seedable_querylength(Shortread_quality_string(queryseq5),querylength5);
    mincoverage_querylength_3 =
      segment_quality_seedable_querylength(Shortread_quality_string(queryseq3),querylength3);
    mincoverage_filter_5 = (int) rint(user_mincoverage_filter_float * (double) mincoverage_querylength_5);
    mincoverage_filter_3 = (int) rint(user_mincoverage_filter_float * (double) mincoverage_querylength_3);
  } else {
    mincoverage_querylength_5 = querylength5;
    mincoverage_querylength_3 = querylength3;
    mincoverage_filter_5 = mincoverage_filter_3 = (int) user_mincoverage_filter_float;
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,
	  "[PAIRED_MINCOVERAGE_EFFECTIVE] read5=%s querylength5=%d effective5=%d mincoverage5=%d read3=%s querylength3=%d effective3=%d mincoverage3=%d user_mincoverage=%g\n",
	  Shortread_accession(queryseq5),querylength5,
	  user_mincoverage_filter_float <= 1.0 ? mincoverage_querylength_5 : querylength5,
	  mincoverage_filter_5,
	  Shortread_accession(queryseq3),querylength3,
	  user_mincoverage_filter_float <= 1.0 ? mincoverage_querylength_3 : querylength3,
	  mincoverage_filter_3,user_mincoverage_filter_float);
#endif

  /* Now allocated in Stage1hr_new */
  /* mismatch_positions_alloc_5 = (int *) MALLOC((querylength5+MISMATCH_EXTRA)*sizeof(int)); */
  /* mismatch_positions_alloc_3 = (int *) MALLOC((querylength3+MISMATCH_EXTRA)*sizeof(int)); */

  query5_compress_fwd = Compress_new_fwd(queryuc_ptr_5,querylength5);
  query5_compress_rev = Compress_new_rev(queryuc_ptr_5,querylength5);
  query3_compress_fwd = Compress_new_fwd(queryuc_ptr_3,querylength3);
  query3_compress_rev = Compress_new_rev(queryuc_ptr_3,querylength3);
  CASCADE_STAGE_ADD(cascade_timing,setup_us,"setup",queryseq5,univdiagpool,stage_start_us);



  if (mode == CMET_NONSTRANDED || mode == ATOI_NONSTRANDED || mode == TTOC_NONSTRANDED) {
    /* Not implemented yet */
    fprintf(stderr,"Nonstranded modes not yet implemented\n");
    exit(9);

    pathpairarray = (Pathpair_T *) NULL;
#ifdef TO_FIX
    geneplus_pathpairs = paired_read(&geneplus_abort_pairing_p,&geneplus_hits5,&geneplus_hits3,
				     &geneplus_samechr,&geneplus_conc_transloc,
				     queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				     knownsplicing,knownindels,this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
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
				      knownsplicing,knownindels,this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
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
    found_score_5 = end_score_5 = querylength5;
    found_score_3 = end_score_3 = querylength3;
    found_score_paired = querylength5 + querylength3 + 1; /* Add 1 to enter loop */
    end_score_paired = querylength5 + querylength3 + 1; /* Add 1 to enter loop */
    sufficient_score_5 = querylength5/20;
    sufficient_score_3 = querylength3/20;
    max_nsplices_5 = querylength5/EXON_LENGTH + 1;
    max_nsplices_3 = querylength3/EXON_LENGTH + 1;

    complete5_p = complete3_p = false;

    if (transcriptome_align_p == true) {
      stage_start_us = CASCADE_STAGE_START("tr",queryseq5,univdiagpool);
      Trchain_timing_reset();
      pathpairs = paired_search_tr(&found_score_paired,&found_score_5,&found_score_3,
				   &end_score_paired,&end_score_5,&end_score_3,
				   sufficient_score_5,sufficient_score_3,
				   
				   queryseq5,queryseq3,querylength5,querylength3,/*genestrand*/0,
				   this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
				   this5,this3,knownsplicing,
				   query5_compress_fwd,query5_compress_rev,
				   query3_compress_fwd,query3_compress_rev,

				   nmismatches_allowed_5,nmismatches_allowed_3,
				   nmismatches_filter_5,nmismatches_filter_3,
				   mincoverage_filter_5,mincoverage_filter_3,

				   novel_univdiags_alloc,novel_diagonals_alloc,
				   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				   localdb,localdb_nmismatches_allowed,

				   trdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool);
      CASCADE_STAGE_ADD(cascade_timing,tr_search_us,"tr",queryseq5,univdiagpool,stage_start_us);
      if (cascade_timing != NULL) {
	cascade_timing->tr_chain_k13_us += Trchain_timing_k13_us();
	cascade_timing->tr_chain_k11_us += Trchain_timing_k11_us();
	cascade_timing->tr_chain_fm_us += Trchain_timing_fm_us();
	cascade_timing->tr_chain_fm_trdiag_us += Trchain_timing_fm_trdiag_us();
	cascade_timing->tr_chain_fm_lchain_us += Trchain_timing_fm_lchain_us();
	Trchain_fm_counters_get(&(cascade_timing->tr_chain_fm_mate5),
				 &(cascade_timing->tr_chain_fm_mate3));
      }
    }

    if (genome_align_p == true &&
	(method_goal = paired_gen_method_goal()) != METHOD_INIT) {
      last_method_5 = last_method_3 = METHOD_INIT;
#if 0
      /* NEW */
      /* Strange: Calling this can result in missing an EXACT1 alignment */
      pathpairs =
	paired_search_exact(&complete5_p,&complete3_p,
			    &found_score_paired,&found_score_5,&found_score_3,
			    &end_score_paired,&end_score_5,&end_score_3,
			    &last_method_5,&last_method_3,
			    &any_imperfect_ends5_p,&any_imperfect_ends3_p,
				      
			    &unextended_pathpairs,pathpairs,this5,this3,queryseq5,queryseq3,
			    queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
			    knownsplicing,knownindels,this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			    query5_compress_fwd,query5_compress_rev,
			    query3_compress_fwd,query3_compress_rev,
			    sufficient_score_5,sufficient_score_3,/*genestrand*/0,
				      
			    nmismatches_filter_5,nmismatches_filter_3,
			    mincoverage_filter_5,mincoverage_filter_3,
				      
			    repetitive_ef64,auxinfopool,pathstorepool,univdiagpool,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,
			    pathpool,transcriptpool,vectorpool,hitlistpool,
			    spliceendsgen5,spliceendsgen3,/*method_goal*/KMER_EXACT1,pass);
#endif

      /* Comparing found_score_paired against found_score_5 +
	 found_score_3 tells us whether we have missed any good
	 ends */
      /* Exit if found_score_5 <= sufficient_score && found_score_3 <= sufficient_score &&
	 found_score_paired <= found_score_5 + found_score_3 */
      while ((found_score_paired > found_score_5 + found_score_3 ||
	      found_score_5 > sufficient_score_5 || found_score_3 > sufficient_score_3) &&
	     (paired_gen_method_available_p(last_method_5,method_goal) ||
	      paired_gen_method_available_p(last_method_3,method_goal))) {
	pathpairs =
	  paired_search_univdiagonals(&complete5_p,&complete3_p,
				      &found_score_paired,&found_score_5,&found_score_3,
				      &end_score_paired,&end_score_5,&end_score_3,
				      &last_method_5,&last_method_3,
				      &any_imperfect_ends5_p,&any_imperfect_ends3_p,
				      
				      &unresolved_pathpairs,pathpairs,this5,this3,queryseq5,queryseq3,
				      queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
				      knownsplicing,knownindels,this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
				      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				      query5_compress_fwd,query5_compress_rev,
				      query3_compress_fwd,query3_compress_rev,
				      sufficient_score_5,sufficient_score_3,
				      max_nsplices_5,max_nsplices_3,/*genestrand*/0,
				      
				      nmismatches_filter_5,nmismatches_filter_3,
				      mincoverage_filter_5,mincoverage_filter_3,
				      
				      repetitive_ef64,auxinfopool,pathstorepool,
				      univdiagpool,segmenteltpool,intlistpool,uintlistpool,univcoordlistpool,listpool,
				      pathpool,transcriptpool,vectorpool,hitlistpool,
				      spliceendsgen5,spliceendsgen3,cascade_timing,method_goal);
	debug(printf("pathpairs %p (length %d)\n",pathpairs,List_length(pathpairs)));
	debug(printf("found_score_paired %d, found_score_5 %d, found_score_3 %d\n",
		     found_score_paired,found_score_5,found_score_3));

      }

      /* TODO: If found_score_paired == found_score_5 + found_score_3, then we can filter pathpairs now */
    }
  }


  /* Want to continue searching until we found a resolved pathpair,
     but can now consider the unresolved ones.  Pathpair_eval_and_sort
     will call Pathpair_resolve */

#if 0
  /* With Segmentelt_T procedures, no longer returning exhaustive arrays */
  if (pathpairs == NULL) {
    /* Univdiagonals solved paths vs Exhaustive */
    /* Does yield results, since we might have been too strict on complete_p for the mates */
    pathpairs = paths5_mates(&found_score_paired,&found_score_3,&end_score_paired,&end_score_3,&unresolved_pathpairs,
			     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_3,queryrc3,querylength3,
			     knownsplicing,knownindels,this3->mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			     this5,this3,query3_compress_fwd,query3_compress_rev,

			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,

			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen3,/*only_complete5_p*/complete5_p);
    
    pathpairs = paths3_mates(&found_score_paired,&found_score_5,&end_score_paired,&end_score_5,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_5,queryrc5,querylength5,
			     knownsplicing,knownindels,this5->mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			     this5,this3,query5_compress_fwd,query5_compress_rev,

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
#endif


#if 0
  if (pathpairs == NULL) {
    /* Univdiagonals vs Exhaustive.  This procedure is too
       time-consuming.  Better to go to paired_search_exhaustive,
       which takes counts into account */
    pathpairs =
      univdiagonals5_mates(&found_score_paired,&found_score_5,&found_score_3,
			   &end_score_paired,&end_score_5,&end_score_3,
			   &unresolved_pathpairs,
				  
			   pathpairs,queryseq5,queryseq3,
			   queryuc_ptr_5,queryrc5,querylength5,
			   queryuc_ptr_3,queryrc3,querylength3,
			   knownsplicing,knownindels,
			   this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
			   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			   this5,this3,
			   query5_compress_fwd,query5_compress_rev,
			   query3_compress_fwd,query3_compress_rev,

			   nmismatches_filter_5,nmismatches_filter_3,
			   mincoverage_filter_5,mincoverage_filter_3,

			   univdiagpool,auxinfopool,intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			   spliceendsgen5,spliceendsgen3);
    
    pathpairs =
      univdiagonals3_mates(&found_score_paired,&found_score_5,&found_score_3,
			   &end_score_paired,&end_score_5,&end_score_3,
			   &unresolved_pathpairs,
				  
			   pathpairs,queryseq5,queryseq3,
			   queryuc_ptr_5,queryrc5,querylength5,
			   queryuc_ptr_3,queryrc3,querylength3,
			   knownsplicing,knownindels,
			   this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
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


#if 0
  /* With Segmentelt_T procedures, no longer computing exhaustive arrays */
  if (pathpairs == NULL) {
    debug(printf("Beginning paired_search_exhaustive\n"));
    /* Exhaustive vs Exhaustive.  Takes counts into account, so should
       be better than univdiagonals_mates procedures */
    pathpairs = paired_search_exhaustive(&found_score_paired,&found_score_5,&found_score_3,
					 &end_score_paired,&end_score_5,&end_score_3,
					 &unresolved_pathpairs_from_exhaustive,
					 queryseq5,queryuc_ptr_5,queryrc5,
					 queryseq3,queryuc_ptr_3,queryrc3,
					 querylength5,querylength3,

					 query5_compress_fwd,query5_compress_rev,
					 query3_compress_fwd,query3_compress_rev,
					 this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					 this5,this3,knownsplicing,knownindels,

					 nmismatches_filter_5,nmismatches_filter_3,
					 mincoverage_filter_5,mincoverage_filter_3,

					 auxinfopool,pathstorepool,univdiagpool,intlistpool,uintlistpool,
					 univcoordlistpool,listpool,pathpool,vectorpool,
					 transcriptpool,hitlistpool,
					 spliceendsgen5,spliceendsgen3);
    debug(printf("Done with paired_search_exhaustive\n"));
  }
#endif

  /* Exhaustive should succeed, unless one querylength < index1part */
  debug(printf("After exhaustive step (if run), have %d pathpairs, %d unresolved pathpairs, %d unresolved from exhaustive\n",
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

#ifdef GSNAP_PAIRED_SOLVE_STATS
  paired_solve_stats_print_pools("after-exhaustive",Shortread_accession(queryseq5),
				 listpool,univdiagpool,pathpool,hitlistpool);
#endif

#if 0
  /* With Segmentelt_T procedures, no longer returning exhaustive arrays */
  if (pathpairs == NULL) {
    solved_all_p = true;

    /* Solve all 5' solutions and find mates locally */
    solve_all(&complete5_p,&found_score_5,&end_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
	      max_nsplices_5,this5,knownsplicing,knownindels,this5->mismatch_positions_alloc,
	      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
	      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
	      
	      query5_compress_fwd,query5_compress_rev,/*first_read_p*/true,
		
	      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
	      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
	      spliceendsgen5);

    pathpairs = paths5_mates(&found_score_paired,&found_score_3,&end_score_paired,&end_score_3,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_3,queryrc3,querylength3,
			     knownsplicing,knownindels,this3->mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			     this5,this3,query3_compress_fwd,query3_compress_rev,
				     
			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
				     
			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen3,/*only_complete5_p*/complete5_p);

    /* Solve all 3' solutions and find mates locally */
    solve_all(&complete3_p,&found_score_3,&end_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
	      max_nsplices_3,this3,knownsplicing,knownindels,this3->mismatch_positions_alloc,
	      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
	      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
		
	      query3_compress_fwd,query3_compress_rev,/*first_read_p*/false,
		
	      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
	      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
	      spliceendsgen3);
    
    pathpairs = paths3_mates(&found_score_paired,&found_score_5,&end_score_paired,&end_score_5,&unresolved_pathpairs,
				     
			     pathpairs,queryseq5,queryseq3,
			     queryuc_ptr_5,queryrc5,querylength5,
			     knownsplicing,knownindels,this5->mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			     this5,this3,query5_compress_fwd,query5_compress_rev,
			     
			     nmismatches_filter_5,nmismatches_filter_3,
			     mincoverage_filter_5,mincoverage_filter_3,
				     
			     auxinfopool,pathstorepool,univdiagpool,
			     intlistpool,uintlistpool,univcoordlistpool,
			     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
			     /*mate*/spliceendsgen5,/*only_complete3_p*/complete3_p);
    solved_all_p = true;
  }
#endif

  debug(Stage1_list_all_univdiagonals(this5));
  debug(Stage1_list_all_univdiagonals(this3));

  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE)) {
    /* Previously had a call to paired_search_anchored, but now we have mates procedures */

    /* Short transcript/terminal anchors can be too weak to make the
       mate's global univdiagonal list.  Use those anchors to perform a
       bounded local mate rescue before scoring pathpairs. */
    if (paired_preeval_rescue_skip_p == true &&
	pathpairs != NULL && complete5_p == true && complete3_p == true &&
	found_score_paired <= found_score_5 + found_score_3) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
      fprintf(stderr,
	      "paired_solve_stats\tmate_rescue\tstage=pre-eval\tskipped=score-complete\tbegin_pathpairs=%d\tbegin_unresolved=%d\tcomplete5=%d\tcomplete3=%d\tfound_pair=%d\tfound5=%d\tfound3=%d\n",
	      List_length(pathpairs),List_length(unresolved_pathpairs),
	      complete5_p,complete3_p,found_score_paired,found_score_5,found_score_3);
#endif

    } else if (paired_preeval_rescue_skip_p == true &&
	       pathpairs != NULL && complete5_p == true && complete3_p == true &&
	       pathpairs_have_resolved_best_pair_p(pathpairs,found_score_paired) == true) {
#ifdef GSNAP_PAIRED_SOLVE_STATS
      fprintf(stderr,
	      "paired_solve_stats\tmate_rescue\tstage=pre-eval\tskipped=resolved-best-pair\tbegin_pathpairs=%d\tbegin_unresolved=%d\tcomplete5=%d\tcomplete3=%d\tfound_pair=%d\tfound5=%d\tfound3=%d\n",
	      List_length(pathpairs),List_length(unresolved_pathpairs),
	      complete5_p,complete3_p,found_score_paired,found_score_5,found_score_3);
#endif

    } else {
      stage_start_us = CASCADE_STAGE_START("preMR",queryseq5,univdiagpool);
      pathpairs = paired_local_mate_rescue("pre-eval",
					   &found_score_paired,&found_score_5,&found_score_3,
					   &end_score_paired,&end_score_5,&end_score_3,
					   &unresolved_pathpairs,pathpairs,
					   &local_mate_singlepaths5,&local_mate_singlepaths3,
					   queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,querylength5,
					   queryuc_ptr_3,queryrc3,querylength3,
					   knownsplicing,knownindels,
					   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					   this5,this3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   nmismatches_filter_5,nmismatches_filter_3,
					   mincoverage_filter_5,mincoverage_filter_3,
					   auxinfopool,pathstorepool,univdiagpool,
					   intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					   spliceendsgen5,spliceendsgen3,complete5_p,complete3_p,
					   exact_full_overlap_p,
					   paired_preeval_ungapped_mate_p,
					   /*reliable_anchors_only_p*/false);
      CASCADE_STAGE_ADD(cascade_timing,pre_eval_rescue_us,"preMR",queryseq5,univdiagpool,stage_start_us);
    }

#ifdef GSNAP_PAIRED_SOLVE_STATS
    paired_solve_stats_print_pools("after-pre-eval",Shortread_accession(queryseq5),
				   listpool,univdiagpool,pathpool,hitlistpool);
#endif

    debug(printf("Have %d pathpairs, %d unresolved pathpairs, complete5p %d, complete3p %d\n",
		 List_length(pathpairs),List_length(unresolved_pathpairs),
		 complete5_p,complete3_p));

    stage_start_us = CASCADE_STAGE_START("preGC",queryseq5,univdiagpool);
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
    CASCADE_STAGE_ADD(cascade_timing,pre_eval_cleanup_us,"preGC",queryseq5,univdiagpool,stage_start_us);
  }
    
  /* Terminal paired-end flow evaluates existing candidate pathpairs first,
     then optional anchor-first rescue, then full post-solve rescue/fallback. */
  while (1) {
    if (pathpairs == NULL) {
      if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_SOLVE_ALL) == false &&
	  Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE) == false) {
	stage_start_us = CASCADE_STAGE_START("fbCollect",queryseq5,univdiagpool);
	collect_single_end_paths(&sense_paths5,&antisense_paths5,
				 &sense_paths3,&antisense_paths3,
				 this5,this3,hitlistpool);
	CASCADE_STAGE_ADD(cascade_timing,fallback_collect_us,"fbCollect",queryseq5,univdiagpool,stage_start_us);
	break;
      }

      if (paired_solve_all_anchor_first_p == true &&
	  anchor_first_solve_attempted_p == false &&
	  solved_all_p == false &&
	  Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_SOLVE_ALL) == true &&
	  Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE) == true) {
	anchor_first_solve_attempted_p = true;

	if (stage1_all_univdiagonal_count(this3) <= stage1_all_univdiagonal_count(this5)) {
	  debug(printf("Calling anchor-first solve_all on 3' end\n"));
	  stage_start_us = CASCADE_STAGE_START_ROLE("afSolve","mate3",queryseq5,univdiagpool);
	  solve_all(&complete3_p,&found_score_3,&end_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
		    max_nsplices_3,this3,knownsplicing,knownindels,this3->mismatch_positions_alloc,
		    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		    query3_compress_fwd,query3_compress_rev,/*first_read_p*/false,

		    auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		    pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		    spliceendsgen3);
	  CASCADE_STAGE_ADD_ROLE(cascade_timing,anchor_first_solve_us,"afSolve","mate3",queryseq5,univdiagpool,stage_start_us);
	  debug(printf("Done with anchor-first solve_all on 3' end\n"));
	  solved_all_3_p = true;

	} else {
	  debug(printf("Calling anchor-first solve_all on 5' end\n"));
	  stage_start_us = CASCADE_STAGE_START_ROLE("afSolve","mate5",queryseq5,univdiagpool);
	  solve_all(&complete5_p,&found_score_5,&end_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
		    max_nsplices_5,this5,knownsplicing,knownindels,this5->mismatch_positions_alloc,
		    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		    query5_compress_fwd,query5_compress_rev,/*first_read_p*/true,

		    auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		    pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		    spliceendsgen5);
	  CASCADE_STAGE_ADD_ROLE(cascade_timing,anchor_first_solve_us,"afSolve","mate5",queryseq5,univdiagpool,stage_start_us);
	  debug(printf("Done with anchor-first solve_all on 5' end\n"));
	  solved_all_5_p = true;
	}
	solved_all_p = solved_all_5_p && solved_all_3_p;

	stage_start_us = CASCADE_STAGE_START("afMR",queryseq5,univdiagpool);
	pathpairs = paired_local_mate_rescue("anchor-first-solve_all",
					     &found_score_paired,&found_score_5,&found_score_3,
					     &end_score_paired,&end_score_5,&end_score_3,
					     &unresolved_pathpairs,pathpairs,
					     &local_mate_singlepaths5,&local_mate_singlepaths3,
					     queryseq5,queryseq3,
					     queryuc_ptr_5,queryrc5,querylength5,
					     queryuc_ptr_3,queryrc3,querylength3,
					     knownsplicing,knownindels,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
					     this5,this3,
					     query5_compress_fwd,query5_compress_rev,
					     query3_compress_fwd,query3_compress_rev,
					     nmismatches_filter_5,nmismatches_filter_3,
					     mincoverage_filter_5,mincoverage_filter_3,
					     auxinfopool,pathstorepool,univdiagpool,
					     intlistpool,uintlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					     spliceendsgen5,spliceendsgen3,complete5_p,complete3_p,
					     exact_full_overlap_p,
					     /*allow_ungapped_mate_p*/true,
					     /*reliable_anchors_only_p*/false);

	if (pathpairs != NULL) {
	  pathpairs = keep_local_mate_pathpairs(pathpairs,intlistpool,univcoordlistpool,
						listpool,pathpool,transcriptpool,hitlistpool);
	  Pathpair_gc(&unresolved_pathpairs,intlistpool,univcoordlistpool,listpool,pathpool,
		      transcriptpool,hitlistpool);
	}
	CASCADE_STAGE_ADD(cascade_timing,anchor_first_rescue_us,"afMR",queryseq5,univdiagpool,stage_start_us);
      }

      if (pathpairs == NULL) {
	post_solve_fallback_entered_p = true;

	if (solved_all_p == true) {
	  debug(printf("Already solved all\n"));
	  /* paired_search5_mates and paired_search3_mates were called previously after solve_all */

	} else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_SOLVE_ALL)) {
	  if (solved_all_5_p == false) {
	    debug(printf("Calling solve_all on 5' end\n"));
	    stage_start_us = CASCADE_STAGE_START_ROLE("psSolve","mate5",queryseq5,univdiagpool);
	    solve_all(&complete5_p,&found_score_5,&end_score_5,queryseq5,queryuc_ptr_5,queryrc5,querylength5,
		      max_nsplices_5,this5,knownsplicing,knownindels,this5->mismatch_positions_alloc,
		      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		      query5_compress_fwd,query5_compress_rev,/*first_read_p*/true,

		      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		      spliceendsgen5);
	    CASCADE_STAGE_ADD_ROLE(cascade_timing,postsolve_solve_us,"psSolve","mate5",queryseq5,univdiagpool,stage_start_us);
	    debug(printf("Done with solve_all on 5' end\n"));
	    solved_all_5_p = true;
	  }

	  if (solved_all_3_p == false) {
	    debug(printf("Calling solve_all on 3' end\n"));
	    stage_start_us = CASCADE_STAGE_START_ROLE("psSolve","mate3",queryseq5,univdiagpool);
	    solve_all(&complete3_p,&found_score_3,&end_score_3,queryseq3,queryuc_ptr_3,queryrc3,querylength3,
		      max_nsplices_3,this3,knownsplicing,knownindels,this3->mismatch_positions_alloc,
		      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		      query3_compress_fwd,query3_compress_rev,/*first_read_p*/false,

		      auxinfopool,intlistpool,uintlistpool,univcoordlistpool,listpool,
		      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
		      spliceendsgen3);
	    CASCADE_STAGE_ADD_ROLE(cascade_timing,postsolve_solve_us,"psSolve","mate3",queryseq5,univdiagpool,stage_start_us);
	    debug(printf("Done with solve_all on 3' end\n"));
	    solved_all_3_p = true;
	  }

	  solved_all_p = solved_all_5_p && solved_all_3_p;
	}

	if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE)) {
	  if (paired_postsolve_reliable_anchors_first_p == true) {
	    stage_start_us = CASCADE_STAGE_START("psMRr",queryseq5,univdiagpool);
	    pathpairs = paired_local_mate_rescue("post-solve_all-reliable",
						 &found_score_paired,&found_score_5,&found_score_3,
						 &end_score_paired,&end_score_5,&end_score_3,
						 &unresolved_pathpairs,pathpairs,
						 &local_mate_singlepaths5,&local_mate_singlepaths3,
						 queryseq5,queryseq3,
						 queryuc_ptr_5,queryrc5,querylength5,
						 queryuc_ptr_3,queryrc3,querylength3,
						 knownsplicing,knownindels,
						 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
						 this5,this3,
						 query5_compress_fwd,query5_compress_rev,
						 query3_compress_fwd,query3_compress_rev,
						 nmismatches_filter_5,nmismatches_filter_3,
						 mincoverage_filter_5,mincoverage_filter_3,
						 auxinfopool,pathstorepool,univdiagpool,
						 intlistpool,uintlistpool,univcoordlistpool,
						 listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
						 spliceendsgen5,spliceendsgen3,complete5_p,complete3_p,
						 exact_full_overlap_p,
						 /*allow_ungapped_mate_p*/true,
						 /*reliable_anchors_only_p*/true);
	    CASCADE_STAGE_ADD(cascade_timing,postsolve_rescue_reliable_us,"psMRr",queryseq5,univdiagpool,stage_start_us);
	  }

	  if (pathpairs == NULL || found_score_paired > found_score_5 + found_score_3) {
	    stage_start_us = CASCADE_STAGE_START("psMR",queryseq5,univdiagpool);
	    pathpairs = paired_local_mate_rescue("post-solve_all",
						 &found_score_paired,&found_score_5,&found_score_3,
						 &end_score_paired,&end_score_5,&end_score_3,
						 &unresolved_pathpairs,pathpairs,
						 &local_mate_singlepaths5,&local_mate_singlepaths3,
						 queryseq5,queryseq3,
						 queryuc_ptr_5,queryrc5,querylength5,
						 queryuc_ptr_3,queryrc3,querylength3,
						 knownsplicing,knownindels,
						 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
						 this5,this3,
						 query5_compress_fwd,query5_compress_rev,
						 query3_compress_fwd,query3_compress_rev,
						 nmismatches_filter_5,nmismatches_filter_3,
						 mincoverage_filter_5,mincoverage_filter_3,
						 auxinfopool,pathstorepool,univdiagpool,
						 intlistpool,uintlistpool,univcoordlistpool,
						 listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
						 spliceendsgen5,spliceendsgen3,complete5_p,complete3_p,
						 exact_full_overlap_p,
						 /*allow_ungapped_mate_p*/true,
						 /*reliable_anchors_only_p*/false);
	    CASCADE_STAGE_ADD(cascade_timing,postsolve_rescue_us,"psMR",queryseq5,univdiagpool,stage_start_us);
	  }
	}

	if (pathpairs == NULL) {
	  stage_start_us = CASCADE_STAGE_START("fbCollect",queryseq5,univdiagpool);
	  collect_single_end_paths(&sense_paths5,&antisense_paths5,
				   &sense_paths3,&antisense_paths3,
				   this5,this3,hitlistpool);
	  CASCADE_STAGE_ADD(cascade_timing,fallback_collect_us,"fbCollect",queryseq5,univdiagpool,stage_start_us);
	}
      }

#if 0
      if (splicingp == true) {
	/* Inner path fusions require further work, and is too time-consuming */
	if ((pathpairs = find_inner_fusions(&found_score_5,&found_score_3,&end_score_5,&end_score_3,

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
#ifdef GSNAP_PAIRED_SOLVE_STATS
      double eval_t0 = paired_solve_stats_now();
      int eval_input_npaths = List_length(pathpairs);

      fprintf(stderr,
	      "paired_solve_stats\tpathpair_eval\tbegin_pathpairs=%d\tfound_pair=%d\tfound5=%d\tfound3=%d\tend_pair=%d\tend5=%d\tend3=%d\n",
	      eval_input_npaths,found_score_paired,found_score_5,found_score_3,
	      end_score_paired,end_score_5,end_score_3);
#endif
      debug(printf("Starting Pathpair_eval_and_sort\n"));
      pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
      stage_start_us = CASCADE_STAGE_START("eval",queryseq5,univdiagpool);
      pathpairarray = Pathpair_eval_and_sort(&found_score_5,&found_score_3,&end_score_5,&end_score_3,
					     &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					     &rejected_singlepaths5,&rejected_singlepaths3,
					     pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					     query5_compress_fwd,query5_compress_rev,
					     query3_compress_fwd,query3_compress_rev,
					     queryseq5,queryseq3,
					     queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					     /*quality_string_5*/Shortread_quality_string(queryseq5),
					     /*quality_string_3*/Shortread_quality_string(queryseq3),

					     this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

					     knownsplicing,knownindels,
					     nmismatches_filter_5,nmismatches_filter_3,
					     mincoverage_filter_5,mincoverage_filter_3,
					     querylength5,querylength3,
					     univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					     &terminal_anchor_pathpairs_deferred,
					     &terminal_anchor_paths5_deferred,
					     &terminal_anchor_paths3_deferred,
					     /*filterp*/true,
					     /*allow_postfallback_terminal_anchor_rescue_p*/false,
					     /*allow_near_exact_terminal_anchor_rescue_p*/false);
      CASCADE_STAGE_ADD(cascade_timing,pathpair_eval_us,"eval",queryseq5,univdiagpool,stage_start_us);
      debug(printf("Done with Pathpair_eval_and_sort\n"));
#ifdef GSNAP_PAIRED_SOLVE_STATS
      fprintf(stderr,
	      "paired_solve_stats\tpathpair_eval\telapsed=%.6f\tinput=%d\tprimary=%d\taltloc=%d\tfirst_absmq=%d\tsecond_absmq=%d\tfound5=%d\tfound3=%d\tend5=%d\tend3=%d\n",
	      paired_solve_stats_now() - eval_t0,eval_input_npaths,
	      *npaths_primary,*npaths_altloc,*first_absmq,*second_absmq,
	      found_score_5,found_score_3,end_score_5,end_score_3);
#endif
    }

    if (pathpairs == NULL) {
      /* We never generated pathpairarray, and we already collected singlepaths */
      break;

    } else if ((*npaths_primary) + (*npaths_altloc) > 0) {
      *final_pairtype = CONCORDANT;
      break;

    } else {
      FREE_OUT(pathpairarray);

      Hitlistpool_free_list(&pathpairs,hitlistpool
			    hitlistpool_trace(__FILE__,__LINE__));
      pathpairs = (List_T) NULL;

      if (post_solve_fallback_entered_p == false &&
	  (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_SOLVE_ALL) ||
	   Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_MATE_RESCUE))) {
	debug(printf("Pathpairs were rejected; entering post-solve fallback\n"));
	continue;

      } else {
	debug(printf("Pathpairs were rejected; collecting single-end fallback paths\n"));
	stage_start_us = CASCADE_STAGE_START("fbCollect",queryseq5,univdiagpool);
	collect_single_end_paths(&sense_paths5,&antisense_paths5,
				 &sense_paths3,&antisense_paths3,
				 this5,this3,hitlistpool);
	CASCADE_STAGE_ADD(cascade_timing,fallback_collect_us,"fbCollect",queryseq5,univdiagpool,stage_start_us);
	break;
      }
    }
  }

  if ((*npaths_primary) + (*npaths_altloc) == 0) {
    debug(printf("Calling consolidate_results\n"));

    stage_start_us = CASCADE_STAGE_START("fbCons",queryseq5,univdiagpool);
    pathpairarray =
      consolidate_results(&found_score_5,&found_score_3,&end_score_5,&end_score_3,&(*final_pairtype),
			  &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
			  &(*patharray5),&(*npaths5_primary),&(*npaths5_altloc),&(*first_absmq5),&(*second_absmq5),
			  &(*patharray3),&(*npaths3_primary),&(*npaths3_altloc),&(*first_absmq3),&(*second_absmq3),
			  sense_paths5,antisense_paths5,sense_paths3,antisense_paths3,
			  rejected_singlepaths5,rejected_singlepaths3,
			  local_mate_singlepaths5,local_mate_singlepaths3,
			    
			  queryseq5,queryseq3,
			  queryuc_ptr_5,queryrc5,querylength5,queryuc_ptr_3,queryrc3,querylength3,
			    
			  knownsplicing,knownindels,
			  this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
			  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			  localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			  this5,this3,query5_compress_fwd,query5_compress_rev,query3_compress_fwd,query3_compress_rev,
			  nmismatches_filter_5,nmismatches_filter_3,mincoverage_filter_5,mincoverage_filter_3,
			  
			  intlistpool,uintlistpool,univcoordlistpool,listpool,univdiagpool,pathpool,vectorpool,
			  hitlistpool,transcriptpool,spliceendsgen5,spliceendsgen3);
    CASCADE_STAGE_ADD(cascade_timing,fallback_consolidate_us,"fbCons",queryseq5,univdiagpool,stage_start_us);
    rejected_singlepaths5 = rejected_singlepaths3 = (List_T) NULL;
    local_mate_singlepaths5 = local_mate_singlepaths3 = (List_T) NULL;

    Hitlistpool_free_list(&sense_paths5,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&sense_paths3,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&antisense_paths5,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
    Hitlistpool_free_list(&antisense_paths3,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }

  if ((*npaths_primary) + (*npaths_altloc) == 0 &&
      terminal_anchor_pathpairs_deferred != NULL &&
      (((*npaths5_primary) + (*npaths5_altloc) == 0 &&
	(*npaths3_primary) + (*npaths3_altloc) > 0) ||
       ((*npaths3_primary) + (*npaths3_altloc) == 0 &&
	(*npaths5_primary) + (*npaths5_altloc) > 0))) {
    /* Use relaxed terminal-anchor pairs only when consolidation still leaves one mate absent. */
    stage_start_us = CASCADE_STAGE_START("term",queryseq5,univdiagpool);
    pathpairs = terminal_anchor_pathpairs_deferred;
    terminal_anchor_pathpairs_deferred = (List_T) NULL;

    *npaths_primary = List_length(pathpairs);
    *npaths_altloc = 0;
    pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
    pathpairarray = Pathpair_eval_and_sort(&found_score_5,&found_score_3,&end_score_5,&end_score_3,
					   &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					   /*singlepaths5_from_rejected*/NULL,/*singlepaths3_from_rejected*/NULL,
					   pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					   query5_compress_fwd,query5_compress_rev,
					   query3_compress_fwd,query3_compress_rev,
					   queryseq5,queryseq3,
					   queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					   /*quality_string_5*/Shortread_quality_string(queryseq5),
					   /*quality_string_3*/Shortread_quality_string(queryseq3),
					   this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					   localdb_tally,localdb_overflowp,max_nblocks,
					   localdb_nmismatches_allowed,

					   knownsplicing,knownindels,
					   nmismatches_filter_5,nmismatches_filter_3,
					   mincoverage_filter_5,mincoverage_filter_3,
					   querylength5,querylength3,
					   univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					   /*terminal_anchor_pathpairs_deferred*/NULL,
					   /*terminal_anchor_paths5_deferred*/NULL,
					   /*terminal_anchor_paths3_deferred*/NULL,
					   /*filterp*/true,
					   /*allow_postfallback_terminal_anchor_rescue_p*/true,
					   /*allow_near_exact_terminal_anchor_rescue_p*/false);
    Hitlistpool_free_list(&pathpairs,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));

    if ((*npaths_primary) + (*npaths_altloc) > 0) {
      *final_pairtype = CONCORDANT;

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
    }
    CASCADE_STAGE_ADD(cascade_timing,terminal_anchor_eval_us,"term",queryseq5,univdiagpool,stage_start_us);
  }

  if ((*npaths_primary) + (*npaths_altloc) == 0 &&
      (*npaths5_primary) + (*npaths5_altloc) == 0 &&
      (*npaths3_primary) + (*npaths3_altloc) > 0 &&
      terminal_anchor_paths5_deferred != NULL) {
    /* Replay near-exact 5' terminal anchors only after consolidation leaves mate 5 absent. */
    stage_start_us = CASCADE_STAGE_START("term",queryseq5,univdiagpool);
    pathpairs = replay_near_exact_terminal_anchor_pathpairs(&unresolved_pathpairs,/*replay5p*/true,
							    *patharray3,*npaths3_primary + *npaths3_altloc,
							    terminal_anchor_paths5_deferred,
							    queryseq5,queryseq3,
							    nmismatches_filter_5,nmismatches_filter_3,
							    mincoverage_filter_5,mincoverage_filter_3,
							    intlistpool,univcoordlistpool,listpool,
							    pathpool,vectorpool,transcriptpool,hitlistpool);

    if (pathpairs != NULL) {
      *npaths_primary = List_length(pathpairs);
      *npaths_altloc = 0;
      pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
      pathpairarray = Pathpair_eval_and_sort(&found_score_5,&found_score_3,&end_score_5,&end_score_3,
					     &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					     /*singlepaths5_from_rejected*/NULL,/*singlepaths3_from_rejected*/NULL,
					     pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					     query5_compress_fwd,query5_compress_rev,
					     query3_compress_fwd,query3_compress_rev,
					     queryseq5,queryseq3,
					     queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					     /*quality_string_5*/Shortread_quality_string(queryseq5),
					     /*quality_string_3*/Shortread_quality_string(queryseq3),
					     this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,

					     knownsplicing,knownindels,
					     nmismatches_filter_5,nmismatches_filter_3,
					     mincoverage_filter_5,mincoverage_filter_3,
					     querylength5,querylength3,
					     univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					     /*terminal_anchor_pathpairs_deferred*/NULL,
					     /*terminal_anchor_paths5_deferred*/NULL,
					     /*terminal_anchor_paths3_deferred*/NULL,
					     /*filterp*/true,
					     /*allow_postfallback_terminal_anchor_rescue_p*/false,
					     /*allow_near_exact_terminal_anchor_rescue_p*/true);
      Hitlistpool_free_list(&pathpairs,hitlistpool
			    hitlistpool_trace(__FILE__,__LINE__));

      if ((*npaths_primary) + (*npaths_altloc) > 0) {
	*final_pairtype = CONCORDANT;

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
      }
    }
    CASCADE_STAGE_ADD(cascade_timing,terminal_anchor_eval_us,"term",queryseq5,univdiagpool,stage_start_us);
  }

  if ((*npaths_primary) + (*npaths_altloc) == 0 &&
      (*npaths5_primary) + (*npaths5_altloc) > 0 &&
      (*npaths3_primary) + (*npaths3_altloc) == 0 &&
      terminal_anchor_paths3_deferred != NULL) {
    /* Replay near-exact 3' terminal anchors only after consolidation leaves mate 3 absent. */
    stage_start_us = CASCADE_STAGE_START("term",queryseq5,univdiagpool);
    pathpairs = replay_near_exact_terminal_anchor_pathpairs(&unresolved_pathpairs,/*replay5p*/false,
							    *patharray5,*npaths5_primary + *npaths5_altloc,
							    terminal_anchor_paths3_deferred,
							    queryseq5,queryseq3,
							    nmismatches_filter_5,nmismatches_filter_3,
							    mincoverage_filter_5,mincoverage_filter_3,
							    intlistpool,univcoordlistpool,listpool,
							    pathpool,vectorpool,transcriptpool,hitlistpool);

    if (pathpairs != NULL) {
      *npaths_primary = List_length(pathpairs);
      *npaths_altloc = 0;
      pathpairarray = (Pathpair_T *) List_to_array_out(pathpairs,NULL);
      pathpairarray = Pathpair_eval_and_sort(&found_score_5,&found_score_3,&end_score_5,&end_score_3,
					     &(*npaths_primary),&(*npaths_altloc),&(*first_absmq),&(*second_absmq),
					     /*singlepaths5_from_rejected*/NULL,/*singlepaths3_from_rejected*/NULL,
					     pathpairarray,/*npaths*/List_length(pathpairs),this5,this3,
					     query5_compress_fwd,query5_compress_rev,
					     query3_compress_fwd,query3_compress_rev,
					     queryseq5,queryseq3,
					     queryuc_ptr_5,queryrc5,queryuc_ptr_3,queryrc3,
					     /*quality_string_5*/Shortread_quality_string(queryseq5),
					     /*quality_string_3*/Shortread_quality_string(queryseq3),
					     this5->mismatch_positions_alloc,this3->mismatch_positions_alloc,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,
					     localdb_nmismatches_allowed,

					     knownsplicing,knownindels,
					     nmismatches_filter_5,nmismatches_filter_3,
					     mincoverage_filter_5,mincoverage_filter_3,
					     querylength5,querylength3,
					     univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,vectorpool,hitlistpool,
					     /*terminal_anchor_pathpairs_deferred*/NULL,
					     /*terminal_anchor_paths5_deferred*/NULL,
					     /*terminal_anchor_paths3_deferred*/NULL,
					     /*filterp*/true,
					     /*allow_postfallback_terminal_anchor_rescue_p*/false,
					     /*allow_near_exact_terminal_anchor_rescue_p*/true);
      Hitlistpool_free_list(&pathpairs,hitlistpool
			    hitlistpool_trace(__FILE__,__LINE__));

      if ((*npaths_primary) + (*npaths_altloc) > 0) {
	*final_pairtype = CONCORDANT;

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
      }
    }
    CASCADE_STAGE_ADD(cascade_timing,terminal_anchor_eval_us,"term",queryseq5,univdiagpool,stage_start_us);
  }

  stage_start_us = CASCADE_STAGE_START("cleanup",queryseq5,univdiagpool);
  Pathpair_gc(&terminal_anchor_pathpairs_deferred,intlistpool,univcoordlistpool,listpool,pathpool,
	      transcriptpool,hitlistpool);
  Path_gc(&terminal_anchor_paths5_deferred,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);
  Path_gc(&terminal_anchor_paths3_deferred,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);

  Path_gc(&rejected_singlepaths5,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);
  Path_gc(&rejected_singlepaths3,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);
  Path_gc(&local_mate_singlepaths5,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);
  Path_gc(&local_mate_singlepaths3,intlistpool,univcoordlistpool,listpool,
	  pathpool,transcriptpool,hitlistpool);

  Compress_free(&query3_compress_rev);
  Compress_free(&query3_compress_fwd);
  Compress_free(&query5_compress_rev);
  Compress_free(&query5_compress_fwd);
  CASCADE_STAGE_ADD(cascade_timing,cleanup_us,"cleanup",queryseq5,univdiagpool,stage_start_us);

  /* FREE(queryrc3); -- Taken from Shortread */
  /* FREE(queryrc5); -- Taken from Shortread */

  debug(printf("Returning with final_pairtype %s\n",Pairtype_string(*final_pairtype)));
  debug(print_pathpairarray_contents(pathpairarray,/*n*/(*npaths_primary) + (*npaths_altloc)));

  return pathpairarray;
}


void
Stage1hr_paired_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       Genomebits_T genomebits_in, EF64_T chromosome_ef64_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       int max_deletionlen, int max_insertlength, Chrpos_T shortsplicedist, bool splicingp_in,
		       int maxpaths_search_in, int maxpaths_report_in,
			       bool *circularp_in, int pairmax_linear_in, int pairmax_circular_in,
			       int paired_segment_exact_max_mismatch_divisor_in,
			       bool paired_mate_rescue_gapped_p_in,
			       bool paired_preeval_ungapped_mate_p_in,
			       bool paired_preeval_rescue_skip_p_in,
			       bool paired_segment_stream_chains_p_in,
			       bool paired_solve_all_anchor_first_p_in,
			       bool paired_postsolve_reliable_anchors_first_p_in,
			       bool paired_seg1_exhaustive_concordance_p_in,
			       bool paired_duplicate_auxinfo_frontier_p_in,
			       unsigned int paired_mate_minpure_floor_in,
			       bool paired_transcript_completion_p_in,
			       bool paired_tr_adjacent_indel_repair_p_in,
			       bool tr_chain_search_p_in,
			       double se_segment_min_kmer_avg_quality_in,
			       int se_segment_max_kmer_bad_bases_in,
			       int se_segment_bad_base_quality_in,
			       int segment_max_avg_kmer_positions_in,
			       Align_cascade_T align_cascade_in) {

  mode = mode_in;
  index1part = index1part_in;
  index1interval = index1interval_in;
  index1part_tr = index1part_tr_in;

  transcriptome = transcriptome_in;
  genome_align_p = genome_align_p_in;
  transcriptome_align_p = transcriptome_align_p_in;

  genomebits = genomebits_in;
  localdb = (Localdb_T) NULL;

  chromosome_ef64 = chromosome_ef64_in;

  user_nmismatches_filter_float = user_nmismatches_filter_float_in;
  user_mincoverage_filter_float = user_mincoverage_filter_float_in;

  concordance_distance = (Chrpos_T) max_insertlength + shortsplicedist;
  paired_max_insertlength = max_insertlength;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  splicingp = splicingp_in;
  maxpaths_search = maxpaths_search_in;
  maxpaths_report = maxpaths_report_in;

  circularp = circularp_in;
  pairmax_linear = pairmax_linear_in;
  pairmax_circular = pairmax_circular_in;
  paired_segment_exact_max_mismatch_divisor = paired_segment_exact_max_mismatch_divisor_in;
  paired_mate_rescue_gapped_p = paired_mate_rescue_gapped_p_in;
  paired_preeval_ungapped_mate_p = paired_preeval_ungapped_mate_p_in;
  paired_preeval_rescue_skip_p = paired_preeval_rescue_skip_p_in;
  paired_segment_stream_chains_p = paired_segment_stream_chains_p_in;
  paired_solve_all_anchor_first_p = paired_solve_all_anchor_first_p_in;
  paired_postsolve_reliable_anchors_first_p = paired_postsolve_reliable_anchors_first_p_in;
  paired_seg1_exhaustive_concordance_p = paired_seg1_exhaustive_concordance_p_in;
  paired_duplicate_auxinfo_frontier_p = paired_duplicate_auxinfo_frontier_p_in;
  paired_mate_minpure_floor = paired_mate_minpure_floor_in;
  paired_transcript_completion_p = paired_transcript_completion_p_in;
  paired_tr_adjacent_indel_repair_p = paired_tr_adjacent_indel_repair_p_in;
  tr_chain_search_p = tr_chain_search_p_in;
  se_segment_min_kmer_avg_quality = se_segment_min_kmer_avg_quality_in;
  se_segment_max_kmer_bad_bases = se_segment_max_kmer_bad_bases_in;
  se_segment_bad_base_quality = se_segment_bad_base_quality_in;
  segment_max_avg_kmer_positions = segment_max_avg_kmer_positions_in;
  align_cascade = align_cascade_in;

  return;
}


void
Stage1hr_paired_pass2_setup (int max_insertlength, Chrpos_T shortsplicedist) {

  concordance_distance = (Chrpos_T) max_insertlength + shortsplicedist;
  paired_max_insertlength = max_insertlength;

  return;
}

void
Stage1hr_paired_localdb_setup (Localdb_T localdb_in) {
  localdb = localdb_in;
  return;
}
