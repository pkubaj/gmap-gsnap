#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path.h"
#include "path-solve.h"
#include "debug_groups.h"
#include "parameters.h"
#include "path-eval.h"
#include "path-trim.h"
#include "accel-profile-stats.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "mem.h"

#include "genomebits_count.h"
#include "genomebits_indel.h"
#include "genomebits_trim.h"
#include "splice.h"
#include "indel.h"
#include "spliceends.h"
#include "altsplice.h"

#include "intron.h"

#include "univdiagdef.h"
#include "junction.h"

#include "sedgesort.h"
#ifdef LARGE_GENOMES
#include "merge-uint8.h"
#include "uint8list.h"
#include "uint8table.h"
#else
#include "uintlist.h"
#endif

#include "repair.h"
#include "transcript-remap.h"

/* EXT-anchor reuse machinery (lookup/store + shadow/copy-check telemetry).
   Compiled only for the explicit reuse build, the substitution proof, or a
   stats build measuring reuse opportunity.  Off by default, so production
   builds carry no per-chain shadow-key overhead and never consult the cache. */
#if defined(GSNAP_EXT_ANCHOR_REUSE_ENABLE) || defined(GSNAP_EXT_ANCHOR_SUBSTITUTION_PROOF) || defined(GSNAP_ACCEL_PROFILE_STATS)
#define GSNAP_EXT_ANCHOR_PROOF 1
#endif


/* One of these appears to slow down program significantly, but for
   RNA velocity, we need to at least remap paths */
#define TRY_REPAIRS 1
#define REMAP_PATHS 1

/* Cannot measure exon length since we run into the start or end of the read */
/* #define MIN_EXONLEN 20 */
#define MIN_INTRONLEN 9

/* #define MAX_DEPTH_MIDDLE 5	-- Was 3, but that can be limiting with multiple indels */
/* #define MAX_DEPTH_LOCAL 2 */

#define PROB_SLOP 0.2


/* Creation and freeing of paths */
#if defined(DEBUG0) || defined(DEBUG_PAIRED_DEEP) || defined(DEBUG_SINGLE_DEEP) || defined(DEBUG_PATH_CORE)
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Decision on whether to consider path */
#if defined(DEBUG1) || defined(DEBUG_PAIRED_DEEP) || defined(DEBUG_SINGLE_DEEP) || defined(DEBUG_PATH_CORE)
#define debug1(x) x
#else
#define debug1(x)
#endif

/* endpoints_acceptable_p */
#if defined(DEBUG2) || defined(DEBUG_PAIRED_DEEP) || defined(DEBUG_SINGLE_DEEP) || defined(DEBUG_PATH_CORE)
#define debug2(x) x
#else
#define debug2(x)
#endif

/* qstart and qend resolve */
#if defined(DEBUG3) || defined(DEBUG_PAIRED_DEEP) || defined(DEBUG_SINGLE_DEEP) || defined(DEBUG_PATH_CORE)
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Path_nmatches */
#if defined(DEBUG7) || defined(DEBUG_PATH_CORE)
#define debug7(x) x
#else
#define debug7(x)
#endif

/* Path_solve_from_univdiagonal */
#if defined(DEBUG6) || defined(DEBUG_PATH_CORE)
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Path_extend */
#if defined(DEBUG11) || defined(DEBUG_PATH_CORE)
#define debug11(x) x
#else
#define debug11(x)
#endif

/* Path_solve_from_chain */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* best_path_genome */
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
#define debug13(x) x
#else
#define debug13(x)
#endif

/* Path_solve_junctions, Path_add_junctions */
#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
#define debug14(x) x
#else
#define debug14(x)
#endif

/* Location of path */
#if defined(DEBUG20) || defined(DEBUG_PATH_CORE)
#define debug20(x) x
#else
#define debug20(x)
#endif

#define SEGMENT_CACHE_SIZE 4096
#define SEGMENT_CACHE_MASK (SEGMENT_CACHE_SIZE - 1)
#define EMITTED_PATH_CACHE_SIZE 65536
#define EMITTED_PATH_CACHE_MASK (EMITTED_PATH_CACHE_SIZE - 1)
/* Single-anchor path caching is only an optimization; large path lists are safer to recompute. */
#define SEGMENT_CACHE_MAX_CACHED_PATHS 32
#define EXT_ANCHOR_RESULT_CACHE_SIZE 2

#ifdef GSNAP_EXT_ANCHOR_PROOF
static unsigned long long ext_anchor_lookup_count = 0ULL;
static unsigned long long ext_anchor_hit_count = 0ULL;
static unsigned long long ext_anchor_first_seen_count = 0ULL;
static unsigned long long ext_anchor_second_seen_count = 0ULL;
static unsigned long long ext_anchor_promoted_count = 0ULL;
static unsigned long long ext_anchor_table_reject_count = 0ULL;
static unsigned long long ext_anchor_large_result_reject_count = 0ULL;
static bool ext_anchor_report_registered_p = false;
static bool ext_anchor_report_enabled_p = false;

#define EXT_ANCHOR_COUNT(counter) do { \
    if (ext_anchor_report_enabled_p == true) { \
      __atomic_fetch_add(&(counter),1ULL,__ATOMIC_RELAXED); \
    } \
  } while (0)

static void
ext_anchor_reuse_report (void) {
  fprintf(stderr,"ext_anchor_reuse_summary\tlookups=%llu\thits=%llu\tfirst_seen=%llu\t"
	  "second_seen=%llu\tpromoted=%llu\ttable_rejects=%llu\tlarge_result_rejects=%llu\n",
	  __atomic_load_n(&ext_anchor_lookup_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_hit_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_first_seen_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_second_seen_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_promoted_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_table_reject_count,__ATOMIC_RELAXED),
	  __atomic_load_n(&ext_anchor_large_result_reject_count,__ATOMIC_RELAXED));
}
#endif

struct Path_segment_cache_entry_T {
  bool occupiedp;
  int qstart;
  int qend;
  int nmismatches;
  int querylength;
  int localdb_nmismatches_allowed;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  Univcoord_T univdiagonal;
  bool plusp;
  int genestrand;
  bool first_read_p;
  bool find_splices_p;
  List_T sense_paths;
  List_T antisense_paths;
};

struct Path_segment_qend_cache_entry_T {
  bool occupiedp;
  int head_qstart;
  int head_qend;
  int head_nmismatches;
  Univcoord_T head_univdiagonal;
  int next_qstart;
  int next_qend;
  int next_nmismatches;
  Univcoord_T next_univdiagonal;
  int querylength;
  int localdb_nmismatches_allowed;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  bool plusp;
  int genestrand;
  bool first_read_p;
  bool find_splices_p;
};

struct Path_segment_qstart_cache_entry_T {
  bool occupiedp;
  int head_qstart;
  int head_qend;
  int head_nmismatches;
  Univcoord_T head_univdiagonal;
  int next_qstart;
  int next_qend;
  int next_nmismatches;
  Univcoord_T next_univdiagonal;
  int querylength;
  int localdb_nmismatches_allowed;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  bool plusp;
  int genestrand;
  bool first_read_p;
  bool find_splices_p;
};

struct Path_qend_local_cache_entry_T {
  bool occupiedp;
  int depth;
  int querylength;
  int max_nblocks;
  int localdb_nmismatches_allowed;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  bool plusp;
  int genestrand;
  int sensedir;
  bool innerp;
  bool find_splices_p;
  bool salvagep;
  Path_T key_path;
  List_T result_paths;
};

struct Path_emitted_path_entry_T {
  bool occupiedp;
  Path_T key_path;
};

#ifdef GSNAP_EXT_ANCHOR_PROOF
struct Path_ext_anchor_result_entry_T {
  bool occupiedp;
  struct Accel_profile_ext_anchor_key_T key;
  List_T sense_paths;
  List_T antisense_paths;
};

struct Path_ext_anchor_cache_T {
  struct Path_ext_anchor_result_entry_T entries[EXT_ANCHOR_RESULT_CACHE_SIZE];
};
#endif

struct Path_segment_cache_T {
  bool emitted_path_dedup_p;
  int qstart_nentries;
  struct Path_segment_cache_entry_T entries[SEGMENT_CACHE_SIZE];
  int qend_nentries;
  struct Path_segment_qend_cache_entry_T qend_entries[SEGMENT_CACHE_SIZE];
  int qstart_pair_nentries;
  struct Path_segment_qstart_cache_entry_T qstart_pair_entries[SEGMENT_CACHE_SIZE];
  int qend_local_nentries;
  struct Path_qend_local_cache_entry_T qend_local_entries[SEGMENT_CACHE_SIZE];
  int emitted_path_nentries;
  struct Path_emitted_path_entry_T emitted_path_entries[EMITTED_PATH_CACHE_SIZE];
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
  unsigned long long qend_local_lookup_attempts;
  unsigned long long qend_local_hits;
  unsigned long long qend_local_empty_hits;
  unsigned long long qend_local_nonempty_hits;
  unsigned long long qend_local_empty_misses;
  unsigned long long qend_local_full_misses;
  unsigned long long qend_local_inserts;
  unsigned long long qend_local_skip_uncacheable_lookup;
  unsigned long long qend_local_skip_uncacheable_insert;
  unsigned long long qend_local_skip_large_result;
  unsigned long long qend_local_skip_full;
  unsigned long long emitted_path_marks;
  unsigned long long emitted_path_skips;
  unsigned long long emitted_path_skip_uncacheable;
  unsigned long long emitted_path_skip_full;
#endif
};


static unsigned long long
segment_cache_mix (unsigned long long h, unsigned long long x) {
  h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

static unsigned int
segment_cache_hash (Univdiag_T univdiag, int querylength, int localdb_nmismatches_allowed,
		   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		   bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) univdiag->qstart);
  h = segment_cache_mix(h,(unsigned long long) univdiag->qend);
  h = segment_cache_mix(h,(unsigned long long) univdiag->nmismatches);
  h = segment_cache_mix(h,(unsigned long long) univdiag->univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) querylength);
  h = segment_cache_mix(h,(unsigned long long) localdb_nmismatches_allowed);
  h = segment_cache_mix(h,(unsigned long long) chrnum);
  h = segment_cache_mix(h,(unsigned long long) chroffset);
  h = segment_cache_mix(h,(unsigned long long) chrhigh);
  h = segment_cache_mix(h,(unsigned long long) plusp);
  h = segment_cache_mix(h,(unsigned long long) genestrand);
  h = segment_cache_mix(h,(unsigned long long) first_read_p);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);

  return (unsigned int) (h & SEGMENT_CACHE_MASK);
}

static bool
segment_cache_entry_matches_p (struct Path_segment_cache_entry_T *entry, Univdiag_T univdiag,
			      int querylength, int localdb_nmismatches_allowed,
			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			      bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  return entry->qstart == univdiag->qstart &&
    entry->qend == univdiag->qend &&
    entry->nmismatches == univdiag->nmismatches &&
    entry->querylength == querylength &&
    entry->localdb_nmismatches_allowed == localdb_nmismatches_allowed &&
    entry->chrnum == chrnum &&
    entry->chroffset == chroffset &&
    entry->chrhigh == chrhigh &&
    entry->univdiagonal == univdiag->univdiagonal &&
    entry->plusp == plusp &&
    entry->genestrand == genestrand &&
    entry->first_read_p == first_read_p &&
    entry->find_splices_p == find_splices_p;
}

static List_T
segment_cache_copy_paths (List_T old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T new = NULL, p;

  for (p = old; p != NULL; p = List_next(p)) {
    new = Hitlist_push(new,hitlistpool,
		       (void *) Path_copy((Path_T) List_head(p),intlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool,hitlistpool)
		       hitlistpool_trace(__FILE__,__LINE__));
  }

  return List_reverse(new);
}

static List_T
solver_result_copy_paths (List_T old, Intlistpool_T intlistpool,
			  Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool,
			  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			  Transcriptpool_T transcriptpool) {
  List_T new = NULL, p;
  for (p = old; p != NULL; p = List_next(p)) {
    new = Hitlist_push(new,hitlistpool,
	(void *) Path_copy_solver_result((Path_T) List_head(p),intlistpool,
	  univcoordlistpool,listpool,pathpool,vectorpool,transcriptpool,hitlistpool)
	hitlistpool_trace(__FILE__,__LINE__));
  }
  return List_reverse(new);
}

static bool
segment_cache_small_path_list_p (List_T sense_paths, List_T antisense_paths) {
  List_T p;
  int npaths = 0;

  for (p = sense_paths; p != NULL; p = List_next(p)) {
    if (++npaths > SEGMENT_CACHE_MAX_CACHED_PATHS) {
      return false;
    }
  }

  for (p = antisense_paths; p != NULL; p = List_next(p)) {
    if (++npaths > SEGMENT_CACHE_MAX_CACHED_PATHS) {
      return false;
    }
  }

  return true;
}

static unsigned long long
qend_local_intlist_hash (Intlist_T list) {
  Intlist_T p;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = Intlist_next(p)) {
    hash = segment_cache_mix(hash,(unsigned long long) Intlist_head(p));
    n++;
  }

  return segment_cache_mix(hash,(unsigned long long) n);
}

static unsigned long long
qend_local_univcoordlist_hash (Univcoordlist_T list) {
  Univcoordlist_T p;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = Univcoordlist_next(p)) {
    hash = segment_cache_mix(hash,(unsigned long long) Univcoordlist_head(p));
    n++;
  }

  return segment_cache_mix(hash,(unsigned long long) n);
}

static unsigned long long
qend_local_junction_list_hash (List_T list) {
  List_T p;
  Junction_T junction;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (junction == (Junction_T) NULL) {
      hash = segment_cache_mix(hash,0xabcdefULL);
    } else {
      hash = segment_cache_mix(hash,(unsigned long long) junction->type);
      hash = segment_cache_mix(hash,(unsigned long long) junction->nindels);
      if (junction->type == SPLICE_JUNCTION || junction->type == AMB_JUNCTION) {
	hash = segment_cache_mix(hash,(unsigned long long) junction->splice_distance);
	hash = segment_cache_mix(hash,(unsigned long long) junction->sensedir);
      } else if (junction->type == CHIMERA_JUNCTION) {
	hash = segment_cache_mix(hash,(unsigned long long) junction->donor1);
	hash = segment_cache_mix(hash,(unsigned long long) junction->donor2);
	hash = segment_cache_mix(hash,(unsigned long long) junction->acceptor1);
	hash = segment_cache_mix(hash,(unsigned long long) junction->acceptor2);
      }
    }
    n++;
  }

  return segment_cache_mix(hash,(unsigned long long) n);
}

static bool
qend_local_junction_equal_p (Junction_T x, Junction_T y) {
  if (x == y) {
    return true;
  } else if (x == (Junction_T) NULL || y == (Junction_T) NULL) {
    return false;
  } else if (x->type != y->type || x->nindels != y->nindels) {
    return false;
  } else if (x->type == SPLICE_JUNCTION || x->type == AMB_JUNCTION) {
    return x->splice_distance == y->splice_distance && x->sensedir == y->sensedir &&
      x->donor_prob == y->donor_prob && x->acceptor_prob == y->acceptor_prob;
  } else if (x->type == CHIMERA_JUNCTION) {
    return x->donor1 == y->donor1 && x->donor2 == y->donor2 &&
      x->acceptor1 == y->acceptor1 && x->acceptor2 == y->acceptor2 &&
      x->donor_prob == y->donor_prob && x->acceptor_prob == y->acceptor_prob;
  } else {
    return true;
  }
}

static bool
qend_local_junction_list_equal_p (List_T x, List_T y) {
  while (x != (List_T) NULL && y != (List_T) NULL) {
    if (qend_local_junction_equal_p((Junction_T) List_head(x),(Junction_T) List_head(y)) == false) {
      return false;
    } else {
      x = List_next(x);
      y = List_next(y);
    }
  }

  return (x == (List_T) NULL && y == (List_T) NULL) ? true : false;
}

static bool
qend_local_cacheable_path_p (Path_T path) {
  if (path == (Path_T) NULL) {
    return false;
  } else if (path->qstart_alts != NULL || path->qend_alts != NULL) {
    return false;
  } else if (path->circular_endpoints != NULL) {
    return false;
  } else if (path->fusion_querystart_junction != NULL || path->fusion_queryend_junction != NULL ||
	     path->fusion_endpoints != NULL) {
    return false;
  } else if (path->transcripts != NULL || path->invalid_transcripts != NULL ||
	     path->fusion_transcripts != NULL || path->fusion_invalid_transcripts != NULL ||
	     path->genomic_diff != NULL || path->transcriptome_method_p == true) {
    return false;
  } else {
    return true;
  }
}

#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
static bool
qend_local_cache_stats_print_p (unsigned long long count) {
  return (count <= 16ULL || (count & 0x3ffULL) == 0ULL) ? true : false;
}

static void
qend_local_cache_stats_report (Path_segment_cache_T cache, const char *event, Path_T path, int depth) {
  fprintf(stderr,
	  "[QEND_LOCAL_CACHE] event=%s attempts=%llu hits=%llu empty_hits=%llu nonempty_hits=%llu inserts=%llu empty_misses=%llu "
	  "full_misses=%llu skip_uncache_lookup=%llu skip_uncache_insert=%llu "
	  "skip_large=%llu skip_full=%llu path_marks=%llu path_skips=%llu "
	  "path_uncacheable=%llu path_full=%llu nentries=%d path_nentries=%d depth=%d",
	  event,cache->qend_local_lookup_attempts,cache->qend_local_hits,
	  cache->qend_local_empty_hits,cache->qend_local_nonempty_hits,
	  cache->qend_local_inserts,cache->qend_local_empty_misses,
	  cache->qend_local_full_misses,cache->qend_local_skip_uncacheable_lookup,
	  cache->qend_local_skip_uncacheable_insert,cache->qend_local_skip_large_result,
	  cache->qend_local_skip_full,cache->emitted_path_marks,
	  cache->emitted_path_skips,cache->emitted_path_skip_uncacheable,
	  cache->emitted_path_skip_full,cache->qend_local_nentries,
	  cache->emitted_path_nentries,depth);

  if (path != (Path_T) NULL) {
    fprintf(stderr," qstart=%d qend=%d main=%u method=%d sensedir=%d\n",
	    Path_qstart(path),Path_qend(path),path->main_univdiagonal,path->method,path->sensedir);
  } else {
    fprintf(stderr,"\n");
  }

  return;
}

static void
qend_local_cache_stats_maybe_report (Path_segment_cache_T cache, Path_T path, int depth) {
  if ((cache->qend_local_lookup_attempts & 0x3ffULL) == 0ULL) {
    qend_local_cache_stats_report(cache,"periodic",path,depth);
  }

  return;
}
#endif

static unsigned int
qend_local_cache_hash (Path_T path, int depth, int querylength, int max_nblocks,
		       int localdb_nmismatches_allowed,
		       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		       bool plusp, int genestrand, int sensedir, bool innerp,
		       bool find_splices_p, bool salvagep) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) path->plusp);
  h = segment_cache_mix(h,(unsigned long long) path->genestrand);
  h = segment_cache_mix(h,(unsigned long long) path->sensedir);
  h = segment_cache_mix(h,(unsigned long long) path->querylength);
  h = segment_cache_mix(h,(unsigned long long) path->chrnum);
  h = segment_cache_mix(h,(unsigned long long) path->chroffset);
  h = segment_cache_mix(h,(unsigned long long) path->chrhigh);
  h = segment_cache_mix(h,(unsigned long long) path->main_univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) path->method);
  h = segment_cache_mix(h,(unsigned long long) path->splice5p);
  h = segment_cache_mix(h,(unsigned long long) path->splice3p);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype5);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype3);
  h = segment_cache_mix(h,qend_local_intlist_hash(path->endpoints));
  h = segment_cache_mix(h,qend_local_univcoordlist_hash(path->univdiagonals));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->nmismatches));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->ref_nmismatches));
  h = segment_cache_mix(h,qend_local_junction_list_hash(path->junctions));

  h = segment_cache_mix(h,(unsigned long long) depth);
  h = segment_cache_mix(h,(unsigned long long) querylength);
  h = segment_cache_mix(h,(unsigned long long) max_nblocks);
  h = segment_cache_mix(h,(unsigned long long) localdb_nmismatches_allowed);
  h = segment_cache_mix(h,(unsigned long long) chrnum);
  h = segment_cache_mix(h,(unsigned long long) chroffset);
  h = segment_cache_mix(h,(unsigned long long) chrhigh);
  h = segment_cache_mix(h,(unsigned long long) plusp);
  h = segment_cache_mix(h,(unsigned long long) genestrand);
  h = segment_cache_mix(h,(unsigned long long) sensedir);
  h = segment_cache_mix(h,(unsigned long long) innerp);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);
  h = segment_cache_mix(h,(unsigned long long) salvagep);

  return (unsigned int) (h & EMITTED_PATH_CACHE_MASK);
}

static bool
qend_local_path_equal_p (Path_T x, Path_T y) {
  if (x == y) {
    return true;
  } else if (x == (Path_T) NULL || y == (Path_T) NULL) {
    return false;
  } else if (x->plusp == y->plusp && x->genestrand == y->genestrand &&
	     x->sensedir == y->sensedir && x->querylength == y->querylength &&
	     x->chrnum == y->chrnum && x->chroffset == y->chroffset &&
	     x->chrhigh == y->chrhigh &&
	     x->main_univdiagonal == y->main_univdiagonal &&
	     x->method == y->method &&
	     x->splice5p == y->splice5p && x->splice3p == y->splice3p &&
	     x->splicetype5 == y->splicetype5 && x->splicetype3 == y->splicetype3 &&
	     x->ambig_prob_5 == y->ambig_prob_5 && x->ambig_prob_3 == y->ambig_prob_3 &&
	     Intlist_equal(x->endpoints,y->endpoints) == true &&
	     Univcoordlist_equal(x->univdiagonals,y->univdiagonals) == true &&
	     Intlist_equal(x->nmismatches,y->nmismatches) == true &&
	     Intlist_equal(x->ref_nmismatches,y->ref_nmismatches) == true &&
	     qend_local_junction_list_equal_p(x->junctions,y->junctions) == true) {
    return true;
  } else {
    return false;
  }
}

#ifdef GSNAP_ACCEL_PROFILE_STATS
static unsigned long long
pathsolve_path_fingerprint (Path_T path) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) path->plusp);
  h = segment_cache_mix(h,(unsigned long long) path->genestrand);
  h = segment_cache_mix(h,(unsigned long long) path->sensedir);
  h = segment_cache_mix(h,(unsigned long long) path->querylength);
  h = segment_cache_mix(h,(unsigned long long) path->chrnum);
  h = segment_cache_mix(h,(unsigned long long) path->chroffset);
  h = segment_cache_mix(h,(unsigned long long) path->chrhigh);
  h = segment_cache_mix(h,(unsigned long long) path->main_univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) path->method);
  h = segment_cache_mix(h,(unsigned long long) path->splice5p);
  h = segment_cache_mix(h,(unsigned long long) path->splice3p);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype5);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype3);
  h = segment_cache_mix(h,qend_local_intlist_hash(path->endpoints));
  h = segment_cache_mix(h,qend_local_univcoordlist_hash(path->univdiagonals));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->nmismatches));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->ref_nmismatches));
  h = segment_cache_mix(h,qend_local_junction_list_hash(path->junctions));
  return h;
}

static unsigned long long
pathsolve_frontier_fingerprint (Path_T path, const char *side, int depth,
				bool innerp, bool find_splices_p, bool salvagep) {
  unsigned long long h = 1469598103934665603ULL;
  Junction_T junction = path->junctions == NULL ? NULL : (Junction_T) List_head(path->junctions);

  h = segment_cache_mix(h,(unsigned long long) path->plusp);
  h = segment_cache_mix(h,(unsigned long long) path->genestrand);
  h = segment_cache_mix(h,(unsigned long long) path->sensedir);
  h = segment_cache_mix(h,(unsigned long long) path->querylength);
  h = segment_cache_mix(h,(unsigned long long) path->chrnum);
  h = segment_cache_mix(h,(unsigned long long) path->chroffset);
  h = segment_cache_mix(h,(unsigned long long) path->chrhigh);
  h = segment_cache_mix(h,(unsigned long long) Intlist_head(path->endpoints));
  h = segment_cache_mix(h,(unsigned long long) Univcoordlist_head(path->univdiagonals));
  h = segment_cache_mix(h,(unsigned long long) Intlist_head(path->nmismatches));
  h = segment_cache_mix(h,(unsigned long long) Intlist_head(path->ref_nmismatches));
  if (junction == NULL) {
    h = segment_cache_mix(h,0ULL);
  } else {
    h = segment_cache_mix(h,(unsigned long long) junction->type);
    h = segment_cache_mix(h,(unsigned long long) junction->nindels);
    h = segment_cache_mix(h,(unsigned long long) junction->splice_distance);
    h = segment_cache_mix(h,(unsigned long long) junction->sensedir);
  }
  h = segment_cache_mix(h,(unsigned long long) (side != NULL && !strcmp(side,"qstart")));
  h = segment_cache_mix(h,(unsigned long long) depth);
  h = segment_cache_mix(h,(unsigned long long) innerp);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);
  h = segment_cache_mix(h,(unsigned long long) salvagep);
  return h;
}

static unsigned long long
spliceends_call_fingerprint (const char *side, int sensedir, Univcoord_T univdiagonal,
			     int querylength, int boundary, int exon_origin,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, int localdb_nmismatches_allowed,
			     bool innerp, bool salvagep, int max_nblocks,
			     Stage1_T stage1, Knownsplicing_T knownsplicing,
			     Spliceendsgen_T spliceendsgen, Compress_T query_compress,
			     char *queryptr, Genomebits_T genomebits,
			     Genomebits_T genomebits_alt, bool find_splices_p) {
  unsigned long long h = 1469598103934665603ULL;
  h = segment_cache_mix(h,(unsigned long long) (side != NULL && !strcmp(side,"qstart")));
  h = segment_cache_mix(h,(unsigned long long) sensedir);
  h = segment_cache_mix(h,(unsigned long long) univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) querylength);
  h = segment_cache_mix(h,(unsigned long long) boundary);
  h = segment_cache_mix(h,(unsigned long long) exon_origin);
  h = segment_cache_mix(h,(unsigned long long) chrnum);
  h = segment_cache_mix(h,(unsigned long long) chroffset);
  h = segment_cache_mix(h,(unsigned long long) chrhigh);
  h = segment_cache_mix(h,(unsigned long long) plusp);
  h = segment_cache_mix(h,(unsigned long long) genestrand);
  h = segment_cache_mix(h,(unsigned long long) localdb_nmismatches_allowed);
  h = segment_cache_mix(h,(unsigned long long) innerp);
  h = segment_cache_mix(h,(unsigned long long) salvagep);
  h = segment_cache_mix(h,(unsigned long long) max_nblocks);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) stage1);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) knownsplicing);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) spliceendsgen);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) query_compress);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) queryptr);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) genomebits);
  h = segment_cache_mix(h,(unsigned long long) (uintptr_t) genomebits_alt);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);
  return h;
}

static unsigned long long
spliceends_result_fingerprint (Spliceends_T spliceends, int nspliceends,
			       int nosplice_trimpos, int farsplice_trimpos,
			       int nosplice_nmismatches, int farsplice_nmismatches,
			       bool splicep, Splicetype_T splicetype,
			       double ambig_prob) {
  unsigned long long h = 1469598103934665603ULL;
  union { double d; unsigned long long u; } bits;
  int i;
  h = segment_cache_mix(h,(unsigned long long) nspliceends);
  h = segment_cache_mix(h,(unsigned long long) nosplice_trimpos);
  h = segment_cache_mix(h,(unsigned long long) farsplice_trimpos);
  h = segment_cache_mix(h,(unsigned long long) nosplice_nmismatches);
  h = segment_cache_mix(h,(unsigned long long) farsplice_nmismatches);
  h = segment_cache_mix(h,(unsigned long long) splicep);
  h = segment_cache_mix(h,(unsigned long long) splicetype);
  bits.d = ambig_prob;
  h = segment_cache_mix(h,bits.u);
  if (spliceends != NULL) {
    h = segment_cache_mix(h,(unsigned long long) spliceends->boundedp);
    h = segment_cache_mix(h,(unsigned long long) spliceends->nspliceends);
    h = segment_cache_mix(h,(unsigned long long) spliceends->splicetype);
    h = segment_cache_mix(h,(unsigned long long) spliceends->sensedir);
    for (i = 0; i < nspliceends; i++) {
      h = segment_cache_mix(h,(unsigned long long) spliceends->splice_qpos[i]);
      h = segment_cache_mix(h,(unsigned long long) spliceends->medial_nmismatches[i]);
      h = segment_cache_mix(h,(unsigned long long) spliceends->partners[i]);
      h = segment_cache_mix(h,(unsigned long long) spliceends->distal_univdiags[i]->univdiagonal);
      h = segment_cache_mix(h,(unsigned long long) spliceends->distal_univdiags[i]->qstart);
      h = segment_cache_mix(h,(unsigned long long) spliceends->distal_univdiags[i]->qend);
      h = segment_cache_mix(h,(unsigned long long) spliceends->distal_univdiags[i]->nmismatches);
      bits.d = spliceends->medial_probs[i];
      h = segment_cache_mix(h,bits.u);
      bits.d = spliceends->distal_probs[i];
      h = segment_cache_mix(h,bits.u);
      h = segment_cache_mix(h,(unsigned long long) spliceends->medial_support[i]);
      h = segment_cache_mix(h,(unsigned long long) spliceends->distal_support[i]);
    }
  }
  return h;
}

static void
profile_pathsolve_result (unsigned long long state_id, const char *side, int depth,
			  unsigned long long input_path_hash,
			  unsigned long long frontier_hash,
			  List_T paths, List_T old_head) {
  unsigned long long h = 1469598103934665603ULL;
  List_T p;
  int npaths = 0;

  if (state_id == 0ULL) return;
  for (p = paths; p != old_head && p != NULL; p = List_next(p)) {
    h = segment_cache_mix(h,pathsolve_path_fingerprint((Path_T) List_head(p)));
    npaths++;
  }
  /* If old_head is not reached, ownership/list topology changed unexpectedly;
     encode that fact rather than silently treating a partial walk as exact. */
  if (p != old_head) h = segment_cache_mix(h,0xbad0bad0bad0bad0ULL);
  Accel_profile_pathsolve_result(state_id,side,depth,input_path_hash,frontier_hash,npaths,h);
}
#endif

static unsigned int
emitted_path_hash (Path_T path) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) path->plusp);
  h = segment_cache_mix(h,(unsigned long long) path->genestrand);
  h = segment_cache_mix(h,(unsigned long long) path->sensedir);
  h = segment_cache_mix(h,(unsigned long long) path->querylength);
  h = segment_cache_mix(h,(unsigned long long) path->chrnum);
  h = segment_cache_mix(h,(unsigned long long) path->chroffset);
  h = segment_cache_mix(h,(unsigned long long) path->chrhigh);
  h = segment_cache_mix(h,(unsigned long long) path->main_univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) path->method);
  h = segment_cache_mix(h,(unsigned long long) path->splice5p);
  h = segment_cache_mix(h,(unsigned long long) path->splice3p);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype5);
  h = segment_cache_mix(h,(unsigned long long) path->splicetype3);
  h = segment_cache_mix(h,qend_local_intlist_hash(path->endpoints));
  h = segment_cache_mix(h,qend_local_univcoordlist_hash(path->univdiagonals));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->nmismatches));
  h = segment_cache_mix(h,qend_local_intlist_hash(path->ref_nmismatches));
  h = segment_cache_mix(h,qend_local_junction_list_hash(path->junctions));

  return (unsigned int) (h & EMITTED_PATH_CACHE_MASK);
}

static bool
emitted_path_seen_or_mark (Path_segment_cache_T cache, Path_T path,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_emitted_path_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return false;
  } else if (qend_local_cacheable_path_p(path) == false) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    cache->emitted_path_skip_uncacheable++;
    if (qend_local_cache_stats_print_p(cache->emitted_path_skip_uncacheable) == true) {
      qend_local_cache_stats_report(cache,"path-uncacheable",path,/*depth*/-1);
    }
#endif
    return false;
  }

  hash = emitted_path_hash(path);
  for (i = 0; i < EMITTED_PATH_CACHE_SIZE; i++) {
    entry = &(cache->emitted_path_entries[(hash + i) & EMITTED_PATH_CACHE_MASK]);

    if (entry->occupiedp == false) {
      if (cache->emitted_path_nentries >= EMITTED_PATH_CACHE_SIZE) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
	cache->emitted_path_skip_full++;
	if (qend_local_cache_stats_print_p(cache->emitted_path_skip_full) == true) {
	  qend_local_cache_stats_report(cache,"path-full-table",path,/*depth*/-1);
	}
#endif
	return false;
      }

      entry->occupiedp = true;
      entry->key_path =
	Path_copy(path,intlistpool,univcoordlistpool,
		  listpool,pathpool,vectorpool,transcriptpool,hitlistpool);
      cache->emitted_path_nentries++;
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
      cache->emitted_path_marks++;
      if (qend_local_cache_stats_print_p(cache->emitted_path_marks) == true) {
	qend_local_cache_stats_report(cache,"path-mark",path,/*depth*/-1);
      }
#endif
      return false;

    } else if (qend_local_path_equal_p(entry->key_path,path) == true) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
      cache->emitted_path_skips++;
      if (qend_local_cache_stats_print_p(cache->emitted_path_skips) == true) {
	qend_local_cache_stats_report(cache,"path-skip",path,/*depth*/-1);
      }
#endif
      return true;
    }
  }

#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
  cache->emitted_path_skip_full++;
  if (qend_local_cache_stats_print_p(cache->emitted_path_skip_full) == true) {
    qend_local_cache_stats_report(cache,"path-full-table",path,/*depth*/-1);
  }
#endif

  return false;
}

static List_T
filter_emitted_paths (List_T paths, Path_segment_cache_T cache,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		      Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T kept = (List_T) NULL, p;
  Path_T path;

  if (cache == (Path_segment_cache_T) NULL || paths == (List_T) NULL) {
    return paths;
  }

  for (p = paths; p != (List_T) NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (emitted_path_seen_or_mark(cache,path,intlistpool,univcoordlistpool,
				  listpool,pathpool,vectorpool,hitlistpool,transcriptpool) == true) {
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    } else {
      kept = Hitlist_push(kept,hitlistpool,(void *) path
			  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  Hitlistpool_free_list(&paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  return List_reverse(kept);
}

static bool
qend_local_cache_entry_matches_p (struct Path_qend_local_cache_entry_T *entry,
				  Path_T path, int depth, int querylength, int max_nblocks,
				  int localdb_nmismatches_allowed,
				  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				  bool plusp, int genestrand, int sensedir, bool innerp,
				  bool find_splices_p, bool salvagep) {
  return entry->depth == depth &&
    entry->querylength == querylength &&
    entry->max_nblocks == max_nblocks &&
    entry->localdb_nmismatches_allowed == localdb_nmismatches_allowed &&
    entry->chrnum == chrnum &&
    entry->chroffset == chroffset &&
    entry->chrhigh == chrhigh &&
    entry->plusp == plusp &&
    entry->genestrand == genestrand &&
    entry->sensedir == sensedir &&
    entry->innerp == innerp &&
    entry->find_splices_p == find_splices_p &&
    entry->salvagep == salvagep &&
    qend_local_path_equal_p(entry->key_path,path) == true;
}


static bool
qend_local_cache_lookup (List_T *qend_paths, Path_segment_cache_T cache, Path_T path,
			 int depth, int querylength, int max_nblocks, int localdb_nmismatches_allowed,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool plusp, int genestrand, int sensedir, bool innerp,
			 bool find_splices_p, bool salvagep,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_qend_local_cache_entry_T *entry;
  List_T copies, p;
  Path_T result_path;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return false;
  } else if (qend_local_cacheable_path_p(path) == false) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    cache->qend_local_skip_uncacheable_lookup++;
    if (qend_local_cache_stats_print_p(cache->qend_local_skip_uncacheable_lookup) == true) {
      qend_local_cache_stats_report(cache,"lookup-uncacheable",path,depth);
    }
#endif
    return false;
  }

  #ifdef GSNAP_QEND_LOCAL_CACHE_STATS
  cache->qend_local_lookup_attempts++;
  qend_local_cache_stats_maybe_report(cache,path,depth);
#endif

  hash = qend_local_cache_hash(path,depth,querylength,max_nblocks,localdb_nmismatches_allowed,
			       chrnum,chroffset,chrhigh,plusp,genestrand,
			       sensedir,innerp,find_splices_p,salvagep);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qend_local_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == false) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
      cache->qend_local_empty_misses++;
#endif
      return false;
    } else if (qend_local_cache_entry_matches_p(entry,path,depth,querylength,max_nblocks,
						localdb_nmismatches_allowed,chrnum,
						chroffset,chrhigh,plusp,genestrand,
						sensedir,innerp,find_splices_p,salvagep) == true) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
      cache->qend_local_hits++;
      if (entry->result_paths == NULL) {
	cache->qend_local_empty_hits++;
      } else {
	cache->qend_local_nonempty_hits++;
      }
      if (qend_local_cache_stats_print_p(cache->qend_local_hits) == true) {
	qend_local_cache_stats_report(cache,"hit",path,depth);
      }
#endif
      if (cache->emitted_path_dedup_p == true) {
	copies = (List_T) NULL;
	for (p = entry->result_paths; p != (List_T) NULL; p = List_next(p)) {
	  result_path = (Path_T) List_head(p);
	  if (emitted_path_seen_or_mark(cache,result_path,intlistpool,univcoordlistpool,
					listpool,pathpool,vectorpool,hitlistpool,transcriptpool) == false) {
	    copies = Hitlist_push(copies,hitlistpool,
				  (void *) Path_copy(result_path,intlistpool,univcoordlistpool,
						     listpool,pathpool,vectorpool,transcriptpool,hitlistpool)
				  hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
	copies = List_reverse(copies);
      } else {
	copies = segment_cache_copy_paths(entry->result_paths,intlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
      }
      *qend_paths = List_append(copies,*qend_paths);
      return true;
    }
  }

#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
  cache->qend_local_full_misses++;
  if (qend_local_cache_stats_print_p(cache->qend_local_full_misses) == true) {
    qend_local_cache_stats_report(cache,"lookup-full-table",path,depth);
  }
#endif

  return false;
}

static bool
qend_local_added_paths_small_p (List_T current, List_T old_head) {
  List_T p;
  int npaths = 0;

  for (p = current; p != old_head; p = List_next(p)) {
    if (p == (List_T) NULL || ++npaths > SEGMENT_CACHE_MAX_CACHED_PATHS) {
      return false;
    }
  }

  return true;
}

static List_T
qend_local_copy_added_paths (List_T current, List_T old_head,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			     Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T copies = NULL, p;

  for (p = current; p != old_head; p = List_next(p)) {
    copies = Hitlist_push(copies,hitlistpool,
			  (void *) Path_copy((Path_T) List_head(p),intlistpool,univcoordlistpool,
					     listpool,pathpool,vectorpool,transcriptpool,hitlistpool)
			  hitlistpool_trace(__FILE__,__LINE__));
  }

  return List_reverse(copies);
}

static void
qend_local_cache_insert (Path_segment_cache_T cache, Path_T path,
			 List_T old_qend_paths, List_T current_qend_paths,
			 int depth, int querylength, int max_nblocks, int localdb_nmismatches_allowed,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool plusp, int genestrand, int sensedir, bool innerp,
			 bool find_splices_p, bool salvagep,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_qend_local_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return;
  } else if (cache->qend_local_nentries >= SEGMENT_CACHE_SIZE) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    cache->qend_local_skip_full++;
    if (qend_local_cache_stats_print_p(cache->qend_local_skip_full) == true) {
      qend_local_cache_stats_report(cache,"insert-full-table",path,depth);
    }
#endif
    return;
  } else if (qend_local_cacheable_path_p(path) == false) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    cache->qend_local_skip_uncacheable_insert++;
    if (qend_local_cache_stats_print_p(cache->qend_local_skip_uncacheable_insert) == true) {
      qend_local_cache_stats_report(cache,"insert-uncacheable",path,depth);
    }
#endif
    return;
  } else if (qend_local_added_paths_small_p(current_qend_paths,old_qend_paths) == false) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    cache->qend_local_skip_large_result++;
    if (qend_local_cache_stats_print_p(cache->qend_local_skip_large_result) == true) {
      qend_local_cache_stats_report(cache,"insert-large-result",path,depth);
    }
#endif
    return;
  }

  hash = qend_local_cache_hash(path,depth,querylength,max_nblocks,localdb_nmismatches_allowed,
			       chrnum,chroffset,chrhigh,plusp,genestrand,
			       sensedir,innerp,find_splices_p,salvagep);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qend_local_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == true) {
      if (qend_local_cache_entry_matches_p(entry,path,depth,querylength,max_nblocks,
					   localdb_nmismatches_allowed,chrnum,
					   chroffset,chrhigh,plusp,genestrand,
					   sensedir,innerp,find_splices_p,salvagep) == true) {
	return;
      }
    } else {
      entry->occupiedp = true;
      entry->depth = depth;
      entry->querylength = querylength;
      entry->max_nblocks = max_nblocks;
      entry->localdb_nmismatches_allowed = localdb_nmismatches_allowed;
      entry->chrnum = chrnum;
      entry->chroffset = chroffset;
      entry->chrhigh = chrhigh;
      entry->plusp = plusp;
      entry->genestrand = genestrand;
      entry->sensedir = sensedir;
      entry->innerp = innerp;
      entry->find_splices_p = find_splices_p;
      entry->salvagep = salvagep;
      entry->key_path =
	Path_copy(path,intlistpool,univcoordlistpool,
		  listpool,pathpool,vectorpool,transcriptpool,hitlistpool);
      entry->result_paths =
	qend_local_copy_added_paths(current_qend_paths,old_qend_paths,
				    intlistpool,univcoordlistpool,listpool,pathpool,
				    vectorpool,hitlistpool,transcriptpool);
      cache->qend_local_nentries++;
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
      cache->qend_local_inserts++;
      if (qend_local_cache_stats_print_p(cache->qend_local_inserts) == true) {
	qend_local_cache_stats_report(cache,"insert",path,depth);
      }
#endif
      return;
    }
  }

  return;
}

Path_segment_cache_T
Path_segment_cache_new (bool emitted_path_dedup_p) {
  Path_segment_cache_T new = (Path_segment_cache_T) CALLOC(1,sizeof(struct Path_segment_cache_T));

  new->emitted_path_dedup_p = emitted_path_dedup_p;
  return new;
}

Path_ext_anchor_cache_T
Path_ext_anchor_cache_new (void) {
#ifdef GSNAP_EXT_ANCHOR_PROOF
  ext_anchor_report_enabled_p = getenv("GSNAP_EXT_ANCHOR_REUSE_REPORT") != NULL;
  if (ext_anchor_report_enabled_p == true &&
      ext_anchor_report_registered_p == false) {
    atexit(ext_anchor_reuse_report);
    ext_anchor_report_registered_p = true;
  }
  return (Path_ext_anchor_cache_T) CALLOC(1,sizeof(struct Path_ext_anchor_cache_T));
#else
  return (Path_ext_anchor_cache_T) NULL;
#endif
}

void
Path_ext_anchor_cache_free (Path_ext_anchor_cache_T *old,
			    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool,
			    Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
#ifdef GSNAP_EXT_ANCHOR_PROOF
  int i;
  if (*old != NULL) {
    for (i = 0; i < EXT_ANCHOR_RESULT_CACHE_SIZE; i++) {
      if ((*old)->entries[i].occupiedp == true) {
	Path_gc(&(*old)->entries[i].sense_paths,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
	Path_gc(&(*old)->entries[i].antisense_paths,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
    FREE(*old);
    *old = (Path_ext_anchor_cache_T) NULL;
  }
#else
  (void) old; (void) intlistpool; (void) univcoordlistpool; (void) listpool;
  (void) pathpool; (void) transcriptpool; (void) hitlistpool;
#endif
}

void
Path_segment_cache_free (Path_segment_cache_T *old,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool,
			 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  int i;

  if (*old) {
#ifdef GSNAP_QEND_LOCAL_CACHE_STATS
    qend_local_cache_stats_report(*old,"free",(Path_T) NULL,/*depth*/-1);
#endif
    for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
      if ((*old)->entries[i].occupiedp == true) {
	Path_gc(&(*old)->entries[i].sense_paths,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
	Path_gc(&(*old)->entries[i].antisense_paths,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      }
      if ((*old)->qend_local_entries[i].occupiedp == true) {
	Path_free(&(*old)->qend_local_entries[i].key_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	Path_gc(&(*old)->qend_local_entries[i].result_paths,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
    for (i = 0; i < EMITTED_PATH_CACHE_SIZE; i++) {
      if ((*old)->emitted_path_entries[i].occupiedp == true) {
	Path_free(&(*old)->emitted_path_entries[i].key_path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
    }
    FREE(*old);
    *old = (Path_segment_cache_T) NULL;
  }

  return;
}

static bool
segment_cache_lookup (List_T *qstart_sense_paths, List_T *qstart_antisense_paths,
		     Path_segment_cache_T cache, Univdiag_T univdiag,
		     int querylength, int localdb_nmismatches_allowed,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     bool plusp, int genestrand, bool first_read_p, bool find_splices_p,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		     Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_segment_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return false;
  }

  hash = segment_cache_hash(univdiag,querylength,localdb_nmismatches_allowed,
			   chrnum,chroffset,chrhigh,plusp,genestrand,
			   first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == false) {
      return false;
    } else if (segment_cache_entry_matches_p(entry,univdiag,querylength,
					    localdb_nmismatches_allowed,chrnum,
					    chroffset,chrhigh,plusp,genestrand,
					    first_read_p,find_splices_p) == true) {
      *qstart_sense_paths =
	segment_cache_copy_paths(entry->sense_paths,intlistpool,univcoordlistpool,
				listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
      *qstart_antisense_paths =
	segment_cache_copy_paths(entry->antisense_paths,intlistpool,univcoordlistpool,
				listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
      return true;
    }
  }

  return false;
}

static void
segment_cache_insert (Path_segment_cache_T cache, Univdiag_T univdiag,
		     List_T qstart_sense_paths, List_T qstart_antisense_paths,
		     int querylength, int localdb_nmismatches_allowed,
		     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     bool plusp, int genestrand, bool first_read_p, bool find_splices_p,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		     Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_segment_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL || cache->qstart_nentries >= SEGMENT_CACHE_SIZE) {
    return;
  } else if (segment_cache_small_path_list_p(qstart_sense_paths,qstart_antisense_paths) == false) {
    return;
  }

  hash = segment_cache_hash(univdiag,querylength,localdb_nmismatches_allowed,
			   chrnum,chroffset,chrhigh,plusp,genestrand,
			   first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == true) {
      if (segment_cache_entry_matches_p(entry,univdiag,querylength,
				       localdb_nmismatches_allowed,chrnum,
				       chroffset,chrhigh,plusp,genestrand,
				       first_read_p,find_splices_p) == true) {
	return;
      }
    } else {
      entry->occupiedp = true;
      entry->qstart = univdiag->qstart;
      entry->qend = univdiag->qend;
      entry->nmismatches = univdiag->nmismatches;
      entry->querylength = querylength;
      entry->localdb_nmismatches_allowed = localdb_nmismatches_allowed;
      entry->chrnum = chrnum;
      entry->chroffset = chroffset;
      entry->chrhigh = chrhigh;
      entry->univdiagonal = univdiag->univdiagonal;
      entry->plusp = plusp;
      entry->genestrand = genestrand;
      entry->first_read_p = first_read_p;
      entry->find_splices_p = find_splices_p;
      entry->sense_paths =
	segment_cache_copy_paths(qstart_sense_paths,intlistpool,univcoordlistpool,
				listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
      entry->antisense_paths =
	segment_cache_copy_paths(qstart_antisense_paths,intlistpool,univcoordlistpool,
				listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
      cache->qstart_nentries++;
      return;
    }
  }

  return;
}

static unsigned int
segment_qend_cache_hash (Univdiag_T head, Univdiag_T next, int querylength,
			 int localdb_nmismatches_allowed, Chrnum_T chrnum,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool plusp, int genestrand, bool first_read_p,
			 bool find_splices_p) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) head->qstart);
  h = segment_cache_mix(h,(unsigned long long) head->qend);
  h = segment_cache_mix(h,(unsigned long long) head->nmismatches);
  h = segment_cache_mix(h,(unsigned long long) head->univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) next->qstart);
  h = segment_cache_mix(h,(unsigned long long) next->qend);
  h = segment_cache_mix(h,(unsigned long long) next->nmismatches);
  h = segment_cache_mix(h,(unsigned long long) next->univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) querylength);
  h = segment_cache_mix(h,(unsigned long long) localdb_nmismatches_allowed);
  h = segment_cache_mix(h,(unsigned long long) chrnum);
  h = segment_cache_mix(h,(unsigned long long) chroffset);
  h = segment_cache_mix(h,(unsigned long long) chrhigh);
  h = segment_cache_mix(h,(unsigned long long) plusp);
  h = segment_cache_mix(h,(unsigned long long) genestrand);
  h = segment_cache_mix(h,(unsigned long long) first_read_p);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);

  return (unsigned int) (h & SEGMENT_CACHE_MASK);
}

static bool
segment_qend_cache_entry_matches_p (struct Path_segment_qend_cache_entry_T *entry,
				    Univdiag_T head, Univdiag_T next,
				    int querylength, int localdb_nmismatches_allowed,
				    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				    bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  return entry->head_qstart == head->qstart &&
    entry->head_qend == head->qend &&
    entry->head_nmismatches == head->nmismatches &&
    entry->head_univdiagonal == head->univdiagonal &&
    entry->next_qstart == next->qstart &&
    entry->next_qend == next->qend &&
    entry->next_nmismatches == next->nmismatches &&
    entry->next_univdiagonal == next->univdiagonal &&
    entry->querylength == querylength &&
    entry->localdb_nmismatches_allowed == localdb_nmismatches_allowed &&
    entry->chrnum == chrnum &&
    entry->chroffset == chroffset &&
    entry->chrhigh == chrhigh &&
    entry->plusp == plusp &&
    entry->genestrand == genestrand &&
    entry->first_read_p == first_read_p &&
    entry->find_splices_p == find_splices_p;
}

static bool
segment_qend_cache_contains (Path_segment_cache_T cache, Univdiag_T head, Univdiag_T next,
			     int querylength, int localdb_nmismatches_allowed,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, bool first_read_p,
			     bool find_splices_p) {
  struct Path_segment_qend_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return false;
  }

  hash = segment_qend_cache_hash(head,next,querylength,localdb_nmismatches_allowed,
				 chrnum,chroffset,chrhigh,plusp,genestrand,
				 first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qend_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == false) {
      return false;
    } else if (segment_qend_cache_entry_matches_p(entry,head,next,querylength,
						  localdb_nmismatches_allowed,chrnum,
						  chroffset,chrhigh,plusp,genestrand,
						  first_read_p,find_splices_p) == true) {
      return true;
    }
  }

  return false;
}

bool
Path_segment_cache_qend_known_p (Path_segment_cache_T cache,
				 int head_qstart, int head_qend, int head_nmismatches,
				 Univcoord_T head_univdiagonal,
				 int next_qstart, int next_qend, int next_nmismatches,
				 Univcoord_T next_univdiagonal,
				 int querylength, int localdb_nmismatches_allowed,
				 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				 bool plusp, int genestrand, bool first_read_p,
				 bool find_splices_p) {
  struct Univdiag_T head, next;

  head.univdiagonal = head_univdiagonal;
  head.qstart = head_qstart;
  head.qend = head_qend;
  head.nmismatches = head_nmismatches;

  next.univdiagonal = next_univdiagonal;
  next.qstart = next_qstart;
  next.qend = next_qend;
  next.nmismatches = next_nmismatches;

  return segment_qend_cache_contains(cache,&head,&next,querylength,
				     localdb_nmismatches_allowed,chrnum,
				     chroffset,chrhigh,plusp,genestrand,
				     first_read_p,find_splices_p);
}

static void
segment_qend_cache_insert (Path_segment_cache_T cache, Univdiag_T head, Univdiag_T next,
			   int querylength, int localdb_nmismatches_allowed,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  struct Path_segment_qend_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL || cache->qend_nentries >= SEGMENT_CACHE_SIZE) {
    return;
  }

  hash = segment_qend_cache_hash(head,next,querylength,localdb_nmismatches_allowed,
				 chrnum,chroffset,chrhigh,plusp,genestrand,
				 first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qend_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == true) {
      if (segment_qend_cache_entry_matches_p(entry,head,next,querylength,
					     localdb_nmismatches_allowed,chrnum,
					     chroffset,chrhigh,plusp,genestrand,
					     first_read_p,find_splices_p) == true) {
	return;
      }
    } else {
      entry->occupiedp = true;
      entry->head_qstart = head->qstart;
      entry->head_qend = head->qend;
      entry->head_nmismatches = head->nmismatches;
      entry->head_univdiagonal = head->univdiagonal;
      entry->next_qstart = next->qstart;
      entry->next_qend = next->qend;
      entry->next_nmismatches = next->nmismatches;
      entry->next_univdiagonal = next->univdiagonal;
      entry->querylength = querylength;
      entry->localdb_nmismatches_allowed = localdb_nmismatches_allowed;
      entry->chrnum = chrnum;
      entry->chroffset = chroffset;
      entry->chrhigh = chrhigh;
      entry->plusp = plusp;
      entry->genestrand = genestrand;
      entry->first_read_p = first_read_p;
      entry->find_splices_p = find_splices_p;
      cache->qend_nentries++;
      return;
    }
  }

  return;
}


static unsigned int
segment_qstart_cache_hash (Univdiag_T head, Univdiag_T next, int querylength,
			   int localdb_nmismatches_allowed, Chrnum_T chrnum,
			   Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand, bool first_read_p,
			   bool find_splices_p) {
  unsigned long long h = 1469598103934665603ULL;

  h = segment_cache_mix(h,(unsigned long long) head->qstart);
  h = segment_cache_mix(h,(unsigned long long) head->qend);
  h = segment_cache_mix(h,(unsigned long long) head->nmismatches);
  h = segment_cache_mix(h,(unsigned long long) head->univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) next->qstart);
  h = segment_cache_mix(h,(unsigned long long) next->qend);
  h = segment_cache_mix(h,(unsigned long long) next->nmismatches);
  h = segment_cache_mix(h,(unsigned long long) next->univdiagonal);
  h = segment_cache_mix(h,(unsigned long long) querylength);
  h = segment_cache_mix(h,(unsigned long long) localdb_nmismatches_allowed);
  h = segment_cache_mix(h,(unsigned long long) chrnum);
  h = segment_cache_mix(h,(unsigned long long) chroffset);
  h = segment_cache_mix(h,(unsigned long long) chrhigh);
  h = segment_cache_mix(h,(unsigned long long) plusp);
  h = segment_cache_mix(h,(unsigned long long) genestrand);
  h = segment_cache_mix(h,(unsigned long long) first_read_p);
  h = segment_cache_mix(h,(unsigned long long) find_splices_p);

  return (unsigned int) (h & SEGMENT_CACHE_MASK);
}

static bool
segment_qstart_cache_entry_matches_p (struct Path_segment_qstart_cache_entry_T *entry,
				      Univdiag_T head, Univdiag_T next,
				      int querylength, int localdb_nmismatches_allowed,
				      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				      bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  return entry->head_qstart == head->qstart &&
    entry->head_qend == head->qend &&
    entry->head_nmismatches == head->nmismatches &&
    entry->head_univdiagonal == head->univdiagonal &&
    entry->next_qstart == next->qstart &&
    entry->next_qend == next->qend &&
    entry->next_nmismatches == next->nmismatches &&
    entry->next_univdiagonal == next->univdiagonal &&
    entry->querylength == querylength &&
    entry->localdb_nmismatches_allowed == localdb_nmismatches_allowed &&
    entry->chrnum == chrnum &&
    entry->chroffset == chroffset &&
    entry->chrhigh == chrhigh &&
    entry->plusp == plusp &&
    entry->genestrand == genestrand &&
    entry->first_read_p == first_read_p &&
    entry->find_splices_p == find_splices_p;
}

static bool
segment_qstart_cache_contains (Path_segment_cache_T cache, Univdiag_T head, Univdiag_T next,
			       int querylength, int localdb_nmismatches_allowed,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			       bool plusp, int genestrand, bool first_read_p,
			       bool find_splices_p) {
  struct Path_segment_qstart_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL) {
    return false;
  }

  hash = segment_qstart_cache_hash(head,next,querylength,localdb_nmismatches_allowed,
				   chrnum,chroffset,chrhigh,plusp,genestrand,
				   first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qstart_pair_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == false) {
      return false;
    } else if (segment_qstart_cache_entry_matches_p(entry,head,next,querylength,
						    localdb_nmismatches_allowed,chrnum,
						    chroffset,chrhigh,plusp,genestrand,
						    first_read_p,find_splices_p) == true) {
      return true;
    }
  }

  return false;
}

bool
Path_segment_cache_qstart_known_p (Path_segment_cache_T cache,
				   int head_qstart, int head_qend, int head_nmismatches,
				   Univcoord_T head_univdiagonal,
				   int next_qstart, int next_qend, int next_nmismatches,
				   Univcoord_T next_univdiagonal,
				   int querylength, int localdb_nmismatches_allowed,
				   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				   bool plusp, int genestrand, bool first_read_p,
				   bool find_splices_p) {
  struct Univdiag_T head, next;

  head.univdiagonal = head_univdiagonal;
  head.qstart = head_qstart;
  head.qend = head_qend;
  head.nmismatches = head_nmismatches;

  next.univdiagonal = next_univdiagonal;
  next.qstart = next_qstart;
  next.qend = next_qend;
  next.nmismatches = next_nmismatches;

  return segment_qstart_cache_contains(cache,&head,&next,querylength,
				       localdb_nmismatches_allowed,chrnum,
				       chroffset,chrhigh,plusp,genestrand,
				       first_read_p,find_splices_p);
}

static void
segment_qstart_cache_insert (Path_segment_cache_T cache, Univdiag_T head, Univdiag_T next,
			     int querylength, int localdb_nmismatches_allowed,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, bool first_read_p, bool find_splices_p) {
  struct Path_segment_qstart_cache_entry_T *entry;
  unsigned int hash;
  int i;

  if (cache == (Path_segment_cache_T) NULL || cache->qstart_pair_nentries >= SEGMENT_CACHE_SIZE) {
    return;
  }

  hash = segment_qstart_cache_hash(head,next,querylength,localdb_nmismatches_allowed,
				   chrnum,chroffset,chrhigh,plusp,genestrand,
				   first_read_p,find_splices_p);
  for (i = 0; i < SEGMENT_CACHE_SIZE; i++) {
    entry = &(cache->qstart_pair_entries[(hash + i) & SEGMENT_CACHE_MASK]);
    if (entry->occupiedp == true) {
      if (segment_qstart_cache_entry_matches_p(entry,head,next,querylength,
					       localdb_nmismatches_allowed,chrnum,
					       chroffset,chrhigh,plusp,genestrand,
					       first_read_p,find_splices_p) == true) {
	return;
      }
    } else {
      entry->occupiedp = true;
      entry->head_qstart = head->qstart;
      entry->head_qend = head->qend;
      entry->head_nmismatches = head->nmismatches;
      entry->head_univdiagonal = head->univdiagonal;
      entry->next_qstart = next->qstart;
      entry->next_qend = next->qend;
      entry->next_nmismatches = next->nmismatches;
      entry->next_univdiagonal = next->univdiagonal;
      entry->querylength = querylength;
      entry->localdb_nmismatches_allowed = localdb_nmismatches_allowed;
      entry->chrnum = chrnum;
      entry->chroffset = chroffset;
      entry->chrhigh = chrhigh;
      entry->plusp = plusp;
      entry->genestrand = genestrand;
      entry->first_read_p = first_read_p;
      entry->find_splices_p = find_splices_p;
      cache->qstart_pair_nentries++;
      return;
    }
  }

  return;
}



static bool *circularp;

static Transcriptome_T transcriptome;
static EF64_T chromosome_ef64;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Univcoord_T genomelength;

static int index1part;

static Localdb_T localdb;

static Chrpos_T min_intronlength;

/* Splicing */
static bool splicingp;
/* static bool novelsplicingp; */

static int max_insertionlen;
static int max_deletionlen;
static Chrpos_T shortsplicedist;
static int n_splice_paths = 1;
static int max_recursive_ends = 2; /* Needs to be 2 or more, because the case of 1 is already checked */
static int min_support_indel = 6;
static int indel_support_mismatch_penalty = 3;
static bool materialize_equiv_splice_alts_p = false;

static inline bool
sufficient_indel_support_p (int support, int nmismatches) {
  return support - indel_support_mismatch_penalty*nmismatches >= min_support_indel;
}

/* For splice plus indel */
/* static int max_splice_deletionlen = 3; */
/* static int max_splice_insertionlen = 3; */


#define T Path_T


/* Used with knownindels */
static T
attach_indel_qstart_simple (int adj, T path, int indel_pos,
			    Univcoord_T univdiagonal, int querylength, int try_sensedir,
			    bool plusp, int genestrand, int *mismatch_positions_alloc,
			    Stage1_T stage1, Knownsplicing_T knownsplicing, Spliceendsgen_T spliceendsgen,
			    Compress_T query_compress, char *queryptr,
			    Genomebits_T genomebits, Genomebits_T genomebits_alt,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, bool find_splices_p, 

			    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			    unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,

			    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Univdiagpool_T univdiagpool,
			    Vectorpool_T vectorpool) {
  T newpath;

  Univcoord_T distal_univdiagonal;
  int trimpos, nmismatches_to_trimpos, exon_origin;
  int found_sensedir;
  bool splice5p;
  Splicetype_T splicetype5;
  double ambig_prob_5;

  /* Subtract adj to get low diagonal for qstart, but add adj to get high diagonal for qend */
  distal_univdiagonal = univdiagonal - adj;

  /* Becauses splices are common, it is more likely that we have a
     splice distal to an indel than an indel distal to a splice */
  exon_origin = Path_exon_origin(path);
  debug13(printf("Calling Spliceends_qstart_trim with indel_pos %d\n",indel_pos));
  splice5p = Spliceends_qstart_trim(&trimpos,&nmismatches_to_trimpos,
				    &found_sensedir,&splicetype5,&ambig_prob_5,
				    stage1,knownsplicing,try_sensedir,
				    distal_univdiagonal,querylength,
				    /*pos3*/indel_pos,exon_origin,chrnum,chroffset,chrhigh,
				    plusp,genestrand,mismatch_positions_alloc,
				    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				    localdb_tally,localdb_overflowp,max_nblocks,
				    univdiagpool,vectorpool,spliceendsgen,
				    query_compress,queryptr,
				    genomebits,genomebits_alt,find_splices_p);
  debug13(printf("(1) Spliceends_qstart_trim returns trimpos %d and %d nmismatches, splice5 prob %f\n",
		 trimpos,nmismatches_to_trimpos,ambig_prob_5));

  if (trimpos == indel_pos) {
    debug13(printf("New indel does not have a good distal segment\n"));
    return (T) NULL;

  } else {
    newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    if (indel_pos == Intlist_head(newpath->endpoints)) {
      /* No change, so can keep nmismatches */
    } else if (indel_pos > Intlist_head(newpath->endpoints) && Intlist_head(newpath->nmismatches) == 0) {
      /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
      Intlist_head_set(newpath->endpoints,indel_pos);
    } else {
      /* Need to re-compute nmismatches */
      Intlist_head_set(newpath->endpoints,indel_pos);
      Intlist_head_set(newpath->nmismatches,-1); /* From previous endpoint to indel_pos */
      Intlist_head_set(newpath->ref_nmismatches,-1); /* From previous endpoint to indel_pos */
    }
    newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos
					  intlistpool_trace(__FILE__,__LINE__));

    newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_to_trimpos
					    intlistpool_trace(__FILE__,__LINE__)); /* from indel_pos to trimpos */
    newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_to_trimpos
						intlistpool_trace(__FILE__,__LINE__)); /* from indel_pos to trimpos */

    newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,distal_univdiagonal
						    univcoordlistpool_trace(__FILE__,__LINE__));

    if (adj < 0) {
      /* Start insertion */
      debug13(printf("Start insertion.  New diagonal is %u\n",distal_univdiagonal));
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_insertion(/*nindels*/-adj,pathpool)
					 listpool_trace(__FILE__,__LINE__));
    } else {
      /* Start deletion */
      debug13(printf("Start deletion.  New diagonal is %u\n",distal_univdiagonal));

      /* deletionpos = (univdiagonal - querylength) + indel_pos; */
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_deletion(/*nindels*/adj,pathpool)
					 listpool_trace(__FILE__,__LINE__));
    }

    debug13(Path_print(newpath));
    debug13(printf("\n"));

    assert(newpath != path);
    return newpath;
  }
}


/* Uses parts of attach_unknown_qstart */
/* Returns a newpath without modifying or deleting path */
static T
attach_indel_qstart (T path, Univcoord_T low_univdiagonal, int low_qstart,
		     Univcoord_T chroffset, Univcoord_T chrhigh, int querylength, Indelinfo_T indelinfo,
		     Compress_T query_compress, bool plusp, int genestrand,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool) {
  T newpath = NULL;
  Univcoord_T univdiagonal;
  int qend;
  int nindels, indel_pos;
  int supporti, supportj;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j;
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  int qstart;
#endif
  
  /* Do not need to call Spliceends_qstart_trim, because
     Genomebits_indel_solve_low has already computed trimpos, provided
     as low_qstart */

  univdiagonal = Univcoordlist_head(path->univdiagonals);

  /* Assume that left+qend gives a coordinate within genome */
  qend = Intlist_second_value(path->endpoints) /*+ ninserts*/;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  if (path->junctions == NULL) {
    /* ninserts = 0; */
  } else {
    /* ninserts = Junction_ninserts(List_head(path->junctions)); */
  }

  qstart = Intlist_head(path->endpoints) /*+ ninserts*/;
  printf("Entering attach_indel_qstart with low_univdiagonal %u with low_qstart %d and univdiagonal %u %d..%d (diff %d)\n",
	 low_univdiagonal,low_qstart,univdiagonal,qstart,qend,univdiagonal - low_univdiagonal);
#endif

  if (low_univdiagonal > univdiagonal + max_insertionlen) {
    /* Impossible */
    debug13(printf("Impossible\n"));

  } else if (low_univdiagonal + low_qstart < chroffset + querylength) {
    debug13(printf("Extends beyond start of chromosome: low_univdiagonal %u - querylength %d + low_qstart %d vs chroffset %u\n",
		   low_univdiagonal,querylength,low_qstart,chroffset));

  } else if (low_univdiagonal > univdiagonal) {
    /* (A) Insertion */
    nindels = low_univdiagonal - univdiagonal;
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal*/low_univdiagonal,/*indels*/+nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/low_qstart,/*pos3*/qend,querylength,
						    indelinfo,plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("(1) Insertion fails\n"));
      
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - (indel_pos + nindels);
      debug13(printf("(1) supporti %d - %d, supportj %d - (%d + %d)\n",
		     indel_pos,low_qstart,qend,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(1) Not enough support for indel: supporti %d and mismatches %d\n",
		       supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(1) Not enough support for indel: supportj %d and mismatches %d\n",
		       supportj,nmismatches_j));
      } else {
	newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	debug13(printf("(3) attach_unknown_qstart is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart
					      intlistpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,nmismatches_j);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_i
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(1) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }
    
  } else if (low_univdiagonal + max_deletionlen >= univdiagonal) {
    /* (B) Deletion (or short intron) */
    nindels = univdiagonal - low_univdiagonal;
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/low_univdiagonal,/*indels*/-nindels,chrhigh,
						   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/low_qstart,/*pos3*/qend,querylength,
						   indelinfo,plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("Deletion or short intron fails\n"));
	  
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - indel_pos;
      debug13(printf("(2) supporti %d - %d, supportj %d - %d\n",
		     indel_pos,low_qstart,qend,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(2) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(2) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	    
	assert(nindels >= 0);
	/* deletionpos = (low_univdiagonal - querylength) + indel_pos; */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	debug13(printf("(4) attach_unknown_qstart is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart
					      intlistpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,nmismatches_j);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_i
						    intlistpool_trace(__FILE__,__LINE__));
      
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("Deletion or short splice in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }
  }

  assert(newpath != path);
  return newpath;
}


#if 0
/* Using attach_unknown_qend instead, which can handle spliceindels */
/* Returns a newpath without modifying or deleting path */
static T
attach_splice_qstart (T path, int splice_qpos, int trimpos, int nmismatches_to_trimpos, Splicetype_T splicetype5,
		      bool plusp, int try_sensedir, double medial_prob, double distal_prob,
		      Univcoord_T medial_univdiagonal, Univcoord_T distal_univdiagonal,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, Pathpool_T pathpool, Univdiagpool_T univdiagpool,
		      Vectorpool_T vectorpool) {
  T newpath;
  /* int exon_origin; */
  bool splice5p = false;
  int found_sensedir = try_sensedir;
  /* Splicetype_T splicetype5; */
  double ambig_prob_5 = medial_prob;
  Chrpos_T splice_distance;

  /* exon_origin = Path_exon_origin(path); */
  debug13(printf("attach_splice_qstart with splice_qpos %d, trimpos %d\n",splice_qpos,trimpos));
#if 0
  splice5p = Spliceends_qstart_trim(&trimpos,&nmismatches_to_trimpos,&found_sensedir,&splicetype5,&ambig_prob_5,
				    stage1,knownsplicing,try_sensedir,
				    distal_univdiagonal,querylength,
				    /*pos3*/splice_qpos,exon_origin,chrnum,chroffset,
				    plusp,genestrand,mismatch_positions_alloc,
				    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				    localdb_tally,localdb_overflowp,max_nblocks,
				    univdiagpool,vectorpool,spliceendsgen,
				    query_compress,queryptr,
				    genomebits,genomebits_alt,find_splices_p);
  debug13(printf("(2) Spliceends_qstart_trim returns trimpos %d and %d nmismatches, splice5 prob %f\n",
		 trimpos,nmismatches_to_trimpos,ambig_prob_5));
#endif

  if (trimpos == splice_qpos) {
    debug13(printf("New splice does not have a good distal segment with trimpos %d == splice_qpos %d\n",
		   trimpos,splice_qpos));
    return (T) NULL;

  } else {
    newpath = Path_copy_5(path,splice5p,splicetype5,ambig_prob_5,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    debug13(printf("(1) attach_splice_qstart (splice) is modifying path %p\n",newpath));
    debug13(printf("Adding an intron on top of qstart diagonal with splice pos %d: %u to %u\n",
		   splice_qpos,medial_univdiagonal,distal_univdiagonal));
    debug13(printf("Before addition\n"));
    debug13(Path_print(newpath));
  
    if (splice_qpos == Intlist_head(newpath->endpoints)) {
      /* No change, so can keep nmismatches */
    } else if (splice_qpos > Intlist_head(newpath->endpoints) && Intlist_head(newpath->nmismatches) == 0) {
      /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
      Intlist_head_set(newpath->endpoints,splice_qpos);
    } else {
      /* Need to re-compute nmismatches */
      Intlist_head_set(newpath->endpoints,splice_qpos);
      Intlist_head_set(newpath->nmismatches,-1); /* From previous endpoint to splice_qpos */
      Intlist_head_set(newpath->ref_nmismatches,-1); /* From previous endpoint to splice_qpos */
    }
    newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos
					  intlistpool_trace(__FILE__,__LINE__));

    newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_to_trimpos
					    intlistpool_trace(__FILE__,__LINE__)); /* from splice_qpos to trimpos */
    newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_to_trimpos
						intlistpool_trace(__FILE__,__LINE__)); /* from splice_qpos to trimpos */
	  
    newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,distal_univdiagonal
						    univcoordlistpool_trace(__FILE__,__LINE__));

    splice_distance = (Chrpos_T) (medial_univdiagonal - distal_univdiagonal);
    if (plusp) {
      if (found_sensedir == SENSE_FORWARD) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/distal_prob,/*acceptor_prob*/medial_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else if (found_sensedir == SENSE_ANTI) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_ANTI,
									/*donor_prob*/medial_prob,/*acceptor_prob*/distal_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	fprintf(stderr,"Unexpected sensedir %d\n",found_sensedir);
	abort();
      }
    } else {
      if (found_sensedir == SENSE_ANTI) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_ANTI,
									/*donor_prob*/distal_prob,/*acceptor_prob*/medial_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else if (found_sensedir == SENSE_FORWARD) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/medial_prob,/*acceptor_prob*/distal_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	fprintf(stderr,"Unexpected sensedir %d\n",found_sensedir);
	abort();
      }
    }
	  
    newpath->splice5p = false;
    newpath->splicetype5 = NO_SPLICE;
    newpath->ambig_prob_5 = 0.0;
  
    debug13(printf("After addition\n"));
    debug13(Path_print(newpath));
  
    assert(newpath != path);
    return newpath;
  }
}
#endif


#if 0
static bool
extend_qstart_p (int *distal_matchlengths, int npartners, int common_splice_qpos) {
  int i;

  for (i = 0; i < npartners; i++) {
    if (common_splice_qpos - distal_matchlengths[i] <= 1) {
      return true;
    }
  }

  return false;
}
#endif



/* Implements the addition of a diagonal to the start of a path.
   Returns a Path_T object, or NULL if the diagonal cannot
   be attached.  All Path_T objects are copies of the original path */
/* Assume that left+low_qstart and left+low_qend give coordinates within genome */

/* Returns newpath, which is different from path */
static T
attach_unknown_qstart (T path, Univcoord_T low_univdiagonal, int low_qstart,
		       Univcoord_T chroffset, Univcoord_T chrhigh, int querylength,
		       Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		       Stage1_T stage1, Compress_T query_compress, char *queryptr,
		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		       Localdb_T localdb, int localdb_nmismatches_allowed,
		       bool plusp, int genestrand, bool innerp,
		       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, int try_sensedir, int middle_rank) {
  T newpath = NULL;
  Univcoord_T univdiagonal, low_left;
  int qend;

  Univcoordlist_T new_univdiagonals;
  Intlist_T new_nmismatches, new_ref_nmismatches, new_endpoints;
  List_T new_junctions;
  int ignore;

  int nindels, indel_pos;
  int supporti, supportj;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j, nmismatches, ref_nmismatches;
  int cont_trimpos, cont_nmismatches;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  int qstart;
#endif
  
#if 0
  if (low_univdiagonal + low_qstart < (Univcoord_T) querylength) {
    /* low segment is before the beginning of the genome */
    return (T) NULL;
  } else { 
    univdiagonal = Univcoordlist_head(path->univdiagonals);
  }
#else
  univdiagonal = Univcoordlist_head(path->univdiagonals);
#endif

  /* Assume that left+qend gives a coordinate within genome */
  qend = Intlist_second_value(path->endpoints) /*+ ninserts*/;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  if (path->junctions == NULL) {
    /* ninserts = 0; */
  } else {
    /* ninserts = Junction_ninserts(List_head(path->junctions)); */
  }

  qstart = Intlist_head(path->endpoints) /*+ ninserts*/;
  printf("Entering attach_unknown_qstart with try_sensedir %d, low_univdiagonal %u and univdiagonal %u %d..%d (diff %d)\n",
	 try_sensedir,low_univdiagonal - chroffset,univdiagonal - chroffset,qstart,qend,univdiagonal - low_univdiagonal);
#endif

  /* Compute continuation along univdiagonal for baseline */
  cont_trimpos = Genomebits_trim_qstart(&cont_nmismatches,query_compress,genomebits,
					univdiagonal,querylength,
					/*pos5*/0,/*pos3*/qend,plusp,genestrand);

  if (low_qstart >= qend) {
    debug13(printf("Does not add to start of path: low_qstart %d >= qend %d\n",low_qstart,qend));

  } else if (low_univdiagonal + low_qstart < chroffset + querylength) {
    debug13(printf("Extends beyond start of chromosome: low_univdiagonal %u - querylength %d + low_qstart %u vs chroffset %u\n",
		   low_univdiagonal,querylength,low_qstart,chroffset));

  } else if (low_univdiagonal == univdiagonal) {
    if (low_qstart >= Intlist_head(path->endpoints)) {
      debug13(printf("Mismatch fails, since new endpoint %d >= old endpoint %d\n",low_qstart,Intlist_head(path->endpoints)));
    } else {
      /* Mismatch: Revise the endpoint */
      debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(path->endpoints),low_qstart));
      newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
      
      /* Determine nmismatches */
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,querylength,/*pos5*/low_qstart,/*pos3*/qend,
							  plusp,genestrand);
      /* debug13(printf("Counting mismatches from %d to %d => %d (%d ref)\n",
	 low_qstart,qend,nmismatches,ref_nmismatches)); */

      debug13(printf("(2) attach_unknown_qstart is modifying path %p\n",newpath));
      Intlist_head_set(newpath->nmismatches,nmismatches);
      Intlist_head_set(newpath->ref_nmismatches,nmismatches);
      Intlist_head_set(newpath->endpoints,low_qstart);
    }
    
#if 0
  } else if (low_univdiagonal > univdiagonal + max_insertionlen) {
    /* Possible with segment search */
    debug13(printf("Impossible\n"));
#endif

  } else if (low_univdiagonal > univdiagonal) {
    /* (A) Insertion */
    nindels = low_univdiagonal - univdiagonal;
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal*/low_univdiagonal,/*indels*/+nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/low_qstart,/*pos3*/qend,querylength,
						    indelinfo,plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("(2) Insertion fails\n"));
      
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - (indel_pos + nindels);
      debug13(printf("(3) supporti %d - %d, supportj %d - (%d + %d)\n",
		     indel_pos,low_qstart,qend,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(3) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(3) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	debug13(printf("(3) attach_unknown_qstart is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart
					      intlistpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,nmismatches_j);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_i
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(2) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }
    
  } else if (low_univdiagonal + max_deletionlen >= univdiagonal) {
    /* (B) Deletion (or short intron) */
    nindels = univdiagonal - low_univdiagonal;
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/low_univdiagonal,/*indels*/-nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/low_qstart,/*pos3*/qend,querylength,
						   indelinfo,plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("Deletion or short intron fails\n"));
	  
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - indel_pos;
      debug13(printf("(4) supporti %d - %d, supportj %d - %d\n",
		     indel_pos,low_qstart,qend,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(4) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(4) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	    
	assert(nindels >= 0);
	/* deletionpos = (low_univdiagonal - querylength) + indel_pos; */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	      
	debug13(printf("(4) attach_unknown_qstart is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,low_qstart
					      intlistpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(newpath->nmismatches,nmismatches_j);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_i
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("Deletion or short splice in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* No splicing on circular chromosomes */
    debug13(printf("No splicing on circular chromosomes\n"));
#endif

  } else if (low_univdiagonal + shortsplicedist < univdiagonal) {
    /* Gap is too large for a splice */
    debug13(printf("Gap from %u to %u is too large for a splice\n",low_univdiagonal,univdiagonal));

  } else if (splicingp == false) {
    /* Unable to try splicing */

  } else {
    /* (C) Splice with or without indel */
    low_left = low_univdiagonal - (Univcoord_T) querylength;
    /* left = univdiagonal - (Univcoord_T) querylength; */

    /* Previously filled spliceinfo with knownsplicing information
       on donors and antiacceptors for segmenti and segmentj, now
       done in Splice_resolve */
    
    debug13(printf("BEFORE call to Splice_resolve_qstart, %d..%d\n",low_qstart,qend));
    debug13(Path_print(path));

    if (Splice_resolve_qstart(&new_univdiagonals,&new_nmismatches,&new_ref_nmismatches,&new_endpoints,&new_junctions,
			      /*univdiagonal_i*/low_univdiagonal,/*univdiagonal_j*/univdiagonal,
			      stage1,query_compress,queryptr,plusp,chroffset,chrhigh,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,localdb,
			      localdb_nmismatches_allowed,
			      /*pos5*/low_qstart,/*pos3*/qend,querylength,cont_trimpos,cont_nmismatches,
			      indelinfo,spliceinfo,knownsplicing,univcoordlistpool,intlistpool,listpool,pathpool,
			      /*sensedir*/try_sensedir,genestrand,/*trim5p*/true,/*trim3p*/false,innerp,
			      /*allow_indels_p*/true,/*find_middle_p*/true,middle_rank) == 0) {
      debug13(printf("Splice_resolve_qstart: fails\n"));
      newpath = (T) NULL;

    } else {
      newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

      newpath->univdiagonals = Univcoordlist_append(new_univdiagonals,newpath->univdiagonals);
      newpath->nmismatches = Intlistpool_pop(newpath->nmismatches,intlistpool,&ignore
					     intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlist_append(new_nmismatches,newpath->nmismatches);
      newpath->ref_nmismatches = Intlistpool_pop(newpath->ref_nmismatches,intlistpool,&ignore
						 intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlist_append(new_ref_nmismatches,newpath->ref_nmismatches);
      newpath->endpoints = Intlistpool_pop(newpath->endpoints,intlistpool,&ignore
					   intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlist_append(new_endpoints,newpath->endpoints);
      newpath->junctions = List_append(new_junctions,newpath->junctions);
	
      debug13(printf("AFTER call to Splice_resolve_qstart\n"));
      debug13(Path_print(newpath));
    }
  }

  debug13(printf("attach_unknown_qstart is returning newpath %p\n",newpath));
  debug13(Path_print(newpath));

  assert(newpath != path);
  return newpath;
}


/* Adds a diagonal to the start of a path, as found by
   compute_qstart_local.  Calls attach_unknown_qstart to modify the path,
   potentially leading to multiple paths */
static List_T
multiadd_splice_qstarts (T path, Univcoord_T univdiagonal, Splicetype_T splicetype,
			 bool boundedp, int splice_qpos, int anchor_qpos, int medial_nmismatches,
			 double medial_prob, int medial_support, Univcoord_T *distal_splice_positions,
			 int *distal_support, Univdiag_T *distal_univdiags, double *distal_probs, int npartners,

			 Stage1_T stage1, Knownsplicing_T knownsplicing, Compress_T query_compress, char *queryptr,

			 Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			 unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			 Localdb_T localdb, int localdb_nmismatches_allowed,
			 bool plusp, int querylength, bool innerp,

			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Hitlistpool_T hitlistpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			 Vectorpool_T vectorpool, int found_sensedir) {
  List_T newpaths;
  Altsplice_T qstart_alts;
  T newpath, nextpath;
  Univcoord_T distal_univdiagonal;
  int best_endpoint, distal_trimpos;
  int n_good_partners, goodi, i, middle_rank;

  debug13(printf("Entered multiadd_splice_qstarts with %d partners at univdiagonal %u, common_splice_qpos %d.  Current path: ",
		 npartners,univdiagonal,splice_qpos));
  debug13(Path_print(path));
  debug13(printf("\n"));

  if (npartners == 1) {
    /* Needed to find spliceindel, where distal_prob is not known yet */
    debug13(printf("Only one partner: position %u\n",distal_splice_positions[0]));
    n_good_partners = 1;
    goodi = 0;

  } else {
    n_good_partners = 0;
    for (i = 0; i < npartners; i++) {
      if (distal_univdiags[i]->univdiagonal < univdiagonal &&
	  Splice_accept_p(splice_qpos,querylength,distal_probs[i],medial_prob,
			  /*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			  /*univdiagonal_j*/univdiagonal,
			  distal_support[i],distal_univdiags[i]->nmismatches,
			  medial_support,medial_nmismatches,
			  plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
			  innerp) == true) {
	debug13(printf("position %d:%d is good\n",i,distal_splice_positions[i]));
	n_good_partners += 1;
	goodi = i;
      }
    }
    debug13(printf("%d partners are good\n",n_good_partners));
  }

  if (materialize_equiv_splice_alts_p == true && n_good_partners > 1) {
    bool have_equiv = false, equivalentp = true;
    Univcoord_T equiv_univdiagonal = 0;
    int equiv_qstart = 0, equiv_qend = 0, equiv_nmismatches = 0, equiv_support = 0;
    double best_prob = -1.0;

    for (i = 0; i < npartners; i++) {
      if (distal_univdiags[i]->univdiagonal < univdiagonal &&
	  Splice_accept_p(splice_qpos,querylength,distal_probs[i],medial_prob,
			  /*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			  /*univdiagonal_j*/univdiagonal,
			  distal_support[i],distal_univdiags[i]->nmismatches,
			  medial_support,medial_nmismatches,
			  plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
			  innerp) == true) {
	if (have_equiv == false) {
	  have_equiv = true;
	  equiv_univdiagonal = distal_univdiags[i]->univdiagonal;
	  equiv_qstart = distal_univdiags[i]->qstart;
	  equiv_qend = distal_univdiags[i]->qend;
	  equiv_nmismatches = distal_univdiags[i]->nmismatches;
	  equiv_support = distal_support[i];
	} else if (distal_univdiags[i]->univdiagonal != equiv_univdiagonal ||
		   distal_univdiags[i]->qstart != equiv_qstart ||
		   distal_univdiags[i]->qend != equiv_qend ||
		   distal_univdiags[i]->nmismatches != equiv_nmismatches ||
		   distal_support[i] != equiv_support) {
	  equivalentp = false;
	}

	if (distal_probs[i] > best_prob) {
	  best_prob = distal_probs[i];
	  goodi = i;
	}
      }
    }

    if (have_equiv == true && equivalentp == true) {
      debug13(printf("Collapsing %d equivalent qstart splice partners at univdiagonal %u, %d..%d\n",
		     n_good_partners,equiv_univdiagonal,equiv_qstart,equiv_qend));
      n_good_partners = 1;
    }
  }

  if (n_good_partners == 0) {
    /* Multiple partners, but none good */
    return (List_T) NULL;

  } else if (n_good_partners == 1 /*|| extend_qstart_p(distal_matchlengths,npartners,common_splice_qpos) == false*/) {

    newpaths = (List_T) NULL;
    assert(distal_splice_positions[goodi] != 0);
#if 0
    distal_univdiagonal = distal_splice_positions[goodi] - splice_qpos + querylength;
#else
    distal_univdiagonal = distal_univdiags[goodi]->univdiagonal;
    distal_trimpos = distal_univdiags[goodi]->qstart;
#endif
    if (distal_univdiagonal == univdiagonal) {
      /* Not a splice */
      debug13(printf("Continuing the alignment to the start\n"));
      newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
      Intlist_head_set(newpath->endpoints,distal_trimpos);
      Intlist_head_set(newpath->nmismatches,-1);
      Intlist_head_set(newpath->ref_nmismatches,-1);
      
      debug13(Path_print(newpath));
      newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
			      hitlistpool_trace(__FILE__,__LINE__));
      
#if 0
    } else if (0) {
      /* Previously called attach_splice_qstart, but that did not handle spliceindels */
#endif

    } else {
      for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	if ((newpath = attach_unknown_qstart(path,/*low_univdiagonal*/distal_univdiagonal,
						  /*low_qstart*/distal_trimpos,
						  path->chroffset,path->chrhigh,querylength,
						  stage1->indelinfo,stage1->spliceinfo,knownsplicing,
						  stage1,query_compress,queryptr,
						  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						  localdb_tally,localdb_overflowp,max_nblocks,localdb,
						  localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
						  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
						  vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,middle_rank)) == NULL) {
	  break;
	} else {
	  assert(newpath != path);
	  newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
				  hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      newpaths = List_reverse(newpaths);
    }
      
    return newpaths;

  } else if (innerp == false &&
	     (distal_univdiagonal =
	      Altsplice_qstart_best(&distal_trimpos,splice_qpos,querylength,
				    /*medial_univdiagonal*/univdiagonal,
				    medial_prob,medial_support,medial_nmismatches,distal_support,
				    distal_univdiags,distal_probs,npartners,
				    plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
				    innerp)) > 0) {
    debug13(printf("Found a best outer splice\n"));

    newpaths = (List_T) NULL;
    for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
      if ((newpath = attach_unknown_qstart(path,/*low_univdiagonal*/distal_univdiagonal,
					       /*low_qstart*/distal_trimpos,
					       path->chroffset,path->chrhigh,querylength,
					       stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					       stage1,query_compress,queryptr,
					       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					       localdb_tally,localdb_overflowp,max_nblocks,localdb,
					       localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					       vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,middle_rank)) == NULL) {
	break;
      } else {
	assert(newpath != path);
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
    return List_reverse(newpaths);


  } else if ((qstart_alts =
	      Altsplice_qstart_new(&best_endpoint,boundedp,splice_qpos,anchor_qpos,querylength,
				   /*medial_univdiagonal*/univdiagonal,medial_nmismatches,
				   medial_prob,medial_support,distal_splice_positions,
				   distal_support,distal_univdiags,distal_probs,
				   /*order_medial_qpos*/NULL,npartners,
				   plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,innerp,
				   pathpool,vectorpool,/*sort_bydistal_p*/false)) == NULL) {
    /* Multiple partners do not work */
    return (List_T) NULL;
    
  } else {
    /* Multiple partners extend to the start, so create an alt substring */
    /* Altsplice_qstart_new needs to filter for good partners */
    debug13(printf("Creating an alt substring\n"));
    newpath = Path_copy_5(path,/*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    debug13(printf("Combining qstart alts for path:\n"));
    /* No need to sort because we are using a common splice qpos */
    newpath->qstart_alts = qstart_alts;

    Intlist_head_set(newpath->endpoints,best_endpoint);
    Intlist_head_set(newpath->nmismatches,-1);
    Intlist_head_set(newpath->ref_nmismatches,-1);
    
    debug13(printf("(2) Resulting path: "));
    debug13(Path_print(newpath));
    debug13(printf("add_qstart_local is returning path %p as newpaths\n",newpath));
    
    assert(newpath != path);

    if (innerp == true) {
      /* Allow Altsplice_resolve procedures to select later */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((goodi = Altsplice_select_qstart(newpath->qstart_alts,querylength,
						query_compress,plusp,/*genestrand*/0)) < 0) {
      /* Keep outer altsplice */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((nextpath = attach_unknown_qstart(path,/*low_univdiagonal*/distal_univdiags[goodi]->univdiagonal,
						 /*low_qstart*/distal_univdiags[goodi]->qstart,
						 path->chroffset,path->chrhigh,querylength,
						 stage1->indelinfo,stage1->spliceinfo,knownsplicing,
						 stage1,query_compress,queryptr,
						 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						 localdb_tally,localdb_overflowp,max_nblocks,localdb,
						 localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
						 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
						 vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,/*middle_rank*/0)) == NULL) {
      /* Keep outer altsplice */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else {
      /* Use solved outer splice */
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      debug13(Path_print(nextpath));
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(nextpath)
			  hitlistpool_trace(__FILE__,__LINE__));
    }
  }
}


#if 0
static void
check_for_ascending_values (Univcoord_T *diagonals, int ndiagonals) {
  int i;

  for (i = 0; i < ndiagonals - 1; i++) {
    if (diagonals[i+1] < diagonals[i]) {
      /* abort(); */
    }
  }
  return;
}
#endif


/* Used with knownindels */
static T
attach_indel_qend_simple (int adj, T path, int indel_pos,
			  Univcoord_T univdiagonal, int querylength, int try_sensedir,
			  bool plusp, int genestrand, int *mismatch_positions_alloc,
			  Stage1_T stage1, Knownsplicing_T knownsplicing, Spliceendsgen_T spliceendsgen,
			  Compress_T query_compress, char *queryptr,
			  Genomebits_T genomebits, Genomebits_T genomebits_alt,
			  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, bool find_splices_p,

			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,

			  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Univdiagpool_T univdiagpool,
			  Vectorpool_T vectorpool) {
  T newpath;

  Univcoord_T distal_univdiagonal;
  int trimpos, nmismatches_to_trimpos, exon_origin;
  int found_sensedir;
  bool splice3p;
  Splicetype_T splicetype3;
  double ambig_prob_3;

  /* Subtract adj to get low diagonal for qstart, but add adj to get high diagonal for qend */
  distal_univdiagonal = univdiagonal + adj;

  /* Becauses splices are common, it is more likely that we have a
     splice distal to an indel than an indel distal to a splice */
  exon_origin = Path_exon_origin(path);
  debug13(printf("Calling Spliceends_qend_trim with indel_pos %d\n",indel_pos));
  splice3p = Spliceends_qend_trim(&trimpos,&nmismatches_to_trimpos,
				  &found_sensedir,&splicetype3,&ambig_prob_3,
				  stage1,knownsplicing,try_sensedir,
				  distal_univdiagonal,querylength,
				  /*pos5*/indel_pos,exon_origin,chrnum,chroffset,chrhigh,
				  plusp,genestrand,mismatch_positions_alloc,
				  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				  localdb_tally,localdb_overflowp,max_nblocks,
				  univdiagpool,vectorpool,spliceendsgen,
				  query_compress,queryptr,
				  genomebits,genomebits_alt,find_splices_p);
  debug13(printf("(1) Spliceends_qend_trim returns trimpos %d and %d nmismatches, splice3 prob %f\n",
		 trimpos,nmismatches_to_trimpos,ambig_prob_3));

  if (trimpos == indel_pos) {
    debug13(printf("New indel does not have a good distal segment\n"));
    return (T) NULL;

  } else {
    newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    if (indel_pos == Intlist_head(newpath->endpoints)) {
      /* No change, so can keep nmismatches */
    } else if (indel_pos > Intlist_head(newpath->endpoints) && Intlist_head(newpath->nmismatches) == 0) {
      /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
      Intlist_head_set(newpath->endpoints,indel_pos);
    } else {
      /* Need to re-compute nmismatches */
      Intlist_head_set(newpath->endpoints,indel_pos);
      Intlist_head_set(newpath->nmismatches,-1); /* From previous endpoint to indel_pos */
      Intlist_head_set(newpath->ref_nmismatches,-1); /* From previous endpoint to indel_pos */
    }
    newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos
					  intlistpool_trace(__FILE__,__LINE__));

    newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_to_trimpos
					    intlistpool_trace(__FILE__,__LINE__)); /* From indel_pos to trimpos */
    newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_to_trimpos
						intlistpool_trace(__FILE__,__LINE__));/* From indel_pos to trimpos */

    newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,distal_univdiagonal
						    univcoordlistpool_trace(__FILE__,__LINE__));


    /* Subtract adj to get low diagonal, but add adj to get high diagonal */
    /* For qend, push i first, then push j.  For qend, add adj. */
    if (adj < 0) {
      /* End insertion */
      debug13(printf("End insertion.  New diagonal is %u\n",distal_univdiagonal));
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_insertion(/*nindels*/-adj,pathpool)
					 listpool_trace(__FILE__,__LINE__));
    } else {
      /* End deletion */
      debug13(printf("End deletion.  New diagonal is %u\n",distal_univdiagonal));
      /* deletionpos = (univdiagonal - querylength) + indel_pos; */
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_deletion(/*nindels*/adj,pathpool)
					 listpool_trace(__FILE__,__LINE__));
    }

    debug13(Path_print(newpath));
    debug13(printf("\n"));

    assert(newpath != path);
    return newpath;
  }
}



/* Uses code from attach_unknown_qend */
static T
attach_indel_qend (T path, Univcoord_T high_univdiagonal, int high_qend,
		   Univcoord_T chrhigh, int querylength, Indelinfo_T indelinfo,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool) {
  T newpath = NULL;
  Univcoord_T univdiagonal;
  int qstart, ninserts;
  int nindels, indel_pos;
  int supporti, supportj;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j;
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  int qend;
#endif


  /* Do not need to call Spliceends_qend_trim, because
     Genomebits_indel_solve_high has already computed trimpos,
     provided as high_qend*/

  univdiagonal = Univcoordlist_head(path->univdiagonals);

#if 0
  ninserts = Junction_total_ninserts(path->junctions);
#else
  if (path->junctions == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts(List_head(path->junctions));
  }
#endif

  /* Assume that left+qstart gives a coordinate within genome */
  qstart = Intlist_second_value(path->endpoints) + ninserts;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  qend = Intlist_head(path->endpoints) /*+ ninserts*/;
  printf("Entering attach_indel_qend with univdiagonal %u %d..%d and high_univdiagonal %u with high_qend %d (diff %d)\n",
	 univdiagonal,qstart,qend,high_univdiagonal,high_univdiagonal - univdiagonal,high_qend);
#endif

  if (high_univdiagonal + max_insertionlen < univdiagonal) {
    /* Impossible */
    debug13(printf("Impossible\n"));

  } else if (high_univdiagonal + high_qend >= chrhigh + querylength) {
    debug13(printf("Extends beyond end of chromosome: high_univdiagonal %u - querylength %d + high_qend %d vs chrhigh %u\n",
		   high_univdiagonal,querylength,high_qend,chrhigh));

  } else if (high_univdiagonal < univdiagonal) {
    /* (A) Insertion */
    nindels = univdiagonal - high_univdiagonal;
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    univdiagonal,/*indels*/+nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/qstart,/*pos3*/high_qend,querylength,
						    indelinfo,plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("(3) Insertion fails\n"));

    } else {
      supporti = indel_pos - qstart;
      supportj = high_qend - (indel_pos + nindels);
      debug13(printf("(5) supporti %d - %d, supportj %d - (%d + %d)\n",
		     indel_pos,qstart,high_qend,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(5) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(5) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	debug13(printf("(3) attach_indel_qend is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend
					    intlistpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,nmismatches_i);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_j
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(3) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
	
	debug13(Path_print(newpath));
      }
    }

  } else if (high_univdiagonal <= univdiagonal + max_deletionlen) {
    /* (B) Deletion (or short intron) */
    nindels = high_univdiagonal - univdiagonal;
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/univdiagonal,/*indels*/-nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/qstart,/*pos3*/high_qend,querylength,
						   indelinfo,plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("Deletion or short intron fails\n"));
      
    } else {
      supporti = indel_pos - qstart;
      supportj = high_qend - indel_pos;
      debug13(printf("(6) supporti %d - %d, supportj %d - %d\n",
		     indel_pos,qstart,high_qend,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(6) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(6) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	assert(nindels >= 0);
	/* deletionpos = (univdiagonal - querylength) + indel_pos; */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	debug13(printf("(4) attach_indel_qend is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend
					      intlistpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,nmismatches_i);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_j
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(1) Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,nmismatches_i,nmismatches_j));

	debug13(Path_print(newpath));
      }
    }
  }

  assert(newpath != path);
  return newpath;
}


#if 0
/* Using attach_unknown_qend instead, which can handle spliceindels */
/* Returns a new path without modifying deleting path */
static T
attach_splice_qend (T path, int splice_qpos, int trimpos, int nmismatches_to_trimpos, Splicetype_T splicetype3,
		    bool plusp, int try_sensedir, double medial_prob, double distal_prob,
		    Univcoord_T medial_univdiagonal, Univcoord_T distal_univdiagonal,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Univdiagpool_T univdiagpool,
		    Vectorpool_T vectorpool) {
  T newpath;
  /* int exon_origin; */
  bool splice3p = false;
  int found_sensedir = try_sensedir;
  /* Splicetype_T splicetype3; */
  double ambig_prob_3 = medial_prob;
  Chrpos_T splice_distance;

  /* exon_origin = Path_exon_origin(path); */
  debug13(printf("Calling attach_splice_qend with splice_qpos %d, trimpos %d\n",splice_qpos,trimpos));
#if 0
  splice3p = Spliceends_qend_trim(&trimpos,&nmismatches_to_trimpos,
				  &found_sensedir,&splicetype3,&ambig_prob_3,
				  stage1,knownsplicing,try_sensedir,
				  distal_univdiagonal,querylength,
				  /*pos5*/splice_qpos,exon_origin,chrnum,chroffset,chrhigh,
				  plusp,genestrand,mismatch_positions_alloc,
				  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				  localdb_tally,localdb_overflowp,max_nblocks,
				  univdiagpool,vectorpool,spliceendsgen,
				  query_compress,queryptr,
				  genomebits,genomebits_alt,find_splices_p);
  debug13(printf("(2) Spliceends_qend_trim returns trimpos %d and %d nmismatches, splice3 prob %f\n",
		 trimpos,nmismatches_to_trimpos,ambig_prob_3));
#endif

  if (trimpos == splice_qpos) {
    debug13(printf("New splice does not have a good distal segment with trimpos %d == splice_qpos %d\n",
		   trimpos,splice_qpos));
    return (T) NULL;

  } else {
    newpath = Path_copy_3(path,splice3p,splicetype3,ambig_prob_3,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    debug13(printf("(1) attach_splice_qend (splice) is modifying path %p\n",path));
    debug13(printf("Adding an intron on top of qend diagonal with splice pos %d: %u to %u\n",
		   splice_qpos,medial_univdiagonal,distal_univdiagonal));
    debug13(printf("Before addition\n"));
    debug13(Path_print(newpath));

    if (splice_qpos == Intlist_head(newpath->endpoints)) {
      /* No change, so can keep nmismatches */
    } else if (splice_qpos < Intlist_head(newpath->endpoints) && Intlist_head(newpath->nmismatches) == 0) {
      /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
      Intlist_head_set(newpath->endpoints,splice_qpos);
    } else {
      /* Need to re-compute nmismatches */
      Intlist_head_set(newpath->endpoints,splice_qpos);
      Intlist_head_set(newpath->nmismatches,-1); /* From previous endpoint to splice_qpos */
      Intlist_head_set(newpath->ref_nmismatches,-1); /* From previous endpoint to splice_qpos */
    }
    newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos
					  intlistpool_trace(__FILE__,__LINE__));

    newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_to_trimpos
					    intlistpool_trace(__FILE__,__LINE__)); /* From splice_qpos to trimpos */
    newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_to_trimpos
						intlistpool_trace(__FILE__,__LINE__)); /* From splice_qpos to trimpos */
       
    newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,distal_univdiagonal
						    univcoordlistpool_trace(__FILE__,__LINE__));

    splice_distance = (Chrpos_T) (distal_univdiagonal - medial_univdiagonal);
    if (plusp) {
      if (found_sensedir == SENSE_FORWARD) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/medial_prob,/*acceptor_prob*/distal_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else if (found_sensedir == SENSE_ANTI) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_ANTI,
									/*donor_prob*/distal_prob,/*acceptor_prob*/medial_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	fprintf(stderr,"Unexpected sensedir %d\n",found_sensedir);
	abort();
      }
    } else {
      if (found_sensedir == SENSE_ANTI) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_ANTI,
									/*donor_prob*/medial_prob,/*acceptor_prob*/distal_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else if (found_sensedir == SENSE_FORWARD) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
									/*donor_prob*/distal_prob,/*acceptor_prob*/medial_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	fprintf(stderr,"Unexpected sensedir %d\n",found_sensedir);
	abort();
      }
    }    
	  
    newpath->splice3p = false;
    newpath->splicetype3 = NO_SPLICE;
    newpath->ambig_prob_3 = 0.0;

    debug13(printf("After addition\n"));
    debug13(Path_print(newpath));

    assert(newpath != path);
    return newpath;
  }
}
#endif


#if 0
static bool
extend_qend_p (int *distal_matchlengths, int npartners, int common_splice_qpos, int querylength) {
  int i;

  for (i = 0; i < npartners; i++) {
    if (common_splice_qpos + distal_matchlengths[i] >= querylength - 1) {
      return true;
    }
  }

  return false;
}
#endif


/* Implements the addition of a diagonal to the end of a path.
   Returns a Path_T object, or NULL if the diagonal cannot
   be attached.  All Path_T objects are copies of the original path */

/* Returns newpath, which is different from path */
static T
attach_unknown_qend (T path, Univcoord_T high_univdiagonal, int high_qend,
		     Univcoord_T chroffset, Univcoord_T chrhigh, int querylength,
		     Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     Stage1_T stage1, Compress_T query_compress, char *queryptr,
		     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		     unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		     Localdb_T localdb, int localdb_nmismatches_allowed,
		     bool plusp, int genestrand, bool innerp,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		     Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, int try_sensedir, int middle_rank) {
  T newpath = NULL;
  Univcoord_T univdiagonal, left;
  int qstart, ninserts;

  Univcoordlist_T new_univdiagonals;
  Intlist_T new_nmismatches, new_ref_nmismatches, new_endpoints;
  List_T new_junctions;
  int ignore;

  int nindels, indel_pos;
  int supporti, supportj;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j, nmismatches, ref_nmismatches;
  int cont_trimpos, cont_nmismatches;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  int qend;
#endif


#if 0
  if (high_univdiagonal + high_qstart < (Univcoord_T) querylength) {
    /* high segment is before the beginning of the genome */
    return (T) NULL;
  } else {
    univdiagonal = Univcoordlist_head(path->univdiagonals);
  }
#else
  univdiagonal = Univcoordlist_head(path->univdiagonals);
#endif

  /* Should not need to check for high_univdiagonal >= genomelength, since that should be true */


#if 0
  ninserts = Junction_total_ninserts(path->junctions);
#else
  if (path->junctions == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts(List_head(path->junctions));
  }
#endif

  /* Assume that left+qstart gives a coordinate within genome */
  qstart = Intlist_second_value(path->endpoints) + ninserts;

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  qend = Intlist_head(path->endpoints) /*+ ninserts*/;
  printf("Entering attach_unknown_qend with try_sensedir %d, univdiagonal %u %d..%d and high_univdiagonal %u (diff %d)\n",
	 try_sensedir,univdiagonal - chroffset,qstart,qend,high_univdiagonal - chroffset,high_univdiagonal - univdiagonal);
#endif

  /* Compute continuation along univdiagonal for baseline */
  cont_trimpos = Genomebits_trim_qend(&cont_nmismatches,query_compress,genomebits,
				      univdiagonal,querylength,
				      /*pos5*/qstart,/*pos3*/querylength,plusp,genestrand);

  if (qstart >= high_qend) {
    debug13(printf("Does not add to end of path: qstart %d >= high_qend %d\n",qstart,high_qend));
    
  } else if (high_univdiagonal + high_qend >= chrhigh + querylength) {
    debug13(printf("Extends beyond end of chromosome: high_univdiagonal %u - querylength %d + high_qend %d vs chrhigh %u\n",
		   high_univdiagonal,querylength,high_qend,chrhigh));

  } else if (high_univdiagonal == univdiagonal) {
    if (high_qend <= Intlist_head(path->endpoints)) {
      debug13(printf("Mismatch fails, since new endpoint %d <= old endpoint %d\n",high_qend,Intlist_head(path->endpoints)));

    } else {
      /* Mismatch: Revise the endpoint */
      debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(path->endpoints),high_qend));
      newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
      
      /* Determine nmismatches */
      if (path->junctions == NULL) {
	ninserts = 0;
      } else {
	ninserts = Junction_ninserts(List_head(path->junctions));
      }
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,querylength,
							  /*pos5*/Intlist_second_value(newpath->endpoints) + ninserts,
							  /*pos3*/high_qend,plusp,genestrand);
      /* debug13(printf("Counting mismatches from %d to %d => %d (%d ref)\n",
	 Intlist_head(Intlist_next(newpath->endpoints)),high_qend,nmismatches,ref_nmismatches)); */

      debug13(printf("(2) attach_unknown_qend is modifying path %p\n",newpath));
      Intlist_head_set(newpath->nmismatches,nmismatches);
      Intlist_head_set(newpath->ref_nmismatches,nmismatches);
      Intlist_head_set(newpath->endpoints,high_qend);
    }

#if 0
  } else if (high_univdiagonal + max_insertionlen < univdiagonal) {
    /* Possible with segment search */
    debug13(printf("Impossible\n"));
#endif

  } else if (high_univdiagonal < univdiagonal) {
    /* (A) Insertion */
    nindels = univdiagonal - high_univdiagonal;
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    univdiagonal,/*indels*/+nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/qstart,/*pos3*/high_qend,querylength,
						    indelinfo,plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("(4) Insertion fails\n"));

    } else {
      supporti = indel_pos - qstart;
      supportj = high_qend - (indel_pos + nindels);
      debug13(printf("(7) supporti %d - %d, supportj %d - (%d + %d)\n",
		     indel_pos,qstart,high_qend,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(7) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(7) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	debug13(printf("(3) attach_unknown_qend is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend
					      intlistpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,nmismatches_i);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_j
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(4) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
	
	debug13(Path_print(newpath));
      }
    }

  } else if (high_univdiagonal <= univdiagonal + max_deletionlen) {
    /* (B) Deletion (or short intron) */
    nindels = high_univdiagonal - univdiagonal;
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/univdiagonal,/*indels*/-nindels,chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/qstart,/*pos3*/high_qend,querylength,
						   indelinfo,plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug13(printf("Deletion or short intron fails\n"));
      
    } else {
      supporti = indel_pos - qstart;
      supportj = high_qend - indel_pos;
      debug13(printf("(8) supporti %d - %d, supportj %d - %d\n",
		     indel_pos,qstart,high_qend,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug13(printf("(8) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug13(printf("(8) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
	
	assert(nindels >= 0);
	/* deletionpos = (univdiagonal - querylength) + indel_pos; */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
	  
	debug13(printf("(4) attach_unknown_qend is modifying path %p\n",newpath));
	Intlist_head_set(newpath->endpoints,indel_pos);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,high_qend
					      intlistpool_trace(__FILE__,__LINE__));
	  
	/* For qend, push i first, then push j */
	Intlist_head_set(newpath->nmismatches,nmismatches_i);
	Intlist_head_set(newpath->ref_nmismatches,nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_j
						    intlistpool_trace(__FILE__,__LINE__));
	
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	debug13(printf("(2) Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       qstart,high_qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
	
	debug13(Path_print(newpath));
      }
    }
    
#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* No splicing on circular chromosomes */
    debug13(printf("No splicing on circular chromosomes\n"));
#endif

  } else if (splicingp == false) {
    /* Unable to try splicing */

  } else if (high_univdiagonal > univdiagonal + shortsplicedist) {
    /* Gap is too large for a splice */
    debug13(printf("Gap from %u to %u is too large for a splice\n",univdiagonal,high_univdiagonal));

  } else {
    /* (C) Splice with or without indel */
    left = univdiagonal - (Univcoord_T) querylength;
    /* high_left = high_univdiagonal - (Univcoord_T) querylength; */

    /* Previously filled spliceinfo with knownsplicing information on donors and antiacceptors for segmenti, now done in Splice_resolve_sense */
    /* Previously filled spliceinfo with knownsplicing information on acceptors and antidonors for segmentj, now done in Splice_resolve_sense */

    debug13(printf("BEFORE call to Splice_resolve_qend, %d..%d\n",qstart,high_qend));
    debug13(Path_print(path));

    if (Splice_resolve_qend(&new_univdiagonals,&new_nmismatches,&new_ref_nmismatches,&new_endpoints,&new_junctions,
			    /*univdiagonal_i*/univdiagonal,/*univdiagonal_j*/high_univdiagonal,
			    stage1,query_compress,queryptr,plusp,chroffset,chrhigh,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb,
			    localdb_nmismatches_allowed,
			    /*pos5*/qstart,/*pos3*/high_qend,querylength,cont_trimpos,cont_nmismatches,
			    indelinfo,spliceinfo,knownsplicing,univcoordlistpool,intlistpool,listpool,pathpool,
			    /*sensedir*/try_sensedir,genestrand,/*trim5p*/false,/*trim3p*/true,innerp,
			    /*allow_indels_p*/true,/*find_middle_p*/true,middle_rank) == 0) {
      debug13(printf("Splice_resolve_qend: fails\n"));
      newpath = (T) NULL;

    } else {
      newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

      newpath->univdiagonals = Univcoordlist_append(new_univdiagonals,newpath->univdiagonals);
      newpath->nmismatches = Intlistpool_pop(newpath->nmismatches,intlistpool,&ignore
					     intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlist_append(new_nmismatches,newpath->nmismatches);
      newpath->ref_nmismatches = Intlistpool_pop(newpath->ref_nmismatches,intlistpool,&ignore
						 intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlist_append(new_ref_nmismatches,newpath->ref_nmismatches);
      newpath->endpoints = Intlistpool_pop(newpath->endpoints,intlistpool,&ignore
					   intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlist_append(new_endpoints,newpath->endpoints);
      newpath->junctions = List_append(new_junctions,newpath->junctions);
	
      debug13(printf("AFTER call to Splice_resolve_qend\n"));
      debug13(Path_print(newpath));
    }
  }

  debug13(printf("attach_unknown_qend is returning newpath %p\n",newpath));
  debug13(Path_print(newpath));

  assert(newpath != path);
  return newpath;
}


/* Adds a diagonal to the end of a path, as found by
   compute_qend_local.  Calls attach_unknown_qend to modify the path,
   potentially leading to multiple paths */
static List_T
multiadd_splice_qends (T path, Univcoord_T univdiagonal, Splicetype_T splicetype,
		       bool boundedp, int splice_qpos, int anchor_qpos, int medial_nmismatches,
		       double medial_prob, double medial_support, Univcoord_T *distal_splice_positions,
		       int *distal_support, Univdiag_T *distal_univdiags, double *distal_probs, int npartners,

		       Stage1_T stage1, Knownsplicing_T knownsplicing, Compress_T query_compress, char *queryptr,

		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		       Localdb_T localdb, int localdb_nmismatches_allowed,
		       bool plusp, int querylength, bool innerp,

		       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Hitlistpool_T hitlistpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Vectorpool_T vectorpool, int found_sensedir) {
  List_T newpaths;
  Altsplice_T qend_alts;
  T newpath, nextpath;
  Univcoord_T distal_univdiagonal;
  int best_endpoint, distal_trimpos;
  int n_good_partners, goodi, i, middle_rank;

  debug13(printf("Entered multiadd_splice_qends with %d partners at univdiagonal %u, common_splice_qpos %d.  Current path: ",
		 npartners,univdiagonal,splice_qpos));
  debug13(Path_print(path));
  debug13(printf("\n"));

  if (npartners == 1) {
    /* Needed to find spliceindel, where distal_prob is not known yet */
    debug13(printf("Only one partner: position %u\n",distal_splice_positions[0]));
    n_good_partners = 1;
    goodi = 0;
  } else {
    n_good_partners = 0;
    for (i = 0; i < npartners; i++) {
      if (univdiagonal < distal_univdiags[i]->univdiagonal &&
	  Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			  /*univdiagonal_i*/univdiagonal,
			  /*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			  medial_support,medial_nmismatches,
			  distal_support[i],distal_univdiags[i]->nmismatches,
			  plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
			  innerp) == true) {
	debug13(printf("position %d:%d is good\n",i,distal_splice_positions[i]));
	n_good_partners += 1;
	goodi = i;
      }
    }
    debug13(printf("%d partners are good\n",n_good_partners));
  }

  if (materialize_equiv_splice_alts_p == true && n_good_partners > 1) {
    bool have_equiv = false, equivalentp = true;
    Univcoord_T equiv_univdiagonal = 0;
    int equiv_qstart = 0, equiv_qend = 0, equiv_nmismatches = 0, equiv_support = 0;
    double best_prob = -1.0;

    for (i = 0; i < npartners; i++) {
      if (univdiagonal < distal_univdiags[i]->univdiagonal &&
	  Splice_accept_p(splice_qpos,querylength,medial_prob,distal_probs[i],
			  /*univdiagonal_i*/univdiagonal,
			  /*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			  medial_support,medial_nmismatches,
			  distal_support[i],distal_univdiags[i]->nmismatches,
			  plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
			  innerp) == true) {
	if (have_equiv == false) {
	  have_equiv = true;
	  equiv_univdiagonal = distal_univdiags[i]->univdiagonal;
	  equiv_qstart = distal_univdiags[i]->qstart;
	  equiv_qend = distal_univdiags[i]->qend;
	  equiv_nmismatches = distal_univdiags[i]->nmismatches;
	  equiv_support = distal_support[i];
	} else if (distal_univdiags[i]->univdiagonal != equiv_univdiagonal ||
		   distal_univdiags[i]->qstart != equiv_qstart ||
		   distal_univdiags[i]->qend != equiv_qend ||
		   distal_univdiags[i]->nmismatches != equiv_nmismatches ||
		   distal_support[i] != equiv_support) {
	  equivalentp = false;
	}

	if (distal_probs[i] > best_prob) {
	  best_prob = distal_probs[i];
	  goodi = i;
	}
      }
    }

    if (have_equiv == true && equivalentp == true) {
      debug13(printf("Collapsing %d equivalent qend splice partners at univdiagonal %u, %d..%d\n",
		     n_good_partners,equiv_univdiagonal,equiv_qstart,equiv_qend));
      n_good_partners = 1;
    }
  }


  if (n_good_partners == 0) {
    /* Multiple partners, but none good */
    return (List_T) NULL;

  } else if (n_good_partners == 1 /*|| extend_qend_p(distal_matchlengths,npartners,common_splice_qpos,querylength) == false*/) {
    /* attach_splice_qend may need to filter for GOOD_DISTAL_PROB */
    newpaths = (List_T) NULL;
    assert(distal_splice_positions[goodi] != 0);
#if 0
    distal_univdiagonal = distal_splice_positions[goodi] - splice_qpos + querylength;    
#else
    distal_univdiagonal = distal_univdiags[goodi]->univdiagonal;
    distal_trimpos = distal_univdiags[goodi]->qend;
#endif
    if (distal_univdiagonal == univdiagonal) {
      /* Not a splice */
      debug13(printf("Continuing the alignment to the end\n"));
      newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);
      Intlist_head_set(newpath->endpoints,distal_univdiags[goodi]->qend);
      Intlist_head_set(newpath->nmismatches,-1);
      Intlist_head_set(newpath->ref_nmismatches,-1);
      
      debug13(Path_print(newpath));
      newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
			      hitlistpool_trace(__FILE__,__LINE__));
      
#if 0
    } else if (0) {
      /* Previously called attach_splice_qend, but that did not handle spliceindels */
#endif

    } else {
      for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	if ((newpath = attach_unknown_qend(path,/*high_univdiagonal*/distal_univdiagonal,/*high_qend*/distal_trimpos,
						path->chroffset,path->chrhigh,querylength,
						stage1->indelinfo,stage1->spliceinfo,knownsplicing,
						stage1,query_compress,queryptr,
						novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						localdb_tally,localdb_overflowp,max_nblocks,localdb,
						localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
						intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
						vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,middle_rank)) == NULL) {
	  break;
	} else {
	  assert(newpath != path);
	  newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
				  hitlistpool_trace(__FILE__,__LINE__));
	}
      }
      newpaths = List_reverse(newpaths);
    }

    return newpaths;

  } else if (innerp == false &&
	     (distal_univdiagonal =
	      Altsplice_qend_best(&distal_trimpos,splice_qpos,querylength,
				  /*medial_univdiagonal*/univdiagonal,
				  medial_prob,medial_support,medial_nmismatches,distal_support,
				  distal_univdiags,distal_probs,npartners,
				  plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,
				  innerp)) > 0) {
    debug13(printf("Found a best outer splice\n"));

    newpaths = (List_T) NULL;
    for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
      if ((newpath = attach_unknown_qend(path,/*high_univdiagonal*/distal_univdiagonal,/*high_qend*/distal_trimpos,
					     path->chroffset,path->chrhigh,querylength,
					     stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					     stage1,query_compress,queryptr,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,localdb,
					     localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
					     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,middle_rank)) == NULL) {
	break;
      } else {
	assert(newpath != path);
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
    return List_reverse(newpaths);


  } else if ((qend_alts =
	      Altsplice_qend_new(&best_endpoint,boundedp,splice_qpos,anchor_qpos,querylength,
				 /*medial_univdiagonal*/univdiagonal,medial_nmismatches,
				 medial_prob,medial_support,distal_splice_positions,
				 distal_support,distal_univdiags,distal_probs,
				 /*order_medial_qpos*/NULL,npartners,
				 plusp,/*sense_forward_p*/found_sensedir == SENSE_FORWARD,innerp,
				 pathpool,vectorpool,/*sort_bydistal_p*/false)) == NULL) {
    /* Multiple partners do not work */
    return (List_T) NULL;
    
  } else {
    /* Multiple partners extend to the start, so create an alt substring */
    /* Altsplice_qend_new needs to filter for good partners */
    debug13(printf("Creating an alt substring\n"));
    newpath = Path_copy_3(path,/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool);

    debug13(printf("Combining qend alts for path:\n"));
    /* No need to sort because we are using a common splice qpos */
    newpath->qend_alts = qend_alts;
    Intlist_head_set(newpath->endpoints,best_endpoint);
    Intlist_head_set(newpath->nmismatches,-1);
    Intlist_head_set(newpath->ref_nmismatches,-1);

    debug13(printf("(4) Resulting path: "));
    debug13(Path_print(newpath));
    debug13(printf("add_qend_local is returning path %p as newpaths\n",newpath));

    assert(newpath != path);

    if (innerp == true) {
      /* Allow Altsplice_resolve procedures to select later */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((goodi = Altsplice_select_qend(newpath->qend_alts,querylength,
					      query_compress,plusp,/*genestrand*/0)) < 0) {
      /* Keep outer altsplice */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else if ((nextpath = attach_unknown_qend(path,/*high_univdiagonal*/distal_univdiags[goodi]->univdiagonal,
					       /*high_qend*/distal_univdiags[goodi]->qend,
					       path->chroffset,path->chrhigh,querylength,
					       stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					       stage1,query_compress,queryptr,
					       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					       localdb_tally,localdb_overflowp,max_nblocks,localdb,
					       localdb_nmismatches_allowed,plusp,/*genestrand*/0,innerp,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					       vectorpool,hitlistpool,/*try_sensedir*/found_sensedir,/*middle_rank*/0)) == NULL) {
      /* Keep outer altsplice */
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(newpath)
			  hitlistpool_trace(__FILE__,__LINE__));

    } else {
      /* Use solved outer splice */
      Path_free(&newpath,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      debug13(Path_print(nextpath));
      return Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(nextpath)
			  hitlistpool_trace(__FILE__,__LINE__));
    }
  }
}


/* Sometimes merging of left and right paths can result in anomalies */
static bool
endpoints_acceptable_p (Intlist_T endpoints, List_T junctions) {
  Intlist_T p;
  List_T q;
  Junction_T junction, next_junction;
  int last_endpoint;

  debug2(printf("Evaluating endpoints for acceptability: %s\n",Intlist_to_string(endpoints)));

  /* last_endpoint = 0; */
  /* Skip first endpoint */
  for (p = Intlist_next(endpoints), q = junctions; Intlist_next(p) != NULL; p = Intlist_next(p), q = List_next(q)) {
    last_endpoint = Intlist_head(p);
    junction = (Junction_T) List_head(q);
    if (last_endpoint + Junction_ninserts(junction) > Intlist_second_value(p)) {
      debug2(printf("Endpoint %d + %d > %d, so unacceptable\n",
		    last_endpoint,Junction_ninserts(junction),Intlist_second_value(p)));
      return false;
    } else if (last_endpoint + Junction_ninserts(junction) == Intlist_second_value(p) &&
	       (Parameters_fix_splice_boundary_indels_p == false ||
		q == NULL || List_next(q) == NULL ||
		(next_junction = (Junction_T) List_head(List_next(q))) == NULL ||
		!(((Junction_type(junction) == INS_JUNCTION || Junction_type(junction) == DEL_JUNCTION) &&
		   Junction_type(next_junction) == SPLICE_JUNCTION) ||
		  (Junction_type(junction) == SPLICE_JUNCTION &&
		   (Junction_type(next_junction) == INS_JUNCTION || Junction_type(next_junction) == DEL_JUNCTION))))) {
      debug2(printf("Endpoint %d + %d == %d without an enabled indel/splice boundary, so unacceptable\n",
		    last_endpoint,Junction_ninserts(junction),Intlist_second_value(p)));
      return false;
    } else {
      debug2(printf("Endpoint %d + %d <= %d, so acceptable\n",
		    last_endpoint,Junction_ninserts(junction),Intlist_second_value(p)));
    }
  }

  return true;
}


static bool
indel_splice_boundary_p (Junction_T junction, Junction_T next_junction) {
  Junctiontype_T type, next_type;

  if (Parameters_fix_splice_boundary_indels_p == false) {
    return false;
  } else if (junction == NULL || next_junction == NULL) {
    return false;
  }

  type = Junction_type(junction);
  next_type = Junction_type(next_junction);
  return (((type == INS_JUNCTION || type == DEL_JUNCTION) && next_type == SPLICE_JUNCTION) ||
          (type == SPLICE_JUNCTION && (next_type == INS_JUNCTION || next_type == DEL_JUNCTION))) ? true : false;
}

static bool
endpoints_monotonic_p (Intlist_T endpoints, List_T junctions) {
  int prev_endpoint;
  Intlist_T q;
  List_T j;
  Junction_T junction, next_junction;

  prev_endpoint = Intlist_head(endpoints);
  for (q = Intlist_next(endpoints), j = junctions; q != NULL; q = Intlist_next(q), j = (j == NULL) ? NULL : List_next(j)) {
    if (Intlist_head(q) < prev_endpoint) {
      return false;
    } else if (Intlist_head(q) == prev_endpoint) {
      junction = (j == NULL) ? NULL : (Junction_T) List_head(j);
      next_junction = (j == NULL || List_next(j) == NULL) ? NULL : (Junction_T) List_head(List_next(j));
      if (indel_splice_boundary_p(junction,next_junction) == false) {
        return false;
      }
    }
    prev_endpoint = Intlist_head(q);
  }

  return true;
}


#if 0
/* Ignores values of -1 (unknown) */
static int
preliminary_score_within_trims (Intlist_T nmismatches) {
  int score = 0;
  Intlist_T p;

  for (p = nmismatches; p != NULL; p = Intlist_next(p)) {
    if (Intlist_head(p) >= 0) {
      score += Intlist_head(p);
    }
  }

  return score;
}
#endif


/* Modified from Path_exon_origin */
static int
compute_exon_origin (Intlist_T endpoints, List_T junctions) {
  int exon_origin;
  Intlist_T p = endpoints;
  List_T j = junctions;

  p = Intlist_next(p);
  exon_origin = Intlist_head(p);

  while (j != NULL && Junction_type((Junction_T) List_head(j)) != SPLICE_JUNCTION) {
    p = Intlist_next(p);
    exon_origin = Intlist_head(p);

    j = List_next(j);
  }

  return exon_origin;
}


/* Always solves against plus strand of genome.  Just provide either
   queryuc/query_compress_fwd (coords measured from beginning of
   sequence) or queryrc/query_compress_rev (coords measured from end
   of sequence).  All coordinates measured from low end.
   Sense/antisense is with respect to the plus strand.  But to
   interface with Stage3end_new_substrings command, need to flip
   coordinates for case where queryrc aligns to plus strand. */

/* chrnum is fixed from middle_diagonal */
static List_T
combine_leftright_paths (int *found_score, int *end_score, Univcoord_T main_univdiagonal,
			 List_T qstart_paths, List_T qend_paths,
			 Compress_T query_compress, Compress_T query_compress_fwd,
			 Compress_T query_compress_rev, char *queryptr, int querylength,
			 bool plusp, bool first_read_p, int genestrand, int sensedir,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 int *mismatch_positions_alloc, Stage1_T stage1,
			 Knownsplicing_T knownsplicing, Spliceendsgen_T spliceendsgen,

			 Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			 unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,

			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Univdiagpool_T univdiagpool,
			 Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, Method_T method,
			 bool find_splices_p) {
  List_T paths = NULL, best_qstart_paths, best_qend_paths;
  List_T a, b;
  T path, qstart_path, qend_path;

  /* int best_nmatches, best_ref_nmatches, nmatches; */
  /* double best_splice_prob; */

  Intlist_T endpoints, q;
  Univcoordlist_T univdiagonals, u;
  Univcoord_T univdiagonal_trimmed;
  Intlist_T nmismatches, ref_nmismatches, s, r;
  List_T junctions, j;
  Junction_T junction;

  int qstart1, qend1, qstart2, qend2, ninserts1, ninserts2;
  int middle_nmismatches, middle_ref_nmismatches;

  int found_sensedir;
  int trimpos, ninserts, exon_origin, qstart, qend, pos5, pos3;
  int nmismatches_to_trimpos;
  bool cassettep, splice5p, splice3p;
  Splicetype_T splicetype5, splicetype3;
  double ambig_prob_5, ambig_prob_3;



#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  printf("\n");
  printf("*** %s: Entered combine_leftright_paths with %d qstart paths and %d qend paths\n",
	 sensedir == SENSE_FORWARD ? "sense" : "antisense",List_length(qstart_paths),List_length(qend_paths));
#endif

  /* assert(qstart_paths != NULL && qend_paths != NULL); */

#ifdef PRUNE_PATHS
  /* Find best qstart paths */
  if (List_length(qstart_paths) == 1) {
    best_qstart_paths = Hitlist_copy(qstart_paths,hitlistpool);

  } else {
    best_qstart_paths = (List_T) NULL;
    best_nmatches = 0;
    best_ref_nmatches = 0;
    best_splice_prob = 0.0;
    for (a = qstart_paths; a != NULL; a = List_next(a)) {
      qstart_path = (T) List_head(a);
      nmatches = Path_eval_nmatches(&ignore_found_score,qstart_path,
				    query_compress_fwd,query_compress_rev);
      if (nmatches < best_nmatches) {
	/* Worse than current best */
	debug13(printf("=> (X) Worse than current best by nmatches\n"));

      } else if (nmatches > best_nmatches) {
	/* Better than current best */
	Hitlistpool_free_list(&best_qstart_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qstart_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(qstart_path)
					 hitlistpool_trace(__FILE__,__LINE__));
	best_nmatches = nmatches;
	best_ref_nmatches = qstart_path->ref_nmatches;
	best_splice_prob = qstart_path->junction_splice_prob;

      } else if (qstart_path->ref_nmatches < best_ref_nmatches) {
	/* Worse than current best */
	debug13(printf("=> (Y) Worse than current best by ref_nmatches\n"));

      } else if (qstart_path->ref_nmatches > best_ref_nmatches) {
	/* Better than current best */
	Hitlistpool_free_list(&best_qstart_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qstart_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(qstart_path)
					 hitlistpool_trace(__FILE__,__LINE__));
	/* best_nmatches = nmatches; */
	best_ref_nmatches = qstart_path->ref_nmatches;
	best_splice_prob = qstart_path->junction_splice_prob;

      } else if (qstart_path->junction_splice_prob < best_splice_prob) {
	/* Worse than current best */
	debug13(printf("=> (Z) Worse than current best by splice_prob\n"));

      } else if (qstart_path->junction_splice_prob > best_splice_prob) {
	Hitlistpool_free_list(&best_qstart_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qstart_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_fwd(qstart_path)
					 hitlistpool_trace(__FILE__,__LINE__));
	/* best_nmatches = nmatches; */
	/* best_ref_nmatches = qstart_path->ref_nmatches; */
	best_splice_prob = qstart_path->junction_splice_prob;

      } else {
	/* Same as current best */
	best_qstart_paths = Hitlist_push(best_qstart_paths,hitlistpool,(void *) Path_expect_fwd(qstart_path)
					 hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }
#else
  best_qstart_paths = qstart_paths;
#endif

    
#ifdef PRUNE_PATHS
  /* Find best qend paths */
  if (List_length(qend_paths) == 1) {
    best_qend_paths = Hitlist_copy(qend_paths,hitlistpool);

  } else {
    best_qend_paths = (List_T) NULL;
    best_nmatches = 0;
    best_ref_nmatches = 0;
    best_splice_prob = 0.0;

    /* Don't use splice prob, since that can introduce bad splices */
    for (b = qend_paths; b != NULL; b = List_next(b)) {
      qend_path = (T) List_head(b);

      Path_reverse(qend_path,/*expect_fwd_p*/true);
      nmatches = Path_eval_nmatches(&ignore_found_score,qend_path,
				    query_compress_fwd,query_compress_rev);
      Path_reverse(qend_path,/*expect_fwd_p*/false);

      if (nmatches < best_nmatches) {
	/* Worse than current best */
	debug13(printf("=> (X) Worse than current best by nmatches\n"));

      } else if (nmatches > best_nmatches) {
	/* Better than current best */
	debug13(printf("=> (X) Better than current best by nmatches\n"));
	Hitlistpool_free_list(&best_qend_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qend_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(qend_path)
				       hitlistpool_trace(__FILE__,__LINE__));
	best_nmatches = nmatches;
	best_ref_nmatches = qend_path->ref_nmatches;
	best_splice_prob = qend_path->junction_splice_prob;

      } else if (qend_path->ref_nmatches < best_ref_nmatches) {
	/* Worse than current best */
	debug13(printf("=> (Y) Worse than current best by ref_nmatches\n"));

      } else if (qend_path->ref_nmatches > best_ref_nmatches) {
	/* Better than current best */
	debug13(printf("=> (X) Better than current best by nmatches\n"));
	Hitlistpool_free_list(&best_qend_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qend_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(qend_path)
				       hitlistpool_trace(__FILE__,__LINE__));
	/* best_nmatches = nmatches; */
	best_ref_nmatches = qend_path->ref_nmatches;
	best_splice_prob = qend_path->junction_splice_prob;

      } else if (qend_path->junction_splice_prob < best_splice_prob) {
	/* Worse than current best */
	debug13(printf("=> (Z) Worse than current best by splice_prob\n"));

      } else if (qend_path->junction_splice_prob > best_splice_prob) {
	Hitlistpool_free_list(&best_qend_paths,hitlistpool
			      hitlistpool_trace(__FILE__,__LINE__));
	best_qend_paths = Hitlist_push(NULL,hitlistpool,(void *) Path_expect_rev(qend_path)
				       hitlistpool_trace(__FILE__,__LINE__));
	/* best_nmatches = nmatches; */
	/* best_ref_nmatches = qend_path->ref_nmatches; */
	best_splice_prob = qend_path->junction_splice_prob;

      } else {
	/* Same as current best */
	best_qend_paths = Hitlist_push(best_qend_paths,hitlistpool,(void *) Path_expect_rev(qend_path)
				       hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }
#else
  best_qend_paths = qend_paths;
#endif

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  printf("Now have %d best qstart paths and %d best qend paths\n",
	 List_length(best_qstart_paths),List_length(best_qend_paths));
  for (a = best_qstart_paths; a != NULL; a = List_next(a)) {
    qstart_path = (T) List_head(a);
    printf("Qstart: ");
    Path_print(qstart_path);
  }
  for (b = best_qend_paths; b != NULL; b = List_next(b)) {
    qend_path = (T) List_head(b);
    printf("Qend: ");
    Path_print(qend_path);
  }
#endif

  for (a = best_qstart_paths; a != NULL; a = List_next(a)) {
    qstart_path = (T) List_head(a);
    assert(qstart_path->sensedir == sensedir);

    for (b = best_qend_paths; b != NULL; b = List_next(b)) {
      qend_path = (T) List_head(b);
      assert(qend_path->sensedir == sensedir);
	
      debug13(printf("++ Qstart/left path %p: ",qstart_path));
      debug13(Path_print(qstart_path));
      debug13(printf("++ Qend/right path %p: ",qend_path));
      debug13(Path_print(qend_path));
      debug13(printf("\n"));

      /* Combine qstart_path with qend_path */
      /* If either list is NULL, must have obtained an unacceptable result */
	
      /* Example:
	 qstart_path->endpoints:        59(qstart1), 77(qend1).
	 qend_path->endpoints: 100, 79, 68(qend2),   57(qstart2).
	   
	 endpoints1:           59                77
	 endpoints2:        57               68, 79, 100
	   
	 Desired result: 57, 68, 79, 100.
	   
	 Could use middle_univdiagonal to find the original, common
	 segment, but it should be the last segment in each case,
	 since we pushed results on top for qstart and for qend */

      /* endpoints = (Intlist_T) NULL; -- Initialized with first push */
      univdiagonals = (Univcoordlist_T) NULL;
      nmismatches = (Intlist_T) NULL;
      ref_nmismatches = (Intlist_T) NULL;
      junctions = (List_T) NULL;
	
      qend1 = Intlist_last_value(qstart_path->endpoints);
      qstart1 = Intlist_penultimate_value(qstart_path->endpoints);
	
      qstart2 = Intlist_last_value(qend_path->endpoints);
      qend2 = Intlist_penultimate_value(qend_path->endpoints);
	
      q = qstart_path->endpoints;
      u = qstart_path->univdiagonals;
      s = qstart_path->nmismatches;
      r = qstart_path->ref_nmismatches;
      j = qstart_path->junctions;
	
      junction = (Junction_T) NULL;
      endpoints = Intlistpool_push(NULL,intlistpool,Intlist_head(q)
				   intlistpool_trace(__FILE__,__LINE__));
      q = Intlist_next(q);
      ninserts1 = 0;
      while (j != NULL) {
	endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(q)
				     intlistpool_trace(__FILE__,__LINE__));
	univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(u)
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(s)
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(r)
					   intlistpool_trace(__FILE__,__LINE__));
	junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(j),pathpool)
				  listpool_trace(__FILE__,__LINE__));
	junction = (Junction_T) List_head(junctions);
	ninserts1 = Junction_ninserts(junction);
	  
	q = Intlist_next(q);
	u = Univcoordlist_next(u);
	s = Intlist_next(s);
	r = Intlist_next(r);
	j = List_next(j);
      }
	
      /* Reached middle_univdiagonal of qstart1 */
      qstart1 = Intlist_head(endpoints);
      qend1 = Intlist_head(q);
	
      Path_reverse(qend_path,/*expect_fwd_p*/true);
	
      qstart2 = Intlist_head(qend_path->endpoints);
      qend2 = Intlist_second_value(qend_path->endpoints);
      ninserts2 = 0;
      if (qend2 < qstart1) {
	/* No overlap: Apparent cassette difference */
	debug13(printf("++ Combined path has a middle cassette difference, due to lack of overlap between %d..%d and %d..%d, where %d < %d\n",
		       qstart1,qend1,qstart2,qend2,qend2,qstart1));
#if 0
	Intlistpool_free_list(&endpoints,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				    univcoordlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&ref_nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Junction_list_gc(&junctions,listpool,pathpool);
#else
	/* The endpoints now specify the cassette difference, but should re-extended later by calls to Spliceends_qstart_trim and Spliceends_qend_trim */
	cassettep = true;
	endpoints = Intlistpool_push(endpoints,intlistpool,qend2
				     intlistpool_trace(__FILE__,__LINE__));
	univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(u)
					       univcoordlistpool_trace(__FILE__,__LINE__));

	nmismatches = Intlistpool_push(nmismatches,intlistpool,-1
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,-1
					   intlistpool_trace(__FILE__,__LINE__));
#endif

      } else {
	cassettep = false;
	endpoints = Intlistpool_push(endpoints,intlistpool,qend2
				     intlistpool_trace(__FILE__,__LINE__));
	univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(u)
					       univcoordlistpool_trace(__FILE__,__LINE__));

	if (qstart1 + ninserts1 == qstart2 + ninserts2 && qend1 == qend2) {
	  /* Take nmismatches from either qstart_path or qend_path (if available) */
	  if ((middle_nmismatches = Intlist_head(s)) == -1) {
	    middle_nmismatches = Intlist_head(qend_path->nmismatches);
	  }
	  if ((middle_ref_nmismatches = Intlist_head(r)) == -1) {
	    middle_ref_nmismatches = Intlist_head(qend_path->ref_nmismatches);
	  }
	  
	} else if (qstart1 + ninserts1 == qstart2 + ninserts2) {
	  /* Take nmismatches from either qend_path */
	  middle_nmismatches = Intlist_head(qend_path->nmismatches);
	  middle_ref_nmismatches = Intlist_head(qend_path->ref_nmismatches);
	  
	} else if (qend1 == qend2) {
	  /* Take nmismatches from either qstart_path */
	  middle_nmismatches = Intlist_head(s);
	  middle_ref_nmismatches = Intlist_head(r);
	  
	} else {
	  middle_nmismatches = -1;
	  middle_ref_nmismatches = -1;
	}
	nmismatches = Intlistpool_push(nmismatches,intlistpool,middle_nmismatches
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,middle_ref_nmismatches
					   intlistpool_trace(__FILE__,__LINE__));
      }

      q = Intlist_next(qend_path->endpoints);
      u = qend_path->univdiagonals;
      s = qend_path->nmismatches;
      r = qend_path->ref_nmismatches;
      j = qend_path->junctions;
	
      while (j != NULL) {
	junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(j),pathpool)
				  listpool_trace(__FILE__,__LINE__));
	junction = (Junction_T) List_head(junctions);
	ninserts2 = Junction_ninserts(junction);
	
	q = Intlist_next(q);
	u = Univcoordlist_next(u);
	s = Intlist_next(s);
	r = Intlist_next(r);
	endpoints = Intlistpool_push(endpoints,intlistpool,Intlist_head(q)
				     intlistpool_trace(__FILE__,__LINE__));
	univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,Univcoordlist_head(u)
					       univcoordlistpool_trace(__FILE__,__LINE__));
	nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(s)
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(r)
					   intlistpool_trace(__FILE__,__LINE__));
	j = List_next(j);
      }

      endpoints = Intlist_reverse(endpoints);
      junctions = List_reverse(junctions);
      if (endpoints_acceptable_p(endpoints,junctions) == false) {
	debug13(printf("++ Combined path not possible, due to unacceptable endpoints\n"));
	Intlistpool_free_list(&endpoints,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				    univcoordlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&ref_nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Junction_list_gc(&junctions,listpool,pathpool);
	
      } else {
	univdiagonals = Univcoordlist_reverse(univdiagonals);
	nmismatches = Intlist_reverse(nmismatches);
	ref_nmismatches = Intlist_reverse(ref_nmismatches);

	pos5 = Intlist_penultimate_value(endpoints); /* For qend trimming */
	pos3 = Intlist_second_value(endpoints);      /* For qstart trimming */
	univdiagonal_trimmed = 0;

	/* qend_path->splice3p == true || qend_path->qend_alts != NULL || (qend = Intlist_last_value(endpoints)) == querylength */
	qend = Intlist_last_value(endpoints);
	if (cassettep == false) {
	  /* Should have found trimming from Path_solve procedures */
	  splice3p = qend_path->splice3p;
	  splicetype3 = qend_path->splicetype3;
	  ambig_prob_3 = qend_path->ambig_prob_3;

	} else {
	  endpoints = Intlist_reverse(endpoints);
	  univdiagonals = Univcoordlist_reverse(univdiagonals);
	  junctions = List_reverse(junctions);
	  if (junctions == NULL) {
	    ninserts = 0;
	  } else {
	    ninserts = Junction_ninserts((Junction_T) List_head(junctions));
	  }
	    
	  debug13(printf("endpoints: %s\n",Intlist_to_string(endpoints)));
	  debug13(printf("nmismatches: %s\n",Intlist_to_string(nmismatches)));

	  exon_origin = compute_exon_origin(endpoints,junctions);
	  debug13(printf("Calling Spliceends_qend_trim with pos5 %d\n",pos5 + ninserts));
	  splice3p = Spliceends_qend_trim(&trimpos,&nmismatches_to_trimpos,
					  &found_sensedir,&splicetype3,&ambig_prob_3,
					  stage1,knownsplicing,sensedir,
					  Univcoordlist_head(univdiagonals),querylength,
					  pos5 + ninserts,
					  exon_origin,chrnum,chroffset,chrhigh,
					  plusp,genestrand,mismatch_positions_alloc,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,
					  univdiagpool,vectorpool,spliceendsgen,
					  query_compress,queryptr,
					  genomebits,genomebits_alt,find_splices_p);
	  debug13(printf("(3) Spliceends_qend_trim returns trimpos %d and %d nmismatches, splice3 prob %f\n",
			 trimpos,nmismatches_to_trimpos,ambig_prob_3));

	  if (trimpos == pos5) {
	    debug13(printf("Segment is not great, since nosplice_trimpos erases segment, so not trimming\n"));
	      
	  } else if (trimpos == qend) {
	    /* No change, so can keep nmismatches */
	    debug13(printf("No change, so keeping nmismatches as %d\n",Intlist_last_value(nmismatches)));
	      
	  } else if (Univcoordlist_head(univdiagonals) - querylength + trimpos < chroffset) {
	    debug13(printf("Attempt to trim beyond start of chromosome (%u < %u), so not trimming\n",
			   Univcoordlist_head(univdiagonals) - querylength + trimpos,chroffset));
	      
	  } else if (trimpos <= pos5 + ninserts) {
	    debug13(printf("Attempt to trim before an insertion (%d <= %d+%d), so not trimming\n",
			   trimpos,pos5,ninserts));
	  } else {
	    /* Change both endpoints later */
	    /* Use new nmismatches_to_trimpos */
	    debug13(printf("combine_leftright_paths changing 3' endpoint from %d to %d\n",Intlist_head(endpoints),trimpos));
	    Altsplice_free(&qend_path->qend_alts,pathpool);
	    qend_path->qend_alts = (Altsplice_T) NULL;

	    Intlist_head_set(endpoints,trimpos);
	    univdiagonal_trimmed = Univcoordlist_head(univdiagonals);
	      
	    nmismatches = Intlist_reverse(nmismatches);
	    ref_nmismatches = Intlist_reverse(ref_nmismatches);
	    if (cassettep == true) {
	      Intlist_head_set(nmismatches,-1);
	      Intlist_head_set(ref_nmismatches,-1);
	    } else {
	      Intlist_head_set(nmismatches,nmismatches_to_trimpos);
	      Intlist_head_set(ref_nmismatches,nmismatches_to_trimpos);
	    }
	    nmismatches = Intlist_reverse(nmismatches);
	    ref_nmismatches = Intlist_reverse(ref_nmismatches);
	  }
	    
	  endpoints = Intlist_reverse(endpoints);
	  univdiagonals = Univcoordlist_reverse(univdiagonals);
	  junctions = List_reverse(junctions);
	}
	  
	/*qstart_path->splice5p == true || qstart_path->qstart_alts != NULL || (qstart = Intlist_head(endpoints)) == 0 */
	qstart = Intlist_head(endpoints);
	if (cassettep == false) {
	  /* Should have found endpoints from Path_solve procedures */
	  splice5p = qstart_path->splice5p;
	  splicetype5 = qstart_path->splicetype5;
	  ambig_prob_5 = qstart_path->ambig_prob_5;

	} else {
	  debug13(printf("endpoints: %s\n",Intlist_to_string(endpoints)));
	  debug13(printf("nmismatches: %s\n",Intlist_to_string(nmismatches)));
	    
	  exon_origin = compute_exon_origin(endpoints,junctions);
	  debug13(printf("Calling Spliceends_qstart_trim with pos3 %d\n",pos3));
	  splice5p = Spliceends_qstart_trim(&trimpos,&nmismatches_to_trimpos,&found_sensedir,&splicetype5,&ambig_prob_5,
					    stage1,knownsplicing,sensedir,
					    Univcoordlist_head(univdiagonals),querylength,pos3,
					    exon_origin,chrnum,chroffset,chrhigh,
					    plusp,genestrand,mismatch_positions_alloc,
					    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					    localdb_tally,localdb_overflowp,max_nblocks,
					    univdiagpool,vectorpool,spliceendsgen,
					    query_compress,queryptr,
					    genomebits,genomebits_alt,find_splices_p);
	  debug13(printf("(3) Spliceends_qstart_trim returns trimpos %d and %d nmismatches, splice5 prob %f\n",
			 trimpos,nmismatches_to_trimpos,ambig_prob_5));

	  if (trimpos == pos3) {
	    debug13(printf("Segment is not great, since nosplice_trimpos erases segment, so not trimming\n"));
	      
	  } else if (trimpos == qstart) {
	    /* No change, so can keep nmismatches */
	    debug13(printf("No change, so keeping nmismatches as %d\n",Intlist_head(nmismatches)));
	      
	  } else if (Univcoordlist_head(univdiagonals) - querylength + trimpos >= chrhigh) {
	    debug13(printf("Attempt to trim beyond end of chromosome (%u > %u), so not trimming\n",
			   Univcoordlist_head(univdiagonals) - querylength + trimpos,chrhigh));
	  } else {
	    /* Use new nmismatches_to_trimpos */
	    debug13(printf("Changing 5' endpoint from %d to %d\n",Intlist_head(endpoints),trimpos));
	    Altsplice_free(&qstart_path->qstart_alts,pathpool);
	    qstart_path->qstart_alts = (Altsplice_T) NULL;

	    Intlist_head_set(endpoints,trimpos);
	    if (cassettep == true) {
	      Intlist_head_set(nmismatches,-1);
	      Intlist_head_set(ref_nmismatches,-1);
	    } else if (Univcoordlist_head(univdiagonals) == univdiagonal_trimmed) {
	      /* Trimmed both ends with different parameters, so nmismatches_to_trimpos is not valid */
	      Intlist_head_set(nmismatches,-1);
	      Intlist_head_set(ref_nmismatches,-1);
	    } else {
	      Intlist_head_set(nmismatches,nmismatches_to_trimpos);
	      Intlist_head_set(ref_nmismatches,nmismatches_to_trimpos);
	    }
	  }
	}
	  
	debug13(printf("Before Path_create, endpoints are %s\n",Intlist_to_string(endpoints)));
	if (endpoints_monotonic_p(endpoints,junctions) == false) {
	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Junction_list_gc(&junctions,listpool,pathpool);

	} else {
	  path = Path_create(main_univdiagonal,
			     endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,
			     plusp,first_read_p,genestrand,sensedir,querylength,
			     method,chrnum,chroffset,chrhigh,
			     splice5p,splicetype5,ambig_prob_5,
			     splice3p,splicetype3,ambig_prob_3,
			     qstart_path->qstart_alts,qend_path->qend_alts,
			     query_compress_fwd,query_compress_rev,pathpool,vectorpool);
	  Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
	  debug13(printf("++ Combined path %p: ",path));
	  debug13(Path_print(path));
	  debug13(printf("\n"));

	  debug13(printf("combine_leftright_paths yields result with score %d\n",path->score_within_trims));
	  paths = Hitlist_push(paths,hitlistpool,(void *) Path_expect_fwd(path)
			       hitlistpool_trace(__FILE__,__LINE__));
	}
      }

      /* Undo reversal */
      Path_reverse(qend_path,/*expect_fwd_p*/false);
      debug13(printf("Done with qend path\n"));
    }
    debug13(printf("Done with qstart path\n"));
  }
 
#ifdef PRUNE_PATHS
  Hitlistpool_free_list(&best_qend_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&best_qstart_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
#endif

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  printf("\n");
  printf("*** Exiting combine_leftright_paths.  Paths now has length %d\n",List_length(paths));
#endif

  return paths;
}


#if 0
static void
check_for_descending_qend (int prev_qend, Univcoord_T prev_univdiagonal, List_T qstart_diagonals) {
  List_T p;
  Univdiag_T diagonal;

  for (p = qstart_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(qstart_diagonals);
    if (diagonal->qend > prev_qend) {
      abort();
    } else if (diagonal->qend == prev_qend && diagonal->univdiagonal > prev_univdiagonal) {
      abort();
    }
    prev_qend = diagonal->qend;
    prev_univdiagonal = diagonal->univdiagonal;
  }
  return;
}
#endif

#if 0
static void
check_for_ascending_qstart (int prev_qstart, Univcoord_T prev_univdiagonal, List_T qend_diagonals) {
  List_T p;
  Univdiag_T diagonal;

  for (p = qend_diagonals; p != NULL; p = List_next(p)) {
    diagonal = (Univdiag_T) List_head(qend_diagonals);
    if (diagonal->qstart < prev_qstart) {
      abort();
    } else if (diagonal->qstart == prev_qstart && diagonal->univdiagonal < prev_univdiagonal) {
      abort();
    }
    prev_qstart = diagonal->qstart;
    prev_univdiagonal = diagonal->univdiagonal;
  }
  return;
}
#endif
  

/* Recursively adds diagonals to the start of a path, using either knownsplicing or localdb */
/* Note: is is important to trim ends before calling Localdb_get,
   because the presence of a splice site indicates the boundary point
   for searching the distal end */

/* Either returns NULL, which means caller should use the given path (or a copy)
   as the result, or returns newpaths, which of which is different
   from path, where the caller needs to free path, if desired */

/* unextended_paths are all fwd.  complete_paths are all fwd */
static void
compute_qstart_local (List_T *qstart_paths,
		      int depth, Path_T path, char *queryptr, int querylength,
		      int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		      Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		      Compress_T query_compress, Chrnum_T chrnum, Univcoord_T chroffset,
		      Univcoord_T chrhigh, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
		      Hitlistpool_T hitlistpool, Spliceendsgen_T spliceendsgen,
		      int localdb_nmismatches_allowed, bool plusp, int genestrand,
		      int sensedir, bool innerp, bool find_splices_p, bool salvagep);

static void
compute_qstart_local_work (List_T *qstart_paths,
		      int depth, Path_T path, char *queryptr, int querylength,
		      int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		      Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels, Compress_T query_compress, 
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		      Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		      Spliceendsgen_T spliceendsgen, int localdb_nmismatches_allowed,
		      bool plusp, int genestrand, int sensedir, bool innerp, bool find_splices_p, bool salvagep) {
  List_T newpaths = NULL, p;
  T splice_path, end_path, newpath, parent;
  Univcoord_T indel_univdiagonal, univdiagonal;
  int qstart, pos3, exon_origin, indel_pos, trimpos, adj;
  int orig_qstart, child_qstart;
  int nmismatches_i;

  bool splice5p;
  Splicetype_T splicetype5;
  int nspliceends, i, j;
  Spliceends_T spliceends;
  int nosplice_trimpos, medial_splice_trimpos_0;
  double ambig_prob_5;
  int nosplice_nmismatches, medial_splice_nmismatches_0;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long spliceends_key_hash, spliceends_result_hash;
#endif


  debug13(printf("Entering compute_qstart_local, salvage %d, at depth %d with path:\n",salvagep,depth));
  debug13(Path_print(path));
  debug13(printf("\n"));
  Path_expect_fwd(path);

  orig_qstart = Path_qstart(path);

  if (path->qstart_alts != NULL) {
    /* Not possible to add localdb because of the ambiguity */
    debug13(printf("compute_qstart_local is complete because of qstart_alts\n"));
    path->completep = Path_completep(path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(1) Putting into complete qstart paths: ")); debug20(Path_print(newpath));
    *qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else if ((qstart = Intlist_head(path->endpoints)) == 0) {
    /* Already at start */
    debug13(printf("compute_qstart_local is complete because of qstart == 0\n"));
    path->completep = Path_completep(path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(2) Putting into complete qstart paths: ")); debug20(Path_print(newpath));
    *qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else if (depth > max_recursive_ends) {
    /* Too much recursive depth */
    debug13(printf("compute_qstart_local is unextended because depth %d > max %d\n",
		   depth,max_recursive_ends));
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(3) Putting into unextended qstart paths: ")); debug20(Path_print(newpath));
    assert(newpath->sensedir == sensedir);
    *qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else if (splicingp == false) {
    /* Not able to try splicing */
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(4) Putting into unextended qstart paths: ")); debug20(Path_print(newpath));
    assert(newpath->sensedir == sensedir);
    *qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
				 hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else {
    /* Attempt to extend */
    debug13(printf("compute_qstart_local extension state: depth %d, orig_qstart %d\n",
		   depth,orig_qstart));

    /* pos3 is the end of the current segment.  qstart is the
       approximate start of the current segment.  pos5 (0) to qstart
       is where we look for a new segment. */
    pos3 = Intlist_second_value(path->endpoints);
    univdiagonal = Univcoordlist_head(path->univdiagonals);
    /* left = univdiagonal - (Univcoord_T) querylength; */
    
    exon_origin = Path_exon_origin(path);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    spliceends_key_hash = spliceends_call_fingerprint(
      "qstart",sensedir,univdiagonal,querylength,pos3,exon_origin,
      chrnum,chroffset,chrhigh,plusp,genestrand,localdb_nmismatches_allowed,
      innerp,salvagep,max_nblocks,stage1,knownsplicing,spliceendsgen,
      query_compress,queryptr,genomebits,genomebits_alt,find_splices_p);
#endif
    debug13(printf("(SPLICE) Calling Spliceends_trimmed_qstarts with univdiagonal %u, pos3 %d, exon_origin %d\n",
		   univdiagonal,pos3,exon_origin));
    nspliceends = Spliceends_trimmed_qstarts(&spliceends,&nosplice_trimpos,&medial_splice_trimpos_0,
					     &nosplice_nmismatches,&medial_splice_nmismatches_0,
					     &splice5p,&splicetype5,&ambig_prob_5,
					     sensedir,univdiagonal,querylength,
					     /*qstart:0,*//*qend*/pos3,exon_origin,
					     chrnum,chroffset,chrhigh,plusp,genestrand,
					     localdb_nmismatches_allowed,innerp,salvagep,
					     mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
					     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					     stage1,knownsplicing,univdiagpool,vectorpool,spliceendsgen,
					     query_compress,queryptr,genomebits,genomebits_alt,find_splices_p);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    spliceends_result_hash = spliceends_result_fingerprint(
      spliceends,nspliceends,nosplice_trimpos,medial_splice_trimpos_0,
      nosplice_nmismatches,medial_splice_nmismatches_0,splice5p,splicetype5,
      ambig_prob_5);
    Accel_profile_spliceends_result("qstart",spliceends_key_hash,nspliceends,
				    spliceends_result_hash);
#endif
    /* path->spliceends_qstart = pos3; */
    debug13(printf("Spliceends_trimmed_qstarts returning %d spliceends, nosplice %d (%d nmismatches), medial splice %d (%d nmismatches)\n",
		   nspliceends,nosplice_trimpos,nosplice_nmismatches,medial_splice_trimpos_0,medial_splice_nmismatches_0));
    
    if (nspliceends == 0) {
      debug13(printf("=> No trimming\n"));
      
    } else if (spliceends == NULL) {
      debug13(printf("=> No partners found\n"));

    } else if (nosplice_trimpos == 0) {
      debug13(printf("=> Actually, not solving for splice because nosplice_trimpos is %d\n",nosplice_trimpos));
      
    } else {
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
      printf("=> Got %d spliceends\n",nspliceends);
      for (int i = 0; i < nspliceends; i++) {
	printf("splice_qpos %d => partner %u, %d..%d, splicedist %u, mismatches %d, prob %f, trimpos %d\n",
	       spliceends->splice_qpos[i],spliceends->partners[i],
	       spliceends->distal_univdiags[i]->qstart,spliceends->distal_univdiags[i]->qend,
	       univdiagonal - (spliceends->partners[i] - spliceends->splice_qpos[i] + querylength),
	       spliceends->distal_univdiags[i]->nmismatches,spliceends->distal_probs[i],
	       spliceends->distal_univdiags[i]->qstart);
      }
#endif
      
      i = 0;
      while (i < nspliceends) {
	/* Assumes that spliceends are sorted primarily by splice_qpos */
	/* ambig_prob_5 = spliceends->medial_probs[i]; */
	qstart = spliceends->splice_qpos[i]; /* was low_qstart */
	j = i + 1;
	while (j < nspliceends && spliceends->splice_qpos[j] == qstart) {
	  j++;
	}
	
	/* Handle a cluster of solutions for qstart */
	debug13(printf("Splice qpos: %d, medial prob %f\n",qstart,ambig_prob_5));
	newpaths = List_append(newpaths,
			       multiadd_splice_qstarts(path,univdiagonal,spliceends->splicetype,
						       spliceends->boundedp,/*splice_qpos*/qstart,/*anchor_qpos*/pos3,
						       spliceends->medial_nmismatches[i],spliceends->medial_probs[i],
						       spliceends->medial_support[i],
						       &(spliceends->partners[i]),&(spliceends->distal_support[i]),
						       &(spliceends->distal_univdiags[i]),&(spliceends->distal_probs[i]),
						       /*npartners*/(j - i),
						       stage1,knownsplicing,query_compress,queryptr,
						       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						       localdb_tally,localdb_overflowp,max_nblocks,
						       localdb,localdb_nmismatches_allowed,plusp,querylength,innerp,
						       intlistpool,univcoordlistpool,listpool,
						       hitlistpool,pathpool,transcriptpool,vectorpool,sensedir));
	i = j;
      }
    }

    /* Spliceends_free(&spliceends); */
    Spliceendsgen_return(spliceendsgen,&spliceends);

    if (newpaths != NULL) {
      debug13(printf("compute_qstart_local has %d splice child paths; still trying start indel\n",
		     List_length(newpaths)));
    }

    if (1 || newpaths == NULL) {
      /* Always try an indel, since splice may be incorrect */
      /* Splice end without partner, or an end indel */
      /* Try start indel */
      debug13(printf("Path before trying start indel:\n"));
      debug13(Path_print(path));
      debug13(printf("\n"));
      
      /* Computed above */
      /* pos3 = Intlist_head(Intlist_next(path->endpoints)); */
      /* univdiagonal = Univcoordlist_head(path->univdiagonals); */
      /* left = univdiagonal - (Univcoord_T) querylength; */
      
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
      if (nosplice_trimpos != 0) {
	printf("(INDEL) Calling Genomebits_indel_solve_low and Spliceends_indel_qstart with %d..%d\n",0,nosplice_trimpos);
      }
#endif
      
      debug13(printf("nosplice_trimpos at qstart is %d\n",nosplice_trimpos));
      if (nosplice_trimpos == 0) {
	/* Already extended to the start */
	debug13(printf("Actually, not solving for indel because nosplice_trimpos is %d\n",nosplice_trimpos));
	
      } else if (knownindels != NULL &&
		 (adj = Knownindels_find_highest(&indel_pos,knownindels,univdiagonal,querylength,
						 /*pos5*/0,/*pos3*/nosplice_trimpos+1)) != 0 &&
		 (newpath = attach_indel_qstart_simple(adj,path,indel_pos,univdiagonal,querylength,sensedir,
						       plusp,genestrand,mismatch_positions_alloc,
						       stage1,knownsplicing,spliceendsgen,query_compress,queryptr,
						       genomebits,genomebits_alt,chrnum,chroffset,chrhigh,find_splices_p,
						       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						       localdb_tally,localdb_overflowp,max_nblocks,
						       intlistpool,univcoordlistpool,listpool,
						       pathpool,univdiagpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d\n",adj));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
				hitlistpool_trace(__FILE__,__LINE__));
	
      } else if ((adj = Genomebits_indel_solve_low(&trimpos,&nmismatches_i,
						   univdiagonal,querylength,/*pos5*/0,/*pos3*/nosplice_trimpos,
						   query_compress,mismatch_positions_alloc,
						   genomebits,genomebits_alt,plusp,genestrand)) != 0 &&
		 (newpath = attach_indel_qstart(path,/*low_diagonal*/univdiagonal - adj,/*low_qstart*/trimpos,
						chroffset,chrhigh,querylength,stage1->indelinfo,
						query_compress,plusp,genestrand,
						intlistpool,univcoordlistpool,listpool,pathpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d\n",adj));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
				hitlistpool_trace(__FILE__,__LINE__));
	
      } else if ((indel_univdiagonal =
		  Spliceends_indel_qstart(nosplice_trimpos,univdiagonal,querylength,chroffset,chrhigh,
					  plusp,genestrand,localdb_nmismatches_allowed,
					  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					  localdb_tally,localdb_overflowp,max_nblocks,stage1,
					  univdiagpool,query_compress,queryptr)) != 0 &&
		 (newpath = attach_indel_qstart(path,/*low_diagonal*/indel_univdiagonal,/*low_qstart*/0,
						chroffset,chrhigh,querylength,stage1->indelinfo,
						query_compress,plusp,genestrand,
						intlistpool,univcoordlistpool,listpool,pathpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d succeeds\n",univdiagonal - indel_univdiagonal));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_fwd(newpath)
				hitlistpool_trace(__FILE__,__LINE__));
      } else {
	debug13(printf("=> Indel fails\n"));
      }
    }
    /* Done with attempt to extend */
      
    if (newpaths == NULL) {
      /* Terminate this path with splice end and/or mismatches to the end.  No recursion. */

      /* Splice end */
      splice_path = (T) NULL;
      if (medial_splice_trimpos_0 == -1) {
	/* No splice end */
      } else if (0 && medial_splice_trimpos_0 == Intlist_head(path->endpoints)) {
	/* No change, but want to set splice5p so this path is considered complete and undergoes repair */
      } else if (ambig_prob_5 > 0.80) {
	/* Change endpoint to the splice trimpos */
	/* assert(nspliceends == 1); */
	splice_path = Path_copy(path,intlistpool,univcoordlistpool,listpool,
				pathpool,vectorpool,transcriptpool,hitlistpool);

	debug13(printf("(1) Changing endpoint from %d",Intlist_head(path->endpoints)));
	if (medial_splice_trimpos_0 > Intlist_head(path->endpoints) && Intlist_head(path->nmismatches) == 0) {
	  /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
	  Intlist_head_set(splice_path->endpoints,medial_splice_trimpos_0);
	} else {
	  Intlist_head_set(splice_path->endpoints,medial_splice_trimpos_0);
	  Intlist_head_set(splice_path->nmismatches,medial_splice_nmismatches_0);
	  Intlist_head_set(splice_path->ref_nmismatches,medial_splice_nmismatches_0);
	}
	splice_path->splice5p = splice5p;
	splice_path->splicetype5 = splicetype5;
	splice_path->ambig_prob_5 = ambig_prob_5;
	debug13(printf(" to splice end %d (splice5p %d) with %d mismatches\n",
		       medial_splice_trimpos_0,splice5p,Intlist_head(splice_path->nmismatches)));
	debug13(Path_print(splice_path));

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	if (Intlist_head(splice_path->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
			 Intlist_head(splice_path->endpoints)));
	  debug13(Path_print(splice_path));
	  debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(splice_path));
	  splice_path->completep = Path_completep(splice_path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	  *complete_qstart_paths = Hitlist_push(*complete_qstart_paths,hitlistpool,(void *) Path_expect_fwd(splice_path)
						hitlistpool_trace(__FILE__,__LINE__));
	} else {	
	  debug13(printf("splice path from compute_qstart_local is unextended because of qstart %d => putting into unextended\n",
			 Intlist_head(splice_path->endpoints)));
	  debug20(printf("(5) Putting into unextended qstart paths: ")); debug20(Path_print(splice_path));
	  assert(splice_path->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,(void *) Path_expect_fwd(splice_path)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
		       Intlist_head(splice_path->endpoints)));
	debug13(Path_print(splice_path));
	debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(splice_path));
	splice_path->completep = Path_completep(splice_path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(splice_path)
				     hitlistpool_trace(__FILE__,__LINE__));
#endif

	
      }

      /* Extend to end */
      end_path = (T) NULL;
      if (nosplice_trimpos == Intlist_head(path->endpoints)) {
	/* No change */
      } else if (nosplice_trimpos == Intlist_second_value(path->endpoints)) {
	debug13(printf("Segment is not great, since nosplice_trimpos erases segment, so not trimming\n"));
      
      } else if (Univcoordlist_head(path->univdiagonals) - querylength + nosplice_trimpos >= chrhigh) {
	debug13(printf("Attempt to trim beyond end of chromosome (%u > %u), so not trimming\n",
		       Univcoordlist_head(path->univdiagonals) - querylength + trimpos,chrhigh));
      } else {
	/* Change endpoint to the non-splice trimpos */
	end_path = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool);

	debug13(printf("(2) Changing endpoint from %d",Intlist_head(path->endpoints)));
	if (nosplice_trimpos > Intlist_head(path->endpoints) && Intlist_head(path->nmismatches) == 0) {
	  /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
	  Intlist_head_set(end_path->endpoints,nosplice_trimpos);
	} else {
	  Intlist_head_set(end_path->endpoints,nosplice_trimpos);
	  Intlist_head_set(end_path->nmismatches,nosplice_nmismatches);
	  Intlist_head_set(end_path->ref_nmismatches,nosplice_nmismatches);
	}
	end_path->splice5p = false;
	end_path->splicetype5 = NO_SPLICE;
	end_path->ambig_prob_5 = 0.0;
	debug13(printf(" to nosplice end %d (splice5p %d) with %d mismatches\n",
		       nosplice_trimpos,splice5p,Intlist_head(end_path->endpoints)));
	debug13(Path_print(end_path));

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	/* Previously put into complete always, but this makes it hard to find fusions */
	if (Intlist_head(end_path->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
			 Intlist_head(end_path->endpoints)));
	  debug13(Path_print(end_path));
	  debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(end_path));
	  end_path->completep = true;
	  *complete_qstart_paths = Hitlist_push(*complete_qstart_paths,hitlistpool,(void *) Path_expect_fwd(end_path)
						hitlistpool_trace(__FILE__,__LINE__));

	} else {	
	  debug13(printf("splice path from compute_qstart_local is unextended because of qstart %d => putting into unextended\n",
			 Intlist_head(end_path->endpoints)));
	  debug20(printf("(5) Putting into unextended qstart paths: ")); debug20(Path_print(end_path));
	  assert(end_path->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,(void *) Path_expect_fwd(end_path)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
		       Intlist_head(end_path->endpoints)));
	debug13(Path_print(end_path));
	debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(end_path));
	end_path->completep = Path_completep(end_path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(end_path)
				     hitlistpool_trace(__FILE__,__LINE__));
#endif


      }

      if (splice_path == NULL && end_path == NULL) {
	newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			    pathpool,vectorpool,transcriptpool,hitlistpool);

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	/* Previously put into complete always, but this makes it hard to find fusions */
	if (Intlist_head(newpath->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
			 Intlist_head(newpath->endpoints)));
	  debug13(Path_print(newpath));
	  debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(newpath));
	  newpath->completep = true;
	  *complete_qstart_paths = Hitlist_push(*complete_qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
						hitlistpool_trace(__FILE__,__LINE__));

	} else {	
	  debug13(printf("splice path from compute_qstart_local is unextended because of qstart %d => putting into unextended\n",
			 Intlist_head(newpath->endpoints)));
	  debug20(printf("(5) Putting into unextended qstart paths: ")); debug20(Path_print(newpath));
	  assert(newpath->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("splice path from compute_qstart_local is reasonably complete because of qstart %d => putting into complete\n",
		       Intlist_head(newpath->endpoints)));
	debug13(Path_print(newpath));
	debug20(printf("(4) Putting into complete qstart paths: ")); debug20(Path_print(newpath));
	newpath->completep = Path_completep(newpath,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qstart_paths = Hitlist_push(*qstart_paths,hitlistpool,(void *) Path_expect_fwd(newpath)
				     hitlistpool_trace(__FILE__,__LINE__));
#endif
      }

      return;

    } else {
      /* Extensions were added, so recurse */
      debug13(printf("compute_qstart_local recursing over %d child paths from depth %d, orig_qstart %d\n",
		     List_length(newpaths),depth,orig_qstart));
      Accel_profile_pathsolve_children("qstart",depth,List_length(newpaths));
      for (p = newpaths; p != NULL; p = List_next(p)) {
	parent = (T) List_head(p);
	child_qstart = Path_qstart(parent);
	if (child_qstart >= orig_qstart) {
	  debug13(printf("Skipping compute_qstart_local non-progress child: depth %d -> %d, child qstart %d >= orig qstart %d\n",
			 depth,depth+1,child_qstart,orig_qstart));
	} else {
	  debug13(printf("compute_qstart_local progress child: depth %d -> %d, child qstart %d < orig qstart %d\n",
			 depth,depth+1,child_qstart,orig_qstart));
	  compute_qstart_local(&(*qstart_paths),depth+1,parent,queryptr,querylength,
			       mismatch_positions_alloc,novel_univdiags_alloc,
			       novel_diagonals_alloc,localdb_alloc,
			       localdb_tally,localdb_overflowp,max_nblocks,
			       stage1,knownsplicing,knownindels,
			       query_compress,chrnum,chroffset,chrhigh,intlistpool,
			       univcoordlistpool,listpool,pathpool,transcriptpool,
			       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			       localdb_nmismatches_allowed,plusp,genestrand,
			       sensedir,innerp,find_splices_p,salvagep);
	}
	Path_free(&parent,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      return;
    }
  }
}


static void
compute_qstart_local (List_T *qstart_paths,
		      int depth, Path_T path, char *queryptr, int querylength,
		      int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		      Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		      Compress_T query_compress, Chrnum_T chrnum, Univcoord_T chroffset,
		      Univcoord_T chrhigh, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
		      Hitlistpool_T hitlistpool, Spliceendsgen_T spliceendsgen,
		      int localdb_nmismatches_allowed, bool plusp, int genestrand,
		      int sensedir, bool innerp, bool find_splices_p, bool salvagep) {
#ifdef GSNAP_ACCEL_PROFILE_STATS
  List_T old_qstart_paths = *qstart_paths;
  unsigned long long state_id = Accel_profile_pathsolve_call("qstart",depth);
  unsigned long long input_path_hash = pathsolve_path_fingerprint(path);
  unsigned long long frontier_hash = pathsolve_frontier_fingerprint(
    path,"qstart",depth,innerp,find_splices_p,salvagep);
#endif

  compute_qstart_local_work(qstart_paths,depth,path,queryptr,querylength,
			    mismatch_positions_alloc,novel_univdiags_alloc,
			    novel_diagonals_alloc,localdb_alloc,localdb_tally,
			    localdb_overflowp,max_nblocks,stage1,knownsplicing,
			    knownindels,query_compress,chrnum,chroffset,chrhigh,
			    intlistpool,univcoordlistpool,listpool,pathpool,
			    transcriptpool,univdiagpool,vectorpool,hitlistpool,
			    spliceendsgen,localdb_nmismatches_allowed,plusp,
			    genestrand,sensedir,innerp,find_splices_p,salvagep);
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_pathsolve_result(state_id,"qstart",depth,input_path_hash,frontier_hash,
			   *qstart_paths,old_qstart_paths);
#endif
}

/* Recursively adds diagonals to the end of a path */
/* Either returns NULL, which means caller should use the given path (or a copy)
   as the result, or returns newpaths in qend_paths, which of which is different
   from path, where the caller needs to free path, if desired */

/* unextended_paths are all fwd.  complete_paths are all rev */
static void
compute_qend_local (List_T *qend_paths, int depth, T path, char *queryptr, int querylength,
		    int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels, Compress_T query_compress, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		    Spliceendsgen_T spliceendsgen, int localdb_nmismatches_allowed,
		    bool plusp, int genestrand, int sensedir, bool innerp, bool find_splices_p,
		    bool salvagep, Path_segment_cache_T segment_cache,
		    bool path_qend_local_cache_p);

static void
compute_qend_local_work (List_T *qend_paths, int depth, T path, char *queryptr, int querylength,
			 int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
			 Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			 int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			 Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels, Compress_T query_compress,
			 Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			 Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			 Spliceendsgen_T spliceendsgen, int localdb_nmismatches_allowed,
			 bool plusp, int genestrand, int sensedir, bool innerp, bool find_splices_p,
			 bool salvagep, Path_segment_cache_T segment_cache,
			 bool path_qend_local_cache_p) {
  List_T newpaths = NULL, p;
  T splice_path, end_path, newpath, parent;
  Univcoord_T indel_univdiagonal, univdiagonal;
  int qend, pos5, exon_origin, indel_pos, trimpos, adj;
  int orig_qend, child_qend;
  int nmismatches_j;

  bool splice3p;
  Splicetype_T splicetype3;
  int nspliceends, i, j;
  Spliceends_T spliceends;
  int nosplice_trimpos, medial_splice_trimpos_0; /* Most proximal medial_splice_trimpos */
  double ambig_prob_3;
  int nosplice_nmismatches, medial_splice_nmismatches_0;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long spliceends_key_hash, spliceends_result_hash;
#endif


  debug13(printf("Entering compute_qend_local, salvage %d, at depth %d with path:\n",salvagep,depth));
  debug13(Path_print(path));
  debug13(printf("\n"));
  Path_expect_rev(path);

  orig_qend = Path_qend(path);

  if (path->qend_alts != NULL) {
    /* Not possible to add localdb because of the ambiguity */
    debug13(printf("compute_qend_local is complete because of qend_alts\n"));
    path->completep = Path_completep(path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(6) Putting into complete qend paths: ")); debug20(Path_print(newpath));
    *qend_paths = Hitlist_push(*qend_paths,hitlistpool,(void *) Path_expect_rev(newpath)
			       hitlistpool_trace(__FILE__,__LINE__));
    return;
	
  } else if ((qend = Intlist_head(path->endpoints)) == querylength) {
    /* Already at end */
    debug13(printf("compute_qend_local is complete because of qend == querylength\n"));
    path->completep = Path_completep(path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(7) Putting into complete qend paths: ")); debug20(Path_print(newpath));
    *qend_paths = Hitlist_push(*qend_paths,hitlistpool,(void *) Path_expect_rev(newpath)
			       hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else if (depth > max_recursive_ends) {
    /* Too much recursive depth */
    debug13(printf("compute_qend_local is unextended because depth %d > max %d\n",
		   depth,max_recursive_ends));
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(8) Putting into unextended qend paths: ")); debug20(Path_print(newpath));
    assert(newpath->sensedir == sensedir);
    *qend_paths = Hitlist_push(*qend_paths,hitlistpool,
			       (void *) Path_expect_rev(newpath)
			       hitlistpool_trace(__FILE__,__LINE__));
    return;

  } else if (splicingp == false) {
    /* Not able to try splicing */
    newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			pathpool,vectorpool,transcriptpool,hitlistpool);
    debug20(printf("(8) Putting into unextended qend paths: ")); debug20(Path_print(newpath));
    assert(newpath->sensedir == sensedir);
    *qend_paths = Hitlist_push(*qend_paths,hitlistpool,
			       (void *) Path_expect_rev(newpath)
			       hitlistpool_trace(__FILE__,__LINE__));
  } else {
    /* Attempt to extend */
    debug13(printf("compute_qend_local extension state: depth %d, orig_qend %d\n",
		   depth,orig_qend));

    /* pos5 is the start of the current segment.  qend is the
       approximate end of the current segment.  qend to pos3
       (querylength) is where we look for a new segment. */
    pos5 = Intlist_second_value(path->endpoints);
    if (path->junctions != NULL) {
      pos5 += Junction_ninserts((Junction_T) List_head(path->junctions));
    }
    univdiagonal = Univcoordlist_head(path->univdiagonals);
    /* left = univdiagonal - (Univcoord_T) querylength; */

    exon_origin = Path_exon_origin(path);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    spliceends_key_hash = spliceends_call_fingerprint(
      "qend",sensedir,univdiagonal,querylength,pos5,exon_origin,
      chrnum,chroffset,chrhigh,plusp,genestrand,localdb_nmismatches_allowed,
      innerp,salvagep,max_nblocks,stage1,knownsplicing,spliceendsgen,
      query_compress,queryptr,genomebits,genomebits_alt,find_splices_p);
#endif
    debug13(printf("(SPLICE) Calling Spliceends_trimmed_qends with univdiagonal %u, pos5 %d, exon_origin %d\n",
		   univdiagonal,pos5,exon_origin));
    nspliceends = Spliceends_trimmed_qends(&spliceends,&nosplice_trimpos,&medial_splice_trimpos_0,
					   &nosplice_nmismatches,&medial_splice_nmismatches_0,
					   &splice3p,&splicetype3,&ambig_prob_3,
					   sensedir,univdiagonal,querylength,
					   /*qstart*/pos5,/*qend:querylength,*/exon_origin,chrnum,chroffset,chrhigh,
					   plusp,genestrand,localdb_nmismatches_allowed,innerp,salvagep,
					   mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,stage1,
					   knownsplicing,univdiagpool,vectorpool,spliceendsgen,
					   query_compress,queryptr,genomebits,genomebits_alt,find_splices_p);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    spliceends_result_hash = spliceends_result_fingerprint(
      spliceends,nspliceends,nosplice_trimpos,medial_splice_trimpos_0,
      nosplice_nmismatches,medial_splice_nmismatches_0,splice3p,splicetype3,
      ambig_prob_3);
    Accel_profile_spliceends_result("qend",spliceends_key_hash,nspliceends,
				    spliceends_result_hash);
#endif
    /* path->spliceends_qend = pos5; */
    debug13(printf("Spliceends_trimmed_qends returning %d spliceends, nosplice %d (%d nmismatches), medial_splice %d (%d nmismatches)\n",
		   nspliceends,nosplice_trimpos,nosplice_nmismatches,medial_splice_trimpos_0,medial_splice_nmismatches_0));
    
    if (nspliceends == 0) {
      /* No trimming.  Means that trimpos < pos5, so a bad segment */
      debug13(printf("=> No trimming\n"));
      
    } else if (spliceends == NULL) {
      debug13(printf("=> No partners found\n"));
      
    } else if (nosplice_trimpos == querylength) {
      debug13(printf("=> Actually, not solving for splice because nosplice_trimpos is %d\n",nosplice_trimpos));
      
    } else {
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
      printf("=> Got %d spliceends\n",nspliceends);
      for (int i = 0; i < nspliceends; i++) {
	printf("splice_qpos %d => partner %u, %d..%d, splicedist %u, mismatches %d, prob %f, trimpos %d\n",
	       spliceends->splice_qpos[i],spliceends->partners[i],
	       spliceends->distal_univdiags[i]->qstart,spliceends->distal_univdiags[i]->qend,
	       (spliceends->partners[i] - spliceends->splice_qpos[i] + querylength) - univdiagonal,
	       spliceends->distal_univdiags[i]->nmismatches,spliceends->distal_probs[i],
	       spliceends->distal_univdiags[i]->qend);
      }
#endif

      i = 0;
      while (i < nspliceends) {
	/* Assumes that spliceends are sorted primarily by splice_qpos */
	/* ambig_prob_3 = spliceends->medial_probs[i]; */
	qend = spliceends->splice_qpos[i]; /* was high_qend */
	j = i + 1;
	while (j < nspliceends && spliceends->splice_qpos[j] == qend) {
	  j++;
	}
	
	/* Handle a cluster of solutions for qend */
	debug13(printf("Splice qpos: %d, medial prob %f\n",qend,ambig_prob_3));
	newpaths = List_append(newpaths,
			       multiadd_splice_qends(path,univdiagonal,spliceends->splicetype,
						     spliceends->boundedp,/*splice_qpos*/qend,/*anchor_qpos*/pos5,
						     spliceends->medial_nmismatches[i],spliceends->medial_probs[i],
						     spliceends->medial_support[i],
						     &(spliceends->partners[i]),&(spliceends->distal_support[i]),
						     &(spliceends->distal_univdiags[i]),&(spliceends->distal_probs[i]),
						     /*npartners*/(j - i),
						     stage1,knownsplicing,query_compress,queryptr,
						     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						     localdb_tally,localdb_overflowp,max_nblocks,
						     localdb,localdb_nmismatches_allowed,plusp,querylength,innerp,
						     intlistpool,univcoordlistpool,listpool,
						     hitlistpool,pathpool,transcriptpool,vectorpool,sensedir));
	i = j;
      }
    }

    /* Spliceends_free(&spliceends); */
    Spliceendsgen_return(spliceendsgen,&spliceends);

    if (newpaths != NULL) {
      debug13(printf("compute_qend_local has %d splice child paths; still trying end indel\n",
		     List_length(newpaths)));
    }

    if (1 || newpaths == NULL) {
      /* Always try an indel, since splice may be incorrect */
      /* Splice end without partner, or an end indel */
      /* Try end indel */
      debug13(printf("Path before trying end indel:\n"));
      debug13(Path_print(path));
      debug13(printf("\n"));

      /* Computed above */
      /* pos5 = Intlist_head(Intlist_next(path->endpoints)); */
      /* univdiagonal = Univcoordlist_head(path->univdiagonals); */
      /* left = univdiagonal - (Univcoord_T) querylength; */
      
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
      if (nosplice_trimpos != querylength) {
	printf("(INDEL) Calling Genomebits_indel_solve_high and Spliceends_indel_qend with %d..%d\n",nosplice_trimpos,querylength);
      }
#endif
  
      debug13(printf("nosplice_trimpos at qend is %d\n",nosplice_trimpos));
      if (nosplice_trimpos == querylength) {
	/* Already extended to the end */
	debug13(printf("Actually, not solving for indel because nosplice_trimpos is %d\n",nosplice_trimpos));
	
      } else if (knownindels != NULL &&
		 (adj = Knownindels_find_lowest(&indel_pos,knownindels,univdiagonal,querylength,
						/*pos5*/nosplice_trimpos-1,/*pos3*/querylength)) != 0 &&
		 (newpath = attach_indel_qend_simple(adj,path,indel_pos,univdiagonal,querylength,sensedir,
						     plusp,genestrand,mismatch_positions_alloc,
						     stage1,knownsplicing,spliceendsgen,query_compress,queryptr,
						     genomebits,genomebits_alt,chrnum,chroffset,chrhigh,find_splices_p,
						     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						     localdb_tally,localdb_overflowp,max_nblocks,
						     intlistpool,univcoordlistpool,listpool,
						     pathpool,univdiagpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d\n",adj));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
				hitlistpool_trace(__FILE__,__LINE__));
	
      } else if ((adj = Genomebits_indel_solve_high(&trimpos,&nmismatches_j,
						    univdiagonal,querylength,/*pos5*/nosplice_trimpos,/*pos3*/querylength,
						    query_compress,mismatch_positions_alloc,
						    genomebits,genomebits_alt,plusp,genestrand)) != 0 &&
		  (newpath = attach_indel_qend(path,/*high_diagonal*/univdiagonal + adj,/*high_qend*/trimpos,
					       chrhigh,querylength,stage1->indelinfo,
					       query_compress,plusp,genestrand,
					       intlistpool,univcoordlistpool,listpool,pathpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d\n",adj));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
				hitlistpool_trace(__FILE__,__LINE__));

      } else if ((indel_univdiagonal =
		  Spliceends_indel_qend(nosplice_trimpos,univdiagonal,querylength,chroffset,chrhigh,
					plusp,genestrand,localdb_nmismatches_allowed,
					novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					localdb_tally,localdb_overflowp,max_nblocks,stage1,
					univdiagpool,query_compress,queryptr)) != 0 &&
		 (newpath = attach_indel_qend(path,/*high_diagonal*/indel_univdiagonal,/*high_qend*/querylength,
					      chrhigh,querylength,stage1->indelinfo,
					      query_compress,plusp,genestrand,
					      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool)) != NULL) {
	debug13(printf("=> Indel of adj %d succeeds\n",indel_univdiagonal - univdiagonal));
	newpaths = Hitlist_push(newpaths,hitlistpool,(void *) Path_expect_rev(newpath)
				hitlistpool_trace(__FILE__,__LINE__));
	
      } else {
	debug13(printf("=> Indel fails\n"));
      }
    }
    /* Done with attempt to extend */
    
    if (newpaths == NULL) {
      /* Terminate this path with splice end and/or mismatches to the end.  No recursion. */

      /* Splice end */
      splice_path = (T) NULL;
      if (medial_splice_trimpos_0 == -1) {
	/* No splice end */
      } else if (0 && medial_splice_trimpos_0 == Intlist_head(path->endpoints)) {
	/* No change, but want to set splice3p so this path is considered complete and undergoes repair */
      } else if (path->junctions != NULL &&
		 Intlist_second_value(path->endpoints) +
		 Junction_ninserts((Junction_T) List_head(path->junctions)) >= medial_splice_trimpos_0) {
	debug13(printf("(3) End indel %d + ins:%d vs changing endpoint from %d to splice end %d\n",
		       Intlist_second_value(path->endpoints),Junction_ninserts((Junction_T) List_head(path->junctions)),
		       Intlist_head(path->endpoints),medial_splice_trimpos_0));
      } else if (ambig_prob_3 > 0.80) {
	/* Change endpoint to the splice trimpos */
	/* assert(nspliceends == 1); */
	splice_path = Path_copy(path,intlistpool,univcoordlistpool,listpool,
				pathpool,vectorpool,transcriptpool,hitlistpool);
	debug13(printf("(3) Changing endpoint from %d",Intlist_head(path->endpoints)));
	if (medial_splice_trimpos_0 < Intlist_head(path->endpoints) && Intlist_head(path->nmismatches) == 0) {
	  /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
	  Intlist_head_set(splice_path->endpoints,medial_splice_trimpos_0);
	} else {
	  Intlist_head_set(splice_path->endpoints,medial_splice_trimpos_0);
	  Intlist_head_set(splice_path->nmismatches,medial_splice_nmismatches_0);
	  Intlist_head_set(splice_path->ref_nmismatches,medial_splice_nmismatches_0);
	}
	splice_path->splice3p = splice3p;
	splice_path->splicetype3 = splicetype3;
	splice_path->ambig_prob_3 = ambig_prob_3;
	debug13(printf(" to splice end %d (splice3p %d) with %d mismatches\n",
		       medial_splice_trimpos_0,splice3p,Intlist_head(splice_path->nmismatches)));
	debug13(Path_print(splice_path));

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	/* Previously put into complete always, but this makes it hard to find fusions */
	if (querylength - Intlist_head(splice_path->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
			 Intlist_head(splice_path->endpoints)));
	  debug13(Path_print(splice_path));
	  debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(splice_path));
	  splice_path->completep = true;
	  *complete_qend_paths = Hitlist_push(*complete_qend_paths,hitlistpool,(void *) Path_expect_rev(splice_path)
					      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  debug13(printf("compute_qend_local is unextended because of qend %d => putting into unextended\n",
			 Intlist_head(splice_path->endpoints)));
	  debug20(printf("(10) Putting into unextended qend paths: ")); debug20(Path_print(splice_path));
	  assert(splice_path->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,
					   (void *) Path_reverse(splice_path,/*expect_fwd_p*/true)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
		       Intlist_head(splice_path->endpoints)));
	debug13(Path_print(splice_path));
	debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(splice_path));
	splice_path->completep = Path_completep(splice_path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qend_paths = Hitlist_push(*qend_paths,hitlistpool,(void *) Path_expect_rev(splice_path)
				   hitlistpool_trace(__FILE__,__LINE__));
#endif
      }
	
      /* Extend to end */
      end_path = (T) NULL;
      if (nosplice_trimpos == Intlist_head(path->endpoints)) {
	/* No change */
      } else if (nosplice_trimpos == Intlist_second_value(path->endpoints)) {
	debug13(printf("Segment is not great, since nosplice_trimpos erases segment, so not trimming\n"));
	
      } else if (Univcoordlist_head(path->univdiagonals) - querylength + nosplice_trimpos < chroffset) {
	debug13(printf("Attempt to trim beyond start of chromosome (%u < %u), so not trimming\n",
		       Univcoordlist_head(path->univdiagonals) - querylength + trimpos,chroffset));

      } else if (path->junctions != NULL &&
		 Intlist_second_value(path->endpoints) +
		 Junction_ninserts((Junction_T) List_head(path->junctions)) >= nosplice_trimpos) {
	debug13(printf("(4) End indel %d + ins:%d vs changing endpoint from %d to splice end %d\n",
		       Intlist_second_value(path->endpoints),Junction_ninserts((Junction_T) List_head(path->junctions)),
		       Intlist_head(path->endpoints),nosplice_trimpos));
      } else {
	/* Change endpoint to the non-splice trimpos */
	end_path = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool);

	debug13(printf("(4) Changing endpoint from %d",Intlist_head(path->endpoints)));
	if (nosplice_trimpos < Intlist_head(path->endpoints) && Intlist_head(path->nmismatches) == 0) {
	  /* Shorter segment in region with no nmismatches, so can keep nmismatches being 0 */
	  Intlist_head_set(end_path->endpoints,nosplice_trimpos);
	} else {
	  Intlist_head_set(end_path->endpoints,nosplice_trimpos);
	  Intlist_head_set(end_path->nmismatches,nosplice_nmismatches);
	  Intlist_head_set(end_path->ref_nmismatches,nosplice_nmismatches);
	}
	end_path->splice3p = false;
	end_path->splicetype3 = NO_SPLICE;
	end_path->ambig_prob_3 = 0.0;
	debug13(printf(" to nosplice end %d (splice3p %d) with %d mismatches\n",
		       nosplice_trimpos,splice3p,Intlist_head(end_path->nmismatches)));
	debug13(Path_print(end_path));

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	/* Previously put into complete always, but this makes it hard to find fusions */
	if (querylength - Intlist_head(end_path->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
			 Intlist_head(end_path->endpoints)));
	  debug13(Path_print(end_path));
	  debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(end_path));
	  end_path->completep = true;
	  *complete_qend_paths = Hitlist_push(*complete_qend_paths,hitlistpool,(void *) Path_expect_rev(end_path)
					      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  debug13(printf("compute_qend_local is unextended because of qend %d => putting into unextended\n",
			 Intlist_head(end_path->endpoints)));
	  debug20(printf("(10) Putting into unextended qend paths: ")); debug20(Path_print(end_path));
	  assert(end_path->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,
					   (void *) Path_reverse(end_path,/*expect_fwd_p*/true)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
		       Intlist_head(end_path->endpoints)));
	debug13(Path_print(end_path));
	debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(end_path));
	end_path->completep = Path_completep(end_path,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qend_paths = Hitlist_push(*qend_paths,hitlistpool,(void *) Path_expect_rev(end_path)
				   hitlistpool_trace(__FILE__,__LINE__));
#endif
      }

      if (splice_path == NULL && end_path == NULL) {
	newpath = Path_copy(path,intlistpool,univcoordlistpool,listpool,
			    pathpool,vectorpool,transcriptpool,hitlistpool);

#if 0
	/* No longer filtering based on 8 bp at ends, which misses some reads and misses the chance to combine qstart and qend */
	/* Previously put into complete always, but this makes it hard to find fusions */
	if (querylength - Intlist_head(newpath->endpoints) <= 8) {
	  /* Reasonably complete */
	  debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
			 Intlist_head(newpath->endpoints)));
	  debug13(Path_print(newpath));
	  debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(newpath));
	  newpath->completep = true;
	  *complete_qend_paths = Hitlist_push(*complete_qend_paths,hitlistpool,(void *) Path_expect_rev(newpath)
					      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  debug13(printf("compute_qend_local is unextended because of qend %d => putting into unextended\n",
			 Intlist_head(newpath->endpoints)));
	  debug20(printf("(10) Putting into unextended qend paths: ")); debug20(Path_print(newpath));
	  assert(newpath->sensedir == sensedir);
	  *unextended_paths = Hitlist_push(*unextended_paths,hitlistpool,
					   (void *) Path_reverse(newpath,/*expect_fwd_p*/true)
					   hitlistpool_trace(__FILE__,__LINE__));
	}
#else
	debug13(printf("compute_qend_local is reasonably complete because of qend %d => putting into complete\n",
		       Intlist_head(newpath->endpoints)));
	debug13(Path_print(newpath));
	debug20(printf("(9) Putting into complete qend paths: ")); debug20(Path_print(newpath));
	newpath->completep = Path_completep(newpath,/*endtrim_allowed*/8,/*allow_ambig_p*/true);
	*qend_paths = Hitlist_push(*qend_paths,hitlistpool,(void *) Path_expect_rev(newpath)
				   hitlistpool_trace(__FILE__,__LINE__));
#endif
      }

      return;

    } else {
      /* Extensions were added, so recurse */
      debug13(printf("compute_qend_local recursing over %d child paths from depth %d, orig_qend %d\n",
		     List_length(newpaths),depth,orig_qend));
      Accel_profile_pathsolve_children("qend",depth,List_length(newpaths));
      for (p = newpaths; p != NULL; p = List_next(p)) {
	parent = (T) List_head(p);
	child_qend = Path_qend(parent);
	if (child_qend <= orig_qend) {
	  debug13(printf("Skipping compute_qend_local non-progress child: depth %d -> %d, child qend %d <= orig qend %d\n",
			 depth,depth+1,child_qend,orig_qend));
	} else {
	  debug13(printf("compute_qend_local progress child: depth %d -> %d, child qend %d > orig qend %d\n",
			 depth,depth+1,child_qend,orig_qend));
	  compute_qend_local(&(*qend_paths),depth+1,parent,queryptr,querylength,
			     mismatch_positions_alloc,novel_univdiags_alloc,
			     novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,
			     stage1,knownsplicing,knownindels,
			     query_compress,chrnum,chroffset,chrhigh,intlistpool,
			     univcoordlistpool,listpool,pathpool,transcriptpool,
			     univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			     localdb_nmismatches_allowed,plusp,genestrand,
			     sensedir,innerp,find_splices_p,salvagep,segment_cache,
			     path_qend_local_cache_p);
	}
	Path_free(&parent,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      return;
    }
  }
}

static void
compute_qend_local (List_T *qend_paths, int depth, T path, char *queryptr, int querylength,
		    int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels, Compress_T query_compress,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		    Spliceendsgen_T spliceendsgen, int localdb_nmismatches_allowed,
		    bool plusp, int genestrand, int sensedir, bool innerp, bool find_splices_p,
		    bool salvagep, Path_segment_cache_T segment_cache,
		    bool path_qend_local_cache_p) {
  List_T old_qend_paths = *qend_paths;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long state_id = Accel_profile_pathsolve_call("qend",depth);
  unsigned long long input_path_hash = pathsolve_path_fingerprint(path);
  unsigned long long frontier_hash = pathsolve_frontier_fingerprint(
    path,"qend",depth,innerp,find_splices_p,salvagep);
#endif
  Path_segment_cache_T qend_local_cache =
    path_qend_local_cache_p == true ? segment_cache : (Path_segment_cache_T) NULL;

  if (qend_local_cache_lookup(&(*qend_paths),qend_local_cache,path,depth,querylength,
			      max_nblocks,localdb_nmismatches_allowed,chrnum,chroffset,chrhigh,
			      plusp,genestrand,sensedir,innerp,find_splices_p,salvagep,
			      intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			      hitlistpool,transcriptpool) == true) {
    Accel_profile_pathsolve_cache("qend_local",/*hitp*/true);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    profile_pathsolve_result(state_id,"qend",depth,input_path_hash,frontier_hash,
			     *qend_paths,old_qend_paths);
#endif
    return;
  }
  Accel_profile_pathsolve_cache("qend_local",/*hitp*/false);

  compute_qend_local_work(&(*qend_paths),depth,path,queryptr,querylength,
			  mismatch_positions_alloc,novel_univdiags_alloc,
			  novel_diagonals_alloc,localdb_alloc,
			  localdb_tally,localdb_overflowp,max_nblocks,
			  stage1,knownsplicing,knownindels,query_compress,
			  chrnum,chroffset,chrhigh,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
			  hitlistpool,spliceendsgen,localdb_nmismatches_allowed,
			  plusp,genestrand,sensedir,innerp,find_splices_p,salvagep,
			  segment_cache,path_qend_local_cache_p);

  qend_local_cache_insert(qend_local_cache,path,old_qend_paths,*qend_paths,depth,querylength,
			  max_nblocks,localdb_nmismatches_allowed,chrnum,chroffset,chrhigh,
			  plusp,genestrand,sensedir,innerp,find_splices_p,salvagep,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			  hitlistpool,transcriptpool);
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_pathsolve_result(state_id,"qend",depth,input_path_hash,frontier_hash,
			   *qend_paths,old_qend_paths);
#endif
  return;
}


#ifdef REMAP_PATHS
static List_T
try_repairs (int *found_score, int *end_score, Path_T path,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     Shortread_T queryseq, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	     Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
	     Transcriptpool_T transcriptpool, Pathpool_T pathpool,
	     Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T repaired_paths = NULL;
  T newpath;
  List_T repairs, r;
  Repair_T repair;
  int desired_genestrand;
  Compress_T query_compress;

#if 0
  if (path->sensedir == SENSE_FORWARD) {
    desired_genestrand = (path->plusp == true) ? +1 : -1;
  } else {
    desired_genestrand = (path->plusp == true) ? -1 : +1;
  }
#endif

  if (path->plusp == true) {
    query_compress = query_compress_fwd;
    desired_genestrand = (path->sensedir == SENSE_FORWARD) ? +1 : -1;
  } else {
    query_compress = query_compress_rev;
    desired_genestrand = (path->sensedir == SENSE_FORWARD) ? -1 : +1;
  }

  repairs = Transcript_remap_all(&path->transcripts,&path->invalid_transcripts,
				 path->endpoints,path->univdiagonals,path->junctions,
				 queryseq,path->querylength,path->plusp,
				 path->chrnum,path->chroffset,path->chrhigh,
				 uintlistpool,listpool,transcriptpool,desired_genestrand,
				 /*extend_qstart_p*/true,/*extend_qend_p*/true,/*repairp*/true);

#ifdef TRY_REPAIRS
  for (r = repairs; r != NULL; r = List_next(r)) {
    repair = (Repair_T) List_head(r);
    if ((newpath = Repair_path(&(*found_score),&(*end_score),repair,path,path->sensedir,
			       query_compress,query_compress_fwd,query_compress_rev,
			       intlistpool,univcoordlistpool,listpool,
			       pathpool,vectorpool,transcriptpool,hitlistpool)) == NULL) {
      Repair_free(&repair,listpool,transcriptpool,/*free_transcripts_p*/true);
    } else {
      repaired_paths = Hitlist_push(repaired_paths,hitlistpool,(void *) newpath
				    hitlistpool_trace(__FILE__,__LINE__));
      /* Repair transcripts got transferred onto newpath->transcripts */
      Repair_free(&repair,listpool,transcriptpool,/*free_transcripts_p*/false);
    } 
  }
#else
  for (r = repairs; r != NULL; r = List_next(r)) {
    repair = (Repair_T) List_head(r);
    Repair_free(&repair,listpool,transcriptpool,/*free_transcripts_p*/true);
  }
#endif

  Listpool_free_list(&repairs,listpool
		     listpool_trace(__FILE__,__LINE__));

  return repaired_paths;
}
#endif


#if 0
/* Used for SEGMENT1 if a Path_solve_univdiagonal does not succeed */
/* Calls combine_leftright_paths, which calls Path_eval_nmatches */
void
Path_solve_from_diagonals (int *found_score, int *end_score,

			   List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			   List_T *sense_paths, List_T *antisense_paths,

			   Univcoord_T middle_univdiagonal, Auxinfo_T auxinfo,

			   Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			   int localdb_nmismatches_allowed,

			   Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand,
			   bool paired_end_p, bool first_read_p,

			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			   Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p) {

  T path, newpath;

  int middle_diagonal_qstart = auxinfo->qstart;
  int middle_diagonal_qend = auxinfo->qend;
  /* int middle_nmismatches = auxinfo->nmismatches; */
  List_T qend_univdiags = auxinfo->right_univdiags;
  List_T qstart_univdiags = auxinfo->left_univdiags;

#ifdef TRY_REPAIRS
  List_T repaired_paths;
#endif

  List_T qstart_sense_paths = NULL, qstart_antisense_paths = NULL,
    qend_sense_paths = NULL, qend_antisense_paths = NULL, combined_sense_paths, combined_antisense_paths, p;
  Univdiag_T univdiag;
  bool qstart_innerp, qend_innerp;
  bool sense_completep, antisense_completep;
  int middle_rank;


  debug13(printf("Entered Path_solve_from_diagonals from method %s, find_splices_p %d, first_read_p %d, with middle_diagonal %u (#%d %u) %d..%d, %d qstart diagonals, and %d qend diagonals\n",
		 Method_string(method),find_splices_p,first_read_p,middle_univdiagonal,
		 chrnum,middle_univdiagonal - chroffset,middle_diagonal_qstart,middle_diagonal_qend,
		 List_length(qstart_univdiags),List_length(qend_univdiags)));

#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  printf("Start: first_read_p %d, plusp %d\n",first_read_p,plusp);
  printf("Start: %d complete sense paths\n",List_length(*sense_paths));
  printf("Start: %d complete antisense paths\n",List_length(*antisense_paths));
  printf("Start: %d unextended sense paths\n",List_length(*unextended_sense_paths));
  printf("Start: %d unextended antisense paths\n",List_length(*unextended_antisense_paths));

  for (p = qstart_univdiags; p != NULL; p = List_next(p)) {
    printf("qstart diagonal %u %d %d..%d\n",
	   ((Univdiag_T) List_head(p))->univdiagonal,
	   middle_univdiagonal - ((Univdiag_T) List_head(p))->univdiagonal,
	   ((Univdiag_T) List_head(p))->qstart,
	   ((Univdiag_T) List_head(p))->qend);
  }
  for (p = qend_univdiags; p != NULL; p = List_next(p)) {
    printf("qend diagonal %u %d %d..%d\n",
	   ((Univdiag_T) List_head(p))->univdiagonal,
	   ((Univdiag_T) List_head(p))->univdiagonal - middle_univdiagonal,
	   ((Univdiag_T) List_head(p))->qstart,
	   ((Univdiag_T) List_head(p))->qend);
  }
#endif


#if 0
  check_for_descending_qend(middle_diagonal_qend,middle_univdiagonal,qstart_diagonals);
  check_for_ascending_qstart(middle_diagonal_qstart,middle_univdiagonal,qend_diagonals);
#endif

  if (middle_univdiagonal < (Univcoord_T) querylength) {
    return;
  } else if (paired_end_p == false) {
    qstart_innerp = qend_innerp = false;
  } else if (first_read_p == plusp) {
    qstart_innerp = false;
    qend_innerp = true;
  } else {
    qstart_innerp = true;
    qend_innerp = false;
  }


  /* Qstart */
  path = Path_new_for_qstart_extension(middle_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
				       /*middle_nmismatches*/-1,plusp,first_read_p,genestrand,
				       /*sensedir*/SENSE_FORWARD,querylength,
				       method,chrnum,chroffset,chrhigh,
				       /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				       query_compress_fwd,query_compress_rev,
				       intlistpool,univcoordlistpool,pathpool);

  for (p = qstart_univdiags; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    debug13(printf("Processing qstart univdiag %u %u %d..%d\n",
		   univdiag->univdiagonal,middle_univdiagonal - univdiag->univdiagonal,univdiag->qstart,univdiag->qend));
    /* Previously checked for a large gap between middle_diagonal_qstart and univdiag->qend */

    if (univdiag->univdiagonal < chroffset + querylength) {
      /* Skip, because univdiagonal is in an earlier chromosome */

    } else {
      for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	if ((newpath = attach_unknown_qstart(path,/*low_univdiagonal*/univdiag->univdiagonal,
					       /*low_qstart*/univdiag->qstart,
					       chroffset,chrhigh,querylength,
					       stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					       stage1,query_compress,queryptr,
					       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					       localdb_tally,localdb_overflowp,max_nblocks,
					       localdb,localdb_nmismatches_allowed,plusp,genestrand,qstart_innerp,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					       vectorpool,hitlistpool,/*try_sensedir*/SENSE_FORWARD,middle_rank)) == NULL) {
	  break;
	} else {
	  compute_qstart_local(&qstart_sense_paths,/*depth*/0,newpath,queryptr,querylength,
			       mismatch_positions_alloc,novel_univdiags_alloc,
			       novel_diagonals_alloc,localdb_alloc,
			       localdb_tally,localdb_overflowp,max_nblocks,
			       stage1,knownsplicing,knownindels,
			       query_compress,chrnum,chroffset,chrhigh,intlistpool,
			       univcoordlistpool,listpool,pathpool,transcriptpool,
			       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			       localdb_nmismatches_allowed,plusp,genestrand,
			       /*try_sensedir*/SENSE_FORWARD,qstart_innerp,find_splices_p,/*salvagep*/false);
	  Path_free(&newpath,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }

  if (1 || qstart_sense_paths == NULL) {
    /* Always try without a qstart_univdiag, in case the given univdiag is bad */
    debug13(printf("Processing without a qstart sense univdiag\n"));
    compute_qstart_local(&qstart_sense_paths,/*depth*/0,path,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_FORWARD,qstart_innerp,find_splices_p,/*salvagep*/false);
  }
  Path_free(&path,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);

  if (splicingp == true) {
    path = Path_new_for_qstart_extension(middle_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
					 /*middle_nmismatches*/-1,plusp,first_read_p,genestrand,
					 /*sensedir*/SENSE_ANTI,querylength,
					 method,chrnum,chroffset,chrhigh,
					 /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
					 query_compress_fwd,query_compress_rev,
					 intlistpool,univcoordlistpool,pathpool);
    for (p = qstart_univdiags; p != NULL; p = List_next(p)) {
      univdiag = (Univdiag_T) List_head(p);
      debug13(printf("Processing qstart univdiag %u %u %d..%d\n",
		     univdiag->univdiagonal,middle_univdiagonal - univdiag->univdiagonal,univdiag->qstart,univdiag->qend));
      /* Previously checked for a large gap between middle_diagonal_qstart and univdiag->qend */

      if (univdiag->univdiagonal < chroffset + querylength) {
	/* Skip, because univdiagonal is in an earlier chromosome */

      } else {
	for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	  if ((newpath = attach_unknown_qstart(path,/*low_univdiagonal*/univdiag->univdiagonal,
						 /*low_qstart*/univdiag->qstart,
						 chroffset,chrhigh,querylength,
						 stage1->indelinfo,stage1->spliceinfo,knownsplicing,
						 stage1,query_compress,queryptr,
						 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						 localdb_tally,localdb_overflowp,max_nblocks,
						 localdb,localdb_nmismatches_allowed,plusp,genestrand,qstart_innerp,
						 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
						 vectorpool,hitlistpool,/*try_sensedir*/SENSE_ANTI,middle_rank)) == NULL) {
	    break;
	  } else {
	    compute_qstart_local(&qstart_antisense_paths,/*depth*/0,newpath,queryptr,querylength,
				 mismatch_positions_alloc,novel_univdiags_alloc,
				 novel_diagonals_alloc,localdb_alloc,
				 localdb_tally,localdb_overflowp,max_nblocks,
				 stage1,knownsplicing,knownindels,
				 query_compress,chrnum,chroffset,chrhigh,intlistpool,
				 univcoordlistpool,listpool,pathpool,transcriptpool,
				 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				 localdb_nmismatches_allowed,plusp,genestrand,
				 /*try_sensedir*/SENSE_ANTI,qstart_innerp,find_splices_p,/*salvagep*/false);
	    Path_free(&newpath,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}
      }
    }

    if (1 || qstart_antisense_paths == NULL) {
      /* Always try without a qstart_univdiag, in case the given univdiag is bad */
      debug13(printf("Processing without a qstart antisense univdiag\n"));
      compute_qstart_local(&qstart_antisense_paths,/*depth*/0,path,queryptr,querylength,
			   mismatch_positions_alloc,novel_univdiags_alloc,
			   novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   stage1,knownsplicing,knownindels,
			   query_compress,chrnum,chroffset,chrhigh,intlistpool,
			   univcoordlistpool,listpool,pathpool,transcriptpool,
			   univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			   localdb_nmismatches_allowed,plusp,genestrand,
			   /*try_sensedir*/SENSE_ANTI,qstart_innerp,find_splices_p,/*salvagep*/false);
    }
    Path_free(&path,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }


  /* Qend */
  path = Path_new_for_qend_extension(middle_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
				     /*middle_nmismatches*/-1,plusp,first_read_p,genestrand,
				     /*sensedir*/SENSE_FORWARD,querylength,
				     method,chrnum,chroffset,chrhigh,
				     /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				     query_compress_fwd,query_compress_rev,
				     intlistpool,univcoordlistpool,pathpool);

  for (p = qend_univdiags; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    debug13(printf("Processing qend univdiag %u %u %d..%d\n",
		 univdiag->univdiagonal,univdiag->univdiagonal - middle_univdiagonal,univdiag->qstart,univdiag->qend));
    /* Previously checked for a large gap between middle_diagonal_qend and univdiag->qstart */

    if (univdiag->univdiagonal >= chrhigh) {
      /* Skip, because univdiagonal is in a later chromosome */
    
    } else {
      for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	if ((newpath = attach_unknown_qend(path,/*high_univdiagonal*/univdiag->univdiagonal,
					     /*high_qend*/univdiag->qend,
					     chroffset,chrhigh,querylength,
					     stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					     stage1,query_compress,queryptr,
					     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					     localdb_tally,localdb_overflowp,max_nblocks,
					     localdb,localdb_nmismatches_allowed,plusp,genestrand,qend_innerp,
					     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					     vectorpool,hitlistpool,/*try_sensedir*/SENSE_FORWARD,middle_rank)) == NULL) {
	  break;
	} else {
	  compute_qend_local(&qend_sense_paths,/*depth*/0,newpath,queryptr,querylength,
			     mismatch_positions_alloc,novel_univdiags_alloc,
			     novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,
			     stage1,knownsplicing,knownindels,
			     query_compress,chrnum,chroffset,chrhigh,intlistpool,
			     univcoordlistpool,listpool,pathpool,transcriptpool,
			     univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			     localdb_nmismatches_allowed,plusp,genestrand,
			     /*try_sensedir*/SENSE_FORWARD,qend_innerp,find_splices_p,/*salvagep*/false,
				     /*segment_cache*/(Path_segment_cache_T) NULL,
				     /*path_qend_local_cache_p*/false);
	  Path_free(&newpath,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
      }
    }
  }
      
  if (1 || qend_sense_paths == NULL) {
    /* Always try without a qend_univdiag, in case the given univdiag is bad */
    debug13(printf("Processing without a qend sense univdiag\n"));
    compute_qend_local(&qend_sense_paths,/*depth*/0,path,queryptr,querylength,
		       mismatch_positions_alloc,novel_univdiags_alloc,
		       novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       stage1,knownsplicing,knownindels,
		       query_compress,chrnum,chroffset,chrhigh,intlistpool,
		       univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       localdb_nmismatches_allowed,plusp,genestrand,
		       /*try_sensedir*/SENSE_FORWARD,qend_innerp,find_splices_p,/*salvagep*/false,
			       /*segment_cache*/(Path_segment_cache_T) NULL,
			       /*path_qend_local_cache_p*/false);
  }
  Path_free(&path,intlistpool,univcoordlistpool,
	    listpool,pathpool,transcriptpool,hitlistpool);

  if (splicingp == true) {
    path = Path_new_for_qend_extension(middle_univdiagonal,middle_diagonal_qstart,middle_diagonal_qend,
				       /*middle_nmismatches*/-1,plusp,first_read_p,genestrand,
				       /*sensedir*/SENSE_ANTI,querylength,
				       method,chrnum,chroffset,chrhigh,
				       /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				       query_compress_fwd,query_compress_rev,
				       intlistpool,univcoordlistpool,pathpool);
    for (p = qend_univdiags; p != NULL; p = List_next(p)) {
      univdiag = (Univdiag_T) List_head(p);
      debug13(printf("Processing qend univdiag %u %u %d..%d\n",
		   univdiag->univdiagonal,univdiag->univdiagonal - middle_univdiagonal,univdiag->qstart,univdiag->qend));
      /* Previously checked for a large gap between middle_diagonal_qend and univdiag->qstart */

      if (univdiag->univdiagonal >= chrhigh) {
	/* Skip, because univdiagonal is in a later chromosome */
    
      } else {
	for (middle_rank = 0; middle_rank < n_splice_paths; middle_rank++) {
	  if ((newpath = attach_unknown_qend(path,/*high_univdiagonal*/univdiag->univdiagonal,
					       /*high_qend*/univdiag->qend,
					       chroffset,chrhigh,querylength,
					       stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					       stage1,query_compress,queryptr,
					       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					       localdb_tally,localdb_overflowp,max_nblocks,
					       localdb,localdb_nmismatches_allowed,plusp,genestrand,qend_innerp,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					       vectorpool,hitlistpool,/*try_sensedir*/SENSE_ANTI,middle_rank)) == NULL) {
	    break;
	  } else {
	    compute_qend_local(&qend_antisense_paths,/*depth*/0,newpath,queryptr,querylength,
			       mismatch_positions_alloc,novel_univdiags_alloc,
			       novel_diagonals_alloc,localdb_alloc,
			       localdb_tally,localdb_overflowp,max_nblocks,
			       stage1,knownsplicing,knownindels,
			       query_compress,chrnum,chroffset,chrhigh,intlistpool,
			       univcoordlistpool,listpool,pathpool,transcriptpool,
			       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			       localdb_nmismatches_allowed,plusp,genestrand,
			       /*try_sensedir*/SENSE_ANTI,qend_innerp,find_splices_p,/*salvagep*/false,
			       /*segment_cache*/(Path_segment_cache_T) NULL,
			       /*path_qend_local_cache_p*/false);
	    Path_free(&newpath,intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
	  }
	}
      }
    }

    if (1 || qend_antisense_paths == NULL) {
      /* Always try without a qend_univdiag, in case the univdiag is bad */
      debug13(printf("Processing without a qend antisense univdiag\n"));
      compute_qend_local(&qend_antisense_paths,/*depth*/0,path,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_ANTI,qend_innerp,find_splices_p,/*salvagep*/false,
				       /*segment_cache*/(Path_segment_cache_T) NULL,
				       /*path_qend_local_cache_p*/false);
    }
    Path_free(&path,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }

  /* Qstart sense.  Also handles the non-splicing case. */
#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  for (p = qstart_sense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    printf("** Complete qstart path, sense %p: ",path);
    Path_print(path);
  }

  if (splicingp == true) {
    for (p = qstart_antisense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      printf("** Complete qstart path, antisense %p: ",path);
      Path_print(path);
    }
  }

  for (p = qend_sense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    printf("** Complete qend path, sense %p: ",path);
    Path_print(path);
  }

  if (splicingp == true) {
    for (p = qend_antisense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      printf("** Complete qend path, antisense %p: ",path);
      Path_print(path);
      printf("\n");
    }
  }
#endif

  if (qstart_sense_paths != NULL && qend_sense_paths != NULL) {
    combined_sense_paths = combine_leftright_paths(&(*found_score),&(*end_score),/*main_univdiagonal*/middle_univdiagonal,
						   qstart_sense_paths,qend_sense_paths,
						   query_compress,query_compress_fwd,query_compress_rev,
						   queryptr,querylength,plusp,first_read_p,genestrand,
						   /*sensedir*/SENSE_FORWARD,chrnum,chroffset,chrhigh,
						   mismatch_positions_alloc,stage1,knownsplicing,spliceendsgen,
						   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						   localdb_tally,localdb_overflowp,max_nblocks,
						   intlistpool,univcoordlistpool,listpool,pathpool,
						   univdiagpool,vectorpool,hitlistpool,method,find_splices_p);

    Path_gc(&qstart_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    Path_gc(&qend_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

  } else if (qstart_sense_paths != NULL) {
    combined_sense_paths = qstart_sense_paths;

  } else if (qend_sense_paths != NULL) {
    combined_sense_paths = qend_sense_paths;
    for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
      Path_reverse((Path_T) List_head(p),/*expect_fwd_p*/true);
    }

  } else {
    combined_sense_paths = (List_T) NULL;
  }


  if (qstart_antisense_paths != NULL && qend_antisense_paths != NULL) {
    combined_antisense_paths = combine_leftright_paths(&(*found_score),&(*end_score),/*main_univdiagonal*/middle_univdiagonal,
						       qstart_antisense_paths,qend_antisense_paths,
						       query_compress,query_compress_fwd,query_compress_rev,
						       queryptr,querylength,plusp,first_read_p,genestrand,
						       /*sensedir*/SENSE_ANTI,chrnum,chroffset,chrhigh,
						       mismatch_positions_alloc,stage1,knownsplicing,spliceendsgen,
						       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						       localdb_tally,localdb_overflowp,max_nblocks,
						       intlistpool,univcoordlistpool,listpool,pathpool,
						       univdiagpool,vectorpool,hitlistpool,method,find_splices_p);
    Path_gc(&qstart_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    Path_gc(&qend_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    
  } else if (qstart_antisense_paths != NULL) {
    combined_antisense_paths = qstart_antisense_paths;

  } else if (qend_antisense_paths != NULL) {
    combined_antisense_paths = qend_antisense_paths;
    for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
      Path_reverse((Path_T) List_head(p),/*expect_fwd_p*/true);
    }
    
  } else {
    combined_antisense_paths = (List_T) NULL;
  }

  debug13(printf("After combining qstart and qend, combined_sense_paths: %d, combined_antisense_paths %d\n",
		 List_length(combined_sense_paths),List_length(combined_antisense_paths)));


  /* Want to treat sense and antisense equally in terms of complete
     versus unextended paths.  So we set endtrim_allowed to be 0 and
     set allow_ambig_p to be true.  We perform the same action on all
     paths found */
  /* Previously had allow_ambig_p to be true, but this falsely puts
     paths into complete, and prevents us from finding fusions from the
     unextended paths */
  
  sense_completep = false;
  for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    /* Cannot trust ambig */
    if (Path_unextendedp(path,/*endtrim_allowed*/0,/*allow_ambig_p*/false) == false) {
      debug13(printf("Sense path is complete: ")); debug13(Path_print(path));
      sense_completep = true;
    }
  }

  antisense_completep = false;
  for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    /* Cannot trust ambig */
    if (Path_unextendedp(path,/*endtrim_allowed*/0,/*allow_ambig_p*/false) == false) {
      debug13(printf("Antisense path is complete: ")); debug13(Path_print(path));
      antisense_completep = true;
    }
  }
  debug13(printf("sense_completep %d, antisense_completep %d\n",sense_completep,antisense_completep));

    
  if (sense_completep == false) {
    /* Put into unextended in case the complete antisense path over-reaches the other end of a paired-end read */
    *unextended_sense_paths = List_append(combined_sense_paths,*unextended_sense_paths);
  } else {
#ifdef REMAP_PATHS
    for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (transcriptome == NULL) {
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
#if 0
      } else if (path->transcriptome_method_p == true) {
	/* Not possible */
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
#endif

      } else if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					       path,query_compress_fwd,query_compress_rev,queryseq,
					       intlistpool,uintlistpool,univcoordlistpool,listpool,
					       transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	*sense_paths = List_append(repaired_paths,*sense_paths);
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);

      } else {
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    Hitlistpool_free_list(&combined_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by Hitlist_push */
#else
    *sense_paths = List_append(combined_sense_paths,*sense_paths);
#endif
  }


  if (antisense_completep == false) {
    /* Put into unextended in case the complete sense path over-reaches the other end of a paired-end read */
    *unextended_antisense_paths = List_append(combined_antisense_paths,*unextended_antisense_paths);

  } else {
#ifdef REMAP_PATHS
    for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (transcriptome == NULL) {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));

#if 0
      } else if (path->transcriptome_method_p(path) == true) {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));
#endif

      } else if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					       path,query_compress_fwd,query_compress_rev,queryseq,
					       intlistpool,uintlistpool,univcoordlistpool,listpool,
					       transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	*antisense_paths = List_append(repaired_paths,*antisense_paths);
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);

      } else {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    Hitlistpool_free_list(&combined_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by Hitlist_push */
#else
    *antisense_paths = List_append(combined_antisense_paths,*antisense_paths);
#endif
  }


#if defined(DEBUG13) || defined(DEBUG_PATH_CORE)
  printf("Result: first_read_p %d, plusp %d\n",first_read_p,plusp);
  printf("Result: %d complete sense paths\n",List_length(*sense_paths));
  printf("Result: %d complete antisense paths\n",List_length(*antisense_paths));
  printf("Result: %d unextended sense paths\n",List_length(*unextended_sense_paths));
  printf("Result: %d unextended antisense paths\n",List_length(*unextended_antisense_paths));
#endif

#if 0  
  /* Already freed or used above */
  Path_gc(&qstart_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qend_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qstart_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qend_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#endif

  return;
}
#endif


/* Used for KMER_EXACT1 */
/* Calls Genome_count_mismatches_substring, which sets nmismatches correctly */
bool
Path_solve_exact (int *found_score, int *end_score,
		  List_T *sense_paths, List_T *antisense_paths,
		  Univcoord_T univdiagonal, Auxinfo_T auxinfo, int querylength,
		  int max_nmismatches,
		  bool plusp, bool first_read_p, int genestrand,
		  Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Shortread_T queryseq, char *queryuc_ptr, char *queryrc,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool, Method_T method) {

  bool perfect_ends_p = true;

#ifdef REMAP_PATHS
  List_T repaired_paths;
#endif

  T sense_path, antisense_path;
  int nmismatches, ref_nmismatches;

  /* auxinfo stores kmer_querystart and kmer_queryend in qstart and qend */
  int kmer_querystart = auxinfo->qstart;
  int kmer_queryend = auxinfo->qend;

  /* Determine nmismatches */
  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
						      univdiagonal,querylength,/*pos5*/0,/*pos3*/querylength,
						      plusp,genestrand);
  /* debug13(printf("Counting mismatches from %d to %d => %d (%d ref)\n",
     pos5,pos3,nmismatches,ref_nmismatches)); */

  if (max_nmismatches >= 0 && nmismatches > max_nmismatches) {
    return false;
  }
      
#if 0
  if (nmismatches > nmismatches_allowed) {
    debug1(printf("Path_solve_exact: Not considering because nmismatches %d > nmismatches_allowed %d\n",
		  nmismatches,nmismatches_allowed));
    return true;		/* Caller is looking for cases of imperfect ends */

  } else {
#endif
    /* Since pos5 is 0 and pos3 is querylength, found_score is nmismatches */
    sense_path = Path_new_exact(univdiagonal,/*pos5*/0,/*pos3*/querylength,nmismatches,ref_nmismatches,
				plusp,first_read_p,genestrand,querylength,/*found_score*/nmismatches,
				chrnum,chroffset,chrhigh,query_compress_fwd,query_compress_rev,
				intlistpool,univcoordlistpool,pathpool,method);
    if (/*found_score*/nmismatches < *found_score) {
      *found_score = nmismatches;
    }
    *end_score = 0;

#ifdef CHECK_ASSERTIONS
    Path_eval_nmatches(&(*found_score),&(*end_score),sense_path,query_compress_fwd,query_compress_rev);
#endif
    perfect_ends_p = Path_eval_perfect_ends_p(sense_path,query_compress_fwd,queryuc_ptr,
					      query_compress_rev,queryrc,
					      kmer_querystart,kmer_queryend,pathpool);

    if (splicingp == false) {
      sense_path->completep = true;
      *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) sense_path
				  hitlistpool_trace(__FILE__,__LINE__));

    } else {
      antisense_path = Path_copy(sense_path,intlistpool,univcoordlistpool,listpool,
				 pathpool,vectorpool,transcriptpool,hitlistpool);
      antisense_path->sensedir = SENSE_ANTI;

#ifdef REMAP_PATHS
      if (transcriptome == NULL) {
	sense_path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) sense_path
				    hitlistpool_trace(__FILE__,__LINE__));
	  
	antisense_path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) antisense_path
					hitlistpool_trace(__FILE__,__LINE__));

      } else {
	if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					  sense_path,query_compress_fwd,query_compress_rev,queryseq,
					  intlistpool,uintlistpool,univcoordlistpool,listpool,
					  transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	  *sense_paths = List_append(repaired_paths,*sense_paths);
	  Path_free(&sense_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else {
	  sense_path->completep = true;
	  *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) sense_path
				      hitlistpool_trace(__FILE__,__LINE__));
	}

	if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					  antisense_path,query_compress_fwd,query_compress_rev,queryseq,
					  intlistpool,uintlistpool,univcoordlistpool,listpool,
					  transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	  *antisense_paths = List_append(repaired_paths,*antisense_paths);
	  Path_free(&antisense_path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);

	} else {
	  antisense_path->completep = true;
	  *antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) antisense_path
					  hitlistpool_trace(__FILE__,__LINE__));
	}
      }
#else
      sense_path->completep = true;
      *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) sense_path
				  hitlistpool_trace(__FILE__,__LINE__));
	  
      antisense_path->completep = true;
      *antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) antisense_path
				      hitlistpool_trace(__FILE__,__LINE__));
#endif
    }

    return perfect_ends_p;
#if 0
  }
#endif
}


int
Path_solve_count_mismatches_substring (int *ref_nmismatches, Compress_T query_compress,
				       Univcoord_T univdiagonal, int querylength,
				       int pos5, int pos3, bool plusp, int genestrand) {
  return Genomebits_count_mismatches_substring(&(*ref_nmismatches),genomebits,genomebits_alt,
					       query_compress,univdiagonal,querylength,
					       pos5,pos3,plusp,genestrand);
}


#if 0
/* kmer_size should generally be index1part, but should be 4 for localdb matches */
/* Calls either Path_solve_exact or Path_solve_from_diagonals, which both set nmismatches correctly */
void
Path_solve_from_univdiagonal (int *found_score, int *end_score,

			      List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			      List_T *sense_paths, List_T *antisense_paths,
			      
			      Univcoord_T univdiagonal, Auxinfo_T auxinfo,
			      Shortread_T queryseq, char *queryptr, Compress_T query_compress,
			      Compress_T query_compress_fwd, Compress_T query_compress_rev,
			      bool plusp, int querylength, int *mismatch_positions_alloc,
			      
			      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			      int localdb_nmismatches_allowed,

			      Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			      bool paired_end_p, bool first_read_p,

			      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			      Univcoordlistpool_T univcoordlistpool,
			      Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			      Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p) {

  T sense_path, antisense_path;

  int score, nmismatches, ref_nmismatches;
  int qstart, qend;

#ifdef INDIVIDUAL_CHRINFO
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
#endif

  qstart = auxinfo->qstart;
  qend = auxinfo->qend;

  debug6(printf("Entered Path_solve_from_univdiagonal with univdiagonal %u => %d..%d\n",
		univdiagonal,qstart,qend));
  assert(qstart < qend);

#ifdef INDIVIDUAL_CHRINFO
  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,univdiagonal - querylength,univdiagonal);
#endif

#if 0
  if (0 && qstart <= 8 && querylength - qend <= 8) {
    /* This prevents us from extending ends of the local mate */
    nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							univdiagonal,querylength,
							/*pos5*/qstart,/*pos3*/qend,
							plusp,/*genestrand*/0);
    score = qstart + (querylength - qend) + nmismatches;
    if (score < *found_score) {
      *found_score = score;
    }

    sense_path = Path_new_exact(univdiagonal,qstart,qend,nmismatches,/*ref_nmismatches*/nmismatches,
				plusp,first_read_p,/*genestrand*/0,querylength,/*found_score*/score,
				chrnum,chroffset,chrhigh,query_compress_fwd,query_compress_rev,
				intlistpool,univcoordlistpool,pathpool,method);
    sense_path->completep = true;
    *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) sense_path
				hitlistpool_trace(__FILE__,__LINE__));

    if (splicingp == true) {
      antisense_path = Path_new_exact(univdiagonal,qstart,qend,nmismatches,/*ref_nmismatches*/nmismatches,
				      plusp,first_read_p,/*genestrand*/0,querylength,/*found_score*/score,
				      chrnum,chroffset,chrhigh,query_compress_fwd,query_compress_rev,
				      intlistpool,univcoordlistpool,pathpool,method);
      antisense_path->sensedir = SENSE_ANTI;
      antisense_path->completep = true;
      *antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) antisense_path
				      hitlistpool_trace(__FILE__,__LINE__));
    }

  } else if (find_splices_p == false) {
    /* This also prevents us from extending ends of the local mate */
    sense_path = Path_new_exact(univdiagonal,qstart,qend,/*nmismatches*/-1,/*ref_nmismatches*/-1,
				plusp,first_read_p,/*genestrand*/0,querylength,/*found_score*/querylength,
				chrnum,chroffset,chrhigh,query_compress_fwd,query_compress_rev,
				intlistpool,univcoordlistpool,pathpool,method);
    *unextended_sense_paths = Hitlist_push(*unextended_sense_paths,hitlistpool,(void *) sense_path
					   hitlistpool_trace(__FILE__,__LINE__));

    if (splicingp == true) {
      antisense_path = Path_new_exact(univdiagonal,qstart,qend,/*nmismatches*/-1,/*ref_nmismatches*/-1,
				      plusp,first_read_p,/*genestrand*/0,querylength,/*found_score*/querylength,
				      chrnum,chroffset,chrhigh,query_compress_fwd,query_compress_rev,
				      intlistpool,univcoordlistpool,pathpool,method);
      antisense_path->sensedir = SENSE_ANTI;
      *unextended_antisense_paths = Hitlist_push(*unextended_antisense_paths,hitlistpool,(void *) antisense_path
						 hitlistpool_trace(__FILE__,__LINE__));
    }

  } else {
#endif

    /* Original auxinfo for KMER_PREVALENT used 0,querylength */
    auxinfo->qstart = qstart;
    auxinfo->qend = qend;
    auxinfo->nmismatches = -1;

    Path_solve_from_diagonals(&(*found_score),&(*end_score),

			      &(*unextended_sense_paths),&(*unextended_antisense_paths),
			      &(*sense_paths),&(*antisense_paths),

			      /*middle_univdiagonal*/univdiagonal,auxinfo,
			      
			      queryseq,queryptr,querylength,mismatch_positions_alloc,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,
			      localdb_nmismatches_allowed,

			      stage1,knownsplicing,knownindels,
			      query_compress,query_compress_fwd,query_compress_rev,
			      chrnum,chroffset,chrhigh,plusp,/*genestrand*/0,
			      paired_end_p,first_read_p,
			      intlistpool,uintlistpool,univcoordlistpool,
			      listpool,pathpool,transcriptpool,
			      univdiagpool,vectorpool,hitlistpool,spliceendsgen,method,
			      /*find_splices_p*/false);
#if 0
  }
#endif

  return;
}
#endif



/* Used for EXT (extension search) and SEGMENT1 (Kmer_segment) */
static void
solve_chain_qstart (List_T *qstart_sense_paths, List_T *qstart_antisense_paths,

		    List_T left_chain, char *queryptr, int querylength, int *mismatch_positions_alloc,

		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		    int localdb_nmismatches_allowed,

		    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
		    Path_segment_cache_T segment_cache, bool segment_prefix_skip_p,
		    bool skip_initial_local_p) {

  T path_sense, path_antisense, newpath;

  Univdiag_T middle_univdiag, univdiag;
  Univcoord_T middle_univdiagonal;

  List_T p;
  bool qstart_innerp, qend_innerp, first_univdiag_p;


  debug12(printf("Entered solve_chain_qstart from method %s, find_splices_p %d, first_read_p %d\n",
		 Method_string(method),find_splices_p,first_read_p));

#ifdef DEBUG12
  for (p = left_chain; p != NULL; p = List_next(p)) {
    printf("qstart diagonal %u %d %d..%d\n",
	   ((Univdiag_T) List_head(p))->univdiagonal,
	   middle_univdiagonal - ((Univdiag_T) List_head(p))->univdiagonal,
	   ((Univdiag_T) List_head(p))->qstart,
	   ((Univdiag_T) List_head(p))->qend);
  }
#endif

  if (paired_end_p == false) {
    qstart_innerp = qend_innerp = false;
  } else if (first_read_p == plusp) {
    qstart_innerp = false;
    qend_innerp = true;
  } else {
    qstart_innerp = true;
    qend_innerp = false;
  }

  assert(left_chain != NULL);
  middle_univdiag = (Univdiag_T) List_head(left_chain);
  middle_univdiagonal = middle_univdiag->univdiagonal;

  debug12(printf("Qstart: Starting with univdiag %u %d..%d\n",
		 middle_univdiag->univdiagonal,middle_univdiag->qstart,middle_univdiag->qend));

  /* Skip the first element of the chain, which is the middle univdiagonal */
  p = List_next(left_chain);
  if (segment_prefix_skip_p == true && p != NULL && method == SEGMENT1 && paired_end_p == false) {
    univdiag = (Univdiag_T) List_head(p);
    if (univdiag->univdiagonal >= chroffset + querylength &&
	segment_qstart_cache_contains(segment_cache,middle_univdiag,univdiag,
				      querylength,localdb_nmismatches_allowed,
				      chrnum,chroffset,chrhigh,plusp,genestrand,
				      first_read_p,find_splices_p) == true) {
      return;
    }
  }

  path_sense =
    Path_new_for_qstart_extension(middle_univdiagonal,middle_univdiag->qstart,middle_univdiag->qend,
				  /*nmismatches*/-1,plusp,first_read_p,genestrand,
				  /*sensedir*/SENSE_FORWARD,querylength,
				  method,chrnum,chroffset,chrhigh,
				  /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				  query_compress_fwd,query_compress_rev,
				  intlistpool,univcoordlistpool,pathpool);

  if (splicingp == false) {
    path_antisense = (Path_T) NULL;
  } else {
    path_antisense =
      Path_new_for_qstart_extension(middle_univdiagonal,middle_univdiag->qstart,middle_univdiag->qend,
				    /*nmismatches*/-1,plusp,first_read_p,genestrand,
				    /*sensedir*/SENSE_ANTI,querylength,
				    method,chrnum,chroffset,chrhigh,
				    /*splice5p*/false,/*splicetype5*/NO_SPLICE,/*ambig_prob_5*/0.0,
				    query_compress_fwd,query_compress_rev,
				    intlistpool,univcoordlistpool,pathpool);
  }

  first_univdiag_p = true;
  while (p != NULL && (path_sense != (Path_T) NULL || path_antisense != (Path_T) NULL)) {
    /* First, ignore the univdiag */
    if (skip_initial_local_p == false && path_sense != (Path_T) NULL) {
      compute_qstart_local(&(*qstart_sense_paths),/*depth*/0,path_sense,queryptr,querylength,
			   mismatch_positions_alloc,novel_univdiags_alloc,
			   novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   stage1,knownsplicing,knownindels,
			   query_compress,chrnum,chroffset,chrhigh,intlistpool,
			   univcoordlistpool,listpool,pathpool,transcriptpool,
			   univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			   localdb_nmismatches_allowed,plusp,genestrand,
			   /*try_sensedir*/SENSE_FORWARD,qstart_innerp,find_splices_p,/*salvagep*/false);
    }

    if (skip_initial_local_p == false && path_antisense != (Path_T) NULL) {
      compute_qstart_local(&(*qstart_antisense_paths),/*depth*/0,path_antisense,queryptr,querylength,
			   mismatch_positions_alloc,novel_univdiags_alloc,
			   novel_diagonals_alloc,localdb_alloc,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   stage1,knownsplicing,knownindels,
			   query_compress,chrnum,chroffset,chrhigh,intlistpool,
			   univcoordlistpool,listpool,pathpool,transcriptpool,
			   univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			   localdb_nmismatches_allowed,plusp,genestrand,
			   /*try_sensedir*/SENSE_ANTI,qstart_innerp,find_splices_p,/*salvagep*/false);
    }


    /* Second, process the univdiag */
    univdiag = (Univdiag_T) List_head(p);
    debug12(printf("Processing univdiag %u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend));
    
    if (univdiag->univdiagonal < chroffset + querylength) {
      /* Skip, because univdiagonal is in an earlier chromosome */

    } else {
      if (path_sense != (Path_T) NULL) {
	newpath = attach_unknown_qstart(path_sense,/*low_univdiagonal*/univdiag->univdiagonal,
					/*low_qstart*/univdiag->qstart,chroffset,chrhigh,querylength,
					stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					stage1,query_compress,queryptr,
					novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					localdb_tally,localdb_overflowp,max_nblocks,
					localdb,localdb_nmismatches_allowed,plusp,genestrand,qstart_innerp,
					intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					vectorpool,hitlistpool,/*try_sensedir*/SENSE_FORWARD,/*middle_rank*/0);
	Path_free(&path_sense,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	path_sense = newpath;
      }
	
      if (path_antisense != (Path_T) NULL) {
	newpath = attach_unknown_qstart(path_antisense,/*low_univdiagonal*/univdiag->univdiagonal,
					/*low_qstart*/univdiag->qstart,chroffset,chrhigh,querylength,
					stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					stage1,query_compress,queryptr,
					novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					localdb_tally,localdb_overflowp,max_nblocks,
					localdb,localdb_nmismatches_allowed,plusp,genestrand,qstart_innerp,
					intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
					vectorpool,hitlistpool,/*try_sensedir*/SENSE_ANTI,/*middle_rank*/0);
	Path_free(&path_antisense,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	path_antisense = newpath;
      }
    }

    if (segment_prefix_skip_p == true &&
	first_univdiag_p == true && method == SEGMENT1 && paired_end_p == false &&
	path_sense == (Path_T) NULL && path_antisense == (Path_T) NULL &&
	univdiag->univdiagonal >= chroffset + querylength) {
      segment_qstart_cache_insert(segment_cache,middle_univdiag,univdiag,
				  querylength,localdb_nmismatches_allowed,
				  chrnum,chroffset,chrhigh,plusp,genestrand,
				  first_read_p,find_splices_p);
    }

    first_univdiag_p = false;
    skip_initial_local_p = false;
    p = List_next(p);
  }

  /* Finish the chains */
  if (skip_initial_local_p == false && path_sense != (Path_T) NULL) {
    compute_qstart_local(&(*qstart_sense_paths),/*depth*/0,path_sense,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_FORWARD,qstart_innerp,find_splices_p,/*salvagep*/false);
    Path_free(&path_sense,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }

  if (skip_initial_local_p == false && path_antisense != (Path_T) NULL) {
    compute_qstart_local(&(*qstart_antisense_paths),/*depth*/0,path_antisense,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_ANTI,qstart_innerp,find_splices_p,/*salvagep*/false);
    Path_free(&path_antisense,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }

  return;
}


static void
solve_chain_qend (List_T *qend_sense_paths, List_T *qend_antisense_paths,

		  List_T right_chain, char *queryptr, int querylength, int *mismatch_positions_alloc,

		  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		  int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		  int localdb_nmismatches_allowed,

		  Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		  Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		  bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		  Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
		  Path_segment_cache_T segment_cache, bool segment_prefix_skip_p,
		  bool path_qend_local_cache_p, bool skip_initial_local_p) {

  T path_sense, path_antisense, newpath;

  Univdiag_T middle_univdiag, univdiag;
  Univcoord_T middle_univdiagonal;

  List_T p;
  bool qstart_innerp, qend_innerp, first_univdiag_p;



  debug12(printf("Entered solve_chain_qend from method %s, find_splices_p %d, first_read_p %d\n",
		 Method_string(method),find_splices_p,first_read_p));

#ifdef DEBUG12
  for (p = right_chain; p != NULL; p = List_next(p)) {
    printf("qend diagonal %u %d %d..%d\n",
	   ((Univdiag_T) List_head(p))->univdiagonal,
	   ((Univdiag_T) List_head(p))->univdiagonal - middle_univdiagonal,
	   ((Univdiag_T) List_head(p))->qstart,
	   ((Univdiag_T) List_head(p))->qend);
  }
#endif

  if (paired_end_p == false) {
    qstart_innerp = qend_innerp = false;
  } else if (first_read_p == plusp) {
    qstart_innerp = false;
    qend_innerp = true;
  } else {
    qstart_innerp = true;
    qend_innerp = false;
  }

  assert(right_chain != NULL);
  middle_univdiag = (Univdiag_T) List_head(right_chain);
  middle_univdiagonal = middle_univdiag->univdiagonal;
  debug12(printf("Qend: Starting with univdiag %u %d..%d\n",
		 middle_univdiag->univdiagonal,middle_univdiag->qstart,middle_univdiag->qend));

  /* Skip the first element of the chain, which is the middle univdiagonal */
  p = List_next(right_chain);
  if (segment_prefix_skip_p == true && p != NULL && method == SEGMENT1 && paired_end_p == false) {
    univdiag = (Univdiag_T) List_head(p);
    if (univdiag->univdiagonal < chrhigh &&
	segment_qend_cache_contains(segment_cache,middle_univdiag,univdiag,
				    querylength,localdb_nmismatches_allowed,
				    chrnum,chroffset,chrhigh,plusp,genestrand,
				    first_read_p,find_splices_p) == true) {
      return;
    }
  }

  path_sense =
    Path_new_for_qend_extension(middle_univdiagonal,middle_univdiag->qstart,middle_univdiag->qend,
				/*nmismatches*/-1,plusp,first_read_p,genestrand,
				/*sensedir*/SENSE_FORWARD,querylength,
				method,chrnum,chroffset,chrhigh,
				/*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				query_compress_fwd,query_compress_rev,
				intlistpool,univcoordlistpool,pathpool);

  if (splicingp == false) {
    path_antisense = (Path_T) NULL;
  } else {
    path_antisense =
      Path_new_for_qend_extension(middle_univdiagonal,middle_univdiag->qstart,middle_univdiag->qend,
				  /*nmismatches*/-1,plusp,first_read_p,genestrand,
				  /*sensedir*/SENSE_ANTI,querylength,
				  method,chrnum,chroffset,chrhigh,
				  /*splice3p*/false,/*splicetype3*/NO_SPLICE,/*ambig_prob_3*/0.0,
				  query_compress_fwd,query_compress_rev,
				  intlistpool,univcoordlistpool,pathpool);
  }


  first_univdiag_p = true;
  while (p != NULL && (path_sense != (Path_T) NULL || path_antisense != (Path_T) NULL)) {
    /* First, ignore the univdiag */
    if (skip_initial_local_p == false && path_sense != (Path_T) NULL) {
      compute_qend_local(&(*qend_sense_paths),/*depth*/0,path_sense,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_FORWARD,qend_innerp,find_splices_p,/*salvagep*/false,
			 segment_cache,path_qend_local_cache_p);
    }

    if (skip_initial_local_p == false && path_antisense != (Path_T) NULL) {
      compute_qend_local(&(*qend_antisense_paths),/*depth*/0,path_antisense,queryptr,querylength,
			 mismatch_positions_alloc,novel_univdiags_alloc,
			 novel_diagonals_alloc,localdb_alloc,
			 localdb_tally,localdb_overflowp,max_nblocks,
			 stage1,knownsplicing,knownindels,
			 query_compress,chrnum,chroffset,chrhigh,intlistpool,
			 univcoordlistpool,listpool,pathpool,transcriptpool,
			 univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			 localdb_nmismatches_allowed,plusp,genestrand,
			 /*try_sensedir*/SENSE_ANTI,qend_innerp,find_splices_p,/*salvagep*/false,
			 segment_cache,path_qend_local_cache_p);
    }

    /* Second, process the univdiag */
    univdiag = (Univdiag_T) List_head(p);
    debug12(printf("Processing univdiag %u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend));

    if (univdiag->univdiagonal >= chrhigh) {
      /* Skip, because univdiagonal is in a later chromosome */
    
    } else {
      if (path_sense != (Path_T) NULL) {
	newpath = attach_unknown_qend(path_sense,/*high_univdiagonal*/univdiag->univdiagonal,
				      /*high_qend*/univdiag->qend,chroffset,chrhigh,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,
				      stage1,query_compress,queryptr,
				      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				      localdb_tally,localdb_overflowp,max_nblocks,
				      localdb,localdb_nmismatches_allowed,plusp,genestrand,qend_innerp,
				      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				      vectorpool,hitlistpool,/*try_sensedir*/SENSE_FORWARD,/*middle_rank*/0);
	Path_free(&path_sense,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	path_sense = newpath;
      }

      if (path_antisense != (Path_T) NULL) {
	newpath = attach_unknown_qend(path_antisense,/*high_univdiagonal*/univdiag->univdiagonal,
				      /*high_qend*/univdiag->qend,chroffset,chrhigh,querylength,
				      stage1->indelinfo,stage1->spliceinfo,knownsplicing,
				      stage1,query_compress,queryptr,
				      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				      localdb_tally,localdb_overflowp,max_nblocks,
				      localdb,localdb_nmismatches_allowed,plusp,genestrand,qend_innerp,
				      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				      vectorpool,hitlistpool,/*try_sensedir*/SENSE_ANTI,/*middle_rank*/0);
	Path_free(&path_antisense,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	path_antisense = newpath;
      }
    }

    if (segment_prefix_skip_p == true &&
	first_univdiag_p == true && method == SEGMENT1 && paired_end_p == false &&
	path_sense == (Path_T) NULL && path_antisense == (Path_T) NULL &&
	univdiag->univdiagonal < chrhigh) {
      segment_qend_cache_insert(segment_cache,middle_univdiag,univdiag,
				querylength,localdb_nmismatches_allowed,
				chrnum,chroffset,chrhigh,plusp,genestrand,
				first_read_p,find_splices_p);
    }

    first_univdiag_p = false;
    skip_initial_local_p = false;
    p = List_next(p);
  }

  /* Finish the chains */
  if (skip_initial_local_p == false && path_sense != (Path_T) NULL) {
    compute_qend_local(&(*qend_sense_paths),/*depth*/0,path_sense,queryptr,querylength,
		       mismatch_positions_alloc,novel_univdiags_alloc,
		       novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       stage1,knownsplicing,knownindels,
		       query_compress,chrnum,chroffset,chrhigh,intlistpool,
		       univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       localdb_nmismatches_allowed,plusp,genestrand,
		       /*try_sensedir*/SENSE_FORWARD,qend_innerp,find_splices_p,/*salvagep*/false,
		       segment_cache,path_qend_local_cache_p);
    Path_free(&path_sense,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }

  if (skip_initial_local_p == false && path_antisense != (Path_T) NULL) {
    compute_qend_local(&(*qend_antisense_paths),/*depth*/0,path_antisense,queryptr,querylength,
		       mismatch_positions_alloc,novel_univdiags_alloc,
		       novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       stage1,knownsplicing,knownindels,
		       query_compress,chrnum,chroffset,chrhigh,intlistpool,
		       univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       localdb_nmismatches_allowed,plusp,genestrand,
		       /*try_sensedir*/SENSE_ANTI,qend_innerp,find_splices_p,/*salvagep*/false,
		       segment_cache,path_qend_local_cache_p);
    Path_free(&path_antisense,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }

  return;
}

static void
solve_combine (int *found_score, int *end_score,

	       List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
	       List_T *sense_paths, List_T *antisense_paths,

	       List_T qstart_sense_paths, List_T qstart_antisense_paths,
	       List_T qend_sense_paths, List_T qend_antisense_paths,
	       bool qstart_borrowed_p, bool qend_borrowed_p,

	       Univcoord_T middle_univdiagonal,
	       Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

	       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	       int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,

	       Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
	       Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
	       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
	       bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
	       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
	       Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
	       Path_segment_cache_T segment_cache, bool emitted_path_dedup_p) {
  
  Path_T path;

#ifdef TRY_REPAIRS
  List_T repaired_paths;
#endif

  List_T combined_sense_paths, combined_antisense_paths, p;
  bool sense_completep, antisense_completep;


#ifdef DEBUG12
  printf("Start: paired_end_p %d, first_read_p %d, plusp %d\n",paired_end_p,first_read_p,plusp);
  printf("Start: %d complete sense paths\n",List_length(*sense_paths));
  printf("Start: %d complete antisense paths\n",List_length(*antisense_paths));
  printf("Start: %d unextended sense paths\n",List_length(*unextended_sense_paths));
  printf("Start: %d unextended antisense paths\n",List_length(*unextended_antisense_paths));
#endif


  /* Code below is taken from Path_solve_from_diagonals */
  if (qstart_sense_paths != NULL && qend_sense_paths != NULL) {
    combined_sense_paths = combine_leftright_paths(&(*found_score),&(*end_score),/*main_univdiagonal*/middle_univdiagonal,
						   qstart_sense_paths,qend_sense_paths,
						   query_compress,query_compress_fwd,query_compress_rev,
						   queryptr,querylength,plusp,first_read_p,genestrand,
						   /*sensedir*/SENSE_FORWARD,chrnum,chroffset,chrhigh,
						   mismatch_positions_alloc,stage1,knownsplicing,spliceendsgen,
						   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						   localdb_tally,localdb_overflowp,max_nblocks,
						   intlistpool,univcoordlistpool,listpool,pathpool,
						   univdiagpool,vectorpool,hitlistpool,method,find_splices_p);

    if (qstart_borrowed_p == false) {
      Path_gc(&qstart_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    }
    if (qend_borrowed_p == false) {
      Path_gc(&qend_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    }

  } else if (qstart_sense_paths != NULL) {
    /* Preserve legacy duplicate multiplicity.  Equal paths participate in a
       tied qsort before maxpaths reporting, so silently dropping a duplicate
       can change which genomic records survive the reporting cap. */
    combined_sense_paths = qstart_borrowed_p == true ?
      segment_cache_copy_paths(qstart_sense_paths,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) :
      qstart_sense_paths;

  } else if (qend_sense_paths != NULL) {
    combined_sense_paths = qend_borrowed_p == true ?
      segment_cache_copy_paths(qend_sense_paths,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) :
      qend_sense_paths;
    for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
      Path_reverse((Path_T) List_head(p),/*expect_fwd_p*/true);
    }

  } else {
    combined_sense_paths = (List_T) NULL;
  }


  if (qstart_antisense_paths != NULL && qend_antisense_paths != NULL) {
    combined_antisense_paths = combine_leftright_paths(&(*found_score),&(*end_score),/*main_univdiagonal*/middle_univdiagonal,
						       qstart_antisense_paths,qend_antisense_paths,
						       query_compress,query_compress_fwd,query_compress_rev,
						       queryptr,querylength,plusp,first_read_p,genestrand,
						       /*sensedir*/SENSE_ANTI,chrnum,chroffset,chrhigh,
						       mismatch_positions_alloc,stage1,knownsplicing,spliceendsgen,
						       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
						       localdb_tally,localdb_overflowp,max_nblocks,
						       intlistpool,univcoordlistpool,listpool,pathpool,
						       univdiagpool,vectorpool,hitlistpool,method,find_splices_p);
    if (qstart_borrowed_p == false) {
      Path_gc(&qstart_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    }
    if (qend_borrowed_p == false) {
      Path_gc(&qend_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
    }
    
  } else if (qstart_antisense_paths != NULL) {
    combined_antisense_paths = qstart_borrowed_p == true ?
      segment_cache_copy_paths(qstart_antisense_paths,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) :
      qstart_antisense_paths;

  } else if (qend_antisense_paths != NULL) {
    combined_antisense_paths = qend_borrowed_p == true ?
      segment_cache_copy_paths(qend_antisense_paths,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) :
      qend_antisense_paths;
    for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
      Path_reverse((Path_T) List_head(p),/*expect_fwd_p*/true);
    }
    
  } else {
    combined_antisense_paths = (List_T) NULL;
  }

  if (emitted_path_dedup_p == true) {
    combined_sense_paths =
      filter_emitted_paths(combined_sense_paths,segment_cache,
			   intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,hitlistpool,transcriptpool);
    combined_antisense_paths =
      filter_emitted_paths(combined_antisense_paths,segment_cache,
			   intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,hitlistpool,transcriptpool);
  }

  debug12(printf("After combining qstart and qend, combined_sense_paths: %d, combined_antisense_paths %d\n",
		 List_length(combined_sense_paths),List_length(combined_antisense_paths)));


  /* Want to treat sense and antisense equally in terms of complete
     versus unextended paths.  So we set endtrim_allowed to be 0 and
     set allow_ambig_p to be true.  We perform the same action on all
     paths found */
  /* Previously had allow_ambig_p to be true, but this falsely puts
     paths into complete, and prevents us from finding fusions from the
     unextended paths */
  
  sense_completep = false;
  for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    /* Cannot trust ambig */
    if (Path_unextendedp(path,/*endtrim_allowed*/0,/*allow_ambig_p*/false) == false) {
      debug12(printf("Sense path is complete: ")); debug12(Path_print(path));
      sense_completep = true;
    }
  }

  antisense_completep = false;
  for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    /* Cannot trust ambig */
    if (Path_unextendedp(path,/*endtrim_allowed*/0,/*allow_ambig_p*/false) == false) {
      debug12(printf("Antisense path is complete: ")); debug12(Path_print(path));
      antisense_completep = true;
    }
  }
  debug12(printf("sense_completep %d, antisense_completep %d\n",sense_completep,antisense_completep));

    
  if (sense_completep == false) {
    /* Put into unextended in case the complete antisense path over-reaches the other end of a paired-end read */
    *unextended_sense_paths = List_append(combined_sense_paths,*unextended_sense_paths);
  } else {
#ifdef REMAP_PATHS
    for (p = combined_sense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (transcriptome == NULL) {
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
#if 0
      } else if (path->transcriptome_method_p == true) {
	/* Not possible */
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
#endif

      } else if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					       path,query_compress_fwd,query_compress_rev,queryseq,
					       intlistpool,uintlistpool,univcoordlistpool,listpool,
					       transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	*sense_paths = List_append(repaired_paths,*sense_paths);
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);

      } else {
	path->completep = true;
	*sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) Path_expect_fwd(path)
				    hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    Hitlistpool_free_list(&combined_sense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by Hitlist_push */
#else
    *sense_paths = List_append(combined_sense_paths,*sense_paths);
#endif
  }


  if (antisense_completep == false) {
    /* Put into unextended in case the complete sense path over-reaches the other end of a paired-end read */
    *unextended_antisense_paths = List_append(combined_antisense_paths,*unextended_antisense_paths);

  } else {
#ifdef REMAP_PATHS
    for (p = combined_antisense_paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (transcriptome == NULL) {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));

#if 0
      } else if (path->transcriptome_method_p(path) == true) {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));
#endif

      } else if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					       path,query_compress_fwd,query_compress_rev,queryseq,
					       intlistpool,uintlistpool,univcoordlistpool,listpool,
					       transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {
	*antisense_paths = List_append(repaired_paths,*antisense_paths);
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);

      } else {
	path->completep = true;
	*antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) Path_expect_fwd(path)
					hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    Hitlistpool_free_list(&combined_antisense_paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__)); /* Allocated by Hitlist_push */
#else
    *antisense_paths = List_append(combined_antisense_paths,*antisense_paths);
#endif
  }


#ifdef DEBUG12
  printf("Result: first_read_p %d, plusp %d\n",first_read_p,plusp);
  printf("Result: %d complete sense paths\n",List_length(*sense_paths));
  printf("Result: %d complete antisense paths\n",List_length(*antisense_paths));
  printf("Result: %d unextended sense paths\n",List_length(*unextended_sense_paths));
  printf("Result: %d unextended antisense paths\n",List_length(*unextended_antisense_paths));
#endif

#if 0  
  /* Already freed or used above */
  Path_gc(&qstart_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qend_sense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qstart_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&qend_antisense_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#endif

  return;
}


/* A chain is a list of Univdiag_T objects */
/* This procedure converts a chain to paths */
/* For single-end, computing as left chains.  For paired-end, compute low end as left chain and high end as right chain */

#ifdef GSNAP_EXT_ANCHOR_PROOF
static bool
profile_ext_anchor_shadow (struct Accel_profile_ext_anchor_key_T *key_out,
			   Univdiag_T head, unsigned char solve_side,
			   char *queryptr, int querylength, int max_nblocks,
			   int localdb_nmismatches_allowed,
			   Stage1_T stage1, Knownsplicing_T knownsplicing,
			   Knownindels_T knownindels,
			   Compress_T query_compress, Compress_T query_compress_fwd,
			   Compress_T query_compress_rev,
			   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, int genestrand, bool paired_end_p,
			   bool first_read_p, Spliceendsgen_T spliceendsgen,
			   Method_T method, bool find_splices_p) {
  struct Accel_profile_ext_anchor_key_T key;

  if (paired_end_p == false || (method != EXT && method != SEGMENT1)) {
    return false;
  }
  memset(&key,0,sizeof(key));
  key.queryptr_identity = (uintptr_t) queryptr;
  key.query_compress_identity = (uintptr_t) query_compress;
  key.query_compress_fwd_identity = (uintptr_t) query_compress_fwd;
  key.query_compress_rev_identity = (uintptr_t) query_compress_rev;
  key.stage1_identity = (uintptr_t) stage1;
  key.knownsplicing_identity = (uintptr_t) knownsplicing;
  key.knownindels_identity = (uintptr_t) knownindels;
  key.spliceendsgen_identity = (uintptr_t) spliceendsgen;
  key.head_univdiagonal = (uint64_t) head->univdiagonal;
  key.chroffset = (uint64_t) chroffset;
  key.chrhigh = (uint64_t) chrhigh;
  key.chrnum = (int) chrnum;
  key.head_qstart = head->qstart;
  key.head_qend = head->qend;
  key.querylength = querylength;
  key.max_nblocks = max_nblocks;
  key.localdb_nmismatches_allowed = localdb_nmismatches_allowed;
  key.genestrand = genestrand;
  key.method = (int) method;
  key.solve_side = solve_side;
  key.plusp = plusp ? 1U : 0U;
  key.paired_end_p = paired_end_p ? 1U : 0U;
  key.first_read_p = first_read_p ? 1U : 0U;
  key.find_splices_p = find_splices_p ? 1U : 0U;
  *key_out = key;
  return Accel_profile_ext_anchor_shadow(&key);
}

/* First ownership-proof subset: paths whose complete owned graph consists of
   the base Path_T plus endpoint/diagonal/mismatch/junction lists.  Complex
   graphs are refused and counted until their equality contracts are audited. */
static unsigned int
ext_anchor_path_complex_mask (Path_T path) {
  unsigned int mask = 0U;
  if (path->genomic_diff != NULL) mask |= 1U;
  if (path->qstart_alts != NULL || path->qend_alts != NULL) mask |= 2U;
  if (path->circular_endpoints != NULL) mask |= 4U;
  if (path->fusion_querystart_junction != NULL || path->fusion_queryend_junction != NULL) {
    mask |= 8U;
    if (path->fusion_alts != NULL) mask |= 32U;
  }
  if (path->transcripts != NULL || path->invalid_transcripts != NULL ||
      path->fusion_transcripts != NULL || path->fusion_invalid_transcripts != NULL) mask |= 16U;
  return mask;
}

static const char *ext_anchor_path_mismatch_reason = NULL;
static Path_T ext_anchor_mismatch_x = NULL, ext_anchor_mismatch_y = NULL;

static bool
ext_anchor_simple_path_equal_p (Path_T x, Path_T y) {
  if (x->nmatches != y->nmatches || x->ref_nmatches != y->ref_nmatches ||
      x->junction_splice_prob != y->junction_splice_prob ||
      x->total_splice_prob != y->total_splice_prob ||
      x->found_score != y->found_score || x->end_score != y->end_score ||
      x->score_within_trims != y->score_within_trims) {
    ext_anchor_path_mismatch_reason = "scores"; return false;
  } else if (x->plusp != y->plusp || x->genestrand != y->genestrand ||
             x->sensedir != y->sensedir || x->querylength != y->querylength ||
             x->chrnum != y->chrnum || x->chroffset != y->chroffset ||
             x->chrhigh != y->chrhigh || x->main_univdiagonal != y->main_univdiagonal) {
    ext_anchor_path_mismatch_reason = "coordinates"; return false;
  } else if (Intlist_equal(x->endpoints,y->endpoints) == false) {
    ext_anchor_path_mismatch_reason = "endpoints"; return false;
  } else if (Univcoordlist_equal(x->univdiagonals,y->univdiagonals) == false) {
    ext_anchor_path_mismatch_reason = "univdiagonals"; return false;
  } else if (Intlist_equal(x->nmismatches,y->nmismatches) == false) {
    ext_anchor_path_mismatch_reason = "nmismatches"; return false;
  } else if (Intlist_equal(x->ref_nmismatches,y->ref_nmismatches) == false) {
    ext_anchor_path_mismatch_reason = "ref_nmismatches"; return false;
  } else if (qend_local_junction_list_equal_p(x->junctions,y->junctions) == false) {
    ext_anchor_path_mismatch_reason = "junctions";
    ext_anchor_mismatch_x = x; ext_anchor_mismatch_y = y; return false;
  } else if (x->splice5p != y->splice5p || x->splice3p != y->splice3p ||
             x->splicetype5 != y->splicetype5 || x->splicetype3 != y->splicetype3 ||
             x->ambig_prob_5 != y->ambig_prob_5 || x->ambig_prob_3 != y->ambig_prob_3) {
    ext_anchor_path_mismatch_reason = "splice_ends"; return false;
  } else if (x->completep != y->completep || x->childp != y->childp ||
             x->extendedp != y->extendedp || x->method != y->method ||
             x->transcriptome_method_p != y->transcriptome_method_p ||
             x->query_compress_fwd != y->query_compress_fwd ||
             x->query_compress_rev != y->query_compress_rev) {
    ext_anchor_path_mismatch_reason = "flags_or_compress"; return false;
  } else if (x->fusion_querystart_junction == NULL && x->fusion_queryend_junction == NULL &&
	     y->fusion_querystart_junction == NULL && y->fusion_queryend_junction == NULL) {
    return true;
  }

  if (qend_local_junction_equal_p(x->fusion_querystart_junction,
				y->fusion_querystart_junction) == true &&
    qend_local_junction_equal_p(x->fusion_queryend_junction,
				y->fusion_queryend_junction) == true &&
    x->fusion_chrnum == y->fusion_chrnum && x->fusion_chroffset == y->fusion_chroffset &&
    x->fusion_chrhigh == y->fusion_chrhigh && x->fusion_plusp == y->fusion_plusp &&
    Intlist_equal(x->fusion_endpoints,y->fusion_endpoints) == true &&
    Univcoordlist_equal(x->fusion_univdiagonals,y->fusion_univdiagonals) == true &&
    Intlist_equal(x->fusion_nmismatches,y->fusion_nmismatches) == true &&
    Intlist_equal(x->fusion_ref_nmismatches,y->fusion_ref_nmismatches) == true &&
    qend_local_junction_list_equal_p(x->fusion_junctions,y->fusion_junctions) == true &&
    x->fusion_splicep == y->fusion_splicep &&
    x->fusion_splicetype == y->fusion_splicetype &&
      x->fusion_ambig_prob == y->fusion_ambig_prob) {
    return true;
  } else {
    ext_anchor_path_mismatch_reason = "fusion"; return false;
  }
}

static bool
ext_anchor_simple_list_equal_p (List_T x, List_T y, int *npaths, int *nnodes,
				bool *supportedp, unsigned int *complex_mask) {
  Path_T xpath, ypath;

  while (x != NULL && y != NULL) {
    xpath = (Path_T) List_head(x);
    ypath = (Path_T) List_head(y);
    (*npaths)++;
    *nnodes += 1 + Intlist_length(xpath->endpoints) +
      Univcoordlist_length(xpath->univdiagonals) +
      Intlist_length(xpath->nmismatches) + Intlist_length(xpath->ref_nmismatches) +
      List_length(xpath->junctions);
    if (xpath->fusion_querystart_junction != NULL || xpath->fusion_queryend_junction != NULL) {
      *nnodes += Intlist_length(xpath->fusion_endpoints) +
	Univcoordlist_length(xpath->fusion_univdiagonals) +
	Intlist_length(xpath->fusion_nmismatches) +
	Intlist_length(xpath->fusion_ref_nmismatches) + List_length(xpath->fusion_junctions);
    }
    *complex_mask |= ext_anchor_path_complex_mask(xpath);
    *complex_mask |= ext_anchor_path_complex_mask(ypath);
    if ((*complex_mask & ~8U) != 0U) {
      *supportedp = false;
    } else if (ext_anchor_simple_path_equal_p(xpath,ypath) == false) {
      return false;
    }
    x = List_next(x);
    y = List_next(y);
  }
  if (x == NULL && y == NULL) {
    return true;
  } else {
    ext_anchor_path_mismatch_reason = "list_length";
    return false;
  }
}

static bool
ext_anchor_copy_check (List_T fresh_sense_paths, List_T fresh_antisense_paths,
		       List_T canonical_sense_paths, List_T canonical_antisense_paths,
		       List_T *sense_copy_out, List_T *antisense_copy_out,
		       unsigned long long position_gap,
		       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		       Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T sense_copy, antisense_copy;
  bool supportedp = true, equalp;
  unsigned int complex_mask = 0U;
  int npaths = 0, nnodes = 0;

  ext_anchor_path_mismatch_reason = NULL;
  ext_anchor_mismatch_x = ext_anchor_mismatch_y = NULL;

  sense_copy = NULL;
  for (List_T p = canonical_sense_paths; p != NULL; p = List_next(p)) {
    sense_copy = Hitlist_push(sense_copy,hitlistpool,
			      (void *) Path_copy_solver_result((Path_T) List_head(p),intlistpool,
					      univcoordlistpool,listpool,pathpool,vectorpool,
					      transcriptpool,hitlistpool)
			      hitlistpool_trace(__FILE__,__LINE__));
  }
  sense_copy = List_reverse(sense_copy);
  antisense_copy = NULL;
  for (List_T p = canonical_antisense_paths; p != NULL; p = List_next(p)) {
    antisense_copy = Hitlist_push(antisense_copy,hitlistpool,
				  (void *) Path_copy_solver_result((Path_T) List_head(p),intlistpool,
						  univcoordlistpool,listpool,pathpool,vectorpool,
						  transcriptpool,hitlistpool)
				  hitlistpool_trace(__FILE__,__LINE__));
  }
  antisense_copy = List_reverse(antisense_copy);
  equalp = ext_anchor_simple_list_equal_p(fresh_sense_paths,sense_copy,&npaths,&nnodes,
					 &supportedp,&complex_mask) &&
    ext_anchor_simple_list_equal_p(fresh_antisense_paths,antisense_copy,&npaths,&nnodes,
					  &supportedp,&complex_mask);
  Accel_profile_ext_anchor_copy_check(npaths,nnodes,supportedp,equalp,complex_mask,position_gap);
  *sense_copy_out = sense_copy;
  *antisense_copy_out = antisense_copy;
  return supportedp == true && equalp == true;
}

static bool
ext_anchor_result_store (Path_ext_anchor_cache_T cache,
			 const struct Accel_profile_ext_anchor_key_T *key,
			 List_T fresh_sense_paths, List_T fresh_antisense_paths,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_ext_anchor_result_entry_T *entry;
  unsigned int side;

  if (cache == NULL || key == NULL) return false;
  side = (unsigned int) key->solve_side;
  assert(side < EXT_ANCHOR_RESULT_CACHE_SIZE);
  entry = &cache->entries[side];
  if (entry->occupiedp == true) {
    /* The cache is scoped to one exact EXT head and one immutable solver
       context.  A second distinct key in the same side is a programming
       error, never a reason to recompute silently. */
    assert(memcmp(&entry->key,key,sizeof(*key)) == 0);
    return false;
  }

  entry->occupiedp = true;
  entry->key = *key;
  entry->sense_paths =
    solver_result_copy_paths(fresh_sense_paths,intlistpool,univcoordlistpool,
			     listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
  entry->antisense_paths =
    solver_result_copy_paths(fresh_antisense_paths,intlistpool,univcoordlistpool,
			     listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
  EXT_ANCHOR_COUNT(ext_anchor_first_seen_count);
  return true;
}

static bool
ext_anchor_result_lookup (List_T *sense_paths, List_T *antisense_paths,
			  Path_ext_anchor_cache_T cache,
			  const struct Accel_profile_ext_anchor_key_T *key,
			  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  struct Path_ext_anchor_result_entry_T *entry;
  unsigned int side;

  if (cache == NULL || key == NULL) return false;
  EXT_ANCHOR_COUNT(ext_anchor_lookup_count);
  side = (unsigned int) key->solve_side;
  assert(side < EXT_ANCHOR_RESULT_CACHE_SIZE);
  entry = &cache->entries[side];
  if (entry->occupiedp == false) return false;
  assert(memcmp(&entry->key,key,sizeof(*key)) == 0);

  /* The combiner has explicit borrowed-input handling and never mutates or
     frees these canonical graphs.  Tail solvers use ext_anchor_result_replay,
     which returns owned copies. */
  *sense_paths = entry->sense_paths;
  *antisense_paths = entry->antisense_paths;
  EXT_ANCHOR_COUNT(ext_anchor_hit_count);
  return true;
}

#ifdef GSNAP_CHAIN_ROOT_REPLAY
static bool
ext_anchor_result_replay (List_T *sense_paths, List_T *antisense_paths,
			  Path_ext_anchor_cache_T cache,
			  const struct Accel_profile_ext_anchor_key_T *key,
			  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			  Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T canonical_sense, canonical_antisense;
  if (ext_anchor_result_lookup(&canonical_sense,&canonical_antisense,cache,key,
			       intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			       hitlistpool,transcriptpool) == false) return false;
  *sense_paths =
    solver_result_copy_paths(canonical_sense,intlistpool,univcoordlistpool,
			     listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
  *antisense_paths =
    solver_result_copy_paths(canonical_antisense,intlistpool,univcoordlistpool,
			     listpool,pathpool,vectorpool,hitlistpool,transcriptpool);
  return true;
}

static void
ext_anchor_result_admit (Path_ext_anchor_cache_T cache,
			 const struct Accel_profile_ext_anchor_key_T *key,
			 List_T sense_paths, List_T antisense_paths,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  (void) ext_anchor_result_store(cache,key,sense_paths,antisense_paths,
				 intlistpool,univcoordlistpool,listpool,pathpool,
				 vectorpool,hitlistpool,transcriptpool);
}
#endif

#endif

/* A left chain has higher univdiagonal at head */
static void
path_solve_from_left_chain (int *found_score, int *end_score,

			    List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			    List_T *sense_paths, List_T *antisense_paths,

			    List_T left_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			    int localdb_nmismatches_allowed,

			    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
			    Path_ext_anchor_cache_T ext_anchor_cache,
			    Path_segment_cache_T segment_cache,
			    bool segment_prefix_skip_p,
			    bool segment_qstart_anchor_cache_p,
			    bool path_qend_local_cache_p,
			    bool emitted_path_dedup_p) {
  
  List_T qstart_sense_paths = NULL, qstart_antisense_paths = NULL,
    qend_sense_paths = NULL, qend_antisense_paths = NULL, p;
  List_T right_chain;
#ifdef GSNAP_CHAIN_ROOT_REPLAY
  List_T root_chain;
#endif
  Univdiag_T right_univdiag;
  Univcoord_T middle_univdiagonal;
  T path;

  bool qstart_innerp, qend_innerp;
  bool ext_anchor_reused_p = false;
#ifdef GSNAP_CHAIN_ROOT_REPLAY
  bool chain_root_materialized_p = false;
#endif
#ifdef GSNAP_EXT_ANCHOR_PROOF
  struct Accel_profile_ext_anchor_key_T ext_shadow_key;
#endif
  int found_sensedir;
  int trimpos, exon_origin, qstart, qend, pos5, pos3;
  int nmatches, nmismatches_to_trimpos;
  bool splice5p, splice3p;
  Splicetype_T splicetype5, splicetype3;
  double ambig_prob_5, ambig_prob_3;

  /* Create right chain and solve => qend paths */
  right_univdiag = (Univdiag_T) List_head(left_chain);
  middle_univdiagonal = right_univdiag->univdiagonal;

  p = List_next(left_chain);
  if (segment_prefix_skip_p == true && p != NULL && method == SEGMENT1 && paired_end_p == false) {
    Univdiag_T next_univdiag = (Univdiag_T) List_head(p);

    if (next_univdiag->univdiagonal >= chroffset + querylength &&
	segment_qstart_cache_contains(segment_cache,right_univdiag,next_univdiag,
				      querylength,localdb_nmismatches_allowed,
				      chrnum,chroffset,chrhigh,plusp,genestrand,
				      first_read_p,find_splices_p) == true) {
      /* A previous chain with the same head+next proved that next cannot
	 attach, so the remaining tail cannot affect qstart or the combined path. */
      return;
    }
  }

  right_chain = Univdiagpool_push(NULL,univdiagpool,right_univdiag->qstart,right_univdiag->qend,
				  right_univdiag->nmismatches,middle_univdiagonal
			  univdiagpool_trace(__FILE__,__LINE__));
#ifdef GSNAP_EXT_ANCHOR_PROOF

  (void) profile_ext_anchor_shadow(&ext_shadow_key,
				   right_univdiag,/*solve_side:qend*/1U,
			    queryptr,querylength,max_nblocks,localdb_nmismatches_allowed,
			    stage1,knownsplicing,knownindels,
			    query_compress,query_compress_fwd,query_compress_rev,
			    chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,
			    first_read_p,spliceendsgen,method,find_splices_p);
  ext_anchor_reused_p = ext_anchor_result_lookup(&qend_sense_paths,&qend_antisense_paths,
					   ext_anchor_cache,&ext_shadow_key,
					   intlistpool,univcoordlistpool,listpool,pathpool,
					   vectorpool,hitlistpool,transcriptpool);
#endif
  if (ext_anchor_reused_p == false) {
    solve_chain_qend(&qend_sense_paths,&qend_antisense_paths,right_chain,
		   queryptr,querylength,mismatch_positions_alloc,
		   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		   localdb_tally,localdb_overflowp,max_nblocks,
		   localdb_nmismatches_allowed,

		   stage1,knownsplicing,knownindels,
		   query_compress,query_compress_fwd,query_compress_rev,
		   chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		   intlistpool,univcoordlistpool,listpool,
		   pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		   method,find_splices_p,segment_cache,segment_prefix_skip_p,
		   path_qend_local_cache_p,/*skip_initial_local_p*/false);
  }
#ifdef GSNAP_EXT_ANCHOR_PROOF
  Accel_profile_ext_anchor_shadow_end();
  if (ext_anchor_reused_p == false) {
    (void) ext_anchor_result_store(ext_anchor_cache,&ext_shadow_key,
				   qend_sense_paths,qend_antisense_paths,
				   intlistpool,univcoordlistpool,listpool,pathpool,
				   vectorpool,hitlistpool,transcriptpool);
  }
#endif
  Univdiagpool_gc(&right_chain,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__));

#ifdef GSNAP_CHAIN_ROOT_REPLAY
  root_chain = Univdiagpool_push(NULL,univdiagpool,right_univdiag->qstart,
				right_univdiag->qend,right_univdiag->nmismatches,
				middle_univdiagonal
				univdiagpool_trace(__FILE__,__LINE__));
  (void) profile_ext_anchor_shadow(&ext_shadow_key,right_univdiag,
				   /*solve_side:qstart*/0U,queryptr,querylength,max_nblocks,
				   localdb_nmismatches_allowed,stage1,knownsplicing,knownindels,
				   query_compress,query_compress_fwd,query_compress_rev,
				   chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,
				   first_read_p,spliceendsgen,method,find_splices_p);
  if (ext_anchor_result_replay(&qstart_sense_paths,&qstart_antisense_paths,
			       ext_anchor_cache,&ext_shadow_key,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) == false) {
    solve_chain_qstart(&qstart_sense_paths,&qstart_antisense_paths,root_chain,
		       queryptr,querylength,mismatch_positions_alloc,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
		       stage1,knownsplicing,knownindels,
		       query_compress,query_compress_fwd,query_compress_rev,
		       chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,method,find_splices_p,
		       segment_cache,segment_prefix_skip_p,/*skip_initial_local_p*/false);
    ext_anchor_result_admit(ext_anchor_cache,&ext_shadow_key,
			    qstart_sense_paths,qstart_antisense_paths,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			    hitlistpool,transcriptpool);
  }
  Univdiagpool_gc(&root_chain,univdiagpool univdiagpool_trace(__FILE__,__LINE__));
  chain_root_materialized_p = true;
#endif

  /* Solve given left chain => qstart paths */
  solve_chain_qstart(&qstart_sense_paths,&qstart_antisense_paths,left_chain,
		     queryptr,querylength,mismatch_positions_alloc,
		     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		     localdb_tally,localdb_overflowp,max_nblocks,
		     localdb_nmismatches_allowed,

		     stage1,knownsplicing,knownindels,
		     query_compress,query_compress_fwd,query_compress_rev,
		     chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		     intlistpool,univcoordlistpool,listpool,
		     pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		     method,find_splices_p,segment_cache,segment_prefix_skip_p,
		     /*skip_initial_local_p*/
#ifdef GSNAP_CHAIN_ROOT_REPLAY
		     chain_root_materialized_p
#else
		     false
#endif
		     );

  /* Combine qstart and qend paths */
  solve_combine(&(*found_score),&(*end_score),

		&(*unextended_sense_paths),&(*unextended_antisense_paths),
		&(*sense_paths),&(*antisense_paths),

		qstart_sense_paths,qstart_antisense_paths,
		qend_sense_paths,qend_antisense_paths,
		/*qstart_borrowed_p*/false,/*qend_borrowed_p*/ext_anchor_reused_p,

		middle_univdiagonal,queryseq,queryptr,querylength,mismatch_positions_alloc,
		novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		localdb_tally,localdb_overflowp,max_nblocks,

		stage1,knownsplicing,knownindels,
		query_compress,query_compress_fwd,query_compress_rev,
		chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		intlistpool,uintlistpool,univcoordlistpool,listpool,
		pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		method,find_splices_p,segment_cache,emitted_path_dedup_p);
  return;
}

void
Path_solve_from_left_chain (int *found_score, int *end_score,

			    List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			    List_T *sense_paths, List_T *antisense_paths,

			    List_T left_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			    int localdb_nmismatches_allowed,

			    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p) {

  path_solve_from_left_chain(&(*found_score),&(*end_score),
			     &(*unextended_sense_paths),&(*unextended_antisense_paths),
			     &(*sense_paths),&(*antisense_paths),
			     left_chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,
			     localdb_nmismatches_allowed,
			     stage1,knownsplicing,knownindels,
			     query_compress,query_compress_fwd,query_compress_rev,
			     chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
				     intlistpool,uintlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen,method,find_splices_p,
				     /*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,
				     /*segment_cache*/(Path_segment_cache_T) NULL,
				     /*segment_prefix_skip_p*/false,
				     /*segment_qstart_anchor_cache_p*/false,
				     /*path_qend_local_cache_p*/false,
				     /*emitted_path_dedup_p*/false);
  return;
}

void
Path_solve_from_left_chain_cached (int *found_score, int *end_score,

				   List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
				   List_T *sense_paths, List_T *antisense_paths,

				   List_T left_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

				   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
				   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
				   int localdb_nmismatches_allowed,

				   Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
				   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				   Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				   bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
				   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
					   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
					   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
				   Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
				   Path_ext_anchor_cache_T ext_anchor_cache,
				   Path_segment_cache_T segment_cache,
					   bool segment_prefix_skip_p,
					   bool segment_qstart_anchor_cache_p,
					   bool path_qend_local_cache_p,
					   bool emitted_path_dedup_p) {

  path_solve_from_left_chain(&(*found_score),&(*end_score),
			     &(*unextended_sense_paths),&(*unextended_antisense_paths),
			     &(*sense_paths),&(*antisense_paths),
			     left_chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     localdb_tally,localdb_overflowp,max_nblocks,
			     localdb_nmismatches_allowed,
			     stage1,knownsplicing,knownindels,
			     query_compress,query_compress_fwd,query_compress_rev,
			     chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
				     intlistpool,uintlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				     spliceendsgen,method,find_splices_p,ext_anchor_cache,segment_cache,
				     segment_prefix_skip_p,segment_qstart_anchor_cache_p,
				     path_qend_local_cache_p,emitted_path_dedup_p);
  return;
}


/* A right chain has lower univdiagonal at head */
static void
path_solve_from_right_chain (int *found_score, int *end_score,

			     List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			     List_T *sense_paths, List_T *antisense_paths,

			     List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			     int localdb_nmismatches_allowed,

			     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
				     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			     Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
			     Path_ext_anchor_cache_T ext_anchor_cache,
			     Path_segment_cache_T segment_cache,
				     bool segment_prefix_skip_p,
				     bool segment_qstart_anchor_cache_p,
				     bool path_qend_local_cache_p,
				     bool emitted_path_dedup_p) {
  
  List_T qstart_sense_paths = NULL, qstart_antisense_paths = NULL,
    qend_sense_paths = NULL, qend_antisense_paths = NULL, p;
  List_T left_chain;
#ifdef GSNAP_CHAIN_ROOT_REPLAY
  List_T root_chain;
#endif
  Univdiag_T left_univdiag;
  Univcoord_T middle_univdiagonal;
  T path;

  bool qstart_innerp, qend_innerp;
  bool ext_anchor_reused_p = false;
#ifdef GSNAP_CHAIN_ROOT_REPLAY
  bool chain_root_materialized_p = false;
#endif
#ifdef GSNAP_EXT_ANCHOR_PROOF
  struct Accel_profile_ext_anchor_key_T ext_shadow_key;
#endif
  int found_sensedir;
  int trimpos, exon_origin, qstart, qend, pos5, pos3;
  int nmatches, nmismatches_to_trimpos;
  bool splice5p, splice3p;
  Splicetype_T splicetype5, splicetype3;
  double ambig_prob_5, ambig_prob_3;

  /* Create left chain and solve => qstart paths */
  left_univdiag = (Univdiag_T) List_head(right_chain);
  middle_univdiagonal = left_univdiag->univdiagonal;

  p = List_next(right_chain);
  if (segment_prefix_skip_p == true && p != NULL && method == SEGMENT1 && paired_end_p == false) {
    Univdiag_T next_univdiag = (Univdiag_T) List_head(p);

    if (next_univdiag->univdiagonal < chrhigh &&
	segment_qend_cache_contains(segment_cache,left_univdiag,next_univdiag,
				    querylength,localdb_nmismatches_allowed,
				    chrnum,chroffset,chrhigh,plusp,genestrand,
				    first_read_p,find_splices_p) == true) {
      /* A previous chain with the same head+next proved that next cannot
	 attach, so the remaining tail cannot affect qend or the combined path. */
      return;
    }
  }

  if (segment_qstart_anchor_cache_p == true &&
      method == SEGMENT1 && paired_end_p == false &&
      segment_cache_lookup(&qstart_sense_paths,&qstart_antisense_paths,segment_cache,left_univdiag,
			  querylength,localdb_nmismatches_allowed,chrnum,chroffset,chrhigh,
			  plusp,genestrand,first_read_p,find_splices_p,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			  hitlistpool,transcriptpool) == true) {
    /* Reused an exact single-anchor qstart solve. */
    Accel_profile_pathsolve_cache("qstart_anchor",/*hitp*/true);

  } else {
    if (segment_qstart_anchor_cache_p == true && method == SEGMENT1 && paired_end_p == false) {
      Accel_profile_pathsolve_cache("qstart_anchor",/*hitp*/false);
    }
    left_chain = Univdiagpool_push(NULL,univdiagpool,left_univdiag->qstart,left_univdiag->qend,
				   left_univdiag->nmismatches,middle_univdiagonal
				   univdiagpool_trace(__FILE__,__LINE__));
#ifdef GSNAP_EXT_ANCHOR_PROOF
    (void) profile_ext_anchor_shadow(&ext_shadow_key,
				     left_univdiag,/*solve_side:qstart*/0U,
			      queryptr,querylength,max_nblocks,localdb_nmismatches_allowed,
			      stage1,knownsplicing,knownindels,
			      query_compress,query_compress_fwd,query_compress_rev,
			      chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,
			      first_read_p,spliceendsgen,method,find_splices_p);
    ext_anchor_reused_p = ext_anchor_result_lookup(&qstart_sense_paths,&qstart_antisense_paths,
					     ext_anchor_cache,&ext_shadow_key,
					     intlistpool,univcoordlistpool,listpool,pathpool,
					     vectorpool,hitlistpool,transcriptpool);
#endif
    if (ext_anchor_reused_p == false) {
      solve_chain_qstart(&qstart_sense_paths,&qstart_antisense_paths,left_chain,
		       queryptr,querylength,mismatch_positions_alloc,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       localdb_nmismatches_allowed,

		       stage1,knownsplicing,knownindels,
		       query_compress,query_compress_fwd,query_compress_rev,
		       chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		       intlistpool,univcoordlistpool,listpool,
		       pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       method,find_splices_p,/*segment_cache*/(Path_segment_cache_T) NULL,
		       /*segment_prefix_skip_p*/false,/*skip_initial_local_p*/false);
    }
#ifdef GSNAP_EXT_ANCHOR_PROOF
    Accel_profile_ext_anchor_shadow_end();
    if (ext_anchor_reused_p == false) {
      (void) ext_anchor_result_store(ext_anchor_cache,&ext_shadow_key,
				     qstart_sense_paths,qstart_antisense_paths,
				     intlistpool,univcoordlistpool,listpool,pathpool,
				     vectorpool,hitlistpool,transcriptpool);
    }
#endif
    Univdiagpool_gc(&left_chain,univdiagpool
		    univdiagpool_trace(__FILE__,__LINE__));

    if (segment_qstart_anchor_cache_p == true &&
	method == SEGMENT1 && paired_end_p == false) {
      segment_cache_insert(segment_cache,left_univdiag,qstart_sense_paths,qstart_antisense_paths,
			  querylength,localdb_nmismatches_allowed,chrnum,chroffset,chrhigh,
			  plusp,genestrand,first_read_p,find_splices_p,
			  intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			  hitlistpool,transcriptpool);
    }
  }

#ifdef GSNAP_CHAIN_ROOT_REPLAY
  root_chain = Univdiagpool_push(NULL,univdiagpool,left_univdiag->qstart,
				left_univdiag->qend,left_univdiag->nmismatches,
				middle_univdiagonal
				univdiagpool_trace(__FILE__,__LINE__));
  (void) profile_ext_anchor_shadow(&ext_shadow_key,left_univdiag,
				   /*solve_side:qend*/1U,queryptr,querylength,max_nblocks,
				   localdb_nmismatches_allowed,stage1,knownsplicing,knownindels,
				   query_compress,query_compress_fwd,query_compress_rev,
				   chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,
				   first_read_p,spliceendsgen,method,find_splices_p);
  if (ext_anchor_result_replay(&qend_sense_paths,&qend_antisense_paths,
			       ext_anchor_cache,&ext_shadow_key,intlistpool,univcoordlistpool,
			       listpool,pathpool,vectorpool,hitlistpool,transcriptpool) == false) {
    solve_chain_qend(&qend_sense_paths,&qend_antisense_paths,root_chain,
		     queryptr,querylength,mismatch_positions_alloc,
		     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
		     stage1,knownsplicing,knownindels,
		     query_compress,query_compress_fwd,query_compress_rev,
		     chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
		     univdiagpool,vectorpool,hitlistpool,spliceendsgen,method,find_splices_p,
		     segment_cache,segment_prefix_skip_p,path_qend_local_cache_p,
		     /*skip_initial_local_p*/false);
    ext_anchor_result_admit(ext_anchor_cache,&ext_shadow_key,
			    qend_sense_paths,qend_antisense_paths,
			    intlistpool,univcoordlistpool,listpool,pathpool,vectorpool,
			    hitlistpool,transcriptpool);
  }
  Univdiagpool_gc(&root_chain,univdiagpool univdiagpool_trace(__FILE__,__LINE__));
  chain_root_materialized_p = true;
#endif

  /* Solve given right chain => qend paths */
  solve_chain_qend(&qend_sense_paths,&qend_antisense_paths,right_chain,
		   queryptr,querylength,mismatch_positions_alloc,
		   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		   localdb_tally,localdb_overflowp,max_nblocks,
		   localdb_nmismatches_allowed,

		   stage1,knownsplicing,knownindels,
		   query_compress,query_compress_fwd,query_compress_rev,
		   chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		   intlistpool,univcoordlistpool,listpool,
		   pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		   method,find_splices_p,segment_cache,segment_prefix_skip_p,
		   path_qend_local_cache_p,/*skip_initial_local_p*/
#ifdef GSNAP_CHAIN_ROOT_REPLAY
		   chain_root_materialized_p
#else
		   false
#endif
		   );

  /* Combine qstart and qend paths */
  solve_combine(&(*found_score),&(*end_score),

		&(*unextended_sense_paths),&(*unextended_antisense_paths),
		&(*sense_paths),&(*antisense_paths),

		qstart_sense_paths,qstart_antisense_paths,
		qend_sense_paths,qend_antisense_paths,
		/*qstart_borrowed_p*/ext_anchor_reused_p,/*qend_borrowed_p*/false,

		middle_univdiagonal,queryseq,queryptr,querylength,mismatch_positions_alloc,
		novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		localdb_tally,localdb_overflowp,max_nblocks,

		stage1,knownsplicing,knownindels,
		query_compress,query_compress_fwd,query_compress_rev,
		chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
		intlistpool,uintlistpool,univcoordlistpool,listpool,
		pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		method,find_splices_p,segment_cache,emitted_path_dedup_p);
  return;
}


void
Path_solve_from_right_chain (int *found_score, int *end_score,

			     List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
			     List_T *sense_paths, List_T *antisense_paths,

			     List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

			     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
			     int localdb_nmismatches_allowed,

			     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
			     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
			     Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p) {

  path_solve_from_right_chain(&(*found_score),&(*end_score),
			      &(*unextended_sense_paths),&(*unextended_antisense_paths),
			      &(*sense_paths),&(*antisense_paths),
			      right_chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,
			      localdb_nmismatches_allowed,
			      stage1,knownsplicing,knownindels,
			      query_compress,query_compress_fwd,query_compress_rev,
			      chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
				      intlistpool,uintlistpool,univcoordlistpool,listpool,
				      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				      spliceendsgen,method,find_splices_p,
				      /*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,
				      /*segment_cache*/(Path_segment_cache_T) NULL,
				      /*segment_prefix_skip_p*/false,
				      /*segment_qstart_anchor_cache_p*/false,
				      /*path_qend_local_cache_p*/false,
				      /*emitted_path_dedup_p*/false);
  return;
}


void
Path_solve_from_right_chain_cached (int *found_score, int *end_score,

				    List_T *unextended_sense_paths, List_T *unextended_antisense_paths,
				    List_T *sense_paths, List_T *antisense_paths,

				    List_T right_chain, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc,

				    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
				    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
				    int localdb_nmismatches_allowed,

				    Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
				    Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
				    bool plusp, int genestrand, bool paired_end_p, bool first_read_p,
				    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
					    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
				    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
				    Spliceendsgen_T spliceendsgen, Method_T method, bool find_splices_p,
				    Path_ext_anchor_cache_T ext_anchor_cache,
				    Path_segment_cache_T segment_cache,
					    bool segment_prefix_skip_p,
					    bool segment_qstart_anchor_cache_p,
					    bool path_qend_local_cache_p,
					    bool emitted_path_dedup_p) {

  path_solve_from_right_chain(&(*found_score),&(*end_score),
			      &(*unextended_sense_paths),&(*unextended_antisense_paths),
			      &(*sense_paths),&(*antisense_paths),
			      right_chain,queryseq,queryptr,querylength,mismatch_positions_alloc,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,
			      localdb_nmismatches_allowed,
			      stage1,knownsplicing,knownindels,
			      query_compress,query_compress_fwd,query_compress_rev,
			      chrnum,chroffset,chrhigh,plusp,genestrand,paired_end_p,first_read_p,
				      intlistpool,uintlistpool,univcoordlistpool,listpool,
				      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,
				      spliceendsgen,method,find_splices_p,ext_anchor_cache,segment_cache,
				      segment_prefix_skip_p,segment_qstart_anchor_cache_p,
				      path_qend_local_cache_p,emitted_path_dedup_p);
  return;
}



#ifdef REMAP_PATHS
/* Called only by Path_extend, so can set childp to be true here */
static List_T
remap_paths (int *found_score, int *end_score, List_T paths,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     Shortread_T queryseq, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	     Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
	     Transcriptpool_T transcriptpool, Vectorpool_T vectorpool,
	     Pathpool_T pathpool, Hitlistpool_T hitlistpool) {

  List_T result_paths = NULL, repaired_paths, p;
  T path;

  for (p = paths; p != NULL; p = List_next(p)) {
    path = (T) List_head(p);
    path->nmatches = -1;
    Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);

    if (transcriptome == NULL) {
      path->childp = true;
      result_paths = Hitlist_push(result_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    } else if ((repaired_paths = try_repairs(&(*found_score),&(*end_score),
					     path,query_compress_fwd,query_compress_rev,queryseq,
					     intlistpool,uintlistpool,univcoordlistpool,listpool,
					     transcriptpool,pathpool,vectorpool,hitlistpool)) != NULL) {

      result_paths = List_append(repaired_paths,result_paths);
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    } else {
      path->childp = true;
      result_paths = Hitlist_push(result_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  Hitlistpool_free_list(&paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));

  /* assert(result_paths != NULL); -- incorrect assumption */
  return result_paths;
}
#endif


/* Either unextended_paths or complete_paths is non-NULL */
/* Sets completep field for parent, but not for children, which are carried with parent */
/* Returns complete paths.  Calls Path_eval_nmatches */
List_T
Path_extend (int *found_score, int *end_score, List_T *global_unextended_paths,
	     T original_path, Shortread_T queryseq, char *queryptr, int querylength, int *mismatch_positions_alloc, 
	     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
	     int localdb_nmismatches_allowed, 

	     Stage1_T stage1, Knownsplicing_T knownsplicing, Knownindels_T knownindels,
	     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     int genestrand, bool paired_end_p, bool lowp,
	     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
	     Spliceendsgen_T spliceendsgen, bool extend_qstart_p, bool extend_qend_p) {

  List_T unextended_paths = NULL, complete_paths = NULL;

  List_T qstart_paths = NULL, qend_paths = NULL, p;
  T path;
  bool qstart_innerp, qend_innerp;
  Chrnum_T chrnum = original_path->chrnum;
  Univcoord_T chroffset = original_path->chroffset, chrhigh = original_path->chrhigh;

#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long profile_path_hash = 1469598103934665603ULL;
  Univcoordlist_T profile_u;
  Intlist_T profile_i;
  List_T profile_j;
  unsigned long long profile_value;
  int profile_k;

#define PROFILE_PATH_HASH_VALUE(value) do { \
    profile_value = (unsigned long long) (value); \
    for (profile_k = 0; profile_k < 8; profile_k++) { \
      profile_path_hash ^= profile_value & 0xffU; \
      profile_path_hash *= 1099511628211ULL; \
      profile_value >>= 8; \
    } \
  } while (0)
  PROFILE_PATH_HASH_VALUE(Path_effective_sensedir(original_path));
  for (profile_u = original_path->univdiagonals; profile_u != NULL;
	 profile_u = Univcoordlist_next(profile_u)) {
    PROFILE_PATH_HASH_VALUE(Univcoordlist_head(profile_u));
  }
  for (profile_i = original_path->endpoints; profile_i != NULL; profile_i = Intlist_next(profile_i)) {
    PROFILE_PATH_HASH_VALUE(Intlist_head(profile_i));
  }
  for (profile_j = original_path->junctions; profile_j != NULL; profile_j = List_next(profile_j)) {
    PROFILE_PATH_HASH_VALUE(Junction_nindels((Junction_T) List_head(profile_j)));
  }
  Accel_profile_path_extend_begin(profile_path_hash,(int) original_path->chrnum,
				  Intlist_head(original_path->endpoints),
				  Intlist_last_value(original_path->endpoints),
				  original_path->plusp,Path_effective_sensedir(original_path));
#undef PROFILE_PATH_HASH_VALUE
#endif


  debug13(printf("Entered Path_extend, sensedir %s, with path\n",
		 original_path->sensedir == SENSE_FORWARD ? "fwd" : "anti"));
  debug13(Path_print(original_path));

  assert(original_path->transcriptome_method_p == false);
  assert(original_path->extendedp == false);

  Path_expect_fwd(original_path);
  
  original_path->extendedp = true;

  if (paired_end_p == false) {
    qstart_innerp = qend_innerp = false;
  } else if (lowp) {
    qstart_innerp = false;
    qend_innerp = true;
  } else {
    qstart_innerp = true;
    qend_innerp = false;
  }

  assert(extend_qstart_p == true);
  assert(extend_qend_p == true);

  /* Extend qstart with salvagep being true */
  /* qstart_paths should be fwd */
  compute_qstart_local(&qstart_paths,/*depth*/0,original_path,queryptr,querylength,
		       mismatch_positions_alloc,novel_univdiags_alloc,
		       novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       stage1,knownsplicing,knownindels,
		       query_compress,chrnum,chroffset,chrhigh,intlistpool,
		       univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       localdb_nmismatches_allowed,original_path->plusp,genestrand,
		       original_path->sensedir,qstart_innerp,/*find_splices_p*/true,/*salvagep*/true);
  
  /* Extend qend with salvagep being true */
  for (p = qstart_paths; p != NULL; p = List_next(p)) {
    path = Path_reverse((Path_T) List_head(p),/*expect_fwd_p*/false);
    compute_qend_local(&qend_paths,/*depth*/0,path,queryptr,querylength,
		       mismatch_positions_alloc,novel_univdiags_alloc,
		       novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,
		       stage1,knownsplicing,knownindels,
		       query_compress,chrnum,chroffset,chrhigh,intlistpool,
		       univcoordlistpool,listpool,pathpool,transcriptpool,
		       univdiagpool,vectorpool,hitlistpool,spliceendsgen,
		       localdb_nmismatches_allowed,original_path->plusp,genestrand,
		       original_path->sensedir,qend_innerp,/*find_splices_p*/true,/*salvagep*/true,
			       /*segment_cache*/(Path_segment_cache_T) NULL,
			       /*path_qend_local_cache_p*/false);
  }
  Path_gc(&qstart_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  /* Hitlistpool_free_list(&qstart_paths,hitlistpool 
     hitlistpool_trace(__FILE__,__LINE__)); */

  complete_paths = unextended_paths = (List_T) NULL;
  for (p = qend_paths; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    path = Path_reverse(path,/*expect_fwd_p*/true);
    path->nmatches = -1;
    Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
    if (Intlist_head(path->endpoints) <= 8 && querylength - Intlist_last_value(path->endpoints) <= 8) {
      complete_paths = Hitlist_push(complete_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else {
      unextended_paths = Hitlist_push(unextended_paths,hitlistpool,(void *) path
				      hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  if (complete_paths != NULL) {
    debug11(printf("complete paths is not NULL\n"));
    original_path->completep = true;
    Path_gc(&unextended_paths,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

#ifdef REMAP_PATHS
    complete_paths = remap_paths(&(*found_score),&(*end_score),complete_paths,
				 query_compress_fwd,query_compress_rev,queryseq,
				 intlistpool,uintlistpool,univcoordlistpool,
				 listpool,transcriptpool,vectorpool,pathpool,hitlistpool);
#endif

    debug11(printf("(3) Path_extend returning %d complete paths\n",
		   List_length(complete_paths)));

    *global_unextended_paths = (List_T) NULL;
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_path_extend_end();
#endif
    return complete_paths;

  } else {
    debug11(printf("complete paths is NULL\n"));
#ifdef REMAP_PATHS
    unextended_paths = remap_paths(&(*found_score),&(*end_score),unextended_paths,
				   query_compress_fwd,query_compress_rev,queryseq,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,transcriptpool,vectorpool,pathpool,hitlistpool);
#endif

    debug11(printf("(4) Path_extend returning %d unextended paths\n",
		   List_length(unextended_paths)));

    *global_unextended_paths = unextended_paths;
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_path_extend_end();
#endif
    return (List_T) NULL;
  }
}


/* Modifies path and does not create a copy */
/* innerp is true by definition (resolve means to determine inner splices) */
void
Path_qstart_resolve (int *found_score, int *end_score, T path,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal, char *queryptr, int querylength,
		     Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		     int localdb_nmismatches_allowed,
		     
		     Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		     Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     int genestrand, 
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		     Univdiagpool_T univdiagpool, Vectorpool_T vectorpool) {

  Univcoord_T new_univdiagonal, indel_univdiagonal;
  int nosplice_trimpos, low_qstart, qend;
  int supporti, supportj;

  Univcoordlist_T univdiagonals;
  Intlist_T nmismatches, ref_nmismatches, endpoints;
  List_T junctions;
  int ignore;

  Univcoord_T *new_univdiagonals, nearest_univdiagonal, univdiagonal_i,
    univdiagonal, middle_univdiagonal;
  Univcoord_T *distal_positions;
  Univdiag_T *distal_univdiags;

  int *order;
  int *medial_nmismatches, *medial_qpositions;
  double *medial_probs, *distal_probs;
  int *medial_support, *distal_support;
  int n_medial_partners, npartners, ndiagonals, goodi, i, j;
  bool sense_forward_p, multiple_splice_qpos_p;
  Junction_T junction;

  Chrpos_T splice_distance;
  int type;
  int cont_trimpos, local_trimpos, trimpos5, trimpos3, qstart, pos3;
  int best_endpoint, medial_qpos, splice_qpos;
  int nindels, indel_pos;
  int cont_nmismatches, cont_ref_nmismatches, local_nmismatches, resolve_nmismatches,
    nmismatches_i, nmismatches_j, nmismatches_indel;
  int ref_nmismatches_i, ref_nmismatches_j, ref_nmismatches_indel;
  double best_prob, donor_prob, acceptor_prob;


  debug3(printf("Entered Path_qstart_resolve, sensedir %d, univdiagonal %u..%u and path\n",
		path->sensedir,low_univdiagonal,high_univdiagonal));
  debug3(Path_print(path));

  Path_expect_fwd(path);

  univdiagonal = Univcoordlist_head(path->univdiagonals);
  debug3(printf("Original univdiagonal: %u\n",univdiagonal));

  if (high_univdiagonal > univdiagonal) {
    debug3(printf("Changing high_univdiagonal to be %u\n",univdiagonal));
    high_univdiagonal = univdiagonal;
  }

  if (high_univdiagonal <= low_univdiagonal) {
    /* Ends overlap, so no way to resolve */
    debug3(printf("Interval %u..%u is negative, so not resolving\n",
		  low_univdiagonal,high_univdiagonal));

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[path->chrnum] == true) {
    /* No splicing on circular chromosomes */
    debug3(printf("No splicing on circular chromosomes\n"));
#endif

  } else {
    debug3(printf("Interval %u..%u is positive, so resolving\n",
		  low_univdiagonal,high_univdiagonal));

    qstart = Intlist_head(path->endpoints);
    pos3 = Intlist_second_value(path->endpoints);

#if 0
    /* This can lead to poor results */
    /* Compute continuation to 0 for baseline */
    cont_trimpos = 0;
    cont_nmismatches =
      Genomebits_count_mismatches_substring(&cont_ref_nmismatches,genomebits,genomebits_alt,
					    query_compress,univdiagonal,querylength,
					    /*pos5*/0,pos3,path->plusp,genestrand);
#else
    cont_trimpos = Genomebits_trim_qstart(&cont_nmismatches,query_compress,genomebits,
					  univdiagonal,querylength,
					  /*pos5*/0,pos3,path->plusp,genestrand);
#endif

    new_univdiagonals =
      Spliceends_qstart_resolve(&ndiagonals,/*pos3*/qstart,querylength,
				low_univdiagonal,high_univdiagonal,
				query_compress,queryptr,path->plusp,genestrand,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,stage1,
				localdb_nmismatches_allowed);
    
    if (ndiagonals == 0) {
      debug3(printf("Spliceends_qstart_resolve returns nothing\n"));
      local_trimpos = Genomebits_trim_qstart(&local_nmismatches,query_compress,genomebits,
					     univdiagonal,querylength,
					     /*pos5*/0,pos3,path->plusp,genestrand);
      debug3(printf("(1) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
      if (local_trimpos >= pos3) {
	debug3(printf("Need to trim exon completely\n"));
	Path_trim_qstart_n(/*noutside*/local_trimpos - Intlist_head(path->endpoints),path,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

      } else {
	Intlist_head_set(path->endpoints,local_trimpos);
	Intlist_head_set(path->nmismatches,local_nmismatches);
	Intlist_head_set(path->ref_nmismatches,local_nmismatches);

	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
      }

    } else if ((nearest_univdiagonal = new_univdiagonals[0]) == univdiagonal) {
      /* Continuation of univdiagonal */
      debug3(printf("Continuation of univdiagonal now with %d nmismatches\n",local_nmismatches));
      local_trimpos = Genomebits_trim_qstart(&local_nmismatches,query_compress,genomebits,
					     nearest_univdiagonal,querylength,
					     /*pos5*/0,pos3,path->plusp,genestrand);
      debug3(printf("(2) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
      if (local_trimpos >= pos3) {
	debug3(printf("Need to trim exon completely\n"));
	Path_trim_qstart_n(/*noutside*/local_trimpos - Intlist_head(path->endpoints),path,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

      } else {
	Intlist_head_set(path->endpoints,local_trimpos);
#if 0
	if (Intlist_head(path->nmismatches) >= 0) {
	  Intlist_head_incr(path->nmismatches,local_nmismatches);
	  Intlist_head_incr(path->ref_nmismatches,local_nmismatches);
	}
#else
	Intlist_head_set(path->nmismatches,local_nmismatches);
	Intlist_head_set(path->ref_nmismatches,local_nmismatches);
#endif

	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
      }

    } else if (nearest_univdiagonal > univdiagonal) {
      /* Insertion.  Skip */

    } else if (nearest_univdiagonal + max_deletionlen >= univdiagonal) {
      /* Deletion (or short intron) */
      nindels = univdiagonal - nearest_univdiagonal;
      trimpos5 = Genomebits_trim_qstart(&resolve_nmismatches,query_compress,genomebits,
					/*univdiagonal*/nearest_univdiagonal,querylength,
					/*pos5*/0,pos3,path->plusp,genestrand);

      if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						     &ref_nmismatches_i,&ref_nmismatches_j,
						     /*univdiagonal_i*/nearest_univdiagonal,/*indels*/-nindels,
						     Path_chrhigh(path),
						     /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						     /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						     /*pos5*/trimpos5,pos3,querylength,
						     stage1->indelinfo,path->plusp,genestrand,
						     /*want_lowest_coordinate_p*/true)) <= 0) {
	debug3(printf("Deletion or short intron fails\n"));
	
      } else if ((pos3 - trimpos5) - nmismatches_i - nmismatches_j <= (pos3 - cont_trimpos) - cont_nmismatches) {
	debug3(printf("Deletion or short intron yields fewer matches than continuation: (%d - %d) - %d - %d vs (%d - %d) - %d\n",
		      pos3,trimpos5,nmismatches_i,nmismatches_j,pos3,cont_trimpos,cont_nmismatches));

      } else {
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos5 /* was 0 */
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* deletionpos = nearest_univdiagonal + indel_pos; */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = nearest_univdiagonal;
	
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
      }
      
    } else if (splicingp == false) {
      /* Unable to try splicing */

    } else if (ndiagonals == 1) {

#if 0
      new_univdiagonal = 0;
      /* Splice_nomiddle call type 1 */
      if ((splice_qpos = Splice_nomiddle(&type,&trimpos5,&trimpos3,&supporti,&supportj,
					 &nindels,&indel_pos,&nmismatches_i,&nmismatches_j,&nmismatches_indel,
					 &ref_nmismatches_i,&ref_nmismatches_j,&ref_nmismatches_indel,
					 &donor_prob,&acceptor_prob,
					 /*univdiagonal_i*/nearest_univdiagonal,/*univdiagonal_j*/univdiagonal,
					 query_compress,path->plusp,Path_chroffset(path),Path_chrhigh(path),
					 /*pos5*/0,pos3,querylength,stage1->indelinfo,stage1->spliceinfo,knownsplicing,
					 intlistpool,/*sense_forward_p*/(path->sensedir == SENSE_FORWARD) ? true : false,
					 genestrand,/*trim5p*/true,/*trim3p*/false,/*innerp*/true,
					 /*allow_indels_p*/true)) < 0) {
	debug3(printf("Single splice finds no solution\n"));
	local_trimpos = Genomebits_trim_qstart(&local_nmismatches,query_compress,
					       /*bits*/genomebits,univdiagonal,querylength,
					       /*pos5*/0,pos3,path->plusp,genestrand);
	debug3(printf("(3) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos >= pos3) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qstart_n(/*noutside*/local_trimpos - Intlist_head(path->endpoints),path,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	
	  path->splice5p = false;
	  path->splicetype5 = NO_SPLICE;
	  path->ambig_prob_5 = 0.0;
	}

      } else if ((pos3 - trimpos5) - nmismatches_i - nmismatches_j <= (pos3 - cont_trimpos) - cont_nmismatches) {
	debug3(printf("Deletion or short intron yields fewer matches than continuation: (%d - %d) - %d - %d vs (%d - %d) - %d\n",
		      pos3,trimpos5,nmismatches_i,nmismatches_j,pos3,cont_trimpos,cont_nmismatches));

      } else if (nindels == 0) {
	/* Single splice succeeds */
	debug3(printf("Single splice succeeds with trimpos5 %d\n",trimpos5));
	Intlist_head_set(path->endpoints,splice_qpos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos5   /* was 0 */
					   intlistpool_trace(__FILE__,__LINE__));
      
	splice_distance = univdiagonal - nearest_univdiagonal;
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = nearest_univdiagonal;
      
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;

	Transcript_list_gc(&path->transcripts,listpool,transcriptpool);
	Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool);

      } else if (splice_qpos < indel_pos) {
	/* Single spliceindel succeeds */
	/* Push indel (based on left) then splice.  splice is distal, indel is medial */
	debug13(printf("splice_qpos %d < indel_pos %d\n",splice_qpos,indel_pos));
	
	middle_univdiagonal = univdiagonal - nindels; /* nindels = univdiagonal - middle_univdiagonal */
	splice_distance = middle_univdiagonal - nearest_univdiagonal;
	
	debug3(printf("Splice_nomiddle (qstart, sense): %d indels at %d, mismatches %d, then splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
		      nindels,indel_pos,nmismatches_indel,
		      /*low_qstart*/0,qend,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob));

	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,splice_qpos
					   intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos5    /* was low_qstart */
					   intlistpool_trace(__FILE__,__LINE__));

	/* Indel first */
	if (nindels < 0) {
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_insertion(-nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	} else {
	  /* deletionpos = low_left + splice_distance + indel_pos; -- qstart; add splice_distance if indel first */
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_deletion(nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	}
	
	/* Splice second */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_indel
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_indel
						 intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,middle_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	
      } else {
	/* Single spliceindel succeeds */
	/* Push splice then indel (based on nearest_univdiagonal).  indel is distal, splice is medial */
	debug13(printf("splice_qpos %d >= indel_pos %d\n",splice_qpos,indel_pos));
	
	middle_univdiagonal = nearest_univdiagonal + nindels; /* nindels = middle_univdiagonal - nearest_univdiagonal */
	splice_distance = univdiagonal - middle_univdiagonal;
	
	debug3(printf("Splice_nomiddle (qstart, sense): splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f, then %d indels at %d, nmismatches %d\n",
		      /*low_qstart*/0,qend,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob,
		      nindels,indel_pos,nmismatches_indel));

	Intlist_head_set(path->endpoints,splice_qpos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,indel_pos
					   intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos5   /* was low_qstart */
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* Splice first */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* Indel second */
	if (nindels < 0) {
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_insertion(-nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	} else {
	  /* deletionpos = low_left + indel_pos; -- qstart; do not add splice_distance if indel second */
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_deletion(nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	}
	
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_indel
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_indel
						 intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,middle_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
      }

#else
      if (Splice_resolve_qstart(&univdiagonals,&nmismatches,&ref_nmismatches,&endpoints,&junctions,
				/*univdiagonal_i*/nearest_univdiagonal,/*univdiagonal_j*/univdiagonal,
				stage1,query_compress,queryptr,
				path->plusp,Path_chroffset(path),Path_chrhigh(path),
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,localdb,
				localdb_nmismatches_allowed,
				/*pos5*/0,pos3,querylength,cont_trimpos,cont_nmismatches,
				indelinfo,stage1->spliceinfo,knownsplicing,univcoordlistpool,intlistpool,listpool,pathpool,
				path->sensedir,genestrand,/*trim5p*/true,/*trim3p*/false,/*innerp*/true,
				/*allow_indels_p*/true,/*find_middle_p*/false,/*middle_rank*/0) == 0) {


	debug3(printf("Single splice finds no solution\n"));
	local_trimpos = Genomebits_trim_qstart(&local_nmismatches,query_compress,
					       /*bits*/genomebits,univdiagonal,querylength,
					       /*pos5*/0,pos3,path->plusp,genestrand);
	debug3(printf("(3) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos >= pos3) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qstart_n(/*noutside*/local_trimpos - Intlist_head(path->endpoints),path,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	
	  path->splice5p = false;
	  path->splicetype5 = NO_SPLICE;
	  path->ambig_prob_5 = 0.0;
	}

#if 0
      } else if ((pos3 - trimpos5) - nmismatches_i - nmismatches_j <= (pos3 - cont_trimpos) - cont_nmismatches) {
	debug3(printf("Single splice yields fewer matches than continuation: (%d - %d) - %d - %d vs (%d - %d) - %d\n",
		      pos3,trimpos5,nmismatches_i,nmismatches_j,pos3,cont_trimpos,cont_nmismatches));
#endif
	
      } else {
	/* No need to copy path */
	debug3(printf("Adding results from Splice_resolve_qstart\n"));

	path->univdiagonals = Univcoordlist_append(univdiagonals,path->univdiagonals);
	path->nmismatches = Intlistpool_pop(path->nmismatches,intlistpool,&ignore
					       intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlist_append(nmismatches,path->nmismatches);
	path->ref_nmismatches = Intlistpool_pop(path->ref_nmismatches,intlistpool,&ignore
						   intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlist_append(ref_nmismatches,path->ref_nmismatches);
	path->endpoints = Intlistpool_pop(path->endpoints,intlistpool,&ignore
					     intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlist_append(endpoints,path->endpoints);
	path->junctions = List_append(junctions,path->junctions);
	
	debug3(printf("AFTER call to Splice_resolve_qstart\n"));
	debug3(Path_print(path));
      }
#endif

    } else {
      /* Multiple splices */
      debug3(printf("(1) Multiple diagonals: %d\n",ndiagonals));
      npartners = 0;

      medial_qpositions = Vectorpool_new_intvector(vectorpool,ndiagonals+1); /* Need +1 for Sedgesort_order */
      distal_positions = Vectorpool_new_univcoordvector(vectorpool,ndiagonals);

      distal_support = Vectorpool_new_intvector(vectorpool,ndiagonals);
      distal_univdiags = Univdiag_new_array(ndiagonals);

      medial_probs = Vectorpool_new_doublevector(vectorpool,ndiagonals);
      distal_probs = Vectorpool_new_doublevector(vectorpool,ndiagonals);
      medial_support = Vectorpool_new_intvector(vectorpool,ndiagonals);
      medial_nmismatches = Vectorpool_new_intvector(vectorpool,ndiagonals);

      sense_forward_p = (path->sensedir == SENSE_FORWARD) ? true : false;

      multiple_splice_qpos_p = false;
      for (i = 0; i < ndiagonals; i++) {
	univdiagonal_i = new_univdiagonals[i];

#if 0
	debug3(printf("Calling Splice_nomiddle with %u and %u\n",univdiagonal_i,univdiagonal));
	/* Splice_nomiddle call type 2 */
	if ((splice_qpos = Splice_nomiddle(&type,&trimpos5,&trimpos3,&supporti,&supportj,
					   &nindels,&indel_pos,&nmismatches_i,&nmismatches_j,&nmismatches_indel,
					   &ref_nmismatches_i,&ref_nmismatches_j,&ref_nmismatches_indel,
					   &donor_prob,&acceptor_prob,
					   univdiagonal_i,/*univdiagonal_j*/univdiagonal,
					   query_compress,path->plusp,Path_chroffset(path),Path_chrhigh(path),
					   /*pos5*/0,pos3,querylength,
					   stage1->indelinfo,stage1->spliceinfo,knownsplicing,intlistpool,
					   sense_forward_p,genestrand,/*trim5p*/true,/*trim3p*/false,/*innerp*/true,
					   /*allow_indels_p*/false)) < 0) {
	  /* Skip */

	} else if (npartners == 0) {
	  medial_qpos = medial_qpositions[0] = splice_qpos;
	  medial_support[0] = supportj;
	  medial_nmismatches[0] = nmismatches_j;
	  distal_positions[0] = univdiagonal_i - querylength + splice_qpos;

	  distal_support[0] = supporti;
	  distal_univdiags[0]->univdiagonal = univdiagonal_i;
	  distal_univdiags[0]->qstart = trimpos5;
	  distal_univdiags[0]->qend = splice_qpos;
	  distal_univdiags[0]->nmismatches = nmismatches_i;

	  if (path->plusp == sense_forward_p) {
	    medial_probs[0] = acceptor_prob;
	    distal_probs[0] = donor_prob;
	  } else {
	    medial_probs[0] = donor_prob;
	    distal_probs[0] = acceptor_prob;
	  }
	  npartners++;
	  
	} else {
	  if (splice_qpos != medial_qpos) {
	    multiple_splice_qpos_p = true;
	  }

	  medial_qpositions[npartners] = splice_qpos;
	  medial_support[npartners] = supportj;
	  medial_nmismatches[npartners] = nmismatches_j;
	  distal_positions[npartners] = univdiagonal_i - querylength + splice_qpos;

	  distal_support[npartners] = supporti;
	  distal_univdiags[npartners]->univdiagonal = univdiagonal_i;
	  distal_univdiags[npartners]->qstart = trimpos5;
	  distal_univdiags[npartners]->qend = splice_qpos;
	  distal_univdiags[npartners]->nmismatches = nmismatches_i;

	  if (path->plusp == sense_forward_p) {
	    medial_probs[npartners] = acceptor_prob;
	    distal_probs[npartners] = donor_prob;
	  } else {
	    medial_probs[npartners] = donor_prob;
	    distal_probs[npartners] = acceptor_prob;
	  }
	  npartners++;
	}
	debug3(printf("%u %u => %d, %f, %f\n",univdiagonal_i,univdiagonal - univdiagonal_i,
		      splice_qpos,donor_prob,acceptor_prob));
      }

#else
	/* allow_indels_p is false because Altsplice_T cannot handle indels */
	debug3(printf("Calling Splice_result_qstart with %u and %u\n",univdiagonal_i,univdiagonal));
	if (Splice_resolve_qstart(&univdiagonals,&nmismatches,&ref_nmismatches,&endpoints,&junctions,
				  univdiagonal_i,/*univdiagonal_j*/univdiagonal,
				  stage1,query_compress,queryptr,
				  path->plusp,Path_chroffset(path),Path_chrhigh(path),
				  
				  novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				  localdb_tally,localdb_overflowp,max_nblocks,localdb,
				  localdb_nmismatches_allowed,

				  /*pos5*/0,pos3,querylength,cont_trimpos,cont_nmismatches,
				  stage1->indelinfo,stage1->spliceinfo,knownsplicing,
				  univcoordlistpool,intlistpool,listpool,pathpool,
				  path->sensedir,genestrand,/*trim5p*/true,/*trim3p*/false,/*innerp*/true,
				  /*allow_indels_p*/false,/*find_middle_p*/false,/*middle_rank*/0) == 0) {
	  /* Skip */

	} else {
	  splice_qpos = Intlist_second_value(endpoints);
	  if (npartners == 0) {
	    /* Do not set multiple_splice_qpos_p */
	    medial_qpos = splice_qpos;
	  } else if (splice_qpos != medial_qpos) {
	    multiple_splice_qpos_p = true;
	  }

	  debug3(printf("endpoints %s => supporti %d, supportj %d\n",
			Intlist_to_string(endpoints),Intlist_second_value(endpoints) - Intlist_head(endpoints),
			pos3 - Intlist_second_value(endpoints)));

	  medial_qpositions[npartners] = splice_qpos;
	  medial_support[npartners] = /*supportj*/ pos3 - Intlist_second_value(endpoints);
	  medial_nmismatches[npartners] = /*nmismatches_j*/ Intlist_second_value(nmismatches);
	  distal_positions[npartners] = univdiagonal_i - querylength + splice_qpos;

	  distal_support[npartners] = /*supporti*/ Intlist_second_value(endpoints) - Intlist_head(endpoints);
	  distal_univdiags[npartners]->univdiagonal = univdiagonal_i;
	  distal_univdiags[npartners]->qstart = /*trimpos5*/ Intlist_head(endpoints);
	  distal_univdiags[npartners]->qend = splice_qpos;
	  distal_univdiags[npartners]->nmismatches = /*nmismatches_i*/ Intlist_head(nmismatches);

	  junction = (Junction_T) List_head(junctions);
	  /* good_outer_accept_p = junction->outer_accept_p; -- Used if npartners == 1 */
	  /* assert(Junction_type(junction) == SPLICE_JUNCTION); -- Could have been converted to a deletion */
	  debug3(Junction_print(junction));

	  if (path->plusp == sense_forward_p) {
	    medial_probs[npartners] = junction->acceptor_prob;
	    distal_probs[npartners] = junction->donor_prob;
	  } else {
	    medial_probs[npartners] = junction->donor_prob;
	    distal_probs[npartners] = junction->acceptor_prob;
	  }

	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Junction_list_gc(&junctions,listpool,pathpool);

	  debug3(printf("%u %u => %d, %f, %f\n",univdiagonal_i,univdiagonal - univdiagonal_i,
			splice_qpos,junction->donor_prob,junction->acceptor_prob));
	  npartners++;
	}
      }
#endif

      if (npartners == 0) {
	debug3(printf("(1) Multiple splices find no good solution\n"));
	local_trimpos = Genomebits_trim_qstart(&local_nmismatches,query_compress,
					       /*bits*/genomebits,univdiagonal,querylength,
					       /*pos5*/0,pos3,path->plusp,genestrand);
	debug3(printf("(4) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos >= pos3) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qstart_n(/*noutside*/local_trimpos - Intlist_head(path->endpoints),path,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);

	  path->splice5p = false;
	  path->splicetype5 = NO_SPLICE;
	  path->ambig_prob_5 = 0.0;
	}

      } else if (npartners == 1) {
	/* Single valid splice */
	debug3(printf("(1) Multiple splices yield one valid solution\n"));

	Intlist_head_set(path->endpoints,medial_qpositions[0]);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,distal_univdiags[0]->qstart /* was 0 */
					   intlistpool_trace(__FILE__,__LINE__));
      
	univdiagonal_i = distal_positions[0] - medial_qpositions[0] + querylength;  /* univdiagonal; */
	splice_distance = univdiagonal - univdiagonal_i;
	if (path->plusp == sense_forward_p) {
	  donor_prob = distal_probs[0];
	  acceptor_prob = medial_probs[0];
	} else {
	  donor_prob = medial_probs[0];
	  acceptor_prob = distal_probs[0];
	}
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
      
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,medial_nmismatches[0]);
	Intlist_head_set(path->ref_nmismatches,medial_nmismatches[0]);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,distal_univdiags[0]->nmismatches
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,distal_univdiags[0]->nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
      
	debug3(printf("Pushing univdiagonal %u\n",univdiagonal_i));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,univdiagonal_i
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = univdiagonal_i;
      
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
	
      } else if (multiple_splice_qpos_p == false &&
		 (path->qstart_alts =
		  Altsplice_qstart_new(&best_endpoint,/*boundedp*/false,
				       medial_qpos,/*anchor_qpos*/pos3,querylength,
				       /*medial_univdiagonal*/univdiagonal,medial_nmismatches[0],
				       medial_probs[0],medial_support[0],distal_positions,
				       distal_support,distal_univdiags,distal_probs,
				       /*order_medial_qpos*/NULL,npartners,
				       path->plusp,/*sense_forward_p*/path->sensedir == SENSE_FORWARD,
				       /*innerp*/true,pathpool,vectorpool,/*sort_bydistal_p*/false)) != NULL) {

	debug3(printf("(1) %d partners with same medial qpos => altsplice\n",npartners));
	
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
	
	Intlist_head_set(path->endpoints,best_endpoint);
	Intlist_head_set(path->nmismatches,-1);
	Intlist_head_set(path->ref_nmismatches,-1);
	
      } else {
	/* Determine the best medial_qpos, based on the total prob */
	best_prob = 0.0;
	for (i = 0; i < npartners; i++) {
	  if (medial_probs[i] + distal_probs[i] > best_prob) {
	    best_prob = medial_probs[i] + distal_probs[i];
	    medial_qpos = medial_qpositions[i];
	  }
	}
	debug3(printf("best medial qpos is %d with prob %f\n",medial_qpos,best_prob));

	n_medial_partners = 0;
	for (i = 0; i < npartners; i++) {
#if 0
	  /* Need to use distal_support when a distal indel is found by Spliceends */
	  if (distal_univdiags[i]->univdiagonal < univdiagonal &&
	      Splice_accept_p(splice_qpos,querylength,distal_probs[i],medial_probs[i],
			      /*univdiagonal_i*/distal_univdiags[i]->univdiagonal,
			      /*univdiagonal_j*/univdiagonal,
			      distal_support[i],distal_univdiags[i]->nmismatches,
			      medial_support[i],medial_nmismatches[i],
			      path->plusp,sense_forward_p,/*innerp*/true) == true) {
	    /* Skip */
	  } else if (medial_qpositions[i] == medial_qpos) {
	    n_good_partners += 1;
	    goodi = i;
	  }
#else
	  if (medial_qpositions[i] == medial_qpos) {
	    n_medial_partners += 1;
	    goodi = i;		/* Used in case n_medial_partners == 1 */
	  }
#endif
	}

	debug3(printf("Determined that medial_qpos %d is best, with %d partners\n",
		      medial_qpos,n_medial_partners));
	
	if (n_medial_partners == 0) {
	  /* All good distal probs are paired with a low medial prob.  Leave path as is. */

	} else if (n_medial_partners == 1) {
	  Intlist_head_set(path->endpoints,medial_qpositions[goodi]);
	  path->endpoints = Intlistpool_push(path->endpoints,intlistpool,distal_univdiags[goodi]->qstart /* was 0 */
					     intlistpool_trace(__FILE__,__LINE__));
	  
	  univdiagonal_i = distal_positions[goodi] - medial_qpositions[goodi] + querylength;  /* univdiagonal; */
	  splice_distance = univdiagonal - univdiagonal_i;
	  if (path->plusp == sense_forward_p) {
	    donor_prob = distal_probs[goodi];
	    acceptor_prob = medial_probs[goodi];
	  } else {
	    donor_prob = medial_probs[goodi];
	    acceptor_prob = distal_probs[goodi];
	  }
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_splice(splice_distance,path->sensedir,
								       donor_prob,acceptor_prob,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	  
	  /* For qstart, push j first, then push i */
	  Intlist_head_set(path->nmismatches,medial_nmismatches[goodi]);
	  Intlist_head_set(path->ref_nmismatches,medial_nmismatches[goodi]);
	  path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,distal_univdiags[goodi]->nmismatches
					       intlistpool_trace(__FILE__,__LINE__));
	  path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,distal_univdiags[goodi]->nmismatches
						   intlistpool_trace(__FILE__,__LINE__));
	  
	  debug3(printf("Pushing univdiagonal %u\n",univdiagonal_i));
	  path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,univdiagonal_i
						       univcoordlistpool_trace(__FILE__,__LINE__));
	  new_univdiagonal = univdiagonal_i;
      
	  path->splice5p = false;
	  path->splicetype5 = NO_SPLICE;
	  path->ambig_prob_5 = 0.0;
	  
	} else {
	  order = Sedgesort_order_int(medial_qpositions,npartners);

	  i = 0;
	  while (i < npartners && medial_qpositions[order[i]] < medial_qpos) {
	    i++;
	  }
	  j = i + 1;
	  while (j < npartners && medial_qpositions[order[j]] == medial_qpos) {
	    j++;
	  }
	  assert(j - i == n_medial_partners);
	  
	  path->splice5p = false;
	  path->splicetype5 = NO_SPLICE;
	  path->ambig_prob_5 = 0.0;
	
	  if ((path->qstart_alts =
	       Altsplice_qstart_new(&best_endpoint,/*boundedp*/false,
				    medial_qpos,/*anchor_qpos*/pos3,querylength,
				    /*medial_univdiagonal*/univdiagonal,medial_nmismatches[goodi],
				    medial_probs[goodi],medial_support[goodi],
				    distal_positions,distal_support,distal_univdiags,distal_probs,
				    /*order_medial_qpos*/&(order[i]),/*npartners*/n_medial_partners,
				    path->plusp,/*sense_forward_p*/path->sensedir == SENSE_FORWARD,
				    /*innerp*/true,pathpool,vectorpool,/*sort_bydistal_p*/false)) == NULL) {
	    Intlist_head_set(path->endpoints,medial_qpos);
	    Intlist_head_set(path->nmismatches,medial_nmismatches[goodi]);
	    Intlist_head_set(path->ref_nmismatches,medial_nmismatches[goodi]);

	  } else {
	    Intlist_head_set(path->endpoints,best_endpoint);
	    Intlist_head_set(path->nmismatches,-1);
	    Intlist_head_set(path->ref_nmismatches,-1);
	  }

	  FREE(order);
	}
      }

      /* Vectorpool_free_univcoordvector(distal_positions,vectorpool); */
      /* Vectorpool_free_doublevector(distal_probs,vectorpool); */
      Univdiag_free_array(&distal_univdiags);
    }
  }
  debug3(Path_print(path));


  /* If we have added a new univdiagonal (splice), see if there is an indel beyond that */
  /* Code taken from attach_indel_qstart */

  /* Note: new_univdiagonal is not defined with Splice_resolve_qstart */
  new_univdiagonal = Univcoordlist_head(path->univdiagonals);
  debug3(printf("Comparing new_univdiagonal %u with univdiagonal %u\n",new_univdiagonal,univdiagonal));

  if (new_univdiagonal == univdiagonal) { /* Previously checked for new_univdiagonal == 0 */
    /* Skip */
  } else if ((nosplice_trimpos = Intlist_head(path->endpoints)) == 0) {
    /* Skip */
  } else if ((indel_univdiagonal =
	      Spliceends_indel_qstart(nosplice_trimpos,/*univdiagonal*/new_univdiagonal,
				      querylength,path->chroffset,path->chrhigh,path->plusp,genestrand,
				      localdb_nmismatches_allowed,novel_univdiags_alloc,novel_diagonals_alloc,
				      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,stage1,
				      univdiagpool,query_compress,queryptr)) == 0) {
    /* Skip */

  } else if ((low_qstart = Genomebits_trim_qstart(&nmismatches_i,query_compress,
						  genomebits,indel_univdiagonal,querylength,
						  /*pos5*/0,/*pos3*/nosplice_trimpos,
						  path->plusp,genestrand)) >= nosplice_trimpos) {
    debug3(printf("No distal segment before qstart indel: low_qstart %d >= nosplice_trimpos %d\n",
		  low_qstart,nosplice_trimpos));

  } else if (indel_univdiagonal > new_univdiagonal) {
    /* Insertion */
    nindels = indel_univdiagonal - new_univdiagonal;
    qend = Intlist_second_value(path->endpoints);
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal_i*/indel_univdiagonal,/*indels*/+nindels,path->chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/low_qstart,/*pos3*/qend,querylength,
						    indelinfo,path->plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug3(printf("Insertion fails\n"));
      
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - (indel_pos + nindels);
      debug3(printf("(1) supporti %d - %d, supportj %d - (%d + %d)\n",
		    indel_pos,low_qstart,qend,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug3(printf("(1) Not enough support for indel: supporti %d and mismatches %d\n",
		       supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug3(printf("(1) Not enough support for indel: supportj %d and mismatches %d\n",
		       supportj,nmismatches_j));
      } else {
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,low_qstart
					   intlistpool_trace(__FILE__,__LINE__));
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_insertion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,indel_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	debug3(printf("(5) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		       low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }

  } else if (indel_univdiagonal + max_deletionlen >= new_univdiagonal) {
    /* Deletion (or short intron) */
    nindels = new_univdiagonal - indel_univdiagonal;
    qend = Intlist_second_value(path->endpoints);
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/indel_univdiagonal,/*indels*/-nindels,path->chrhigh,
						   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/low_qstart,/*pos3*/qend,querylength,
						   indelinfo,path->plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug3(printf("Deletion or short intron fails\n"));
	  
    } else {
      supporti = indel_pos - low_qstart;
      supportj = qend - indel_pos;
      debug3(printf("(2) supporti %d - %d, supportj %d - %d\n",
		    indel_pos,low_qstart,qend,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug3(printf("(2) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug3(printf("(2) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	assert(nindels >= 0);
	/* deletionpos = (indel_univdiagonal - querylength) + indel_pos; */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,low_qstart
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* For qstart, push j first, then push i */
	Intlist_head_set(path->nmismatches,nmismatches_j);
	Intlist_head_set(path->ref_nmismatches,nmismatches_j);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_i
						 intlistpool_trace(__FILE__,__LINE__));
      
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,indel_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	debug3(printf("(3) Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		      low_qstart,qend,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }
  }


#if 0
  path->nmatches = -1;
Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
#else
  path->nmatches = -1;		/* Calls Path_eval_nmatches */
Path_trim_chrbounds(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev,
  		      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
		    
#endif

  debug3(printf("Result of Path_qstart_resolve: ")); debug3(Path_print(path));

  return;
}


/* Modifies path and does not create a copy */
/* innerp is true by definition (resolve means to determine inner splices) */
void
Path_qend_resolve (int *found_score, int *end_score, T path,
		   Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		   char *queryptr, int querylength,
		   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		   int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		   int localdb_nmismatches_allowed,
		   
		   Stage1_T stage1, Knownsplicing_T knownsplicing, Indelinfo_T indelinfo,
		   Compress_T query_compress, Compress_T query_compress_fwd, Compress_T query_compress_rev, int genestrand,
		   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		   Univdiagpool_T univdiagpool, Vectorpool_T vectorpool) {

  Univcoord_T new_univdiagonal, indel_univdiagonal;
  int nosplice_trimpos, qstart;
  int supporti, supportj;

  Univcoordlist_T univdiagonals;
  Intlist_T nmismatches, ref_nmismatches, endpoints;
  List_T junctions;
  int ignore;

  Univcoord_T *new_univdiagonals, nearest_univdiagonal, univdiagonal_j,
    univdiagonal, middle_univdiagonal;
  Univcoord_T *distal_positions;
  Univdiag_T *distal_univdiags;

  int *order;
  int *medial_nmismatches, *medial_qpositions;
  double *medial_probs, *distal_probs;
  int *medial_support, *distal_support;
  int n_medial_partners, npartners, ndiagonals, goodi, i, j;
  bool sense_forward_p, multiple_splice_qpos_p;
  Junction_T junction;

  Chrpos_T splice_distance;
  int type;
  int cont_trimpos, local_trimpos, trimpos5, trimpos3, qend, pos5;
  int best_endpoint, medial_qpos, splice_qpos;
  int nindels, indel_pos;
  int cont_nmismatches, cont_ref_nmismatches, local_nmismatches, resolve_nmismatches,
    nmismatches_i, nmismatches_j, nmismatches_indel;
  int ref_nmismatches_i, ref_nmismatches_j, ref_nmismatches_indel;
  double best_prob, donor_prob, acceptor_prob;

  int ninserts;			/* Needed for Path_qend_resolve but not Path_qstart_resolve */


  debug3(printf("Entered Path_qend_resolve, sensedir %d, univdiagonal %u..%u and path\n",
		path->sensedir,low_univdiagonal,high_univdiagonal));
  debug3(Path_print(path));

  Path_expect_fwd(path);

  univdiagonal = Univcoordlist_last_value(path->univdiagonals);
  debug3(printf("Original univdiagonal: %u\n",univdiagonal));

  if (low_univdiagonal < univdiagonal) {
    debug3(printf("Changing low_univdiagonal to be %u\n",univdiagonal));
    low_univdiagonal = univdiagonal;
  }

  if (high_univdiagonal <= low_univdiagonal) {
    /* Ends overlap, so no way to resolve */
    debug3(printf("Interval %u..%u is negative, so not resolving\n",
		  low_univdiagonal,high_univdiagonal));

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[path->chrnum] == true) {
    /* No splicing on circular chromosomes */
    debug3(printf("No splicing on circular chromosomes\n"));
#endif

  } else {
    debug3(printf("Interval %u..%u is positive, so resolving\n",
		  low_univdiagonal,high_univdiagonal));

    qend = Intlist_last_value(path->endpoints);
    pos5 = Intlist_penultimate_value(path->endpoints);
    if (path->junctions != NULL) {
      pos5 += Junction_ninserts((Junction_T) List_last_value(path->junctions,NULL));
    }

#if 0
    /* This can lead to poor results */
    /* Compute continuation to querylength for baseline */
    cont_trimpos = querylength;
    cont_nmismatches =
      Genomebits_count_mismatches_substring(&cont_ref_nmismatches,genomebits,genomebits_alt,
					    query_compress,univdiagonal,querylength,
					    pos5,/*pos3*/querylength,path->plusp,genestrand);
#else
    cont_trimpos = Genomebits_trim_qend(&cont_nmismatches,query_compress,genomebits,
					univdiagonal,querylength,
					pos5,/*pos3*/querylength,path->plusp,genestrand);
#endif

    new_univdiagonals =
      Spliceends_qend_resolve(&ndiagonals,/*pos5*/qend,querylength,
			      low_univdiagonal,high_univdiagonal,
			      query_compress,queryptr,path->plusp,genestrand,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,stage1,
			      localdb_nmismatches_allowed);

    if (ndiagonals == 0) {
      debug3(printf("Spliceends_qend_resolve returns nothing\n"));
      local_trimpos = Genomebits_trim_qend(&local_nmismatches,query_compress,genomebits,
					   univdiagonal,querylength,
					   pos5,/*pos3*/querylength,path->plusp,genestrand);
      debug3(printf("(5) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
      if (local_trimpos <= pos5) {
	debug3(printf("Need to trim exon completely\n"));
	Path_trim_qend_n(/*noutside*/Intlist_last_value(path->endpoints) - local_trimpos,path,
			 query_compress_fwd,query_compress_rev,
			 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
	
      } else {
	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,local_trimpos);
	Intlist_head_set(path->nmismatches,local_nmismatches);
	Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	Path_reverse(path,/*expect_fwd_p*/true);
      
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
      }

    } else if ((nearest_univdiagonal = new_univdiagonals[0]) == univdiagonal) {
      /* Continuation of univdiagonal */
      debug3(printf("Continuation of univdiagonal now with %d nmismatches\n",local_nmismatches));
      local_trimpos = Genomebits_trim_qend(&local_nmismatches,query_compress,genomebits,
					   nearest_univdiagonal,querylength,
					   pos5,/*pos3*/querylength,path->plusp,genestrand);

      debug3(printf("(6) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
      if (local_trimpos <= pos5) {
	debug3(printf("Need to trim exon completely\n"));
	Path_trim_qend_n(/*noutside*/Intlist_last_value(path->endpoints) - local_trimpos,path,
			 query_compress_fwd,query_compress_rev,
			 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

      } else {
	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,local_trimpos);
#if 0
	if (Intlist_head(path->nmismatches) >= 0) {
	  Intlist_head_incr(path->nmismatches,local_nmismatches);
	  Intlist_head_incr(path->ref_nmismatches,local_nmismatches);
	}
#else
	Intlist_head_set(path->nmismatches,local_nmismatches);
	Intlist_head_set(path->ref_nmismatches,local_nmismatches);
#endif
	Path_reverse(path,/*expect_fwd_p*/true);
	
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
      }
      
    } else if (univdiagonal > nearest_univdiagonal) {
      /* Insertion.  Skip */
      
    } else if (univdiagonal + max_deletionlen >= nearest_univdiagonal) {
      /* Deletion (or short intron) */
      nindels = nearest_univdiagonal - univdiagonal;
      trimpos3 = Genomebits_trim_qend(&resolve_nmismatches,query_compress,genomebits,
				      /*univdiagonal*/nearest_univdiagonal,querylength,
				      pos5,/*pos3*/querylength,path->plusp,genestrand);

      if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						     &ref_nmismatches_i,&ref_nmismatches_j,
						     /*univdiagonal_i*/univdiagonal,/*indels*/-nindels,
						     Path_chrhigh(path),
						     /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						     /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						     /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						     pos5,/*pos3*/trimpos3,querylength,
						     stage1->indelinfo,path->plusp,genestrand,
						     /*want_lowest_coordinate_p*/true)) <= 0) {
	debug3(printf("Deletion or short intron fails\n"));
	
      } else if (trimpos3 - nmismatches_i - nmismatches_j <= cont_trimpos - cont_nmismatches) {
	debug3(printf("Deletion or short intron yields fewer matches than continuation: %d - %d - %d  vs %d - %d\n",
		      trimpos3,nmismatches_i,nmismatches_j,cont_trimpos,cont_nmismatches));

      } else {
	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos3 /* was querylength */
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* deletionpos = univdiagonal + indel_pos; */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = nearest_univdiagonal;
	
	Path_reverse(path,/*expect_fwd_p*/true);
	
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
      }
      
    } else if (splicingp == false) {
      /* Unable to try splicing */

    } else if (ndiagonals == 1) {

#if 0
      new_univdiagonal = 0;
      /* Splice_nomiddle call type 1 */
      if ((splice_qpos = Splice_nomiddle(&type,&trimpos5,&trimpos3,&supporti,&supportj,
					 &nindels,&indel_pos,&nmismatches_i,&nmismatches_j,&nmismatches_indel,
					 &ref_nmismatches_i,&ref_nmismatches_j,&ref_nmismatches_indel,
					 &donor_prob,&acceptor_prob,
					 /*univdiagonal_i*/univdiagonal,/*univdiagonal_j*/nearest_univdiagonal,
					 query_compress,path->plusp,Path_chroffset(path),Path_chrhigh(path),
					 pos5,/*pos3*/querylength,querylength,
					 stage1->indelinfo,stage1->spliceinfo,knownsplicing,intlistpool,
					 /*sense_forward_p*/(path->sensedir == SENSE_FORWARD) ? true : false,
					 genestrand,/*trim5p*/false,/*trim3p*/true,/*innerp*/true,
					 /*allow_indels_p*/true)) < 0) {

	debug3(printf("Single splice finds no solution\n"));
	local_trimpos = Genomebits_trim_qend(&local_nmismatches,query_compress,
					     /*bits*/genomebits,univdiagonal,querylength,
					     pos5,/*pos3*/querylength,path->plusp,genestrand);
	debug3(printf("(7) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos <= pos5) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qend_n(/*noutside*/Intlist_last_value(path->endpoints) - local_trimpos,path,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Path_reverse(path,/*expect_fwd_p*/false);
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	  Path_reverse(path,/*expect_fwd_p*/true);
	
	  path->splice3p = false;
	  path->splicetype3 = NO_SPLICE;
	  path->ambig_prob_3 = 0.0;
	}
	
      } else if (trimpos3 - nmismatches_i - nmismatches_j <= cont_trimpos - cont_nmismatches) {
	debug3(printf("Single splice yields fewer matches than continuation: %d - %d - %d  vs %d - %d\n",
		      trimpos3,nmismatches_i,nmismatches_j,cont_trimpos,cont_nmismatches));

      } else if (nindels == 0) {
	/* Single splice succeeds */
	debug3(printf("Single splice succeeds with trimpos3 %d\n",trimpos3));
	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,splice_qpos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos3 /* was querylength */
					   intlistpool_trace(__FILE__,__LINE__));
	
	splice_distance = nearest_univdiagonal - univdiagonal;
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = nearest_univdiagonal;
	Path_reverse(path,/*expect_fwd_p*/true);
	
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
	
	Transcript_list_gc(&path->transcripts,listpool,transcriptpool);
	Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool);
	
      } else if (indel_pos < splice_qpos) {
	/* Single spliceindel succeeds */
	/* Push indel (based on left) then splice.  indel is medial, splice is distal. */
	debug13(printf("indel_pos %d < splice_qpos %d\n",indel_pos,splice_qpos));
	
	middle_univdiagonal = univdiagonal + nindels; /* nindels = middle_univdiagonal - univdiagonal */
	splice_distance = nearest_univdiagonal - middle_univdiagonal;
	
	debug3(printf("Splice_nomiddle (qend, sense): %d indels at %d, nmismatches %d, then splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
		      nindels,indel_pos,nmismatches_indel,
		      qstart,/*high_qend*/querylength,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob));

	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,splice_qpos
					   intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos3   /* was high_qend */
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* Indel first */
	if (nindels < 0) {
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_insertion(-nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	} else {
	  /* deletionpos = left + indel_pos; -- qend; do not add splice_distance if indel first */
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_deletion(nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	}
	  
	/* Splice second */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_indel
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_indel
						 intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,middle_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	Path_reverse(path,/*expect_fwd_p*/true);
	
      } else {
	/* Single spliceindel succeeds */
	/* Push splice then indel (based on nearest_univdiagonal).  splice is medial, indel is distal. */
	debug3(printf("indel_pos %d >= splice_qpos %d\n",indel_pos,splice_qpos));
	
	middle_univdiagonal = nearest_univdiagonal - nindels; /* nindels = nearest_univdiagonal - middle_univdiagonal */
	splice_distance = middle_univdiagonal - univdiagonal;
	
	debug3(printf("Splice_nomiddle (qend, sense): splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f, then %d indels at %d\n",
		      qstart,/*high_qend*/querylength,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob,
		      nindels,indel_pos));

	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,splice_qpos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,indel_pos
					   intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,trimpos3   /* was high_qend */
					   intlistpool_trace(__FILE__,__LINE__));
	
	/* Splice first */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* Indel second */
	if (nindels < 0) {
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_insertion(-nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	} else {
	  /* deletionpos = left + splice_distance + indel_pos; -- qend; add splice_distance if indel second */
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_deletion(nindels,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	}
	
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,ref_nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_indel
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_indel
						 intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,middle_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,nearest_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	Path_reverse(path,/*expect_fwd_p*/true);
      }

#else
	
      if (Splice_resolve_qend(&univdiagonals,&nmismatches,&ref_nmismatches,&endpoints,&junctions,
			      /*univdiagonal_i*/univdiagonal,/*univdiagonal_j*/nearest_univdiagonal,
			      stage1,query_compress,queryptr,
			      path->plusp,Path_chroffset(path),Path_chrhigh(path),
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      localdb_tally,localdb_overflowp,max_nblocks,localdb,
			      localdb_nmismatches_allowed,
			      pos5,/*pos3*/querylength,querylength,cont_trimpos,cont_nmismatches,
			      indelinfo,stage1->spliceinfo,knownsplicing,univcoordlistpool,
			      intlistpool,listpool,pathpool,
			      path->sensedir,genestrand,/*trim5p*/false,/*trim3p*/true,/*innerp*/true,
			      /*allow_indels_p*/true,/*find_middle_p*/false,/*middle_rank*/0) == 0) {

	debug3(printf("Single splice finds no solution\n"));
	local_trimpos = Genomebits_trim_qend(&local_nmismatches,query_compress,
					     /*bits*/genomebits,univdiagonal,querylength,
					     pos5,/*pos3*/querylength,path->plusp,genestrand);
	debug3(printf("(7) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos <= pos5) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qend_n(/*noutside*/Intlist_last_value(path->endpoints) - local_trimpos,path,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Path_reverse(path,/*expect_fwd_p*/false);
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	  Path_reverse(path,/*expect_fwd_p*/true);
	
	  path->splice3p = false;
	  path->splicetype3 = NO_SPLICE;
	  path->ambig_prob_3 = 0.0;
	}
	
#if 0
      } else if (trimpos3 - nmismatches_i - nmismatches_j <= cont_trimpos - cont_nmismatches) {
	debug3(printf("Single splice yields fewer matches than continuation: %d - %d - %d  vs %d - %d\n",
		      trimpos3,nmismatches_i,nmismatches_j,cont_trimpos,cont_nmismatches));
#endif

      } else {
	/* No need to copy path */
	debug3(printf("Adding results from Splice_resolve_qend\n"));

	Path_reverse(path,/*expect_fwd_p*/false);
	path->univdiagonals = Univcoordlist_append(univdiagonals,path->univdiagonals);
	path->nmismatches = Intlistpool_pop(path->nmismatches,intlistpool,&ignore
					       intlistpool_trace(__FILE__,__LINE__));
	path->nmismatches = Intlist_append(nmismatches,path->nmismatches);
	path->ref_nmismatches = Intlistpool_pop(path->ref_nmismatches,intlistpool,&ignore
						   intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlist_append(ref_nmismatches,path->ref_nmismatches);
	path->endpoints = Intlistpool_pop(path->endpoints,intlistpool,&ignore
					     intlistpool_trace(__FILE__,__LINE__));
	path->endpoints = Intlist_append(endpoints,path->endpoints);
	path->junctions = List_append(junctions,path->junctions);
	Path_reverse(path,/*expect_fwd_p*/true);
	
	debug3(printf("AFTER call to Splice_resolve_qend\n"));
	debug3(Path_print(path));
      }
#endif
      
    } else {
      /* Multiple splices */
      debug3(printf("(2) Multiple diagonals: %d\n",ndiagonals));
      npartners = 0;

      medial_qpositions = Vectorpool_new_intvector(vectorpool,ndiagonals+1); /* Need +1 for Sedgesort_order */
      distal_positions = Vectorpool_new_univcoordvector(vectorpool,ndiagonals);

      distal_support = Vectorpool_new_intvector(vectorpool,ndiagonals);
      distal_univdiags = Univdiag_new_array(ndiagonals);

      medial_probs = Vectorpool_new_doublevector(vectorpool,ndiagonals);
      distal_probs = Vectorpool_new_doublevector(vectorpool,ndiagonals);
      medial_support = Vectorpool_new_intvector(vectorpool,ndiagonals);
      medial_nmismatches = Vectorpool_new_intvector(vectorpool,ndiagonals);

      sense_forward_p = (path->sensedir == SENSE_FORWARD) ? true : false;

      multiple_splice_qpos_p = false;
      for (i = 0; i < ndiagonals; i++) {
	univdiagonal_j = new_univdiagonals[i];

#if 0
	debug3(printf("Calling Splice_nomiddle with %u and %u\n",univdiagonal,univdiagonal_j));
	/* Splice_nomiddle call type 2 */
	if ((splice_qpos = Splice_nomiddle(&type,&trimpos5,&trimpos3,&supporti,&supportj,
					   &nindels,&indel_pos,&nmismatches_i,&nmismatches_j,&nmismatches_indel,
					   &ref_nmismatches_i,&ref_nmismatches_j,&ref_nmismatches_indel,
					   &donor_prob,&acceptor_prob,
					   /*univdiagonal_i*/univdiagonal,univdiagonal_j,
					   query_compress,path->plusp,Path_chroffset(path),Path_chrhigh(path),
					   pos5,/*pos3*/querylength,querylength,
					   stage1->indelinfo,stage1->spliceinfo,knownsplicing,intlistpool,
					   sense_forward_p,genestrand,/*trimp5*/false,/*trim3p*/true,/*innerp*/true,
					   /*allow_indels_p*/false)) < 0) {
	  /* Skip */

	} else if (npartners == 0) {
	  medial_qpos = medial_qpositions[0] = splice_qpos;
	  medial_support[0] = supporti;
	  medial_nmismatches[0] = nmismatches_i;
	  distal_positions[0] = univdiagonal_j - querylength + splice_qpos;
	  
	  distal_support[0] = supportj;
	  distal_univdiags[0]->univdiagonal = univdiagonal_j;
	  distal_univdiags[0]->qstart = splice_qpos;
	  distal_univdiags[0]->qend = trimpos3;
	  distal_univdiags[0]->nmismatches = nmismatches_j;

	  if (path->plusp == sense_forward_p) { /* was sense_forward_p == true */
	    medial_probs[0] = donor_prob;
	    distal_probs[0] = acceptor_prob;
	  } else {
	    medial_probs[0] = acceptor_prob;
	    distal_probs[0] = donor_prob;
	  }
	  npartners++;
	  
	} else {
	  if (splice_qpos != medial_qpos) {
	    multiple_splice_qpos_p = true;
	  }

	  medial_qpositions[npartners] = splice_qpos;
	  medial_support[npartners] = supporti;
	  medial_nmismatches[npartners] = nmismatches_i;
	  distal_positions[npartners] = univdiagonal_j - querylength + splice_qpos;

	  distal_support[npartners] = supportj;
	  distal_univdiags[npartners]->univdiagonal = univdiagonal_j;
	  distal_univdiags[npartners]->qstart = splice_qpos;
	  distal_univdiags[npartners]->qend = trimpos3;
	  distal_univdiags[npartners]->nmismatches = nmismatches_j;

	  if (path->plusp == sense_forward_p) {
	    medial_probs[npartners] = donor_prob;
	    distal_probs[npartners] = acceptor_prob;
	  } else {
	    medial_probs[npartners] = acceptor_prob;
	    distal_probs[npartners] = donor_prob;
	  }
	  npartners++;
	}
	
#else
	/* allow_indels_p is false because Altsplice_T cannot handle indels */
	debug3(printf("Calling Splice_resolve_qend with %u and %u\n",univdiagonal,univdiagonal_j));
	if (Splice_resolve_qend(&univdiagonals,&nmismatches,&ref_nmismatches,&endpoints,&junctions,
				/*univdiagonal_i*/univdiagonal,univdiagonal_j,
				stage1,query_compress,queryptr,
				path->plusp,Path_chroffset(path),Path_chrhigh(path),

				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,localdb,
				localdb_nmismatches_allowed,

				pos5,/*pos3*/querylength,querylength,cont_trimpos,cont_nmismatches,
				stage1->indelinfo,stage1->spliceinfo,knownsplicing,
				univcoordlistpool,intlistpool,listpool,pathpool,
				path->sensedir,genestrand,/*trim5p*/false,/*trim3p*/true,/*innerp*/true,
				/*allow_indels_p*/false,/*find_middle_p*/false,/*middle_rank*/0) == 0) {
	  /* Skip */

	} else {
	  splice_qpos = Intlist_second_value(endpoints);
	  if (npartners == 0) {
	    /* Do not set multiple_splice_qpos_p */
	    medial_qpos = splice_qpos;
	  } else if (splice_qpos != medial_qpos) {
	    multiple_splice_qpos_p = true;
	  }

	  debug3(printf("endpoints %s => supporti %d, supportj %d\n",
			Intlist_to_string(endpoints),Intlist_second_value(endpoints) - pos5,
			Intlist_head(endpoints) - Intlist_second_value(endpoints)));

	  medial_qpositions[npartners] = splice_qpos;
	  medial_support[npartners] = /*supporti*/ Intlist_second_value(endpoints) - pos5;
	  medial_nmismatches[npartners] = /*nmismatches_i*/ Intlist_second_value(nmismatches);
	  distal_positions[npartners] = univdiagonal_j - querylength + splice_qpos;

	  distal_support[npartners] = /*supportj*/ Intlist_head(endpoints) - Intlist_second_value(endpoints);
	  distal_univdiags[npartners]->univdiagonal = univdiagonal_j;
	  distal_univdiags[npartners]->qstart = splice_qpos;
	  distal_univdiags[npartners]->qend = /*trimpos3*/ Intlist_head(endpoints);
	  distal_univdiags[npartners]->nmismatches = /*nmismatches_j*/ Intlist_head(nmismatches);

	  junction = (Junction_T) List_head(junctions);
	  /* good_outer_accept_p = junction->outer_accept_p; -- Used if npartners == 1 */
	  /* assert(Junction_type(junction) == SPLICE_JUNCTION); -- Could have been converted to a deletion */
	  debug3(Junction_print(junction));

	  if (path->plusp == sense_forward_p) {
	    medial_probs[npartners] = junction->donor_prob;
	    distal_probs[npartners] = junction->acceptor_prob;
	  } else {
	    medial_probs[npartners] = junction->acceptor_prob;
	    distal_probs[npartners] = junction->donor_prob;
	  }

	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__));
	  Junction_list_gc(&junctions,listpool,pathpool);

	  debug3(printf("%u %u => %d, %f, %f\n",univdiagonal_j,univdiagonal_j - univdiagonal,
			splice_qpos,junction->donor_prob,junction->acceptor_prob));
	  npartners++;
	}
#endif
      }

      if (npartners == 0) {
	debug3(printf("(2) Multiple splices find no solution\n"));
	local_trimpos = Genomebits_trim_qend(&local_nmismatches,query_compress,
					     /*bits*/genomebits,univdiagonal,querylength,
					     pos5,/*pos3*/querylength,path->plusp,genestrand);
	debug3(printf("(8) Have local trimpos %d with %d nmismatches\n",local_trimpos,local_nmismatches));
	if (local_trimpos <= pos5) {
	  debug3(printf("Need to trim exon completely\n"));
	  Path_trim_qend_n(/*noutside*/Intlist_last_value(path->endpoints) - local_trimpos,path,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);

	} else {
	  Path_reverse(path,/*expect_fwd_p*/false);
	  Intlist_head_set(path->endpoints,local_trimpos);
	  Intlist_head_set(path->nmismatches,local_nmismatches);
	  Intlist_head_set(path->ref_nmismatches,local_nmismatches);
	  Path_reverse(path,/*expect_fwd_p*/true);
	
	  path->splice3p = false;
	  path->splicetype3 = NO_SPLICE;
	  path->ambig_prob_3 = 0.0;
	}

      } else if (npartners == 1) {
	/* Single valid splice */
	debug3(printf("(2) Multiple splices yield one valid solution\n"));

	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,medial_qpositions[0]);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,distal_univdiags[0]->qend /* was querylength */
					   intlistpool_trace(__FILE__,__LINE__));

	univdiagonal_j = distal_positions[0] - medial_qpositions[0] + querylength; /* univdiagonal */
	splice_distance = univdiagonal_j - univdiagonal;
	if (path->plusp == sense_forward_p) { /* was sense_forward_p == true */
	  donor_prob = medial_probs[0];
	  acceptor_prob = distal_probs[0];
	} else {
	  donor_prob = distal_probs[0];
	  acceptor_prob = medial_probs[0];
	}
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_splice(splice_distance,path->sensedir,
								     donor_prob,acceptor_prob,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,medial_nmismatches[0]);
	Intlist_head_set(path->ref_nmismatches,medial_nmismatches[0]);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,distal_univdiags[0]->nmismatches
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,distal_univdiags[0]->nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
	
	debug3(printf("Pushing univdiagonal %u\n",univdiagonal_j));
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,univdiagonal_j
						     univcoordlistpool_trace(__FILE__,__LINE__));
	new_univdiagonal = univdiagonal_j;
	Path_reverse(path,/*expect_fwd_p*/true);
	
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;

      } else if (multiple_splice_qpos_p == false &&
		 (path->qend_alts =
		  Altsplice_qend_new(&best_endpoint,/*boundedp*/false,
				     medial_qpos,/*anchor_qpos*/pos5,querylength,
				     /*medial_univdiagonal*/univdiagonal,medial_nmismatches[0],
				     medial_probs[0],medial_support[0],distal_positions,
				     distal_support,distal_univdiags,distal_probs,
				     /*order_medial_qpos*/NULL,npartners,
				     path->plusp,/*sense_forward_p*/path->sensedir == SENSE_FORWARD,
				     /*innerp*/true,pathpool,vectorpool,/*sort_bydistal_p*/false)) != NULL) {
	debug3(printf("(2) %d partners with same medial qpos => altsplice\n",npartners));
	Path_reverse(path,/*expect_fwd_p*/false);

	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;

	Intlist_head_set(path->endpoints,best_endpoint);
	Intlist_head_set(path->nmismatches,-1);
	Intlist_head_set(path->ref_nmismatches,-1);

	Path_reverse(path,/*expect_fwd_p*/true);

      } else {
	/* Determine the best medial_qpos, based on the total prob */
	best_prob = 0.0;
	for (i = 0; i < npartners; i++) {
	  if (medial_probs[i] + distal_probs[i] > best_prob) {
	    best_prob = medial_probs[i] + distal_probs[i];
	    medial_qpos = medial_qpositions[i];
	  }
	}
	debug3(printf("best medial qpos is %d with prob %f\n",medial_qpos,best_prob));

	n_medial_partners = 0;
	for (i = 0; i < npartners; i++) {
#if 0
	  /* Need to use distal_support when a distal indel is found by Spliceends */
	  if (univdiagonal < distal_univdiags[i]->univdiagonal &&
	      Splice_accept_p(splice_qpos,querylength,medial_probs[i],distal_probs[i],
			      /*univdiagonal_i*/univdiagonal,
			      /*univdiagonal_j*/distal_univdiags[i]->univdiagonal,
			      medial_support[i],medial_nmismatches[i],
			      distal_support[i],distal_univdiags[i]->nmismatches,
			      path->plusp,sense_forward_p,/*innerp*/true) == true) {
	    /* Skip */
	  } else if (medial_qpositions[i] == medial_qpos) {
	    n_good_partners += 1;
	    goodi = i;
	  }
#else
	  if (medial_qpositions[i] == medial_qpos) {
	    n_medial_partners += 1;
	    goodi = i;		/* Used in case n_medial_partners == 1 */
	  }
#endif
	}

	debug3(printf("Determined that medial_qpos %d is best, with %d partners\n",
		      medial_qpos,n_medial_partners));
	
	if (n_medial_partners == 0) {
	  /* All good distal probs are paired with a low medial prob.  Leave path as is. */

	} else if (n_medial_partners == 1) {
	  Path_reverse(path,/*expect_fwd_p*/false);
	  Intlist_head_set(path->endpoints,medial_qpositions[goodi]);
	  path->endpoints = Intlistpool_push(path->endpoints,intlistpool,distal_univdiags[goodi]->qend /* was querylength */
					     intlistpool_trace(__FILE__,__LINE__));

	  univdiagonal_j = distal_positions[goodi] - medial_qpositions[goodi] + querylength; /* univdiagonal */
	  splice_distance = univdiagonal_j - univdiagonal;
	  if (path->plusp == sense_forward_p) { /* was sense_forward_p == true */
	    donor_prob = medial_probs[goodi];
	    acceptor_prob = distal_probs[goodi];
	  } else {
	    donor_prob = distal_probs[goodi];
	    acceptor_prob = medial_probs[goodi];
	  }
	  path->junctions = Listpool_push(path->junctions,listpool,
					  (void *) Junction_new_splice(splice_distance,path->sensedir,
								       donor_prob,acceptor_prob,pathpool)
					  listpool_trace(__FILE__,__LINE__));
	
	  /* For qend, push i first, then push j */
	  Intlist_head_set(path->nmismatches,medial_nmismatches[goodi]);
	  Intlist_head_set(path->ref_nmismatches,medial_nmismatches[goodi]);
	  path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,distal_univdiags[goodi]->nmismatches
					       intlistpool_trace(__FILE__,__LINE__));
	  path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,distal_univdiags[goodi]->nmismatches
						   intlistpool_trace(__FILE__,__LINE__));
	
	  debug3(printf("Pushing univdiagonal %u\n",univdiagonal_j));
	  path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,univdiagonal_j
						       univcoordlistpool_trace(__FILE__,__LINE__));
	  new_univdiagonal = univdiagonal_j;
	  Path_reverse(path,/*expect_fwd_p*/true);
	
	  path->splice3p = false;
	  path->splicetype3 = NO_SPLICE;
	  path->ambig_prob_3 = 0.0;

	} else {
	  order = Sedgesort_order_int(medial_qpositions,npartners);

	  i = 0;
	  while (i < npartners && medial_qpositions[order[i]] < medial_qpos) {
	    i++;
	  }
	  j = i + 1;
	  while (j < npartners && medial_qpositions[order[j]] == medial_qpos) {
	    j++;
	  }
	  assert(j - i == n_medial_partners);

	  Path_reverse(path,/*expect_fwd_p*/false);

	  path->splice3p = false;
	  path->splicetype3 = NO_SPLICE;
	  path->ambig_prob_3 = 0.0;

	  if ((path->qend_alts =
	       Altsplice_qend_new(&best_endpoint,/*boundedp*/false,
				  medial_qpos,/*anchor_qpos*/pos5,querylength,
				  /*medial_univdiagonal*/univdiagonal,medial_nmismatches[goodi],
				  medial_probs[goodi],medial_support[goodi],
				  distal_positions,distal_support,distal_univdiags,distal_probs,
				  /*order_medial_qpos*/&(order[i]),/*npartners*/n_medial_partners,
				  path->plusp,/*sense_forward_p*/path->sensedir == SENSE_FORWARD,
				  /*innerp*/true,pathpool,vectorpool,/*sort_bydistal_p*/false)) == NULL) {
	    Intlist_head_set(path->endpoints,medial_qpos);
	    Intlist_head_set(path->nmismatches,medial_nmismatches[goodi]);
	    Intlist_head_set(path->ref_nmismatches,medial_nmismatches[goodi]);

	  } else {
	    Intlist_head_set(path->endpoints,best_endpoint);
	    Intlist_head_set(path->nmismatches,-1);
	    Intlist_head_set(path->ref_nmismatches,-1);
	  }
	  
	  Path_reverse(path,/*expect_fwd_p*/true);
	  FREE(order);
	}
      }

      /* Vectorpool_free_univcoordvector(distal_positions,vectorpool); */
      /* Vectorpool_free_doublevector(distal_probs,vectorpool); */
      Univdiag_free_array(&distal_univdiags);
    }
  }
  debug3(Path_print(path));


  /* If we have added a new univdiagonal (splice), see if there is an indel beyond that */
  /* Code taken from attach_indel_qend */

  /* Note: new_univdiagonal is not defined with Splice_resolve_qend */
  new_univdiagonal = Univcoordlist_last_value(path->univdiagonals);
  debug3(printf("Comparing new_univdiagonal %u with univdiagonal %u\n",new_univdiagonal,univdiagonal));

  if (new_univdiagonal == univdiagonal) { /* Previously checked for new_univdiagonal == 0 */
    /* Skip: No extension */
  } else if ((nosplice_trimpos = Intlist_last_value(path->endpoints)) == querylength) {
    /* Skip */
  } else if ((indel_univdiagonal =
	      Spliceends_indel_qend(nosplice_trimpos,/*univdiagonal*/new_univdiagonal,
				    querylength,path->chroffset,path->chrhigh,path->plusp,genestrand,
				    localdb_nmismatches_allowed,novel_univdiags_alloc,novel_diagonals_alloc,
				    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,stage1,
				    univdiagpool,query_compress,queryptr)) == 0) {
    /* Skip */

  } else if (indel_univdiagonal < new_univdiagonal) {
    /* Insertion */
    nindels = new_univdiagonal - indel_univdiagonal;
    if (path->junctions == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_last_value(path->junctions,NULL));
    }
    qstart = Intlist_penultimate_value(path->endpoints) + ninserts;

    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal_i*/new_univdiagonal,/*indels*/+nindels,path->chrhigh,
						    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						    /*pos5*/qstart,/*pos3*/querylength,querylength,
						    indelinfo,path->plusp,genestrand,
						    /*want_lowest_coordinate_p*/true)) <= 0) {
      debug3(printf("Insertion fails\n"));

    } else {
      supporti = indel_pos - qstart;
      supportj = /*high_qend*/querylength - (indel_pos + nindels);
      debug3(printf("(5) supporti %d - %d, supportj %d - (%d + %d)\n",
		    indel_pos,qstart,/*high_qend*/querylength,indel_pos,nindels));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug3(printf("(5) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug3(printf("(5) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	Path_reverse(path,/*expect_fwd_p*/false);
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,/*high_qend*/querylength
					   intlistpool_trace(__FILE__,__LINE__));
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_insertion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,indel_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	Path_reverse(path,/*expect_fwd_p*/true);

	debug3(printf("(6) Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		      qstart,/*high_qend*/querylength,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }

  } else if (indel_univdiagonal <= new_univdiagonal + max_deletionlen) {
    /* Deletion (or short intron) */
    nindels = indel_univdiagonal - new_univdiagonal;
    if (path->junctions == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_last_value(path->junctions,NULL));
    }
    qstart = Intlist_penultimate_value(path->endpoints) + ninserts;

    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/new_univdiagonal,/*indels*/-nindels,path->chrhigh,
						   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
						   /*pos5*/qstart,/*pos3*/querylength,querylength,
						   indelinfo,path->plusp,genestrand,
						   /*want_lowest_coordinate_p*/true)) <= 0) {
      debug3(printf("Deletion or short intron fails\n"));
      
    } else {
      supporti = indel_pos - qstart;
      supportj = /*high_qend*/querylength - indel_pos;
      debug3(printf("(6) supporti %d - %d, supportj %d - %d\n",
		    indel_pos,qstart,/*high_qend*/querylength,indel_pos));
      if (sufficient_indel_support_p(supporti,nmismatches_i) == false) {
	debug3(printf("(6) Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (sufficient_indel_support_p(supportj,nmismatches_j) == false) {
	debug3(printf("(6) Not enough support for indel: supportj %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	assert(nindels >= 0);
	Path_reverse(path,/*expect_fwd_p*/false);

	/* deletionpos = (new_univdiagonal - querylength) + indel_pos; */
	path->junctions = Listpool_push(path->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
	Intlist_head_set(path->endpoints,indel_pos);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,/*high_qend*/querylength
					      intlistpool_trace(__FILE__,__LINE__));
	
	/* For qend, push i first, then push j */
	Intlist_head_set(path->nmismatches,nmismatches_i);
	Intlist_head_set(path->ref_nmismatches,nmismatches_i);
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,nmismatches_j
						 intlistpool_trace(__FILE__,__LINE__));
	
	path->univdiagonals = Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,indel_univdiagonal
						     univcoordlistpool_trace(__FILE__,__LINE__));
	Path_reverse(path,/*expect_fwd_p*/true);

	debug3(printf("(4) Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		      qstart,/*high_qend*/querylength,indel_pos,nindels,nmismatches_i,nmismatches_j));
      }
    }
  }

#if 0
  path->nmatches = -1;
  Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
#else
  path->nmatches = -1;		/* Calls Path_eval_nmatches */
  Path_trim_chrbounds(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev,
		      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool);
#endif

  debug3(printf("Result of Path_qend_resolve: ")); debug3(Path_print(path));

  return;
}


/* Solves unsolved junctions */
/* Called by Trpath_convert procedures */
/* Set innerp to be true (looser criteria), because we are allowing for introns specified by known transcripts */

#if 0
/* Uses Splice_nomiddle */
T
Path_solve_junctions (int *found_score, int *end_score, T this, int sensedir, int genestrand,
		      Compress_T query_compress,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      Shortread_T queryseq, char *queryptr, int querylength,

		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Localdb_T localdb, int localdb_nmismatches_allowed,

		      Stage1_T stage1, Knownsplicing_T knownsplicing,
		      Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  
  int nmismatches_i, nmismatches_j, nmismatches_indel,
    ref_nmismatches_i, ref_nmismatches_j, ref_nmismatches_indel;
  int nindels;
  double donor_prob, acceptor_prob;

  int type;
  int trimpos5, trimpos3, supporti, supportj;
  int splice_qpos, indel_qpos, prev_ninserts, qstart, qend;
  Univcoord_T univdiagonal_i, univdiagonal_j, univdiagonal_k,
    segmenti_left, segmentj_left;
  Chrpos_T splice_distance;

  int ignore;
  Univcoordlist_T univdiagonals;
  Intlist_T nmismatches, ref_nmismatches, endpoints;
  List_T junctions;

  Junction_T junction;
  List_T j;
  Intlist_T q, r, s;

  List_T jnext;	/* junctions */
  Intlist_T qmid, qnext; /* endpoints */
  Intlist_T rnext, snext; /* ref_nmismatches, nmismatches */
  Univcoordlist_T p, pnext;	/* univdiagonals */

  bool all_solved_p = true;
  bool sense_forward_p = (this->sensedir == SENSE_FORWARD) ? true : false;

  List_T invalid_transcripts, t;
  Transcript_T transcript;
  bool validp;


  debug14(printf("Entered Path_solve_junctions with path\n"));
  debug14(Path_print(this));

  prev_ninserts = 0;
  j = this->junctions;
  p = this->univdiagonals;
  q = this->endpoints;
  r = this->ref_nmismatches;
  s = this->nmismatches;

  while (j != NULL) {
#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
    printf("Endpoints: %s\n",Intlist_to_string(q));
    printf("Univdiagonals: %s\n",Univcoordlist_to_string(p));
    printf("Mismatches: %s\n",Intlist_to_string(s));
    printf("Junctions: ");
    Junction_print_list(j);
    printf("\n");
#endif

    assert(Univcoordlist_length(/*univdiagonals*/p) == Intlist_length(/*endpoints*/q) - 1);
    assert(Intlist_length(/*nmismatches*/s) == Intlist_length(/*endpoints*/q) - 1);
    assert(Intlist_length(/*ref_nmismatches*/r) == Intlist_length(/*endpoints*/q) - 1);
    assert(List_length(/*junctions*/j) == Intlist_length(/*endpoints*/q) - 2);

    junction = (Junction_T) List_head(j);
    debug14(Junction_print(junction));

    if (
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	Junction_type(junction) != UNSOLVED_JUNCTION
#else
	junction != JUNCTION_UNSOLVED
#endif
	) {

      prev_ninserts = Junction_ninserts(junction);
      j = List_next(j);
      p = Univcoordlist_next(p);
      q = Intlist_next(q);
      r = Intlist_next(r);
      s = Intlist_next(s);

    } else {
      debug14(printf("Solving unsolved junction\n"));

      univdiagonal_i = Univcoordlist_head(p);
      univdiagonal_j = Univcoordlist_head(Univcoordlist_next(p));
      segmenti_left = univdiagonal_i - querylength;
      segmentj_left = univdiagonal_j - querylength;

      qstart = Intlist_head(q) + prev_ninserts;	/* The endpoint before the junction */
      qend = Intlist_head(Intlist_next(Intlist_next(q))); /* The endpoint after the junction */

      if (univdiagonal_i > univdiagonal_j + max_insertionlen) {
	/* Impossible */

      } else if (univdiagonal_i > univdiagonal_j) {
	nindels = univdiagonal_i - univdiagonal_j;
	if ((indel_qpos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
							&ref_nmismatches_i,&ref_nmismatches_j,
							univdiagonal_i,/*indels*/+nindels,this->chrhigh,
							/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							/*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							/*pos5*/qstart,/*pos3*/qend,querylength,
							 stage1->indelinfo,this->plusp,genestrand,
							/*want_lowest_coordinate_p*/true)) <= 0) {
	  debug14(printf("Insertion fails\n"));

	  /* Should advance pointers, if this code were used */

	} else {
	  /* No need to change univdiagonals */
	  p = Univcoordlist_next(p);

	  Intlist_head_set(q->rest,indel_qpos);
	  q = Intlist_next(q);

#ifdef ALLOCATE_UNSOLVED_JUNCTION
	  junction = (Junction_T) List_head(j);
	  Junction_free(&junction,pathpool);
#endif
	  List_head_set(j,(void *) Junction_new_insertion(nindels,pathpool));
	  prev_ninserts = nindels;
	  j = List_next(j);

	  Intlist_head_set(s,nmismatches_i);
	  Intlist_head_set(s->rest,nmismatches_j);
	  s = Intlist_next(s);

	  Intlist_head_set(r,ref_nmismatches_i);
	  Intlist_head_set(r->rest,ref_nmismatches_j);
	  r = Intlist_next(r);
	}

      } else if (univdiagonal_i + max_deletionlen >= univdiagonal_j) {
	/* Deletion (or short intron */
	nindels = univdiagonal_j - univdiagonal_i;
	if ((indel_qpos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
							&ref_nmismatches_i,&ref_nmismatches_j,
							univdiagonal_i,/*indels*/-nindels,this->chrhigh,
							/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							/*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							/*pos5*/qstart,/*pos3*/qend,querylength,
							stage1->indelinfo,this->plusp,genestrand,
							/*want_lowest_coordinate_p*/true)) <= 0) {
	  debug14(printf("Deletion or short intron fails\n"));

	  /* Should advance pointers, if this code were used */

	} else {
	  /* No need to change univdiagonals */
	  p = Univcoordlist_next(p);

	  Intlist_head_set(q->rest,indel_qpos);
	  q = Intlist_next(q);

#ifdef ALLOCATE_UNSOLVED_JUNCTION
	  junction = (Junction_T) List_head(j);
	  Junction_free(&junction,pathpool);
#endif
	  /* deletionpos = (univdiagonal_i - querylength) + indel_qpos; */
	  List_head_set(j,(void *) Junction_new_deletion(nindels,pathpool));
	  prev_ninserts = 0;
	  j = List_next(j);

	  Intlist_head_set(s,nmismatches_i);
	  Intlist_head_set(s->rest,nmismatches_j);
	  s = Intlist_next(s);

	  Intlist_head_set(r,ref_nmismatches_i);
	  Intlist_head_set(r->rest,ref_nmismatches_j);
	  r = Intlist_next(r);
	}

	/* Splice_nomiddle call type 3 */
      } else if ((splice_qpos = Splice_nomiddle(&type,&trimpos5,&trimpos3,&supporti,&supportj,&nindels,&indel_qpos,
						&nmismatches_i,&nmismatches_j,&nmismatches_indel,
						&ref_nmismatches_i,&ref_nmismatches_j,
						&ref_nmismatches_indel,&donor_prob,&acceptor_prob,
						univdiagonal_i,univdiagonal_j,
						query_compress,this->plusp,this->chroffset,this->chrhigh,
						/*pos5*/qstart,/*pos3*/qend,querylength,
						stage1->indelinfo,stage1->spliceinfo,knownsplicing,intlistpool,
						sense_forward_p,genestrand,/*trim5p*/false,/*trim3p*/false,/*innerp*/true,
						/*allow_indels_p*/true)) < 0) {

#if 0
	/* Does not need to handle middle exons */
	if (middle_univdiagonal != 0) {
	  debug14(printf("Splice_resolve (qstart, sense): middle univdiagonal %u from %d to %d\n",
			 middle_univdiagonal,splice_qpos_i,splice_qpos_j));

	  splice_distance_j = univdiagonal_j - middle_univdiagonal;
	  splice_distance_i = middle_univdiagonal - univdiagonal_i;
	  
	  jnext = List_next(j);
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	  junction = (Junction_T) List_head(j);
	  Junction_free(&junction,pathpool);
#endif
	  List_head_set(j,(void *) Junction_new_splice(splice_distance_i,sensedir,
						       donor1_prob,acceptor1_prob,
						       pathpool)
			listpool_trace(__FILE__,__LINE__));
	  j->rest = Listpool_push(jnext,listpool,
				  (void *) Junction_new_splice(splice_distance_j,sensedir,
							       donor2_prob,acceptor2_prob,
							       pathpool)
				  listpool_trace(__FILE__,__LINE__));
	  prev_ninserts = 0;
	  j = jnext;
	  
	  pnext = Univcoordlist_next(p);
	  p->rest = Univcoordlistpool_push(pnext,univcoordlistpool,middle_univdiagonal
					   univcoordlistpool_trace(__FILE__,__LINE__));
	  p = pnext;
	  
	  /* This is the correct order for splice_qpos_j and splice_qpos_i */
	  qnext = Intlist_next(q);
	  Intlist_head_set(q->rest,splice_qpos_j);
	  q->rest = Intlistpool_push(qnext,intlistpool,splice_qpos_i
				     intlistpool_trace(__FILE__,__LINE__));
	  q = qnext;
	  
	  snext = Intlist_next(s);
	  Intlist_head_set(s,nmismatches_i);
	  s->rest = Intlistpool_push(snext,intlistpool,nmismatches_middle
				     intlistpool_trace(__FILE__,__LINE__));
	  Intlist_head_set(snext,nmismatches_j);
	  s = snext;

	  rnext = Intlist_next(r);
	  Intlist_head_set(r,ref_nmismatches_i);
	  r->rest = Intlistpool_push(rnext,intlistpool,ref_nmismatches_middle
				     intlistpool_trace(__FILE__,__LINE__));
	  Intlist_head_set(rnext,ref_nmismatches_j);
	  r = rnext;
	}
#endif

	all_solved_p = false;

	/* prev_ninserts = 0; -- Doesn't matter */
	j = List_next(j);
	p = Univcoordlist_next(p);
	q = Intlist_next(q);
	r = Intlist_next(r);
	s = Intlist_next(s);

      } else if (nindels == 0) {
	/* Splice only */
	splice_distance = univdiagonal_j - univdiagonal_i;
	debug14(printf("Splice_resolve (qstart, sense): splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
		       qstart,qend,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob));

	/* No need to change univdiagonals */
	p = Univcoordlist_next(p);

	Intlist_head_set(q->rest,splice_qpos);
	q = Intlist_next(q);

#ifdef ALLOCATE_UNSOLVED_JUNCTION
	junction = (Junction_T) List_head(j);
	Junction_free(&junction,pathpool);
#endif
	List_head_set(j,(void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
						     donor_prob,acceptor_prob,pathpool));
	prev_ninserts = 0;
	j = List_next(j);

	Intlist_head_set(s,nmismatches_i);
	Intlist_head_set(s->rest,nmismatches_j);
	s = Intlist_next(s);

	Intlist_head_set(r,ref_nmismatches_i);
	Intlist_head_set(r->rest,ref_nmismatches_j);
	r = Intlist_next(r);

      } else if (indel_qpos < splice_qpos) {
	debug14(printf("Splice_resolve (qstart, sense) with indel: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
		       qstart,qend,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob));

	debug14(printf("Before fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));
	
	/* Indel on segmenti: revise indel (based on segmenti) then splice */
	debug14(printf("Indel on segmenti: indel_qpos %d, splice_qpos %d\n",indel_qpos,splice_qpos));
	univdiagonal_k = univdiagonal_i + nindels; /* nindels = univdiagonal_k - univdiagonal_i */
	splice_distance = univdiagonal_j - univdiagonal_k;
	
	jnext = List_next(j);
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	junction = (Junction_T) List_head(j);
	Junction_free(&junction,pathpool);
#endif
	if (nindels < 0) {
	  List_head_set(j,(void *) Junction_new_insertion(-nindels,pathpool));
	} else {
	  /* deletionpos = segmenti_left + indel_qpos; */
	  List_head_set(j,(void *) Junction_new_deletion(nindels,pathpool));
	}

	j->rest = Listpool_push(jnext,listpool,
				(void *) Junction_new_splice(splice_distance,sensedir,
							     donor_prob,acceptor_prob,pathpool)
				listpool_trace(__FILE__,__LINE__));
	prev_ninserts = 0;
	j = jnext;
	  
	pnext = Univcoordlist_next(p);
	p->rest = Univcoordlistpool_push(pnext,univcoordlistpool,univdiagonal_k
					 univcoordlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	p = pnext;
	
	qmid = Intlist_next(q);
	qnext = Intlist_next(qmid);
	Intlist_head_set(qmid,indel_qpos);
	q = qmid->rest = Intlistpool_push(qnext,intlistpool,splice_qpos
					  intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	
	snext = Intlist_next(s);
	Intlist_head_set(s,nmismatches_i);
	s->rest = Intlistpool_push(snext,intlistpool,nmismatches_indel
				   intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	Intlist_head_set(snext,nmismatches_j);
	s = snext;

	rnext = Intlist_next(r);
	Intlist_head_set(r,ref_nmismatches_i);
	r->rest = Intlistpool_push(rnext,intlistpool,ref_nmismatches_indel
				   intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	Intlist_head_set(rnext,ref_nmismatches_j);
	r = rnext;
	
	debug14(printf("After fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));
	  
      } else {
	debug14(printf("Splice_resolve (qstart, sense) with indel: splice_qpos in range %d..%d is %d with distance %u, mismatches %d+%d, and probs %f and %f\n",
		       qstart,qend,splice_qpos,splice_distance,nmismatches_i,nmismatches_j,donor_prob,acceptor_prob));

	debug14(printf("Before fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));

	/* Indel on segmentj: revise splice then indel (based on segmentj) */
	debug14(printf("Indel on segmentj: splice_qpos %d, indel_qpos %d\n",splice_qpos,indel_qpos));
	univdiagonal_k = univdiagonal_j - nindels; /* nindels = univdiagonal_j - univdiagonal_k */
	splice_distance = univdiagonal_k - univdiagonal_i;
	
	jnext = List_next(j);
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	junction = (Junction_T) List_head(j);
	Junction_free(&junction,pathpool);
#endif
	List_head_set(j,(void *) Junction_new_splice(splice_distance,sensedir,donor_prob,acceptor_prob,pathpool));

	if (nindels < 0) {
	  j->rest = Listpool_push(jnext,listpool,(void *) Junction_new_insertion(-nindels,pathpool)
				  listpool_trace(__FILE__,__LINE__));
	  prev_ninserts = -nindels;
	} else {
	  /* deletionpos = segmentj_left + indel_qpos; */
	  j->rest = Listpool_push(jnext,listpool,(void *) Junction_new_deletion(nindels,pathpool)
				  listpool_trace(__FILE__,__LINE__));
	  prev_ninserts = 0;
	}
	j = jnext;
	  
	pnext = Univcoordlist_next(p);
	p->rest = Univcoordlistpool_push(pnext,univcoordlistpool,univdiagonal_k
					 univcoordlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	p = pnext;
	  
	qmid = Intlist_next(q);
	qnext = Intlist_next(qmid);
	Intlist_head_set(qmid,splice_qpos);
	q = qmid->rest = Intlistpool_push(qnext,intlistpool,indel_qpos
					  intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	
	snext = Intlist_next(s);
	Intlist_head_set(s,nmismatches_i);
	s->rest = Intlistpool_push(snext,intlistpool,nmismatches_indel
				   intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	Intlist_head_set(snext,nmismatches_j);
	s = snext;

	rnext = Intlist_next(r);
	Intlist_head_set(r,ref_nmismatches_i);
	r->rest = Intlistpool_push(rnext,intlistpool,ref_nmismatches_indel
				   intlistpool_trace(__FILE__,__LINE__)); /* Insert elt and advance */
	Intlist_head_set(rnext,ref_nmismatches_j);
	r = rnext;
	
	debug14(printf("After fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));
      }
    }
  }

#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
  printf("Endpoints: %s\n",Intlist_to_string(q));
  printf("Univdiagonals: %s\n",Univcoordlist_to_string(p));
  printf("Mismatches: %s\n",Intlist_to_string(s));
  printf("Junctions: ");
  Junction_print_list(j);
  printf("\n");
#endif

  assert(Univcoordlist_length(/*univdiagonals*/p) == Intlist_length(/*endpoints*/q) - 1);
  assert(Intlist_length(/*nmismatches*/s) == Intlist_length(/*endpoints*/q) - 1);
  assert(Intlist_length(/*ref_nmismatches*/r) == Intlist_length(/*endpoints*/q) - 1);
  assert(List_length(/*junctions*/j) == Intlist_length(/*endpoints*/q) - 2);


#if 0
  if (all_solved_p == false) {
    debug14(printf("Not all junctions were solved\n"));
  } else {
    /* Need to make copies of endpoints, junctions, univdiagonals, nmismatches, and ref_nmismatches */
    endpoints = Intlistpool_copy(this->endpoints,intlistpool);
    junctions = Listpool_copy(this->junctions,listpool);
  
    if (endpoints_acceptable_p(endpoints,junctions) == false) {
      debug14(printf("Endpoints were not acceptable\n"));
    } else {
      nmismatches = Intlistpool_copy(this->nmismatches,intlistpool);
      ref_nmismatches = Intlistpool_copy(this->ref_nmismatches,intlistpool);
      univdiagonals = Univcoordlistpool_copy(this->univdiagonals,univcoordlistpool);
  
#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
      printf("\nCombined ");
      Path_print(this);
#endif
    }
  }
#endif
	
  if (all_solved_p == false) {
    debug14(printf("Not all junctions were solved\n"));
    return (T) NULL;

  } else {
    /* assert(List_length(this->invalid_transcripts) == 1); */

    invalid_transcripts = this->invalid_transcripts;
    this->invalid_transcripts = (List_T) NULL;

    validp = false;
    for (t = invalid_transcripts; t != NULL; t = List_next(t)) {
      transcript = (Transcript_T) List_head(t);
      if (Transcript_remap_invalid(transcript,this,transcriptome,queryseq,
				   uintlistpool,listpool,transcriptpool) == true) {
	validp = true;
      }
    }
    if (validp == true) {
      this->transcriptome_method_p = true;
    } else {
      this->transcriptome_method_p = false;
    }
    Transcript_list_gc(&invalid_transcripts,listpool,transcriptpool);


    debug14(printf("After Path_solve_junctions:\n"));
    debug14(Path_print(this));
    debug14(printf("\n"));

    this->nmatches = -1;
    Path_eval_nmatches(&(*found_score),&(*end_score),this,query_compress_fwd,query_compress_rev);

    debug14(printf("Exiting Path_solve_junctions with path\n"));
    debug14(Path_print(this));

    return this;
  }
}


#else


/* Solves unsolved junctions */
/* Called by Trpath_convert procedures */
/* Set innerp to be true (looser criteria), because we are allowing for introns specified by known transcripts */

/* Uses Splice_resolve_qstart */
T
Path_solve_junctions (int *found_score, int *end_score, T this, int sensedir, int genestrand,
		      Compress_T query_compress,
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      Shortread_T queryseq, char *queryptr, int querylength,

		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      Localdb_T localdb, int localdb_nmismatches_allowed,

		      Stage1_T stage1, Knownsplicing_T knownsplicing,
		      Uintlistpool_T uintlistpool, Intlistpool_T intlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  
  int nmismatches_i, nmismatches_j, nmismatches_indel,
    ref_nmismatches_i, ref_nmismatches_j, ref_nmismatches_indel;
  int nindels;
  double donor_prob, acceptor_prob;

  int type;
  int trimpos5, trimpos3, supporti, supportj;
  int splice_qpos, indel_qpos, prev_ninserts, qstart, qend;
  Univcoord_T univdiagonal_i, univdiagonal_j, univdiagonal_k,
    segmenti_left, segmentj_left;
  Chrpos_T splice_distance;

  int ignore;
  Univcoordlist_T univdiagonals;
  Intlist_T nmismatches, ref_nmismatches, endpoints;
  List_T junctions;

  Junction_T junction;
  List_T j;
  Intlist_T q, r, s;

  Junction_T splice_junction, indel_junction;
  List_T exit_list, first_list, last_list, *jprev;	/* junctions */
  Intlist_T entry_intlist, exit_intlist, first_intlist, last_intlist,
    *rprev, *sprev; /* endpoints, ref_nmismatches, nmismatches */
  Univcoordlist_T exit_univcoordlist,
    first_univcoordlist, last_univcoordlist, p, *pprev; /* univdiagonals */

  bool all_solved_p = true;
  bool sense_forward_p = (this->sensedir == SENSE_FORWARD) ? true : false;

  List_T invalid_transcripts, t;
  Transcript_T transcript;
  bool validp;


  debug14(printf("Entered Path_solve_junctions with path\n"));
  debug14(Path_print(this));

  prev_ninserts = 0;
  jprev = &this->junctions;
  pprev = &this->univdiagonals;
  q = this->endpoints;
  rprev = &this->ref_nmismatches;
  sprev = &this->nmismatches;

  while (*jprev != NULL) {
#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
    printf("Endpoints: %s\n",Intlist_to_string(q));
    printf("Univdiagonals: %s\n",Univcoordlist_to_string(*pprev));
    printf("Mismatches: %s\n",Intlist_to_string(*sprev));
    printf("Junctions: ");
    Junction_print_list(*jprev);
    printf("\n");
#endif

    assert(Univcoordlist_length(/*univdiagonals*/*pprev) == Intlist_length(/*endpoints*/q) - 1);
    assert(Intlist_length(/*nmismatches*/*sprev) == Intlist_length(/*endpoints*/q) - 1);
    assert(Intlist_length(/*ref_nmismatches*/*rprev) == Intlist_length(/*endpoints*/q) - 1);
    assert(List_length(/*junctions*/*jprev) == Intlist_length(/*endpoints*/q) - 2);

    junction = (Junction_T) List_head(*jprev);
    debug14(Junction_print(junction));

    if (
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	Junction_type(junction) != UNSOLVED_JUNCTION
#else
	junction != JUNCTION_UNSOLVED
#endif
	) {

      prev_ninserts = Junction_ninserts(junction);
      jprev = &(*jprev)->rest;	/* j = List_next(j); */
      pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */
      q = Intlist_next(q);
      rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
      sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

    } else {
      debug14(printf("Solving unsolved junction\n"));

      univdiagonal_i = Univcoordlist_head(*pprev);
      univdiagonal_j = Univcoordlist_head(Univcoordlist_next(*pprev));
      /* segmenti_left = univdiagonal_i - querylength; */
      /* segmentj_left = univdiagonal_j - querylength; */

      qstart = Intlist_head(q) + prev_ninserts;	/* The endpoint before the junction */
      qend = Intlist_head(Intlist_next(Intlist_next(q))); /* The endpoint after the junction */

#if 0
      if (univdiagonal_i > univdiagonal_j + max_insertionlen) {
	/* Possible with segment search */
	all_solved_p = false;

	prev_ninserts = 0;
	jprev = &(*jprev)->rest; /* j = List_next(j); */
	pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */
	q = Intlist_next(q);
	rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	sprev = &(*sprev)->rest; /* s = Intlist_next(s); */
      }
#endif

      if (univdiagonal_i > univdiagonal_j) {
	nindels = univdiagonal_i - univdiagonal_j;
	if ((indel_qpos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
							&ref_nmismatches_i,&ref_nmismatches_j,
							univdiagonal_i,/*indels*/+nindels,this->chrhigh,
							/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							/*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							/*pos5*/qstart,/*pos3*/qend,querylength,
							 stage1->indelinfo,this->plusp,genestrand,
							/*want_lowest_coordinate_p*/true)) <= 0) {
	  debug14(printf("Insertion fails\n"));

	  /* Advance pointers */
	  all_solved_p = false;
	  prev_ninserts = 0;

	  jprev = &(*jprev)->rest; /* j = List_next(j); */
	  pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */
	  q = Intlist_next(q);
	  rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	  sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

	} else {
	  /* No need to change univdiagonals */
	  pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */

	  Intlist_head_set(q->rest,indel_qpos);
	  q = Intlist_next(q);

#ifdef ALLOCATE_UNSOLVED_JUNCTION
	  junction = (Junction_T) List_head(j);
	  Junction_free(&junction,pathpool);
#endif
	  List_head_set(*jprev,(void *) Junction_new_insertion(nindels,pathpool));
	  prev_ninserts = nindels;
	  jprev = &(*jprev)->rest; /* j = List_next(j); */

	  Intlist_head_set(*sprev,nmismatches_i);
	  Intlist_head_set((*sprev)->rest,nmismatches_j);
	  sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

	  Intlist_head_set(*rprev,ref_nmismatches_i);
	  Intlist_head_set((*rprev)->rest,ref_nmismatches_j);
	  rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	}

      } else if (univdiagonal_i + max_deletionlen >= univdiagonal_j) {
	/* Deletion (or short intron */
	nindels = univdiagonal_j - univdiagonal_i;
	if ((indel_qpos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
							&ref_nmismatches_i,&ref_nmismatches_j,
							univdiagonal_i,/*indels*/-nindels,this->chrhigh,
							/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							/*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
							/*pos5*/qstart,/*pos3*/qend,querylength,
							stage1->indelinfo,this->plusp,genestrand,
							/*want_lowest_coordinate_p*/true)) <= 0) {
	  debug14(printf("Deletion or short intron fails\n"));

	  /* Advance pointers */
	  all_solved_p = false;
	  prev_ninserts = 0;

	  jprev = &(*jprev)->rest; /* j = List_next(j); */
	  pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */
	  q = Intlist_next(q);
	  rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	  sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

	} else {
	  /* No need to change univdiagonals */
	  pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */

	  Intlist_head_set(q->rest,indel_qpos);
	  q = Intlist_next(q);

#ifdef ALLOCATE_UNSOLVED_JUNCTION
	  junction = (Junction_T) List_head(j);
	  Junction_free(&junction,pathpool);
#endif
	  /* deletionpos = (univdiagonal_i - querylength) + indel_qpos; */
	  List_head_set(*jprev,(void *) Junction_new_deletion(nindels,pathpool));
	  prev_ninserts = 0;
	  jprev = &(*jprev)->rest; /* j = List_next(j); */

	  Intlist_head_set(*sprev,nmismatches_i);
	  Intlist_head_set((*sprev)->rest,nmismatches_j);
	  sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

	  Intlist_head_set(*rprev,ref_nmismatches_i);
	  Intlist_head_set((*rprev)->rest,ref_nmismatches_j);
	  rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	}

      } else if (Splice_resolve_qstart(&univdiagonals,&nmismatches,&ref_nmismatches,&endpoints,&junctions,
				       univdiagonal_i,univdiagonal_j,
				       stage1,query_compress,queryptr,
				       this->plusp,Path_chroffset(this),Path_chrhigh(this),

				       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				       localdb_tally,localdb_overflowp,max_nblocks,localdb,
				       localdb_nmismatches_allowed,

				       /*pos5*/qstart,/*pos3*/qend,querylength,
				       /*cont_trimpos*/-1,/*cont_nmismatches*/-1,
				       stage1->indelinfo,stage1->spliceinfo,knownsplicing,
				       univcoordlistpool,intlistpool,listpool,pathpool,
				       this->sensedir,genestrand,/*trim5p*/false,/*trim3p*/false,
				       /*innerp*/true,/*allow_indels_p*/true,/*find_middle_p*/false,/*middle_rank*/0) == 0) {
	all_solved_p = false;

	/* prev_ninserts = 0; -- Doesn't matter */
	jprev = &(*jprev)->rest; /* j = List_next(j); */
	pprev = &(*pprev)->rest; /* p = Univcoordlist_next(p); */
	q = Intlist_next(q);
	rprev = &(*rprev)->rest; /* r = Intlist_next(r); */
	sprev = &(*sprev)->rest; /* s = Intlist_next(s); */

      } else {
	debug14(printf("Before fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));

#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
	printf("Splice_resolve_qstart returns\n");
	printf("Univdiagonals: %s\n",Univcoordlist_to_string(univdiagonals));
	printf("Endpoints: %s\n",Intlist_to_string(endpoints));
	printf("Mismatches: %s\n",Intlist_to_string(nmismatches));
	Junction_print_list(junctions);
	printf("\n");
#endif

	/* Junctions: Needs to be (nsegments - 1).  Splice_resolve
	   provides (n_new_segments).  Replace unsolved junction with
	   the new ones */
	/* entry_list is a pointer, and handled below */
	exit_list = List_next(*jprev); /* Skips unsolved junction */

	debug14(printf("exit list is "));
	debug14(Junction_print_list(exit_list));
	debug14(printf("\n"));
#ifdef ALLOCATE_UNSOLVED_JUNCTION
	junction = (Junction_T) List_head(*jprev);
	Junction_free(&junction,pathpool);
#endif

	first_list = junctions;
	last_list = List_last_elt(junctions,NULL);
	prev_ninserts = Junction_ninserts((Junction_T) List_head(last_list));

	*jprev = first_list; /* Beginning of insert */
	last_list->rest = exit_list; /* End of insert */
	jprev = &last_list->rest;


	/* Univdiagonals: Needs to be (nsegments).  Splice_resolve
	   provides (n_new_segments).  Insert new univdiagonals after
	   skipping the first */
	exit_univcoordlist = Univcoordlist_next(*pprev);
	debug14(printf("exit univcoordlist is %s\n",Univcoordlist_to_string(exit_univcoordlist)));
	
	first_univcoordlist = univdiagonals;
	last_univcoordlist = Univcoordlist_last_elt(univdiagonals);
	debug14(printf("(1) last univcoordlist is %s\n",Univcoordlist_to_string(last_univcoordlist)));

	*pprev = first_univcoordlist;			 /* Beginning of insert */
	last_univcoordlist->rest = exit_univcoordlist; /* End of insert */
	debug14(printf("(2) last univcoordlist is %s\n",Univcoordlist_to_string(last_univcoordlist)));
	debug14(printf("(3) univdiagonals is %s\n",Univcoordlist_to_string(this->univdiagonals)));
	pprev = &last_univcoordlist->rest;
	

	/* Endpoints: Needs to be (nsegments + 1).  Splice_resolve
	   provides (n_new_segments + 1).  Skip the existing endpoint,
	   and insert the new ones after skipping the first  */

	entry_intlist = q;
	exit_intlist = Intlist_next(Intlist_next(q)); /* Skip the existing endpoint */

	first_intlist = Intlist_next(endpoints);
	last_intlist = Intlist_last_elt(endpoints);
	q = last_intlist;

	entry_intlist->rest = first_intlist; /* Beginning of insert */
	last_intlist->rest = exit_intlist; /* End of insert */
	
	/* Nmismatches and ref_nmismatches: Needs to be (nsegments).
	   Splice_resolve provides (n_new_segments + 1).  Overwrite
	   first and last values with the end values and insert the
	   middle ones */

	/* Number of nmismatches equals (nsegments + 1).  Overwriting first and last value for nmismatches. */
	exit_intlist = Intlist_next(*sprev);
	debug14(printf("exit intlist is %s\n",Intlist_to_string(exit_intlist)));

	Intlist_head_set(*sprev,Intlist_head(nmismatches));
	Intlist_head_set(exit_intlist,Intlist_last_value(nmismatches));

	first_intlist = nmismatches;
	last_intlist = Intlist_penultimate_elt(nmismatches);
	debug14(printf("(1) last intlist is %s\n",Intlist_to_string(last_intlist)));

	*sprev = first_intlist; /* Beginning of insert */
	last_intlist->rest = exit_intlist; /* End of insert */
	debug14(printf("(2) last intlist is %s\n",Intlist_to_string(last_intlist)));
	debug14(printf("(3) intlist is %s\n",Intlist_to_string(this->nmismatches)));
	sprev = &last_intlist->rest;


	exit_intlist = Intlist_next(*rprev);
	Intlist_head_set(*rprev,Intlist_head(ref_nmismatches));
	Intlist_head_set(exit_intlist,Intlist_last_value(ref_nmismatches));

	first_intlist = ref_nmismatches;
	last_intlist = Intlist_penultimate_elt(ref_nmismatches);

	*rprev = first_intlist; /* Beginning of insert */
	last_intlist->rest = exit_intlist; /* End of insert */
	rprev = &last_intlist->rest;
	
	debug14(printf("After fill:\n"));
	debug14(Path_print(this));
	debug14(printf("\n"));
      }
    }
  }

#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
  printf("Endpoints: %s\n",Intlist_to_string(q));
  printf("Univdiagonals: %s\n",Univcoordlist_to_string(*pprev));
  printf("Mismatches: %s\n",Intlist_to_string(*sprev));
  printf("Junctions: ");
  Junction_print_list(*jprev);
  printf("\n");
#endif

  assert(Univcoordlist_length(/*univdiagonals*/*pprev) == Intlist_length(/*endpoints*/q) - 1);
  assert(Intlist_length(/*nmismatches*/*sprev) == Intlist_length(/*endpoints*/q) - 1);
  assert(Intlist_length(/*ref_nmismatches*/*rprev) == Intlist_length(/*endpoints*/q) - 1);
  assert(List_length(/*junctions*/*jprev) == Intlist_length(/*endpoints*/q) - 2);


#if 0
  if (all_solved_p == false) {
    debug14(printf("Not all junctions were solved\n"));
  } else {
    /* Need to make copies of endpoints, junctions, univdiagonals, nmismatches, and ref_nmismatches */
    endpoints = Intlistpool_copy(this->endpoints,intlistpool);
    junctions = Listpool_copy(this->junctions,listpool);
  
    if (endpoints_acceptable_p(endpoints,junctions) == false) {
      debug14(printf("Endpoints were not acceptable\n"));
    } else {
      nmismatches = Intlistpool_copy(this->nmismatches,intlistpool);
      ref_nmismatches = Intlistpool_copy(this->ref_nmismatches,intlistpool);
      univdiagonals = Univcoordlistpool_copy(this->univdiagonals,univcoordlistpool);
  
#if defined(DEBUG14) || defined(DEBUG_PATH_CORE)
      printf("\nCombined ");
      Path_print(this);
#endif
    }
  }
#endif
	
  if (all_solved_p == false) {
    debug14(printf("Not all junctions were solved\n"));
    /* Caller will free the path */
    return (T) NULL;

  } else {
    /* assert(List_length(this->invalid_transcripts) == 1); */

    invalid_transcripts = this->invalid_transcripts;
    this->invalid_transcripts = (List_T) NULL;

    validp = false;
    for (t = invalid_transcripts; t != NULL; t = List_next(t)) {
      transcript = (Transcript_T) List_head(t);
      if (Transcript_remap_invalid(transcript,this,transcriptome,queryseq,
				   uintlistpool,listpool,transcriptpool) == true) {
	validp = true;
      }
    }
    if (validp == true) {
      this->transcriptome_method_p = true;
    } else {
      this->transcriptome_method_p = false;
    }
    Transcript_list_gc(&invalid_transcripts,listpool,transcriptpool);


    debug14(printf("After Path_solve_junctions:\n"));
    debug14(Path_print(this));
    debug14(printf("\n"));

    this->nmatches = -1;
    Path_eval_nmatches(&(*found_score),&(*end_score),this,query_compress_fwd,query_compress_rev);

    debug14(printf("Exiting Path_solve_junctions with path\n"));
    debug14(Path_print(this));

    return this;
  }
}

#endif


void
Path_solve_setup (bool *circularp_in, Transcriptome_T transcriptome_in, EF64_T chromosome_ef64_in,
		  Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in, Univcoord_T genomelength_in,
		  int index1part_in, int min_intronlength_in,
		  int max_insertionlen_in, int max_deletionlen_in, int shortsplicedist_in,
		  bool novelsplicingp, bool knownsplicingp, int n_splice_paths_in,
		  int max_recursive_ends_in, int min_support_indel_in,
		  int indel_support_mismatch_penalty_in,
		  bool materialize_equiv_splice_alts_p_in) {

  circularp = circularp_in;

  transcriptome = transcriptome_in;
  chromosome_ef64 = chromosome_ef64_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  genomelength = genomelength_in;

  index1part = index1part_in;

  localdb = (Localdb_T) NULL;

  min_intronlength = min_intronlength_in;

  max_insertionlen = max_insertionlen_in;
  max_deletionlen = max_deletionlen_in;
  shortsplicedist = shortsplicedist_in;
  n_splice_paths = n_splice_paths_in > 0 ? n_splice_paths_in : 1;
  max_recursive_ends = max_recursive_ends_in;
  min_support_indel = min_support_indel_in;
  indel_support_mismatch_penalty = indel_support_mismatch_penalty_in;
  materialize_equiv_splice_alts_p = materialize_equiv_splice_alts_p_in;

  if (novelsplicingp == true || knownsplicingp == true) {
    splicingp = true;
  } else {
    splicingp = false;
  }

  return;
}


void
Path_solve_localdb_setup (Localdb_T localdb_in) {
  localdb = localdb_in;
  return;
}
