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
#include "align-cascade.h"
#include "debug_groups.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For rint */
#include <string.h>		/* For memset */
#if defined(GSNAP_DEBUG_SEGMENT_STREAM_STATS) || defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
#include <sys/time.h>
#endif

#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"

#include "list.h"
#include "compress.h"
#include "sedgesort.h"
#include "univdiagdef.h"

#include "genomebits_mismatches.h" /* For MISMATCH_EXTRA */
#include "orderstat.h"
#include "transcriptome-search.h"
#include "tr-extension-search.h"
#include "tr-chain-search.h"
#include "kmer-search.h"
#include "extension-search.h"
/* #include "merge-search.h" */
/* #include "segment-search.h" */

#include "trpath-solve.h"
#include "trpath-convert.h"

#include "auxinfo.h"
#include "path.h"
#include "path-solve.h"
#include "path-fusion.h"

#include "transcript-remap.h"
#include "transcript-velocity.h"
#include "path-eval.h"


/* Also defined in stage1hr-paired.c */
#define EXON_LENGTH 50		/* Only 5% of exons are less than 50 bp */

#define MIN_SIZELIMIT 100
#define MAX_HITS_EXACT 100	/* Excessive exact paths indicate a repetitive sequence */
#define MAX_END_POSITIONS 100	/* For END method.  Could be different for DNA-seq compared with RNA-seq */
#define SEGMENT_CHAIN_SEEN_SIZE 65536
#define SEGMENT_CHAIN_SEEN_MASK (SEGMENT_CHAIN_SEEN_SIZE - 1)


static Mode_T mode;
static int index1part;
static int index1interval;
static int index1part_tr;

static Transcriptome_T transcriptome;
static bool transcriptome_align_p;
static bool genome_align_p;
static Align_cascade_T align_cascade = ALIGN_CASCADE_DEFAULT;
static bool refine_partial_tr_p;
static bool tr_chain_search_p;

static double user_nmismatches_filter_float;
static double user_mincoverage_filter_float;
static bool split_paths_p;
static int split_path_max_overlap;
static int split_path_min_segment;
static double split_path_min_combined_coverage_float;

static bool splicingp;
static bool se_ext_dedup_exact_chains_p;
static bool se_ext_equiv_prev_prune_p;
static bool se_keep_trpath_candidates_p;
static bool se_segment_stream_chains_p;
static bool se_segment_dedup_exact_chains_p;
static bool se_segment_prefix_skip_p;
static bool se_segment_qstart_anchor_cache_p;
static bool se_path_qend_local_cache_p;
static bool se_emitted_path_dedup_p;
static bool se_segment_equiv_prev_prune_p;
static int se_segment_cache_threshold;
static int se_segment_max_window;
static double se_segment_min_kmer_avg_quality = -1.0;
static int se_segment_max_kmer_bad_bases = -1;
static int se_segment_bad_base_quality = 10;
static int segment_max_avg_kmer_positions = 0;
static double se_max_mismatches_filter_float;
static double se_max_softclip_filter_float;

struct Segment_chain_key_elt_T {
  int qstart;
  int qend;
  int nmismatches;
  Univcoord_T univdiagonal;
};

struct Segment_chain_seen_entry_T {
  bool occupiedp;
  unsigned long long hash;
  bool leftp;
  bool plusp;
  int genestrand;
  int querylength;
  int max_nblocks;
  int localdb_nmismatches_allowed;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  int chain_length;
  struct Segment_chain_key_elt_T *chain_elts;
};

struct Segment_chain_seen_T {
  int nentries;
  struct Segment_chain_seen_entry_T entries[SEGMENT_CHAIN_SEEN_SIZE];
  unsigned long long lookups;
  unsigned long long marks;
  unsigned long long skips;
  unsigned long long full_misses;
  unsigned long long collisions;
  unsigned long long null_chains;
};

#if defined(GSNAP_DEBUG_SEGMENT_STREAM_STATS) || defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
static unsigned long long
stage1_single_now_us (void) {
  struct timeval tv;

  gettimeofday(&tv,(struct timezone *) NULL);
  return (unsigned long long) tv.tv_sec * 1000000ULL + (unsigned long long) tv.tv_usec;
}
#endif

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
#define segment_stream_now_us stage1_single_now_us
#endif


#if defined(DEBUG) || defined(DEBUG_SINGLE_FLOW)
#define debug(x) x
#else
#define debug(x)
#endif

/* Deep-tier single-end narration: per-stage path counts, found_score transitions, and branch decisions. */
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
#define debug1(x) x
#else
#define debug1(x)
#endif

static bool
single_ext_univdiag_equal_p (Univdiag_T x, Univdiag_T y) {
  if (x == y) {
    return true;
  } else if (x == (Univdiag_T) NULL || y == (Univdiag_T) NULL) {
    return false;
  } else if (x->qstart == y->qstart && x->qend == y->qend &&
	     x->nmismatches == y->nmismatches &&
	     x->univdiagonal == y->univdiagonal) {
    return true;
  } else {
    return false;
  }
}

static bool
single_ext_chain_equal_p (List_T x, List_T y) {
  while (x != (List_T) NULL && y != (List_T) NULL) {
    if (single_ext_univdiag_equal_p((Univdiag_T) List_head(x),
				    (Univdiag_T) List_head(y)) == false) {
      return false;
    } else {
      x = List_next(x);
      y = List_next(y);
    }
  }

  return (x == (List_T) NULL && y == (List_T) NULL) ? true : false;
}

static bool
single_ext_chain_seen_p (List_T seen_chains, List_T chain) {
  List_T p;

  for (p = seen_chains; p != NULL; p = List_next(p)) {
    if (single_ext_chain_equal_p((List_T) List_head(p),chain) == true) {
      return true;
    }
  }

  return false;
}

#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
static unsigned long long
single_debug_mix_u64 (unsigned long long hash, unsigned long long value) {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
  return hash;
}

static Univdiag_T
single_debug_chain_last (List_T chain) {
  List_T p;
  Univdiag_T last = (Univdiag_T) NULL;

  for (p = chain; p != NULL; p = List_next(p)) {
    last = (Univdiag_T) List_head(p);
  }

  return last;
}

static bool
single_debug_univdiag_equal_p (Univdiag_T x, Univdiag_T y) {
  if (x == y) {
    return true;
  } else if (x == (Univdiag_T) NULL || y == (Univdiag_T) NULL) {
    return false;
  } else if (x->qstart == y->qstart && x->qend == y->qend &&
	     x->nmismatches == y->nmismatches &&
	     x->univdiagonal == y->univdiagonal) {
    return true;
  } else {
    return false;
  }
}

static bool
single_debug_chain_equal_p (List_T x, List_T y) {
  Univdiag_T univdiag_x, univdiag_y;

  while (x != (List_T) NULL && y != (List_T) NULL) {
    univdiag_x = (Univdiag_T) List_head(x);
    univdiag_y = (Univdiag_T) List_head(y);
    if (single_debug_univdiag_equal_p(univdiag_x,univdiag_y) == false) {
      return false;
    } else {
      x = List_next(x);
      y = List_next(y);
    }
  }

  return (x == (List_T) NULL && y == (List_T) NULL) ? true : false;
}

static unsigned long long
single_debug_chain_hash (List_T chain) {
  List_T p;
  Univdiag_T univdiag;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    if (univdiag == (Univdiag_T) NULL) {
      hash = single_debug_mix_u64(hash,0xfeedfaceULL);
    } else {
      hash = single_debug_mix_u64(hash,(unsigned long long) univdiag->qstart);
      hash = single_debug_mix_u64(hash,(unsigned long long) univdiag->qend);
      hash = single_debug_mix_u64(hash,(unsigned long long) univdiag->nmismatches);
      hash = single_debug_mix_u64(hash,(unsigned long long) univdiag->univdiagonal);
    }
    n++;
  }

  return single_debug_mix_u64(hash,(unsigned long long) n);
}

static int
single_debug_chain_duplicate_rank (List_T chain_list, List_T current_cell) {
  List_T p, chain;
  int rank = 1;

  if (current_cell == (List_T) NULL) {
    return 0;
  } else {
    chain = (List_T) List_head(current_cell);
    for (p = chain_list; p != NULL && p != current_cell; p = List_next(p)) {
      if (single_debug_chain_equal_p((List_T) List_head(p),chain) == true) {
	return rank;
      } else {
	rank++;
      }
    }
    return 0;
  }
}

static int
single_debug_chain_unique_count (List_T chain_list) {
  List_T p;
  int unique = 0;

  for (p = chain_list; p != NULL; p = List_next(p)) {
    if (single_debug_chain_duplicate_rank(chain_list,p) == 0) {
      unique++;
    }
  }

  return unique;
}

static void
single_debug_print_univdiag (FILE *fp, Univdiag_T univdiag) {
  if (univdiag == (Univdiag_T) NULL) {
    fprintf(fp,"null");
  } else {
    fprintf(fp,"%d..%d/%u",univdiag->qstart,univdiag->qend,univdiag->univdiagonal);
  }
  return;
}

static void
single_debug_print_chain_full (FILE *fp, List_T chain) {
  List_T p;
  Univdiag_T univdiag;
  int i = 0;

  fprintf(fp," chain_full=");
  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    if (i++ > 0) {
      fprintf(fp,";");
    }
    single_debug_print_univdiag(fp,univdiag);
    if (univdiag != (Univdiag_T) NULL) {
      fprintf(fp,"/nm=%d",univdiag->nmismatches);
    }
  }
  return;
}

static unsigned long long
single_debug_intlist_hash (Intlist_T list) {
  Intlist_T p;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = Intlist_next(p)) {
    hash = single_debug_mix_u64(hash,(unsigned long long) Intlist_head(p));
    n++;
  }

  return single_debug_mix_u64(hash,(unsigned long long) n);
}

static unsigned long long
single_debug_univcoordlist_hash (Univcoordlist_T list) {
  Univcoordlist_T p;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = Univcoordlist_next(p)) {
    hash = single_debug_mix_u64(hash,(unsigned long long) Univcoordlist_head(p));
    n++;
  }

  return single_debug_mix_u64(hash,(unsigned long long) n);
}

static bool
single_debug_junction_equal_p (Junction_T x, Junction_T y) {
  if (x == y) {
    return true;
  } else if (x == (Junction_T) NULL || y == (Junction_T) NULL) {
    return false;
  } else if (x->type == y->type && x->nindels == y->nindels &&
	     x->splice_distance == y->splice_distance &&
	     x->sensedir == y->sensedir &&
	     x->donor1 == y->donor1 && x->donor2 == y->donor2 &&
	     x->acceptor1 == y->acceptor1 && x->acceptor2 == y->acceptor2 &&
	     x->donor_prob == y->donor_prob &&
	     x->acceptor_prob == y->acceptor_prob) {
    return true;
  } else {
    return false;
  }
}

static bool
single_debug_junction_list_equal_p (List_T x, List_T y) {
  while (x != (List_T) NULL && y != (List_T) NULL) {
    if (single_debug_junction_equal_p((Junction_T) List_head(x),(Junction_T) List_head(y)) == false) {
      return false;
    } else {
      x = List_next(x);
      y = List_next(y);
    }
  }

  return (x == (List_T) NULL && y == (List_T) NULL) ? true : false;
}

static unsigned long long
single_debug_junction_list_hash (List_T list) {
  List_T p;
  Junction_T junction;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  for (p = list; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (junction == (Junction_T) NULL) {
      hash = single_debug_mix_u64(hash,0xabcdefULL);
    } else {
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->type);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->nindels);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->splice_distance);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->sensedir);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->donor1);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->donor2);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->acceptor1);
      hash = single_debug_mix_u64(hash,(unsigned long long) junction->acceptor2);
    }
    n++;
  }

  return single_debug_mix_u64(hash,(unsigned long long) n);
}

static unsigned long long
single_debug_path_hash (Path_T path) {
  unsigned long long hash = 1469598103934665603ULL;

  if (path == (Path_T) NULL) {
    return 0ULL;
  } else {
    hash = single_debug_mix_u64(hash,(unsigned long long) path->plusp);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->genestrand);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->sensedir);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->querylength);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->chrnum);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->chroffset);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->chrhigh);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->main_univdiagonal);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->found_score);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->end_score);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->score_within_trims);
    hash = single_debug_mix_u64(hash,(unsigned long long) path->method);
    hash = single_debug_mix_u64(hash,single_debug_intlist_hash(path->endpoints));
    hash = single_debug_mix_u64(hash,single_debug_univcoordlist_hash(path->univdiagonals));
    hash = single_debug_mix_u64(hash,single_debug_intlist_hash(path->nmismatches));
    hash = single_debug_mix_u64(hash,single_debug_intlist_hash(path->ref_nmismatches));
    hash = single_debug_mix_u64(hash,single_debug_junction_list_hash(path->junctions));
    return hash;
  }
}

static bool
single_debug_path_equal_p (Path_T x, Path_T y) {
  if (x == y) {
    return true;
  } else if (x == (Path_T) NULL || y == (Path_T) NULL) {
    return false;
  } else if (x->plusp == y->plusp && x->genestrand == y->genestrand &&
	     x->sensedir == y->sensedir && x->querylength == y->querylength &&
	     x->chrnum == y->chrnum && x->chroffset == y->chroffset &&
	     x->chrhigh == y->chrhigh &&
	     x->main_univdiagonal == y->main_univdiagonal &&
	     x->found_score == y->found_score && x->end_score == y->end_score &&
	     x->score_within_trims == y->score_within_trims &&
	     x->method == y->method &&
	     Intlist_equal(x->endpoints,y->endpoints) == true &&
	     Univcoordlist_equal(x->univdiagonals,y->univdiagonals) == true &&
	     Intlist_equal(x->nmismatches,y->nmismatches) == true &&
	     Intlist_equal(x->ref_nmismatches,y->ref_nmismatches) == true &&
	     single_debug_junction_list_equal_p(x->junctions,y->junctions) == true) {
    return true;
  } else {
    return false;
  }
}

static int
single_debug_path_duplicate_rank_after_new (List_T path_list, int nnew, Path_T path) {
  List_T p;
  int skip, rank;

  for (skip = 0, p = path_list; p != NULL && skip < nnew; p = List_next(p), skip++) {
    /* Skip paths added by this solve call. */
  }
  for (rank = nnew + 1; p != NULL; p = List_next(p), rank++) {
    if (single_debug_path_equal_p((Path_T) List_head(p),path) == true) {
      return rank;
    }
  }

  return 0;
}

static void
single_debug_print_junction_list (FILE *fp, List_T list) {
  List_T p;
  Junction_T junction;
  int i = 0;

  fprintf(fp,"[");
  for (p = list; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (i++ > 0) {
      fprintf(fp,";");
    }
    if (junction == (Junction_T) NULL) {
      fprintf(fp,"UNSOLVED");
    } else {
      fprintf(fp,"%d/indel=%d/dist=%u/sense=%d",
	      junction->type,junction->nindels,junction->splice_distance,junction->sensedir);
    }
  }
  fprintf(fp,"]");
  return;
}

static void
single_debug_print_path_added (const char *dir, int nsplices, int chaini,
			       const char *bucket, int rank, int older_dup_rank,
			       Path_T path) {
  char *endpoints, *univdiagonals, *nmismatches, *ref_nmismatches;

  if (path == (Path_T) NULL) {
    fprintf(stderr,
	    "[EXT_PATH_ADDED] dir=%s nsplices=%d chain=%d bucket=%s rank=%d path=null\n",
	    dir,nsplices,chaini,bucket,rank);
  } else {
    endpoints = Intlist_to_string(path->endpoints);
    univdiagonals = Univcoordlist_to_string(path->univdiagonals);
    nmismatches = Intlist_to_string(path->nmismatches);
    ref_nmismatches = Intlist_to_string(path->ref_nmismatches);

    fprintf(stderr,
	    "[EXT_PATH_ADDED] dir=%s nsplices=%d chain=%d bucket=%s rank=%d "
	    "path=%p exact_dup_in_older_rank=%d hash=%016llx plusp=%d sensedir=%d "
	    "method=%d chr=%d:%u..%u main=%u score found=%d end=%d within=%d "
	    "endpoints=%s univdiags=%s nm=%s refnm=%s junctions=",
	    dir,nsplices,chaini,bucket,rank,(void *) path,older_dup_rank,
	    single_debug_path_hash(path),path->plusp,path->sensedir,path->method,
	    path->chrnum,path->chroffset,path->chrhigh,path->main_univdiagonal,
	    path->found_score,path->end_score,path->score_within_trims,
	    endpoints,univdiagonals,nmismatches,ref_nmismatches);
    single_debug_print_junction_list(stderr,path->junctions);
    fprintf(stderr,"\n");

    FREE(endpoints);
    FREE(univdiagonals);
    FREE(nmismatches);
    FREE(ref_nmismatches);
  }

  return;
}

static void
single_debug_print_new_paths (const char *dir, int nsplices, int chaini,
			      const char *bucket, int before_count, List_T paths) {
  List_T p;
  Path_T path;
  int after_count = List_length(paths);
  int nnew = after_count - before_count;
  int rank;

  if (nnew <= 0) {
    return;
  } else {
    fprintf(stderr,
	    "[EXT_PATH_DELTA] dir=%s nsplices=%d chain=%d bucket=%s before=%d after=%d nnew=%d\n",
	    dir,nsplices,chaini,bucket,before_count,after_count,nnew);
    for (rank = 1, p = paths; p != NULL && rank <= nnew; p = List_next(p), rank++) {
      path = (Path_T) List_head(p);
      single_debug_print_path_added(dir,nsplices,chaini,bucket,rank,
				    single_debug_path_duplicate_rank_after_new(paths,nnew,path),
				    path);
    }
  }

  return;
}

static void
single_debug_print_chain_catalog (const char *dir, int nsplices, int candidatei,
				  Univcoord_T anchor_univdiagonal, Auxinfo_T auxinfo,
				  List_T chain_list) {
  List_T p, chain;
  int chaini, exact_dup_of;
  int chain_count = List_length(chain_list);
  int unique_count = single_debug_chain_unique_count(chain_list);

  fprintf(stderr,
	  "[EXT_CHAIN_CATALOG_BEGIN] dir=%s nsplices=%d candidate=%d anchor=%u "
	  "aux=%d..%d chain_count=%d exact_unique=%d exact_duplicates=%d\n",
	  dir,nsplices,candidatei,anchor_univdiagonal,auxinfo->qstart,auxinfo->qend,
	  chain_count,unique_count,chain_count - unique_count);
  for (chaini = 1, p = chain_list; p != NULL; p = List_next(p), chaini++) {
    chain = (List_T) List_head(p);
    exact_dup_of = single_debug_chain_duplicate_rank(chain_list,p);
    fprintf(stderr,
	    "[EXT_CHAIN_CATALOG] dir=%s nsplices=%d candidate=%d chain=%d "
	    "chain_cell=%p chain_head_cell=%p chain_hash=%016llx exact_dup_of=%d len=%d",
	    dir,nsplices,candidatei,chaini,(void *) p,(void *) chain,
	    single_debug_chain_hash(chain),exact_dup_of,List_length(chain));
    single_debug_print_chain_full(stderr,chain);
    fprintf(stderr,"\n");
  }
  fprintf(stderr,
	  "[EXT_CHAIN_CATALOG_END] dir=%s nsplices=%d candidate=%d anchor=%u\n",
	  dir,nsplices,candidatei,anchor_univdiagonal);
  fflush(stderr);

  return;
}

static void
single_debug_print_chain_summary (const char *tag, const char *dir, int nsplices,
				  int chaini, int candidatei, int ncandidates,
				  Univcoord_T anchor_univdiagonal, Auxinfo_T auxinfo,
				  Pathstore_T pathstore, List_T chain_list, List_T chain_cell,
				  int found_score, int end_score,
				  List_T unextended_sense_paths,
				  List_T unextended_antisense_paths,
				  List_T sense_paths, List_T antisense_paths) {
  Univdiag_T head, next, last;
  List_T chain = (List_T) List_head(chain_cell);
  int exact_dup_of = single_debug_chain_duplicate_rank(chain_list,chain_cell);
  int chain_length = List_length(chain);

  head = chain ? (Univdiag_T) List_head(chain) : (Univdiag_T) NULL;
  next = chain && List_next(chain) ? (Univdiag_T) List_head(List_next(chain)) : (Univdiag_T) NULL;
  last = single_debug_chain_last(chain);

  fprintf(stderr,
	  "[EXT_CHAIN_%s] dir=%s nsplices=%d chain=%d candidate=%d/%d "
	  "anchor=%u aux=%d..%d elti=%d chr=%d:%u..%u chain_len=%d "
	  "chain_cell=%p chain_head_cell=%p chain_hash=%016llx exact_dup_of=%d "
	  "found=%d end=%d paths us=%d ua=%d s=%d a=%d head=",
	  tag,dir,nsplices,chaini,candidatei,ncandidates,
	  anchor_univdiagonal,auxinfo->qstart,auxinfo->qend,auxinfo->elti,
	  pathstore->chrnum,pathstore->chroffset,pathstore->chrhigh,chain_length,
	  (void *) chain_cell,(void *) chain,single_debug_chain_hash(chain),exact_dup_of,
	  found_score,end_score,
	  List_length(unextended_sense_paths),List_length(unextended_antisense_paths),
	  List_length(sense_paths),List_length(antisense_paths));
  single_debug_print_univdiag(stderr,head);
  fprintf(stderr," next=");
  single_debug_print_univdiag(stderr,next);
  fprintf(stderr," last=");
  single_debug_print_univdiag(stderr,last);
  single_debug_print_chain_full(stderr,chain);
  fprintf(stderr,"\n");
  fflush(stderr);

  return;
}

static void
single_debug_print_chain_result (const char *dir, int nsplices, int chaini,
				 unsigned long long elapsed_us,
				 int found_before, int end_before,
				 int found_after, int end_after,
				 int us_before, int ua_before, int s_before, int a_before,
				 List_T unextended_sense_paths,
				 List_T unextended_antisense_paths,
				 List_T sense_paths, List_T antisense_paths) {
  int us_after = List_length(unextended_sense_paths);
  int ua_after = List_length(unextended_antisense_paths);
  int s_after = List_length(sense_paths);
  int a_after = List_length(antisense_paths);

  fprintf(stderr,
	  "[EXT_CHAIN_END] dir=%s nsplices=%d chain=%d elapsed_us=%llu "
	  "found=%d->%d end=%d->%d paths us=%d->%d ua=%d->%d s=%d->%d a=%d->%d\n",
	  dir,nsplices,chaini,elapsed_us,found_before,found_after,end_before,end_after,
	  us_before,us_after,ua_before,ua_after,s_before,s_after,a_before,a_after);
  single_debug_print_new_paths(dir,nsplices,chaini,"unextended_sense",us_before,unextended_sense_paths);
  single_debug_print_new_paths(dir,nsplices,chaini,"unextended_antisense",ua_before,unextended_antisense_paths);
  single_debug_print_new_paths(dir,nsplices,chaini,"sense",s_before,sense_paths);
  single_debug_print_new_paths(dir,nsplices,chaini,"antisense",a_before,antisense_paths);
  fflush(stderr);

  return;
}
#endif

#define T Stage1_T

#if 0
static int
determine_sizelimit (T this, int querylength) {
  int cutoff, *set, count;
  int n;
  int query_lastpos, querypos;

  assert(querylength >= index1part);

  query_lastpos = querylength - index1part;
  set = (int *) MALLOC(2*(query_lastpos+1)*sizeof(int));
  n = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->validp[querypos] == true) {
      set[n++] = count = this->plus_npositions[querypos];
      set[n++] = count = this->minus_npositions[querypos];
    }
  }

  if (n < 5) {
    cutoff = MIN_SIZELIMIT;
  } else if ((cutoff = Orderstat_int_pct_inplace(set,n,/*pct*/0.60)) < MIN_SIZELIMIT) {
    cutoff = MIN_SIZELIMIT;
  }
  FREE(set);

  return cutoff;
}
#endif


#if defined(DEBUG) || defined(DEBUG_SINGLE_FLOW)
static void
list_paths (List_T list, char *destination, bool expected_sensedir) {
  List_T p;
  Path_T path;

  for (p = list; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    printf("Destination %s: ",destination);
    Path_print(path);
    assert(path->sensedir == expected_sensedir);
  }

  return;
}
#endif


Method_T
single_read_trdiagonal_method_goal (void) {
  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_PREV)) {
    return TR_PREVALENT;
  } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXACT)) {
    return TR_EXACT1;
  } else {
    return METHOD_INIT;
  }
}


bool
single_read_trdiagonal_method_available_p (Method_T last_method, Method_T method_goal) {
  return
    (last_method < TR_EXACT1 && method_goal >= TR_EXACT1 &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXACT)) ||
    (last_method < TR_PREVALENT && method_goal >= TR_PREVALENT &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_PREV));
}


Method_T
single_read_tr_method_goal (void) {
  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_CHAIN)) {
    return TR_CHAIN;
  } else if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXT)) {
    return TR_EXT;
  } else {
    return METHOD_INIT;
  }
}


bool
single_read_tr_method_available_p (Method_T last_method, Method_T method_goal) {
  return
    (last_method < TR_EXT && method_goal >= TR_EXT &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXT)) ||
    (last_method < TR_CHAIN && method_goal >= TR_CHAIN &&
     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_CHAIN));
}


/* Benchmarks show slight speed improvement with TR_EXACT2 before TR_ANYPAIR */
Method_T
single_read_next_method_trdiagonal (Method_T last_method, T this, int querylength,
				    Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    bool first_read_p) {

  if (last_method < TR_EXACT1 &&
      Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXACT)) {
    /* 1.1.  TR_EXACT1 */
    debug(printf("%s Read: 1.1.  Running Transcriptome_exact1 for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_init_end_tr(this,querylength);
    Transcriptome_exact1(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			 &this->_sense_trdiagonals,&this->n_sense_trdiagonals,
			 &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			 &this->_antisense_trdiagonals,&this->n_antisense_trdiagonals,
			 this,querylength);
    return TR_EXACT1;

#if 0
  } else if (last_method < TR_EXACT2) {
    /* Can lead to ends that don't match the transcript */
    Stage1_trdiagonals_gc(this);
    
    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_init_end2_positions_tr(this,querylength);
    Transcriptome_exact2(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			 &this->_sense_trdiagonals,&this->n_sense_trdiagonals,
			 &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			 &this->_antisense_trdiagonals,&this->n_antisense_trdiagonals,
			 this,querylength);
    return TR_EXACT2;
#endif

#if 0
  } else if (last_method < TR_ANYPAIR) {
    /* Can lead to incorrect read loci */
    debug(printf("%s Read: 1.2.  Running Transcriptome_anypair for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_trdiagonals_gc(this);

    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_fill_all_positions_tr(this,querylength);

    Transcriptome_anypair(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			  &this->_sense_trdiagonals,&this->sense_tstarts,&this->sense_tends,
			  &this->n_sense_trdiagonals,
			  &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			  &this->_antisense_trdiagonals,&this->antisense_tstarts,&this->antisense_tends,
			  &this->n_antisense_trdiagonals,this,querylength);

    return TR_ANYPAIR;
#endif

  } else if (last_method < TR_PREVALENT &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_PREV)) {
    debug(printf("%s Read: 1.3.  Running Transcriptome_prevalent for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_trdiagonals_gc(this);

    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_fill_all_positions_tr(this,querylength);

    Transcriptome_prevalent(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			    &this->_sense_trdiagonals,&this->sense_tstarts,&this->sense_tends,
			    &this->n_sense_trdiagonals,
			    &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			    &this->_antisense_trdiagonals,&this->antisense_tstarts,&this->antisense_tends,
			    &this->n_antisense_trdiagonals,
			    this,querylength,query_compress_fwd,query_compress_rev);
    
    return TR_PREVALENT;

  } else {
    fprintf(stderr,"No enabled trdiagonal method after %s\n",Method_string(last_method));
    abort();
  }
}


Method_T
single_read_next_method_tr (int *found_score, int *end_score, int sufficient_score, Method_T last_method,

			    List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
			    List_T *sense_trpaths, List_T *antisense_trpaths,

			    T this, int genestrand, int querylength,
			    int *mismatch_positions_alloc,
			    char *queryuc_ptr, char *queryrc,
			    Compress_T query_compress_fwd, Compress_T query_compress_rev,

			    int read_nmismatches_allowed,
			 
			    Trdiagpool_T trdiagpool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			    Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			    Hitlistpool_T hitlistpool, bool paired_end_p, bool first_read_p,
			    bool appendp) {

  if (appendp == false) {
    *partial_sense_trpaths = *partial_antisense_trpaths = (List_T) NULL;
    *sense_trpaths = *antisense_trpaths = (List_T) NULL;
  }

  if (last_method < TR_CHAIN &&
      Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_CHAIN)) {
    debug(printf("%s Read: 1.4.  Running Trchain_search for trpaths\n",
		 first_read_p ? "5'" : "3'"));
    Trchain_search(&(*found_score),&(*end_score),sufficient_score,
		   &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		   &(*sense_trpaths),&(*antisense_trpaths),

		   this,querylength,mismatch_positions_alloc,
		   queryuc_ptr,queryrc,
		   query_compress_fwd,query_compress_rev,paired_end_p,first_read_p,

		   trdiagpool,intlistpool,uintlistpool,listpool,
		   trpathpool,pathpool,hitlistpool,
		   genestrand,/*method*/TR_CHAIN);
    debug(printf(" => %d sense and %d antisense trpaths\n",
		 List_length(*sense_trpaths),List_length(*antisense_trpaths)));
    return TR_CHAIN;

  } else if (last_method < TR_EXT &&
	     Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_TR_EXT)) {

#if 0
    /* Supplanted by new Tr_extension_search,which also uses
       Genomebits_trim procedures, but still offers some speed
       advantage */
    /* 1.2.  TR_END */
    debug(printf("%s Read: 1.2.  Running Transcriptome_search_end\n",
		 first_read_p ? "5'" : "3'"));
    Transcriptome_search_end(&(*found_score),&(*end_score),&(*sense_trpaths),&(*antisense_trpaths),

			     this,querylength,query_compress_fwd,query_compress_rev, 
			     intlistpool,uintlistpool,trpathpool,hitlistpool,
			     /*method*/TR_END);
    return TR_END;

  } else if (last_method == TR_END) {
#endif

    /* 1.3.  TR_EXT */
    debug(printf("%s Read: 1.4.  Running Tr_extension_search for trpaths\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_fill_all_oligos_tr(this,querylength);
    Tr_extension_search(&(*found_score),&(*end_score),sufficient_score,
			&(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			&(*sense_trpaths),&(*antisense_trpaths),

			this,querylength,mismatch_positions_alloc,
			query_compress_fwd,query_compress_rev, 
			paired_end_p,first_read_p,

			trdiagpool,intlistpool,uintlistpool,listpool,
			trpathpool,pathpool,hitlistpool,
			read_nmismatches_allowed,genestrand,/*method*/TR_EXT);
    debug(printf(" => %d sense and %d antisense trpaths\n",
		 List_length(*sense_trpaths),List_length(*antisense_trpaths)));
    return TR_EXT;
    
  } else {
    fprintf(stderr,"No enabled trpath method after %s\n",Method_string(last_method));
    abort();
  }
}


static void
single_read_extend (int *found_score, int *end_score, T this,

		    List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		    List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		    List_T *unextended_sense_paths_gplus, List_T *unextended_sense_paths_gminus,
		    List_T *unextended_antisense_paths_gplus, List_T *unextended_antisense_paths_gminus,

		    Shortread_T queryseq, char *queryuc_ptr, char *queryrc, int querylength,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels,

		    int *mismatch_positions_alloc,
		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks, int localdb_nmismatches_allowed,

		    Intlistpool_T intlistpool,
		    Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Spliceendsgen_T spliceendsgen) {

  List_T filtered_paths, unextended_paths, complete_paths, p;
  Path_T path;

  /* Attempt to extend paths */
  debug(printf("Entered single_read_extend\n"));
  debug1(printf("Entered single_read_extend: found_score=%d end_score=%d  unextended sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d  complete sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d\n",
		*found_score,*end_score,
		List_length(*unextended_sense_paths_gplus), List_length(*unextended_sense_paths_gminus),
		List_length(*unextended_antisense_paths_gplus), List_length(*unextended_antisense_paths_gminus),
		List_length(*sense_paths_gplus), List_length(*sense_paths_gminus),
		List_length(*antisense_paths_gplus), List_length(*antisense_paths_gminus)));

#if 0
  /* Path_filter appears to select incorrect splice lengths */
  debug(printf("Before Path_filter: %d unextended_sense_paths_gplus\n",
	       List_length(this->unextended_sense_paths_gplus)));
  debug(printf("Before Path_filter: %d unextended_antisense_paths_gplus\n",
	       List_length(this->unextended_antisense_paths_gplus)));
  debug(printf("Before Path_filter: %d unextended_sense_paths_gminus\n",
	       List_length(this->unextended_sense_paths_gminus)));
  debug(printf("Before Path_filter: %d unextended_antisense_paths_gminus\n",
	       List_length(this->unextended_antisense_paths_gminus)));

  this->unextended_sense_paths_gplus = Path_filter(this->unextended_sense_paths_gplus,
						   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_sense_paths_gminus = Path_filter(this->unextended_sense_paths_gminus,
						    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_antisense_paths_gplus = Path_filter(this->unextended_antisense_paths_gplus,
						       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_antisense_paths_gminus = Path_filter(this->unextended_antisense_paths_gminus,
							intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#elif 0
  /* Path_unique gets rid of sense/antisense, needed for finding fusions */
  *unextended_sense_paths_gplus = Path_unique(*unextended_sense_paths_gplus,
					      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_sense_paths_gminus = Path_unique(*unextended_sense_paths_gminus,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_antisense_paths_gplus = Path_unique(*unextended_antisense_paths_gplus,
						  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_antisense_paths_gminus = Path_unique(*unextended_antisense_paths_gminus,
						   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#endif
      
  debug(printf("Have %d unextended_sense_paths_gplus\n",List_length(*unextended_sense_paths_gplus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&(*end_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryuc_ptr,querylength,this->mismatch_positions_alloc,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			    this,knownsplicing,knownindels,/*query_compress*/query_compress_fwd,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    /*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_sense_paths_gplus\n"));
      path->completep = true;
      *sense_paths_gplus = List_append(complete_paths,*sense_paths_gplus);

    } else {
      /* TODO: Path_gc(&unextended_paths); */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_sense_paths_gplus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_sense_paths_gplus = filtered_paths;
  debug1(printf("single_read_extend after sense_g+ block: complete=%d kept_unextended=%d found_score=%d end_score=%d\n",
		List_length(*sense_paths_gplus),List_length(*unextended_sense_paths_gplus),*found_score,*end_score));


  debug(printf("Have %d unextended_sense_paths_gminus\n",List_length(*unextended_sense_paths_gminus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths =
		Path_extend(&(*found_score),&(*end_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryrc,querylength,this->mismatch_positions_alloc,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			    this,knownsplicing,knownindels,/*query_compress*/query_compress_rev,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    /*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_sense_paths_gminus\n"));
      path->completep = true;
      *sense_paths_gminus = List_append(complete_paths,*sense_paths_gminus);

    } else {
      /* TODO: Path_gc(&unextended_paths); */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_sense_paths_gminus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_sense_paths_gminus = filtered_paths;
  debug1(printf("single_read_extend after sense_g- block: complete=%d kept_unextended=%d found_score=%d end_score=%d\n",
		List_length(*sense_paths_gminus),List_length(*unextended_sense_paths_gminus),*found_score,*end_score));


  debug(printf("Have %d unextended_antisense_paths_gplus\n",List_length(*unextended_antisense_paths_gplus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&(*end_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryuc_ptr,querylength,this->mismatch_positions_alloc,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			    this,knownsplicing,knownindels,/*query_compress*/query_compress_fwd,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    /*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_antisense_paths_gplus\n"));
      path->completep = true;
      *antisense_paths_gplus = List_append(complete_paths,*antisense_paths_gplus);
      
    } else {
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_antisense_paths_gplus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_antisense_paths_gplus = filtered_paths;
  debug1(printf("single_read_extend after anti_g+ block: complete=%d kept_unextended=%d found_score=%d end_score=%d\n",
		List_length(*antisense_paths_gplus),List_length(*unextended_antisense_paths_gplus),*found_score,*end_score));


  debug(printf("Have %d unextended_antisense_paths_gminus\n",List_length(*unextended_antisense_paths_gminus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&(*end_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryrc,querylength,this->mismatch_positions_alloc,
			    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

			    this,knownsplicing,knownindels,/*query_compress*/query_compress_rev,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    /*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_antisense_paths_gminus\n"));
      path->completep = true;
      *antisense_paths_gminus = List_append(complete_paths,*antisense_paths_gminus);
    } else {
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_antisense_paths_gminus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_antisense_paths_gminus = filtered_paths;
  debug1(printf("single_read_extend after anti_g- block: complete=%d kept_unextended=%d found_score=%d end_score=%d\n",
		List_length(*antisense_paths_gminus),List_length(*unextended_antisense_paths_gminus),*found_score,*end_score));
  debug1(printf("Leaving single_read_extend: found_score=%d end_score=%d  complete sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d  unextended sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d\n",
		*found_score,*end_score,
		List_length(*sense_paths_gplus), List_length(*sense_paths_gminus),
		List_length(*antisense_paths_gplus), List_length(*antisense_paths_gminus),
		List_length(*unextended_sense_paths_gplus), List_length(*unextended_sense_paths_gminus),
		List_length(*unextended_antisense_paths_gplus), List_length(*unextended_antisense_paths_gminus)));

  return;
}



#ifdef FUSION
static List_T
single_read_fusion (int *found_score, int *end_score, T this, int querylength,

		    List_T unextended_sense_paths_gplus, List_T unextended_sense_paths_gminus,
		    List_T unextended_antisense_paths_gplus, List_T unextended_antisense_paths_gminus,

		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Shortread_T queryseq, Knownsplicing_T knownsplicing, int nmismatches_allowed,

		    Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
		    Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T paths = NULL, p;
  Path_T path;

  debug1(printf("Entered single_read_fusion: found_score=%d end_score=%d  unextended sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d\n",
		*found_score,*end_score,
		List_length(unextended_sense_paths_gplus), List_length(unextended_sense_paths_gminus),
		List_length(unextended_antisense_paths_gplus), List_length(unextended_antisense_paths_gminus)));

  /* Look for possible fusions (or combinations of existing paths) */
  /* Find fusions.  Use code similar to finding outer fusions in Pathpair_eval_and_sort */

  debug(printf("Have %d unextended_sense_paths_gplus\n",List_length(unextended_sense_paths_gplus)));
  for (p = unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_plus(&(*found_score),&(*end_score),/*main*/path,this,

						      query_compress_fwd,query_compress_rev,
						      queryseq,querylength,knownsplicing,/*genestrand*/0,
						      nmismatches_allowed,
						      intlistpool,uintlistpool,univcoordlistpool,
						      listpool,univdiagpool,pathpool,vectorpool,
						      transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_plus(&(*found_score),&(*end_score),/*main*/path,this,

						    query_compress_fwd,query_compress_rev,
						    queryseq,querylength,knownsplicing,/*genestrand*/0,
						    nmismatches_allowed,
						    intlistpool,uintlistpool,univcoordlistpool,
						    listpool,univdiagpool,pathpool,vectorpool,
						    transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_antisense_paths_gplus\n",List_length(unextended_antisense_paths_gplus)));
  for (p = unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_plus(&(*found_score),&(*end_score),/*main*/path,this,

						      query_compress_fwd,query_compress_rev,
						      queryseq,querylength,knownsplicing,/*genestrand*/0,
						      nmismatches_allowed,
						      intlistpool,uintlistpool,univcoordlistpool,
						      listpool,univdiagpool,pathpool,vectorpool,
						      transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_plus(&(*found_score),&(*end_score),/*main*/path,this,

						    query_compress_fwd,query_compress_rev,
						    queryseq,querylength,knownsplicing,/*genestrand*/0,
						    nmismatches_allowed,
						    intlistpool,uintlistpool,univcoordlistpool,
						    listpool,univdiagpool,pathpool,vectorpool,
						    transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_sense_paths_gminus\n",List_length(unextended_sense_paths_gminus)));
  for (p = unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_minus(&(*found_score),&(*end_score),/*main*/path,this,

						       query_compress_fwd,query_compress_rev,
						       queryseq,querylength,knownsplicing,/*genestrand*/0,
						       nmismatches_allowed,
						       intlistpool,uintlistpool,univcoordlistpool,
						       listpool,univdiagpool,pathpool,vectorpool,
						       transcriptpool,hitlistpool),paths);
    }

    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_minus(&(*found_score),&(*end_score),/*main*/path,this,

						     query_compress_fwd,query_compress_rev,
						     queryseq,querylength,knownsplicing,/*genestrand*/0,
						     nmismatches_allowed,
						     intlistpool,uintlistpool,univcoordlistpool,
						     listpool,univdiagpool,pathpool,vectorpool,
						     transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_antisense_paths_gminus\n",List_length(unextended_antisense_paths_gminus)));
  for (p = unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_minus(&(*found_score),&(*end_score),/*main*/path,this,

						       query_compress_fwd,query_compress_rev,
						       queryseq,querylength,knownsplicing,/*genestrand*/0,
						       nmismatches_allowed,
						       intlistpool,uintlistpool,univcoordlistpool,
						       listpool,univdiagpool,pathpool,vectorpool,
						       transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_minus(&(*found_score),&(*end_score),/*main*/path,this,

						     query_compress_fwd,query_compress_rev,
						     queryseq,querylength,knownsplicing,/*genestrand*/0,
						     nmismatches_allowed,
						     intlistpool,uintlistpool,univcoordlistpool,
						     listpool,univdiagpool,pathpool,vectorpool,
						     transcriptpool,hitlistpool),paths);
    }
  }

  debug1(printf("Leaving single_read_fusion: %d fusion paths produced, found_score=%d end_score=%d\n",
		List_length(paths),*found_score,*end_score));

  return paths;
}
#endif


static void
convert_trpaths (int *found_score, int *end_score,

		 List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		 List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		 List_T sense_trpaths, List_T antisense_trpaths, bool first_read_p,

		 T this, Knownsplicing_T knownsplicing,
		 Shortread_T queryseq, int querylength,
		 Compress_T query_compress_fwd, Compress_T query_compress_rev,

		 Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		 unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		 int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

		 Uintlistpool_T uintlistpool, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  *sense_paths_gplus = *sense_paths_gminus = (List_T) NULL;
  Trpath_convert_sense(&(*found_score),&(*end_score),
		       &(*sense_paths_gplus),&(*sense_paths_gminus),
		       /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,
		       sense_trpaths,first_read_p,
		       queryseq,querylength,this,knownsplicing,
		       query_compress_fwd,query_compress_rev,

		       novel_univdiags_alloc,novel_diagonals_alloc,
		       localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
		       localdb,localdb_nmismatches_allowed,

		       intlistpool,uintlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool,
		       /*crossover_sense_p*/false);

  *antisense_paths_gplus = *antisense_paths_gminus = (List_T) NULL;
  Trpath_convert_antisense(&(*found_score),&(*end_score),
			   &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			   /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			   antisense_trpaths,first_read_p,
			   queryseq,querylength,this,knownsplicing,
			   query_compress_fwd,query_compress_rev,

			   novel_univdiags_alloc,novel_diagonals_alloc,
			   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			   localdb,localdb_nmismatches_allowed,

			   intlistpool,uintlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool,
			   /*crossover_sense_p*/false);

  *sense_paths_gplus = Path_consolidate(*sense_paths_gplus,queryseq,
					query_compress_fwd,query_compress_rev,
					uintlistpool,intlistpool,univcoordlistpool,
					listpool,pathpool,transcriptpool,hitlistpool);
  *sense_paths_gminus = Path_consolidate(*sense_paths_gminus,queryseq,
					query_compress_fwd,query_compress_rev,
					uintlistpool,intlistpool,univcoordlistpool,
					listpool,pathpool,transcriptpool,hitlistpool);
  *antisense_paths_gplus = Path_consolidate(*antisense_paths_gplus,queryseq,
					   query_compress_fwd,query_compress_rev,
					   uintlistpool,intlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,hitlistpool);
  *antisense_paths_gminus = Path_consolidate(*antisense_paths_gminus,queryseq,
					    query_compress_fwd,query_compress_rev,
					    uintlistpool,intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,hitlistpool);

  return;
}
    

static bool
no_converted_paths_p (List_T sense_paths_gplus, List_T sense_paths_gminus,
		      List_T antisense_paths_gplus, List_T antisense_paths_gminus) {
  return (sense_paths_gplus == NULL && sense_paths_gminus == NULL &&
	  antisense_paths_gplus == NULL && antisense_paths_gminus == NULL);
}


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
segment_apply_quality_filter (T this, Shortread_T queryseq, int querylength) {
  char *quality_string;
  int query_lastpos, querypos, nfiltered = 0;

  if (segment_quality_filter_enabled_p() == false) {
    return;
  } else if ((quality_string = Shortread_quality_string(queryseq)) == NULL) {
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


static bool
single_read_tr_paths (int *found_score, int *end_score, int sufficient_score, Method_T *last_method,

		      List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		      T this, int genestrand,
		      
		      Shortread_T queryseq, int querylength,
		      Knownsplicing_T knownsplicing, 
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,
		      int nmismatches_allowed,
			 
		      int *mismatch_positions_alloc, 
		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		      int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

		      Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		      Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Trpath_T trpath;
  int tstart, tend;
  int index;

  Trnum_T trnum;
  Chrnum_T chrnum;
  int transcript_genestrand;
  bool partialp;
  Method_T method_goal;


  /* Univcoord_T *univdiagonals_gplus, *univdiagonals_gminus; */
  /* int nunivdiagonals_gplus, nunivdiagonals_gminus, i; */


  *sense_paths_gplus = NULL;
  *sense_paths_gminus = NULL;
  *antisense_paths_gplus = NULL;
  *antisense_paths_gminus = NULL;
    
  debug1(printf("Entered single_read_tr_paths: found_score=%d end_score=%d sufficient_score=%d transcriptome_align_p=%d last_method=%d\n",
		*found_score,*end_score,sufficient_score,transcriptome_align_p,*last_method));

  if (transcriptome_align_p == true) {
    /* A.  Transcriptome search using univdiagonals */
    method_goal = single_read_trdiagonal_method_goal();
    while (tr_chain_search_p == false &&
	   *found_score > sufficient_score &&
	   single_read_trdiagonal_method_available_p(*last_method,method_goal)) {
      /* Append results to lists in Stage1_T object */
      *last_method = single_read_next_method_trdiagonal(*last_method,this,querylength,
							query_compress_fwd,query_compress_rev,
							/*first_read_p*/true);
      tstart = 0;
      tend = querylength;
      for (index = 0; index < this->n_sense_trdiagonals; index++) {
	if (this->sense_tstarts != NULL) {
	  tstart = this->sense_tstarts[index];
	  tend = this->sense_tends[index];
	}

	trnum = this->sense_trnums[index];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&(*end_score),&partialp,sufficient_score,
						   /*trdiagonal*/this->_sense_trdiagonals[index],
						   tstart,tend,trnum,/*troffset*/this->sense_troffsets[index],
						   /*trhigh*/this->sense_trhighs[index],
						   /*query_compress_tr*/query_compress_fwd,/*tplusp*/true,querylength,
						   mismatch_positions_alloc,chrnum,transcript_genestrand,
						   /*paired_end_p*/false,/*first_read_p*/true,this->indelinfo,
						   intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						   *last_method)) != NULL) {
	  if (partialp == true) {
	    this->partial_sense_trpaths = Hitlist_push(this->partial_sense_trpaths,hitlistpool,(void *) trpath
					       hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this->sense_trpaths = Hitlist_push(this->sense_trpaths,hitlistpool,(void *) trpath
					       hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
      }

      tstart = 0;
      tend = querylength;
      for (index = 0; index < this->n_antisense_trdiagonals; index++) {
	if (this->antisense_tstarts != NULL) {
	  tstart = this->antisense_tstarts[index];
	  tend = this->antisense_tends[index];
	}

	trnum = this->antisense_trnums[index];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&(*end_score),&partialp,sufficient_score,
						   /*trdiagonal*/this->_antisense_trdiagonals[index],
						   tstart,tend,trnum,/*troffset*/this->antisense_troffsets[index],
						   /*trhigh*/this->antisense_trhighs[index],
						   /*query_compress_tr*/query_compress_rev,/*tplusp*/false,querylength,
						   mismatch_positions_alloc,chrnum,transcript_genestrand,
						   /*paired_end_p*/false,/*first_read_p*/true,this->indelinfo,
						   intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						   *last_method)) != NULL) {
	  if (partialp == true) {
	    this->partial_antisense_trpaths = Hitlist_push(this->partial_antisense_trpaths,hitlistpool,(void *) trpath
							   hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this->antisense_trpaths = Hitlist_push(this->antisense_trpaths,hitlistpool,(void *) trpath
						   hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
      }
    }

    /* Should be handled by Stage1hr_free */
    /* FREE_ALIGN(this->_sense_trdiagonals); */
    /* FREE_ALIGN(this->_antisense_trdiagonals); */

    debug1(printf("single_read_tr_paths Section A done: found_score=%d end_score=%d sense_trpaths=%d antisense_trpaths=%d partial_sense=%d partial_antisense=%d\n",
		  *found_score,*end_score,List_length(this->sense_trpaths),List_length(this->antisense_trpaths),
		  List_length(this->partial_sense_trpaths),List_length(this->partial_antisense_trpaths)));

    if (*found_score <= sufficient_score) {
      convert_trpaths(&(*found_score),&(*end_score),
		      
		      &(*sense_paths_gplus),&(*sense_paths_gminus),
		      &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		      this->sense_trpaths,this->antisense_trpaths,/*first_read_p*/true,
		      this,knownsplicing,queryseq,querylength,
		      query_compress_fwd,query_compress_rev,

		      novel_univdiags_alloc,novel_diagonals_alloc,
		      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
		      localdb,localdb_nmismatches_allowed,

		      uintlistpool,intlistpool,univcoordlistpool,listpool,
		      pathpool,vectorpool,transcriptpool,hitlistpool);
      return true;
    }



    /* B.  Transcriptome search using trpaths */
    method_goal = single_read_tr_method_goal();
    while (*found_score > sufficient_score &&
	   single_read_tr_method_available_p(*last_method,method_goal)) {
      /* Append results to lists in Stage1_T object */
      *last_method = single_read_next_method_tr(&(*found_score),&(*end_score),sufficient_score,*last_method,
					     
						&this->partial_sense_trpaths,&this->partial_antisense_trpaths,
						&this->sense_trpaths,&this->antisense_trpaths,
					     
						this,genestrand,querylength,
						mismatch_positions_alloc,
						Shortread_queryuc_ptr(queryseq),Shortread_queryrc(queryseq),
						query_compress_fwd,query_compress_rev,
						nmismatches_allowed,
						trdiagpool,intlistpool,uintlistpool,
						listpool,trpathpool,pathpool,hitlistpool,
						/*paired_end_p*/false,/*first_read_p*/true,/*appendp*/true);
    }

    debug1(printf("single_read_tr_paths Section B done: found_score=%d end_score=%d sense_trpaths=%d antisense_trpaths=%d partial_sense=%d partial_antisense=%d refine_partial_tr_p=%d\n",
		  *found_score,*end_score,List_length(this->sense_trpaths),List_length(this->antisense_trpaths),
		  List_length(this->partial_sense_trpaths),List_length(this->partial_antisense_trpaths),refine_partial_tr_p));

    if (refine_partial_tr_p == true || *found_score <= sufficient_score) {
      convert_trpaths(&(*found_score),&(*end_score),

		      &(*sense_paths_gplus),&(*sense_paths_gminus),
		      &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		      this->sense_trpaths,this->antisense_trpaths,/*first_read_p*/true,
		      this,knownsplicing,queryseq,querylength,
		      query_compress_fwd,query_compress_rev,

		      novel_univdiags_alloc,novel_diagonals_alloc,
		      localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
		      localdb,localdb_nmismatches_allowed,

		      uintlistpool,intlistpool,univcoordlistpool,listpool,
		      pathpool,vectorpool,transcriptpool,hitlistpool);
      return true;
    }
    

    /* C.  Prep for genome search.  Convert all partial and complete trpaths to paths */
    debug1(printf("single_read_tr_paths Section C: converting all trpaths (sense_full=%d sense_partial=%d anti_full=%d anti_partial=%d) for genome refinement\n",
		  List_length(this->sense_trpaths),List_length(this->partial_sense_trpaths),
		  List_length(this->antisense_trpaths),List_length(this->partial_antisense_trpaths)));
    *sense_paths_gplus = *sense_paths_gminus = (List_T) NULL;
    Trpath_convert_sense(&(*found_score),&(*end_score),
			 &(*sense_paths_gplus),&(*sense_paths_gminus),
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 this->partial_sense_trpaths,/*first_read_p*/true,
			 queryseq,querylength,this,knownsplicing,
			 query_compress_fwd,query_compress_rev,

			 novel_univdiags_alloc,novel_diagonals_alloc,
			 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			 localdb,localdb_nmismatches_allowed,

			 intlistpool,uintlistpool,univcoordlistpool,listpool,
			 pathpool,vectorpool,transcriptpool,hitlistpool,
			 /*crossover_sense_p*/false);

    Trpath_convert_sense(&(*found_score),&(*end_score),
			 &(*sense_paths_gplus),&(*sense_paths_gminus),
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 this->sense_trpaths,/*first_read_p*/true,
			 queryseq,querylength,this,knownsplicing,
			 query_compress_fwd,query_compress_rev,

			 novel_univdiags_alloc,novel_diagonals_alloc,
			 localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			 localdb,localdb_nmismatches_allowed,

			 intlistpool,uintlistpool,univcoordlistpool,listpool,
			 pathpool,vectorpool,transcriptpool,hitlistpool,
			 /*crossover_sense_p*/false);

    *antisense_paths_gplus = *antisense_paths_gminus = (List_T) NULL;
    Trpath_convert_antisense(&(*found_score),&(*end_score),
			     &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     this->partial_antisense_trpaths,/*first_read_p*/true,
			     queryseq,querylength,this,knownsplicing,
			     query_compress_fwd,query_compress_rev,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool,
			     /*crossover_sense_p*/false);

    Trpath_convert_antisense(&(*found_score),&(*end_score),
			     &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     this->antisense_trpaths,/*first_read_p*/true,
			     queryseq,querylength,this,knownsplicing,
			     query_compress_fwd,query_compress_rev,

			     novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool,
			     /*crossover_sense_p*/false);


    if (refine_partial_tr_p == true || *found_score <= sufficient_score ||
	(se_keep_trpath_candidates_p == true &&
	 no_converted_paths_p(*sense_paths_gplus,*sense_paths_gminus,
			      *antisense_paths_gplus,*antisense_paths_gminus) == false)) {
      *sense_paths_gplus = Path_consolidate(*sense_paths_gplus,queryseq,
					    query_compress_fwd,query_compress_rev,
					    uintlistpool,intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,hitlistpool);
      *sense_paths_gminus = Path_consolidate(*sense_paths_gminus,queryseq,
					     query_compress_fwd,query_compress_rev,
					     uintlistpool,intlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,hitlistpool);
      *antisense_paths_gplus = Path_consolidate(*antisense_paths_gplus,queryseq,
						query_compress_fwd,query_compress_rev,
						uintlistpool,intlistpool,univcoordlistpool,
						listpool,pathpool,transcriptpool,hitlistpool);
      *antisense_paths_gminus = Path_consolidate(*antisense_paths_gminus,queryseq,
						 query_compress_fwd,query_compress_rev,
						 uintlistpool,intlistpool,univcoordlistpool,
						 listpool,pathpool,transcriptpool,hitlistpool);
      debug1(printf("single_read_tr_paths Section C: post-consolidation sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d\n",
		    List_length(*sense_paths_gplus),List_length(*sense_paths_gminus),
		    List_length(*antisense_paths_gplus),List_length(*antisense_paths_gminus)));
      if (refine_partial_tr_p == true || *found_score <= sufficient_score) {
	return true;
      } else {
	debug1(printf("single_read_tr_paths Section C: preserving converted trpaths for genome fallback\n"));
	return false;
      }
    }

    /* Previously solved unsolved paths here, but now Trpath_convert procedures solve them */
  }

  return false;
}


struct Segment_chain_apply_T {
  int *found_score;
  int *end_score;
  List_T *unextended_sense_paths;
  List_T *unextended_antisense_paths;
  List_T *sense_paths;
  List_T *antisense_paths;

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
  int sufficient_score;

  Stage1_T stage1;
  Knownsplicing_T knownsplicing;
  Knownindels_T knownindels;
  Compress_T query_compress_fwd;
  Compress_T query_compress_rev;
  Compress_T query_compress;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;
  bool plusp;
  int genestrand;
  bool leftp;

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

  Path_segment_cache_T *segment_cache;
  struct Segment_chain_seen_T *segment_chain_seen;
  int *segment_chaini;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long applied_chains;
  unsigned long long skipped_exact_chains;
  unsigned long long skipped_subtrees;
  unsigned long long total_solve_us;
  unsigned long long max_solve_us;
#endif
};

static unsigned long long
segment_chain_seen_mix (unsigned long long hash, unsigned long long value) {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
  return hash;
}

static struct Segment_chain_seen_T *
segment_chain_seen_new (void) {
  return (struct Segment_chain_seen_T *) CALLOC(1,sizeof(struct Segment_chain_seen_T));
}

static void
segment_chain_seen_free (struct Segment_chain_seen_T **old) {
  int i;

  if (*old != (struct Segment_chain_seen_T *) NULL) {
    for (i = 0; i < SEGMENT_CHAIN_SEEN_SIZE; i++) {
      if ((*old)->entries[i].occupiedp == true) {
	FREE((*old)->entries[i].chain_elts);
      }
    }
    FREE(*old);
    *old = (struct Segment_chain_seen_T *) NULL;
  }

  return;
}

static unsigned long long
segment_chain_seen_hash (int *chain_length, struct Segment_chain_apply_T *this, List_T chain) {
  List_T p;
  Univdiag_T univdiag;
  unsigned long long hash = 1469598103934665603ULL;
  int n = 0;

  hash = segment_chain_seen_mix(hash,(unsigned long long) this->leftp);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->plusp);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->genestrand);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->querylength);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->max_nblocks);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->localdb_nmismatches_allowed);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->chrnum);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->chroffset);
  hash = segment_chain_seen_mix(hash,(unsigned long long) this->chrhigh);

  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    hash = segment_chain_seen_mix(hash,(unsigned long long) univdiag->qstart);
    hash = segment_chain_seen_mix(hash,(unsigned long long) univdiag->qend);
    hash = segment_chain_seen_mix(hash,(unsigned long long) univdiag->nmismatches);
    hash = segment_chain_seen_mix(hash,(unsigned long long) univdiag->univdiagonal);
    n++;
  }

  *chain_length = n;
  return segment_chain_seen_mix(hash,(unsigned long long) n);
}

static bool
segment_chain_seen_entry_matches_p (struct Segment_chain_seen_entry_T *entry,
				    struct Segment_chain_apply_T *this, List_T chain,
				    unsigned long long hash, int chain_length) {
  List_T p;
  Univdiag_T univdiag;
  int i;

  if (entry->hash != hash ||
      entry->leftp != this->leftp ||
      entry->plusp != this->plusp ||
      entry->genestrand != this->genestrand ||
      entry->querylength != this->querylength ||
      entry->max_nblocks != this->max_nblocks ||
      entry->localdb_nmismatches_allowed != this->localdb_nmismatches_allowed ||
      entry->chrnum != this->chrnum ||
      entry->chroffset != this->chroffset ||
      entry->chrhigh != this->chrhigh ||
      entry->chain_length != chain_length) {
    return false;
  }

  for (i = 0, p = chain; i < chain_length && p != NULL; i++, p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    if (entry->chain_elts[i].qstart != univdiag->qstart ||
	entry->chain_elts[i].qend != univdiag->qend ||
	entry->chain_elts[i].nmismatches != univdiag->nmismatches ||
	entry->chain_elts[i].univdiagonal != univdiag->univdiagonal) {
      return false;
    }
  }

  return (i == chain_length && p == (List_T) NULL) ? true : false;
}

static void
segment_chain_seen_entry_fill (struct Segment_chain_seen_entry_T *entry,
			       struct Segment_chain_apply_T *this, List_T chain,
			       unsigned long long hash, int chain_length) {
  List_T p;
  Univdiag_T univdiag;
  int i;

  entry->occupiedp = true;
  entry->hash = hash;
  entry->leftp = this->leftp;
  entry->plusp = this->plusp;
  entry->genestrand = this->genestrand;
  entry->querylength = this->querylength;
  entry->max_nblocks = this->max_nblocks;
  entry->localdb_nmismatches_allowed = this->localdb_nmismatches_allowed;
  entry->chrnum = this->chrnum;
  entry->chroffset = this->chroffset;
  entry->chrhigh = this->chrhigh;
  entry->chain_length = chain_length;
  entry->chain_elts = (struct Segment_chain_key_elt_T *)
    MALLOC(chain_length * sizeof(struct Segment_chain_key_elt_T));

  for (i = 0, p = chain; i < chain_length && p != NULL; i++, p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    entry->chain_elts[i].qstart = univdiag->qstart;
    entry->chain_elts[i].qend = univdiag->qend;
    entry->chain_elts[i].nmismatches = univdiag->nmismatches;
    entry->chain_elts[i].univdiagonal = univdiag->univdiagonal;
  }

  return;
}

static bool
segment_chain_seen_or_mark (struct Segment_chain_seen_T *seen,
			    struct Segment_chain_apply_T *this, List_T chain) {
  struct Segment_chain_seen_entry_T *entry;
  unsigned long long hash;
  unsigned int hash_index;
  int chain_length, i;

  if (seen == (struct Segment_chain_seen_T *) NULL) {
    return false;
  } else if (chain == (List_T) NULL) {
    seen->null_chains++;
    return false;
  } else if (seen->nentries >= SEGMENT_CHAIN_SEEN_SIZE) {
    seen->full_misses++;
    return false;
  }

  seen->lookups++;
  hash = segment_chain_seen_hash(&chain_length,this,chain);
  if (chain_length == 0) {
    seen->null_chains++;
    return false;
  }

  hash_index = (unsigned int) (hash & SEGMENT_CHAIN_SEEN_MASK);
  for (i = 0; i < SEGMENT_CHAIN_SEEN_SIZE; i++) {
    entry = &(seen->entries[(hash_index + i) & SEGMENT_CHAIN_SEEN_MASK]);
    if (entry->occupiedp == false) {
      segment_chain_seen_entry_fill(entry,this,chain,hash,chain_length);
      seen->nentries++;
      seen->marks++;
      return false;
    } else if (segment_chain_seen_entry_matches_p(entry,this,chain,hash,chain_length) == true) {
      seen->skips++;
      return true;
    } else {
      seen->collisions++;
    }
  }

  seen->full_misses++;
  return false;
}

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
static void
segment_chain_seen_report (struct Segment_chain_seen_T *seen) {
  if (seen != (struct Segment_chain_seen_T *) NULL) {
    fprintf(stderr,
	    "[SEGMENT_CHAIN_SEEN] lookups=%llu marks=%llu skips=%llu "
	    "full_misses=%llu collisions=%llu null_chains=%llu nentries=%d size=%d\n",
	    seen->lookups,seen->marks,seen->skips,seen->full_misses,
	    seen->collisions,seen->null_chains,seen->nentries,SEGMENT_CHAIN_SEEN_SIZE);
  }

  return;
}
#endif

static bool
segment_chain_skip_chain_p (int head_qstart, int head_qend, Univcoord_T head_univdiagonal,
			    int next_qstart, int next_qend, Univcoord_T next_univdiagonal,
			    void *data) {
  struct Segment_chain_apply_T *this = (struct Segment_chain_apply_T *) data;

  if (se_segment_prefix_skip_p == false) {
    return false;
  } else if (*(this->segment_cache) == (Path_segment_cache_T) NULL) {
    return false;
  } else if (this->leftp == true) {
    if (Path_segment_cache_qstart_known_p(*(this->segment_cache),
					  head_qstart,head_qend,/*head_nmismatches*/-1,
					  head_univdiagonal,
					  next_qstart,next_qend,/*next_nmismatches*/-1,
					  next_univdiagonal,
					  this->querylength,this->localdb_nmismatches_allowed,
					  this->chrnum,this->chroffset,this->chrhigh,
					  this->plusp,this->genestrand,
					  /*first_read_p*/true,/*find_splices_p*/true) == true) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      this->skipped_subtrees++;
#endif
      return true;
    } else {
      return false;
    }
  } else {
    if (Path_segment_cache_qend_known_p(*(this->segment_cache),
					head_qstart,head_qend,/*head_nmismatches*/-1,
					head_univdiagonal,
					next_qstart,next_qend,/*next_nmismatches*/-1,
					next_univdiagonal,
					this->querylength,this->localdb_nmismatches_allowed,
					this->chrnum,this->chroffset,this->chrhigh,
					this->plusp,this->genestrand,
					/*first_read_p*/true,/*find_splices_p*/true) == true) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      this->skipped_subtrees++;
#endif
      return true;
    } else {
      return false;
    }
  }
}

static void
segment_chain_apply_chain (List_T chain, void *data) {
  struct Segment_chain_apply_T *this = (struct Segment_chain_apply_T *) data;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  int before_found_score = *(this->found_score), before_end_score = *(this->end_score);
  unsigned long long start_us = segment_stream_now_us(), elapsed_us;
#endif

  (*this->segment_chaini)++;
  if (se_segment_dedup_exact_chains_p == true &&
      segment_chain_seen_or_mark(this->segment_chain_seen,this,chain) == true) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    this->skipped_exact_chains++;
#endif
    return;
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  this->applied_chains++;
#endif

  if ((se_segment_prefix_skip_p == true ||
       se_segment_qstart_anchor_cache_p == true ||
       se_path_qend_local_cache_p == true ||
       se_emitted_path_dedup_p == true) &&
      *this->segment_chaini > se_segment_cache_threshold) {
    if (*(this->segment_cache) == (Path_segment_cache_T) NULL) {
      *(this->segment_cache) = Path_segment_cache_new(se_emitted_path_dedup_p);
    }
    if (this->leftp == true) {
      Path_solve_from_left_chain_cached(this->found_score,this->end_score,

					this->unextended_sense_paths,this->unextended_antisense_paths,
					this->sense_paths,this->antisense_paths,

					chain,this->queryseq,this->queryptr,this->querylength,
					this->mismatch_positions_alloc,
					this->novel_univdiags_alloc,this->novel_diagonals_alloc,
					this->localdb_alloc,this->localdb_tally,this->localdb_overflowp,
					this->max_nblocks,this->localdb_nmismatches_allowed,

					this->stage1,this->knownsplicing,this->knownindels,
					this->query_compress,
					this->query_compress_fwd,this->query_compress_rev,
					this->chrnum,this->chroffset,this->chrhigh,
					this->plusp,this->genestrand,
					/*paired_end_p*/false,/*first_read_p*/true,
					this->intlistpool,this->uintlistpool,this->univcoordlistpool,
					this->listpool,this->pathpool,this->transcriptpool,
					this->univdiagpool,this->vectorpool,this->hitlistpool,
					this->spliceendsgen,
					/*method*/SEGMENT1,/*find_splices_p*/true,
					/*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,
					*(this->segment_cache),
					se_segment_prefix_skip_p,
					se_segment_qstart_anchor_cache_p,
					se_path_qend_local_cache_p,
					se_emitted_path_dedup_p);
    } else {
      Path_solve_from_right_chain_cached(this->found_score,this->end_score,

					 this->unextended_sense_paths,this->unextended_antisense_paths,
					 this->sense_paths,this->antisense_paths,

					 chain,this->queryseq,this->queryptr,this->querylength,
					 this->mismatch_positions_alloc,
					 this->novel_univdiags_alloc,this->novel_diagonals_alloc,
					 this->localdb_alloc,this->localdb_tally,this->localdb_overflowp,
					 this->max_nblocks,this->localdb_nmismatches_allowed,

					 this->stage1,this->knownsplicing,this->knownindels,
					 this->query_compress,
					 this->query_compress_fwd,this->query_compress_rev,
					 this->chrnum,this->chroffset,this->chrhigh,
					 this->plusp,this->genestrand,
					 /*paired_end_p*/false,/*first_read_p*/true,
					 this->intlistpool,this->uintlistpool,this->univcoordlistpool,
					 this->listpool,this->pathpool,this->transcriptpool,
					 this->univdiagpool,this->vectorpool,this->hitlistpool,
					 this->spliceendsgen,
					 /*method*/SEGMENT1,/*find_splices_p*/true,
					 /*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,
					 *(this->segment_cache),
					 se_segment_prefix_skip_p,
					 se_segment_qstart_anchor_cache_p,
					 se_path_qend_local_cache_p,
					 se_emitted_path_dedup_p);
    }
  } else {
    if (this->leftp == true) {
      Path_solve_from_left_chain(this->found_score,this->end_score,

				 this->unextended_sense_paths,this->unextended_antisense_paths,
				 this->sense_paths,this->antisense_paths,

				 chain,this->queryseq,this->queryptr,this->querylength,
				 this->mismatch_positions_alloc,
				 this->novel_univdiags_alloc,this->novel_diagonals_alloc,
				 this->localdb_alloc,this->localdb_tally,this->localdb_overflowp,
				 this->max_nblocks,this->localdb_nmismatches_allowed,

				 this->stage1,this->knownsplicing,this->knownindels,
				 this->query_compress,
				 this->query_compress_fwd,this->query_compress_rev,
				 this->chrnum,this->chroffset,this->chrhigh,
				 this->plusp,this->genestrand,
				 /*paired_end_p*/false,/*first_read_p*/true,
				 this->intlistpool,this->uintlistpool,this->univcoordlistpool,
				 this->listpool,this->pathpool,this->transcriptpool,
				 this->univdiagpool,this->vectorpool,this->hitlistpool,
				 this->spliceendsgen,
				 /*method*/SEGMENT1,/*find_splices_p*/true);
    } else {
      Path_solve_from_right_chain(this->found_score,this->end_score,

				  this->unextended_sense_paths,this->unextended_antisense_paths,
				  this->sense_paths,this->antisense_paths,

				  chain,this->queryseq,this->queryptr,this->querylength,
				  this->mismatch_positions_alloc,
				  this->novel_univdiags_alloc,this->novel_diagonals_alloc,
				  this->localdb_alloc,this->localdb_tally,this->localdb_overflowp,
				  this->max_nblocks,this->localdb_nmismatches_allowed,

				  this->stage1,this->knownsplicing,this->knownindels,
				  this->query_compress,
				  this->query_compress_fwd,this->query_compress_rev,
				  this->chrnum,this->chroffset,this->chrhigh,
				  this->plusp,this->genestrand,
				  /*paired_end_p*/false,/*first_read_p*/true,
				  this->intlistpool,this->uintlistpool,this->univcoordlistpool,
				  this->listpool,this->pathpool,this->transcriptpool,
				  this->univdiagpool,this->vectorpool,this->hitlistpool,
				  this->spliceendsgen,
				  /*method*/SEGMENT1,/*find_splices_p*/true);
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  elapsed_us = segment_stream_now_us() - start_us;
  this->total_solve_us += elapsed_us;
  if (elapsed_us > this->max_solve_us) {
    this->max_solve_us = elapsed_us;
  }
  if (*(this->found_score) != before_found_score || *(this->end_score) != before_end_score) {
    fprintf(stderr,
	    "[SEGMENT_STREAM] chain=%llu elapsed_us=%llu found_score=%d->%d end_score=%d->%d sufficient=%d\n",
	    this->applied_chains,elapsed_us,before_found_score,*(this->found_score),
	    before_end_score,*(this->end_score),this->sufficient_score);
  } else if (elapsed_us >= 100000ULL) {
    fprintf(stderr,
	    "[SEGMENT_STREAM] slow_chain=%llu elapsed_us=%llu found_score=%d end_score=%d sufficient=%d\n",
	    this->applied_chains,elapsed_us,*(this->found_score),
	    *(this->end_score),this->sufficient_score);
  }
#endif

  return;
}


/* For single-end reads, we will just call
   Path_solve_from_univdiagonals on most prevalent univdiagonals.
   For paired-end reads, we will return all univdiagonals (and most
   prevalent) and try to anchor */

/* TODO: Need to return Auxinfo_T arrays, which have the paths */
static bool
single_read_gen_paths (int *found_score, int *end_score, Method_T *last_method,

		       List_T *unextended_sense_paths_gplus, List_T *unextended_sense_paths_gminus,
		       List_T *unextended_antisense_paths_gplus, List_T *unextended_antisense_paths_gminus,

		       List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		       List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		       T this, int genestrand,
		 
		       Shortread_T queryseq, char *queryuc_ptr, char *queryrc, int querylength,
		       Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,

		       int *mismatch_positions_alloc, 
		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		       int max_nblocks,

		       int localdb_nmismatches_allowed, EF64_T repetitive_ef64,
		       
		       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		       Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		       Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		       Spliceendsgen_T spliceendsgen) {
  
  int sufficient_score = querylength/20;
  int max_nsplices = querylength/EXON_LENGTH + 1;

  bool any_imperfect_ends_p = false;
  bool gen_exact_p = Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXACT);
  bool gen_ext_p = Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_EXT);
  bool gen_segment_p = Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_SEGMENT);
  int kmer_querystart, kmer_queryend;
  int ii, jj, i, j;

  int **indices_by_nsplices_gplus, **indices_by_nsplices_gminus, *indices;
  int *nindices_by_nsplices_gplus, *nindices_by_nsplices_gminus;
  int nsplices;

  int nunivdiagonals_gplus, nunivdiagonals_gminus;
  int total_npositions_plus, total_npositions_minus;

  Univcoord_T *_univdiagonals_gplus, *_univdiagonals_gminus, univdiagonal;
  Auxinfo_T *auxinfo_gplus = NULL, *auxinfo_gminus = NULL, auxinfo;
  Pathstore_T *pathstores_gplus = NULL, *pathstores_gminus = NULL, pathstore;

  List_T chains, chain, p;
  List_T ext_seen_chains;
  int ext_chain_rank;
  Path_segment_cache_T segment_cache = NULL;
  struct Segment_chain_seen_T *segment_chain_seen = (struct Segment_chain_seen_T *) NULL;
  Singleelt_T *eltarray;
  int nexhaustive;
  int segment_chaini = 0;
  int extension_chaini = 0;
  int *order_gplus = NULL, *order_gminus = NULL;
  int *segment_counts_gplus = NULL, *segment_counts_gminus = NULL;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long stage_start_us;
#endif

#ifdef INDIVIDUAL_CHRINFO
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
#endif

  
  *sense_paths_gplus = NULL;
  *sense_paths_gminus = NULL;
  *antisense_paths_gplus = NULL;
  *antisense_paths_gminus = NULL;

  if (gen_exact_p == false && gen_ext_p == false && gen_segment_p == false) {
    return false;
  }

  /* 1. Exact */
  if (gen_exact_p == true) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  stage_start_us = segment_stream_now_us();
#endif
  debug(printf("Single Read: 1.  Running Kmer_exact1\n"));
  Stage1_init_end_gen(&kmer_querystart,&kmer_queryend,this,querylength,genestrand);
  Kmer_exact1(&_univdiagonals_gplus,&auxinfo_gplus,
	      &pathstores_gplus,&nunivdiagonals_gplus,
	      &_univdiagonals_gminus,&auxinfo_gminus,
	      &pathstores_gminus,&nunivdiagonals_gminus,
	      this,kmer_querystart,kmer_queryend,querylength,
	      auxinfopool,pathstorepool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf(" => %d plus and %d minus univdiagonals\n",nunivdiagonals_gplus,nunivdiagonals_gminus));
  debug(printf("Kmer_exact1 returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = KMER_EXACT1;
  any_imperfect_ends_p = false;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

    if (univdiagonal < (Univcoord_T) querylength) {
      /* Skip */
    } else {
#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength,univdiagonal);
#endif
      if (Path_solve_exact(&(*found_score),&(*end_score),
			   
			   &(*sense_paths_gplus),&(*antisense_paths_gplus),

			   univdiagonal,auxinfo,querylength,
			   /*max_nmismatches*/-1,
			   /*plusp*/true,/*first_read_p*/true,genestrand,
			   /*query_compress*/query_compress_fwd,
			   query_compress_fwd,query_compress_rev,
			   queryseq,queryuc_ptr,queryrc,
			   /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			   /*chrhigh*/pathstore->chrhigh,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			   /*method*/KMER_EXACT1) == false) {
	debug(printf("Imperfect end\n"));
	any_imperfect_ends_p = true;
      }
    }
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

    if (univdiagonal < (Univcoord_T) querylength) {
      /* Skip */
    } else {
#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength,univdiagonal);
#endif
      if (Path_solve_exact(&(*found_score),&(*end_score),

			   &(*sense_paths_gminus),&(*antisense_paths_gminus),

			   univdiagonal,auxinfo,querylength,
			   /*max_nmismatches*/-1,
			   /*plusp*/false,/*first_read_p*/true,genestrand,
			   /*query_compress*/query_compress_rev,
			   query_compress_fwd,query_compress_rev,
			   queryseq,queryuc_ptr,queryrc,
			   /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			   /*chrhigh*/pathstore->chrhigh,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			   /*method*/KMER_EXACT1) == false) {
	debug(printf("Imperfect end\n"));
	any_imperfect_ends_p = true;
      }
    }
  }
    
  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool,listpool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool,listpool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[STAGE1_GEN_TIMING] exact_us=%llu found_score=%d end_score=%d\n",
	  segment_stream_now_us() - stage_start_us,*found_score,*end_score);
#endif


  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    return true;
  }
  }


  /* 2. Extension search */
  if (gen_ext_p == false) {
    if (gen_segment_p == true) {
      if (this->all_oligos_gen_filledp == false) {
	Stage1_fill_all_oligos_gen(this,querylength,genestrand);
      }
      goto segment_search;
    } else {
      return true;
    }
  }
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  stage_start_us = segment_stream_now_us();
#endif
  debug(printf("Single Read: 2.  Running Extension search\n"));
  Stage1_fill_all_oligos_gen(this,querylength,genestrand);
  Extension_search(&indices_by_nsplices_gplus,&nindices_by_nsplices_gplus,
		   &indices_by_nsplices_gminus,&nindices_by_nsplices_gminus,

		   &_univdiagonals_gplus,&auxinfo_gplus,
		   &pathstores_gplus,&nunivdiagonals_gplus,
		   &_univdiagonals_gminus,&auxinfo_gminus,
		   &pathstores_gminus,&nunivdiagonals_gminus,

		   this,query_compress_fwd,query_compress_rev,
		   querylength,max_nsplices,/*first_read_p*/true,/*retain_terminal_p*/false,
		   se_ext_equiv_prev_prune_p,
		   segmenteltpool,auxinfopool,pathstorepool,
		   univcoordlistpool,univdiagpool,intlistpool,listpool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf("Extension_search returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = EXT;

  if (se_path_qend_local_cache_p == true || se_emitted_path_dedup_p == true) {
    segment_cache = Path_segment_cache_new(se_emitted_path_dedup_p);
  }

  nsplices = 0;
  while (nsplices <= max_nsplices && *found_score > sufficient_score) {
    /* Finish the full plus/minus bucket once started.  Stopping inside
       the current nsplices bucket can skip a later same-bucket chain that
       improves a merely sufficient alignment. */
    /* plus */
    indices = indices_by_nsplices_gplus[nsplices];
    for (i = 0; i < nindices_by_nsplices_gplus[nsplices]; i++) {
      ii = indices[i];
      univdiagonal = _univdiagonals_gplus[ii];
      debug(printf("Extension plus univdiagonal %u\n",univdiagonal));
      auxinfo = auxinfo_gplus[ii];
      pathstore = pathstores_gplus[ii];

#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength + auxinfo->qstart,
			   univdiagonal - querylength + auxinfo->qend);
#endif

      assert(auxinfo->chaindir == CHAIN_LEFT);
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
      single_debug_print_chain_catalog("plus",nsplices,i + 1,univdiagonal,auxinfo,
				       auxinfo->chains_by_nsplices[nsplices]);
#endif
      ext_seen_chains = (List_T) NULL;
      ext_chain_rank = 0;
      for (p = auxinfo->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
	chain = (List_T) List_head(p);
	if (se_ext_dedup_exact_chains_p == true) {
	  ext_chain_rank++;
	  if (single_ext_chain_seen_p(ext_seen_chains,chain) == true) {
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	    fprintf(stderr,
		    "[EXT_CHAIN_SKIP_DUP] dir=plus nsplices=%d candidate=%d chain_rank=%d "
		    "chain_hash=%016llx found=%d sufficient=%d\n",
		    nsplices,i + 1,ext_chain_rank,single_debug_chain_hash(chain),
		    *found_score,sufficient_score);
	    fflush(stderr);
#endif
	    continue;
	  } else {
	    ext_seen_chains = List_push(ext_seen_chains,(void *) chain);
	  }
	}
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	{
	  int chain_id = ++extension_chaini;
	  int found_before = *found_score, end_before = *end_score;
	  int us_before = List_length(*unextended_sense_paths_gplus);
	  int ua_before = List_length(*unextended_antisense_paths_gplus);
	  int s_before = List_length(*sense_paths_gplus);
	  int a_before = List_length(*antisense_paths_gplus);
	  unsigned long long start_us;

	  single_debug_print_chain_summary("BEGIN","plus",nsplices,chain_id,
					   i + 1,nindices_by_nsplices_gplus[nsplices],
					   univdiagonal,auxinfo,pathstore,
					   auxinfo->chains_by_nsplices[nsplices],p,
					   found_before,end_before,
					   *unextended_sense_paths_gplus,
					   *unextended_antisense_paths_gplus,
					   *sense_paths_gplus,*antisense_paths_gplus);
	  start_us = stage1_single_now_us();
#endif
	Path_solve_from_left_chain_cached(&(*found_score),&(*end_score),
			   
				   &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
				   &(*sense_paths_gplus),&(*antisense_paths_gplus),
				   chain,queryseq,/*queryptr*/queryuc_ptr,querylength,mismatch_positions_alloc,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				   /*stage1*/this,knownsplicing,knownindels,
				   /*query_compress*/query_compress_fwd,query_compress_fwd,query_compress_rev,
				   /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				   /*chrhigh*/pathstore->chrhigh,/*plusp*/true,genestrand,
				   /*paired_end_p*/false,/*first_read_p*/true,
				   intlistpool,uintlistpool,univcoordlistpool,listpool,
				   pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				   /*method*/EXT,/*find_splices_p*/true,
				   /*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,segment_cache,
				   /*segment_prefix_skip_p*/false,
				   /*segment_qstart_anchor_cache_p*/false,
				   se_path_qend_local_cache_p,
				   se_emitted_path_dedup_p);
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	  single_debug_print_chain_result("plus",nsplices,chain_id,
					  stage1_single_now_us() - start_us,
					  found_before,end_before,*found_score,*end_score,
					  us_before,ua_before,s_before,a_before,
					  *unextended_sense_paths_gplus,
					  *unextended_antisense_paths_gplus,
					  *sense_paths_gplus,*antisense_paths_gplus);
	}
#endif
      }
      if (se_ext_dedup_exact_chains_p == true) {
	List_free(&ext_seen_chains);
      }
    }

    /* minus */
    indices = indices_by_nsplices_gminus[nsplices];
    for (i = 0; i < nindices_by_nsplices_gminus[nsplices]; i++) {
      ii = indices[i];
      univdiagonal = _univdiagonals_gminus[ii];
      debug(printf("Extension minus univdiagonal %u\n",univdiagonal));
      auxinfo = auxinfo_gminus[ii];
      pathstore = pathstores_gminus[ii];

#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength + auxinfo->qstart,
			   univdiagonal - querylength + auxinfo->qend);
#endif

      assert(auxinfo->chaindir == CHAIN_RIGHT);
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
      single_debug_print_chain_catalog("minus",nsplices,i + 1,univdiagonal,auxinfo,
				       auxinfo->chains_by_nsplices[nsplices]);
#endif
      ext_seen_chains = (List_T) NULL;
      ext_chain_rank = 0;
      for (p = auxinfo->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
	chain = (List_T) List_head(p);
	if (se_ext_dedup_exact_chains_p == true) {
	  ext_chain_rank++;
	  if (single_ext_chain_seen_p(ext_seen_chains,chain) == true) {
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	    fprintf(stderr,
		    "[EXT_CHAIN_SKIP_DUP] dir=minus nsplices=%d candidate=%d chain_rank=%d "
		    "chain_hash=%016llx found=%d sufficient=%d\n",
		    nsplices,i + 1,ext_chain_rank,single_debug_chain_hash(chain),
		    *found_score,sufficient_score);
	    fflush(stderr);
#endif
	    continue;
	  } else {
	    ext_seen_chains = List_push(ext_seen_chains,(void *) chain);
	  }
	}
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	{
	  int chain_id = ++extension_chaini;
	  int found_before = *found_score, end_before = *end_score;
	  int us_before = List_length(*unextended_sense_paths_gminus);
	  int ua_before = List_length(*unextended_antisense_paths_gminus);
	  int s_before = List_length(*sense_paths_gminus);
	  int a_before = List_length(*antisense_paths_gminus);
	  unsigned long long start_us;

	  single_debug_print_chain_summary("BEGIN","minus",nsplices,chain_id,
					   i + 1,nindices_by_nsplices_gminus[nsplices],
					   univdiagonal,auxinfo,pathstore,
					   auxinfo->chains_by_nsplices[nsplices],p,
					   found_before,end_before,
					   *unextended_sense_paths_gminus,
					   *unextended_antisense_paths_gminus,
					   *sense_paths_gminus,*antisense_paths_gminus);
	  start_us = stage1_single_now_us();
#endif
	Path_solve_from_right_chain_cached(&(*found_score),&(*end_score),

				    &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
				    &(*sense_paths_gminus),&(*antisense_paths_gminus),
				    chain,queryseq,/*queryptr*/queryrc,querylength,mismatch_positions_alloc,
				    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				    localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,
			      
				    /*stage1*/this,knownsplicing,knownindels,
				    /*query_compress*/query_compress_rev,query_compress_fwd,query_compress_rev,
				    /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				    /*chrhigh*/pathstore->chrhigh,/*plusp*/false,genestrand,
				    /*paired_end_p*/false,/*first_read_p*/true,
				    intlistpool,uintlistpool,univcoordlistpool,listpool,
				    pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				    /*method*/EXT,/*find_splices_p*/true,
				    /*ext_anchor_cache*/(Path_ext_anchor_cache_T) NULL,segment_cache,
				    /*segment_prefix_skip_p*/false,
				    /*segment_qstart_anchor_cache_p*/false,
				    se_path_qend_local_cache_p,
				    se_emitted_path_dedup_p);
#if defined(DEBUG1) || defined(DEBUG_SINGLE_DEEP)
	  single_debug_print_chain_result("minus",nsplices,chain_id,
					  stage1_single_now_us() - start_us,
					  found_before,end_before,*found_score,*end_score,
					  us_before,ua_before,s_before,a_before,
					  *unextended_sense_paths_gminus,
					  *unextended_antisense_paths_gminus,
					  *sense_paths_gminus,*antisense_paths_gminus);
	}
#endif
      }
      if (se_ext_dedup_exact_chains_p == true) {
	List_free(&ext_seen_chains);
      }
    }

    debug(printf("At extension search %d splices (%d plus and %d minus), found score %d vs sufficient score %d, end_score %d\n",
		 nsplices,nindices_by_nsplices_gplus[nsplices],nindices_by_nsplices_gminus[nsplices],
		 *found_score,sufficient_score,*end_score));
    nsplices++;
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    FREE(indices_by_nsplices_gminus[nsplices]);
    FREE(indices_by_nsplices_gplus[nsplices]);
  }
  FREE(indices_by_nsplices_gminus);
  FREE(indices_by_nsplices_gplus);
  FREE(nindices_by_nsplices_gminus);
  FREE(nindices_by_nsplices_gplus);


  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool,listpool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool,listpool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[STAGE1_GEN_TIMING] extension_us=%llu found_score=%d end_score=%d\n",
	  segment_stream_now_us() - stage_start_us,*found_score,*end_score);
#endif

  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    Path_segment_cache_free(&segment_cache,intlistpool,univcoordlistpool,
			    listpool,pathpool,transcriptpool,hitlistpool);
    segment_chain_seen_free(&segment_chain_seen);
    return true;
  }

  if (gen_segment_p == false) {
    Path_segment_cache_free(&segment_cache,intlistpool,univcoordlistpool,
			    listpool,pathpool,transcriptpool,hitlistpool);
    segment_chain_seen_free(&segment_chain_seen);
    goto final_extend;
  }

  /* 3. Segment search */
 segment_search:
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  stage_start_us = segment_stream_now_us();
#endif
  debug(printf("Single Read: 3.  Running Kmer segment\n"));
  if (this->all_oligos_gen_filledp == false) {
    Stage1_fill_all_oligos_gen(this,querylength,genestrand);
  }
  assert(this->all_oligos_gen_filledp == true);
  segment_apply_quality_filter(this,queryseq,querylength);
  Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				this,querylength,genestrand);
  Stage1_filter_segment_kmer_positions_gen(this,querylength,segment_max_avg_kmer_positions,"single");
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  Stage1_print_segment_kmer_positions_gen(this,querylength,"single");
#endif
  Kmer_segment(&_univdiagonals_gplus,&auxinfo_gplus
	       ,&pathstores_gplus,&nunivdiagonals_gplus,
	       &_univdiagonals_gminus,&auxinfo_gminus,
	       &pathstores_gminus,&nunivdiagonals_gminus,

	       this,querylength,max_nsplices,/*first_read_p*/true,
	       Shortread_accession(queryseq),"single",repetitive_ef64,
	       univdiagpool,listpool,univcoordlistpool,auxinfopool,pathstorepool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf("Kmer_segment returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = SEGMENT1;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,
	  "[STAGE1_GEN_TIMING] segment_prep_us=%llu plus=%d minus=%d total_pos_plus=%d total_pos_minus=%d\n",
	  segment_stream_now_us() - stage_start_us,nunivdiagonals_gplus,nunivdiagonals_gminus,
	  total_npositions_plus,total_npositions_minus);
  stage_start_us = segment_stream_now_us();
#endif

  if (nunivdiagonals_gplus == 0) {
    i = -1;
  } else {
    segment_counts_gplus = (int *) MALLOC((nunivdiagonals_gplus + 1) * sizeof(int));
    for (i = 0; i < nunivdiagonals_gplus; i++) {
      auxinfo = auxinfo_gplus[i];
      segment_counts_gplus[i] = this->exhaustive_counts_gplus[auxinfo->elti];
    }
    order_gplus = Sedgesort_order_int(segment_counts_gplus,nunivdiagonals_gplus);
    i = nunivdiagonals_gplus - 1;
  }

  if (nunivdiagonals_gminus == 0) {
    j = -1;
  } else {
    segment_counts_gminus = (int *) MALLOC((nunivdiagonals_gminus + 1) * sizeof(int));
    for (j = 0; j < nunivdiagonals_gminus; j++) {
      auxinfo = auxinfo_gminus[j];
      segment_counts_gminus[j] = this->exhaustive_counts_gminus[auxinfo->elti];
    }
    order_gminus = Sedgesort_order_int(segment_counts_gminus,nunivdiagonals_gminus);
    j = nunivdiagonals_gminus - 1;
  }
  if (se_segment_stream_chains_p == true &&
      se_segment_dedup_exact_chains_p == true) {
    segment_chain_seen = segment_chain_seen_new();
  }
  while ((i >= 0 || j >= 0) && *found_score > sufficient_score) {
    if (j < 0 ||
	(i >= 0 && segment_counts_gplus[order_gplus[i]] > segment_counts_gminus[order_gminus[j]])) {
      /* plus */
      eltarray = this->exhaustive_eltarray_gplus;
      nexhaustive = this->nexhaustive_gplus;
      ii = order_gplus[i];
      auxinfo = auxinfo_gplus[ii];
      pathstore = pathstores_gplus[ii];

      if (se_segment_stream_chains_p == false) {
	chains = Kmer_solve_eltarray_left(eltarray,/*elti*/auxinfo->elti,nexhaustive,
					  max_nsplices,univdiagpool,listpool);
	for (p = chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_left_chain(&(*found_score),&(*end_score),

				     &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
				     &(*sense_paths_gplus),&(*antisense_paths_gplus),

				     chain,queryseq,/*queryptr*/queryuc_ptr,querylength,mismatch_positions_alloc,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				     this,knownsplicing,knownindels,
				     /*query_compress*/query_compress_fwd,query_compress_fwd,query_compress_rev,
				     /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				     /*chrhigh*/pathstore->chrhigh,/*plusp*/true,genestrand,
				     /*paired_end_p*/false,/*first_read_p*/true,
				     intlistpool,uintlistpool,univcoordlistpool,listpool,
				     pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				     /*method*/SEGMENT1,/*find_splices_p*/true);
	}

      } else {
	struct Segment_chain_apply_T segment_apply;

	segment_apply.found_score = &(*found_score);
	segment_apply.end_score = &(*end_score);
	segment_apply.unextended_sense_paths = &(*unextended_sense_paths_gplus);
	segment_apply.unextended_antisense_paths = &(*unextended_antisense_paths_gplus);
	segment_apply.sense_paths = &(*sense_paths_gplus);
	segment_apply.antisense_paths = &(*antisense_paths_gplus);

	segment_apply.queryseq = queryseq;
	segment_apply.queryptr = queryuc_ptr;
	segment_apply.querylength = querylength;
	segment_apply.mismatch_positions_alloc = mismatch_positions_alloc;
	segment_apply.novel_univdiags_alloc = novel_univdiags_alloc;
	segment_apply.novel_diagonals_alloc = novel_diagonals_alloc;
	segment_apply.localdb_alloc = localdb_alloc;
	segment_apply.localdb_tally = localdb_tally;
	segment_apply.localdb_overflowp = localdb_overflowp;
	segment_apply.max_nblocks = max_nblocks;
	segment_apply.localdb_nmismatches_allowed = localdb_nmismatches_allowed;
	segment_apply.sufficient_score = sufficient_score;

	segment_apply.stage1 = this;
	segment_apply.knownsplicing = knownsplicing;
	segment_apply.knownindels = knownindels;
	segment_apply.query_compress_fwd = query_compress_fwd;
	segment_apply.query_compress_rev = query_compress_rev;
	segment_apply.query_compress = query_compress_fwd;
	segment_apply.chrnum = pathstore->chrnum;
	segment_apply.chroffset = pathstore->chroffset;
	segment_apply.chrhigh = pathstore->chrhigh;
	segment_apply.plusp = true;
	segment_apply.genestrand = genestrand;
	segment_apply.leftp = true;

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

	segment_apply.segment_cache = &segment_cache;
	segment_apply.segment_chain_seen = segment_chain_seen;
	segment_apply.segment_chaini = &segment_chaini;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
	segment_apply.applied_chains = 0U;
	segment_apply.skipped_exact_chains = 0U;
	segment_apply.skipped_subtrees = 0U;
	segment_apply.total_solve_us = 0U;
	segment_apply.max_solve_us = 0U;
#endif

	Kmer_apply_eltarray_left(eltarray,/*elti*/auxinfo->elti,nexhaustive,max_nsplices,
				 listpool,segment_chain_skip_chain_p,se_segment_equiv_prev_prune_p,
				 se_segment_max_window,
				 /*materialized_order_p*/false,
				 segment_chain_apply_chain,(void *) &segment_apply);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
	fprintf(stderr,
		"[SEGMENT_STREAM_LEFT] elti=%d applied_chains=%llu skipped_exact=%llu skipped_subtrees=%llu "
		"total_solve_us=%llu max_solve_us=%llu total_chains=%d\n",
		auxinfo->elti,segment_apply.applied_chains,
		segment_apply.skipped_exact_chains,segment_apply.skipped_subtrees,
		segment_apply.total_solve_us,segment_apply.max_solve_us,
		segment_chaini);
#endif
      }
      i--;
      
    } else {
      /* minus */
      eltarray = this->exhaustive_eltarray_gminus;
      nexhaustive = this->nexhaustive_gminus;
      jj = order_gminus[j];
      auxinfo = auxinfo_gminus[jj];
      pathstore = pathstores_gminus[jj];

      if (se_segment_stream_chains_p == false) {
	chains = Kmer_solve_eltarray_right(eltarray,/*elti*/auxinfo->elti,nexhaustive,
					   max_nsplices,univdiagpool,listpool);
	for (p = chains; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Path_solve_from_right_chain(&(*found_score),&(*end_score),

				      &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
				      &(*sense_paths_gminus),&(*antisense_paths_gminus),

				      chain,queryseq,/*queryptr*/queryrc,querylength,mismatch_positions_alloc,
				      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				      localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				      this,knownsplicing,knownindels,
				      /*query_compress*/query_compress_rev,query_compress_fwd,query_compress_rev,
				      /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
				      /*chrhigh*/pathstore->chrhigh,/*plusp*/false,genestrand,
				      /*paired_end_p*/false,/*first_read_p*/true,
				      intlistpool,uintlistpool,univcoordlistpool,listpool,
				      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				      /*method*/SEGMENT1,/*find_splices_p*/true);
	}

      } else {
	struct Segment_chain_apply_T segment_apply;

	segment_apply.found_score = &(*found_score);
	segment_apply.end_score = &(*end_score);
	segment_apply.unextended_sense_paths = &(*unextended_sense_paths_gminus);
	segment_apply.unextended_antisense_paths = &(*unextended_antisense_paths_gminus);
	segment_apply.sense_paths = &(*sense_paths_gminus);
	segment_apply.antisense_paths = &(*antisense_paths_gminus);

	segment_apply.queryseq = queryseq;
	segment_apply.queryptr = queryrc;
	segment_apply.querylength = querylength;
	segment_apply.mismatch_positions_alloc = mismatch_positions_alloc;
	segment_apply.novel_univdiags_alloc = novel_univdiags_alloc;
	segment_apply.novel_diagonals_alloc = novel_diagonals_alloc;
	segment_apply.localdb_alloc = localdb_alloc;
	segment_apply.localdb_tally = localdb_tally;
	segment_apply.localdb_overflowp = localdb_overflowp;
	segment_apply.max_nblocks = max_nblocks;
	segment_apply.localdb_nmismatches_allowed = localdb_nmismatches_allowed;
	segment_apply.sufficient_score = sufficient_score;

	segment_apply.stage1 = this;
	segment_apply.knownsplicing = knownsplicing;
	segment_apply.knownindels = knownindels;
	segment_apply.query_compress_fwd = query_compress_fwd;
	segment_apply.query_compress_rev = query_compress_rev;
	segment_apply.query_compress = query_compress_rev;
	segment_apply.chrnum = pathstore->chrnum;
	segment_apply.chroffset = pathstore->chroffset;
	segment_apply.chrhigh = pathstore->chrhigh;
	segment_apply.plusp = false;
	segment_apply.genestrand = genestrand;
	segment_apply.leftp = false;

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

	segment_apply.segment_cache = &segment_cache;
	segment_apply.segment_chain_seen = segment_chain_seen;
	segment_apply.segment_chaini = &segment_chaini;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
	segment_apply.applied_chains = 0U;
	segment_apply.skipped_exact_chains = 0U;
	segment_apply.skipped_subtrees = 0U;
	segment_apply.total_solve_us = 0U;
	segment_apply.max_solve_us = 0U;
#endif

	Kmer_apply_eltarray_right(eltarray,/*elti*/auxinfo->elti,nexhaustive,max_nsplices,
				  listpool,segment_chain_skip_chain_p,se_segment_equiv_prev_prune_p,
				  se_segment_max_window,
				  /*materialized_order_p*/false,
				  segment_chain_apply_chain,(void *) &segment_apply);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
	fprintf(stderr,
		"[SEGMENT_STREAM_RIGHT] elti=%d applied_chains=%llu skipped_exact=%llu skipped_subtrees=%llu "
		"total_solve_us=%llu max_solve_us=%llu total_chains=%d\n",
		auxinfo->elti,segment_apply.applied_chains,
		segment_apply.skipped_exact_chains,segment_apply.skipped_subtrees,
		segment_apply.total_solve_us,segment_apply.max_solve_us,
		segment_chaini);
#endif
      }
      j--;
    }

    debug(printf("At segment search i %d, j %d, found score %d vs sufficient score %d, end_score %d\n",
		 i,j,*found_score,sufficient_score,*end_score));
  }
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[STAGE1_GEN_TIMING] segment_solve_us=%llu found_score=%d end_score=%d\n",
	  segment_stream_now_us() - stage_start_us,*found_score,*end_score);
  segment_chain_seen_report(segment_chain_seen);
#endif

  Path_segment_cache_free(&segment_cache,intlistpool,univcoordlistpool,
			  listpool,pathpool,transcriptpool,hitlistpool);
  segment_chain_seen_free(&segment_chain_seen);

  if (order_gminus != NULL) {
    FREE(order_gminus);
  }
  if (order_gplus != NULL) {
    FREE(order_gplus);
  }
  if (segment_counts_gminus != NULL) {
    FREE(segment_counts_gminus);
  }
  if (segment_counts_gplus != NULL) {
    FREE(segment_counts_gplus);
  }

  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool,listpool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool,listpool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);

#if 0
  /* 3. Prevalent (merging).  Equivalent of segment search */
  debug(printf("Single Read: 3.  Running Kmer prevalent\n"));
  assert(this->all_oligos_gen_filledp == true); /* From Extension_search */
  Stage1_fill_all_positions_gen(this,querylength,genestrand);
  Kmer_prevalent(&_univdiagonals_gplus,&auxinfo_gplus,&nunivdiagonals_gplus,
		 &_univdiagonals_gminus,&auxinfo_gminus,&nunivdiagonals_gminus,
		 this,querylength);
  debug(printf(" => %d plus and %d minus univdiagonals\n",nunivdiagonals_gplus,nunivdiagonals_gminus));
  debug(printf("Kmer_prevalent returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = KMER_PREVALENT;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

    Path_solve_from_univdiagonal(&(*found_score),&(*end_score),

				 &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
				 &(*sense_paths_gplus),&(*antisense_paths_gplus),

				 univdiagonal,auxinfo,queryseq,/*queryptr*/queryuc_ptr,
				 /*query_compress*/query_compress_fwd,
				 query_compress_fwd,query_compress_rev,
				 /*plusp*/true,querylength,mismatch_positions_alloc,
				   
				 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				 this,knownsplicing,knownindels,
				 /*paired_end_p*/false,/*first_read_p*/true,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				 /*method*/KMER_PREVALENT,/*find_splices_p*/true);
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

    Path_solve_from_univdiagonal(&(*found_score),&(*end_score),

				 &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
				 &(*sense_paths_gminus),&(*antisense_paths_gminus),

				 univdiagonal,auxinfo,queryseq,/*queryptr*/queryrc,
				 /*query_compress*/query_compress_rev,
				 query_compress_fwd,query_compress_rev,
				 /*plusp*/false,querylength,mismatch_positions_alloc,
				 
				 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				 localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

				 this,knownsplicing,knownindels,
				 /*paired_end_p*/false,/*first_read_p*/true,
				 
				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				 /*method*/KMER_PREVALENT,/*find_splices_p*/true);
  }

  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool,listpool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool,listpool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
#endif



  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    return true;

  } else if (gen_ext_p == true) {
    /* For single-end reads, perform extension, since no paired-end read
       to help with extension */
  final_extend:
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    stage_start_us = segment_stream_now_us();
#endif
    debug(printf("Single Read: 4.  Performing extension\n"));

    single_read_extend(&(*found_score),&(*end_score),this,

		       &(*sense_paths_gplus),&(*sense_paths_gminus),
		       &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		       &(*unextended_sense_paths_gplus),&(*unextended_sense_paths_gminus),
		       &(*unextended_antisense_paths_gplus),&(*unextended_antisense_paths_gminus),

		       queryseq,queryuc_ptr,queryrc,querylength,
		       query_compress_fwd,query_compress_rev,
		       knownsplicing,knownindels,mismatch_positions_alloc,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       localdb_tally,localdb_overflowp,max_nblocks,localdb_nmismatches_allowed,

		       intlistpool,uintlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
		       hitlistpool,spliceendsgen);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    fprintf(stderr,"[STAGE1_GEN_TIMING] final_extend_us=%llu found_score=%d end_score=%d\n",
	    segment_stream_now_us() - stage_start_us,*found_score,*end_score);
#endif
  } else {
    /* GENOME-EXTN controls final Path_extend cleanup of unextended paths. */
  }

  if (Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_FUSION)) {
    return *found_score <= sufficient_score ? true : false;
  } else {
    /* Current default behavior always processes complete paths and skips fusion. */
    return true;
  }
}


Path_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    T this, Shortread_T queryseq, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Localdb_T localdb,

		    Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		    unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		    int max_nblocks, int localdb_nmismatches_allowed,

		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool, Segmenteltpool_T segmenteltpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
		    Spliceendsgen_T spliceendsgen, bool single_cell_p, bool first_read_p,
		    Pass_T pass) {
  Path_T *patharray;
  List_T paths = NULL;

  bool any_imperfect_ends_p;

  List_T unextended_sense_paths_gplus = NULL, unextended_sense_paths_gminus = NULL,
    unextended_antisense_paths_gplus = NULL, unextended_antisense_paths_gminus = NULL;

  List_T sense_paths_gplus, sense_paths_gminus,
    antisense_paths_gplus, antisense_paths_gminus;
  List_T tr_sense_paths_gplus = NULL, tr_sense_paths_gminus = NULL,
    tr_antisense_paths_gplus = NULL, tr_antisense_paths_gminus = NULL;

  /* int *mismatch_positions_alloc; */

  int nmismatches_filter, mincoverage_filter, max_softclip_filter;
  int split_path_min_combined_coverage;
  int nmismatches_allowed;

  int sufficient_score, max_nsplices, querylength;
  char *queryuc_ptr, *queryrc;
  Compress_T query_compress_fwd, query_compress_rev;

  int found_score, end_score;
  Method_T last_method;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long stage_start_us;
#endif
  
#if 0
  bool first_read_p;
  if (single_cell_p == true) {
    first_read_p = false;
  } else {
    first_read_p = true;
  }
#endif

  if ((querylength = Shortread_fulllength(queryseq)) < index1part + index1interval - 1) {
    *npaths_primary = *npaths_altloc = 0;
    return (Path_T *) NULL;
  } else {
    queryuc_ptr = Shortread_queryuc_ptr(queryseq);
    queryrc = Shortread_queryrc(queryseq);
  }

  /* nmismatches_allowed means nmismatches_search and is not specified
     by the user.  The user-specified value for -m represents
     nmismatches_filter */
  /* TODO: make this dependent upon the defect_rate */
  nmismatches_allowed = querylength/20; /* was querylength/index1part */

  if (user_nmismatches_filter_float < 0.0) {
    /* Not specified, so don't filter */
    nmismatches_filter = querylength;
  } else if (user_nmismatches_filter_float < 1.0) {
    nmismatches_filter = (int) rint(user_nmismatches_filter_float * (double) querylength);
  } else {
    nmismatches_filter = (int) user_nmismatches_filter_float;
  }
  if (se_max_mismatches_filter_float >= 0.0) {
    int se_nmismatches_filter;

    if (se_max_mismatches_filter_float < 1.0) {
      se_nmismatches_filter = (int) rint(se_max_mismatches_filter_float * (double) querylength);
    } else {
      se_nmismatches_filter = (int) se_max_mismatches_filter_float;
    }
    if (se_nmismatches_filter < nmismatches_filter) {
      nmismatches_filter = se_nmismatches_filter;
    }
  }

  if (user_mincoverage_filter_float <= 0.0) {
    mincoverage_filter = 0;
  } else if (user_mincoverage_filter_float <= 1.0) {
    /* Assuming that --min-coverage=1 must mean 1.0 and not a coverage of 1 bp */
    mincoverage_filter = (int) rint(user_mincoverage_filter_float * (double) querylength);
  } else {
    mincoverage_filter = (int) user_mincoverage_filter_float;
  }
  if (se_max_softclip_filter_float < 0.0) {
    max_softclip_filter = -1;
  } else if (se_max_softclip_filter_float < 1.0) {
    max_softclip_filter = (int) rint(se_max_softclip_filter_float * (double) querylength);
  } else {
    max_softclip_filter = (int) se_max_softclip_filter_float;
  }

  if (split_path_min_combined_coverage_float <= 0.0) {
    split_path_min_combined_coverage = 0;
  } else if (split_path_min_combined_coverage_float <= 1.0) {
    split_path_min_combined_coverage =
      (int) rint(split_path_min_combined_coverage_float * (double) querylength);
  } else {
    split_path_min_combined_coverage = (int) split_path_min_combined_coverage_float;
  }

#if 0
  if (max_insertionlen > querylength) {
    max_insertionlen = querylength;
  }
#endif

  /* Now allocated in Stage1hr_new */
  /* mismatch_positions_alloc = (int *) MALLOC((querylength+MISMATCH_EXTRA)*sizeof(int)); */

  query_compress_fwd = Compress_new_fwd(queryuc_ptr,querylength);
  query_compress_rev = Compress_new_rev(queryuc_ptr,querylength);

  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    found_score = end_score = querylength;
    sufficient_score = querylength/20;
    max_nsplices = querylength/EXON_LENGTH + 1;
    last_method = METHOD_INIT;

    if (transcriptome_align_p == true &&
	single_read_tr_paths(&found_score,&end_score,sufficient_score,&last_method,

			     &sense_paths_gplus,&sense_paths_gminus,
			     &antisense_paths_gplus,&antisense_paths_gminus,
		       
			     this,/*genestrand*/0,queryseq,querylength,
			     knownsplicing,
			     query_compress_fwd,query_compress_rev,nmismatches_allowed,

			     this->mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
			     localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
			     localdb,localdb_nmismatches_allowed,

			     trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
			     listpool,trpathpool,pathpool,vectorpool,transcriptpool,
			     hitlistpool) == true &&
	no_converted_paths_p(sense_paths_gplus,sense_paths_gminus,
			     antisense_paths_gplus,antisense_paths_gminus) == false) {

      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else if (genome_align_p == false) {
      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else {
      if (transcriptome_align_p == true &&
	  no_converted_paths_p(sense_paths_gplus,sense_paths_gminus,
			       antisense_paths_gplus,antisense_paths_gminus)) {
	debug1(printf("Transcriptome search produced no converted paths; resetting scores for genomic fallback\n"));
	found_score = end_score = querylength;
	last_method = METHOD_INIT;
      }
      if (se_keep_trpath_candidates_p == true &&
	  no_converted_paths_p(sense_paths_gplus,sense_paths_gminus,
			       antisense_paths_gplus,antisense_paths_gminus) == false) {
	debug1(printf("Preserving converted transcriptome candidates across genomic fallback: sense_g+=%d sense_g-=%d anti_g+=%d anti_g-=%d\n",
		      List_length(sense_paths_gplus),List_length(sense_paths_gminus),
		      List_length(antisense_paths_gplus),List_length(antisense_paths_gminus)));
	tr_sense_paths_gplus = sense_paths_gplus;
	tr_sense_paths_gminus = sense_paths_gminus;
	tr_antisense_paths_gplus = antisense_paths_gplus;
	tr_antisense_paths_gminus = antisense_paths_gminus;
      }

      if (single_read_gen_paths(&found_score,&end_score,&last_method,
			       
				&unextended_sense_paths_gplus,&unextended_sense_paths_gminus,
				&unextended_antisense_paths_gplus,&unextended_antisense_paths_gminus,

				&sense_paths_gplus,&sense_paths_gminus,
				&antisense_paths_gplus,&antisense_paths_gminus,

				this,/*genestrand*/0,queryseq,queryuc_ptr,queryrc,querylength,
				knownsplicing,knownindels,query_compress_fwd,query_compress_rev,

				this->mismatch_positions_alloc,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				localdb_tally,localdb_overflowp,max_nblocks,
				nmismatches_allowed,repetitive_ef64,

				auxinfopool,pathstorepool,univdiagpool,segmenteltpool,
				intlistpool,uintlistpool,univcoordlistpool,
				listpool,pathpool,transcriptpool,vectorpool,
				hitlistpool,spliceendsgen) == true ||
	  Align_cascade_enabled_p(align_cascade,ALIGN_CASCADE_GEN_FUSION) == false) {

	if (se_keep_trpath_candidates_p == true) {
	  sense_paths_gplus = List_append(tr_sense_paths_gplus,sense_paths_gplus);
	  sense_paths_gminus = List_append(tr_sense_paths_gminus,sense_paths_gminus);
	  antisense_paths_gplus = List_append(tr_antisense_paths_gplus,antisense_paths_gplus);
	  antisense_paths_gminus = List_append(tr_antisense_paths_gminus,antisense_paths_gminus);
	  tr_sense_paths_gplus = tr_sense_paths_gminus = (List_T) NULL;
	  tr_antisense_paths_gplus = tr_antisense_paths_gminus = (List_T) NULL;
	}

	paths = List_append(sense_paths_gplus,
			    List_append(sense_paths_gminus,
					List_append(antisense_paths_gplus,
						    antisense_paths_gminus)));

      } else if (splicingp == false) {
	paths = List_append(sense_paths_gplus,
			    List_append(sense_paths_gminus,
					List_append(antisense_paths_gplus,
						    antisense_paths_gminus)));

      } else {
#ifdef FUSION
	/* Perform fusions, if single-end and splicingp is true */
	debug(printf("Single Read: 6.  Performing fusions\n"));

	debug(Stage1_list_extension(this));

	paths = single_read_fusion(&found_score,&end_score,this,querylength,

				   unextended_sense_paths_gplus,unextended_sense_paths_gminus,
				   unextended_antisense_paths_gplus,unextended_antisense_paths_gminus,

				   query_compress_fwd,query_compress_rev,
				   queryseq,knownsplicing,nmismatches_allowed,
				   univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
#endif
      }
    }
  }

  if (paths != NULL) {
    Path_gc(&unextended_sense_paths_gplus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_sense_paths_gminus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_antisense_paths_gplus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_antisense_paths_gminus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);

  } else {
    /* As last resort, use unextended paths */
     /* Should have called single_read_extend, which means no further extensions are possible */
    paths = List_append(unextended_sense_paths_gplus,
			List_append(unextended_sense_paths_gminus,
				    List_append(unextended_antisense_paths_gplus,
						unextended_antisense_paths_gminus)));
    /* unextended_sense_paths_gplus = (List_T) NULL; */
    /* unextended_sense_paths_gminus = (List_T) NULL; */
    /* unextended_antisense_paths_gplus = (List_T) NULL;*/
    /* unextended_antisense_paths_gminus = (List_T) NULL; */
  }
    
  if (paths == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    patharray = (Path_T *) NULL;

  } else {
    patharray = (Path_T *) List_to_array_out(paths,NULL);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    stage_start_us = segment_stream_now_us();
#endif
    patharray = Path_eval_and_sort(&(*npaths_primary),&(*npaths_altloc),
				   &(*first_absmq),&(*second_absmq),patharray,
				   /*npaths*/List_length(paths),
				   query_compress_fwd,query_compress_rev,queryuc_ptr,queryrc,
				   Shortread_quality_string(queryseq),nmismatches_filter,mincoverage_filter,
				   max_softclip_filter,
				   split_paths_p,split_path_max_overlap,split_path_min_segment,
				   split_path_min_combined_coverage,
				   intlistpool,univcoordlistpool,listpool,
				   pathpool,transcriptpool,hitlistpool,/*filterp*/true);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    fprintf(stderr,
	    "[STAGE1_SINGLE_TIMING] path_eval_us=%llu npaths_primary=%d npaths_altloc=%d\n",
	    segment_stream_now_us() - stage_start_us,*npaths_primary,*npaths_altloc);
#endif


    Hitlistpool_free_list(&paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }
  
  Compress_free(&query_compress_fwd);
  Compress_free(&query_compress_rev);
  /* FREE(mismatch_positions_alloc); */

  /* FREE(queryrc); -- Now taken from Shortread */

  return patharray;
}


void
Stage1hr_single_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       bool refine_partial_tr_p_in, bool tr_chain_search_p_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       bool split_paths_p_in, int split_path_max_overlap_in,
		       int split_path_min_segment_in, double split_path_min_combined_coverage_float_in,
		       bool splicingp_in, bool se_ext_dedup_exact_chains_p_in,
		       bool se_ext_equiv_prev_prune_p_in,
		       bool se_keep_trpath_candidates_p_in,
		       bool se_segment_stream_chains_p_in,
		       bool se_segment_dedup_exact_chains_p_in,
		       bool se_segment_prefix_skip_p_in,
		       bool se_segment_qstart_anchor_cache_p_in,
		       bool se_path_qend_local_cache_p_in,
		       bool se_emitted_path_dedup_p_in,
		       bool se_segment_equiv_prev_prune_p_in,
		       int se_segment_cache_threshold_in,
		       int se_segment_max_window_in,
		       double se_segment_min_kmer_avg_quality_in,
		       int se_segment_max_kmer_bad_bases_in,
		       int se_segment_bad_base_quality_in,
		       int segment_max_avg_kmer_positions_in,
		       double se_max_mismatches_filter_float_in,
		       double se_max_softclip_filter_float_in,
		       Align_cascade_T align_cascade_in) {

  mode = mode_in;
  index1part = index1part_in;
  index1interval = index1interval_in;
  index1part_tr = index1part_tr_in;

  transcriptome = transcriptome_in;
  genome_align_p = genome_align_p_in;
  transcriptome_align_p = transcriptome_align_p_in;
  align_cascade = align_cascade_in;
  refine_partial_tr_p = refine_partial_tr_p_in;
  tr_chain_search_p = tr_chain_search_p_in;

  user_nmismatches_filter_float = user_nmismatches_filter_float_in;
  user_mincoverage_filter_float = user_mincoverage_filter_float_in;
  split_paths_p = split_paths_p_in;
  split_path_max_overlap = split_path_max_overlap_in;
  split_path_min_segment = split_path_min_segment_in;
  split_path_min_combined_coverage_float = split_path_min_combined_coverage_float_in;

  splicingp = splicingp_in;
  se_ext_dedup_exact_chains_p = se_ext_dedup_exact_chains_p_in;
  se_ext_equiv_prev_prune_p = se_ext_equiv_prev_prune_p_in;
  se_keep_trpath_candidates_p = se_keep_trpath_candidates_p_in;
  se_segment_stream_chains_p = se_segment_stream_chains_p_in;
  se_segment_dedup_exact_chains_p = se_segment_dedup_exact_chains_p_in;
  se_segment_prefix_skip_p = se_segment_prefix_skip_p_in;
  se_segment_qstart_anchor_cache_p = se_segment_qstart_anchor_cache_p_in;
  se_path_qend_local_cache_p = se_path_qend_local_cache_p_in;
  se_emitted_path_dedup_p = se_emitted_path_dedup_p_in;
  se_segment_equiv_prev_prune_p = se_segment_equiv_prev_prune_p_in;
  se_segment_cache_threshold = se_segment_cache_threshold_in;
  se_segment_max_window = se_segment_max_window_in;
  se_segment_min_kmer_avg_quality = se_segment_min_kmer_avg_quality_in;
  se_segment_max_kmer_bad_bases = se_segment_max_kmer_bad_bases_in;
  se_segment_bad_base_quality = se_segment_bad_base_quality_in;
  segment_max_avg_kmer_positions = segment_max_avg_kmer_positions_in;
  se_max_mismatches_filter_float = se_max_mismatches_filter_float_in;
  se_max_softclip_filter_float = se_max_softclip_filter_float_in;

  return;
}
