#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "singleelt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
#include <sys/time.h>
#endif

#include "assert.h"
#include "univdiagdef.h"
#include "auxinfo.h"
#include "sedgesort.h"


#define SUFFICIENT_WIDTH 30

/* backtrack */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Singleelt_left_chains and Singleelt_right_chains */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
static unsigned long long
singleelt_now_us (void) {
  struct timeval tv;

  gettimeofday(&tv,(struct timezone *) NULL);
  return (unsigned long long) tv.tv_sec * 1000000ULL + (unsigned long long) tv.tv_usec;
}
#endif

#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
#define chainstats(x) x
#define CHAINPAIR_TABLE_SIZE 65536
#define CHAINPAIR_TABLE_MASK (CHAINPAIR_TABLE_SIZE - 1)
#define CHAINPAIR_TOP_N 12

struct Chainpair_count_T {
  bool occupiedp;
  int head_qstart;
  int head_qend;
  Univcoord_T head_univdiagonal;
  int next_qstart;
  int next_qend;
  Univcoord_T next_univdiagonal;
  unsigned long long count;
};

static unsigned long long
chainstats_mix (unsigned long long h, unsigned long long x) {
  h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
  return h;
}

static unsigned int
chainpair_hash_values (int head_qstart, int head_qend, Univcoord_T head_univdiagonal,
		       int next_qstart, int next_qend, Univcoord_T next_univdiagonal) {
  unsigned long long h = 1469598103934665603ULL;

  h = chainstats_mix(h,(unsigned long long) head_qstart);
  h = chainstats_mix(h,(unsigned long long) head_qend);
  h = chainstats_mix(h,(unsigned long long) head_univdiagonal);
  h = chainstats_mix(h,(unsigned long long) next_qstart);
  h = chainstats_mix(h,(unsigned long long) next_qend);
  h = chainstats_mix(h,(unsigned long long) next_univdiagonal);

  return (unsigned int) (h & CHAINPAIR_TABLE_MASK);
}

static bool
chainpair_matches_values_p (struct Chainpair_count_T *entry,
			    int head_qstart, int head_qend, Univcoord_T head_univdiagonal,
			    int next_qstart, int next_qend, Univcoord_T next_univdiagonal) {
  return entry->head_qstart == head_qstart &&
    entry->head_qend == head_qend &&
    entry->head_univdiagonal == head_univdiagonal &&
    entry->next_qstart == next_qstart &&
    entry->next_qend == next_qend &&
    entry->next_univdiagonal == next_univdiagonal;
}

static void
chainpair_add_values (struct Chainpair_count_T *table, int *nunique,
		      int head_qstart, int head_qend, Univcoord_T head_univdiagonal,
		      int next_qstart, int next_qend, Univcoord_T next_univdiagonal) {
  struct Chainpair_count_T *entry;
  unsigned int hash;
  int i;

  hash = chainpair_hash_values(head_qstart,head_qend,head_univdiagonal,
			       next_qstart,next_qend,next_univdiagonal);
  for (i = 0; i < CHAINPAIR_TABLE_SIZE; i++) {
    entry = &(table[(hash + i) & CHAINPAIR_TABLE_MASK]);
    if (entry->occupiedp == false) {
      entry->occupiedp = true;
      entry->head_qstart = head_qstart;
      entry->head_qend = head_qend;
      entry->head_univdiagonal = head_univdiagonal;
      entry->next_qstart = next_qstart;
      entry->next_qend = next_qend;
      entry->next_univdiagonal = next_univdiagonal;
      entry->count = 1U;
      *nunique += 1;
      return;
    } else if (chainpair_matches_values_p(entry,head_qstart,head_qend,head_univdiagonal,
					  next_qstart,next_qend,next_univdiagonal) == true) {
      entry->count++;
      return;
    }
  }

  return;
}

static void
chainpair_add (struct Chainpair_count_T *table, int *nunique, Univdiag_T head, Univdiag_T next) {
  chainpair_add_values(table,nunique,head->qstart,head->qend,head->univdiagonal,
		       next->qstart,next->qend,next->univdiagonal);
  return;
}

static unsigned long long
chainstats_record_chains (struct Chainpair_count_T *table, int *nunique,
			  unsigned long long *singleton_chains, List_T chains) {
  List_T p, chain;
  Univdiag_T head, next;
  unsigned long long nchains = 0U;

  for (p = chains; p != NULL; p = List_next(p)) {
    chain = (List_T) List_head(p);
    nchains++;
    head = (Univdiag_T) List_head(chain);
    if (List_next(chain) == NULL) {
      *singleton_chains += 1;
    } else {
      next = (Univdiag_T) List_second_value(chain);
      chainpair_add(table,&(*nunique),head,next);
    }
  }

  return nchains;
}

static void
chainstats_print_top_pairs (const char *tag, struct Chainpair_count_T *table) {
  struct Chainpair_count_T *entry, *top[CHAINPAIR_TOP_N];
  int i, j, k;

  for (i = 0; i < CHAINPAIR_TOP_N; i++) {
    top[i] = (struct Chainpair_count_T *) NULL;
  }

  for (i = 0; i < CHAINPAIR_TABLE_SIZE; i++) {
    entry = &(table[i]);
    if (entry->occupiedp == true) {
      for (j = 0; j < CHAINPAIR_TOP_N; j++) {
	if (top[j] == (struct Chainpair_count_T *) NULL || entry->count > top[j]->count) {
	  for (k = CHAINPAIR_TOP_N - 1; k > j; k--) {
	    top[k] = top[k - 1];
	  }
	  top[j] = entry;
	  break;
	}
      }
    }
  }

  for (i = 0; i < CHAINPAIR_TOP_N && top[i] != (struct Chainpair_count_T *) NULL; i++) {
    fprintf(stderr,
	    "%s top_pair rank=%d count=%llu "
	    "head=%llu:%d..%d next=%llu:%d..%d\n",
	    tag,i + 1,top[i]->count,
	    (unsigned long long) top[i]->head_univdiagonal,
	    top[i]->head_qstart,top[i]->head_qend,
	    (unsigned long long) top[i]->next_univdiagonal,
	    top[i]->next_qstart,top[i]->next_qend);
  }

  return;
}

static void
chainstats_summarize_prev_links (const char *tag, Singleelt_T *eltarray, int narray) {
  Singleelt_T elt, max_elt = (Singleelt_T) NULL;
  int i, nprev, multi_prev = 0, max_prev = 0;
  unsigned long long total_prev_links = 0U;

  for (i = 0; i < narray; i++) {
    elt = eltarray[i];
    nprev = List_length(elt->prev_elt_list);
    total_prev_links += (unsigned long long) nprev;
    if (nprev > 1) {
      multi_prev++;
    }
    if (nprev > max_prev) {
      max_prev = nprev;
      max_elt = elt;
    }
  }

  fprintf(stderr,
	  "%s dp narray=%d total_prev_links=%llu "
	  "multi_prev_elts=%d max_prev=%d",
	  tag,narray,total_prev_links,multi_prev,max_prev);
  if (max_elt != (Singleelt_T) NULL) {
    fprintf(stderr," max_prev_q=%d..%d max_prev_diag=%llu score=%d nsplices=%d",
	    max_elt->qstart,max_elt->qend,(unsigned long long) max_elt->univdiagonal,
	    max_elt->score,max_elt->nsplices);
  }
  fprintf(stderr,"\n");

  return;
}

static int
chainstats_find_elt_index (Singleelt_T *eltarray, int narray, Singleelt_T elt) {
  int i;

  for (i = 0; i < narray; i++) {
    if (eltarray[i] == elt) {
      return i;
    }
  }

  return -1;
}

static unsigned long long
chainstats_add_saturating (unsigned long long x, unsigned long long y) {
  unsigned long long max = ~0ULL;

  if (x > max - y) {
    return max;
  } else {
    return x + y;
  }
}

static unsigned long long
chainstats_count_chains_for_elt (Singleelt_T elt, Singleelt_T *eltarray, int narray,
				 unsigned long long *memo, bool *seen) {
  List_T p;
  Singleelt_T prev;
  unsigned long long total = 0U;
  int index;

  if ((index = chainstats_find_elt_index(eltarray,narray,elt)) < 0) {
    return 0U;
  } else if (seen[index] == true) {
    return memo[index];
  }

  seen[index] = true;
  if (elt->prev_elt_list == NULL) {
    memo[index] = 1U;
  } else {
    for (p = elt->prev_elt_list; p != NULL; p = List_next(p)) {
      prev = (Singleelt_T) List_head(p);
      total = chainstats_add_saturating(total,
		chainstats_count_chains_for_elt(prev,eltarray,narray,memo,seen));
    }
    memo[index] = total;
  }

  return memo[index];
}

static unsigned long long
chainstats_summarize_chain_counts (const char *tag, Singleelt_T *eltarray, int narray,
				   int max_nsplices, int *global_score_by_nsplices,
				   List_T *global_elt_list_by_nsplices) {
  Singleelt_T elt;
  List_T a;
  unsigned long long *memo, count, total_final_chains = 0U;
  bool *seen, *reported;
  int nsplices, nfinals, final_count_total = 0;
  int i, rank, best_i;

  memo = (unsigned long long *) CALLOC(narray,sizeof(unsigned long long));
  seen = (bool *) CALLOC(narray,sizeof(bool));

  for (i = 0; i < narray; i++) {
    chainstats_count_chains_for_elt(eltarray[i],eltarray,narray,memo,seen);
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    if ((nfinals = List_length(global_elt_list_by_nsplices[nsplices])) > 0) {
      final_count_total += nfinals;
      fprintf(stderr,
	      "%s global nsplices=%d score=%d finals=%d\n",
	      tag,nsplices,global_score_by_nsplices[nsplices],nfinals);
    }
    for (a = global_elt_list_by_nsplices[nsplices]; a != NULL; a = List_next(a)) {
      elt = (Singleelt_T) List_head(a);
      count = chainstats_count_chains_for_elt(elt,eltarray,narray,memo,seen);
      total_final_chains = chainstats_add_saturating(total_final_chains,count);
      fprintf(stderr,
	      "%s final_est nsplices=%d final=%d..%d diag=%llu score=%d prevs=%d chains=%llu\n",
	      tag,nsplices,elt->qstart,elt->qend,(unsigned long long) elt->univdiagonal,
	      elt->score,List_length(elt->prev_elt_list),count);
    }
  }

  fprintf(stderr,"%s predicted finals=%d total_chains=%llu\n",
	  tag,final_count_total,total_final_chains);

  if (total_final_chains >= 1000U || narray >= 64) {
    reported = (bool *) CALLOC(narray,sizeof(bool));
    for (rank = 1; rank <= CHAINPAIR_TOP_N; rank++) {
      best_i = -1;
      for (i = 0; i < narray; i++) {
	if (reported[i] == false && memo[i] > 1U &&
	    (best_i < 0 || memo[i] > memo[best_i])) {
	  best_i = i;
	}
      }
      if (best_i < 0) {
	break;
      }
      reported[best_i] = true;
      elt = eltarray[best_i];
      fprintf(stderr,
	      "%s top_elt rank=%d chains=%llu q=%d..%d diag=%llu score=%d nsplices=%d prevs=%d\n",
	      tag,rank,memo[best_i],elt->qstart,elt->qend,
	      (unsigned long long) elt->univdiagonal,elt->score,elt->nsplices,
	      List_length(elt->prev_elt_list));
    }
    FREE(reported);
  }

  FREE(seen);
  FREE(memo);

  return total_final_chains;
}

struct Backtrack_apply_stats_T {
  struct Chainpair_count_T *chainpair_table;
  int chainpair_unique;
  unsigned long long leaves;
  unsigned long long skipped_prefixes;
  unsigned long long singleton_chains;
  unsigned long long recursive_calls;
  int max_depth;
};

static void
chainstats_record_apply_leaf (struct Backtrack_apply_stats_T *stats, Singleelt_T *elt_stack,
			      int depth) {
  Singleelt_T head, next;

  stats->leaves++;
  if (depth > stats->max_depth) {
    stats->max_depth = depth;
  }

  if (depth == 0) {
    stats->singleton_chains++;
  } else {
    head = elt_stack[0];
    next = elt_stack[1];
    chainpair_add_values(stats->chainpair_table,&stats->chainpair_unique,
			 head->qstart,head->qend,head->univdiagonal,
			 next->qstart,next->qend,next->univdiagonal);
  }

  return;
}
#else
#define chainstats(x)
#endif



static int negative_gap_distance;
static int max_deletionlen;
static Univcoord_T positive_gap_distance;


#define T Singleelt_T

void
Singleelt_setup (int max_insertionlen, Chrpos_T shortsplicedist, int max_deletionlen_in) {

  negative_gap_distance = max_insertionlen;
  positive_gap_distance = (Univcoord_T) ((shortsplicedist > (Chrpos_T) max_deletionlen_in) ? shortsplicedist : (Chrpos_T) max_deletionlen_in);
  max_deletionlen = max_deletionlen_in;

  return;
}


T *
Singleelt_new_array (int n) {
  T *array;
  struct T *space;
  int i;

  array = (T *) MALLOC(n*sizeof(T));
  space = (struct T *) MALLOC(n*sizeof(struct T));

  array[0] = space;
  for (i = 1; i < n; i++) {
    array[i] = &(array[i-1][1]);
  }
#ifdef GSNAP_DIAGONAL_ROW_STATS
  for (i = 0; i < n; i++) {
    array[i]->diagstats_minfreq_min = 0U;
    array[i]->diagstats_minfreq_max = 0U;
    array[i]->diagstats_minfreq_sum = 0ULL;
    array[i]->diagstats_minfreq_n = 0;
    array[i]->diagstats_kmer_hits_min = 0U;
    array[i]->diagstats_kmer_hits_max = 0U;
    array[i]->diagstats_kmer_hits_sum = 0ULL;
    array[i]->diagstats_kmer_hits_n = 0;
  }
#endif

  return array;
}


void
Singleelt_gc_array (T *array, int n, Listpool_T listpool) {
  int i;
  struct T *space;
  T elt;

  if (array) {
    space = array[0];
    for (i = 0; i < n; i++) {
      elt = array[i];
      Listpool_free_list(&elt->prev_elt_list,listpool
			 listpool_trace(__FILE__,__LINE__));
    }

    FREE(space);
    FREE(array);
  }

  return;
}


void
Singleelt_print (T this) {
  printf("%d..%d with univdiagonal %u\n",
	 this->qstart,this->qend,this->univdiagonal);
  return;
}

static void
Singleelt_print_array (T *eltarray, int narray) {
  int i;
  T elt;

  for (i = 0; i < narray; i++) {
    printf("%d: ",i);
    elt = eltarray[i];
    Singleelt_print(elt);
    printf("\n");
  }

  printf("\n");

  return;
}


void
Singleelt_print_wlinks (T this) {
  List_T p;
  T prev;

  printf("%d..%d with univdiagonal %u\n",this->qstart,this->qend,this->univdiagonal);
  for (p = this->prev_elt_list; p != NULL; p = List_next(p)) {
    prev = (T) List_head(p);

    printf("    %d..%d %u score:%d on %d splices <- %u dist:%u score:%d on %d splices\n",
	   prev->qstart,prev->qend,prev->univdiagonal,prev->score,prev->nsplices,
	   this->univdiagonal,this->univdiagonal - prev->univdiagonal,
	   this->score,this->nsplices);
  }

  return;
}

#if 0
static void
Singleelt_print_array_wlinks (T *eltarray, int narray) {
  int i;
  T elt;

  for (i = 0; i < narray; i++) {
    printf(">%d: ",i);
    elt = eltarray[i];
    Singleelt_print_wlinks(elt);
    printf(" %u",elt->univdiagonal);
    printf("\n");
  }

  printf("\n");

  return;
}
#endif


/* Induces a partial order, where subsumed segments have no comparison */
int
Singleelt_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  int x_center, y_center;

  if (x->qstart < y->qstart && x->qend < y->qend) {
    return -1;
  } else if (y->qend > x->qend && y->qstart > x->qstart) {
    return +1;
  } else {
    /* Subsumption.  Order based on midpoints */
    x_center = (x->qstart + x->qend)/2;
    y_center = (y->qstart + y->qend)/2;
    if (x_center < y_center) {
      return -1;
    } else if (y_center < x_center) {
      return +1;
    } else {
      return 0;
    }
  }
}


static bool
significant_overlap_p (int *overlap, T p, T q) {
  int extent;
  int qstart_min, qstart_max, qend_min, qend_max;

  if (q->qend < p->qstart) {
    *overlap = 0;
    return false;

  } else if (p->qend < q->qstart) {
    *overlap = 0;
    return false;

  } else {
    qstart_min = q->qstart; qstart_max = p->qstart;
    if (p->qend < q->qend) {
      qend_min = p->qend; qend_max = q->qend;
    } else {
      qend_min = q->qend; qend_max = p->qend;
    }
	
    *overlap = qend_min - qstart_max;
    extent = qend_max - qstart_min;
	
    if (2*(*overlap) > extent) {
      /* Skip: Mostly overlap */
      return true;
    } else {
      return false;
    }
  }
}


static List_T
backtrack (T q, Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T chains = NULL, children, p;
  List_T chain;
  List_T prev_elt_list, a;
  T prev_elt;

  debug8(printf("Entered backtrack with %d..%d %u\n",
		q->qstart,q->qend,q->univdiagonal));
  prev_elt_list = q->prev_elt_list;

  if (prev_elt_list == NULL) {
    debug8(printf("prev_elt_list is NULL, so terminating recursion\n"));
    /* Starting elt, so terminate recursion and return a list with a single chain containing a single univdiag */
    chain = Univdiagpool_push(NULL,univdiagpool,
			      q->qstart,q->qend,/*nmismatches*/-1,
			      /*univdiagonal*/q->univdiagonal
			      univdiagpool_trace(__FILE__,__LINE__));
    return Listpool_push(NULL,listpool,(void *) chain
			 listpool_trace(__FILE__,__LINE__));

  } else {
    /* Recurse */
    debug8(printf("recursing on %d hits\n",List_length(prev_elt_list)));
    for (a = prev_elt_list; a != NULL; a = List_next(a)) {
      prev_elt = List_head(a);

      children = backtrack(prev_elt,univdiagpool,listpool);
      debug8(printf("hit returned %d children\n",List_length(children)));
      for (p = children; p != NULL; p = List_next(p)) {
	chain = (List_T) List_head(p);
	chain = Univdiagpool_push(chain,univdiagpool,
				  q->qstart,q->qend,/*nmismatches*/-1,
				  /*univdiagonal*/q->univdiagonal
				  univdiagpool_trace(__FILE__,__LINE__));
	chains = Listpool_push(chains,listpool,(void *) chain
			       listpool_trace(__FILE__,__LINE__));
      }
      Listpool_free_list(&children,listpool
			 listpool_trace(__FILE__,__LINE__));
    }

    return chains;
  }
}

static void
chain_stack_apply (Singleelt_chain_apply_T apply_chain, void *data,
		   T *elt_stack, int depth, struct Univdiag_T *univdiag_stack,
		   struct List_T *list_stack) {
  int i;
  T elt;

  for (i = 0; i <= depth; i++) {
    elt = elt_stack[i];
    univdiag_stack[i].univdiagonal = elt->univdiagonal;
    univdiag_stack[i].qstart = elt->qstart;
    univdiag_stack[i].qend = elt->qend;
    univdiag_stack[i].nmismatches = -1;

    list_stack[i].first = (void *) &(univdiag_stack[i]);
    list_stack[i].rest = (i < depth) ? &(list_stack[i + 1]) : (List_T) NULL;
  }

  apply_chain((List_T) &(list_stack[0]),data);
  return;
}

static void
backtrack_apply (T q, T *elt_stack, int depth,
		 Singleelt_chain_skip_T skip_chain_p,
		 Singleelt_chain_apply_T apply_chain, void *data,
		 struct Univdiag_T *univdiag_stack, struct List_T *list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
		 , struct Backtrack_apply_stats_T *chainstats
#endif
		 ) {
  List_T prev_elt_list, a;
  T prev_elt, head, next;

  elt_stack[depth] = q;

#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
    chainstats->recursive_calls++;
  }
#endif

  if (depth >= 1 && skip_chain_p != (Singleelt_chain_skip_T) NULL) {
    head = elt_stack[0];
    next = elt_stack[1];
    if (skip_chain_p(head->qstart,head->qend,head->univdiagonal,
		     next->qstart,next->qend,next->univdiagonal,data) == true) {
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
      if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
	chainstats->skipped_prefixes++;
      }
#endif
      return;
    }
  }

  prev_elt_list = q->prev_elt_list;
  if (prev_elt_list == NULL) {
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
    if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
      chainstats_record_apply_leaf(chainstats,elt_stack,depth);
    }
#endif
    chain_stack_apply(apply_chain,data,elt_stack,depth,univdiag_stack,list_stack);

  } else {
    for (a = prev_elt_list; a != NULL; a = List_next(a)) {
      prev_elt = List_head(a);
      backtrack_apply(prev_elt,elt_stack,depth + 1,
		      skip_chain_p,apply_chain,data,univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
		      ,chainstats
#endif
		      );
    }
  }

  return;
}

/* Emits the same chain order as backtrack(), without materializing chains.
   backtrack() prepends every returned child chain.  Therefore, at a node q,
   its output B(q) is the reverse of the concatenated B(p) outputs for its
   predecessor list.  Emitting B(q) requires visiting predecessors in reverse
   order while emitting reverse(B(p)); emitting reverse(B(q)) requires visiting
   predecessors in forward order while emitting B(p).  The traversal mode thus
   alternates at each depth.

   A predecessor list belongs to its Singleelt node.  We reverse it temporarily
   for the B(q) traversal and restore it before returning, which avoids both the
   per-chain pool allocation in backtrack() and an unbounded recursive walk over
   the predecessor list itself. */
static void
backtrack_apply_materialized_order (T q, T *elt_stack, int depth,
				    bool emit_B_p,
				    Singleelt_chain_skip_T skip_chain_p,
				    Singleelt_chain_apply_T apply_chain, void *data,
				    struct Univdiag_T *univdiag_stack, struct List_T *list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
				    , struct Backtrack_apply_stats_T *chainstats
#endif
				    ) {
  List_T prev_elt_list, a;
  T prev_elt, head, next;

  elt_stack[depth] = q;

#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
    chainstats->recursive_calls++;
  }
#endif

  if (depth >= 1 && skip_chain_p != (Singleelt_chain_skip_T) NULL) {
    head = elt_stack[0];
    next = elt_stack[1];
    if (skip_chain_p(head->qstart,head->qend,head->univdiagonal,
		     next->qstart,next->qend,next->univdiagonal,data) == true) {
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
      if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
	chainstats->skipped_prefixes++;
      }
#endif
      return;
    }
  }

  if (q->prev_elt_list == NULL) {
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
    if (chainstats != (struct Backtrack_apply_stats_T *) NULL) {
      chainstats_record_apply_leaf(chainstats,elt_stack,depth);
    }
#endif
    chain_stack_apply(apply_chain,data,elt_stack,depth,univdiag_stack,list_stack);

  } else if (emit_B_p == true) {
    /* B(q): reverse predecessor order, emitting reverse(B(p)). */
    prev_elt_list = List_reverse(q->prev_elt_list);
    q->prev_elt_list = prev_elt_list;
    for (a = prev_elt_list; a != NULL; a = List_next(a)) {
      prev_elt = List_head(a);
      backtrack_apply_materialized_order(prev_elt,elt_stack,depth + 1,
					 /*emit_B_p*/false,
					 skip_chain_p,apply_chain,data,
					 univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
					 ,chainstats
#endif
					 );
    }
    q->prev_elt_list = List_reverse(prev_elt_list);

  } else {
    /* reverse(B(q)): forward predecessor order, emitting B(p). */
    for (a = q->prev_elt_list; a != NULL; a = List_next(a)) {
      prev_elt = List_head(a);
      backtrack_apply_materialized_order(prev_elt,elt_stack,depth + 1,
					 /*emit_B_p*/true,
					 skip_chain_p,apply_chain,data,
					 univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
					 ,chainstats
#endif
					 );
    }
  }

  return;
}



/* Dynamic programming */
static int
compute_score (int *nsplices, Univcoord_T lower_univdiagonal, Univcoord_T higher_univdiagonal, int overlap,
	       int lower_qend, int higher_qstart, int prev_score, int prev_nsplices, int this_width) {
  int ninsertions;

  if (lower_univdiagonal == higher_univdiagonal) {
    /* Include all nucleotides in between */
    *nsplices = prev_nsplices;
    return prev_score + (higher_qstart - lower_qend) + this_width;

  } else if (lower_univdiagonal > higher_univdiagonal) {
    /* Penalize for ninsertions */
    *nsplices = prev_nsplices;
    ninsertions = (int) (lower_univdiagonal - higher_univdiagonal);
    return prev_score - overlap - ninsertions + this_width;

  } else if (higher_univdiagonal <= lower_univdiagonal + max_deletionlen) {
    *nsplices = prev_nsplices;
    return prev_score - overlap + this_width;

  } else if (higher_univdiagonal <= lower_univdiagonal + positive_gap_distance) {
    *nsplices = prev_nsplices + 1;
    return prev_score - overlap + this_width;

  } else {
    /* Distance is too long.  Can happen with segment search, which allows chains of positive_gap_distance */
    *nsplices = prev_nsplices;
    return 0;
  }
}
    

static bool
look_leftward (bool *extendp, bool *nosplicep,
	       T q, T p, int max_nsplices, Listpool_T listpool) {
  int score, overlap;
  int nsplices, prev_nsplices;

  /* printf("%d..%d looking back at %d..%d\n",q->qstart,q->qend,p->qstart,p->qend); */
  *extendp = *nosplicep = false;

  if (p->qstart >= q->qstart && p->qend <= q->qend) {
    /* q subsumes p */
    return false;

  } else if (p->qstart <= q->qstart && p->qend >= q->qend) {
    /* p subsumes q */
    return false;

  } else if (significant_overlap_p(&overlap,p,q) == true) {
    /* Skip */
    return true;

  } else {
    prev_nsplices = p->nsplices;
      
    debug9(printf("For %d..%d:%u (higher) looking leftward to %d..%d:%u (lower)",
		  q->qstart,q->qend,q->univdiagonal,p->qstart,p->qend,p->univdiagonal));
    if ((score = compute_score(&nsplices,/*lower_univdiagonal*/p->univdiagonal,
			       /*higher_univdiagonal*/q->univdiagonal,
			       overlap,/*lower_qend*/p->qend,/*higher_qstart*/q->qstart,
			       /*prev_score*/p->score,prev_nsplices,
			       /*this_width*/q->qend - q->qstart)) < q->score) {
      /* Skip: worse score */
      debug9(printf(" => worse score\n"));
      return false;

    } else if (nsplices > max_nsplices) {
      /* Skip: too many splices */
      debug9(printf(" => too many splices\n"));
      return false;
	
    } else if (score > q->score) {
      /* Reset: better score */
      debug9(printf(" => better score\n"));
      
      Listpool_free_list(&q->prev_elt_list,listpool
			 listpool_trace(__FILE__,__LINE__));

      q->prev_elt_list = Listpool_push(NULL,listpool,(void *) p
				       listpool_trace(__FILE__,__LINE__));
      q->nsplices = nsplices;
      q->score = score;
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;

    } else if (nsplices < q->nsplices) {
      /* Reset: equal score but fewer splices */
      debug9(printf(" => equal score but fewer splices\n"));
      
      Listpool_free_list(&q->prev_elt_list,listpool
			 listpool_trace(__FILE__,__LINE__));
      q->prev_elt_list = Listpool_push(NULL,listpool,(void *) p
				       listpool_trace(__FILE__,__LINE__));
      q->nsplices = nsplices;
      /* q->score = score; */
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;

    } else if (nsplices == q->nsplices) {
      /* Add: equal score and nsplices */
      debug9(printf(" => equal score and nsplices\n"));
      
      q->prev_elt_list = Listpool_push(q->prev_elt_list,listpool,(void *) p
				       listpool_trace(__FILE__,__LINE__));
      /* q->nsplices = nsplices; */
      /* q->score = score; */
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;
	
    } else {
      /* Skip: equal score and more splices */
      debug9(printf(" => equal score and more splices\n"));
      return false;
    }
  }
}


/* Uses simple elts, with only one univdiagonal each */
List_T
Singleelt_left_chains (T *eltarray, int narray, int max_nsplices,
		       Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T results;
  List_T chains;

  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  
  int elti, eltj;
  T p = (T) NULL, q, p_init;
  int nsplices;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;


  debug8(printf("Entering Singleelt_left_chains\n"));
  debug8(Singleelt_print_array_wlinks(eltarray,narray));
  debug9(printf("Entering Singleelt_left_chains\n"));
  debug9(Singleelt_print_array_wlinks(eltarray,narray));


  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));

  /* Quadruple loop: Looking at pairs of elts, p and q.  Within each
     pair, if they are non-overlapping, find links (using approximate
     intersection) between univdiagonals in p and those in q.  Since
     intersection is linear, complexity is O(narray * narray *
     ndiagonals).  Also, with sufficiency check, we could reduce
     complexity to O(narray * ndiagonals). */

  for (eltj = 0; eltj < narray; eltj++) { /*  */
    q = eltarray[eltj];

    /* Initialize links */
    q->prev_elt_list = (List_T) NULL;
    q->score = q->qend - q->qstart;
    q->nsplices = 0;


    /* Compute links */
    sufficientp = false;

    elti = eltj - 1;		/*  */
    while (elti >= 0 && sufficientp == false) { /*  */
      p = eltarray[elti];

      if (look_leftward(&extendp,&nosplicep,q,p,max_nsplices,listpool) == false) {
	/* Skip */

      } else if (extendp == false) {
	/* Not sufficient */

      } else if (p->qend - p->qstart >= SUFFICIENT_WIDTH) {
	sufficientp = true;

      } else if (nosplicep == true) {
	sufficientp = true;
      }

      elti--;			/*  */
    }

    /* Continue while p overlaps original p */
    p_init = p;
    while (elti >= 0 && significant_overlap_p(&overlap,eltarray[elti],p_init) == true) { /*  */
      p = eltarray[elti];
      look_leftward(&extendp,&nosplicep,q,p,max_nsplices,listpool);
      elti--;			/*  */
    }

    /* Revise global scores */
    nsplices = q->nsplices;
    if ((score = q->score) < global_score_by_nsplices[nsplices]) {
      /* Skip: worse score */
	  
    } else if (score > global_score_by_nsplices[nsplices]) {
      /* Reset: better score */
      Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			 listpool_trace(__FILE__,__LINE__));
      global_score_by_nsplices[nsplices] = score;
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(NULL,listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));

    } else {
      /* Add: equal score */
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));
    }
  }


  /* Compute chains by backtracking */
  results = (List_T) NULL;
  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices]; a != NULL; a = List_next(a)) {
      q = (T) List_head(a);
      chains = backtrack(q,univdiagpool,listpool);

#ifdef DEBUG9
      for (List_T c = chains; c != NULL; c = List_next(c)) {
	List_T chain = (List_T) List_head(c);
	Univdiag_chain_print_left(chain);
      }
#endif

      results = List_append(results,chains);
    }
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Listpool_free_list(&global_elt_list_by_nsplices[nsplices],listpool
		       listpool_trace(__FILE__,__LINE__));
  }

  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

  return results;
}

/* Organized by rightmost (highest) univdiagonal */
/* q (higher) looks back at p (lower) */
/* low univdiagonal (tail) <- high univdiagonal (head) */

static bool
singleelt_same_genomic_interval_p (T x, T y) {
  return x->univdiagonal + (Univcoord_T) x->qstart ==
    y->univdiagonal + (Univcoord_T) y->qstart &&
    x->univdiagonal + (Univcoord_T) x->qend ==
    y->univdiagonal + (Univcoord_T) y->qend;
}

static bool
singleelt_prev_list_has_same_genomic_interval_p (List_T prev_elt_list, T elt) {
  List_T p;

  for (p = prev_elt_list; p != NULL; p = List_next(p)) {
    if (singleelt_same_genomic_interval_p((T) List_head(p),elt) == true) {
      return true;
    }
  }

  return false;
}


static bool
look_rightward (bool *extendp, bool *nosplicep,
		T q, T p, int max_nsplices, Listpool_T listpool,
		bool dedup_equiv_prevp, unsigned long long *equiv_prev_skips) {
  int score, overlap;
  int nsplices, prev_nsplices;

  /* printf("%d..%d looking back at %d..%d\n",q->qstart,q->qend,p->qstart,p->qend); */
  *extendp = *nosplicep = false;

  if (p->qstart >= q->qstart && p->qend <= q->qend) {
    /* q subsumes p */
    return false;

  } else if (p->qstart <= q->qstart && p->qend >= q->qend) {
    /* p subsumes q */
    return false;

  } else if (significant_overlap_p(&overlap,p,q) == true) {
    /* Skip */
    return true;

  } else {
    prev_nsplices = p->nsplices;
      
    debug9(printf("For %d..%d:%u (lower) looking rightward to %d..%d:%u (higher)",
		  q->qstart,q->qend,q->univdiagonal,p->qstart,p->qend,p->univdiagonal));
    if ((score = compute_score(&nsplices,/*lower_univdiagonal*/q->univdiagonal,
			       /*higher_univdiagonal*/p->univdiagonal,
			       overlap,/*lower_qend*/q->qend,/*higher_qstart*/p->qstart,
			       /*prev_score*/p->score,prev_nsplices,
			       /*this_width*/q->qend - q->qstart)) < q->score) {
      /* Skip: worse score */
      debug9(printf(" => worse score\n"));
      return false;

    } else if (nsplices > max_nsplices) {
      /* Skip: too many splices */
      debug9(printf(" => too many splices\n"));
      return false;
	
    } else if (score > q->score) {
      /* Reset: better score */
      debug9(printf(" => better score\n"));
      
      Listpool_free_list(&q->prev_elt_list,listpool
			 listpool_trace(__FILE__,__LINE__));
	
      q->prev_elt_list = Listpool_push(NULL,listpool,(void *) p
				       listpool_trace(__FILE__,__LINE__));
      q->nsplices = nsplices;
      q->score = score;
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;

    } else if (nsplices < q->nsplices) {
      /* Reset: equal score but fewer splices */
      debug9(printf(" => equal score but fewer splices\n"));
      
      Listpool_free_list(&q->prev_elt_list,listpool
			 listpool_trace(__FILE__,__LINE__));
      q->prev_elt_list = Listpool_push(NULL,listpool,(void *) p
				       listpool_trace(__FILE__,__LINE__));
      q->nsplices = nsplices;
      /* q->score = score; */
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;

    } else if (nsplices == q->nsplices) {
      /* Add: equal score and nsplices */
      debug9(printf(" => equal score and nsplices\n"));

      if (dedup_equiv_prevp == true &&
	  singleelt_prev_list_has_same_genomic_interval_p(q->prev_elt_list,p) == true) {
	if (equiv_prev_skips != (unsigned long long *) NULL) {
	  *equiv_prev_skips += 1U;
	}
      } else {
	q->prev_elt_list = Listpool_push(q->prev_elt_list,listpool,(void *) p
					 listpool_trace(__FILE__,__LINE__));
      }
      /* q->nsplices = nsplices; */
      /* q->score = score; */
      *extendp = true;
      if (nsplices == prev_nsplices) {
	*nosplicep = true;
      }
      return true;
	
    } else {
      /* Skip: equal score and more splices */
      debug9(printf(" => equal score and more splices\n"));
      return false;
    }
  }
}


/* Uses simple elts, with only one univdiagonal each */
List_T
Singleelt_right_chains (T *eltarray, int narray, int max_nsplices,
			Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T results;
  List_T chains;
  
  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  
  int elti, eltj;
  T p = (T) NULL, q, p_init;
  int nsplices;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;

#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  struct Chainpair_count_T *chainpair_table;
  unsigned long long total_chains = 0U, singleton_chains = 0U, final_chains;
  int chainpair_unique = 0, final_count_total = 0, nfinals;
#endif


  debug8(printf("Entering Singleelt_right_chains\n"));
  debug8(Singleelt_print_array_wlinks(eltarray,narray));
  debug9(printf("Entering Singleelt_right_chains\n"));
  debug9(Singleelt_print_array_wlinks(eltarray,narray));


  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  chainpair_table = (struct Chainpair_count_T *)
    CALLOC(CHAINPAIR_TABLE_SIZE,sizeof(struct Chainpair_count_T));
#endif

  /* Quadruple loop: Looking at pairs of elts, p and q.  Within each
     pair, if they are non-overlapping, find links (using approximate
     intersection) between univdiagonals in p and those in q.  Since
     intersection is linear, complexity is O(narray * narray *
     ndiagonals).  Also, with sufficiency check, we could reduce
     complexity to O(narray * ndiagonals). */

  for (eltj = narray - 1; eltj >= 0; eltj--) { /*  */
    q = eltarray[eltj];

    /* Initialize links */
    q->prev_elt_list = (List_T) NULL;
    q->score = q->qend - q->qstart;
    q->nsplices = 0;


    /* Compute links */
    sufficientp = false;

    elti = eltj + 1;		/*  */
    while (elti < narray && sufficientp == false) { /*  */
      p = eltarray[elti];

      if (look_rightward(&extendp,&nosplicep,q,p,max_nsplices,listpool,
			 /*dedup_equiv_prevp*/false,
			 /*equiv_prev_skips*/(unsigned long long *) NULL) == false) {
	/* Skip */

      } else if (extendp == false) {
	/* Not sufficient */

      } else if (p->qend - p->qstart >= SUFFICIENT_WIDTH) {
	sufficientp = true;

      } else if (nosplicep == true) {
	sufficientp = true;
      }

      elti++;			/*  */
    }

    /* Continue while p overlaps original p */
    p_init = p;
    while (elti < narray && significant_overlap_p(&overlap,eltarray[elti],p_init) == true) { /*  */
      p = eltarray[elti];
      look_rightward(&extendp,&nosplicep,q,p,max_nsplices,listpool,
		     /*dedup_equiv_prevp*/false,
		     /*equiv_prev_skips*/(unsigned long long *) NULL);
      elti++;			/*  */
    }

    /* Revise global scores */
    nsplices = q->nsplices;
    if ((score = q->score) < global_score_by_nsplices[nsplices]) {
      /* Skip: worse score */
	  
    } else if (score > global_score_by_nsplices[nsplices]) {
      /* Reset: better score */
      Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			 listpool_trace(__FILE__,__LINE__));
      global_score_by_nsplices[nsplices] = score;
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(NULL,listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));

    } else {
      /* Add: equal score */
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));
    }
  }


#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  chainstats_summarize_prev_links("[CHAINSTAT_RIGHT]",eltarray,narray);
  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    if ((nfinals = List_length(global_elt_list_by_nsplices[nsplices])) > 0) {
      final_count_total += nfinals;
      fprintf(stderr,
	      "[CHAINSTAT_RIGHT] global nsplices=%d score=%d finals=%d\n",
	      nsplices,global_score_by_nsplices[nsplices],nfinals);
    }
  }
#endif

  /* Compute chains by backtracking, and add to univdiagonals and auxinfo lists */
  results = (List_T) NULL;
  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices]; a != NULL; a = List_next(a)) {
      q = (T) List_head(a);
      chains = backtrack(q,univdiagpool,listpool);
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
      final_chains = chainstats_record_chains(chainpair_table,&chainpair_unique,
					      &singleton_chains,chains);
      total_chains += final_chains;
      fprintf(stderr,
	      "[CHAINSTAT_RIGHT] final nsplices=%d final=%d..%d diag=%llu "
	      "score=%d prevs=%d chains=%llu\n",
	      nsplices,q->qstart,q->qend,(unsigned long long) q->univdiagonal,
	      q->score,List_length(q->prev_elt_list),final_chains);
#endif

#ifdef DEBUG9
      for (List_T c = chains; c != NULL; c = List_next(c)) {
	List_T chain = (List_T) List_head(c);
	Univdiag_chain_print_left(chain);
      }
#endif

      results = List_append(results,chains);
    }
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Listpool_free_list(&global_elt_list_by_nsplices[nsplices],listpool
		       listpool_trace(__FILE__,__LINE__));
  }

  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  fprintf(stderr,
	  "[CHAINSTAT_RIGHT] summary finals=%d total_chains=%llu "
	  "unique_head_next=%d singleton_chains=%llu\n",
	  final_count_total,total_chains,chainpair_unique,singleton_chains);
  chainstats_print_top_pairs("[CHAINSTAT_RIGHT]",chainpair_table);
  FREE(chainpair_table);
#endif

  return results;
}

void
Singleelt_left_chains_apply (T *eltarray, int narray, int max_nsplices,
			     Listpool_T listpool,
			     Singleelt_chain_skip_T skip_chain_p, bool equiv_prev_prune_p,
			     bool materialized_order_p,
			     Singleelt_chain_apply_T apply_chain, void *data) {
  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  T *elt_stack;
  struct Univdiag_T *univdiag_stack;
  struct List_T *list_stack;

  int elti, eltj;
  T p = (T) NULL, q, p_init;
  int nsplices;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long start_us, dp_us, apply_us, lookups = 0U, overlap_lookups = 0U;
#endif

  if (narray == 0) {
    return;
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  start_us = singleelt_now_us();
#endif

  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));

  for (eltj = 0; eltj < narray; eltj++) {
    q = eltarray[eltj];

    q->prev_elt_list = (List_T) NULL;
    q->score = q->qend - q->qstart;
    q->nsplices = 0;

    sufficientp = false;

    elti = eltj - 1;
    while (elti >= 0 && sufficientp == false) {
      p = eltarray[elti];
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      lookups++;
#endif

      if (look_leftward(&extendp,&nosplicep,q,p,max_nsplices,listpool) == false) {
	/* Skip */

      } else if (extendp == false) {
	/* Not sufficient */

      } else if (p->qend - p->qstart >= SUFFICIENT_WIDTH) {
	sufficientp = true;

      } else if (nosplicep == true) {
	sufficientp = true;
      }

      elti--;
    }

    p_init = p;
    while (elti >= 0 && significant_overlap_p(&overlap,eltarray[elti],p_init) == true) {
      p = eltarray[elti];
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      lookups++;
      overlap_lookups++;
#endif
      look_leftward(&extendp,&nosplicep,q,p,max_nsplices,listpool);
      elti--;
    }

    nsplices = q->nsplices;
    if ((score = q->score) < global_score_by_nsplices[nsplices]) {
      /* Skip: worse score */

    } else if (score > global_score_by_nsplices[nsplices]) {
      Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			 listpool_trace(__FILE__,__LINE__));
      global_score_by_nsplices[nsplices] = score;
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(NULL,listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));

    } else {
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  dp_us = singleelt_now_us() - start_us;
  start_us = singleelt_now_us();
#endif

  elt_stack = (T *) MALLOC(narray * sizeof(T));
  univdiag_stack = (struct Univdiag_T *) MALLOC(narray * sizeof(struct Univdiag_T));
  list_stack = (struct List_T *) MALLOC(narray * sizeof(struct List_T));

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices]; a != NULL; a = List_next(a)) {
      q = (T) List_head(a);
      if (materialized_order_p == true) {
	backtrack_apply_materialized_order(q,elt_stack,/*depth*/0,/*emit_B_p*/true,
					 skip_chain_p,apply_chain,data,
					 univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
					 ,(struct Backtrack_apply_stats_T *) NULL
#endif
					 );
      } else {
	backtrack_apply(q,elt_stack,/*depth*/0,skip_chain_p,apply_chain,data,
			univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
			,(struct Backtrack_apply_stats_T *) NULL
#endif
			);
      }
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  apply_us = singleelt_now_us() - start_us;
  fprintf(stderr,
	  "[SINGLEELT_APPLY_LEFT] narray=%d lookups=%llu overlap_lookups=%llu dp_us=%llu apply_us=%llu\n",
	  narray,lookups,overlap_lookups,dp_us,apply_us);
#endif

  FREE(list_stack);
  FREE(univdiag_stack);
  FREE(elt_stack);

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Listpool_free_list(&global_elt_list_by_nsplices[nsplices],listpool
		       listpool_trace(__FILE__,__LINE__));
  }

  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

  return;
}

void
Singleelt_right_chains_apply (T *eltarray, int narray, int max_nsplices,
			      Listpool_T listpool,
			      Singleelt_chain_skip_T skip_chain_p, bool equiv_prev_prune_p,
			      bool materialized_order_p,
			      Singleelt_chain_apply_T apply_chain, void *data) {
  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  T *elt_stack;
  struct Univdiag_T *univdiag_stack;
  struct List_T *list_stack;

  int elti, eltj;
  T p = (T) NULL, q, p_init;
  int nsplices;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long start_us, dp_us, apply_us, lookups = 0U, overlap_lookups = 0U;
  unsigned long long equiv_prev_skips = 0U;
#endif
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  struct Backtrack_apply_stats_T apply_stats;
  unsigned long long predicted_total_chains = 0U, leaves_before, skipped_before;
#endif

  if (narray == 0) {
    return;
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  start_us = singleelt_now_us();
#endif

  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  memset(&apply_stats,0,sizeof(apply_stats));
  apply_stats.chainpair_table = (struct Chainpair_count_T *)
    CALLOC(CHAINPAIR_TABLE_SIZE,sizeof(struct Chainpair_count_T));
#endif

  for (eltj = narray - 1; eltj >= 0; eltj--) {
    q = eltarray[eltj];

    q->prev_elt_list = (List_T) NULL;
    q->score = q->qend - q->qstart;
    q->nsplices = 0;

    sufficientp = false;

    elti = eltj + 1;
    while (elti < narray && sufficientp == false) {
      p = eltarray[elti];
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      lookups++;
#endif

      if (look_rightward(&extendp,&nosplicep,q,p,max_nsplices,listpool,
			 /*dedup_equiv_prevp*/equiv_prev_prune_p,
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
			 &equiv_prev_skips
#else
			 /*equiv_prev_skips*/(unsigned long long *) NULL
#endif
			 ) == false) {
	/* Skip */

      } else if (extendp == false) {
	/* Not sufficient */

      } else if (p->qend - p->qstart >= SUFFICIENT_WIDTH) {
	sufficientp = true;

      } else if (nosplicep == true) {
	sufficientp = true;
      }

      elti++;
    }

    p_init = p;
    while (elti < narray && significant_overlap_p(&overlap,eltarray[elti],p_init) == true) {
      p = eltarray[elti];
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
      lookups++;
      overlap_lookups++;
#endif
      look_rightward(&extendp,&nosplicep,q,p,max_nsplices,listpool,
		     /*dedup_equiv_prevp*/equiv_prev_prune_p,
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
		     &equiv_prev_skips
#else
		     /*equiv_prev_skips*/(unsigned long long *) NULL
#endif
		     );
      elti++;
    }

    nsplices = q->nsplices;
    if ((score = q->score) < global_score_by_nsplices[nsplices]) {
      /* Skip: worse score */

    } else if (score > global_score_by_nsplices[nsplices]) {
      Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			 listpool_trace(__FILE__,__LINE__));
      global_score_by_nsplices[nsplices] = score;
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(NULL,listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));

    } else {
      global_elt_list_by_nsplices[nsplices] =
	Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
		      listpool_trace(__FILE__,__LINE__));
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  dp_us = singleelt_now_us() - start_us;
  start_us = singleelt_now_us();
#endif
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  chainstats_summarize_prev_links("[CHAINSTAT_RIGHT_APPLY]",eltarray,narray);
  predicted_total_chains =
    chainstats_summarize_chain_counts("[CHAINSTAT_RIGHT_APPLY]",eltarray,narray,
				      max_nsplices,global_score_by_nsplices,
				      global_elt_list_by_nsplices);
#endif

  elt_stack = (T *) MALLOC(narray * sizeof(T));
  univdiag_stack = (struct Univdiag_T *) MALLOC(narray * sizeof(struct Univdiag_T));
  list_stack = (struct List_T *) MALLOC(narray * sizeof(struct List_T));

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices]; a != NULL; a = List_next(a)) {
      q = (T) List_head(a);
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
      leaves_before = apply_stats.leaves;
      skipped_before = apply_stats.skipped_prefixes;
#endif
      if (materialized_order_p == true) {
	backtrack_apply_materialized_order(q,elt_stack,/*depth*/0,/*emit_B_p*/true,
					 skip_chain_p,apply_chain,data,
					 univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
					 ,&apply_stats
#endif
					 );
      } else {
	backtrack_apply(q,elt_stack,/*depth*/0,skip_chain_p,apply_chain,data,
			univdiag_stack,list_stack
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
			,&apply_stats
#endif
			);
      }
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
      fprintf(stderr,
	      "[CHAINSTAT_RIGHT_APPLY] final_apply nsplices=%d final=%d..%d diag=%llu "
	      "score=%d prevs=%d applied=%llu skipped_prefixes=%llu\n",
	      nsplices,q->qstart,q->qend,(unsigned long long) q->univdiagonal,
	      q->score,List_length(q->prev_elt_list),
	      apply_stats.leaves - leaves_before,
	      apply_stats.skipped_prefixes - skipped_before);
#endif
    }
  }

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  apply_us = singleelt_now_us() - start_us;
  fprintf(stderr,
	  "[SINGLEELT_APPLY] narray=%d lookups=%llu overlap_lookups=%llu "
	  "equiv_prev_skips=%llu dp_us=%llu apply_us=%llu\n",
	  narray,lookups,overlap_lookups,equiv_prev_skips,dp_us,apply_us);
#endif
#ifdef GSNAP_DEBUG_SINGLEELT_CHAIN_STATS
  fprintf(stderr,
	  "[CHAINSTAT_RIGHT_APPLY] summary predicted_total_chains=%llu "
	  "applied_chains=%llu skipped_prefixes=%llu unique_head_next=%d "
	  "singleton_chains=%llu recursive_calls=%llu max_depth=%d\n",
	  predicted_total_chains,apply_stats.leaves,apply_stats.skipped_prefixes,
	  apply_stats.chainpair_unique,apply_stats.singleton_chains,
	  apply_stats.recursive_calls,apply_stats.max_depth);
  chainstats_print_top_pairs("[CHAINSTAT_RIGHT_APPLY]",apply_stats.chainpair_table);
  FREE(apply_stats.chainpair_table);
#endif

  FREE(list_stack);
  FREE(univdiag_stack);
  FREE(elt_stack);

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Listpool_free_list(&global_elt_list_by_nsplices[nsplices],listpool
		       listpool_trace(__FILE__,__LINE__));
  }

  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

  return;
}
