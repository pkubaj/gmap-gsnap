#ifndef SINGLEELT_INCLUDED
#define SINGLEELT_INCLUDED


typedef struct Singleelt_T *Singleelt_T;

#include "univcoord.h"
#include "list.h"


struct Singleelt_T {
  int qstart;
  int qend;
  Univcoord_T univdiagonal;

  /* For dynamic programming. Allows equi-optimal paths */
  List_T prev_elt_list;
  int score;
  int nsplices;

#ifdef GSNAP_DIAGONAL_ROW_STATS
  /* Raw observables for the exact query k-mers that generated this
     exhaustive diagonal.  Populated while the position streams are mapped
     back onto the unique diagonal array. */
  unsigned int diagstats_minfreq_min;
  unsigned int diagstats_minfreq_max;
  unsigned long long diagstats_minfreq_sum;
  int diagstats_minfreq_n;
  unsigned int diagstats_kmer_hits_min;
  unsigned int diagstats_kmer_hits_max;
  unsigned long long diagstats_kmer_hits_sum;
  int diagstats_kmer_hits_n;
#endif
};


#include "univdiagpool.h"
#include "auxinfopool.h"
#include "listpool.h"

#define T Singleelt_T

extern void
Singleelt_setup (int max_insertionlen, Chrpos_T shortsplicedist, int max_deletionlen_in);

extern T *
Singleelt_new_array (int n);

extern void
Singleelt_gc_array (T *array, int n, Listpool_T listpool);

extern void
Singleelt_print (T this);

extern void
Singleelt_dump (T this);

extern int
Singleelt_cmp (const void *a, const void *b);

extern List_T
Singleelt_left_chains (T *eltarray, int narray, int max_nsplices,
		       Univdiagpool_T univdiagpool, Listpool_T listpool);

extern List_T
Singleelt_right_chains (T *eltarray, int narray, int max_nsplices,
			Univdiagpool_T univdiagpool, Listpool_T listpool);

typedef bool (*Singleelt_chain_skip_T)(int head_qstart, int head_qend, Univcoord_T head_univdiagonal,
				       int next_qstart, int next_qend, Univcoord_T next_univdiagonal,
				       void *data);
typedef void (*Singleelt_chain_apply_T)(List_T chain, void *data);

extern void
Singleelt_left_chains_apply (T *eltarray, int narray, int max_nsplices,
			     Listpool_T listpool,
			     Singleelt_chain_skip_T skip_chain_p,
			     bool equiv_prev_prune_p,
			     bool materialized_order_p,
			     Singleelt_chain_apply_T apply_chain, void *data);

extern void
Singleelt_right_chains_apply (T *eltarray, int narray, int max_nsplices,
			      Listpool_T listpool,
			      Singleelt_chain_skip_T skip_chain_p,
			      bool equiv_prev_prune_p,
			      bool materialized_order_p,
			      Singleelt_chain_apply_T apply_chain, void *data);

#undef T
#endif
