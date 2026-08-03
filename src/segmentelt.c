#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "segmentelt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset */

#include "assert.h"
#include "univdiagdef.h"
#include "sedgesort.h"


#ifdef LARGE_GENOMES
#include "intersect-approx-indices-uint8.h"
#else
#include "intersect-approx-indices-uint4.h"
#endif


#define SUFFICIENT_WIDTH SEGMENTELT_DEFAULT_TERMINAL_ANCHOR_MIN_SCORE

/* backtrack */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Segmentelt_smooth_left and Segmentelt_smooth_right */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Segmentelt_smooth, details */
#ifdef DEBUG9A
#define debug9a(x) x
#else
#define debug9a(x)
#endif

/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif


static int negative_gap_distance;
static int max_deletionlen;
static Univcoord_T positive_gap_distance;
static int terminal_anchor_min_score;
static int terminal_anchor_slop;


#define T Segmentelt_T

extern void
Segmentelt_setup (int max_insertionlen, Chrpos_T shortsplicedist, int max_deletionlen_in,
		  int terminal_anchor_min_score_in, int terminal_anchor_slop_in) {

  negative_gap_distance = max_insertionlen;
  positive_gap_distance = (Univcoord_T) ((shortsplicedist > (Chrpos_T) max_deletionlen_in) ? shortsplicedist : (Chrpos_T) max_deletionlen_in);
  max_deletionlen = max_deletionlen_in;
  terminal_anchor_min_score = terminal_anchor_min_score_in;
  terminal_anchor_slop = terminal_anchor_slop_in;

  return;
}

void
Segmentelt_free (T *old, Segmenteltpool_T segmenteltpool,
		 Intlistpool_T intlistpool, Listpool_T listpool) {
  int i;

  FREE((*old)->univdiagonals);

  for (i = 0; i < (*old)->ndiagonals; i++) {
    Listpool_free_list(&((*old)->prev_elt_lists[i]),listpool
		       listpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&((*old)->prev_diagonali_lists[i]),intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
  }

  FREE((*old)->nsplices);
  FREE((*old)->scores);
  FREE((*old)->prev_elt_lists);
  FREE((*old)->prev_diagonali_lists);

  Segmenteltpool_free_segmentelt(&(*old),segmenteltpool
				 segmenteltpool_trace(__FILE__,__LINE__));
}

void
Segmentelt_gc (List_T set, Segmenteltpool_T segmenteltpool,
	       Intlistpool_T intlistpool, Listpool_T listpool) {
  List_T p;
  T elt;

  for (p = set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    Segmentelt_free(&elt,segmenteltpool,intlistpool,listpool);
  }

  Listpool_free_list(&set,listpool
		     listpool_trace(__FILE__,__LINE__));
  return;
}


T
Segmentelt_new (int qstart, int qend, Univcoord_T *univdiagonals,
		int ndiagonals, Segmenteltpool_T segmenteltpool) {
  T new;

  new = Segmenteltpool_new_segmentelt(segmenteltpool
				      segmenteltpool_trace(__FILE__,__LINE__));

  new->qstart = qstart;
  new->qend = qend;
  new->univdiagonals = univdiagonals;
  new->ndiagonals = ndiagonals;
  
  new->prev_elt_lists = (List_T *) CALLOC(ndiagonals,sizeof(List_T));
  new->prev_diagonali_lists = (Intlist_T *) CALLOC(ndiagonals,sizeof(Intlist_T));
  new->scores = (int *) CALLOC(ndiagonals,sizeof(int));
  new->nsplices = (int *) CALLOC(ndiagonals,sizeof(int));

  return new;
}


void
Segmentelt_print (T this) {
  printf("%d..%d with %d univdiagonals\n",
	 this->qstart,this->qend,this->ndiagonals);
  return;
}

static void
Segmentelt_print_array (T *eltarray, int narray) {
  int i;
  T elt;

  for (i = 0; i < narray; i++) {
    printf("%d: ",i);
    elt = eltarray[i];
    Segmentelt_print(elt);
#if 0
    for (k = 0; k < elt->ndiagonals; k++) {
      printf(" %u",elt->univdiagonals[k]);
    }
#endif
    printf("\n");
  }

  printf("\n");

  return;
}


void
Segmentelt_print_wlinks (T this) {
  List_T p;
  Intlist_T q;
  T prev;
  int i, j;

  printf("%d..%d with %d univdiagonals\n",this->qstart,this->qend,this->ndiagonals);
  for (i = 0; i < this->ndiagonals; i++) {
    for (p = this->prev_elt_lists[i], q = this->prev_diagonali_lists[i]; p != NULL;
	 p = List_next(p), q = Intlist_next(q)) {
      prev = (T) List_head(p);
      j = Intlist_head(q);

      printf("    %d..%d #%d:%u score:%d on %d splices <- #%d:%u dist:%u score:%d on %d splices\n",
	     prev->qstart,prev->qend,j,prev->univdiagonals[j],prev->scores[j],prev->nsplices[j],
	     i,this->univdiagonals[i],
	     this->univdiagonals[i] - prev->univdiagonals[j],this->scores[i],this->nsplices[i]);
    }
  }

  return;
}

static void
Segmentelt_print_array_wlinks (T *eltarray, int narray) {
  int i, k;
  T elt;

  for (i = 0; i < narray; i++) {
    printf(">%d: ",i);
    elt = eltarray[i];
    Segmentelt_print_wlinks(elt);
    for (k = 0; k < elt->ndiagonals; k++) {
      printf(" %u",elt->univdiagonals[k]);
    }
    printf("\n");
  }

  printf("\n");

  return;
}


void
Segmentelt_dump (T this) {
  printf("Elt %p with bounds %d..%d and %d univdiags:\n",
	 this,this->qstart,this->qend,this->ndiagonals);
#if 0
  for (k = 0; k < this->ndiagonals; k++) {
    printf(" %u\n",this->univdiagonals[k]);
  }
#endif
  printf("\n");


  return;
}


/* Induces a partial order, where subsumed segments have no comparison */
int
Segmentelt_cmp (const void *a, const void *b) {
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
backtrack (T q, int diagonali, Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T chains = NULL, children, p;
  List_T chain;
  List_T prev_elt_list, a;
  Intlist_T prev_diagonali_list, b;
  T prev_elt;
  int prev_diagonali;

  debug8(printf("Entered backtrack with %d..%d, diagonali %d (%d total)\n",
		q->qstart,q->qend,diagonali,q->ndiagonals));
  prev_elt_list = q->prev_elt_lists[diagonali];
  prev_diagonali_list = q->prev_diagonali_lists[diagonali];

  if (prev_elt_list == NULL) {
    debug8(printf("prev_elt_list is NULL, so terminating recursion\n"));
    /* Starting elt, so terminate recursion and return a list with a single chain containing a single univdiag */
    chain = Univdiagpool_push(NULL,univdiagpool,
			      q->qstart,q->qend,/*nmismatches*/-1,
			      /*univdiagonal*/q->univdiagonals[diagonali]
			      univdiagpool_trace(__FILE__,__LINE__));
    return Listpool_push(NULL,listpool,(void *) chain
			 listpool_trace(__FILE__,__LINE__));

  } else {
    /* Recurse */
    debug8(printf("recursing on %d hits\n",List_length(prev_elt_list)));
    for (a = prev_elt_list, b = prev_diagonali_list; a != NULL; a = List_next(a), b = Intlist_next(b)) {
      prev_elt = List_head(a);
      prev_diagonali = Intlist_head(b);

      children = backtrack(prev_elt,prev_diagonali,univdiagpool,listpool);
      debug8(printf("hit returned %d children\n",List_length(children)));
      for (p = children; p != NULL; p = List_next(p)) {
	chain = (List_T) List_head(p);
	chain = Univdiagpool_push(chain,univdiagpool,
				  q->qstart,q->qend,/*nmismatches*/-1,
				  /*univdiagonal*/q->univdiagonals[diagonali]
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


static bool
left_terminal_anchor_p (T q, int score, int querylength) {
  return (score >= terminal_anchor_min_score && q->qend >= querylength - terminal_anchor_slop);
}

static bool
right_terminal_anchor_p (T q, int score) {
  return (score >= terminal_anchor_min_score && q->qstart <= terminal_anchor_slop);
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

  } else {
    assert(higher_univdiagonal <= lower_univdiagonal + positive_gap_distance);
    *nsplices = prev_nsplices + 1;
    return prev_score - overlap + this_width;
  }
}
    

/* Organized by rightmost (highest) univdiagonal */
/* q (higher) looks back at p (lower) */
/* low univdiagonal (tail) <- high univdiagonal (head) */
static bool
look_leftward (bool *extendp, bool *nosplicep,
	       T q, T p, int querylength, int max_nsplices,
	       bool prune_equal_score_p, unsigned long long *n_equal_prev_pruned,
	       Intlistpool_T intlistpool, Listpool_T listpool) {
  int *indices;
  int nindices;

  int p_diagonali, q_diagonali;
  int k, i;
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
    return false;

  } else {
#ifdef LARGE_GENOMES
    indices = Intersect_approx_indices_uint8(&nindices,
					     p->univdiagonals,p->ndiagonals,/*diagterm1*/-querylength,
					     q->univdiagonals,q->ndiagonals,/*diagterm2*/0,
					     /*below_slop*/negative_gap_distance,/*above_slop*/positive_gap_distance);
#else
    indices = Intersect_approx_indices_uint4(&nindices,
					     p->univdiagonals,p->ndiagonals,/*diagterm1*/-querylength,
					     q->univdiagonals,q->ndiagonals,/*diagterm2*/0,
					     /*below_slop*/negative_gap_distance,/*above_slop*/positive_gap_distance);
#endif
    for (k = i = 0; k < nindices; k++, i += 2) {
      p_diagonali = indices[i]; /* lower */
      q_diagonali = indices[i+1]; /* higher */
      prev_nsplices = p->nsplices[p_diagonali];
      
      debug9(printf("For %d..%d:%u (higher) looking leftward to %d..%d:%u (lower)",
		    q->qstart,q->qend,q->univdiagonals[q_diagonali],
		    p->qstart,p->qend,p->univdiagonals[p_diagonali]));
      if ((score = compute_score(&nsplices,/*lower_univdiagonal*/p->univdiagonals[p_diagonali],
				 /*higher_univdiagonal*/q->univdiagonals[q_diagonali],
				 overlap,/*lower_qend*/p->qend,/*higher_qstart*/q->qstart,
				 /*prev_score*/p->scores[p_diagonali],prev_nsplices,
				 /*this_width*/q->qend - q->qstart)) < q->scores[q_diagonali]) {
	/* Skip: worse score */
	debug9(printf(" => worse score\n"));

      } else if (nsplices > max_nsplices) {
	/* Skip: too many splices */
	debug9(printf(" => too many splices\n"));
	
      } else if (score > q->scores[q_diagonali]) {
	/* Reset: better score */
	debug9(printf(" => better score\n"));
	
	Listpool_free_list(&q->prev_elt_lists[q_diagonali],listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&q->prev_diagonali_lists[q_diagonali],intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	
	q->prev_elt_lists[q_diagonali] = Listpool_push(NULL,listpool,(void *) p
						       listpool_trace(__FILE__,__LINE__));
	q->prev_diagonali_lists[q_diagonali] =
	  Intlistpool_push(NULL,intlistpool,p_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	q->nsplices[q_diagonali] = nsplices;
	q->scores[q_diagonali] = score;
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}

      } else if (nsplices < q->nsplices[q_diagonali]) {
	/* Reset: equal score but fewer splices */
	debug9(printf(" => equal score but fewer splices\n"));
	
	Listpool_free_list(&q->prev_elt_lists[q_diagonali],listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&q->prev_diagonali_lists[q_diagonali],intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	
	q->prev_elt_lists[q_diagonali] = Listpool_push(NULL,listpool,(void *) p
						       listpool_trace(__FILE__,__LINE__));
	q->prev_diagonali_lists[q_diagonali] =
	  Intlistpool_push(NULL,intlistpool,p_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	q->nsplices[q_diagonali] = nsplices;
	/* q->scores[q_diagonali] = score; */
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}

      } else if (nsplices == q->nsplices[q_diagonali]) {
	/* Add: equal score and nsplices */
	debug9(printf(" => equal score and nsplices\n"));

	if (prune_equal_score_p == true && q->prev_elt_lists[q_diagonali] != NULL) {
	  *n_equal_prev_pruned += 1;
	} else {
	  q->prev_elt_lists[q_diagonali] = Listpool_push(q->prev_elt_lists[q_diagonali],listpool,(void *) p
							 listpool_trace(__FILE__,__LINE__));
	  q->prev_diagonali_lists[q_diagonali] =
	    Intlistpool_push(q->prev_diagonali_lists[q_diagonali],intlistpool,p_diagonali
			     intlistpool_trace(__FILE__,__LINE__));
	}
	/* q->nsplices[q_diagonali] = nsplices; */
	/* q->scores[q_diagonali] = score; */
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}
	
      } else {
	/* Skip: equal score and more splices */
	debug9(printf(" => equal score and more splices\n"));
      }
    }

    FREE(indices);
    return true;
  }
}


int **
Segmentelt_left_chains (int **nindices_by_nsplices, int max_nsplices,
			Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores, int *nunivdiagonals,
			T *eltarray, int narray, int querylength, bool retain_terminal_p,
			bool prune_equal_score_p,
			Univcoordlistpool_T univcoordlistpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Listpool_T listpool) {
  int **indices_by_nsplices;
  Univcoordlist_T univdiagonals_list;
  List_T auxinfo_list;
  Univcoord_T *unordered_univdiagonals, univdiagonal;
  Auxinfo_T *unordered_auxinfo, info_i, info_j;

  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  Intlist_T *global_diagonali_list_by_nsplices, b;
  List_T *terminal_elt_list_by_nsplices;
  Intlist_T *terminal_diagonali_list_by_nsplices;
  
#ifdef DEBUG9
  List_T chain, c;
#endif
  List_T chains = NULL;
  int elti, eltj, k, i, j, ii;
  int q_diagonali;
  T p, q, p_init;
  int nsplices, max_level, n;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;

  int *order;
  unsigned long long n_equal_prev_pruned = 0U;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long total_ndiagonals = 0U;
#endif


  debug8(printf("Entering Segmentelt_left_chains\n"));
  debug8(Segmentelt_print_array(eltarray,narray));
  debug9(printf("Entering Segmentelt_left_chains\n"));
  debug9(Segmentelt_print_array(eltarray,narray));

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  for (i = 0; i < narray; i++) {
    total_ndiagonals += (unsigned long long) eltarray[i]->ndiagonals;
  }
  fprintf(stderr,"[SEGMENTELT_LEFT_START] narray=%d total_ndiagonals=%llu prune_equal=%d\n",
	  narray,total_ndiagonals,prune_equal_score_p);
#endif


  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
  global_diagonali_list_by_nsplices = (Intlist_T *) CALLOC(max_nsplices + 1,sizeof(Intlist_T));
  terminal_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
  terminal_diagonali_list_by_nsplices = (Intlist_T *) CALLOC(max_nsplices + 1,sizeof(Intlist_T));

  /* Quadruple loop: Looking at pairs of elts, p and q.  Within each
     pair, if they are non-overlapping, find links (using approximate
     intersection) between univdiagonals in p and those in q.  Since
     intersection is linear, complexity is O(narray * narray *
     ndiagonals).  Also, with sufficiency check, we could reduce
     complexity to O(narray * ndiagonals). */

  for (eltj = 0; eltj < narray; eltj++) { /*  */
    q = eltarray[eltj];

    /* Initialize links */
    for (q_diagonali = 0; q_diagonali < q->ndiagonals; q_diagonali++) {
      q->prev_elt_lists[q_diagonali] = (List_T) NULL;
      q->prev_diagonali_lists[q_diagonali] = (Intlist_T) NULL;
      q->scores[q_diagonali] = q->qend - q->qstart;
      q->nsplices[q_diagonali] = 0;
    }

    /* Compute links */
    sufficientp = false;

    elti = eltj - 1;		/*  */
    while (elti >= 0 && sufficientp == false) { /*  */
      p = eltarray[elti];

      if (look_leftward(&extendp,&nosplicep,
			q,p,querylength,max_nsplices,
			prune_equal_score_p,&n_equal_prev_pruned,
			intlistpool,listpool) == false) {
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
      look_leftward(&extendp,&nosplicep,
		    q,p,querylength,max_nsplices,
		    prune_equal_score_p,&n_equal_prev_pruned,intlistpool,listpool);
      elti--;			/*  */
    }

    /* Revise global scores */
    for (q_diagonali = 0; q_diagonali < q->ndiagonals; q_diagonali++) {
      nsplices = q->nsplices[q_diagonali];

      if ((score = q->scores[q_diagonali]) < global_score_by_nsplices[nsplices]) {
	/* Skip: worse score */
	  
      } else if (score > global_score_by_nsplices[nsplices]) {
	/* Reset: better score */
	Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&(global_diagonali_list_by_nsplices[nsplices]),intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	global_score_by_nsplices[nsplices] = score;
	global_elt_list_by_nsplices[nsplices] =
	  Listpool_push(NULL,listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	global_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(NULL,intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	
      } else {
	/* Add: equal score */
	global_elt_list_by_nsplices[nsplices] =
	  Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	global_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(global_diagonali_list_by_nsplices[nsplices],intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
      }

      if (retain_terminal_p == true && left_terminal_anchor_p(q,score,querylength) == true) {
	terminal_elt_list_by_nsplices[nsplices] =
	  Listpool_push(terminal_elt_list_by_nsplices[nsplices],listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	terminal_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(terminal_diagonali_list_by_nsplices[nsplices],intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
      }
    }
  }


  /* Compute chains by backtracking, and store in univdiagonals and auxinfo */
  univdiagonals_list = (Univcoordlist_T) NULL;
  auxinfo_list = (List_T) NULL;

  *nindices_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices], b = global_diagonali_list_by_nsplices[nsplices];
	 a != NULL; a = List_next(a), b = Intlist_next(b)) {
      q = (T) List_head(a);
      q_diagonali = Intlist_head(b);
      chains = backtrack(q,q_diagonali,univdiagpool,listpool);

      /* Store chains in Auxinfo_T object */
      univdiagonals_list = Univcoordlistpool_push(univdiagonals_list,univcoordlistpool,
						  q->univdiagonals[q_diagonali]
						  univcoordlistpool_trace(__FILE__,__LINE__));
      info_i = Auxinfo_new_chains(/*method*/EXT,q->qstart,q->qend,
				  chains,/*chaindir*/CHAIN_LEFT,nsplices,max_nsplices,
				  auxinfopool);
      auxinfo_list = Listpool_push(auxinfo_list,listpool,(void *) info_i
				   listpool_trace(__FILE__,__LINE__));
      (*nindices_by_nsplices)[nsplices] += 1;
    }

    /* Paired concordance can promote terminal anchors that are not the
       best single-read chains.  Global-best terminal anchors are already
       present in the global list above. */
    for (a = terminal_elt_list_by_nsplices[nsplices], b = terminal_diagonali_list_by_nsplices[nsplices];
	 a != NULL; a = List_next(a), b = Intlist_next(b)) {
      q = (T) List_head(a);
      q_diagonali = Intlist_head(b);
      if (q->scores[q_diagonali] < global_score_by_nsplices[nsplices]) {
	chains = backtrack(q,q_diagonali,univdiagpool,listpool);

	univdiagonals_list = Univcoordlistpool_push(univdiagonals_list,univcoordlistpool,
						    q->univdiagonals[q_diagonali]
						    univcoordlistpool_trace(__FILE__,__LINE__));
	info_i = Auxinfo_new_chains(/*method*/EXT,q->qstart,q->qend,
				    chains,/*chaindir*/CHAIN_LEFT,nsplices,max_nsplices,
				    auxinfopool);
	auxinfo_list = Listpool_push(auxinfo_list,listpool,(void *) info_i
				     listpool_trace(__FILE__,__LINE__));
	(*nindices_by_nsplices)[nsplices] += 1;
      }
    }
  }


  /* Sort results by univdiagonal, eliminate duplicate univdiagonals, and store indices */
  indices_by_nsplices = (int **) CALLOC(max_nsplices + 1,sizeof(int *));
  if ((*nunivdiagonals = Univcoordlist_length(univdiagonals_list)) == 0) {

    *_univdiagonals = (Univcoord_T *) NULL;
    *auxinfo = (Auxinfo_T *) NULL;
    *pathstores = (Pathstore_T *) NULL;

  } else {
    /* Create indices */
    max_level = 0;
    for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
      if ((n = (*nindices_by_nsplices)[nsplices]) > max_level) {
	max_level = n;
      }
    }

    for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
      indices_by_nsplices[nsplices] = (int *) MALLOC(max_level * sizeof(int));
    }
    memset(*nindices_by_nsplices,0,nsplices * sizeof(int));

    /* Sort */
    *_univdiagonals = (Univcoord_T *) MALLOC_ALIGN((*nunivdiagonals)*sizeof(Univcoord_T));
    *auxinfo = (Auxinfo_T *) MALLOC((*nunivdiagonals)*sizeof(Auxinfo_T));
    *pathstores = (Pathstore_T *) MALLOC((*nunivdiagonals)*sizeof(Pathstore_T));

    unordered_univdiagonals = Univcoordlist_to_array(univdiagonals_list,/*end*/0);
    unordered_auxinfo = (Auxinfo_T *) List_to_array(auxinfo_list,NULL); /* Returns n+1 elts */

#ifdef LARGE_GENOMES
    order = Sedgesort_order_uint8(unordered_univdiagonals,*nunivdiagonals);
#else
    order = Sedgesort_order_uint4(unordered_univdiagonals,*nunivdiagonals);
#endif

    /* Eliminate duplicate univdiagonals and merge chains */
    k = 0;
    i = 0;
    while (i < *nunivdiagonals) {
      ii = order[i];
      univdiagonal = unordered_univdiagonals[ii];
      info_i = unordered_auxinfo[ii];

      j = i + 1;
      while (j < *nunivdiagonals && unordered_univdiagonals[order[j]] == univdiagonal) {
	info_j = unordered_auxinfo[order[j]];
	nsplices = info_j->nsplices;
	info_i->chains_by_nsplices[nsplices] = List_append(info_i->chains_by_nsplices[nsplices],
							   info_j->chains_by_nsplices[nsplices]);
	Auxinfo_free(&info_j,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);
	j++;
      }
      
#ifdef DEBUG9
      printf("%d..%d %u\n",info_i->qstart,info_i->qend,univdiagonal);
      for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
	for (c = info_i->chains_by_nsplices[nsplices]; c != NULL; c = List_next(c)) {
	  chain = (List_T) List_head(c);
	  Univdiag_chain_print_left(chain);
	}
      }
#endif

      (*_univdiagonals)[k] = univdiagonal;
      (*auxinfo)[k] = info_i;
      (*pathstores)[k] = Pathstore_new(pathstorepool);

      for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
	if (info_i->chains_by_nsplices[nsplices] != (List_T) NULL) {
	  n = (*nindices_by_nsplices)[nsplices]++;
	  debug9(printf("left nsplices %d, index %d => %d:%u\n",nsplices,n,k,univdiagonal));
	  indices_by_nsplices[nsplices][n] = k;
	}
      }

      k++;
      i = j;
    }

    *nunivdiagonals = k;
    
    FREE(order);
    FREE(unordered_auxinfo);
    FREE(unordered_univdiagonals);
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Intlistpool_free_list(&(terminal_diagonali_list_by_nsplices[nsplices]),intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(terminal_elt_list_by_nsplices[nsplices]),listpool
		       listpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&(global_diagonali_list_by_nsplices[nsplices]),intlistpool
			  intlistpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
		       listpool_trace(__FILE__,__LINE__));
  }
  FREE(terminal_diagonali_list_by_nsplices);
  FREE(terminal_elt_list_by_nsplices);
  FREE(global_diagonali_list_by_nsplices);
  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[SEGMENTELT_LEFT_DONE] narray=%d nunivdiagonals=%d equal_prev_pruned=%llu prune_equal=%d\n",
	  narray,*nunivdiagonals,n_equal_prev_pruned,prune_equal_score_p);
#endif

  return indices_by_nsplices;
}


/* Organized by leftmost (lowest) univdiagonal */
/* q (lower) looks back at p (higher) */
/* low univdiagonal (head) -> high univdiagonal (tail) */
static bool
look_rightward (bool *extendp, bool *nosplicep,
	       T q, T p, int querylength, int max_nsplices,
	       bool prune_equal_score_p, unsigned long long *n_equal_prev_pruned,
	       Intlistpool_T intlistpool, Listpool_T listpool) {
  int *indices;
  int nindices;

  int p_diagonali, q_diagonali;
  int k, i;
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
    return false;

  } else {
#ifdef LARGE_GENOMES
    indices = Intersect_approx_indices_uint8(&nindices,
					     q->univdiagonals,q->ndiagonals,/*diagterm1*/-querylength,
					     p->univdiagonals,p->ndiagonals,/*diagterm2*/0,
					     /*below_slop*/negative_gap_distance,/*above_slop*/positive_gap_distance);
#else
    indices = Intersect_approx_indices_uint4(&nindices,
					     q->univdiagonals,q->ndiagonals,/*diagterm1*/-querylength,
					     p->univdiagonals,p->ndiagonals,/*diagterm2*/0,
					     /*below_slop*/negative_gap_distance,/*above_slop*/positive_gap_distance);
#endif
    for (k = i = 0; k < nindices; k++, i += 2) {
      q_diagonali = indices[i]; /* lower */
      p_diagonali = indices[i+1]; /* higher */
      prev_nsplices = p->nsplices[p_diagonali];

      debug9(printf("For %d..%d:%u (lower) looking rightward to %d..%d:%u (higher)",
		    q->qstart,q->qend,q->univdiagonals[q_diagonali],
		    p->qstart,p->qend,p->univdiagonals[p_diagonali]));
      if ((score = compute_score(&nsplices,/*lower_univdiagonal*/q->univdiagonals[q_diagonali],
				 /*higher_univdiagonal*/p->univdiagonals[p_diagonali],
				 overlap,/*lower_qend*/q->qend,/*higher_qstart*/p->qstart,
				 /*prev_score*/p->scores[p_diagonali],prev_nsplices,
				 /*this_width*/q->qend - q->qstart)) < q->scores[q_diagonali]) {
	/* Skip: worse score */
	debug9(printf(" => worse score\n"));
	
      } else if (nsplices > max_nsplices) {
	/* Skip: too many splices */
	debug9(printf(" => too many splices\n"));
	
      } else if (score > q->scores[q_diagonali]) {
	/* Reset: better score */
	debug9(printf(" => better score\n"));
	
	Listpool_free_list(&(q->prev_elt_lists[q_diagonali]),listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&(q->prev_diagonali_lists[q_diagonali]),intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	
	q->prev_elt_lists[q_diagonali] = Listpool_push(NULL,listpool,(void *) p
						       listpool_trace(__FILE__,__LINE__));
	q->prev_diagonali_lists[q_diagonali] =
	  Intlistpool_push(NULL,intlistpool,p_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	q->nsplices[q_diagonali] = nsplices;
	q->scores[q_diagonali] = score;
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}
	
      } else if (nsplices < q->nsplices[q_diagonali]) {
	/* Reset: equal score but fewer splices */
	debug9(printf(" => equal score but fewer splices\n"));
	
	Listpool_free_list(&(q->prev_elt_lists[q_diagonali]),listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&(q->prev_diagonali_lists[q_diagonali]),intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	
	q->prev_elt_lists[q_diagonali] = Listpool_push(NULL,listpool,(void *) p
						       listpool_trace(__FILE__,__LINE__));
	q->prev_diagonali_lists[q_diagonali] =
	  Intlistpool_push(NULL,intlistpool,p_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	q->nsplices[q_diagonali] = nsplices;
	/* q->scores[q_diagonali] = score; */
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}
	
      } else if (nsplices == q->nsplices[q_diagonali]) {
	/* Add: equal score and nsplices */
	debug9(printf(" => equal score and nsplices\n"));

	if (prune_equal_score_p == true && q->prev_elt_lists[q_diagonali] != NULL) {
	  *n_equal_prev_pruned += 1;
	} else {
	  q->prev_elt_lists[q_diagonali] = Listpool_push(q->prev_elt_lists[q_diagonali],listpool,(void *) p
							 listpool_trace(__FILE__,__LINE__));
	  q->prev_diagonali_lists[q_diagonali] =
	    Intlistpool_push(q->prev_diagonali_lists[q_diagonali],intlistpool,p_diagonali
			     intlistpool_trace(__FILE__,__LINE__));
	}
	/* q->nsplices[q_diagonali] = nsplices; */
	/* q->scores[q_diagonali] = score; */
	*extendp = true;
	if (nsplices == prev_nsplices) {
	  *nosplicep = true;
	}
	
      } else {
	/* Skip: equal score and more splices */
	debug9(printf(" => equal score and more splices\n"));
      }
    }
    
    FREE(indices);
    return true;
  }
}


int **
Segmentelt_right_chains (int **nindices_by_nsplices, int max_nsplices,
			 Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores, int *nunivdiagonals,
			 T *eltarray, int narray, int querylength, bool retain_terminal_p,
			 bool prune_equal_score_p,
			 Univcoordlistpool_T univcoordlistpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			 Univdiagpool_T univdiagpool, Intlistpool_T intlistpool, Listpool_T listpool) {
  int **indices_by_nsplices;
  Univcoordlist_T univdiagonals_list;
  List_T auxinfo_list;
  Univcoord_T *unordered_univdiagonals, univdiagonal;
  Auxinfo_T *unordered_auxinfo, info_i, info_j;

  int *global_score_by_nsplices;
  List_T *global_elt_list_by_nsplices, a;
  Intlist_T *global_diagonali_list_by_nsplices, b;
  List_T *terminal_elt_list_by_nsplices;
  Intlist_T *terminal_diagonali_list_by_nsplices;
  
#ifdef DEBUG9
  List_T chain, c;
#endif
  List_T chains = NULL;
  int elti, eltj, k, i, j, ii;
  int q_diagonali;
  T p, q, p_init;
  int nsplices, max_level, n;

  int overlap, score;
  bool sufficientp, extendp, nosplicep;

  int *order;
  unsigned long long n_equal_prev_pruned = 0U;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long total_ndiagonals = 0U;
#endif


  debug8(printf("Entering Segmentelt_right_chains\n"));
  debug8(Segmentelt_print_array(eltarray,narray));
  debug9(printf("Entering Segmentelt_right_chains\n"));
  debug9(Segmentelt_print_array(eltarray,narray));

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  for (i = 0; i < narray; i++) {
    total_ndiagonals += (unsigned long long) eltarray[i]->ndiagonals;
  }
  fprintf(stderr,"[SEGMENTELT_RIGHT_START] narray=%d total_ndiagonals=%llu prune_equal=%d\n",
	  narray,total_ndiagonals,prune_equal_score_p);
#endif


  global_score_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  global_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
  global_diagonali_list_by_nsplices = (Intlist_T *) CALLOC(max_nsplices + 1,sizeof(Intlist_T));
  terminal_elt_list_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
  terminal_diagonali_list_by_nsplices = (Intlist_T *) CALLOC(max_nsplices + 1,sizeof(Intlist_T));

  /* Quadruple loop: Looking at pairs of elts, p and q.  Within each
     pair, if they are non-overlapping, find links (using approximate
     intersection) between univdiagonals in p and those in q.  Since
     intersection is linear, complexity is O(narray * narray *
     ndiagonals).  Also, with sufficiency check, we could reduce
     complexity to O(narray * ndiagonals). */

  for (eltj = narray - 1; eltj >= 0; eltj--) { /*  */
    q = eltarray[eltj];

    /* Initialize links */
    for (q_diagonali = 0; q_diagonali < q->ndiagonals; q_diagonali++) {
      q->prev_elt_lists[q_diagonali] = (List_T) NULL;
      q->prev_diagonali_lists[q_diagonali] = (Intlist_T) NULL;
      q->scores[q_diagonali] = q->qend - q->qstart;
      q->nsplices[q_diagonali] = 0;
    }

    /* Compute links */
    sufficientp = false;

    elti = eltj + 1;		/*  */
    while (elti < narray && sufficientp == false) { /*  */
      p = eltarray[elti];

      if (look_rightward(&extendp,&nosplicep,
			 q,p,querylength,max_nsplices,
			 prune_equal_score_p,&n_equal_prev_pruned,
			 intlistpool,listpool) == false) {
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
    while (elti < narray && significant_overlap_p(&overlap,eltarray[elti],p_init) == true) { /*  */      p = eltarray[elti];
      look_rightward(&extendp,&nosplicep,
		     q,p,querylength,max_nsplices,
		     prune_equal_score_p,&n_equal_prev_pruned,intlistpool,listpool);
      elti++;			/*  */
    }

    /* Revise global scores */
    for (q_diagonali = 0; q_diagonali < q->ndiagonals; q_diagonali++) {
      nsplices = q->nsplices[q_diagonali];

      if ((score = q->scores[q_diagonali]) < global_score_by_nsplices[nsplices]) {
	/* Skip: worse score */
	  
      } else if (score > global_score_by_nsplices[nsplices]) {
	/* Reset: better score */
	Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
			   listpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&(global_diagonali_list_by_nsplices[nsplices]),intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	global_score_by_nsplices[nsplices] = score;
	global_elt_list_by_nsplices[nsplices] =
	  Listpool_push(NULL,listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	global_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(NULL,intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
	
      } else {
	/* Add: equal score */
	global_elt_list_by_nsplices[nsplices] =
	  Listpool_push(global_elt_list_by_nsplices[nsplices],listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	global_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(global_diagonali_list_by_nsplices[nsplices],intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
      }

      if (retain_terminal_p == true && right_terminal_anchor_p(q,score) == true) {
	terminal_elt_list_by_nsplices[nsplices] =
	  Listpool_push(terminal_elt_list_by_nsplices[nsplices],listpool,(void *) q
			listpool_trace(__FILE__,__LINE__));
	terminal_diagonali_list_by_nsplices[nsplices] =
	  Intlistpool_push(terminal_diagonali_list_by_nsplices[nsplices],intlistpool,q_diagonali
			   intlistpool_trace(__FILE__,__LINE__));
      }
    }
  }


  /* Compute chains by backtracking, and store in univdiagonals and auxinfo */
  univdiagonals_list = (Univcoordlist_T) NULL;
  auxinfo_list = (List_T) NULL;

  *nindices_by_nsplices = (int *) CALLOC(max_nsplices + 1,sizeof(int));
  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    for (a = global_elt_list_by_nsplices[nsplices], b = global_diagonali_list_by_nsplices[nsplices];
	 a != NULL; a = List_next(a), b = Intlist_next(b)) {
      q = (T) List_head(a);
      q_diagonali = Intlist_head(b);
      chains = backtrack(q,q_diagonali,univdiagpool,listpool);

      /* Store chains in Auxinfo_T object */
      univdiagonals_list = Univcoordlistpool_push(univdiagonals_list,univcoordlistpool,
						  q->univdiagonals[q_diagonali]
						  univcoordlistpool_trace(__FILE__,__LINE__));
      info_i = Auxinfo_new_chains(/*method*/EXT,q->qstart,q->qend,
				  chains,/*chaindir*/CHAIN_RIGHT,nsplices,max_nsplices,
				  auxinfopool);
      auxinfo_list = Listpool_push(auxinfo_list,listpool,(void *) info_i
				   listpool_trace(__FILE__,__LINE__));
      (*nindices_by_nsplices)[nsplices] += 1;
    }

    /* Paired concordance can promote terminal anchors that are not the
       best single-read chains.  Global-best terminal anchors are already
       present in the global list above. */
    for (a = terminal_elt_list_by_nsplices[nsplices], b = terminal_diagonali_list_by_nsplices[nsplices];
	 a != NULL; a = List_next(a), b = Intlist_next(b)) {
      q = (T) List_head(a);
      q_diagonali = Intlist_head(b);
      if (q->scores[q_diagonali] < global_score_by_nsplices[nsplices]) {
	chains = backtrack(q,q_diagonali,univdiagpool,listpool);

	univdiagonals_list = Univcoordlistpool_push(univdiagonals_list,univcoordlistpool,
						    q->univdiagonals[q_diagonali]
						    univcoordlistpool_trace(__FILE__,__LINE__));
	info_i = Auxinfo_new_chains(/*method*/EXT,q->qstart,q->qend,
				    chains,/*chaindir*/CHAIN_RIGHT,nsplices,max_nsplices,
				    auxinfopool);
	auxinfo_list = Listpool_push(auxinfo_list,listpool,(void *) info_i
				     listpool_trace(__FILE__,__LINE__));
	(*nindices_by_nsplices)[nsplices] += 1;
      }
    }
  }


  /* Sort results by univdiagonal, eliminate duplicates, and store indices */
  indices_by_nsplices = (int **) CALLOC(max_nsplices + 1,sizeof(int *));
  if ((*nunivdiagonals = Univcoordlist_length(univdiagonals_list)) == 0) {
    *_univdiagonals = (Univcoord_T *) NULL;
    *auxinfo = (Auxinfo_T *) NULL;
    *pathstores = (Pathstore_T *) NULL;

  } else {
    /* Create indices */
    max_level = 0;
    for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
      if ((n = (*nindices_by_nsplices)[nsplices]) > max_level) {
	max_level = n;
      }
    }

    for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
      indices_by_nsplices[nsplices] = (int *) MALLOC(max_level * sizeof(int));
    }
    memset(*nindices_by_nsplices,0,nsplices * sizeof(int));

    /* Sort */
    *_univdiagonals = (Univcoord_T *) MALLOC_ALIGN((*nunivdiagonals)*sizeof(Univcoord_T));
    *auxinfo = (Auxinfo_T *) MALLOC((*nunivdiagonals)*sizeof(Auxinfo_T));
    *pathstores = (Pathstore_T *) MALLOC((*nunivdiagonals)*sizeof(Pathstore_T));

    unordered_univdiagonals = Univcoordlist_to_array(univdiagonals_list,/*end*/0);
    unordered_auxinfo = (Auxinfo_T *) List_to_array(auxinfo_list,NULL); /* Returns n+1 elts */

#ifdef LARGE_GENOMES
    order = Sedgesort_order_uint8(unordered_univdiagonals,*nunivdiagonals);
#else
    order = Sedgesort_order_uint4(unordered_univdiagonals,*nunivdiagonals);
#endif

    /* Eliminate duplicate univdiagonals and merge chains */
    k = 0;
    i = 0;
    while (i < *nunivdiagonals) {
      ii = order[i];
      univdiagonal = unordered_univdiagonals[ii];
      info_i = unordered_auxinfo[ii];

      j = i + 1;
      while (j < *nunivdiagonals && unordered_univdiagonals[order[j]] == univdiagonal) {
	info_j = unordered_auxinfo[order[j]];
	nsplices = info_j->nsplices;
	info_i->chains_by_nsplices[nsplices] = List_append(info_i->chains_by_nsplices[nsplices],
							   info_j->chains_by_nsplices[nsplices]);
	Auxinfo_free(&info_j,univdiagpool,auxinfopool,listpool,/*free_chains_p*/false);
	j++;
      }
      
#ifdef DEBUG9
      printf("%d..%d %u\n",info_i->qstart,info_i->qend,univdiagonal);
      for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
	for (c = info_i->chains_by_nsplices[nsplices]; c != NULL; c = List_next(c)) {
	  chain = (List_T) List_head(c);
	  Univdiag_chain_print_right(chain);
	}
      }
#endif

      (*_univdiagonals)[k] = univdiagonal;
      (*auxinfo)[k] = info_i;
      (*pathstores)[k] = Pathstore_new(pathstorepool);

      for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
	if (info_i->chains_by_nsplices[nsplices] != (List_T) NULL) {
	  n = (*nindices_by_nsplices)[nsplices]++;
	  debug9(printf("right nsplices %d, index %d => %d:%u\n",nsplices,n,k,univdiagonal));
	  indices_by_nsplices[nsplices][n] = k;
	}
      }

      k++;
      i = j;
    }

    *nunivdiagonals = k;
    
    FREE(order);
    FREE(unordered_auxinfo);
    FREE(unordered_univdiagonals);
  }

  for (nsplices = 0; nsplices <= max_nsplices; nsplices++) {
    Intlistpool_free_list(&(terminal_diagonali_list_by_nsplices[nsplices]),intlistpool
		     intlistpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(terminal_elt_list_by_nsplices[nsplices]),listpool
		       listpool_trace(__FILE__,__LINE__));
    Intlistpool_free_list(&(global_diagonali_list_by_nsplices[nsplices]),intlistpool
		     intlistpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(global_elt_list_by_nsplices[nsplices]),listpool
		       listpool_trace(__FILE__,__LINE__));
  }
  FREE(terminal_diagonali_list_by_nsplices);
  FREE(terminal_elt_list_by_nsplices);
  FREE(global_diagonali_list_by_nsplices);
  FREE(global_elt_list_by_nsplices);
  FREE(global_score_by_nsplices);

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[SEGMENTELT_RIGHT_DONE] narray=%d nunivdiagonals=%d equal_prev_pruned=%llu prune_equal=%d\n",
	  narray,*nunivdiagonals,n_equal_prev_pruned,prune_equal_score_p);
#endif

  return indices_by_nsplices;
}
