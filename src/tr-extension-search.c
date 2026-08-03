#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "tr-extension-search.h"

#ifdef WORDS_BIGENDIAN
#define CONVERT(x) Bigendian_convert_uint(x)
#include "bigendian.h"
#else
#define CONVERT(x) (x)
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>		/* For munmap */

#include "mem.h"
#include "bool.h"
#include "assert.h"
#include "access.h"
#include "types.h"
#include "oligo.h"

#include "genomebits_consec.h"
#include "genomebits_trim.h"

#include "trdiagdef.h"
#include "sedgesort.h"

#include "trpath-solve.h"

#include "simd.h"

#if 0
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_AVX512
#include <immintrin.h>
#endif
#endif


/* Conservative uses Genomebits_consecutive procedures.  Aggressive uses Genomebits_trim procedures. */
/* #define CONSERVATIVE 1 */

#define TRIM_AT_TRANSCRIPTOME_BOUNDS 1
#define MAX_TOTAL_NPOSITIONS 10000


/* #define CHECK_OLIGOS 1 */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* binary_search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Sorting of diagonals */
#ifdef DEBUG12
#define debug12(x) x
#else
#define debug12(x)
#endif

/* filter_elts */
#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif


static Trcoord_T transcriptomelength;
static Transcriptome_T transcriptome;
static EF64_T transcript_ef64;

static Genomebits_T transcriptomebits;
static Indexdb_T tr_indexdb;

static int index1part_tr;
static int leftreadshift_tr;
static Oligospace_T oligobase_mask_tr;

static int max_insertionlen;
static int max_deletionlen;

/* Some limit is needed to prevent GSNAP from running very slowly */
/* Was MAX_HITS_FOR_BEST_ELT 1000 */

static int maxpaths_search;


void
Tr_extension_search_setup (Transcriptome_T transcriptome_in, Trcoord_T transcriptomelength_in, EF64_T transcript_ef64_in,
			   Genomebits_T transcriptomebits_in, Indexdb_T tr_indexdb_in,
			   int index1part_tr_in, int max_insertionlen_in, int max_deletionlen_in,
			   int maxpaths_search_in) {

  transcriptome = transcriptome_in;
  transcriptomelength = transcriptomelength_in;
  transcript_ef64 = transcript_ef64_in;

  tr_indexdb = tr_indexdb_in;

  transcriptomebits = transcriptomebits_in;

  index1part_tr = index1part_tr_in;

#ifdef HAVE_64_BIT
  leftreadshift_tr = 64 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#else
  leftreadshift_tr = 32 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#endif

  max_insertionlen = max_insertionlen_in;
  max_deletionlen = max_deletionlen_in;

  maxpaths_search = maxpaths_search_in;

  return;
}


/* Simplified version of Spanningelt_T */
#define T Tr_elt_T

static void
Tr_elt_free (T *old, Trdiagpool_T trdiagpool) {
  int i;

  for (i = 0; i < (*old)->n_all_trdiags; i++) {
    Trdiagpool_free_trdiag(&(*old)->all_trdiags[i],trdiagpool
			   trdiagpool_trace(__FILE__,__LINE__));
  }
  FREE((*old)->all_trdiags);
  FREE(*old);
  return;
}

void
Tr_elt_gc (List_T *set, Listpool_T listpool, Trdiagpool_T trdiagpool) {
  List_T p;
  T elt;

  for (p = *set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    Tr_elt_free(&elt,trdiagpool);
  }
  Listpool_free_list(&(*set),listpool
		     listpool_trace(__FILE__,__LINE__)); /* allocated by Listpool_push */
  return;
}


#if 0
static int
nt_querylength (char *query, int querylength) {
  int i;
  char c;

  i = 0;
  while (i < querylength && ((c = query[i]) == 'A' || c == 'C' || c == 'G' || c == 'T')) {
    i++;
  }

  return i;
}
#endif

#ifdef CHECK_OLIGOS
static Oligospace_T
nt_oligo (char *query, int indexsize) {
  Oligospace_T oligo = 0U;
  int i;

  for (i = 0; i < indexsize; i++) {
    oligo *= 4;
    
    printf("%c",query[i]);
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default:
      fprintf(stderr,"Saw N in nt_oligo\n");
      abort();
    }
  }
  printf("\n");

  return oligo;
}

#define LOW_TWO_BITS 0x3


static char *
oligo_nt (UINT4 oligo, int oligosize) {
  char *nt = MALLOC((oligosize+1)*sizeof(char));
  int i, j;
  UINT4 lowbits;

  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case RIGHT_A: nt[j] = 'A'; break;
    case RIGHT_C: nt[j] = 'C'; break;
    case RIGHT_G: nt[j] = 'G'; break;
    case RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }
  nt[oligosize] = '\0';

  return nt;
}

#endif


/* query is a substring of the original, starting with queryoffset */
/* Performs the same function as Sarray_lookup, which returns the diagonals */
static T
Elt_read_queryfwd (UINT4 *positions, int n, int diagterm, int querylength, int querystart,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Trdiagpool_T trdiagpool) {
  T new;
  int max_nmatches;
  Trdiag_T *out, *best_trdiags, *trdiags;
  Trcoord_T trdiagonal;

  int best_trimpos, trimpos, nmismatches;

  int *nmatches, pos3;
  int i;
  

  debug(printf("Got %d positions at querystart %d, plusp %d\n",n,querystart,plusp));
  if (n == 0) {
    return (T) NULL;
  } else if (querystart >= querylength - index1part_tr) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = querystart;

    nmatches = (int *) MALLOC(n*sizeof(int));
    trdiags = (Trdiag_T *) MALLOC(n*sizeof(Trdiag_T));
    for (i = 0; i < n; i++) {
      trdiagonal = positions[i] + diagterm;

#ifdef TRIM_AT_TRANSCRIPTOME_BOUNDS
      pos3 = (trdiagonal <= transcriptomelength) ? querylength : (int) (transcriptomelength - trdiagonal + querylength);
#else
      pos3 = querylength;
#endif

#ifdef CONSERVATIVE
      /* Use Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part_tr +
	Genomebits_consecutive_matches_rightward(transcriptomebits,query_compress,trdiagonal,querylength,
						 /*pos5*/querystart+index1part_tr,pos3,
						 plusp,genestrand);
      trdiags[i] = Trdiagpool_new_trdiag(trdiagpool,querystart,/*qend*/querystart+nmatches[i],
					 /*nmismatches*/0,trdiagonal
					 trdiagpool_trace(__FILE__,__LINE__));
      trimpos = querystart + nmatches[i];
      debug(printf("plusp %d, trdiagonal is %u, querystart is %d + nmatches %d => queryend is %d\n",
		   plusp,trdiagonal,querystart,nmatches[i],querystart + nmatches[i]));

#else
      /* Try Genomebits_trim_qend, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qend(&nmismatches,query_compress,
				     transcriptomebits,(Univcoord_T) trdiagonal,querylength,/*pos5*/querystart,pos3,
				     plusp,genestrand);
      debug(printf("qend conservative: At trdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		   trdiagonal,querystart,pos3,nmatches[i],querystart+nmatches[i]));
      debug(printf("qend aggressive: At trdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		   trdiagonal,querystart,pos3,trimpos,nmismatches,(trimpos - querystart) - nmismatches,trimpos));
      nmatches[i] = (trimpos - querystart) - nmismatches;
      trdiags[i] = Trdiagpool_new_trdiag(trdiagpool,querystart,/*qend*/trimpos,
					 nmismatches,trdiagonal
					 trdiagpool_trace(__FILE__,__LINE__));
#endif

      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
	best_trimpos = trimpos;
      }
    }

    out = best_trdiags = MALLOC(n*sizeof(Trdiag_T));
    for (i = 0; i < n; i++) {
      /* Accept diagonals within a certain threshold of the maximum */
      if (nmatches[i] >= max_nmatches /*- 1 (index1interval) + 1*/) {
	*out++ = trdiags[i];
      } else {
	Trdiagpool_free_trdiag(&(trdiags[i]),trdiagpool
			       trdiagpool_trace(__FILE__,__LINE__));
      }
    }
    FREE(trdiags);
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->min_qstart = querystart;
    /* new->max_qend = querystart + max_nmatches; */
    new->max_qend = best_trimpos;
    new->nmatches = max_nmatches;

    assert(new->max_qend <= (int) querylength);
    debug(printf("Making a queryfwd Elt %p with querystart %d, max_nmatches %d => max_queryend %d\n",
		 new,new->min_qstart,new->nmatches,new->max_qend));

    new->all_trdiags = new->trdiags = best_trdiags;
    new->n_all_trdiags = new->ntrdiags = out - best_trdiags;

#ifdef DEBUG
    printf("Trdiags:");
    for (i = 0; i < new->ntrdiags; i++) {
      printf("(%p) %u %d..%d",
	     new,new->trdiags[i]->trdiagonal,new->trdiags[i]->qstart,new->trdiags[i]->qend);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_trdiags;

    return new;
  }
}


static T
Elt_read_queryrev (UINT4 *positions, int n, int diagterm, int queryend, int querylength,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Trdiagpool_T trdiagpool) {
  T new;
  int max_nmatches;
  Trdiag_T *out, *best_trdiags, *trdiags;
  Trcoord_T trdiagonal;

  int best_trimpos, trimpos, nmismatches;

  int *nmatches, pos5;
  int i;
  

  debug(printf("Got %d positions at queryend %d, plusp %d\n",n,queryend,plusp));
  if (n == 0) {
    return (T) NULL;
  } else if (queryend < index1part_tr) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = queryend;

    nmatches = (int *) MALLOC(n*sizeof(int));
    trdiags = (Trdiag_T *) MALLOC(n*sizeof(Trdiag_T));
    for (i = 0; i < n; i++) {
      trdiagonal = positions[i] + diagterm;

#ifdef TRIM_AT_TRANSCRIPTOME_BOUNDS
      pos5 = (trdiagonal >= 0 + (Trcoord_T) querylength) ? 0 : (int) (querylength - trdiagonal);
#else
      pos5 = 0;
#endif
      
#ifdef CONSERVATIVE
      /* Use Genome_consecutive_matches_leftward, which is conservative */
      nmatches[i] = index1part_tr +
	Genomebits_consecutive_matches_leftward(transcriptomebits,query_compress,trdiagonal,querylength,
						pos5,/*pos3*/queryend - index1part_tr,plusp,genestrand);
      trdiags[i] = Trdiagpool_new_trdiag(trdiagpool,/*qstart*/queryend - nmatches[i],queryend,
					       /*nmismatches*/0,trdiagonal
					       trdiagpool_trace(__FILE__,__LINE__));
      trimpos = queryend - nmatches[i];
      debug(printf("plusp %d, trdiagonal is %u, queryend is %d - nmatches %d => querystart is %d\n",
		   plusp,trdiagonal,queryend,nmatches[i],queryend - nmatches[i]));
#else
      /* Try Genomebits_trim_qstart, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qstart(&nmismatches,query_compress,
				       transcriptomebits,(Univcoord_T) trdiagonal,querylength,pos5,/*pos3*/queryend,
				       plusp,genestrand);
      debug(printf("qstart conservative: At trdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		   trdiagonal,pos5,queryend,nmatches[i],queryend-nmatches[i]));
      debug(printf("qstart aggressive: At trdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		   trdiagonal,pos5,queryend,trimpos,nmismatches,(queryend - trimpos) - nmismatches,trimpos));
      nmatches[i] = (queryend - trimpos) - nmismatches;
      trdiags[i] = Trdiagpool_new_trdiag(trdiagpool,/*qstart*/trimpos,queryend,
					 nmismatches,trdiagonal
					 trdiagpool_trace(__FILE__,__LINE__));
#endif

      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
	best_trimpos = trimpos;
      }
    }

    out = best_trdiags = MALLOC(n*sizeof(Trdiag_T));
    for (i = 0; i < n; i++) {
      /* Accept diagonals within a certain threshold of the maximum */
      if (nmatches[i] >= max_nmatches /*- 1 (index1interval) + 1*/) {
	*out++ = trdiags[i];
      } else {
	Trdiagpool_free_trdiag(&(trdiags[i]),trdiagpool
			       trdiagpool_trace(__FILE__,__LINE__));
      }
    }
    FREE(trdiags);
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->max_qend = queryend;
    /* new->min_qstart = queryend - max_nmatches; */
    new->min_qstart = best_trimpos;
    new->nmatches = max_nmatches;

    assert(new->min_qstart >= 0);
    debug(printf("Making a queryrev Elt %p with queryend %d, max_nmatches %d => min_querystart %d\n",
		 new,new->max_qend,new->nmatches,new->min_qstart));

    new->all_trdiags = new->trdiags = best_trdiags;
    new->n_all_trdiags = new->ntrdiags = out - best_trdiags;

#ifdef DEBUG
    printf("Trdiags:");
    for (i = 0; i < new->ntrdiags; i++) {
      printf(" %u %d..%d",
	     new->trdiags[i]->trdiagonal,new->trdiags[i]->qstart,new->trdiags[i]->qend);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_trdiags;

    return new;
  }
}


static int
binary_search_trdiag (int lowi, int highi, Trdiag_T *trdiags, Trcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,trdiags[lowi]->trdiagonal,middlei,trdiags[middlei]->trdiagonal,
		   highi-1,trdiags[highi-1]->trdiagonal,goal));
    if (goal < trdiags[middlei]->trdiagonal) {
      highi = middlei;
    } else if (goal > trdiags[middlei]->trdiagonal) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


static void
Elt_filter_trdiags (T this, Trcoord_T low, Trcoord_T high) {
  int lowi, highi;
#ifdef DEBUG
  int i;
#endif

  debug(printf("Entered Elt_filter_trdiags on %d..%d with %d trdiags with low %u and high %u, nmatches %d\n",
	       this->min_qstart,this->max_qend,this->n_all_trdiags,low,high,this->nmatches));

  /* low_adj and high_adj are inclusive */
  lowi = binary_search_trdiag(/*lowi*/0,/*highi*/this->n_all_trdiags,this->all_trdiags,/*goal*/low);
  highi = binary_search_trdiag(lowi,/*highi*/this->n_all_trdiags,this->all_trdiags,/*goal*/high + 1) - 1;
  if ((this->ntrdiags = highi - lowi + 1) == 0) {
    this->trdiags = (Trdiag_T *) NULL;

  } else {
    this->trdiags = &(this->all_trdiags[lowi]);
  }

#ifdef DEBUG
  printf("Setting lowi %d and highi %d\n",lowi,highi);
  for (i = lowi; i <= highi; i++) {
    printf("  %u %d..%d\n",
	   this->all_trdiags[i]->trdiagonal,this->all_trdiags[i]->qstart,this->all_trdiags[i]->qend);
  }
#endif

  return;
}


#ifdef DEBUG
static void
Tr_elt_dump (T elt) {
  int k;

  printf("Elt %p with max querybounds %d..%d and %d trdiags:\n",
	 elt,elt->min_qstart,elt->max_qend,elt->ntrdiags);
  for (k = 0; k < elt->ntrdiags; k++) {
    printf("  (%p) %u %d..%d\n",
	   elt->trdiags[k],elt->trdiags[k]->trdiagonal,elt->trdiags[k]->qstart,elt->trdiags[k]->qend);
  }
  printf("\n");

  return;
}

static void
Tr_elt_dump_set (List_T set) {
  List_T p;

  for (p = set; p != NULL; p = List_next(p)) {
    Tr_elt_dump((T) List_head(p));
  }

  return;
}
#endif


#ifdef TRIM_AT_TRANSCRIPT_BOUNDS
/* lowbound was troffset; highbound was trhigh */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))
#else
#define add_bounded(x,plusterm) (x + (plusterm))
#define subtract_bounded(x,minusterm) ((x < (minusterm)) ? 0 : x - (minusterm))
#endif


static bool
elt_startp (T elt, int middle_qstart, int middle_qend) {
  if (elt->min_qstart >= middle_qstart && elt->max_qend <= middle_qend) {
    debug13(printf("Not allowing left elt that is subsumed by middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if (elt->max_qend >= middle_qend) {
    debug13(printf("Not allowing left elt that extends right of middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if ((elt->max_qend - middle_qstart) > (middle_qend - middle_qstart) / 2) {
    debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if ((elt->max_qend - middle_qstart) > (elt->max_qend - elt->min_qstart) / 2) {
    debug13(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else {
    return true;
  }
}


static bool
elt_endp (T elt, int middle_qstart, int middle_qend) {
  if (elt->min_qstart >= middle_qstart && elt->max_qend <= middle_qend) {
    debug13(printf("Not allowing right elt that is subsumed by middle elt: qpos %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if (elt->min_qstart <= middle_qstart) {
    debug13(printf("Not allowing right elt that extends left of middle elt: qpos %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if ((middle_qend - elt->min_qstart) > (middle_qend - middle_qstart) / 2) {
    debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if ((middle_qend - elt->min_qstart) > (elt->max_qend - elt->min_qstart) / 2) {
    debug13(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else {
    return true;
  }
}


/* Caller needs to call Tr_elt_gc(&plus_set) and Tr_elt_gc(&minus_set) */
static void
get_elt_sets_queryfwd (List_T *plus_set, List_T *minus_set, List_T *best_plus_elts, List_T *best_minus_elts,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int read_nmismatches_allowed, int genestrand, Listpool_T listpool,
		       Trdiagpool_T trdiagpool) {
  T elt;
  List_T p;
  int best_plus_nmatches, best_minus_nmatches;
  int plus_qpos, minus_qpos;

  int niter_plus, niter_minus;

  UINT4 *positions;
  int total_npositions = 0, npositions;

  int query_lastpos = querylength - index1part_tr;
  int queryoffset, querypos_rc;


  debug(printf("\nStarting get_elt_sets_queryfwd with querylength %d and nmismatches_allowed %d, genestrand %d\n",
	       querylength,read_nmismatches_allowed,genestrand));

#if 0
  /* Allow calls to Extension_search_queryfwd in addition to other methods */
  *paths_gplus = *paths_gminus = (List_T) NULL;
#endif

  if (read_nmismatches_allowed < 0) {
    read_nmismatches_allowed = 0;
#if 0
  } else {
    /* It is possible that this makes GSNAP too slow */
    read_nmismatches_allowed = querylength;
#endif
  }

  /* I.  Race from plus and minus start to end.  Compute best_plus_nmatches and best_plus_nmatches */
  best_plus_nmatches = best_minus_nmatches = 0;
  plus_qpos = minus_qpos = 0;
  niter_plus = niter_minus = 0;

  while (total_npositions <= MAX_TOTAL_NPOSITIONS &&
	 niter_plus <= read_nmismatches_allowed && niter_minus <= read_nmismatches_allowed &&
	 plus_qpos < query_lastpos && minus_qpos < query_lastpos) {
    if ((queryoffset = plus_qpos) >= query_lastpos) {
      /* Skip */
      plus_qpos += 1 /*index1interval*/;
    } else if (stage1->tr_validp[queryoffset] == false) {
      /* Skip */
      plus_qpos += 1 /*index1interval*/;
    } else {
      if (stage1->tr_plus_retrievedp[queryoffset] == true) {
	positions = stage1->tr_plus_positions[queryoffset];
	npositions = stage1->tr_plus_npositions[queryoffset];
      } else {
	assert(stage1->tr_plus_positions[queryoffset] == NULL);
	npositions = stage1->tr_plus_npositions[queryoffset] =
	  Indexdb_ptr(&stage1->tr_plus_positions[queryoffset],tr_indexdb,
		      stage1->tr_forward_oligos[queryoffset]);
	positions = stage1->tr_plus_positions[queryoffset];
	stage1->tr_plus_retrievedp[queryoffset] = true;
      }
      total_npositions += npositions;


      if (total_npositions > MAX_TOTAL_NPOSITIONS) {
	/* Skip */
	plus_qpos += 1 /*index1interval*/;
      } else if ((elt = Elt_read_queryfwd(positions,npositions,/*diagterm*/querylength - queryoffset,
					  querylength,/*querystart*/queryoffset,
					  query_compress_fwd,/*tplusp*/true,genestrand,
					  trdiagpool)) == NULL) {
	plus_qpos += 1 /*index1interval*/;
      } else {
	if (elt->nmatches > best_plus_nmatches && elt->n_all_trdiags <= maxpaths_search) {
	  if (elt->ntrdiags > 0) {
	    /* Could be 0 if there are too many trdiags */
	    best_plus_nmatches = elt->nmatches;
	  }
	}
	debug(printf("Pushing elt %p onto plus set\n",elt));
	*plus_set = Listpool_push(*plus_set,listpool,(void *) elt
				  listpool_trace(__FILE__,__LINE__));
	plus_qpos += elt->nmatches;
	niter_plus++;
      }
    }


    if ((queryoffset = minus_qpos) >= query_lastpos) {
      /* Skip */
      minus_qpos += 1 /*index1interval*/;
    } else if (stage1->tr_validp[(querypos_rc = query_lastpos - queryoffset)] == false) {
      /* Skip */
      minus_qpos += 1 /*index1interval*/;
    } else {
      if (stage1->tr_minus_retrievedp[querypos_rc] == true) {
	positions = stage1->tr_minus_positions[querypos_rc];
	npositions = stage1->tr_minus_npositions[querypos_rc];
      } else {
	assert(stage1->tr_minus_positions[querypos_rc] == NULL);
	npositions = stage1->tr_minus_npositions[querypos_rc] =
	  Indexdb_ptr(&stage1->tr_minus_positions[querypos_rc],tr_indexdb,
		      stage1->tr_revcomp_oligos[querypos_rc]);
	positions = stage1->tr_minus_positions[querypos_rc];
	stage1->tr_minus_retrievedp[querypos_rc] = true;
      }
      total_npositions += npositions;
	

      if (total_npositions > MAX_TOTAL_NPOSITIONS) {
	/* Skip */
	minus_qpos += 1 /*index1interval*/;
      } else if ((elt = Elt_read_queryfwd(positions,npositions,/*diagterm*/querylength - queryoffset,
					  querylength,/*querystart*/queryoffset,
					  query_compress_rev,/*tplusp*/false,genestrand,
					  trdiagpool)) == NULL) {
	minus_qpos += 1 /*index1interval*/;
      } else {
	if (elt->nmatches > best_minus_nmatches && elt->n_all_trdiags <= maxpaths_search) {
	  if (elt->ntrdiags > 0) {
	    /* Could be 0 if there are too many trdiags */
	    best_minus_nmatches = elt->nmatches;
	  }
	}
	debug(printf("Pushing elt %p onto minus set\n",elt));
	*minus_set = Listpool_push(*minus_set,listpool,(void *) elt
				   listpool_trace(__FILE__,__LINE__));
	minus_qpos += elt->nmatches;
	niter_minus++;
      }
    }
	
#ifdef DEBUG
    printf("\n");
    printf("plus_qpos %d\n",plus_qpos);
    printf("minus_qpos %d\n",minus_qpos);
    printf("\n");
#endif

#ifdef CONSERVATIVE
    /* Skip the presumed mismatch */
    plus_qpos += 1;
    minus_qpos += 1;
#endif
  }


  /* II.  Fill best_plus_elts and best_minus_elts */
  *best_plus_elts = *best_minus_elts = (List_T) NULL;
  for (p = *plus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_plus_nmatches - 3 && elt->n_all_trdiags <= maxpaths_search) {
      if (elt->ntrdiags > 0) {
	/* Could be 0 if there are too many trdiagonals */
	debug(printf("Pushing elt %p onto best_plus_elts\n",elt));
	*best_plus_elts = Listpool_push(*best_plus_elts,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (p = *minus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_minus_nmatches - 3 && elt->n_all_trdiags <= maxpaths_search) {
      if (elt->ntrdiags > 0) {
	/* Could be 0 if there are too many trdiags */
	debug(printf("Pushing elt %p onto best_minus_elts\n",elt));
	*best_minus_elts = Listpool_push(*best_minus_elts,listpool,(void *) elt
					 listpool_trace(__FILE__,__LINE__));
      }
    }
  }


  *plus_set = List_reverse(*plus_set);
  *minus_set = List_reverse(*minus_set);

#ifdef DEBUG
  printf("queryfwd plus set:\n");
  Tr_elt_dump_set(*plus_set);
  printf("\n");
  printf("queryfwd minus set:\n");
  Tr_elt_dump_set(*minus_set);
  printf("\n");
#endif

  if (minus_qpos >= query_lastpos && plus_qpos >= query_lastpos) {
    debug(printf("QUERYFWD: both sides won (minus_qpos %d, plus_qpos %d), so use both sides\n",
		 minus_qpos,plus_qpos));
    debug(printf("QUERYFWD PLUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYFWD MINUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_minus_elts));
    debug(printf("\n"));

  } else if (minus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD PLUS: minus side won (minus_qpos %d), so skip plus side\n",
		 minus_qpos));
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;

  } else if (plus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD MINUS: plus side won (plus_qpos %d), so skip minus side\n",
		 plus_qpos));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_minus_elts = (List_T) NULL;
    
  } else {
#if 0
    /* This code also has a high cost in speed */
    debug(printf("QUERYFWD: both sides lost (minus_qpos %d, plus_qpos %d), so use both sides\n",
		 minus_qpos,plus_qpos));
    debug(printf("QUERYFWD PLUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYFWD MINUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_minus_elts));
    debug(printf("\n"));
#else
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;
    *best_minus_elts = (List_T) NULL;
#endif
  }

  return;
}


/* Note: A second try, starting at queryrev, may solve only a few more
   cases presented to this method (i.e., cases that cannot be solved
   by kmer-end search).  It therefore may be best to push these cases
   to the next method. */

/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryrev (List_T *plus_set, List_T *minus_set, List_T *best_plus_elts, List_T *best_minus_elts,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int read_nmismatches_allowed, int genestrand, Listpool_T listpool,
		       Trdiagpool_T trdiagpool) {
  T elt;
  List_T p;
  int best_plus_nmatches, best_minus_nmatches;
  int plus_qpos, minus_qpos;

  int niter_plus, niter_minus;

  UINT4 *positions;
  int total_npositions = 0, npositions;

  int query_lastpos = querylength - index1part_tr;
  int queryoffset, querypos_rc;

  debug(printf("\nStarting get_elt_sets_queryrev with querylength %d and nmismatches_allowed %d, genestrand %d\n",
	       querylength,read_nmismatches_allowed,genestrand));

#if 0
  /* Allow calls to Extension_search_queryrev in addition to other methods */
  *paths_gplus = *paths_gminus = (List_T) NULL;
#endif

  if (read_nmismatches_allowed < 0) {
    read_nmismatches_allowed = 0;
#if 0
  } else {
    /* It is possible that this makes GSNAP too slow */
    read_nmismatches_allowed = querylength;
#endif
  }

  /* I.  Race from plus and minus end to start.  Compute best_plus_nmatches and best_minus_nmatches */
  best_plus_nmatches = best_minus_nmatches = 0;
  plus_qpos = minus_qpos = query_lastpos;
  niter_plus = niter_minus = 0;

  while (total_npositions <= MAX_TOTAL_NPOSITIONS &&
	 niter_plus <= read_nmismatches_allowed && niter_minus <= read_nmismatches_allowed &&
	 plus_qpos > 0 && minus_qpos > 0) {
    if ((queryoffset = plus_qpos) < 0) {
      /* Skip */
      plus_qpos -= 1 /*index1interval*/;
    } else if (stage1->tr_validp[queryoffset] == false) {
      /* Skip */
      plus_qpos -= 1 /*index1interval*/;
    } else {
      if (stage1->tr_plus_retrievedp[queryoffset] == true) {
	positions = stage1->tr_plus_positions[queryoffset];
	npositions = stage1->tr_plus_npositions[queryoffset];
      } else {
	assert(stage1->tr_plus_positions[queryoffset] == NULL);
	npositions = stage1->tr_plus_npositions[queryoffset] =
	  Indexdb_ptr(&stage1->tr_plus_positions[queryoffset],tr_indexdb,
		      stage1->tr_forward_oligos[queryoffset]);
	positions = stage1->tr_plus_positions[queryoffset];
	stage1->tr_plus_retrievedp[queryoffset] = true;
      }
      total_npositions += npositions;
	

      if (total_npositions > MAX_TOTAL_NPOSITIONS) {
	/* Skip */
	plus_qpos -= 1 /*index1interval*/;
      } else if ((elt = Elt_read_queryrev(positions,npositions,/*diagterm*/querylength - queryoffset,
					  /*queryend*/queryoffset + index1part_tr,querylength,
					  query_compress_fwd,/*tplusp*/true,genestrand,
					  trdiagpool)) == NULL) {
	plus_qpos -= 1 /*index1interval*/;
      } else {
	if (elt->nmatches > best_plus_nmatches && elt->n_all_trdiags <= maxpaths_search) {
	  if (elt->ntrdiags > 0) {
	    /* Could be 0 if there are too many trdiags */
	    best_plus_nmatches = elt->nmatches;
	  }
	}
	debug(printf("Pushing elt %p onto plus_set\n",elt));
	*plus_set = Listpool_push(*plus_set,listpool,(void *) elt
				  listpool_trace(__FILE__,__LINE__));
	plus_qpos -= elt->nmatches;
	niter_plus++;
      }
    }


    if ((queryoffset = minus_qpos) < 0) {
      /* Skip */
      minus_qpos -= 1 /*index1interval*/;
    } else if (stage1->tr_validp[(querypos_rc = query_lastpos - queryoffset)] == false) {
      /* Skip */
      minus_qpos -= 1 /*index1interval*/;
    } else {
      if (stage1->tr_minus_retrievedp[querypos_rc] == true) {
	positions = stage1->tr_minus_positions[querypos_rc];
	npositions = stage1->tr_minus_npositions[querypos_rc];
      } else {
	assert(stage1->tr_minus_positions[querypos_rc] == NULL);
	npositions = stage1->tr_minus_npositions[querypos_rc] =
	  Indexdb_ptr(&stage1->tr_minus_positions[querypos_rc],tr_indexdb,
		      stage1->tr_revcomp_oligos[querypos_rc]);
	positions = stage1->tr_minus_positions[querypos_rc];
	stage1->tr_minus_retrievedp[querypos_rc] = true;
      }
      total_npositions += npositions;
	

      if (total_npositions > MAX_TOTAL_NPOSITIONS) {
	/* Skip */
	minus_qpos -= 1 /*index1interval*/;
      } else if ((elt = Elt_read_queryrev(positions,npositions,/*diagterm*/querylength - queryoffset,
					  /*queryend*/queryoffset + index1part_tr,querylength,
					  query_compress_rev,/*tplusp*/false,genestrand,
					  trdiagpool)) == NULL) {
	minus_qpos -= 1 /*index1interval*/;
      } else {
	if (elt->nmatches > best_minus_nmatches && elt->n_all_trdiags <= maxpaths_search) {
	  if (elt->ntrdiags > 0) {
	    /* Could be 0 if there are too many trdiags */
	    best_minus_nmatches = elt->nmatches;
	  }
	}
	debug(printf("Pushing elt %p onto minus_set\n",elt));
	*minus_set = Listpool_push(*minus_set,listpool,(void *) elt
				   listpool_trace(__FILE__,__LINE__));
	minus_qpos -= elt->nmatches;
	niter_minus++;
      }
    }

#ifdef DEBUG
    printf("\n");
    printf("plus_qpos %d\n",plus_qpos);
    printf("minus_qpos %d\n",minus_qpos);
    printf("\n");
#endif

#ifdef CONSERVATIVE
    /* Skip the presumed mismatch */
    plus_qpos -= 1;
    minus_qpos -= 1;
#endif
  }


  /* II.  Fill best_plus_elts and best_minus_elts */
  *best_plus_elts = *best_minus_elts = (List_T) NULL;
  for (p = *plus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_plus_nmatches - 3 && elt->n_all_trdiags <= maxpaths_search) {
      if (elt->ntrdiags > 0) {
	/* Could be 0 if there are too many trdiags */
	debug(printf("Pushing elt %p onto best_plus_elts\n",elt));
	*best_plus_elts = Listpool_push(*best_plus_elts,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (p = *minus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_minus_nmatches - 3 && elt->n_all_trdiags <= maxpaths_search) {
      if (elt->ntrdiags > 0) {
	/* Could be 0 if there are too many trdiags */
	debug(printf("Pushing elt %p onto best_minus_elts\n",elt));
	*best_minus_elts = Listpool_push(*best_minus_elts,listpool,(void *) elt
					 listpool_trace(__FILE__,__LINE__));
      }
    }
  }
  

#if 0
  /* Not needed for queryrev */
  *plus_set = List_reverse(*plus_set);
  *minus_set = List_reverse(*minus_set);
#endif

#ifdef DEBUG
  printf("queryrev plus set:\n");
  Tr_elt_dump_set(*plus_set);
  printf("\n");
  printf("queryrev minus set:\n");
  Tr_elt_dump_set(*minus_set);
  printf("\n");
#endif

  if (minus_qpos <= 0 && plus_qpos <= 0) {
    debug(printf("QUERYREV: both sides won (minus_qpos %d, plus_qpos %d), so use both sides\n",
		 minus_qpos,plus_qpos));
    debug(printf("QUERYREV PLUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYREV MINUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_minus_elts));
    debug(printf("\n"));

  } else if (minus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV PLUS: minus side won (minus_qpos %d), so skip plus side\n",
		 minus_qpos));
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;

  } else if (plus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV MINUS: plus side won (plus_qpos %d), so skip minus side\n",
		 plus_qpos));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_minus_elts = (List_T) NULL;

  } else {
#if 0
    /* This code also has a high cost in speed */
    debug(printf("QUERYREV: both sides lost (minus_qpos %d, plus_qpos %d), so use both sides\n",
		 minus_qpos,plus_qpos));
    debug(printf("QUERYREV PLUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYREV MINUS HAS BEST ELTS:\n"));
    debug(Tr_elt_dump_set(*best_minus_elts));
    debug(printf("\n"));
#else
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;
    *best_minus_elts = (List_T) NULL;
#endif
  }

  return;
}


static void
process_seed (int *found_score, int *end_score, int sufficient_score,
	      List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
	      List_T *sense_trpaths, List_T *antisense_trpaths,

	      Trcoord_T middle_trdiagonal, int middle_qstart, int middle_qend, int middle_nmismatches,
	      List_T queryfwd_set, List_T queryrev_set, T queryfwd_best_elt, T queryrev_best_elt,
	      Stage1_T stage1, int querylength, int *mismatch_positions_alloc, 

	      Compress_T query_compress_fwd, Compress_T query_compress_rev, bool tplusp,
	      bool paired_end_p, bool first_read_p,
	      Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, 
	      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
	      Hitlistpool_T hitlistpool, Method_T method) {

  Trpath_T trpath;
  bool partialp;

  List_T left_trdiags = NULL, right_trdiags = NULL, p;
  Trdiag_T *left_trdiag_array, *right_trdiag_array, qstart_trdiag, qend_trdiag;
  Trcoord_T trdiagonal, start_low, start_high, end_low, end_high;
  int qstart, qend;
  Tr_elt_T elt;
  int nleft, nright, i, j, k;

  Trnum_T trnum;
  Trcoord_T troffset, trhigh;

  Chrnum_T chrnum;
  int transcript_genestrand;

  /* Previously computed middle_qstart and middle_qend are the maximum
     bounds for the entire elt.  But now, in order to set
     path->nmismatches to be 0, we are computing the minimum bounds */
  /* When we were computing maximum bounds for the entire elt, for an individual trdiagonal,
     they might go past the chromosome boundaries */

  /* left = middle_trdiagonal - querylength; */
  debug(printf("Entering process_seed with middle_trdiagonal %u\n",middle_trdiagonal));

#ifdef TRIM_AT_TRANSCRIPTOME_BOUNDS
  middle_qstart = (middle_trdiagonal + middle_qstart >= 0 + (Trcoord_T) querylength) ? middle_qstart : (int) (querylength - middle_trdiagonal);
  middle_qend = (middle_trdiagonal + middle_qend <= transcriptomelength + (Trcoord_T) querylength) ? middle_qend : (int) (transcriptomelength - middle_trdiagonal + querylength);
#endif

  assert(middle_trdiagonal - querylength + middle_qend <= transcriptomelength);

  trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,
		     middle_trdiagonal - querylength + middle_qstart,
		     middle_trdiagonal - querylength + middle_qend);
  
#ifdef TRIM_AT_TRANSCRIPT_BOUNDS
  middle_qstart = (middle_trdiagonal + middle_qstart >= troffset + (Trcoord_T) querylength) ? middle_qstart : (int) (troffset - left);
  middle_qend = (middle_trdiagonal + middle_qend <= trhigh + (Trcoord_T) querylength) ? middle_qend : (int) (trhigh - left);
#endif
  debug(printf("PROCESS SEED at %u, qstart %d, qend %d\n",
	       middle_trdiagonal,middle_qstart,middle_qend));


  start_low = subtract_bounded(middle_trdiagonal,/*minusterm*/(Trcoord_T) max_deletionlen);
  start_high = add_bounded(middle_trdiagonal,/*plusterm*/(Trcoord_T) max_insertionlen);

  end_low = subtract_bounded(middle_trdiagonal,/*minusterm*/(Trcoord_T) max_insertionlen);
  end_high = add_bounded(middle_trdiagonal,/*plusterm*/(Trcoord_T) max_deletionlen);

  debug(printf("Computing start %u..%u and end %u..%u\n",start_low,start_high,end_low,end_high));


  for (p = queryfwd_set; p != NULL; p = List_next(p)) {
    elt = (Tr_elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_trdiags(elt,start_low,start_high);
	for (k = 0; k < elt->ntrdiags; k++) {
	  left_trdiags = Trdiagpool_push_existing(left_trdiags,trdiagpool,elt->trdiags[k]
						  trdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing %u %d..%d onto left of query fwd trdiagonal %u\n",
		       elt->trdiags[k]->trdiagonal,elt->trdiags[k]->qstart,elt->trdiags[k]->qend,
		       middle_trdiagonal));
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_trdiags(elt,end_low,end_high);
	for (k = 0; k < elt->ntrdiags; k++) {
	  right_trdiags = Trdiagpool_push_existing(right_trdiags,trdiagpool,elt->trdiags[k]
						   trdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing %u %d..%d onto right of query fwd trdiagonal %u\n",
		       elt->trdiags[k]->trdiagonal,elt->trdiags[k]->qstart,elt->trdiags[k]->qend,
		       middle_trdiagonal));
	}
      }
    }
  }
  
  for (p = queryrev_set; p != NULL; p = List_next(p)) {
    elt = (Tr_elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_trdiags(elt,start_low,start_high);
	for (k = 0; k < elt->ntrdiags; k++) {
	  left_trdiags = Trdiagpool_push_existing(left_trdiags,trdiagpool,elt->trdiags[k]
						  trdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing %u %d..%d onto left of query rev trdiagonal %u\n",
		       elt->trdiags[k]->trdiagonal,elt->trdiags[k]->qstart,elt->trdiags[k]->qend,
		       middle_trdiagonal));
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_trdiags(elt,end_low,end_high);
	for (k = 0; k < elt->ntrdiags; k++) {
	  right_trdiags = Trdiagpool_push_existing(right_trdiags,trdiagpool,elt->trdiags[k]
						   trdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing %u %d..%d onto right of query rev trdiagonal %u\n",
		       elt->trdiags[k]->trdiagonal,elt->trdiags[k]->qstart,elt->trdiags[k]->qend,
		       middle_trdiagonal));
	}
      }
    }
  }

  /* At this point, left_trdiags and right_trdiags are a subset of those in Tr_elt_T objects */
  /* Below, the contents of left_trdiags and right_trdiags changes */

  if (left_trdiags != NULL) {
    /* Sort the left trdiagonals and unique them */
    left_trdiag_array = (Trdiag_T *) List_to_array_n(&nleft,left_trdiags);
    qsort(left_trdiag_array,nleft,sizeof(Trdiag_T),Trdiag_diagonal_cmp);
    Trdiagpool_free_list(&left_trdiags,trdiagpool
			 trdiagpool_trace(__FILE__,__LINE__)); /* allocated by Trdiagpool_push */

    /* Combine left trdiagonals to expand their qstart and qend range */
    left_trdiags = (List_T) NULL;
    i = 0;
    while (i < nleft) {
      trdiagonal = left_trdiag_array[i]->trdiagonal;
      qstart = left_trdiag_array[i]->qstart;
      qend = left_trdiag_array[i]->qend;
      
      j = i+1;
      while (j < nleft && left_trdiag_array[j]->trdiagonal == trdiagonal) {
	debug(printf("At left diagonal %u, combining %d..%d with %d..%d\n",
		     trdiagonal,qstart,qend,left_trdiag_array[j]->qstart,left_trdiag_array[j]->qend));
	if (left_trdiag_array[j]->qstart < qstart) {
	  qstart = left_trdiag_array[j]->qstart;
	}
	if (left_trdiag_array[j]->qend > qend) {
	  qend = left_trdiag_array[j]->qend;
	}
	j++;
      }
      
      if (qstart == left_trdiag_array[i]->qstart && qend == left_trdiag_array[i]->qend) {
#if 0
	/* This shortcut confuses existing Trdiag_T objects with newly created ones */
	left_trdiags = Trdiagpool_push_existing(left_trdiags,trdiagpool,left_trdiag_array[i]
						trdiagpool_trace(__FILE__,__LINE__));
#else
	left_trdiags = Trdiagpool_push(left_trdiags,trdiagpool,qstart,qend,
				       left_trdiag_array[i]->nmismatches,trdiagonal
				       trdiagpool_trace(__FILE__,__LINE__));
#endif
      } else {
	left_trdiags = Trdiagpool_push(left_trdiags,trdiagpool,qstart,qend,/*nmismatches*/-1,trdiagonal
				       trdiagpool_trace(__FILE__,__LINE__));
      }	

      i = j;
    }
    FREE(left_trdiag_array);
  }

  if (right_trdiags != NULL) {
    /* Sort the right trdiagonals and unique them */
    right_trdiag_array = (Trdiag_T *) List_to_array_n(&nright,right_trdiags);
    qsort(right_trdiag_array,nright,sizeof(Trdiag_T),Trdiag_diagonal_cmp);
    Trdiagpool_free_list(&right_trdiags,trdiagpool
			 trdiagpool_trace(__FILE__,__LINE__)); /* allocated by Trdiagpool_push */
    
    debug(printf("right trdiagonals is not NULL => %d\n",nright));

    /* Combine right trdiagonals to expand their qstart and qend range */
    right_trdiags = (List_T) NULL;
    i = 0;
    while (i < nright) {
      trdiagonal = right_trdiag_array[i]->trdiagonal;
      qstart = right_trdiag_array[i]->qstart;
      qend = right_trdiag_array[i]->qend;
      
      j = i+1;
      while (j < nright && right_trdiag_array[j]->trdiagonal == trdiagonal) {
	debug(printf("At right diagonal %u, combining %d..%d with %d..%d\n",
		     trdiagonal,qstart,qend,right_trdiag_array[j]->qstart,right_trdiag_array[j]->qend));
	if (right_trdiag_array[j]->qstart < qstart) {
	  qstart = right_trdiag_array[j]->qstart;
	}
	if (right_trdiag_array[j]->qend > qend) {
	  qend = right_trdiag_array[j]->qend;
	}
	j++;
      }
      
      if (qstart == right_trdiag_array[i]->qstart && qend == right_trdiag_array[i]->qend) {
#if 0
	/* This shortcut confuses existing Trdiag_T objects with newly created ones */
	right_trdiags = Trdiagpool_push_existing(right_trdiags,trdiagpool,right_trdiag_array[i]
						 trdiagpool_trace(__FILE__,__LINE__));
#else
	right_trdiags = Trdiagpool_push(right_trdiags,trdiagpool,qstart,qend,
					right_trdiag_array[i]->nmismatches,trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
#endif	
      } else {
	right_trdiags = Trdiagpool_push(right_trdiags,trdiagpool,qstart,qend,/*nmismatches*/-1,trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
      }	

      i = j;
    }
    FREE(right_trdiag_array);
  }

  /* At this point, left_trdiags and right_trdiags are new Trdiag_T objects, not in any Tr_elt_T object */


#ifdef DEBUG
  Trdiag_T trdiag;
  printf("Calling Trpath_solve with middle diagonal %u [%u], %d..%d, and %d left and %d right diagonals\n",
	 middle_trdiagonal,middle_trdiagonal - troffset,middle_qstart,middle_qend,
	 List_length(left_trdiags),List_length(right_trdiags));
  for (p = left_trdiags; p != NULL; p = List_next(p)) {
    trdiag = (Trdiag_T) List_head(p);
    printf("Left %u, %d..%d\n",trdiag->trdiagonal,trdiag->qstart,trdiag->qend);
  }
  for (p = right_trdiags; p != NULL; p = List_next(p)) {
    trdiag = (Trdiag_T) List_head(p);
    printf("Right %u, %d..%d\n",trdiag->trdiagonal,trdiag->qstart,trdiag->qend);
  }
#endif

  if (List_length(left_trdiags) == 1) {
    qstart_trdiag = List_head(left_trdiags);
  } else {
    qstart_trdiag = (Trdiag_T) NULL;
  }

  if (List_length(right_trdiags) == 1) {
    qend_trdiag = List_head(right_trdiags);
  } else {
    qend_trdiag = (Trdiag_T) NULL;
  }


  if (tplusp == true) {
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath = Trpath_solve_from_diagonals(&(*found_score),&(*end_score),&partialp,sufficient_score,
					      middle_trdiagonal,middle_qstart,middle_qend,
					      middle_nmismatches,qstart_trdiag,qend_trdiag,/*tplusp*/true,
					      querylength,/*query_compress_tr*/query_compress_fwd,
					      mismatch_positions_alloc,trnum,troffset,trhigh,
					      chrnum,transcript_genestrand,paired_end_p,first_read_p,
					      stage1->indelinfo,intlistpool,uintlistpool,listpool,
					      trpathpool,pathpool,method)) != NULL) {
      if (partialp == true) {
	/* Currently partialp always set to false */
	*partial_sense_trpaths = Hitlist_push(*partial_sense_trpaths,hitlistpool,(void *) trpath
					      hitlistpool_trace(__FILE__,__LINE__));
      } else {
	*sense_trpaths = Hitlist_push(*sense_trpaths,hitlistpool,(void *) trpath
				      hitlistpool_trace(__FILE__,__LINE__));
      }
    }

  } else {
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    if ((trpath = Trpath_solve_from_diagonals(&(*found_score),&(*end_score),&partialp,sufficient_score,
					      middle_trdiagonal,middle_qstart,middle_qend,
					      middle_nmismatches,qstart_trdiag,qend_trdiag,/*tplusp*/false,
					      querylength,/*query_compress_tr*/query_compress_rev,
					      mismatch_positions_alloc,trnum,troffset,trhigh,
					      chrnum,transcript_genestrand,paired_end_p,first_read_p,
					      stage1->indelinfo,intlistpool,uintlistpool,listpool,
					      trpathpool,pathpool,method)) != NULL) {
      if (partialp == true) {
	/* Currently partialp always set to false */
	*partial_antisense_trpaths = Hitlist_push(*partial_antisense_trpaths,hitlistpool,(void *) trpath
						  hitlistpool_trace(__FILE__,__LINE__));
      } else {
	*antisense_trpaths = Hitlist_push(*antisense_trpaths,hitlistpool,(void *) trpath
					  hitlistpool_trace(__FILE__,__LINE__));
      }
    }
  }

  Trdiagpool_gc(&right_trdiags,trdiagpool
		trdiagpool_trace(__FILE__,__LINE__)); /* allocated by Trdiagpool_push */
  Trdiagpool_gc(&left_trdiags,trdiagpool
		trdiagpool_trace(__FILE__,__LINE__)); /* allocated by Trdiagpool_push */

  return;
}


#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b


static void
extend_seeds_union (int *found_score, int *end_score, int sufficient_score,
		    List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
		    List_T *sense_trpaths, List_T *antisense_trpaths,
		    
		    List_T queryfwd_best_elts, List_T queryrev_best_elts,
		    List_T queryfwd_set, List_T queryrev_set,
		    
		    Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
		    
		    Compress_T query_compress_fwd, Compress_T query_compress_rev, bool tplusp, 
		    bool paired_end_p, bool first_read_p,
		    Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
		    Hitlistpool_T hitlistpool, Method_T method) {
  
  List_T trdiags = NULL, queryfwd_elt_list = NULL, queryrev_elt_list = NULL;
  Uintlist_T trdiagonals = NULL; /* To determine order and avoid duplicates */
  Trcoord_T *array;

  int *order;
  Trdiag_T *trdiag_array, queryfwd_trdiag, queryrev_trdiag, trdiag;

  int qstart, qend;
  int queryfwd_nseeds, queryrev_nseeds, n, i, j;
  List_T p, q;
  T *queryfwd_elt_array, *queryrev_elt_array, queryfwd_elt, queryrev_elt;

  
  debug(printf("Entering extend_seeds_union with %d queryfwd_best_elts and %d queryrev_best_elts\n",
	       List_length(queryfwd_best_elts),List_length(queryrev_best_elts)));

  /* Create a union of seeds */
  for (p = queryfwd_best_elts; p != NULL; p = List_next(p)) {
    queryfwd_elt = (T) List_head(p);
    queryfwd_nseeds = queryfwd_elt->n_all_trdiags;

    for (q = queryrev_best_elts; q != NULL; q = List_next(q)) {
      queryrev_elt = (T) List_head(q);
      queryrev_nseeds = queryrev_elt->n_all_trdiags;
    
#ifdef DEBUG
      printf("queryfwd_elt: ");
      Tr_elt_dump(queryfwd_elt);
      printf("\n");
      printf("queryrev_elt: ");
      Tr_elt_dump(queryrev_elt);
      printf("\n\n");
#endif

      i = j = 0;
      while (i < queryfwd_nseeds && j < queryrev_nseeds) {
	queryfwd_trdiag = queryfwd_elt->all_trdiags[i];
	queryrev_trdiag = queryrev_elt->all_trdiags[j];

	if (queryfwd_trdiag->trdiagonal == queryrev_trdiag->trdiagonal) {
	  /* Combine the queryfwd and queryrev trdiags */
	  debug(printf("At middle diagonal %u, combining #%d %d..%d with #%d %d..%d\n",
		       queryfwd_trdiag->trdiagonal,
		       i,queryfwd_trdiag->qstart,queryfwd_trdiag->qend,
		       j,queryrev_trdiag->qstart,queryrev_trdiag->qend));
	  qstart = min(queryfwd_trdiag->qstart,queryrev_trdiag->qstart);
	  qend = max(queryfwd_trdiag->qend,queryrev_trdiag->qend);
	  trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryfwd_trdiag->trdiagonal
					  uintlistpool_trace(__FILE__,__LINE__));
	  trdiags = Trdiagpool_push(trdiags,trdiagpool,qstart,qend,
					/*nmismatches*/-1,queryfwd_trdiag->trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));

	  i++; j++;

	} else if (queryfwd_trdiag->trdiagonal < queryrev_trdiag->trdiagonal) {
	  trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryfwd_trdiag->trdiagonal
					  uintlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Trdiag_T objects with those in Tr_elt_T objects */
	  trdiags = Trdiagpool_push_existing(trdiags,trdiagpool,queryfwd_trdiag
						 trdiagpool_trace(__FILE__,__LINE__));
#else
	  trdiags = Trdiagpool_push(trdiags,trdiagpool,queryfwd_trdiag->qstart,queryfwd_trdiag->qend,
					queryfwd_trdiag->nmismatches,queryfwd_trdiag->trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
#endif
	  i++;

	} else {
	  trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryrev_trdiag->trdiagonal
					  uintlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Trdiag_T objects with those in Tr_elt_T objects */
	  trdiags = Trdiagpool_push_existing(trdiags,trdiagpool,queryrev_trdiag
						 trdiagpool_trace(__FILE__,__LINE__));
#else
	  trdiags = Trdiagpool_push(trdiags,trdiagpool,queryrev_trdiag->qstart,queryrev_trdiag->qend,
					queryrev_trdiag->nmismatches,queryrev_trdiag->trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
#endif
	  j++;
	}

	queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) queryfwd_elt
					  listpool_trace(__FILE__,__LINE__));
	queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) queryrev_elt
					  listpool_trace(__FILE__,__LINE__));
      }

      if (i < queryfwd_nseeds) {
	/* qstart = queryfwd_elt->min_qstart; */
	/* qend = queryfwd_elt->max_qend; */
	while (i < queryfwd_nseeds) {
	  queryfwd_trdiag = queryfwd_elt->all_trdiags[i++];
	  
	  /* Allow trdiagonals from queryrev_elt */
	  trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryfwd_trdiag->trdiagonal
					  uintlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Trdiag_T objects with those in Tr_elt_T objects */
	  trdiags = Trdiagpool_push_existing(trdiags,trdiagpool,queryfwd_trdiag
						 trdiagpool_trace(__FILE__,__LINE__));
#else
	  trdiags = Trdiagpool_push(trdiags,trdiagpool,queryfwd_trdiag->qstart,queryfwd_trdiag->qend,
					queryfwd_trdiag->nmismatches,queryfwd_trdiag->trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
#endif

	  queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) queryfwd_elt
					    listpool_trace(__FILE__,__LINE__));
	  queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) NULL
					    listpool_trace(__FILE__,__LINE__));
	}
      }

      if (j < queryrev_nseeds) {
	/* qstart = queryrev_elt->min_qstart; */
	/* qend = queryrev_elt->max_qend; */
	while (j < queryrev_nseeds) {
	  queryrev_trdiag = queryrev_elt->all_trdiags[j++];
	  
	  /* Allow trdiagonals from queryfwd_elt */
	  trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryrev_trdiag->trdiagonal
					  uintlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Trdiag_T objects with those in Tr_elt_T objects */
	  trdiags = Trdiagpool_push_existing(trdiags,trdiagpool,queryrev_trdiag
						 trdiagpool_trace(__FILE__,__LINE__));
#else
	  trdiags = Trdiagpool_push(trdiags,trdiagpool,queryrev_trdiag->qstart,queryrev_trdiag->qend,
					queryrev_trdiag->nmismatches,queryrev_trdiag->trdiagonal
					trdiagpool_trace(__FILE__,__LINE__));
#endif
	  queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) NULL
					    listpool_trace(__FILE__,__LINE__));
	  queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) queryrev_elt
					    listpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }
      
  /* trdiags contains newly created Trdiag_T objects, not in Tr_elt_T objects */


  /* Want to sort three arrays simultaneously: trdiags,
     queryfwd_elts, and queryrev_elts.  Use order on the trdiagonals
     to do this. */
  n = Uintlist_length(trdiagonals);
  array = Uintlist_to_array(trdiagonals,/*end*/0);
  trdiag_array = (Trdiag_T *) List_to_array(trdiags,NULL);
  queryfwd_elt_array = (T *) List_to_array(queryfwd_elt_list,NULL);
  queryrev_elt_array = (T *) List_to_array(queryrev_elt_list,NULL);

  Listpool_free_list(&queryrev_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Listpool_free_list(&queryfwd_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Trdiagpool_free_list(&trdiags,trdiagpool
			 trdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Trdiagpool_push */
  Uintlistpool_free_list(&trdiagonals,uintlistpool
			 uintlistpool_trace(__FILE__,__LINE__)); /* Allocated by Uintlistpool_push */

  order = Sedgesort_order_uint4(array,n);
  
  /* trdiag_array contains newly created Trdiag_T objects, not in Tr_elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    trdiag = trdiag_array[order[i]];
    process_seed(&(*found_score),&(*end_score),sufficient_score,
		 &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		 &(*sense_trpaths),&(*antisense_trpaths),

		 trdiag->trdiagonal,trdiag->qstart,trdiag->qend,trdiag->nmismatches,
		 queryfwd_set,queryrev_set,
		 /*queryfwd_best_elt*/queryfwd_elt_array[order[i]],
		 /*queryrev_best_elt*/queryrev_elt_array[order[i]],
		 stage1,querylength,mismatch_positions_alloc,
		 query_compress_fwd,query_compress_rev,tplusp,
		 paired_end_p,first_read_p,		 
		 trdiagpool,intlistpool,uintlistpool,
		 listpool,trpathpool,pathpool,hitlistpool,method);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Trdiagpool_free_trdiag(&(trdiag_array[i]),trdiagpool
			   trdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(trdiag_array);

  FREE(order);
  FREE(queryrev_elt_array);
  FREE(queryfwd_elt_array);
  FREE(array);

  return;
}


static void
extend_seeds_queryfwd (int *found_score, int *end_score, int sufficient_score,
		       List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
		       List_T *sense_trpaths, List_T *antisense_trpaths,

		       List_T queryfwd_best_elts, List_T queryfwd_set,
		       
		       Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
		       
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, bool tplusp, 
		       bool paired_end_p, bool first_read_p,
		       Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
		       Hitlistpool_T hitlistpool, Method_T method) {
  
  List_T trdiags = NULL, queryfwd_elt_list = NULL;
  Uintlist_T trdiagonals = NULL; /* To determine order and avoid duplicates */
  Trcoord_T *array;

  int *order;
  Trdiag_T *trdiag_array, queryfwd_trdiag, trdiag;

  int queryfwd_nseeds, n, i, j;
  List_T p;
  T *queryfwd_elt_array, queryfwd_elt;

  
  debug(printf("Entering use_seeds_queryfwd with %d queryfwd_best_elts\n",
	       List_length(queryfwd_best_elts)));

  for (p = queryfwd_best_elts; p != NULL; p = List_next(p)) {
    queryfwd_elt = (T) List_head(p);
    queryfwd_nseeds = queryfwd_elt->n_all_trdiags;

#ifdef DEBUG
    printf("queryfwd_elt: ");
    Tr_elt_dump(queryfwd_elt);
    printf("\n");
#endif

    i = 0;
    while (i < queryfwd_nseeds) {
      queryfwd_trdiag = queryfwd_elt->all_trdiags[i++];
      trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryfwd_trdiag->trdiagonal
				      uintlistpool_trace(__FILE__,__LINE__));
      trdiags = Trdiagpool_push(trdiags,trdiagpool,queryfwd_trdiag->qstart,queryfwd_trdiag->qend,
				queryfwd_trdiag->nmismatches,queryfwd_trdiag->trdiagonal
				trdiagpool_trace(__FILE__,__LINE__));
      queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) queryfwd_elt
					listpool_trace(__FILE__,__LINE__));
    }
  }
      
  /* trdiags contains newly created Trdiag_T objects, not in Tr_elt_T objects */

  /* Want to sort two arrays simultaneously: trdiags and queryfwd_elts.
     Use order on the trdiagonals to do this. */
  n = Uintlist_length(trdiagonals);
  array = Uintlist_to_array(trdiagonals,/*end*/0);
  trdiag_array = (Trdiag_T *) List_to_array(trdiags,NULL);
  queryfwd_elt_array = (T *) List_to_array(queryfwd_elt_list,NULL);

  Listpool_free_list(&queryfwd_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Trdiagpool_free_list(&trdiags,trdiagpool
			 trdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Trdiagpool_push */
  Uintlistpool_free_list(&trdiagonals,uintlistpool
			 uintlistpool_trace(__FILE__,__LINE__)); /* Allocated by Uintlistpool_push */

  order = Sedgesort_order_uint4(array,n);
  
  /* trdiag_array contains newly created Trdiag_T objects, not in Tr_elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    trdiag = trdiag_array[order[i]];
    process_seed(&(*found_score),&(*end_score),sufficient_score,
		 &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		 &(*sense_trpaths),&(*antisense_trpaths),

		 trdiag->trdiagonal,trdiag->qstart,trdiag->qend,trdiag->nmismatches,
		 queryfwd_set,/*queryrev_set*/NULL,
		 /*queryfwd_best_elt*/queryfwd_elt_array[order[i]],
		 /*queryrev_best_elt*/NULL,
		 stage1,querylength,mismatch_positions_alloc,
		 query_compress_fwd,query_compress_rev,tplusp,
		 paired_end_p,first_read_p,		 
		 trdiagpool,intlistpool,uintlistpool,
		 listpool,trpathpool,pathpool,hitlistpool,method);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Trdiagpool_free_trdiag(&(trdiag_array[i]),trdiagpool
			   trdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(trdiag_array);

  FREE(order);
  FREE(queryfwd_elt_array);
  FREE(array);

  return;
}


static void
extend_seeds_queryrev (int *found_score, int *end_score, int sufficient_score,
		       List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
		       List_T *sense_trpaths, List_T *antisense_trpaths,
		    
		       List_T queryrev_best_elts, List_T queryrev_set,
		       
		       Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
		       
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, bool tplusp, 
		       bool paired_end_p, bool first_read_p,
		       Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
		       Hitlistpool_T hitlistpool, Method_T method) {
  
  List_T trdiags = NULL, queryrev_elt_list = NULL;
  Uintlist_T trdiagonals = NULL; /* To determine order and avoid duplicates */
  Trcoord_T *array;

  int *order;
  Trdiag_T *trdiag_array, queryrev_trdiag, trdiag;

  int queryrev_nseeds, n, i, j;
  List_T q;
  T *queryrev_elt_array, queryrev_elt;

  
  debug(printf("Entering extend_seeds_queryrev with %d queryrev_best_elts\n",
	       List_length(queryrev_best_elts)));

  for (q = queryrev_best_elts; q != NULL; q = List_next(q)) {
    queryrev_elt = (T) List_head(q);
    queryrev_nseeds = queryrev_elt->n_all_trdiags;
    
#ifdef DEBUG
    printf("queryrev_elt: ");
    Tr_elt_dump(queryrev_elt);
    printf("\n\n");
#endif

    j = 0;
    while (j < queryrev_nseeds) {
      queryrev_trdiag = queryrev_elt->all_trdiags[j++];
	  
      trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,queryrev_trdiag->trdiagonal
				      uintlistpool_trace(__FILE__,__LINE__));
      trdiags = Trdiagpool_push(trdiags,trdiagpool,queryrev_trdiag->qstart,queryrev_trdiag->qend,
				queryrev_trdiag->nmismatches,queryrev_trdiag->trdiagonal
				trdiagpool_trace(__FILE__,__LINE__));
      queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) queryrev_elt
					listpool_trace(__FILE__,__LINE__));
    }
  }
      
  /* trdiags contains newly created Trdiag_T objects, not in Tr_elt_T objects */

  /* Want to sort two arrays simultaneously: trdiags and
     queryrev_elts.  Use order on the trdiagonals to do this. */
  n = Uintlist_length(trdiagonals);
  array = Uintlist_to_array(trdiagonals,/*end*/0);
  trdiag_array = (Trdiag_T *) List_to_array(trdiags,NULL);
  queryrev_elt_array = (T *) List_to_array(queryrev_elt_list,NULL);

  Listpool_free_list(&queryrev_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Trdiagpool_free_list(&trdiags,trdiagpool
			 trdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Trdiagpool_push */
  Uintlistpool_free_list(&trdiagonals,uintlistpool
			 uintlistpool_trace(__FILE__,__LINE__)); /* Allocated by Uintlistpool_push */

  order = Sedgesort_order_uint4(array,n);
  
  /* trdiag_array contains newly created Trdiag_T objects, not in Tr_elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    trdiag = trdiag_array[order[i]];
    process_seed(&(*found_score),&(*end_score),sufficient_score,
		 &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		 &(*sense_trpaths),&(*antisense_trpaths),

		 trdiag->trdiagonal,trdiag->qstart,trdiag->qend,trdiag->nmismatches,
		 /*queryfwd_set*/NULL,queryrev_set,
		 /*queryfwd_best_elt*/NULL,
		 /*queryrev_best_elt*/queryrev_elt_array[order[i]],
		 stage1,querylength,mismatch_positions_alloc,
		 query_compress_fwd,query_compress_rev,tplusp,
		 paired_end_p,first_read_p,		 
		 trdiagpool,intlistpool,uintlistpool,
		 listpool,trpathpool,pathpool,hitlistpool,method);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Trdiagpool_free_trdiag(&(trdiag_array[i]),trdiagpool
			   trdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(trdiag_array);

  FREE(order);
  FREE(queryrev_elt_array);
  FREE(array);

  return;
}


void
Tr_extension_search (int *found_score, int *end_score, int sufficient_score,
		     List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
		     List_T *sense_trpaths, List_T *antisense_trpaths,

		     Stage1_T stage1, int querylength, int *mismatch_positions_alloc,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     bool paired_end_p, bool first_read_p,

		     Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
		     Hitlistpool_T hitlistpool,

		     int read_nmismatches_allowed, int genestrand, Method_T method) {

  List_T queryfwd_best_plus_elts, queryrev_best_plus_elts, queryfwd_best_minus_elts, queryrev_best_minus_elts;


  get_elt_sets_queryfwd(&stage1->tr_queryfwd_plus_set,&stage1->tr_queryfwd_minus_set,
			&queryfwd_best_plus_elts,&queryfwd_best_minus_elts,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			read_nmismatches_allowed,genestrand,listpool,trdiagpool);
  
  get_elt_sets_queryrev(&stage1->tr_queryrev_plus_set,&stage1->tr_queryrev_minus_set,
			&queryrev_best_plus_elts,&queryrev_best_minus_elts,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			read_nmismatches_allowed,genestrand,listpool,trdiagpool);

  if (queryfwd_best_plus_elts == NULL && queryrev_best_plus_elts == NULL) {
    /* Skip */
  } else if (queryrev_best_plus_elts == NULL) {
    extend_seeds_queryfwd(&(*found_score),&(*end_score),sufficient_score,
			  &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			  &(*sense_trpaths),&(*antisense_trpaths),

			  queryfwd_best_plus_elts,stage1->tr_queryfwd_plus_set,
			  stage1,querylength,mismatch_positions_alloc,
			  query_compress_fwd,query_compress_rev,/*tplusp*/true,
			  paired_end_p,first_read_p,
			  trdiagpool,intlistpool,uintlistpool,listpool,
			  trpathpool,pathpool,hitlistpool,method);
  } else if (queryfwd_best_plus_elts == NULL) {
    extend_seeds_queryrev(&(*found_score),&(*end_score),sufficient_score,
			  &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			  &(*sense_trpaths),&(*antisense_trpaths),

			  queryrev_best_plus_elts,stage1->tr_queryrev_plus_set,
			  stage1,querylength,mismatch_positions_alloc,
			  query_compress_fwd,query_compress_rev,/*tplusp*/true,
			  paired_end_p,first_read_p,
			  trdiagpool,intlistpool,uintlistpool,listpool,
			  trpathpool,pathpool,hitlistpool,method);
  } else {
    extend_seeds_union(&(*found_score),&(*end_score),sufficient_score,
		       &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		       &(*sense_trpaths),&(*antisense_trpaths),

		       queryfwd_best_plus_elts,queryrev_best_plus_elts,
		       stage1->tr_queryfwd_plus_set,stage1->tr_queryrev_plus_set,
		       stage1,querylength,mismatch_positions_alloc,
		       query_compress_fwd,query_compress_rev,/*tplusp*/true,
		       paired_end_p,first_read_p,
		       trdiagpool,intlistpool,uintlistpool,listpool,
		       trpathpool,pathpool,hitlistpool,method);
  }

  if (queryfwd_best_minus_elts == NULL && queryrev_best_minus_elts == NULL) {
    /* Skip */
  } else if (queryrev_best_minus_elts == NULL) {
    extend_seeds_queryfwd(&(*found_score),&(*end_score),sufficient_score,
			  &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			  &(*sense_trpaths),&(*antisense_trpaths),

			  queryfwd_best_minus_elts,stage1->tr_queryfwd_minus_set,
			  stage1,querylength,mismatch_positions_alloc,
			  query_compress_fwd,query_compress_rev,/*tplusp*/false,
			  paired_end_p,first_read_p,
			  trdiagpool,intlistpool,uintlistpool,listpool,
			  trpathpool,pathpool,hitlistpool,method);
  } else if (queryfwd_best_minus_elts == NULL) {
    extend_seeds_queryrev(&(*found_score),&(*end_score),sufficient_score,
			  &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			  &(*sense_trpaths),&(*antisense_trpaths),

			  queryrev_best_minus_elts,stage1->tr_queryrev_minus_set,
			  stage1,querylength,mismatch_positions_alloc,
			  query_compress_fwd,query_compress_rev,/*tplusp*/false,
			  paired_end_p,first_read_p,
			  trdiagpool,intlistpool,uintlistpool,listpool,
			  trpathpool,pathpool,hitlistpool,method);
  } else {
    extend_seeds_union(&(*found_score),&(*end_score),sufficient_score,
		       &(*partial_sense_trpaths),&(*partial_antisense_trpaths),
		       &(*sense_trpaths),&(*antisense_trpaths),

		       queryfwd_best_minus_elts,queryrev_best_minus_elts,
		       stage1->tr_queryfwd_minus_set,stage1->tr_queryrev_minus_set,
		       stage1,querylength,mismatch_positions_alloc,
		       query_compress_fwd,query_compress_rev,/*tplusp*/false,
		       paired_end_p,first_read_p,
		       trdiagpool,intlistpool,uintlistpool,listpool,
		       trpathpool,pathpool,hitlistpool,method);
  }

  /* Allocated by listpool */
  Listpool_free_list(&queryfwd_best_plus_elts,listpool
		     listpool_trace(__FILE__,__LINE__));
  Listpool_free_list(&queryrev_best_plus_elts,listpool
		     listpool_trace(__FILE__,__LINE__));
  Listpool_free_list(&queryfwd_best_minus_elts,listpool
		     listpool_trace(__FILE__,__LINE__));
  Listpool_free_list(&queryrev_best_minus_elts,listpool
		     listpool_trace(__FILE__,__LINE__));

  return;
}


