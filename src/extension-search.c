static char rcsid[] = "$Id: 832a78249c406a4a68734ba5e5e687817efac2c0 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "extension-search.h"


/***********************************************************************
 *   Definitions:
 *
 *   For plus:
 *
 *   ACGT .................... ACGT
 *   queryoffset, distance from beginning of read
 *   | ------> | 
 *
 *
 *   For minus:
 *
 *   ACGT .................... ACGT
 *                   queryoffset, distance from end of read (actually query_lastpos)
 *                   | <------ |
 *                            query_lastpos = querylength - index1part
 *
 *
 *   In both cases, diagterm = querylength - queryoffset
 *   For plus, diagterm = querylength - querypos
 *   For minus, diagterm = querylength - (query_lastpos - querypos_rc)
 *                       = querypos_rc + index1part
 *
 *   In Stage1_T object:
 *        forward_oligo
 *        at [querypos], where querypos = queryoffset
 * 
 *                revcomp_oligo
 *                at [querypos_rc], where querypos_rc = query_lastpos - queryoffset
 *
 *   queryoffset should be the same as qpos
 *
 ************************************************************************/

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
#include "univcoord.h"
#include "oligo.h"

#include "list.h"
#include "maxent_hr.h"

#include "genomebits_consec.h"
#include "genomebits_trim.h"

#include "univdiagdef.h"
#include "univdiag.h"
#include "univdiagpool.h"
#include "sedgesort.h"

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
/* Conservative appears to be 20% faster than aggressive */
#define CONSERVATIVE 1

#define TRIM_AT_GENOME_BOUNDS 1
#define MAX_NPOSITIONS 1000
#define MAX_TOTAL_NPOSITIONS 3000
#define MAX_INDEX1INTERVAL 6


/* #define CHECK_OLIGOS 1 */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Details */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
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


static Mode_T mode;

static Univcoord_T genomelength;
static int circular_typeint;
static bool *circularp;
static EF64_T chromosome_ef64;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Indexdb_T indexdb;
static Indexdb_T indexdb_nonstd;

static Chrpos_T negative_gap_distance; /* typically max_insertionlen */
static Chrpos_T positive_gap_distance;


static char conversion_fwd[128];
static char conversion_rev[128];

static int index1part;
static int index1interval;
static int leftreadshift;
static Oligospace_T oligobase_mask;

/* Some limit is needed to prevent GSNAP from running very slowly */
/* Was MAX_HITS_FOR_BEST_ELT 1000 */

static int maxpaths_search;


#ifdef LARGE_GENOMES
#define GETPOS(high,low) (((Univcoord_T) high << 32) + low)
#endif


void
Extension_search_setup (Mode_T mode_in,
			Univcoord_T genomelength_in, int circular_typeint_in, bool *circularp_in, EF64_T chromosome_ef64_in,
			Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
			Indexdb_T indexdb_in, Indexdb_T indexdb_nonstd_in,
			int index1part_in, int index1interval_in, int maxpaths_search_in,
			int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist) {
  int i;

  mode = mode_in;

  genomelength = genomelength_in;
  circular_typeint = circular_typeint_in;
  circularp = circularp_in;
  chromosome_ef64 = chromosome_ef64_in;

  indexdb = indexdb_in;
  indexdb_nonstd = indexdb_nonstd_in;

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  for (i = 0; i < 128; i++) {
    conversion_fwd[i] = i;
    conversion_rev[i] = i;
  }
  if (mode == STANDARD) {
    /* Don't change conversion */
  } else if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    conversion_fwd['C'] = 'T';	/* CT */
    conversion_rev['G'] = 'A';	/* GA */
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    conversion_fwd['A'] = 'G';	/* AG */
    conversion_rev['T'] = 'C';	/* TC */
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    conversion_fwd['T'] = 'C';	/* TC */
    conversion_rev['A'] = 'G';	/* AG */
  }

  index1part = index1part_in;
  index1interval = index1interval_in;

  negative_gap_distance = (Chrpos_T) max_insertionlen;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;


#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
#endif

  maxpaths_search = maxpaths_search_in;

  return;
}


/* Simplified version of Spanningelt_T */
#define T Elt_T


#if 0
/* Works only for queryfwd.  Use Elt_read functions instead */
static T
Elt_new (int querypos, int nmatches, Univcoord_T *all_diagonals, int n_all_diagonals) {
  T new = (T) MALLOC(sizeof(*new));
#ifdef DEBUG
  int i;
#endif

  new->min_qstart = querypos;
  new->max_qend = querypos + nmatches - 1;
  new->nmatches = nmatches;
  debug(printf("Making an Elt with querystart %d, nmatches %d => queryend %d\n",
	       new->qstart,new->nmatches,new->qend));

  new->all_diagonals = all_diagonals;
  new->n_all_diagonals = n_all_diagonals;

  new->diagonals = &(new->all_diagonals[0]);
  new->ndiagonals = n_all_diagonals;

#ifdef DEBUG
  printf("Diagonals:");
  for (i = 0; i < n_all_diagonals; i++) {
    printf(" %llu",all_diagonals[i]);
  }
  printf("\n");
#endif

  new->lowi = 0;
  new->highi = n_all_diagonals;

  return new;
}
#endif


static void
Elt_free (T *old, Univdiagpool_T univdiagpool) {
  int i;

  for (i = 0; i < (*old)->n_all_univdiags; i++) {
    Univdiagpool_free_univdiag(&(*old)->all_univdiags[i],univdiagpool
			       univdiagpool_trace(__FILE__,__LINE__));
  }
  FREE((*old)->all_univdiags);
  FREE(*old);
  return;
}

void
Elt_gc (List_T *set, Listpool_T listpool, Univdiagpool_T univdiagpool) {
  List_T p;
  T elt;

  for (p = *set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    Elt_free(&elt,univdiagpool);
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
/* When plusp is false, need to read from querypos_rc */
static T
Elt_read_queryfwd (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int querylength, int querystart,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Univdiagpool_T univdiagpool) {
  T new;
  int max_nmatches;
  Univdiag_T *out, *best_univdiags, *univdiags;
  Univcoord_T univdiagonal;

  int best_trimpos, trimpos;

  int *nmatches, pos3;
  int i;
  

  debug9(printf("Got %d positions at querystart %d, plusp %d\n",n,querystart,plusp));
  if (n == 0) {
    return (T) NULL;
  } else if (querystart >= querylength - index1part) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = querystart;

    nmatches = (int *) MALLOC(n*sizeof(int));
    univdiags = (Univdiag_T *) MALLOC(n*sizeof(Univdiag_T));
    for (i = 0; i < n; i++) {
#ifdef LARGE_GENOMES
      univdiagonal = GETPOS(positions_high[i],positions[i]) + diagterm;
#else
      univdiagonal = positions[i] + diagterm;
#endif


#ifdef TRIM_AT_GENOME_BOUNDS
      pos3 = (univdiagonal <= genomelength) ? querylength : (int) (genomelength - univdiagonal + querylength);
#else
      pos3 = querylength;
#endif

#ifdef CONSERVATIVE
      /* Use Genome_consecutive_matches_rightward, which is conservative */
      nmatches[i] = index1part +
	Genomebits_consecutive_matches_rightward(genomebits,query_compress,univdiagonal,querylength,
						 /*pos5*/querystart+index1part,pos3,
						 plusp,genestrand);
      univdiags[i] = Univdiagpool_new_univdiag(univdiagpool,querystart,/*qend*/querystart+nmatches[i],
					       /*nmismatches*/0,univdiagonal
					       univdiagpool_trace(__FILE__,__LINE__));
      trimpos = querystart + nmatches[i];
      debug9(printf("rightward plusp %d, univdiagonal is %u, querystart is %d + nmatches %d => queryend is %d\n",
		    plusp,univdiagonal,querystart,nmatches[i],querystart + nmatches[i]));

#else
      /* Try Genomebits_trim_qend, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qend(&nmismatches,query_compress,
				     genomebits,univdiagonal,querylength,/*pos5*/querystart,pos3,
				     plusp,genestrand);
      debug9(printf("qend conservative: At univdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		    univdiagonal,querystart,pos3,nmatches[i],querystart+nmatches[i]));
      debug9(printf("qend aggressive: At univdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		    univdiagonal,querystart,pos3,trimpos,nmismatches,(trimpos - querystart) - nmismatches,trimpos));
      nmatches[i] = (trimpos - querystart) - nmismatches;
      univdiags[i] = Univdiagpool_new_univdiag(univdiagpool,querystart,/*qend*/trimpos,
					       nmismatches,univdiagonal
					       univdiagpool_trace(__FILE__,__LINE__));
#endif

      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
	best_trimpos = trimpos;
      }
    }

    out = best_univdiags = MALLOC(n*sizeof(Univdiag_T));
    for (i = 0; i < n; i++) {
      /* Accept diagonals within a certain threshold of the maximum */
      if (nmatches[i] >= max_nmatches - index1interval + 1) {
	*out++ = univdiags[i];
      } else {
	Univdiagpool_free_univdiag(&(univdiags[i]),univdiagpool
				   univdiagpool_trace(__FILE__,__LINE__));
      }
    }
    FREE(univdiags);
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->min_qstart = querystart;
    /* new->max_qend = querystart + max_nmatches; */
    new->max_qend = best_trimpos;
    new->nmatches = max_nmatches;

    assert(new->max_qend <= (int) querylength);
    debug(printf("Making a queryfwd Elt %p with querystart %d, max_nmatches %d => max_queryend %d\n",
		 new,new->min_qstart,new->nmatches,new->max_qend));

    new->all_univdiags = new->univdiags = best_univdiags;
    new->n_all_univdiags = new->nunivdiags = out - best_univdiags;

#ifdef DEBUG
    printf("Univdiags:");
    for (i = 0; i < new->nunivdiags; i++) {
      printf("(%p) %u %d..%d",
	     new,new->univdiags[i]->univdiagonal,new->univdiags[i]->qstart,new->univdiags[i]->qend);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_univdiags;

    return new;
  }
}


static T
Elt_read_queryrev (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int queryend, int querylength,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Univdiagpool_T univdiagpool) {
  T new;
  int max_nmatches;
  Univdiag_T *out, *best_univdiags, *univdiags;
  Univcoord_T univdiagonal;

  int best_trimpos, trimpos;

  int *nmatches, pos5;
  int i;
  

  debug9(printf("Got %d positions at queryend %d, plusp %d\n",n,queryend,plusp));
  if (n == 0) {
    return (T) NULL;
  } else if (queryend < index1part) {
    return (T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = queryend;

    nmatches = (int *) MALLOC(n*sizeof(int));
    univdiags = (Univdiag_T *) MALLOC(n*sizeof(Univdiag_T));
    for (i = 0; i < n; i++) {
#ifdef LARGE_GENOMES
      univdiagonal = GETPOS(positions_high[i],positions[i]) + diagterm;
#else
      univdiagonal = positions[i] + diagterm;
#endif

#ifdef TRIM_AT_GENOME_BOUNDS
      pos5 = (univdiagonal >= 0 + (Univcoord_T) querylength) ? 0 : (int) (querylength - univdiagonal);
#else
      pos5 = 0;
#endif
      
#ifdef CONSERVATIVE
      /* Use Genome_consecutive_matches_leftward, which is conservative */
      nmatches[i] = index1part +
	Genomebits_consecutive_matches_leftward(genomebits,query_compress,univdiagonal,querylength,
						pos5,/*pos3*/queryend - index1part,plusp,genestrand);
      univdiags[i] = Univdiagpool_new_univdiag(univdiagpool,/*qstart*/queryend - nmatches[i],queryend,
					       /*nmismatches*/0,univdiagonal
					       univdiagpool_trace(__FILE__,__LINE__));
      trimpos = queryend - nmatches[i];
      debug9(printf("leftward plusp %d, univdiagonal is %u, queryend is %d - nmatches %d => querystart is %d\n",
		    plusp,univdiagonal,queryend,nmatches[i],queryend - nmatches[i]));

#else
      /* Try Genomebits_trim_qstart, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qstart(&nmismatches,query_compress,
				       genomebits,univdiagonal,querylength,pos5,/*pos3*/queryend,
				       plusp,genestrand);
      debug9(printf("qstart conservative: At univdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		    univdiagonal,pos5,queryend,nmatches[i],queryend-nmatches[i]));
      debug9(printf("qstart aggressive: At univdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		    univdiagonal,pos5,queryend,trimpos,nmismatches,(queryend - trimpos) - nmismatches,trimpos));
      nmatches[i] = (queryend - trimpos) - nmismatches;
      univdiags[i] = Univdiagpool_new_univdiag(univdiagpool,/*qstart*/trimpos,queryend,
					       nmismatches,univdiagonal
					       univdiagpool_trace(__FILE__,__LINE__));
#endif

      if (nmatches[i] > max_nmatches) {
	max_nmatches = nmatches[i];
	best_trimpos = trimpos;
      }
    }

    out = best_univdiags = MALLOC(n*sizeof(Univdiag_T));
    for (i = 0; i < n; i++) {
      /* Accept diagonals within a certain threshold of the maximum */
      if (nmatches[i] >= max_nmatches - index1interval + 1) {
	*out++ = univdiags[i];
      } else {
	Univdiagpool_free_univdiag(&(univdiags[i]),univdiagpool
				   univdiagpool_trace(__FILE__,__LINE__));
      }
    }
    FREE(univdiags);
    FREE(nmatches);
   

    new = (T) MALLOC(sizeof(*new));
    new->max_qend = queryend;
    /* new->min_qstart = queryend - max_nmatches; */
    new->min_qstart = best_trimpos;
    new->nmatches = max_nmatches;

    assert(new->min_qstart >= 0);
    debug(printf("Making a queryrev Elt %p with queryend %d, max_nmatches %d => min_querystart %d\n",
		 new,new->max_qend,new->nmatches,new->min_qstart));

    new->all_univdiags = new->univdiags = best_univdiags;
    new->n_all_univdiags = new->nunivdiags = out - best_univdiags;

#ifdef DEBUG
    printf("Univdiagonals:");
    for (i = 0; i < new->nunivdiags; i++) {
      printf(" %u %d..%d",
	     new->univdiags[i]->univdiagonal,new->univdiags[i]->qstart,new->univdiags[i]->qend);
    }
    printf("\n");
#endif

    new->lowi = 0;
    new->highi = new->n_all_univdiags;

    return new;
  }
}


static int
binary_search_univdiag (int lowi, int highi, Univdiag_T *univdiags, Univcoord_T goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,univdiags[lowi]->univdiagonal,middlei,univdiags[middlei]->univdiagonal,
		   highi-1,univdiags[highi-1]->univdiagonal,goal));
    if (goal < univdiags[middlei]->univdiagonal) {
      highi = middlei;
    } else if (goal > univdiags[middlei]->univdiagonal) {
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
Elt_filter_univdiags (T this, Univcoord_T low, Univcoord_T high) {
  int lowi, highi;
#ifdef DEBUG
  int i;
#endif

  debug(printf("Entered Elt_filter_univdiags on %d..%d with %d univdiags with low %u and high %u, nmatches %d\n",
	       this->min_qstart,this->max_qend,this->n_all_univdiags,low,high,this->nmatches));
#ifdef DEBUG
  for (i = 0; i < this->n_all_univdiags; i++) {
    printf("Input: %u\n",this->all_univdiags[i]->univdiagonal);
  }
#endif

  /* low_adj and high_adj are inclusive */
  lowi = binary_search_univdiag(/*lowi*/0,/*highi*/this->n_all_univdiags,this->all_univdiags,/*goal*/low);
  highi = binary_search_univdiag(lowi,/*highi*/this->n_all_univdiags,this->all_univdiags,/*goal*/high + 1) - 1;
  if ((this->nunivdiags = highi - lowi + 1) == 0) {
    this->univdiags = (Univdiag_T *) NULL;

  } else {
    this->univdiags = &(this->all_univdiags[lowi]);
  }

#ifdef DEBUG
  printf("Setting lowi %d and highi %d\n",lowi,highi);
  for (i = lowi; i <= highi; i++) {
    printf("  %u %d..%d\n",
	   this->all_univdiags[i]->univdiagonal,this->all_univdiags[i]->qstart,this->all_univdiags[i]->qend);
  }
#endif

  return;
}


#ifdef DEBUG
static void
Elt_dump (T elt) {
  int k;

  printf("Elt %p with max querybounds %d..%d and %d univdiags:\n",
	 elt,elt->min_qstart,elt->max_qend,elt->nunivdiags);
  for (k = 0; k < elt->nunivdiags; k++) {
    printf("  (%p) %u %d..%d\n",
	   elt->univdiags[k],elt->univdiags[k]->univdiagonal,elt->univdiags[k]->qstart,elt->univdiags[k]->qend);
  }
  printf("\n");

  return;
}

static void
Elt_dump_set (List_T set) {
  List_T p;

  for (p = set; p != NULL; p = List_next(p)) {
    Elt_dump((T) List_head(p));
  }

  return;
}
#endif


#ifdef TRIM_AT_CHROMOSOME_BOUNDS
/* lowbound was chroffset; highbound was chrhigh */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))
#else
#define add_bounded(x,plusterm) (x + (plusterm))
#define subtract_bounded(x,minusterm) ((x < (minusterm)) ? 0 : x - (minusterm))
#endif


static bool
elt_startp (T elt, int middle_qstart, int middle_qend) {
  if (elt->min_qstart >= middle_qstart && elt->max_qend <= middle_qend) {
    debug(printf("Not allowing left elt that is subsumed by middle elt: q %d..%d\n",
		 elt->min_qstart,elt->max_qend));
    return false;
  } else if (elt->max_qend >= middle_qend) {
    debug(printf("Not allowing left elt that extends right of middle elt: q %d..%d\n",
		 elt->min_qstart,elt->max_qend));
    return false;
#if 0
  } else if ((elt->max_qend - middle_qstart) > (middle_qend - middle_qstart) / 2) {
    debug(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		 elt->min_qstart,elt->max_qend));
    return false;
#endif
#if 0
  } else if ((elt->max_qend - middle_qstart) > (elt->max_qend - elt->min_qstart) / 2) {
    debug(printf("Not allowing left elt that mainly overlaps middle elt: q %d..%d\n",
		 elt->min_qstart,elt->max_qend));
    return false;
#endif
  } else {
    return true;
  }
}


static bool
elt_endp (T elt, int middle_qstart, int middle_qend) {
  if (elt->min_qstart >= middle_qstart && elt->max_qend <= middle_qend) {
    debug(printf("Not allowing right elt that is subsumed by middle elt: qpos %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
  } else if (elt->min_qstart <= middle_qstart) {
    debug(printf("Not allowing right elt that extends left of middle elt: qpos %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
#if 0
  } else if ((middle_qend - elt->min_qstart) > (middle_qend - middle_qstart) / 2) {
    debug(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
#endif
#if 0
  } else if ((middle_qend - elt->min_qstart) > (elt->max_qend - elt->min_qstart) / 2) {
    debug(printf("Not allowing right elt that mainly overlaps middle elt: q %d..%d\n",
		   elt->min_qstart,elt->max_qend));
    return false;
#endif
  } else {
    return true;
  }
}


/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryfwd (List_T *plus_set, List_T *minus_set, List_T *best_plus_elts, List_T *best_minus_elts,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int max_iter, int genestrand, Listpool_T listpool,
		       Univdiagpool_T univdiagpool) {
  T elt;
  List_T p;
  int best_plus_nmatches, best_minus_nmatches;
  int max_plus_qpos, max_minus_qpos, plus_qpos[MAX_INDEX1INTERVAL], minus_qpos[MAX_INDEX1INTERVAL];

  int mod;
  int niter_plus, niter_minus;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
#endif
  UINT4 *positions;
  int total_npositions_plus = 0, total_npositions_minus = 0, mod_npositions, npositions;

  int query_lastpos = querylength - index1part;
  int queryoffset, querypos_rc;


  debug(printf("\nStarting get_elt_sets_queryfwd with querylength %d and max_iter %d, genestrand %d\n",
	       querylength,max_iter,genestrand));

#if 0
  /* Allow calls to Extension_search_queryfwd in addition to other methods */
  *paths_gplus = *paths_gminus = (List_T) NULL;
#endif


  /* I.  Race from plus and minus start to end.  Compute best_plus_nmatches and best_plus_nmatches */
  best_plus_nmatches = best_minus_nmatches = 0;
  max_plus_qpos = max_minus_qpos = 0;
  plus_qpos[0] = minus_qpos[0] = 0;
  plus_qpos[1] = minus_qpos[1] = 1;
  plus_qpos[2] = minus_qpos[2] = 2;
  niter_plus = niter_minus = 0;

  while ((total_npositions_plus <= MAX_TOTAL_NPOSITIONS ||
	  total_npositions_minus <= MAX_TOTAL_NPOSITIONS) &&
	 niter_plus <= max_iter && niter_minus <= max_iter &&
	 max_plus_qpos < query_lastpos && max_minus_qpos < query_lastpos) {

    if (total_npositions_plus <= MAX_TOTAL_NPOSITIONS) {
      mod_npositions = 0;
      for (mod = 0; mod < index1interval; mod++) {
	debug(printf("queryfwd plus mod %d",mod));
	if ((queryoffset = plus_qpos[mod]) >= query_lastpos) {
	  /* Skip */
	  debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	  
	} else if (stage1->validp[queryoffset] == false) {
	  /* Skip */
	  debug(printf(" => queryoffset %d not valid\n",queryoffset));
	  /* plus_qpos[mod] += index1interval; */
	  
	} else if (stage1->plus_retrievedp[queryoffset] == true) {
	  mod_npositions += stage1->plus_npositions[queryoffset];
	  debug(printf(" => saw %d positions\n",stage1->plus_npositions[queryoffset]));
	    
	} else {
	  assert(stage1->plus_positions[queryoffset] == NULL);
#ifdef LARGE_GENOMES
	  mod_npositions += stage1->plus_npositions[queryoffset] =
	    Indexdb_largeptr(&stage1->plus_positions_high[queryoffset],&stage1->plus_positions[queryoffset],
			     indexdb,stage1->forward_oligos[queryoffset]);
#else
	  mod_npositions += stage1->plus_npositions[queryoffset] =
	    Indexdb_ptr(&stage1->plus_positions[queryoffset],indexdb,
			stage1->forward_oligos[queryoffset]);
#endif
	  stage1->plus_retrievedp[queryoffset] = true;
	  debug(printf(" => retrieved %d positions\n",stage1->plus_npositions[queryoffset]));
	}
      }
	  
      debug(printf("mod positions plus: %d\n",mod_npositions));
      if (mod_npositions > MAX_NPOSITIONS) {
	/* Skip all */
	for (mod = 0; mod < index1interval; mod++) {
	  plus_qpos[mod] += index1interval;
	}

      } else {
	for (mod = 0; mod < index1interval; mod++) {
	  debug(printf("queryfwd plus mod %d",mod));
	  if ((queryoffset = plus_qpos[mod]) >= query_lastpos) {
	    /* Skip */
	    debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	  
	  } else if (stage1->validp[queryoffset] == false) {
	    /* Skip */
	    debug(printf(" => queryoffset %d not valid\n",queryoffset));
	    plus_qpos[mod] += index1interval;
	  
	  } else {
	    assert(stage1->plus_retrievedp[queryoffset] == true);
#ifdef LARGE_GENOMES
	    positions_high = stage1->plus_positions_high[queryoffset];
#endif
	    positions = stage1->plus_positions[queryoffset];
	    npositions = stage1->plus_npositions[queryoffset];
	    debug(printf(" => saw %d positions\n",npositions));

	    if ((elt = Elt_read_queryfwd(
#ifdef LARGE_GENOMES
					 positions_high,
#endif
					 positions,npositions,/*diagterm*/querylength - queryoffset,
					 querylength,/*querystart*/queryoffset,
					 query_compress_fwd,/*plusp*/true,genestrand,
					 univdiagpool)) == NULL) {
	      plus_qpos[mod] += index1interval;
	    
	    } else {
	      if (elt->nmatches > best_plus_nmatches && elt->n_all_univdiags <= maxpaths_search &&
		  elt->nunivdiags > 0) { /* Could be 0 if there are too many trdiags */
		best_plus_nmatches = elt->nmatches;
	      }

	      debug(printf("Pushing elt %p onto plus set\n",elt));
	      *plus_set = Listpool_push(*plus_set,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
	      plus_qpos[mod] += elt->nmatches;
	    }
	  }

	  if (plus_qpos[mod] > max_plus_qpos) {
	    max_plus_qpos = plus_qpos[mod];
	  }
	}

	total_npositions_plus += mod_npositions;
	niter_plus++;

      }
    }

    if (total_npositions_minus <= MAX_TOTAL_NPOSITIONS) {
      mod_npositions = 0;
      for (mod = 0; mod < index1interval; mod++) {
	debug(printf("queryfwd minus mod %d",mod));
	if ((queryoffset = minus_qpos[mod]) >= query_lastpos) {
	  /* Skip */
	  debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	  
	} else if (stage1->validp[(querypos_rc = query_lastpos - queryoffset)] == false) {
	  /* Skip */
	  debug(printf(" => querypos_rc %d not valid\n",querypos_rc));
	  /* minus_qpos[mod] += index1interval; */
	  
	} else if (stage1->minus_retrievedp[querypos_rc] == true) {
	  mod_npositions += stage1->minus_npositions[querypos_rc];
	  debug(printf(" => saw %d positions\n",stage1->minus_npositions[querypos_rc]));
	    
	} else {
	  assert(stage1->minus_positions[querypos_rc] == NULL);
#ifdef LARGE_GENOMES
	  mod_npositions += stage1->minus_npositions[querypos_rc] =
	    Indexdb_largeptr(&stage1->minus_positions_high[querypos_rc],&stage1->minus_positions[querypos_rc],
			     indexdb_nonstd,stage1->revcomp_oligos[querypos_rc]);
#else
	  mod_npositions = stage1->minus_npositions[querypos_rc] =
	    Indexdb_ptr(&stage1->minus_positions[querypos_rc],indexdb_nonstd,
			stage1->revcomp_oligos[querypos_rc]);
#endif
	  stage1->minus_retrievedp[querypos_rc] = true;
	  debug(printf(" => retrieved %d positions\n",stage1->minus_npositions[querypos_rc]));
	}
      }
	  
      debug(printf("mod positions minus: %d\n",mod_npositions));
      if (mod_npositions > MAX_NPOSITIONS) {
	/* Skip all */
	for (mod = 0; mod < index1interval; mod++) {
	  minus_qpos[mod] += index1interval;
	}

      } else {
	for (mod = 0; mod < index1interval; mod++) {
	  debug(printf("queryfwd minus mod %d",mod));
	  if ((queryoffset = minus_qpos[mod]) >= query_lastpos) {
	    /* Skip */
	    debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	    
	  } else if (stage1->validp[(querypos_rc = query_lastpos - queryoffset)] == false) {
	    /* Skip */
	    debug(printf(" => querypos_rc %d not valid\n",querypos_rc));
	    minus_qpos[mod] += index1interval;
	  
	  } else {
	    assert(stage1->minus_retrievedp[querypos_rc] == true);
#ifdef LARGE_GENOMES
	    positions_high = stage1->minus_positions_high[querypos_rc];
#endif
	    positions = stage1->minus_positions[querypos_rc];
	    npositions = stage1->minus_npositions[querypos_rc];
	    debug(printf(" => saw %d positions\n",npositions));
	    
	    if ((elt = Elt_read_queryfwd(
#ifdef LARGE_GENOMES
					 positions_high,
#endif
					 positions,npositions,/*diagterm*/querylength - queryoffset,
					 querylength,/*querystart*/queryoffset,
					 query_compress_rev,/*plusp*/false,genestrand,
					 univdiagpool)) == NULL) {
	      minus_qpos[mod] += index1interval;
	    
	    } else {
	      if (elt->nmatches > best_minus_nmatches && elt->n_all_univdiags <= maxpaths_search &&
		  elt->nunivdiags > 0) { /* Could be 0 if there are too many univdiags */
		best_minus_nmatches = elt->nmatches;
	      }

	      debug(printf("Pushing elt %p onto minus set\n",elt));
	      *minus_set = Listpool_push(*minus_set,listpool,(void *) elt
					 listpool_trace(__FILE__,__LINE__));
	      minus_qpos[mod] += elt->nmatches;
	    }
	  }

	  if (minus_qpos[mod] > max_minus_qpos) {
	    max_minus_qpos = minus_qpos[mod];
	  }
	}

	total_npositions_minus += mod_npositions;
	niter_minus++;
      }
    }

#ifdef DEBUG
    printf("\nStatus:\n");
    for (mod = 0; mod < index1interval; mod++) {
      printf("mod %d, plus_qpos %d\n",mod,plus_qpos[mod]);
    }
    printf("max_plus_qpos %d\n",max_plus_qpos);
    printf("\n");
    for (mod = 0; mod < index1interval; mod++) {
      printf("mod %d, minus_qpos %d\n",mod,minus_qpos[mod]);
    }
    printf("max_minus_qpos %d\n",max_minus_qpos);
    printf("\n");
#endif

#ifdef CONSERVATIVE
    /* Skip the presumed mismatch */
    max_plus_qpos += 1;
    max_minus_qpos += 1;
#endif

    plus_qpos[0] = max_plus_qpos;
    plus_qpos[1] = max_plus_qpos + 1;
    plus_qpos[2] = max_plus_qpos + 2;

    minus_qpos[0] = max_minus_qpos;
    minus_qpos[1] = max_minus_qpos + 1;
    minus_qpos[2] = max_minus_qpos + 2;
  }


  /* II.  Fill best_plus_elts and best_minus_elts */
  *best_plus_elts = *best_minus_elts = (List_T) NULL;
  for (p = *plus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_plus_nmatches - 3 && elt->n_all_univdiags <= maxpaths_search) {
      if (elt->nunivdiags > 0) {
	/* Could be 0 if there are too many univdiagonals */
	debug(printf("Pushing elt %p onto best_plus_elts\n",elt));
	*best_plus_elts = Listpool_push(*best_plus_elts,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (p = *minus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_minus_nmatches - 3 && elt->n_all_univdiags <= maxpaths_search) {
      if (elt->nunivdiags > 0) {
	/* Could be 0 if there are too many univdiags */
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
  Elt_dump_set(*plus_set);
  printf("\n");
  printf("queryfwd minus set:\n");
  Elt_dump_set(*minus_set);
  printf("\n");
#endif

#ifdef DEBUG
  printf("QUERYFWD PLUS HAS BEST ELTS:\n");
  Elt_dump_set(*best_plus_elts);
  printf("\n");

  printf("QUERYFWD MINUS HAS BEST ELTS:\n");
  Elt_dump_set(*best_minus_elts);
  printf("\n");
#endif

#ifdef PICK_SIDE
  if (max_minus_qpos >= query_lastpos && max_plus_qpos >= query_lastpos) {
    debug(printf("QUERYFWD: both sides won (max_minus_qpos %d, max_plus_qpos %d), so use both sides\n",
		 max_minus_qpos,max_plus_qpos));
    debug(printf("QUERYFWD PLUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYFWD MINUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_minus_elts));
    debug(printf("\n"));

  } else if (max_minus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD PLUS: minus side won (max_minus_qpos %d), so skip plus side\n",
		 max_minus_qpos));
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;

  } else if (max_plus_qpos >= query_lastpos) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYFWD MINUS: plus side won (max_plus_qpos %d), so skip minus side\n",
		 max_plus_qpos));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_minus_elts = (List_T) NULL;
    
  } else {
#if 0
    /* This code also has a high cost in speed */
    debug(printf("QUERYFWD: both sides lost (max_minus_qpos %d, max_plus_qpos %d), so use both sides\n",
		 max_minus_qpos,max_plus_qpos));
    debug(printf("QUERYFWD PLUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYFWD MINUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_minus_elts));
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
#endif

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
		       int max_iter, int genestrand, Listpool_T listpool,
		       Univdiagpool_T univdiagpool) {
  T elt;
  List_T p;
  int best_plus_nmatches, best_minus_nmatches;
  int min_plus_qpos, min_minus_qpos, plus_qpos[MAX_INDEX1INTERVAL], minus_qpos[MAX_INDEX1INTERVAL];

  int mod;
  int niter_plus, niter_minus;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
#endif
  UINT4 *positions;
  int total_npositions_plus = 0, total_npositions_minus = 0, npositions;

  int query_lastpos = querylength - index1part;
  int queryoffset, querypos_rc;

  debug(printf("\nStarting get_elt_sets_queryrev with querylength %d and max_iter %d, genestrand %d\n",
	       querylength,max_iter,genestrand));

#if 0
  /* Allow calls to Extension_search_queryrev in addition to other methods */
  *paths_gplus = *paths_gminus = (List_T) NULL;
#endif


  /* I.  Race from plus and minus end to start.  Compute best_plus_nmatches and best_minus_nmatches */
  best_plus_nmatches = best_minus_nmatches = 0;
  min_plus_qpos = min_minus_qpos = query_lastpos;
  plus_qpos[0] = minus_qpos[0] = query_lastpos;
  plus_qpos[1] = minus_qpos[1] = query_lastpos-1;
  plus_qpos[2] = minus_qpos[2] = query_lastpos-2;
  niter_plus = niter_minus = 0;

  while ((total_npositions_plus <= MAX_TOTAL_NPOSITIONS ||
	  total_npositions_minus <= MAX_TOTAL_NPOSITIONS) &&
	 niter_plus <= max_iter && niter_minus <= max_iter &&
	 min_plus_qpos > 0 && min_minus_qpos > 0) {

    if (total_npositions_plus <= MAX_TOTAL_NPOSITIONS) {
      for (mod = 0; mod < index1interval; mod++) {
	debug(printf("queryrev plus mod %d",mod));
	if ((queryoffset = plus_qpos[mod]) < 0) {
	  /* Skip */
	  debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	  plus_qpos[mod] -= index1interval;
	  
	} else if (stage1->validp[queryoffset] == false) {
	  /* Skip */
	  debug(printf(" => queryoffset %d not valid\n",queryoffset));
	  plus_qpos[mod] -= index1interval;
	  
	} else {
	  if (stage1->plus_retrievedp[queryoffset] == true) {
#ifdef LARGE_GENOMES
	    positions_high = stage1->plus_positions_high[queryoffset];
#endif
	    positions = stage1->plus_positions[queryoffset];
	    npositions = stage1->plus_npositions[queryoffset];
	    debug(printf(" => saw %d positions\n",npositions));
	    
	  } else {
	    assert(stage1->plus_positions[queryoffset] == NULL);
#ifdef LARGE_GENOMES
	    npositions = stage1->plus_npositions[queryoffset] =
	      Indexdb_largeptr(&stage1->plus_positions_high[queryoffset],&stage1->plus_positions[queryoffset],
			       indexdb,stage1->forward_oligos[queryoffset]);
	    positions_high = stage1->plus_positions_high[queryoffset];
#else
	    npositions = stage1->plus_npositions[queryoffset] =
	      Indexdb_ptr(&stage1->plus_positions[queryoffset],indexdb,
			  stage1->forward_oligos[queryoffset]);
#endif
	    positions = stage1->plus_positions[queryoffset];
	    stage1->plus_retrievedp[queryoffset] = true;
	    debug(printf(" => retrieved %d positions\n",npositions));
	  }
	  
	  
	  /* debug9(Stage1_dump(stage1,querylength)); */
	  if (npositions > MAX_NPOSITIONS) {
	    /* Skip */
	    plus_qpos[mod] -= index1interval;
	    
	  } else if ((elt = Elt_read_queryrev(
#ifdef LARGE_GENOMES
					      positions_high,
#endif
					      positions,npositions,/*diagterm*/querylength - queryoffset,
					      /*queryend*/queryoffset + index1part,querylength,
					      query_compress_fwd,/*plusp*/true,genestrand,
					      univdiagpool)) == NULL) {
	    total_npositions_plus += npositions;
	    plus_qpos[mod] -= index1interval;
	    
	  } else {
	    if (elt->nmatches > best_plus_nmatches && elt->n_all_univdiags <= maxpaths_search) {
	      if (elt->nunivdiags > 0) {
		/* Could be 0 if there are too many univdiags */
		best_plus_nmatches = elt->nmatches;
	      }
	    }
	    debug(printf("Pushing elt %p onto plus_set\n",elt));
	    *plus_set = Listpool_push(*plus_set,listpool,(void *) elt
				      listpool_trace(__FILE__,__LINE__));
	    plus_qpos[mod] -= elt->nmatches;
	    total_npositions_plus += npositions;
	    niter_plus++;
	  }
	}
	
	if (plus_qpos[mod] < min_plus_qpos) {
	  min_plus_qpos = plus_qpos[mod];
	}
      }
    }


    if (total_npositions_minus <= MAX_TOTAL_NPOSITIONS) {
      for (mod = 0; mod < index1interval; mod++) {
	debug(printf("queryrev minus mod %d",mod));
	if ((queryoffset = minus_qpos[mod]) < 0) {
	  /* Skip */
	  debug(printf(" => queryoffset %d greater than query_lastpos\n",queryoffset));
	  minus_qpos[mod] -= index1interval;
	  
	} else if (stage1->validp[(querypos_rc = query_lastpos - queryoffset)] == false) {
	  /* Skip */
	  debug(printf(" => querypos_rc %d not valid\n",querypos_rc));
	  minus_qpos[mod] -= index1interval;
	  
	} else {
	  if (stage1->minus_retrievedp[querypos_rc] == true) {
#ifdef LARGE_GENOMES
	    positions_high = stage1->minus_positions_high[querypos_rc];
#endif
	    positions = stage1->minus_positions[querypos_rc];
	    npositions = stage1->minus_npositions[querypos_rc];
	    debug(printf(" => saw %d positions\n",npositions));
	    
	  } else {
	    assert(stage1->minus_positions[querypos_rc] == NULL);
#ifdef LARGE_GENOMES
	    npositions = stage1->minus_npositions[querypos_rc] =
	      Indexdb_largeptr(&stage1->minus_positions_high[querypos_rc],&stage1->minus_positions[querypos_rc],
			       indexdb_nonstd,stage1->revcomp_oligos[querypos_rc]);
	    positions_high = stage1->minus_positions_high[querypos_rc];
#else
	    npositions = stage1->minus_npositions[querypos_rc] =
	      Indexdb_ptr(&stage1->minus_positions[querypos_rc],indexdb_nonstd,
			  stage1->revcomp_oligos[querypos_rc]);
#endif
	    positions = stage1->minus_positions[querypos_rc];
	    stage1->minus_retrievedp[querypos_rc] = true;
	    debug(printf(" => retrieved %d positions\n",npositions));
	  }
	  
	  
	  /* debug9(Stage1_dump(stage1,querylength)); */
	  if (npositions > MAX_NPOSITIONS) {
	    /* Skip */
	    minus_qpos[mod] -= index1interval;
	    
	  } else if ((elt = Elt_read_queryrev(
#ifdef LARGE_GENOMES
					      positions_high,
#endif
					      positions,npositions,/*diagterm*/querylength - queryoffset,
					      /*queryend*/queryoffset + index1part,querylength,
					      query_compress_rev,/*plusp*/false,genestrand,
					      univdiagpool)) == NULL) {
	    total_npositions_minus += npositions;
	    minus_qpos[mod] -= index1interval;
	    
	  } else {
	    if (elt->nmatches > best_minus_nmatches && elt->n_all_univdiags <= maxpaths_search) {
	      if (elt->nunivdiags > 0) {
		/* Could be 0 if there are too many univdiags */
		best_minus_nmatches = elt->nmatches;
	      }
	    }
	    debug(printf("Pushing elt %p onto minus_set\n",elt));
	    *minus_set = Listpool_push(*minus_set,listpool,(void *) elt
				       listpool_trace(__FILE__,__LINE__));
	    minus_qpos[mod] -= elt->nmatches;
	    total_npositions_minus += npositions;
	    niter_minus++;
	  }
	}
	
	if (minus_qpos[mod] < min_minus_qpos) {
	  min_minus_qpos = minus_qpos[mod];
	}
      }
    }

#ifdef DEBUG
    printf("\nStatus:\n");
    for (mod = 0; mod < index1interval; mod++) {
      printf("mod %d, plus_qpos %d\n",mod,plus_qpos[mod]);
      printf("mod %d, minus_qpos %d\n",mod,minus_qpos[mod]);
    }
    printf("min_plus_qpos %d\n",min_plus_qpos);
    printf("min_minus_qpos %d\n",min_minus_qpos);
    printf("\n");
#endif

#ifdef CONSERVATIVE
    /* Skip the presumed mismatch */
    min_plus_qpos -= 1;
    min_minus_qpos -= 1;
#endif

    plus_qpos[0] = min_plus_qpos;
    plus_qpos[1] = min_plus_qpos - 1;
    plus_qpos[2] = min_plus_qpos - 2;

    minus_qpos[0] = min_minus_qpos;
    minus_qpos[1] = min_minus_qpos - 1;
    minus_qpos[2] = min_minus_qpos - 2;
  }


  /* II.  Fill best_plus_elts and best_minus_elts */
  *best_plus_elts = *best_minus_elts = (List_T) NULL;
  for (p = *plus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_plus_nmatches - 3 && elt->n_all_univdiags <= maxpaths_search) {
      if (elt->nunivdiags > 0) {
	/* Could be 0 if there are too many univdiags */
	debug(printf("Pushing elt %p onto best_plus_elts\n",elt));
	*best_plus_elts = Listpool_push(*best_plus_elts,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
      }
    }
  }

  for (p = *minus_set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    if (elt->nmatches > best_minus_nmatches - 3 && elt->n_all_univdiags <= maxpaths_search) {
      if (elt->nunivdiags > 0) {
	/* Could be 0 if there are too many univdiags */
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
  Elt_dump_set(*plus_set);
  printf("\n");
  printf("queryrev minus set:\n");
  Elt_dump_set(*minus_set);
  printf("\n");
#endif

#ifdef DEBUG
  printf("QUERYREV PLUS HAS BEST ELTS:\n");
  Elt_dump_set(*best_plus_elts);
  printf("\n");

  printf("QUERYREV MINUS HAS BEST ELTS:\n");
  Elt_dump_set(*best_minus_elts);
  printf("\n");
#endif

#ifdef PICK_SIDE
  if (min_minus_qpos <= 0 && min_plus_qpos <= 0) {
    debug(printf("QUERYREV: both sides won (min_minus_qpos %d, min_plus_qpos %d), so use both sides\n",
		 min_minus_qpos,min_plus_qpos));
    debug(printf("QUERYREV PLUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYREV MINUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_minus_elts));
    debug(printf("\n"));

  } else if (min_minus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV PLUS: minus side won (min_minus_qpos %d), so skip plus side\n",
		 min_minus_qpos));
    Listpool_free_list(&(*best_plus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_plus_elts = (List_T) NULL;

  } else if (min_plus_qpos <= 0) {
    /* This branch is critical for high speed (about 4x faster),
       although we could miss some hits that need to be solved by
       another method */
    debug(printf("QUERYREV MINUS: plus side won (min_plus_qpos %d), so skip minus side\n",
		 min_plus_qpos));
    Listpool_free_list(&(*best_minus_elts),listpool
		       listpool_trace(__FILE__,__LINE__));
    *best_minus_elts = (List_T) NULL;

  } else {
#if 0
    /* This code also has a high cost in speed */
    debug(printf("QUERYREV: both sides lost (min_minus_qpos %d, min_plus_qpos %d), so use both sides\n",
		 min_minus_qpos,min_plus_qpos));
    debug(printf("QUERYREV PLUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_plus_elts));
    debug(printf("\n"));

    debug(printf("QUERYREV MINUS HAS BEST ELTS:\n"));
    debug(Elt_dump_set(*best_minus_elts));
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
#endif

  return;
}


static void
process_seed (Univcoordlist_T *univdiagonals_list, List_T *auxinfo_list,

	      Univcoord_T middle_univdiagonal, int middle_qstart, int middle_qend, int middle_nmismatches,
	      List_T queryfwd_set, List_T queryrev_set, T queryfwd_best_elt, T queryrev_best_elt,
	      int querylength, Auxinfopool_T auxinfopool, Univdiagpool_T univdiagpool,
	      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool) {

  Auxinfo_T auxinfo;

  List_T left_univdiags = NULL, right_univdiags = NULL, p;
  Univdiag_T *left_univdiag_array, *right_univdiag_array;
  Univcoord_T univdiagonal, start_low, start_high, end_low, end_high;
  int qstart, qend;
  Elt_T elt;
  int nleft, nright, i, j, k;

  /* Chrnum_T chrnum; */
  Univcoord_T chroffset, chrhigh;


  /* Previously computed middle_qstart and middle_qend are the maximum
     bounds for the entire elt.  But now, in order to set
     path->nmismatches to be 0, we are computing the minimum bounds */
  /* When we were computing maximum bounds for the entire elt, for an individual univdiagonal,
     they might go past the chromosome boundaries */

  /* left = middle_univdiagonal - querylength; */
  debug(printf("Entering process_seed with middle_univdiagonal %u\n",middle_univdiagonal));

#ifdef TRIM_AT_GENOME_BOUNDS
  middle_qstart = (middle_univdiagonal + middle_qstart >= 0 + (Univcoord_T) querylength) ? middle_qstart : (int) (querylength - middle_univdiagonal);
  middle_qend = (middle_univdiagonal + middle_qend <= genomelength + (Univcoord_T) querylength) ? middle_qend : (int) (genomelength - middle_univdiagonal + querylength);
#endif

#if 0
  /* No longer works now that we have sets of best elts */
  /* Need to update chrnum, chroffset, chrhigh, chrlength through multiple calls to process_seed */
  *chrhigh = Univ_IIT_update_chrnum(&(*chrnum),&(*chroffset),*chrhigh,&(*chrlength),chromosome_iit,
				    /*low*/left + middle_qstart,/*high*/left + middle_qend,circular_typeint);
#endif
#ifdef USE_CHROMOSOME_IIT
  /* chrnum = */ Univ_IIT_get_chrnum(&chroffset,&chrhigh,&chrlength,chromosome_iit,
			       /*low*/left + middle_qstart,/*high*/left + middle_qend,circular_typeint);
#else
  /* chrnum = */ EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			     middle_univdiagonal - querylength + middle_qstart,
			     middle_univdiagonal - querylength + middle_qend);
#endif

  
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  middle_qstart = (middle_univdiagonal + middle_qstart >= chroffset + (Univcoord_T) querylength) ? middle_qstart : (int) (chroffset - left);
  middle_qend = (middle_univdiagonal + middle_qend <= chrhigh + (Univcoord_T) querylength) ? middle_qend : (int) (chrhigh - left);
#endif
  debug(printf("PROCESS SEED at %u, qstart %d, qend %d\n",
	       middle_univdiagonal,middle_qstart,middle_qend));


  start_low = subtract_bounded(middle_univdiagonal,/*minusterm*/(Univcoord_T) positive_gap_distance);
  start_high = add_bounded(middle_univdiagonal,/*plusterm*/(Univcoord_T) negative_gap_distance);

  end_low = subtract_bounded(middle_univdiagonal,/*minusterm*/(Univcoord_T) negative_gap_distance);
  end_high = add_bounded(middle_univdiagonal,/*plusterm*/(Univcoord_T) positive_gap_distance);

  debug(printf("Computing start %u..%u and end %u..%u\n",start_low,start_high,end_low,end_high));


  for (p = queryfwd_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_univdiags(elt,start_low,start_high);
	for (k = 0; k < elt->nunivdiags; k++) {
	  left_univdiags = Univdiagpool_push_existing(left_univdiags,univdiagpool,elt->univdiags[k]
						      univdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing univdiag %u %d..%d onto left of query fwd\n",
		       elt->univdiags[k]->univdiagonal,elt->univdiags[k]->qstart,elt->univdiags[k]->qend));
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_univdiags(elt,end_low,end_high);
	for (k = 0; k < elt->nunivdiags; k++) {
	  right_univdiags = Univdiagpool_push_existing(right_univdiags,univdiagpool,elt->univdiags[k]
						       univdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing univdiag %u %d..%d onto right of query fwd\n",
		       elt->univdiags[k]->univdiagonal,elt->univdiags[k]->qstart,elt->univdiags[k]->qend));
	}
      }
    }
  }
  
  for (p = queryrev_set; p != NULL; p = List_next(p)) {
    elt = (Elt_T) List_head(p);
    if (elt != queryfwd_best_elt && elt != queryrev_best_elt) {
      if (elt_startp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_univdiags(elt,start_low,start_high);
	for (k = 0; k < elt->nunivdiags; k++) {
	  left_univdiags = Univdiagpool_push_existing(left_univdiags,univdiagpool,elt->univdiags[k]
						      univdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing univdiag %u %d..%d onto left of query rev\n",
		       elt->univdiags[k]->univdiagonal,elt->univdiags[k]->qstart,elt->univdiags[k]->qend));
	}
      }
      
      if (elt_endp(elt,middle_qstart,middle_qend) == true) {
	Elt_filter_univdiags(elt,end_low,end_high);
	for (k = 0; k < elt->nunivdiags; k++) {
	  right_univdiags = Univdiagpool_push_existing(right_univdiags,univdiagpool,elt->univdiags[k]
						       univdiagpool_trace(__FILE__,__LINE__));
	  debug(printf("Pushing univdiag %u %d..%d onto right of query rev\n",
		       elt->univdiags[k]->univdiagonal,elt->univdiags[k]->qstart,elt->univdiags[k]->qend));
	}
      }
    }
  }

  /* At this point, left_univdiags and right_univdiags are a subset of those in Elt_T objects */
  /* Below, the contents of left_univdiags and right_univdiags changes */

  if (left_univdiags != NULL) {
    /* Sort the left univdiagonals and unique them */
    left_univdiag_array = (Univdiag_T *) List_to_array_n(&nleft,left_univdiags);
    qsort(left_univdiag_array,nleft,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
    Univdiagpool_free_list(&left_univdiags,univdiagpool
			   univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */

    /* Combine left univdiagonals to expand their qstart and qend range */
    left_univdiags = (List_T) NULL;
    i = 0;
    while (i < nleft) {
      univdiagonal = left_univdiag_array[i]->univdiagonal;
      qstart = left_univdiag_array[i]->qstart;
      qend = left_univdiag_array[i]->qend;
      
      j = i+1;
      while (j < nleft && left_univdiag_array[j]->univdiagonal == univdiagonal) {
	debug(printf("At left diagonal %u, combining %d..%d with %d..%d\n",
		     univdiagonal,qstart,qend,left_univdiag_array[j]->qstart,left_univdiag_array[j]->qend));
	if (left_univdiag_array[j]->qstart < qstart) {
	  qstart = left_univdiag_array[j]->qstart;
	}
	if (left_univdiag_array[j]->qend > qend) {
	  qend = left_univdiag_array[j]->qend;
	}
	j++;
      }
      
      if (qstart == left_univdiag_array[i]->qstart && qend == left_univdiag_array[i]->qend) {
#if 0
	/* This shortcut confuses existing Univdiag_T objects with newly created ones */
	left_univdiags = Univdiagpool_push_existing(left_univdiags,univdiagpool,left_univdiag_array[i]
						    univdiagpool_trace(__FILE__,__LINE__));
#else
	left_univdiags = Univdiagpool_push(left_univdiags,univdiagpool,qstart,qend,
					   left_univdiag_array[i]->nmismatches,univdiagonal
					   univdiagpool_trace(__FILE__,__LINE__));
#endif
      } else {
	left_univdiags = Univdiagpool_push(left_univdiags,univdiagpool,qstart,qend,/*nmismatches*/-1,univdiagonal
					   univdiagpool_trace(__FILE__,__LINE__));
      }	

      i = j;
    }
    FREE(left_univdiag_array);
  }

  if (right_univdiags != NULL) {
    /* Sort the right univdiagonals and unique them */
    right_univdiag_array = (Univdiag_T *) List_to_array_n(&nright,right_univdiags);
    qsort(right_univdiag_array,nright,sizeof(Univdiag_T),Univdiag_diagonal_cmp);
    Univdiagpool_free_list(&right_univdiags,univdiagpool
			   univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */
    
    debug(printf("right univdiagonals is not NULL => %d\n",nright));

    /* Combine right univdiagonals to expand their qstart and qend range */
    right_univdiags = (List_T) NULL;
    i = 0;
    while (i < nright) {
      univdiagonal = right_univdiag_array[i]->univdiagonal;
      qstart = right_univdiag_array[i]->qstart;
      qend = right_univdiag_array[i]->qend;
      
      j = i+1;
      while (j < nright && right_univdiag_array[j]->univdiagonal == univdiagonal) {
	debug(printf("At right diagonal %u, combining %d..%d with %d..%d\n",
		     univdiagonal,qstart,qend,right_univdiag_array[j]->qstart,right_univdiag_array[j]->qend));
	if (right_univdiag_array[j]->qstart < qstart) {
	  qstart = right_univdiag_array[j]->qstart;
	}
	if (right_univdiag_array[j]->qend > qend) {
	  qend = right_univdiag_array[j]->qend;
	}
	j++;
      }
      
      if (qstart == right_univdiag_array[i]->qstart && qend == right_univdiag_array[i]->qend) {
#if 0
	/* This shortcut confuses existing Univdiag_T objects with newly created ones */
	right_univdiags = Univdiagpool_push_existing(right_univdiags,univdiagpool,right_univdiag_array[i]
						     univdiagpool_trace(__FILE__,__LINE__));
#else
	right_univdiags = Univdiagpool_push(right_univdiags,univdiagpool,qstart,qend,
					    right_univdiag_array[i]->nmismatches,univdiagonal
					    univdiagpool_trace(__FILE__,__LINE__));
#endif	
      } else {
	right_univdiags = Univdiagpool_push(right_univdiags,univdiagpool,qstart,qend,/*nmismatches*/-1,univdiagonal
					    univdiagpool_trace(__FILE__,__LINE__));
      }	

      i = j;
    }
    FREE(right_univdiag_array);
  }

  /* At this point, left_univdiags and right_univdiags are new Univdiag_T objects, not in any Elt_T object */


#if 0

#ifdef DEBUG
  Univdiag_T univdiag;
  printf("Calling Path_solve_from_diagonals with middle diagonal %u [%u], %d..%d, and %d left and %d right diagonals.  chroffset %u, chrhigh %u\n",
	 middle_univdiagonal,middle_univdiagonal - chroffset,middle_qstart,middle_qend,
	 List_length(left_univdiags),List_length(right_univdiags),chroffset,chrhigh);
  for (p = left_univdiags; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    printf("Left %u, %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
  }
  for (p = right_univdiags; p != NULL; p = List_next(p)) {
    univdiag = (Univdiag_T) List_head(p);
    printf("Right %u, %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
  }
#endif

  Path_solve_from_diagonals(&(*found_score),
			    &(*unextended_sense_paths),&(*unextended_antisense_paths),
			    &(*sense_paths),&(*antisense_paths),
			    middle_univdiagonal,middle_qstart,middle_qend,middle_nmismatches,
			    /*qend_univdiags*/right_univdiags,/*qstart_univdiags*/left_univdiags,
			    queryseq,queryptr,querylength,
			    mismatch_positions_alloc,novel_diagonals_alloc,localdb_alloc,
			    stage1,knownsplicing,knownindels,
			    query_compress,query_compress_fwd,query_compress_rev,
			    chrnum,chroffset,chrhigh,plusp,genestrand,
			    nmismatches_allowed,paired_end_p,first_read_p,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,
			    pathpool,transcriptpool,vectorpool,hitlistpool,spliceendsgen,method,
			    find_splices_p);

  debug(printf("Path_solve_from_diagonals now at %d + %d unextended, %d sense paths and %d antisense paths\n\n",
	       List_length(*unextended_sense_paths),List_length(*unextended_sense_paths),
	       List_length(*sense_paths),List_length(*antisense_paths)));
  Univdiagpool_gc(&right_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */
  Univdiagpool_gc(&left_univdiags,univdiagpool
		  univdiagpool_trace(__FILE__,__LINE__)); /* allocated by Univdiagpool_push */

#else
  
  /* New algorithm: Instead of calling Path_solve_from_diagonals here,
     try return univdiagonals.  But also need indices to
     qend_univdiags and qstart_univdiags */
  debug(printf("Pushing middle univdiagonal %u\n",middle_univdiagonal));
  *univdiagonals_list = Univcoordlistpool_push(*univdiagonals_list,univcoordlistpool,middle_univdiagonal
					       univcoordlistpool_trace(__FILE__,__LINE__));
  auxinfo = Auxinfo_new_univdiags(/*method*/EXT,middle_qstart,middle_qend,middle_nmismatches,
				  right_univdiags,left_univdiags,auxinfopool);
  *auxinfo_list = Listpool_push(*auxinfo_list,listpool,(void *) auxinfo
				listpool_trace(__FILE__,__LINE__));
#endif

  return;
}


#define min(a,b) (a < b) ? a : b
#define max(a,b) (a > b) ? a : b


static void
extend_seeds_union (Univcoordlist_T *univdiagonals_list, List_T *auxinfo_list,

		    List_T queryfwd_best_elts, List_T queryrev_best_elts,
		    List_T queryfwd_set, List_T queryrev_set, int querylength,

		    Auxinfopool_T auxinfopool, Univdiagpool_T univdiagpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool) {
  
  List_T univdiags = NULL, queryfwd_elt_list = NULL, queryrev_elt_list = NULL;
  Univcoordlist_T univdiagonals = NULL; /* To determine order and avoid duplicates */
  Univcoord_T *array;

  int *order;
  Univdiag_T *univdiag_array, queryfwd_univdiag, queryrev_univdiag, univdiag;

  int qstart, qend;
  int queryfwd_nseeds, queryrev_nseeds, n, i, j;
  List_T p, q;
  T *queryfwd_elt_array, *queryrev_elt_array, queryfwd_elt, queryrev_elt;

  
  debug(printf("Entering extend_seeds_union with %d queryfwd_best_elts and %d queryrev_best_elts\n",
	       List_length(queryfwd_best_elts),List_length(queryrev_best_elts)));

  /* Create a union of seeds */
  for (p = queryfwd_best_elts; p != NULL; p = List_next(p)) {
    queryfwd_elt = (T) List_head(p);
    queryfwd_nseeds = queryfwd_elt->n_all_univdiags;

    for (q = queryrev_best_elts; q != NULL; q = List_next(q)) {
      queryrev_elt = (T) List_head(q);
      queryrev_nseeds = queryrev_elt->n_all_univdiags;
    
#ifdef DEBUG
      printf("queryfwd_elt: ");
      Elt_dump(queryfwd_elt);
      printf("\n");
      printf("queryrev_elt: ");
      Elt_dump(queryrev_elt);
      printf("\n\n");
#endif

      i = j = 0;
      while (i < queryfwd_nseeds && j < queryrev_nseeds) {
	queryfwd_univdiag = queryfwd_elt->all_univdiags[i];
	queryrev_univdiag = queryrev_elt->all_univdiags[j];
	
	if (queryfwd_univdiag->univdiagonal == queryrev_univdiag->univdiagonal) {
	  /* Combine the queryfwd and queryrev univdiags */
	  debug(printf("At middle diagonal %u, combining %d..%d with %d..%d\n",
		       queryfwd_univdiag->univdiagonal,
		       queryfwd_univdiag->qstart,queryfwd_univdiag->qend,
		       queryrev_univdiag->qstart,queryrev_univdiag->qend));
	  qstart = min(queryfwd_univdiag->qstart,queryrev_univdiag->qstart);
	  qend = max(queryfwd_univdiag->qend,queryrev_univdiag->qend);
	  univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryfwd_univdiag->univdiagonal
						 univcoordlistpool_trace(__FILE__,__LINE__));
	  univdiags = Univdiagpool_push(univdiags,univdiagpool,qstart,qend,
					/*nmismatches*/-1,queryfwd_univdiag->univdiagonal
					univdiagpool_trace(__FILE__,__LINE__));
	  
	  i++; j++;
	  
	} else if (queryfwd_univdiag->univdiagonal < queryrev_univdiag->univdiagonal) {
	  univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryfwd_univdiag->univdiagonal
						 univcoordlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Univdiag_T objects with those in Elt_T objects */
	  univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,queryfwd_univdiag
						 univdiagpool_trace(__FILE__,__LINE__));
#else
	  univdiags = Univdiagpool_push(univdiags,univdiagpool,queryfwd_univdiag->qstart,queryfwd_univdiag->qend,
					queryfwd_univdiag->nmismatches,queryfwd_univdiag->univdiagonal
					univdiagpool_trace(__FILE__,__LINE__));
#endif
	  i++;
	  
	} else {
	  univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryrev_univdiag->univdiagonal
						 univcoordlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Univdiag_T objects with those in Elt_T objects */
	  univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,queryrev_univdiag
						 univdiagpool_trace(__FILE__,__LINE__));
#else
	  univdiags = Univdiagpool_push(univdiags,univdiagpool,queryrev_univdiag->qstart,queryrev_univdiag->qend,
					queryrev_univdiag->nmismatches,queryrev_univdiag->univdiagonal
					univdiagpool_trace(__FILE__,__LINE__));
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
	  queryfwd_univdiag = queryfwd_elt->all_univdiags[i++];
	  
	  /* Allow univdiagonals from queryrev_elt */
	  univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryfwd_univdiag->univdiagonal
						 univcoordlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Univdiag_T objects with those in Elt_T objects */
	  univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,queryfwd_univdiag
						 univdiagpool_trace(__FILE__,__LINE__));
#else
	  univdiags = Univdiagpool_push(univdiags,univdiagpool,queryfwd_univdiag->qstart,queryfwd_univdiag->qend,
					queryfwd_univdiag->nmismatches,queryfwd_univdiag->univdiagonal
					univdiagpool_trace(__FILE__,__LINE__));
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
	  queryrev_univdiag = queryrev_elt->all_univdiags[j++];
	  
	  /* Allow univdiagonals from queryfwd_elt */
	  univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryrev_univdiag->univdiagonal
						 univcoordlistpool_trace(__FILE__,__LINE__));
#if 0
	  /* This shortcut confuses newly created Univdiag_T objects with those in Elt_T objects */
	  univdiags = Univdiagpool_push_existing(univdiags,univdiagpool,queryrev_univdiag
						 univdiagpool_trace(__FILE__,__LINE__));
#else
	  univdiags = Univdiagpool_push(univdiags,univdiagpool,queryrev_univdiag->qstart,queryrev_univdiag->qend,
					queryrev_univdiag->nmismatches,queryrev_univdiag->univdiagonal
					univdiagpool_trace(__FILE__,__LINE__));
#endif
	  queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) NULL
					    listpool_trace(__FILE__,__LINE__));
	  queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) queryrev_elt
					    listpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }
      
  /* univdiags contains newly created Univdiag_T objects, not in Elt_T objects */


  /* Want to sort three arrays simultaneously: univdiags,
     queryfwd_elts, and queryrev_elts.  Use order on the univdiagonals
     to do this. */
  n = Univcoordlist_length(univdiagonals);
  array = Univcoordlist_to_array(univdiagonals,/*end*/0);
  univdiag_array = (Univdiag_T *) List_to_array(univdiags,NULL);
  queryfwd_elt_array = (T *) List_to_array(queryfwd_elt_list,NULL);
  queryrev_elt_array = (T *) List_to_array(queryrev_elt_list,NULL);

  Listpool_free_list(&queryrev_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Listpool_free_list(&queryfwd_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Univdiagpool_free_list(&univdiags,univdiagpool
			 univdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Univdiagpool_push */
  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_push */

#ifdef LARGE_GENOMES
  order = Sedgesort_order_uint8(array,n);
#else
  order = Sedgesort_order_uint4(array,n);
#endif
  
  /* univdiag_array contains newly created Univdiag_T objects, not in Elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    univdiag = univdiag_array[order[i]];
    process_seed(&(*univdiagonals_list),&(*auxinfo_list),

		 univdiag->univdiagonal,univdiag->qstart,univdiag->qend,univdiag->nmismatches,
		 queryfwd_set,queryrev_set,
		 /*queryfwd_best_elt*/queryfwd_elt_array[order[i]],
		 /*queryrev_best_elt*/queryrev_elt_array[order[i]],
		 querylength,auxinfopool,univdiagpool,univcoordlistpool,listpool);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Univdiagpool_free_univdiag(&(univdiag_array[i]),univdiagpool
			       univdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(univdiag_array);

  FREE(order);
  FREE(queryrev_elt_array);
  FREE(queryfwd_elt_array);
  FREE(array);

  return;
}


static void
extend_seeds_queryfwd (Univcoordlist_T *univdiagonals_list, List_T *auxinfo_list,
		       List_T queryfwd_best_elts, List_T queryfwd_set, int querylength,

		       Auxinfopool_T auxinfopool, Univdiagpool_T univdiagpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool) {
  
  List_T univdiags = NULL, queryfwd_elt_list = NULL;
  Univcoordlist_T univdiagonals = NULL; /* To determine order and avoid duplicates */
  Univcoord_T *array;

  int *order;
  Univdiag_T *univdiag_array, queryfwd_univdiag, univdiag;

  int queryfwd_nseeds, n, i, j;
  List_T p;
  T *queryfwd_elt_array, queryfwd_elt;

  
  debug(printf("Entering extend_seeds_queryfwd with %d queryfwd_best_elts\n",
	       List_length(queryfwd_best_elts)));

  /* Handle only queryfwd_best_elts */
  for (p = queryfwd_best_elts; p != NULL; p = List_next(p)) {
    queryfwd_elt = (T) List_head(p);
    queryfwd_nseeds = queryfwd_elt->n_all_univdiags;

    i = 0;
    while (i < queryfwd_nseeds) {
      queryfwd_univdiag = queryfwd_elt->all_univdiags[i++];
	    
      /* Allow univdiagonals from queryrev_elt */
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryfwd_univdiag->univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));
      univdiags = Univdiagpool_push(univdiags,univdiagpool,queryfwd_univdiag->qstart,queryfwd_univdiag->qend,
				    queryfwd_univdiag->nmismatches,queryfwd_univdiag->univdiagonal
				    univdiagpool_trace(__FILE__,__LINE__));
	    
      queryfwd_elt_list = Listpool_push(queryfwd_elt_list,listpool,(void *) queryfwd_elt
					listpool_trace(__FILE__,__LINE__));
    }
  }
      
  /* univdiags contains newly created Univdiag_T objects, not in Elt_T objects */


  /* Want to sort two arrays simultaneously: univdiags and
    queryfwd_elts.  Use order on the univdiagonals to do this. */
  n = Univcoordlist_length(univdiagonals);
  array = Univcoordlist_to_array(univdiagonals,/*end*/0);
  univdiag_array = (Univdiag_T *) List_to_array(univdiags,NULL);
  queryfwd_elt_array = (T *) List_to_array(queryfwd_elt_list,NULL);

  Listpool_free_list(&queryfwd_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Univdiagpool_free_list(&univdiags,univdiagpool
			 univdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Univdiagpool_push */
  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_push */

#ifdef LARGE_GENOMES
  order = Sedgesort_order_uint8(array,n);
#else
  order = Sedgesort_order_uint4(array,n);
#endif
  
  /* univdiag_array contains newly created Univdiag_T objects, not in Elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    univdiag = univdiag_array[order[i]];
    process_seed(&(*univdiagonals_list),&(*auxinfo_list),

		 univdiag->univdiagonal,univdiag->qstart,univdiag->qend,univdiag->nmismatches,
		 queryfwd_set,/*queryrev_set*/NULL,
		 /*queryfwd_best_elt*/queryfwd_elt_array[order[i]],
		 /*queryrev_best_elt*/NULL,
		 querylength,auxinfopool,univdiagpool,univcoordlistpool,listpool);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Univdiagpool_free_univdiag(&(univdiag_array[i]),univdiagpool
			       univdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(univdiag_array);

  FREE(order);
  FREE(queryfwd_elt_array);
  FREE(array);

  return;
}


static void
extend_seeds_queryrev (Univcoordlist_T *univdiagonals_list, List_T *auxinfo_list,
		       List_T queryrev_best_elts, List_T queryrev_set, int querylength,

		       Auxinfopool_T auxinfopool, Univdiagpool_T univdiagpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool) {
  
  List_T univdiags = NULL, queryrev_elt_list = NULL;
  Univcoordlist_T univdiagonals = NULL; /* To determine order and avoid duplicates */
  Univcoord_T *array;

  int *order;
  Univdiag_T *univdiag_array, queryrev_univdiag, univdiag;

  int queryrev_nseeds, n, i, j;
  List_T q;
  T *queryrev_elt_array, queryrev_elt;

  
  debug(printf("Entering extend_seeds_queryrev with %d queryrev_best_elts\n",
	       List_length(queryrev_best_elts)));

  /* Handle only queryrev_best_elts */
  for (q = queryrev_best_elts; q != NULL; q = List_next(q)) {
    queryrev_elt = (T) List_head(q);
    queryrev_nseeds = queryrev_elt->n_all_univdiags;

    j = 0;
    while (j < queryrev_nseeds) {
      queryrev_univdiag = queryrev_elt->all_univdiags[j++];
	    
      /* Allow univdiagonals from queryfwd_elt */
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,queryrev_univdiag->univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));
      univdiags = Univdiagpool_push(univdiags,univdiagpool,queryrev_univdiag->qstart,queryrev_univdiag->qend,
				    queryrev_univdiag->nmismatches,queryrev_univdiag->univdiagonal
				    univdiagpool_trace(__FILE__,__LINE__));
      queryrev_elt_list = Listpool_push(queryrev_elt_list,listpool,(void *) queryrev_elt
					listpool_trace(__FILE__,__LINE__));
    }
  }
      
  /* univdiags contains newly created Univdiag_T objects, not in Elt_T objects */


  /* Want to sort two arrays simultaneously: univdiags and
     queryrev_elts.  Use order on the univdiagonals to do this. */
  n = Univcoordlist_length(univdiagonals);
  array = Univcoordlist_to_array(univdiagonals,/*end*/0);
  univdiag_array = (Univdiag_T *) List_to_array(univdiags,NULL);
  queryrev_elt_array = (T *) List_to_array(queryrev_elt_list,NULL);

  Listpool_free_list(&queryrev_elt_list,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  Univdiagpool_free_list(&univdiags,univdiagpool
			 univdiagpool_trace(__FILE__,__LINE__)); /* Allocated by Univdiagpool_push */
  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_push */

#ifdef LARGE_GENOMES
  order = Sedgesort_order_uint8(array,n);
#else
  order = Sedgesort_order_uint4(array,n);
#endif
  
  /* univdiag_array contains newly created Univdiag_T objects, not in Elt_T objects */

  i = 0;
  while (i < n) {
    /* Skip duplicates */
    j = i + 1;
    while (array[order[j]] == array[order[i]]) {
      j++;
    }

    univdiag = univdiag_array[order[i]];
    process_seed(&(*univdiagonals_list),&(*auxinfo_list),

		 univdiag->univdiagonal,univdiag->qstart,univdiag->qend,univdiag->nmismatches,
		 /*queryfwd_set*/NULL,queryrev_set,
		 /*queryfwd_best_elt*/NULL,
		 /*queryrev_best_elt*/queryrev_elt_array[order[i]],

		 querylength,auxinfopool,univdiagpool,univcoordlistpool,listpool);
    i = j;
  }

  for (i = 0; i < n; i++) {
    Univdiagpool_free_univdiag(&(univdiag_array[i]),univdiagpool
			       univdiagpool_trace(__FILE__,__LINE__));
  }
  FREE(univdiag_array);

  FREE(order);
  FREE(queryrev_elt_array);
  FREE(array);

  return;
}


#ifdef CHECK_ASSERTIONS
static void
check_ascending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  if (n > 0) {
    prev_coord = coords[0];
    for (i = 1; i < n; i++) {
      if (coords[i] <= prev_coord) {
	printf("Expecting ascending, but at %d, got %u <= %u\n",
	       i,coords[i],prev_coord);
	abort();
      }
      prev_coord = coords[i];
    }
  }
 
  return;
}
#endif


/* _univdiagonals_gplus and _univdiagonals_gminus should be aligned to match Kmer_exact1 */
void
Extension_search (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
		  Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
		  Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
		  Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,

		  Stage1_T stage1, Compress_T query_compress_fwd, Compress_T query_compress_rev,

		  int querylength, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		  Univdiagpool_T univdiagpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool) {

  Univcoordlist_T univdiagonals_gplus_list = NULL, univdiagonals_gminus_list = NULL;
  List_T auxinfo_gplus_list = NULL, auxinfo_gminus_list = NULL;
  Auxinfo_T auxinfo;

  List_T queryfwd_best_plus_elts, queryrev_best_plus_elts, queryfwd_best_minus_elts, queryrev_best_minus_elts;
  int max_iter = querylength / index1part;
  int i;

  get_elt_sets_queryfwd(&stage1->queryfwd_plus_set,&stage1->queryfwd_minus_set,
			&queryfwd_best_plus_elts,&queryfwd_best_minus_elts,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			max_iter,/*genestrand*/0,listpool,univdiagpool);
  
  get_elt_sets_queryrev(&stage1->queryrev_plus_set,&stage1->queryrev_minus_set,
			&queryrev_best_plus_elts,&queryrev_best_minus_elts,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			max_iter,/*genestrand*/0,listpool,univdiagpool);

  if (queryfwd_best_plus_elts == NULL && queryrev_best_plus_elts == NULL) {
    /* Skip */
  } else if (queryrev_best_plus_elts == NULL) {
    extend_seeds_queryfwd(&univdiagonals_gplus_list,&auxinfo_gplus_list,

			  queryfwd_best_plus_elts,stage1->queryfwd_plus_set,

			  querylength,auxinfopool,
			  univdiagpool,univcoordlistpool,listpool);

  } else if (queryfwd_best_plus_elts == NULL) {
    extend_seeds_queryrev(&univdiagonals_gplus_list,&auxinfo_gplus_list,

			  queryrev_best_plus_elts,stage1->queryrev_plus_set,

			  querylength,auxinfopool,
			  univdiagpool,univcoordlistpool,listpool);

  } else {
    extend_seeds_union(&univdiagonals_gplus_list,&auxinfo_gplus_list,

		       queryfwd_best_plus_elts,queryrev_best_plus_elts,
		       stage1->queryfwd_plus_set,stage1->queryrev_plus_set,
		       
		       querylength,auxinfopool,
		       univdiagpool,univcoordlistpool,listpool);
  }

  if (queryfwd_best_minus_elts == NULL && queryrev_best_minus_elts == NULL) {
    /* Skip */
  } else if (queryrev_best_minus_elts == NULL) {
    extend_seeds_queryfwd(&univdiagonals_gminus_list,&auxinfo_gminus_list,

			  queryfwd_best_minus_elts,stage1->queryfwd_minus_set,

			  querylength,auxinfopool,
			  univdiagpool,univcoordlistpool,listpool);

  } else if (queryfwd_best_minus_elts == NULL) {
    extend_seeds_queryrev(&univdiagonals_gminus_list,&auxinfo_gminus_list,

			  queryrev_best_minus_elts,stage1->queryrev_minus_set,

			  querylength,auxinfopool,
			  univdiagpool,univcoordlistpool,listpool);

  } else {
    extend_seeds_union(&univdiagonals_gminus_list,&auxinfo_gminus_list,

		       queryfwd_best_minus_elts,queryrev_best_minus_elts,

		       stage1->queryfwd_minus_set,stage1->queryrev_minus_set,

		       querylength,auxinfopool,
		       univdiagpool,univcoordlistpool,listpool);
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

  assert(Univcoordlist_length(univdiagonals_gplus_list) == List_length(auxinfo_gplus_list));
  assert(Univcoordlist_length(univdiagonals_gminus_list) == List_length(auxinfo_gminus_list));

  if ((*nunivdiagonals_gplus = Univcoordlist_length(univdiagonals_gplus_list)) == 0) {
    *_univdiagonals_gplus = (Univcoord_T *) NULL;
    *auxinfo_gplus = (Auxinfo_T *) NULL;
    *pathstores_gplus = (Pathstore_T *) NULL;

  } else {
    univdiagonals_gplus_list = Univcoordlist_reverse(univdiagonals_gplus_list);
    auxinfo_gplus_list = List_reverse(auxinfo_gplus_list);

    MALLOC_ALIGN(*_univdiagonals_gplus,(*nunivdiagonals_gplus)*sizeof(UINT4));
    Univcoordlist_fill_array(*_univdiagonals_gplus,univdiagonals_gplus_list);
    *auxinfo_gplus = (Auxinfo_T *) List_to_array(auxinfo_gplus_list,NULL); /* Returns n+1 elts */
    *pathstores_gplus = (Pathstore_T *) MALLOC((*nunivdiagonals_gplus)*sizeof(Pathstore_T));
    for (i = 0; i < *nunivdiagonals_gplus; i++) {
      (*pathstores_gplus)[i] = Pathstore_new(pathstorepool);
    }
  }

  if ((*nunivdiagonals_gminus = Univcoordlist_length(univdiagonals_gminus_list)) == 0) {
    *_univdiagonals_gminus = (Univcoord_T *) NULL;
    *auxinfo_gminus = (Auxinfo_T *) NULL;
    *pathstores_gminus = (Pathstore_T *) NULL;

  } else {
    univdiagonals_gminus_list = Univcoordlist_reverse(univdiagonals_gminus_list);
    auxinfo_gminus_list = List_reverse(auxinfo_gminus_list);

    MALLOC_ALIGN(*_univdiagonals_gminus,(*nunivdiagonals_gminus)*sizeof(UINT4));
    Univcoordlist_fill_array(*_univdiagonals_gminus,univdiagonals_gminus_list);
    *auxinfo_gminus = (Auxinfo_T *) List_to_array(auxinfo_gminus_list,NULL); /* Returns n+1 elts */
    *pathstores_gminus = (Pathstore_T *) MALLOC((*nunivdiagonals_gminus)*sizeof(Pathstore_T));
    for (i = 0; i < *nunivdiagonals_gminus; i++) {
      (*pathstores_gminus)[i] = Pathstore_new(pathstorepool);
    }
  }

#ifdef CHECK_ASSERTIONS
  check_ascending(*_univdiagonals_gplus,*nunivdiagonals_gplus);
  check_ascending(*_univdiagonals_gminus,*nunivdiagonals_gminus);
#endif

  if (*nunivdiagonals_gplus > 0) {
    stage1->extension_gplus = (Univcoord_T *) MALLOC((*nunivdiagonals_gplus)*sizeof(Univcoord_T));
    memcpy(stage1->extension_gplus,*_univdiagonals_gplus,(*nunivdiagonals_gplus)*sizeof(Univcoord_T));

    stage1->extension_qstart_gplus = (int *) MALLOC((*nunivdiagonals_gplus)*sizeof(int));
    stage1->extension_qend_gplus = (int *) MALLOC((*nunivdiagonals_gplus)*sizeof(int));
    for (i = 0; i < *nunivdiagonals_gplus; i++) {
      auxinfo = (*auxinfo_gplus)[i];
      stage1->extension_qstart_gplus[i] = auxinfo->qstart;
      stage1->extension_qend_gplus[i] = auxinfo->qend;
    }

    stage1->nextension_gplus = *nunivdiagonals_gplus;
  }

  if (*nunivdiagonals_gminus > 0) {
    stage1->extension_gminus = (Univcoord_T *) MALLOC((*nunivdiagonals_gminus)*sizeof(Univcoord_T));
    memcpy(stage1->extension_gminus,*_univdiagonals_gminus,(*nunivdiagonals_gminus)*sizeof(Univcoord_T));

    stage1->extension_qstart_gminus = (int *) MALLOC((*nunivdiagonals_gminus)*sizeof(int));
    stage1->extension_qend_gminus = (int *) MALLOC((*nunivdiagonals_gminus)*sizeof(int));
    for (i = 0; i < *nunivdiagonals_gminus; i++) {
      auxinfo = (*auxinfo_gminus)[i];
      stage1->extension_qstart_gminus[i] = auxinfo->qstart;
      stage1->extension_qend_gminus[i] = auxinfo->qend;
    }

    stage1->nextension_gminus = *nunivdiagonals_gminus;
  }

  Univcoordlistpool_free_list(&univdiagonals_gplus_list,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__));
  Listpool_free_list(&auxinfo_gplus_list,listpool
		     listpool_trace(__FILE__,__LINE__));

  Univcoordlistpool_free_list(&univdiagonals_gminus_list,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__));
  Listpool_free_list(&auxinfo_gminus_list,listpool
		     listpool_trace(__FILE__,__LINE__));

  debug(Stage1_dump(stage1,querylength));

  return;
}


