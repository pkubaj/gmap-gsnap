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
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
#include <sys/time.h>
#endif

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

#include "segmentelt.h"

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
#define MAX_CHAIN_OVERLAP 25


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

#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
static unsigned long long
extension_now_us (void) {
  struct timeval tv;

  gettimeofday(&tv,(struct timezone *) NULL);
  return (unsigned long long) tv.tv_sec * 1000000ULL + (unsigned long long) tv.tv_usec;
}
#endif


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


#define T Segmentelt_T

#if 0
static void
Elt_free (T *old) {

  Univdiag_free_array(&(*old)->all_univdiags);
  FREE(*old);

  return;
}

void
Elt_gc (List_T *set, Listpool_T listpool) {
  List_T p;
  T elt;

  for (p = *set; p != NULL; p = List_next(p)) {
    elt = (T) List_head(p);
    Elt_free(&elt);
  }
  Listpool_free_list(&(*set),listpool
		     listpool_trace(__FILE__,__LINE__)); /* allocated by Listpool_push */
  return;
}
#endif


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
static Segmentelt_T
Elt_read_queryfwd (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int querylength, int querystart,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Segmenteltpool_T segmenteltpool) {
  Segmentelt_T elt;
  int max_nmatches;
  Univcoord_T *univdiagonals, univdiagonal;

  int best_trimpos, trimpos;

  int nmatches, pos3;
  int i;
  

  debug9(printf("Got %d positions at querystart %d, plusp %d\n",n,querystart,plusp));
  if (n == 0) {
    return (Segmentelt_T) NULL;
  } else if (querystart >= querylength - index1part) {
    return (Segmentelt_T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = querystart;
    univdiagonals = (Univcoord_T *) MALLOC(n*sizeof(Univcoord_T));

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
      nmatches = index1part +
	Genomebits_consecutive_matches_rightward(genomebits,query_compress,univdiagonal,querylength,
						 /*pos5*/querystart+index1part,pos3,
						 plusp,genestrand);
      univdiagonals[i] = univdiagonal;
      trimpos = querystart + nmatches;
      debug9(printf("rightward plusp %d, univdiagonal is %u, querystart is %d + nmatches %d => queryend is %d\n",
		    plusp,univdiagonal,querystart,nmatches,querystart + nmatches));

#else
      /* Try Genomebits_trim_qend, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qend(&nmismatches,query_compress,
				     genomebits,univdiagonal,querylength,/*pos5*/querystart,pos3,
				     plusp,genestrand);
      debug9(printf("qend conservative: At univdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		    univdiagonal,querystart,pos3,nmatches[i],querystart+nmatches[i]));
      debug9(printf("qend aggressive: At univdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		    univdiagonal,querystart,pos3,trimpos,nmismatches,(trimpos - querystart) - nmismatches,trimpos));
      nmatches = (trimpos - querystart) - nmismatches;
      univdiagonals[i] = univdiagonal;
#endif

      if (nmatches > max_nmatches) {
	max_nmatches = nmatches;
	best_trimpos = trimpos;
      }
    }

    elt = Segmentelt_new(/*qstart*/querystart,/*qend*/best_trimpos,
			 univdiagonals,/*ndiagonals*/n,segmenteltpool);
    debug(printf("Making a queryfwd Elt %p with querystart %d, max_nmatches %d => max_queryend %d\n",
		 elt,elt->qstart,elt->qend));

    return elt;
  }
}


static Segmentelt_T
Elt_read_queryrev (
#ifdef LARGE_GENOMES
		   unsigned char *positions_high,
#endif
		   UINT4 *positions, int n, int diagterm, int queryend, int querylength,
		   Compress_T query_compress, bool plusp, int genestrand,
		   Segmenteltpool_T segmenteltpool) {
  Segmentelt_T elt;
  int max_nmatches;
  Univcoord_T *univdiagonals, univdiagonal;

  int best_trimpos, trimpos;

  int nmatches, pos5;
  int i;
  

  debug9(printf("Got %d positions at queryend %d, plusp %d\n",n,queryend,plusp));
  if (n == 0) {
    return (Segmentelt_T) NULL;
  } else if (queryend < index1part) {
    return (Segmentelt_T) NULL;
  } else {
    max_nmatches = 0;
    best_trimpos = queryend;
    univdiagonals = (Univcoord_T *) MALLOC(n*sizeof(Univcoord_T));

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
      nmatches = index1part +
	Genomebits_consecutive_matches_leftward(genomebits,query_compress,univdiagonal,querylength,
						pos5,/*pos3*/queryend - index1part,plusp,genestrand);
      univdiagonals[i] = univdiagonal;
      trimpos = queryend - nmatches;
      debug9(printf("leftward plusp %d, univdiagonal is %u, queryend is %d - nmatches %d => querystart is %d\n",
		    plusp,univdiagonal,queryend,nmatches,queryend - nmatches));

#else
      /* Try Genomebits_trim_qstart, which can handle mismatches and maximize extension */ 
      trimpos = Genomebits_trim_qstart(&nmismatches,query_compress,
				       genomebits,univdiagonal,querylength,pos5,/*pos3*/queryend,
				       plusp,genestrand);
      debug9(printf("qstart conservative: At univdiagonal %u, %d..%d, consecutive matches leftward yields %d nmatches => pos %d.\n",
		    univdiagonal,pos5,queryend,nmatches[i],queryend-nmatches[i]));
      debug9(printf("qstart aggressive: At univdiagonal %u, %d..%d, trimpos is %d with %d nmismatches (matches %d) => pos %d.\n",
		    univdiagonal,pos5,queryend,trimpos,nmismatches,(queryend - trimpos) - nmismatches,trimpos));
      nmatches = (queryend - trimpos) - nmismatches;
      univdiagonals[i] = univdiagonal;
#endif

      if (nmatches > max_nmatches) {
	max_nmatches = nmatches;
	best_trimpos = trimpos;
      }
    }

    elt = Segmentelt_new(/*qstart*/best_trimpos,/*qend*/queryend,
			 univdiagonals,/*ndiagonals*/n,segmenteltpool);
    debug(printf("Making a queryrev Elt %p with queryend %d, max_nmatches %d => min_querystart %d\n",
		 elt,elt->qend,elt->qstart));

    return elt;
  }
}


#ifdef DEBUG
static void
Elt_dump_set (List_T set) {
  List_T p;

  for (p = set; p != NULL; p = List_next(p)) {
    Segmentelt_dump((Segmentelt_T) List_head(p));
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


/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryfwd (List_T *plus_set, List_T *minus_set,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int max_iter, int genestrand,
		       Segmenteltpool_T segmenteltpool, Listpool_T listpool) {
  Segmentelt_T elt;
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


  /* I.  Race from plus and minus start to end.  Compute best_plus_nmatches and best_minus_nmatches */
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
			     indexdb,stage1->forward_signatures[queryoffset]);
#else
	  mod_npositions += stage1->plus_npositions[queryoffset] =
	    Indexdb_ptr(&stage1->plus_positions[queryoffset],indexdb,
			stage1->forward_signatures[queryoffset]);
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
					 segmenteltpool)) == NULL) {
	      plus_qpos[mod] += index1interval;
	    
	    } else {
#if 0
	      if (elt->nmatches > *best_plus_nmatches && elt->ndiagonals <= maxpaths_search) {
		*best_plus_nmatches = elt->nmatches;
	      }
#endif
	      debug(printf("Pushing elt %p onto plus set\n",elt));
	      *plus_set = Listpool_push(*plus_set,listpool,(void *) elt
					listpool_trace(__FILE__,__LINE__));
	      plus_qpos[mod] += (elt->qend - elt->qstart);
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
			     indexdb_nonstd,stage1->revcomp_signatures[querypos_rc]);
#else
	  mod_npositions = stage1->minus_npositions[querypos_rc] =
	    Indexdb_ptr(&stage1->minus_positions[querypos_rc],indexdb_nonstd,
			stage1->revcomp_signatures[querypos_rc]);
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
					 segmenteltpool)) == NULL) {
	      minus_qpos[mod] += index1interval;
	    
	    } else {
#if 0
	      if (elt->nmatches > *best_minus_nmatches && elt->ndiagonals <= maxpaths_search) {
		*best_minus_nmatches = elt->nmatches;
	      }
#endif
	      debug(printf("Pushing elt %p onto minus set\n",elt));
	      *minus_set = Listpool_push(*minus_set,listpool,(void *) elt
					 listpool_trace(__FILE__,__LINE__));

	      minus_qpos[mod] += (elt->qend - elt->qstart);
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

  return;
}


/* Note: A second try, starting at queryrev, may solve only a few more
   cases presented to this method (i.e., cases that cannot be solved
   by kmer-end search).  It therefore may be best to push these cases
   to the next method. */

/* Caller needs to call Elt_gc(&plus_set) and Elt_gc(&minus_set) */
static void
get_elt_sets_queryrev (List_T *plus_set, List_T *minus_set,
		       Stage1_T stage1, int querylength,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev, 
		       int max_iter, int genestrand,
		       Segmenteltpool_T segmenteltpool, Listpool_T listpool) {
  Segmentelt_T elt;
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
			       indexdb,stage1->forward_signatures[queryoffset]);
	    positions_high = stage1->plus_positions_high[queryoffset];
#else
	    npositions = stage1->plus_npositions[queryoffset] =
	      Indexdb_ptr(&stage1->plus_positions[queryoffset],indexdb,
			  stage1->forward_signatures[queryoffset]);
#endif
	    positions = stage1->plus_positions[queryoffset];
	    stage1->plus_retrievedp[queryoffset] = true;
	    debug(printf(" => retrieved %d positions\n",npositions));
	  }
	  
	  
	  debug(Stage1_dump(stage1,querylength));
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
					      segmenteltpool)) == NULL) {
	    total_npositions_plus += npositions;
	    plus_qpos[mod] -= index1interval;
	    
	  } else {
#if 0
	    if (elt->nmatches > *best_plus_nmatches && elt->ndiagonals <= maxpaths_search) {
	      *best_plus_nmatches = elt->nmatches;
	    }
#endif
	    *plus_set = Listpool_push(*plus_set,listpool,(void *) elt
				      listpool_trace(__FILE__,__LINE__));

	    plus_qpos[mod] -= (elt->qend - elt->qstart);
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
			       indexdb_nonstd,stage1->revcomp_signatures[querypos_rc]);
	    positions_high = stage1->minus_positions_high[querypos_rc];
#else
	    npositions = stage1->minus_npositions[querypos_rc] =
	      Indexdb_ptr(&stage1->minus_positions[querypos_rc],indexdb_nonstd,
			  stage1->revcomp_signatures[querypos_rc]);
#endif
	    positions = stage1->minus_positions[querypos_rc];
	    stage1->minus_retrievedp[querypos_rc] = true;
	    debug(printf(" => retrieved %d positions\n",npositions));
	  }
	  
	  
	  debug(Stage1_dump(stage1,querylength));
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
					      segmenteltpool)) == NULL) {
	    total_npositions_minus += npositions;
	    minus_qpos[mod] -= index1interval;
	    
	  } else {
#if 0
	    if (elt->nmatches > *best_minus_nmatches && elt->ndiagonals <= maxpaths_search) {
	      *best_minus_nmatches = elt->nmatches;
	    }
#endif
	    *minus_set = Listpool_push(*minus_set,listpool,(void *) elt
				       listpool_trace(__FILE__,__LINE__));

	    minus_qpos[mod] -= (elt->qend - elt->qstart);
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
      if (coords[i] < prev_coord) {
	/* Using < and not <=, because now allowing for duplicates */
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
/* May want to store chain in auxinfo, or use auxinfo only for segment search */
void
Extension_search (int ***indices_by_nsplices_gplus, int **nindices_by_nsplices_gplus,
		  int ***indices_by_nsplices_gminus, int **nindices_by_nsplices_gminus,

		  Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
		  Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
		  Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
		  Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
		      
		  Stage1_T stage1, Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  int querylength, int max_nsplices, bool first_read_p, bool retain_terminal_p,
		  bool prune_equal_score_p,

		  Segmenteltpool_T segmenteltpool,
		  Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		  Univcoordlistpool_T univcoordlistpool, Univdiagpool_T univdiagpool,
		  Intlistpool_T intlistpool, Listpool_T listpool) {

  Segmentelt_T *plus_eltarray, *minus_eltarray;
  int plus_narray, minus_narray;
  List_T plus_elt_set = NULL, minus_elt_set = NULL;

  int max_iter = querylength / index1part;
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  unsigned long long stage_start_us;
#endif

  debug(printf("Extension_search starting.  Stage1_T status:\n"));
  debug(Stage1_dump(stage1,querylength));

  /* Previously saved as stage1->extension_plus_set */
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  stage_start_us = extension_now_us();
#endif
  get_elt_sets_queryfwd(&plus_elt_set,&minus_elt_set,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			max_iter,/*genestrand*/0,segmenteltpool,listpool);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[EXTENSION_STATS] queryfwd_us=%llu plus_set=%d minus_set=%d\n",
	  extension_now_us() - stage_start_us,List_length(plus_elt_set),List_length(minus_elt_set));
  stage_start_us = extension_now_us();
#endif
  
  get_elt_sets_queryrev(&plus_elt_set,&minus_elt_set,
			stage1,querylength,query_compress_fwd,query_compress_rev,
			max_iter,/*genestrand*/0,segmenteltpool,listpool);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,"[EXTENSION_STATS] queryrev_us=%llu plus_set=%d minus_set=%d prune_equal=%d\n",
	  extension_now_us() - stage_start_us,List_length(plus_elt_set),List_length(minus_elt_set),
	  prune_equal_score_p);
#endif

  if ((plus_narray = List_length(plus_elt_set)) == 0) {
    *nunivdiagonals_gplus = 0;
    *_univdiagonals_gplus = (Univcoord_T *) NULL;
    *auxinfo_gplus = (Auxinfo_T *) NULL;
    *pathstores_gplus = (Pathstore_T *) NULL;
    *indices_by_nsplices_gplus = (int **) CALLOC(max_nsplices + 1,sizeof(int *));
    *nindices_by_nsplices_gplus = (int *) CALLOC(max_nsplices + 1,sizeof(int));

  } else if (first_read_p == true) {
    /* For first read, plus is low end, so compute left chains */
    plus_eltarray = (Segmentelt_T *) List_to_array(plus_elt_set,NULL);
    qsort(plus_eltarray,plus_narray,sizeof(Segmentelt_T),Segmentelt_cmp);

    *indices_by_nsplices_gplus =
      Segmentelt_left_chains(&(*nindices_by_nsplices_gplus),max_nsplices,
			     &(*_univdiagonals_gplus),&(*auxinfo_gplus),&(*pathstores_gplus),
			     &(*nunivdiagonals_gplus),plus_eltarray,plus_narray,querylength,retain_terminal_p,
			     prune_equal_score_p,
			     univcoordlistpool,auxinfopool,pathstorepool,
			     univdiagpool,intlistpool,listpool);
    FREE(plus_eltarray);
    Segmentelt_gc(plus_elt_set,segmenteltpool,intlistpool,listpool);

  } else {
    plus_eltarray = (Segmentelt_T *) List_to_array(plus_elt_set,NULL);
    qsort(plus_eltarray,plus_narray,sizeof(Segmentelt_T),Segmentelt_cmp);

    *indices_by_nsplices_gplus =
      Segmentelt_right_chains(&(*nindices_by_nsplices_gplus),max_nsplices,
			      &(*_univdiagonals_gplus),&(*auxinfo_gplus),&(*pathstores_gplus),
			      &(*nunivdiagonals_gplus),plus_eltarray,plus_narray,querylength,retain_terminal_p,
			      prune_equal_score_p,
			      univcoordlistpool,auxinfopool,pathstorepool,
			      univdiagpool,intlistpool,listpool);
    FREE(plus_eltarray);
    Segmentelt_gc(plus_elt_set,segmenteltpool,intlistpool,listpool);
  }

  if ((minus_narray = List_length(minus_elt_set)) == 0) {
    *nunivdiagonals_gminus = 0;
    *_univdiagonals_gminus = (Univcoord_T *) NULL;
    *auxinfo_gminus = (Auxinfo_T *) NULL;
    *pathstores_gminus = (Pathstore_T *) NULL;
    *indices_by_nsplices_gminus = (int **) CALLOC(max_nsplices + 1,sizeof(int *));
    *nindices_by_nsplices_gminus = (int *) CALLOC(max_nsplices + 1,sizeof(int));

  } else if (first_read_p == true) {
    /* For first read, minus is high end, so compute right chains */
    minus_eltarray = (Segmentelt_T *) List_to_array(minus_elt_set,NULL);
    qsort(minus_eltarray,minus_narray,sizeof(Segmentelt_T),Segmentelt_cmp);

    *indices_by_nsplices_gminus =
      Segmentelt_right_chains(&(*nindices_by_nsplices_gminus),max_nsplices,
			      &(*_univdiagonals_gminus),&(*auxinfo_gminus),&(*pathstores_gminus),
			      &(*nunivdiagonals_gminus),minus_eltarray,minus_narray,querylength,retain_terminal_p,
			      prune_equal_score_p,
			      univcoordlistpool,auxinfopool,pathstorepool,
			      univdiagpool,intlistpool,listpool);
    FREE(minus_eltarray);
    Segmentelt_gc(minus_elt_set,segmenteltpool,intlistpool,listpool);

  } else {
    minus_eltarray = (Segmentelt_T *) List_to_array(minus_elt_set,NULL);
    qsort(minus_eltarray,minus_narray,sizeof(Segmentelt_T),Segmentelt_cmp);

    *indices_by_nsplices_gminus =
      Segmentelt_left_chains(&(*nindices_by_nsplices_gminus),max_nsplices,
			     &(*_univdiagonals_gminus),&(*auxinfo_gminus),&(*pathstores_gminus),
			     &(*nunivdiagonals_gminus),minus_eltarray,minus_narray,querylength,retain_terminal_p,
			     prune_equal_score_p,
			     univcoordlistpool,auxinfopool,pathstorepool,
			     univdiagpool,intlistpool,listpool);
    FREE(minus_eltarray);
    Segmentelt_gc(minus_elt_set,segmenteltpool,intlistpool,listpool);
  }

  return;
}
