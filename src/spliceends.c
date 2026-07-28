#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "spliceends.h"

#include <stdio.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "sense.h"
#include "genome_sites.h"
#include "genomebits_count.h"
#include "genomebits_mismatches.h"
#include "genomebits_trim.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "univcoord.h"
#include "univdiagdef.h"

#ifdef LARGE_GENOMES
#include "intersect-lower-large.h"
#include "intersect-higher-large.h"
#else
#include "intersect-lower-small.h"
#include "intersect-higher-small.h"
#endif

#include "merge-method.h"
#ifdef USE_HEAP_MERGE
#include "merge-diagonals-heap-uint4.h"
#elif defined(USE_SIMD_MERGE)
#include "merge-diagonals-simd-uint4.h"
#endif

#ifdef LARGE_GENOMES
#ifdef USE_HEAP_MERGE
#include "merge-diagonals-heap-uint8.h"
#elif defined(USE_SIMD_MERGE)
#include "merge-diagonals-simd-uint8.h"
#endif
#endif


/* Trimming at chromosome bounds can cause problems with endpoints, so
   do all trimming later */
/* #define TRIM_AT_CHROMOSOME_BOUNDS 1 */

#define USE_VECTORPOOL 1

#define MAX_SITES 3
#define MAX_OUTER_PARTNERS_PER_SITE 5
#define MAX_INNER_PARTNERS_PER_SITE 30
#define MAX_NSPLICEENDS MAX_SITES*MAX_INNER_PARTNERS_PER_SITE

#define ACCEPTABLE_TRIM 3
#define SUFFICIENT_NMATCHES 8 /* For end indel or resolve */

#define DEFAULT_MEDIAL_SPLICESITE_PROB 0.90
#define DEFAULT_DISTAL_SPLICESITE_PROB 0.90

#define SALVAGE_MEDIAL_SPLICESITE_PROB 0.80
#define SALVAGE_DISTAL_SPLICESITE_PROB 0.80

#define PROB_SLOP 0.05
#define OUTER_PROB_SLOP 0.05

#define MAX_NCONSECUTIVE_CLOSE 6
#define MAX_NCONSECUTIVE_FAR 6	/* Needs to be generous to find splices.  Was 20, but now checking for min_nmismatches */

/* #define END_SPLICESITE_SEARCH_MM 1 -- Amount to search in the trimmed area */
#define END_SPLICESITE_SEARCH 10   /* Amount to search in the matching area (beyond the correct splice site) */

#define MIN_EXON_LENGTH 20	/* Minimum length from the exon_origin before we accept a splice site */
#define MIN_INTRON_LENGTH 9

/* Designed to allow 1 match to offset 1 mismatch.  To handle 2 matches vs 2 mismatches, penalize for multiple mismatches */
#define TRIM_MATCH_SCORE 1
#define TRIM_MISMATCH_SCORE_LAST -1 /* Requires 1 match to compensate */
#define TRIM_MISMATCH_SCORE_MULT -3 /* Requires 3 matches to compensate */


/* Known and novel splicing at ends */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Extension using indexdb */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Resolve */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* Trimming of univdiags */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Trimming nosplice */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Trimming at ends */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Extending using indexdb */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Do need to trim univdiagonals at chromosome bounds, but not pos5 or pos3 */
#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))


static bool *circularp;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Univcoord_T genomelength;

static Indexdb_T indexdb;
static Localdb_T localdb;

static int index1part;
static int index1interval;

static Chrpos_T shortsplicedist;
static bool splicingp;

static int max_insertionlen;
static int max_deletionlen;

static bool allow_soft_clips_p;

static int distal_nmismatches_allowed = 1;

#define T Spliceends_T

void
Spliceends_setup (bool *circularp_in,
		  Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		  Univcoord_T genomelength_in, Indexdb_T indexdb_in,
		  int index1part_in, int index1interval_in,
		  int max_insertionlen_in, int max_deletionlen_in,
		  Chrpos_T shortsplicedist_in, bool allow_soft_clips_p_in,
		  bool novelsplicingp, bool knownsplicingp) {

  circularp = circularp_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  genomelength = genomelength_in;

  indexdb = indexdb_in;
  index1part = index1part_in;
  index1interval = index1interval_in;

  localdb = (Localdb_T) NULL;
  shortsplicedist = shortsplicedist_in;

  max_insertionlen = max_insertionlen_in;
  max_deletionlen = max_deletionlen_in;

  allow_soft_clips_p = allow_soft_clips_p_in;

  if (novelsplicingp == true || knownsplicingp == true) {
    splicingp = true;
  } else {
    splicingp = false;
  }

  return;
}


void
Spliceends_localdb_setup (Localdb_T localdb_in) {
  localdb = localdb_in;
  return;
}


/* was Spliceends_free (T *old, Spliceendspool_T spliceendspool) */

void
Spliceends_free (T *old) {
  if (*old) {
    /* Allocated by Vectorpool_T, and not reclaiming due to variable lengths */

#ifdef USE_VECTORPOOL    
    /* Vectorpool frees its memory only after the query is handled */
#else
    FREE((*old)->splice_qpos);
    FREE((*old)->medial_nmismatches);
    FREE((*old)->partners);
#endif

#if 0
#ifdef USE_VECTORPOOL
    /* Vectorpool frees its memory only after the query is handled */
#else
    FREE((*old)->matchlengths);
    FREE((*old)->distal_lengths);
    FREE((*old)->distal_trimpos);
    FREE((*old)->distal_nmismatches);
#endif
#else
    Univdiag_free_array(&(*old)->distal_univdiags); /* pointers */
#endif

#ifdef USE_VECTORPOOL
    /* Vectorpool frees its memory only after the query is handled */
#else
    FREE((*old)->medial_probs);
    FREE((*old)->distal_probs);

    FREE((*old)->medial_support);
    FREE((*old)->distal_support);

    FREE((*old)->mismatch_positions_left);
    FREE((*old)->mismatch_positions_right);
#endif

#if 0
    Spliceendspool_free_spliceends(&(*old),spliceendspool
				   spliceendspool_trace(__FILE__,__LINE__));
#else
    FREE(*old);
#endif
  }

  return;
}


/* was Spliceends_new (int id, int querylength, Vectorpool_T vectorpool, Spliceendspool_T spliceendspool) */
/* For Spliceendsgen_checkout and Spliceendsgen_return to work, the sizes for each field must be the same */
T
Spliceends_new (int id, int querylength, Vectorpool_T vectorpool) {
#if 0
  T new = Spliceendspool_new_spliceends(spliceendspool
					spliceendspool_trace(__FILE__,__LINE__));
#else
  T new = (T) MALLOC(sizeof(*new));
#endif

  /* Multiply by index1interval for extend procedures, which test various mods */
  int n = index1interval * MAX_NSPLICEENDS;

  new->id = id;		    /* For debugging of Spliceendsgen_T */
  new->checkedout_p = true;	/* Because created by Spliceendsgen_T only when needed */

  new->boundedp = false;

  /* Computed later */
  /* new->nspliceends = nspliceends; */
  /* new->splicetype = splicetype; */
  /* new->sensedir = sensedir; */

  /* MISMATCH_EXTRA defined in genomebits_mismatches.h */

#ifdef USE_VECTORPOOL

  new->mismatch_positions_left = Vectorpool_new_intvector(vectorpool,querylength + MISMATCH_EXTRA);
  new->mismatch_positions_right = Vectorpool_new_intvector(vectorpool,querylength + MISMATCH_EXTRA);

  new->splice_qpos = Vectorpool_new_intvector(vectorpool,n + 1);
  new->medial_nmismatches = Vectorpool_new_intvector(vectorpool,n + 1);
  new->partners = Vectorpool_new_univcoordvector(vectorpool,n + 1); /* n + 1 needed for Sedgesort_order */

  /* new->matchlengths = Vectorpool_new_intvector(vectorpool,n + 1); */
  /* new->distal_lengths = Vectorpool_new_intvector(vectorpool,n + 1); */
  /* new->distal_trimpos = Vectorpool_new_intvector(vectorpool,n + 1); */
  /* new->distal_nmismatches = Vectorpool_new_intvector(vectorpool,n + 1); */

  new->medial_probs = Vectorpool_new_doublevector(vectorpool,n);
  new->distal_probs = Vectorpool_new_doublevector(vectorpool,n);

  new->medial_support = Vectorpool_new_intvector(vectorpool,n);
  new->distal_support = Vectorpool_new_intvector(vectorpool,n);

#else

  new->mismatch_positions_left = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));
  new->mismatch_positions_right = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));

  new->splice_qpos = (int *) MALLOC((n + 1)*sizeof(int));
  new->medial_nmismatches = (int *) MALLOC((n + 1)*sizeof(int));
  new->partners = (Univcoord_T *) MALLOC((n + 1)*sizeof(Univcoord_T));

  /* new->matchlengths = (int *) MALLOC((n + 1)*sizeof(int)); */
  /* new->distal_lengths = (int *) MALLOC((n + 1)*sizeof(int)); */
  /* new->distal_trimpos = (int *) MALLOC((n + 1)*sizeof(int)); */
  /* new->distal_nmismatches = (int *) MALLOC((n + 1)*sizeof(int)); */

  new->medial_probs = (double *) MALLOC(n * sizeof(double));
  new->distal_probs = (double *) MALLOC(n * sizeof(double));

  new->medial_support = (int *) MALLOC(n * sizeof(int));
  new->distal_support = (int *) MALLOC(n * sizeof(int));
#endif
  
  new->distal_univdiags = Univdiag_new_array(n);

  return new;
}


static int
make_unique (Univcoord_T *diagonals, int ndiagonals) {
  int k = 0, i, j;

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    diagonals[k++] = diagonals[i];
    
    i = j;
  }

  return k;
}


int
Spliceends_middle_plus (Univcoord_T **diagonals,
			Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			int qstart, int qend, int querylength,
			Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			Compress_T query_compress, char *queryptr,
			Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			Localdb_T localdb, int localdb_nmismatches_allowed) {

  int ndiagonals, starti, endi, i;
  Univcoord_T *_diagonals;
  int querystart, queryend;
  int nstreams;

  bool sortedp, trimmedp;
  int matchlength, local_nmismatches;
  int total_npositions_plus, total_npositions_minus;


  debug2(printf("Spliceends_middle_plus for univdiagonals %u to %u\n",low_univdiagonal,high_univdiagonal));

  assert(stage1->all_oligos_gen_filledp == true);
  if (stage1->all_positions_gen_filledp == false) {
    Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				  stage1,querylength,/*genestrand*/0);
  }

  querystart = qstart;
  queryend = qend;

  if ((nstreams = (queryend - index1part) - querystart + 1) > 0) {
#ifdef LARGE_GENOMES
    _diagonals = Merge_diagonals_large(&ndiagonals,&(stage1->plus_positions_high[querystart]),
				       &(stage1->plus_positions[querystart]),&(stage1->plus_npositions[querystart]),
				       &(stage1->plus_diagterms[querystart]),nstreams,stage1->mergeinfo);
#else
    _diagonals = Merge_diagonals(&ndiagonals,&(stage1->plus_positions[querystart]),&(stage1->plus_npositions[querystart]),
				 &(stage1->plus_diagterms[querystart]),nstreams,stage1->mergeinfo);
#endif

    if (ndiagonals == 0) {
      FREE_ALIGN(_diagonals);

    } else {
      ndiagonals = make_unique(_diagonals,ndiagonals);

      starti = 0;
      while (starti < ndiagonals && _diagonals[starti] < low_univdiagonal) {
	starti++;
      }

      endi = starti;
      while (endi < ndiagonals && _diagonals[endi] < high_univdiagonal) {
	endi++;
      }

      if ((ndiagonals = endi - starti) <= 0) {
	/* No diagonals in the given range */
	FREE_ALIGN(_diagonals);

      } else {
	*diagonals = (Univcoord_T *) MALLOC(ndiagonals*sizeof(Univcoord_T));
	memcpy(*diagonals,&(_diagonals[starti]),ndiagonals*sizeof(Univcoord_T));
	FREE_ALIGN(_diagonals);
	return ndiagonals;
      }
    }
  }

  if (localdb == NULL) {
    *diagonals = (Univcoord_T *) NULL;
    return 0;

  } else if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				       localdb,diagonals_alloc,localdb_alloc,
				       stage1,streamspace_max_alloc,streamspace_alloc,
				       queryptr,/*pos5*/qstart,/*pos3*/qend,querylength,
				       /*anchor_univdiagonal*/0,low_univdiagonal,high_univdiagonal,
				       query_compress,/*plusp*/true,/*genestrand*/0,genomebits,localdb_nmismatches_allowed,
				       /*extend5p*/true,/*trim5p*/true,/*trim3p*/true)) == 0) {
    *diagonals = (Univcoord_T *) NULL;
    return 0;

  } else {
    *diagonals = (Univcoord_T *) MALLOC(ndiagonals*sizeof(Univcoord_T));
    for (i = 0; i < ndiagonals; i++) {
      (*diagonals)[i] = univdiags_alloc[i]->univdiagonal;
    }
    return ndiagonals;
  }
}


int
Spliceends_middle_minus (Univcoord_T **diagonals,
			 Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			 int qstart, int qend, int querylength,
			 Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			 Compress_T query_compress, char *queryptr,
			 Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			 Localdb_T localdb, int localdb_nmismatches_allowed) {
  
  int ndiagonals, starti, endi, i;
  Univcoord_T *_diagonals;
  int querystart, queryend;
  int nstreams;

  bool sortedp, trimmedp;
  int matchlength, local_nmismatches;
  int total_npositions_plus, total_npositions_minus;


  debug2(printf("Spliceends_middle_minus for univdiagonals %u to %u\n",low_univdiagonal,high_univdiagonal));

  assert(stage1->all_oligos_gen_filledp == true);
  if (stage1->all_positions_gen_filledp == false) {
    Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				  stage1,querylength,/*genestrand*/0);
  }

  querystart = querylength - qend;
  queryend = querylength - qstart;

  if ((nstreams = (queryend - index1part) - querystart + 1) > 0) {
#ifdef LARGE_GENOMES
    _diagonals = Merge_diagonals_large(&ndiagonals,&(stage1->minus_positions_high[querystart]),
				       &(stage1->minus_positions[querystart]),&(stage1->minus_npositions[querystart]),
				       &(stage1->minus_diagterms[querystart]),nstreams,stage1->mergeinfo);
#else
    _diagonals = Merge_diagonals(&ndiagonals,&(stage1->minus_positions[querystart]),&(stage1->minus_npositions[querystart]),
				 &(stage1->minus_diagterms[querystart]),nstreams,stage1->mergeinfo);
#endif

    if (ndiagonals == 0) {
      FREE_ALIGN(_diagonals);

    } else {
      ndiagonals = make_unique(_diagonals,ndiagonals);

      starti = 0;
      while (starti < ndiagonals && _diagonals[starti] < low_univdiagonal) {
	starti++;
      }
    
      endi = starti;
      while (endi < ndiagonals && _diagonals[endi] < high_univdiagonal) {
	endi++;
      }

      if ((ndiagonals = endi - starti) <= 0) {
	/* No diagonals in the given range */
	FREE_ALIGN(_diagonals);

      } else {
	*diagonals = (Univcoord_T *) MALLOC(ndiagonals*sizeof(Univcoord_T));
	memcpy(*diagonals,&(_diagonals[starti]),ndiagonals*sizeof(Univcoord_T));
	FREE_ALIGN(_diagonals);
	return ndiagonals;
      }
    }
  }

  if (localdb == NULL) {
    *diagonals = (Univcoord_T *) NULL;
    return 0;

  } else if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				       localdb,diagonals_alloc,localdb_alloc,
				       stage1,streamspace_max_alloc,streamspace_alloc,
				       queryptr,/*pos5*/qstart,/*pos3*/qend,querylength,
				       /*anchor_univdiagonal*/0,low_univdiagonal,high_univdiagonal,
				       query_compress,/*plusp*/false,/*genestrand*/0,genomebits,localdb_nmismatches_allowed,
				       /*extend5p*/true,/*trim5p*/true,/*trim3p*/true)) == 0) {
    *diagonals = (Univcoord_T *) NULL;
    return 0;

  } else {
    *diagonals = (Univcoord_T *) MALLOC(ndiagonals*sizeof(Univcoord_T));
    for (i = 0; i < ndiagonals; i++) {
      (*diagonals)[i] = univdiags_alloc[i]->univdiagonal;
    }
    return ndiagonals;
  }
}


static inline int
min_spliceend (int *splice_qpos, int nspliceends) {
  int min;
  int i;

  min = splice_qpos[0];
  for (i = 1; i < nspliceends; i++) {
    if (splice_qpos[i] < min) {
      min = splice_qpos[i];
    }
  }
  return min;
}

static inline int
max_spliceend (int *splice_qpos, int nspliceends) {
  int max;
  int i;

  max = splice_qpos[0];
  for (i = 1; i < nspliceends; i++) {
    if (splice_qpos[i] > max) {
      max = splice_qpos[i];
    }
  }
  return max;
}


static int
univcoord_ascending_cmp (const void *x, const void *y) {
  Univcoord_T a = * (Univcoord_T *) x;
  Univcoord_T b = * (Univcoord_T *) y;

  if (a < b) {
    return -1;
  } else if (b < a) {
    return +1;
  } else {
    return 0;
  }
}


static int
extend_trim5 (Univdiag_T **univdiags, int start_qpos, Stage1_T stage1,
	      Univcoord_T univdiagonal, int querylength, bool plusp,
	      int slop, int insertion_slop, Univdiagpool_T univdiagpool) {

  int nmerged, nunivdiags, ndiagonals, npositions, i, j;
  List_T univdiag_list = NULL;
  Univcoord_T *diagonals;
  Univdiag_T univdiag;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
  UINT4 *positions;
#else
  Univcoord_T *positions;
#endif
  int querypos, qpos;


  if (plusp == true) {
    /* plus */
    debug9(printf("extend_trim5_plus at univdiagonal %u\n",univdiagonal));

    /* Check from medial to distal.  Previously greedy, now exhaustive */
    for (qpos = start_qpos; qpos >= index1part; qpos--) {
      querypos = qpos - index1part;
      debug9(printf("Testing kmer at qpos %d => querypos %d\n",qpos,querypos));

      if (stage1->validp[querypos] == false) {
	/* Skip */
      } else {
	if (stage1->plus_retrievedp[querypos] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->plus_positions_high[querypos];
#endif
	  positions = stage1->plus_positions[querypos];
	  npositions = stage1->plus_npositions[querypos];
	} else {
	  assert(stage1->plus_positions[querypos] == NULL);
#ifdef LARGE_GENOMES
	  npositions = stage1->plus_npositions[querypos] =
	    Indexdb_largeptr(&stage1->plus_positions_high[querypos],&stage1->plus_positions[querypos],
			     indexdb,stage1->forward_oligos[querypos]);
	  positions_high = stage1->plus_positions_high[querypos];
#else
	  npositions = stage1->plus_npositions[querypos] =
	    Indexdb_ptr(&stage1->plus_positions[querypos],indexdb,
			stage1->forward_oligos[querypos]);
#endif
	  positions = stage1->plus_positions[querypos];
	  stage1->plus_retrievedp[querypos] = true;
	}
      
	if (npositions > 0) {
	  diagonals = (Univcoord_T *) MALLOC(npositions*sizeof(Univcoord_T));
	  ndiagonals = Intersect_lower(diagonals,
#ifdef LARGE_GENOMES
				       positions_high,
#endif
				       positions,npositions,
				       /*diagterm1, plus*/querylength - querypos,
				       /*set2*/&univdiagonal,/*length2*/1,
				       slop,insertion_slop);
	  for (i = 0; i < ndiagonals; i++) {
	    univdiag_list = Univdiagpool_push(univdiag_list,univdiagpool,/*qstart*/qpos - index1part,/*qend*/qpos,
					      /*nmismatches*/0,diagonals[i]
					      univdiagpool_trace(__FILE__,__LINE__));
	  }
	  FREE(diagonals);
	}
      }
    }
    
  } else {
    /* minus */
    debug9(printf("extend_trim5_minus at univdiagonal %u\n",univdiagonal));

    /* Check from medial to distal */
    for (qpos = start_qpos; qpos >= index1part; qpos--) {
      querypos = querylength - qpos;
      debug9(printf("Testing kmer at qpos %d => querypos %d\n",qpos,querypos));

      if (stage1->validp[querypos] == false) {
	/* Skip */
      } else {
	if (stage1->minus_retrievedp[querypos] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->minus_positions_high[querypos];
#endif
	  positions = stage1->minus_positions[querypos];
	  npositions = stage1->minus_npositions[querypos];
	} else {
	  assert(stage1->minus_positions[querypos] == NULL);
#ifdef LARGE_GENOMES
	  npositions = stage1->minus_npositions[querypos] =
	    Indexdb_largeptr(&stage1->minus_positions_high[querypos],&stage1->minus_positions[querypos],
			     indexdb,stage1->revcomp_oligos[querypos]);
	  positions_high = stage1->minus_positions_high[querypos];
#else
	  npositions = stage1->minus_npositions[querypos] =
	    Indexdb_ptr(&stage1->minus_positions[querypos],indexdb,
			stage1->revcomp_oligos[querypos]);
#endif
	  positions = stage1->minus_positions[querypos];
	  stage1->minus_retrievedp[querypos] = true;
	}
      
	if (npositions > 0) {
	  diagonals = (Univcoord_T *) MALLOC(npositions*sizeof(Univcoord_T));
	  ndiagonals = Intersect_lower(diagonals,
#ifdef LARGE_GENOMES
				       positions_high,
#endif
				       positions,npositions,
				       /*diagterm1, minus*/querypos + index1part,
				       /*set2*/&univdiagonal,/*length2*/1,
				       slop,insertion_slop);
	  for (i = 0; i < ndiagonals; i++) {
	    univdiag_list = Univdiagpool_push(univdiag_list,univdiagpool,/*qstart*/qpos - index1part,/*qend*/qpos,
					      /*nmismatches*/0,diagonals[i]
					      univdiagpool_trace(__FILE__,__LINE__));
	  }
	  FREE(diagonals);
	}
      }
    }
  }

  if (univdiag_list == 0) {
    *univdiags = (Univdiag_T *) NULL;
    return 0;

  } else {
    nunivdiags = List_length(univdiag_list);
    *univdiags = Univdiag_new_array(nunivdiags);
    Univdiag_transfer_list_to_array(*univdiags,univdiag_list);
    Univdiagpool_free_list(&univdiag_list,univdiagpool
			   univdiagpool_trace(__FILE__,__LINE__));
    if (nunivdiags == 1) {
      return 1;

    } else {
      /* Merge univdiags.  Output in descending order */
      qsort((*univdiags)[0],nunivdiags,sizeof(struct Univdiag_T),Univdiag_struct_diagonal_rev_cmp);
      
#ifdef DEBUG9
      printf("Before merging: %d univdiags\n",nunivdiags);
      for (i = 0; i < nunivdiags; i++) {
	univdiag = (*univdiags)[i];
	printf("%u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
      }
      printf("\n");
#endif
    
      nmerged = 0; i = 0;
      while (i < nunivdiags) {
	univdiag = (*univdiags)[i];

	j = i + 1;
	while (j < nunivdiags && (*univdiags)[j]->univdiagonal == univdiag->univdiagonal) {
	  univdiag->qend = (*univdiags)[j]->qend;
	  j++;
	}

	Univdiag_transfer((*univdiags)[nmerged++],univdiag);
	i = j;
      }

#ifdef DEBUG9
      printf("After merging: returning %d univdiags\n",nmerged);
      for (i = 0; i < nmerged; i++) {
	univdiag = (*univdiags)[i];
	printf("%u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
      }
      printf("\n");
#endif

      return nmerged;
    }
  }
}


static Univdiag_T
novel_trim5_indel (Univcoord_T outbound_genomicpos, Univcoord_T univdiagonal, int querylength,
		   Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
		   Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		   Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		   bool plusp, int genestrand, int localdb_nmismatches_allowed) {

  Univdiag_T best_univdiag = NULL;
  Univcoord_T deletion_univdiagonal, insertion_univdiagonal;
  bool sortedp, trimmedp;

  int local_nmismatches, nmatches, matchlength;
  int indel_qpos;
  int ndiagonals, i;
  int best_adj, adj;
  Univcoord_T left = univdiagonal - querylength;

  /* Test for deletion */
  indel_qpos = outbound_genomicpos - left;

  debug9(printf("Start with univdiagonal %u\n",univdiagonal));
  deletion_univdiagonal = subtract_bounded(univdiagonal,max_deletionlen,chroffset);
  debug9(printf("Subtract %d (bounded by %u) to yield (low) deletion_univdiagonal %u\n",
		max_deletionlen,chroffset,deletion_univdiagonal));

  debug9(printf("Testing indel_qpos %d for deletion\n",indel_qpos));
  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				localdb,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,
				queryptr,/*pos5*/0,/*pos3*/indel_qpos,querylength,
				/*anchor_univdiagonal*/univdiagonal,
				/*low_univdiagonal*/deletion_univdiagonal,/*high_univdiagonal*/univdiagonal,
				query_compress,plusp,genestrand,genomebits,localdb_nmismatches_allowed,
				/*extend5p*/true,/*trim5p*/true,/*trim3p*/false)) == 0) {
    /* Skip */

  } else if ((nmatches = matchlength - local_nmismatches) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*local_nmismatches) {
    debug9(printf("(1) Got %d localdb diagonals with matchlength %d and local_nmismatches %d => Skipping\n",
		  ndiagonals,matchlength,local_nmismatches));
    /* Skip */

  } else {
    debug9(printf("(2) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		  ndiagonals,matchlength,local_nmismatches));

    best_adj = 0;
    for (i = ndiagonals - 1; i >= 0; i--) {
      debug9(printf("indel_univdiagonal %u, adj %d\n",
		    univdiags_alloc[i]->univdiagonal,univdiagonal - univdiags_alloc[i]->univdiagonal));
      /* Deletion */
      if (univdiags_alloc[i]->univdiagonal >= univdiagonal) {
	/* Skip */
      } else if ((adj = univdiagonal - univdiags_alloc[i]->univdiagonal) >= matchlength) {
	debug9(printf("adj >= matchlength\n"));
      } else if (best_adj == 0 || adj < best_adj) {
	best_univdiag = univdiags_alloc[i];
	best_adj = adj;
      }
    }
    if (best_univdiag != NULL) {
      debug9(printf("Returning deletion %u %d..%d\n",
		    best_univdiag->unidiagonal,best_univdiag->qstart,best_univdiag->qend));
      return best_univdiag;
    }
  }
      

  /* Test for insertion */
  insertion_univdiagonal = add_bounded(univdiagonal,max_insertionlen,chrhigh);
  debug9(printf("Add %d (bounded by %u) to yield (high) insertion_univdiagonal %u\n",
		max_insertionlen,chrhigh,insertion_univdiagonal));

  debug9(printf("Testing indel_qpos %d for insertion\n",indel_qpos));
  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				localdb,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,
				queryptr,/*pos5*/0,/*pos3*/indel_qpos,querylength,
				/*anchor_univdiagonal*/univdiagonal,
				/*low_univdiagonal*/univdiagonal,/*high_univdiagonal*/insertion_univdiagonal,
				query_compress,plusp,genestrand,genomebits,localdb_nmismatches_allowed,
				/*extend5p*/true,/*trim5p*/true,/*trim3p*/true)) == 0) {
    /* Skip */
    
  } else if ((nmatches = matchlength - local_nmismatches) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*local_nmismatches) {
    debug9(printf("(3) Got %d localdb diagonals with matchlength %d and local_nmismatches %d => Skipping\n",
		  ndiagonals,matchlength,local_nmismatches));
    /* Skip */
    
  } else {
    debug9(printf("(4) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		  ndiagonals,matchlength,local_nmismatches));
    
    best_adj = 0;
    for (i = ndiagonals - 1; i >= 0; i--) {
      debug9(printf("indel_univdiagonal %u, adj %d\n",
		    univdiags_alloc[i]->univdiagonal,univdiagonal - univdiags_alloc[i]->univdiagonal));
      /* Insertion */
      if (univdiags_alloc[i]->univdiagonal <= univdiagonal) {
	/* Skip */
      } else if ((adj = univdiags_alloc[i]->univdiagonal - univdiagonal) >= matchlength) {
	debug9(printf("adj >= matchlength\n"));
      } else if (best_adj == 0 || adj < best_adj) {
	best_univdiag = univdiags_alloc[i];
	best_adj = adj;
      }
    }
    if (best_univdiag != NULL) {
      debug9(printf("Returning insertion %u %d..%d\n",
		    best_univdiag->unidiagonal,best_univdiag->qstart,best_univdiag->qend));
      return best_univdiag;
    }
  }

  return (Univdiag_T) NULL;
}


static void
reverse_int_inplace (int *values, int starti, int endi) {
  int temp;
  int i, j, n = endi - starti;

  for (i = starti, j = endi-1; i < starti + n/2; i++, j--) {
    temp = values[i];
    values[i] = values[j];
    values[j] = temp;
  }

#if 0
  if (i == j) {
    values[i] = values[j];
  }
#endif

  return;
}


static void
reverse_double_inplace (double *values, int starti, int endi) {
  double temp;
  int i, j, n = endi - starti;

  for (i = starti, j = endi-1; i < starti + n/2; i++, j--) {
    temp = values[i];
    values[i] = values[j];
    values[j] = temp;
  }

#if 0
  if (i == j) {
    values[i] = values[j];
  }
#endif

  return;
}


static void
reverse_univcoord_inplace (Univcoord_T *coords, int starti, int endi) {
  Univcoord_T temp;
  int i, j, n = endi - starti;;

  for (i = starti, j = endi-1; i < starti + n/2; i++, j--) {
    temp = coords[i];
    coords[i] = coords[j];
    coords[j] = temp;
  }

#if 0
  if (i == j) {
    coords[i] = coords[j];
  }
#endif

  return;
}


static void
Spliceends_reverse (T this, int starti, int endi) {

  if (endi > starti) {
    reverse_int_inplace(this->splice_qpos,starti,endi);
    reverse_int_inplace(this->medial_nmismatches,starti,endi);
    reverse_univcoord_inplace(this->partners,starti,endi);

    /* reverse_int_inplace(this->matchlengths,starti,endi); */
    /* reverse_int_inplace(this->distal_lengths,starti,endi); */
    /* reverse_int_inplace(this->distal_trimpos,starti,endi); */
    /* reverse_int_inplace(this->distal_nmismatches,starti,endi); */

    Univdiag_reverse_inplace(this->distal_univdiags,starti,endi);

    reverse_double_inplace(this->medial_probs,starti,endi);
    reverse_double_inplace(this->distal_probs,starti,endi);

    reverse_int_inplace(this->medial_support,starti,endi);
    reverse_int_inplace(this->distal_support,starti,endi);
  }

  return;
}


static void
append_int (int *dest, int *values, int n) {
  int i;

  for (i = 0; i < n; i++) {
    *dest++ = *values++;
  }

  return;
}

static void
append_double (double *dest, double *values, int n) {
  int i;

  for (i = 0; i < n; i++) {
    *dest++ = *values++;
  }

  return;
}

static void
append_univcoord (Univcoord_T *dest, Univcoord_T *values, int n) {
  int i;

  for (i = 0; i < n; i++) {
    *dest++ = *values++;
  }

  return;
}

void
append_univdiags (Univdiag_T *dest, Univdiag_T *values, int n) {
  int i;

  for (i = 0; i < n; i++) {
    Univdiag_transfer(*dest++,*values++);
  }

  return;
}

static void
Spliceends_combine (T this_outward, int k_outward, T this_inward, int k_inward) {

  Spliceends_reverse(this_outward,/*starti*/0,/*endi*/k_outward);
  
  append_int(/*dest*/&(this_outward->splice_qpos[k_outward]),this_inward->splice_qpos,k_inward);
  append_int(/*dest*/&(this_outward->medial_nmismatches[k_outward]),this_inward->medial_nmismatches,k_inward);
  append_univcoord(/*dest*/&(this_outward->partners[k_outward]),this_inward->partners,k_inward);

#if 0
  append_int(/*dest*/&(this_outward->matchlengths[k_outward]),this_inward->matchlengths,k_inward);
  append_int(/*dest*/&(this_outward->distal_lengths[k_outward]),this_inward->distal_lengths,k_inward);
  append_int(/*dest*/&(this_outward->distal_trimpos[k_outward]),this_inward->distal_trimpos,k_inward);
  append_int(/*dest*/&(this_outward->distal_nmismatches[k_outward]),this_inward->distal_nmismatches,k_inward);
#endif

  append_univdiags(/*dest*/&(this_outward->distal_univdiags[k_outward]),this_inward->distal_univdiags,k_inward);

  append_double(/*dest*/&(this_outward->medial_probs[k_outward]),this_inward->medial_probs,k_inward);
  append_double(/*dest*/&(this_outward->distal_probs[k_outward]),this_inward->distal_probs,k_inward);

  append_int(/*dest*/&(this_outward->medial_support[k_outward]),this_inward->medial_support,k_inward);
  append_int(/*dest*/&(this_outward->distal_support[k_outward]),this_inward->distal_support,k_inward);

  return;
}


/* allocp enters with a value for novel_diagonals */
/* Always alloc, so we can use univdiags_alloc to find distal indels */
static int
merge_known_novel (Univdiag_T **univdiags, bool **knownp,
		   Univdiag_T *known_univdiags, int known_nunivdiags,
		   Univdiag_T *novel_univdiags, int novel_nunivdiags) {
  int i, j, k;


  if (known_nunivdiags == 0 && novel_nunivdiags == 0) {
    *univdiags = (Univdiag_T *) NULL;
    *knownp = (bool *) NULL;
    return 0;

#if 0
  } else if (known_nunivdiags == 0) {
    *univdiags = novel_univdiags;
    *knownp = (bool *) CALLOC(novel_nunivdiags,sizeof(bool));
    /* Use the value of allocp upon entry (novel_univdiags via indexdb or localdb) */
    return novel_nunivdiags;

  } else if (novel_nunivdiags == 0) {
    *univdiags = known_univdiags;
    *knownp = (bool *) MALLOC(known_nunivdiags * sizeof(bool));
    for (k = 0; k < known_nunivdiags; k++) {
      (*knownp)[k] = true;
    }
    *allocp = false;
    return known_nunivdiags;
#endif

  } else {
    i = j = 0;
    k = 0;

    /* *univdiags = (Univdiag_T *) MALLOC((known_nunivdiags + novel_nunivdiags) * sizeof(Univdiag_T)); */
    *univdiags = Univdiag_new_array(known_nunivdiags + novel_nunivdiags);
    *knownp = (bool *) MALLOC((known_nunivdiags + novel_nunivdiags) * sizeof(bool));

    while (i < known_nunivdiags && j < novel_nunivdiags) {
      if (known_univdiags[i]->univdiagonal < novel_univdiags[j]->univdiagonal) {
	Univdiag_transfer((*univdiags)[k],known_univdiags[i]);
	(*knownp)[k++] = true;
	i++;
      } else if (novel_univdiags[j]->univdiagonal < known_univdiags[i]->univdiagonal) {
	Univdiag_transfer((*univdiags)[k],novel_univdiags[j]);
	(*knownp)[k++] = false;
	j++;
      } else {
	Univdiag_transfer((*univdiags)[k],known_univdiags[i]);
	(*knownp)[k++] = true;
	i++; j++;
      }
    }

    while (i < known_nunivdiags) {
      Univdiag_transfer((*univdiags)[k],known_univdiags[i]);
      (*knownp)[k++] = true;
      i++;
    }

    while (j < novel_nunivdiags) {
      Univdiag_transfer((*univdiags)[k],novel_univdiags[j]);
      (*knownp)[k++] = false;
      j++;
    }

#if 0
    if (*allocp == true) {
      Univdiag_free_array(&novel_univdiags);
    }

    *allocp = true;
#endif
  }

  return k;
}


static int
solve_trim5 (int *min_nmismatches, double *max_prob, double *max_medial_prob, bool *max_medial_setp,
	     bool *partnerp, bool *boundedp, int *nsites, T this, int k, bool *sitep,
	     double (*Medial_prob_fcn)(Univcoord_T,Univcoord_T),
	     double (*Distal_prob_fcn)(Univcoord_T,Univcoord_T),

	     Univcoord_T genomicpos, int pos5, int splice_qpos, int qend, int querylength, int mismatchi,
	     Univcoord_T univdiagonal, Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,

	     Univdiag_T *known_univdiags_alloc, int known_nunivdiags,
			
	     Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
	     Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
	     Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
	     Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
	     bool plusp, int genestrand, int localdb_nmismatches_allowed,
	     bool search_localdb_p, bool innerp, int max_npartners, bool inwardp,
	     bool medial_salvagep) {

  int kstart = k;
  Univdiag_T *univdiags, *novel_univdiags, indel_univdiag;
  int nunivdiags, novel_nunivdiags = 0;
  bool *knownp, sortedp, trimmedp;
  bool allocp = false; 		/* For extend_trim5 */
  Chrpos_T positive_gap_distance;

  Univcoord_T distal_genomicpos;
  int local_nmismatches, matchlength, trimpos;
  double medial_prob, max_distal_prob, distal_prob;
  int nmismatches;
  int npartners, best_partneri, i;

  Univcoord_T left = univdiagonal - querylength;

#ifdef DEBUG1
  medial_prob = Medial_prob_fcn(genomicpos,chroffset);
  printf("5' %s: %u %u %d site:%d %f medial\n",
	 plusp ? "plus" : "minus",genomicpos,genomicpos-chroffset,genomicpos-left,
	 sitep[genomicpos - left],medial_prob);
#endif

  if (known_nunivdiags > 0) {
    medial_prob = 1.0;
    *max_medial_setp = true;	/* was not set before */
  } else if (sitep[genomicpos-left] == false) {
    medial_prob = 0.0;
    *max_medial_setp = false;
  } else if ((medial_prob = Medial_prob_fcn(genomicpos,chroffset)) > *max_medial_prob) {
    *max_medial_prob = medial_prob;
    *max_medial_setp = true;
  } else {
    *max_medial_setp = false;
  }

  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  if (medial_prob >= medial_splicesite_prob || medial_salvagep == true) {
    debug1(printf("Testing splice_qpos %d with mismatchi %d\n",splice_qpos,mismatchi));
    if (0 && innerp == false && splice_qpos < 4) {
      /* Forces program to find a bad splicing solution */
      /* Skip */

    } else if (stage1 != NULL &&
	       (novel_nunivdiags = extend_trim5(&novel_univdiags,splice_qpos,stage1,
						univdiagonal,querylength,plusp,
						/*slop*/positive_gap_distance,
						/*insertion_slop*/0,univdiagpool)) > 0) {
      allocp = true;
      debug1(printf("(1) Got %d indexdb univdiags\n",novel_nunivdiags));
      /* No need to sort.  Iterate from end */
	
    } else if (search_localdb_p == true &&
	       (novel_nunivdiags = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
					       localdb,diagonals_alloc,localdb_alloc,
					       stage1,streamspace_max_alloc,streamspace_alloc,
					       queryptr,pos5,/*pos3*/splice_qpos,querylength,
					       /*anchor_univdiagonal*/univdiagonal,
					       low_univdiagonal,high_univdiagonal,
					       query_compress,plusp,genestrand,genomebits,
					       localdb_nmismatches_allowed,/*extend5p*/true,
					       /*trim5p*/true,/*trim3p*/false)) > 0) {
      debug1(printf("(5) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		    novel_nunivdiags,matchlength,local_nmismatches));
      if (novel_nunivdiags == 1) {
	/* No need to sort */
      } else if (sortedp == true) {
	/* No need to sort.  Iterate from end */
      } else {
	qsort(univdiags_alloc[0],novel_nunivdiags,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
      }
      novel_univdiags = univdiags_alloc;
    }

    if ((nunivdiags = merge_known_novel(&univdiags,&knownp,
					known_univdiags_alloc,known_nunivdiags,
					novel_univdiags,novel_nunivdiags)) == 0) {
      /* Mark as a trimming position without a partner */
      debug1(printf("Got no novel univdiags\n"));
      this->splice_qpos[k] = splice_qpos;
      this->medial_nmismatches[k] = mismatchi;
      this->partners[k] = 0;

      this->medial_probs[k] = medial_prob;
      this->distal_probs[k] = 0.0;
      this->medial_support[k] = (qend - splice_qpos);  /*- 4*mismatchi for adj_support */
      this->distal_support[k] = 0;
      k++;

      (*nsites)++;
	  
    } else {
      /* Process starting from the closest diagonal */
      best_partneri = -1; max_distal_prob = 0.0;
      for (i = nunivdiags - 1; i >= 0; i--) {
	distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos;
	if (knownp[i] == true) {
	  distal_prob = 1.0;
	} else {
	  distal_prob = Distal_prob_fcn(distal_genomicpos,chroffset);
	}
	/* debug1(printf("Evaluating diagonal %u, distal_genomicpos %u, prob %f\n",
	   univdiags[i]->univdiagonal,distal_genomicpos,distal_prob)); */

	if (distal_prob > max_distal_prob) {
	  max_distal_prob = distal_prob;
	  best_partneri = i;
	}
      }

      npartners = 0;
      if (max_distal_prob >= distal_splicesite_prob || medial_salvagep == true) {
	i = nunivdiags - 1;
	while (npartners < max_npartners && i >= 0) {
	  distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos;
	  if (knownp[i] == true) {
	    distal_prob = 1.0;
	  } else {
	    distal_prob = Distal_prob_fcn(distal_genomicpos,chroffset);
	  }

	  if (distal_prob >= max_distal_prob - PROB_SLOP) {

	    debug5(printf("Univdiag %u, %d..%d\n",univdiags[i]->univdiagonal,univdiags[i]->qstart,univdiags[i]->qend));
	    if ((trimpos = Genomebits_trim_qstart(&(univdiags[i]->nmismatches),query_compress,
						  genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
						  pos5,/*pos3*/splice_qpos,plusp,genestrand)) <= univdiags[i]->qstart) {
	      /* Extends to splice_qpos */
	      debug5(printf("(1) At univdiagonal %u, trimpos %d <= qstart %d (preserved), so extends to splice_qpos %d\n",
			    univdiags[i]->univdiagonal,trimpos,univdiags[i]->qstart,splice_qpos));
	      univdiags[i]->qend = splice_qpos;
	      univdiags[i]->qstart = trimpos;
	    } else {
	      /* Cassette mismatch, perhaps due to a spliceindel */
	      debug5(printf("(1) At univdiagonal %u, trimpos %d > qstart %d, so cassette mismatch.  ",
			    univdiags[i]->univdiagonal,trimpos,univdiags[i]->qstart));
	      trimpos = univdiags[i]->qstart =
		Genomebits_trim_qstart(&(univdiags[i]->nmismatches),query_compress,
				       genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
				       pos5,/*pos3*/univdiags[i]->qend,plusp,genestrand);
	      debug5(printf("New bounds are %d..%d\n",univdiags[i]->qstart,univdiags[i]->qend));
	    }
		 
	    debug1(printf("For univdiagonal %u from %d to %d, distal trimpos is %d\n",
			  univdiags[i]->univdiagonal,pos5,univdiags[i]->qend,trimpos));
	  
	    if (/*univdiags[i]->qstart*/trimpos == univdiags[i]->qend) {
	      /* Skip: Not a good partner */
	      debug1(printf("trimpos %d == qend %d (length zero), so not a good partner\n",
			    trimpos,univdiags[i]->qend));
	      
	    } else if (splice_qpos - pos5 < 8 || univdiags[i]->qend - univdiags[i]->qstart >= 8) {
	      /* Adequate extension */
	      this->splice_qpos[k] = splice_qpos;
	      this->medial_nmismatches[k] = mismatchi;
	      this->partners[k] = distal_genomicpos;

	      Univdiag_transfer(this->distal_univdiags[k],univdiags[i]);

	      this->medial_probs[k] = medial_prob;
	      this->distal_probs[k] = distal_prob;
	      this->medial_support[k] = (qend - splice_qpos); /*- 4*mismatchi for adj_support*/
	      this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart);
	      /*- 4*univdiags[i]->nmismatches for adj_support*/
	    
	      if ((nmismatches = mismatchi + univdiags[i]->nmismatches + trimpos) < *min_nmismatches) {
		*min_nmismatches = nmismatches;
		*max_prob = medial_prob + distal_prob;
	      } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
		*max_prob = medial_prob + distal_prob;
	      }
	      debug1(printf("(1) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			  univdiags[i]->univdiagonal,this->partners[k],univdiagonal - univdiags[i]->univdiagonal,
			  this->medial_probs[k],this->distal_probs[k],
			  this->splice_qpos[k],this->medial_nmismatches[k],
			  univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	      k++;
	  
	      npartners++;
	      *partnerp = true;

	    } else if (localdb == NULL) {
	      /* Skip */

	    } else if ((indel_univdiag =
			novel_trim5_indel(/*genomicpos*/univdiags[i]->univdiagonal - querylength + univdiags[i]->qstart,
					  univdiags[i]->univdiagonal,querylength,
					  query_compress,queryptr,chroffset,chrhigh,
					  univdiags_alloc,diagonals_alloc,localdb_alloc,
					  stage1,streamspace_max_alloc,streamspace_alloc,
					  plusp,genestrand,localdb_nmismatches_allowed)) == (Univdiag_T) NULL) {
	      /* Skip: Non-extension without distal indel */
	      debug1(printf("non-extension without distal indel\n"));

	    } else {
	      /* Provide indel univdiag instead so we find a spliceindel */
	      debug1(printf("found indel %u %d..%d\n",
			    indel_univdiag->univdiagonal,indel_univdiag->qstart,indel_univdiag->qend));

	      this->splice_qpos[k] = splice_qpos;
	      this->medial_nmismatches[k] = mismatchi;
	      this->partners[k] = distal_genomicpos;

	      Univdiag_transfer(this->distal_univdiags[k],indel_univdiag);

	      this->medial_probs[k] = medial_prob;
	      this->distal_probs[k] = distal_prob;
	      this->medial_support[k] = (qend - splice_qpos); /*- 4*mismatchi for adj_support*/
	      this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart) +
		(indel_univdiag->qend - indel_univdiag->qstart);
	      /*- 4*univdiags[i]->nmismatches for adj_support*/
	      if (this->distal_support[k] > splice_qpos) {
		this->distal_support[k] = splice_qpos;
	      }
	    
	      if ((nmismatches = mismatchi + univdiags[i]->nmismatches + trimpos) < *min_nmismatches) {
		*min_nmismatches = nmismatches;
		*max_prob = medial_prob + distal_prob;
	      } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
		*max_prob = medial_prob + distal_prob;
	      }
	      debug1(printf("(1) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			    this->distal_univdiags[k]->univdiagonal,this->partners[k],
			    univdiagonal - this->distal_univdiags[k]->univdiagonal,
			    this->medial_probs[k],this->distal_probs[k],
			    this->splice_qpos[k],this->medial_nmismatches[k],
			    this->distal_univdiags[k]->nmismatches,
			    this->distal_univdiags[k]->qstart,this->distal_univdiags[k]->qend));
	      k++;
	  
	      npartners++;
	      *partnerp = true;
	    }
	  }
	  
	  i--;
	}

	if (i >= 0) {
	  /* Limited by max_npartners */
	  debug1(printf("BOUNDED\n"));
	  *boundedp = true;

	} else if (medial_salvagep == true) {
	  /* Do not salvage both ends */
	  
	} else if (npartners == 0 && max_distal_prob > 0.5) {
	  /* Salvage: use best partneri */
	  i = best_partneri;
	  distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos; /* retrieve */
	  
	  debug5(printf("Univdiag %u, %d..%d\n",univdiags[i]->univdiagonal,univdiags[i]->qstart,univdiags[i]->qend));
	  if ((trimpos = Genomebits_trim_qstart(&(univdiags[i]->nmismatches),query_compress,
						genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
						pos5,/*pos3*/splice_qpos,plusp,genestrand)) <= univdiags[i]->qstart) {
	    /* Extends to splice_qpos */
	    debug5(printf("(2) At univdiagonal %u, trimpos %d <= qstart %d (preserved), so extends to splice_qpos %d\n",
			  univdiags[i]->univdiagonal,trimpos,univdiags[i]->qstart,splice_qpos));
	    univdiags[i]->qend = splice_qpos;
	    univdiags[i]->qstart = trimpos;
	  } else {
	    /* Cassette mismatch, perhaps due to a spliceindel */
	    debug5(printf("(2) At univdiagonal %u, trimpos %d > qstart %d, so cassette mismatch.  ",
			  univdiags[i]->univdiagonal,trimpos,univdiags[i]->qstart));
	    trimpos = univdiags[i]->qstart =
	      Genomebits_trim_qstart(&(univdiags[i]->nmismatches),query_compress,
				     genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
				     pos5,/*pos3*/univdiags[i]->qend,plusp,genestrand);
	    debug5(printf("New bounds are %d..%d\n",univdiags[i]->qstart,univdiags[i]->qend));
	  }

	  debug1(printf("For univdiagonal %u from %d to %d, distal trimpos is %d\n",
			univdiags[i]->univdiagonal,pos5,univdiags[i]->qend,univdiags[i]->qstart));

	  if (/*univdiags[i]->qstart*/trimpos == univdiags[i]->qend) {
	    /* Skip: Still not a good partner */
	    
	  } else if ((nmismatches = mismatchi + univdiags[i]->nmismatches + trimpos) > 1) {
	    /* Have a stricter standard for salvage cases */
	    
	  } else if (splice_qpos - pos5 < 8 || univdiags[i]->qend - univdiags[i]->qstart >= 8) {
	    /* Adequate extension */
	    this->splice_qpos[k] = splice_qpos;
	    this->medial_nmismatches[k] = mismatchi;
	    this->partners[k] = distal_genomicpos;

	    Univdiag_transfer(this->distal_univdiags[k],univdiags[i]);

	    this->medial_probs[k] = medial_prob;
	    this->distal_probs[k] = max_distal_prob; /* retrieve */
	    this->medial_support[k] = (qend - splice_qpos); /*- 4*mismatchi for adj_support*/
	    this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart);
	    /*- 4*univdiags[i]->nmismatches for adj_support*/
	  
	    if (nmismatches < *min_nmismatches) {
	      *min_nmismatches = nmismatches;
	      *max_prob = medial_prob + distal_prob;
	    } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
	      *max_prob = medial_prob + distal_prob;
	    }
	    debug1(printf("(2) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d (salvage)\n",
			  univdiags[i]->univdiagonal,this->partners[k],univdiagonal - univdiags[i]->univdiagonal,
			  this->medial_probs[k],this->distal_probs[k],
			  this->splice_qpos[k],this->medial_nmismatches[k],
			  univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	    k++;
	  
	    /* npartners++; */
	    *partnerp = true;

	  } else if (localdb == NULL) {
	    /* Skip */

	  } else if ((indel_univdiag =
		      novel_trim5_indel(/*genomicpos*/univdiags[i]->univdiagonal - querylength + univdiags[i]->qstart,
					univdiags[i]->univdiagonal,querylength,
					query_compress,queryptr,chroffset,chrhigh,
					univdiags_alloc,diagonals_alloc,localdb_alloc,
					stage1,streamspace_max_alloc,streamspace_alloc,
					plusp,genestrand,localdb_nmismatches_allowed)) == (Univdiag_T) NULL) {
	    /* Skip: Non-extension without distal indel */
	    debug1(printf("non-extension without distal indel\n"));

	  } else {
	    /* Provide indel univdiag instead so we find a spliceindel */
	    debug1(printf("found indel %u %d..%d\n",
			  indel_univdiag->univdiagonal,indel_univdiag->qstart,indel_univdiag->qend));

	    this->splice_qpos[k] = splice_qpos;
	    this->medial_nmismatches[k] = mismatchi;
	    this->partners[k] = distal_genomicpos;

	    Univdiag_transfer(this->distal_univdiags[k],indel_univdiag);

	    this->medial_probs[k] = medial_prob;
	    this->distal_probs[k] = max_distal_prob; /* retrieve */
	    this->medial_support[k] = (qend - splice_qpos); /*- 4*mismatchi for adj_support*/
	    this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart) +
	      (indel_univdiag->qend - indel_univdiag->qstart);
	    /*- 4*univdiags[i]->nmismatches for adj_support*/
	    if (this->distal_support[k] > splice_qpos) {
	      this->distal_support[k] = splice_qpos;
	    }
	  
	    if (nmismatches < *min_nmismatches) {
	      *min_nmismatches = nmismatches;
	      *max_prob = medial_prob + distal_prob;
	    } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
	      *max_prob = medial_prob + distal_prob;
	    }
	    debug1(printf("(2) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d (salvage)\n",
			  this->distal_univdiags[k]->univdiagonal,this->partners[k],
			  univdiagonal - this->distal_univdiags[k]->univdiagonal,
			  this->medial_probs[k],this->distal_probs[k],
			  this->splice_qpos[k],this->medial_nmismatches[k],
			  this->distal_univdiags[k]->nmismatches,
			  this->distal_univdiags[k]->qstart,this->distal_univdiags[k]->qend));
	    k++;
	  
	    /* npartners++; */
	    *partnerp = true;

	  }
	}
      }

      if (inwardp == false) {
	Spliceends_reverse(this,kstart,k);
      }

      Univdiag_free_array(&univdiags);
      FREE(knownp);
      (*nsites)++;
    }

    if (allocp == true) {
      Univdiag_free_array(&novel_univdiags);
    }
  }

  return k;
}


static int
spliceends_trim5_sense (double *max_prob, bool *partnerp, T this_outward, T this_inward, int qend,
			Univcoord_T inbound_genomicpos, Univcoord_T outbound_genomicpos, Univcoord_T middle_genomicpos,
			int *mismatch_positions, int total_nmismatches,
			Univcoord_T univdiagonal, int querylength,
			Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
			Knownsplicing_T knownsplicing, int max_nconsecutive,
			Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
			bool plusp, int genestrand, int localdb_nmismatches_allowed,
			bool search_localdb_p, bool innerp) {
  int nspliceends = 0, k, k_inward = 0, k_outward = 0;
  Chrpos_T splice_chrbound, positive_gap_distance;

  Univcoord_T *endpoints_inward, *endpoints_outward;
  uint64_t low_rank_inward = 0, high_rank_inward = 0, rank_inward = 0,
    low_rank_outward = 0, high_rank_outward = 0, rank_outward = 0;
  int n_inward, n_outward;

  Univdiag_T *known_univdiags_alloc = NULL, univdiag;
  Univcoord_T best_genomicpos, genomicpos_inward, genomicpos_outward, low_univdiagonal, high_univdiagonal;
  int known_nunivdiags;

  int pos5, best_splice_qpos, splice_qpos_inward, splice_qpos_outward, nconsecutive;
  int best_mismatchi, mismatchi_inward, mismatchi_outward;
  int nsites = 0, i;
  Univcoord_T left = univdiagonal - querylength;
  double max_medial_prob = 0.0;
  int min_nmismatches = querylength;
  bool boundedp = false, max_medial_setp, best_inwardp;

  int max_npartners = (innerp == true) ? MAX_INNER_PARTNERS_PER_SITE : MAX_OUTER_PARTNERS_PER_SITE;
  bool *sitep;

  *max_prob = 0.0;
  *partnerp = false;

  assert(outbound_genomicpos < inbound_genomicpos);

  /* Prevent splicing past chromosome bounds */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  if (univdiagonal < chroffset + positive_gap_distance) {
    splice_chrbound = univdiagonal - chroffset;
  } else {
    splice_chrbound = positive_gap_distance;
  }

  debug1(printf("Start with univdiagonal %u\n",univdiagonal));
  high_univdiagonal = univdiagonal - MIN_INTRON_LENGTH; /* Subtract so we don't return a continuation */
  debug1(printf("Subtract %d to yield high_univdiagonal %u\n",MIN_INTRON_LENGTH,high_univdiagonal));
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos5 = (high_univdiagonal >= chroffset + (Univcoord_T) querylength) ? 0 : (int) (chroffset - left);
#else
  pos5 = 0;
#endif
  low_univdiagonal = subtract_bounded(high_univdiagonal,splice_chrbound,chroffset);
  debug1(printf("Subtract %d (bounded by %u) to yield low_univdiagonal %u\n",
		splice_chrbound,chroffset,low_univdiagonal));
  debug1(printf("  (1) splice_chrbound %u, low %u, high %u\n",
		splice_chrbound,low_univdiagonal,high_univdiagonal));

#ifdef DEBUG1
  printf("(1) mismatch positions:");
  for (i = 0; i <= total_nmismatches; i++) {
    printf(" %d",mismatch_positions[i]);
  }
  printf("\n");
#endif


  if (low_univdiagonal > high_univdiagonal) {
    /* Skip */

  } else if (plusp) {
    debug1(printf("spliceends_trim5_sense, plus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));
    
    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_acceptors(&low_rank_inward,&high_rank_inward,
						 knownsplicing,univdiagonal,querylength,
						 middle_genomicpos - left,inbound_genomicpos - left);
      rank_inward = low_rank_inward;
#ifdef DEBUG1
      printf("(1) Knownsplicing acceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos,inbound_genomicpos,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank] - endpoints_inward[2*rank+1]));
      }
#endif

      endpoints_outward = Knownsplicing_acceptors(&low_rank_outward,&high_rank_outward,
						  knownsplicing,univdiagonal,querylength,
						  outbound_genomicpos - left,middle_genomicpos - left);
      rank_outward = high_rank_outward;
#ifdef DEBUG1
      printf("(1) Knownsplicing acceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     outbound_genomicpos,middle_genomicpos,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank] - endpoints_outward[2*rank+1]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_acceptor_ag_sites(sitep,left,/*pos5*/outbound_genomicpos - left,/*pos3*/inbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] < splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos - 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] < splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward < inbound_genomicpos || genomicpos_outward >= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward < inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward < high_rank_inward && endpoints_inward[2*rank_inward] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*rank_inward + 1] + querylength - splice_qpos_inward;
	  univdiag->qend = splice_qpos_inward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }

	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward++;
	}
	
	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_donor_prob,

			       genomicpos_inward,pos5,splice_qpos_inward,qend,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
			       
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);
      
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
	splice_qpos_inward++;
	genomicpos_inward++;
      }

      if (genomicpos_outward >= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward > low_rank_outward && endpoints_outward[2*(rank_outward - 1)] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*(rank_outward - 1) + 1] + querylength - splice_qpos_outward;
	  univdiag->qend = splice_qpos_outward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }
	  
	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward--;
	}

	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_donor_prob,

				genomicpos_outward,pos5,splice_qpos_outward,qend,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	splice_qpos_outward--;
	genomicpos_outward--;
	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_donor_prob,

			       best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,/*known_nunivdiags*/0,
			       
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_donor_prob,

				best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,/*known_nunivdiags*/0,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);

  } else {
    debug1(printf("spliceends_trim5_sense, minus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_antidonors(&low_rank_inward,&high_rank_inward,
						  knownsplicing,univdiagonal,querylength,
						  middle_genomicpos - left,inbound_genomicpos - left);
      rank_inward = low_rank_inward;
#ifdef DEBUG1
      printf("(2) Knownsplicing antidonors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos,inbound_genomicpos,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank] - endpoints_inward[2*rank+1]));
      }
#endif

      endpoints_outward = Knownsplicing_antidonors(&low_rank_outward,&high_rank_outward,
						   knownsplicing,univdiagonal,querylength,
						   outbound_genomicpos - left,middle_genomicpos - left);
      rank_outward = high_rank_outward;
#ifdef DEBUG1
      printf("(2) Knownsplicing antidonors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     outbound_genomicpos,inbound_genomicpos,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank] - endpoints_outward[2*rank+1]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }
  
    Genome_fill_antidonor_gtgc_sites(sitep,left,/*pos5*/outbound_genomicpos - left,/*pos3*/inbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] < splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos - 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] < splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward < inbound_genomicpos || genomicpos_outward >= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward < inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward < high_rank_inward && endpoints_inward[2*rank_inward] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*rank_inward + 1] + querylength - splice_qpos_inward;
	  univdiag->qend = splice_qpos_inward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }
	  
	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward++;
	}

	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

			       genomicpos_inward,pos5,splice_qpos_inward,qend,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
	splice_qpos_inward++;
	genomicpos_inward++;
      }

      if (genomicpos_outward >= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward > low_rank_outward &&
	       endpoints_outward[2*(rank_outward - 1)] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*(rank_outward - 1) + 1] + querylength - splice_qpos_outward;
	  univdiag->qend = splice_qpos_outward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }
	  
	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward--;
	}

	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

				genomicpos_outward,pos5,splice_qpos_outward,qend,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	splice_qpos_outward--;
	genomicpos_outward--;
	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

			       best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

				best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);
  }

  debug1(printf("\n"));

  /* Combine inward and outward */
  debug1(printf("k_inward = %d, k_outward = %d, partnerp %d, max_prob %f, max_medial_prob %f\n",
		k_inward,k_outward,*partnerp,*max_prob,max_medial_prob));
  Spliceends_combine(this_outward,k_outward,this_inward,k_inward);
  this_outward->boundedp = boundedp;

  if ((k = k_inward + k_outward) == 0) {
    /* Skip */
    debug1(printf("k is zero\n"));
    nspliceends = 0;

  } else if (*partnerp == false) {
    /* min_nmismatches never set */
    debug1(printf("No partner, so use prob only\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_probs[i] > max_medial_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i]; /* Should be 0 */

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	/* No partner */
	/* Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]); */

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else if (innerp == false) {
    /* Keep only the best spliceends */
    debug1(printf("outer, so keeping only the best\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_nmismatches[i] + this_outward->distal_univdiags[i]->nmismatches + this_outward->distal_univdiags[i]->qstart <= min_nmismatches &&
	  this_outward->medial_probs[i] + this_outward->distal_probs[i] > *max_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i];*/

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else {
    /* Keep all partner ends for later resolving of ambiguity */
    debug1(printf("inner with partners, so keeping all\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->partners[i] != 0) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }
  }

#ifdef DEBUG1
  printf("spliceends_trim5_sense yielded %d spliceends\n",nspliceends);
  for (k = 0; k < nspliceends; k++) {
    printf("%u %u %d\n",this_outward->partners[k],this_outward->partners[k] - chroffset,this_outward->splice_qpos[k]);
  }
#endif

  return nspliceends;
}


static int
spliceends_trim5_anti (double *max_prob, bool *partnerp, T this_outward, T this_inward, int qend,
		       Univcoord_T inbound_genomicpos, Univcoord_T outbound_genomicpos, Univcoord_T middle_genomicpos,
		       int *mismatch_positions, int total_nmismatches,
		       Univcoord_T univdiagonal, int querylength,
		       Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
		       Knownsplicing_T knownsplicing, int max_nconsecutive,
		       Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		       Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
		       bool plusp, int genestrand, int localdb_nmismatches_allowed,
		       bool search_localdb_p, bool innerp) {
  int nspliceends = 0, k, k_inward = 0, k_outward = 0;
  Chrpos_T splice_chrbound, positive_gap_distance;

  Univcoord_T *endpoints_inward, *endpoints_outward;
  uint64_t low_rank_inward = 0, high_rank_inward = 0, rank_inward = 0,
    low_rank_outward = 0, high_rank_outward = 0, rank_outward = 0;
  int n_inward, n_outward;

  Univdiag_T *known_univdiags_alloc = NULL, univdiag;
  Univcoord_T best_genomicpos, genomicpos_inward, genomicpos_outward, low_univdiagonal, high_univdiagonal;
  int known_nunivdiags;

  int pos5, best_splice_qpos, splice_qpos_inward, splice_qpos_outward, nconsecutive;
  int best_mismatchi, mismatchi_inward, mismatchi_outward;
  int nsites = 0, i;
  Univcoord_T left = univdiagonal - querylength;
  double max_medial_prob = 0.0;
  int min_nmismatches = querylength;
  bool boundedp = false, max_medial_setp, best_inwardp;

  int max_npartners = (innerp == true) ? MAX_INNER_PARTNERS_PER_SITE : MAX_OUTER_PARTNERS_PER_SITE;
  bool *sitep;

  *max_prob = 0.0;
  *partnerp = false;

  assert(outbound_genomicpos < inbound_genomicpos);

  /* Prevent splicing past chromosome bounds */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  if (univdiagonal < chroffset + positive_gap_distance) {
    splice_chrbound = univdiagonal - chroffset;
  } else {
    splice_chrbound = positive_gap_distance;
  }

  debug1(printf("Start with univdiagonal %u\n",univdiagonal));
  high_univdiagonal = univdiagonal - MIN_INTRON_LENGTH; /* Subtract so we don't return a continuation */
  debug1(printf("Subtract %d to yield high_univdiagonal %u\n",MIN_INTRON_LENGTH,high_univdiagonal));
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos5 = (high_univdiagonal >= chroffset + (Univcoord_T) querylength) ? 0 : (int) (chroffset - left);
#else
  pos5 = 0;
#endif
  low_univdiagonal = subtract_bounded(high_univdiagonal,splice_chrbound,chroffset);
  debug1(printf("Subtract %d (bounded by %u) to yield low_univdiagonal %u\n",
		splice_chrbound,chroffset,low_univdiagonal));
  debug1(printf("  (2) splice_chrbound %u, low %u, high %u\n",
		splice_chrbound,low_univdiagonal,high_univdiagonal));

#ifdef DEBUG1
  printf("(2) mismatch positions:");
  for (i = 0; i <= total_nmismatches; i++) {
    printf(" %d",mismatch_positions[i]);
  }
  printf("\n");
#endif


  if (low_univdiagonal > high_univdiagonal) {
    /* Skip */

  } else if (plusp) {
    debug1(printf("spliceends_trim5_anti, plus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));
 
    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_antidonors(&low_rank_inward,&high_rank_inward,knownsplicing,univdiagonal,querylength,
					   middle_genomicpos - left,inbound_genomicpos - left);
      rank_inward = low_rank_inward;
#ifdef DEBUG1
      printf("(3) Knownsplicing antidonors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos,inbound_genomicpos,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank] - endpoints_inward[2*rank+1]));
      }
#endif

      endpoints_outward = Knownsplicing_antidonors(&low_rank_outward,&high_rank_outward,knownsplicing,univdiagonal,querylength,
						   outbound_genomicpos - left,middle_genomicpos - left);
      rank_outward = high_rank_outward;

#ifdef DEBUG1
      printf("(3) Knownsplicing antidonors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     outbound_genomicpos,middle_genomicpos,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank] - endpoints_outward[2*rank+1]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_antidonor_gtgc_sites(sitep,left,/*pos5*/outbound_genomicpos - left,/*pos3*/inbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] < splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos - 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] < splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward < inbound_genomicpos || genomicpos_outward >= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward < inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward < high_rank_inward && endpoints_inward[2*rank_inward] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*rank_inward + 1] + querylength - splice_qpos_inward;
	  univdiag->qend = splice_qpos_inward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }
	  
	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward++;
	}

	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

			       genomicpos_inward,pos5,splice_qpos_inward,qend,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
	splice_qpos_inward++;
	genomicpos_inward++;
      }

      if (genomicpos_outward >= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward > low_rank_outward && endpoints_outward[2*(rank_outward - 1)] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*(rank_outward - 1) + 1] + querylength - splice_qpos_outward;
	  univdiag->qend = splice_qpos_outward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }

	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward--;
	}

	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

				genomicpos_outward,pos5,splice_qpos_outward,qend,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	splice_qpos_outward--;
	genomicpos_outward--;
	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

			       best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antidonor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antiacceptor_prob,

				best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);

  } else {
    debug1(printf("spliceends_trim5_anti, minus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_acceptors(&low_rank_inward,&high_rank_inward,
						 knownsplicing,univdiagonal,querylength,
						 middle_genomicpos - left,inbound_genomicpos - left);
      rank_inward = low_rank_inward;

#ifdef DEBUG1
      printf("(4) Knownsplicing acceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos,inbound_genomicpos,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank] - endpoints_inward[2*rank+1]));
      }
#endif

      endpoints_outward = Knownsplicing_acceptors(&low_rank_outward,&high_rank_outward,knownsplicing,univdiagonal,querylength,
						  outbound_genomicpos - left,middle_genomicpos - left);
      rank_outward = high_rank_outward;

#ifdef DEBUG1
      printf("(4) Knownsplicing acceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     outbound_genomicpos,middle_genomicpos,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank] - endpoints_outward[2*rank+1]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_acceptor_ag_sites(sitep,left,/*pos5*/outbound_genomicpos - left,/*pos3*/inbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] < splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos - 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] < splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward < inbound_genomicpos || genomicpos_outward >= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward < inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward < high_rank_inward && endpoints_inward[2*rank_inward] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*rank_inward + 1] + querylength - splice_qpos_inward;
	  univdiag->qend = splice_qpos_inward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }

	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward++;
	}

	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_donor_prob,

			       genomicpos_inward,pos5,splice_qpos_inward,qend,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
	splice_qpos_inward++;
	genomicpos_inward++;
      }

      if (genomicpos_outward >= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward > low_rank_outward && endpoints_outward[2*(rank_outward - 1)] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*(rank_outward - 1) + 1] + querylength - splice_qpos_outward;
	  univdiag->qend = splice_qpos_outward;
	  univdiag->qstart = Genomebits_trim_qstart(&univdiag->nmismatches,query_compress,genomebits,
						    univdiag->univdiagonal,querylength,
						    /*pos5*/0,/*pos3*/univdiag->qend,plusp,genestrand);
	  if (univdiag->qstart < univdiag->qend) {
	    known_nunivdiags++;
	  }

	  debug1(printf("5' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward--;
	}

	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_donor_prob,
				
				genomicpos_outward,pos5,splice_qpos_outward,qend,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);
	
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	splice_qpos_outward--;
	genomicpos_outward--;
	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_donor_prob,

			       best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim5(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_acceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_donor_prob,
				
				best_genomicpos,pos5,best_splice_qpos,qend,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
				
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);
  }

  debug1(printf("\n"));

  /* Combine inward and outward */
  debug1(printf("k_inward = %d, k_outward = %d, partnerp %d, max_prob %f, max_medial_prob %f\n",
		k_inward,k_outward,*partnerp,*max_prob,max_medial_prob));
  Spliceends_combine(this_outward,k_outward,this_inward,k_inward);
  this_outward->boundedp = boundedp;

  if ((k = k_inward + k_outward) == 0) {
    /* Skip */
    debug1(printf("k is zero\n"));
    nspliceends = 0;

  } else if (*partnerp == false) {
    /* min_nmismatches never set */
    debug1(printf("No partner, so use prob only\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_probs[i] > max_medial_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	/* No partner */
	/* Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]); */

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else if (innerp == false) {
    /* Keep only the best spliceends */
    debug1(printf("outer, so keeping only the best\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_nmismatches[i] + this_outward->distal_univdiags[i]->nmismatches + this_outward->distal_univdiags[i]->qstart <= min_nmismatches &&
	  this_outward->medial_probs[i] + this_outward->distal_probs[i] > *max_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else {
    /* Keep all partner ends for later resolving of ambiguity */
    debug1(printf("inner with partners, so keeping all\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->partners[i] != 0) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }
  }

#ifdef DEBUG1
  printf("spliceends_trim5_anti yielded %d spliceends\n",nspliceends);
  for (k = 0; k < nspliceends; k++) {
    printf("%u %u %d\n",this_outward->partners[k],this_outward->partners[k] - chroffset,this_outward->splice_qpos[k]);
  }
#endif

  return nspliceends;
}



/* TODO: Consider whether there is a conflict between the regular prob
   and mismatch prob, and if so, set splicedir to be SENSE_NULL */
static T
trim_5 (bool *partnerp, Compress_T query_compress, char *queryptr,
	Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
	Knownsplicing_T knownsplicing, int try_sensedir,
	Univcoord_T univdiagonal, int querylength,
	int pos5, int qstart, int qend, int exon_origin, Univcoord_T chroffset, Univcoord_T chrhigh,
	int *mismatch_positions, int total_nmismatches,
	Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
	int max_nconsecutive, bool plusp, int genestrand, int localdb_nmismatches_allowed,
	bool search_localdb_p, bool innerp, bool salvagep) {
  T new, new_inward;
  int found_sensedir;
  int qdistal;
  int mismatchi;
  int nspliceends = 0;
  Univcoord_T outbound_genomicpos, middle_genomicpos, inbound_genomicpos;
  double max_sense_prob = 0.0, max_antisense_prob = 0.0;
  double medial_splicesite_prob, distal_splicesite_prob;
  Univcoord_T left = univdiagonal - querylength;


  debug1(printf("\nEntered trim_5 with try_sensedir %d, qstart %d to qend %d\n",
		 try_sensedir,qstart,qend));
  assert(try_sensedir == SENSE_FORWARD || try_sensedir == SENSE_ANTI);

  if (salvagep == true) {
    medial_splicesite_prob = SALVAGE_MEDIAL_SPLICESITE_PROB;
    distal_splicesite_prob = SALVAGE_DISTAL_SPLICESITE_PROB;
  } else {
    medial_splicesite_prob = DEFAULT_MEDIAL_SPLICESITE_PROB;
    distal_splicesite_prob = DEFAULT_DISTAL_SPLICESITE_PROB;
  }

  /* Search from qstart - END_SPLICESITE_SEARCH_MM to qstart +
     END_SPLICESITE_SEARCH, but not past qend */

  middle_genomicpos = left + qstart;

  debug1(
	 printf("%d mismatches from qend %d down to pos5 %d:",total_nmismatches,qend,pos5);
	 for (int i = 0; i <= total_nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );

  /* middle_genomicpos is where the current trim is.  Previously
     subtracted END_SPLICESITE_SEARCH_MM, but now allowing some number
     of mismatches distally to get to outbound_genomicpos */
  mismatchi = 0;
  while (mismatchi < total_nmismatches && mismatch_positions[mismatchi] >= qstart) {
    mismatchi++;
  }

  if (mismatchi + distal_nmismatches_allowed >= total_nmismatches) {
    qdistal = pos5;
    debug1(printf("qdistal = pos5 %d\n",qdistal));
  } else {
    qdistal = mismatch_positions[mismatchi + distal_nmismatches_allowed];
    debug1(printf("qdistal = mismatch position %d\n",qdistal));
  }
  if (qdistal == 0) {
    /* Don't need to evaluate at the end of the read */
    qdistal = 1;
  }

#if 0
  if (middle_genomicpos < left + END_SPLICESITE_SEARCH_MM) {
    outbound_genomicpos = left;
  } else {
    outbound_genomicpos = middle_genomicpos - END_SPLICESITE_SEARCH_MM;
  }
#else
  outbound_genomicpos = left + qdistal;
#endif

  if (middle_genomicpos < outbound_genomicpos) {
    middle_genomicpos = outbound_genomicpos;
  }

  if ((inbound_genomicpos = middle_genomicpos + END_SPLICESITE_SEARCH) > left + qend) {
    inbound_genomicpos = left + qend;
  }

  if (left + exon_origin < MIN_EXON_LENGTH) {
    /* At beginning of genome, so don't subtract MIN_EXON_LENGTH */
  } else if (inbound_genomicpos > left + exon_origin - MIN_EXON_LENGTH) {
    inbound_genomicpos = left + exon_origin - MIN_EXON_LENGTH;
  }

  debug1(printf("\n1 Set end points for 5' trim to be %u..%u..%u\n",
		outbound_genomicpos,middle_genomicpos,inbound_genomicpos));

  if (outbound_genomicpos >= inbound_genomicpos) {
    debug1(printf("Got no spliceends\n"));
    return (T) NULL;
  } else {
#if 0
    new = Spliceends_new(/*id*/0,querylength,vectorpool,spliceendspool);
    new_inward = Spliceends_new(/*id*/0,querylength,vectorpool,spliceendspool);
#else
    new = Spliceendsgen_checkout(spliceendsgen,querylength,vectorpool);
    new_inward = Spliceendsgen_checkout(spliceendsgen,querylength,vectorpool);
#endif
  }

  if (stage1 != NULL && stage1->all_oligos_gen_filledp == false) {
    Stage1_fill_all_oligos_gen(stage1,querylength,genestrand);
  }

  if (try_sensedir == SENSE_FORWARD) {
    nspliceends = spliceends_trim5_sense(&max_sense_prob,&(*partnerp),new,new_inward,qend,
					 inbound_genomicpos,outbound_genomicpos,middle_genomicpos,
					 mismatch_positions,total_nmismatches,univdiagonal,querylength,
					 query_compress,queryptr,chroffset,chrhigh,
					 knownsplicing,max_nconsecutive,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 stage1,streamspace_max_alloc,streamspace_alloc,
					 univdiagpool,medial_splicesite_prob,distal_splicesite_prob,
					 plusp,genestrand,localdb_nmismatches_allowed,
					 search_localdb_p,innerp);
    found_sensedir = SENSE_FORWARD;

  } else if (try_sensedir == SENSE_ANTI) {
    nspliceends = spliceends_trim5_anti(&max_antisense_prob,&(*partnerp),new,new_inward,qend,
					inbound_genomicpos,outbound_genomicpos,middle_genomicpos,
					mismatch_positions,total_nmismatches,univdiagonal,querylength,
					query_compress,queryptr,chroffset,chrhigh,
					knownsplicing,max_nconsecutive,
					novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					stage1,streamspace_max_alloc,streamspace_alloc,
					univdiagpool,medial_splicesite_prob,distal_splicesite_prob,
					plusp,genestrand,localdb_nmismatches_allowed,
					search_localdb_p,innerp);
    found_sensedir = SENSE_ANTI;

  } else {
    /* SENSE_NULL */
    fprintf(stderr,"try_sensedir is neither SENSE_FORWARD nor SENSE_ANTI\n");
    abort();
  }


  Spliceendsgen_return(spliceendsgen,&new_inward);
  if (nspliceends == 0) {
    debug1(printf("Got no spliceends\n"));
    /* Spliceends_free(&new,spliceendspool); */
    Spliceendsgen_return(spliceendsgen,&new);
    return (T) NULL;

  } else {
    new->nspliceends = nspliceends;
    new->sensedir = found_sensedir;
    if (found_sensedir == SENSE_FORWARD) {
      if (plusp) {
	new->splicetype = ACCEPTOR;
      } else {
	new->splicetype = ANTIDONOR;
      }

    } else {
      /* SENSE_ANTI */
      if (plusp) {
	new->splicetype = ANTIDONOR;
      } else {
	new->splicetype = ACCEPTOR;
      }
    }

#ifdef DEBUG1
    printf("trim_5 got %d spliceends\n",nspliceends);
    for (int i = 0; i < nspliceends; i++) {
      if (new->partners[i] == 0) {
	printf("%d (no partner)\n",new->splice_qpos[i]);
      } else {
	printf("%u %d..%d\n",new->distal_univdiags[i]->univdiagonal,new->distal_univdiags[i]->qstart,new->distal_univdiags[i]->qend);
      }
    }
#endif

    /* Sorted primarily by medial_qpos, and secondarily by descending distal positions, which Path_solve depends on */
    return new;
  }
}


/* Note: this procedure is called repeatedly with the same univdiagonal */
static int
extend_trim3 (Univdiag_T **univdiags, int start_qpos, Stage1_T stage1,
	      Univcoord_T univdiagonal, int querylength, bool plusp,
	      int slop, int insertion_slop, Univdiagpool_T univdiagpool) {

  int nmerged, nunivdiags, ndiagonals, npositions, i, j;
  List_T univdiag_list = NULL;
  Univcoord_T *diagonals;
  Univdiag_T univdiag;

#ifdef LARGE_GENOMES
  unsigned char *positions_high;
  UINT4 *positions;
#else
  Univcoord_T *positions;
#endif
  int querypos, qpos;

  if (plusp == true) {
    /* plus */
    debug2(printf("extend_trim3_plus, at univdiagonal %u\n",univdiagonal));

    /* Check from medial to distal.  Previously greedy, now exhaustive */
    for (qpos = start_qpos; qpos <= querylength - index1part; qpos++) {
      querypos = qpos;
      debug2(printf("Testing kmer at qpos %d => querypos %d\n",qpos,querypos));

      if (stage1->validp[querypos] == false) {
	/* Skip */
      } else {
	if (stage1->plus_retrievedp[querypos] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->plus_positions_high[querypos];
#endif
	  positions = stage1->plus_positions[querypos];
	  npositions = stage1->plus_npositions[querypos];
	} else {
	  assert(stage1->plus_positions[querypos] == NULL);
#ifdef LARGE_GENOMES
	  npositions = stage1->plus_npositions[querypos] =
	    Indexdb_largeptr(&stage1->plus_positions_high[querypos],&stage1->plus_positions[querypos],
			     indexdb,stage1->forward_oligos[querypos]);
	  positions_high = stage1->plus_positions_high[querypos];
#else
	  npositions = stage1->plus_npositions[querypos] =
	    Indexdb_ptr(&stage1->plus_positions[querypos],indexdb,
			stage1->forward_oligos[querypos]);
#endif
	  positions = stage1->plus_positions[querypos];
	  stage1->plus_retrievedp[querypos] = true;
	}
      
	if (npositions > 0) {
	  diagonals = (Univcoord_T *) MALLOC(npositions*sizeof(Univcoord_T));
	  ndiagonals = Intersect_higher(diagonals,
#ifdef LARGE_GENOMES
					positions_high,
#endif
					positions,npositions,
					/*diagterm1, plus*/querylength - querypos,
					/*set2*/&univdiagonal,/*length2*/1,
					slop,insertion_slop);
	  for (i = 0; i < ndiagonals; i++) {
	    univdiag_list = Univdiagpool_push(univdiag_list,univdiagpool,/*qstart*/qpos,/*qend*/qpos + index1part,
					      /*nmismatches*/0,diagonals[i]
					      univdiagpool_trace(__FILE__,__LINE__));
	  }
	  FREE(diagonals);
	}
      }
    }

  } else {
    /* minus */
    debug2(printf("extend_trim3_minus at univdiagonal %u\n",univdiagonal));

    /* Check from medial to distal */
    for (qpos = start_qpos; qpos <= querylength - index1part; qpos++) {
      querypos = (querylength - index1part) - qpos;
      debug2(printf("Testing kmer at qpos %d => querypos %d\n",qpos,querypos));

      if (stage1->validp[querypos] == false) {
	/* Skip */
      } else {
	if (stage1->minus_retrievedp[querypos] == true) {
#ifdef LARGE_GENOMES
	  positions_high = stage1->minus_positions_high[querypos];
#endif
	  positions = stage1->minus_positions[querypos];
	  npositions = stage1->minus_npositions[querypos];
	} else {
	  assert(stage1->minus_positions[querypos] == NULL);
#ifdef LARGE_GENOMES
	  npositions = stage1->minus_npositions[querypos] =
	    Indexdb_largeptr(&stage1->minus_positions_high[querypos],&stage1->minus_positions[querypos],
			     indexdb,stage1->revcomp_oligos[querypos]);
	  positions_high = stage1->minus_positions_high[querypos];
#else
	  npositions = stage1->minus_npositions[querypos] =
	    Indexdb_ptr(&stage1->minus_positions[querypos],indexdb,
			stage1->revcomp_oligos[querypos]);
#endif
	  positions = stage1->minus_positions[querypos];
	  stage1->minus_retrievedp[querypos] = true;
	}
      
	if (npositions > 0) {
	  diagonals = (Univcoord_T *) MALLOC(npositions*sizeof(Univcoord_T));
	  ndiagonals = Intersect_higher(diagonals,
#ifdef LARGE_GENOMES
					positions_high,
#endif
					positions,npositions,
					/*diagterm1, minus*/querypos + index1part,
					/*set2*/&univdiagonal,/*length2*/1,
					slop,insertion_slop);
	  for (i = 0; i < ndiagonals; i++) {
	    univdiag_list = Univdiagpool_push(univdiag_list,univdiagpool,/*qstart*/qpos,/*qend*/qpos + index1part,
					      /*nmismatches*/0,diagonals[i]
					      univdiagpool_trace(__FILE__,__LINE__));
	  }
	  FREE(diagonals);
	}
      }
    }
  }

  if (univdiag_list == 0) {
    *univdiags = (Univdiag_T *) NULL;
    return 0;

  } else {
    nunivdiags = List_length(univdiag_list);
    *univdiags = Univdiag_new_array(nunivdiags);
    Univdiag_transfer_list_to_array(*univdiags,univdiag_list);
    Univdiagpool_free_list(&univdiag_list,univdiagpool
			   univdiagpool_trace(__FILE__,__LINE__));

    if (nunivdiags == 1) {
      return 1;

    } else {
      /* Merge univdiags.  Output in ascending order */
      qsort((*univdiags)[0],nunivdiags,sizeof(struct Univdiag_T),Univdiag_struct_diagonal_fwd_cmp);
      
#ifdef DEBUG9
      printf("Before merging: %d univdiags\n",nunivdiags);
      for (i = 0; i < nunivdiags; i++) {
	univdiag = (*univdiags)[i];
	printf("%u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
      }
      printf("\n");
#endif
    
      nmerged = 0; i = 0;
      while (i < nunivdiags) {
	univdiag = (*univdiags)[i];

	j = i + 1;
	while (j < nunivdiags && (*univdiags)[j]->univdiagonal == univdiag->univdiagonal) {
	  univdiag->qend = (*univdiags)[j]->qend;
	  j++;
	}

	Univdiag_transfer((*univdiags)[nmerged++],univdiag);
	i = j;
      }

#ifdef DEBUG9
      printf("After merging: returning %d univdiags\n",nmerged);
      for (i = 0; i < nmerged; i++) {
	univdiag = (*univdiags)[i];
	printf("%u %d..%d\n",univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
      }
      printf("\n");
#endif

      return nmerged;
    }
  }

}


static Univdiag_T
novel_trim3_indel (Univcoord_T outbound_genomicpos, Univcoord_T univdiagonal, int querylength,
		   Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
		   Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		   Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		   bool plusp, int genestrand, int localdb_nmismatches_allowed) {

  Univdiag_T best_univdiag = NULL;
  Univcoord_T deletion_univdiagonal, insertion_univdiagonal;
  bool sortedp, trimmedp;

  int local_nmismatches, nmatches, matchlength;
  int indel_qpos;
  int ndiagonals, i;
  int best_adj, adj;
  Univcoord_T left = univdiagonal - querylength;

  /* Test for deletion */
  indel_qpos = outbound_genomicpos - left;

  debug9(printf("Start with univdiagonal %u\n",univdiagonal));
  deletion_univdiagonal = add_bounded(univdiagonal,max_deletionlen,chrhigh);
  debug9(printf("Add %d (bounded by %u) to yield (high) deletion_univdiagonal %u\n",
		max_deletionlen,chrhigh,deletion_univdiagonal));

  debug9(printf("Testing indel_qpos %d for deletion\n",indel_qpos));
  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				localdb,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,
				queryptr,/*pos5*/indel_qpos,/*pos3*/querylength,querylength,
				/*anchor_univdiagonal*/univdiagonal,
				/*low_univdiagonal*/univdiagonal,/*high_univdiagonal*/deletion_univdiagonal,
				query_compress,plusp,genestrand,genomebits,localdb_nmismatches_allowed,
				/*extend5p*/false,/*trim5p*/false,/*trim3p*/true)) == 0) {
    /* Skip */

  } else if ((nmatches = matchlength - local_nmismatches) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*local_nmismatches) {
    debug9(printf("(6) Got %d localdb diagonals with matchlength %d and local_nmismatches %d => Skipping\n",
		  ndiagonals,matchlength,local_nmismatches));
    /* Skip */
      
  } else {
    debug9(printf("(7) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		  ndiagonals,matchlength,local_nmismatches));

    best_adj = 0;
    for (i = 0; i < ndiagonals; i++) {
      debug9(printf("indel_univdiagonal %u, adj %d\n",
		    univdiags_alloc[i]->univdiagonal,univdiags_alloc[i]->univdiagonal - univdiagonal));
      /* Deletion */
      if (univdiags_alloc[i]->univdiagonal <= univdiagonal) {
	/* Skip */
      } else if ((adj = univdiags_alloc[i]->univdiagonal - univdiagonal) >= matchlength) {
	debug9(printf("adj >= matchlength\n"));
      } else if (best_adj == 0 || adj < best_adj) {
	best_univdiag = univdiags_alloc[i];
	best_adj = adj;
      }
    }
    if (best_univdiag != NULL) {
      debug9(printf("Returning deletion %u %d..%d\n",
		    best_univdiag->unidiagonal,best_univdiag->qstart,best_univdiag->qend));
      return best_univdiag;
    }
  }


  /* Test for insertion */
  insertion_univdiagonal = subtract_bounded(univdiagonal,max_insertionlen,chroffset);
  debug9(printf("Subtract %d (bounded by %u) to yield (low) insertion_univdiagonal %u\n",
		max_insertionlen,chroffset,insertion_univdiagonal));

  debug9(printf("Testing indel_qpos %d for insertion\n",indel_qpos));
  if ((ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
				localdb,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,queryptr,
				/*pos5*/indel_qpos,/*pos3*/querylength,querylength,
				/*anchor_univdiagonal*/univdiagonal,
				/*low_univdiagonal*/insertion_univdiagonal,/*high_univdiagonal*/univdiagonal,
				query_compress,plusp,genestrand,genomebits,localdb_nmismatches_allowed,
				/*extend5p*/false,/*trim5p*/true,/*trim3p*/true)) == 0) {
    /* Skip */

  } else if ((nmatches = matchlength - local_nmismatches) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*local_nmismatches) {
    debug9(printf("(8) Got %d localdb diagonals with matchlength %d and local_nmismatches %d => Skipping\n",
		  ndiagonals,matchlength,local_nmismatches));
    /* Skip */
      
  } else {
    debug9(printf("(9) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		  ndiagonals,matchlength,local_nmismatches));

    best_adj = 0;
    for (i = 0; i < ndiagonals; i++) {
      debug9(printf("indel_univdiagonal %u, adj %d\n",
		    univdiags_alloc[i]->univdiagonal,univdiags_alloc[i]->univdiagonal - univdiagonal));
      /* Insertion */
      if (univdiags_alloc[i]->univdiagonal >= univdiagonal) {
	/* Skip */
      } else if ((adj = univdiagonal - univdiags_alloc[i]->univdiagonal) >= matchlength) {
	debug9(printf("adj >= matchlength\n"));
      } else if (best_adj == 0 || adj < best_adj) {
	best_univdiag = univdiags_alloc[i];
	best_adj = adj;
      }
    }
    
    if (best_univdiag != NULL) {
      debug9(printf("Returning insertion %u %d..%d\n",
		    best_univdiag->unidiagonal,best_univdiag->qstart,best_univdiag->qend));
      return best_univdiag;
    }
  }
  
  return (Univdiag_T) NULL;
}


static int
solve_trim3 (int *min_nmismatches, double *max_prob, double *max_medial_prob, bool *max_medial_setp,
	     bool *partnerp, bool *boundedp, int *nsites, T this, int k, bool *sitep,
	     double (*Medial_prob_fcn)(Univcoord_T,Univcoord_T),
	     double (*Distal_prob_fcn)(Univcoord_T,Univcoord_T),

	     Univcoord_T genomicpos, int pos3, int splice_qpos, int qstart, int querylength, int mismatchi,
	     Univcoord_T univdiagonal, Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,

	     Univdiag_T *known_univdiags_alloc, int known_nunivdiags,
			
	     Compress_T query_compress, char *queryptr, Univcoord_T chroffset, Univcoord_T chrhigh,
	     Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
	     Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
	     Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
	     bool plusp, int genestrand, int localdb_nmismatches_allowed,
	     bool search_localdb_p, bool innerp, int max_npartners, bool inwardp,
	     bool medial_salvagep) {

  int kstart = k;
  Univdiag_T *univdiags, *novel_univdiags, indel_univdiag;
  int nunivdiags, novel_nunivdiags = 0;
  bool *knownp, sortedp, trimmedp;
  bool allocp = false;		/* For extend_trim3 */
  Chrpos_T positive_gap_distance;

  Univcoord_T distal_genomicpos;
  int local_nmismatches, matchlength, trimpos;
  double medial_prob, max_distal_prob, distal_prob;
  int nmismatches;
  int npartners, best_partneri, i;

  Univcoord_T left = univdiagonal - querylength;

#ifdef DEBUG1
  medial_prob = Medial_prob_fcn(genomicpos,chroffset);
  printf("3' %s: %u %u %d site:%d %f medial\n",
	 plusp ? "plus" : "minus",genomicpos,genomicpos-chroffset,genomicpos-left,
	 sitep[genomicpos-left],medial_prob);
#endif

  if (known_nunivdiags > 0) {
    medial_prob = 1.0;
    *max_medial_setp = true;	/* was not set before */
  } else if (sitep[genomicpos-left] == false) {
    medial_prob = 0.0;
    *max_medial_setp = false;
  } else if ((medial_prob = Medial_prob_fcn(genomicpos,chroffset)) > *max_medial_prob) {
    *max_medial_prob = medial_prob;
    *max_medial_setp = true;
  } else {
    *max_medial_setp = false;
  }

  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  if (medial_prob >= medial_splicesite_prob || medial_salvagep == true) {
    debug1(printf("Testing splice_qpos %d with mismatchi %d\n",splice_qpos,mismatchi));
    if (0 && innerp == false && querylength - splice_qpos < 4) {
      /* Forces program to find a bad splicing solution */
      /* Skip */

    } else if (stage1 != NULL &&
	       (novel_nunivdiags = extend_trim3(&novel_univdiags,splice_qpos,stage1,
						univdiagonal,querylength,plusp,
						/*slop*/positive_gap_distance,
						/*insertion_slop*/0,univdiagpool)) > 0) {
      allocp = true;
      debug1(printf("(2) Got %d indexdb univdiags\n",novel_nunivdiags));
      /* No need to sort.  Iterate from start */

    } else if (search_localdb_p == true &&
	       (novel_nunivdiags = Localdb_get(&sortedp,&trimmedp,&matchlength,&local_nmismatches,univdiags_alloc,
					       localdb,diagonals_alloc,localdb_alloc,
					       stage1,streamspace_max_alloc,streamspace_alloc,
					       queryptr,/*pos5*/splice_qpos,pos3,querylength,
					       /*anchor_univdiagonal*/univdiagonal,
					       low_univdiagonal,high_univdiagonal,
					       query_compress,plusp,genestrand,genomebits,
					       localdb_nmismatches_allowed,/*extend5p*/false,
					       /*trim5p*/false,/*trim3p*/true)) > 0) {
      debug1(printf("(10) Got %d localdb diagonals with matchlength %d and local_nmismatches %d\n",
		    novel_nunivdiags,matchlength,local_nmismatches));
      if (novel_nunivdiags == 1) {
	/* No need to sort */
      } else if (sortedp == true) {
	/* No need to sort */
      } else {
	qsort(univdiags_alloc[0],novel_nunivdiags,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
      }
      novel_univdiags = univdiags_alloc;
    }
	
    if ((nunivdiags = merge_known_novel(&univdiags,&knownp,
					known_univdiags_alloc,known_nunivdiags,
					novel_univdiags,novel_nunivdiags)) == 0) {
      /* Mark as a trimming position without a partner */
      debug1(printf("Got no novel diagonals\n"));
      this->splice_qpos[k] = splice_qpos;
      this->medial_nmismatches[k] = mismatchi;
      this->partners[k] = 0.0;

      this->medial_probs[k] = medial_prob;
      this->distal_probs[k] = 0.0;
      this->medial_support[k] = (splice_qpos - qstart); /*- 4*mismatchi for adj_support*/
      this->distal_support[k] = 0;
      k++;

      (*nsites)++;
    
    } else {
      /* Process starting from the closest diagonal */
      best_partneri = -1; max_distal_prob = 0.0;
      for (i = 0; i < nunivdiags; i++) {
	distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos;
	if (knownp[i] == true) {
	  distal_prob = 1.0;
	} else {
	  distal_prob = Distal_prob_fcn(distal_genomicpos,chroffset);
	}
	/* debug1(printf("Evaluating diagonal %u, distal_genomicpos %u, prob %f\n",
	   univdiags[i]->univdiagonal,distal_genomicpos,distal_prob)); */

	if (distal_prob > max_distal_prob) {
	  max_distal_prob = distal_prob;
	  best_partneri = i;
	}
      }

      npartners = 0;
      if (max_distal_prob >= distal_splicesite_prob || medial_salvagep == true) {
	i = 0;
	while (npartners < max_npartners && i < nunivdiags) {
	  distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos;
	  if (knownp[i] == true) {
	    distal_prob = 1.0;
	  } else {
	    distal_prob = Distal_prob_fcn(distal_genomicpos,chroffset);
	  }

	  if (distal_prob >= max_distal_prob - PROB_SLOP) {

	    debug5(printf("Univdiag %u, %d..%d\n",univdiags[i]->univdiagonal,univdiags[i]->qstart,univdiags[i]->qend));
	    if ((trimpos = Genomebits_trim_qend(&(univdiags[i]->nmismatches),query_compress,
						genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
						/*pos5*/splice_qpos,pos3,plusp,genestrand)) >= univdiags[i]->qend) {
	      /* Extends to splice_qpos */
	      debug5(printf("(3) At univdiagonal %u, trimpos %d >= qend %d (preserved), so extends to splice_qpos %d\n",
			    univdiags[i]->univdiagonal,trimpos,univdiags[i]->qend,splice_qpos));
	      univdiags[i]->qstart = splice_qpos;
	      univdiags[i]->qend = trimpos;
	    } else {
	      /* Cassette mismatch, perhaps due to a spliceindel */
	      debug5(printf("(3) At univdiagonal %u, trimpos %d < qend %d, so cassette mismatch.  ",
			    univdiags[i]->univdiagonal,trimpos,univdiags[i]->qend));
	      trimpos = univdiags[i]->qend =
		Genomebits_trim_qend(&(univdiags[i]->nmismatches),query_compress,
				     genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
				     /*pos5*/univdiags[i]->qstart,pos3,plusp,genestrand);	      
	      debug5(printf("New bounds are %d..%d\n",univdiags[i]->qstart,univdiags[i]->qend));
	    }

	    debug1(printf("For univdiagonal %u from %d to %d, distal trimpos is %d\n",
			  univdiags[i]->univdiagonal,univdiags[i]->qstart,pos3,trimpos));
	  
	    if (/*univdiags[i]->qend*/trimpos == univdiags[i]->qstart) {
	      /* Skip: Not a good partner */
	      debug1(printf("trimpos %d == qstart %d (length zero), so not a good partner\n",
			    trimpos,univdiags[i]->qstart));
	      
	    } else if (pos3 - splice_qpos < 8 || univdiags[i]->qend - univdiags[i]->qstart >= 8) {
	      /* Adequate extension */
	      this->splice_qpos[k] = splice_qpos;
	      this->medial_nmismatches[k] = mismatchi;
	      this->partners[k] = distal_genomicpos;

	      Univdiag_transfer(this->distal_univdiags[k],univdiags[i]);

	      this->medial_probs[k] = medial_prob;
	      this->distal_probs[k] = distal_prob;
	      this->medial_support[k] = (splice_qpos - qstart); /*- 4*mismatchi for adj_support*/
	      this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart);
	      /*- 4*univdiags[i]->nmismatches for adj_support*/
	      
	      if ((nmismatches = mismatchi + univdiags[i]->nmismatches + (querylength - trimpos)) < *min_nmismatches) {
		*min_nmismatches = nmismatches;
		*max_prob = medial_prob + distal_prob;
	      } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
		*max_prob = medial_prob + distal_prob;
	      }
	      debug1(printf("(3) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			    univdiags[i]->univdiagonal,this->partners[k],univdiags[i]->univdiagonal - univdiagonal,
			    this->medial_probs[k],this->distal_probs[k],
			    this->splice_qpos[k],this->medial_nmismatches[k],
			    univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	      k++;
	  
	      npartners++;
	      *partnerp = true;

	    } else if (localdb == NULL) {
	      /* Skip */

	    } else if ((indel_univdiag =
			novel_trim3_indel(/*genomicpos*/univdiags[i]->univdiagonal - querylength + univdiags[i]->qend,
					  univdiags[i]->univdiagonal,querylength,
					  query_compress,queryptr,chroffset,chrhigh,
					  univdiags_alloc,diagonals_alloc,localdb_alloc,
					  stage1,streamspace_max_alloc,streamspace_alloc,
					  plusp,genestrand,localdb_nmismatches_allowed)) == (Univdiag_T) NULL) {
	      /* Skip: Non-extension without distal indel */
	      debug1(printf("non-extension without distal indel\n"));

	    } else {
	      /* Provide indel univdiag instead so we find a spliceindel */
	      debug1(printf("found indel %u %d..%d\n",
			    indel_univdiag->univdiagonal,indel_univdiag->qstart,indel_univdiag->qend));

	      this->splice_qpos[k] = splice_qpos;
	      this->medial_nmismatches[k] = mismatchi;
	      this->partners[k] = distal_genomicpos;

	      Univdiag_transfer(this->distal_univdiags[k],indel_univdiag);

	      this->medial_probs[k] = medial_prob;
	      this->distal_probs[k] = distal_prob;
	      this->medial_support[k] = (splice_qpos - qstart); /*- 4*mismatchi for adj_support*/
	      this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart) +
		(indel_univdiag->qend - indel_univdiag->qstart);
	      /*- 4*univdiags[i]->nmismatches for adj_support*/
	      if (this->distal_support[k] > querylength - splice_qpos) {
		this->distal_support[k] = querylength - splice_qpos;
	      }
	      
	      if ((nmismatches = mismatchi + univdiags[i]->nmismatches + (querylength - trimpos)) < *min_nmismatches) {
		*min_nmismatches = nmismatches;
		*max_prob = medial_prob + distal_prob;
	      } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
		*max_prob = medial_prob + distal_prob;
	      }
	      debug1(printf("(3) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			    univdiags[i]->univdiagonal,this->partners[k],univdiags[i]->univdiagonal - univdiagonal,
			    this->medial_probs[k],this->distal_probs[k],
			    this->splice_qpos[k],this->medial_nmismatches[k],
			    univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	      k++;
	  
	      npartners++;
	      *partnerp = true;
	    }
	  }
	  
	  i++;
	}
      }
      
      if (i < nunivdiags) {
	/* Limited by max_npartners */
	debug1(printf("BOUNDED\n"));
	*boundedp = true;

      } else if (medial_salvagep == true) {
	/* Do not salvage both ends */

      } else if (npartners == 0 && max_distal_prob > 0.5) {
	/* Salvage: use best partneri */
	i = best_partneri;
	distal_genomicpos = univdiags[i]->univdiagonal - querylength + splice_qpos; /* retrieve */

	debug5(printf("Univdiag %u, %d..%d\n",univdiags[i]->univdiagonal,univdiags[i]->qstart,univdiags[i]->qend));
	if ((trimpos = Genomebits_trim_qend(&(univdiags[i]->nmismatches),query_compress,
					    genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
					    /*pos5*/splice_qpos,pos3,plusp,genestrand)) >= univdiags[i]->qend) {
	  /* Extends to splice_qpos */
	  debug5(printf("(4) At univdiagonal %u, trimpos %d >= qend %d (preserved), so extends to splice_qpos %d\n",
			univdiags[i]->univdiagonal,trimpos,univdiags[i]->qend,splice_qpos));
	  univdiags[i]->qstart = splice_qpos;
	  univdiags[i]->qend = trimpos;
	} else {
	  /* Cassette mismatch, perhaps due to a spliceindel */
	  debug5(printf("(4) At univdiagonal %u, trimpos %d < qend %d, so cassette mismatch.  ",
			univdiags[i]->univdiagonal,trimpos,univdiags[i]->qend));
	  trimpos = univdiags[i]->qend =
	    Genomebits_trim_qend(&(univdiags[i]->nmismatches),query_compress,
				 genomebits,/*univdiagonal*/univdiags[i]->univdiagonal,querylength,
				 /*pos5*/univdiags[i]->qstart,pos3,plusp,genestrand);
	  debug5(printf("New bounds are %d..%d\n",univdiags[i]->qstart,univdiags[i]->qend));
	}
	  
	debug1(printf("For univdiagonal %u from %d to %d, distal trimpos is %d\n",
		      univdiags[i]->univdiagonal,univdiags[i]->qstart,pos3,trimpos));
	  
	if (/*univdiags[i]->qend*/trimpos == univdiags[i]->qstart) {
	  /* Skip: Still not a good partner */

	} else if ((nmismatches = mismatchi + univdiags[i]->nmismatches + (querylength - trimpos)) > 1) {
	  /* Have a stricter standard for salvage cases */
  
	} else if (pos3 - splice_qpos < 8 || univdiags[i]->qend - univdiags[i]->qstart >= 8) {
	  /* Adequate extension */
	  this->splice_qpos[k] = splice_qpos;
	  this->medial_nmismatches[k] = mismatchi;
	  this->partners[k] = distal_genomicpos;

	  Univdiag_transfer(this->distal_univdiags[k],univdiags[i]);

	  this->medial_probs[k] = medial_prob;
	  this->distal_probs[k] = max_distal_prob; /* retrieve */
	  this->medial_support[k] = (splice_qpos - qstart); /*- 4*mismatchi for adj_support*/
	  this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart);
	  /*- 4*univdiags[i]->nmismatches for adj_support*/

	  if (nmismatches < *min_nmismatches) {
	    *min_nmismatches = nmismatches;
	    *max_prob = medial_prob + distal_prob;
	  } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
	    *max_prob = medial_prob + distal_prob;
	  }
	  debug1(printf("(4) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			univdiags[i]->univdiagonal,this->partners[k],univdiags[i]->univdiagonal - univdiagonal,
			this->medial_probs[k],this->distal_probs[k],
			this->splice_qpos[k],this->medial_nmismatches[k],
			univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	  k++;
	  
	  /* npartners++; */
	  *partnerp = true;

	} else if (localdb == NULL) {
	  /* Skip */

	} else if ((indel_univdiag =
		    novel_trim3_indel(/*genomicpos*/univdiags[i]->univdiagonal - querylength + univdiags[i]->qend,
				      univdiags[i]->univdiagonal,querylength,
				      query_compress,queryptr,chroffset,chrhigh,
				      univdiags_alloc,diagonals_alloc,localdb_alloc,
				      stage1,streamspace_max_alloc,streamspace_alloc,
				      plusp,genestrand,localdb_nmismatches_allowed)) == (Univdiag_T) NULL) {
	  /* Skip: Non-extension without distal indel */
	  debug1(printf("non-extension without distal indel\n"));

	} else {
	  /* Provide indel univdiag instead so we find a spliceindel */
	  debug1(printf("found indel %u %d..%d\n",
			indel_univdiag->univdiagonal,indel_univdiag->qstart,indel_univdiag->qend));

	  this->splice_qpos[k] = splice_qpos;
	  this->medial_nmismatches[k] = mismatchi;
	  this->partners[k] = distal_genomicpos;

	  Univdiag_transfer(this->distal_univdiags[k],indel_univdiag);

	  this->medial_probs[k] = medial_prob;
	  this->distal_probs[k] = max_distal_prob; /* retrieve */
	  this->medial_support[k] = (splice_qpos - qstart); /*- 4*mismatchi for adj_support*/
	  this->distal_support[k] = (univdiags[i]->qend - univdiags[i]->qstart) +
	    (indel_univdiag->qend - indel_univdiag->qstart);
	  /*- 4*univdiags[i]->nmismatches for adj_support*/
	  if (this->distal_support[k] > querylength - splice_qpos) {
	    this->distal_support[k] = querylength - splice_qpos;
	  }

	  if (nmismatches < *min_nmismatches) {
	    *min_nmismatches = nmismatches;
	    *max_prob = medial_prob + distal_prob;
	  } else if (nmismatches == *min_nmismatches && medial_prob + distal_prob > *max_prob) {
	    *max_prob = medial_prob + distal_prob;
	  }
	  debug1(printf("(4) %u %u dist:%u %f,%f qpos:%d mismatches:%d+%d distal:%d..%d\n",
			univdiags[i]->univdiagonal,this->partners[k],univdiags[i]->univdiagonal - univdiagonal,
			this->medial_probs[k],this->distal_probs[k],
			this->splice_qpos[k],this->medial_nmismatches[k],
			univdiags[i]->nmismatches,univdiags[i]->qstart,univdiags[i]->qend));
	  k++;
	  
	  /* npartners++; */
	  *partnerp = true;
	}
      }
      
      if (inwardp == false) {
	Spliceends_reverse(this,kstart,k);
      }

      Univdiag_free_array(&univdiags);
      FREE(knownp);
      (*nsites)++;
    }

    if (allocp == true) {
      Univdiag_free_array(&novel_univdiags);
    }
  }

  return k;
}



static int
spliceends_trim3_sense (double *max_prob, bool *partnerp, T this_outward, T this_inward, int qstart,
			Univcoord_T inbound_genomicpos, Univcoord_T outbound_genomicpos, Univcoord_T middle_genomicpos,
			int *mismatch_positions, int total_nmismatches,
			Univcoord_T univdiagonal, int querylength,
			Compress_T query_compress, char *queryptr,
			Univcoord_T chroffset, Univcoord_T chrhigh,
			Knownsplicing_T knownsplicing, int max_nconsecutive,
			Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
			Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
			bool plusp, int genestrand, int localdb_nmismatches_allowed,
			bool search_localdb_p, bool innerp) {
  int nspliceends = 0, k, k_inward = 0, k_outward = 0;
  Chrpos_T splice_chrbound, positive_gap_distance;

  Univcoord_T *endpoints_inward, *endpoints_outward;
  uint64_t low_rank_inward = 0, high_rank_inward = 0, rank_inward = 0,
    low_rank_outward = 0, high_rank_outward = 0, rank_outward = 0;
  int n_inward, n_outward;

  Univdiag_T *known_univdiags_alloc = NULL, univdiag;
  Univcoord_T best_genomicpos, genomicpos_inward, genomicpos_outward, low_univdiagonal, high_univdiagonal;
  int known_nunivdiags;

  int pos3, best_splice_qpos, splice_qpos_inward, splice_qpos_outward, nconsecutive;
  int best_mismatchi, mismatchi_inward, mismatchi_outward;
  int nsites = 0, i;
  Univcoord_T left = univdiagonal - querylength;
  double max_medial_prob = 0.0;
  int min_nmismatches = querylength;
  bool boundedp = false, max_medial_setp, best_inwardp;

  int max_npartners = (innerp == true) ? MAX_INNER_PARTNERS_PER_SITE : MAX_OUTER_PARTNERS_PER_SITE;
  bool *sitep;

  *max_prob = 0.0;
  *partnerp = false;

  assert(outbound_genomicpos > inbound_genomicpos);

  /* Prevent splicing past chromosome bounds */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  if (univdiagonal + positive_gap_distance >= chrhigh) {
    splice_chrbound = chrhigh - univdiagonal;
  } else {
    splice_chrbound = positive_gap_distance;
  }

  debug1(printf("Start with univdiagonal %u\n",univdiagonal));
  low_univdiagonal = univdiagonal + MIN_INTRON_LENGTH; /* Add so we don't return a continuation */
  debug1(printf("Add %d to yield low_univdiagonal %u\n",MIN_INTRON_LENGTH,low_univdiagonal));
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos3 = (low_univdiagonal <= chrhigh) ? querylength : (int) (chrhigh - left);
#else
  pos3 = querylength;
#endif
  high_univdiagonal = add_bounded(low_univdiagonal,splice_chrbound,chrhigh);
  debug1(printf("Add %d (bounded by %u) to yield high_univdiagonal %u\n",
		splice_chrbound,chrhigh,high_univdiagonal));
  debug1(printf("  (3) splice_chrbound %u, low %u, high %u\n",
		splice_chrbound,low_univdiagonal,high_univdiagonal));

#ifdef DEBUG1
  printf("(3) mismatch positions:");
  for (i = 0; i <= total_nmismatches; i++) {
    printf(" %d",mismatch_positions[i]);
  }
  printf("\n");
#endif


  if (low_univdiagonal > high_univdiagonal) {
    /* Skip */

  } else if (plusp) {
    debug1(printf("spliceends_trim3_sense, plus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_donors(&low_rank_inward,&high_rank_inward,
					      knownsplicing,univdiagonal,querylength,
					      (inbound_genomicpos + 1) - left,(middle_genomicpos + 1) - left);
      rank_inward = high_rank_inward;

#ifdef DEBUG1
      printf("(5) Knownsplicing donors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     inbound_genomicpos + 1,middle_genomicpos + 1,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank+1] - endpoints_inward[2*rank]));
      }
#endif

      endpoints_outward = Knownsplicing_donors(&low_rank_outward,&high_rank_outward,
					       knownsplicing,univdiagonal,querylength,
					       (middle_genomicpos + 1) - left,(outbound_genomicpos + 1) - left);
      rank_outward = low_rank_outward;

#ifdef DEBUG1
      printf("(5) Knownsplicing donors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos + 1,outbound_genomicpos + 1,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank+1] - endpoints_outward[2*rank]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_donor_gtgc_sites(sitep,left,/*pos5*/inbound_genomicpos - left,/*pos3*/outbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] >= splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos + 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] >= splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward > inbound_genomicpos || genomicpos_outward <= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward > inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward > low_rank_inward && endpoints_inward[2*(rank_inward - 1)] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*(rank_inward - 1) + 1] + querylength - splice_qpos_inward;
	  univdiag->qstart = splice_qpos_inward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward--;
	}

	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_donor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			       
			       genomicpos_inward,pos3,splice_qpos_inward,qstart,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	splice_qpos_inward--;
	genomicpos_inward--;
	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
      }

      if (genomicpos_outward <= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward < high_rank_outward && endpoints_outward[2*rank_outward] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*rank_outward + 1] + querylength - splice_qpos_outward;
	  univdiag->qstart = splice_qpos_outward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward++;
	}
	
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_donor_prob,
				/*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			
				genomicpos_outward,pos3,splice_qpos_outward,qstart,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);
	
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
	splice_qpos_outward++;
	genomicpos_outward++;
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_donor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			       
			       best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_donor_prob,
				/*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			
				best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);

  } else {
    debug1(printf("spliceends_trim3_sense, minus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_antiacceptors(&low_rank_inward,&high_rank_inward,
						     knownsplicing,univdiagonal,querylength,
						     (inbound_genomicpos + 1) - left,(middle_genomicpos + 1) - left);
      rank_inward = high_rank_inward;

#ifdef DEBUG1
      printf("(6) Knownsplicing antiacceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     inbound_genomicpos + 1,middle_genomicpos + 1,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank+1] - endpoints_inward[2*rank]));
      }
#endif

      endpoints_outward = Knownsplicing_antiacceptors(&low_rank_outward,&high_rank_outward,
						      knownsplicing,univdiagonal,querylength,
						      (middle_genomicpos + 1) - left,(outbound_genomicpos + 1) - left);
      rank_outward = low_rank_outward;

#ifdef DEBUG1
      printf("(5) Knownsplicing antiacceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos + 1,outbound_genomicpos + 1,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank+1] - endpoints_outward[2*rank]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_antiacceptor_ag_sites(sitep,left,/*pos5*/inbound_genomicpos - left,/*pos3*/outbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] >= splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos + 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] >= splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward > inbound_genomicpos || genomicpos_outward <= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward > inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward > low_rank_inward && endpoints_inward[2*(rank_inward - 1)] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*(rank_inward - 1) + 1] + querylength - splice_qpos_inward;
	  univdiag->qstart = splice_qpos_inward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward--;
	}

	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antidonor_prob,

			       genomicpos_inward,pos3,splice_qpos_inward,qstart,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	splice_qpos_inward--;
	genomicpos_inward--;
	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
      }

      if (genomicpos_outward <= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward < high_rank_outward && endpoints_outward[2*rank_outward] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*rank_outward + 1] + querylength - splice_qpos_outward;
	  univdiag->qstart = splice_qpos_outward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward++;
	}

	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antidonor_prob,
			
				genomicpos_outward,pos3,splice_qpos_outward,qstart,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);
	
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
	splice_qpos_outward++;
	genomicpos_outward++;
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antidonor_prob,

			       best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antidonor_prob,
			
				best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);
  }

  debug1(printf("\n"));

  /* Combine inward and outward */
  debug1(printf("k_inward = %d, k_outward = %d, partnerp %d, max_prob %f, max_medial_prob %f\n",
		k_inward,k_outward,*partnerp,*max_prob,max_medial_prob));
  Spliceends_combine(this_outward,k_outward,this_inward,k_inward);
  this_outward->boundedp = boundedp;

  if ((k = k_inward + k_outward) == 0) {
    /* Skip */
    debug1(printf("k is zero\n"));
    nspliceends = 0;

  } else if (*partnerp == false) {
    /* min_nmismatches never set */
    debug1(printf("No partner, so use prob only\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_probs[i] > max_medial_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	/* No partner */
	/* Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]); */

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else if (innerp == false) {
    /* Keep only the best spliceends */
    debug1(printf("outer, so keeping only the best\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_nmismatches[i] + this_outward->distal_univdiags[i]->nmismatches + (querylength - this_outward->distal_univdiags[i]->qend) <= min_nmismatches &&
	  this_outward->medial_probs[i] + this_outward->distal_probs[i] > *max_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else {
    /* Keep all partner ends for later resolving of ambiguity */
    debug1(printf("inner with partners, so keeping all\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->partners[i] != 0) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }
  }

#ifdef DEBUG1
  printf("spliceends_trim3_sense yielded %d spliceends\n",nspliceends);
  for (k = 0; k < nspliceends; k++) {
    printf("%u %u %d\n",this_outward->partners[k],this_outward->partners[k] - chroffset,this_outward->splice_qpos[k]);
  }
#endif

  return nspliceends;
}


static int
spliceends_trim3_anti (double *max_prob, bool *partnerp, T this_outward, T this_inward, int qstart,
		       Univcoord_T inbound_genomicpos, Univcoord_T outbound_genomicpos, Univcoord_T middle_genomicpos,
		       int *mismatch_positions, int total_nmismatches,
		       Univcoord_T univdiagonal, int querylength,
		       Compress_T query_compress, char *queryptr,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       Knownsplicing_T knownsplicing, int max_nconsecutive,
		       Univdiag_T *univdiags_alloc, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		       Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       Univdiagpool_T univdiagpool, double medial_splicesite_prob, double distal_splicesite_prob,
		       bool plusp, int genestrand, int localdb_nmismatches_allowed,
		       bool search_localdb_p, bool innerp) {
  int nspliceends = 0, k, k_inward = 0, k_outward = 0;
  Chrpos_T splice_chrbound, positive_gap_distance;

  Univcoord_T *endpoints_inward, *endpoints_outward;
  uint64_t low_rank_inward = 0, high_rank_inward = 0, rank_inward = 0,
    low_rank_outward = 0, high_rank_outward = 0, rank_outward = 0;
  int n_inward, n_outward;

  Univdiag_T *known_univdiags_alloc = NULL, univdiag;
  Univcoord_T best_genomicpos, genomicpos_inward, genomicpos_outward, low_univdiagonal, high_univdiagonal;
  int known_nunivdiags;

  int pos3, best_splice_qpos, splice_qpos_inward, splice_qpos_outward, nconsecutive;
  int best_mismatchi, mismatchi_inward, mismatchi_outward;
  int nsites = 0, i;
  Univcoord_T left = univdiagonal - querylength;
  double max_medial_prob = 0.0;
  int min_nmismatches = querylength;
  bool boundedp = false, max_medial_setp, best_inwardp;

  int max_npartners = (innerp == true) ? MAX_INNER_PARTNERS_PER_SITE : MAX_OUTER_PARTNERS_PER_SITE;
  bool *sitep;

  *max_prob = 0.0;
  *partnerp = false;

  assert(outbound_genomicpos > inbound_genomicpos);

  /* Prevent splicing past chromosome bounds */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  if (univdiagonal + positive_gap_distance >= chrhigh) {
    splice_chrbound = chrhigh - univdiagonal;
  } else {
    splice_chrbound = positive_gap_distance;
  }

  debug1(printf("Start with univdiagonal %u\n",univdiagonal));
  low_univdiagonal = univdiagonal + MIN_INTRON_LENGTH; /* Add so we don't return a continuation */
  debug1(printf("Add %d to yield low_univdiagonal %u\n",MIN_INTRON_LENGTH,low_univdiagonal));
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos3 = (low_univdiagonal <= chrhigh) ? querylength : (int) (chrhigh - left);
#else
  pos3 = querylength;
#endif
  high_univdiagonal = add_bounded(low_univdiagonal,splice_chrbound,chrhigh);
  debug1(printf("Add %d (bounded by %u) to yield high_univdiagonal %u\n",
		splice_chrbound,chrhigh,high_univdiagonal));
  debug1(printf("  (4) splice_chrbound %u, low %u, high %u\n",
		splice_chrbound,low_univdiagonal,high_univdiagonal));

#ifdef DEBUG1
  printf("(4) mismatch positions:");
  for (i = 0; i <= total_nmismatches; i++) {
    printf(" %d",mismatch_positions[i]);
  }
  printf("\n");
#endif

  if (low_univdiagonal > high_univdiagonal) {
    /* Skip */

  } else if (plusp) {
    debug1(printf("spliceends_trim3_anti, plus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_antiacceptors(&low_rank_inward,&high_rank_inward,
						     knownsplicing,univdiagonal,querylength,
						     (inbound_genomicpos + 1) - left,(middle_genomicpos + 1) - left);
      rank_inward = high_rank_inward;

#ifdef DEBUG1
      printf("(7) Knownsplicing antiacceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     inbound_genomicpos + 1,middle_genomicpos + 1,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank+1] - endpoints_inward[2*rank]));
      }
#endif

      endpoints_outward = Knownsplicing_antiacceptors(&low_rank_outward,&high_rank_outward,
						      knownsplicing,univdiagonal,querylength,
						      (middle_genomicpos + 1) - left,(outbound_genomicpos + 1) - left);
      rank_outward = low_rank_outward;

#ifdef DEBUG1
	printf("(7) Knownsplicing antiacceptors at %u..%u yields low_rank %lu to high_rank %lu\n",
	       middle_genomicpos + 1,outbound_genomicpos + 1,low_rank_outward,high_rank_outward);
	for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	  printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
		 rank,/*splice_pos*/endpoints_outward[2*rank] - left,
		 endpoints_outward[2*rank],endpoints_outward[2*rank+1],
		 (Chrpos_T) (endpoints_outward[2*rank+1] - endpoints_outward[2*rank]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_antiacceptor_ag_sites(sitep,left,/*pos5*/inbound_genomicpos - left,/*pos3*/outbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] >= splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos + 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] >= splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward > inbound_genomicpos || genomicpos_outward <= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward > inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward > low_rank_inward && endpoints_inward[2*(rank_inward - 1)] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*(rank_inward - 1) + 1] + querylength - splice_qpos_inward;
	  univdiag->qstart = splice_qpos_inward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward--;
	}

	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antidonor_prob,

			       genomicpos_inward,pos3,splice_qpos_inward,qstart,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	splice_qpos_inward--;
	genomicpos_inward--;
	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
      }

      if (genomicpos_outward <= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward < high_rank_outward && endpoints_outward[2*rank_outward] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*rank_outward + 1] + querylength - splice_qpos_outward;
	  univdiag->qstart = splice_qpos_outward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward++;
	}

	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antidonor_prob,
			
				genomicpos_outward,pos3,splice_qpos_outward,qstart,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);
	
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
	splice_qpos_outward++;
	genomicpos_outward++;
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_antidonor_prob,

			       best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
		      
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_antiacceptor_prob,
				/*Distal_prob_fcn*/Maxent_hr_antidonor_prob,
			
				best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }
	
    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);

  } else {
    debug1(printf("spliceends_trim3_anti, minus\n"));
    sitep = (bool *) CALLOC(querylength,sizeof(bool));

    if (knownsplicing != NULL) {
      endpoints_inward = Knownsplicing_donors(&low_rank_inward,&high_rank_inward,
					      knownsplicing,univdiagonal,querylength,
					      (inbound_genomicpos + 1) - left,(middle_genomicpos + 1) - left);
      rank_inward = high_rank_inward;

#ifdef DEBUG1
      printf("(8) Knownsplicing donors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     inbound_genomicpos + 1,middle_genomicpos + 1,low_rank_inward,high_rank_inward);
      for (uint64_t rank = low_rank_inward; rank < high_rank_inward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_inward[2*rank] - left,
	       endpoints_inward[2*rank],endpoints_inward[2*rank+1],
	       (Chrpos_T) (endpoints_inward[2*rank+1] - endpoints_inward[2*rank]));
      }
#endif
      
      endpoints_outward = Knownsplicing_donors(&low_rank_outward,&high_rank_outward,
					       knownsplicing,univdiagonal,querylength,
					       (middle_genomicpos + 1) - left,(outbound_genomicpos + 1) - left);
      rank_outward = low_rank_outward;

#ifdef DEBUG1
      printf("(8) Knownsplicing donors at %u..%u yields low_rank %lu to high_rank %lu\n",
	     middle_genomicpos + 1,outbound_genomicpos + 1,low_rank_outward,high_rank_outward);
      for (uint64_t rank = low_rank_outward; rank < high_rank_outward; rank++) { /* For qstart, want lower qpos first */
	printf("Rank #%lu at qpos %d: %u..%u (splice distance %u)\n",
	       rank,/*splice_pos*/endpoints_outward[2*rank] - left,
	       endpoints_outward[2*rank],endpoints_outward[2*rank+1],
	       (Chrpos_T) (endpoints_outward[2*rank+1] - endpoints_outward[2*rank]));
      }
#endif

      n_inward = (int) (high_rank_inward - low_rank_inward);
      n_outward = (int) (high_rank_outward - low_rank_outward);
      if (n_inward == 0 && n_outward == 0) {
	/* Skip */
      } else if (n_inward >= n_outward) {
	known_univdiags_alloc = Univdiag_new_array(n_inward);
      } else {
	known_univdiags_alloc = Univdiag_new_array(n_outward);
      }
    }

    Genome_fill_donor_gtgc_sites(sitep,left,/*pos5*/inbound_genomicpos - left,/*pos3*/outbound_genomicpos - left + 1);

    genomicpos_inward = middle_genomicpos;
    splice_qpos_inward = genomicpos_inward - left;
    mismatchi_inward = total_nmismatches;
    while (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] >= splice_qpos_inward) {
      mismatchi_inward--;
    }

    genomicpos_outward = middle_genomicpos + 1;
    splice_qpos_outward = genomicpos_outward - left;
    mismatchi_outward = total_nmismatches;
    while (mismatchi_outward - 1 >= 0 && mismatch_positions[mismatchi_outward - 1] >= splice_qpos_outward) {
      mismatchi_outward--;
    }


    /* Iterate */
    nconsecutive = 0;
    while ((genomicpos_inward > inbound_genomicpos || genomicpos_outward <= outbound_genomicpos) &&
	   nsites < MAX_SITES /*&& *partnerp == false -- too greedy*/) {
      if (genomicpos_inward > inbound_genomicpos) {
	/* Inward */
	known_nunivdiags = 0;
	while (rank_inward > low_rank_inward && endpoints_inward[2*(rank_inward - 1)] == genomicpos_inward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_inward[2*(rank_inward - 1) + 1] + querylength - splice_qpos_inward;
	  univdiag->qstart = splice_qpos_inward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_inward,genomicpos_inward-chroffset,genomicpos_inward-left));
	  rank_inward--;
	}

	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_donor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_acceptor_prob,

			       genomicpos_inward,pos3,splice_qpos_inward,qstart,querylength,mismatchi_inward,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
			       
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/false);

	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_inward;
	  best_splice_qpos = splice_qpos_inward;
	  best_mismatchi = mismatchi_inward;
	  best_inwardp = true;
	}

	splice_qpos_inward--;
	genomicpos_inward--;
	if (mismatchi_inward - 1 >= 0 && mismatch_positions[mismatchi_inward - 1] == splice_qpos_inward) {
	  nconsecutive = 0;
	  mismatchi_inward--;
	} else {
	  nconsecutive++;
	}
      }

      if (genomicpos_outward <= outbound_genomicpos) {
	/* Outward */
	known_nunivdiags = 0;
	while (rank_outward < high_rank_outward && endpoints_outward[2*rank_outward] == genomicpos_outward) {
	  univdiag = known_univdiags_alloc[known_nunivdiags];
	  univdiag->univdiagonal = endpoints_outward[2*rank_outward + 1] + querylength - splice_qpos_outward;
	  univdiag->qstart = splice_qpos_outward;
	  univdiag->qend = Genomebits_trim_qend(&univdiag->nmismatches,query_compress,genomebits,
						univdiag->univdiagonal,querylength,
						/*pos5*/univdiag->qstart,/*pos3*/querylength,plusp,genestrand);
	  if (univdiag->qend > univdiag->qstart) {
	    known_nunivdiags++;
	  }

	  debug1(printf("3' %s: %u %u %d known medial\n",
			plusp ? "plus" : "minus",genomicpos_outward,genomicpos_outward-chroffset,genomicpos_outward-left));
	  rank_outward++;
	}
	
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_donor_prob,
				/*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			
				genomicpos_outward,pos3,splice_qpos_outward,qstart,querylength,mismatchi_outward,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/false);
	
	if (max_medial_setp == true) {
	  best_genomicpos = genomicpos_outward;
	  best_splice_qpos = splice_qpos_outward;
	  best_mismatchi = mismatchi_outward;
	  best_inwardp = false;
	}

	if (mismatchi_outward < total_nmismatches && mismatch_positions[mismatchi_outward] == splice_qpos_outward) {
	  mismatchi_outward++;
	}
	splice_qpos_outward++;
	genomicpos_outward++;
      }
    }

    debug1(printf("k_inward and k_outward are 0.  max_medial_prob is %f\n",max_medial_prob));
    if (k_inward == 0 && k_outward == 0 && max_medial_prob > 0.0) {
      /* Salvage medial */
      debug1(printf("Salvage medial\n"));
      if (best_inwardp == true) {
	k_inward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
			       &(*partnerp),&boundedp,&nsites,this_inward,k_inward,sitep,
			       /*Medial_prob_fcn*/Maxent_hr_donor_prob,
			       /*Distal_prob_fcn*/Maxent_hr_acceptor_prob,

			       best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
			       univdiagonal,low_univdiagonal,high_univdiagonal,
			       known_univdiags_alloc,known_nunivdiags,
			       
			       query_compress,queryptr,chroffset,chrhigh,
			       univdiags_alloc,diagonals_alloc,localdb_alloc,
			       stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
			       medial_splicesite_prob,distal_splicesite_prob,
			       plusp,genestrand,localdb_nmismatches_allowed,
			       search_localdb_p,innerp,max_npartners,/*inwardp*/true,
			       /*medial_salvagep*/true);
      } else {
	k_outward = solve_trim3(&min_nmismatches,&(*max_prob),&max_medial_prob,&max_medial_setp,
				&(*partnerp),&boundedp,&nsites,this_outward,k_outward,sitep,
				/*Medial_prob_fcn*/Maxent_hr_donor_prob,
				/*Distal_prob_fcn*/Maxent_hr_acceptor_prob,
			
				best_genomicpos,pos3,best_splice_qpos,qstart,querylength,best_mismatchi,
				univdiagonal,low_univdiagonal,high_univdiagonal,
				known_univdiags_alloc,known_nunivdiags,
			
				query_compress,queryptr,chroffset,chrhigh,
				univdiags_alloc,diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,univdiagpool,
				medial_splicesite_prob,distal_splicesite_prob,
				plusp,genestrand,localdb_nmismatches_allowed,
				search_localdb_p,innerp,max_npartners,/*inwardp*/false,
				/*medial_salvagep*/true);
      }
    }

    if (knownsplicing != NULL && (n_inward > 0 || n_outward > 0)) {
      Univdiag_free_array(&known_univdiags_alloc);
    }

    FREE(sitep);
  }

  debug1(printf("\n"));

  /* Combine inward and outward */
  debug1(printf("k_inward = %d, k_outward = %d, partnerp %d, max_prob %f, max_medial_prob %f\n",
		k_inward,k_outward,*partnerp,*max_prob,max_medial_prob));
  Spliceends_combine(this_outward,k_outward,this_inward,k_inward);
  this_outward->boundedp = boundedp;

  if ((k = k_inward + k_outward) == 0) {
    /* Skip */
    debug1(printf("k is zero\n"));
    nspliceends = 0;

  } else if (*partnerp == false) {
    /* min_nmismatches never set */
    debug1(printf("No partner, so use prob only\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_probs[i] > max_medial_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	/* No partner */
	/* Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]); */

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else if (innerp == false) {
    /* Keep only the best spliceends */
    debug1(printf("outer, so keeping only the best\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->medial_nmismatches[i] + this_outward->distal_univdiags[i]->nmismatches + (querylength - this_outward->distal_univdiags[i]->qend) <= min_nmismatches &&
	  this_outward->medial_probs[i] + this_outward->distal_probs[i] > *max_prob - OUTER_PROB_SLOP) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];
	
	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }

  } else {
    /* Keep all partner ends for later resolving of ambiguity */
    debug1(printf("inner with partners, so keeping all\n"));
    for (i = 0; i < k; i++) {
      if (this_outward->partners[i] != 0) {
	this_outward->splice_qpos[nspliceends] = this_outward->splice_qpos[i];
	this_outward->medial_nmismatches[nspliceends] = this_outward->medial_nmismatches[i];
	this_outward->partners[nspliceends] = this_outward->partners[i];

	/* this_outward->matchlengths[nspliceends] = this_outward->matchlengths[i]; */
	/* this_outward->distal_lengths[nspliceends] = this_outward->distal_lengths[i]; */
	/* this_outward->distal_trimpos[nspliceends] = this_outward->distal_trimpos[i]; */
	/* this_outward->distal_nmismatches[nspliceends] = this_outward->distal_nmismatches[i]; */

	Univdiag_transfer(this_outward->distal_univdiags[nspliceends],this_outward->distal_univdiags[i]);

	this_outward->medial_probs[nspliceends] = this_outward->medial_probs[i];
	this_outward->distal_probs[nspliceends] = this_outward->distal_probs[i];
	this_outward->medial_support[nspliceends] = this_outward->medial_support[i];
	this_outward->distal_support[nspliceends] = this_outward->distal_support[i];
	nspliceends++;
      }
    }
  }

#ifdef DEBUG1
  printf("spliceends_trim3_anti yielded %d spliceends\n",nspliceends);
  for (k = 0; k < nspliceends; k++) {
    printf("%u %u %d\n",this_outward->partners[k],this_outward->partners[k] - chroffset,this_outward->splice_qpos[k]);
  }
#endif

  return nspliceends;
}


static T
trim_3 (bool *partnerp, Compress_T query_compress, char *queryptr,
	Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
	Knownsplicing_T knownsplicing, int try_sensedir,
	Univcoord_T univdiagonal, int querylength,
	int qstart, int qend, int pos3, int exon_origin,
	Univcoord_T chroffset, Univcoord_T chrhigh,
	int *mismatch_positions, int total_nmismatches,
	Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
	Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
	int max_nconsecutive, bool plusp, int genestrand, int localdb_nmismatches_allowed,
	bool search_localdb_p, bool innerp, bool salvagep) {
  T new, new_inward;
  int found_sensedir;
  int qdistal;
  int mismatchi;
  int nspliceends = 0;
  Univcoord_T outbound_genomicpos, middle_genomicpos, inbound_genomicpos;
  double max_sense_prob = 0.0, max_antisense_prob = 0.0;
  double medial_splicesite_prob, distal_splicesite_prob;
  Univcoord_T left = univdiagonal - querylength;


  debug1(printf("\nEntered trim_3 with try_sensedir %d, qstart %d to qend %d\n",
		 try_sensedir,qstart,qend));
  assert(try_sensedir == SENSE_FORWARD || try_sensedir == SENSE_ANTI);

  if (salvagep == true) {
    medial_splicesite_prob = SALVAGE_MEDIAL_SPLICESITE_PROB;
    distal_splicesite_prob = SALVAGE_DISTAL_SPLICESITE_PROB;
  } else {
    medial_splicesite_prob = DEFAULT_MEDIAL_SPLICESITE_PROB;
    distal_splicesite_prob = DEFAULT_DISTAL_SPLICESITE_PROB;
  }

  /* Search from qend - END_SPLICESITE_SEARCH, but not past qstart, to
     qend + END_SPLICESITE_SEARCH_MM */
  middle_genomicpos = left + qend;

  debug1(
	 printf("%d mismatches from qstart %d up to pos3 %d:",total_nmismatches,qstart,pos3);
	 for (int i = 0; i <= total_nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );


  /* middle_genomicpos is where the current trim is.  Previously
     added END_SPLICESITE_SEARCH_MM, but now allowing some number
     of mismatches distally to get to outbound_genomicpos */
  mismatchi = 0;
  while (mismatchi < total_nmismatches && mismatch_positions[mismatchi] < qend) {
    mismatchi++;
  }

  if (mismatchi + distal_nmismatches_allowed >= total_nmismatches) {
    qdistal = pos3;
    debug1(printf("qdistal = pos3 %d\n",qdistal));
  } else {
    qdistal = mismatch_positions[mismatchi + distal_nmismatches_allowed];
    debug1(printf("qdistal = mismatch position %d\n",qdistal));
  }
  if (qdistal == querylength) {
    /* Don't need to evaluate at the end of the read */
    qdistal = querylength - 1;
  }

#if 0
  if ((outbound_genomicpos = middle_genomicpos + END_SPLICESITE_SEARCH_MM) > left + querylength) {
    outbound_genomicpos = left + querylength;
  }
#else
  outbound_genomicpos = left + qdistal;
#endif
  if (middle_genomicpos > outbound_genomicpos) {
    middle_genomicpos = outbound_genomicpos;
  }

  if (middle_genomicpos < left + qstart + END_SPLICESITE_SEARCH) {
    inbound_genomicpos = left + qstart;
  } else {
    inbound_genomicpos = middle_genomicpos - END_SPLICESITE_SEARCH;
  }

#ifdef TRIM_AT_GENOME_BOUNDS
  if (left + qstart + MIN_EXON_LENGTH >= genomelength) {
    /* At end of genome, so don't add MIN_EXON_LENGTH */
  } else if (inbound_genomicpos < left + exon_origin + MIN_EXON_LENGTH) {
    inbound_genomicpos = left + exon_origin + MIN_EXON_LENGTH;
  }
#else
  if (inbound_genomicpos < left + exon_origin + MIN_EXON_LENGTH) {
    inbound_genomicpos = left + exon_origin + MIN_EXON_LENGTH;
  }
#endif

  debug1(printf("\n1 Set end points for 3' trim to be %u..%u..%u\n",
		 outbound_genomicpos,middle_genomicpos,inbound_genomicpos));

  if (outbound_genomicpos <= inbound_genomicpos) {
    /* Skip */
    debug1(printf("Got no spliceends\n"));
    return (T) NULL;
  } else {
#if 0
    new = Spliceends_new(/*id*/0,querylength,vectorpool,spliceendspool);
    new_inward = Spliceends_new(/*id*/0,querylength,vectorpool,spliceendspool);
#else
    new = Spliceendsgen_checkout(spliceendsgen,querylength,vectorpool);
    new_inward = Spliceendsgen_checkout(spliceendsgen,querylength,vectorpool);
#endif
  }

  if (stage1 != NULL && stage1->all_oligos_gen_filledp == false) {
    Stage1_fill_all_oligos_gen(stage1,querylength,genestrand);
  }

  if (try_sensedir == SENSE_FORWARD) {
    nspliceends = spliceends_trim3_sense(&max_sense_prob,&(*partnerp),new,new_inward,qstart,
					 inbound_genomicpos,outbound_genomicpos,middle_genomicpos,
					 mismatch_positions,total_nmismatches,univdiagonal,querylength,
					 query_compress,queryptr,chroffset,chrhigh,
					 knownsplicing,max_nconsecutive,
					 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					 stage1,streamspace_max_alloc,streamspace_alloc,
					 univdiagpool,medial_splicesite_prob,distal_splicesite_prob,
					 plusp,genestrand,localdb_nmismatches_allowed,
					 search_localdb_p,innerp);
    found_sensedir = SENSE_FORWARD;

  } else if (try_sensedir == SENSE_ANTI) {
    nspliceends = spliceends_trim3_anti(&max_antisense_prob,&(*partnerp),new,new_inward,qstart,
					inbound_genomicpos,outbound_genomicpos,middle_genomicpos,
					mismatch_positions,total_nmismatches,univdiagonal,querylength,
					query_compress,queryptr,chroffset,chrhigh,
					knownsplicing,max_nconsecutive,
					novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
					stage1,streamspace_max_alloc,streamspace_alloc,
					univdiagpool,medial_splicesite_prob,distal_splicesite_prob,
					plusp,genestrand,localdb_nmismatches_allowed,
					search_localdb_p,innerp);
    found_sensedir = SENSE_ANTI;

  } else {
    /* SENSE_NULL */
    fprintf(stderr,"try_sensedir is neither SENSE_FORWARD nor SENSE_ANTI\n");
    abort();
  }


  Spliceendsgen_return(spliceendsgen,&new_inward);
  if (nspliceends == 0) {
    debug1(printf("Got no spliceends\n"));
    /* Spliceends_free(&new,spliceendspool); */
    Spliceendsgen_return(spliceendsgen,&new);
    return (T) NULL;

  } else {
    new->nspliceends = nspliceends;
    new->sensedir = found_sensedir;
    if (found_sensedir == SENSE_FORWARD) {
      if (plusp) {
	new->splicetype = DONOR;
      } else {
	new->splicetype = ANTIACCEPTOR;
      }
    } else {
      /* SENSE_ANTI */
      if (plusp) {
	new->splicetype = ANTIACCEPTOR;
      } else {
	new->splicetype = DONOR;
      }
    }

#ifdef DEBUG1
    printf("trim_3 got %d spliceends\n",nspliceends);
    for (int i = 0; i < nspliceends; i++) {
      if (new->partners[i] == 0) {
	printf("%d (no partner)\n",new->splice_qpos[i]);
      } else {
	printf("%u %d..%d\n",new->distal_univdiags[i]->univdiagonal,new->distal_univdiags[i]->qstart,new->distal_univdiags[i]->qend);
      }
    }
#endif

    /* Sorted primarily by medial_qpos, and secondarily by ascending distal positions, which Path_solve depends on */
    return new;
  }
}



/* Taken from Substring_trim_qstart_nosplice */
/* Want to return the most distal good region */

#if 0
/* Assumes that caller has done */
total_nmismatches = Genomebits_mismatches_fromright_for_trim(mismatch_positions,/*max_mismatches*/alignlength,
							     /*ome*/bits,/*ome_alt*/bits_alt,query_compress,
							     univdiagonal,querylength,pos5,pos3,plusp,genestrand);
#endif

/* Called only by Splice_resolve and Splice_nomiddle */
int
Spliceends_trim_qstart_nosplice (int *nmismatches_to_trimpos, int *mismatch_positions, int total_nmismatches,
				 int pos5, int pos3) {
  int max_score, score;
  int trimpos = pos5, pos, prevpos, i;
  bool donep = false;

  debug6(printf("Entered trim_qstart_nosplice with pos5 %d, pos3 %d, mismatch_scores %d/%d, match_score %d\n",
		pos5,pos3,TRIM_MISMATCH_SCORE_LAST,TRIM_MISMATCH_SCORE_MULT,TRIM_MATCH_SCORE));
  debug6(printf("%d mismatches:",total_nmismatches));
  debug6(
	 for (i = 0; i <= total_nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );

  if (allow_soft_clips_p == false) {
    /* Report mismatches and do not soft clip */
    *nmismatches_to_trimpos = -1;
    debug6(printf("Returning 0\n"));
    return 0;
  } else if (total_nmismatches == 0) {
    *nmismatches_to_trimpos = 0;
    debug6(printf("Returning %d\n",pos5));
    return pos5;
  }

  /* pos3 | mismatch_positions | (pos5 - 1) */
  prevpos = pos3;
  trimpos = pos = mismatch_positions[0];
  /* Don't add mismatch initially because we stop before the mismatch */
  max_score = score = (prevpos - pos - 1)*TRIM_MATCH_SCORE /*+ TRIM_MISMATCH_SCORE_MULT*/;
  *nmismatches_to_trimpos = 0;
  debug6(printf("Initialize trimpos to be %d with 0 nmismatches and score %d\n",trimpos,score));
  prevpos = pos;

  i = 1;
  while (donep == false && i < total_nmismatches) {
    pos = mismatch_positions[i];
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
    debug6(printf("pos %d, score %d",pos,score));
    if (score >= max_score) {
      debug6(printf(" **"));
      trimpos = pos;
      *nmismatches_to_trimpos = i;
      max_score = score;
#if 0
    } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
      /* No redemption possible with local stretch */
      debug6(printf(" redemption: %d => terminate",pos + 1 - pos5));
      donep = true;
#endif
    } else if (score < 0) {
      /* Want local stretch */
      score = 0;
    }
    debug6(printf("\n"));
    prevpos = pos;
    i++;
  }
    
  if (donep == true) {
    /* No further computation */

  } else if (*nmismatches_to_trimpos == total_nmismatches - 1) {
    /* If last mismatch compensated for previous, then take the last
       segment, regardless of whether it compensates for the last
       mismatch */
    debug6(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,total_nmismatches));
    trimpos = pos5 - 1;
    *nmismatches_to_trimpos += 1;

  } else {
    /* See if last segment compensates */
    pos = pos5 - 1;
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
    debug6(printf("pos %d, score %d",pos,score));
    if (score >= max_score) {
      debug6(printf(" **"));
      trimpos = pos;
      *nmismatches_to_trimpos = i;
      /* max_score = score; */
#if 0
    } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
      /* No redemption needed for final iteration */
      debug6(printf(" redemption: %d => terminate",pos + 1 - pos5));
      donep = true;
#endif
    }
    debug6(printf("\n"));
    /* prevpos = pos; */
  }

  debug6(printf("Returning %d\n",trimpos + 1));
  return trimpos + 1;		/* One position after the mismatch for qstart */

#if 0
  if (trimpos != 0) { /* trimpos + 1 != 1 */
    debug6(printf("Final qstart pos %d => trimpos %d, nmismatches_to_trimpos %d\n",pos,trimpos+1,*nmismatches_to_trimpos));
    return trimpos + 1;		/* One position after the mismatch */
  } else {
    /* For DNA-seq or RNA-seq, if still within chromosome, accept the initial mismatch at the beginning of the read */
    debug6(printf("Backing up 1 bp to start => trimpos %d, nmismatches_to_trimpos %d+1\n",0,*nmismatches_to_trimpos));
    *nmismatches_to_trimpos += 1;
    return 0;			/* trimpos + 1 - 1 */
  }
#endif
}



Univcoord_T
Spliceends_indel_qstart (int nosplice_trimpos, 
			 Univcoord_T univdiagonal, int querylength,
			 Univcoord_T chroffset, Univcoord_T chrhigh,
			 bool plusp, int genestrand, int localdb_nmismatches_allowed,
			 Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			 unsigned short *localdb_alloc, Stage1_T stage1,
			 int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			 Univdiagpool_T univdiagpool, Compress_T query_compress, char *queryptr) {

  Univcoord_T indel_univdiagonal;
  Univdiag_T *univdiags, indel_univdiag;
  int nunivdiags;

  if (stage1->all_oligos_gen_filledp == false) {
    Stage1_fill_all_oligos_gen(stage1,querylength,genestrand);
  }

  debug9(printf("Entering Spliceends_indel_qstart\n"));
  if ((nunivdiags = extend_trim5(&univdiags,/*qpos*/nosplice_trimpos,
				 stage1,univdiagonal,querylength,plusp,
				 /*slop*/max_deletionlen,
				 /*insertion_slop*/max_insertionlen,univdiagpool)) > 0) {
    indel_univdiagonal = univdiags[nunivdiags - 1]->univdiagonal;
    Univdiag_free_array(&univdiags);
    debug9(printf("Returning %u based on indexdb\n",indel_univdiagonal));
    return indel_univdiagonal;

  } else if (localdb != NULL &&
	     (indel_univdiag =
	      novel_trim5_indel(/*outbound_genomicpos*/univdiagonal - querylength + nosplice_trimpos,
				univdiagonal,querylength,query_compress,queryptr,chroffset,chrhigh,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,
				plusp,genestrand,localdb_nmismatches_allowed)) != (Univdiag_T) NULL) {
    debug9(printf("Returning %u based on localdb\n",indel_univdiag->univdiagonal));
    return indel_univdiag->univdiagonal;
    
  } else {
    debug9(printf("Returning 0\n"));
    return 0;
  }
}


/* Returns number of spliceends found.  If a partner is found, then returns spliceendsgen */
/* TODO: No need to return nosplice_trimpos or nosplice_nmismatches */
/* find_splices_p can be false when we are solving for a local mate */
int
Spliceends_trimmed_qstarts (T *new, int *nosplice_trimpos, int *farsplice_trimpos,
			    int *nosplice_nmismatches, int *farsplice_nmismatches,
			    bool *splice5p, Splicetype_T *splicetype5, double *ambig_prob_5,
			    int try_sensedir, Univcoord_T univdiagonal, int querylength,
			    int qend, int exon_origin,
			    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			    bool plusp, int genestrand, int localdb_nmismatches_allowed, bool innerp, bool salvagep,
			    int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
			    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			    Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			    Knownsplicing_T knownsplicing, Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
			    Spliceendsgen_T spliceendsgen, Compress_T query_compress, char *queryptr,
			    Genomebits_T genomebits, Genomebits_T genomebits_alt, bool find_splices_p) {
  int pos5;
  int pos5_nmismatches;
  bool partnerp;
  bool search_localdb_p;

  search_localdb_p = (localdb == NULL) ? false : true;

  debug8(printf("\n***Entered Spliceends_trimmed_qstarts with univdiagonal %u [%u], qend %d..%d, plusp %d, try_sensedir %d, salvagep %d\n",
		univdiagonal,univdiagonal-chroffset,0,qend,plusp,try_sensedir,salvagep));

  /* left = univdiagonal - (Univcoord_T) querylength; */

  *new = (T) NULL;
  *farsplice_trimpos = -1;
  *farsplice_nmismatches = 0;
  *splice5p = false; *splicetype5 = NO_SPLICE; *ambig_prob_5 = 0.0;

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos5 = (univdiagonal >= chroffset + (Univcoord_T) querylength) ? 0 : (int) (chroffset - left);
#else
  pos5 = 0;
#endif

  if (pos5 >= qend) {
    debug8(printf("trimpos %d >= qend %d, so returning -1\n",pos5,qend));
    /* No spliceends.  No trimming. */
    *nosplice_trimpos = pos5;
    return 0;

  } else {
    /* Note: This procedure allows for 1 mismatch at the start of the read, so we don't need to handle it here */
    *nosplice_trimpos = Genomebits_trim_qstart(&(*nosplice_nmismatches),query_compress,
					       /*bits*/genomebits,univdiagonal,querylength,
					       pos5,/*pos3*/qend,plusp,genestrand);
    debug8(printf("nosplice: trimpos %d (relative to %d)\n",*nosplice_trimpos,pos5));
  }

  if (*nosplice_trimpos >= qend) {
    debug8(printf("trimpos %d >= qend %d, so returning -1\n",*nosplice_trimpos,qend));
    /* No spliceends.  No trimming. */
    return 0;

  } else if (splicingp == false || find_splices_p == false) {
    /* No spliceends.  Found a trim. */
    debug8(printf("Keeping given trim because find_splices_p is false\n"));
    return 1;

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* No spliceends.  Found a trim. */
    debug8(printf("Keeping given trim because chrnum %d is circular\n",chrnum));
    return 1;
#endif

  } else {
    /* Want this because it is equivalent to having peelback */
    debug8(printf("Finding splice based on trim_5 on trimpos %d to qend %d\n",*nosplice_trimpos,qend));

    /* Previously called by Spliceends_qstart_nosplice */
    pos5_nmismatches =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_alloc,/*max_mismatches*/qend - pos5,
					       /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					       univdiagonal,querylength,pos5,/*pos3*/qend,plusp,genestrand);
    if ((*new = trim_5(&partnerp,query_compress,queryptr,
		       stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
		       univdiagonal,querylength,pos5,/*qstart*/*nosplice_trimpos,qend,
		       exon_origin,chroffset,chrhigh,
		       mismatch_positions_alloc,/*total_nmismatches*/pos5_nmismatches,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       univdiagpool,vectorpool,spliceendsgen,
		       MAX_NCONSECUTIVE_FAR,plusp,genestrand,localdb_nmismatches_allowed,
		       search_localdb_p,innerp,salvagep)) == NULL) {
      /* No spliceends.  Found a trim. */
      debug8(printf("spliceends is NULL\n"));
      return 1;

    } else {
      *farsplice_trimpos = (*new)->splice_qpos[0];
      *farsplice_nmismatches = (*new)->medial_nmismatches[0];

      *splice5p = true; *splicetype5 = (*new)->splicetype; *ambig_prob_5 = (*new)->medial_probs[0];
      if (partnerp == true) {
	debug8(printf("found spliceends with partner\n"));
	return (*new)->nspliceends;
      } else {
	debug8(printf("found spliceends without partner\n"));
	/* Spliceends_free(&(*new),spliceendspool); */
	Spliceendsgen_return(spliceendsgen,&(*new));
	assert(*new == NULL);
	return 1;
      }
    }
  }
}


/* TODO: Find the best splice site near the nosplice trimpos,
   regardless of probability.  Then these splice sites can be used in
   two-pass mode to use this site, regardless of probability */

/* Cannot return -1 for trimpos, because caller needs to use its value */
/* find_splices_p can be false when we are solving for a local mate */
bool
Spliceends_qstart_trim (int *trimpos, int *nmismatches_to_trimpos,
			int *found_sensedir, Splicetype_T *splicetype, double *ambig_prob_qstart,
			Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			Knownsplicing_T knownsplicing, int try_sensedir,
			Univcoord_T univdiagonal, int querylength, int pos3, int exon_origin,
			Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			bool plusp, int genestrand, int *mismatch_positions_alloc,

			Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			unsigned short *localdb_alloc, 

			Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
			Compress_T query_compress, char *queryptr,
			Genomebits_T genomebits, Genomebits_T genomebits_alt,
			bool find_splices_p) {
  T spliceends;
  int pos5;
  int total_nmismatches;
  bool partnerp;


  *found_sensedir = try_sensedir;
  /* left = univdiagonal - (Univcoord_T) querylength; */

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos5 = (univdiagonal >= chroffset + (Univcoord_T) querylength) ? 0 : (int) (chroffset - left);
#else
  pos5 = 0;
#endif
  debug8(printf("\n***Entered Spliceends_qstart_trim with univdiagonal %u, pos5 %d, pos3 %d, plusp %d, try_sensedir %d\n",
		univdiagonal,pos5,pos3,plusp,try_sensedir));

  if (pos5 >= pos3) {
    debug8(printf("pos5 %d >= pos3 %d, so returning pos3\n",pos5,pos3));
    *trimpos = pos3;
    *splicetype = NO_SPLICE;
    *ambig_prob_qstart = 0.0;
    return false;

  } else {
    /* Note: This procedure allows for 1 mismatch at the start of the read, so we don't need to handle it here */
    *trimpos = Genomebits_trim_qstart(&(*nmismatches_to_trimpos),query_compress,
				      /*bits*/genomebits,univdiagonal,querylength,
				      pos5,pos3,plusp,genestrand);
    debug8(printf("nosplice trimpos %d (relative to %d)\n",*trimpos,pos5));
  }

  if (*trimpos >= pos3) {
    debug8(printf("trimpos %d >= pos3 %d, so returning pos3\n",*trimpos,pos3));
    *trimpos = pos3;
    *splicetype = NO_SPLICE;
    *ambig_prob_qstart = 0.0;
    return false;

  } else if (*trimpos == pos5) {
    /* Use trimpos, which extends to the start, which could be longer than the caller had */
    debug8(printf("Using trimpos %d, which extends to the start, pos5 %d\n",*trimpos,pos5));
    *splicetype = NO_SPLICE;
    *ambig_prob_qstart = 0.0;
    return false;

  } else if (splicingp == false || find_splices_p == false) {
    /* Keep given trim */
    debug8(printf("Keeping given trim because find_splices_p is false\n"));
    *splicetype = NO_SPLICE;
    *ambig_prob_qstart = 0.0;
    return false;

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* Keep given trim */
    debug8(printf("Keeping given trim because chrnum %d is circular\n",chrnum));
    *splicetype = NO_SPLICE;
    *ambig_prob_qstart = 0.0;
    return false;
#endif

#if 0
  } else if (0 && *trimpos <= pos5 + ACCEPTABLE_TRIM) {
    debug8(printf("Accepting splice based on trim_5 on pos5 %d to pos3 %d\n",pos5,pos3));
    if ((spliceends = trim_5(&partnerp,query_compress,queryptr,
			     stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
			     unifdiagonal,querylength,pos5,qstart,/*qend*/pos3,
			     exon_origin,chroffset,chrhigh,
			     mismatch_positions_alloc,total_nmismatches,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     univdiagpool,vectorpool,spliceendsgen,
			     MAX_NCONSECUTIVE_CLOSE,plusp,genestrand,
			     /*localdb_nmismatches_allowed*/querylength,
			     /*search_localdb_p*/false,/*innerp*/false,/*salvagep*/false)) != NULL) {
      /* TODO: Make sure this is the farthest one */
      *found_sensedir = spliceends->sensedir;
      *trimpos = spliceends->splice_qpos[0];
      *nmismatches_to_trimpos = spliceends->medial_nmismatches[0];
      *splicetype = spliceends->splicetype;
      *ambig_prob_qstart = spliceends->medial_probs[0];
      /* Spliceends_free(&spliceends,spliceendspool); */
      Spliceendsgen_return(spliceendsgen,&spliceends);
      debug8(printf("found spliceends at trimpos %d with %d nmismatches\n",*trimpos,*nmismatches_to_trimpos));
      return true;

    } else {
      debug8(printf("spliceends is NULL, but multiple mismatches at the end => trimpos %d\n",*trimpos));
      *splicetype = NO_SPLICE;
      *ambig_prob_qstart = 0.0;
      return false;
    }
#endif
    
  } else {
    debug8(printf("Finding splice based on trim_5 on trimpos %d to pos3 %d\n",*trimpos,pos3));

    /* Previously called by Spliceends_qstart_nosplice */
    total_nmismatches =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_alloc,/*max_mismatches*/pos3 - pos5,
					       /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					       univdiagonal,querylength,pos5,pos3,plusp,genestrand);
    if ((spliceends = trim_5(&partnerp,query_compress,queryptr,
			     stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
			     univdiagonal,querylength,pos5,/*qstart*/(*trimpos),/*qend*/pos3,
			     exon_origin,chroffset,chrhigh,
			     mismatch_positions_alloc,total_nmismatches,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     univdiagpool,vectorpool,spliceendsgen,
			     MAX_NCONSECUTIVE_FAR,plusp,genestrand,
			     /*localdb_nmismatches_allowed*/querylength,
			     /*search_localdb_p*/false,/*innerp*/false,/*salvagep*/false)) != NULL) {
      /* TODO: Make sure this is the farthest one */
      *found_sensedir = spliceends->sensedir;
      *trimpos = spliceends->splice_qpos[0];
      *nmismatches_to_trimpos = spliceends->medial_nmismatches[0];
      *splicetype = spliceends->splicetype;
      *ambig_prob_qstart = spliceends->medial_probs[0];
      /* Spliceends_free(&spliceends,spliceendspool); */
      Spliceendsgen_return(spliceendsgen,&spliceends);
      debug8(printf("found spliceends at trimpos %d with %d nmismatches\n",*trimpos,*nmismatches_to_trimpos));
      return true;
      
    } else {
      /* Keep given trim */
      debug8(printf("spliceends is NULL\n"));
      *splicetype = NO_SPLICE;
      *ambig_prob_qstart = 0.0;
      return false;
    }
  }
}


/* Taken from Substring_trim_qend_nosplice */
/* Want to return the most distal good region */

#if 0
/* Assumes that caller has done */
total_nmismatches = Genomebits_mismatches_fromleft_for_trim(mismatch_positions,/*max_mismatches*/alignlength,
							    /*ome*/bits,/*ome_alt*/bits_alt,query_compress,
							    univdiagonal,querylength,pos5,pos3,plusp,genestrand);
#endif

/* Called only by Splice_resolve and Splice_nomiddle */
int
Spliceends_trim_qend_nosplice (int *nmismatches_to_trimpos, int *mismatch_positions, int total_nmismatches,
			       int pos5, int pos3, int querylength) {
  int max_score, score;
  int trimpos = pos3, pos, prevpos, i;
  bool donep = false;

  debug6(printf("Entered trim_qend_nosplice with pos5 %d, pos3 %d, mismatch_scores %d/%d, match_score %d\n",
		pos5,pos3,TRIM_MISMATCH_SCORE_LAST,TRIM_MISMATCH_SCORE_MULT,TRIM_MATCH_SCORE));
  debug6(printf("%d mismatches:",total_nmismatches));
  debug6(
	 for (i = 0; i <= total_nmismatches; i++) {
	   printf(" %d",mismatch_positions[i]);
	 }
	 printf("\n");
	 );


  if (allow_soft_clips_p == false) {
    /* Report mismatches and do not soft clip */
    *nmismatches_to_trimpos = -1;
    debug6(printf("Returning %d\n",querylength));
    return querylength;
  } else if (total_nmismatches == 0) {
    *nmismatches_to_trimpos = 0;
    debug6(printf("Returning %d\n",pos3));
    return pos3;
  }

  /* (pos5 - 1) | mismatch_positions | pos3 */
  prevpos = pos5 - 1;
  trimpos = pos = mismatch_positions[0];
  /* Don't add mismatch initially because we stop before the mismatch */
  max_score = score = (pos - prevpos - 1)*TRIM_MATCH_SCORE /*+ TRIM_MISMATCH_SCORE_MULT*/;
  *nmismatches_to_trimpos = 0;
  debug6(printf("Initialize trimpos to be %d with 0 nmismatches and score %d\n",trimpos,score));
  prevpos = pos;

  i = 1;
  while (donep == false && i < total_nmismatches) {
    pos = mismatch_positions[i];
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug6(printf("pos %d, score %d",pos,score));
    if (score >= max_score) {
      debug6(printf(" **"));
      trimpos = pos;
      *nmismatches_to_trimpos = i;
      max_score = score;
#if 0
    } else if (score + /*redemption*/(pos3 - pos) < 0) {
      /* No redemption possible with local stretch */
      debug6(printf(" redemption: %d => terminate",pos3 - pos));
      donep = true;
#endif
    } else if (score < 0) {
      /* Want local stretch */
      score = 0;
    }
    debug6(printf("\n"));
    prevpos = pos;
    i++;
  }
  
  if (donep == true) {
    /* No further computation */

  } else if (*nmismatches_to_trimpos == total_nmismatches - 1) {
    /* If last mismatch compensated for previous, then take the last
       segment, regardless of whether it compensates for the last
       mismatch */
    debug6(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,total_nmismatches));
    trimpos = pos3;
    *nmismatches_to_trimpos += 1;

  } else {
    /* See if last segment compensates */
    pos = pos3;
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug6(printf("pos %d, score %d",pos,score));
    if (score >= max_score) {
      debug6(printf(" **"));
      trimpos = pos;
      *nmismatches_to_trimpos = i;
      /* max_score = score; */
#if 0
    } else if (score + /*redemption*/(pos3 - pos) < 0) {
      /* No redemption needed for final iteration */
      debug6(printf(" redemption: %d => terminate",pos3 - pos));
      donep = true;
#endif
    }
    debug6(printf("\n"));
    /* prevpos = pos; */
  }

  debug6(printf("Returning %d\n",trimpos));
  return trimpos;		/* At the mismatch for qend */
  
#if 0
  if (trimpos != querylength - 1) {
    debug6(printf("Final qend pos %d => trimpos %d, nmismatches_to_trimpos %d\n",pos,trimpos,*nmismatches_to_trimpos));
    return trimpos;		/* qend is outside the region */
  } else {
    /* For DNA-seq or RNA-seq, if still within the chromosome, accept the final mismatch at end of read */
    debug6(printf("Advancing 1 bp to end => trimpos %d, nmismatches_to_trimpos %d+1\n",querylength,*nmismatches_to_trimpos));
    *nmismatches_to_trimpos += 1;
    return querylength;
  }
#endif
}


Univcoord_T
Spliceends_indel_qend (int nosplice_trimpos, 
		       Univcoord_T univdiagonal, int querylength,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       bool plusp, int genestrand, int localdb_nmismatches_allowed,
		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc, Stage1_T stage1,
		       int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       Univdiagpool_T univdiagpool, Compress_T query_compress, char *queryptr) {

  Univcoord_T indel_univdiagonal;
  Univdiag_T *univdiags, indel_univdiag;
  int nunivdiags;

  if (stage1->all_oligos_gen_filledp == false) {
    Stage1_fill_all_oligos_gen(stage1,querylength,genestrand);
  }

  debug9(printf("Entering Spliceends_indel_qend\n"));
  if ((nunivdiags = extend_trim3(&univdiags,/*qpos*/nosplice_trimpos,
				 stage1,univdiagonal,querylength,plusp,
				 /*slop*/max_deletionlen,
				 /*insertion_slop*/max_insertionlen,univdiagpool)) > 0) {
    indel_univdiagonal = univdiags[0]->univdiagonal;
    Univdiag_free_array(&univdiags);
    debug9(printf("Returning %u based on indexdb\n",indel_univdiagonal));
    return indel_univdiagonal;

  } else if (localdb != NULL &&
	     (indel_univdiag =
	      novel_trim3_indel(/*outbound_genomicpos*/univdiagonal - querylength + nosplice_trimpos,
				univdiagonal,querylength,query_compress,queryptr,chroffset,chrhigh,
				novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				stage1,streamspace_max_alloc,streamspace_alloc,
				plusp,genestrand,localdb_nmismatches_allowed)) != (Univdiag_T) NULL) {
    debug9(printf("Returning %u based on localdb\n",indel_univdiag->univdiagonal));
    return indel_univdiag->univdiagonal;
      
  } else {
    debug9(printf("Returning 0\n"));
    return 0;
  }
}


/* Returns number of spliceends found.  If a partner is found, then returns spliceendsgen */
/* find_splices_p can be false when we are solving for a local mate */
int
Spliceends_trimmed_qends (T *new, int *nosplice_trimpos, int *farsplice_trimpos,
			  int *nosplice_nmismatches, int *farsplice_nmismatches,
			  bool *splice3p, Splicetype_T *splicetype3, double *ambig_prob_3,
			  int try_sensedir, Univcoord_T univdiagonal, int querylength,
			  int qstart, int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  bool plusp, int genestrand, int localdb_nmismatches_allowed, bool innerp, bool salvagep,
			  int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
			  Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			  Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			  Knownsplicing_T knownsplicing, Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
			  Spliceendsgen_T spliceendsgen, Compress_T query_compress, char *queryptr,
			  Genomebits_T genomebits, Genomebits_T genomebits_alt, bool find_splices_p) {
  int pos3;
  int pos3_nmismatches;
  bool partnerp;
  bool search_localdb_p;

  search_localdb_p = (localdb == NULL) ? false : true;

  debug8(printf("\n***Entered Spliceends_trimmed_qends with univdiagonal %u [%u], qstart %d..%d, plusp %d, try_sensedir %d, salvagep %d\n",
		univdiagonal,univdiagonal-chroffset,qstart,querylength,plusp,try_sensedir,salvagep));

  /* left = univdiagonal - (Univcoord_T) querylength; */

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos3 = (univdiagonal <= chrhigh) ? querylength : (int) (chrhigh - left);
#else
  pos3 = querylength;
#endif

  *new = (T) NULL;
  *farsplice_trimpos = -1;
  *farsplice_nmismatches = 0;
  *splice3p = false; *splicetype3 = NO_SPLICE; *ambig_prob_3 = 0.0;

  if (pos3 <= qstart) {
    debug8(printf("pos3 %d <= qstart %d, so returning -1\n",pos3,qstart));
    /* No spliceends.  No trimming. */
    *nosplice_trimpos = pos3;
    return 0;

  } else {
    /* Note: This procedure allows for 1 mismatch at the end of the read, so we don't need to handle it here */
    *nosplice_trimpos = Genomebits_trim_qend(&(*nosplice_nmismatches),query_compress,
					     /*bits*/genomebits,univdiagonal,querylength,
					     /*pos5*/qstart,pos3,plusp,genestrand);
    debug8(printf("nosplice: trimpos %d (relative to %d)\n",*nosplice_trimpos,pos3));
  }

  if (*nosplice_trimpos <= qstart) {
    debug8(printf("trimpos %d <= qstart %d, so returning -1\n",*nosplice_trimpos,qstart));
    /* No spliceends.  No trimming. */
    return 0;

  } else if (splicingp == false || find_splices_p == false) {
    /* No spliceends.  Found a trim. */
    debug8(printf("Keeping given trim because find_splices_p is false\n"));
    return 1;

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* No spliceends.  Found a trim. */
    debug8(printf("Keeping given trim because chrnum %d is circular\n",chrnum));
    return 1;
#endif

  } else {
    /* Want this because it is equivalent to peelback */
    debug8(printf("Finding splice based on trim_3 from qstart %d to trimpos %d\n",qstart,*nosplice_trimpos));

    /* Previously called by Spliceends_trim_qend_nosplice */
    pos3_nmismatches =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_alloc,/*max_mismatches*/pos3 - qstart,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					      univdiagonal,querylength,/*pos5*/qstart,pos3,plusp,genestrand);
    if ((*new = trim_3(&partnerp,query_compress,queryptr,
		       stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
		       univdiagonal,querylength,
		       qstart,/*qend*/*nosplice_trimpos,pos3,exon_origin,
		       chroffset,chrhigh,mismatch_positions_alloc,/*total_nmismatches*/pos3_nmismatches,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       univdiagpool,vectorpool,spliceendsgen,
		       MAX_NCONSECUTIVE_FAR,plusp,genestrand,localdb_nmismatches_allowed,
		       search_localdb_p,innerp,salvagep)) == NULL) {
      /* No spliceends.  Found a trim. */
      debug8(printf("spliceends is NULL\n"));
      return 1;

    } else {
      *farsplice_trimpos = (*new)->splice_qpos[0];
      *farsplice_nmismatches = (*new)->medial_nmismatches[0];

      *splice3p = true; *splicetype3 = (*new)->splicetype; *ambig_prob_3 = (*new)->medial_probs[0];
      if (partnerp == true) {
	debug8(printf("found spliceends with partner\n"));
	return (*new)->nspliceends;
      } else {
	debug8(printf("found spliceends without partner\n"));
	/* Spliceends_free(&(*new),spliceendspool); */
	Spliceendsgen_return(spliceendsgen,&(*new));
	assert(*new == NULL);
	return 1;
      }
    }
  }
}


/* TODO: Find the best splice site near the nosplice trimpos,
   regardless of probability.  Then these splice sites can be used in
   two-pass mode to use this site, regardless of probability */

/* Cannot return -1 for trimpos, because caller needs to use its value */
/* find_splices_p can be false when we are solving for a local mate */
bool
Spliceends_qend_trim (int *trimpos, int *nmismatches_to_trimpos,
		      int *found_sensedir, Splicetype_T *splicetype, double *ambig_prob_qend,
		      Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		      Knownsplicing_T knownsplicing, int try_sensedir,
		      Univcoord_T univdiagonal, int querylength, int pos5,
		      int exon_origin, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      bool plusp, int genestrand, int *mismatch_positions_alloc,

		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, 

		      Univdiagpool_T univdiagpool, Vectorpool_T vectorpool, Spliceendsgen_T spliceendsgen,
		      Compress_T query_compress, char *queryptr,
		      Genomebits_T genomebits, Genomebits_T genomebits_alt, bool find_splices_p) {
  T spliceends;
  int pos3;
  int total_nmismatches;
  bool partnerp;

  *found_sensedir = try_sensedir;

  /* left = univdiagonal - (Univcoord_T) querylength; */

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  pos3 = (univdiagonal <= chrhigh) ? querylength : (int) (chrhigh - left);
#else
  pos3 = querylength;
#endif
  debug8(printf("\n***Entered Spliceends_qend_trim with univdiagonal %u, pos5 %d, pos3 %d, plusp %d, try_sensedir %d\n",
		univdiagonal,pos5,pos3,plusp,try_sensedir));

  if (pos3 <= pos5) {
    debug8(printf("trimpos %d <= pos5 %d, so returning pos5\n",*trimpos,pos5));
    *trimpos = pos5;
    *splicetype = NO_SPLICE;
    *ambig_prob_qend = 0.0;
    return false;

  } else {
    /* Note: This procedure allows for 1 mismatch at the start of the read, so we don't need to handle it here */
    *trimpos = Genomebits_trim_qend(&(*nmismatches_to_trimpos),query_compress,
				    /*bits*/genomebits,univdiagonal,querylength,
				    pos5,pos3,plusp,genestrand);
    debug8(printf("nosplice trimpos %d (relative to %d)\n",*trimpos,pos3));
  }

  if (*trimpos <= pos5) {
    debug8(printf("trimpos %d <= pos5 %d, so returning -1\n",*trimpos,pos5));
    *trimpos = pos5;
    *splicetype = NO_SPLICE;
    *ambig_prob_qend = 0.0;
    return false;

  } else if (*trimpos == pos3) {
    /* Use trimpos, which extends to the end, which could be longer than the caller had */
    debug8(printf("Using trimpos %d, which extends to the end, pos3 %d\n",*trimpos,pos3));
    *splicetype = NO_SPLICE;
    *ambig_prob_qend = 0.0;
    return false;

  } else if (splicingp == false || find_splices_p == false) {
    /* Keep given trim */
    debug8(printf("Keeping given trim because find_splices_p is false\n"));
    *splicetype = NO_SPLICE;
    *ambig_prob_qend = 0.0;
    return false;

#ifdef DISALLOW_CIRCULAR_SPLICING
  } else if (circularp[chrnum] == true) {
    /* Keep given trim */
    debug8(printf("Keeping given trim because chrnum %d is circular\n",chrnum));
    *splicetype = NO_SPLICE;
    *ambig_prob_qend = 0.0;
    return false;
#endif

#if 0
  } else if (0 && *trimpos >= pos3 - ACCEPTABLE_TRIM) {
    debug8(printf("Accepting splice based on trim_3 from pos5 %d to pos3 %d\n",pos5,pos3));
    if ((spliceends = trim_3(&partnerp,query_compress,queryptr,
			     stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
			     univdiagonal,querylength,/*qstart*/pos5,qend,pos3,exon_origin,
			     chroffset,chrhigh,mismatch_positions_alloc,total_nmismatches,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     univdiagpool,vectorpool,spliceendsgen,
			     MAX_NCONSECUTIVE_CLOSE,plusp,genestrand,
			     /*localdb_nmismatches_allowed*/0,/*search_localdb_p*/false,
			     /*innerp*/false,/*salvagep*/false)) != NULL) {
      /* TODO: Make sure this is the farthest one */
      *found_sensedir = spliceends->sensedir;
      *trimpos = spliceends->splice_qpos[0];
      *nmismatches_to_trimpos = spliceends->medial_nmismatches[0];
      *splicetype = spliceends->splicetype;
      *ambig_prob_qend = spliceends->medial_probs[0];
      /* Spliceends_free(&spliceends,spliceendspool); */
      Spliceendsgen_return(spliceendsgen,&spliceends);
      debug8(printf("found spliceends at trimpos %d and %d nmismatches\n",*trimpos,*nmismatches_to_trimpos));
      return true;

    } else {
      debug8(printf("spliceends is NULL, but multiple mismatches at the end => trimpos %d\n",*trimpos));
      *splicetype = NO_SPLICE;
      *ambig_prob_qend = 0.0;
      return false;
    }
#endif

  } else {
    debug8(printf("Finding splice based on trim_3 from pos5 %d to trimpos %d\n",pos5,*trimpos));

    /* Previously called by Spliceends_trim_qend_nosplice */
    total_nmismatches =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_alloc,/*max_mismatches*/pos3 - pos5,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,query_compress,
					      univdiagonal,querylength,pos5,pos3,plusp,genestrand);
    if ((spliceends = trim_3(&partnerp,query_compress,queryptr,
			     stage1,streamspace_max_alloc,streamspace_alloc,knownsplicing,try_sensedir,
			     univdiagonal,querylength,
			     /*qstart*/pos5,/*qend*/*trimpos,pos3,exon_origin,
			     chroffset,chrhigh,mismatch_positions_alloc,total_nmismatches,
			     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			     univdiagpool,vectorpool,spliceendsgen,
			     MAX_NCONSECUTIVE_FAR,plusp,genestrand,
			     /*localdb_nmismatches_allowed*/0,/*search_localdb_p*/false,
			     /*innerp*/false,/*salvagep*/false)) != NULL) {
      /* TODO: Make sure this is the farthest one */
      *found_sensedir = spliceends->sensedir;
      *trimpos = spliceends->splice_qpos[0];
      *nmismatches_to_trimpos = spliceends->medial_nmismatches[0];
      *splicetype = spliceends->splicetype;
      *ambig_prob_qend = spliceends->medial_probs[0];
      /* Spliceends_free(&spliceends,spliceendspool); */
      Spliceendsgen_return(spliceendsgen,&spliceends);
      debug8(printf("found spliceends at trimpos %d and %d nmismatches\n",*trimpos,*nmismatches_to_trimpos));
      return true;

    } else {
      /* Keep given trim */
      debug8(printf("spliceends is NULL\n"));
      *splicetype = NO_SPLICE;
      *ambig_prob_qend = 0.0;
      return false;
    }
  }
}


static int
univcoord_descending_cmp (const void *x, const void *y) {
  Univcoord_T a = * (Univcoord_T *) x;
  Univcoord_T b = * (Univcoord_T *) y;

  if (a > b) {
    return -1;
  } else if (b > a) {
    return +1;
  } else {
    return 0;
  }
}



/* For resolve procedure, could require *ndiagonals == 1 */
/* Returns diagonals in descending order, same as trim_5 */
Univcoord_T *
Spliceends_qstart_resolve (int *ndiagonals, int *local_nmismatches, int pos3, int querylength,
			   Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			   Compress_T query_compress, char *queryptr, bool plusp, int genestrand,
			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			   Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			   int localdb_nmismatches_allowed) {
  int i;
  int nmatches, matchlength;
  bool sortedp, trimmedp;


  /* Not sure whether to use nmismatches_allowed instead of 0 */
  /* TODO: Use Boyer-Moore instead */
  /* Called only after Path_extend, so we have already tried salvage
     mode over a large region.  Here are relying on a small region
     near the anchor end, so we can require few mismatches */
  /* Want nmismatches_allowed based on the fragment, not the entire
     read.  This value is ignored for salvage mode anyway. */

  /* Setting salvagep to be true results in very slow speed */
  debug3(printf("\nEntered Spliceends_qstart_resolve with querystart 0 to pos3 %d\n",pos3));
  /* assert(high_univdiagonal - low_univdiagonal <= positive_gap_distance); */

  if (localdb == NULL) {
    /* No localdb present */
    *ndiagonals = 0;
    return (Univcoord_T *) NULL;

  } else if ((*ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&(*local_nmismatches),
					novel_univdiags_alloc,localdb,novel_diagonals_alloc,localdb_alloc,
					stage1,streamspace_max_alloc,streamspace_alloc,
					queryptr,/*pos5*/0,pos3,querylength,
					/*anchor_univdiagonal*/high_univdiagonal,
					low_univdiagonal,high_univdiagonal - /*prevent continuation*/1,
					query_compress,plusp,genestrand,genomebits,
					localdb_nmismatches_allowed,/*extend5p*/true,
					/*trim5p*/true,/*trim3p*/false)) == 0) {
    debug3(printf("no diagonals\n"));
    return (Univcoord_T *) NULL;

  } else if ((nmatches = matchlength - (*local_nmismatches)) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*(*local_nmismatches)) {
    debug3(printf("matchlength %d and local_nmismatches %d => Skipping\n",
		  matchlength,*local_nmismatches));
    *ndiagonals = 0;
    return (Univcoord_T *) NULL;

  } else {
    /* Transfer from univdiags to diagonals */
    for (i = 0; i < *ndiagonals; i++) {
      novel_diagonals_alloc[i] = novel_univdiags_alloc[i]->univdiagonal;
    }

    if (*ndiagonals == 1) {
      debug3(printf("Returning single univdiagonal %u with matchlength %d and local_nmismatches %d\n\n",
		    novel_diagonals_alloc[0],matchlength,*local_nmismatches));
      return novel_diagonals_alloc;

    } else {
      /* Highest (last) is closest to pathL */
      debug3(printf("multiple ndiagonals %d with matchlength %d and local_nmismatches %d.  Sort ascending\n",
		    *ndiagonals,matchlength,*local_nmismatches));
      if (sortedp == true) {
	reverse_univcoord_inplace(novel_diagonals_alloc,/*starti*/0,/*endi*/*ndiagonals);
      } else {
	qsort(novel_diagonals_alloc,*ndiagonals,sizeof(Univcoord_T),univcoord_descending_cmp);
      }
      return novel_diagonals_alloc;
    }
  }
}


/* Returns diagonals in ascending order, same as trim_3 */
Univcoord_T *
Spliceends_qend_resolve (int *ndiagonals, int *local_nmismatches, int pos5, int querylength,
			 Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			 Compress_T query_compress, char *queryptr, bool plusp, int genestrand,
			 Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
			 Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
			 int localdb_nmismatches_allowed) {
  int i;
  int nmatches, matchlength;
  bool sortedp, trimmedp;

  /* Not sure whether to use nmismatches_allowed instead of 0 */
  /* TODO: Use Boyer-Moore instead */
  /* Called only after Path_extend, so we have already tried salvage
     mode over a large region.  Here are relying on a small region
     near the anchor end, so we can require few mismatches */
  /* Want nmismatches_allowed based on the fragment, not the entire
     read.  This value is ignored for salvage mode anyway. */

  /* Setting salvagep to be true results in very slow speed */
  debug3(printf("\nEntered Spliceends_qend_resolve with pos5 %d to querylength %d\n",pos5,querylength));
  /* assert(high_univdiagonal - low_univdiagonal <= positive_gap_distance); */

  if (localdb == NULL) {
    /* No localdb present */
    *ndiagonals = 0;
    return (Univcoord_T *) NULL;

  } else if ((*ndiagonals = Localdb_get(&sortedp,&trimmedp,&matchlength,&(*local_nmismatches),
					novel_univdiags_alloc,localdb,novel_diagonals_alloc,localdb_alloc,
					stage1,streamspace_max_alloc,streamspace_alloc,
					queryptr,pos5,/*pos3*/querylength,querylength,
					/*anchor_univdiagonal*/low_univdiagonal,
					low_univdiagonal + /*prevent continuation*/1,high_univdiagonal,
					query_compress,plusp,genestrand,genomebits,
					localdb_nmismatches_allowed,/*extend5p*/false,
					/*trim5p*/false,/*trim3p*/true)) == 0) {
    debug3(printf("no diagonals\n"));
    return (Univcoord_T *) NULL;

  } else if ((nmatches = matchlength - (*local_nmismatches)) < SUFFICIENT_NMATCHES &&
	     nmatches < 3*(*local_nmismatches)) {
    debug3(printf("matchlength %d and local_nmismatches %d => Skipping\n",
		  matchlength,*local_nmismatches));
    *ndiagonals = 0;
    return (Univcoord_T *) NULL;

  } else {
    /* Transfer from univdiags to diagonals */
    for (i = 0; i < *ndiagonals; i++) {
      novel_diagonals_alloc[i] = novel_univdiags_alloc[i]->univdiagonal;
    }
    
    if (*ndiagonals == 1) {
      debug3(printf("Returning single univdiagonal %u with matchlength %d and local_nmismatches %d\n\n",
		    novel_diagonals_alloc[0],matchlength,*local_nmismatches));
      return novel_diagonals_alloc;
      
    } else {
      /* Lowest (last) is closest to pathH */
      debug3(printf("multiple ndiagonals %d with matchlength %d and local_nmismatches %d.  Sort descending\n",
		    *ndiagonals,matchlength,*local_nmismatches));
      if (sortedp == false) {
	qsort(novel_diagonals_alloc,*ndiagonals,sizeof(Univcoord_T),univcoord_ascending_cmp);
      }
      return novel_diagonals_alloc;
    }
  }
}

