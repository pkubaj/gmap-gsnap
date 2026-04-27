static char rcsid[] = "$Id: d1776cb9ba96b349f8139e43b5a0d03e7c04f452 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "kmer-search.h"

#include <string.h>		/* For strlen */
#include "assert.h"
#include "mem.h"
#include "types.h"
#include "chrnum.h"
#include "reader.h"
#include "oligo.h"

#ifdef LARGE_GENOMES
#include "intersect-large.h"
#include "intersect-approx-uint8.h"
#include "intersect-approx-indices-uint8.h"
#include "intersect-indices2-large.h"
#else
#include "intersect-small.h"
#include "intersect-approx-uint4.h"
#include "intersect-approx-indices-uint4.h"
#include "intersect-indices2-small.h"
#endif

#include "transcript.h"
#include "popcount.h"
#include "junction.h"

#include "univdiag.h"
#include "univdiagdef.h"
#include "genomebits_count.h"
#include "genomebits_kmer.h"
#include "splice.h"
#include "indel.h"
#include "intron.h"
#include "maxent_hr.h"
#include "sedgesort.h"
#include "path-solve.h"


#ifndef LARGE_GENOMES
#include "merge-diagonals-simd-uint4.h"
#elif !defined(HAVE_AVX512) && !defined(HAVE_AVX2)
#include "merge-diagonals-heap.h" /* For Merge_diagonals_large */
#include "merge-diagonals-simd-uint4.h"
#else
#include "merge-diagonals-simd-uint8.h" /* For Merge_diagonals_large */
#include "merge-diagonals-simd-uint4.h"
#endif

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


#ifdef LARGE_GENOMES
#define GETPOS(high,low) (((UINT8) high << 32) + low)
#endif


#define MAX_UNIVDIAGS 5


/* General flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* most_prevalent_univcoord */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Kmer_anypair */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Kmer_segment, Kmer_prevalent, Kmer_search_complete */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Kmer_anchored */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif



/* Merging faster than count table */
#define USE_MERGE 1

#define MAX_NEIGHBORS 3		/* Cannot be 0 */
#define SUBOPT 3

#define LONG_END 6
#define ALLOWED_END_MISMATCHES 2 /* For long ends */
#define ALLOWED_TRANSCRIPTOME_TRIM 3

#define MAX_ACCUMULATED 100	/* For Kmer_anypair */


static int index1part;
static int index1interval;
static Indexdb_T indexdb;

static EF64_T chromosome_ef64;
static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Univcoord_T genomelength;

static Chrpos_T negative_gap_distance; /* typically max_insertionlen */
static Chrpos_T positive_gap_distance;

static bool splicingp;


/* Indexdb: oligo + diagterm -> streams of diagonals */
/* Merge: an array of diagonals with duplicates */
/* Path: genomic endpoints + gaps + trnums */

/* All calls to Substring_new are for transcriptome.  May need to make call to Univ_IIT_update_chrnum */

/* Genome */
/* Ultrafast: check ends only */
/* find_local_sets: merged->diagonals -> middle_diagonal, left_diagonals, right_diagonals */
/* Algorithm 1a: find_best_path_genome -> complete_path (list of Univdiag_T) */
/* Algorithm 1b: Kmer_search_genome (solve_via_segments_genome): complete_path -> hits */
/* Algorithm 2: Stage3end_complete_path_run_gmap: complete_path -> hits */


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


/* kmer_querystart and kmer_queryend may not go to the ends because of invalid or repetitive oligos */
/* However, we are still using 0 and querylength for auxinfo */
/* _univdiagonals_gplus and _univdiagonals_gminus are aligned */
int
Kmer_exact1 (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	     Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	     Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	     Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
	     Stage1_T stage1, int kmer_querystart, int kmer_queryend, int querylength,
	     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {

  Univcoord_T **mod_univdiagonals_gplus, **mod_univdiagonals_gminus;
  int *mod_nunivdiagonals_gplus, *mod_nunivdiagonals_gminus;
  /* int query_lastpos = querylength - index1part; */
  int mod5, mod3, adj;
  int querypos5, querypos3;

  int i;


  debug(printf("Entered Kmer_exact1 with kmer_querystart %d and kmer_queryend %d\n",
	       kmer_querystart,kmer_queryend));
  debug(Stage1_dump(stage1,querylength));

  if (kmer_queryend < kmer_querystart) {
    *nunivdiagonals_gplus = *nunivdiagonals_gminus = 0;
    *_univdiagonals_gplus = *_univdiagonals_gminus = (Univcoord_T *) NULL;
    return 0;

  } else {
    mod_univdiagonals_gplus = (Univcoord_T **) MALLOC(index1interval*sizeof(Univcoord_T *));
    mod_nunivdiagonals_gplus = (int *) MALLOC((index1interval + 1)*sizeof(int)); /* Add 1 needed for Merge_diagonals */
    mod_univdiagonals_gminus = (Univcoord_T **) MALLOC(index1interval*sizeof(Univcoord_T *));
    mod_nunivdiagonals_gminus = (int *) MALLOC((index1interval + 1)*sizeof(int)); /* Add 1 needed for Merge_diagonals */
  }

  for (mod5 = 0; mod5 < index1interval; mod5++) {
    querypos5 = kmer_querystart + mod5;
    adj = (kmer_queryend - kmer_querystart) % index1interval;
    mod3 = (index1interval + adj - mod5) % index1interval;
    querypos3 = kmer_queryend - mod3;
    debug(printf("mod5 %d => querypos5 %d => mod3 %d => querypos3 %d\n",mod5,querypos5,mod3,querypos3));
    if (querypos3 < querypos5) {
      /* Skip */
      mod_nunivdiagonals_gplus[mod5] = 0;
      mod_univdiagonals_gplus[mod5] = (Univcoord_T *) NULL;
      mod_nunivdiagonals_gminus[mod3] = 0;
      mod_univdiagonals_gminus[mod3] = (Univcoord_T *) NULL;
      
    } else {
      assert((querypos5 % index1interval) == (querypos3 % index1interval));

      /* gplus */
#ifdef LARGE_GENOMES
      mod_univdiagonals_gplus[mod5] =
	Intersect_large(&(mod_nunivdiagonals_gplus[mod5]),
			stage1->plus_positions_high[querypos5],stage1->plus_positions[querypos5],
			stage1->plus_npositions[querypos5],/*diagterm1*/querylength - querypos5,
			stage1->plus_positions_high[querypos3],stage1->plus_positions[querypos3],
			stage1->plus_npositions[querypos3],/*diagterm2*/querylength - querypos3);
#else
      mod_univdiagonals_gplus[mod5] =
	Intersect_small(&(mod_nunivdiagonals_gplus[mod5]),
			stage1->plus_positions[querypos5],stage1->plus_npositions[querypos5],
			/*diagterm1*/querylength - querypos5,
			stage1->plus_positions[querypos3],stage1->plus_npositions[querypos3],
			/*diagterm2*/querylength - querypos3,/*alignp*/false);
#endif
      
      debug(printf("gplus mod5 %d, mod3 %d: diagterm5 %d, diagterm3 %d, %d univdiagonals\n",
		   mod5,mod3,querylength - querypos5,querylength - querypos3,
		   mod_nunivdiagonals_gplus[mod5]));

      /* gminus */
#ifdef LARGE_GENOMES
      mod_univdiagonals_gminus[mod3] =
	Intersect_large(&(mod_nunivdiagonals_gminus[mod3]),
			stage1->minus_positions_high[querypos5],stage1->minus_positions[querypos5],
			stage1->minus_npositions[querypos5],/*diagterm1*/querypos5 + index1part,
			stage1->minus_positions_high[querypos3],stage1->minus_positions[querypos3],
			stage1->minus_npositions[querypos3],/*diagterm2*/querypos3 + index1part);
#else
      mod_univdiagonals_gminus[mod3] =
	Intersect_small(&(mod_nunivdiagonals_gminus[mod3]),
			stage1->minus_positions[querypos5],stage1->minus_npositions[querypos5],
			/*diagterm1*/querypos5 + index1part,
			stage1->minus_positions[querypos3],stage1->minus_npositions[querypos3],
			/*diagterm2*/querypos3 + index1part,/*alignp*/false);
#endif
      
      debug(printf("gminus mod5 %d, mod3 %d: diagterm5 %d, diagterm3 %d, %d univdiagonals\n",
		   mod5,mod3,querypos5 + index1part,querypos3 + index1part,
		   mod_nunivdiagonals_gminus[mod3]));
    }
  }
    
#ifdef LARGE_GENOMES  
  *_univdiagonals_gplus = Merge_diagonals_uint8(&(*nunivdiagonals_gplus),
						mod_univdiagonals_gplus,mod_nunivdiagonals_gplus,
						/*nstreams*/index1interval,stage1->mergeinfo);
  *_univdiagonals_gminus = Merge_diagonals_uint8(&(*nunivdiagonals_gminus),
						 mod_univdiagonals_gminus,mod_nunivdiagonals_gminus,
						 /*nstreams*/index1interval,stage1->mergeinfo);
#else
  *_univdiagonals_gplus = Merge_diagonals_uint4(&(*nunivdiagonals_gplus),
						mod_univdiagonals_gplus,mod_nunivdiagonals_gplus,
						/*nstreams*/index1interval,stage1->mergeinfo);
  *_univdiagonals_gminus = Merge_diagonals_uint4(&(*nunivdiagonals_gminus),
						 mod_univdiagonals_gminus,mod_nunivdiagonals_gminus,
						 /*nstreams*/index1interval,stage1->mergeinfo);
#endif
  
  /* Can have duplicates among the various mods */
  *nunivdiagonals_gplus = make_unique(*_univdiagonals_gplus,*nunivdiagonals_gplus);
  *nunivdiagonals_gminus = make_unique(*_univdiagonals_gminus,*nunivdiagonals_gminus);

  /* auxinfo */
  if (*nunivdiagonals_gplus == 0) {
    *auxinfo_gplus = (Auxinfo_T *) NULL;
    *pathstores_gplus = (Pathstore_T *) NULL;
  } else {
    *auxinfo_gplus = MALLOC((*nunivdiagonals_gplus)*sizeof(Auxinfo_T));
    *pathstores_gplus = MALLOC((*nunivdiagonals_gplus)*sizeof(Pathstore_T));
    for (i = 0; i < *nunivdiagonals_gplus; i++) {
      (*auxinfo_gplus)[i] = Auxinfo_new(KMER_EXACT1,/*qstart*/kmer_querystart,/*qend*/kmer_queryend,
					auxinfopool);
      (*pathstores_gplus)[i] = Pathstore_new(pathstorepool);
    }
  }

  if (*nunivdiagonals_gminus == 0) {
    *auxinfo_gminus = (Auxinfo_T *) NULL;
    *pathstores_gminus = (Pathstore_T *) NULL;
  } else {
    *auxinfo_gminus = MALLOC((*nunivdiagonals_gminus)*sizeof(Auxinfo_T));
    *pathstores_gminus = MALLOC((*nunivdiagonals_gminus)*sizeof(Pathstore_T));
    for (i = 0; i < *nunivdiagonals_gminus; i++) {
      (*auxinfo_gminus)[i] = Auxinfo_new(KMER_EXACT1,/*qstart*/kmer_querystart,/*qend*/kmer_queryend,
					 auxinfopool);
      (*pathstores_gminus)[i] = Pathstore_new(pathstorepool);
    }
  }

  for (mod3 = 0; mod3 < index1interval; mod3++) {
    FREE(mod_univdiagonals_gminus[mod3]);
  }
  FREE(mod_univdiagonals_gminus);
  FREE(mod_nunivdiagonals_gminus);

  for (mod5 = 0; mod5 < index1interval; mod5++) {
    FREE(mod_univdiagonals_gplus[mod5]);
  }
  FREE(mod_univdiagonals_gplus);
  FREE(mod_nunivdiagonals_gplus);

#ifdef DEBUG
  printf("Kmer_exact1 returning %d plus and %d minus hits\n",
	 (*nunivdiagonals_gplus),(*nunivdiagonals_gminus));

  for (int i = 0; i < *nunivdiagonals_gplus; i++) {
    printf("plus %u\n",(*_univdiagonals_gplus)[i]);
  }
  for (int i = 0; i < *nunivdiagonals_gminus; i++) {
    printf("minus %u\n",(*_univdiagonals_gminus)[i]);
  }
#endif

  return (*nunivdiagonals_gplus) + (*nunivdiagonals_gminus);
}


#if 0
static int
keep_duplicates (Univcoord_T *diagonals, int ndiagonals) {
  int k = 0, i, j;

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    if (j - i > 1) {
      diagonals[k++] = diagonals[i];
    }
    
    i = j;
  }

  return k;
}
#endif


/* We need to limit the number of univdiags, or else we could get
   hundreds and the solution by Path_solve_diagonals will be very
   time-consuming */
Auxinfo_T
Kmer_compute_auxinfo_univdiags (Univcoord_T main_univdiagonal, int qstart, int qend, int i,
				Univcoord_T *exhaustive, int *qstarts, int *qends, int nexhaustive,
				Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
				Method_T method) {
  List_T left_univdiags = NULL, right_univdiags = NULL;
  int nmismatches;
  int jlow, jhigh, j, k;

  /* left_univdiags */
  jlow = i - 1;
  while (jlow >= 0 && exhaustive[jlow] + positive_gap_distance > main_univdiagonal) {
    jlow--;
  }
  jlow++;
      
  jhigh = i + 1;
  while (jhigh < nexhaustive && exhaustive[jhigh] < main_univdiagonal + negative_gap_distance) {
    jhigh++;
  }

  /* Should traverse in order from medial to distal, meaning downstream */
  k = 0;
  j = jhigh - 1;
  while (j >= jlow && k < MAX_UNIVDIAGS) {
    if (qstarts[j] < qstart && qends[j] < qend) {
      nmismatches = (qends[j] - qstarts[j] <= 2*index1part) ? 0 : -1;
      left_univdiags = Univdiagpool_push(left_univdiags,univdiagpool,
					 qstarts[j],qends[j],nmismatches,exhaustive[j]
					 univdiagpool_trace(__FILE__,__LINE__));
      k++;
    }
    j--;
  }
      
  /* right_univdiags */
  right_univdiags = (List_T) NULL;
  jlow = i - 1;
  while (jlow >= 0 && exhaustive[jlow] + negative_gap_distance > main_univdiagonal) {
    jlow--;
  }
  jlow++;
  
  jhigh = i + 1;
  while (jhigh < nexhaustive && exhaustive[jhigh] < main_univdiagonal + positive_gap_distance) {
    jhigh++;
  }

  /* Should traverse in order from medial to distal, meaning upstream */
  k = 0;
  j = jlow;
  while (j < jhigh && k < MAX_UNIVDIAGS) {
    if (qstarts[j] > qstart && qends[j] > qend) {
      nmismatches = (qends[j] - qstarts[j] <= 2*index1part) ? 0 : -1;
      right_univdiags = Univdiagpool_push(right_univdiags,univdiagpool,
					  qstarts[j],qends[j],nmismatches,exhaustive[j]
					  univdiagpool_trace(__FILE__,__LINE__));
      k++;
    }
    j++;
  }

  nmismatches = (qend - qstart <= 2*index1part) ? 0 : -1;
  return Auxinfo_new_univdiags(method,qstart,qend,nmismatches,right_univdiags,left_univdiags,
			       auxinfopool);
}


/* _univdiagonals needs to be aligned for Kmer_segment to match Kmer_exact1 */
static int
compute_univdiagonals_auxinfo (Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores,
			       Univcoord_T *exhaustive, int *qstarts, int *qends,
			       int *counts, int nexhaustive, Univdiagpool_T univdiagpool,
			       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
			       Method_T method) {
  int k, i;
  
  int max_width = 0, width, qstart, qend;
  int max_count = 0;
  Univcoord_T main_univdiagonal;


  /* Find max_width and max_count */
  for (i = 0; i < nexhaustive; i++) {
    if ((width = qends[i] - qstarts[i]) > max_width) {
      max_width = width;
    }
    if (counts[i] > max_count) {
      max_count = counts[i];
    }
  }
  max_width -= (index1interval - 1);
  debug9(printf("max_width is %d\n",max_width));
  debug9(printf("max_count is %d\n",max_count));

  /* Compute univdiagonals and auxinfo */
  MALLOC_ALIGN(*_univdiagonals,nexhaustive*sizeof(Univcoord_T));
  *auxinfo = (Auxinfo_T *) MALLOC(nexhaustive*sizeof(Auxinfo_T));
  *pathstores = (Pathstore_T *) MALLOC(nexhaustive*sizeof(Pathstore_T));

  k = 0;
  for (i = 0; i < nexhaustive; i++) {
    if (qends[i] - qstarts[i] >= max_width) {
      main_univdiagonal = (*_univdiagonals)[k] = exhaustive[i];
      debug9(printf("%u wide\n",main_univdiagonal));
      qstart = qstarts[i];
      qend = qends[i];
      (*auxinfo)[k] = Kmer_compute_auxinfo_univdiags(main_univdiagonal,qstart,qend,i,
						     exhaustive,qstarts,qends,nexhaustive,
						     univdiagpool,auxinfopool,method);
      (*pathstores)[k] = Pathstore_new(pathstorepool);
      k++;

    } else if (counts[i] >= max_count) {
      main_univdiagonal = (*_univdiagonals)[k] = exhaustive[i];
      debug9(printf("%u prevalent\n",main_univdiagonal));
      qstart = qstarts[i];
      qend = qends[i];
      (*auxinfo)[k] = Kmer_compute_auxinfo_univdiags(main_univdiagonal,qstart,qend,i,
						     exhaustive,qstarts,qends,nexhaustive,
						     univdiagpool,auxinfopool,method);
      (*pathstores)[k] = Pathstore_new(pathstorepool);
      k++;
    }
  }

  return k;
}


static int
get_records_gplus (Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores,
		   Stage1_T stage1, int querylength, int total_npositions_plus,
		   Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {

  int nunivdiagonals, nexhaustive, ndiagonals, nindices;
  int query_lastpos, querypos;

  Univcoord_T *_diagonals;	/* aligned */
  int *qstarts, *qends;
  int *counts;
  int *indices, index2, k;

#ifndef BRANCHED_CODE
  int **all_indices, *space;
  int *all_nindices;
  int qend;
#endif


  query_lastpos = querylength - index1part;

  /* Compute exhaustive univdiagonals, aligned */
#ifdef LARGE_GENOMES
  _diagonals = Merge_diagonals_large(&ndiagonals,stage1->plus_positions_high,
				     stage1->plus_positions,stage1->plus_npositions,
				     stage1->plus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#else
  _diagonals = Merge_diagonals(&ndiagonals,stage1->plus_positions,stage1->plus_npositions,
			       stage1->plus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#endif

  nexhaustive = make_unique(_diagonals,ndiagonals); /* overwrites _diagonals, aligned */


  /* Compute qstarts and qends */
  qstarts = (int *) MALLOC(nexhaustive*sizeof(int));
#ifdef BRANCHED_CODE
  qends = (int *) CALLOC(nexhaustive,sizeof(int));
#else
  qends = (int *) MALLOC(nexhaustive*sizeof(int));
#endif
  counts = (int *) CALLOC(nexhaustive,sizeof(int));
  
#ifdef BRANCHED_CODE
  /* Code with branches inside a loop */
  /* For plus, iterate from 0 to query_lastpos, so we update qends in increasing order */
  indices = MALLOC(max_npositions * sizeof(int));
  
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#ifdef LARGE_GENOMES
    nindices = Intersect_indices2_large(indices,stage1->plus_positions_high[querypos],
					stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
					/*diagterm1, plus*/stage1->plus_diagterms[querypos],
					_diagonals,nexhaustive,/*diagterm2*/0);
#else
    nindices = Intersect_indices2_small(indices,stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
					/*diagterm1, plus*/stage1->plus_diagterms[querypos],
					_diagonals,nexhaustive,/*diagterm2*/0);
#endif
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      if (qends[index2] == 0) {
	qstarts[index2] = querypos; /* plus */
      }
      qends[index2] = querypos + index1part; /* plus */
      counts[index2] += 1;
    }
  }

  FREE(indices);

#else
  /* Code without branches.  Requires two iterations, one to update qends and one to update qstarts */
  /* To avoid duplicate calculations, we save the indices */

  /* Allocate a 2-dimensional array */
  all_indices = (int **) MALLOC((query_lastpos + 2) * sizeof(int *)); /* Allocate 1 extra to write the end pointer in the loop */
  all_indices[0] = space = (int *) MALLOC(total_npositions_plus * sizeof(int));
  all_nindices = (int *) MALLOC((query_lastpos + 1) * sizeof(int));


  /* Loop 1: Update qends */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#ifdef LARGE_GENOMES
    nindices = all_nindices[querypos] =
      Intersect_indices2_large(all_indices[querypos],stage1->plus_positions_high[querypos],
			      stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
			      /*diagterm1, plus*/stage1->plus_diagterms[querypos],
			      _diagonals,nexhaustive,/*diagterm2*/0);
#else
    nindices = all_nindices[querypos] =
      Intersect_indices2_small(all_indices[querypos],
			      stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
			      /*diagterm1, plus*/stage1->plus_diagterms[querypos],
			      _diagonals,nexhaustive,/*diagterm2*/0);
#endif
    indices = all_indices[querypos];
    qend = querypos + index1part; /* plus */
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      qends[index2] = qend;
      counts[index2] += 1;
    }
    all_indices[querypos+1] = &(indices[nindices]); /* Advance pointer */
  }
  
  /* Loop 2: Update qstarts */
  for (querypos = query_lastpos; querypos >= 0; querypos--) {
    nindices = all_nindices[querypos];
    indices = all_indices[querypos];
    /* qstart = querypos; -- plus */
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      qstarts[index2] = querypos; /* plus */
      /* counts[index2] += 1; -- Already counted in loop 1 */
    }
  }

  FREE(all_nindices);
  FREE(all_indices);
  FREE(space);
#endif

  nunivdiagonals =
    compute_univdiagonals_auxinfo(&(*_univdiagonals),&(*auxinfo),&(*pathstores),
				  /*exhaustive*/_diagonals,qstarts,qends,counts,nexhaustive,
				  univdiagpool,auxinfopool,pathstorepool,/*method*/SEGMENT1);

  stage1->exhaustive_gplus = _diagonals;
  stage1->exhaustive_qstart_gplus = qstarts;
  stage1->exhaustive_qend_gplus = qends;
  stage1->exhaustive_counts_gplus = counts;
  stage1->nexhaustive_gplus = nexhaustive;

  return nunivdiagonals;
}


static int
get_records_gminus (Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo, Pathstore_T **pathstores,
		    Stage1_T stage1, int querylength, int total_npositions_minus,
		    Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {

  int nunivdiagonals, nexhaustive, ndiagonals, nindices;
  int query_lastpos, querypos;

  Univcoord_T *_diagonals;	/* aligned */
  int *qstarts, *qends;
  int *counts;
  int *indices, index2, k;

#ifndef BRANCHED_CODE
  int **all_indices, *space;
  int *all_nindices;
  int qstart, qend;
#endif


  query_lastpos = querylength - index1part;

  /* Compute exhaustive univdiagonals, aligned */
#ifdef LARGE_GENOMES
  _diagonals = Merge_diagonals_large(&ndiagonals,stage1->minus_positions_high,
				    stage1->minus_positions,stage1->minus_npositions,
				    stage1->minus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#else
  _diagonals = Merge_diagonals(&ndiagonals,stage1->minus_positions,stage1->minus_npositions,
			      stage1->minus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#endif
  nexhaustive = make_unique(_diagonals,ndiagonals); /* overwrites _diagonals, aligned */


  /* Compute qstarts and qends */
  qstarts = (int *) MALLOC(nexhaustive*sizeof(int));
#ifdef BRANCHED_CODE
  qends = (int *) CALLOC(nexhaustive,sizeof(int));
#else
  qends = (int *) MALLOC(nexhaustive * sizeof(int));
#endif
  counts = (int *) CALLOC(nexhaustive,sizeof(int));

#ifdef BRANCHED_CODE
  /* Code with branches inside a loop */
  /* For minus, iterate from query_lastpos to 0, so we update qends in increasing order */
  indices = MALLOC(max_npositions * sizeof(int));
  
  for (querypos = query_lastpos; querypos >= 0; querypos--) {
#ifdef LARGE_GENOMES
    nindices = Intersect_indices2_large(indices,stage1->minus_positions_high[querypos],
					stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
					/*diagterm1, minus*/stage1->minus_diagterms[querypos],
					_diagonals,nexhaustive,/*diagterm2*/0);
#else
    nindices = Intersect_indices2_small(indices,stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
					/*diagterm1*/stage1->minus_diagterms[querypos],
					_diagonals,nexhaustive,/*diagterm2*/0);
#endif
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      if (qends[index2] == 0) {
	qstarts[index2] = query_lastpos - querypos; /* minus */
      }
      qends[index2] = querylength - querypos; /* minus */
      counts[index2] += 1;
    }
  }
  FREE(indices);

#else
  /* Code without branches.  Requires two iterations, one to update qends and one to update qstarts */
  /* To avoid duplicate calculations, we save the indices */

  /* Allocate a 2-dimensional array */
  all_indices = (int **) MALLOC((query_lastpos + 2) * sizeof(int *)); /* Allocate 1 extra to write the end pointer in the loop */
  all_indices[0] = space = (int *) MALLOC(total_npositions_minus * sizeof(int));
  all_nindices = (int *) MALLOC((query_lastpos + 1) * sizeof(int));


  /* Loop 1: Update qstarts */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
#ifdef LARGE_GENOMES
    nindices = all_nindices[querypos] =
      Intersect_indices2_large(all_indices[querypos],stage1->minus_positions_high[querypos],
			       stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
			       /*diagterm1, minus*/stage1->minus_diagterms[querypos],
			       _diagonals,nexhaustive,/*diagterm2*/0);
#else
    nindices = all_nindices[querypos] =
      Intersect_indices2_small(all_indices[querypos],
			       stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
			       /*diagterm1*/stage1->minus_diagterms[querypos],
			       _diagonals,nexhaustive,/*diagterm2*/0);
#endif
    indices = all_indices[querypos];
    qstart = query_lastpos - querypos; /* minus */
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      qstarts[index2] = qstart;
      counts[index2] += 1;
    }
    all_indices[querypos+1] = &(indices[nindices]); /* Advance pointer */
  }

  /* Loop 2: Update qends */
  for (querypos = query_lastpos; querypos >= 0; querypos--) {
    nindices = all_nindices[querypos];
    indices = all_indices[querypos];
    qend = querylength - querypos; /* minus */
    for (k = 0; k < nindices; k++) {
      index2 = indices[k];	/* Index of exhaustive */
      qends[index2] = qend;
      /* counts[index2] += 1; -- Already counted in loop 1 */
    }
  }

  FREE(all_nindices);
  FREE(all_indices);
  FREE(space);
#endif

  nunivdiagonals =
    compute_univdiagonals_auxinfo(&(*_univdiagonals),&(*auxinfo),&(*pathstores),
				  /*exhaustive*/_diagonals,qstarts,qends,counts,nexhaustive,
				  univdiagpool,auxinfopool,pathstorepool,/*method*/SEGMENT1);

  stage1->exhaustive_gminus = _diagonals;
  stage1->exhaustive_qstart_gminus = qstarts;
  stage1->exhaustive_qend_gminus = qends;
  stage1->exhaustive_counts_gminus = counts;
  stage1->nexhaustive_gminus = nexhaustive;

  return nunivdiagonals;
}


/* Needed before we call get_records_gplus and get_records_gminus, or
   else paired_search_exhaustive will yield too many intersections */
static void
clear_repetitive_oligos (int *total_npositions_plus, int *total_npositions_minus,
			 Stage1_T stage1, int querylength, EF64_T repetitive_ef64) {
  int querypos;
  int query_lastpos = querylength - index1part;

  *total_npositions_plus = *total_npositions_minus = 0;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (stage1->validp[querypos] == false) {
      /* Skip */
    } else {
      if (EF64_presentp(stage1->forward_oligos[querypos],repetitive_ef64) == true) {
	/* Set npositions to be 0 */
	stage1->plus_npositions[querypos] = 0;
      } else if (stage1->plus_npositions[querypos] > 100000) {
	stage1->plus_npositions[querypos] = 0;
      } else {
	*total_npositions_plus += stage1->plus_npositions[querypos];
      }
      if (EF64_presentp(stage1->revcomp_oligos[querypos],repetitive_ef64) == true) {
	/* Set npositions to be 0 */
	stage1->minus_npositions[querypos] = 0;
      } else if (stage1->minus_npositions[querypos] > 100000) {
	stage1->minus_npositions[querypos] = 0;
      } else {
	*total_npositions_minus += stage1->minus_npositions[querypos];
      }
    }
  }

  return;
}


/* Rewrite of Segment_search */
/* _univdiagonals_gplus and _univdiagonals_gminus must be aligned to match Kmer_exact1 */

int
Kmer_segment (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	      Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	      Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	      Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,
	      Stage1_T stage1, int querylength, EF64_T repetitive_ef64,
	      Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {

  int total_npositions_plus, total_npositions_minus;

  debug9(printf("Entered Kmer_segment\n"));
  /* debug9(Stage1_dump(stage1,querylength)); */
  assert(stage1->all_positions_gen_filledp == true);

  clear_repetitive_oligos(&total_npositions_plus,&total_npositions_minus,
			  stage1,querylength,repetitive_ef64);

  if (total_npositions_plus == 0) {
    *_univdiagonals_gplus = (Univcoord_T *) NULL;
    *auxinfo_gplus = (Auxinfo_T *) NULL;
    *nunivdiagonals_gplus = 0;
  } else {
    *nunivdiagonals_gplus =
      get_records_gplus(&(*_univdiagonals_gplus),&(*auxinfo_gplus),&(*pathstores_gplus),
			stage1,querylength,total_npositions_plus,univdiagpool,auxinfopool,
			pathstorepool);
  }

  if (total_npositions_minus == 0) {
    *_univdiagonals_gminus = (Univcoord_T *) NULL;
    *auxinfo_gminus = (Auxinfo_T *) NULL;
    *nunivdiagonals_gminus = 0;
  } else {
    *nunivdiagonals_gminus =
      get_records_gminus(&(*_univdiagonals_gminus),&(*auxinfo_gminus),&(*pathstores_gminus),
			 stage1,querylength,total_npositions_minus,univdiagpool,auxinfopool,
			 pathstorepool);
  }

#ifdef DEBUG9
  printf("Kmer_segment returning %d plus and %d minus univdiagonals\n",
	 (*nunivdiagonals_gplus),(*nunivdiagonals_gminus));

  for (int i = 0; i < *nunivdiagonals_gplus; i++) {
    printf("plus %u\n",(*_univdiagonals_gplus)[i]);
  }
  for (int i = 0; i < *nunivdiagonals_gminus; i++) {
    printf("minus %u\n",(*_univdiagonals_gminus)[i]);
  }
#endif

  return (*nunivdiagonals_gplus) + (*nunivdiagonals_gminus);
}



void
Kmer_search_setup (int index1part_in, int index1interval_in,
		   Indexdb_T indexdb_in, EF64_T chromosome_ef64_in,
		   Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		   Univcoord_T genomelength_in, 
		   int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist,
		   bool splicingp_in) {

  index1part = index1part_in;
  index1interval = index1interval_in;
  indexdb = indexdb_in;

  chromosome_ef64 = chromosome_ef64_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  genomelength = genomelength_in;

  negative_gap_distance = (Chrpos_T) max_insertionlen;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  
  splicingp = splicingp_in;

  return;
}


