#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "kmer-search.h"

#include <stdio.h>
#include <limits.h>
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
#include "singleelt.h"
#include "minpure.h"

#include "genomebits_count.h"
#include "genomebits_kmer.h"
#include "splice.h"
#include "indel.h"
#include "intron.h"
#include "maxent_hr.h"
#include "sedgesort.h"
#include "path-solve.h"
#include "accel-profile-stats.h"


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

#define SEGMENT_OVERLAP 3


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
static Minpure_Table_T genome_minpure = NULL;
#ifdef GSNAP_DIAGONAL_ROW_STATS
static bool diagstats_header_printed_p = false;
#endif
static int index1interval;

static EF64_T chromosome_ef64;
static Univ_IIT_T chromosome_iit;
static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;
static Univcoord_T genomelength;

static Chrpos_T negative_gap_distance; /* typically max_insertionlen */
static Chrpos_T positive_gap_distance;

static bool splicingp;
static bool minpure_diagonal_debug_p = false;
static bool minpure_diagonal_debug_header_printed_p = false;

static unsigned int
minpure_specificity_from_frequency (UINT4 minfreq) {
  unsigned int denom, specificity;

  denom = minfreq == 0U ? 1U : (unsigned int) minfreq;
  specificity = STAGE1_MINPURE_SPECIFICITY_SCALE/denom;
  if (specificity == 0U) {
    specificity = 1U;
  }

  return specificity;
}

static void
minpure_specificity_add (unsigned int *specificity_sums, int *specificity_ns,
			 unsigned int *specificity_mins, unsigned int *specificity_maxs,
			 unsigned int *frequency_mins, unsigned int *frequency_maxs,
			 unsigned long long *frequency_sums,
			 int i, UINT4 minfreq) {
  unsigned int frequency, specificity;
  int n;

  frequency = (unsigned int) minfreq;
  specificity = minpure_specificity_from_frequency(minfreq);

  n = specificity_ns[i];

  if (UINT_MAX - specificity_sums[i] < specificity) {
    specificity_sums[i] = UINT_MAX;
  } else {
    specificity_sums[i] += specificity;
  }

  if (specificity_mins != NULL && (n == 0 || specificity < specificity_mins[i])) {
    specificity_mins[i] = specificity;
  }
  if (specificity_maxs != NULL && (n == 0 || specificity > specificity_maxs[i])) {
    specificity_maxs[i] = specificity;
  }
  if (frequency_mins != NULL && (n == 0 || frequency < frequency_mins[i])) {
    frequency_mins[i] = frequency;
  }
  if (frequency_maxs != NULL && (n == 0 || frequency > frequency_maxs[i])) {
    frequency_maxs[i] = frequency;
  }
  if (frequency_sums != NULL) {
    frequency_sums[i] += (unsigned long long) frequency;
  }

  specificity_ns[i]++;

  return;
}

static void
minpure_specificity_finalize (unsigned int *specificity_sums, int *specificity_ns, int n) {
  int i;

  if (specificity_sums == NULL) {
    return;
  }
  for (i = 0; i < n; i++) {
    if (specificity_ns[i] > 0) {
      specificity_sums[i] /= (unsigned int) specificity_ns[i];
    }
  }

  return;
}

static Segment_diagonal_algorithm_T segment_diagonal_algorithm = SEGMENT_DIAGONAL_ALGORITHM_WIDTH_OR_COUNT;
static int adaptive_hybrid_floor_per_kmer = 2;
static int adaptive_hybrid_budget_per_kmer = 12;
static int adaptive_hybrid_hard_max = 0;
static int adaptive_hybrid_weak_width_per_kmer = 6;
static int adaptive_hybrid_weak_support_divisor = 8;
static double adaptive_hybrid_strong_support_fraction = 0.10;
static double stats_weak_support_fraction = 0.05;
static bool segment_repetitive_kmer_filter_p = true;
static int segment_max_kmer_positions = 100000;


static const char *
segment_diagonal_algorithm_name (void) {
  switch (segment_diagonal_algorithm) {
  case SEGMENT_DIAGONAL_ALGORITHM_WIDTH_OR_COUNT:
    return "width-or-count";
  case SEGMENT_DIAGONAL_ALGORITHM_ADAPTIVE_HYBRID:
    return "adaptive-hybrid";
  default:
    abort();
  }
}


#ifdef GSNAP_ACCEL_PROFILE_STATS
static int
int_ascending_cmp (const void *a, const void *b) {
  int x = * (const int *) a;
  int y = * (const int *) b;

  return x < y ? -1 : x > y ? 1 : 0;
}

static int
double_ascending_cmp (const void *a, const void *b) {
  double x = * (const double *) a;
  double y = * (const double *) b;

  return x < y ? -1 : x > y ? 1 : 0;
}

static int
uint8_ascending_cmp (const void *a, const void *b) {
  unsigned long long x = * (const unsigned long long *) a;
  unsigned long long y = * (const unsigned long long *) b;

  return x < y ? -1 : x > y ? 1 : 0;
}

static int
quantile_index (int n, int percentile) {
  return (int) (((long long) (n - 1) * (long long) percentile)/100LL);
}

static void
segment_diagonal_report_distribution (Univcoord_T *exhaustive, Singleelt_T *eltarray,
				      int *counts, int nexhaustive,
				      int min_width, int max_count,
				      const char *label, const char *strand) {
  struct Accel_profile_segment_distribution_T d;
  int *widths, *support_counts;
  double *densities;
  unsigned long long *gaps;
  Univcoord_T previous = 0U;
  Singleelt_T elt;
  int i, k = 0, ngaps = 0, width;
  bool width_selected_p, count_selected_p, have_previous_p = false;

  memset(&d,0,sizeof(d));
  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    if (width >= min_width || counts[i] >= max_count) {
      d.candidates++;
    }
  }
  if (d.candidates == 0) {
    Accel_profile_segment_distribution(segment_diagonal_algorithm_name(),label,strand,&d);
    return;
  }

  widths = (int *) MALLOC(d.candidates * sizeof(int));
  support_counts = (int *) MALLOC(d.candidates * sizeof(int));
  densities = (double *) MALLOC(d.candidates * sizeof(double));
  gaps = d.candidates > 1 ?
    (unsigned long long *) MALLOC((d.candidates - 1) * sizeof(unsigned long long)) : NULL;

  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    width_selected_p = width >= min_width;
    count_selected_p = counts[i] >= max_count;
    if (width_selected_p == true || count_selected_p == true) {
      widths[k] = width;
      support_counts[k] = counts[i];
      densities[k] = width <= 0 ? 0.0 : (double) counts[i]/(double) width;
      if (width_selected_p == true && count_selected_p == true) {
	d.width_and_count++;
      } else if (width_selected_p == true) {
	d.width_only++;
      } else {
	d.count_only++;
      }
      if (have_previous_p == true) {
	gaps[ngaps++] = exhaustive[i] >= previous ?
	  (unsigned long long) (exhaustive[i] - previous) :
	  (unsigned long long) (previous - exhaustive[i]);
      }
      previous = exhaustive[i];
      have_previous_p = true;
      k++;
    }
  }

  qsort(widths,d.candidates,sizeof(int),int_ascending_cmp);
  qsort(support_counts,d.candidates,sizeof(int),int_ascending_cmp);
  qsort(densities,d.candidates,sizeof(double),double_ascending_cmp);
  d.width_min = widths[0];
  d.width_p50 = widths[quantile_index(d.candidates,50)];
  d.width_p90 = widths[quantile_index(d.candidates,90)];
  d.width_p95 = widths[quantile_index(d.candidates,95)];
  d.width_p99 = widths[quantile_index(d.candidates,99)];
  d.width_max = widths[d.candidates - 1];
  d.count_min = support_counts[0];
  d.count_p50 = support_counts[quantile_index(d.candidates,50)];
  d.count_p90 = support_counts[quantile_index(d.candidates,90)];
  d.count_p95 = support_counts[quantile_index(d.candidates,95)];
  d.count_p99 = support_counts[quantile_index(d.candidates,99)];
  d.count_max = support_counts[d.candidates - 1];
  d.density_min = densities[0];
  d.density_p50 = densities[quantile_index(d.candidates,50)];
  d.density_p90 = densities[quantile_index(d.candidates,90)];
  d.density_p95 = densities[quantile_index(d.candidates,95)];
  d.density_p99 = densities[quantile_index(d.candidates,99)];
  d.density_max = densities[d.candidates - 1];

  if (ngaps > 0) {
    qsort(gaps,ngaps,sizeof(unsigned long long),uint8_ascending_cmp);
    d.gap_p50 = gaps[quantile_index(ngaps,50)];
    d.gap_p90 = gaps[quantile_index(ngaps,90)];
    d.gap_p95 = gaps[quantile_index(ngaps,95)];
    d.gap_p99 = gaps[quantile_index(ngaps,99)];
    d.gap_max = gaps[ngaps - 1];
  }

  Accel_profile_segment_distribution(segment_diagonal_algorithm_name(),label,strand,&d);
  FREE(gaps);
  FREE(densities);
  FREE(support_counts);
  FREE(widths);
  return;
}
#endif


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
make_unique (int **counts, Univcoord_T *diagonals, int ndiagonals) {
  int k = 0, i, j;

  *counts = (int *) MALLOC((ndiagonals + 1) * sizeof(int));

  i = 0;
  while (i < ndiagonals) {
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    diagonals[k] = diagonals[i];
    (*counts)[k++] = j - i;
    
    i = j;
  }

  return k;
}


#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
#define SEGMENT_DIAGONAL_TOP_N 8

static void
segment_diagonal_print_merge_stats (const char *label, const char *strand,
				    int raw_ndiagonals, int nexhaustive,
				    int *counts) {
  int i, count;
  int max_count = 0, count_ge2 = 0, count_ge5 = 0, count_ge10 = 0;
  int count_ge50 = 0, count_ge100 = 0, count_ge1000 = 0;
  unsigned long long total_counts = 0U;

  for (i = 0; i < nexhaustive; i++) {
    count = counts[i];
    total_counts += (unsigned long long) count;
    if (count > max_count) {
      max_count = count;
    }
    if (count >= 2) {
      count_ge2++;
    }
    if (count >= 5) {
      count_ge5++;
    }
    if (count >= 10) {
      count_ge10++;
    }
    if (count >= 50) {
      count_ge50++;
    }
    if (count >= 100) {
      count_ge100++;
    }
    if (count >= 1000) {
      count_ge1000++;
    }
  }

  fprintf(stderr,
	  "[SEGMENT_DIAGONAL_MERGE] label=%s strand=%s raw_diagonals=%d "
	  "unique_diagonals=%d total_counts=%llu avg_count=%.2f max_count=%d "
	  "count_ge2=%d count_ge5=%d count_ge10=%d count_ge50=%d "
	  "count_ge100=%d count_ge1000=%d\n",
	  label,strand,raw_ndiagonals,nexhaustive,total_counts,
	  nexhaustive == 0 ? 0.0 : (double) total_counts/(double) nexhaustive,
	  max_count,count_ge2,count_ge5,count_ge10,count_ge50,count_ge100,count_ge1000);

  return;
}


static void
segment_diagonal_record_top (int *top_indices, int *top_counts, int *top_widths,
			     int elti, int count, int width) {
  int rank, j;

  for (rank = 0; rank < SEGMENT_DIAGONAL_TOP_N; rank++) {
    if (count > top_counts[rank] ||
	(count == top_counts[rank] && width > top_widths[rank])) {
      for (j = SEGMENT_DIAGONAL_TOP_N - 1; j > rank; j--) {
	top_indices[j] = top_indices[j - 1];
	top_counts[j] = top_counts[j - 1];
	top_widths[j] = top_widths[j - 1];
      }
      top_indices[rank] = elti;
      top_counts[rank] = count;
      top_widths[rank] = width;
      return;
    }
  }

  return;
}


static void
segment_diagonal_print_select_stats (const char *label, const char *strand,
				     Univcoord_T *exhaustive, Singleelt_T *eltarray,
				     int *counts, int nexhaustive, int selected,
				     int min_width, int max_width, int max_count) {
  Singleelt_T elt;
  int i, width, rank, selected_by_width = 0, selected_by_count = 0;
  int count_ge2 = 0, count_ge5 = 0, count_ge10 = 0, count_ge50 = 0;
  int top_indices[SEGMENT_DIAGONAL_TOP_N];
  int top_counts[SEGMENT_DIAGONAL_TOP_N];
  int top_widths[SEGMENT_DIAGONAL_TOP_N];
  unsigned long long total_width = 0U;

  for (rank = 0; rank < SEGMENT_DIAGONAL_TOP_N; rank++) {
    top_indices[rank] = -1;
    top_counts[rank] = 0;
    top_widths[rank] = 0;
  }

  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    total_width += (unsigned long long) width;
    if (width >= min_width) {
      selected_by_width++;
    }
    if (counts[i] >= max_count) {
      selected_by_count++;
    }
    if (counts[i] >= 2) {
      count_ge2++;
    }
    if (counts[i] >= 5) {
      count_ge5++;
    }
    if (counts[i] >= 10) {
      count_ge10++;
    }
    if (counts[i] >= 50) {
      count_ge50++;
    }
    segment_diagonal_record_top(top_indices,top_counts,top_widths,i,counts[i],width);
  }

  fprintf(stderr,
	  "[SEGMENT_DIAGONAL_SELECT] label=%s strand=%s exhaustive=%d selected=%d "
	  "min_width=%d max_width=%d avg_width=%.2f max_count=%d "
	  "selected_by_width=%d selected_by_count=%d count_ge2=%d count_ge5=%d "
	  "count_ge10=%d count_ge50=%d\n",
	  label,strand,nexhaustive,selected,min_width,max_width,
	  nexhaustive == 0 ? 0.0 : (double) total_width/(double) nexhaustive,
	  max_count,selected_by_width,selected_by_count,
	  count_ge2,count_ge5,count_ge10,count_ge50);

  for (rank = 0; rank < SEGMENT_DIAGONAL_TOP_N; rank++) {
    i = top_indices[rank];
    if (i >= 0) {
      elt = eltarray[i];
      width = elt->qend - elt->qstart;
      fprintf(stderr,
	      "[SEGMENT_DIAGONAL_TOP] label=%s strand=%s rank=%d elti=%d "
	      "univdiagonal=%u qstart=%d qend=%d width=%d count=%d selected=%d\n",
	      label,strand,rank + 1,i,exhaustive[i],elt->qstart,elt->qend,width,
	      counts[i],(width >= min_width || counts[i] >= max_count) ? 1 : 0);
    }
  }

  return;
}
#endif


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
  int *counts_gplus, *counts_gminus;
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
    
  /* Only three streams, so SIMD merge might be faster than heap merge */
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
  if (*nunivdiagonals_gplus > 0) {
    *nunivdiagonals_gplus = make_unique(&counts_gplus,*_univdiagonals_gplus,*nunivdiagonals_gplus);
    FREE(counts_gplus);
  }
  if (*nunivdiagonals_gminus > 0) {
    *nunivdiagonals_gminus = make_unique(&counts_gminus,*_univdiagonals_gminus,*nunivdiagonals_gminus);
    FREE(counts_gminus);
  }

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


static int
store_exhaustive_gplus (Stage1_T stage1, int querylength, const char *label) {
  int nexhaustive;
  int ndiagonals;
  int query_lastpos;

  query_lastpos = querylength - index1part;

  /* Compute exhaustive univdiagonals, aligned */
#ifdef LARGE_GENOMES
  stage1->exhaustive_gplus = Merge_diagonals_large(&ndiagonals,stage1->plus_positions_high,
						   stage1->plus_positions,stage1->plus_npositions,
						   stage1->plus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#else
  stage1->exhaustive_gplus = Merge_diagonals(&ndiagonals,stage1->plus_positions,stage1->plus_npositions,
					     stage1->plus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#endif

  nexhaustive = make_unique(&stage1->exhaustive_counts_gplus,stage1->exhaustive_gplus,ndiagonals);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  segment_diagonal_print_merge_stats(label,"plus",ndiagonals,nexhaustive,
				     stage1->exhaustive_counts_gplus);
#endif

  return nexhaustive;
}


static int
store_exhaustive_gminus (Stage1_T stage1, int querylength, const char *label) {
  int nexhaustive;
  int ndiagonals;
  int query_lastpos;

  query_lastpos = querylength - index1part;

  /* Compute exhaustive univdiagonals, aligned */
#ifdef LARGE_GENOMES
  stage1->exhaustive_gminus = Merge_diagonals_large(&ndiagonals,stage1->minus_positions_high,
						   stage1->minus_positions,stage1->minus_npositions,
						   stage1->minus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#else
  stage1->exhaustive_gminus = Merge_diagonals(&ndiagonals,stage1->minus_positions,stage1->minus_npositions,
					     stage1->minus_diagterms,/*nstreams*/query_lastpos+1,stage1->mergeinfo);
#endif

  nexhaustive = make_unique(&stage1->exhaustive_counts_gminus,stage1->exhaustive_gminus,ndiagonals);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  segment_diagonal_print_merge_stats(label,"minus",ndiagonals,nexhaustive,
				     stage1->exhaustive_counts_gminus);
#endif

  return nexhaustive;
}


static Singleelt_T *
get_eltarray_gplus (Univcoord_T *exhaustive, int nexhaustive,
		    Stage1_T stage1, int querylength, int total_npositions_plus) {
  Singleelt_T *eltarray, elt;
  /* struct Singleelt_T *eltspace; */
  int nindices;
  int query_lastpos, querypos;
  int *indices, exhaustive_index, k;

  int **all_indices, *space;
  int *all_nindices;
  int qend;
  int *minpure_ns = (int *) NULL;
  UINT4 minfreq;
  bool minfreq_valid_p;

  query_lastpos = querylength - index1part;

  /* Compute qstarts and qends => univdiag_array */
  eltarray = Singleelt_new_array(nexhaustive);
  /* eltspace = eltarray[0]; -- Handle for freeing eltarray.  Need to save before we do qsorts */

  /* Not sure if we need to save these */
  /* qstarts = (int *) MALLOC(nexhaustive * sizeof(int)); */
  /* qends = (int *) MALLOC(nexhaustive * sizeof(int)); */


  /* Code without branches.  Requires two iterations, one to update qends and one to update qstarts */
  /* To avoid duplicate calculations, we save the indices */

  /* Allocate a 2-dimensional array */
  all_indices = (int **) MALLOC((query_lastpos + 2) * sizeof(int *)); /* Allocate 1 extra to write the end pointer in the loop */
  all_indices[0] = space = (int *) MALLOC(total_npositions_plus * sizeof(int));
  all_nindices = (int *) MALLOC((query_lastpos + 1) * sizeof(int));
  if (genome_minpure != NULL) {
    stage1->exhaustive_minpure_specificity_gplus =
      (unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
    if (minpure_diagonal_debug_p == true) {
      stage1->exhaustive_minpure_observations_gplus =
	(int *) CALLOC(nexhaustive,sizeof(int));
      stage1->exhaustive_minpure_specificity_min_gplus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_specificity_max_gplus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_min_gplus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_max_gplus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_sum_gplus =
	(unsigned long long *) CALLOC(nexhaustive,sizeof(unsigned long long));
      minpure_ns = stage1->exhaustive_minpure_observations_gplus;
    } else {
      minpure_ns = (int *) CALLOC(nexhaustive,sizeof(int));
    }
  }

  /* Loop 1: Map nonunique stage1 positions onto unique univdiagonals (exhaustive) */
  /* Update qends; finds maximum qend for each univdiagonal */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    /* Match univdiagonals at querypos with the master list (exhaustive), which is unique */
#ifdef LARGE_GENOMES
    nindices = all_nindices[querypos] =
      Intersect_indices2_large(all_indices[querypos],stage1->plus_positions_high[querypos],
			       stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
			       /*diagterm1, plus*/stage1->plus_diagterms[querypos],
			       exhaustive,nexhaustive,/*diagterm2*/0);
#else
    nindices = all_nindices[querypos] =
      Intersect_indices2_small(all_indices[querypos],
			       stage1->plus_positions[querypos],stage1->plus_npositions[querypos],
			       /*diagterm1, plus*/stage1->plus_diagterms[querypos],
			       exhaustive,nexhaustive,/*diagterm2*/0);
#endif
    indices = all_indices[querypos];
    qend = querypos + index1part; /* plus */
#ifdef GSNAP_DIAGONAL_ROW_STATS
    minfreq_valid_p = genome_minpure != NULL &&
      Minpure_table_lookup(&minfreq,genome_minpure,
			   stage1->forward_oligos[querypos],index1part) == true;
    for (k = 0; k < nindices; k++) {
      exhaustive_index = indices[k];	/* Index of exhaustive */
      elt = eltarray[exhaustive_index];
      elt->univdiagonal = exhaustive[exhaustive_index];
      elt->qend /*= qends[exhaustive_index]*/ = qend;
      if (minfreq_valid_p == true && stage1->exhaustive_minpure_specificity_gplus != NULL) {
	minpure_specificity_add(stage1->exhaustive_minpure_specificity_gplus,
				minpure_ns,
				stage1->exhaustive_minpure_specificity_min_gplus,
				stage1->exhaustive_minpure_specificity_max_gplus,
				stage1->exhaustive_minpure_frequency_min_gplus,
				stage1->exhaustive_minpure_frequency_max_gplus,
				stage1->exhaustive_minpure_frequency_sum_gplus,
				exhaustive_index,minfreq);
      }
#ifdef GSNAP_DIAGONAL_ROW_STATS
      if (elt->diagstats_kmer_hits_n == 0 ||
	  (unsigned int) stage1->plus_npositions[querypos] < elt->diagstats_kmer_hits_min) {
	elt->diagstats_kmer_hits_min = (unsigned int) stage1->plus_npositions[querypos];
      }
      if ((unsigned int) stage1->plus_npositions[querypos] > elt->diagstats_kmer_hits_max) {
	elt->diagstats_kmer_hits_max = (unsigned int) stage1->plus_npositions[querypos];
      }
      elt->diagstats_kmer_hits_sum +=
	(unsigned long long) stage1->plus_npositions[querypos];
      elt->diagstats_kmer_hits_n++;
      if (minfreq_valid_p == true) {
	if (elt->diagstats_minfreq_n == 0 || minfreq < elt->diagstats_minfreq_min) {
	  elt->diagstats_minfreq_min = minfreq;
	}
	if (minfreq > elt->diagstats_minfreq_max) {
	  elt->diagstats_minfreq_max = minfreq;
	}
	elt->diagstats_minfreq_sum += (unsigned long long) minfreq;
	elt->diagstats_minfreq_n++;
      }
#endif
    }
#else
    if (stage1->exhaustive_minpure_specificity_gplus != NULL &&
	Minpure_table_lookup(&minfreq,genome_minpure,
			     stage1->forward_oligos[querypos],index1part) == true) {
      for (k = 0; k < nindices; k++) {
	exhaustive_index = indices[k];	/* Index of exhaustive */
	elt = eltarray[exhaustive_index];
	elt->univdiagonal = exhaustive[exhaustive_index];
	elt->qend /*= qends[exhaustive_index]*/ = qend;
	minpure_specificity_add(stage1->exhaustive_minpure_specificity_gplus,
				minpure_ns,
				stage1->exhaustive_minpure_specificity_min_gplus,
				stage1->exhaustive_minpure_specificity_max_gplus,
				stage1->exhaustive_minpure_frequency_min_gplus,
				stage1->exhaustive_minpure_frequency_max_gplus,
				stage1->exhaustive_minpure_frequency_sum_gplus,
				exhaustive_index,minfreq);
      }
    } else {
      for (k = 0; k < nindices; k++) {
	exhaustive_index = indices[k];	/* Index of exhaustive */
	elt = eltarray[exhaustive_index];
	elt->univdiagonal = exhaustive[exhaustive_index];
	elt->qend /*= qends[exhaustive_index]*/ = qend;
      }
    }
#endif
    all_indices[querypos+1] = &(indices[nindices]); /* Advance pointer */
  }
  
  /* Loop 2: Update qstarts; finds minimum qstart for each univdiagonal */
  for (querypos = query_lastpos; querypos >= 0; querypos--) {
    nindices = all_nindices[querypos];
    indices = all_indices[querypos];
    /* qstart = querypos; -- plus */
    for (k = 0; k < nindices; k++) {
      exhaustive_index = indices[k];	/* Index of exhaustive */
      elt = eltarray[exhaustive_index];
      elt->qstart /*= qstarts[exhaustive_index]*/ = querypos;
    }
  }

  minpure_specificity_finalize(stage1->exhaustive_minpure_specificity_gplus,
			       minpure_ns,nexhaustive);
  if (minpure_ns != NULL && minpure_ns != stage1->exhaustive_minpure_observations_gplus) {
    FREE(minpure_ns);
  }

  FREE(all_nindices);
  FREE(all_indices);
  FREE(space);

#ifdef DEBUG9
  for (k = 0; k < nexhaustive; k++) {
    printf("%d: %d..%d %u\n",k,eltarray[k]->qstart,eltarray[k]->qend,eltarray[k]->univdiagonal);
  }
#endif

  return eltarray;
}


/* Without elts, do we need to keep track of qstarts and qends, or to have the sorted? */
/* exhaustive, qstarts, qends, and counts are sorted by ascending
   univcoord, but univdiag_array is sorted by qstart..qend */
static Singleelt_T *
get_eltarray_gminus (Univcoord_T *exhaustive, int nexhaustive,
		     Stage1_T stage1, int querylength, int total_npositions_minus) {
  Singleelt_T *eltarray, elt;
  /* struct Singleelt_T *eltspace; */
  int nindices;
  int query_lastpos, querypos;
  int *indices, exhaustive_index, k;

  int **all_indices, *space;
  int *all_nindices;
  int qstart, qend;
  int *minpure_ns = (int *) NULL;
  UINT4 minfreq;
  bool minfreq_valid_p;

  query_lastpos = querylength - index1part;

  /* Compute qstarts and qends => univdiag_array */
  eltarray = Singleelt_new_array(nexhaustive);
  /* eltspace = eltarray[0]; -- Handle for freeing eltarray.  Need to save before we do qsorts */

  /* Not sure if we need to save these */
  /* qstarts = (int *) MALLOC(nexhaustive *sizeof(int)); */
  /* qends = (int *) MALLOC(nexhaustive * sizeof(int)); */


  /* Code without branches.  Requires two iterations, one to update qends and one to update qstarts */
  /* To avoid duplicate calculations, we save the indices */

  /* Allocate a 2-dimensional array */
  all_indices = (int **) MALLOC((query_lastpos + 2) * sizeof(int *)); /* Allocate 1 extra to write the end pointer in the loop */
  all_indices[0] = space = (int *) MALLOC(total_npositions_minus * sizeof(int));
  all_nindices = (int *) MALLOC((query_lastpos + 1) * sizeof(int));
  if (genome_minpure != NULL) {
    stage1->exhaustive_minpure_specificity_gminus =
      (unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
    if (minpure_diagonal_debug_p == true) {
      stage1->exhaustive_minpure_observations_gminus =
	(int *) CALLOC(nexhaustive,sizeof(int));
      stage1->exhaustive_minpure_specificity_min_gminus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_specificity_max_gminus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_min_gminus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_max_gminus =
	(unsigned int *) CALLOC(nexhaustive,sizeof(unsigned int));
      stage1->exhaustive_minpure_frequency_sum_gminus =
	(unsigned long long *) CALLOC(nexhaustive,sizeof(unsigned long long));
      minpure_ns = stage1->exhaustive_minpure_observations_gminus;
    } else {
      minpure_ns = (int *) CALLOC(nexhaustive,sizeof(int));
    }
  }

  /* Loop 1: Map nonunique stage1 positions onto unique univdiagonals (exhaustive) */
  /* Update qstarts; finds minimum qstart for each univdiagonal */
  /* Ascending querypos corresponds to descending qstart */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    /* Match univdiagonals at querypos with the master list (exhaustive), which is unique */
#ifdef LARGE_GENOMES
    nindices = all_nindices[querypos] =
      Intersect_indices2_large(all_indices[querypos],stage1->minus_positions_high[querypos],
			       stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
			       /*diagterm1, minus*/stage1->minus_diagterms[querypos],
			       exhaustive,nexhaustive,/*diagterm2*/0);
#else
    nindices = all_nindices[querypos] =
      Intersect_indices2_small(all_indices[querypos],
			       stage1->minus_positions[querypos],stage1->minus_npositions[querypos],
			       /*diagterm1*/stage1->minus_diagterms[querypos],
			       exhaustive,nexhaustive,/*diagterm2*/0);
#endif
    indices = all_indices[querypos];
    qstart = query_lastpos - querypos; /* minus */
#ifdef GSNAP_DIAGONAL_ROW_STATS
    minfreq_valid_p = genome_minpure != NULL &&
      Minpure_table_lookup(&minfreq,genome_minpure,
			   stage1->revcomp_oligos[querypos],index1part) == true;
    for (k = 0; k < nindices; k++) {
      exhaustive_index = indices[k];	/* Index of exhaustive */
      elt = eltarray[exhaustive_index];
      elt->univdiagonal = exhaustive[exhaustive_index];
      elt->qstart /*= qstarts[exhaustive_index]*/ = qstart;
      if (minfreq_valid_p == true && stage1->exhaustive_minpure_specificity_gminus != NULL) {
	minpure_specificity_add(stage1->exhaustive_minpure_specificity_gminus,
				minpure_ns,
				stage1->exhaustive_minpure_specificity_min_gminus,
				stage1->exhaustive_minpure_specificity_max_gminus,
				stage1->exhaustive_minpure_frequency_min_gminus,
				stage1->exhaustive_minpure_frequency_max_gminus,
				stage1->exhaustive_minpure_frequency_sum_gminus,
				exhaustive_index,minfreq);
      }
#ifdef GSNAP_DIAGONAL_ROW_STATS
      if (elt->diagstats_kmer_hits_n == 0 ||
	  (unsigned int) stage1->minus_npositions[querypos] < elt->diagstats_kmer_hits_min) {
	elt->diagstats_kmer_hits_min = (unsigned int) stage1->minus_npositions[querypos];
      }
      if ((unsigned int) stage1->minus_npositions[querypos] > elt->diagstats_kmer_hits_max) {
	elt->diagstats_kmer_hits_max = (unsigned int) stage1->minus_npositions[querypos];
      }
      elt->diagstats_kmer_hits_sum +=
	(unsigned long long) stage1->minus_npositions[querypos];
      elt->diagstats_kmer_hits_n++;
      if (minfreq_valid_p == true) {
	if (elt->diagstats_minfreq_n == 0 || minfreq < elt->diagstats_minfreq_min) {
	  elt->diagstats_minfreq_min = minfreq;
	}
	if (minfreq > elt->diagstats_minfreq_max) {
	  elt->diagstats_minfreq_max = minfreq;
	}
	elt->diagstats_minfreq_sum += (unsigned long long) minfreq;
	elt->diagstats_minfreq_n++;
      }
#endif
    }
#else
    if (stage1->exhaustive_minpure_specificity_gminus != NULL &&
	Minpure_table_lookup(&minfreq,genome_minpure,
			     stage1->revcomp_oligos[querypos],index1part) == true) {
      for (k = 0; k < nindices; k++) {
	exhaustive_index = indices[k];	/* Index of exhaustive */
	elt = eltarray[exhaustive_index];
	elt->univdiagonal = exhaustive[exhaustive_index];
	elt->qstart /*= qstarts[exhaustive_index]*/ = qstart;
	minpure_specificity_add(stage1->exhaustive_minpure_specificity_gminus,
				minpure_ns,
				stage1->exhaustive_minpure_specificity_min_gminus,
				stage1->exhaustive_minpure_specificity_max_gminus,
				stage1->exhaustive_minpure_frequency_min_gminus,
				stage1->exhaustive_minpure_frequency_max_gminus,
				stage1->exhaustive_minpure_frequency_sum_gminus,
				exhaustive_index,minfreq);
      }
    } else {
      for (k = 0; k < nindices; k++) {
	exhaustive_index = indices[k];	/* Index of exhaustive */
	elt = eltarray[exhaustive_index];
	elt->univdiagonal = exhaustive[exhaustive_index];
	elt->qstart /*= qstarts[exhaustive_index]*/ = qstart;
      }
    }
#endif
    all_indices[querypos+1] = &(indices[nindices]); /* Advance pointer */
  }

  /* Loop 2: Update qends; finds maximum qend for each univdiagonal */
  /* Descending querypos corresponds to ascending qend */
  for (querypos = query_lastpos; querypos >= 0; querypos--) {
    nindices = all_nindices[querypos];
    indices = all_indices[querypos];
    qend = querylength - querypos; /* minus */
    for (k = 0; k < nindices; k++) {
      exhaustive_index = indices[k];	/* Index of exhaustive */
      elt = eltarray[exhaustive_index];
      elt->qend /*= qends[exhaustive_index]*/ = qend;
    }
  }

  minpure_specificity_finalize(stage1->exhaustive_minpure_specificity_gminus,
			       minpure_ns,nexhaustive);
  if (minpure_ns != NULL && minpure_ns != stage1->exhaustive_minpure_observations_gminus) {
    FREE(minpure_ns);
  }

  FREE(all_nindices);
  FREE(all_indices);
  FREE(space);


#ifdef DEBUG9
  for (k = 0; k < nexhaustive; k++) {
    printf("%d: %d..%d %u\n",k,eltarray[k]->qstart,eltarray[k]->qend,eltarray[k]->univdiagonal);
  }
#endif

  return eltarray;
}



struct Segment_diagonal_budget_T {
  int elti;
  int width;
  int count;
  double score;
  Univcoord_T univdiagonal;
};


static int
segment_effective_kmers (Stage1_T stage1, int querylength, bool plusp) {
  int querypos, effective_kmers = 0;
  int query_lastpos = querylength - index1part;

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (stage1->validp[querypos] == false) {
      /* Skip */
    } else if (plusp == true && stage1->plus_npositions[querypos] > 0) {
      effective_kmers++;
    } else if (plusp == false && stage1->minus_npositions[querypos] > 0) {
      effective_kmers++;
    }
  }

  return effective_kmers;
}


static int
Segment_diagonal_cmp_width (const void *a, const void *b) {
  const struct Segment_diagonal_budget_T *x = (const struct Segment_diagonal_budget_T *) a;
  const struct Segment_diagonal_budget_T *y = (const struct Segment_diagonal_budget_T *) b;

  if (x->width < y->width) {
    return 1;
  } else if (x->width > y->width) {
    return -1;
  } else if (x->count < y->count) {
    return 1;
  } else if (x->count > y->count) {
    return -1;
  } else if (x->score < y->score) {
    return 1;
  } else if (x->score > y->score) {
    return -1;
  } else if (x->univdiagonal > y->univdiagonal) {
    return 1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else {
    return x->elti - y->elti;
  }
}


static int
Segment_diagonal_cmp_count (const void *a, const void *b) {
  const struct Segment_diagonal_budget_T *x = (const struct Segment_diagonal_budget_T *) a;
  const struct Segment_diagonal_budget_T *y = (const struct Segment_diagonal_budget_T *) b;

  if (x->count < y->count) {
    return 1;
  } else if (x->count > y->count) {
    return -1;
  } else if (x->width < y->width) {
    return 1;
  } else if (x->width > y->width) {
    return -1;
  } else if (x->score < y->score) {
    return 1;
  } else if (x->score > y->score) {
    return -1;
  } else if (x->univdiagonal > y->univdiagonal) {
    return 1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else {
    return x->elti - y->elti;
  }
}


static int
Segment_diagonal_cmp_score (const void *a, const void *b) {
  const struct Segment_diagonal_budget_T *x = (const struct Segment_diagonal_budget_T *) a;
  const struct Segment_diagonal_budget_T *y = (const struct Segment_diagonal_budget_T *) b;

  if (x->score < y->score) {
    return 1;
  } else if (x->score > y->score) {
    return -1;
  } else if (x->width < y->width) {
    return 1;
  } else if (x->width > y->width) {
    return -1;
  } else if (x->count < y->count) {
    return 1;
  } else if (x->count > y->count) {
    return -1;
  } else if (x->univdiagonal > y->univdiagonal) {
    return 1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else {
    return x->elti - y->elti;
  }
}


static int
Segment_diagonal_cmp_coord (const void *a, const void *b) {
  const struct Segment_diagonal_budget_T *x = (const struct Segment_diagonal_budget_T *) a;
  const struct Segment_diagonal_budget_T *y = (const struct Segment_diagonal_budget_T *) b;

  if (x->univdiagonal > y->univdiagonal) {
    return 1;
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else {
    return x->elti - y->elti;
  }
}


static int
segment_diagonal_mark_top (bool *keepp, struct Segment_diagonal_budget_T *candidates,
			   int ncandidates, int ntop, int nkept, int budget) {
  int i, elti, limit;

  if (ntop > ncandidates) {
    limit = ncandidates;
  } else {
    limit = ntop;
  }

  for (i = 0; i < limit && nkept < budget; i++) {
    elti = candidates[i].elti;
    if (keepp[elti] == false) {
      keepp[elti] = true;
      nkept++;
    }
  }

  return nkept;
}


static int
segment_diagonal_mark_coordinate_spread (bool *keepp,
					 struct Segment_diagonal_budget_T *candidates,
					 int ncandidates, int nkept, int budget) {
  int i, rank, j, target, remaining = 0, to_add;

  if (nkept >= budget) {
    return nkept;
  }

  for (i = 0; i < ncandidates; i++) {
    if (keepp[candidates[i].elti] == false) {
      remaining++;
    }
  }

  to_add = budget - nkept;
  if (to_add > remaining) {
    to_add = remaining;
  }
  if (to_add <= 0) {
    return nkept;
  }

  qsort(candidates,ncandidates,sizeof(struct Segment_diagonal_budget_T),
	Segment_diagonal_cmp_coord);

  j = 0;
  rank = 0;
  for (i = 0; i < ncandidates && j < to_add; i++) {
    if (keepp[candidates[i].elti] == false) {
      target = (int) (((long long) j * (long long) remaining)/(long long) to_add);
      if (rank == target) {
	keepp[candidates[i].elti] = true;
	nkept++;
	j++;
      }
      rank++;
    }
  }

  return nkept;
}


static int
segment_diagonal_apply_adaptive_hybrid (bool *selectedp, Univcoord_T *exhaustive,
			       Singleelt_T *eltarray, int *counts, int nexhaustive,
			       int nselected, int selected_by_width,
			       int selected_by_count, int querylength,
			       int effective_kmers, int min_width, int max_width,
			       int max_count, const char *label, const char *strand) {
  Singleelt_T elt;
  struct Segment_diagonal_budget_T *candidates, *candidate;
  int i, k, width, count, effective, floor_keep, budget, selected_after;
  int keep_top_width, keep_top_count, keep_top_score, strong_support = 0, weak_support = 0;
  int weak_width_fanout, triggered, nkept;
  double width_fraction, max_width_fraction, width_norm, count_norm;
  double density, max_density = 0.0, density_norm, support_fraction;
  double selected_per_kmer, selected_by_width_per_kmer, selected_by_count_per_kmer;
  double strong_support_per_kmer, weak_support_per_kmer, weak_fanout_limit, max_count_fraction;
  bool *keepp;

  if (nselected <= 0) {
    return nselected;
  }

  if (effective_kmers <= 0) {
    effective = 1;
  } else {
    effective = effective_kmers;
  }

  candidates = (struct Segment_diagonal_budget_T *)
    MALLOC(nselected * sizeof(struct Segment_diagonal_budget_T));

  k = 0;
  for (i = 0; i < nexhaustive; i++) {
    if (selectedp[i] == true) {
      elt = eltarray[i];
      width = elt->qend - elt->qstart;
      count = counts[i];
      density = width <= 0 ? 0.0 : (double) count/(double) width;
      if (density > max_density) {
	max_density = density;
      }

      candidate = &(candidates[k++]);
      candidate->elti = i;
      candidate->width = width;
      candidate->count = count;
      candidate->score = 0.0;
      candidate->univdiagonal = exhaustive[i];
    }
  }

  max_width_fraction = querylength <= 0 ? 0.0 : (double) max_width/(double) querylength;
  max_count_fraction = (double) max_count/(double) effective;
  for (k = 0; k < nselected; k++) {
    candidate = &(candidates[k]);
    width_fraction = querylength <= 0 ? 0.0 : (double) candidate->width/(double) querylength;
    width_norm = max_width_fraction <= 0.0 ? 0.0 : width_fraction/max_width_fraction;
    count_norm = max_count <= 0 ? 0.0 : (double) candidate->count/(double) max_count;
    density = candidate->width <= 0 ? 0.0 : (double) candidate->count/(double) candidate->width;
    density_norm = max_density <= 0.0 ? 0.0 : density/max_density;
    candidate->score = 0.60 * width_norm + 0.30 * count_norm + 0.10 * density_norm;

    support_fraction = (double) candidate->count/(double) effective;
    if (support_fraction >= adaptive_hybrid_strong_support_fraction) {
      strong_support++;
    }
    if (support_fraction >= stats_weak_support_fraction) {
      weak_support++;
    }
  }

  floor_keep = adaptive_hybrid_floor_per_kmer * effective;
  if (floor_keep < 1) {
    floor_keep = 1;
  }

  budget = adaptive_hybrid_budget_per_kmer * effective;
  if (budget < floor_keep) {
    budget = floor_keep;
  }
  if (adaptive_hybrid_hard_max > 0 && budget > adaptive_hybrid_hard_max) {
    budget = adaptive_hybrid_hard_max;
  }
  if (budget < 1) {
    budget = 1;
  }
  if (budget > nselected) {
    budget = nselected;
  }
  if (floor_keep > budget) {
    floor_keep = budget;
  }

  selected_per_kmer = (double) nselected/(double) effective;
  selected_by_width_per_kmer = (double) selected_by_width/(double) effective;
  selected_by_count_per_kmer = (double) selected_by_count/(double) effective;
  strong_support_per_kmer = (double) strong_support/(double) effective;
  weak_support_per_kmer = (double) weak_support/(double) effective;
  weak_fanout_limit = selected_by_width_per_kmer /
    (double) adaptive_hybrid_weak_support_divisor;

  weak_width_fanout =
    selected_by_width_per_kmer > (double) adaptive_hybrid_weak_width_per_kmer &&
    selected_by_width >= 4 * selected_by_count &&
    strong_support_per_kmer < weak_fanout_limit ? 1 : 0;

  triggered =
    nselected > budget ||
    ((double) nselected > 0.05 * (double) nexhaustive &&
     selected_per_kmer > (double) adaptive_hybrid_weak_width_per_kmer) ||
    weak_width_fanout == 1 ? 1 : 0;

  selected_after = triggered == 1 ? budget : nselected;

#ifdef GSNAP_DEBUG_SEGMENT_DIAGONAL_ADAPTIVE_HYBRID
  fprintf(stderr,
	  "[SEGMENT_DIAGONAL_ADAPTIVE_HYBRID] label=%s strand=%s triggered=%d "
	  "effective_kmers=%d exhaustive=%d selected_before=%d selected_after=%d "
	  "dropped=%d budget=%d floor_keep=%d hard_max=%d selected_per_kmer=%.2f "
	  "selected_by_width=%d selected_by_count=%d selected_by_width_per_kmer=%.2f "
	  "selected_by_count_per_kmer=%.2f strong_support=%d strong_support_per_kmer=%.2f "
	  "weak_support=%d weak_support_per_kmer=%.2f "
	  "weak_width_fanout=%d min_width=%d max_width=%d max_count=%d max_count_fraction=%.3f\n",
	  label,strand,triggered,
	  effective_kmers,nexhaustive,nselected,selected_after,
	  nselected - selected_after,budget,floor_keep,adaptive_hybrid_hard_max,
	  selected_per_kmer,selected_by_width,selected_by_count,
	  selected_by_width_per_kmer,selected_by_count_per_kmer,
	  strong_support,strong_support_per_kmer,weak_support,weak_support_per_kmer,
	  weak_width_fanout,
	  min_width,max_width,max_count,max_count_fraction);
#endif

  if (triggered == 0) {
    FREE(candidates);
    return nselected;
  }

  keepp = (bool *) CALLOC(nexhaustive,sizeof(bool));
  keep_top_width = weak_width_fanout == 1 ? budget / 8 : budget / 4;
  keep_top_count = weak_width_fanout == 1 ? budget / 3 : budget / 4;
  keep_top_score = budget - keep_top_width - keep_top_count;
  if (keep_top_score < 0) {
    keep_top_score = 0;
  }

  nkept = 0;
  qsort(candidates,nselected,sizeof(struct Segment_diagonal_budget_T),
	Segment_diagonal_cmp_width);
  nkept = segment_diagonal_mark_top(keepp,candidates,nselected,keep_top_width,nkept,budget);

  qsort(candidates,nselected,sizeof(struct Segment_diagonal_budget_T),
	Segment_diagonal_cmp_count);
  nkept = segment_diagonal_mark_top(keepp,candidates,nselected,keep_top_count,nkept,budget);

  qsort(candidates,nselected,sizeof(struct Segment_diagonal_budget_T),
	Segment_diagonal_cmp_score);
  nkept = segment_diagonal_mark_top(keepp,candidates,nselected,keep_top_score,nkept,budget);
  if (nkept < floor_keep) {
    nkept = segment_diagonal_mark_top(keepp,candidates,nselected,nselected,nkept,floor_keep);
  }
  if (nkept < budget) {
    nkept = segment_diagonal_mark_coordinate_spread(keepp,candidates,nselected,nkept,budget);
  }
  if (nkept < budget) {
    qsort(candidates,nselected,sizeof(struct Segment_diagonal_budget_T),
	  Segment_diagonal_cmp_score);
    nkept = segment_diagonal_mark_top(keepp,candidates,nselected,nselected,nkept,budget);
  }

  for (i = 0; i < nexhaustive; i++) {
    selectedp[i] = keepp[i];
  }

  FREE(keepp);
  FREE(candidates);

  return nkept;
}


static int
select_segment_elts (Univcoord_T **_univdiagonals, Auxinfo_T **auxinfo,
		     Pathstore_T **pathstores, Univcoord_T *exhaustive,
		     Singleelt_T *eltarray, int *counts, int nexhaustive,
		     int querylength, int effective_kmers,
		     const char *label, const char *strand,
		     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {
  Singleelt_T elt;
  bool *selectedp;
  int i, k, width, max_width = 0, min_width, max_count = 0, nbest = 0, selected_before;
  int selected_by_width = 0, selected_by_count = 0;

  if (nexhaustive == 0) {
    *_univdiagonals = (Univcoord_T *) NULL;
    *auxinfo = (Auxinfo_T *) NULL;
    *pathstores = (Pathstore_T *) NULL;
    return 0;
  }

  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    if ((width = elt->qend - elt->qstart) > max_width) {
      max_width = width;
    }
    if (counts[i] > max_count) {
      max_count = counts[i];
    }
  }

  min_width = max_width - (index1interval - 1);
  if (min_width < 0) {
    min_width = 0;
  }

#ifdef GSNAP_ACCEL_PROFILE_STATS
  segment_diagonal_report_distribution(exhaustive,eltarray,counts,nexhaustive,
				       min_width,max_count,label,strand);
#endif

  switch (segment_diagonal_algorithm) {
  case SEGMENT_DIAGONAL_ALGORITHM_WIDTH_OR_COUNT:
    for (i = 0; i < nexhaustive; i++) {
      elt = eltarray[i];
      width = elt->qend - elt->qstart;
      if (width >= min_width || counts[i] >= max_count) {
	nbest++;
	if (width >= min_width) {
	  selected_by_width++;
	}
	if (counts[i] >= max_count) {
	  selected_by_count++;
	}
      }
    }

    *_univdiagonals = (Univcoord_T *) MALLOC_ALIGN(nbest * sizeof(Univcoord_T));
    *auxinfo = (Auxinfo_T *) MALLOC(nbest * sizeof(Auxinfo_T));
    *pathstores = (Pathstore_T *) MALLOC(nbest * sizeof(Pathstore_T));

    k = 0;
    for (i = 0; i < nexhaustive; i++) {
      elt = eltarray[i];
      width = elt->qend - elt->qstart;
      if (width >= min_width || counts[i] >= max_count) {
	(*_univdiagonals)[k] = exhaustive[i];
	(*auxinfo)[k] = Auxinfo_new_elt(/*method*/SEGMENT1,elt->qstart,elt->qend,/*elti*/i,
					auxinfopool);
	(*pathstores)[k] = Pathstore_new(pathstorepool);
	k++;
      }
    }

    debug9(printf("Selected %d segment candidates from %d exhaustive diagonals (min width %d, max count %d)\n",
		  nbest,nexhaustive,min_width,max_count));
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    segment_diagonal_print_select_stats(label,strand,exhaustive,eltarray,counts,nexhaustive,
					nbest,min_width,max_width,max_count);
#endif
    Accel_profile_segment_select(segment_diagonal_algorithm_name(),label,strand,
				 nexhaustive,/*selected_before*/nbest,
				 /*selected_after*/nbest,selected_by_width,selected_by_count,
				 min_width,max_width,max_count,effective_kmers);

    return nbest;

  case SEGMENT_DIAGONAL_ALGORITHM_ADAPTIVE_HYBRID:
    break;

  default:
    abort();
  }

  selectedp = (bool *) CALLOC(nexhaustive,sizeof(bool));
  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    if (width >= min_width || counts[i] >= max_count) {
      selectedp[i] = true;
      nbest++;
      if (width >= min_width) {
	selected_by_width++;
      }
      if (counts[i] >= max_count) {
	selected_by_count++;
      }
    }
  }

  selected_before = nbest;
  nbest = segment_diagonal_apply_adaptive_hybrid(selectedp,exhaustive,eltarray,counts,nexhaustive,
					nbest,selected_by_width,selected_by_count,
					querylength,effective_kmers,min_width,max_width,
					max_count,label,strand);

  *_univdiagonals = (Univcoord_T *) MALLOC_ALIGN(nbest * sizeof(Univcoord_T));
  *auxinfo = (Auxinfo_T *) MALLOC(nbest * sizeof(Auxinfo_T));
  *pathstores = (Pathstore_T *) MALLOC(nbest * sizeof(Pathstore_T));

  k = 0;
  for (i = 0; i < nexhaustive; i++) {
    if (selectedp[i] == true) {
      elt = eltarray[i];
      (*_univdiagonals)[k] = exhaustive[i];
      (*auxinfo)[k] = Auxinfo_new_elt(/*method*/SEGMENT1,elt->qstart,elt->qend,/*elti*/i,
				      auxinfopool);
      (*pathstores)[k] = Pathstore_new(pathstorepool);
      k++;
    }
  }

  debug9(printf("Selected %d segment candidates from %d exhaustive diagonals (min width %d, max count %d)\n",
		nbest,nexhaustive,min_width,max_count));
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  segment_diagonal_print_select_stats(label,strand,exhaustive,eltarray,counts,nexhaustive,
				      nbest,min_width,max_width,max_count);
#endif
  Accel_profile_segment_select(segment_diagonal_algorithm_name(),label,strand,
			       nexhaustive,selected_before,/*selected_after*/nbest,
			       selected_by_width,selected_by_count,min_width,max_width,max_count,
			       effective_kmers);

  FREE(selectedp);

  return nbest;
}


/* Needed before we call get_eltarray_plus and get_eltarray_gminus, or
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
      if (segment_repetitive_kmer_filter_p == true && repetitive_ef64 != NULL &&
	  EF64_presentp(stage1->forward_oligos[querypos],repetitive_ef64) == true) {
	/* Set npositions to be 0 */
	stage1->plus_npositions[querypos] = 0;
      } else if (segment_max_kmer_positions > 0 &&
		 stage1->plus_npositions[querypos] > segment_max_kmer_positions) {
	stage1->plus_npositions[querypos] = 0;
      } else {
	*total_npositions_plus += stage1->plus_npositions[querypos];
      }
      if (segment_repetitive_kmer_filter_p == true && repetitive_ef64 != NULL &&
	  EF64_presentp(stage1->revcomp_oligos[querypos],repetitive_ef64) == true) {
	/* Set npositions to be 0 */
	stage1->minus_npositions[querypos] = 0;
      } else if (segment_max_kmer_positions > 0 &&
		 stage1->minus_npositions[querypos] > segment_max_kmer_positions) {
	stage1->minus_npositions[querypos] = 0;
      } else {
	*total_npositions_minus += stage1->minus_npositions[querypos];
      }
    }
  }

  return;
}


List_T
Kmer_solve_eltarray_left (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			  Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T chains;
  Singleelt_T *copy;
  int i, jlow, jhigh, n;
  Univcoord_T last_univdiagonal, orig_univdiagonal, univdiagonal;
  Univcoord_T max_total_distance = max_nsplices * positive_gap_distance;

  i = elti;
  orig_univdiagonal = eltarray[i]->univdiagonal;

  jlow = i - 1;
  last_univdiagonal = orig_univdiagonal;
  while (jlow >= 0 &&
	 last_univdiagonal - (univdiagonal = eltarray[jlow]->univdiagonal) <= positive_gap_distance &&
	 orig_univdiagonal - univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jlow--;
  }

  jhigh = i + 1;
  last_univdiagonal = orig_univdiagonal;
  while (jhigh < nexhaustive &&
	 (univdiagonal = eltarray[jhigh]->univdiagonal) - last_univdiagonal <= positive_gap_distance &&
	 univdiagonal - orig_univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jhigh++;
  }

  jlow++;
  n = jhigh - jlow;

#if 0
  for (int k = jlow; k < jhigh; k++) {
    printf("%d [%d..%d]: %d..%d %u\n",
	   k,jlow,jhigh - 1,eltarray[k]->qstart,eltarray[k]->qend,eltarray[k]->univdiagonal);
  }
#endif

  /* Need to copy eltarray so we don't disturb the original */
  copy = (Singleelt_T *) MALLOC(n * sizeof(Singleelt_T));
  memcpy(copy,&(eltarray[jlow]),n * sizeof(Singleelt_T));
  qsort(copy,n,sizeof(Singleelt_T),Singleelt_cmp);

  chains = Singleelt_left_chains(copy,n,max_nsplices,univdiagpool,listpool);

  FREE(copy);

  return chains;
}


List_T
Kmer_solve_eltarray_right (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			   Univdiagpool_T univdiagpool, Listpool_T listpool) {
  List_T chains;
  Singleelt_T *copy;
  int i, jlow, jhigh, n;
  Univcoord_T last_univdiagonal, orig_univdiagonal, univdiagonal;
  Univcoord_T max_total_distance = max_nsplices * positive_gap_distance;

  i = elti;
  orig_univdiagonal = eltarray[i]->univdiagonal;

  jlow = i - 1;
  last_univdiagonal = orig_univdiagonal;
  while (jlow >= 0 &&
	 last_univdiagonal - (univdiagonal = eltarray[jlow]->univdiagonal) <= positive_gap_distance &&
	 orig_univdiagonal - univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jlow--;
  }

  jhigh = i + 1;
  last_univdiagonal = orig_univdiagonal;
  while (jhigh < nexhaustive &&
	 (univdiagonal = eltarray[jhigh]->univdiagonal) - last_univdiagonal <= positive_gap_distance &&
	 univdiagonal - orig_univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jhigh++;
  }

  jlow++;
  n = jhigh - jlow;

#if 0
  for (int k = jlow; k < jhigh; k++) {
    printf("%d [%d..%d]: %d..%d %u\n",
	   k,jlow,jhigh - 1,eltarray[k]->qstart,eltarray[k]->qend,eltarray[k]->univdiagonal);
  }
#endif

  /* Need to copy eltarray so we don't disturb the original */
  copy = (Singleelt_T *) MALLOC(n * sizeof(Singleelt_T));
  memcpy(copy,&(eltarray[jlow]),n * sizeof(Singleelt_T));
  qsort(copy,n,sizeof(Singleelt_T),Singleelt_cmp);

  chains = Singleelt_right_chains(copy,n,max_nsplices,univdiagpool,listpool);

  FREE(copy);

  return chains;
}

void
Kmer_apply_eltarray_left (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			  Listpool_T listpool,
			  Singleelt_chain_skip_T skip_chain_p,
			  bool equiv_prev_prune_p, int max_chain_window,
			  bool materialized_order_p,
			  Singleelt_chain_apply_T apply_chain, void *data) {
  Singleelt_T *copy;
  int i, jlow, jhigh, n;
  Univcoord_T last_univdiagonal, orig_univdiagonal, univdiagonal;
  Univcoord_T max_total_distance = max_nsplices * positive_gap_distance;

  i = elti;
  orig_univdiagonal = eltarray[i]->univdiagonal;

  jlow = i - 1;
  last_univdiagonal = orig_univdiagonal;
  while (jlow >= 0 &&
	 last_univdiagonal - (univdiagonal = eltarray[jlow]->univdiagonal) <= positive_gap_distance &&
	 orig_univdiagonal - univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jlow--;
  }

  jhigh = i + 1;
  last_univdiagonal = orig_univdiagonal;
  while (jhigh < nexhaustive &&
	 (univdiagonal = eltarray[jhigh]->univdiagonal) - last_univdiagonal <= positive_gap_distance &&
	 univdiagonal - orig_univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jhigh++;
  }

  jlow++;
  n = jhigh - jlow;

  if (max_chain_window > 0 && n > max_chain_window) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    fprintf(stderr,
	    "[KMER_APPLY_LEFT_SKIP_WINDOW] elti=%d narray=%d max_chain_window=%d "
	    "jlow=%d jhigh=%d\n",
	    elti,n,max_chain_window,jlow,jhigh);
#endif
    return;
  }

  copy = (Singleelt_T *) MALLOC(n * sizeof(Singleelt_T));
  memcpy(copy,&(eltarray[jlow]),n * sizeof(Singleelt_T));
  qsort(copy,n,sizeof(Singleelt_T),Singleelt_cmp);

  Singleelt_left_chains_apply(copy,n,max_nsplices,listpool,
			      skip_chain_p,equiv_prev_prune_p,materialized_order_p,
			      apply_chain,data);

  FREE(copy);

  return;
}

void
Kmer_apply_eltarray_right (Singleelt_T *eltarray, int elti, int nexhaustive, int max_nsplices,
			   Listpool_T listpool,
			   Singleelt_chain_skip_T skip_chain_p,
			   bool equiv_prev_prune_p, int max_chain_window,
			   bool materialized_order_p,
			   Singleelt_chain_apply_T apply_chain, void *data) {
  Singleelt_T *copy;
  int i, jlow, jhigh, n;
  Univcoord_T last_univdiagonal, orig_univdiagonal, univdiagonal;
  Univcoord_T max_total_distance = max_nsplices * positive_gap_distance;

  i = elti;
  orig_univdiagonal = eltarray[i]->univdiagonal;

  jlow = i - 1;
  last_univdiagonal = orig_univdiagonal;
  while (jlow >= 0 &&
	 last_univdiagonal - (univdiagonal = eltarray[jlow]->univdiagonal) <= positive_gap_distance &&
	 orig_univdiagonal - univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jlow--;
  }

  jhigh = i + 1;
  last_univdiagonal = orig_univdiagonal;
  while (jhigh < nexhaustive &&
	 (univdiagonal = eltarray[jhigh]->univdiagonal) - last_univdiagonal <= positive_gap_distance &&
	 univdiagonal - orig_univdiagonal <= max_total_distance) {
    last_univdiagonal = univdiagonal;
    jhigh++;
  }

  jlow++;
  n = jhigh - jlow;

  if (max_chain_window > 0 && n > max_chain_window) {
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
    fprintf(stderr,
	    "[KMER_APPLY_RIGHT_SKIP_WINDOW] elti=%d narray=%d max_chain_window=%d "
	    "jlow=%d jhigh=%d\n",
	    elti,n,max_chain_window,jlow,jhigh);
#endif
    return;
  }

  copy = (Singleelt_T *) MALLOC(n * sizeof(Singleelt_T));
  memcpy(copy,&(eltarray[jlow]),n * sizeof(Singleelt_T));
  qsort(copy,n,sizeof(Singleelt_T),Singleelt_cmp);

  Singleelt_right_chains_apply(copy,n,max_nsplices,listpool,
			       skip_chain_p,equiv_prev_prune_p,materialized_order_p,
			       apply_chain,data);

  FREE(copy);

  return;
}


/* Optional genome minimizer-frequency sidecar, consulted only by the
   paired mate-rescue specificity policy and by the compile-gated per-diagonal
   stats emitter below.  NULL when no sidecar is present. */
void
Kmer_search_setup_minpure (Minpure_Table_T genome_minpure_in) {
  genome_minpure = genome_minpure_in;
  return;
}

static Univcoord_T
diagonal_genomic_coord (Univcoord_T univdiagonal, int querylength, int querypos) {
  int offset;

  if (querypos < 0) {
    querypos = 0;
  } else if (querypos > querylength) {
    querypos = querylength;
  }

  offset = querylength - querypos;
  if (univdiagonal < (Univcoord_T) offset) {
    return 0U;
  } else {
    return univdiagonal - (Univcoord_T) offset;
  }
}

static void
segment_minpure_diagonal_emit_debug (Univcoord_T *exhaustive, Singleelt_T *eltarray,
				     int *counts, int nexhaustive, int querylength,
				     bool first_read_p, bool plusp,
				     const char *read_name, const char *label,
				     int *observations,
				     unsigned int *frequency_mins,
				     unsigned int *frequency_maxs,
				     unsigned long long *frequency_sums,
				     unsigned int *specificity_avgs,
				     unsigned int *specificity_mins,
				     unsigned int *specificity_maxs) {
  Singleelt_T elt;
  Univcoord_T diagonal_left, diagonal_high, segment_start, segment_end;
  Univcoord_T chroffset, chrhigh;
  Chrnum_T chrnum;
  Chrpos_T chrpos0;
  char *chr;
  const char *strand = plusp == true ? "plus" : "minus";
  bool allocp;
  int i, width, n;
  double frequency_mean;

  if (minpure_diagonal_debug_p == false || nexhaustive <= 0) {
    return;
  }

  flockfile(stderr);
  if (minpure_diagonal_debug_header_printed_p == false) {
    fprintf(stderr,
	    "MINPURE_DIAG_HEADER\tread_name\tmate\tsearch_label\tstrand\tdiagonal_index"
	    "\tunivdiagonal\tchrnum\tchr\tchrpos0\tchrpos1"
	    "\tdiagonal_left_univcoord\tsegment_start_univcoord\tsegment_end_univcoord"
	    "\tquery_length\tqstart\tqend\twidth\tsupport_count"
	    "\tpure_minimizer_available\tpure_minimizer_observations"
	    "\tpure_minimizer_frequency_min\tpure_minimizer_frequency_mean"
	    "\tpure_minimizer_frequency_max"
	    "\tpure_minimizer_specificity_score_min"
	    "\tpure_minimizer_specificity_score_mean"
	    "\tpure_minimizer_specificity_score_max\n");
    minpure_diagonal_debug_header_printed_p = true;
  }

  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    n = observations != NULL ? observations[i] : 0;
    frequency_mean = n > 0 && frequency_sums != NULL ?
      (double) frequency_sums[i]/(double) n : 0.0;

    diagonal_left = diagonal_genomic_coord(exhaustive[i],querylength,/*querypos*/0);
    diagonal_high = exhaustive[i] > diagonal_left ? exhaustive[i] : diagonal_left + 1U;
    segment_start = diagonal_genomic_coord(exhaustive[i],querylength,elt->qstart);
    segment_end = diagonal_genomic_coord(exhaustive[i],querylength,elt->qend);

    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 diagonal_left,diagonal_high);
    if (chromosome_iit != NULL && chrnum > 0) {
      chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
    } else {
      chr = (char *) "";
      allocp = false;
    }
    chrpos0 = segment_start >= chroffset ? (Chrpos_T) (segment_start - chroffset) : 0U;

    fprintf(stderr,
	    "MINPURE_DIAG\t%s\t%d\t%s\t%s\t%d"
	    "\t%llu\t%d\t%s\t%u\t%u"
	    "\t%llu\t%llu\t%llu"
	    "\t%d\t%d\t%d\t%d\t%d"
	    "\t%d\t%d"
	    "\t%u\t%.2f\t%u"
	    "\t%u\t%u\t%u\n",
	    read_name != NULL ? read_name : "",first_read_p ? 1 : 2,label,strand,i,
	    (unsigned long long) exhaustive[i],chrnum,chr,
	    chrpos0,chrpos0 + 1U,
	    (unsigned long long) diagonal_left,
	    (unsigned long long) segment_start,
	    (unsigned long long) segment_end,
	    querylength,elt->qstart,elt->qend,width,counts[i],
	    genome_minpure != NULL ? 1 : 0,n,
	    n > 0 && frequency_mins != NULL ? frequency_mins[i] : 0U,
	    frequency_mean,
	    n > 0 && frequency_maxs != NULL ? frequency_maxs[i] : 0U,
	    n > 0 && specificity_mins != NULL ? specificity_mins[i] : 0U,
	    n > 0 && specificity_avgs != NULL ? specificity_avgs[i] : 0U,
	    n > 0 && specificity_maxs != NULL ? specificity_maxs[i] : 0U);

    if (allocp == true) {
      FREE(chr);
    }
  }

  funlockfile(stderr);
  return;
}

#ifdef GSNAP_DIAGONAL_ROW_STATS
/* Emit one TSV record for every raw exhaustive diagonal.  Minimizer frequency
   is aggregated only over the exact query k-mer occurrences that generated
   the diagonal, as recorded while mapping position streams to exhaustive[]. */
static void
segment_diagonal_emit_stats (Univcoord_T *exhaustive, Singleelt_T *eltarray,
			     int *counts, int nexhaustive, int querylength,
			     bool first_read_p, bool plusp,
			     const char *read_name, const char *label) {
  Singleelt_T elt;
  int i, width, kmer_span;
  unsigned long long prev_delta, next_delta;
  const char *strand = plusp == true ? "plus" : "minus";

  if (nexhaustive <= 0) {
    return;
  }

  flockfile(stderr);
  if (diagstats_header_printed_p == false) {
    fprintf(stderr,
	    "DIAGROW_HEADER\tread_name\tmate\tsearch_label\tstrand\tdiagonal_index"
	    "\tunivdiagonal\tprevious_diagonal_delta\tnext_diagonal_delta"
	    "\tquery_length\tindex_k\tqstart\tqend\twidth\tkmer_span"
	    "\tsupport_count\tsupport_density"
	    "\tsupporting_kmer_index_hits_min\tsupporting_kmer_index_hits_mean"
	    "\tsupporting_kmer_index_hits_max\tpure_minimizer_available\tpure_minimizer_m"
	    "\tpure_minimizer_observations\tpure_minimizer_context_frequency_min"
	    "\tpure_minimizer_context_frequency_mean"
	    "\tpure_minimizer_context_frequency_max\n");
    diagstats_header_printed_p = true;
  }

  for (i = 0; i < nexhaustive; i++) {
    elt = eltarray[i];
    width = elt->qend - elt->qstart;
    kmer_span = width >= index1part ? width - index1part + 1 : 0;
    prev_delta = i > 0 ?
      (unsigned long long) (exhaustive[i] - exhaustive[i-1]) : 0ULL;
    next_delta = i + 1 < nexhaustive ?
      (unsigned long long) (exhaustive[i+1] - exhaustive[i]) : 0ULL;
    fprintf(stderr,
	    "DIAGROW\t%s\t%d\t%s\t%s\t%d\t%llu\t%llu\t%llu"
	    "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.8f"
	    "\t%u\t%.2f\t%u\t%d\t%d\t%d\t%u\t%.2f\t%u\n",
	    read_name != NULL ? read_name : "",first_read_p ? 1 : 2,label,strand,i,
	    (unsigned long long) exhaustive[i],prev_delta,next_delta,
	    querylength,index1part,elt->qstart,elt->qend,width,kmer_span,counts[i],
	    kmer_span > 0 ? (double) counts[i]/(double) kmer_span : 0.0,
	    elt->diagstats_kmer_hits_n > 0 ? elt->diagstats_kmer_hits_min : 0U,
	    elt->diagstats_kmer_hits_n > 0 ?
	    (double) elt->diagstats_kmer_hits_sum/(double) elt->diagstats_kmer_hits_n : 0.0,
	    elt->diagstats_kmer_hits_max,
	    genome_minpure != NULL ? 1 : 0,
	    genome_minpure != NULL ? Minpure_table_m(genome_minpure) : 0,
	    elt->diagstats_minfreq_n,
	    elt->diagstats_minfreq_n > 0 ? elt->diagstats_minfreq_min : 0U,
	    elt->diagstats_minfreq_n > 0 ?
	    (double) elt->diagstats_minfreq_sum/(double) elt->diagstats_minfreq_n : 0.0,
	    elt->diagstats_minfreq_max);
  }
  funlockfile(stderr);
  return;
}
#endif


/* Rewrite of Segment_search */
/* _univdiagonals_gplus and _univdiagonals_gminus must be aligned to match Kmer_exact1 */

/* For segment1, _univdiagonals are the best univdiagonals.  For segment2, they are all univdiagonals */
/* Kmer_segment also assigns stage1->exhaustive_gplus and
   stage1->exhaustive_gminus, but they must be separate from
   univdiagonals */
void
Kmer_segment (Univcoord_T **_univdiagonals_gplus, Auxinfo_T **auxinfo_gplus,
	      Pathstore_T **pathstores_gplus, int *nunivdiagonals_gplus,
	      Univcoord_T **_univdiagonals_gminus, Auxinfo_T **auxinfo_gminus,
	      Pathstore_T **pathstores_gminus, int *nunivdiagonals_gminus,

	      Stage1_T stage1, int querylength, int max_nsplices, bool first_read_p,
	      const char *read_name, const char *label,

	      EF64_T repetitive_ef64, Univdiagpool_T univdiagpool,
	      Listpool_T listpool,
	      Univcoordlistpool_T univcoordlistpool,
	      Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool) {

  Singleelt_T *eltarray;
  int total_npositions_plus, total_npositions_minus;
  int effective_kmers_plus, effective_kmers_minus;
  debug9(printf("Entered Kmer_segment\n"));
  debug9(Stage1_dump(stage1,querylength));
  assert(stage1->all_positions_gen_filledp == true);

  /* total_npositions_plus and total_npositions_minus have duplicates */
  clear_repetitive_oligos(&total_npositions_plus,&total_npositions_minus,
			  stage1,querylength,repetitive_ef64);
  effective_kmers_plus = segment_effective_kmers(stage1,querylength,/*plusp*/true);
  effective_kmers_minus = segment_effective_kmers(stage1,querylength,/*plusp*/false);
  Accel_profile_segment_input(first_read_p,querylength,total_npositions_plus,total_npositions_minus,
			      effective_kmers_plus,effective_kmers_minus);
#ifdef GSNAP_DEBUG_SEGMENT_STREAM_STATS
  fprintf(stderr,
	  "[SEGMENT_DIAGONAL_INPUT] label=%s stage1=%p first_read=%d querylength=%d "
	  "max_nsplices=%d removed_by_position_filter=%d total_plus_after_repeat_clear=%d "
	  "total_minus_after_repeat_clear=%d effective_kmers_plus=%d effective_kmers_minus=%d\n",
	  label,stage1,first_read_p ? 1 : 0,querylength,max_nsplices,
	  stage1->segment_positions_filter_removed,
	  total_npositions_plus,total_npositions_minus,effective_kmers_plus,effective_kmers_minus);
#endif

  if (total_npositions_plus == 0) {
    *nunivdiagonals_gplus = stage1->nexhaustive_gplus = 0;
    *_univdiagonals_gplus = stage1->exhaustive_gplus = (Univcoord_T *) NULL;
    *auxinfo_gplus = (Auxinfo_T *) NULL;
    *pathstores_gplus = (Pathstore_T *) NULL;

  } else {
    /* Sets stage1->exhaustive_gplus and stage1->exhaustive_counts_gplus */
    stage1->nexhaustive_gplus = store_exhaustive_gplus(stage1,querylength,label);

    stage1->exhaustive_eltarray_gplus = eltarray =
      get_eltarray_gplus(stage1->exhaustive_gplus,stage1->nexhaustive_gplus,
			 stage1,querylength,total_npositions_plus);
#ifdef GSNAP_DIAGONAL_ROW_STATS
    segment_diagonal_emit_stats(stage1->exhaustive_gplus,eltarray,
				stage1->exhaustive_counts_gplus,stage1->nexhaustive_gplus,
				querylength,first_read_p,/*plusp*/true,read_name,label);
#endif
    segment_minpure_diagonal_emit_debug(stage1->exhaustive_gplus,eltarray,
					stage1->exhaustive_counts_gplus,stage1->nexhaustive_gplus,
					querylength,first_read_p,/*plusp*/true,read_name,label,
					stage1->exhaustive_minpure_observations_gplus,
					stage1->exhaustive_minpure_frequency_min_gplus,
					stage1->exhaustive_minpure_frequency_max_gplus,
					stage1->exhaustive_minpure_frequency_sum_gplus,
					stage1->exhaustive_minpure_specificity_gplus,
					stage1->exhaustive_minpure_specificity_min_gplus,
					stage1->exhaustive_minpure_specificity_max_gplus);

    *nunivdiagonals_gplus =
      select_segment_elts(&(*_univdiagonals_gplus),&(*auxinfo_gplus),&(*pathstores_gplus),
			  stage1->exhaustive_gplus,eltarray,stage1->exhaustive_counts_gplus,
			  stage1->nexhaustive_gplus,querylength,effective_kmers_plus,
			  label,"plus",auxinfopool,pathstorepool);
  }


  if (total_npositions_minus == 0) {
    *nunivdiagonals_gminus = stage1->nexhaustive_gminus = 0;
    *_univdiagonals_gminus = stage1->exhaustive_gminus = (Univcoord_T *) NULL;
    *auxinfo_gminus = (Auxinfo_T *) NULL;
    *pathstores_gminus = (Pathstore_T *) NULL;

  } else {
    /* Sets stage1->exhaustive_gminus and stage1->exhaustive_counts_gminus */
    stage1->nexhaustive_gminus = store_exhaustive_gminus(stage1,querylength,label);

    stage1->exhaustive_eltarray_gminus = eltarray =
      get_eltarray_gminus(stage1->exhaustive_gminus,stage1->nexhaustive_gminus,
			  stage1,querylength,total_npositions_minus);
#ifdef GSNAP_DIAGONAL_ROW_STATS
    segment_diagonal_emit_stats(stage1->exhaustive_gminus,eltarray,
				stage1->exhaustive_counts_gminus,stage1->nexhaustive_gminus,
				querylength,first_read_p,/*plusp*/false,read_name,label);
#endif
    segment_minpure_diagonal_emit_debug(stage1->exhaustive_gminus,eltarray,
					stage1->exhaustive_counts_gminus,stage1->nexhaustive_gminus,
					querylength,first_read_p,/*plusp*/false,read_name,label,
					stage1->exhaustive_minpure_observations_gminus,
					stage1->exhaustive_minpure_frequency_min_gminus,
					stage1->exhaustive_minpure_frequency_max_gminus,
					stage1->exhaustive_minpure_frequency_sum_gminus,
					stage1->exhaustive_minpure_specificity_gminus,
					stage1->exhaustive_minpure_specificity_min_gminus,
					stage1->exhaustive_minpure_specificity_max_gminus);

    *nunivdiagonals_gminus =
      select_segment_elts(&(*_univdiagonals_gminus),&(*auxinfo_gminus),&(*pathstores_gminus),
			  stage1->exhaustive_gminus,eltarray,stage1->exhaustive_counts_gminus,
			  stage1->nexhaustive_gminus,querylength,effective_kmers_minus,
			  label,"minus",auxinfopool,pathstorepool);
  }

  debug9(printf("Exiting Kmer_segment\n"));

  return;
}



void
Kmer_search_setup (int index1part_in, int index1interval_in, EF64_T chromosome_ef64_in,
		   Univ_IIT_T chromosome_iit_in,
		   Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
		   Univcoord_T genomelength_in, 
		   int max_insertionlen, int max_deletionlen, Chrpos_T shortsplicedist,
		   bool splicingp_in, Segment_diagonal_algorithm_T segment_diagonal_algorithm_in,
		   int adaptive_hybrid_floor_per_kmer_in,
		   int adaptive_hybrid_budget_per_kmer_in,
		   int adaptive_hybrid_hard_max_in,
		   int adaptive_hybrid_weak_width_per_kmer_in,
		   int adaptive_hybrid_weak_support_divisor_in,
		   double adaptive_hybrid_strong_support_fraction_in,
		   double stats_weak_support_fraction_in,
		   bool segment_repetitive_kmer_filter_p_in,
		   int segment_max_kmer_positions_in,
		   bool minpure_diagonal_debug_p_in) {

  index1part = index1part_in;
  index1interval = index1interval_in;

  chromosome_ef64 = chromosome_ef64_in;
  chromosome_iit = chromosome_iit_in;
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  genomelength = genomelength_in;

  negative_gap_distance = (Chrpos_T) max_insertionlen;
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  
  splicingp = splicingp_in;

  segment_diagonal_algorithm = segment_diagonal_algorithm_in;
  adaptive_hybrid_floor_per_kmer = adaptive_hybrid_floor_per_kmer_in;
  adaptive_hybrid_budget_per_kmer = adaptive_hybrid_budget_per_kmer_in;
  adaptive_hybrid_hard_max = adaptive_hybrid_hard_max_in;
  adaptive_hybrid_weak_width_per_kmer = adaptive_hybrid_weak_width_per_kmer_in;
  adaptive_hybrid_weak_support_divisor = adaptive_hybrid_weak_support_divisor_in;
  adaptive_hybrid_strong_support_fraction = adaptive_hybrid_strong_support_fraction_in;
  stats_weak_support_fraction = stats_weak_support_fraction_in;
  segment_repetitive_kmer_filter_p = segment_repetitive_kmer_filter_p_in;
  segment_max_kmer_positions = segment_max_kmer_positions_in;
  minpure_diagonal_debug_p = minpure_diagonal_debug_p_in;

  if (adaptive_hybrid_floor_per_kmer < 1) {
    adaptive_hybrid_floor_per_kmer = 1;
  }
  if (adaptive_hybrid_budget_per_kmer < adaptive_hybrid_floor_per_kmer) {
    adaptive_hybrid_budget_per_kmer = adaptive_hybrid_floor_per_kmer;
  }
  if (adaptive_hybrid_hard_max < 0) {
    adaptive_hybrid_hard_max = 0;
  }
  if (adaptive_hybrid_weak_width_per_kmer < 1) {
    adaptive_hybrid_weak_width_per_kmer = 1;
  }
  if (adaptive_hybrid_weak_support_divisor < 1) {
    adaptive_hybrid_weak_support_divisor = 1;
  }
  if (adaptive_hybrid_strong_support_fraction < 0.0) {
    adaptive_hybrid_strong_support_fraction = 0.0;
  }
  if (stats_weak_support_fraction < 0.0) {
    stats_weak_support_fraction = 0.0;
  }
  if (segment_max_kmer_positions < 0) {
    segment_max_kmer_positions = 0;
  }

  return;
}
