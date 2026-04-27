static char rcsid[] = "$Id$";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Code is modified from the SIMDCompressionAndIntersection package,
   by Daniel Lemire, Leonid Boytsov, and Nathan Kurz, file
   intersection.cpp.  */

/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 */

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "bool.h"
#include "intersect-concordance-uint4.h"

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#endif


#define PICK_BEST_DISTANCE 1
#define MAX_HITS_PER_RARE 5

#define EPI16_MAX 0x8000
#define EPI32_MAX 0x80000000


#ifdef HAVE_SSE4_1
#define LOAD_SI128 _mm_lddqu_si128
#else
#define LOAD_SI128 _mm_loadu_si128
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* initialize_match_values */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Pick best distance */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif



#if defined(HAVE_AVX2) && defined(DEBUG)
static void
print_vector_signed (__m256i x) {
  printf("%10d %10d %10d %10d %10d %10d %10d %10d\n",
	 _mm256_extract_epi32(x,0) + EPI32_MAX,
	 _mm256_extract_epi32(x,1) + EPI32_MAX,
	 _mm256_extract_epi32(x,2) + EPI32_MAX,
	 _mm256_extract_epi32(x,3) + EPI32_MAX,
	 _mm256_extract_epi32(x,4) + EPI32_MAX,
	 _mm256_extract_epi32(x,5) + EPI32_MAX,
	 _mm256_extract_epi32(x,6) + EPI32_MAX,
	 _mm256_extract_epi32(x,7) + EPI32_MAX);
  return;
}

static void
print_vector_unsigned (__m256i x) {
  printf("%10u %10u %10u %10u %10u %10u %10u %10u\n",
	 _mm256_extract_epi32(x,0),
	 _mm256_extract_epi32(x,1),
	 _mm256_extract_epi32(x,2),
	 _mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),
	 _mm256_extract_epi32(x,5),
	 _mm256_extract_epi32(x,6),
	 _mm256_extract_epi32(x,7));
  return;
}

static void
print_vector_hex (__m256i x) {
  printf("%08X %08X %08X %08X %08X %08X %08X %08X\n",
	 _mm256_extract_epi32(x,0),
	 _mm256_extract_epi32(x,1),
	 _mm256_extract_epi32(x,2),
	 _mm256_extract_epi32(x,3),
	 _mm256_extract_epi32(x,4),
	 _mm256_extract_epi32(x,5),
	 _mm256_extract_epi32(x,6),
	 _mm256_extract_epi32(x,7));
  return;
}
#endif


#ifdef HAVE_AVX512
#define NELTS 16
#define MAX_IDX 65536		/* 2^16 */
#else
#define NELTS 8
#define MAX_IDX 256		/* 2^8 */
#endif

static int match_start[MAX_IDX];
static int match_n[MAX_IDX];

static void
initialize_match_values () {
  int high_bit, low_bit, bit;
  int idx;

  for (idx = 0; idx < MAX_IDX; idx++) {
    match_start[idx] = -1;
    match_n[idx] = 0;
  }

  for (high_bit = 0; high_bit < NELTS; high_bit++) {
    for (low_bit = 0; low_bit <= high_bit; low_bit++) {
      idx = 0;
      for (bit = low_bit; bit <= high_bit; bit++) {
	idx |= (1 << bit);
      }
      match_start[idx] = low_bit;
      match_n[idx] = high_bit - low_bit + 1;
      debug1(printf("%d %d %d\n",idx,match_start[idx],match_n[idx]));
    }
  }
      
  return;
}



/* Fast scalar scheme designed by N. Kurz. */
/* Delta is required to be added to A to become equivalent with B */
/* Desired result should be A + diagtermA (== B + diagtermB) */
static int *
scalar (int *indices, int **init_indices, int **end_indices,
	const int A_index, const int B_index,
	const unsigned int *A, const unsigned int *endA,
	const unsigned int *B, const unsigned int *endB,
	const unsigned int A_slop_neg, const unsigned int A_slop_pos,
	bool A_first_p) {

  const unsigned int *initA = A;
  const unsigned int *initB = B;

  const unsigned int *B_save;
  int allocation, nindices;
  int *more_indices;

  debug(printf("Entered scalar with %lu and %lu entries\n",endA - A,endB - B));

  while (A < endA) {
    /* Advance B until (A - A_slop_neg) */
    while (B < endB && (*B) + A_slop_neg < (*A)) {
      B++;
    }

    B_save = B;
    debug(printf("Advancing B from %u up to %u\n",*B,(*A) + A_slop_pos));
    while (B < endB && (*B) <= (*A) + A_slop_pos) {
      debug(printf("(1) Generating indices for %u %u\n",(*A),(*B)));
      if (indices >= (*end_indices) /* same as indices + 2 > (*end_indices) */) {
	allocation = (*end_indices) - (*init_indices);
	more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
	
	nindices = indices - (*init_indices);
	memcpy(more_indices,*init_indices,nindices * sizeof(int));

	FREE(*init_indices);
	*init_indices = more_indices;
	*end_indices = &(more_indices[2 * allocation]);
	indices = &(more_indices[nindices]);
      }

      if (A_first_p == true) {
	*indices++ = A_index + (A - initA);
	*indices++ = B_index + (B - initB);
      } else {
	*indices++ = B_index + (B - initB);
	*indices++ = A_index + (A - initA);
      }
      assert(indices <= *end_indices);
	  
      B++;
    }
    B = B_save;

    A++;
  }

  return indices;
}


/* Performs scalar against a single value of A */
static inline int *
scalar_one_set (int *indices, int **indices_for_rare, int **init_indices, int **end_indices,
		const int A_index, const int B_index,
		const unsigned int *A, const unsigned int *B, const unsigned int *endB,
		const unsigned int A_slop_neg, const unsigned int A_slop_pos,
		bool A_first_p) {

  const unsigned int *initA = A;
  const unsigned int *initB = B;
  int allocation, nindices;
  int *more_indices;
  int nindices_for_rare;

  debug(printf("Entered scalar_one_set with one entry for A and %lu entries for B\n",endB - B));

  debug(printf("Advancing B from %u up to %u\n",*B,(*A) - A_slop_neg));
  while (B < endB && (*B) + A_slop_neg < (*A)) {
    B++;
  }

  while (B < endB && (*B) <= (*A) + A_slop_pos) {
    debug(printf("(2) Generating indices for %u %u: ",(*A),(*B)));

    if (indices >= (*end_indices) /* same as indices + 2 > (*end_indices) */) {
      allocation = (*end_indices) - (*init_indices);
      more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
	
      nindices = indices - (*init_indices);
      nindices_for_rare = (*indices_for_rare) - (*init_indices);
      memcpy(more_indices,*init_indices,nindices * sizeof(int));

      FREE(*init_indices);
      *init_indices = more_indices;
      *end_indices = &(more_indices[2 * allocation]);
      indices = &(more_indices[nindices]);
      *indices_for_rare = &(more_indices[nindices_for_rare]);
    }

    debug(printf("%d %d\n",A_index + (A - initA),B_index + (B - initB)));
    if (A_first_p == true) {
      *indices++ = A_index + (A - initA);
      *indices++ = B_index + (B - initB);
    } else {
      *indices++ = B_index + (B - initB);
      *indices++ = A_index + (A - initA);
    }
    assert(indices <= (*end_indices));

    B++;
  }

  return indices;
}


#ifdef __GNUC__
#define COMPILER_LIKELY(x) __builtin_expect((x), 1)
#define COMPILER_RARELY(x) __builtin_expect((x), 0)
#else
#define COMPILER_LIKELY(x) x
#define COMPILER_RARELY(x) x
#endif

/**
 * Intersections scheme designed by N. Kurz that works very
 * well when intersecting an array with another where the density
 * differential is small (between 2 to 10).
 *
 * It assumes that lenRare <= lenFreq.
 *
 * Note that this is not symmetric: flipping the rare and freq pointers
 * as well as lenRare and lenFreq could lead to significant performance
 * differences.
 *
 * The matchOut pointer can safely be equal to the rare pointer.
 *
 */

#ifdef HAVE_SSE2
static int
binary_search_uint4 (int lowi, int highi, const unsigned int *positions, unsigned int goal) {
  int middlei;

  debug20(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug20(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug20(printf("binary search returns %d\n",highi));
  return highi;
}
#endif


#ifdef HAVE_SSE2
static int *
v1 (int *indices, int **init_indices, int **end_indices,
    const int rare_index, const int freq_index,
    const unsigned int *rare, int lenRare,
    const unsigned int *freq, int lenFreq,
    const unsigned int rare_slop_neg, const unsigned int rare_slop_pos,
    bool rare_first_p) {

  int allocation, nindices;
  int *more_indices;
  int *indices_for_rare;
  int npairs, pairi, nindices_for_rare;
  /* int best_pairi, k; */
  /* unsigned int min_distance, distance, freq_value; */

  const unsigned int *initRare = rare;
  const unsigned int *initFreq = freq;
  const unsigned int *endRare = &(rare[lenRare]);
  const unsigned int *endFreq = &(freq[lenFreq]);

  const int kRareSpace = 0;	/* 0. type was uint64_t */

  unsigned int valRare;
  unsigned int goal;
  int npositions_rare, npositions_freq, j;

  int i;
  unsigned int rare_slop_neg_plus_1, rare_slop_pos_plus_1; /* Because only cmpgt is available, not cmpge */

  rare_slop_neg_plus_1 = rare_slop_neg + 1; 
  rare_slop_pos_plus_1 = rare_slop_pos + 1;

#ifdef HAVE_AVX512
  const int kFreqSpace = 2 * 8 * (0 + 1) - 1; /* 15. type was uint64_t */
  const int kProbe = (0 + 1) * 2 * 8; /* 16. type was uint64_t */

  __m512i Rare_lowbound, Rare_highbound;
  __m512i F, F_save;
  __m512i _epi32_offset;
  __mmask16 idx, M_above, M_below;
#elif defined(HAVE_AVX2)
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7. type was uint64_t */
  const int kProbe = (0 + 1) * 2 * 4; /* 8. type was uint64_t */

  __m256i Rare_lowbound, Rare_highbound;
  __m256i F, F_save;
  __m256i _epi32_offset;
  __m256i M, M_above, M_below;
  int idx;
#else
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7. type was uint64_t */
  const int kProbe = (0 + 1) * 2 * 4; /* 8. type was uint64_t */

  __m128i Rare_lowbound, Rare_highbound;
  __m128i F0, F1, F0_save, F1_save;
  __m128i _epi32_offset;
  __m128i M, M_above, M_below;
  int idx;
#endif

  const unsigned int *stopRare = endRare - kRareSpace;
  const unsigned int *stopFreq = endFreq - kFreqSpace;

#ifdef HAVE_AVX512
  _epi32_offset = _mm512_set1_epi32(EPI32_MAX);
#elif defined(HAVE_AVX2)
  _epi32_offset = _mm256_set1_epi32(EPI32_MAX);
#else
  _epi32_offset = _mm_set1_epi32(EPI32_MAX);
#endif


  assert(/*lenRare*/(endRare - initRare) <= /*lenFreq*/(endFreq - initFreq));

  /* valRare = (*rare); -- For comparison, now computed inside condition and in main loop */

  /* Avoid issues where valRare - rare_slop_neg_plus_1 < 0 */
  while (COMPILER_RARELY(rare < stopRare &&
			 (valRare = (*rare)) <= rare_slop_neg)) {
    indices_for_rare = indices;
    indices = scalar_one_set(indices, &indices_for_rare, &(*init_indices), &(*end_indices),
			     rare_index + (rare - initRare), freq_index /*+ (freq - initFreq)*/,
			     rare, freq, endFreq,
			     rare_slop_neg, rare_slop_pos, rare_first_p);
    rare++;
    /* valRare = (*rare); -- For comparison, now computed inside condition and in main loop */
  }

  if (COMPILER_LIKELY(rare < stopRare)) {
    valRare = (*rare);

    /* Advance freq */
    npositions_rare = endRare - rare;
    npositions_freq = endFreq - freq;

    if (50 * npositions_rare > npositions_freq) {
      /* Use linear */
      while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + rare_slop_neg < valRare) {
	debug(printf("Advancing freq by %d because %u + slop < %u\n",
		     kProbe,freq[kProbe - 1],valRare));
	freq += kProbe;
      }

    } else {
      /* Use galloping search */
      goal = valRare - rare_slop_neg;
      j = 1;
      while (j < npositions_freq && freq[j] < goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions_freq) {
	j = binary_search_uint4(j >> 1,npositions_freq,freq,goal);
      } else {
	j = binary_search_uint4(j >> 1,j,freq,goal);
      }
      freq += j;
    }
  }

  if (freq < stopFreq) {
#ifdef HAVE_AVX512
    Rare_lowbound = _mm512_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
    Rare_highbound = _mm512_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);

    F = _mm512_loadu_si512((const __m512i *)(freq));
    F = _mm512_sub_epi32(F, _epi32_offset);
#elif defined(HAVE_AVX2)
    Rare_lowbound = _mm256_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
    Rare_highbound = _mm256_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);

    F = _mm256_loadu_si256((const __m256i *)(freq));
    F = _mm256_sub_epi32(F, _epi32_offset);
#else
    Rare_lowbound = _mm_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
    Rare_highbound = _mm_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);

    F0 = LOAD_SI128((const __m128i *)(freq));
    F0 = _mm_sub_epi32(F0, _epi32_offset);
    F1 = LOAD_SI128((const __m128i *)(freq + 4));
    F1 = _mm_sub_epi32(F1, _epi32_offset);
#endif
  }

  while (rare < stopRare && freq < stopFreq) {
    valRare = (*rare);

    debug(printf("rare is at item %lu with value %u\n",rare - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",freq - initFreq,freq[0],freq[kProbe - 1]));

    const unsigned int *freq_save = freq;
    indices_for_rare = indices;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi32_mask(F, Rare_lowbound); /* 16 comparisons */
    M_below = _mm512_cmplt_epi32_mask(F, Rare_highbound);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F_save = F;

    M_above = _mm256_cmpgt_epi32(F, Rare_lowbound); /* 8 comparisons */
    M_below = _mm256_cmpgt_epi32(Rare_highbound, F); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);

    idx = _mm256_movemask_ps((__m256) M);
#else
    F0_save = F0;
    F1_save = F1;

    M_above = _mm_cmpgt_epi32(F1, Rare_lowbound); /* 4 comparisons */
    M_below = _mm_cmplt_epi32(F1, Rare_highbound);
    M = _mm_and_si128(M_above,M_below);
    idx = _mm_movemask_ps((__m128) M);

    M_above = _mm_cmpgt_epi32(F0, Rare_lowbound); /* 4 comparisons */
    M_below = _mm_cmplt_epi32(F0, Rare_highbound);
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 4) + _mm_movemask_ps((__m128) M);
#endif

    debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

    if (indices + 2*match_n[idx] > (*end_indices)) {
      allocation = (*end_indices) - (*init_indices);
      more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
      
      nindices = indices - (*init_indices);
      nindices_for_rare = indices_for_rare - (*init_indices);

      memcpy(more_indices,*init_indices,nindices * sizeof(int));
      
      FREE(*init_indices);
      *init_indices = more_indices;
      *end_indices = &(more_indices[2 * allocation]);
      indices = &(more_indices[nindices]);
      indices_for_rare = &(more_indices[nindices_for_rare]);
      assert(indices + 2*match_n[idx] <= (*end_indices));
    }

    const unsigned int *freq_ptr = freq + match_start[idx];
    for (i = 0; i < match_n[idx]; i++) {
      debug(printf("(3) Generating output rare %u, freq %u\n",valRare,*freq_ptr));

      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq_ptr++) - initFreq;
      } else {
	*indices++ = freq_index + (freq_ptr++) - initFreq;
	*indices++ = rare_index + (rare - initRare);
      }
    }

    /* Subsequent iterations on freq */
    while (freq + kProbe < stopFreq && /*minFreq*/freq[kProbe] <= valRare + rare_slop_pos) {
      freq += kProbe;
      
#ifdef HAVE_AVX512
      F = _mm512_loadu_si512((const __m512i *)(freq));
      F = _mm512_sub_epi32(F, _epi32_offset);
    
      M_above = _mm512_cmpgt_epi32_mask(F, Rare_lowbound); /* 16 comparisons */
      M_below = _mm512_cmplt_epi32_mask(F, Rare_highbound);
      idx = M_above & M_below;
#elif defined(HAVE_AVX2)
      F = _mm256_loadu_si256((const __m256i *)(freq));
      F = _mm256_sub_epi32(F, _epi32_offset);
    
      M_above = _mm256_cmpgt_epi32(F, Rare_lowbound); /* 8 comparisons */
      M_below = _mm256_cmpgt_epi32(Rare_highbound, F); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);

      idx = _mm256_movemask_ps((__m256) M);
#else
      F0 = LOAD_SI128((const __m128i *)(freq));
      F0 = _mm_sub_epi32(F0, _epi32_offset);
      F1 = LOAD_SI128((const __m128i *)(freq + 4));
      F1 = _mm_sub_epi32(F1, _epi32_offset);

      M_above = _mm_cmpgt_epi32(F1, Rare_lowbound); /* 4 comparisons */
      M_below = _mm_cmplt_epi32(F1, Rare_highbound);
      M = _mm_and_si128(M_above,M_below);
      idx = _mm_movemask_ps((__m128) M);

      M_above = _mm_cmpgt_epi32(F0, Rare_lowbound); /* 4 comparisons */
      M_below = _mm_cmplt_epi32(F0, Rare_highbound);
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 4) + _mm_movemask_ps((__m128) M);
#endif

      debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

      if (indices + 2*match_n[idx] > (*end_indices)) {
	allocation = (*end_indices) - (*init_indices);
	more_indices = (int *) MALLOC(2 * allocation * sizeof(int));
      
	nindices = indices - (*init_indices);
	nindices_for_rare = indices_for_rare - (*init_indices);

	memcpy(more_indices,*init_indices,nindices * sizeof(int));
      
	FREE(*init_indices);
	*init_indices = more_indices;
	*end_indices = &(more_indices[2 * allocation]);
	indices = &(more_indices[nindices]);
	indices_for_rare = &(more_indices[nindices_for_rare]);
	assert(indices + 2*match_n[idx] <= (*end_indices));
      }

      const unsigned int *freq_ptr = freq + match_start[idx];
      for (i = 0; i < match_n[idx]; i++) {
	debug(printf("(4) Generating output rare %u, freq %u\n",valRare,*freq_ptr));

	if (rare_first_p == true) {
	  *indices++ = rare_index + (rare - initRare);
	  *indices++ = freq_index + (freq_ptr++) - initFreq;
	} else {
	  *indices++ = freq_index + (freq_ptr++) - initFreq;
	  *indices++ = rare_index + (rare - initRare);
	}
      }
    }

    if (freq + kProbe >= stopFreq) {
      freq += kProbe;
      indices = scalar_one_set(indices, &indices_for_rare, &(*init_indices), &(*end_indices),
			       rare_index + (rare - initRare), freq_index + (freq - initFreq),
			       rare, freq, endFreq,
			       rare_slop_neg, rare_slop_pos, rare_first_p);
    }
    /* End of iterations on freq */

#ifdef PICK_BEST_DISTANCE
    debug2(printf("rare %u => %d pairs of indices\n",*rare,(int) (indices - indices_for_rare)/2));
    if ((npairs = (indices - indices_for_rare)/2) > MAX_HITS_PER_RARE) {
      if (rare_slop_neg == 0) {
	/* Take the first MAX_HITS */
	indices = &(indices_for_rare[2*MAX_HITS_PER_RARE]);
      } else if (rare_slop_pos == 0) {
	/* Take the last MAX_HITS */
	indices = indices_for_rare;
	for (pairi = npairs - MAX_HITS_PER_RARE; pairi < npairs; pairi++) {
	  *indices++ = indices_for_rare[2*pairi];
	  *indices++ = indices_for_rare[2*pairi + 1];
	}
      } else {
	fprintf(stderr,"Expecting one slop to be 0\n");
	abort();
      }
    }
#endif

    if (++rare < stopRare) {
      valRare = (*rare); /* For comparison */

#ifdef HAVE_AVX512
      Rare_lowbound = _mm512_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
      Rare_highbound = _mm512_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_lowbound = _mm256_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
      Rare_highbound = _mm256_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);
      F = F_save;
#else
      Rare_lowbound = _mm_set1_epi32(valRare - rare_slop_neg_plus_1 - EPI32_MAX);
      Rare_highbound = _mm_set1_epi32(valRare + rare_slop_pos_plus_1 - EPI32_MAX);
      F0 = F0_save;
      F1 = F1_save;
#endif

      /* Advance freq */
      freq = freq_save;

      npositions_rare = endRare - rare;
      npositions_freq = endFreq - freq;

      if (50 * npositions_rare > npositions_freq) {
	/* Use linear */
	while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + rare_slop_neg < valRare) {
	  debug(printf("Advancing freq by %d because %u + slop < %u\n",kProbe,freq[kProbe - 1],valRare));
	  freq += kProbe;
	}

      } else {
	/* Use galloping search */
	goal = valRare - rare_slop_neg;
	j = 1;
	while (j < npositions_freq && freq[j] < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions_freq) {
	  j = binary_search_uint4(j >> 1,npositions_freq,freq,goal);
	} else {
	  j = binary_search_uint4(j >> 1,j,freq,goal);
	}
	freq += j;
      }

      if (freq < stopFreq && freq != freq_save) {
#ifdef HAVE_AVX512
	F = _mm512_loadu_si512((const __m512i *)(freq));
	F = _mm512_sub_epi32(F, _epi32_offset);
#elif defined(HAVE_AVX2)
	F = _mm256_loadu_si256((const __m256i *)(freq));
	F = _mm256_sub_epi32(F, _epi32_offset);
#else
	F0 = LOAD_SI128((const __m128i *)(freq));
	F0 = _mm_sub_epi32(F0, _epi32_offset);
	F1 = LOAD_SI128((const __m128i *)(freq + 4));
	F1 = _mm_sub_epi32(F1, _epi32_offset);
#endif
      }

      debug(printf("Advancing rare to %lu to compare with stop %lu\n",
		   rare - initRare,stopRare - initRare));
      debug(printf("Restoring freq to %lu\n",freq - initFreq));
    }
  }

  return scalar(indices, &(*init_indices), &(*end_indices),
		rare_index + (rare - initRare), freq_index + (freq - initFreq),
		rare, endRare, freq, endFreq,
		rare_slop_neg, rare_slop_pos, rare_first_p);
}
#endif


/**
 * Our main heuristic.
 *
 * The out pointer can be set1 if length1<=length2,
 * or else it can be set2 if length1>length2.
 */
#ifdef HAVE_SSE2
int
Intersect_concordance (int **result,
		       const unsigned int *set1, const int length1,
		       const unsigned int *set2, const int length2,
		       const unsigned int slop1_neg, const unsigned int slop1_pos,
		       const unsigned int slop2_neg, const unsigned int slop2_pos) {
  int *indices, *init_indices, *end_indices;
  int allocation, nindices;

  debug(printf("Entered Intersect_concordance with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%d: %u\n",i,set1[i]);
  }
  printf("\n");
  for (i = 0; i < length2; i++) {
    printf("%d: %u\n",i,set2[i]);
  }
  printf("\n");
#endif    

  if ((length1 == 0) || (length2 == 0)) {
    debug(printf("Exiting Intersect_concordance with 0 items\n"));
    *result = (int *) NULL;
    return 0;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      if ((allocation = 2 * length1) < 2*NELTS) {
	/* Needs to be greater than 2*match_n */
	allocation = 2*NELTS;
      }
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);

      indices = v1(indices, &init_indices, &end_indices,
		   /*rare_index*/0, /*freq_index*/0,
		   set1, length1, set2, length2,
		   slop1_neg, slop1_pos, /*rare_first_p*/true);
    } else {
      if ((allocation = 2 * length2) < 2*NELTS) {
	/* Needs to be greater than 2*match_n */
	allocation = 2*NELTS;
      }
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);

      indices = v1(indices, &init_indices, &end_indices,
		   /*rare_index*/0, /*freq_index*/0,
		   set2, length2, set1, length1,
		   slop2_neg, slop2_pos, /*rare_first_p*/false);
    }
    assert(indices <= end_indices);
  }

  *result = init_indices;
  nindices = indices - init_indices;

#ifdef DEBUG
  int k;
  printf("Returning %d pairs\n",nindices/2);
  for (i = 0, k = 0; i < nindices / 2; i++, k += 2) {
    printf("%d %d\n",(*result)[k],(*result)[k+1]);
  }
#endif

  return (nindices / 2);
}

#else

int
Intersect_concordance (int **result,
		       const unsigned int *set1, const int length1,
		       const unsigned int *set2, const int length2,
		       const unsigned int slop1_neg, const unsigned int slop1_pos,
		       const unsigned int slop2_neg, const unsigned int slop2_pos) {

  int *indices, *init_indices, *end_indices;
  int allocation, nindices;

  debug(printf("Entered Intersect_concordance with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%d: %u\n",i,set1[i]);
  }
  printf("\n");
  for (i = 0; i < length2; i++) {
    printf("%d: %u\n",i,set2[i]);
  }
  printf("\n");
#endif    

  if ((length1 == 0) || (length2 == 0)) {
    *result = (int *) NULL;
    return 0;

  } else {
    if (length1 <= length2) {
      allocation = 2 * length1;
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);

      indices = scalar(indices, &init_indices, &end_indices,
		       /*rare_index*/0, /*freq_index*/0,
		       set1, /*endA*/&(set1[length1]), set2, /*endB*/&(set2[length2]),
		       slop1_neg, slop1_pos, /*A_first_p*/true);
    } else {
      allocation = 2 * length2;
      indices = init_indices = (int *) MALLOC(allocation * sizeof(int));
      end_indices = &(indices[allocation]);

      indices = scalar(indices, &init_indices, &end_indices,
		       /*rare_index*/0, /*freq_index*/0,
		       set2, /*endA*/&(set2[length2]), set1, /*endB*/&(set1[length1]),
		       slop2_neg, slop2_pos, /*A_first_p*/false);
    }
  }

  *result = init_indices;
  nindices = indices - init_indices;

#ifdef DEBUG
  int k;
  printf("Returning %d pairs\n",nindices/2);
  for (i = 0, k = 0; i < nindices / 2; i++, k += 2) {
    printf("%d %d\n",(*result)[k],(*result)[k+1]);
  }
#endif

  return (nindices / 2);
}
#endif


void
Intersect_concordance_setup () {
  initialize_match_values();
  return;
}  
