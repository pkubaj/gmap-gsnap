static char rcsid[] = "$Id: 3a30aadcabfb30a9b76cfbbefc2f5e6c250f1e0c $";
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

/* Minimum requirement for SIMD is SSE4.1 for _mm_cvtepi8_epi64 */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "genomebits.h"
#include "intersect-indices2-large.h"

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

/* binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


#define NELTS 8			/* 8 64-mers in 512 bits */
#define MAX_IDX 256		/* 2^8 */

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


#define GETPOS(high,low) (((UINT8) high << 32) + low)
#define GETPOS_I(high,low,i) (((UINT8) high[i] << 32) + low[i])


#ifndef HAVE_SSE4_1
static int
__frogadvanceUntil_large(const unsigned char *array_high, const UINT4 *array_low, const int pos,
			 const int length, const Univcoord_T min) {
  int lower = pos + 1;

  /* special handling for a possibly common sequential case */
  if ((lower >= length) || (GETPOS_I(array_high,array_low,lower) >= min)) {
    return lower;
  }

  int spansize = 1; /* could set larger */
  /* bootstrap an upper limit */

  while ((lower + spansize < length) && (GETPOS_I(array_high,array_low,lower + spansize) < min)) {
    spansize *= 2;
  }
  int upper = (lower + spansize < length) ? lower + spansize : length - 1;

  if (GETPOS_I(array_high,array_low,upper) < min) { /* means array has no item >= min */
    return length;
  }

  /* we know that the next-smallest span was too small */
  lower += (spansize / 2);

  /* else begin binary search */
  int mid = 0;
  while (lower + 1 != upper) {
    /* mid = (lower + upper) / 2; */
    mid = lower + (upper - lower) / 2;
    if (GETPOS_I(array_high,array_low,mid) == min) {
      return mid;
    } else if (GETPOS_I(array_high,array_low,mid) < min) {
      lower = mid;
    } else {
      upper = mid;
    }
  }

  return upper;
}


/* smallset_first_p is false */
static int
galloping_intersection_anchorA (int *indices, const int smallset_index, const int largeset_index,
				const Univcoord_T *smallset, const int smalllength,
				const unsigned char *largeset_high, const UINT4 *largeset_low, const int largelength,
				int delta_small, int diagterm_small) {
  if (smalllength == 0) return 0;
  const int *initout = indices;
  int k1 = 0, k2 = 0;

  while (1) {
    if (GETPOS_I(largeset_high,largeset_low,k1) < smallset[k2] + delta_small) {
      k1 = __frogadvanceUntil_large(largeset_high, largeset_low, k1, largelength, smallset[k2] + delta_small);
      if (k1 == largelength) {
        break;
      }
    }
  midpoint:
    if (smallset[k2] + delta_small < GETPOS_I(largeset_high,largeset_low,k1)) {
      ++k2;
      if (k2 == smalllength) {
        break;
      }
    } else {
      *indices++ = largeset_index + k1;
      *indices++ = smallset_index + k2;
      ++k2;
      if (k2 == smalllength) {
        break;
      }
      k1 = __frogadvanceUntil_large(largeset_high, largeset_low, k1, largelength, smallset[k2] + delta_small);
      if (k1 == largelength) {
        break;
      }
      goto midpoint;
    }
  }

  return indices - initout;
}


static int
__frogadvanceUntil_univcoord(const Univcoord_T *array, const int pos,
			     const int length, const Univcoord_T min) {
  int lower = pos + 1;

  /* special handling for a possibly common sequential case */
  if ((lower >= length) || (array[lower] >= min)) {
    return lower;
  }

  int spansize = 1; /* could set larger */
  /* bootstrap an upper limit */

  while ((lower + spansize < length) && (array[lower + spansize] < min)) {
    spansize *= 2;
  }
  int upper = (lower + spansize < length) ? lower + spansize : length - 1;

  if (array[upper] < min) { /* means array has no item >= min */
    return length;
  }

  /* we know that the next-smallest span was too small */
  lower += (spansize / 2);

  /* else begin binary search */
  int mid = 0;
  while (lower + 1 != upper) {
    /* mid = (lower + upper) / 2; */
    mid = lower + (upper - lower) / 2;
    if (array[mid] == min) {
      return mid;
    } else if (array[mid] < min) {
      lower = mid;
    } else {
      upper = mid;
    }
  }

  return upper;
}


/* smallset_first_p is true */
static int
galloping_intersection_anchorB (int *indices, const int smallset_index, const int largeset_index,
				const unsigned char *smallset_high, const UINT4 *smallset_low, const int smalllength,
				const Univcoord_T *largeset, const int largelength,
				int delta_small, int diagterm_small) {
  if (smalllength == 0) return 0;
  const int *initout = indices;
  int k1 = 0, k2 = 0;

  while (1) {
    if (largeset[k1] < GETPOS_I(smallset_high,smallset_low,k2) + delta_small) {
      k1 = __frogadvanceUntil_univcoord(largeset, k1, largelength, GETPOS_I(smallset_high,smallset_low,k2) + delta_small);
      if (k1 == largelength) {
        break;
      }
    }
  midpoint:
    if (GETPOS_I(smallset_high,smallset_low,k2) + delta_small < largeset[k1]) {
      ++k2;
      if (k2 == smalllength) {
        break;
      }
    } else {
      *indices++ = smallset_index + k2;
      *indices++ = largeset_index + k1;
      ++k2;
      if (k2 == smalllength) {
        break;
      }
      k1 = __frogadvanceUntil_univcoord(largeset, k1, largelength, GETPOS_I(smallset_high,smallset_low,k2) + delta_small);
      if (k1 == largelength) {
        break;
      }
      goto midpoint;
    }
  }

  return indices - initout;
}
#endif



/* Fast scalar scheme designed by N. Kurz. */
/* Delta is required to be added to A to become equivalent with B */
/* Desired result should be A + diagtermA (== B + diagtermB) */
/* A_first_p is false */
static int
scalar_anchorA (int *indices, const int A_index, const int B_index,
		const Univcoord_T *A, const int lenA,
		const unsigned char *B_high, const UINT4 *B_low, const int lenB,
		int deltaA, int diagtermA) {
  const int *initout = indices;
  if (lenA == 0 || lenB == 0) return 0;

  const Univcoord_T *initA = A;
  const UINT4 *initB = B_low;
  const Univcoord_T *endA = A + lenA;
  const UINT4 *endB = B_low + lenB;

  while (1) {
    while ((*A) + deltaA < GETPOS(*B_high,*B_low)) {
    SKIP_FIRST_COMPARE:
      if (++A == endA) {
        return (indices - initout);
      }
    }
    while ((*A) + deltaA > GETPOS(*B_high,*B_low)) {
      ++B_high;
      if (++B_low == endB) {
        return (indices - initout);
      }
    }
    if ((*A) + deltaA == GETPOS(*B_high,*B_low)) {
      /* *indices++ = B_index + (B_low - initB); -- want only indices2 */
      *indices++ = A_index + (A - initA);
      ++B_high;
      if (++A == endA || ++B_low == endB) {
        return (indices - initout);
      }
    } else {
      goto SKIP_FIRST_COMPARE;
    }
  }

  return (indices - initout); 	/* Not reached */
}


/* A_first_p is true */
static int
scalar_anchorB (int *indices, const int A_index, const int B_index,
		const unsigned char *A_high, const UINT4 *A_low, const int lenA,
		const Univcoord_T *B, const int lenB, int deltaA, int diagtermA) {
  const int *initout = indices;
  if (lenA == 0 || lenB == 0) return 0;

  const UINT4 *initA = A_low;
  const Univcoord_T *initB = B;
  const UINT4 *endA = A_low + lenA;
  const Univcoord_T *endB = B + lenB;

  while (1) {
    while (GETPOS(*A_high,*A_low) + deltaA < *B) {
    SKIP_FIRST_COMPARE:
      ++A_high;
      if (++A_low == endA) {
        return (indices - initout);
      }
    }
    while (GETPOS(*A_high,*A_low) + deltaA > *B) {
      if (++B == endB) {
        return (indices - initout);
      }
    }
    if (GETPOS(*A_high,*A_low) + deltaA == *B) {
      /* *indices++ = A_index + (A_low - initA); -- want only indices2 */
      *indices++ = B_index + (B - initB);
      ++A_high;
      if (++A_low == endA || ++B == endB) {
        return (indices - initout);
      }
    } else {
      goto SKIP_FIRST_COMPARE;
    }
  }

  return (indices - initout); 	/* Not reached */
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

#if defined(HAVE_SSE4_1)
/* rare_first_p is false */
static int
v1_anchorA (int *indices, const int rare_index, const int freq_index,
	    const Univcoord_T *rare, int lenRare,
	    const unsigned char *freq_high, const UINT4 *freq_low, int lenFreq,
	    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  Univcoord_T valRare;
  Univcoord_T maxProbe;		/* was uint64_t */
  Univcoord_T maxFreq;		/* was uint64_t */
  const int kRareSpace = 0;	/* was uint64_t */

  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 8 64-mers in 512 bits */

#ifdef HAVE_AVX512
  __m512i Rare;
  __m128i H;
  __m512i F;
  __mmask8 idx;
#elif defined(HAVE_AVX2)
  __m256i Rare;
  __m128i H;
  __m256i F0, F1;
  __m256i M;
  unsigned int idx;
#else
  __m128i Rare;
  __m128i H;
  __m128i F0, F1, F2, F3;
  __m128i M;
  unsigned int idx;
#endif
  int bit;

  const Univcoord_T *stopRare = &rare[lenRare] - kRareSpace;
  const UINT4 *stopFreq = &freq_low[lenFreq] - kFreqSpace;


  if (COMPILER_RARELY((rare >= stopRare) || (freq_low >= stopFreq)))
    goto FINISH_SCALAR;

  valRare = (*rare) + delta_rare; /* For comparison */
  /* outRare = (*rare) + diagterm_rare; -- For output */

#ifdef HAVE_AVX512
  Rare = _mm512_set1_epi64(valRare);
#elif defined(HAVE_AVX2)
  Rare = _mm256_set1_epi64x(valRare);
#else
  Rare = _mm_set1_epi64((__m64) valRare);
#endif

  maxFreq = GETPOS_I(freq_high,freq_low,2 * 4 - 1);
#ifdef HAVE_AVX512
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
		       _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
#elif defined(HAVE_AVX2)
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,4);
  F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
#else
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,2);
  F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 2))));
  H = _mm_bslli_si128(H,2);
  F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
  H = _mm_bslli_si128(H,2);
  F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 6))));
#endif

  if (COMPILER_RARELY(maxFreq < valRare))
    goto ADVANCE_FREQ;

 ADVANCE_RARE:
  do {
    rare += 1;
    if (COMPILER_RARELY(rare >= stopRare)) {
      rare -= 1;
      goto FINISH_SCALAR;
    }

    valRare = (*rare) + delta_rare; /* for next iteration */
    /* outRare = (*rare) + diagterm_rare; -- for next iteration */

#ifdef HAVE_AVX512
    idx = _mm512_cmpeq_epi64_mask(F,Rare);

    Rare = _mm512_set1_epi64(valRare);
#elif defined(HAVE_AVX2)
    M = _mm256_cmpeq_epi64(F1, Rare);
    idx = _mm256_movemask_pd((__m256d) M);
    M = _mm256_cmpeq_epi64(F0, Rare);
    idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);

    Rare = _mm256_set1_epi64x(valRare);
#else
    M = _mm_cmpeq_epi64(F3, Rare);
    idx = _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F2, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F1, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F0, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    Rare = _mm_set1_epi64((__m64) valRare);
#endif

    if ((bit = match_start[idx]) >= 0) {
      assert(match_n[idx] == 1);
      /* *indices++ = freq_index + (freq_low + bit - initFreq); -- want only indices2 */
      *indices++ = rare_index + ((rare - 1) - initRare); /* rare was advanced by 1 */
    }

  } while (maxFreq >= valRare);

 ADVANCE_FREQ:
  do {
    const int kProbe = (0 + 1) * 2 * 4; /* 8 64-mers in 512 bits */
    const unsigned char *probeFreq_high = freq_high + kProbe;
    const UINT4 *probeFreq_low = freq_low + kProbe;

    if (COMPILER_RARELY(probeFreq_low >= stopFreq)) {
      goto FINISH_SCALAR;
    }
    maxProbe = GETPOS_I(freq_high,freq_low,(0 + 2) * 2 * 4 - 1);

    freq_high = probeFreq_high;
    freq_low = probeFreq_low;

  } while (maxProbe < valRare);

  maxFreq = maxProbe;

#ifdef HAVE_AVX512
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
		       _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
#elif defined(HAVE_AVX2)
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,4);
  F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
#else
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,2);
  F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 2))));
  H = _mm_bslli_si128(H,2);
  F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
  H = _mm_bslli_si128(H,2);
  F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 6))));
#endif

  goto ADVANCE_RARE;

  int count;
 FINISH_SCALAR:
  count = indices - initout;

  lenFreq = stopFreq + kFreqSpace - freq_low;
  lenRare = stopRare + kRareSpace - rare;

  int tail = scalar_anchorA(indices, rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
			    rare, lenRare, freq_high, freq_low, lenFreq,
			    delta_rare, diagterm_rare);

  return count + tail;
}


/* rare_first_p is true */
static int
v1_anchorB (int *indices, const int rare_index, const int freq_index,
	    const unsigned char *rare_high, const UINT4 *rare_low, int lenRare,
	    const Univcoord_T *freq, int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  Univcoord_T valRare;
  Univcoord_T maxProbe;		/* was uint64_t */
  Univcoord_T maxFreq;		/* was uint64_t */
  const int kRareSpace = 0;	/* was uint64_t */

  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 8 64-mers in 512 bits */

#ifdef HAVE_AVX512
  __m512i Rare;
  __m512i F;
  __mmask8 idx;
#elif defined(HAVE_AVX2)
  __m256i Rare;
  __m256i F0, F1;
  __m256i M;
  unsigned int idx;
#else
  __m128i Rare;
  __m128i F0, F1, F2, F3;
  __m128i M;
  unsigned int idx;
#endif
  int bit;

  const UINT4 *stopRare = &rare_low[lenRare] - kRareSpace;
  const Univcoord_T *stopFreq = &freq[lenFreq] - kFreqSpace;


  if (COMPILER_RARELY((rare_low >= stopRare) || (freq >= stopFreq)))
    goto FINISH_SCALAR;

  valRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* For comparison */
  /* outRare = (*rare) + diagterm_rare; -- For output */

#ifdef HAVE_AVX512
  Rare = _mm512_set1_epi64(valRare);
#elif defined(HAVE_AVX2)
  Rare = _mm256_set1_epi64x(valRare);
#else
  Rare = _mm_set1_epi64((__m64) valRare);
#endif

  maxFreq = freq[2 * 4 - 1];
#ifdef HAVE_AVX512
  F = _mm512_loadu_si512((const __m512i *)(freq));
#elif defined(HAVE_AVX2)
  F0 = _mm256_loadu_si256((const __m256i *)(freq));
  F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
#else
  maxFreq = freq[2 * 4 - 1];
  F0 = _mm_lddqu_si128((const __m128i *)(freq));
  F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
  F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
  F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
#endif

  if (COMPILER_RARELY(maxFreq < valRare))
    goto ADVANCE_FREQ;

 ADVANCE_RARE:
  do {
    rare_high += 1;
    rare_low += 1;
    if (COMPILER_RARELY(rare_low >= stopRare)) {
      rare_high -= 1;
      rare_low -= 1;
      goto FINISH_SCALAR;
    }

    valRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* for next iteration */
    /* outRare = (*rare) + diagterm_rare; -- for next iteration */

#ifdef HAVE_AVX512
    idx = _mm512_cmpeq_epi64_mask(F,Rare);

    Rare = _mm512_set1_epi64(valRare);
#elif defined(HAVE_AVX2)
    M = _mm256_cmpeq_epi64(F1, Rare);
    idx = _mm256_movemask_pd((__m256d) M);
    M = _mm256_cmpeq_epi64(F0, Rare);
    idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);

    Rare = _mm256_set1_epi64x(valRare);
#else
    M = _mm_cmpeq_epi64(F3, Rare);
    idx = _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F2, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F1, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi64(F0, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    Rare = _mm_set1_epi64((__m64) valRare);
#endif

    if ((bit = match_start[idx]) >= 0) {
      assert(match_n[idx] == 1);
      /* *indices++ = rare_index + ((rare_low - 1) - initRare); -- rare was advanced by 1; want only indices2 */
      *indices++ = freq_index + (freq + bit - initFreq);
    }

  } while (maxFreq >= valRare);

 ADVANCE_FREQ:
  do {
    const int kProbe = (0 + 1) * 2 * 4; /* 8 64-mers in 512 bits */
    const Univcoord_T *probeFreq = freq + kProbe;

    if (COMPILER_RARELY(probeFreq >= stopFreq)) {
      goto FINISH_SCALAR;
    }
    maxProbe = freq[(0 + 2) * 2 * 4 - 1];

    freq = probeFreq;

  } while (maxProbe < valRare);

  maxFreq = maxProbe;

#ifdef HAVE_AVX512
  F = _mm512_loadu_si512((const __m512i *)(freq));
#elif defined(HAVE_AVX2)
  F0 = _mm256_loadu_si256((const __m256i *)(freq));
  F0 = _mm256_loadu_si256((const __m256i *)(freq + 4));
#else
  F0 = _mm_lddqu_si128((const __m128i *)(freq));
  F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
  F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
  F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
#endif

  goto ADVANCE_RARE;

  int count;
 FINISH_SCALAR:
  count = indices - initout;

  lenFreq = stopFreq + kFreqSpace - freq;
  lenRare = stopRare + kRareSpace - rare_low;

  int tail = scalar_anchorB(indices, rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
			    rare_high, rare_low, lenRare, freq, lenFreq, delta_rare, diagterm_rare);

  return count + tail;
}
#endif


/**
 * This intersection function is similar to v1, but is faster when
 * the difference between lenRare and lenFreq is large, but not too large.

 * It assumes that lenRare <= lenFreq.
 *
 * Note that this is not symmetric: flipping the rare and freq pointers
 * as well as lenRare and lenFreq could lead to significant performance
 * differences.
 *
 * This function DOES NOT use inline assembly instructions. Just intrinsics.
 */

#ifdef HAVE_AVX512
/* Processes 32 vectors, each vector having 8 unsigned ints => 256 unsigned ints */
/* rare_first_p is false */
static int
v3_anchorA (int *indices, const int rare_index, const int freq_index,
	    const Univcoord_T *rare, const int lenRare,
	    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  /* 4 masks * 8 bits/mask = 32 bits */
  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < (*rare) + delta_rare) {
    freq_high += veclen * 32;
    freq_low += veclen * 32;
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

    const vec Match = _mm512_set1_epi64(matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1, F2, F3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	if (GETPOS_I(freq_high,freq_low,veclen * 3 + vecmax) >= matchRare) {
	  base = 0 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3)));

	} else {
	  base = 1 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5)));
	
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7)));
	}

      } else {
	if (GETPOS_I(freq_high,freq_low,veclen * 11 + vecmax) >= matchRare) {
	  base = 2 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11)));

	} else {
	  base = 3 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15)));
	}
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	if (GETPOS_I(freq_high,freq_low,veclen * 19 + vecmax) >= matchRare) {
	  base = 4 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 8*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1 + 16)));
	  
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 9*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2 + 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3 + 16)));

	} else {
	  base = 5 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 10*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5 + 16)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 11*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6 + 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7 + 16)));
	}

      } else {
	if (GETPOS_I(freq_high,freq_low,veclen * 27 + vecmax) >= matchRare) {
	  base = 6 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 12*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9 + 16)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 13*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10+ 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11+ 16)));

	} else {
	  base = 7 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 14*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12+ 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13+ 16)));
	  
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 15*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14+ 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15+ 16)));
	}
      }
    }

    idx = _mm512_cmpeq_epi64_mask(F3, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F2, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F1, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F0, Match);

    if (idx == 0) {
      /* Skip */
    } else {
      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices,
		   rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}


/* rare_first_p is true */
static int
v3_anchorB (int *indices, const int rare_index, const int freq_index,
	    const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- nextRare */
    const int indexRare = rare_index + (rare_low - initRare);
    const int indexFreq = freq_index + (freq - initFreq);

    const vec Match = _mm512_set1_epi64(matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	if (freq[veclen * 3 + vecmax] >= matchRare) {
	  base = 0 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 3), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 2), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 1), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 0), Match);
	} else {
	  base = 1 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 7), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 6), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 5), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 4), Match);
	}
      } else {
	if (freq[veclen * 11 + vecmax] >= matchRare) {
	  base = 2 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 11), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 10), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 9), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 8), Match);
	} else {
	  base = 3 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 15), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 14), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 13), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 12), Match);
	}
      }
    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	if (freq[veclen * 19 + vecmax] >= matchRare) {
	  base = 4 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 3 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 2 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 1 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 0 + 16), Match);
	} else {
	  base = 5 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 7 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 6 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 5 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 4 + 16), Match);
	}
      } else {
	if (freq[veclen * 27 + vecmax] >= matchRare) {
	  base = 6 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 11 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 10+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 9 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 8 + 16), Match);
	} else {
	  base = 7 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 15+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 14+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 13+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 12+ 16), Match);
	}
      }
    }

    if (idx == 0) {
      /* Skip */
    } else {
      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}

#elif defined(HAVE_AVX2)
/* Processes 32 vectors, each vector having 4 unsigned ints => 128 unsigned ints */
/* rare_first_p is false */
static int
v3_anchorA (int *indices, const int rare_index, const int freq_index,
	    const Univcoord_T *rare, const int lenRare,
	    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < (*rare) + delta_rare) {
    freq_high += veclen * 32;
    freq_low += veclen * 32;
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

    const vec Match = _mm256_set1_epi64x(matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1, F2, F3, F4, F5, F6, F7;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	base = 0 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));

      } else {
	base = 1 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	base = 2 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));

      } else {
	base = 3 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
      }
    }

    vec Q, Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Q0 = _mm256_cmpeq_epi64(F0, Match);
    Q1 = _mm256_cmpeq_epi64(F1, Match);
    Q2 = _mm256_cmpeq_epi64(F2, Match);
    Q3 = _mm256_cmpeq_epi64(F3, Match);
    Q4 = _mm256_cmpeq_epi64(F4, Match);
    Q5 = _mm256_cmpeq_epi64(F5, Match);
    Q6 = _mm256_cmpeq_epi64(F6, Match);
    Q7 = _mm256_cmpeq_epi64(F7, Match);

    Q = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3)),
			_mm256_or_si256(_mm256_or_si256(Q4, Q5),_mm256_or_si256(Q6, Q7)));

    if (_mm256_testz_si256(Q, Q)) {
      /* Skip */
    } else {
      idx = _mm256_movemask_pd((__m256d) Q7);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q6);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q5);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q4);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q3);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q2);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q1);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q0);

      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices,
		   rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}


/* rare_first_p is true */
static int
v3_anchorB (int *indices, const int rare_index, const int freq_index,
	    const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- nextRare */

    const vec Match = _mm256_set1_epi64x(matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    vec M, M0, M1, M2, M3, M4, M5, M6, M7;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	base = 0 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 0), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 1), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 2), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 3), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 4), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 5), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 6), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 7), Match);

      } else {
	base = 1 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 8), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 9), Match);
	M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15), Match);
      }

    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	base = 2 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 0 + 16), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 1 + 16), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 2 + 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 3 + 16), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 4 + 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 5 + 16), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 6 + 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 7 + 16), Match);

      } else {
	base = 3 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 8 + 16), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 9 + 16), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10 + 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11 + 16), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12 + 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13 + 16), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14 + 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15 + 16), Match);
      }
    }

    M = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(M0, M1),_mm256_or_si256(M2, M3)),
			_mm256_or_si256(_mm256_or_si256(M4, M5),_mm256_or_si256(M6, M7)));
    if (_mm256_testz_si256(M, M)) {
      /* Skip */
    } else {
      idx = _mm256_movemask_pd((__m256d) M7);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M6);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M5);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M4);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M3);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M2);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M1);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M0);

      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE4_1)
/* Processes 32 vectors, each vector having 2 unsigned ints => 64 unsigned ints */
/* rare_first_p is false */
static int
v3_anchorA (int *indices, const int rare_index, const int freq_index,
	    const Univcoord_T *rare, const int lenRare,
	    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  /* Measuring only 16 bits, but procedures expect 32 bits */
  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < (*rare) + delta_rare) {
    freq_high += veclen * 32;
    freq_low += veclen * 32;
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

    const vec Match = _mm_set1_epi64((__m64) matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
        goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1, F2, F3, F4, F5, F6, F7;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	base = 0 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));

      } else {
	base = 1 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	base = 2 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));

      } else {
	base = 3 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
      }
    }

    vec Q, Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Q0 = _mm_cmpeq_epi64(F0, Match);
    Q1 = _mm_cmpeq_epi64(F1, Match);
    Q2 = _mm_cmpeq_epi64(F2, Match);
    Q3 = _mm_cmpeq_epi64(F3, Match);
    Q4 = _mm_cmpeq_epi64(F4, Match);
    Q5 = _mm_cmpeq_epi64(F5, Match);
    Q6 = _mm_cmpeq_epi64(F6, Match);
    Q7 = _mm_cmpeq_epi64(F7, Match);

    Q = _mm_or_si128(_mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3)),
		     _mm_or_si128(_mm_or_si128(Q4, Q5), _mm_or_si128(Q6, Q7)));

    if (_mm_testz_si128(Q, Q)) {
      /* Skip */
    } else {
      idx = _mm_movemask_pd((__m128d) Q7);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q6);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q5);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q4);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q3);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q2);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q1);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q0);

      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices,
		   rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}


/* rare_first_p is true */
static int
v3_anchorB (int *indices, const int rare_index, const int freq_index,
	    const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  /* Measuring only 16 bits, but procedures expect 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

    const vec Match = _mm_set1_epi64((__m64) matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
        goto FINISH_SCALAR;
      }
    }

    vec M, M0, M1, M2, M3, M4, M5, M6, M7;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	base = 0 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 0), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 1), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 2), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 3), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 4), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 5), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 6), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 7), Match);

      } else {
	base = 1 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 8), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 9), Match);
	M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 10), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 11), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 12), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 13), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 14), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 15), Match);
      }
    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	base = 2 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 0 + 16), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 1 + 16), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 2 + 16), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 3 + 16), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 4 + 16), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 5 + 16), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 6 + 16), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 7 + 16), Match);

      } else {
	base = 3 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 8 + 16), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 9 + 16), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 10 + 16), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 11 + 16), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 12 + 16), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 13 + 16), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 14 + 16), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 15 + 16), Match);
      }
    }

    M = _mm_or_si128(_mm_or_si128(_mm_or_si128(M0, M1),_mm_or_si128(M2, M3)),
			_mm_or_si128(_mm_or_si128(M4, M5),_mm_or_si128(M6, M7)));

    if (_mm_testz_si128(M, M)) {
      /* Skip */
    } else {
      idx = _mm_movemask_pd((__m128d) M7);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M6);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M5);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M4);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M3);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M2);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M1);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M0);

      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}
#endif

/**
 * This is the SIMD galloping function. This intersection function works well
 * when lenRare and lenFreq have vastly different values.
 *
 * It assumes that lenRare <= lenFreq.
 *
 * Note that this is not symmetric: flipping the rare and freq pointers
 * as well as lenRare and lenFreq could lead to significant performance
 * differences.
 *
 * This function DOES NOT use assembly. It only relies on intrinsics.
 */

#ifdef HAVE_AVX512
/* rare_first_p is false */
static int
SIMDgalloping_anchorA (int *indices, const int rare_index, const int freq_index,
		       const Univcoord_T *rare, const int lenRare,
		       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
		       int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  /* 4 masks * 8 bits/mask = 32 bits */
  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare;  -- nextRare; */
    const int indexRare = rare_index + (rare - initRare);
    const int indexFreq = freq_index + (freq_low - initFreq);

    const vec Match = _mm512_set1_epi64(matchRare);

    if (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq_low + veclen * 32 > stopFreq) {
	freq_high += veclen * 32;
	freq_low += veclen * 32;
	goto FINISH_SCALAR;
      }
      while (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) { // if no match possible
	if (freq_low + veclen * (2 * offset) * 32 <= stopFreq) {
	  offset *= 2;
	} else if (freq_low + veclen * (offset + 1) * 32 <= stopFreq) {
	  offset = (int) ((stopFreq - freq_low) / (veclen * 32));
	  /* offset += 1; */
	  if (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) {
	    freq_high += veclen * offset * 32;
	    freq_low += veclen * offset * 32;
	    goto FINISH_SCALAR;
	  } else {
	    break;
	  }
	} else {
	  freq_high += veclen * offset * 32;
	  freq_low += veclen * offset * 32;
	  goto FINISH_SCALAR;
	}
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
	/* const int mid = (lower + offset) / 2; */
	const int mid = lower + ((offset - lower) / 2);
	if (GETPOS_I(freq_high,freq_low,veclen * mid * 32 + veclen * 31 + vecmax) < matchRare) {
	  lower = mid;
	} else {
	  offset = mid;
	}
      }
      freq_high += veclen * offset * 32;
      freq_low += veclen * offset * 32;
    }

    __m128i H;
    vec F0, F1, F2, F3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	if (GETPOS_I(freq_high,freq_low,veclen * 3 + vecmax) >= matchRare) {
	  base = 0 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3)));

	} else {
	  base = 1 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5)));
	
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7)));
	}

      } else {
	if (GETPOS_I(freq_high,freq_low,veclen * 11 + vecmax) >= matchRare) {
	  base = 2 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11)));

	} else {
	  base = 3 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15)));
	}
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	if (GETPOS_I(freq_high,freq_low,veclen * 19 + vecmax) >= matchRare) {
	  base = 4 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 8*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1 + 16)));
	  
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 9*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2 + 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3 + 16)));

	} else {
	  base = 5 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 10*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5 + 16)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 11*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6 + 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7 + 16)));
	}

      } else {
	if (GETPOS_I(freq_high,freq_low,veclen * 27 + vecmax) >= matchRare) {
	  base = 6 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 12*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8 + 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9 + 16)));

	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 13*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10+ 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11+ 16)));

	} else {
	  base = 7 * 32;
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 14*16)); /* Reads 16 chars */
	  F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12+ 16)));
	  F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13+ 16)));
	  
	  H = _mm_lddqu_si128((const __m128i *)(freq_high + 15*16)); /* Reads 16 chars */
	  F2 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14+ 16)));
	  F3 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
				_mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15+ 16)));
	}
      }
    }

    idx = _mm512_cmpeq_epi64_mask(F3, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F2, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F1, Match);
    idx = (idx << 8) + _mm512_cmpeq_epi64_mask(F0, Match);

    if (idx == 0) {
      /* Skip */
    } else {
      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices, rare_index, freq_index,
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}


/* rare_first_p is true */
static int
SIMDgalloping_anchorB (int *indices, const int rare_index, const int freq_index,
		       const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
		       const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  -- nextRare; */
    const int indexRare = rare_index + (rare_low - initRare);
    const int indexFreq = freq_index + (freq - initFreq);

    const vec Match = _mm512_set1_epi64(matchRare);
    if (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq + veclen  * 32 > stopFreq) {
	freq += veclen * 32;
	goto FINISH_SCALAR;
      }
      while (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) { // if no match possible
	if (freq + veclen * (2 * offset ) * 32 <= stopFreq) {
	  offset *= 2;
	} else if (freq + veclen * (offset + 1) * 32 <= stopFreq) {
	  offset = (int) ((stopFreq - freq ) / (veclen * 32));
	  /* offset += 1; */
	  if (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) {
	    freq += veclen * offset * 32;
	    goto FINISH_SCALAR;
	  } else {
	    break;
	  }
	} else {
	  freq += veclen * offset * 32;
	  goto FINISH_SCALAR;
	}
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
	/* const int mid = (lower + offset) / 2; */
	const int mid = lower + ((offset - lower) / 2);
	if (freq[veclen * mid * 32 + veclen * 31 + vecmax] < matchRare) {
	  lower = mid;
	} else {
	  offset = mid;
	}
      }
      freq += veclen * offset * 32;
    }

    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	if (freq[veclen * 3 + vecmax] >= matchRare) {
	  base = 0 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 3), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 2), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 1), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 0), Match);
	} else {
	  base = 1 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 7), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 6), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 5), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 4), Match);
	}
      } else {
	if (freq[veclen * 11 + vecmax] >= matchRare) {
	  base = 2 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 11), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 10), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 9), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 8), Match);
	} else {
	  base = 3 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 15), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 14), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 13), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 12), Match);
	}
      }
    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	if (freq[veclen * 19 + vecmax] >= matchRare) {
	  base = 4 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 3 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 2 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 1 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 0 + 16), Match);
	} else {
	  base = 5 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 7 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 6 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 5 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 4 + 16), Match);
	}
      } else {
	if (freq[veclen * 27 + vecmax] >= matchRare) {
	  base = 6 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 11 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 10+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 9 + 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 8 + 16), Match);
	} else {
	  base = 7 * 32;
	  idx = _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 15+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 14+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 13+ 16), Match);
	  idx = (idx << 8) + _mm512_cmpeq_epi64_mask(_mm512_loadu_si512((const vec *)(freq) + 12+ 16), Match);
	}
      }
    }

    if (idx == 0) {
      /* Skip */
    } else {
      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}

#elif defined(HAVE_AVX2)
/* rare_first_p is false */
static int
SIMDgalloping_anchorA (int *indices, const int rare_index, const int freq_index,
		       const Univcoord_T *rare, const int lenRare,
		       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
		       int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare;  -- nextRare; */

    const vec Match = _mm256_set1_epi32(matchRare);
    if (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq_low + veclen  * 32 > stopFreq) {
	freq_high += veclen * 32;
	freq_low += veclen * 32;
	goto FINISH_SCALAR;
      }
      while (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) { // if no match possible
	if (freq_low + veclen * (2 * offset ) * 32 <= stopFreq) {
	  offset *= 2;
	} else if (freq_low + veclen * (offset + 1) * 32 <= stopFreq) {
	  offset = (int) ((stopFreq - freq_low) / (veclen * 32));
	  /* offset += 1; */
	  if (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) {
	    freq_high += veclen * offset * 32;
	    freq_low += veclen * offset * 32;
	    goto FINISH_SCALAR;
	  } else {
	    break;
	  }
	} else {
	  freq_high += veclen * offset * 32;
	  freq_low += veclen * offset * 32;
	  goto FINISH_SCALAR;
	}
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
	/* const int mid = (lower + offset) / 2; */
	const int mid = lower + ((offset - lower) / 2);
	if (GETPOS_I(freq_high,freq_low,veclen * mid * 32 + veclen * 31 + vecmax) < matchRare) {
	  lower = mid;
	} else {
	  offset = mid;
	}
      }
      freq_high += veclen * offset * 32;
      freq_low += veclen * offset * 32;
    }

    __m128i H;
    vec F0, F1, F2, F3, F4, F5, F6, F7;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	base = 0 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));

      } else {
	base = 1 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	base = 2 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));

      } else {
	base = 3 * 32;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	F2 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F3 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F4 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F5 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	F6 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F7 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
      }
    }

    vec Q, Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Q0 = _mm256_cmpeq_epi64(F0, Match);
    Q1 = _mm256_cmpeq_epi64(F1, Match);
    Q2 = _mm256_cmpeq_epi64(F2, Match);
    Q3 = _mm256_cmpeq_epi64(F3, Match);
    Q4 = _mm256_cmpeq_epi64(F4, Match);
    Q5 = _mm256_cmpeq_epi64(F5, Match);
    Q6 = _mm256_cmpeq_epi64(F6, Match);
    Q7 = _mm256_cmpeq_epi64(F7, Match);

    Q = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3)),
			_mm256_or_si256(_mm256_or_si256(Q4, Q5),_mm256_or_si256(Q6, Q7)));

    if (_mm256_testz_si256(Q, Q)) {
      /* Skip */
    } else {
      idx = _mm256_movemask_pd((__m256d) Q7);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q6);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q5);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q4);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q3);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q2);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q1);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) Q0);

      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices,
		   rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}

/* rare_first_p is true */
static int
SIMDgalloping_anchorB (int *indices, const int rare_index, const int freq_index,
		       const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
		       const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  -- nextRare; */

    const vec Match = _mm256_set1_epi64x(matchRare);
    if (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq + veclen  * 32 > stopFreq) {
	freq += veclen * 32;
	goto FINISH_SCALAR;
      }
      while (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) { // if no match possible
	if (freq + veclen * (2 * offset ) * 32 <= stopFreq) {
	  offset *= 2;
	} else if (freq + veclen * (offset + 1) * 32 <= stopFreq) {
	  offset = (int) ((stopFreq - freq ) / (veclen * 32));
	  /* offset += 1; */
	  if (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) {
	    freq += veclen * offset * 32;
	    goto FINISH_SCALAR;
	  } else {
	    break;
	  }
	} else {
	  freq += veclen * offset * 32;
	  goto FINISH_SCALAR;
	}
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
	/* const int mid = (lower + offset) / 2; */
	const int mid = lower + ((offset - lower) / 2);
	if (freq[veclen * mid * 32 + veclen * 31 + vecmax] < matchRare) {
	  lower = mid;
	} else {
	  offset = mid;
	}
      }
      freq += veclen * offset * 32;
    }

    vec M, M0, M1, M2, M3, M4, M5, M6, M7;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	base = 0 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 0), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 1), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 2), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 3), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 4), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 5), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 6), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 7), Match);

      } else {
	base = 1 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 8), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 9), Match);
	M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15), Match);
      }

    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	base = 2 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 0 + 16), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 1 + 16), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 2 + 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 3 + 16), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 4 + 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 5 + 16), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 6 + 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 7 + 16), Match);

      } else {
	base = 3 * 32;
        M0 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 8 + 16), Match);
	M1 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 9 + 16), Match);
        M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10 + 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11 + 16), Match);
        M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12 + 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13 + 16), Match);
        M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14 + 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15 + 16), Match);
      }
    }

    M = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(M0, M1),_mm256_or_si256(M2, M3)),
			_mm256_or_si256(_mm256_or_si256(M4, M5),_mm256_or_si256(M6, M7)));
    if (_mm256_testz_si256(M, M)) {
      /* Skip */
    } else {
      idx = _mm256_movemask_pd((__m256d) M7);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M6);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M5);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M4);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M3);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M2);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M1);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M0);

      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE4_1)
/* rare_first_p is false */
static int
SIMDgalloping_anchorA (int *indices, const int rare_index, const int freq_index,
		       const Univcoord_T *rare, const int lenRare,
		       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
		       int delta_rare, int diagterm_rare) {
  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;

  /* Measuring only 16 bits, but procedures expect 32 bits */
  uint32_t idx;
  int base;

  if (freq_low > stopFreq) {
    return scalar_anchorA(indices, rare_index, freq_index,
			  rare, lenRare, freq_high, freq_low, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare; */

    const vec Match = _mm_set1_epi32(matchRare);
    if (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq_low + veclen * 32 > stopFreq) {
        freq_high += veclen * 32;
        freq_low += veclen * 32;
        goto FINISH_SCALAR;
      }
      while (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) { // if no match possible
        if (freq_low + veclen * (2 * offset) * 32 <= stopFreq) {
          offset *= 2;
        } else if (freq_low + veclen * (offset + 1) * 32 <= stopFreq) {
          offset = (int) ((stopFreq - freq_low) / (veclen * 32));
          /* offset += 1; */
          if (GETPOS_I(freq_high,freq_low,veclen * offset * 32 + veclen * 31 + vecmax) < matchRare) {
            freq_high += veclen * offset * 32;
            freq_low += veclen * offset * 32;
            goto FINISH_SCALAR;
          } else {
            break;
          }
        } else {
          freq_high += veclen * offset * 32;
          freq_low += veclen * offset * 32;
          goto FINISH_SCALAR;
        }
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
        /* const int mid = (lower + offset) / 2; */
        const int mid = lower + (offset - lower) / 2;
        if (GETPOS_I(freq_high,freq_low,veclen * mid * 32 + veclen * 31 + vecmax) < matchRare) {
          lower = mid;
	} else {
          offset = mid;
	}
      }
      freq_high += veclen * offset * 32;
      freq_low += veclen * offset * 32;
    }

    __m128i H;
    vec F0, F1, F2, F3, F4, F5, F6, F7;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) >= matchRare) {
	base = 0 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));

      } else {
	base = 1 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
      }

    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) >= matchRare) {
	base = 2 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));

      } else {
	base = 3 * 16;
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
	F4 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F5 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	F6 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F7 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
      }
    }

    vec Q, Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7;
    Q0 = _mm_cmpeq_epi64(F0, Match);
    Q1 = _mm_cmpeq_epi64(F1, Match);
    Q2 = _mm_cmpeq_epi64(F2, Match);
    Q3 = _mm_cmpeq_epi64(F3, Match);
    Q4 = _mm_cmpeq_epi64(F4, Match);
    Q5 = _mm_cmpeq_epi64(F5, Match);
    Q6 = _mm_cmpeq_epi64(F6, Match);
    Q7 = _mm_cmpeq_epi64(F7, Match);

    Q = _mm_or_si128(_mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3)),
		     _mm_or_si128(_mm_or_si128(Q4, Q5), _mm_or_si128(Q6, Q7)));

    if (_mm_testz_si128(Q, Q)) {
      /* Skip */
    } else {
      idx = _mm_movemask_pd((__m128d) Q7);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q6);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q5);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q4);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q3);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q2);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q1);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) Q0);

      /* *indices++ = freq_index + (freq_low + base + count_trailing_zeroes_32(idx)) - initFreq; -- want only indices2 */
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorA(indices,
		   rare_index + (rare - initRare), freq_index + (freq_low - initFreq),
		   rare, stopRare + rarespace - rare,
		   freq_high, freq_low, stopFreq + freqspace - freq_low,
		   delta_rare, diagterm_rare);
}

/* rare_first_p is true */
static int
SIMDgalloping_anchorB (int *indices, const int rare_index, const int freq_index,
		       const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
		       const Univcoord_T *freq, const int lenFreq,
		       int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  /* Measuring only 16 bits, but procedures expect 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar_anchorB(indices, rare_index, freq_index,
			  rare_high, rare_low, lenRare, freq, lenFreq,
			  delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const Univcoord_T matchRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* nextRare; */
    /* const Univcoord_T outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- nextRare; */

    const vec Match = _mm_set1_epi64((__m64) matchRare);
    if (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      int offset = 1;
      if (freq + veclen * 32 > stopFreq) {
        freq += veclen * 32;
        goto FINISH_SCALAR;
      }
      while (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) { // if no match possible
        if (freq + veclen * (2 * offset) * 32 <= stopFreq) {
          offset *= 2;
        } else if (freq + veclen * (offset + 1) * 32 <= stopFreq) {
          offset = (int) ((stopFreq - freq) / (veclen * 32));
          /* offset += 1; */
          if (freq[veclen * offset * 32 + veclen * 31 + vecmax] < matchRare) {
            freq += veclen * offset * 32;
            goto FINISH_SCALAR;
          } else {
            break;
          }
        } else {
          freq += veclen * offset * 32;
          goto FINISH_SCALAR;
        }
      }

      int lower = offset / 2;
      while (lower + 1 != offset) {
        /* const int mid = (lower + offset) / 2; */
        const int mid = lower + (offset - lower) / 2;
        if (freq[veclen * mid * 32 + veclen * 31 + vecmax] < matchRare) {
          lower = mid;
	} else {
          offset = mid;
	}
      }
      freq += veclen * offset * 32;
    }

    vec M, M0, M1, M2, M3, M4, M5, M6, M7;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] >= matchRare) {
	base = 0 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 0), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 1), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 2), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 3), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 4), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 5), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 6), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 7), Match);

      } else {
	base = 1 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 8), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 9), Match);
	M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 10), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 11), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 12), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 13), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 14), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 15), Match);
      }
    } else {
      if (freq[veclen * 23 + vecmax] >= matchRare) {
	base = 2 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 0 + 16), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 1 + 16), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 2 + 16), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 3 + 16), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 4 + 16), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 5 + 16), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 6 + 16), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 7 + 16), Match);

      } else {
	base = 3 * 16;
        M0 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 8 + 16), Match);
	M1 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 9 + 16), Match);
        M2 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 10 + 16), Match);
	M3 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 11 + 16), Match);
        M4 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 12 + 16), Match);
	M5 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 13 + 16), Match);
        M6 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 14 + 16), Match);
	M7 = _mm_cmpeq_epi64(_mm_loadu_si128((const vec *)(freq) + 15 + 16), Match);
      }
    }

    M = _mm_or_si128(_mm_or_si128(_mm_or_si128(M0, M1),_mm_or_si128(M2, M3)),
			_mm_or_si128(_mm_or_si128(M4, M5),_mm_or_si128(M6, M7)));

    if (_mm_testz_si128(M, M)) {
      /* Skip */
    } else {
      idx = _mm_movemask_pd((__m128d) M7);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M6);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M5);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M4);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M3);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M2);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M1);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M0);

      /* *indices++ = rare_index + (rare_low - initRare); -- want only indices2 */
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    }
  }

 FINISH_SCALAR:
  return (indices - initout) +
    scalar_anchorB(indices,
		   rare_index + (rare_low - initRare), freq_index + (freq - initFreq),
		   rare_high, rare_low, stopRare + rarespace - rare_low,
		   freq, stopFreq + freqspace - freq,
		   delta_rare, diagterm_rare);
}
#endif


/**
 * Our main heuristic.
 *
 * The out pointer can be set1 if length1<=length2,
 * or else it can be set2 if length1>length2.
 */
/* set1 + diagterm1 == set2 + diagterm2 */
#ifdef HAVE_SSE4_1
int
Intersect_indices2_large (int *indices,
			  const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
			  const Univcoord_T *set2, const int length2, int diagterm2) {
  int nindices;

  debug(printf("Entered Intersect_indices2_large with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%d %u\n",set1_high[i],set1_low[i]);
  }
  printf("\n");
  for (i = 0; i < length2; i++) {
    printf("%u\n",set2[i]);
  }
  printf("\n");
#endif    

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else if ((1000 * length1 <= length2) || (1000 * length2 <= length1)) {
    debug(printf("Using SIMDgalloping method\n"));
    if (length1 <= length2) {
      nindices = SIMDgalloping_anchorB(indices, /*rare_index*/0, /*freq_index*/0,
				       set1_high, set1_low, length1, set2, length2,
				       /*delta_rare*/(diagterm1 - diagterm2),
				       /*diagterm_rare*/diagterm1);
    } else {
      nindices = SIMDgalloping_anchorA(indices, /*rare_index*/0, /*freq_index*/0,
				       set2, length2, set1_high, set1_low, length1,
				       /*delta_rare*/(diagterm2 - diagterm1),
				       /*diagterm_rare*/diagterm2);
    }

    return nindices;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    debug(printf("Using v3 method\n"));
    if (length1 <= length2) {
      nindices = v3_anchorB(indices, /*rare_index*/0, /*freq_index*/0,
			    set1_high, set1_low, length1, set2, length2,
			    /*delta_rare*/(diagterm1 - diagterm2),
			    /*diagterm_rare*/diagterm1);
    } else {
      nindices = v3_anchorA(indices, /*rare_index*/0, /*freq_index*/0,
			    set2, length2, set1_high, set1_low, length1,
			    /*delta_rare*/(diagterm2 - diagterm1),
			    /*diagterm_rare*/diagterm2);
    }

    return nindices;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      nindices = v1_anchorB(indices, /*rare_index*/0, /*freq_index*/0,
			    set1_high, set1_low, length1, set2, length2,
			    /*delta_rare*/(diagterm1 - diagterm2),
			    /*diagterm_rare*/diagterm1);
    } else {
      nindices = v1_anchorA(indices, /*rare_index*/0, /*freq_index*/0,
			    set2, length2, set1_high, set1_low, length1,
			    /*delta_rare*/(diagterm2 - diagterm1),
			    /*diagterm_rare*/diagterm2);
    }

    return nindices;
  }
}

#else

int
Intersect_indices2_large (int *indices,
			  const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
			  const Univcoord_T *set2, const int length2, int diagterm2) {
  int nindices;

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    if (length1 <= length2) {
      nindices = galloping_intersection_anchorB(indices, /*rare_index*/0, /*freq_index*/0,
						set1_high, set1_low, length1, set2, length2,
						/*delta_rare*/(diagterm1 - diagterm2),
						/*diagterm_rare*/diagterm1);
    } else {
      nindices = galloping_intersection_anchorA(indices, /*rare_index*/0, /*freq_index*/0,
						set2, length2, set1_high, set1_low, length1,
						/*delta_rare*/(diagterm2 - diagterm1),
						/*diagterm_rare*/diagterm2);
    }

    return nindices;

  } else {
    if (length1 <= length2) {
      nindices = scalar_anchorB(indices, /*rare_index*/0, /*freq_index*/0,
				set1_high, set1_low, length1, set2, length2,
				/*deltaA*/(diagterm1 - diagterm2),
				/*diagtermA*/diagterm1);
    } else {
      nindices = scalar_anchorA(indices, /*rare_index*/0, /*freq_index*/0,
				set2, length2, set1_high, set1_low, length1,
				/*deltaA*/(diagterm2 - diagterm1),
				/*diagtermA*/diagterm2);
    }

    return nindices;
  }
}
#endif


void
Intersect_indices2_large_setup () {
  initialize_match_values();
  return;
}  
