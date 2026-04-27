static char rcsid[] = "$Id: 2f085d27c2c58dc320203efd1c62df4bb70f25d8 $";
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

/* Minimum requirement for SIMD is SSE4_2 for _mm_cmpgt_epi64 */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "bool.h"
#include "genomebits.h"
#include "intersect-indices-uint8.h"

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

#define WRITE_BOTH 1


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


#ifndef HAVE_SSE4_2
static int
__frogadvanceUntil(const Univcoord_T *array, const int pos,
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

static int
galloping_intersection (int *indices, const int smallset_index, const int largeset_index,
			const Univcoord_T *smallset, const int smalllength,
			const Univcoord_T *largeset, const int largelength,
			int delta_small, int diagterm_small, bool smallset_first_p) {
  if (smalllength == 0) return 0;
  const int *initout = indices;
  int k1 = 0, k2 = 0;
  while (1) {
    if (largeset[k1] < smallset[k2] + delta_small) {
      k1 = __frogadvanceUntil(largeset, k1, largelength, smallset[k2] + delta_small);
      if (k1 == largelength) {
        break;
      }
    }
  midpoint:
    if (smallset[k2] + delta_small < largeset[k1]) {
      ++k2;
      if (k2 == smalllength) {
        break;
      }
    } else {
      if (smallset_first_p == true) {
	*indices++ = smallset_index + k2;
	*indices++ = largeset_index + k1;
      } else {
	*indices++ = largeset_index + k1;
	*indices++ = smallset_index + k2;
      }
      ++k2;
      if (k2 == smalllength) {
        break;
      }
      k1 = __frogadvanceUntil(largeset, k1, largelength, smallset[k2] + delta_small);
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
static int
scalar (int *indices, const int A_index, const int B_index,
	const Univcoord_T *A, const int lenA,
	const Univcoord_T *B, const int lenB, int deltaA, int diagtermA,
	bool A_first_p) {
  const int *initout = indices;
  if (lenA == 0 || lenB == 0) return 0;

  const Univcoord_T *initA = A;
  const Univcoord_T *initB = B;
  const Univcoord_T *endA = A + lenA;
  const Univcoord_T *endB = B + lenB;

  while (1) {
    while ((*A) + deltaA < *B) {
    SKIP_FIRST_COMPARE:
      if (++A == endA) {
        return (indices - initout);
      }
    }
    while ((*A) + deltaA > *B) {
      if (++B == endB) {
        return (indices - initout);
      }
    }
    if ((*A) + deltaA == *B) {
      if (A_first_p == true) {
	*indices++ = A_index + (A - initA);
	*indices++ = B_index + (B - initB);
      } else {
	*indices++ = B_index + (B - initB);
	*indices++ = A_index + (A - initA);
      }
      if (++A == endA || ++B == endB) {
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

#if defined(HAVE_SSE4_2)
static int
v1 (int *indices, const int rare_index, const int freq_index,
    const Univcoord_T *rare, int lenRare,
    const Univcoord_T *freq, int lenFreq, int delta_rare, int diagterm_rare,
    bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
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

  const Univcoord_T *stopFreq = &freq[lenFreq] - kFreqSpace;
  const Univcoord_T *stopRare = &rare[lenRare] - kRareSpace;


  if (COMPILER_RARELY((rare >= stopRare) || (freq >= stopFreq)))
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

  maxFreq = freq[2 * 4 - 1];
#ifdef HAVE_AVX512
  F = _mm512_loadu_si512((const __m512i *)(freq));
#elif defined(HAVE_AVX2)
  F0 = _mm256_loadu_si256((const __m256i *)(freq));
  F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
#else
  F0 = _mm_lddqu_si128((const __m128i *)(freq));
  F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
  F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
  F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
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
    idx = _mm256_movemask_pd((__m256d) M); /* pd corresponds to 64 bits */
    M = _mm256_cmpeq_epi64(F0, Rare);
    idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);

    Rare = _mm256_set1_epi64x(valRare);
#else
    M = _mm_cmpeq_epi32(F3, Rare);
    idx = _mm_movemask_pd((__m128d) M); /* pd corresponds to 64 bits */
    M = _mm_cmpeq_epi32(F2, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi32(F1, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    M = _mm_cmpeq_epi32(F0, Rare);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    Rare = _mm_set1_epi64((__m64) valRare);
#endif

    if ((bit = match_start[idx]) >= 0) {
      assert(match_n[idx] == 1);
      if (rare_first_p == true) {
	*indices++ = rare_index + ((rare - 1) - initRare); /* rare is advanced by 1 */
	*indices++ = freq_index + (freq + bit - initFreq);
      } else {
	*indices++ = freq_index + (freq + bit - initFreq);
	*indices++ = rare_index + ((rare - 1) - initRare); /* rare is advanced by 1 */
      }
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
  F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
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
  lenRare = stopRare + kRareSpace - rare;

  int tail = scalar(indices, rare_index + (rare - initRare), freq_index + (freq - initFreq),
		    rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);

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
/* Processes 32 vectors, each vector having 8 Univcoord_Ts => 256 Univcoord_Ts */
static int
v3 (int *indices, const int rare_index, const int freq_index,
    const Univcoord_T *rare, const int lenRare,
    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare,
    bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;
  const Univcoord_T *stopRare = rare + lenRare - rarespace;

  /* 4 masks * 8 bits/mask = 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

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
    } else if (rare_first_p == true) {
      *indices++ = rare_index + (rare - initRare);
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    } else {
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
}

#elif defined(HAVE_AVX2)
/* Processes 32 vectors, each vector having 4 Univcoord_Ts => 128 Univcoord_Ts */
static int
v3 (int *indices, const int rare_index, const int freq_index,
    const Univcoord_T *rare, const int lenRare,
    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare,
    bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;
  const Univcoord_T *stopRare = rare + lenRare - rarespace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare */

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
	M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10+ 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11+ 16), Match);
	M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12+ 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13+ 16), Match);
	M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14+ 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15+ 16), Match);
      }
    }

    M = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(M0, M1), _mm256_or_si256(M2, M3)),
			_mm256_or_si256(_mm256_or_si256(M4, M5), _mm256_or_si256(M6, M7)));

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

      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      } else {
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
	*indices++ = rare_index + (rare - initRare);
      }
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
}

#elif defined(HAVE_SSE4_2)
/* Processes 32 vectors, each vector having 2 unsigned ints => 64 unsigned ints */
static int
v3 (int *indices, const int rare_index, const int freq_index,
    const Univcoord_T *rare, const int lenRare,
    const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare,
    bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;
  const Univcoord_T *stopRare = rare + lenRare - rarespace;

  /* Checking only 16 bits, but bit procedures are expecting 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare */
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
      if (freq[veclen * 7 + vecmax] < matchRare) {
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
      if (freq[veclen * 23 + vecmax] < matchRare) {
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

    M = _mm_or_si128(_mm_or_si128(_mm_or_si128(M0, M1), _mm_or_si128(M2, M3)),
		     _mm_or_si128(_mm_or_si128(M4, M5), _mm_or_si128(M6, M7)));

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

      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      } else {
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
	*indices++ = rare_index + (rare - initRare);
      }
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
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
 * The matchOut pointer can safely be equal to the rare pointer.
 *
 * This function DOES NOT use assembly. It only relies on intrinsics.
 */

#ifdef HAVE_AVX512
static int
SIMDgalloping (int *indices, const int rare_index, const int freq_index,
	       const Univcoord_T *rare, const int lenRare,
	       const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare,
	       bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopRare = rare + lenRare - rarespace;
  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;

  /* 4 masks * 8 bits/mask = 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare;  -- nextRare; */

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
    } else if (rare_first_p == true) {
      *indices++ = rare_index + (rare - initRare);
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
    } else {
      *indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      *indices++ = rare_index + (rare - initRare);
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
}

#elif defined(HAVE_AVX2)
static int
SIMDgalloping (int *indices, const int rare_index, const int freq_index,
	       const Univcoord_T *rare, const int lenRare,
	       const Univcoord_T *freq, const int lenFreq, int delta_rare, int diagterm_rare,
	       bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T); /* 8 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;
  const Univcoord_T *stopRare = rare + lenRare - rarespace;

  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare;  /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare;  -- nextRare; */

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
	M2 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 10+ 16), Match);
	M3 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 11+ 16), Match);
	M4 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 12+ 16), Match);
	M5 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 13+ 16), Match);
	M6 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 14+ 16), Match);
	M7 = _mm256_cmpeq_epi64(_mm256_loadu_si256((const vec *)(freq) + 15+ 16), Match);
      }
    }

    M = _mm256_or_si256(_mm256_or_si256(_mm256_or_si256(M0, M1), _mm256_or_si256(M2, M3)),
			_mm256_or_si256(_mm256_or_si256(M4, M5), _mm256_or_si256(M6, M7)));

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

      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      } else {
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
	*indices++ = rare_index + (rare - initRare);
      }
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
}

#elif defined(HAVE_SSE4_2)
static int
SIMDgalloping (int *indices, const int rare_index, const int freq_index,
	       const Univcoord_T *rare, const int lenRare,
	       const Univcoord_T *freq, const int lenFreq,
	       int delta_rare, int diagterm_rare, bool rare_first_p) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const int *initout = indices;
  const Univcoord_T *initRare = rare;
  const Univcoord_T *initFreq = freq;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(Univcoord_T);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const Univcoord_T *stopFreq = freq + lenFreq - freqspace;
  const Univcoord_T *stopRare = rare + lenRare - rarespace;

  /* Checking only 16 bits, but bit procedures are expecting 32 bits */
  uint32_t idx;
  int base;

  if (freq > stopFreq) {
    return scalar(indices, rare_index, freq_index,
		  rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare, rare_first_p);
  }
  for (; rare < stopRare; ++rare) {
    const Univcoord_T matchRare = (*rare) + delta_rare; /* nextRare; */
    /* const Univcoord_T outRare = (*rare) + diagterm_rare; -- nextRare; */

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
      if (freq[veclen * 7 + vecmax] < matchRare) {
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
      if (freq[veclen * 23 + vecmax] < matchRare) {
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

    M = _mm_or_si128(_mm_or_si128(_mm_or_si128(M0, M1), _mm_or_si128(M2, M3)),
		     _mm_or_si128(_mm_or_si128(M4, M5), _mm_or_si128(M6, M7)));

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

      if (rare_first_p == true) {
	*indices++ = rare_index + (rare - initRare);
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
      } else {
	*indices++ = freq_index + (freq + base + count_trailing_zeroes_32(idx)) - initFreq;
	*indices++ = rare_index + (rare - initRare);
      }
    }
  }

 FINISH_SCALAR:
  return (indices - initout) + scalar(indices,
				      rare_index + (rare - initRare), freq_index + (freq - initFreq),
				      rare, stopRare + rarespace - rare,
				      freq, stopFreq + freqspace - freq,
				      delta_rare, diagterm_rare, rare_first_p);
}
#endif


/**
 * Our main heuristic.
 *
 * The out pointer can be set1 if length1<=length2,
 * or else it can be set2 if length1>length2.
 */
/* set1 + diagterm1 == set2 + diagterm2 */
#ifdef HAVE_SSE4_2
int
Intersect_indices_uint8 (int *indices,
			 const Univcoord_T *set1, const int length1, int diagterm1,
			 const Univcoord_T *set2, const int length2, int diagterm2) {
  int nindices;

  debug(printf("Entered Intersect_exact_indices_univcoord8 with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%u\n",set1[i]);
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
      nindices = SIMDgalloping(indices, /*rare_index*/0, /*freq_index*/0,
			       set1, length1, set2, length2,
			       /*delta_rare*/(diagterm1 - diagterm2),
			       /*diagterm_rare*/diagterm1,/*rare_first_p*/true);
    } else {
      nindices = SIMDgalloping(indices, /*rare_index*/0, /*freq_index*/0,
			       set2, length2, set1, length1,
			       /*delta_rare*/(diagterm2 - diagterm1),
			       /*diagterm_rare*/diagterm2,/*rare_first_p*/false);
    }

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    debug(printf("Using v3 method\n"));
    if (length1 <= length2) {
      nindices = v3(indices, /*rare_index*/0, /*freq_index*/0,
		    set1, length1, set2, length2,
		    /*delta_rare*/(diagterm1 - diagterm2),
		    /*diagterm_rare*/diagterm1,/*rare_first_p*/true);
    } else {
      nindices = v3(indices, /*rare_index*/0, /*freq_index*/0,
		    set2, length2, set1, length1,
		    /*delta_rare*/(diagterm2 - diagterm1),
		    /*diagterm_rare*/diagterm2,/*rare_first_p*/false);
    }

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      nindices = v1(indices, /*rare_index*/0, /*freq_index*/0,
		    set1, length1, set2, length2,
		    /*delta_rare*/(diagterm1 - diagterm2),
		    /*diagterm_rare*/diagterm1,/*rare_first_p*/true);
    } else {
      nindices = v1(indices, /*rare_index*/0, /*freq_index*/0,
		    set2, length2, set1, length1,
		    /*delta_rare*/(diagterm2 - diagterm1),
		    /*diagterm_rare*/diagterm2,/*rare_first_p*/false);
    }
  }

#ifdef WRITE_BOTH
  return (nindices / 2);
#else
  return nindices;
#endif
}

#else

int
Intersect_indices_uint8 (int *indices,
			 const Univcoord_T *set1, const int length1, int diagterm1,
			 const Univcoord_T *set2, const int length2, int diagterm2) {
  int nindices;

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    if (length1 <= length2) {
      nindices = galloping_intersection(indices, /*rare_index*/0, /*freq_index*/0,
					set1, length1, set2, length2,
					/*delta_rare*/(diagterm1 - diagterm2),
					/*diagterm_rare*/diagterm1,/*A_first_p*/true);
    } else {
      nindices = galloping_intersection(indices, /*rare_index*/0, /*freq_index*/0,
					set2, length2, set1, length1,
					/*delta_rare*/(diagterm2 - diagterm1),
					/*diagterm_rare*/diagterm2,/*A_first_p*/false);
    }

  } else {
    if (length1 <= length2) {
      nindices = scalar(indices, /*rare_index*/0, /*freq_index*/0,
			set1, length1, set2, length2,
			/*deltaA*/(diagterm1 - diagterm2),
			/*diagtermA*/diagterm1,/*A_first_p*/true);
    } else {
      nindices = scalar(indices, /*rare_index*/0, /*freq_index*/0,
			set2, length2, set1, length1,
			/*deltaA*/(diagterm2 - diagterm1),
			/*diagtermA*/diagterm2,/*A_first_p*/false);
    }
  }

#ifdef WRITE_BOTH
  return (nindices / 2);
#else
  return nindices;
#endif
}

#endif


void
Intersect_indices_uint8_setup () {
  initialize_match_values();
  return;
}  
