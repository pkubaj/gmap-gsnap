static char rcsid[] = "$Id: 7a81987e19a8bad50026f6e7ee6b473d8039cc36 $";
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

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "intersect-small.h"

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

/* binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Check SIMD against non-SIMD */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif



#ifndef HAVE_SSE2
static int
__frogadvanceUntil(const unsigned int *array, const int pos,
		   const int length, const unsigned int min) {
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
galloping_intersection (unsigned int *out, const unsigned int *smallset, const int smalllength,
			const unsigned int *largeset, const int largelength,
			int delta_small, int diagterm_small) {
  if (smalllength == 0) return 0;
  const unsigned int *initout = out;
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
      *out++ = smallset[k2] + diagterm_small;
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

  return out - initout;
}
#endif



/* Fast scalar scheme designed by N. Kurz. */
/* Delta is required to be added to A to become equivalent with B */
/* Desired result should be A + diagtermA (== B + diagtermB) */
static int
scalar (unsigned int *out, const unsigned int *A, const int lenA,
	const unsigned int *B, const int lenB, int deltaA, int diagtermA) {
  const unsigned int *initout = out;
  if (lenA == 0 || lenB == 0) return 0;

  const unsigned int *endA = A + lenA;
  const unsigned int *endB = B + lenB;

  while (1) {
    while ((*A) + deltaA < *B) {
    SKIP_FIRST_COMPARE:
      if (++A == endA) {
        return (out - initout);
      }
    }
    while ((*A) + deltaA > *B) {
      if (++B == endB) {
        return (out - initout);
      }
    }
    if ((*A) + deltaA == *B) {
      *out++ = (*A) + diagtermA;
      if (++A == endA || ++B == endB) {
        return (out - initout);
      }
    } else {
      goto SKIP_FIRST_COMPARE;
    }
  }

  return (out - initout); 	/* Not reached */
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

#if defined(HAVE_SSE2)
static int
v1 (unsigned int *matchOut, const unsigned int *rare, int lenRare,
    const unsigned int *freq, int lenFreq, int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const unsigned int *matchOrig = matchOut;

  unsigned int valRare, outRare;
  unsigned int maxProbe;		/* was uint64_t */
  unsigned int maxFreq;		/* was uint64_t */
  const int kRareSpace = 0;	/* was uint64_t */

#ifdef HAVE_AVX512
  const int kFreqSpace = 2 * 8 * (0 + 1) - 1; /* 16 32-mers in 512 bits */

  __m512i Rare;
  __m512i F;
  __mmask16 cmp_mask;
#elif defined(HAVE_AVX2)
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 8 32-mers in 256 bits */

  __m256i Rare;
  __m256i F;
  __m256i CMP;
#else
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 8 32-mers in 128+128 bits */

  __m128i Rare;
  __m128i F0, F1;
  __m128i CMP0, CMP1;
#endif

  const unsigned int *stopFreq = &freq[lenFreq] - kFreqSpace;
  const unsigned int *stopRare = &rare[lenRare] - kRareSpace;


  if (COMPILER_RARELY((rare >= stopRare) || (freq >= stopFreq)))
    goto FINISH_SCALAR;

  valRare = (*rare) + delta_rare; /* For comparison */
  outRare = (*rare) + diagterm_rare; /* For output */
#ifdef HAVE_AVX512
  Rare = _mm512_set1_epi32(valRare);
#elif defined(HAVE_AVX2)
  Rare = _mm256_set1_epi32(valRare);
#else
  Rare = _mm_set1_epi32(valRare);
#endif

#ifdef HAVE_AVX512
  maxFreq = freq[2 * 8 - 1];
  F = _mm512_loadu_si512((const __m512i *)(freq));
#elif defined(HAVE_AVX2)
  maxFreq = freq[2 * 4 - 1];
  F = _mm256_loadu_si256((const __m256i *)(freq));
#elif defined(HAVE_SSE4_1)
  maxFreq = freq[2 * 4 - 1];
  F0 = _mm_lddqu_si128((const __m128i *)(freq));
  F1 = _mm_lddqu_si128((const __m128i *)(freq + 4));
#else
  maxFreq = freq[2 * 4 - 1];
  F0 = _mm_loadu_si128((const __m128i *)(freq));
  F1 = _mm_loadu_si128((const __m128i *)(freq + 4));
#endif

  if (COMPILER_RARELY(maxFreq < valRare))
    goto ADVANCE_FREQ;

 ADVANCE_RARE:
  do {
    *matchOut = outRare;
    rare += 1;
    if (COMPILER_RARELY(rare >= stopRare)) {
      rare -= 1;
      goto FINISH_SCALAR;
    }

    valRare = (*rare) + delta_rare; /* for next iteration */
    outRare = (*rare) + diagterm_rare; /* for next iteration */
#ifdef HAVE_AVX512
    cmp_mask = _mm512_cmpeq_epi32_mask(F,Rare);
    Rare = _mm512_set1_epi32(valRare);
    if (cmp_mask != 0) {
      matchOut ++;
    }

#elif defined(HAVE_AVX2)
    CMP = _mm256_cmpeq_epi32(F,Rare);
    Rare = _mm256_set1_epi32(valRare);
    if(_mm256_testz_si256(CMP,CMP) == 0) {
      matchOut ++;
    }

#elif HAVE_SSE4_1
    CMP0 = _mm_cmpeq_epi32(F0, Rare);
    CMP1 = _mm_cmpeq_epi32(F1, Rare);
    Rare = _mm_set1_epi32(valRare);
    CMP0 = _mm_or_si128(CMP0, CMP1);
    if (_mm_testz_si128(CMP0, CMP0) == 0) {
      matchOut++;
    }

#else
    CMP0 = _mm_cmpeq_epi32(F0, Rare);
    CMP1 = _mm_cmpeq_epi32(F1, Rare);
    Rare = _mm_set1_epi32(valRare);
    CMP0 = _mm_or_si128(CMP0, CMP1);
    if (_mm_movemask_epi8(CMP0)) {
      matchOut++;
    }
#endif

  } while (maxFreq >= valRare);

 ADVANCE_FREQ:
  do {
#ifdef HAVE_AVX512
    const int kProbe = (0 + 1) * 2 * 8; /* 16 32-mers in 512 bits */
    const unsigned int *probeFreq = freq + kProbe;

    if (COMPILER_RARELY(probeFreq >= stopFreq)) {
      goto FINISH_SCALAR;
    }
    maxProbe = freq[(0 + 2) * 2 * 8 - 1];
#else
    const int kProbe = (0 + 1) * 2 * 4; /* 8 32-mers in 256 bits */
    const unsigned int *probeFreq = freq + kProbe;

    if (COMPILER_RARELY(probeFreq >= stopFreq)) {
      goto FINISH_SCALAR;
    }
    maxProbe = freq[(0 + 2) * 2 * 4 - 1];
#endif

    freq = probeFreq;

  } while (maxProbe < valRare);

  maxFreq = maxProbe;

#ifdef HAVE_AVX512
  F = _mm512_loadu_si512((const __m512i *)(freq));
#elif defined(HAVE_AVX2)
  F = _mm256_loadu_si256((const __m256i *)(freq));
#elif defined(HAVE_SSE4_1)
  F0 = _mm_lddqu_si128((const __m128i *)(freq));
  F1 = _mm_lddqu_si128((const __m128i *)(freq + 4));
#else
  F0 = _mm_loadu_si128((const __m128i *)(freq));
  F1 = _mm_loadu_si128((const __m128i *)(freq + 4));
#endif

  goto ADVANCE_RARE;

  int count;
 FINISH_SCALAR:
  count = matchOut - matchOrig;

  lenFreq = stopFreq + kFreqSpace - freq;
  lenRare = stopRare + kRareSpace - rare;

  /* Was match_scalar */
  int tail = scalar(matchOut, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);

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
 * The matchOut pointer can safely be equal to the rare pointer.
 *
 * This function DOES NOT use inline assembly instructions. Just intrinsics.
 */

#ifdef HAVE_AVX512
/* Processes 32 vectors, each vector having 16 unsigned ints => 512 unsigned ints */
static int
v3 (unsigned int *out, const unsigned int *rare, const int lenRare,
    const unsigned int *freq, const int lenFreq, int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const unsigned int *initout = out;
  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(unsigned int);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const unsigned int *stopFreq = freq + lenFreq - freqspace;
  const unsigned int *stopRare = rare + lenRare - rarespace;

  if (freq > stopFreq) {
    return scalar(out, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const unsigned int matchRare = (*rare) + delta_rare; /* nextRare */
    const unsigned int outRare = (*rare) + diagterm_rare; /* nextRare */

    const vec Match = _mm512_set1_epi32(matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    __mmask16 Q0, Q1, Q2, Q3;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] < matchRare) {
	Q0 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 8), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 9), Match);
	Q1 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 10), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 11), Match);
	Q2 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 12), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 13), Match);
	Q3 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 14), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 15), Match);

      } else {
	Q0 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 4), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 5), Match);
	Q1 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 6), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 7), Match);
	Q2 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 0), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 1), Match);
	Q3 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 2), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 3), Match);
      }

    } else {
      if (freq[veclen * 23 + vecmax] < matchRare) {
	Q0 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 8 + 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 9 + 16), Match);
	Q1 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 10+ 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 11+ 16), Match);
	Q2 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 12+ 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 13+ 16), Match);
	Q3 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 14+ 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 15+ 16), Match);

      } else {
	Q0 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 4 + 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 5 + 16), Match);
	Q1 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 6 + 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 7 + 16), Match);
	Q2 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 0 + 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 1 + 16), Match);
	Q3 =
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 2 + 16), Match) |
	  _mm512_cmpeq_epi32_mask(_mm512_loadu_si512((const vec *)(freq) + 3 + 16), Match);
      }
    }
    const __mmask16 Q = (Q0 | Q1) | (Q2 | Q3);
    if (Q == 0) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR: return (out - initout) + scalar(out, rare, stopRare + rarespace - rare,
						freq, stopFreq + freqspace - freq,
						delta_rare, diagterm_rare);
}

#elif defined(HAVE_AVX2)
/* Processes 32 vectors, each vector having 4 unsigned ints => 128 unsigned ints */
static int
v3 (unsigned int *out, const unsigned int *rare, const int lenRare,
    const unsigned int *freq, const int lenFreq, int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const unsigned int *initout = out;
  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(unsigned int);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const unsigned int *stopFreq = freq + lenFreq - freqspace;
  const unsigned int *stopRare = rare + lenRare - rarespace;

  if (freq > stopFreq) {
    return scalar(out, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const unsigned int matchRare = (*rare) + delta_rare; /* nextRare */
    const unsigned int outRare = (*rare) + diagterm_rare; /* nextRare */

    const vec Match = _mm256_set1_epi32(matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
	goto FINISH_SCALAR;
      }
    }
    vec Q0, Q1, Q2, Q3;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] < matchRare) {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 8), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 9), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 10), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 11), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 12), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 13), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 14), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 15), Match));
      } else {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 4), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 5), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 6), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 7), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 0), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 1), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 2), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 3), Match));
      }
    } else {
      if (freq[veclen * 23 + vecmax] < matchRare) {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 8 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 9 + 16), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 10+ 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 11+ 16), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 12+ 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 13+ 16), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 14+ 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 15+ 16), Match));
      } else {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 4 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 5 + 16), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 6 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 7 + 16), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 0 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 1 + 16), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 2 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((const vec *)(freq) + 3 + 16), Match));
      }
    }

    const vec Q = _mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3));
    if (_mm256_testz_si256(Q, Q)) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR: return (out - initout) + scalar(out, rare, stopRare + rarespace - rare,
						freq, stopFreq + freqspace - freq,
						delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE2)
/* Processes 32 vectors, each vector having 2 unsigned ints => 64 unsigned ints */
static int
v3 (unsigned int *out, const unsigned int *rare, const int lenRare,
    const unsigned int *freq, const int lenFreq, int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const unsigned int *initout = out;
  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(unsigned int);
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const unsigned int *stopFreq = freq + lenFreq - freqspace;
  const unsigned int *stopRare = rare + lenRare - rarespace;
  if (freq > stopFreq) {
    return scalar(out, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);
  }
  while (freq[veclen * 31 + vecmax] < (*rare) + delta_rare) {
    freq += veclen * 32;
    if (freq > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare < stopRare; ++rare) {
    const unsigned int matchRare = (*rare) + delta_rare;  /* nextRare */
    const unsigned int outRare = (*rare) + diagterm_rare;  /* nextRare */

    const vec Match = _mm_set1_epi32(matchRare);
    while (freq[veclen * 31 + vecmax] < matchRare) { /* if no match possible */
      freq += veclen * 32;	/* advance 32 vectors */
      if (freq > stopFreq) {
        goto FINISH_SCALAR;
      }
    }
    vec Q0, Q1, Q2, Q3;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] < matchRare) {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 8), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 9), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 10), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 11), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 12), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 13), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 14), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 15), Match));
      } else {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 4), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 5), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 6), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 7), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 0), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 1), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 2), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 3), Match));
      }
    } else {
      if (freq[veclen * 23 + vecmax] < matchRare) {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 8 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 9 + 16), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 10 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 11 + 16), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 12 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 13 + 16), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 14 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 15 + 16), Match));
      } else {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 4 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 5 + 16), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 6 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 7 + 16), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 0 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 1 + 16), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 2 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 3 + 16), Match));
      }
    }

    const vec Q = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
#ifdef HAVE_SSE4_1
    if (_mm_testz_si128(Q, Q)) {
    } else {
      *out++ = outRare;
    }
#else
    if (!_mm_movemask_epi8(Q)) {
    } else {
      *out++ = outRare;
    }
#endif
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare, stopRare + rarespace - rare,
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
 * The matchOut pointer can safely be equal to the rare pointer.
 *
 * This function DOES NOT use assembly. It only relies on intrinsics.
 */

#ifdef HAVE_AVX2
static int
SIMDgalloping (unsigned int *out, const unsigned int *rare, const int lenRare,
	       const unsigned int *freq, const int lenFreq, int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const unsigned int *initout = out;
  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(unsigned int); /* 8 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const unsigned int *stopFreq = freq + lenFreq - freqspace;
  const unsigned int *stopRare = rare + lenRare - rarespace;
  if (freq > stopFreq) {
    return scalar(out, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);
  }
  for (; rare < stopRare; ++rare) {
    const unsigned int matchRare = (*rare) + delta_rare;  /* nextRare; */
    const unsigned int outRare = (*rare) + diagterm_rare;  /* nextRare; */
    const vec Match = _mm256_set1_epi32(matchRare);

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
	/* const unsigned int mid = (lower + offset) / 2; */
	const int mid = lower + ((offset - lower) / 2);
	if (freq[veclen * mid * 32 + veclen * 31 + vecmax] < matchRare) {
	  lower = mid;
	} else {
	  offset = mid;
	}
      }
      freq += veclen * offset * 32;
    }
    vec Q0, Q1, Q2, Q3;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] < matchRare) {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 8), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 9), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15), Match));
      } else {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 4), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 5), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 6), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 7), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 0), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 1), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 2), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 3), Match));
      }
    } else {
      if (freq[veclen * 23 + vecmax] < matchRare) {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 8 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 9 + 16), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 10 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 11 + 16), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 12 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 13 + 16), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 14 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 15 + 16), Match));
      } else {
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 4 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 5 + 16), Match));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 6 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 7 + 16), Match));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 0 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 1 + 16), Match));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 2 + 16), Match),
			     _mm256_cmpeq_epi32(_mm256_loadu_si256((vec *) freq + 3 + 16), Match));
      }

    }

    const vec Q = _mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3));
    if (_mm256_testz_si256(Q, Q)) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare, stopRare + rarespace - rare,
				  freq,	stopFreq + freqspace - freq,
				  delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE2)
static int
SIMDgalloping (unsigned int *out, const unsigned int *rare, const int lenRare,
	       const unsigned int *freq, const int lenFreq,
	       int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const unsigned int *initout = out;
  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(unsigned int); /* 4 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const unsigned int *stopFreq = freq + lenFreq - freqspace;
  const unsigned int *stopRare = rare + lenRare - rarespace;
  if (freq > stopFreq) {
    return scalar(out, rare, lenRare, freq, lenFreq, delta_rare, diagterm_rare);
  }
  for (; rare < stopRare; ++rare) {
    const unsigned int matchRare = (*rare) + delta_rare; /* nextRare; */
    const unsigned int outRare = (*rare) + diagterm_rare; /* nextRare; */
    const vec Match = _mm_set1_epi32(matchRare);

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
        /* const unsigned int mid = (lower + offset) / 2; */
        const int mid = lower + (offset - lower) / 2;
        if (freq[veclen * mid * 32 + veclen * 31 + vecmax] < matchRare) {
          lower = mid;
	} else {
          offset = mid;
	}
      }
      freq += veclen * offset * 32;
    }
    vec Q0, Q1, Q2, Q3;
    if (freq[veclen * 15 + vecmax] >= matchRare) {
      if (freq[veclen * 7 + vecmax] < matchRare) {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 8), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 9), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 10), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 11), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 12), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 13), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 14), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 15), Match));
      } else {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 4), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 5), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 6), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 7), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 0), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 1), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 2), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 3), Match));
      }
    } else {
      if (freq[veclen * 23 + vecmax] < matchRare) {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 8 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 9 + 16), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 10 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 11 + 16), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 12 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 13 + 16), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 14 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 15 + 16), Match));
      } else {
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 4 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 5 + 16), Match));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 6 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 7 + 16), Match));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 0 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 1 + 16), Match));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 2 + 16), Match),
			  _mm_cmpeq_epi32(_mm_loadu_si128((const vec *)(freq) + 3 + 16), Match));
      }
    }

    const vec Q = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
#ifdef __SSE4_1__
    if (_mm_testz_si128(Q, Q)) {
    } else {
      *out++ = outRare;
    }
#else
    if (!_mm_movemask_epi8(Q)) {
    } else {
      *out++ = outRare;
    }
#endif
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare, stopRare + rarespace - rare,
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
#ifdef HAVE_SSE2
unsigned int *
Intersect_small (int *ndiagonals,
		 const unsigned int *set1, const int length1, int diagterm1,
		 const unsigned int *set2, const int length2, int diagterm2,
		 bool alignp) {
  unsigned int *diagonals;
#ifdef DEBUG15
  unsigned int *nonsimd_diagonals;
  int nonsimd_ndiagonals, i;
#endif

  debug(printf("Entered Intersect_small with %d and %d items\n",length1,length2));

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
    *ndiagonals = 0;
    return (unsigned int *) NULL;

  } else if ((1000 * length1 <= length2) || (1000 * length2 <= length1)) {
    debug(printf("Using SIMDgalloping method\n"));
    if (length1 <= length2) {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length1 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length1 * sizeof(unsigned int));
      }
      *ndiagonals = SIMDgalloping(diagonals, set1, length1, set2, length2,
				  /*delta_rare*/(diagterm1 - diagterm2),
				  /*diagterm_rare*/diagterm1);
    } else {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length2 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length2 * sizeof(unsigned int));
      }
      *ndiagonals = SIMDgalloping(diagonals, set2, length2, set1, length1,
				  /*delta_rare*/(diagterm2 - diagterm1),
				  /*diagterm_rare*/diagterm2);
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_old(&nonsimd_ndiagonals,
					    set1,length1,diagterm1,
					    set2,length2,diagterm2);
    if (*ndiagonals != nonsimd_ndiagonals) {
      abort();
    } else {
      for (i = 0; i < *ndiagonals; i++) {
	if (diagonals[i] != nonsimd_diagonals[i]) {
	  abort();
	}
      }
    }

#endif
    return diagonals;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    debug(printf("Using v3 method\n"));
    if (length1 <= length2) {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length1 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length1 * sizeof(unsigned int));
      }
      *ndiagonals = v3(diagonals, set1, length1, set2, length2,
		       /*delta_rare*/(diagterm1 - diagterm2),
		       /*diagterm_rare*/diagterm1);
    } else {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length2 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length2 * sizeof(unsigned int));
      }
      *ndiagonals = v3(diagonals, set2, length2, set1, length1,
		       /*delta_rare*/(diagterm2 - diagterm1),
		       /*diagterm_rare*/diagterm2);
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_old(&nonsimd_ndiagonals,
					    set1,length1,diagterm1,
					    set2,length2,diagterm2);
  if (*ndiagonals != nonsimd_ndiagonals) {
    abort();
  } else {
    for (i = 0; i < *ndiagonals; i++) {
      if (diagonals[i] != nonsimd_diagonals[i]) {
	abort();
      }
    }
  }

#endif
    return diagonals;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length1 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length1 * sizeof(unsigned int));
      }
      *ndiagonals = v1(diagonals, set1, length1, set2, length2,
		       /*delta_rare*/(diagterm1 - diagterm2),
		       /*diagterm_rare*/diagterm1);
    } else {
      if (alignp == true) {
	MALLOC_ALIGN(diagonals,length2 * sizeof(unsigned int));
      } else {
	diagonals = (unsigned int *) MALLOC(length2 * sizeof(unsigned int));
      }
      *ndiagonals = v1(diagonals, set2, length2, set1, length1,
		       /*delta_rare*/(diagterm2 - diagterm1),
		       /*diagterm_rare*/diagterm2);
				      
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_old(&nonsimd_ndiagonals,
					    set1,length1,diagterm1,
					    set2,length2,diagterm2);
  if (*ndiagonals != nonsimd_ndiagonals) {
    abort();
  } else {
    for (i = 0; i < *ndiagonals; i++) {
      if (diagonals[i] != nonsimd_diagonals[i]) {
	abort();
      }
    }
  }

#endif
    return diagonals;
  }
}

#else

unsigned int *
Intersect_small (int *ndiagonals,
		 const unsigned int *set1, const int length1, int diagterm1,
		 const unsigned int *set2, const int length2, int diagterm2,
		 bool alignp) {
  unsigned int *diagonals;

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagonals = 0;
    return (unsigned int *) NULL;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    if (length1 <= length2) {
      if (alignp == true) {
      } else {
	diagonals = (unsigned int *) MALLOC(length1 * sizeof(unsigned int));
      }
      *ndiagonals = galloping_intersection(diagonals, set1, length1, set2, length2,
					   /*delta_rare*/(diagterm1 - diagterm2),
					   /*diagterm_rare*/diagterm1);
    } else {
      if (alignp == true) {
      } else {
      diagonals = (unsigned int *) MALLOC(length2 * sizeof(unsigned int));
      }
      *ndiagonals = galloping_intersection(diagonals, set2, length2, set1, length1,
					   /*delta_rare*/(diagterm2 - diagterm1),
					   /*diagterm_rare*/diagterm2);
    }
    return diagonals;

  } else {
    if (length1 <= length2) {
      if (alignp == true) {
      } else {
	diagonals = (unsigned int *) MALLOC(length1 * sizeof(unsigned int));
      }
      *ndiagonals = scalar(diagonals, set1, length1, set2, length2,
			   /*deltaA*/(diagterm1 - diagterm2),
			   /*diagtermA*/diagterm1);
    } else {
      if (alignp == true) {
      } else {
      diagonals = (unsigned int *) MALLOC(length2 * sizeof(unsigned int));
      }
      *ndiagonals = scalar(diagonals, set2, length2, set1, length1,
			   /*deltaA*/(diagterm2 - diagterm1),
			   /*diagtermA*/diagterm2);
    }
    return diagonals;
  }
}
#endif



#ifdef DEBUG15
static int
binary_search_uint4 (int lowi, int highi, unsigned int *positions, unsigned int goal) {
  int middlei;

  debug10(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (goal < positions[middlei]) {
      highi = middlei;
    } else if (goal > positions[middlei]) {
      lowi = middlei + 1;
    } else {
      debug10(printf("binary search returns %d\n",middlei));
      return middlei;
    }
  }

  debug10(printf("binary search returns %d\n",highi));
  return highi;
}


/* Old implementation */
/* LARGE_GENOMES needs this to handle transcriptomes */
unsigned int *
Intersect_exact_old (int *ndiagonals,
		     const unsigned int *positionsa, const int npositionsa, int diagterma,
		     const unsigned int *positionsb, const int npositionsb, int diagtermb) {
  unsigned int *diagonals, local_goal, last_diagonal, this_diagonal;
  unsigned int diagterm, delta;
  const unsigned int *positions0, *positions1;
  int npositions0, npositions1, j;


  *ndiagonals = 0;
  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_exact: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    return (unsigned int *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0 = positionsa;
    positions1 = positionsb;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (unsigned int) diagtermb;	/* local_goal based on larger list */
    delta = (unsigned int) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0 = positionsb;
    positions1 = positionsa;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (unsigned int) diagterma;	/* local_goal based on larger list */
    delta = (unsigned int) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  debug(printf("Intersect_small with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));

  diagonals = (unsigned int *) MALLOC(npositions0 * sizeof(unsigned int));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && (*positions0) < (Univcoord_T) -(diagtermb - diagterma)) {
      ++positions0;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && (*positions1) < (Univcoord_T) -diagterma) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = (*positions0) + delta;
    debug(printf("intersection list 0: %d:%u => local_goal %u\n",npositions0,*positions0,local_goal));
    if (npositions1 > 0 && *positions1 < local_goal) {
      j = 1;
      while (j < npositions1 && positions1[j] < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_uint4(j >> 1,npositions1,positions1,local_goal);
      } else {
	j = binary_search_uint4(j >> 1,j,positions1,local_goal);
      }
      positions1 += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return diagonals;

    } else if ((*positions1) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%u  found\n",npositions1,*positions1));
      if (*ndiagonals == 0) {
	last_diagonal = diagonals[(*ndiagonals)++] = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	last_diagonal = diagonals[(*ndiagonals)++] = this_diagonal;
      }
      ++positions1;
      --npositions1;
    }

    ++positions0;
    --npositions0;
  }
  debug(printf("\n"));

  return diagonals;
}
#endif



