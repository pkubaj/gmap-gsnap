static char rcsid[] = "$Id: 48dcfea3300509a07262e332763ec638e937c014 $";
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

/* Minimum requirement for SIMD is SSE4.1 for _mm_cvtepi8_epi64 and _mm_cmpeq_epi64 */

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "intersect-large.h"

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


#define GETPOS(high,low) (((UINT8) high << 32) + low)
#define GETPOS_I(high,low,i) (((UINT8) high[i] << 32) + low[i])


#ifndef HAVE_SSE4_1
static int
__frogadvanceUntil(const unsigned char *array_high, const UINT4 *array_low, const int pos,
		   const int length, const UINT8 min) {
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

static int
galloping_intersection (UINT8 *out, const unsigned char *smallset_high, const UINT4 *smallset_low, const int smalllength,
			const unsigned char *largeset_high, const UINT4 *largeset_low, const int largelength,
			int delta_small, int diagterm_small) {
  if (smalllength == 0) return 0;
  const UINT8 *initout = out;
  int k1 = 0, k2 = 0;
  while (1) {
    if (GETPOS_I(largeset_high,largeset_low,k1) < GETPOS_I(smallset_high,smallset_low,k2) + delta_small) {
      k1 = __frogadvanceUntil(largeset_high, largeset_low, k1, largelength, GETPOS_I(smallset_high,smallset_low,k2) + delta_small);
      if (k1 == largelength) {
        break;
      }
    }
  midpoint:
    if (GETPOS_I(smallset_high,smallset_low,k2) + delta_small < GETPOS_I(largeset_high,largeset_low,k1)) {
      ++k2;
      if (k2 == smalllength) {
        break;
      }
    } else {
      *out++ = GETPOS_I(smallset_high,smallset_low,k2) + diagterm_small;
      ++k2;
      if (k2 == smalllength) {
        break;
      }
      k1 = __frogadvanceUntil(largeset_high, largeset_low, k1, largelength, GETPOS_I(smallset_high,smallset_low,k2) + delta_small);
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
scalar (UINT8 *out, const unsigned char *A_high, const UINT4 *A_low, const int lenA,
	const unsigned char *B_high, const UINT4 *B_low, const int lenB, int deltaA, int diagtermA) {
  const UINT8 *initout = out;
  if (lenA == 0 || lenB == 0) return 0;

  const UINT4 *endA = A_low + lenA;
  const UINT4 *endB = B_low + lenB;

  while (1) {
    while (GETPOS(*A_high,*A_low) + deltaA < GETPOS(*B_high,*B_low)) {
    SKIP_FIRST_COMPARE:
      ++A_high;
      if (++A_low == endA) {
        return (out - initout);
      }
    }
    while (GETPOS(*A_high,*A_low) + deltaA > GETPOS(*B_high,*B_low)) {
      ++B_high;
      if (++B_low == endB) {
        return (out - initout);
      }
    }
    if (GETPOS(*A_high,*A_low) + deltaA == GETPOS(*B_high,*B_low)) {
      *out++ = GETPOS(*A_high,*A_low) + diagtermA;
      ++A_high;
      if (++A_low == endA) {
        return (out - initout);
      }
      ++B_high;
      if (++B_low == endB) {
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

/* Need SSE4.1 to have _mm_cvtepi8_epi64 */
#if defined(HAVE_SSE4_1)
static int
v1 (UINT8 *matchOut, const unsigned char *rare_high, const UINT4 *rare_low, int lenRare,
    const unsigned char *freq_high, const UINT4 *freq_low, int lenFreq,
    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const UINT8 *matchOrig = matchOut;

  UINT8 valRare, outRare;
  UINT8 maxProbe;		/* was uint64_t */
  UINT8 maxFreq;		/* was uint64_t */
  const int kRareSpace = 0;	/* was uint64_t */

  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 8 64-mers in 512 bits */

#ifdef HAVE_AVX512
  __m512i Rare;
  __m128i H;
  __m512i F;
  __mmask8 cmp_mask;
#elif defined(HAVE_AVX2)
  __m256i Rare;
  __m128i H;
  __m256i F0, F1;
  __m256i CMP0, CMP1;
#else
  __m128i Rare;
  __m128i H;
  __m128i F0, F1, F2, F3;
  __m128i CMP0, CMP1;
#endif

  const UINT4 *stopFreq = &freq_low[lenFreq] - kFreqSpace;
  const UINT4 *stopRare = &rare_low[lenRare] - kRareSpace;


  if (COMPILER_RARELY((rare_low >= stopRare) || (freq_low >= stopFreq)))
    goto FINISH_SCALAR;

  valRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* For comparison */
  outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; /* For output */
#ifdef HAVE_AVX512
  Rare = _mm512_set1_epi64(valRare);
#elif defined(HAVE_AVX2)
  Rare = _mm256_set1_epi64x(valRare);
#else
  Rare = _mm_set1_epi64((__m64) valRare);
#endif

  maxFreq = GETPOS_I(freq_high,freq_low,2 * 4 - 1);
#ifdef HAVE_AVX512
  /* Reads in 16 chars from freq_high, but uses only 8 */
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
		       _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
#elif defined(HAVE_AVX2)
  /* Reads in 16 chars from freq_high, but uses only 8 */
  H = _mm_lddqu_si128((const __m128i *)(freq_high));
  F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,4);
  F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			_mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
#elif defined(HAVE_SSE4_1)
  /* Reads in 16 chars from freq_high, but uses only 8 */
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
#else
  H = _mm_loadu_si128((const __m128i *)(freq_high));
  F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_loadu_si128((const __m128i *)(freq_low))));
  H = _mm_bslli_si128(H,2);
  F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_loadu_si128((const __m128i *)(freq_low + 2))));
  H = _mm_bslli_si128(H,2);
  F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_loadu_si128((const __m128i *)(freq_low + 4))));
  H = _mm_bslli_si128(H,2);
  F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		     _mm_cvtepi32_epi64(_mm_loadu_si128((const __m128i *)(freq_low + 6))));
#endif

  if (COMPILER_RARELY(maxFreq < valRare))
    goto ADVANCE_FREQ;

 ADVANCE_RARE:
  do {
    *matchOut = outRare;
    rare_high += 1;
    rare_low += 1;
    if (COMPILER_RARELY(rare_low >= stopRare)) {
      rare_high -= 1;
      rare_low -= 1;
      goto FINISH_SCALAR;
    }

    valRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* for next iteration */
    outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; /* for next iteration */
#ifdef HAVE_AVX512
    cmp_mask = _mm512_cmpeq_epi64_mask(F,Rare);
    Rare = _mm512_set1_epi64(valRare);
    if (cmp_mask != 0) {
      matchOut ++;
    }

#elif defined(HAVE_AVX2)
    CMP0 = _mm256_cmpeq_epi64(F0, Rare);
    CMP1 = _mm256_cmpeq_epi64(F1, Rare);
    Rare = _mm256_set1_epi64x(valRare);
    CMP0 = _mm256_or_si256(CMP0,CMP1);
    if(_mm256_testz_si256(CMP0,CMP0) == 0) {
      matchOut ++;
    }

#elif HAVE_SSE4_1
    CMP0 = _mm_or_si128(_mm_cmpeq_epi64(F0, Rare),
			_mm_cmpeq_epi64(F1, Rare));
    CMP1 = _mm_or_si128(_mm_cmpeq_epi64(F2, Rare),
			_mm_cmpeq_epi64(F3, Rare));
    Rare = _mm_set1_epi64((__m64) valRare);
    CMP0 = _mm_or_si128(CMP0, CMP1);
    if (_mm_testz_si128(CMP0, CMP0) == 0) {
      matchOut++;
    }
#else
    CMP0 = _mm_or_si128(_mm_cmpeq_epi64(F0, Rare),
			_mm_cmpeq_epi64(F1, Rare));
    CMP1 = _mm_or_si128(_mm_cmpeq_epi64(F2, Rare),
			_mm_cmpeq_epi64(F3, Rare));
    Rare = _mm_set1_epi64((__m64) valRare);
    CMP0 = _mm_or_si128(CMP0, CMP1);
    if (_mm_movemask_epi8(CMP0)) {
      matchOut++;
    }
#endif

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
  count = matchOut - matchOrig;

  lenFreq = stopFreq + kFreqSpace - freq_low;
  lenRare = stopRare + kRareSpace - rare_low;

  /* Was match_scalar */
  int tail = scalar(matchOut, rare_high, rare_low, lenRare,
		    freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);

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
/* Processes 32 vectors, each vector having 8 chars and 8 unsigned ints => 256 chars and 256 unsigned ints */
static int
v3 (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq, int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const UINT8 *initout = out;
  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 8 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;

  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq_high += veclen * 32;	/* advance 32 vectors */
    freq_low += veclen * 32;	/* advance 32 vectors */
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* nextRare */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; /* nextRare */

    const vec Match = _mm512_set1_epi64(matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
	goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1;
    __mmask8 Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
	
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 32 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 12*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9 + 16)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 13*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11+ 16)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 14*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13+ 16)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 15*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15+ 16)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 10*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5 + 16)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 11*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7 + 16)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 8*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1 + 16)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 9*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3 + 16)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
      }
    }
    const __mmask8 Q = (Q0 | Q1) | (Q2 | Q3);
    if (Q == 0) {
    } else {
      *out++ = outRare;
    }
  }

  FINISH_SCALAR:
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low, stopFreq + freqspace - freq_low,
				  delta_rare, diagterm_rare);
}

#elif defined(HAVE_AVX2)
/* Processes 32 vectors, each vector having 4 chars and 4 unsigned ints => 128 chars and 128 unsigned ints */
static int
v3 (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
    int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const UINT8 *initout = out;
  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 4 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq_high += veclen * 32;	/* advance 32 vectors */
    freq_low += veclen * 32;	/* advance 32 vectors */
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  /* nextRare */

    const vec Match = _mm256_set1_epi64x(matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
        goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1;
    vec Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));
      }
    }

    const vec Q = _mm256_or_si256(_mm256_or_si256(Q0, Q1), _mm256_or_si256(Q2, Q3));
    if (_mm256_testz_si256(Q, Q)) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low, stopFreq + freqspace - freq_low,
				  delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE4_1)
/* Processes 32 vectors, each vector having 2 chars and 2 unsigned ints => 64 chars and 64 unsigned ints */
static int
v3 (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
    const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
    int delta_rare, int diagterm_rare) {

  assert(lenRare <= lenFreq);
  if (lenFreq == 0 || lenRare == 0) return 0;

  const UINT8 *initout = out;

  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 2 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < GETPOS(*rare_high,*rare_low) + delta_rare) {
    freq_high += veclen * 32;	/* advance 32 vectors */
    freq_low += veclen * 32;	/* advance 32 vectors */
    if (freq_low > stopFreq) {
      goto FINISH_SCALAR;
    }
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  /* nextRare */

    const vec Match = _mm_set1_epi64((__m64) matchRare);
    while (GETPOS_I(freq_high,freq_low,veclen * 31 + vecmax) < matchRare) { /* if no match possible */
      freq_high += veclen * 32;	/* advance 32 vectors */
      freq_low += veclen * 32;	/* advance 32 vectors */
      if (freq_low > stopFreq) {
        goto FINISH_SCALAR;
      }
    }

    __m128i H;
    vec F0, F1;
    vec Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));
      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));
      }
    }
    const vec Q = _mm_or_si128(_mm_or_si128(Q0, Q1), _mm_or_si128(Q2, Q3));
    if (_mm_testz_si128(Q, Q)) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low, stopFreq + freqspace - freq_low,
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

#ifdef HAVE_AVX512
static int
SIMDgalloping (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	       int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const UINT8 *initout = out;
  typedef __m512i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 8 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare; */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  /* nextRare; */
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
	/* const unsigned int mid = (lower + offset) / 2; */
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
    vec F0, F1;
    __mmask8 Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
	
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 12*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 8 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 9 + 16)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 13*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 10+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 11+ 16)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 14*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 12+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 13+ 16)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 15*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 14+ 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 15+ 16)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 10*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 4 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 5 + 16)));
	Q0 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 11*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 6 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 7 + 16)));
	Q1 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 8*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 0 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 1 + 16)));
	Q2 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 9*16)); /* Reads 16 chars */
	F0 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 2 + 16)));
	F1 = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low) + 3 + 16)));
	Q3 = _mm512_cmpeq_epi64_mask(F0, Match) | _mm512_cmpeq_epi64_mask(F1, Match);
      }

    }
    const __mmask8 Q = (Q0 | Q1) | (Q2 | Q3);
    if (Q == 0) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low, stopFreq + freqspace - freq_low,
				  delta_rare, diagterm_rare);
}

#elif defined(HAVE_AVX2)
static int
SIMDgalloping (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	       int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const UINT8 *initout = out;
  typedef __m256i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 4 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare; */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  /* nextRare; */
    const vec Match = _mm256_set1_epi64x(matchRare);

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
	/* const unsigned int mid = (lower + offset) / 2; */
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
    vec F0, F1;
    vec Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

      } else {
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));
      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 6*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 7*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 5*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
	Q0 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));
	Q1 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	H = _mm_lddqu_si128((const __m128i *)(freq_high + 4*16)); /* Reads 16 chars */
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
	Q2 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));

	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
	Q3 = _mm256_or_si256(_mm256_cmpeq_epi64(F0, Match), _mm256_cmpeq_epi64(F1, Match));
      }

    }
    const vec Q = _mm256_or_si256(_mm256_or_si256(Q0, Q1),_mm256_or_si256(Q2, Q3));
    if (_mm256_testz_si256(Q, Q)) {
    } else {
      *out++ = outRare;
    }
  }

 FINISH_SCALAR:
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low,stopFreq + freqspace - freq_low,
				  delta_rare, diagterm_rare);
}

#elif defined(HAVE_SSE4_1)
static int
SIMDgalloping (UINT8 *out, const unsigned char *rare_high, const UINT4 *rare_low, const int lenRare,
	       const unsigned char *freq_high, const UINT4 *freq_low, const int lenFreq,
	       int delta_rare, int diagterm_rare) {
  if (lenFreq == 0 || lenRare == 0) return 0;
  assert(lenRare <= lenFreq);
  const UINT8 *initout = out;
  typedef __m128i vec;
  const int veclen = sizeof(vec) / sizeof(UINT8); /* 2 */
  const int vecmax = veclen - 1;
  const int freqspace = 32 * veclen;
  const int rarespace = 1;

  const UINT4 *stopFreq = freq_low + lenFreq - freqspace;
  const UINT4 *stopRare = rare_low + lenRare - rarespace;
  if (freq_low > stopFreq) {
    return scalar(out, rare_high, rare_low, lenRare,
		  freq_high, freq_low, lenFreq, delta_rare, diagterm_rare);
  }
  for (; rare_low < stopRare; ++rare_high, ++rare_low) {
    const UINT8 matchRare = GETPOS(*rare_high,*rare_low) + delta_rare;  /* nextRare; */
    const UINT8 outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;  /* nextRare; */
    const vec Match = _mm_set1_epi64((__m64) matchRare);

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
	/* const unsigned int mid = (lower + offset) / 2; */
	const unsigned int mid = lower + ((offset - lower) / 2);
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
    vec F0, F1;
    vec Q0, Q1, Q2, Q3;
    if (GETPOS_I(freq_high,freq_low,veclen * 15 + vecmax) >= matchRare) {
      if (GETPOS_I(freq_high,freq_low,veclen * 7 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 1*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 0*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));
      }
    } else {
      if (GETPOS_I(freq_high,freq_low,veclen * 23 + vecmax) < matchRare) {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 3*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 8 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 9 + 16)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 10+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 11+ 16)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 12+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 13+ 16)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 14+ 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 15+ 16)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

      } else {
	H = _mm_lddqu_si128((const __m128i *)(freq_high + 2*16)); /* Reads 16 chars */
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,8)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 4 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,10)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 5 + 16)));
        Q0 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,12)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 6 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,14)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 7 + 16)));
        Q1 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 0 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,2)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 1 + 16)));
        Q2 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));

	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,4)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 2 + 16)));
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(_mm_bslli_si128(H,6)),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low) + 3 + 16)));
        Q3 = _mm_or_si128(_mm_cmpeq_epi32(F0, Match), _mm_cmpeq_epi32(F1, Match));
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
  return (out - initout) + scalar(out, rare_high, rare_low, stopRare + rarespace - rare_low,
				  freq_high, freq_low, stopFreq + freqspace - freq_low,
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
UINT8 *
Intersect_large (int *ndiagonals,
		 const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		 const unsigned char *set2_high, const UINT4 *set2_low, const int length2, int diagterm2) {
  UINT8 *diagonals;
#ifdef DEBUG15
  UINT8 *nonsimd_diagonals;
  int nonsimd_ndiagonals, i;
#endif

  debug(printf("Entered Intersect_large with %d and %d items\n",length1,length2));

#ifdef DEBUG
  int i;
  for (i = 0; i < length1; i++) {
    printf("%d %u\n",set1_high[i],set1_low[i]);
  }
  printf("\n");
  for (i = 0; i < length2; i++) {
    printf("%d %u\n",set2_high[i],set2_low[i]);
  }
  printf("\n");
#endif    

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagonals = 0;
    return (UINT8 *) NULL;

  } else if ((1000 * length1 <= length2) || (1000 * length2 <= length1)) {
    debug(printf("Using SIMDgalloping method\n"));
    if (length1 <= length2) {
      diagonals = (UINT8 *) MALLOC(length1 * sizeof(UINT8));
      *ndiagonals = SIMDgalloping(diagonals, set1_high, set1_low, length1,
				  set2_high, set2_low, length2,
				  /*delta_rare*/(diagterm1 - diagterm2),
				  /*diagterm_rare*/diagterm1);
    } else {
      diagonals = (UINT8 *) MALLOC(length2 * sizeof(UINT8));
      *ndiagonals = SIMDgalloping(diagonals, set2_high, set2_low, length2,
				  set1_high, set1_low, length1,
				  /*delta_rare*/(diagterm2 - diagterm1),
				  /*diagterm_rare*/diagterm2);
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_large_old(&nonsimd_ndiagonals,
						  set1_high,set1_low,length1,diagterm1,
						  set2_high,set2_low,length2,diagterm2);
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
      diagonals = (UINT8 *) MALLOC(length1 * sizeof(UINT8));
      *ndiagonals = v3(diagonals, set1_high, set1_low, length1,
		       set2_high, set2_low, length2,
		       /*delta_rare*/(diagterm1 - diagterm2),
		       /*diagterm_rare*/diagterm1);
    } else {
      diagonals = (UINT8 *) MALLOC(length2 * sizeof(UINT8));
      *ndiagonals = v3(diagonals, set2_high, set2_low, length2,
		       set1_high, set1_low, length1,
		       /*delta_rare*/(diagterm2 - diagterm1),
		       /*diagterm_rare*/diagterm2);
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_large_old(&nonsimd_ndiagonals,
						  set1_high,set1_low,length1,diagterm1,
						  set2_high,set2_low,,length2,diagterm2);
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
      diagonals = (UINT8 *) MALLOC(length1 * sizeof(UINT8));
      *ndiagonals = v1(diagonals, set1_high, set1_low, length1,
		       set2_high, set2_low, length2,
		       /*delta_rare*/(diagterm1 - diagterm2),
		       /*diagterm_rare*/diagterm1);
    } else {
      diagonals = (UINT8 *) MALLOC(length2 * sizeof(UINT8));
      *ndiagonals = v1(diagonals, set2_high, set2_low, length2,
		       set1_high, set1_low, length1,
		       /*delta_rare*/(diagterm2 - diagterm1),
		       /*diagterm_rare*/diagterm2);
				      
    }
#ifdef DEBUG15
    nonsimd_diagonals = Intersect_exact_large_old(&nonsimd_ndiagonals,
						  set1_high,set1_low,length1,diagterm1,
						  set2_high,set2_low,length2,diagterm2);
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

UINT8 *
Intersect_large (int *ndiagonals,
		 const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		 const unsigned char *set2_high, const UINT4 *set2_low, const int length2, int diagterm2) {
  UINT8 *diagonals;

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagonals = 0;
    return (UINT8 *) NULL;

  } else if ((50 * length1 <= length2) || (50 * length2 <= length1)) {
    if (length1 <= length2) {
      diagonals = (UINT8 *) MALLOC(length1 * sizeof(UINT8));
      *ndiagonals = galloping_intersection(diagonals, set1_high, set1_low, length1,
					   set2_high, set2_low, length2,
					   /*delta_rare*/(diagterm1 - diagterm2),
					   /*diagterm_rare*/diagterm1);
    } else {
      diagonals = (UINT8 *) MALLOC(length2 * sizeof(UINT8));
      *ndiagonals = galloping_intersection(diagonals, set2_high, set2_low, length2,
					   set1_high, set1_low, length1,
					   /*delta_rare*/(diagterm2 - diagterm1),
					   /*diagterm_rare*/diagterm2);
    }
    return diagonals;

  } else {
    if (length1 <= length2) {
      diagonals = (UINT8 *) MALLOC(length1 * sizeof(UINT8));
      *ndiagonals = scalar(diagonals, set1_high, set1_low, length1,
			   set2_high, set2_low, length2,
			   /*deltaA*/(diagterm1 - diagterm2),
			   /*diagtermA*/diagterm1);
    } else {
      diagonals = (UINT8 *) MALLOC(length2 * sizeof(UINT8));
      *ndiagonals = scalar(diagonals, set2_high, set2_low, length2,
			   set1_high, set1_low, length1,
			   /*deltaA*/(diagterm2 - diagterm1),
			   /*diagtermA*/diagterm2);
    }
    return diagonals;
  }
}
#endif



#ifdef DEBUG15
static int
binary_search_large (int lowi, int highi, unsigned char *positions_high, UINT4 *positions_low, UINT8 goal) {
  int middlei;
  UINT8 position;

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    position = ((UINT8) positions_high[middlei] << 32) + positions_low[middlei];
    debug10(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,(positions_high[lowi] << 32) + positions_low[lowi],
		   middlei,position,
		   highi,(positions_high[highi] << 32) + positions_low[highi],goal));
    if (goal < position) {
      highi = middlei;
    } else if (goal > position) {
      lowi = middlei + 1;
    } else {
      return middlei;
    }
  }

  return highi;
}


UINT8 *
Intersect_exact_large_old (int *ndiagonals,
			   unsigned char *positionsa_high, UINT4 *positionsa_low,
			   int npositionsa, int diagterma,
			   unsigned char *positionsb_high, UINT4 *positionsb_low,
			   int npositionsb, int diagtermb) {
  UINT8 *diagonals, local_goal, last_diagonal, this_diagonal;
  UINT8 diagterm, delta;
  unsigned char *positions0_high, *positions1_high;
  UINT4 *positions0_low, *positions1_low;
  int npositions0, npositions1, j;


  *ndiagonals = 0;
  if (npositionsa == 0 || npositionsb == 0) {
    debug(printf("Intersect_exact_large: intersection is null because npositionsa %d or npositionsb %d is zero\n",npositionsa,npositionsb));
    return (UINT8 *) NULL;

  } else if (npositionsa < npositionsb) {
    positions0_high = positionsa_high;
    positions0_low = positionsa_low;
    positions1_high = positionsb_high;
    positions1_low = positionsb_low;
    npositions0 = npositionsa;
    npositions1 = npositionsb;

    diagterm = (UINT8) diagtermb;	/* local_goal based on larger list */
    delta = (UINT8) (diagterma - diagtermb); /* list0 + (diagterm0 - diagterm1) = list1 */

  } else {
    positions0_high = positionsb_high;
    positions0_low = positionsb_low;
    positions1_high = positionsa_high;
    positions1_low = positionsa_low;
    npositions0 = npositionsb;
    npositions1 = npositionsa;

    diagterm = (UINT8) diagterma;	/* local_goal based on larger list */
    delta = (UINT8) (diagtermb - diagterma); /* list0 + (diagterm0 - diagterm1) = list1 */
  }

  debug(printf("Intersect_exact_large with %d positions <= %d positions.  diagterm %d\n",
	       npositions0,npositions1,(int) diagterm));

  diagonals = (UINT8 *) MALLOC(npositions0 * sizeof(UINT8));

#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagtermb - diagterma < 0) {
    while (npositions0 > 0 && GETPOS(*positions0_high,*positions0_low) < (UINT8) -(diagtermb - diagterma)) {
      ++positions0_high;
      ++positions0_low;
      --npositions0;
    }
  }

  if (diagterma < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (UINT8) -diagterma) {
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }
  }
#endif

  while (npositions0 > 0) {
    local_goal = GETPOS(*positions0_high,*positions0_low) + delta;
    debug(printf("intersection list 0: %d:%llu => local_goal %llu\n",
		 npositions0,GETPOS(*positions0_high,*positions0_low),local_goal));
    if (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < local_goal) {
      j = 1;
      while (j < npositions1 && GETPOS(positions1_high[j],positions1_low[j]) < local_goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions1) {
	j = binary_search_large(j >> 1,npositions1,positions1_high,positions1_low,local_goal);
      } else {
	j = binary_search_large(j >> 1,j,positions1_high,positions1_low,local_goal);
      }
      positions1_high += j;
      positions1_low += j;
      npositions1 -= j;
    }

    if (npositions1 <= 0) {
      return diagonals;

    } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
      /* Found local goal.  Save and advance */
      debug(printf("    intersection list 1: %d:%llu  found\n",
		   npositions1,GETPOS(*positions1_high,*positions1_low)));
      if (*ndiagonals == 0) {
	last_diagonal = diagonals[(*ndiagonals)++] = local_goal + diagterm;
      } else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	last_diagonal = diagonals[(*ndiagonals)++] = this_diagonal;
      }
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }

    ++positions0_high;
    ++positions0_low;
    --npositions0;
  }
  debug(printf("\n"));

  return diagonals;
}
#endif



