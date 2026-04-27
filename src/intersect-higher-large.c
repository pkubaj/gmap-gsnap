static char rcsid[] = "$Id: 1ffde2f300921186398c0cfb788042d45b5877b1 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* Code is modified from the SIMDCompressionAndIntersection package,
   by Daniel Lemire, Leonid Boytsov, and Nathan Kurz, file
   intersection.cpp.  */

/* Outputs only a single diagonal */

/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 */

/* Minimum requirement for SIMD is SSE4.2 for _mm_cmpgt_epi64 */

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "bool.h"
#include "intersect-higher-large.h"

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


#define EPI16_MAX 0x8000	     /* 4 nibbles * 4 bits/nibble */
#define EPI32_MAX 0x80000000	     /* 8 nibbles * 4 bits/nibble*/
#define EPI64_MAX 0x8000000000000000 /* 16 nibbles * 4 bits/nibble */
/*                  0123456701234567 */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Intersect_approx_higher_old */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Check SIMD against non-SIMD */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif

/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif


#define GETPOS(high,low) (((UINT8) high << 32) + low)
#define GETPOS_I(high,low,i) (((UINT8) high[i] << 32) + low[i])


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


/* Fast scalar scheme designed by N. Kurz. */
/* Delta is required to be added to A to become equivalent with B */
/* Desired result should be A + diagtermA (== B + diagtermB) */

/* A_first_p is false */
static Univcoord_T *
scalar_anchorA (Univcoord_T *out, Univcoord_T *diagonals,
		const Univcoord_T *A, const Univcoord_T *endA,
		const unsigned char *B_high, const UINT4 *B_low, const UINT4 *endB, int deltaA,
		int diagtermA, int diagtermB,
		const Univcoord_T below_slop, const Univcoord_T above_slop) {

  const unsigned char *B_high_save;
  const UINT4 *B_low_save;
  Univcoord_T B_value;

  debug(printf("Entered scalar with %lu and %lu entries\n",endA - A,endB - B_low));

  while (A < endA) {
    /* Advance B until (A - below_slop) */
    while (B_low < endB && GETPOS(*B_high,*B_low) + below_slop < (*A) + deltaA) {
      B_high++;
      B_low++;
    }

    B_high_save = B_high;
    B_low_save = B_low;
    debug(printf("Advancing B from %u up to %u\n",GETPOS(*B_high,*B_low),(*A) + deltaA + above_slop));
    while (B_low < endB && (B_value = GETPOS(*B_high,*B_low)) <= (*A) + deltaA + above_slop) {
      debug(printf("(1) Generating output %u %u\n",(*A) + diagtermA,GETPOS(*B_high,*B_low) + diagtermB));
      if ((*A) + deltaA == B_value) {
	/* Exclude exact */
      } else if (1 || B_value + diagtermB >= (*A) + diagtermA) { /* higher */
	if (out == diagonals || B_value + diagtermB > out[-1]) {
	  *out++ = B_value + diagtermB;
	}
      }
	  
      B_high++;
      B_low++;
    }
    B_high = B_high_save;
    B_low = B_low_save;

    A++;
  }

  return out;
}


/* A_first_p is true */
static Univcoord_T *
scalar_anchorB (Univcoord_T *out, Univcoord_T *diagonals,
		const unsigned char *A_high, const UINT4 *A_low, const UINT4 *endA,
		const Univcoord_T *B, const Univcoord_T *endB, int deltaA,
		int diagtermA, int diagtermB,
		const Univcoord_T below_slop, const Univcoord_T above_slop) {

  const Univcoord_T *B_save;
  Univcoord_T A_value;

  debug(printf("Entered scalar with %lu and %lu entries\n",endA - A_low,endB - B));

  while (A_low < endA) {
    /* Advance B until (A - below_slop) */
    while (B < endB && (*B) + below_slop < GETPOS(*A_high,*A_low) + deltaA) {
      B++;
    }

    B_save = B;
    debug(printf("Advancing B from %u up to %u\n",*B,GETPOS(*A_high,*A_low) + deltaA + above_slop));
    while (B < endB && (*B) <= (A_value = GETPOS(*A_high,*A_low)) + deltaA + above_slop) {
      debug(printf("(1) Generating output %u %u\n",A_value + diagtermA,(*B) + diagtermB));
      if (A_value + deltaA == (*B)) {
	/* Exclude exact */
      } else if (1 || A_value + diagtermA >= (*B) + diagtermB) { /* higher */
	if (out == diagonals || A_value + diagtermA > out[-1]) {
	  *out++ = A_value + diagtermA;
	}
      }
	  
      B++;
    }
    B = B_save;

    A_high++;
    A_low++;
  }

  return out;
}


/* Performs scalar against a single value of A */
/* A_first_p is false */
static inline Univcoord_T *
scalar_one_set_anchorA (Univcoord_T *out, Univcoord_T *diagonals,
			const Univcoord_T *A,
			const unsigned char *B_high, const UINT4 *B_low, const UINT4 *endB,
			int deltaA, int diagtermA, int diagtermB,
			const Univcoord_T below_slop, const Univcoord_T above_slop) {

  Univcoord_T B_value;

  debug(printf("Entered scalar_one_set with %lu entries\n",endB - B_low));

  debug(printf("Advancing B from %u up to %u\n",GETPOS(*B_high,*B_low),(*A) + deltaA - below_slop));
  while (B_low < endB && GETPOS(*B_high,*B_low) + below_slop < (*A) + deltaA) {
    B_high++;
    B_low++;
  }

  while (B_low < endB && (B_value = GETPOS(*B_high,*B_low)) <= (*A) + deltaA + above_slop) {
    debug(printf("(2) Generating output %u %u\n",(*A) + diagtermA,B_value + diagtermB));
    if ((*A) + deltaA == B_value) {
      /* Exclude exact */
    } else if (1 || B_value + diagtermB >= (*A) + diagtermA) { /* higher */
      if (out == diagonals || B_value + diagtermB > out[-1]) {
	*out++ = B_value + diagtermB;
      }
    }

    B_high++;
    B_low++;
  }

  return out;
}


/* A_first_p is true */
static inline Univcoord_T *
scalar_one_set_anchorB (Univcoord_T *out, Univcoord_T *diagonals,
			const unsigned char *A_high, const UINT4 *A_low,
			const Univcoord_T *B, const Univcoord_T *endB,
			int deltaA, int diagtermA, int diagtermB,
			const Univcoord_T below_slop, const Univcoord_T above_slop) {

  Univcoord_T A_value;

  debug(printf("Entered scalar_one_set with %lu entries\n",endB - B));

  debug(printf("Advancing B from %u up to %u\n",*B,GETPOS(*A_high,*A_low) + deltaA - below_slop));
  while (B < endB && (*B) + below_slop < GETPOS(*A_high,*A_low) + deltaA) {
    B++;
  }

  while (B < endB && (*B) <= (A_value = GETPOS(*A_high,*A_low)) + deltaA + above_slop) {
    debug(printf("(2) Generating output %u %u\n",A_value + diagtermA,(*B) + diagtermB));
    if (A_value + deltaA == (*B)) {
      /* Exclude exact */
    } else if (1 || A_value + diagtermA >= (*B) + diagtermB) { /* higher */
      if (out == diagonals || A_value + diagtermA > out[-1]) {
	*out++ = A_value + diagtermA;
      }
    }

    B++;
  }

  return out;
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

#define NELTS 8
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
      debug(printf("%d %d %d\n",idx,match_start[idx],match_n[idx]));
    }
  }
      
  return;
}


#ifdef HAVE_SSE4_2
static int
binary_search_large (int lowi, int highi,
		     const unsigned char *positions_high, const UINT4 *positions_low, Univcoord_T goal) {
  int middlei;

  debug20(printf("entered binary search with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi < highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,GETPOS_I(positions_high,positions_low,lowi),
		   middlei,GETPOS_I(positions_high,positions_low,middlei),
		   highi,GETPOS_I(positions_high,positions_low,highi),goal));
    if (goal < GETPOS_I(positions_high,positions_low,middlei)) {
      highi = middlei;
    } else if (goal > GETPOS_I(positions_high,positions_low,middlei)) {
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


#ifdef HAVE_SSE4_2
static int
binary_search_univcoord (int lowi, int highi, const Univcoord_T *positions, Univcoord_T goal) {
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


#ifdef HAVE_SSE4_2
/* rare_first_p is false */
static Univcoord_T *
v1_anchorA (Univcoord_T *out, Univcoord_T *diagonals,
	    const Univcoord_T *rare, int lenRare,
	    const unsigned char *freq_high, const UINT4 *freq_low, int lenFreq, int delta_rare,
	    int diagterm_rare, int diagterm_freq,
	    const Univcoord_T below_slop, const Univcoord_T above_slop) {

  const Univcoord_T *initRare = rare;
  const UINT4 *initFreq = freq_low;
  const Univcoord_T *endRare = &(rare[lenRare]);
  const UINT4 *endFreq = &(freq_low[lenFreq]);

  const int kRareSpace = 0;	/* 0. type was uint64_t */

  Univcoord_T valRare, outFreq;
  Univcoord_T outRare;
  Univcoord_T goal;
  int npositions_rare, npositions_freq, j;

  int i;
  Univcoord_T below_slop_plus_1, above_slop_plus_1; /* Because only cmpgt is available, not cmpge */


  below_slop_plus_1 = below_slop + 1; 
  above_slop_plus_1 = above_slop + 1;

  /* 8 comparisons */
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7. type was uint64_t */
  const int kProbe = (0 + 1) * 2 * 4; /* 8. type was uint64_t */

#ifdef HAVE_AVX512
  __m512i Rare_lowbound, Rare_highbound;
  __m128i H;
  __m512i F, F_save;
  __m512i _epi64_offset;
  __mmask16 idx, M_above, M_below;
#elif defined(HAVE_AVX2)
  __m256i Rare_lowbound, Rare_highbound;
  __m128i H;
  __m256i F0, F1, F0_save, F1_save;
  __m256i _epi64_offset;
  __m256i M, M_above, M_below;
  int idx;
#else
  __m128i Rare_lowbound, Rare_highbound;
  __m128i H;
  __m128i F0, F1, F2, F3, F0_save, F1_save, F2_save, F3_save;
  __m128i _epi64_offset;
  __m128i M, M_above, M_below;
  int idx;
#endif

  const Univcoord_T *stopRare = endRare - kRareSpace;
  const UINT4 *stopFreq = endFreq - kFreqSpace;

#ifdef HAVE_AVX512
  _epi64_offset = _mm512_set1_epi64(EPI64_MAX);
#elif defined(HAVE_AVX2)
  _epi64_offset = _mm256_set1_epi64x(EPI64_MAX);
#else
  _epi64_offset = _mm_set1_epi64((__m64) EPI64_MAX);
#endif


  assert(/*lenRare*/(endRare - initRare) <= /*lenFreq*/(endFreq - initFreq));

  /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and in main loop */
  /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */

  /* Avoid issues where valRare - below_slop_plus_1 < 0 */
  while (COMPILER_RARELY(rare < stopRare &&
			 (valRare = (*rare) + delta_rare) <= below_slop)) {
    out = scalar_one_set_anchorA(out, diagonals, rare, freq_high, freq_low, endFreq,
				 delta_rare, diagterm_rare, diagterm_freq,
				 below_slop, above_slop);
    rare++;
    /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and in main loop */
    /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */
  }

  if (COMPILER_LIKELY(rare < stopRare)) {
    valRare = (*rare) + delta_rare;
    outRare = (*rare) + diagterm_rare;

    /* Advance freq */
    npositions_rare = endRare - rare;
    npositions_freq = endFreq - freq_low;

    if (50 * npositions_rare > npositions_freq) {
      /* Use linear */
      while (freq_low < stopFreq && /*maxFreq*/GETPOS_I(freq_high,freq_low,kProbe - 1) + below_slop < valRare) {
	debug(printf("Advancing freq by 8 because %u + slop < %u\n",GETPOS_I(freq_high,freq_low,kProbe - 1),valRare));
	freq_high += kProbe;
	freq_low += kProbe;
      }

    } else {
      /* Use galloping search */
      goal = valRare - below_slop;
      j = 1;
      while (j < npositions_freq && GETPOS_I(freq_high,freq_low,j) < goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions_freq) {
	j = binary_search_large(j >> 1,npositions_freq,freq_high,freq_low,goal);
      } else {
	j = binary_search_large(j >> 1,j,freq_high,freq_low,goal);
      }
      freq_high += j;
      freq_low += j;
    }
  }

  if (freq_low < stopFreq) {
#ifdef HAVE_AVX512
    Rare_lowbound = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_highbound = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);

    H = _mm_lddqu_si128((const __m128i *)(freq_high));
    F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			 _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
    F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
    Rare_lowbound = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_highbound = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);

    H = _mm_lddqu_si128((const __m128i *)(freq_high));
    F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			  _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
    F0 = _mm256_sub_epi64(F0, _epi64_offset);

    H = _mm_bslli_si128(H,4);
    F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			  _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
    F1 = _mm256_sub_epi64(F1, _epi64_offset);
#else
    Rare_lowbound = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
    Rare_highbound = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));

    H = _mm_lddqu_si128((const __m128i *)(freq_high));
    F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		       _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
    F0 = _mm_sub_epi64(F0, _epi64_offset);

    H = _mm_bslli_si128(H,2);
    F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		       _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 2))));
    F1 = _mm_sub_epi64(F1, _epi64_offset);

    H = _mm_bslli_si128(H,2);
    F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		       _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
    F2 = _mm_sub_epi64(F2, _epi64_offset);

    H = _mm_bslli_si128(H,2);
    F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
		       _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 6))));
    F3 = _mm_sub_epi64(F3, _epi64_offset);
#endif
  }

  while (rare < stopRare && freq_low < stopFreq) {
    valRare = (*rare) + delta_rare;
    outRare = (*rare) + diagterm_rare;

    debug(printf("rare is at item %lu with value %u\n",rare - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",
		 freq_low - initFreq,GETPOS_I(freq_high,freq_low,0),GETPOS_I(freq_high,freq_low,7)));

    const unsigned char *freq_high_save = freq_high;
    const UINT4 *freq_low_save = freq_low;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi64_mask(F, Rare_lowbound); /* 8 comparisons */
    M_below = _mm512_cmplt_epi64_mask(F, Rare_highbound);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F0_save = F0;
    F1_save = F1;

    M_above = _mm256_cmpgt_epi64(F1, Rare_lowbound); /* 4 comparisons */
    M_below = _mm256_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);
    idx = _mm256_movemask_pd((__m256d) M);

    M_above = _mm256_cmpgt_epi64(F0, Rare_lowbound); /* 4 comparisons */
    M_below = _mm256_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);
    idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);
#else
    F0_save = F0;
    F1_save = F1;
    F2_save = F2;
    F3_save = F3;

    M_above = _mm_cmpgt_epi64(F3, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F3); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F2, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F2); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F1, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    
    M_above = _mm_cmpgt_epi64(F0, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

    debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

    const unsigned char *freq_high_ptr = freq_high + match_start[idx];
    const UINT4 *freq_low_ptr = freq_low + match_start[idx];

    for (i = 0; i < match_n[idx]; i++) {
      outFreq = GETPOS(*freq_high_ptr++,*freq_low_ptr++) + diagterm_freq;
      debug(printf("(3) Generating output freq %u\n",outFreq));
#if 0
      if (rare_first_p == true) {
	if (out == diagonals || outRare > out[-1]) {
	  *out++ = outRare;
	}
      } else {
	if (out == diagonals || outFreq > out[-1]) {
	  *out++ = outFreq;
	}
      }
#else
      if (outRare == outFreq) {
	/* Exclude exact */
      } else if (out == diagonals || outFreq > out[-1]) {
	*out++ = outFreq;
      }
#endif
    }

    /* Subsequent iterations on freq */
    while (freq_low + kProbe < stopFreq && /*minFreq*/GETPOS_I(freq_high,freq_low,kProbe) <= valRare + above_slop) {
      freq_high += kProbe;
      freq_low += kProbe;
      
#ifdef HAVE_AVX512
      H = _mm_lddqu_si128((const __m128i *)(freq_high));
      F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			   _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
      F = _mm512_sub_epi64(F, _epi64_offset);
    
      M_above = _mm512_cmpgt_epi64_mask(F, Rare_lowbound); /* 8 comparisons */
      M_below = _mm512_cmplt_epi64_mask(F, Rare_highbound);
      idx = M_above & M_below;
#elif defined(HAVE_AVX2)
      H = _mm_lddqu_si128((const __m128i *)(freq_high));
      F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			    _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
      F0 = _mm256_sub_epi64(F0, _epi64_offset);

      H = _mm_bslli_si128(H,4);
      F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			    _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
      F1 = _mm256_sub_epi64(F1, _epi64_offset);

      M_above = _mm256_cmpgt_epi64(F1, Rare_lowbound); /* 4 comparisons */
      M_below = _mm256_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);
      idx = _mm256_movemask_pd((__m256d) M);

      M_above = _mm256_cmpgt_epi64(F0, Rare_lowbound); /* 4 comparisons */
      M_below = _mm256_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);
#else
      H = _mm_lddqu_si128((const __m128i *)(freq_high));
      F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			 _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
      F0 = _mm_sub_epi64(F0, _epi64_offset);

      H = _mm_bslli_si128(H,2);
      F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			 _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 2))));
      F1 = _mm_sub_epi64(F1, _epi64_offset);

      H = _mm_bslli_si128(H,2);
      F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			 _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
      F2 = _mm_sub_epi64(F2, _epi64_offset);

      H = _mm_bslli_si128(H,2);
      F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			 _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 6))));
      F3 = _mm_sub_epi64(F3, _epi64_offset);

      M_above = _mm_cmpgt_epi64(F3, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F3);
      M = _mm_and_si128(M_above,M_below);
      idx = _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F2, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F2);
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F1, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F1);
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
    
      M_above = _mm_cmpgt_epi64(F0, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F0);
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

      debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

      const unsigned char *freq_high_ptr = freq_high + match_start[idx];
      const UINT4 *freq_low_ptr = freq_low + match_start[idx];

      for (i = 0; i < match_n[idx]; i++) {
	outFreq = GETPOS(*freq_high_ptr++,*freq_low_ptr++) + diagterm_freq;
	debug(printf("(4) Generating output freq %u\n",outFreq));
#if 0
	if (rare_first_p == true) {
	  if (out == diagonals || outRare > out[-1]) {
	    *out++ = outRare;
	  }
	} else {
	  if (out == diagonals || outFreq > out[-1]) {
	    *out++ = outFreq;
	  }
	}
#else
	if (outRare == outFreq) {
	  /* Exclude exact */
	} else if (out == diagonals || outFreq > out[-1]) {
	  *out++ = outFreq;
	}
#endif
      }
    }
    /* End of iterations on freq */

    if (freq_low + kProbe >= stopFreq) {
      freq_high += kProbe;
      freq_low += kProbe;
      out = scalar_one_set_anchorA(out, diagonals, rare, freq_high, freq_low, endFreq,
				   delta_rare, diagterm_rare, diagterm_freq,
				   below_slop, above_slop);
    }

    if (++rare < stopRare) {
      valRare = (*rare) + delta_rare; /* For comparison */
      /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */

#ifdef HAVE_AVX512
      Rare_lowbound = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_highbound = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_lowbound = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_highbound = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);
      F0 = F0_save;
      F1 = F1_save;
#else
      Rare_lowbound = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
      Rare_highbound = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));
      F0 = F0_save;
      F1 = F1_save;
      F2 = F2_save;
      F3 = F3_save;
#endif

      /* Advance freq */
      freq_high = freq_high_save;
      freq_low = freq_low_save;

      npositions_rare = endRare - rare;
      npositions_freq = endFreq - freq_low;

      if (50 * npositions_rare > npositions_freq) {
	/* Use linear */
	while (freq_low < stopFreq && /*maxFreq*/GETPOS_I(freq_high,freq_low,kProbe - 1) + below_slop < valRare) {
	  debug(printf("Advancing freq by %d because %u + below_slop < %u\n",
		       kProbe,GETPOS_I(freq_high,freq_low,kProbe - 1),valRare));
	  freq_high += kProbe;
	  freq_low += kProbe;
	}

      } else {
	/* Use galloping search */
	goal = valRare - below_slop;
	j = 1;
	while (j < npositions_freq && GETPOS_I(freq_high,freq_low,j) < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions_freq) {
	  j = binary_search_large(j >> 1,npositions_freq,freq_high,freq_low,goal);
	} else {
	  j = binary_search_large(j >> 1,j,freq_high,freq_low,goal);
	}
	freq_high += j;
	freq_low += j;
      }

      if (freq_low < stopFreq && freq_low != freq_low_save) {
#ifdef HAVE_AVX512
	H = _mm_lddqu_si128((const __m128i *)(freq_high));
	F = _mm512_add_epi64(_mm512_slli_epi64(_mm512_cvtepi8_epi64(H),32),
			     _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i *)(freq_low))));
	F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
	H = _mm_lddqu_si128((const __m128i *)(freq_high));
	F0 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
	F0 = _mm256_sub_epi64(F0, _epi64_offset);

	H = _mm_bslli_si128(H,4);
	F1 = _mm256_add_epi64(_mm256_slli_epi64(_mm256_cvtepi8_epi64(H),32),
			      _mm256_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
	F1 = _mm256_sub_epi64(F1, _epi64_offset);
#else
	H = _mm_lddqu_si128((const __m128i *)(freq_high));
	F0 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low))));
	F0 = _mm_sub_epi64(F0, _epi64_offset);

	H = _mm_bslli_si128(H,2);
	F1 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 2))));
	F1 = _mm_sub_epi64(F1, _epi64_offset);

	H = _mm_bslli_si128(H,2);
	F2 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 4))));
	F2 = _mm_sub_epi64(F2, _epi64_offset);

	H = _mm_bslli_si128(H,2);
	F3 = _mm_add_epi64(_mm_slli_epi64(_mm_cvtepi8_epi64(H),32),
			   _mm_cvtepi32_epi64(_mm_lddqu_si128((const __m128i *)(freq_low + 6))));
	F3 = _mm_sub_epi64(F3, _epi64_offset);
#endif
      }

      debug(printf("Advancing rare to %lu to compare with stop %lu\n",
		   rare - initRare,stopRare - initRare));
      debug(printf("Restoring freq to %lu\n",freq - initFreq));
    }
  }

  return scalar_anchorA(out, diagonals, rare, endRare, freq_high, freq_low, endFreq,
			delta_rare, diagterm_rare, diagterm_freq,
			below_slop, above_slop);
}


/* rare_first_p is true */
static Univcoord_T *
v1_anchorB (Univcoord_T *out, Univcoord_T *diagonals,
	    const unsigned char *rare_high, const UINT4 *rare_low, int lenRare,
	    const Univcoord_T *freq, int lenFreq, int delta_rare,
	    int diagterm_rare, int diagterm_freq,
	    const Univcoord_T below_slop, const Univcoord_T above_slop) {

  const UINT4 *initRare = rare_low;
  const Univcoord_T *initFreq = freq;
  const UINT4 *endRare = &(rare_low[lenRare]);
  const Univcoord_T *endFreq = &(freq[lenFreq]);

  const int kRareSpace = 0;	/* 0. type was uint64_t */

  Univcoord_T valRare, outRare;
  Univcoord_T outFreq;
  Univcoord_T goal;
  int npositions_rare, npositions_freq, j;

  int i;
  Univcoord_T below_slop_plus_1, above_slop_plus_1; /* Because only cmpgt is available, not cmpge */


  below_slop_plus_1 = below_slop + 1; 
  above_slop_plus_1 = above_slop + 1;

  /* 8 comparisons */
  const int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7. type was uint64_t */
  const int kProbe = (0 + 1) * 2 * 4; /* 8. type was uint64_t */

#ifdef HAVE_AVX512
  __m512i Rare_lowbound, Rare_highbound;
  __m512i F, F_save;
  __m512i _epi64_offset;
  __mmask16 idx, M_above, M_below;
#elif defined(HAVE_AVX2)
  __m256i Rare_lowbound, Rare_highbound;
  __m256i F0, F1, F0_save, F1_save;
  __m256i _epi64_offset;
  __m256i M, M_above, M_below;
  unsigned int idx;
#else
  __m128i Rare_lowbound, Rare_highbound;
  __m128i F0, F1, F2, F3, F0_save, F1_save, F2_save, F3_save;
  __m128i _epi64_offset;
  __m128i M, M_above, M_below;
  unsigned int idx;
#endif

  const UINT4 *stopRare = endRare - kRareSpace;
  const Univcoord_T *stopFreq = endFreq - kFreqSpace;

#ifdef HAVE_AVX512
  _epi64_offset = _mm512_set1_epi64(EPI64_MAX);
#elif defined(HAVE_AVX2)
  _epi64_offset = _mm256_set1_epi64x(EPI64_MAX);
#else
  _epi64_offset = _mm_set1_epi64((__m64) EPI64_MAX);
#endif


  assert(/*lenRare*/(endRare - initRare) <= /*lenFreq*/(endFreq - initFreq));

  /* valRare = GETPOS(*rare_high,*rare_low) + delta_rare; -- For comparison, now computed inside condition and in main loop */
  /* outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- For output, not needed here */

  /* Avoid issues where valRare - below_slop_plus_1 < 0 */
  while (COMPILER_RARELY(rare_low < stopRare &&
			 (valRare = GETPOS(*rare_high,*rare_low)) <= below_slop)) {
    out = scalar_one_set_anchorB(out, diagonals, rare_high, rare_low, freq, endFreq,
				 delta_rare, diagterm_rare, diagterm_freq,
				 below_slop, above_slop);
    rare_high++;
    rare_low++;
    /* valRare = GETPOS(*rare_high,*rare_low) + delta_rare; -- For comparison, now computed inside condition and in main loop */
    /* outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- For output, not needed here */
  }

  if (COMPILER_LIKELY(rare_low < stopRare)) {
    valRare = GETPOS(*rare_high,*rare_low) + delta_rare;
    outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;

    /* Advance freq */
    npositions_rare = endRare - rare_low;
    npositions_freq = endFreq - freq;

    if (50 * npositions_rare > npositions_freq) {
      /* Use linear */
      while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	debug(printf("Advancing freq by 8 because %u + below_slop < %u\n",freq[kProbe - 1],valRare));
	freq += kProbe/*8*/;
      }

    } else {
      /* Use galloping search */
      goal = valRare - below_slop;
      j = 1;
      while (j < npositions_freq && freq[j] < goal) {
	j <<= 1;		/* gallop by 2 */
      }
      if (j >= npositions_freq) {
	j = binary_search_univcoord(j >> 1,npositions_freq,freq,goal);
      } else {
	j = binary_search_univcoord(j >> 1,j,freq,goal);
      }
      freq += j;
    }
  }

  if (freq < stopFreq) {
#ifdef HAVE_AVX512
    Rare_lowbound = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_highbound = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);

    F = _mm512_loadu_si512((const __m512i *)(freq));
    F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
    Rare_lowbound = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
    Rare_highbound = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);

    F0 = _mm256_loadu_si256((const __m256i *)(freq));
    F0 = _mm256_sub_epi64(F0, _epi64_offset);
    F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
    F1 = _mm256_sub_epi64(F1, _epi64_offset);
#else
    Rare_lowbound = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
    Rare_highbound = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));

    F0 = _mm_lddqu_si128((const __m128i *)(freq));
    F0 = _mm_sub_epi64(F0, _epi64_offset);
    F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
    F1 = _mm_sub_epi64(F1, _epi64_offset);
    F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
    F2 = _mm_sub_epi64(F2, _epi64_offset);
    F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
    F3 = _mm_sub_epi64(F3, _epi64_offset);
#endif
  }

  while (rare_low < stopRare && freq < stopFreq) {
    valRare = GETPOS(*rare_high,*rare_low) + delta_rare;
    outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare;

    debug(printf("rare is at item %lu with value %u\n",rare_low - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",
		 freq - initFreq,freq[0],freq[kProbe - 1]));

    const Univcoord_T *freq_save = freq;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi64_mask(F, Rare_lowbound); /* 8 comparisons */
    M_below = _mm512_cmplt_epi64_mask(F, Rare_highbound);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F0_save = F0;
    F1_save = F1;

    M_above = _mm256_cmpgt_epi64(F1, Rare_lowbound); /* 4 comparisons */
    M_below = _mm256_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);
    idx = _mm256_movemask_pd((__m256d) M);

    M_above = _mm256_cmpgt_epi64(F0, Rare_lowbound); /* 4 comparisons */
    M_below = _mm256_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);
    idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);
#else
    F0_save = F0;
    F1_save = F1;
    F2_save = F2;
    F3_save = F3;

    M_above = _mm_cmpgt_epi64(F3, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F3); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F2, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F2); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F1, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

    M_above = _mm_cmpgt_epi64(F0, Rare_lowbound); /* 2 comparisons */
    M_below = _mm_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

    debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

    const Univcoord_T *freq_ptr = freq + match_start[idx];

    for (i = 0; i < match_n[idx]; i++) {
      outFreq = *(freq_ptr++) + diagterm_freq;
      debug(printf("(3) Generating output rare %u\n",outRare));
#if 0
      if (rare_first_p == true) {
	if (out == diagonals || outRare > out[-1]) {
	  *out++ = outRare;
	}
      } else {
	if (out == diagonals || outFreq > out[-1]) {
	  *out++ = outFreq;
	}
      }
#else
      if (outRare == outFreq) {
	/* Exclude exact */
      } else if (out == diagonals || outRare > out[-1]) {
	*out++ = outRare;
      }
#endif
    }

    /* Subsequent iterations on freq */
    while (freq + kProbe < stopFreq && /*minFreq*/freq[kProbe] <= valRare + above_slop) {
      freq += kProbe;
      
#ifdef HAVE_AVX512
      F = _mm512_loadu_si512((const __m512i *)(freq));
      F = _mm512_sub_epi64(F, _epi64_offset);
    
      M_above = _mm512_cmpgt_epi64_mask(F, Rare_lowbound);
      M_below = _mm512_cmplt_epi64_mask(F, Rare_highbound);
      idx = M_above & M_below;
#elif defined(HAVE_AVX2)
      F0 = _mm256_loadu_si256((const __m256i *)(freq));
      F0 = _mm256_sub_epi64(F0, _epi64_offset);
      F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
      F1 = _mm256_sub_epi64(F1, _epi64_offset);

      M_above = _mm256_cmpgt_epi64(F1, Rare_lowbound); /* 4 comparisons */
      M_below = _mm256_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);
      idx = _mm256_movemask_pd((__m256d) M);

      M_above = _mm256_cmpgt_epi64(F0, Rare_lowbound); /* 4 comparisons */
      M_below = _mm256_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);
      idx = (idx << 4) + _mm256_movemask_pd((__m256d) M);
#else
      F0 = _mm_lddqu_si128((const __m128i *)(freq));
      F0 = _mm_sub_epi64(F0, _epi64_offset);
      F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
      F1 = _mm_sub_epi64(F1, _epi64_offset);
      F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
      F2 = _mm_sub_epi64(F2, _epi64_offset);
      F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
      F3 = _mm_sub_epi64(F3, _epi64_offset);
      
      M_above = _mm_cmpgt_epi64(F3, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F3); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F2, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F2); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F1, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F1); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);

      M_above = _mm_cmpgt_epi64(F0, Rare_lowbound); /* 2 comparisons */
      M_below = _mm_cmpgt_epi64(Rare_highbound, F0); /* cmplt not available */
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 2) + _mm_movemask_pd((__m128d) M);
#endif

      debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

      const Univcoord_T *freq_ptr = freq + match_start[idx];

      for (i = 0; i < match_n[idx]; i++) {
	outFreq = *(freq_ptr++) + diagterm_freq;
	debug(printf("(4) Generating output rare %u\n",outRare));
#if 0
	if (rare_first_p == true) {
	  if (out == diagonals || outRare > out[-1]) {
	    *out++ = outRare;
	  }
	} else {
	  if (out == diagonals || outFreq > out[-1]) {
	    *out++ = outFreq;
	  }
	}
#else
	if (outRare == outFreq) {
	  /* Exclude exact */
	} else if (out == diagonals || outRare > out[-1]) {
	  *out++ = outRare;
	}
#endif
      }
    }
    /* End of iterations on freq */

    if (freq + kProbe >= stopFreq) {
      freq += kProbe;
      out = scalar_one_set_anchorB(out, diagonals, rare_high, rare_low, freq, endFreq,
				   delta_rare, diagterm_rare, diagterm_freq,
				   below_slop, above_slop);
    }

    ++rare_high;
    if (++rare_low < stopRare) {
      valRare = GETPOS(*rare_high,*rare_low) + delta_rare; /* For comparison */
      /* outRare = GETPOS(*rare_high,*rare_low) + diagterm_rare; -- For output, not needed here */

#ifdef HAVE_AVX512
      Rare_lowbound = _mm512_set1_epi64(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_highbound = _mm512_set1_epi64(valRare + above_slop_plus_1 - EPI64_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_lowbound = _mm256_set1_epi64x(valRare - below_slop_plus_1 - EPI64_MAX);
      Rare_highbound = _mm256_set1_epi64x(valRare + above_slop_plus_1 - EPI64_MAX);
      F0 = F0_save;
      F1 = F1_save;
#else
      Rare_lowbound = _mm_set1_epi64((__m64) (valRare - below_slop_plus_1 - EPI64_MAX));
      Rare_highbound = _mm_set1_epi64((__m64) (valRare + above_slop_plus_1 - EPI64_MAX));
      F0 = F0_save;
      F1 = F1_save;
      F2 = F2_save;
      F3 = F3_save;
#endif

      /* Advance freq */
      freq = freq_save;

      npositions_rare = endRare - rare_low;
      npositions_freq = endFreq - freq;

      if (50 * npositions_rare > npositions_freq) {
	/* Use linear */
	while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	  debug(printf("Advancing freq by %d because %u + slop < %u\n",
		       kProbe,freq[kProbe - 1],valRare));
	  freq += kProbe;
	}

      } else {
	/* Use galloping search */
	goal = valRare - below_slop;
	j = 1;
	while (j < npositions_freq && freq[j] < goal) {
	  j <<= 1;		/* gallop by 2 */
	}
	if (j >= npositions_freq) {
	  j = binary_search_univcoord(j >> 1,npositions_freq,freq,goal);
	} else {
	  j = binary_search_univcoord(j >> 1,j,freq,goal);
	}
	freq += j;
      }

      if (freq < stopFreq && freq != freq_save) {
#ifdef HAVE_AVX512
	F = _mm512_loadu_si512((const __m512i *)(freq));
	F = _mm512_sub_epi64(F, _epi64_offset);
#elif defined(HAVE_AVX2)
	F0 = _mm256_loadu_si256((const __m256i *)(freq));
	F0 = _mm256_sub_epi64(F0, _epi64_offset);
	F1 = _mm256_loadu_si256((const __m256i *)(freq + 4));
	F1 = _mm256_sub_epi64(F1, _epi64_offset);
#else
	F0 = _mm_lddqu_si128((const __m128i *)(freq));
	F0 = _mm_sub_epi64(F0, _epi64_offset);
	F1 = _mm_lddqu_si128((const __m128i *)(freq + 2));
	F1 = _mm_sub_epi64(F1, _epi64_offset);
	F2 = _mm_lddqu_si128((const __m128i *)(freq + 4));
	F2 = _mm_sub_epi64(F2, _epi64_offset);
	F3 = _mm_lddqu_si128((const __m128i *)(freq + 6));
	F3 = _mm_sub_epi64(F3, _epi64_offset);
#endif
      }

      debug(printf("Advancing rare to %lu to compare with stop %lu\n",
		   rare_low - initRare,stopRare - initRare));
      debug(printf("Restoring freq to %lu\n",freq - initFreq));
    }
  }

  return scalar_anchorB(out, diagonals, rare_high, rare_low, endRare, freq, endFreq,
			delta_rare, diagterm_rare, diagterm_freq,
			below_slop, above_slop);
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
Intersect_higher (Univcoord_T *diagonals,
		  const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		  const Univcoord_T *set2, const int length2,
		  const Univcoord_T slop, const Univcoord_T insertion_slop) {
  Univcoord_T *out = diagonals;
#ifdef DEBUG15
  Univcoord_T *nonsimd_diagonals;
  int nonsimd_ndiagonals, k;
#endif

  debug(printf("Entered Intersect_higher with %d and %d items\n",length1,length2));

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
    debug(printf("Exiting Intersect_higher with 0 items\n"));
    return 0;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      out = v1_anchorB(out, diagonals, set1_high, set1_low, length1, set2, length2,
		       /*delta_rare*/+diagterm1,
		       /*diagterm_rare*/diagterm1, /*diagterm_freq*/0,
		       /*below_slop*/slop, /*above_slop*/insertion_slop);
    } else {
      out = v1_anchorA(out, diagonals, set2, length2, set1_high, set1_low, length1,
		       /*delta_rare*/-diagterm1,
		       /*diagterm_rare*/0, /*diagterm_freq*/diagterm1,
		       /*below_slop*/insertion_slop, /*above_slop*/slop);
    }

    debug(printf("Exiting Intersect_approx_higher with %d items\n",(out - diagonals)));
#ifdef DEBUG
    for (i = 0; i < (out - diagonals); i++) {
      printf("%u\n",diagonals[i]);
    }
#endif


#ifdef DEBUG15
    if (length1 <= length2) {
      nonsimd_diagonals = MALLOC(length2*sizeof(unsigned int));
    } else {
      nonsimd_diagonals = MALLOC(length1*sizeof(unsigned int));
    }

    nonsimd_ndiagonals = Intersect_approx_higher_old(nonsimd_diagonals,
						     set1_high,set1_low,length1,diagterm1,
						     set2,length2,slop);
    if ((out - diagonals) != nonsimd_ndiagonals) {
      printf("Number of diagonals %d != %d\n",(out - diagonals),nonsimd_ndiagonals);

      for (k = 0; k < (out - diagonals); k++) {
	printf("%u\n",diagonals[k]);
      }
      printf("\n");
      for (k = 0; k < nonsimd_ndiagonals; k++) {
	printf("%u\n",nonsimd_diagonals[k]);
      }
      printf("\n");

      abort();

    } else {
#if 0
      for (k = 0; k < (out - diagonals); k++) {
	printf("%u\n",diagonals[k]);
      }
      printf("\n");
      for (k = 0; k < nonsimd_ndiagpairs; k++) {
	printf("%u\n",nonsimd_diagonals[k]);
      }
      printf("\n");
#endif

      for (k = 0; k < (out - diagonals); k++) {
	if (diagonals[k] != nonsimd_diagonals[k]) {
	  printf("%u != %u\n",diagonals[k],nonsimd_diagonals[k]);
	  abort();
	}
      }
    }

    FREE(nonsimd_diagonals);
#endif

    return (out - diagonals);
  }
}

#else

int
Intersect_higher (Univcoord_T *diagonals,
		  const unsigned char *set1_high, const UINT4 *set1_low, const int length1, int diagterm1,
		  const Univcoord_T *set2, const int length2,
		  const Univcoord_T slop, const Univcoord_T insertion_slop) {
  Univcoord_T *out = diagonals;

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else {
    if (length1 <= length2) {
      /* A_first_p is true */
      out = scalar_anchorB(out, diagonals, set1_high, set1_low, /*endA*/&(set1_low[length1]),
			   set2, /*endB*/&(set2[length2]),
			   /*deltaA*/+diagterm1,
			   /*diagtermA*/diagterm1, /*diagtermB*/0,
			   /*below_slop*/slop, /*above_slop*/insertion_slop);
    } else {
      /* A_first_p is false */
      out = scalar_anchorA(out, diagonals, set2, /*endA*/&(set2[length2]),
			   set1_high, set1_low, /*endB*/&(set1_low[length1]),
			   /*deltaA*/-diagterm1,
			   /*diagtermA*/0, /*diagtermB*/diagterm1,
			   /*below_slop*/insertion_slop, /*above_slop*/slop);
    }

    return (out - diagonals);
  }
}
#endif


#ifdef DEBUG15
/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
int
Intersect_approx_higher_old (Univcoord_T *diagonals,
			     unsigned char *positions1_high, UINT4 *positions1_low,
			     int npositions1, int diagterm1,
			     Univcoord_T *positions0, int npositions0,
			     Chrpos_T maxdistance) {
  int ndiagonals = 0;
  Univcoord_T local_goal, last_positions0;
  Univcoord_T diagterm, delta;
  Univcoord_T last_diagonal, this_diagonal;
  unsigned char *ptr1_high, *start1;
  UINT4 *ptr1_low;
  int j;
#ifdef DEBUG9
  Univcoord_T *ptr0;
#endif


  diagterm = (Univcoord_T) diagterm1;	/* local_goal based on larger list */
  delta = (Univcoord_T) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1_high;
  /* end1 = &(positions1_high[npositions1]); */

  debug(printf("Intersect_approx_higher with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG9
  for (ptr0 = positions0;  ptr0 < &(positions0[npositions0]); ptr0++) {
    printf("%llu\n",*ptr0);
  }
  printf("\n");

  for (ptr1_high = positions1_high, ptr1_low = positions1_low;
       ptr1_high < &(positions1_high[npositions1]); ptr1_high++, ptr1_low++) {
    printf("%llu + %d\n",GETPOS(*ptr1_high,*ptr1_low),diagterm1);
  }
  printf("\n");
#endif


#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && GETPOS(*positions1_high,*positions1_low) < (Univcoord_T) -diagterm1) {
      ++positions1_high;
      ++positions1_low;
      --npositions1;
    }
  }
#endif

  last_positions0 = (Univcoord_T) -1;
  while (npositions0 > 0) {
    if (*positions0 == last_positions0) {
      /* Skip duplicate in positions0 */
      /* last_positions0 = *positions0 */
      ++positions0;
      --npositions0;
      
    } else {
      local_goal = (*positions0) + delta;
      debug9(printf("intersection list 0: %d:%llu => local_goal %llu.  Searching from npositions1 %d\n",
		   npositions0,*positions0,local_goal,npositions1));
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
#ifdef DEBUG9
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,GETPOS(*positions1_high,*positions1_low));
      }
#endif

      if (npositions1 <= 0) {
	/* Check backwards only */
	debug9(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1_high = &(positions1_high[-1]);
	ptr1_low = &(positions1_low[-1]);
	if (ptr1_high >= start1) {
	  debug9(printf(" prev %d:%llu?",npositions1-1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	    debug9(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug9(printf("\n"));
	
	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
	return ndiagonals;
	
      } else if (GETPOS(*positions1_high,*positions1_low) == local_goal) {
	/* Found local goal.  Save and advance */
	debug9(printf("    intersection list 1: %d:%llu  exact\n",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	/* diagonals[ndiagonals++] = local_goal + diagterm; */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1_high;
	++positions1_low;
	--npositions1;
	
      } else {
	debug9(printf("    intersection list 1: %d:%llu  checking for approximate:",
		     npositions1,GETPOS(*positions1_high,*positions1_low)));
	ptr1_high = &(positions1_high[-1]); /* closest position < local_goal */
	ptr1_low = &(positions1_low[-1]); /* closest position < local_goal */
	if (ptr1_high >= start1) {
	  debug9(printf(" prev %d:%llu?",npositions1+1,GETPOS(*ptr1_high,*ptr1_low)));
	  if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	    debug9(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	    } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}

	ptr1_high = &(positions1_high[0]); /* closest position > local_goal */
	ptr1_low = &(positions1_low[0]); /* closest position > local_goal */
	debug9(printf(" at %d:%llu?",npositions1,GETPOS(*ptr1_high,*ptr1_low)));
	if (GETPOS(*ptr1_high,*ptr1_low) /*(omit for higher) + maxdistance*/ >= local_goal && GETPOS(*ptr1_high,*ptr1_low) <= local_goal + maxdistance) {
	  debug9(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = GETPOS(*ptr1_high,*ptr1_low) + diagterm;
	  } else if ((this_diagonal = GETPOS(*ptr1_high,*ptr1_low) + diagterm) != last_diagonal) {
	    last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	  }
	}
	debug9(printf("\n"));
      }

      last_positions0 = *positions0;
      ++positions0;
      --npositions0;
    }
  }
  debug9(printf("\n"));

  return ndiagonals;
}

#endif


void
Intersect_higher_setup () {
  initialize_match_values();
  return;
}  
