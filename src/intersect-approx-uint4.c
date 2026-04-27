static char rcsid[] = "$Id: 57c05afd123398694aab6279db9137b44f4eb8ff $";
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

/* Minimum requirement for SIMD is SSE2 */

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "bool.h"
#include "intersect-approx-uint4.h"

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


#define INCLUDE_EXACT 1

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

/* Intersect_approx_old */
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
static unsigned int *
scalar (unsigned int *diagonals, unsigned int **init_diagonals, unsigned int **end_diagonals,
	const unsigned int *A, const unsigned int *endA,
	const unsigned int *B, const unsigned int *endB, int deltaA,
	int diagtermA, int diagtermB,
	const unsigned int below_slop, const unsigned int above_slop,
	bool A_first_p) {

  const unsigned int *B_save;
  int allocation, ndiagonals;
  unsigned int *more_diagonals;

  debug(printf("Entered scalar with %lu and %lu entries\n",endA - A,endB - B));

  while (A < endA) {
    /* Advance B until (A - below_slop) */
    while (B < endB && (*B) + below_slop < (*A) + deltaA) {
      B++;
    }

    B_save = B;
    debug(printf("Advancing B from %u up to %u\n",*B,(*A) + deltaA + above_slop));
    while (B < endB && (*B) <= (*A) + deltaA + above_slop) {
      debug(printf("(1) Generating output %u %u\n",(*A) + diagtermA,(*B) + diagtermB));
      if (diagonals >= (*end_diagonals) /* same as diagonals + 2 > (*end_diagonals) */) {
	allocation = (*end_diagonals) - (*init_diagonals);
	more_diagonals = (unsigned int *) MALLOC(2 * allocation * sizeof(unsigned int));

	ndiagonals = diagonals - (*init_diagonals);
	memcpy(more_diagonals,*init_diagonals,ndiagonals * sizeof(unsigned int));

	FREE(*init_diagonals);
	*init_diagonals = more_diagonals;
	*end_diagonals = &(more_diagonals[2 * allocation]);
	diagonals = &(more_diagonals[ndiagonals]);
      }
	
      assert(diagonals + 1 < *end_diagonals);
#ifdef INCLUDE_EXACT
      if (A_first_p == true) {
	*diagonals++ = (*A) + diagtermA;
	*diagonals++ = (*B) + diagtermB;
      } else {
	*diagonals++ = (*B) + diagtermB;
	*diagonals++ = (*A) + diagtermA;
      }
#else
      if ((*A) + deltaA != (*B)) {
	if (A_first_p == true) {
	  *diagonals++ = (*A) + diagtermA;
	  *diagonals++ = (*B) + diagtermB;
	} else {
	  *diagonals++ = (*B) + diagtermB;
	  *diagonals++ = (*A) + diagtermA;
	}
      }
#endif
      assert(diagonals <= (*end_diagonals));

      B++;
    }
    B = B_save;

    A++;
  }

  return diagonals;
}


/* Performs scalar against a single value of A */
static inline unsigned int *
scalar_one_set (unsigned int *diagonals, unsigned int **init_diagonals, unsigned int **end_diagonals,
		const unsigned int *A, const unsigned int *B, const unsigned int *endB,
		int deltaA, int diagtermA, int diagtermB,
		const unsigned int below_slop, const unsigned int above_slop,
		bool A_first_p) {
  int allocation, ndiagonals;
  unsigned int *more_diagonals;

  debug(printf("Entered scalar_one_set with %lu entries\n",endB - B));

  debug(printf("Advancing B from %u up to %u\n",*B,(*A) + deltaA - below_slop));
  while (B < endB && (*B) + below_slop < (*A) + deltaA) {
    B++;
  }

  while (B < endB && (*B) <= (*A) + deltaA + above_slop) {
    debug(printf("(2) Generating output %u %u\n",(*A) + diagtermA,(*B) + diagtermB));

    if (diagonals >= (*end_diagonals) /* same as diagonals + 2 > (*end_diagonals) */) {
      allocation = (*end_diagonals) - (*init_diagonals);
      debug(printf("Reallocating from diagonals %p to end %p (allocation %d)\n",
		   diagonals,*end_diagonals,allocation));
      more_diagonals = (unsigned int *) MALLOC(2 * allocation * sizeof(unsigned int));
      
      ndiagonals = diagonals - (*init_diagonals);
      debug(printf("Copying from diagonals %p to out %p (ndiagonals %d)\n",
		   diagonals,more_diagonals,ndiagonals));
      memcpy(more_diagonals,*init_diagonals,ndiagonals * sizeof(unsigned int));

      FREE(*init_diagonals);
      *init_diagonals = more_diagonals;
      *end_diagonals = &(more_diagonals[2 * allocation]);
      diagonals = &(more_diagonals[ndiagonals]);
    }

    assert(diagonals + 1 < *end_diagonals);
#ifdef INCLUDE_EXACT
    if (A_first_p == true) {
      *diagonals++ = (*A) + diagtermA;
      *diagonals++ = (*B) + diagtermB;
    } else {
      *diagonals++ = (*B) + diagtermB;
      *diagonals++ = (*A) + diagtermA;
    }
#else
    if ((*A) + deltaA != (*B)) {
      if (A_first_p == true) {
	*diagonals++ = (*A) + diagtermA;
	*diagonals++ = (*B) + diagtermB;
      } else {
	*diagonals++ = (*B) + diagtermB;
	*diagonals++ = (*A) + diagtermA;
      }
    }
#endif
    assert(diagonals <= (*end_diagonals));

    B++;
  }

  return diagonals;
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

#ifdef HAVE_AVX512
#define NELTS 16		/* 512 / 32 bits per UINT4 */
#define MAX_IDX 65536		/* 2^16 */
#else
#define NELTS 8			/* 256 / 32 bits per UINT4 */
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
      debug(printf("%d %d %d\n",idx,match_start[idx],match_n[idx]));
      assert(match_n[idx] <= NELTS);
    }
  }
      
  return;
}


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
static unsigned int *
v1 (unsigned int *diagonals, unsigned int **init_diagonals, unsigned int **end_diagonals,
    const unsigned int *rare, int lenRare,
    const unsigned int *freq, int lenFreq, int delta_rare,
    int diagterm_rare, int diagterm_freq,
    const unsigned int below_slop, const unsigned int above_slop, bool rare_first_p) {

  const unsigned int *initRare = rare;
  const unsigned int *initFreq = freq;
  const unsigned int *endRare = &(rare[lenRare]);
  const unsigned int *endFreq = &(freq[lenFreq]);

  const unsigned int kRareSpace = 0;	/* 0. type was uint64_t */

#ifndef INCLUDE_EXACT
  unsigned int valFreq;
#endif
  unsigned int valRare, outRare;
  unsigned int goal;
  int npositions_rare, npositions_freq, j;

  int i;
  const unsigned int below_slop_plus_1 = below_slop + 1; /* Because only cmpgt is available, not cmpge */
  const unsigned int above_slop_plus_1 = above_slop + 1; /* Because only cmpgt is available, not cmpge */
  int allocation, ndiagonals;
  unsigned int *more_diagonals;

#ifdef HAVE_AVX512
  const unsigned int kProbe = (0 + 1) * 2 * 8; /* 16 32-mers in 512 bits */
  const unsigned int kFreqSpace = 2 * 8 * (0 + 1) - 1; /* 15 */

  __m512i Rare_low, Rare_high;
  __m512i F, F_save;
  __m512i _epi32_offset;
  __mmask16 idx, M_above, M_below;
#elif defined(HAVE_AVX2)
  const unsigned int kProbe = (0 + 1) * 2 * 4; /* 8 32-mers in 256 bits */
  const unsigned int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7 */

  __m256i Rare_low, Rare_high;
  __m256i F, F_save;
  __m256i _epi32_offset;
  __m256i M, M_above, M_below;
  int idx;			/* equivalent to __mmask8 */
#else
  const unsigned int kProbe = (0 + 1) * 2 * 4; /* 8 32-mers in 128+128 bits */
  const unsigned int kFreqSpace = 2 * 4 * (0 + 1) - 1; /* 7 */

  __m128i Rare_low, Rare_high;
  __m128i F0, F1, F0_save, F1_save;
  __m128i _epi32_offset;
  __m128i M, M_above, M_below;
  int idx;			/* equivalent to __mmask8 */
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

  /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and inside main loop */
  /* outRare = (*rare) + diagterm_rare; -- For output, now computed inside main loop */

  /* Avoid issues where valRare - below_slop_plus_1 < 0 */
  while (COMPILER_RARELY(rare < stopRare &&
			 (valRare = (*rare) + delta_rare) < below_slop_plus_1)) {
    diagonals = scalar_one_set(diagonals, &(*init_diagonals), &(*end_diagonals),
			       rare, freq, endFreq,
			       delta_rare, diagterm_rare, diagterm_freq,
			       below_slop, above_slop, rare_first_p);
    rare++;
    /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and inside main loop */
    /* outRare = (*rare) + diagterm_rare; -- For output, now computed inside main loop  */
  }

  if (COMPILER_LIKELY(rare < stopRare)) {
    valRare = (*rare) + delta_rare;
    outRare = (*rare) + diagterm_rare; /* For output */

    /* Advance freq */
    npositions_rare = endRare - rare;
    npositions_freq = endFreq - freq;

    if (50 * npositions_rare > npositions_freq) {
      /* Use linear */
      while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	debug(printf("Advancing freq by one vector because %u + below_slop < %u\n",
		     freq[kProbe - 1],valRare));
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
	j = binary_search_uint4(j >> 1,npositions_freq,freq,goal);
      } else {
	j = binary_search_uint4(j >> 1,j,freq,goal);
      }
      freq += j;
    }
  }

  if (freq < stopFreq) {
#ifdef HAVE_AVX512
    Rare_low = _mm512_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_high = _mm512_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F = _mm512_loadu_si512((const __m512i *)(freq));
    F = _mm512_sub_epi32(F, _epi32_offset);
#elif defined(HAVE_AVX2)
    Rare_low = _mm256_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_high = _mm256_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F = _mm256_loadu_si256((const __m256i *)(freq));
    F = _mm256_sub_epi32(F, _epi32_offset);
#else
    Rare_low = _mm_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_high = _mm_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F0 = LOAD_SI128((const __m128i *)(freq));
    F0 = _mm_sub_epi32(F0, _epi32_offset);
    F1 = LOAD_SI128((const __m128i *)(freq + 4));
    F1 = _mm_sub_epi32(F1, _epi32_offset);
#endif
  }

  while (rare < stopRare && freq < stopFreq) {
    debug(printf("rare is at item %lu with value %u\n",rare - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",freq - initFreq,freq[0],freq[kProbe - 1]));

    const unsigned int *freq_save = freq;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi32_mask(F, Rare_low);
    M_below = _mm512_cmplt_epi32_mask(F, Rare_high);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F_save = F;

    M_above = _mm256_cmpgt_epi32(F, Rare_low);
    M_below = _mm256_cmpgt_epi32(Rare_high, F); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);

    idx = _mm256_movemask_ps((__m256) M);
#else
    F0_save = F0;
    F1_save = F1;

    M_above = _mm_cmpgt_epi32(F1, Rare_low);
    M_below = _mm_cmplt_epi32(F1, Rare_high);
    M = _mm_and_si128(M_above,M_below);
    idx = _mm_movemask_ps((__m128) M);

    M_above = _mm_cmpgt_epi32(F0, Rare_low);
    M_below = _mm_cmplt_epi32(F0, Rare_high);
    M = _mm_and_si128(M_above,M_below);
    idx = (idx << 4) + _mm_movemask_ps((__m128) M);
#endif

    debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

    if (diagonals + 2*match_n[idx] > (*end_diagonals)) {
      allocation = (*end_diagonals) - (*init_diagonals);
      debug(printf("Reallocating from diagonals %p to end %p (allocation %d)\n",
		   diagonals,*end_diagonals,allocation));
      more_diagonals = (unsigned int *) MALLOC(2 * allocation * sizeof(unsigned int));
      
      ndiagonals = diagonals - (*init_diagonals);
      debug(printf("Copying from diagonals %p to out %p (ndiagonals %d)\n",
		   diagonals,more_diagonals,ndiagonals));
      memcpy(more_diagonals,*init_diagonals,ndiagonals * sizeof(unsigned int));
      
      FREE(*init_diagonals);
      *init_diagonals = more_diagonals;
      *end_diagonals = &(more_diagonals[2 * allocation]);
      diagonals = &(more_diagonals[ndiagonals]);
      assert(diagonals + 2*match_n[idx] <= (*end_diagonals));
    }

    const unsigned int *freq_ptr = freq + match_start[idx];

    for (i = 0; i < match_n[idx]; i++) {
      debug(printf("(3) Generating output %u %u\n",outRare,(*freq_ptr) + diagterm_freq));

#ifdef INCLUDE_EXACT
      if (rare_first_p == true) {
	*diagonals++ = outRare;
	*diagonals++ = *(freq_ptr++) + diagterm_freq;
      } else {
	*diagonals++ = *(freq_ptr++) + diagterm_freq;
	*diagonals++ = outRare;
      }
#else
      if (valRare != (valFreq = *freq_ptr++)) {
	if (rare_first_p == true) {
	  *diagonals++ = outRare;
	  *diagonals++ = valFreq + diagterm_freq;
	} else {
	  *diagonals++ = valFreq + diagterm_freq;
	  *diagonals++ = outRare;
	}
      }
#endif
    }

    /* Subsequent iterations on freq */
    while (freq + kProbe < stopFreq && /*minFreq*/freq[kProbe] <= valRare + above_slop) {
      freq += kProbe;
      
#ifdef HAVE_AVX512
      F = _mm512_loadu_si512((const __m512i *)(freq));
      F = _mm512_sub_epi32(F, _epi32_offset);
    
      M_above = _mm512_cmpgt_epi32_mask(F, Rare_low);
      M_below = _mm512_cmplt_epi32_mask(F, Rare_high);
      idx = M_above & M_below;
#elif defined(HAVE_AVX2)
      F = _mm256_loadu_si256((const __m256i *)(freq));
      F = _mm256_sub_epi32(F, _epi32_offset);
    
      M_above = _mm256_cmpgt_epi32(F, Rare_low);
      M_below = _mm256_cmpgt_epi32(Rare_high, F); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);

      idx = _mm256_movemask_ps((__m256) M);
#else
      F0 = LOAD_SI128((const __m128i *)(freq));
      F0 = _mm_sub_epi32(F0, _epi32_offset);
      F1 = LOAD_SI128((const __m128i *)(freq + 4));
      F1 = _mm_sub_epi32(F1, _epi32_offset);

      M_above = _mm_cmpgt_epi32(F1, Rare_low);
      M_below = _mm_cmplt_epi32(F1, Rare_high);
      M = _mm_and_si128(M_above,M_below);
      idx = _mm_movemask_ps((__m128) M);

      M_above = _mm_cmpgt_epi32(F0, Rare_low);
      M_below = _mm_cmplt_epi32(F0, Rare_high);
      M = _mm_and_si128(M_above,M_below);
      idx = (idx << 4) + _mm_movemask_ps((__m128) M);
#endif

      debug(printf("Got idx %d => %d x %d\n",idx,match_start[idx],match_n[idx]));

      if (diagonals + 2*match_n[idx] > (*end_diagonals)) {
	allocation = (*end_diagonals) - (*init_diagonals);
	debug(printf("Reallocating from diagonals %p to end %p (allocation %d)\n",
		     diagonals,*end_diagonals,allocation));
	more_diagonals = (unsigned int *) MALLOC(2 * allocation * sizeof(unsigned int));

	ndiagonals = diagonals - (*init_diagonals);
	debug(printf("Copying from diagonals %p to out %p (ndiagonals %d)\n",
		     diagonals,more_diagonals,ndiagonals));
	memcpy(more_diagonals,*init_diagonals,ndiagonals * sizeof(unsigned int));

	FREE(*init_diagonals);
	*init_diagonals = more_diagonals;
	*end_diagonals = &(more_diagonals[2 * allocation]);
	diagonals = &(more_diagonals[ndiagonals]);
	assert(diagonals + 2*match_n[idx] <= (*end_diagonals));
      }

      const unsigned int *freq_ptr = freq + match_start[idx];
      for (i = 0; i < match_n[idx]; i++) {
	debug(printf("(4) Generating output %u %u\n",outRare,(*freq_ptr) + diagterm_freq));
#ifdef INCLUDE_EXACT
	if (rare_first_p == true) {
	  *diagonals++ = outRare;
	  *diagonals++ = *(freq_ptr++) + diagterm_freq;
	} else {
	  *diagonals++ = *(freq_ptr++) + diagterm_freq;
	  *diagonals++ = outRare;
	}
#else
	if (valRare != (valFreq = *freq_ptr++)) {
	  if (rare_first_p == true) {
	    *diagonals++ = outRare;
	    *diagonals++ = valFreq + diagterm_freq;
	  } else {
	    *diagonals++ = valFreq + diagterm_freq;
	    *diagonals++ = outRare;
	  }
	}
#endif
      }
    }
    /* End of iterations on freq */

    if (freq + kProbe >= stopFreq) {
      freq += kProbe;
      diagonals = scalar_one_set(diagonals, &(*init_diagonals), &(*end_diagonals),
				 rare, freq, endFreq,
				 delta_rare, diagterm_rare, diagterm_freq,
				 below_slop, above_slop, rare_first_p);
    }

    if (++rare < stopRare) {
      valRare = (*rare) + delta_rare; /* For comparison */
      outRare = (*rare) + diagterm_rare; /* For output */

#ifdef HAVE_AVX512
      Rare_low = _mm512_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_high = _mm512_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_low = _mm256_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_high = _mm256_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
      F = F_save;
#else
      Rare_low = _mm_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_high = _mm_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
      F0 = F0_save;
      F1 = F1_save;
#endif

      /* Advance freq */
      freq = freq_save;

      npositions_rare = endRare - rare;
      npositions_freq = endFreq - freq;

      if (50 * npositions_rare > npositions_freq) {
	/* Use linear */
	while (freq < stopFreq && /*maxFreq*/freq[kProbe - 1] + below_slop < valRare) {
	  debug(printf("Advancing freq by %d because %u + below_slop < %u\n",
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

  return scalar(diagonals, &(*init_diagonals), &(*end_diagonals),
		rare, endRare, freq, endFreq,
		delta_rare, diagterm_rare, diagterm_freq,
		below_slop, above_slop, rare_first_p);
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
Intersect_approx_uint4 (int *ndiagpairs,
			const unsigned int *set1, const int length1, int diagterm1,
			const unsigned int *set2, const int length2, int diagterm2,
			const unsigned int below_slop, const unsigned int above_slop) {
  unsigned int *diagonals, *init_diagonals, *end_diagonals;
  int allocation;

  debug(printf("Entered Intersect_approx_uint4 with %d and %d items\n",length1,length2));

#ifdef DEBUG0
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
    *ndiagpairs = 0;
    debug(printf("Exiting Intersect_approx with 0 items\n"));
    return (unsigned int *) NULL;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      if ((allocation = 2 * length1) < 2*NELTS) {
	/* Needs to be greater than 2*match_n */
	allocation = 2*NELTS;
      }
      diagonals = init_diagonals = (unsigned int *) MALLOC(allocation * sizeof(unsigned int));
      end_diagonals = &(diagonals[allocation]);
      diagonals = v1(diagonals, &init_diagonals, &end_diagonals,
		     set1, length1, set2, length2,
		     /*delta_rare*/(diagterm1 - diagterm2),
		     /*diagterm_rare*/diagterm1, /*diagterm_freq*/diagterm2,
		     below_slop, above_slop, /*rare_first_p*/true);
    } else {
      if ((allocation = 2 * length2) < 2*NELTS) {
	allocation = 2*NELTS;
      }
      diagonals = init_diagonals = (unsigned int *) MALLOC(allocation * sizeof(unsigned int));
      end_diagonals = &(diagonals[allocation]);
      diagonals = v1(diagonals, &init_diagonals, &end_diagonals,
		     set2, length2, set1, length1,
		     /*delta_rare*/(diagterm2 - diagterm1),
		     /*diagterm_rare*/diagterm2, /*diagterm_freq*/diagterm1,
		     above_slop, below_slop, /*rare_first_p*/false);
    }
    assert(diagonals <= end_diagonals);

    *ndiagpairs = (diagonals - init_diagonals) / 2;

#ifdef DEBUG
    int k;
    printf("Exiting Intersect_approx with %d items\n",*ndiagpairs);
    for (i = 0, k = 0; i < *ndiagpairs; i++, k += 2) {
      printf("%u %u\n",init_diagonals[k],init_diagonals[k+1]);
    }
#endif

    return init_diagonals;
  }
}

#else

unsigned int *
Intersect_approx_uint4 (int *ndiagpairs,
			const unsigned int *set1, const int length1, int diagterm1,
			const unsigned int *set2, const int length2, int diagterm2,
			const unsigned int below_slop, const unsigned int above_slop) {
  unsigned int *diagonals, *init_diagonals, *end_diagonals;
  int allocation;

  if ((length1 == 0) || (length2 == 0)) {
    *ndiagpairs = 0;
    return (unsigned int *) NULL;

  } else {
    if (length1 <= length2) {
      allocation = 2 * length1;
      diagonals = init_diagonals = (unsigned int *) MALLOC(allocation * sizeof(unsigned int));
      end_diagonals = &(diagonals[allocation]);
      diagonals = scalar(diagonals, &init_diagonals, &end_diagonals,
			 set1, /*endA*/&(set1[length1]),
			 set2, /*endB*/&(set2[length2]),
			 /*deltaA*/(diagterm1 - diagterm2),
			 /*diagtermA*/diagterm1, /*diagtermB*/diagterm2,
			 below_slop, above_slop, /*A_first_p*/true);
    } else {
      allocation = 2 * length2;
      diagonals = init_diagonals = (unsigned int *) MALLOC(allocation * sizeof(unsigned int));
      end_diagonals = &(diagonals[allocation]);
      diagonals = scalar(diagonals, &init_diagonals, &end_diagonals,
			 set2, /*endA*/&(set2[length2]),
			 set1, /*endB*/&(set1[length1]),
			 /*deltaA*/(diagterm2 - diagterm1),
			 /*diagtermA*/diagterm2, /*diagtermB*/diagterm1,
			 above_slop, below_slop, /*A_first_p*/false);
    }

    *ndiagpairs = (diagonals - init_diagonals) / 2;
    return init_diagonals;
  }
}
#endif


void
Intersect_approx_uint4_setup () {
  initialize_match_values();
  return;
}  
