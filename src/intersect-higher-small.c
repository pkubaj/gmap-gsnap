static char rcsid[] = "$Id: 0d783f52bdb132d49496468133540e5be9681f2c $";
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

#include "assert.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include "mem.h"
#include "bool.h"
#include "intersect-higher-small.h"

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
static unsigned int *
scalar (unsigned int *out, unsigned int *diagonals,
	const unsigned int *A, const unsigned int *endA,
	const unsigned int *B, const unsigned int *endB, int deltaA,
	int diagtermA, int diagtermB,
	const unsigned int below_slop, const unsigned int above_slop,
	bool A_first_p) {

  const unsigned int *B_save;

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
      if ((*A) + deltaA == (*B)) {
	/* Exclude exact */
      } else if (A_first_p == true) {
	if (1 || (*A) + diagtermA >= (*B) + diagtermB) { /* higher */
	  if (out == diagonals || (*A) + diagtermA > out[-1]) {
	    *out++ = (*A) + diagtermA;
	  }
	}
      } else {
	if (1 || (*B) + diagtermB >= (*A) + diagtermA) { /* higher */
	  if (out == diagonals || (*B) + diagtermB > out[-1]) {
	    *out++ = (*B) + diagtermB;
	  }
	}
      }
	  
      B++;
    }
    B = B_save;

    A++;
  }

  return out;
}


/* Performs scalar against a single value of A */
static inline unsigned int *
scalar_one_set (unsigned int *out, unsigned int *diagonals,
		const unsigned int *A, const unsigned int *B, const unsigned int *endB,
		int deltaA, int diagtermA, int diagtermB,
		const unsigned int below_slop, const unsigned int above_slop,
		bool A_first_p) {

  debug(printf("Entered scalar_one_set with %lu entries\n",endB - B));

  debug(printf("Advancing B from %u up to %u\n",*B,(*A) + deltaA - below_slop));
  while (B < endB && (*B) + below_slop < (*A) + deltaA) {
    B++;
  }

  while (B < endB && (*B) <= (*A) + deltaA + above_slop) {
    debug(printf("(2) Generating output %u %u\n",(*A) + diagtermA,(*B) + diagtermB));

    if ((*A) + deltaA == (*B)) {
      /* Exclude exact */
    } else if (A_first_p == true) {
      if (1 || (*A) + diagtermA >= (*B) + diagtermB) { /* higher */
	if (out == diagonals || (*A) + diagtermA > out[-1]) {
	  *out++ = (*A) + diagtermA;
	}
      }
    } else {
      if (1 || (*B) + diagtermB >= (*A) + diagtermA) { /* higher */
	if (out == diagonals || (*B) + diagtermB > out[-1]) {
	  *out++ = (*B) + diagtermB;
	}
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


#ifdef HAVE_SSE2
static unsigned int *
v1 (unsigned int *out, unsigned int *diagonals,
    const unsigned int *rare, int lenRare,
    const unsigned int *freq, int lenFreq, int delta_rare,
    int diagterm_rare, int diagterm_freq,
    const unsigned int slop, const unsigned int insertion_slop, bool rare_first_p) {

  const unsigned int *initRare = rare;
  const unsigned int *initFreq = freq;
  const unsigned int *endRare = &(rare[lenRare]);
  const unsigned int *endFreq = &(freq[lenFreq]);

  const int kRareSpace = 0;	/* 0. type was uint64_t */

  unsigned int valRare, outRare, outFreq;
  unsigned int goal;
  int npositions_rare, npositions_freq, j;

  int i;
  unsigned int below_slop, above_slop;
  unsigned int below_slop_plus_1, above_slop_plus_1; /* Because only cmpgt is available, not cmpge */

  /* higher */
  if (rare_first_p == true) {
    below_slop = slop;
    above_slop = insertion_slop;
  } else {
    below_slop = insertion_slop;
    above_slop = slop;
  }

  below_slop_plus_1 = below_slop + 1; 
  above_slop_plus_1 = above_slop + 1;

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

  /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and in main loop */
  /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */

  /* Avoid issues where valRare - slop_plus_1 < 0 */
  while (COMPILER_RARELY(rare < stopRare &&
			 (valRare = (*rare) + delta_rare) <= slop)) {
    out = scalar_one_set(out, diagonals, rare, freq, endFreq,
			 delta_rare, diagterm_rare, diagterm_freq,
			 below_slop, above_slop, rare_first_p);
    rare++;
    /* valRare = (*rare) + delta_rare; -- For comparison, now computed inside condition and in main loop */
    /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */
  }

  if (COMPILER_LIKELY(rare < stopRare)) {
    valRare = (*rare) + delta_rare;
    outRare = (*rare) + diagterm_rare;

    /* Advance freq */
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
  }

  if (freq < stopFreq) {
#ifdef HAVE_AVX512
    Rare_lowbound = _mm512_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_highbound = _mm512_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F = _mm512_loadu_si512((const __m512i *)(freq));
    F = _mm512_sub_epi32(F, _epi32_offset);
#elif defined(HAVE_AVX2)
    Rare_lowbound = _mm256_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_highbound = _mm256_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F = _mm256_loadu_si256((const __m256i *)(freq));
    F = _mm256_sub_epi32(F, _epi32_offset);
#else
    Rare_lowbound = _mm_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
    Rare_highbound = _mm_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);

    F0 = LOAD_SI128((const __m128i *)(freq));
    F0 = _mm_sub_epi32(F0, _epi32_offset);
    F1 = LOAD_SI128((const __m128i *)(freq + 4));
    F1 = _mm_sub_epi32(F1, _epi32_offset);
#endif
  }

  while (rare < stopRare && freq < stopFreq) {
    debug(printf("rare is at item %lu with value %u\n",rare - initRare,valRare));
    debug(printf("freq is at item %lu with values %u..%u\n",
		 freq - initFreq,freq[0],freq[kProbe - 1]));

    const unsigned int *freq_save = freq;

    /* First iteration on freq */
#ifdef HAVE_AVX512
    F_save = F;

    M_above = _mm512_cmpgt_epi32_mask(F, Rare_lowbound); /* 16 comparisons */
    M_below = _mm512_cmplt_epi32_mask(F, Rare_highbound);
    idx = M_above & M_below;
#elif defined(HAVE_AVX2)
    F_save = F;

    debug(print_vector_unsigned(Rare_lowbound));
    debug(print_vector_unsigned(F));
    debug(print_vector_unsigned(Rare_highbound));

    M_above = _mm256_cmpgt_epi32(F, Rare_lowbound); /* 8 comparisons */
    M_below = _mm256_cmpgt_epi32(Rare_highbound, F); /* cmplt not available */
    M = _mm256_and_si256(M_above,M_below);

    debug(print_vector_hex(M_above));
    debug(print_vector_hex(M_below));
    debug(print_vector_hex(M));

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

    const unsigned int *freq_ptr = freq + match_start[idx];
    for (i = 0; i < match_n[idx]; i++) {
      outFreq = *(freq_ptr++) + diagterm_freq;
      debug(printf("(3) Generating output rare %u, freq %u\n",outRare,outFreq));
      if (outRare == outFreq) {
	/* Exclude exact */
      } else if (rare_first_p == true) {
	if (out == diagonals || outRare > out[-1]) {
	  *out++ = outRare;
	}
      } else {
	if (out == diagonals || outFreq > out[-1]) {
	  *out++ = outFreq;
	}
      }
    }

    /* Subsequent iterations on freq */
    while (freq + kProbe < stopFreq && /*minFreq*/freq[kProbe] <= valRare + above_slop) {
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
    
      debug(print_vector_unsigned(Rare_lowbound));
      debug(print_vector_unsigned(F));
      debug(print_vector_unsigned(Rare_highbound));

      M_above = _mm256_cmpgt_epi32(F, Rare_lowbound); /* 8 comparisons */
      M_below = _mm256_cmpgt_epi32(Rare_highbound, F); /* cmplt not available */
      M = _mm256_and_si256(M_above,M_below);

      debug(print_vector_hex(M_above));
      debug(print_vector_hex(M_below));
      debug(print_vector_hex(M));

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

      const unsigned int *freq_ptr = freq + match_start[idx];
      for (i = 0; i < match_n[idx]; i++) {
	outFreq = *(freq_ptr++) + diagterm_freq;
	debug(printf("(4) Generating output rare %u, freq %u\n",outRare,outFreq));
	if (outRare == outFreq) {
	  /* Exclude exact */
	} else if (rare_first_p == true) {
	  if (out == diagonals || outRare > out[-1]) {
	    *out++ = outRare;
	  }
	} else {
	  if (out == diagonals || outFreq > out[-1]) {
	    *out++ = outFreq;
	  }
	}
      }
    }
    /* End of iterations on freq */

    if (freq + kProbe >= stopFreq) {
      freq += kProbe;
      out = scalar_one_set(out, diagonals, rare, freq, endFreq,
			   delta_rare, diagterm_rare, diagterm_freq,
			   below_slop, above_slop, rare_first_p);
    }

    if (++rare < stopRare) {
      valRare = (*rare) + delta_rare; /* For comparison */
      /* outRare = (*rare) + diagterm_rare; -- For output, not needed here */

#ifdef HAVE_AVX512
      Rare_lowbound = _mm512_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_highbound = _mm512_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
      F = F_save;
#elif defined(HAVE_AVX2)
      Rare_lowbound = _mm256_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_highbound = _mm256_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
      F = F_save;
#else
      Rare_lowbound = _mm_set1_epi32(valRare - below_slop_plus_1 - EPI32_MAX);
      Rare_highbound = _mm_set1_epi32(valRare + above_slop_plus_1 - EPI32_MAX);
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

  return scalar(out, diagonals, rare, endRare, freq, endFreq,
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
int
Intersect_higher (Univcoord_T *diagonals,
		  const UINT4 *set1, const int length1, int diagterm1,
		  const Univcoord_T *set2, const int length2,
		  const Univcoord_T slop, const Univcoord_T insertion_slop) {
  Univcoord_T *out = diagonals;
#ifdef DEBUG15
  unsigned int *nonsimd_diagonals;
  int nonsimd_ndiagonals, k;
#endif

  debug(printf("Entered Intersect_higher with %d and %d items\n",length1,length2));

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
    debug(printf("Exiting Intersect_higher with 0 items\n"));
    return 0;

  } else {
    debug(printf("Using v1 method\n"));
    if (length1 <= length2) {
      out = v1(out, diagonals, set1, length1, set2, length2,
	       /*delta_rare*/+diagterm1,
	       /*diagterm_rare*/diagterm1, /*diagterm_freq*/0,
	       slop, insertion_slop, /*rare_first_p*/true);
    } else {
      out = v1(out, diagonals, set2, length2, set1, length1,
	       /*delta_rare*/-diagterm1,
	       /*diagterm_rare*/0, /*diagterm_freq*/diagterm1,
	       slop, insertion_slop, /*rare_first_p*/false);
    }

    debug(printf("Exiting Intersect_higher with %d items\n",(int) (out - diagonals)));
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
						     set1,length1,diagterm1,
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
		  const UINT4 *set1, const int length1, int diagterm1,
		  const Univcoord_T *set2, const int length2,
		  const Univcoord_T slop, const Univcoord_T insertion_slop) {
  Univcoord_T *out = diagonals;

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else {
    if (length1 <= length2) {
      out = scalar(out, diagonals, set1, /*endA*/&(set1[length1]),
		   set2, /*endB*/&(set2[length2]),
		   /*deltaA*/+diagterm1,
		   /*diagtermA*/diagterm1, /*diagtermB*/0,
		   /*below_slop*/slop, /*above_slop*/insertion_slop, /*A_first_p*/true);
    } else {
      out = scalar(out, diagonals, set2, /*endA*/&(set2[length2]),
		   set1, /*endB*/&(set1[length1]),
		   /*deltaA*/-diagterm1,
		   /*diagtermA*/0, /*diagtermB*/diagterm1,
		   /*below_slop*/insertion_slop, /*above_slop*/slop, /*A_first_p*/false);
    }

    return (out - diagonals);
  }
}
#endif


#ifdef DEBUG15
/* LARGE_GENOMES still needs this to handle transcriptome */
/* Returns results as pairs of coordinates */
/* diagonals is already allocated by caller */
/* Needs to return diagonals in ascending order.  Need to check if diagonals0 has duplicates. */
/* It appears that this procedure may not be accurate */
int
Intersect_approx_higher_old (unsigned int *diagonals,
			     unsigned int *positions1, int npositions1, int diagterm1,
			     unsigned int *positions0, int npositions0,
			     const unsigned int maxdistance) {
  int ndiagonals = 0;
  unsigned int local_goal, last_positions0;
  unsigned int diagterm, delta;
  unsigned int last_diagonal, this_diagonal;
  UINT4 *ptr1, *start1;
  int j;

  diagterm = (unsigned int) diagterm1;	/* local_goal based on larger list */
  delta = (unsigned int) (-diagterm1); /* list0 + (diagterm0 - diagterm1) = list1 */

  start1 = positions1;
  /* end1 = &(positions1[npositions1]); */

  debug9(printf("Intersect_approx_higher with %d positions against %d positions.  diagterm %d\n",
	       npositions1,npositions0,(int) diagterm));
#ifdef DEBUG9
  for (ptr1 = positions0; ptr1 < &(positions0[npositions0]); ptr1++) {
    printf("%llu\n",*ptr1);
  }
  printf("\n");

  for (ptr1 = positions1; ptr1 < &(positions1[npositions1]); ptr1++) {
    printf("%llu + %d\n",*ptr1,diagterm1);
  }
  printf("\n");
#endif


#if 0
  /* Should not be necessary, since this procedure is called after Indexdb_ptr_with_diagterm */
  if (diagterm1 < 0) {
    while (npositions1 > 0 && (*positions1) < (unsigned int) -diagterm1) {
      ++positions1;
      --npositions1;
    }
  }
#endif

  last_positions0 = (unsigned int) -1;
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
#ifdef DEBUG9
      if (npositions1 > 0) {
	printf("Result of search is npositions1 %d, pointing to %llu\n",npositions1,*positions1);
      }
#endif
      
      if (npositions1 <= 0) {
	/* Check backwards only */
	debug9(printf("    intersection list 1 at end  checking for approximate:"));
	ptr1 = &(positions1[-1]);
	if (ptr1 >= start1) {
	  debug9(printf(" prev %d:%llu?",npositions1-1,*ptr1));
	  if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	    debug9(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	debug9(printf("\n"));
	
	/* Should be able to return at this point.  positions1[n-1] must be < *positions0 */
#ifdef DEBUG9
	printf("Returning %d diagonals\n",ndiagonals);
	for (j = 0; j < ndiagonals; j++) {
	  printf("%u\n",diagonals[j]);
	}
#endif
	return ndiagonals;
	
      } else if ((*positions1) == local_goal) {
	/* Found local goal.  Save and advance */
	debug9(printf("    intersection list 1: %d:%llu  exact\n",npositions1,*positions1));
	/* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	if (ndiagonals == 0) {
	  last_diagonal = diagonals[ndiagonals++] = local_goal + diagterm;
	} else if ((this_diagonal = local_goal + diagterm) != last_diagonal) {
	  last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	}
	++positions1;
	--npositions1;
	
      } else {
	debug9(printf("    intersection list 1: %d:%llu  checking for approximate:",npositions1,*positions1));
	ptr1 = &(positions1[-1]); /* closest position < local_goal */
	if (ptr1 >= start1) {
	  debug9(printf(" prev %d:%llu?",npositions1+1,*ptr1));
	  if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	    debug9(printf(" yes."));
	    /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	    if (ndiagonals == 0) {
	      last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	    } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
	      last_diagonal = diagonals[ndiagonals++] = this_diagonal;
	    }
	  }
	}
	
	ptr1 = &(positions1[0]); /* closest position > local_goal */
	debug9(printf(" at %d:%llu?",npositions1,*ptr1));
	if ((*ptr1) /*(omit for higher) + maxdistance*/ >= local_goal && (*ptr1) <= local_goal + maxdistance) {
	  debug9(printf(" yes."));
	  /* diagonals[ndiagonals++] = local_goal + diagterm; -- For diagpairs */
	  if (ndiagonals == 0) {
	    last_diagonal = diagonals[ndiagonals++] = (*ptr1) + diagterm;
	  } else if ((this_diagonal = (*ptr1) + diagterm) != last_diagonal) {
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

#ifdef DEBUG9
  printf("Returning %d diagonals\n",ndiagonals);
  for (j = 0; j < ndiagonals; j++) {
    printf("%u\n",diagonals[j]);
  }
#endif

  return ndiagonals;
}
#endif


void
Intersect_higher_setup () {
  initialize_match_values();
  return;
}  
