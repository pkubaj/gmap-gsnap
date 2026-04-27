static char rcsid[] = "$Id: 81bf15f0282c82638c0eff30f6ee130862579e99 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>		/* For uintptr_t */
#include "assert.h"
#include "mem.h"
#include "intersect-uint2.h"

#include "simd.h"
#include "popcount.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>		/* For _mm_loadu_si128 and _mm_storeu_si128 */
#endif
#ifdef HAVE_SSSE3
#include <tmmintrin.h>		/* For _mm_shuffle_epi8 */
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>		/* For _mm_extract_epi32 */
#endif
#ifdef HAVE_SSE4_2
#include <nmmintrin.h>		/* For _mm_cmpestrm */
#endif
#endif

#define ALIGN 16		/* 128-bit quantities need to be aligned at 16-byte quantities */

/* _mm_lddqu_si128 is better for unaligned memory */
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

#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


#ifdef DEBUG
static void
print_vector (char *label, __m128i x) {
  printf("%s %10u %10u %10u %10u %10u %10u %10u %10u\n",label,
	 _mm_extract_epi16(x,0),
	 _mm_extract_epi16(x,1),
	 _mm_extract_epi16(x,2),
	 _mm_extract_epi16(x,3),
	 _mm_extract_epi16(x,4),
	 _mm_extract_epi16(x,5),
	 _mm_extract_epi16(x,6),
	 _mm_extract_epi16(x,7));
  return;
}
#endif


int
Intersect_uint2_scalar (Univcoord_T *diagonals, unsigned short **set2_ptr,
			unsigned short *set1, const int length1, int diagterm1,
			unsigned short *set2, const int length2, int diagterm2,
			Univcoord_T region_term) {
  Univcoord_T *out;
  int i, j;
  int value1, value2;

  debug15(printf("Entered Intersect_uint2_scalar\n"));

  if ((length1 == 0) || (length2 == 0)) {
    *set2_ptr = &(set2[0]);
    return 0;

  } else {
    out = diagonals;

    i = j = 0;
    while (i < length1 && j < length2) {
      if ((value1 = set1[i] + diagterm1) < (value2 = set2[j] + diagterm2)) {
	i++;
      } else if (value2 < value1) {
	j++;
      } else {
	debug15(printf("Scalar: equal %d + %d == %d + %d\n",set1[i],diagterm1,set2[j],diagterm2));
	*out++ = region_term + (Univcoord_T) value1;
	i++;
	j++;
      }
    }

    *set2_ptr = &(set2[j]);
    return /*ndiagonals*/ (out - diagonals);
  }
}


/* STTNI implies SSE4_2 */
#if defined(HAVE_SSE4_2) && defined(HAVE_STTNI)

/* 256 128-bit registers */
static const char shuffle_mask16[4096] = 
  {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1,
    8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1,
    6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1,
    10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 10, 11, -1, -1, -1, -1, -1, -1,
    8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1,
    6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1,
    12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, -1, -1, -1, -1, -1, -1,
    8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1,
    6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, -1, -1, -1, -1,
    10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, -1, -1, -1, -1,
    8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1,
    6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, -1, -1,
    14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 14, 15, -1, -1, -1, -1, -1, -1,
    8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1,
    6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, -1, -1, -1, -1,
    10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 14, 15, -1, -1, -1, -1,
    8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1,
    6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, -1, -1,
    12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    4, 5, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 12, 13, 14, 15, -1, -1, -1, -1,
    8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1,
    6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1,
    4, 5, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, -1, -1,
    10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    2, 3, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    4, 5, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 4, 5, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    4, 5, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, -1, -1,
    8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1,
    0, 1, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    2, 3, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    4, 5, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1,
    6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1,
    0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1,
    0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};


/**
 * From Schlegel et al., Fast Sorted-Set Intersection using SIMD Instructions
 *
 * Optimized by D. Lemire on May 3rd 2013
 */
/* Assume that A is aligned, but B is not */
static int
intersect_vector16 (uint16_t *A, uint16_t deltaA, uint16_t **Bptr, uint16_t *B,
		    const int s_a, const int s_b, uint16_t *C) {
  int count = 0;
  int i_a = 0, i_b = 0;
  uint16_t a_adjusted;

  const int st_a = (s_a / 8) * 8; /* stop index: number of 8-element blocks */
  const int st_b = (s_b / 8) * 8; /* stop index: number of 8-element blocks */
  __m128i v_a, v_b, v_delta_a;

  if ((i_a < st_a) && (i_b < st_b)) {
    v_delta_a = _mm_set1_epi16(deltaA);
    v_a = _mm_add_epi16(LOAD_SI128((const __m128i *) &A[i_a]), v_delta_a);
    v_b = LOAD_SI128((const __m128i *) &B[i_b]);
#ifdef DEBUG
    print_vector("v_a",v_a);
    print_vector("v_b",v_b);
    printf("\n");
#endif    

    while (((unsigned short) (A[i_a] + deltaA) == 0) || (B[i_b] == 0)) {
      const __m128i res_v =
          _mm_cmpestrm(v_b, 8, v_a, 8,
                       _SIDD_UWORD_OPS | _SIDD_CMP_EQUAL_ANY | _SIDD_BIT_MASK);

      const int r = _mm_extract_epi32(res_v, 0); /* r must be in 0..255 */
      __m128i p = _mm_shuffle_epi8(v_a, * (__m128i *) &(shuffle_mask16[r*16]));
      _mm_storeu_si128((__m128i *) &C[count], p);

      count += popcount_ones_32(r);
      const uint16_t a_max = A[i_a + 7] + deltaA;
      const uint16_t b_max = B[i_b + 7];
      if (a_max <= b_max) {
        i_a += 8;
        if (i_a == st_a) {
          break;
	}
        v_a = _mm_add_epi16(LOAD_SI128((const __m128i *) &A[i_a]), v_delta_a);
      }

      if (b_max <= a_max) {
        i_b += 8;
        if (i_b == st_b) {
          break;
	}
        v_b = LOAD_SI128((const __m128i *) &B[i_b]);
      }

#ifdef DEBUG
      print_vector("v_a",v_a);
      print_vector("v_b",v_b);
      printf("\n");
#endif    
    }

    if ((i_a < st_a) && (i_b < st_b))
      while (true) {
        const __m128i res_v = _mm_cmpistrm(
            v_b, v_a, _SIDD_UWORD_OPS | _SIDD_CMP_EQUAL_ANY | _SIDD_BIT_MASK);

        const int r = _mm_extract_epi32(res_v, 0); /* r must be in 0..255 */
        __m128i p = _mm_shuffle_epi8(v_a, * (__m128i *) &(shuffle_mask16[r*16]));
        _mm_storeu_si128((__m128i *) &C[count], p);

        count += popcount_ones_32(r);
        const uint16_t a_max = A[i_a + 7] + deltaA;
        const uint16_t b_max = B[i_b + 7];
        if (a_max <= b_max) {
          i_a += 8;
          if (i_a == st_a) {
            break;
	  }
          v_a = _mm_add_epi16(LOAD_SI128((const __m128i *) &A[i_a]), v_delta_a);
        }
        if (b_max <= a_max) {
          i_b += 8;
          if (i_b == st_b) {
            break;
	  }
          v_b = LOAD_SI128((const __m128i *) &B[i_b]);
        }
#ifdef DEBUG
	print_vector("v_a",v_a);
	print_vector("v_b",v_b);
	printf("\n");
#endif    
      }
  }

  /* intersect the tail using scalar intersection */
  while (i_a < s_a && i_b < s_b) {
    if ((a_adjusted = A[i_a] + deltaA) < B[i_b]) {
      i_a++;
    } else if (B[i_b] < a_adjusted) {
      i_b++;
    } else {
      C[count] = a_adjusted;
      count++;
      i_a++;
      i_b++;
    }
  }

  *Bptr = &(B[i_b]);
  return count;
}


#elif defined(HAVE_SSE2)

static int
intersect_vector16 (uint16_t *A, uint16_t deltaA, uint16_t **Bptr, uint16_t *B,
		    const size_t s_a, const size_t s_b,
		    uint16_t *C) {
  uint16_t a_adjusted;

  debug(printf("Entered intersect_vector16 without STTNI with %lu and %lu entries\n",
	       s_a,s_b));

#if 0
  if (s_a > s_b)
    return intersect_vector16(B, A, s_b, s_a, C);
#endif

  size_t count = 0;
  size_t i_a = 0, i_b = 0;
  const size_t st_a = s_a;
  const size_t st_b = (s_b / 8) * 8;
  __m128i v_a, v_b;

  if ((i_a < st_a) && (i_b < st_b)) {
    v_a = _mm_set1_epi16(A[i_a] + deltaA);

    v_b = LOAD_SI128((const __m128i *)&B[i_b]);
    while (B[i_b + 7] < (unsigned short) (A[i_a] + deltaA)) {
      debug(printf("Advancing B because %hu < %hu (%hu + %hd)\n",
		   B[i_b + 7],A[i_a] + deltaA,A[i_a],deltaA));
      i_b += 8;
      if (i_b == st_b) {
	debug(printf("Stopping at end of B because %lu == stop %lu\n",i_b,st_b));
	goto FINISH_SCALAR;
      } else {
	v_b = LOAD_SI128((const __m128i *)&B[i_b]);
      }
    }

    while (true) {
      const __m128i F0 = _mm_cmpeq_epi16(v_a, v_b);
#ifdef DEBUG
      print_vector("v_a",v_a);
      print_vector("v_b",v_b);
      printf("\n");
#endif
#ifdef HAVE_SSE4_1
      if (_mm_testz_si128(F0, F0)) {
#else
      if (!_mm_movemask_epi8(F0)) {
#endif
      } else {
        C[count] = A[i_a] + deltaA;
	debug(printf("Vector output: %u = %u + %u\n",C[count],A[i_a],deltaA));
        count++;
      }

      ++i_a;
      if (i_a == st_a) {
	debug(printf("Stopping at end of A because %lu == stop %lu\n",i_a,st_a));
        goto FINISH_SCALAR;
      } else {
	v_a = _mm_set1_epi16(A[i_a] + deltaA);
      }

      while (B[i_b + 7] < (unsigned short) (A[i_a] + deltaA)) {
	debug(printf("Advancing B because %hu < %hu (%hu + %hd)\n",
		     B[i_b + 7],A[i_a] + deltaA,A[i_a],deltaA));
        i_b += 8;
        if (i_b == st_b) {
	  debug(printf("Stopping at end of B because %lu == stop %lu\n",i_b,st_b));
          goto FINISH_SCALAR;
	} else {
	  v_b = LOAD_SI128((const __m128i *)&B[i_b]);
	}
      }
    }
  }
    FINISH_SCALAR:
  // intersect the tail using scalar intersection
  while (i_a < s_a && i_b < s_b) {
    if ((a_adjusted = A[i_a] + deltaA) < B[i_b]) {
      i_a++;
    } else if (B[i_b] < a_adjusted) {
      i_b++;
    } else {
      C[count] = a_adjusted;
      count++;
      i_a++;
      i_b++;
    }
  }

  *Bptr = &(B[i_b]);
  return count;
}

#endif


 
#ifdef HAVE_SSE2
int
Intersect_uint2 (Univcoord_T *diagonals, unsigned short *localdb_alloc,
		 unsigned short *set1, const int length1, int diagterm1,
		 unsigned short *set2, const int length2, int diagterm2,
		 Univcoord_T region_term) {
  int ndiagonals, ndiagonals_vector, unaligned_length, vector_length, scalar_length;
  Univcoord_T *out;
  unsigned short *aligned, *intersection, *aligned_ptr;
  uint16_t *longset_ptr;
  unsigned short delta, diagterm;
  int i;

#ifdef DEBUG15
  Univcoord_T diagonals_alloc_scalar[65536];
  int ndiagonals_scalar;
#endif

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else if (length1 <= length2) {
    /* set1 is short/aligned; set2 is long */
#ifdef DEBUG
    printf("Set 1 (shorter)\n");
    for (i = 0; i < length1; i++) {
      printf("%u\n",set1[i]);
    }
    
    printf("Set 2 (longer)\n");
    for (i = 0; i < length2; i++) {
      printf("%u\n",set2[i]);
    }
#endif

    if (diagterm1 >= diagterm2) {
      delta = (unsigned short) (diagterm1 - diagterm2); /* set1 + diagterm1 - diagterm2 == set2 */
      debug(printf("delta is +%u\n",delta));

      /* Start with unaligned set1 */
      aligned_ptr = (unsigned short *) (((uintptr_t) set1 + (uintptr_t) (ALIGN - 1)) &~ (uintptr_t) (ALIGN - 1));
      if ((unaligned_length = aligned_ptr - set1) > length1) {
	unaligned_length = length1;
      }
      debug(printf("unaligned length is %d\n",unaligned_length));
      ndiagonals = Intersect_uint2_scalar(diagonals,&longset_ptr,set1,unaligned_length,diagterm1,
					  set2,length2,diagterm2,region_term);
      out = &(diagonals[ndiagonals]);

      /* No overflow by delta: Continue with vector */
      aligned = &(set1[unaligned_length]);
      vector_length = length1 - unaligned_length;
      while (vector_length - 1 >= 0 && aligned[vector_length - 1] >= 65336 - delta) {
	vector_length--;
      }
    
      intersection = &(localdb_alloc[0]); /* Storage for intersecting values */
      ndiagonals +=
	(ndiagonals_vector =
	 intersect_vector16(/*A*/aligned,/*deltaA*/delta,/*Bptr*/&longset_ptr,/*B*/set2,
			    /*lengthA*/vector_length,/*lengthB*/length2,intersection));

      diagterm = (unsigned short) diagterm2; /* result = set2 + diagterm2 */
      for (i = 0; i < ndiagonals_vector; i++) {
	debug15(printf("(Case 1) Intersection %d.  Writing %u + %d + %d\n",
		       intersection[i],region_term,intersection[i],diagterm));
	*out++ = region_term + intersection[i] + diagterm;
      }

      /* Overflow by delta: Finish with scalar */
      ndiagonals += Intersect_uint2_scalar(out,&longset_ptr,
					   &(set1[unaligned_length + vector_length]),
					   length1 - (unaligned_length + vector_length),diagterm1,
					   longset_ptr,length2 - (longset_ptr - set2),diagterm2,region_term);

    } else {
      /* This is really -delta */
      delta = (unsigned short) (diagterm2 - diagterm1); /* set1 + diagterm1 - diagterm2 == set2 */
      debug(printf("delta is -%u\n",delta));

      /* Overflow by delta: Start with scalar */
      i = 0;
      while (i < length1 && set1[i] < delta) {
	i++;
      }
      aligned_ptr = (unsigned short *) (((uintptr_t) (set1 + i) + (uintptr_t) (ALIGN - 1)) &~ (uintptr_t) (ALIGN - 1));
      if ((scalar_length = aligned_ptr - set1) > length1) {
	scalar_length = length1;
      }
      debug(printf("scalar length is %d\n",scalar_length));

      ndiagonals = Intersect_uint2_scalar(diagonals,&longset_ptr,set1,scalar_length,diagterm1,
					  set2,length2,diagterm2,region_term);
      out = &(diagonals[ndiagonals]);

      /* No overflow by delta: Finish with vector */
      aligned = &(set1[scalar_length]);
      intersection = &(localdb_alloc[0]); /* Storage for intersecting values */
      ndiagonals +=
	(ndiagonals_vector =
	 intersect_vector16(/*A*/aligned,/*deltaA*/-delta,/*Bptr*/&longset_ptr,/*B*/longset_ptr,
			    /*lengthA*/length1 - scalar_length,/*lengthB*/length2 - (longset_ptr - set2),
			    intersection));

      diagterm = (unsigned short) diagterm2; /* result = set2 + diagterm2 */
      for (i = 0; i < ndiagonals_vector; i++) {
	debug15(printf("(Case 2) Intersection %d.  Writing %u + %d + %d\n",
		       intersection[i],region_term,intersection[i],diagterm));
	*out++ = region_term + intersection[i] + diagterm;
      }
    }

#ifdef DEBUG15
    ndiagonals_scalar = Intersect_uint2_scalar(diagonals_alloc_scalar,&longset_ptr,
					       set1,length1,diagterm1,set2,length2,diagterm2,region_term);
    if (ndiagonals != ndiagonals_scalar) {
      printf("%d vs %d\n",ndiagonals,ndiagonals_scalar);
      abort();
    } else {
      for (i = 0; i < ndiagonals; i++) {
	if (diagonals[i] != diagonals_alloc_scalar[i]) {
	  for (i = 0; i < ndiagonals; i++) {
	    printf("%u %u\n",diagonals[i],diagonals_alloc_scalar[i]);
	  }
	  abort();
	}
      }
    }
#endif

    return ndiagonals;

  } else {
    /* set2 is short/aligned; set1 is long */
#ifdef DEBUG
    printf("Set 2 (shorter)\n");
    for (i = 0; i < length2; i++) {
      printf("%u\n",set2[i]);
    }

    printf("Set 1 (longer)\n");
    for (i = 0; i < length1; i++) {
      printf("%u\n",set1[i]);
    }
#endif

    if (diagterm2 >= diagterm1) {
      delta = (unsigned short) (diagterm2 - diagterm1); /* set1 == set2 + diagterm2 - diagterm1 */
      debug(printf("delta is +%u\n",delta));
      
      /* Start with unaligned set2 */
      aligned_ptr = (unsigned short *) (((uintptr_t) set2 + (uintptr_t) (ALIGN - 1)) &~ (uintptr_t) (ALIGN - 1));
      if ((unaligned_length = aligned_ptr - set2) > length2) {
	unaligned_length = length2;
      }
      ndiagonals = Intersect_uint2_scalar(diagonals,&longset_ptr,set2,unaligned_length,diagterm2,
					  set1,length1,diagterm1,region_term);
      out = &(diagonals[ndiagonals]);

      /* No overflow by delta: Continue with vector */
      aligned = &(set2[unaligned_length]);
      vector_length = length2 - unaligned_length;
      while (vector_length - 1 >= 0 && aligned[vector_length - 1] >= 65536 - delta) {
	vector_length--;
      }

      intersection = &(localdb_alloc[0]); /* Storage for intersecting values */
      ndiagonals +=
	(ndiagonals_vector =
	 intersect_vector16(/*A*/aligned,/*deltaA*/delta,/*Bptr*/&longset_ptr,/*B*/set1,
			    /*lengthA*/vector_length,/*lengthB*/length1,intersection));

      diagterm = (unsigned short) diagterm1; /* result = set1 + diagterm1 */
      for (i = 0; i < ndiagonals_vector; i++) {
	debug15(printf("(Case 3) Intersection %d.  Writing %u + %d + %d\n",
		       intersection[i],region_term,intersection[i],diagterm));
	*out++ = region_term + intersection[i] + diagterm;
      }

      /* Overflow by delta: Finish with scalar */
      ndiagonals += Intersect_uint2_scalar(out,&longset_ptr,
					   &(set2[unaligned_length + vector_length]),
					   length2 - (unaligned_length + vector_length),diagterm2,
					   longset_ptr,length1 - (longset_ptr - set1),diagterm1,region_term);

    } else {
      /* This is really -delta */
      delta = (unsigned short) (diagterm1 - diagterm2); /* set1 == set2 + diagterm2 - diagterm1 */
      debug(printf("delta is -%u\n",delta));

      /* Overflow by delta: Start with scalar */
      i = 0;
      while (i < length2 && set2[i] < delta) {
	i++;
      }
      aligned_ptr = (unsigned short *) (((uintptr_t) (set2 + i) + (uintptr_t) (ALIGN - 1)) &~ (uintptr_t) (ALIGN - 1));
      if ((scalar_length = aligned_ptr - set2) > length2) {
	scalar_length = length2;
      }
      debug(printf("scalar length is %d\n",scalar_length));

      ndiagonals = Intersect_uint2_scalar(diagonals,&longset_ptr,set2,scalar_length,diagterm2,
					  set1,length1,diagterm1,region_term);
      out = &(diagonals[ndiagonals]);

      /* No overflow by delta: Finish with vector */
      aligned = &(set2[scalar_length]);
      intersection = &(localdb_alloc[0]); /* Storage for intersecting values */
      ndiagonals +=
	(ndiagonals_vector =
	 intersect_vector16(/*A*/aligned,/*deltaA*/-delta,/*Bptr*/&longset_ptr,/*B*/longset_ptr,
			    /*lengthA*/length2 - scalar_length,/*lengthB*/length1 - (longset_ptr - set1),
			    intersection));

      diagterm = (unsigned short) diagterm1; /* result = set1 + diagterm1 */
      for (i = 0; i < ndiagonals_vector; i++) {
	debug15(printf("(Case 4) Intersection %d.  Writing %u + %d + %d\n",
		       intersection[i],region_term,intersection[i],diagterm));
	*out++ = region_term + intersection[i] + diagterm;
      }
    }

#ifdef DEBUG15
    ndiagonals_scalar = Intersect_uint2_scalar(diagonals_alloc_scalar,&longset_ptr,
					       set1,length1,diagterm1,set2,length2,diagterm2,region_term);
    if (ndiagonals != ndiagonals_scalar) {
      printf("%d vs %d\n",ndiagonals,ndiagonals_scalar);
      abort();
    } else {
      for (i = 0; i < ndiagonals; i++) {
	if (diagonals[i] != diagonals_alloc_scalar[i]) {
	  for (i = 0; i < ndiagonals; i++) {
	    printf("%u %u\n",diagonals[i],diagonals_alloc_scalar[i]);
	  }
	  abort();
	}
      }
    }
#endif

    return ndiagonals;
  }
}

#else

/* Non-SIMD version */ 
int
Intersect_uint2 (Univcoord_T *diagonals, unsigned short *localdb_alloc,
		 unsigned short *set1, const int length1, int diagterm1,
		 unsigned short *set2, const int length2, int diagterm2,
		 Univcoord_T region_term) {
  (void)(localdb_alloc);

  uint16_t *longset_ptr;

  if ((length1 == 0) || (length2 == 0)) {
    return 0;

  } else {
    return Intersect_uint2_scalar(diagonals,&longset_ptr,set1,length1,diagterm1,
				  set2,length2,diagterm2,region_term);
  }
}

#endif
