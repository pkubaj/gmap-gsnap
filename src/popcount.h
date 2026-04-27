/* $Id: popcount.h 226315 2023-02-28 18:22:20Z twu $ */
#ifndef POPCOUNT_INCLUDED
#define POPCOUNT_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_BUILTIN_CTZ, HAVE_BUILTIN_POPCOUNT, HAVE_BUILTIN_CLZ */
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_TZCNT) && !defined(HAVE_BUILTIN_CTZ))
extern const int mod_37_bit_position[];
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_POPCNT) && !defined(HAVE_MM_POPCNT) && !defined(HAVE_BUILTIN_POPCOUNT))
extern const int count_bits[];
#endif

#if !defined(HAVE_SSE4_2) || (!defined(HAVE_LZCNT) && !defined(HAVE_BUILTIN_CLZ))
extern const int clz_table[];
#endif

#if !defined(HAVE_SSE4_2)
/* Skip popcnt, which comes after SSE4.2 */
#define popcount_ones_32(diff) (count_bits[diff >> 16] + count_bits[diff & 0xFFFF])
#define popcount_ones_64(diff) (count_bits[diff >> 48] + count_bits[(diff & 0xFFFF00000000) >> 32] + count_bits[(diff & 0xFFFF0000) >> 16] + count_bits[diff & 0xFFFF])
#elif defined(HAVE_POPCNT)
#define popcount_ones_32(diff) (_popcnt32(diff))
#define popcount_ones_64(diff) (_popcnt64(diff))
#define popcount_ones_128(_diff) (_popcnt64(_mm_extract_epi64(_diff,0)) + _popcnt64(_mm_extract_epi64(_diff,1)))
#elif defined(HAVE_MM_POPCNT)
#define popcount_ones_32(diff) (_mm_popcnt_u32(diff))
#define popcount_ones_64(diff) (_mm_popcnt_u64(diff))
#define popcount_ones_128(_diff) (_mm_popcnt_u64(_mm_extract_epi64(_diff,0)) + _mm_popcnt_u64(_mm_extract_epi64(_diff,1)))
#elif defined(HAVE_BUILTIN_POPCOUNT)
#define popcount_ones_32(diff) (__builtin_popcount(diff))
#define popcount_ones_64(diff) (__builtin_popcountll(diff))
#define popcount_ones_128(_diff) (__builtin_popcountll(_mm_extract_epi64(_diff,0)) + __builtin_popcountll(_mm_extract_epi64(_diff,1)))
#else
#define popcount_ones_32(diff) (count_bits[diff >> 16] + count_bits[diff & 0xFFFF])
#define popcount_ones_64(diff) (count_bits[diff >> 48] + count_bits[(diff & 0xFFFF00000000) >> 32] + count_bits[(diff & 0xFFFF0000) >> 16] + count_bits[diff & 0xFFFF])
#endif


#endif

