/* $Id: e8071d604ae3bd32d446e8009efa89669f9a58fe $ */
#ifndef GENOMEBITS_INCLUDED
#define GENOMEBITS_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_UNISTD_H, HAVE_SYS_TYPES_H, HAVE_CADDR_T */
#endif

#include <stdio.h>
#include "simd.h"

#if 0
#if defined(HAVE_SSE2)
#include <emmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip popcnt, which comes after SSE4.2 */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>		/* For _mm_popcnt_u32 */
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt and tzcnt, which come after SSE4.2 */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif
#endif


#include "types.h"
#include "access.h"
#include "bool.h"
#include "univcoord.h"
#include "popcount.h"


typedef struct Genomebits_T *Genomebits_T;
struct Genomebits_T {
  Access_T access;		/* For blocks: ALLOCATED_PRIVATE, ALLOCATED_SHARED, or MMAPPED */

  int high_blocks_shmid;
  key_t high_blocks_key;
  int low_blocks_shmid;
  key_t low_blocks_key;
  int flags_blocks_shmid;
  key_t flags_blocks_key;

  int high_fd;
  int low_fd;
  int flags_fd;
  size_t high_len;			/* filesize */
  size_t low_len;			/* filesize */
  size_t flags_len;			/* filesize */
  /* Univcoord_T genomelength;	-- from chromosome_iit */

  Genomecomp_T *high_blocks;
  Genomecomp_T *low_blocks;
  Genomecomp_T *flags_blocks;
};


/* Need a separate macro that is not in any other *.c file */
#ifdef DEBUGX
#define debugx(x) x
#else
#define debugx(x)
#endif


#define T Genomebits_T
extern void
Genomebits_print (T this, Univcoord_T startpos, Univcoord_T endpos);
extern void
Genomebits_free (T *old);
extern T
Genomebits_new (char *genomesubdir, char *fileroot, char *alt_suffix,
		Access_mode_T access, bool sharedp, bool revcompp);
extern char
Genomebits_get_char (T this, Univcoord_T pos);

#undef T


/* query_high_shifted, query_low_shifted, query_flags_shifted,
   ref_high_ptr, ref_low_ptr, ref_flags_ptr,
   (snp_high_ptr, snp_low_ptr, snp_flags_ptr,)
   plusp, genestrand, query_unk_mismatch_p, genome_unk_mismatch_p */
typedef UINT4 (*Diffproc_32_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				bool, int, bool, bool);
typedef UINT4 (*Diffproc_snp_32_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    bool, int, bool, bool);

typedef UINT8 (*Diffproc_64_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				bool, int, bool, bool);
typedef UINT8 (*Diffproc_snp_64_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				    bool, int, bool, bool);

#ifdef HAVE_SSE2
typedef __m128i (*Diffproc_128_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				   Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				   bool, int, bool, bool);
typedef __m128i (*Diffproc_snp_128_T) (Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				       Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				       Genomecomp_T *, Genomecomp_T *, Genomecomp_T *,
				       bool, int, bool, bool);
#endif


/*                    76543210 */
#define HIGH_BIT_32 0x80000000

/*                    7654321076543210 */
#define HIGH_BIT_64 0x8000000000000000

#define nonzero_p_32(diff) diff
#define nonzero_p_64(diff) diff

#define clear_start_32(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_start_64(diff,startdiscard) (diff & (~0ULL << (startdiscard)))
#define clear_end_32(diff,enddiscard) (diff & ~(~0U << (enddiscard)))
#define clear_end_64(diff,enddiscard) (diff & ~(~0ULL << (enddiscard)))

/* For trimming */
#define set_start_32(diff,startdiscard) (diff | ~(~0U << startdiscard))
#define set_start_64(diff,startdiscard) (diff | ~(~0ULL << startdiscard))
#define set_end_32(diff,enddiscard) (diff | (~0U << enddiscard))
#define set_end_64(diff,enddiscard) (diff | (~0ULL << enddiscard))

/* Same speed: clear_highbit(diff,relpos) (diff - (HIGH_BIT >> relpos)) */
/* Note: xor assumes that bit at relpos was on */
#define clear_highbit_32(diff,relpos) (diff ^ (HIGH_BIT_32 >> relpos))
#define clear_highbit_64(diff,relpos) (diff ^ (HIGH_BIT_64 >> relpos))

/* Slower: clear_lowbit(diff,relpos) diff -= (1 << relpos) */
#define clear_lowbit_32(diff,relpos) (diff & (diff - 1))
#define clear_lowbit_64(diff,relpos) (diff & (diff - 1))


#if defined(HAVE_LZCNT)
#define count_leading_zeroes_32(diff) _lzcnt_u32(diff)
#define count_leading_zeroes_64(diff) _lzcnt_u64(diff)
#elif defined(HAVE_BUILTIN_CLZ)
#define count_leading_zeroes_32(diff) __builtin_clz(diff)
#define count_leading_zeroes_64(diff) __builtin_clzll(diff)
#else
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#define count_leading_zeroes_64(diff) ((diff >> 48) ? clz_table[diff >> 48] : ((diff >> 32) ? 16 + clz_table[diff >> 32] : ((diff >> 16) ? 32 + clz_table[diff >> 16] : 48 + clz_table[diff])))
#endif

#if defined(HAVE_TZCNT)
#define count_trailing_zeroes_32(diff) _tzcnt_u32(diff)
#define count_trailing_zeroes_64(diff) _tzcnt_u64(diff)
#elif defined(HAVE_BUILTIN_CTZ)
#define count_trailing_zeroes_32(diff) __builtin_ctz(diff)
#define count_trailing_zeroes_64(diff) __builtin_ctzll(diff)
#else
/* lowbit = -diff & diff */
#define count_trailing_zeroes_32(diff) mod_37_bit_position[(-diff & diff) % 37]
#define count_trailing_zeroes_64(diff) popcount_ones_64(~diff & (diff - 1)) /* or 64 - popcount_ones_64(diff | ~diff) */
#endif


#ifdef HAVE_SSE2
static inline void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 (_mm_extract_epi16(x,1) << 16) | (_mm_extract_epi16(x,0) & 0x0000FFFF),
	 (_mm_extract_epi16(x,3) << 16) | (_mm_extract_epi16(x,2) & 0x0000FFFF),
	 (_mm_extract_epi16(x,5) << 16) | (_mm_extract_epi16(x,4) & 0x0000FFFF),
	 (_mm_extract_epi16(x,7) << 16) | (_mm_extract_epi16(x,6) & 0x0000FFFF));
  return;
}
#endif

#ifdef HAVE_SSE4_1
#define nonzero_p_128(diff) !_mm_testz_si128(diff,diff)
#else
#define nonzero_p_128(diff) _mm_movemask_epi8(diff);
#endif

#ifdef HAVE_SSE2
#define _BOUND_HIGH _mm_set_epi32(128,96,64,32)
#define _BOUND_LOW _mm_set_epi32(96,64,32,0)

static inline __m128i
clear_start_128 (__m128i _diff, int startdiscard) {
  __m128i _mask, _startdiscard;
#ifdef DEBUGX
  __m128i _result;
#endif

  debugx(printf("Clearing start to startdiscard %d\n",startdiscard));
  debugx(printf("Before: "));
  debugx(print_vector_hex(_diff));

#ifdef SSE2_SLLI_CONST_IMM8
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(startdiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), startdiscard % 32);
#endif
  _startdiscard = _mm_set1_epi32(startdiscard);
  _mask = _mm_or_si128(_mask, _mm_cmplt_epi32(_startdiscard, _BOUND_LOW));
  _mask = _mm_and_si128(_mask, _mm_cmplt_epi32(_startdiscard, _BOUND_HIGH));

#ifdef DEBUGX
  _result = _mm_and_si128(_mask, _diff);
  debugx(printf("After:  "));
  print_vector_hex(_result);
#endif

  return _mm_and_si128(_mask, _diff);
}
#endif


#ifdef HAVE_SSE2
static inline __m128i
clear_end_128 (__m128i _diff, int enddiscard) {
  __m128i _mask, _enddiscard;
#ifdef DEBUGX
  __m128i _result;
#endif

  debugx(printf("Clearing end from enddiscard %d\n",enddiscard));
  debugx(printf("Before: "));
  debugx(print_vector_hex(_diff));

#ifdef SSE2_SLLI_CONST_IMM8
  _mask = _mm_sll_epi32(_mm_set1_epi32(~0U), _mm_setr_epi32(enddiscard % 32,0,0,0));
#else
  _mask = _mm_slli_epi32(_mm_set1_epi32(~0U), enddiscard % 32);
#endif
  _enddiscard = _mm_set1_epi32(enddiscard);
  _mask = _mm_or_si128(_mm_cmplt_epi32(_enddiscard, _BOUND_LOW), _mask);
  _mask = _mm_and_si128(_mm_cmplt_epi32(_enddiscard, _BOUND_HIGH), _mask);

#ifdef DEBUGX
  _result = _mm_andnot_si128(_mask, _diff);
  debugx(printf("After:  "));
  print_vector_hex(_result);
#endif

  return _mm_andnot_si128(_mask, _diff);
}
#endif


#ifdef HAVE_SSE2
static inline __m128i
clear_highbit_128 (__m128i _diff, int leading_zeroes) {
  __m128i _subtract, _relpos;
  int relpos;

  relpos = 127 - leading_zeroes;
  debugx(printf("Clearing high bit at relpos %d\n",relpos));

#ifdef SSE2_SLLI_CONST_IMM8
  _subtract = _mm_sll_epi32(_mm_set1_epi32(1), _mm_setr_epi32(relpos % 32,0,0,0));
#else
  _subtract = _mm_slli_epi32(_mm_set1_epi32(1), relpos % 32);
#endif
  _relpos = _mm_set1_epi32(relpos);
  _subtract = _mm_and_si128(_mm_cmplt_epi32(_relpos, _BOUND_HIGH), _subtract);
  _subtract = _mm_andnot_si128(_mm_cmplt_epi32(_relpos, _BOUND_LOW), _subtract);

  debugx(printf("Subtract: "));
  debugx(print_vector_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm_sub_epi32(_diff, _subtract);
#else
  /* _mm_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm_xor_si128(_diff, _subtract);
#endif
}
#endif


#ifdef HAVE_SSE2
/* relpos is equal to trailing_zeroes */
static inline __m128i
clear_lowbit_128 (__m128i _diff, int relpos) {
  __m128i _subtract, _relpos;

  debugx(printf("Clearing low bit at relpos %d\n",relpos));

#ifdef SSE2_SLLI_CONST_IMM8
  _subtract = _mm_sll_epi32(_mm_set1_epi32(1), _mm_setr_epi32(relpos % 32,0,0,0));
#else
  _subtract = _mm_slli_epi32(_mm_set1_epi32(1), relpos % 32);
#endif
  _relpos = _mm_set1_epi32(relpos);
  _subtract = _mm_and_si128(_mm_cmplt_epi32(_relpos, _BOUND_HIGH), _subtract);
  _subtract = _mm_andnot_si128(_mm_cmplt_epi32(_relpos, _BOUND_LOW), _subtract);

  debugx(printf("Subtract: "));
  debugx(print_vector_hex(_subtract));
#if 0
  /* latency 1, throughput: 0.5 */
  return _mm_sub_epi32(_diff, _subtract);
#else
  /* _mm_xor_si128 also works if all other bits are 0.  latency 1, throughput: 0.33 */
  return _mm_xor_si128(_diff, _subtract);
#endif
}
#endif


#ifdef HAVE_SSE2
static inline int
count_leading_zeroes_128 (__m128i _diff) {
  debugx(printf("Entered count_leading_zeroes_128 with "));
  debugx(print_vector_hex(_diff));

#if defined(HAVE_LZCNT) && defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,1)) != 0) {
    return (int) _lzcnt_u64(x);
  } else {
    return 64 + (int) _lzcnt_u64(_mm_extract_epi64(_diff,0));
  }

#elif defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,1)) != 0) {
    return __builtin_clzll(x);
  } else {
    return 64 + __builtin_clzll(_mm_extract_epi64(_diff,0));
  }

#else
  UINT4 x;

  if ((x = (_mm_extract_epi16(_diff,7) << 16) | (_mm_extract_epi16(_diff,6) & 0x0000FFFF)) != 0) {
    debugx(printf("word 3 is non-empty, so returning %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[3])));
    return __builtin_clz(x);
  } else if ((x = (_mm_extract_epi16(_diff,5) << 16) | (_mm_extract_epi16(_diff,4) & 0x0000FFFF)) != 0) {
    debugx(printf("word 2 is non-empty, so returning 32 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[2])));
    return 32 + __builtin_clz(x);
  } else if ((x = (_mm_extract_epi16(_diff,3) << 16) | (_mm_extract_epi16(_diff,2) & 0x0000FFFF)) != 0) {
    debugx(printf("word 1 is non-empty, so returning 64 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[1])));
    return 64 + __builtin_clz(x);
  } else {
    x = (_mm_extract_epi16(_diff,1) << 16) | (_mm_extract_epi16(_diff,0) & 0x0000FFFF);
    debugx(printf("word 0 is non-empty, so returning 96 + %d vs %d\n",__builtin_clz(x),__builtin_clz(((UINT4 *) &_diff)[0])));
    return 96 + __builtin_clz(x);
  }
#endif
}
#endif


#ifdef HAVE_SSE2
static inline int
count_trailing_zeroes_128 (__m128i _diff) {
  debugx(printf("Entered count_trailing_zeroes_128 with "));
  debugx(print_vector_hex(_diff));

#if defined(HAVE_TZCNT) && defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,0)) != 0) {
    return (int) _tzcnt_u64(x);
  } else {
    return 64 + (int) _tzcnt_u64(_mm_extract_epi64(_diff,1));
  }

#elif defined(HAVE_SSE4_1) && defined(HAVE_MM_EXTRACT_EPI64)
  UINT8 x;

  if ((x = _mm_extract_epi64(_diff,0)) != 0) {
    return __builtin_ctzll(x);
  } else {
    return 64 + __builtin_ctzll(_mm_extract_epi64(_diff,1));
  }

#else
  UINT4 x;

  if ((x = (_mm_extract_epi16(_diff,1) << 16) | (_mm_extract_epi16(_diff,0) & 0x0000FFFF)) != 0) {
    debugx(printf("word 0 is non-empty, so returning %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[0])));
    return __builtin_ctz(x);
  } else if ((x = (_mm_extract_epi16(_diff,3) << 16) | (_mm_extract_epi16(_diff,2) & 0x0000FFFF)) != 0) {
    debugx(printf("word 1 is non-empty, so returning 32 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[1])));
    return 32 + __builtin_ctz(x);
  } else if ((x = (_mm_extract_epi16(_diff,5) << 16) | (_mm_extract_epi16(_diff,4) & 0x0000FFFF)) != 0) {
    debugx(printf("word 2 is non-empty, so returning 64 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[2])));
    return 64 + __builtin_ctz(x);
  } else {
    x = (_mm_extract_epi16(_diff,7) << 16) | (_mm_extract_epi16(_diff,6) & 0x0000FFFF);
    debugx(printf("word 3 is non-empty, so returning 96 + %d vs %d\n",__builtin_ctz(x),__builtin_ctz(((UINT4 *) &_diff)[3])));
    return 96 + __builtin_ctz(x);
  }
#endif
}
#endif


#define UNUSED(x) (void)(x)
#define cast64(x) (* (UINT8 *) x)

static inline UINT4
block_diff_standard_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  UINT4 diff;

#ifdef DEBUGX
  printf("block_diff_standard_32:\n");
  printf("  ref high: %08X\n",*ref_high_ptr);
  printf("query high: %08X\n",*query_high_shifted);

  printf("  ref low:  %08X\n",*ref_low_ptr);
  printf("query low:  %08X\n",*query_low_shifted);

  printf("  ref flag: %08X\n",*ref_flags_ptr);
  printf("query flag: %08X\n",*query_flags_shifted);
#endif

  diff = (*query_high_shifted ^ *ref_high_ptr) | (*query_low_shifted ^ *ref_low_ptr);

  debugx(printf("      diff: %08X\n",diff));


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= *query_flags_shifted;
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(*query_flags_shifted);
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= *ref_flags_ptr;
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(*ref_flags_ptr);
  }

  return diff;
}

static inline UINT4
block_diff_standard_wildcard_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
				 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
				 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
				 Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
				 bool plusp, int genestrand, bool query_unk_mismatch_p,
				 bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  UINT4 diff, snp_noteq_ref, wildcard;

#ifdef DEBUGX
  printf("block_diff_standard_wildcard_32:\n");
  printf("  ref high: %08X\n",*ref_high_ptr);
  printf("query high: %08X\n",*query_high_shifted);

  printf("  ref low:  %08X\n",*ref_low_ptr);
  printf("query low:  %08X\n",*query_low_shifted);
#endif

  /* Taken from block_diff_standard */
  diff = (*query_high_shifted ^ *ref_high_ptr) | (*query_low_shifted ^ *ref_low_ptr);

#ifdef DEBUGX
  printf("      diff: %08X\n",diff);
  printf("\n");
#endif


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= *query_flags_shifted;
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(*query_flags_shifted);
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= *ref_flags_ptr;
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(*ref_flags_ptr);
  }
  /* End of (query ^ ref) */


  /* A difference occurs if the query is different from ref and different from snp */
  /* Already computed difference from ref, so add difference from snp */
  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
#ifdef DEBUGX
  printf("  snp high: %08X\n",*snp_high_ptr);
  printf("  snp low:  %08X\n",*snp_low_ptr);
  printf("\n");

  printf("  ref flags:%08X\n",*ref_flags_ptr);
  printf("query flags:%08X\n",*query_flags_shifted);
  printf("  snp flags:%08X\n",*snp_flags_ptr);
#endif

  diff &= (*query_high_shifted ^ *snp_high_ptr) | (*query_low_shifted ^ *snp_low_ptr);


  /* A wildcard occurs if there is a SNP flag (and not a ref flag) and ref == snp */

  /* Difference between ref and snp */
  snp_noteq_ref = (*ref_high_ptr ^ *snp_high_ptr) | (*ref_low_ptr ^ *snp_low_ptr);
  
  /* Presence of a snp flag and not a ref flag */
  wildcard = ~(*ref_flags_ptr) & (*snp_flags_ptr);

  wildcard = ~snp_noteq_ref & wildcard;

#ifdef DEBUGX
  printf("  wildcard: %08X\n",wildcard);
  printf("    result: %08X\n",~wildcard & diff);
  printf("\n");
#endif

  return ~wildcard & diff;
}


/* Needs to match typedef of Diffproc_snp_128_T */
static inline UINT4
block_diff_standard_masked_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			       Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			       Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			       Genomecomp_T *mask_high_ptr, Genomecomp_T *mask_low_ptr,
			       Genomecomp_T *mask_flags_ptr,
			       bool plusp, int genestrand, bool query_unk_mismatch_p,
			       bool genome_unk_mismatch_p) {
  UNUSED(mask_high_ptr); UNUSED(mask_low_ptr);
  UNUSED(plusp); UNUSED(genestrand);
  UINT4 diff;

#ifdef DEBUGX
  printf("block_diff_standard_masked_32:\n");
  printf("  ref high: %08X\n",*ref_high_ptr);
  printf("query high: %08X\n",*query_high_shifted);

  printf("  ref low:  %08X\n",*ref_low_ptr);
  printf("query low:  %08X\n",*query_low_shifted);

#endif

  /* Taken from block_diff_standard */
  diff = (*query_high_shifted ^ *ref_high_ptr) | (*query_low_shifted ^ *ref_low_ptr);

#ifdef DEBUGX
  printf("      diff: %08X\n",diff);
  printf("\n");
#endif


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= *query_flags_shifted;
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(*query_flags_shifted);
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= *ref_flags_ptr;
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(*ref_flags_ptr);
  }
  /* End of (query ^ ref) */


#ifdef DEBUGX
  printf("  ref flags:%08X\n",*ref_flags_ptr);
  printf("query flags:%08X\n",*query_flags_shifted);
  printf(" mask flags:%08X\n",*mask_flags_ptr);
  printf("    result: %08X\n",~(*mask_flags_ptr) & diff);
  printf("\n");
#endif

  /* A difference occurs if the query is different from ref and is not masked */

  return ~(*mask_flags_ptr) & diff;
}





static inline UINT8
block_diff_standard_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  UINT8 diff;

#ifdef DEBUGX
  printf("block_diff_standard_64:\n");
  printf("  ref high: %016lX\n",cast64(ref_high_ptr));
  printf("query high: %016lX\n",cast64(query_high_shifted));

  printf("  ref low:  %016lX\n",cast64(ref_low_ptr));
  printf("query low:  %016lX\n",cast64(query_low_shifted));
#endif

  diff = (cast64(query_high_shifted) ^ cast64(ref_high_ptr)) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));

  debugx(printf("      diff: %016lX\n",diff));


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= cast64(query_flags_shifted);
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(cast64(query_flags_shifted));
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= cast64(ref_flags_ptr);
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(cast64(ref_flags_ptr));
  }

  return diff;
}

static inline UINT8
block_diff_standard_wildcard_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
				 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
				 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
				 Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
				 bool plusp, int genestrand, bool query_unk_mismatch_p,
				 bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  UINT8 diff, snp_noteq_ref, wildcard;

#ifdef DEBUGX
  printf("block_diff_standard_wildcard_64:\n");
  printf("  ref high: %016lX\n",cast64(ref_high_ptr));
  printf("query high: %016lX\n",cast64(query_high_shifted));

  printf("  ref low:  %016lX\n",cast64(ref_low_ptr));
  printf("query low:  %016lX\n",cast64(query_low_shifted));
#endif

  /* Taken from block_diff_standard */
  diff = (cast64(query_high_shifted) ^ cast64(ref_high_ptr)) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));

#ifdef DEBUGX
  printf("      diff: %016lX\n",diff);
  printf("\n");
#endif


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= cast64(query_flags_shifted);
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(cast64(query_flags_shifted));
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= cast64(ref_flags_ptr);
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(cast64(ref_flags_ptr));
  }
  /* End of (query ^ ref) */


  /* A difference occurs if the query is different from ref and different from snp */
  /* Already computed difference from ref, so add difference from snp */
  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
#ifdef DEBUGX
  printf("  snp high: %016lX\n",cast64(snp_high_ptr));
  printf("  snp low:  %016lX\n",cast64(snp_low_ptr));
  printf("\n");

  printf("  ref flags:%016lX\n",cast64(ref_flags_ptr));
  printf("query flags:%016lX\n",cast64(query_flags_shifted));
  printf("  snp flags:%016lX\n",cast64(snp_flags_ptr));
#endif

  diff &= (cast64(query_high_shifted) ^ cast64(snp_high_ptr)) | (cast64(query_low_shifted) ^ cast64(snp_low_ptr));


  /* A wildcard occurs if there is a SNP flag (and not a ref flag) and ref == snp */

  /* Difference between ref and snp */
  snp_noteq_ref = (cast64(ref_high_ptr) ^ cast64(snp_high_ptr)) | (cast64(ref_low_ptr) ^ cast64(snp_low_ptr));
  
  /* Presence of a snp flag and not a ref flag */
  wildcard = ~(cast64(ref_flags_ptr)) & (cast64(snp_flags_ptr));

  wildcard = ~snp_noteq_ref & wildcard;

#ifdef DEBUGX
  printf("  wildcard: %016lX\n",wildcard);
  printf("    result: %016lX\n",~wildcard & diff);
  printf("\n");
#endif

  return ~wildcard & diff;
}


/* Needs to match typedef of Diffproc_snp_128_T */
static inline UINT8
block_diff_standard_masked_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			       Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			       Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			       Genomecomp_T *mask_high_ptr, Genomecomp_T *mask_low_ptr,
			       Genomecomp_T *mask_flags_ptr,
			       bool plusp, int genestrand, bool query_unk_mismatch_p,
			       bool genome_unk_mismatch_p) {
  UNUSED(mask_high_ptr); UNUSED(mask_low_ptr);
  UNUSED(plusp); UNUSED(genestrand);
  UINT8 diff;

#ifdef DEBUGX
  printf("block_diff_standard_masked_64:\n");
  printf("  ref high: %016lX\n",cast64(ref_high_ptr));
  printf("query high: %016lX\n",cast64(query_high_shifted));

  printf("  ref low:  %016lX\n",cast64(ref_low_ptr));
  printf("query low:  %016lX\n",cast64(query_low_shifted));

#endif

  /* Taken from block_diff_standard */
  diff = (cast64(query_high_shifted) ^ cast64(ref_high_ptr)) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));

#ifdef DEBUGX
  printf("      diff: %016lX\n",diff);
  printf("\n");
#endif


  /* Query Ns */
  if (query_unk_mismatch_p) {
    /* Query: Considering N as a mismatch */
    diff |= cast64(query_flags_shifted);
  } else {
    /* Query: Considering N as a wildcard */
    diff &= ~(cast64(query_flags_shifted));
  }

  /* Genome Ns */
  if (genome_unk_mismatch_p) {
    /* Genome: Considering N as a mismatch */
    diff |= cast64(ref_flags_ptr);
  } else {
    /* Genome: Considering N as a wildcard */
    diff &= ~(cast64(ref_flags_ptr));
  }
  /* End of (query ^ ref) */


#ifdef DEBUGX
  printf("  ref flags:%016lX\n",cast64(ref_flags_ptr));
  printf("query flags:%016lX\n",cast64(query_flags_shifted));
  printf(" mask flags:%016lX\n",cast64(mask_flags_ptr));
  printf("    result: %016lX\n",~(cast64(mask_flags_ptr)) & diff);
  printf("\n");
#endif

  /* A difference occurs if the query is different from ref and is not masked */

  return ~(cast64(mask_flags_ptr)) & diff;
}


#ifdef HAVE_SSE2
static inline __m128i
block_diff_standard_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_p,
			 bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

#ifdef DEBUGX
  printf("block_diff_standard_128:\n");
  printf("  ref high: ");
  print_vector_hex(_ref_high);
  printf("query high: ");
  print_vector_hex(_query_high);

  printf("  ref low:  ");
  print_vector_hex(_ref_low);
  printf("query low:  ");
  print_vector_hex(_query_low);
#endif


  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

#ifdef DEBUGX
  printf("      diff: ");
  print_vector_hex(_diff);
  printf("\n");
#endif

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif


/* wildcard if ref == alt && ref_flag == 0 && alt_flag == 1 */
/* not wildcard if ref != alt || ref_flag == 1 || alt_flag == 0 */
/* diffs are (query ^ ref) & (query ^ alt) & ~wildcard */
/* snp_ptr here is alt_ptr */
#ifdef HAVE_SSE2
static inline __m128i
block_diff_standard_wildcard_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
				  Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
				  Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
				  Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr,
				  Genomecomp_T *snp_flags_ptr,
				  bool plusp, int genestrand, bool query_unk_mismatch_p,
				  bool genome_unk_mismatch_p) {
  UNUSED(plusp); UNUSED(genestrand);
  __m128i _diff, _snp_noteq_ref, _wildcard, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _snp_high, _snp_low, _snp_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

#ifdef DEBUGX
  printf("block_diff_standard_wildcard_128:\n");
  printf("  ref high: ");
  print_vector_hex(_ref_high);
  printf("query high: ");
  print_vector_hex(_query_high);

  printf("  ref low:  ");
  print_vector_hex(_ref_low);
  printf("query low:  ");
  print_vector_hex(_query_low);
#endif

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

#ifdef DEBUGX
  printf("      diff: ");
  print_vector_hex(_diff);
#endif


  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }
  /* End of (query ^ ref) */


  /* A difference occurs if the query is different from ref and different from snp */
  /* Already computed difference from ref, so add difference from snp */
  /* Add (query ^ snp).  Don't need to recompute query flags or use SNP flags. */
  _snp_high = _mm_load_si128((__m128i *) snp_high_ptr);
  _snp_low = _mm_load_si128((__m128i *) snp_low_ptr);
  _snp_flags = _mm_load_si128((__m128i *) snp_flags_ptr);

#ifdef DEBUGX
  printf("  snp high: ");
  print_vector_hex(_snp_high);
  printf("  snp low:  ");
  print_vector_hex(_snp_low);
  printf("\n");

  printf("  ref flags:");
  print_vector_hex(_ref_flags);
  printf("query flags:");
  print_vector_hex(_query_flags);
  printf("  snp flags:");
  print_vector_hex(_snp_flags);
#endif

  _diff = _mm_and_si128(_diff, _mm_or_si128(_mm_xor_si128(_query_high, _snp_high), _mm_xor_si128(_query_low, _snp_low)));


  /* A wildcard occurs if there is a SNP flag (and not a ref flag) and ref == snp */
  _snp_noteq_ref = _mm_or_si128(_mm_xor_si128(_ref_high, _snp_high), _mm_xor_si128(_ref_low, _snp_low));

  _wildcard = _mm_andnot_si128(_ref_flags, _snp_flags); /* snp flag and not a ref flag */
  _wildcard = _mm_andnot_si128(_snp_noteq_ref, _wildcard);


  /* The final difference occurs if there is a diff and it is not a wildcard */
  _diff = _mm_andnot_si128(_wildcard,_diff);


#ifdef DEBUGX
  printf("  wildcard: ");
  print_vector_hex(_wildcard);
  printf("    result: ");
  print_vector_hex(_diff);
  printf("\n");
#endif

  return _diff;
}
#endif


#ifdef HAVE_SSE2
/* Needs to match typedef of Diffproc_snp_128_T */
static inline __m128i
block_diff_standard_masked_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
				Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
				Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
				Genomecomp_T *mask_high_ptr, Genomecomp_T *mask_low_ptr,
				Genomecomp_T *mask_flags_ptr,
				bool plusp, int genestrand, bool query_unk_mismatch_p,
				bool genome_unk_mismatch_p) {
  UNUSED(mask_high_ptr); UNUSED(mask_low_ptr);
  UNUSED(plusp); UNUSED(genestrand);
  __m128i _diff, _query_high, _query_low, _query_flags,
    _ref_high, _ref_low, _ref_flags, _mask_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

#ifdef DEBUGX
  printf("block_diff_standard_masked_128:\n");
  printf("  ref high: ");
  print_vector_hex(_ref_high);
  printf("query high: ");
  print_vector_hex(_query_high);

  printf("  ref low:  ");
  print_vector_hex(_ref_low);
  printf("query low:  ");
  print_vector_hex(_query_low);
#endif

  _diff = _mm_or_si128(_mm_xor_si128(_query_high, _ref_high), _mm_xor_si128(_query_low, _ref_low));

#ifdef DEBUGX
  printf("      diff: ");
  print_vector_hex(_diff);
#endif

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }
  /* End of (query ^ ref) */


  /* A difference occurs if the query is different from ref and is not masked */
  _mask_flags = _mm_load_si128((__m128i *) mask_flags_ptr);

#ifdef DEBUGX
  printf("  ref flags:");
  print_vector_hex(_ref_flags);
  printf("query flags:");
  print_vector_hex(_query_flags);
  printf(" mask flags:");
  print_vector_hex(_mask_flags);
#endif

  _diff = _mm_andnot_si128(_mask_flags,_diff);

#ifdef DEBUGX
  printf("    result: ");
  print_vector_hex(_diff);
  printf("\n");
#endif

  return _diff;
}
#endif



/************************************************************************
 *   CMET
 ************************************************************************/

static inline UINT4
block_diff_metct_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT4 diff;

  /* Mark genome-T to query-C mismatches */
  diff = (~(*query_high_shifted) & *query_low_shifted) & (*ref_high_ptr & *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));

  /* Compare reduced C->T nts  */
  diff |= ((*query_high_shifted | *query_low_shifted) ^ (*ref_high_ptr | *ref_low_ptr)) | (*query_low_shifted ^ *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));


  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= *query_flags_shifted;
  } else {
    diff &= ~(*query_flags_shifted);
  }

  if (genome_unk_mismatch_p) {
    diff |= *ref_flags_ptr;
  } else {
    diff &= ~(*ref_flags_ptr);
  }
  debugx(printf(" => diff %08X\n",diff));

  return diff;
}

static inline UINT4
block_diff_metga_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT4 diff;

  /* Mark genome-A to query-G mismatches */
  diff = (*query_high_shifted & ~(*query_low_shifted)) & ~(*ref_high_ptr | *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));

  /* Compare reduced G->A nts  */
  diff |= ((*query_high_shifted & *query_low_shifted) ^ (*ref_high_ptr & *ref_low_ptr)) | (*query_low_shifted ^ *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));


  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= *query_flags_shifted;
  } else {
    diff &= ~(*query_flags_shifted);
  }

  if (genome_unk_mismatch_p) {
    diff |= (*ref_flags_ptr);
  } else {
    diff &= ~(*ref_flags_ptr);
  }
  debugx(printf(" => diff %08X\n",diff));

  return diff;
}




static inline UINT8
block_diff_metct_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT8 diff;

  /* Mark genome-T to query-C mismatches */
  diff = (~(cast64(query_high_shifted)) & cast64(query_low_shifted)) & (cast64(ref_high_ptr) & cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));

  /* Compare reduced C->T nts  */
  diff |= ((cast64(query_high_shifted) | cast64(query_low_shifted)) ^ (cast64(ref_high_ptr) | cast64(ref_low_ptr))) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));


  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= cast64(query_flags_shifted);
  } else {
    diff &= ~(cast64(query_flags_shifted));
  }

  if (genome_unk_mismatch_p) {
    diff |= cast64(ref_flags_ptr);
  } else {
    diff &= ~(cast64(ref_flags_ptr));
  }
  debugx(printf(" => diff %016lX\n",diff));

  return diff;
}

static inline UINT8
block_diff_metga_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT8 diff;

  /* Mark genome-A to query-G mismatches */
  diff = (cast64(query_high_shifted) & ~(cast64(query_low_shifted))) & ~(cast64(ref_high_ptr) | cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));

  /* Compare reduced G->A nts  */
  diff |= ((cast64(query_high_shifted) & cast64(query_low_shifted)) ^ (cast64(ref_high_ptr) & cast64(ref_low_ptr))) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));


  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= cast64(query_flags_shifted);
  } else {
    diff &= ~(cast64(query_flags_shifted));
  }

  if (genome_unk_mismatch_p) {
    diff |= (cast64(ref_flags_ptr));
  } else {
    diff &= ~(cast64(ref_flags_ptr));
  }
  debugx(printf(" => diff %016lX\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert C to T: high/low (A) 0 0 => new high 0; (C) 0 1 => 1; (G) 1 0 => 1; (T) 1 0 => 1 */
/* new high = high | low */
static inline __m128i
block_diff_metct_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		      Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		      Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

  /* Mark genome-T to query-C mismatches */
  _diff = _mm_and_si128(_mm_andnot_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low));

  /* Compare reduced C->T nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_or_si128(_query_high, _query_low), _mm_or_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif	/* HAVE_SSE2 */


#ifdef HAVE_SSE2
/* Convert G to A: high/low (A) 0 0 => new high 0; (C) 0 1 => 0; (G) 1 0 => 0; (T) 1 0 => 1 */
/* new high = high & low */
static inline __m128i
block_diff_metga_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		      Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		      Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

  /* Mark genome-A to query-G mismatches */
  _diff = _mm_andnot_si128(_query_low, _query_high);
  _diff = _mm_andnot_si128(_ref_high, _diff);
  _diff = _mm_andnot_si128(_ref_low, _diff);

  /* Compare reduced G->A nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_and_si128(_query_high, _query_low), _mm_and_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static inline UINT4
block_diff_cmet_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT4
block_diff_cmet_snp_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

static inline UINT8
block_diff_cmet_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT8
block_diff_cmet_snp_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#ifdef HAVE_SSE2
static inline __m128i
block_diff_cmet_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool plusp, int genestrand, bool query_unk_mismatch_p,
		     bool genome_unk_mismatch_p) {

  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#ifdef HAVE_SSE2
/* Ignores snp ptrs */
static inline __m128i
block_diff_cmet_snp_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			 Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_p,
			 bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_metga_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metct_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_metct_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_metga_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif


/************************************************************************
 *   ATOI
 ************************************************************************/

static inline UINT4
block_diff_a2iag_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT4 diff;

  /* Mark genome-G to query-A mismatches */
  diff = ~(*query_high_shifted | *query_low_shifted) & (*ref_high_ptr & ~(*ref_low_ptr));
  debugx(printf(" => diff %08X\n",diff));

  /* Compare reduced A->G nts  */
  diff |= ((*query_high_shifted | ~(*query_low_shifted)) ^ (*ref_high_ptr | ~(*ref_low_ptr))) | (*query_low_shifted ^ *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));

  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= *query_flags_shifted;
  } else {
    diff &= ~(*query_flags_shifted);
  }

  if (genome_unk_mismatch_p) {
    diff |= (*ref_flags_ptr);
  } else {
    diff &= ~(*ref_flags_ptr);
  }
  debugx(printf(" => diff %08X\n",diff));

  return diff;
}

static inline UINT4
block_diff_a2itc_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT4 diff;

  /* Mark genome-C to query-T mismatches */
  diff = (*query_high_shifted & *query_low_shifted) & (~(*ref_high_ptr) & *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));

  /* Compare reduced T->C nts  */
  diff |= ((*query_high_shifted & ~(*query_low_shifted)) ^ (*ref_high_ptr & ~(*ref_low_ptr))) | (*query_low_shifted ^ *ref_low_ptr);
  debugx(printf(" => diff %08X\n",diff));

  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= *query_flags_shifted;
  } else {
    diff &= ~(*query_flags_shifted);
  }

  if (genome_unk_mismatch_p) {
    diff |= *ref_flags_ptr;
  } else {
    diff &= ~(*ref_flags_ptr);
  }
  debugx(printf(" => diff %08X\n",diff));

  return diff;
}



static inline UINT8
block_diff_a2iag_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT8 diff;

  /* Mark genome-G to query-A mismatches */
  diff = ~(cast64(query_high_shifted) | cast64(query_low_shifted)) & (cast64(ref_high_ptr) & ~(cast64(ref_low_ptr)));
  debugx(printf(" => diff %016lX\n",diff));

  /* Compare reduced A->G nts  */
  diff |= ((cast64(query_high_shifted) | ~(cast64(query_low_shifted))) ^ (cast64(ref_high_ptr) | ~(cast64(ref_low_ptr)))) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));

  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= cast64(query_flags_shifted);
  } else {
    diff &= ~(cast64(query_flags_shifted));
  }

  if (genome_unk_mismatch_p) {
    diff |= (cast64(ref_flags_ptr));
  } else {
    diff &= ~(cast64(ref_flags_ptr));
  }
  debugx(printf(" => diff %016lX\n",diff));

  return diff;
}

static inline UINT8
block_diff_a2itc_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  UINT8 diff;

  /* Mark genome-C to query-T mismatches */
  diff = (cast64(query_high_shifted) & cast64(query_low_shifted)) & (~(cast64(ref_high_ptr)) & cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));

  /* Compare reduced T->C nts  */
  diff |= ((cast64(query_high_shifted) & ~(cast64(query_low_shifted))) ^ (cast64(ref_high_ptr) & ~(cast64(ref_low_ptr)))) | (cast64(query_low_shifted) ^ cast64(ref_low_ptr));
  debugx(printf(" => diff %016lX\n",diff));

  /* Flags: Considering N as a mismatch */
  if (query_unk_mismatch_p) {
    diff |= cast64(query_flags_shifted);
  } else {
    diff &= ~(cast64(query_flags_shifted));
  }

  if (genome_unk_mismatch_p) {
    diff |= cast64(ref_flags_ptr);
  } else {
    diff &= ~(cast64(ref_flags_ptr));
  }
  debugx(printf(" => diff %016lX\n",diff));

  return diff;
}


#ifdef HAVE_SSE2
/* Convert A->G: new high = high | ~low */
static inline __m128i
block_diff_a2iag_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		      Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		      Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

  /* Mark genome-G to query-A mismatches */
  _diff = _mm_andnot_si128(_mm_or_si128(_query_high, _query_low), _mm_andnot_si128(_ref_low, _ref_high));

  /* Compare reduced A->G nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif


#ifdef HAVE_SSE2
/* Convert T->C: new high = high & ~low */
static inline __m128i
block_diff_a2itc_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		      Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		      Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  __m128i _diff, _query_high, _query_low, _query_flags, _ref_high, _ref_low, _ref_flags;

  _query_high = _mm_load_si128((__m128i *) query_high_shifted);
  _query_low = _mm_load_si128((__m128i *) query_low_shifted);
  _ref_high = _mm_load_si128((__m128i *) ref_high_ptr);
  _ref_low = _mm_load_si128((__m128i *) ref_low_ptr);

  /* Mark genome-C to query-T mismatches */
  _diff = _mm_and_si128(_mm_and_si128(_query_high, _query_low), _mm_andnot_si128(_ref_high, _ref_low));

  /* Compare reduced T->C nts  */
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_mm_andnot_si128(_query_low, _query_high), _mm_andnot_si128(_ref_low, _ref_high)));
  _diff = _mm_or_si128(_diff, _mm_xor_si128(_query_low, _ref_low));

  _query_flags = _mm_load_si128((__m128i *) query_flags_shifted);
  if (query_unk_mismatch_p) {
    _diff = _mm_or_si128(_query_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_query_flags, _diff);
  }

  _ref_flags = _mm_load_si128((__m128i *) ref_flags_ptr);
  if (genome_unk_mismatch_p) {
    _diff = _mm_or_si128(_ref_flags, _diff);
  } else {
    _diff = _mm_andnot_si128(_ref_flags, _diff);
  }

  return _diff;
}
#endif


static inline UINT4
block_diff_atoi_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT4
block_diff_atoi_snp_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

static inline UINT8
block_diff_atoi_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT8
block_diff_atoi_snp_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#ifdef HAVE_SSE2
static inline __m128i
block_diff_atoi_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool plusp, int genestrand, bool query_unk_mismatch_p,
		     bool genome_unk_mismatch_p) {

  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#ifdef HAVE_SSE2
/* Ignores snp ptrs */
static inline __m128i
block_diff_atoi_snp_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			 Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_p,
			 bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif


/************************************************************************
 *  TTOC
 ************************************************************************/

static inline UINT4
block_diff_ttoc_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT4
block_diff_ttoc_snp_32 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_32(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

static inline UINT8
block_diff_ttoc_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		    Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		    Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		    bool plusp, int genestrand, bool query_unk_mismatch_p,
		    bool genome_unk_mismatch_p) {
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}

#ifdef GSNAP
/* Ignores snp ptrs */
static inline UINT8
block_diff_ttoc_snp_64 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			bool plusp, int genestrand, bool query_unk_mismatch_p,
			bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_64(query_high_shifted,query_low_shifted,query_flags_shifted,
				 ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				 query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif


#ifdef HAVE_SSE2
static inline __m128i
block_diff_ttoc_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
		     Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
		     Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
		     bool plusp, int genestrand, bool query_unk_mismatch_p,
		     bool genome_unk_mismatch_p) {

  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#ifdef HAVE_SSE2
/* Ignores snp ptrs */
static inline __m128i
block_diff_ttoc_snp_128 (Genomecomp_T *query_high_shifted, Genomecomp_T *query_low_shifted,
			 Genomecomp_T *query_flags_shifted, Genomecomp_T *ref_high_ptr,
			 Genomecomp_T *ref_low_ptr, Genomecomp_T *ref_flags_ptr,
			 Genomecomp_T *snp_high_ptr, Genomecomp_T *snp_low_ptr, Genomecomp_T *snp_flags_ptr,
			 bool plusp, int genestrand, bool query_unk_mismatch_p,
			 bool genome_unk_mismatch_p) {
  UNUSED(snp_high_ptr); UNUSED(snp_low_ptr); UNUSED(snp_flags_ptr);
  if (genestrand == +2) {
    if (plusp) {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  } else {
    if (plusp) {
      return block_diff_a2itc_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    } else {
      return block_diff_a2iag_128(query_high_shifted,query_low_shifted,query_flags_shifted,
				  ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  query_unk_mismatch_p,genome_unk_mismatch_p);
    }
  }
}
#endif

#endif	/* GENOMEBITS_INCLUDED */

