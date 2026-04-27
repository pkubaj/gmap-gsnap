#ifndef SELECT64_COMMON_INCLUDED
#define SELECT64_COMMON_INCLUDED

#include "types.h"

/* Required by select64 */
extern const uint8_t kSelectInByte[2048];


/** Count the number of 1-bits in a word.
 * @param word binary word.
 *
 */
static inline int
nu (uint64_t word) {
  return __builtin_popcountll(word);
}

/** Find the index of the most significant 1-bit in a word.
 * @param word binary word.
 *
 * The Knuth's lambda function is the dual of the rho function.
 *
 * Returns -1 on input zero.
 *
 */
static inline int
lambda_safe (uint64_t word) {
  return word == 0 ? -1 : 63 ^ __builtin_clzll(word);
}


/** Returns the index of the k-th 1-bit in the 64-bit word x.
 * @param x 64-bit word.
 * @param k 0-based rank (`k = 0` returns the position of the first 1-bit).
 *
 * Uses the broadword selection algorithm by Vigna [1], improved by Gog and Petri [2] and Vigna [3].
 * Facebook's Folly implementation [4].
 *
 * [1] Sebastiano Vigna. Broadword Implementation of Rank/Select Queries. WEA, 2008
 *
 * [2] Simon Gog, Matthias Petri. Optimized succinct data structures for massive data. Softw. Pract.
 * Exper., 2014
 *
 * [3] Sebastiano Vigna. MG4J 5.2.1. http://mg4j.di.unimi.it/
 *
 * [4] Facebook Folly library: https://github.com/facebook/folly
 *
 */
static inline uint64_t
select64 (uint64_t x, uint64_t k) {
#ifndef __haswell__
	const uint64_t kOnesStep4 = 0x1111111111111111ULL;
	const uint64_t kOnesStep8 = 0x0101010101010101ULL;
	const uint64_t kLAMBDAsStep8 = 0x80ULL * kOnesStep8;

	uint64_t s = x;
	s = s - ((s & 0xA * kOnesStep4) >> 1);
	s = (s & 0x3 * kOnesStep4) + ((s >> 2) & 0x3 * kOnesStep4);
	s = (s + (s >> 4)) & 0xF * kOnesStep8;
	uint64_t byteSums = s * kOnesStep8;

	uint64_t kStep8 = k * kOnesStep8;
	uint64_t geqKStep8 = (((kStep8 | kLAMBDAsStep8) - byteSums) & kLAMBDAsStep8);
	uint64_t place = nu(geqKStep8) * 8;
	uint64_t byteRank = k - (((byteSums << 8) >> place) & (uint64_t) (0xFF));
	return place + kSelectInByte[((x >> place) & 0xFF) | (byteRank << 8)];
#elif defined(__GNUC__) || defined(__clang__)
	// GCC and Clang won't inline the intrinsics.
	uint64_t result = uint64_t(1) << k;

	asm("pdep %1, %0, %0\n\t"
		"tzcnt %0, %0"
		: "+r"(result)
		: "r"(x));

	return result;
#else
	return _tzcnt_u64(_pdep_u64(1ULL << k, x));
#endif
}

#endif

