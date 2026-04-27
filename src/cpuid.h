/* $Id: cpuid.h 226315 2023-02-28 18:22:20Z twu $ */
#ifndef CPUID_INCLUDED
#define CPUID_INCLUDED
#include "bool.h"

extern void
CPUID_support (bool *arm_support_p,
	       bool *sse2_support_p, bool *ssse3_support_p, bool *sse4_1_support_p, bool *sse4_2_support_p,
	       bool *avx2_support_p, bool *avx512_support_p, bool *avx512bw_support_p);

#endif

