/* $Id: a6a2d0187a740658544b51dfa26a31dd0475f33f $ */
#ifndef SIMD_INCLUDED
#define SIMD_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* configure script determines value of USE_INTEL_INTRINSICS */
/* #define USE_INTEL_INTRINSICS 1 */

#ifdef USE_INTEL_INTRINSICS
/* Intel intrinsic library */

#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#ifdef HAVE_SSE4_2
#include <nmmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#ifdef HAVE_AVX512
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip mm_popcnt */
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt and tzcnt, which come after SSE4.2 */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif

#else
/* SIMD everywhere package */

#define SIMDE_ENABLE_NATIVE_ALIASES 1

#ifdef HAVE_ARM
#include "simde/arm/neon.h"
#include "simde/x86/sse4.2.h"
#else

#ifdef HAVE_SSE2
#include "simde/x86/sse2.h"
#endif
#ifdef HAVE_SSSE3
#include "simde/x86/sse3.h"
#endif
#ifdef HAVE_SSE4_1
#include "simde/x86/sse4.1.h"
#endif
#ifdef HAVE_SSE4_2
#include "simde/x86/sse4.2.h"
#endif
#ifdef HAVE_AVX2
#include "simde/x86/avx2.h"
#endif
#ifdef HAVE_AVX512
#include "simde/x86/avx512.h"
#endif

#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include "simde/x86/sse4.2.h"
#endif

#if !defined(HAVE_SSE4_2)
/* Skip mm_popcnt */
#elif defined(HAVE_MM_POPCNT)
#include "simde/x86/sse4.2.h"
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt and tzcnt, which come after SSE4.2 */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include "simde/x86/sse4.2.h"
#endif

#endif	/* HAVE_ARM */

#endif	/* USE_INTEL_INTRINSICS */

#endif	/* SIMD_INCLUDED */

