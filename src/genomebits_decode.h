/* $Id: eb7cd14cd57a70a8cb07f8be4550daaa71d59d0e $ */
#ifndef GENOMEBITS_DECODE_INCLUDED
#define GENOMEBITS_DECODE_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "types.h"

#ifdef HAVE_SSE2
extern int
Genomebits_decode_trailing_64 (int *out, int nmismatches, uint64_t word, int offset,
			       int max_mismatches);
extern int
Genomebits_decode_trailing_32 (int *out, int nmismatches, uint32_t word, int offset,
			       int max_mismatches);
extern int
Genomebits_decode_leading_64 (int *out, int nmismatches, uint64_t word, int offset,
			      int max_mismatches);
extern int
Genomebits_decode_leading_32 (int *out, int nmismatches, uint32_t word, int offset,
			      int max_mismatches);
#endif

#endif
