/* $Id: b02a0dcd010a3a23b75f24d12c356a2d418a463e $ */
#ifndef COMPRESS_INCLUDED
#define COMPRESS_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_SSE2 */
#endif

#include <stdio.h>
#include "bool.h"
#include "types.h"
#include "genomicpos.h"


#define T Compress_T
typedef struct T *T;

extern void
Compress_free (T *old);
extern void
Compress_print (T this, int nshift, int pos5, int pos3);
extern void
Compress_print_queryseq (T this, int pos5, int pos3);
extern char *
Compress_queryseq (T this, int querylength);
extern T
Compress_new_fwd (char *gbuffer, Chrpos_T length);
extern T
Compress_new_rev (char *gbuffer, Chrpos_T length);
extern void
Compress_shift (Genomecomp_T **query_high_shifted, Genomecomp_T **query_low_shifted,
		Genomecomp_T **query_flags_shifted, T this, int nshift, int initpos);
extern bool
Compress_non_acgt (T this);
extern bool
Compress_fwdp (T this);

#undef T
#endif

