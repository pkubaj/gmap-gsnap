/* $Id: 37510e4a278f8dc1b57e23d436ef843c67e899e6 $ */
#ifndef GENOMEBITS_MISMATCHES_INCLUDED
#define GENOMEBITS_MISMATCHES_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"

#ifdef HAVE_SSE2
/* Could write 8 mismatches after querylength */
#define MISMATCH_EXTRA 8
#else
/* Could write 1 mismatch after querylength */
#define MISMATCH_EXTRA 1
#endif


#define T Genomebits_T

extern int
Genomebits_mismatches_fromleft (int *mismatch_positions, int max_mismatches, T ref, T alt,
				Compress_T query_compress,
				Univcoord_T univdiagonal, int querylength,
				int pos5, int pos3, bool plusp, int genestrand);

extern int
Genomebits_mismatches_fromleft_for_trim (int *mismatch_positions, int max_mismatches, T ref, T alt,
					 Compress_T query_compress,
					 Univcoord_T univdiagonal, int querylength,
					 int pos5, int pos3, bool plusp, int genestrand);

extern int
Genomebits_mismatches_fromright (int *mismatch_positions, int max_mismatches, T ref, T alt,
				 Compress_T query_compress,
				 Univcoord_T univdiagonal, int querylength,
				 int pos5, int pos3, bool plusp, int genestrand);

extern int
Genomebits_mismatches_fromright_for_trim (int *mismatch_positions, int max_mismatches, T ref, T alt,
					  Compress_T query_compress,
					  Univcoord_T univdiagonal, int querylength,
					  int pos5, int pos3, bool plusp, int genestrand);
extern void
Genomebits_mismatches_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			     Mode_T mode, bool maskedp);

#undef T
#endif



