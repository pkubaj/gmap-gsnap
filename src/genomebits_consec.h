/* $Id: a38cb3be7e1a22a9a7475450f61175e998b6802f $ */
#ifndef GENOMEBITS_CONSEC_INCLUDED
#define GENOMEBITS_CONSEC_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"

#define T Genomebits_T

extern int
Genomebits_consecutive_matches_wmm (char *mismatch_char, T ref, Compress_T query_compress,
				    Univcoord_T univdiagonal, int querylength,
				    int pos5, int pos3, bool plusp, int genestrand);
extern int
Genomebits_consecutive_matches_rightward (T ref, Compress_T query_compress,
					  Univcoord_T univdiagonal, int left,
					  int pos5, int pos3, bool plusp, int genestrand);
extern int
Genomebits_consecutive_matches_leftward (T ref, Compress_T query_compress,
					 Univcoord_T univdiagonal, int left,
					 int pos5, int pos3, bool plusp, int genestrand);

extern void
Genomebits_consec_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			 Mode_T mode);

#undef T
#endif

