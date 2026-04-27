/* $Id: 18e31fe30be4ff60bade2ba02ea52ccee85226c3 $ */
#ifndef GENOMEBITS_INDEL_INCLUDED
#define GENOMEBITS_INDEL_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"

#define T Genomebits_T

extern int
Genomebits_indel_solve_high (int *best_trimpos, int *nmismatches_to_trimpos,
			     Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
			     Compress_T query_compress, int *mismatch_positions_alloc,
			     Genomebits_T omebits, Genomebits_T omebits_alt,
			     bool plusp, int genestrand);

extern int
Genomebits_indel_solve_low (int *best_trimpos, int *nmismatches_to_trimpos,
			    Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
			    Compress_T query_compress, int *mismatch_positions_alloc,
			    Genomebits_T omebits, Genomebits_T omebits_alt,
			    bool plusp, int genestrand);

extern void
Genomebits_indel_setup (int max_insertionlen_in, int max_deletionlen_in,
			bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			Mode_T mode, bool maskedp);

#undef T
#endif
