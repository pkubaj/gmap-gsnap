/* $Id: 838649b51cd9dfcf6e77f61b92b05da7140f0d5f $ */
#ifndef GENOMEBITS_COUNT_INCLUDED
#define GENOMEBITS_COUNT_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"

#define T Genomebits_T

extern int
Genomebits_count_mismatches_substring (int *ref_mismatches, T ref, T alt,
				       Compress_T query_compress,
				       Univcoord_T univdiagonal, int querylength,
				       int pos5, int pos3, bool plusp, int genestrand);
extern int
Genomebits_mark_mismatches (int *nmatches_exonic, char *genomic,
			    Compress_T query_compress,
			    Univcoord_T univdiagonal, int querylength,
			    int pos5, int pos3, bool segment_plusp, bool query_plusp,
			    int genestrand);
extern void
Genomebits_count_setup (T ref, T alt,
			bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			Mode_T mode, bool md_report_snps_p_in, bool maskedp_in);

#undef T
#endif
