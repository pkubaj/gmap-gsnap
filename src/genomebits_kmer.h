/* $Id: 133a6af3f2d188bcb78ae2bdb672eed6a04d0dee $ */
#ifndef GENOMEBITS_KMER_INCLUDED
#define GENOMEBITS_KMER_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"

#define T Genomebits_T

extern int
Genomebits_first_kmer_left (int *nmismatches5, T ref, Compress_T query_compress,
			    Univcoord_T univdiagonal, int querylength,
			    int pos5, int pos3, bool plusp, int genestrand,
			    bool query_unk_mismatch_p, int kmer);
extern int
Genomebits_first_kmer_right (int *nmismatches3, T ref, Compress_T query_compress,
			     Univcoord_T univdiagonal, int querylength,
			     int pos5, int pos3, bool plusp, int genestrand,
			     bool query_unk_mismatch_p, int kmer);
extern void
Genomebits_kmer_setup (bool genome_unk_mismatch_p_in, Mode_T mode);

#undef T
#endif


