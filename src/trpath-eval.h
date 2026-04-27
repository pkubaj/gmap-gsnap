#ifndef TRPATH_EVAL_INCLUDED
#define TrPATH_EVAL_INCLUDED

#include "trpath.h"

#include "bool.h"
#include "compress.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"

#include "genomebits.h"


#define T Trpath_T

extern int
Trpath_eval_score (Uintlist_T trdiagonals, Intlist_T tr_endpoints, Intlist_T nmismatches_list, List_T junctions,
		   Compress_T query_compress_tr, int querylength, bool tplusp);

extern void
Trpath_eval_nmatches (int *found_score, T this, int querylength, Compress_T query_compress_tr);

extern void
Trpath_eval_setup (Genomebits_T transcriptome_in);

#undef T
#endif


