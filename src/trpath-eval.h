#ifndef TRPATH_EVAL_INCLUDED
#define TRPATH_EVAL_INCLUDED

#include "trpath.h"

#include "bool.h"
#include "compress.h"

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"

#include "genomebits.h"


#define T Trpath_T

typedef struct Trpath_score_T {
  int matches;
  int mismatches;
  int indels;
  int indel_opens;
  int indel_extends;
  int softclip;
  int coverage;
  bool transcript_discordant_p;
  int found_score;
  int end_score;
  double score;
} Trpath_score_T;

extern int
Trpath_eval_score (Uintlist_T trdiagonals, Intlist_T tr_endpoints, Intlist_T nmismatches_list, List_T junctions,
		   Compress_T query_compress_tr, int querylength, bool tplusp);

extern void
Trpath_score_read (Trpath_score_T *components, T this, int querylength, Compress_T query_compress_tr,
		   double match_reward, double mismatch_penalty,
		   double indel_open_penalty, double indel_extend_penalty,
		   double softclip_penalty);

extern double
Trpath_score_pair (Trpath_score_T *pair_components, Trpath_score_T *components5,
		   Trpath_score_T *components3, T trpath5, T trpath3,
		   int querylength5, int querylength3,
		   Compress_T query5_compress_tr, Compress_T query3_compress_tr,
		   double match_reward, double mismatch_penalty,
		   double indel_open_penalty, double indel_extend_penalty,
		   double softclip_penalty, double transcript_discordance_penalty);

extern void
Trpath_eval_nmatches (int *found_score, int *end_score, T this, int querylength, Compress_T query_compress_tr);

extern void
Trpath_eval_setup (Genomebits_T transcriptome_in);

#undef T
#endif
