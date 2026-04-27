/* $Id: ce331e166dee84bf6ddeed7a869b95c6c6c09000 $ */
#ifndef REPAIR_INCLUDED
#define REPAIR_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct Repair_T *Repair_T;

#include "bool.h"
#include "genomicpos.h"
#include "list.h"

#include "sense.h"
#include "compress.h"
#include "path.h"
#include "transcript.h"
#include "genomebits.h"

#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"


#define T Repair_T
struct T {
  int transcript_genestrand;
  int trstart_overhang;
  int trend_overhang;
  Chrpos_T trstart_splice_distance;
  Chrpos_T trend_splice_distance;
  List_T transcripts;
};


extern void
Repair_free (T *old, Listpool_T listpool, Transcriptpool_T transcriptpool, bool free_transcripts_p);
extern T
Repair_new (int transcript_genestrand,
	    int trstart_overhang, Chrpos_T trstart_splice_distance,
	    int trend_overhang, Chrpos_T trend_splice_distance,
	    Transcript_T transcript, Listpool_T listpool);
extern List_T
Repair_make_unique (List_T repairs, Listpool_T listpool, Transcriptpool_T transcriptpool);
extern Path_T
Repair_path (int *found_score, T this, Path_T oldpath, int sensedir,
	     Compress_T query_compress,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
	     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Repair_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in);

#undef T
#endif

