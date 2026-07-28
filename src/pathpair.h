/* $Id: b7ffa7b67341307ad30875ec94c6789ae4b56199 $ */
#ifndef PATHPAIR_INCLUDED
#define PATHPAIR_INCLUDED

typedef struct Pathpair_T *Pathpair_T;

#include "path.h"

#include "intlistpool.h"
#include "listpool.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "transcriptpool.h"
#include "hitlistpool.h"

#include "bool.h"
#include "univcoord.h"
#include "univdiag.h"
#include "compress.h"
#include "shortread.h"

#include "stage1hr.h"
#include "knownsplicing.h"
#include "knownindels.h"

#include "resulthr.h"
#include "outputtype.h"


#define T Pathpair_T

struct T {
  /* MAPQ scores */
  float mapq_loglik;
  int mapq_score;
  int absmq_score;

#if 0
  int nmismatches;		/* querylength - sum of nmatches */
  int score_eventrim;		/* for storage */
  int alts_status_inside;
#endif

  /* int genestrand; */
  /* int sensedir; */
  /* Chrpos_T insertlength; */
  /* Chrpos_T outerlength; */

  Pairtype_T pairtype;
  int pair_relationship;
  Chrpos_T insertlength;
  Chrpos_T outerlength;

  bool plusp;
  bool transcript_concordant_p;

  Path_T path5;			/* Always a copy from the original */
  Path_T path3;			/* Always a copy from the original */

  /* For concordant plus reads, L=5, H=3; for minus, L=3, H=5 */
  Path_T pathL;			/* Pointer to path5 or path3 */
  Path_T pathH;			/* Pointer to path5 or path3 */

  Shortread_T queryseqL;
  Shortread_T queryseqH;
};


static inline int
Pathpair_found_score (T this) {
  return this->path5->found_score + this->path3->found_score;
}

extern bool
Pathpair_insertlength_knownp (T this);

extern bool
Pathpair_outerlength_knownp (T this);

extern void
Pathpair_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool);

extern void
Pathpair_gc (List_T *list, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool);

extern bool
Pathpair_transcript_intersectp (Path_T pathL, Path_T pathH);

extern bool
Pathpair_overlappingp (T this);

extern Chrpos_T
Pathpair_insertlength (T this);

extern Chrpos_T
Pathpair_outerlength (T this);

extern int
Pathpair_nbadsplices (T this);

static inline int
Pathpair_sensedir (T this) {
  return this->pathL->sensedir;
}

static inline double
Pathpair_total_splice_prob (T this) {
  return this->pathL->total_splice_prob + this->pathH->total_splice_prob;
}

static inline double
Pathpair_avg_junction_splice_prob (T this) {
  return Path_avg_junction_splice_prob(this->pathL) + Path_avg_junction_splice_prob(this->pathH);
}

static inline double
Pathpair_end_splice_prob (T this) {
  return Path_end_splice_prob(this->pathL) + Path_end_splice_prob(this->pathH);
}

extern void
Pathpair_print (T this);

extern T
Pathpair_new_concordant (List_T *unresolved_pathpairs,
			 Path_T pathL, Path_T pathH, Shortread_T queryseqL, Shortread_T queryseqH, bool plusp,

			 int nmismatches_filter_5, int nmismatches_filter_3,
			 int mincoverage_filter_5, int mincoverage_filter_3,

			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			 bool check_inner_p, bool copyLp, bool copyHp);

extern T
Pathpair_new_inner_fusion (Path_T pathL, Path_T pathH, Shortread_T queryseqL, Shortread_T queryseqH, bool plusp,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
			   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool,
			   bool copyLp, bool copyHp);

extern void
Pathpair_resolve (int *found_score_5, int *found_score_3,
		  T this, bool plusp, int genestrand,
		  Compress_T query_compress_L, Compress_T queryL_compress_fwd, Compress_T queryL_compress_rev,
		  Compress_T query_compress_H, Compress_T queryH_compress_fwd, Compress_T queryH_compress_rev,
		  Shortread_T queryseqL, Shortread_T queryseqH, char *queryptrL, char *queryptrH,
		  int querylengthL, int querylengthH, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		  int localdb_nmismatches_allowed,
		  
		  Stage1_T stage1L, Stage1_T stage1H, Knownsplicing_T knownsplicing,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Univdiagpool_T univdiagpool, Vectorpool_T vectorpool);

extern int
Pathpair_interval_cmp (const void *x, const void *y);
extern int
Pathpair_structure_cmp (const void *x, const void *y);
extern int
Pathpair_structure_ignore_sense_cmp (const void *x, const void *y);

extern bool
Pathpair_overlap_p (T x, T y);

extern bool
Pathpair_print_end_alignment (Filestring_T fp, Result_T result, Resulttype_T resulttype,
			      char initchar, bool firstp, 
			      Shortread_T queryseq, Shortread_T headerseq1, Shortread_T headerseq2,
			      int maxpaths, bool quiet_if_excessive_p, bool invertp, int quality_shift,
			      Listpool_T listpool);

extern bool
Pathpair_print_end_m8 (Filestring_T fp, Result_T result, Resulttype_T resulttype,
		       bool firstp, char *accession,
		       int maxpaths, bool quiet_if_excessive_p, bool invertp,
		       Listpool_T listpool);

extern void
Pathpair_setup (Outputtype_T output_type_in, Univcoord_T genomelength_in,
		Chrpos_T shortsplicedist_in);

#undef T
#endif

