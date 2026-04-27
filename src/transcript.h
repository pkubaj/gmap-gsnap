/* $Id: 8a03530984db306995b045782878b36ad5467712 $ */
#ifndef TRANSCRIPT_INCLUDED
#define TRANSCRIPT_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct Transcript_T *Transcript_T;

#include "filestring.h"
#include "list.h"
#include "intlist.h"
#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "outputtype.h"
#include "transcriptome.h"

#include "shortread.h"
#include "path.h"
#include "listpool.h"
#include "transcriptpool.h"

#include "transcript-velocity.h"


#define T Transcript_T
struct T {
  unsigned int num;	      /* Index in transcriptome */
  int genestrand;	/* geneplus/geneminus: > 0 if transcript is on plus strand and < 0 if it is on minus strand */

  int trstart;	      /* Starting coordinate on transcript, corresponding to querystart.  0-based */
  int trend; /* Ending coordinate, corresponding to queryend.  If transcript_plusp is false, trend > trstart */

  int trstart_overhang;		/* Subtract from trstart if repair works */
  int trend_overhang;		/* Add to trend if repair works */

  List_T exons;

  int nexons;			/* Needed for velocity, since a 1-exon gene cannot be unspliced */
  int trlength;

  Velocity_T velocity;
};


static inline unsigned int
Transcript_num (T this) {
  return this->num;
}

static inline int
Transcript_start (T this) {
  return this->trstart;
}

static inline int
Transcript_end (T this) {
  return this->trend;
}

static inline int
Transcript_genestrand (T this) {
  return this->genestrand;
}

extern void
Transcript_free (T *old, Listpool_T listpool, Transcriptpool_T transcriptpool);

extern void
Transcript_list_gc (List_T *list, Listpool_T listpool, Transcriptpool_T transcriptpool);

extern T
Transcript_new (unsigned int num, int transcript_genestrand, int trstart, int trend,
		int trstart_overhang, int trend_overhang,
		List_T exons, int nexons, int trlength, Transcriptpool_T transcriptpool);

extern T
Transcript_copy (T old, Transcriptpool_T transcriptpool, Listpool_T listpool);
extern List_T
Transcript_copy_list (List_T old, Transcriptpool_T transcriptpool, Listpool_T listpool);
extern void
Transcript_list_ascendingp (List_T list);
extern List_T
Transcript_list_sort (List_T transcripts, Listpool_T listpool, Transcriptpool_T transcriptpool);
extern void
Transcript_list_trim (List_T transcripts, int trim5, int trim3, Transcriptome_T transcriptome,
		      Listpool_T listpool, Transcriptpool_T transcriptpool);
extern void
Transcript_print_nums (List_T list);
extern void
Transcript_print_list (Filestring_T fp, List_T transcripts,
		       Univ_IIT_T transcript_iit, char *header);

extern bool
Transcript_intersectp (List_T transcripts5, List_T transcripts3);

extern bool
Transcript_intersection (Path_T path5, Path_T path3,
			 Listpool_T listpool, Transcriptpool_T transcriptpool);
extern int
Transcript_fragment_length (Path_T path5, Path_T path3, Shortread_T queryseq5, Shortread_T queryseq3);

extern void
Transcript_repair_trstart (T this, Listpool_T listpool, Transcriptpool_T transcriptpool);
extern void
Transcript_repair_trend (T this, Listpool_T listpool, Transcriptpool_T transcriptpool);
extern List_T
Transcript_remove_subset (List_T transcripts, List_T subset,
			  Listpool_T listpool, Transcriptpool_T transcriptpool);

extern void
Transcript_setup (int pairmax_transcriptome_in,
		  Outputtype_T output_type_in, Transcriptome_T transcriptome_in,
		  Univ_IIT_T transcript_iit_in);


#undef T
#endif

