/* $Id$ */
#ifndef PATHSTORE_INCLUDED
#define PATHSTORE_INCLUDED

typedef struct Pathstore_T *Pathstore_T;

#include "bool.h"
#include "path.h"
#include "list.h"

#include "ef64.h"
#include "chrnum.h"

#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "pathstorepool.h"

#include "transcriptpool.h"
#include "hitlistpool.h"


#define T Pathstore_T
struct T {
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;

  bool solvedp;
  bool complete_sense_p;		/* tells whether best_sense_paths came from complete (as opposed to unextended) */
  bool complete_antisense_p;		/* tells whether best_antisense_paths came from complete (as opposed to unextended) */

  Intlist_T best_sense_partners;	/* index into other all_univdiagonals/pathstore array */
  Intlist_T best_antisense_partners;	/* index into other all_univdiagonals/pathstore array */

  List_T best_sense_paths;
  List_T best_antisense_paths;

  List_T unextended_sense_paths;
  List_T unextended_antisense_paths;
  List_T complete_sense_paths;
  List_T complete_antisense_paths;

  int mate_count;
  int mate_bestk;
};


extern void
Pathstore_print (T this);

extern void
Pathstore_free (T *old, Pathstorepool_T pathstorepool,
		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, Pathpool_T pathpool,
		Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Pathstore_gc (T *array, int n, Pathstorepool_T pathstorepool,
	      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	      Listpool_T listpool, Pathpool_T pathpool,
	      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);

extern void
Pathstore_set_best_paths (T this, Hitlistpool_T hitlistpool);
extern void
Pathstore_set_best_sense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p);
extern void
Pathstore_set_best_antisense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p);

extern T
Pathstore_merge (T old, T *new, Pathstorepool_T pathstorepool, Intlistpool_T intlistpool);

extern T
Pathstore_new (Pathstorepool_T pathstorepool);

extern T
Pathstore_new_tr (Pathstorepool_T pathstorepool, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh);

extern void
Pathstore_collect_paths (bool *foundp, List_T *sense_paths, List_T *antisense_paths,
		       T *pathstore_array, int nunivdiagonals,
		       Hitlistpool_T hitlistpool);

extern void
Pathstore_assign_chrinfo (Univcoord_T * univdiagonals, T *pathstore, int n, int querylength);

extern void
Pathstore_setup (EF64_T chromosome_ef64_in);

#undef T
#endif


