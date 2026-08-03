#ifndef AUXINFO_INCLUDED
#define AUXINFO_INCLUDED

typedef struct Auxinfo_T *Auxinfo_T;

#include "bool.h"
#include "method.h"
#include "list.h"

#include "univdiagpool.h"
#include "auxinfopool.h"
#include "listpool.h"


/* Stores intermediate information, in addition to the
   middle_univdiagonal, needed to solve a path */

#define CHAIN_NONE 0
#define CHAIN_LEFT -1
#define CHAIN_RIGHT +1

#define T Auxinfo_T
struct T {
  Method_T method;
  int qstart;
  int qend;

  /* From Segmentelt_T procedure, by extension search */
  List_T *chains_by_nsplices;
  int nsplices;			/* Not valid after we merge Auxinfo_T objects */
  int max_nsplices;
  int chaindir;

  /* From Singleelt_T procedure, by segment search */
  int elti;			/* index into stage1->exhaustive_eltarray_{gplus,gminus} */

  struct Auxinfo_T *rest;
};


extern void
Auxinfo_print (T this);

extern T
Auxinfo_pop (T old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	     Listpool_T listpool);

extern void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	      Listpool_T listpool, bool free_chains_p);

extern void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool,
	    Listpool_T listpool);

extern T
Auxinfo_new (Method_T method, int qstart, int qend, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_tr (Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_chains (Method_T method, int qstart, int qend,
		    List_T chains, int chaindir, int nsplices, int max_nsplices,
		    Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_elt (Method_T method, int qstart, int qend, int elti,
		 Auxinfopool_T auxinfopool);

extern T
Auxinfo_append (T new, T old);

#undef T
#endif


