/* $Id: 059e1b94f91de630da3007d6d1ac5f9fa343f889 $ */
#ifndef AUXINFO_INCLUDED
#define AUXINFO_INCLUDED

typedef struct Auxinfo_T *Auxinfo_T;

#include "bool.h"
#include "method.h"
#include "list.h"

#include "univdiagpool.h"
#include "auxinfopool.h"


#define T Auxinfo_T
struct T {
  Method_T method;
  
  int qstart;
  int qend;
  int nmismatches;

  List_T right_univdiags;
  List_T left_univdiags;

  struct Auxinfo_T *rest;
};


extern void
Auxinfo_print (T this);

extern T
Auxinfo_pop (T old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool);

extern void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool);

extern void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new (Method_T method, int qstart, int qend, Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_tr (Auxinfopool_T auxinfopool);

extern T
Auxinfo_new_univdiags (Method_T method, int qstart, int qend, int nmismatches,
		       List_T right_univdiags, List_T left_univdiags, Auxinfopool_T auxinfopool);

extern T
Auxinfo_append (T new, T old);

#undef T
#endif


