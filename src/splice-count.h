/* $Id$ */
#ifndef SPLICE_COUNT_INCLUDED
#define SPLICE_COUNT_INCLUDED

typedef struct Splice_count_T *Splice_count_T;

#include "bool.h"

#define T Splice_count_T
struct T {
  bool knownp;
  int nunique;			/* Number of uniquely mapping reads crossing junction */
  int nmulti;			/* Number of multi-mapping reads crossing junction */
  int support;			/* Maximum alignment overhang */
};


extern T
Splice_count_new (bool knownp);

extern void
Splice_count_free (T *old);

#undef T
#endif

