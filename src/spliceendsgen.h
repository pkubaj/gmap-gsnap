/* $Id: 76d1c8d05e482fb41af5618f67b7a22884eca711 $ */
#ifndef SPLICEENDSGEN_INCLUDED
#define SPLICEENDSGEN_INCLUDED

/* Re-uses memory, as opposed to spliceendspool, and needed to avoid
   getting out-of-memory errors */

typedef struct Spliceendsgen_T *Spliceendsgen_T;

#include "spliceends.h"
#include "vectorpool.h"

#define T Spliceendsgen_T

extern T
Spliceendsgen_new ();
extern void
Spliceendsgen_reset (T this);
extern void
Spliceendsgen_free_memory (T this);
extern void
Spliceendsgen_free (T *old);

extern Spliceends_T
Spliceendsgen_checkout (T this, int querylength, Vectorpool_T vectorpool);
extern void
Spliceendsgen_return (T this, Spliceends_T *old);

#undef T
#endif



