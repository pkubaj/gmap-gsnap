/* $Id: 7d9a2c9e8a1a6d48069af6ef659ba9cf7c076c87 $ */
#ifndef VECTORPOOL_INCLUDED
#define VECTORPOOL_INCLUDED

typedef struct Vectorpool_T *Vectorpool_T;

#include "univcoord.h"

#define T Vectorpool_T

extern void
Vectorpool_reset_memory (T this);
extern void
Vectorpool_free (T *old);
extern T
Vectorpool_new (void);
extern int *
Vectorpool_new_intvector (T this, int nints);
extern unsigned int *
Vectorpool_new_uintvector (T this, int nuints);
extern Univcoord_T *
Vectorpool_new_univcoordvector (T this, int nunivcoords);
extern double *
Vectorpool_new_doublevector (T this, int ndoubles);
extern void
Vectorpool_init (T this);

#undef T
#endif
