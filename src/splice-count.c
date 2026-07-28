static char rcsid[] = "$Id$";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "splice-count.h"

#include "mem.h"

#define T Splice_count_T

T
Splice_count_new (bool knownp) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->knownp = knownp;
  new->nunique = 0;
  new->nmulti = 0;
  new->support = 0;

  return new;
}


void
Splice_count_free (T *old) {

  if (*old) {
    FREE_KEEP(*old);
  }

  return;
}

    
