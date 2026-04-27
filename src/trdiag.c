static char rcsid[] = "$I$";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "trdiag.h"
#include "trdiagdef.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "assert.h"


#define T Trdiag_T


#if 0
/* Now in trdiagpool.c */
T
Trdiag_new (int qstart, int qend, Trcoord_T trdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  assert(qend > qstart);

  new->trdiagonal = trdiagonal;
  new->qstart = qstart;
  new->qend = qend;
  new->nmismatches = 0;
  /* new->nconsecutive = qend - qstart + 1; */

  return new;
}
#endif


#if 0
T
Trdiag_copy (T this) {
  return Trdiag_new(this->qstart,this->qend,this->trdiagonal);
}
#endif


#if 0
List_T
Trdiag_copy_list (List_T list) {
  List_T new = NULL, p;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);

    new = List_push(new,(void *) Trdiag_new(this->qstart,this->qend,this->trdiagonal));
  }
  
  return List_reverse(new);
}
#endif
  


#if 0
T
Trdiag_new_fillin (int qstart, int qend, int indexsize, Trcoord_T trdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->trdiagonal = trdiagonal;
  new->qstart = qstart;
  new->qend = qend + indexsize - 1;
  /* new->nconsecutive = new->qend - qstart + 1; */

  new->intscore = -1;

  return new;
}
#endif


#if 0
/* All Trdiag_T objects and lists now allocated in trdiagpool.h */
void
Trdiag_free (T *old) {
  FREE(*old);
  return;
}
#endif

#if 0
/* All Trdiag_T objects and lists now allocated in trdiagpool.h */
void
Trdiag_gc (List_T *list) {
  T trdiagonal;
  List_T p;

  for (p = *list; p != NULL; p = List_next(p)) {
    trdiagonal = (T) List_head(p);
    FREE(trdiagonal);
  }
  List_free(&(*list));
  return;
}
#endif

#if 0
/* All Trdiag_T objects and lists now allocated in trdiagpool.h */
void
Trdiag_list_gc (List_T *paths) {
  List_T path, p;

  for (p = *paths; p != NULL; p = List_next(p)) {
    path = (List_T) List_head(p);
    Trdiag_gc(&path);
  }
  List_free(&(*paths));
  return;
}
#endif


int
Trdiag_list_length (List_T path) {
  int length = 0;
  T this;
  List_T p;

  for (p = path; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    length += this->qend - this->qstart + 1;
  }
  return length;
}


int
Trdiag_ascending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->qstart < y->qstart) {
    return -1;
  } else if (y->qstart < x->qstart) {
    return +1;
  } else if (x->qend < y->qend) {
    return -1;
  } else if (y->qend < x->qend) {
    return +1;
  } else if (x->trdiagonal < y->trdiagonal) {
    return -1;
  } else if (y->trdiagonal < x->trdiagonal) {
    return +1;
  } else {
    return 0;
  }
}


int
Trdiag_descending_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->qstart > y->qstart) {
    return -1;
  } else if (y->qstart > x->qstart) {
    return +1;
  } else if (x->qend > y->qend) {
    return -1;
  } else if (y->qend > x->qend) {
    return +1;
  } else if (x->trdiagonal > y->trdiagonal) {
    return -1;
  } else if (y->trdiagonal > x->trdiagonal) {
    return +1;
  } else {
    return 0;
  }
}



int
Trdiag_diagonal_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->trdiagonal < y->trdiagonal) {
    return -1;
  } else if (y->trdiagonal < x->trdiagonal) {
    return +1;
  } else if (x->qstart < y->qstart) {
    return -1;
  } else if (y->qstart < x->qstart) {
    return +1;
  } else if (x->qend < y->qend) {
    return -1;
  } else if (y->qend < x->qend) {
    return +1;
  } else {
    return 0;
  }
}
