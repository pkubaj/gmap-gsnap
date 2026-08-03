#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "univdiag.h"
#include "univdiagdef.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mem.h"
#include "assert.h"


#define T Univdiag_T


T
Univdiag_new (int qstart, int qend, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  assert(qend > qstart);

  new->univdiagonal = univdiagonal;
  new->qstart = qstart;
  new->qend = qend;
  new->nmismatches = 0;

  return new;
}


#if 0
T
Univdiag_copy (T this) {
  return Univdiag_new(this->qstart,this->qend,this->univdiagonal);
}
#endif


#if 0
List_T
Univdiag_copy_list (List_T list) {
  List_T new = NULL, p;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);

    new = List_push(new,(void *) Univdiag_new(this->qstart,this->qend,this->univdiagonal));
  }
  
  return List_reverse(new);
}
#endif
  


#if 0
T
Univdiag_new_fillin (int qstart, int qend, int indexsize, Univcoord_T univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->qstart = qstart;
  new->qend = qend + indexsize - 1;
  /* new->nconsecutive = new->qend - qstart + 1; */

  new->intscore = -1;

  return new;
}
#endif


#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_free (T *old) {
  FREE(*old);
  return;
}
#endif

#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_gc (List_T *list) {
  T univdiagonal;
  List_T p;

  for (p = *list; p != NULL; p = List_next(p)) {
    univdiagonal = (T) List_head(p);
    FREE(univdiagonal);
  }
  List_free(&(*list));
  return;
}
#endif

#if 0
/* All Univdiag_T objects and lists now allocated in univdiagpool.h */
void
Univdiag_list_gc (List_T *paths) {
  List_T path, p;

  for (p = *paths; p != NULL; p = List_next(p)) {
    path = (List_T) List_head(p);
    Univdiag_gc(&path);
  }
  List_free(&(*paths));
  return;
}
#endif


void
Univdiag_transfer (T new, T old) {

  new->univdiagonal = old->univdiagonal;
  new->qstart = old->qstart;
  new->qend = old->qend;
  new->nmismatches = old->nmismatches;

  return;
}


void
Univdiag_transfer_list_to_array (T *array, List_T p) {

  while (p != NULL) {
    Univdiag_transfer(*array,(T) List_head(p));
    array++;
    p = List_next(p);
  }

  return;
}


/* Note: Cannot sort the array, or else we will lose the pointer to space */
void
Univdiag_free_array (Univdiag_T **array) {

  if (*array) {
    FREE((*array)[0]);
    FREE(*array);
  }

  return;
}


T *
Univdiag_new_array (int n) {
  T *array;
  struct T *space;
  int i;

  array = (T *) MALLOC(n*sizeof(T));
  space = (struct T *) MALLOC(n*sizeof(struct T));

  if (n > 0) {
    memset(space,0,n*sizeof(struct T));
  }

  array[0] = space;
  for (i = 1; i < n; i++) {
    array[i] = &(array[i-1][1]);
  }

  return array;
}


void
Univdiag_free_array_keep (Univdiag_T **array) {

  if (*array) {
    FREE_KEEP((*array)[0]);
    FREE_KEEP(*array);
  }

  return;
}


T *
Univdiag_new_array_keep (int n) {
  T *array;
  struct T *space;
  int i;

  array = (T *) MALLOC_KEEP(n*sizeof(T));
  space = (struct T *) MALLOC_KEEP(n*sizeof(struct T));

  if (n > 0) {
    memset(space,0,n*sizeof(struct T));
  }

  array[0] = space;
  for (i = 1; i < n; i++) {
    array[i] = &(array[i-1][1]);
  }

  return array;
}


/* Can't swap pointers because we need to keep array[0] to free space */
void
Univdiag_reverse_inplace (T *array, int starti, int endi) {
  struct T temp;
  int i, j, n = endi - starti;;

  for (i = starti, j = endi-1; i < starti + n/2; i++, j--) {
    Univdiag_transfer(&temp,array[i]);
    Univdiag_transfer(array[i],array[j]);
    Univdiag_transfer(array[j],&temp);
  }

#if 0
  if (i == j) {
    coords[i] = coords[j];
  }
#endif

  return;
}


int
Univdiag_list_length (List_T path) {
  int length = 0;
  T this;
  List_T p;

  for (p = path; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    length += this->qend - this->qstart + 1;
  }
  return length;
}


bool
Univdiag_equal (T x, T y) {

  if (x->univdiagonal != y->univdiagonal) {
    return false;
  } else if (x->qstart != y->qstart) {
    return false;
  } else if (x->qend != y->qend) {
    return false;
  } else {
    return true;
  }
}

int
Univdiag_ascending_cmp (const void *a, const void *b) {
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
  } else if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal < x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}


int
Univdiag_descending_cmp (const void *a, const void *b) {
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
  } else if (x->univdiagonal > y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal > x->univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}



#if 0
int
Univdiag_diagonal_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->univdiagonal < y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal < x->univdiagonal) {
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
#endif


int
Univdiag_struct_diagonal_fwd_cmp (const void *a, const void *b) {
  struct T x = * (struct T *) a;
  struct T y = * (struct T *) b;

  if (x.univdiagonal < y.univdiagonal) {
    return -1;
  } else if (y.univdiagonal < x.univdiagonal) {
    return +1;
  } else if (x.qstart < y.qstart) {
    return -1;
  } else if (y.qstart < x.qstart) {
    return +1;
  } else if (x.qend < y.qend) {
    return -1;
  } else if (y.qend < x.qend) {
    return +1;
  } else {
    return 0;
  }
}


#if 0
int
Univdiag_diagonal_rev_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->univdiagonal > y->univdiagonal) {
    return -1;
  } else if (y->univdiagonal > x->univdiagonal) {
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
#endif


int
Univdiag_struct_diagonal_rev_cmp (const void *a, const void *b) {
  struct T x = * (struct T *) a;
  struct T y = * (struct T *) b;

  if (x.univdiagonal > y.univdiagonal) {
    return -1;
  } else if (y.univdiagonal > x.univdiagonal) {
    return +1;
  } else if (x.qstart < y.qstart) {
    return -1;
  } else if (y.qstart < x.qstart) {
    return +1;
  } else if (x.qend < y.qend) {
    return -1;
  } else if (y.qend < x.qend) {
    return +1;
  } else {
    return 0;
  }
}


int
Univdiag_struct_cmp (const void *a, const void *b) {
  struct T x = * (struct T *) a;
  struct T y = * (struct T *) b;

  if (x.univdiagonal < y.univdiagonal) {
    return -1;
  } else if (y.univdiagonal < x.univdiagonal) {
    return +1;
  } else if (x.qstart < y.qstart) {
    return -1;
  } else if (y.qstart < x.qstart) {
    return +1;
  } else if (x.qend < y.qend) {
    return -1;
  } else if (y.qend < x.qend) {
    return +1;
  } else {
    return 0;
  }
}


void
Univdiag_chain_print_left (List_T chain) {
  T univdiag;
  Univcoord_T prev_univdiagonal;
  List_T p;

  univdiag = (T) List_head(chain);
  prev_univdiagonal = univdiag->univdiagonal;

  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (T) List_head(p);
    printf("left %d..%d %u dist:%d\n",
	   univdiag->qstart,univdiag->qend,univdiag->univdiagonal,
	   prev_univdiagonal - univdiag->univdiagonal);
    prev_univdiagonal = univdiag->univdiagonal;
  }
  printf("End of path\n\n");

  return;
}


void
Univdiag_chain_print_right (List_T chain) {
  T univdiag;
  Univcoord_T prev_univdiagonal;
  List_T p;

  univdiag = (T) List_head(chain);
  prev_univdiagonal = univdiag->univdiagonal;

  for (p = chain; p != NULL; p = List_next(p)) {
    univdiag = (T) List_head(p);
    printf("right %d..%d %u dist:%d\n",
	   univdiag->qstart,univdiag->qend,univdiag->univdiagonal,
	   univdiag->univdiagonal - prev_univdiagonal);
    prev_univdiagonal = univdiag->univdiagonal;
  }
  printf("End of path\n\n");

  return;
}

