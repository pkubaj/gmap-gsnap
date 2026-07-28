#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "spliceendsgen.h"
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "list.h"

/* Used to keep a pool of allocated Spliceends_T objects, which are
   used for recursive calls to extend the ends of extend the starts or
   ends of an alignment, and thereby reduce the number of allocations
   and frees for each thread. */

/* Re-uses memory, as opposed to spliceendspool, and needed to avoid
   getting out-of-memory errors */

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#define T Spliceendsgen_T
struct T {
  int nallocated;
  int ncheckedout;
  List_T pool;
};



T
Spliceendsgen_new () {
  T new = (T) MALLOC(sizeof(*new));

  new->nallocated = 0;
  new->ncheckedout = 0;
  new->pool = (List_T) NULL;

  return new;
}


void
Spliceendsgen_reset (T this) {
  this->nallocated = 0;
  this->ncheckedout = 0;
  this->pool = (List_T) NULL;
  return;
}


void
Spliceendsgen_free_memory (T this) {
  Spliceends_T item;
  List_T p;

  debug0(printf("At Spliceendsgen_free_memory, %p has %d checked out\n",
		this,this->ncheckedout));

  assert(this->ncheckedout == 0);

  for (p = this->pool; p != NULL; p = List_next(p)) {
    item = (Spliceends_T) List_head(p);
    debug0(printf("Spliceendsgen_free_memory: freeing spliceends %p (ID %d)\n",item,item->id));
    Spliceends_free(&item);
  }
  List_free(&this->pool);

  return;
}


void
Spliceendsgen_free (T *old) {
  FREE(*old);
  return;
}


/* Use instead of Spliceends_new */
Spliceends_T
Spliceendsgen_checkout (T this, int querylength, Vectorpool_T vectorpool) {
  Spliceends_T new, old;
  List_T p;

#ifdef DEBUG0
  static int call_i = 0;
  printf("%d: %p Spliceendsgen_checkout: ",++call_i,this);
#endif

  if (this->ncheckedout == this->nallocated) {
    new = Spliceends_new(/*id*/this->nallocated++,querylength,vectorpool); /* Already in checkedout state */
    debug0(printf("creating %p (ID %d)\n",new,new->id));
    this->pool = List_push(this->pool,(void *) new);
    this->ncheckedout++;
    return new;

  } else {
    for (p = this->pool; p != NULL; p = List_next(p)) {
      old = (Spliceends_T) List_head(p);
      if (old->checkedout_p == false) {
	old->checkedout_p = true;
	this->ncheckedout++;
	debug0(printf("re-using %p (ID %d)\n",old,old->id));
	return old;
      }
    }

    fprintf(stderr,"Problem with Spliceendsgen_checkout: %d checked out, %d allocated, but cannot find it\n",
	    this->ncheckedout,this->nallocated);
    abort();
    return (Spliceends_T) NULL;
  }
}


/* Use instead of Spliceends_free */
void
Spliceendsgen_return (T this, Spliceends_T *old) {
  Spliceends_T item;
  List_T p;

#ifdef DEBUG0
  static int call_i = 0;
#endif

  if (*old) {
    debug0(printf("%d: %p Spliceendsgen_return: ",++call_i,this));
    for (p = this->pool; p != NULL; p = List_next(p)) {
      item = (Spliceends_T) List_head(p);
      if (item == *old) {
	debug0(printf("returning %p (ID %d)\n",*old,(*old)->id));
	(*old)->checkedout_p = false;
	this->ncheckedout--;

	/* Prevents problems with a subsequent attempt to return this item */
	*old = (Spliceends_T) NULL;

	return;
      }
    }

    fprintf(stderr,"Problem with Spliceendsgen_return: cannot find spliceends %p returned\n",
	    *old);
    abort();
  }

  return;
}
    

