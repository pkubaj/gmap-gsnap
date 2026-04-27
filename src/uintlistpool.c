static char rcsid[] = "$Id: 2db0f48b6ccda3455bc25cdfe28c7452e0fd4204 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "uintlistpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define CHUNKSIZE 16384

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Uintlistpool_T
struct T {
  List_T chunks;
#ifdef UINTLISTPOOL_REUSE
  List_T free_cells;
#else
  struct Uintlist_T *cellptr;
  int cellctr;
#endif
};


#if defined(UINTLISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_chunk (T this, struct Uintlist_T *free_cell) {
  int chunki;
  List_T p;
  struct Uintlist_T *chunk;

  for (p = this->chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Uintlist_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}
  

/* Checks to see if all memory has been returned to free_cells */
static void
check_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->free_cells) == CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->free_cells; p != NULL; p = List_next(p)) {
      chunki = find_chunk(this,(struct Uintlist_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Uintlist_T cells leaked in uintlistpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Uintlistpool_reset_memory (T this) {
  struct Uintlist_T *chunk;

#if defined(UINTLISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_memory(this);
#endif

  while (List_next(this->chunks) != NULL) {
    this->chunks = List_pop_keep(this->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef UINTLISTPOOL_REUSE
  int celli;
  chunk = (struct Uintlist_T *) List_head(this->chunks);
  List_free_keep(&this->free_cells);
  this->free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#else
  this->cellptr = (struct Uintlist_T *) List_head(this->chunks);
  this->cellctr = 0;
#endif

  return;
}

void
Uintlistpool_free (T *old) {
  struct Uintlist_T *chunk;

  while ((*old)->chunks != NULL) {
    (*old)->chunks = List_pop_keep((*old)->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef UINTLISTPOOL_REUSE
  List_free_keep(&(*old)->free_cells);
#endif

  FREE_KEEP(*old);

  return;
}

T
Uintlistpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->chunks = (List_T) NULL;
#ifdef UINTLISTPOOL_REUSE
  new->free_cells = (List_T) NULL;
#else
  new->cellctr = 0;
#endif

  return new;
}


#ifdef UINTLISTPOOL_REUSE
void
Uintlistpool_free_list (Uintlist_T *list, T this
#ifdef UINTLISTPOOL_TRACE
			, const char *file, int line
#endif
			) {
  Uintlist_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->free_cells = List_push_keep(this->free_cells,(void *) prev);
#ifdef UINTLISTPOOL_TRACE
    printf("Uintlistpool: Freed %p -- Uintlistpool_free_list called by %s:%d\n",prev,file,line);
#endif
  }

  *list = (Uintlist_T) NULL;

  return;
}
#endif


#ifdef UINTLISTPOOL_REUSE
Uintlist_T
Uintlistpool_pop (Uintlist_T list, T this, unsigned int *integer
#ifdef UINTLISTPOOL_TRACE
		  , const char *file, int line
#endif
		  ) {
  Uintlist_T head;

  if (list != NULL) {
    this->free_cells = List_push_keep(this->free_cells,(void *) list);
#ifdef UINTLISTPOOL_TRACE
    printf("Uintlistpool: Freed %p -- Uintlistpool_pop called by %s:%d\n",list,file,line);
#endif
    head = list->rest;
    *integer = list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct Uintlist_T *
add_new_chunk (T this) {
  struct Uintlist_T *chunk;

  chunk = (struct Uintlist_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Uintlist_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
#ifdef UINTLISTPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of listcells.  Ptr for chunk %d is %p\n",
	       List_length(this->chunks),chunk));

  return chunk;
}


Uintlist_T
Uintlistpool_push (Uintlist_T list, T this, unsigned int integer
#ifdef UINTLISTPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Uintlist_T listcell;

#ifdef UINTLISTPOOL_REUSE
  if (this->free_cells == (List_T) NULL) {
    add_new_chunk(this);
  }
  this->free_cells = List_pop_keep(this->free_cells,(void **) &listcell);
#else
  if (this->cellctr >= CHUNKSIZE) {
    this->cellptr = add_new_chunk(this);
    this->cellctr = 0;
  }
  this->cellctr += 1;
  listcell = this->cellptr++;
#endif

#ifdef UINTLISTPOOL_TRACE
  printf("Uintlistpool: Allocated %p -- Uintlistpool_push called by %s:%d\n",listcell,file,line);
#endif

  listcell->first = integer;
  listcell->rest = list;

  return listcell;
}

Uintlist_T
Uintlistpool_copy (Uintlist_T source, T this) {
  Uintlist_T dest = NULL;

  while (source != NULL) {
    dest = Uintlistpool_push(dest,this,/*orig*/source->first
			     uintlistpool_trace(__FILE__,__LINE__));
    source = source->rest;
  }
  return Uintlist_reverse(dest);
}


/* Guarantees that a chunk is available */
void
Uintlistpool_init (T this) {
#ifdef UINTLISTPOOL_REUSE
  add_new_chunk(this);
#else
  this->cellptr = add_new_chunk(this);;
  this->cellctr = 0;
#endif

  return;
}

