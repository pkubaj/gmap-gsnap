static char rcsid[] = "$Id: 20a23498d69176ec26f2aabd4447cebe065f1c46 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "intlistpool.h"
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


#define T Intlistpool_T
struct T {
  List_T chunks;
#ifdef INTLISTPOOL_REUSE
  List_T free_cells;
#else
  struct Intlist_T *cellptr;
  int cellctr;
#endif
};


#if defined(INTLISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_chunk (T this, struct Intlist_T *free_cell) {
  int chunki;
  List_T p;
  struct Intlist_T *chunk;

  for (p = this->chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Intlist_T *) List_head(p);
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
      chunki = find_chunk(this,(struct Intlist_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Intlist_T cells leaked in intlistpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Intlistpool_reset_memory (T this) {
  struct Intlist_T *chunk;

#if defined(INTLISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_memory(this);
#endif

  while (List_next(this->chunks) != NULL) {
    this->chunks = List_pop_keep(this->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef INTLISTPOOL_REUSE
  int celli;
  chunk = (struct Intlist_T *) List_head(this->chunks);
  List_free_keep(&this->free_cells);
  this->free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#else
  this->cellptr = (struct Intlist_T *) List_head(this->chunks);
  this->cellctr = 0;
#endif

  return;
}


void
Intlistpool_free (T *old) {
  struct Intlist_T *chunk;

  while ((*old)->chunks != NULL) {
    (*old)->chunks = List_pop_keep((*old)->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef INTLISTPOOL_REUSE
  List_free_keep(&(*old)->free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Intlistpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->chunks = (List_T) NULL;
#ifdef INTLISTPOOL_REUSE
  new->free_cells = (List_T) NULL;
#else
  new->cellctr = 0;
#endif

  return new;
}


#ifdef INTLISTPOOL_REUSE
void
Intlistpool_free_list (Intlist_T *list, T this
#ifdef INTLISTPOOL_TRACE
		       , const char *file, int line
#endif
		       ) {
  Intlist_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->free_cells = List_push_keep(this->free_cells,(void *) prev);
#ifdef INTLISTPOOL_TRACE
    printf("Intlistpool: Freed %p -- Intlistpool_free_list called by %s:%d, contents %d\n",
	   prev,file,line,prev->first);
#endif
  }

  *list = (Intlist_T) NULL;

  return;
}
#endif


#ifdef INTLISTPOOL_REUSE
Intlist_T
Intlistpool_pop (Intlist_T list, T this, int *integer
#ifdef INTLISTPOOL_TRACE
		 , const char *file, int line
#endif
		 ) {
  Intlist_T head;

  if (list != NULL) {
    this->free_cells = List_push_keep(this->free_cells,(void *) list);
#ifdef INTLISTPOOL_TRACE
    printf("Intlistpool: Freed %p -- Intlistpool_pop called by %s:%d, contents %d\n",
	   list,file,line,list->first);
#endif
    head = list->rest;
    *integer = list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct Intlist_T *
add_new_chunk (T this) {
  struct Intlist_T *chunk;

  chunk = (struct Intlist_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Intlist_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
#ifdef INTLISTPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of listcells.  Ptr for chunk %d is %p\n",
	       List_length(this->chunks),chunk));

  return chunk;
}


Intlist_T
Intlistpool_push (Intlist_T list, T this, int integer
#ifdef INTLISTPOOL_TRACE
		  , const char *file, int line
#endif
		  ) {
  Intlist_T listcell;

#ifdef INTLISTPOOL_REUSE
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

#ifdef INTLISTPOOL_TRACE
  printf("Intlistpool: Allocated %p -- Intlistpool_push called by %s:%d, contents %d\n",
	 listcell,file,line,integer);
#endif

  listcell->first = integer;
  listcell->rest = list;

  return listcell;
}


Intlist_T
Intlistpool_copy (Intlist_T source, T this) {
  Intlist_T dest = NULL;

  while (source != NULL) {
    dest = Intlistpool_push(dest,this,/*orig*/source->first
			    intlistpool_trace(__FILE__,__LINE__));
    source = source->rest;
  }
  return Intlist_reverse(dest);
}


/* Guarantees that a chunk is available */
void
Intlistpool_init (T this) {
#ifdef INTLISTPOOL_REUSE
  add_new_chunk(this);
#else
  this->cellptr = add_new_chunk(this);;
  this->cellctr = 0;
#endif

  return;
}

