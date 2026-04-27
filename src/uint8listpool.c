static char rcsid[] = "$Id: db75e77b2b8270031a0066c8659c727ee184ee15 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "uint8listpool.h"
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


#define T Uint8listpool_T
struct T {
  List_T chunks;
#ifdef UINT8LISTPOOL_REUSE
  List_T free_cells;
#else
  struct Uint8list_T *cellptr;
  int cellctr;
#endif
};


#if defined(UINT8LISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_chunk (T this, struct Uint8list_T *free_cell) {
  int chunki;
  List_T p;
  struct Uint8list_T *chunk;

  for (p = this->chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Uint8list_T *) List_head(p);
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
      chunki = find_chunk(this,(struct Uint8list_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Uint8list_T cells leaked in uint8listpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Uint8listpool_reset_memory (T this) {
  struct Uint8list_T *chunk;

#if defined(UINT8LISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_memory(this);
#endif

  while (List_next(this->chunks) != NULL) {
    this->chunks = List_pop_keep(this->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef UINT8LISTPOOL_REUSE
  int celli;
  chunk = (struct Uint8list_T *) List_head(this->chunks);
  List_free_keep(&this->free_cells);
  this->free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#else
  this->cellptr = (struct Uint8list_T *) List_head(this->chunks);
  this->cellctr = 0;
#endif

  return;
}

void
Uint8listpool_free (T *old) {
  struct Uint8list_T *chunk;

  while ((*old)->chunks != NULL) {
    (*old)->chunks = List_pop_keep((*old)->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef UINT8LISTPOOL_REUSE
  List_free_keep(&(*old)->free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Uint8listpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->chunks = (List_T) NULL;
#ifdef UINT8LISTPOOL_REUSE
  new->free_cells = (List_T) NULL;
#else
  new->cellctr = 0;
#endif

  return new;
}


#ifdef UINT8LISTPOOL_REUSE
void
Uint8listpool_free_list (Uint8list_T *list, T this
#ifdef UINT8LISTPOOL_TRACE
			 , const char *file, int line
#endif
			 ) {
  Uint8list_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->free_cells = List_push_keep(this->free_cells,(void *) prev);
#ifdef UINT8LISTPOOL_TRACE
    printf("Uint8listpool: Freed %p -- Uint8listpool_free_list called by %s:%d\n",prev,file,line);
#endif
  }

  *list = (Uint8list_T) NULL;

  return;
}
#endif


#ifdef UINT8LISTPOOL_REUSE
Uint8list_T
Uint8listpool_pop (Uint8list_T list, T this, UINT8 *integer
#ifdef UINT8LISTPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Uint8list_T head;

  if (list != NULL) {
    this->free_cells = List_push_keep(this->free_cells,(void *) list);
#ifdef UINT8LISTPOOL_TRACE
    printf("Uint8listpool: Freed %p -- Uint8listpool_pop called by %s:%d\n",list,file,line);
#endif
    head = list->rest;
    *integer = list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct Uint8list_T *
add_new_chunk (T this) {
  struct Uint8list_T *chunk;

  chunk = (struct Uint8list_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Uint8list_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
#ifdef UINT8LISTPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of listcells.  Ptr for chunk %d is %p\n",
	       List_length(this->chunks),chunk));

  return chunk;
}


Uint8list_T
Uint8listpool_push (Uint8list_T list, T this, UINT8 integer
#ifdef UINT8LISTPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  Uint8list_T listcell;

#ifdef UINT8LISTPOOL_REUSE
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

#ifdef UINT8LISTPOOL_TRACE
  printf("Uint8listpool: Allocated %p -- Uint8listpool_push called by %s:%d\n",listcell,file,line);
#endif

  listcell->first = integer;
  listcell->rest = list;

  return listcell;
}

Uint8list_T
Uint8listpool_copy (Uint8list_T source, T this) {
  Uint8list_T dest = NULL;

  while (source != NULL) {
    dest = Uint8listpool_push(dest,this,/*orig*/source->first);
    source = source->rest;
  }
  return Uint8list_reverse(dest);
}


/* Guarantees that a chunk is available */
void
Uint8listpool_init (T this) {
#ifdef UINT8LISTPOOL_REUSE
  add_new_chunk(this);
#else
  this->cellptr = add_new_chunk(this);
  this->cellctr = 0;
#endif

  return;
}

