static char rcsid[] = "$Id: 3df6f9871309f29a67b0c9b6cccb9890918fb3eb $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "trpathpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "list.h"

#define TRPATH_CHUNKSIZE 1024

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Trpathpool_T
struct T {
  List_T trpath_chunks;
#ifdef TRPATHPOOL_REUSE
  List_T trpath_free_cells;
#else
  struct Trpath_T *trpath_cellptr;
  int trpath_cellctr;
#endif
};


#if defined(TRPATHPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_trpath_chunk (T this, struct Trpath_T *free_cell) {
  int chunki;
  List_T p;
  struct Trpath_T *chunk;

  for (p = this->trpath_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Trpath_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[TRPATH_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}

/* Checks to see if all memory has been returned to free_cells */
static void
check_trpath_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->trpath_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->trpath_free_cells) == TRPATH_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->trpath_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_trpath_chunk(this,(struct Trpath_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < TRPATH_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Trpath_T cells leaked in trpathpool chunk %d\n",
		TRPATH_CHUNKSIZE - nfreecells[chunki],TRPATH_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Trpathpool_reset_memory (T this) {
  struct Trpath_T *trpath_chunk;

#if defined(TRPATHPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_trpath_memory(this);
#endif

  while (List_next(this->trpath_chunks) != NULL) {
    this->trpath_chunks = List_pop_keep(this->trpath_chunks,(void **) &trpath_chunk);
    FREE_KEEP(trpath_chunk);
  }
#ifdef TRPATHPOOL_REUSE
  int celli;
  trpath_chunk = (struct Trpath_T *) List_head(this->trpath_chunks);
  List_free_keep(&this->trpath_free_cells);
  this->trpath_free_cells = (List_T) NULL;
  for (celli = TRPATH_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->trpath_free_cells = List_push_keep(this->trpath_free_cells,(void *) &(trpath_chunk[celli]));
  }
#else
  this->trpath_cellptr = (struct Trpath_T *) List_head(this->trpath_chunks);
  this->trpath_cellctr = 0;
#endif

  return;
}

void
Trpathpool_free (T *old) {
  struct Trpath_T *trpath_chunk;

  while ((*old)->trpath_chunks != NULL) {
    (*old)->trpath_chunks = List_pop_keep((*old)->trpath_chunks,(void **) &trpath_chunk);
    FREE_KEEP(trpath_chunk);
  }
#ifdef TRPATHPOOL_REUSE
  List_free_keep(&(*old)->trpath_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Trpathpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->trpath_chunks = (List_T) NULL;
#ifdef TRPATHPOOL_REUSE
  new->trpath_free_cells = (List_T) NULL;
#else
  new->trpath_cellctr = 0;
#endif

  return new;
}


#ifdef TRPATHPOOL_REUSE
void
Trpathpool_free_trpath (Trpath_T *old, T this
#ifdef TRPATHPOOL_TRACE
			, const char *file, int line
#endif
			) {
  this->trpath_free_cells = List_push_keep(this->trpath_free_cells,(void *) *old);
#ifdef TRPATHPOOL_TRACE
    printf("Trpathpool: Freed %p -- Trpathpool_free_trpath called by %s:%d\n",*old,file,line);
#endif

  *old = (Trpath_T) NULL;
  return;
}  
#endif


static struct Trpath_T *
add_new_trpath_chunk (T this) {
  struct Trpath_T *chunk;

  chunk = (struct Trpath_T *) MALLOC_KEEP(TRPATH_CHUNKSIZE*sizeof(struct Trpath_T));
  this->trpath_chunks = List_push_keep(this->trpath_chunks,(void *) chunk);
#ifdef TRPATHPOOL_REUSE
  int celli;
  for (celli = TRPATH_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->trpath_free_cells = List_push_keep(this->trpath_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of trpath_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->trpath_chunks),chunk));

  return chunk;
}


Trpath_T
Trpathpool_new_trpath (T this
#ifdef TRPATHPOOL_TRACE
		       , const char *file, int line
#endif
		       ) {
  Trpath_T new;

#ifdef TRPATHPOOL_REUSE
  if (this->trpath_free_cells == (List_T) NULL) {
    add_new_trpath_chunk(this);
  }
  this->trpath_free_cells = List_pop_keep(this->trpath_free_cells,(void **) &new);
#else
  if (this->trpath_cellctr >= TRPATH_CHUNKSIZE) {
    this->trpath_cellptr = add_new_trpath_chunk(this);
    this->trpath_cellctr = 0;
  }
  this->trpath_cellctr += 1;
  new = this->trpath_cellptr++;
#endif

#ifdef TRPATHPOOL_TRACE
    printf("Trpathpool: Allocated %p -- Trpathpool_new_trpath called by %s:%d\n",new,file,line);
#endif

  return new;
}  


/* Guarantees that a chunk is available */
void
Trpathpool_init (T this) {
#ifdef TRPATHPOOL_REUSE
  add_new_trpath_chunk(this);
#else
  this->trpath_cellptr = add_new_trpath_chunk(this);
  this->trpath_cellctr = 0;
#endif

  return;
}

