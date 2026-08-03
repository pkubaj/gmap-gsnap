#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "segmenteltpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define SEGMENTELT_CHUNKSIZE 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Segmenteltpool_T
struct T {
  List_T segmentelt_chunks;
#ifdef SEGMENTELTPOOL_REUSE
  List_T segmentelt_free_cells;
#else
  struct Segmentelt_T *segmentelt_cellptr;
  int segmentelt_cellctr;
#endif
};


#if defined(SEGMENTELTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_segmentelt_chunk (T this, struct Segmentelt_T *free_cell) {
  int chunki;
  List_T p;
  struct Segmentelt_T *chunk;

  for (p = this->segmentelt_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Segmentelt_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[SEGMENTELT_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}


/* Checks to see if all memory has been returned to free_cells */
static void
check_segmentelt_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->segmentelt_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->segmentelt_free_cells) == SEGMENTELT_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->segmentelt_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_segmentelt_chunk(this,(struct Segmentelt_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < SEGMENTELT_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Segmentelt_T cells leaked in segmenteltpool chunk %d\n",
		SEGMENTELT_CHUNKSIZE - nfreecells[chunki],SEGMENTELT_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Segmenteltpool_reset_memory (T this) {
  struct Segmentelt_T *segmentelt_chunk;

#if defined(SEGMENTELTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_segmentelt_memory(this);
#endif

  while (List_next(this->segmentelt_chunks) != NULL) {
    this->segmentelt_chunks = List_pop_keep(this->segmentelt_chunks,(void **) &segmentelt_chunk);
    FREE_KEEP(segmentelt_chunk);
  }
#ifdef SEGMENTELTPOOL_REUSE
  int celli;
  segmentelt_chunk = (struct Segmentelt_T *) List_head(this->segmentelt_chunks);
  List_free_keep(&this->segmentelt_free_cells);
  this->segmentelt_free_cells = (List_T) NULL;
  for (celli = SEGMENTELT_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->segmentelt_free_cells = List_push_keep(this->segmentelt_free_cells,(void *) &(segmentelt_chunk[celli]));
  }
#else
  this->segmentelt_cellptr = (struct Segmentelt_T *) List_head(this->segmentelt_chunks);
  this->segmentelt_cellctr = 0;
#endif

  return;
}

void
Segmenteltpool_free (T *old) {
  struct Segmentelt_T *segmentelt_chunk;

  while ((*old)->segmentelt_chunks != NULL) {
    (*old)->segmentelt_chunks = List_pop_keep((*old)->segmentelt_chunks,(void **) &segmentelt_chunk);
    FREE_KEEP(segmentelt_chunk);
  }
#ifdef SEGMENTELTPOOL_REUSE
  List_free_keep(&(*old)->segmentelt_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Segmenteltpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->segmentelt_chunks = (List_T) NULL;
#ifdef SEGMENTELTPOOL_REUSE
  new->segmentelt_free_cells = (List_T) NULL;
#else
  new->segmentelt_cellctr = 0;
#endif

  return new;
}


#ifdef SEGMENTELTPOOL_REUSE
void
Segmenteltpool_free_segmentelt (Segmentelt_T *old, T this
#ifdef SEGMENTELTPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  this->segmentelt_free_cells = List_push_keep(this->segmentelt_free_cells,(void *) *old);
#ifdef SEGMENTELTPOOL_TRACE
  printf("Segmenteltpool/segmentelt: Freed %p -- Segmenteltpool_free_segmentelt called by %s:%d\n",*old,file,line);
#endif

  *old = (Segmentelt_T) NULL;

  return;
}  
#endif


static struct Segmentelt_T *
add_new_segmentelt_chunk (T this) {
  struct Segmentelt_T *chunk;

  chunk = (struct Segmentelt_T *) MALLOC_KEEP(SEGMENTELT_CHUNKSIZE*sizeof(struct Segmentelt_T));
  this->segmentelt_chunks = List_push_keep(this->segmentelt_chunks,(void *) chunk);
#ifdef SEGMENTELTPOOL_REUSE
  int celli;
  for (celli = SEGMENTELT_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->segmentelt_free_cells = List_push_keep(this->segmentelt_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of segmentelt_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->segmentelt_chunks),chunk));

  return chunk;
}

Segmentelt_T
Segmenteltpool_new_segmentelt (T this
#ifdef SEGMENTELTPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Segmentelt_T new;

#ifdef SEGMENTELTPOOL_REUSE
  if (this->segmentelt_free_cells == (List_T) NULL) {
    add_new_segmentelt_chunk(this);
  }
  this->segmentelt_free_cells = List_pop_keep(this->segmentelt_free_cells,(void **) &new);
#else
  if (this->segmentelt_cellctr >= SEGMENTELT_CHUNKSIZE) {
    this->segmentelt_cellptr = add_new_segmentelt_chunk(this);
    this->segmentelt_cellctr = 0;
  }
  this->segmentelt_cellctr += 1;
  new = this->segmentelt_cellptr++;
#endif

#ifdef SEGMENTELTPOOL_TRACE
  printf("Segmenteltpool/segmentelt: Allocated %p -- Segmenteltpool_new_segmentelt called by %s:%d\n",new,file,line);
#endif

  return new;
}  


/* Guarantees that a chunk is available */
void
Segmenteltpool_init (T this) {
#ifdef SEGMENTELTPOOL_REUSE
  add_new_segmentelt_chunk(this);
#else
  this->segmentelt_cellptr = add_new_segmentelt_chunk(this);
  this->segmentelt_cellctr = 0;
#endif

  return;
}

