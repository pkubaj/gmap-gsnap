static char rcsid[] = "$Id: 7cf772cd3c48820771555293cca4dc45a0db1efc $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pathstorepool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define PATHSTORE_CHUNKSIZE 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Pathstorepool_T
struct T {
  List_T pathstore_chunks;
#ifdef PATHSTOREPOOL_REUSE
  List_T pathstore_free_cells;
#else
  struct Pathstore_T *pathstore_cellptr;
  int pathstore_cellctr;
#endif
};


#if defined(PATHSTOREPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_pathstore_chunk (T this, struct Pathstore_T *free_cell) {
  int chunki;
  List_T p;
  struct Pathstore_T *chunk;

  for (p = this->pathstore_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Pathstore_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[PATHSTORE_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}


/* Checks to see if all memory has been returned to free_cells */
static void
check_pathstore_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->pathstore_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->pathstore_free_cells) == PATHSTORE_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->pathstore_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_pathstore_chunk(this,(struct Pathstore_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < PATHSTORE_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Pathstore_T cells leaked in pathstorepool chunk %d\n",
		PATHSTORE_CHUNKSIZE - nfreecells[chunki],PATHSTORE_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Pathstorepool_reset_memory (T this) {
  struct Pathstore_T *pathstore_chunk;

#if defined(PATHSTOREPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_pathstore_memory(this);
#endif

  while (List_next(this->pathstore_chunks) != NULL) {
    this->pathstore_chunks = List_pop_keep(this->pathstore_chunks,(void **) &pathstore_chunk);
    FREE_KEEP(pathstore_chunk);
  }
#ifdef PATHSTOREPOOL_REUSE
  int celli;
  pathstore_chunk = (struct Pathstore_T *) List_head(this->pathstore_chunks);
  List_free_keep(&this->pathstore_free_cells);
  this->pathstore_free_cells = (List_T) NULL;
  for (celli = PATHSTORE_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->pathstore_free_cells = List_push_keep(this->pathstore_free_cells,(void *) &(pathstore_chunk[celli]));
  }
#else
  this->pathstore_cellptr = (struct Pathstore_T *) List_head(this->pathstore_chunks);
  this->pathstore_cellctr = 0;
#endif

  return;
}

void
Pathstorepool_free (T *old) {
  struct Pathstore_T *pathstore_chunk;

  while ((*old)->pathstore_chunks != NULL) {
    (*old)->pathstore_chunks = List_pop_keep((*old)->pathstore_chunks,(void **) &pathstore_chunk);
    FREE_KEEP(pathstore_chunk);
  }
#ifdef PATHSTOREPOOL_REUSE
  List_free_keep(&(*old)->pathstore_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Pathstorepool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->pathstore_chunks = (List_T) NULL;
#ifdef PATHSTOREPOOL_REUSE
  new->pathstore_free_cells = (List_T) NULL;
#else
  new->pathstore_cellctr = 0;
#endif

  return new;
}


#ifdef PATHSTOREPOOL_REUSE
void
Pathstorepool_free_pathstore (Pathstore_T *old, T this
#ifdef PATHSTOREPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  this->pathstore_free_cells = List_push_keep(this->pathstore_free_cells,(void *) *old);
#ifdef PATHSTOREPOOL_TRACE
  printf("Pathstorepool/pathstore: Freed %p -- Pathstorepool_free_pathstore called by %s:%d\n",*old,file,line);
#endif

  *old = (Pathstore_T) NULL;

  return;
}  
#endif


static struct Pathstore_T *
add_new_pathstore_chunk (T this) {
  struct Pathstore_T *chunk;

  chunk = (struct Pathstore_T *) MALLOC_KEEP(PATHSTORE_CHUNKSIZE*sizeof(struct Pathstore_T));
  this->pathstore_chunks = List_push_keep(this->pathstore_chunks,(void *) chunk);
#ifdef PATHSTOREPOOL_REUSE
  int celli;
  for (celli = PATHSTORE_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->pathstore_free_cells = List_push_keep(this->pathstore_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of pathstore_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->pathstore_chunks),chunk));

  return chunk;
}

Pathstore_T
Pathstorepool_new_pathstore (T this
#ifdef PATHSTOREPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Pathstore_T new;

#ifdef PATHSTOREPOOL_REUSE
  if (this->pathstore_free_cells == (List_T) NULL) {
    add_new_pathstore_chunk(this);
  }
  this->pathstore_free_cells = List_pop_keep(this->pathstore_free_cells,(void **) &new);
#else
  if (this->pathstore_cellctr >= PATHSTORE_CHUNKSIZE) {
    this->pathstore_cellptr = add_new_pathstore_chunk(this);
    this->pathstore_cellctr = 0;
  }
  this->pathstore_cellctr += 1;
  new = this->pathstore_cellptr++;
#endif

#ifdef PATHSTOREPOOL_TRACE
  printf("Pathstorepool/pathstore: Allocated %p -- Pathstorepool_new_pathstore called by %s:%d\n",new,file,line);
#endif

  return new;
}  


/* Guarantees that a chunk is available */
void
Pathstorepool_init (T this) {
#ifdef PATHSTOREPOOL_REUSE
  add_new_pathstore_chunk(this);
#else
  this->pathstore_cellptr = add_new_pathstore_chunk(this);
  this->pathstore_cellctr = 0;
#endif

  return;
}

