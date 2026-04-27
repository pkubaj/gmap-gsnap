static char rcsid[] = "$Id: cdc2135617406448003a1647023edeacf888c5f4 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "auxinfopool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define AUXINFO_CHUNKSIZE 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Auxinfopool_T
struct T {
  List_T auxinfo_chunks;
#ifdef AUXINFOPOOL_REUSE
  List_T auxinfo_free_cells;
#else
  struct Auxinfo_T *auxinfo_cellptr;
  int auxinfo_cellctr;
#endif
};


#if defined(AUXINFOPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_auxinfo_chunk (T this, struct Auxinfo_T *free_cell) {
  int chunki;
  List_T p;
  struct Auxinfo_T *chunk;

  for (p = this->auxinfo_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Auxinfo_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[AUXINFO_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}


/* Checks to see if all memory has been returned to free_cells */
static void
check_auxinfo_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->auxinfo_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->auxinfo_free_cells) == AUXINFO_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->auxinfo_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_auxinfo_chunk(this,(struct Auxinfo_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < AUXINFO_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Auxinfo_T cells leaked in auxinfopool chunk %d\n",
		AUXINFO_CHUNKSIZE - nfreecells[chunki],AUXINFO_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Auxinfopool_reset_memory (T this) {
  struct Auxinfo_T *auxinfo_chunk;

#if defined(AUXINFOPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_auxinfo_memory(this);
#endif

  while (List_next(this->auxinfo_chunks) != NULL) {
    this->auxinfo_chunks = List_pop_keep(this->auxinfo_chunks,(void **) &auxinfo_chunk);
    FREE_KEEP(auxinfo_chunk);
  }
#ifdef AUXINFOPOOL_REUSE
  int celli;
  auxinfo_chunk = (struct Auxinfo_T *) List_head(this->auxinfo_chunks);
  List_free_keep(&this->auxinfo_free_cells);
  this->auxinfo_free_cells = (List_T) NULL;
  for (celli = AUXINFO_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->auxinfo_free_cells = List_push_keep(this->auxinfo_free_cells,(void *) &(auxinfo_chunk[celli]));
  }
#else
  this->auxinfo_cellptr = (struct Auxinfo_T *) List_head(this->auxinfo_chunks);
  this->auxinfo_cellctr = 0;
#endif

  return;
}

void
Auxinfopool_free (T *old) {
  struct Auxinfo_T *auxinfo_chunk;

  while ((*old)->auxinfo_chunks != NULL) {
    (*old)->auxinfo_chunks = List_pop_keep((*old)->auxinfo_chunks,(void **) &auxinfo_chunk);
    FREE_KEEP(auxinfo_chunk);
  }
#ifdef AUXINFOPOOL_REUSE
  List_free_keep(&(*old)->auxinfo_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Auxinfopool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->auxinfo_chunks = (List_T) NULL;
#ifdef AUXINFOPOOL_REUSE
  new->auxinfo_free_cells = (List_T) NULL;
#else
  new->auxinfo_cellctr = 0;
#endif

  return new;
}


#ifdef AUXINFOPOOL_REUSE
void
Auxinfopool_free_auxinfo (Auxinfo_T *old, T this
#ifdef AUXINFOPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  this->auxinfo_free_cells = List_push_keep(this->auxinfo_free_cells,(void *) *old);
#ifdef AUXINFOPOOL_TRACE
  printf("Auxinfopool/auxinfo: Freed %p -- Auxinfopool_free_auxinfo called by %s:%d\n",*old,file,line);
#endif

  *old = (Auxinfo_T) NULL;

  return;
}  
#endif


static struct Auxinfo_T *
add_new_auxinfo_chunk (T this) {
  struct Auxinfo_T *chunk;

  chunk = (struct Auxinfo_T *) MALLOC_KEEP(AUXINFO_CHUNKSIZE*sizeof(struct Auxinfo_T));
  this->auxinfo_chunks = List_push_keep(this->auxinfo_chunks,(void *) chunk);
#ifdef AUXINFOPOOL_REUSE
  int celli;
  for (celli = AUXINFO_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->auxinfo_free_cells = List_push_keep(this->auxinfo_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of auxinfo_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->auxinfo_chunks),chunk));

  return chunk;
}

Auxinfo_T
Auxinfopool_new_auxinfo (T this
#ifdef AUXINFOPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Auxinfo_T new;

#ifdef AUXINFOPOOL_REUSE
  if (this->auxinfo_free_cells == (List_T) NULL) {
    add_new_auxinfo_chunk(this);
  }
  this->auxinfo_free_cells = List_pop_keep(this->auxinfo_free_cells,(void **) &new);
#else
  if (this->auxinfo_cellctr >= AUXINFO_CHUNKSIZE) {
    this->auxinfo_cellptr = add_new_auxinfo_chunk(this);
    this->auxinfo_cellctr = 0;
  }
  this->auxinfo_cellctr += 1;
  new = this->auxinfo_cellptr++;
#endif

#ifdef AUXINFOPOOL_TRACE
  printf("Auxinfopool/auxinfo: Allocated %p -- Auxinfopool_new_auxinfo called by %s:%d\n",new,file,line);
#endif

  return new;
}  


/* Guarantees that a chunk is available */
void
Auxinfopool_init (T this) {
#ifdef AUXINFOPOOL_REUSE
  add_new_auxinfo_chunk(this);
#else
  this->auxinfo_cellptr = add_new_auxinfo_chunk(this);
  this->auxinfo_cellctr = 0;
#endif

  return;
}

