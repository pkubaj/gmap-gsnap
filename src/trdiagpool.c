static char rcsid[] = "$Id: 7007f49665a2261fe780dc2387607d38999cb72a $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "trdiagpool.h"
#include "trdiagdef.h"

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "comp.h"


#define CHUNKSIZE 20000

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Trdiagpool_T
struct T {
  List_T trdiag_chunks;
#ifdef TRDIAGPOOL_REUSE
  List_T trdiag_free_cells;
#else
  struct Trdiag_T *trdiag_cellptr;
  int trdiag_cellctr;
#endif

  List_T list_chunks;
#ifdef TRDIAGPOOL_REUSE
  List_T list_free_cells;
#else
  struct List_T *list_cellptr;
  int list_cellctr;
#endif
};


#if defined(TRDIAGPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_trdiag_chunk (T this, struct Trdiag_T *free_cell) {
  int chunki;
  List_T p;
  struct Trdiag_T *chunk;

  for (p = this->trdiag_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Trdiag_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}

static int
find_list_chunk (T this, struct List_T *free_cell) {
  int chunki;
  List_T p;
  struct List_T *chunk;

  for (p = this->list_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct List_T *) List_head(p);
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
check_trdiag_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->trdiag_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->trdiag_free_cells) == CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->trdiag_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_trdiag_chunk(this,(struct Trdiag_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Trdiag_T cells leaked in trdiagpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}

static void
check_list_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->list_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->list_free_cells) == CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->list_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_list_chunk(this,(struct List_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d List_T cells leaked in trdiagpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Trdiagpool_reset_memory (T this) {
  struct Trdiag_T *trdiag_chunk;
  struct List_T *list_chunk;

#if defined(TRDIAGPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_trdiag_memory(this);
  check_list_memory(this);
#endif

  while (List_next(this->trdiag_chunks) != NULL) {
    this->trdiag_chunks = List_pop_keep(this->trdiag_chunks,(void **) &trdiag_chunk);
    FREE_KEEP(trdiag_chunk);
  }
#ifdef TRDIAGPOOL_REUSE
  int celli;
  trdiag_chunk = (struct Trdiag_T *) List_head(this->trdiag_chunks);
  List_free_keep(&this->trdiag_free_cells);
  this->trdiag_free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->trdiag_free_cells = List_push_keep(this->trdiag_free_cells,(void *) &(trdiag_chunk[celli]));
  }
#else
  this->trdiag_cellptr = (struct Trdiag_T *) List_head(this->trdiag_chunks);
  this->trdiag_cellctr = 0;
#endif


  while (List_next(this->list_chunks) != NULL) {
    this->list_chunks = List_pop_keep(this->list_chunks,(void **) &list_chunk);
    FREE_KEEP(list_chunk);
  }
#ifdef TRDIAGPOOL_REUSE
  list_chunk = (struct List_T *) List_head(this->list_chunks);
  List_free_keep(&this->list_free_cells);
  this->list_free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) &(list_chunk[celli]));
  }
#else
  this->list_cellptr = (struct List_T *) List_head(this->list_chunks);
  this->list_cellctr = 0;
#endif

  return;
}


void
Trdiagpool_free (T *old) {
  struct Trdiag_T *trdiag_chunk;
  struct List_T *list_chunk;

  while ((*old)->trdiag_chunks != NULL) {
    (*old)->trdiag_chunks = List_pop_keep((*old)->trdiag_chunks,(void **) &trdiag_chunk);
    FREE_KEEP(trdiag_chunk);
  }
#ifdef TRDIAGPOOL_REUSE
  List_free_keep(&(*old)->trdiag_free_cells);
#endif


  while ((*old)->list_chunks != NULL) {
    (*old)->list_chunks = List_pop_keep((*old)->list_chunks,(void **) &list_chunk);
    FREE_KEEP(list_chunk);
  }
#ifdef TRDIAGPOOL_REUSE
  List_free_keep(&(*old)->list_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Trdiagpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->trdiag_chunks = (List_T) NULL;
#ifdef TRDIAGPOOL_REUSE
  new->trdiag_free_cells = (List_T) NULL;
#else
  new->trdiag_cellctr = 0;
#endif

  new->list_chunks = (List_T) NULL;
#ifdef TRDIAGPOOL_REUSE
  new->list_free_cells = (List_T) NULL;
#else
  new->list_cellctr = 0;
#endif

  return new;
}


#ifdef TRDIAGPOOL_REUSE
void
Trdiagpool_free_trdiag (Trdiag_T *old, T this
#ifdef TRDIAGPOOL_TRACE
			, const char *file, int line
#endif
			) {
  this->trdiag_free_cells = List_push_keep(this->trdiag_free_cells,(void *) *old);
#ifdef TRDIAGPOOL_TRACE
    printf("Trdiagpool/trdiag: Freed %p -- Trdiagpool_free_trdiag called by %s:%d\n",*old,file,line);
#endif

  *old = (Trdiag_T) NULL;
  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
void
Trdiagpool_free_list (List_T *list, T this
#ifdef TRDIAGPOOL_TRACE
		      , const char *file, int line
#endif
		      ) {
  List_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) prev);
#ifdef TRDIAGPOOL_TRACE
    printf("Trdiagpool/list: Freed %p -- Trdiagpool_free_list called by %s:%d\n",prev,file,line);
#endif
  }

  *list = (List_T) NULL;

  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
void
Trdiagpool_gc (List_T *list, T this
#ifdef TRDIAGPOOL_TRACE
	       , const char *file, int line
#endif
	       ) {
  List_T prev;
  Trdiag_T trdiag;

  while ((prev = *list) != NULL) {
    trdiag = (Trdiag_T) prev->first;
    this->trdiag_free_cells = List_push_keep(this->trdiag_free_cells,(void *) trdiag);
#ifdef TRDIAGPOOL_TRACE
    printf("Trdiagpool/trdiag: Freed %p -- Trdiagpool_gc called by %s:%d\n",trdiag,file,line);
#endif

    *list = prev->rest;
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) prev);
#ifdef TRDIAGPOOL_TRACE
    printf("Trdiagpool/list: Freed %p -- Trdiagpool_gc called by %s:%d\n",prev,file,line);
#endif
  }

  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
List_T
Trdiagpool_pop (List_T list, T this, Trdiag_T *x
#ifdef TRDIAGPOOL_TRACE
		, const char *file, int line
#endif
		) {
  List_T head;

  if (list != NULL) {
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) list);
#ifdef TRDIAGPOOL_TRACE
    printf("Trdiagpool/list: Freed %p -- Trdiagpool_pop called by %s:%d\n",list,file,line);
#endif
    head = list->rest;
    *x = (Trdiag_T) list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct Trdiag_T *
add_new_trdiag_chunk (T this) {
  struct Trdiag_T *chunk;

  chunk = (struct Trdiag_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Trdiag_T));
  this->trdiag_chunks = List_push_keep(this->trdiag_chunks,(void *) chunk);
#ifdef TRDIAGPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->trdiag_free_cells = List_push_keep(this->trdiag_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of trdiags.  Ptr for trdiag %d is %p\n",
	       this->ntrdiags,chunk));

  return chunk;
}

static struct List_T *
add_new_list_chunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->list_chunks = List_push_keep(this->list_chunks,(void *) chunk);
#ifdef TRDIAGPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of listcells.  Ptr for listcell %d is %p\n",
	       this->nlistcells,chunk));

  return chunk;
}


List_T
Trdiagpool_push (List_T list, T this, int qstart, int qend, int nmismatches, Trcoord_T trdiagonal
#ifdef TRDIAGPOOL_TRACE
		 , const char *file, int line
#endif
		 ) {
  List_T listcell;
  Trdiag_T trdiag;

#ifdef TRDIAGPOOL_REUSE
  if (this->trdiag_free_cells == (List_T) NULL) {
    add_new_trdiag_chunk(this);
  }
  this->trdiag_free_cells = List_pop_keep(this->trdiag_free_cells,(void **) &trdiag);
#else
  if (this->trdiag_cellctr >= CHUNKSIZE) {
    this->trdiag_cellptr = add_new_trdiag_chunk(this);
    this->trdiag_cellctr = 0;
  }
  this->trdiag_cellctr += 1;
  trdiag = this->trdiag_cellptr++;
#endif

#ifdef TRDIAGPOOL_TRACE
  printf("Trdiagpool/trdiag: Allocated %p -- Trdiagpool_push called by %s:%d, contents %u %d..%d\n",
	 trdiag,file,line,trdiagonal,qstart,qend);
#endif
  trdiag->trdiagonal = trdiagonal;
  trdiag->qstart = qstart;
  trdiag->qend = qend;
  trdiag->nmismatches = nmismatches;

  assert(qstart < qend);


#ifdef TRDIAGPOOL_REUSE
  if (this->list_free_cells == (List_T) NULL) {
    add_new_list_chunk(this);
  }
  this->list_free_cells = List_pop_keep(this->list_free_cells,(void **) &listcell);
#else
  if (this->list_cellctr >= CHUNKSIZE) {
    this->list_cellptr = add_new_list_chunk(this);
    this->list_cellctr = 0;
  }
  this->list_cellctr += 1;
  listcell = this->list_cellptr++;
#endif

#ifdef TRDIAGPOOL_TRACE
  printf("Trdiagpool/list: Allocated %p -- Trdiagpool_push called by %s:%d\n",listcell,file,line);
#endif

  listcell->first = (void *) trdiag;
  listcell->rest = list;

  return listcell;
}


List_T
Trdiagpool_push_existing (List_T list, T this, Trdiag_T trdiag
#ifdef TRDIAGPOOL_TRACE
			  , const char *file, int line
#endif
			  ) {
  List_T listcell;

#ifdef TRDIAGPOOL_REUSE
  if (this->list_free_cells == (List_T) NULL) {
    add_new_list_chunk(this);
  }
  this->list_free_cells = List_pop_keep(this->list_free_cells,(void **) &listcell);
#else
  if (this->list_cellctr >= CHUNKSIZE) {
    this->list_cellptr = add_new_list_chunk(this);
    this->list_cellctr = 0;
  }
  this->list_cellctr += 1;
  listcell = this->list_cellptr++;
#endif

#ifdef TRDIAGPOOL_TRACE
  printf("Trdiagpool/list: Allocated %p -- Trdiagpool_push_existing called by %s:%d, contents %u %d..%d\n",
	 listcell,file,line,trdiag->trdiagonal,trdiag->qstart,trdiag->qend);
#endif

  listcell->first = (void *) trdiag;
  listcell->rest = list;

  return listcell;
}


Trdiag_T
Trdiagpool_new_trdiag (T this, int qstart, int qend, int nmismatches, Trcoord_T trdiagonal
#ifdef TRDIAGPOOL_TRACE
		       , const char *file, int line
#endif	
		       ) {
  Trdiag_T trdiag;

#ifdef TRDIAGPOOL_REUSE
  if (this->trdiag_free_cells == (List_T) NULL) {
    add_new_trdiag_chunk(this);
  }
  this->trdiag_free_cells = List_pop_keep(this->trdiag_free_cells,(void **) &trdiag);
#else
  if (this->trdiag_cellctr >= CHUNKSIZE) {
    this->trdiag_cellptr = add_new_trdiag_chunk(this);
    this->trdiag_cellctr = 0;
  }
  this->trdiag_cellctr += 1;
  trdiag = this->trdiag_cellptr++;
#endif

#ifdef TRDIAGPOOL_TRACE
  printf("Trdiagpool/trdiag: Allocated %p -- Trdiagpool_new_trdiag called by %s:%d, contents %u %d..%d\n",
	 trdiag,file,line,trdiagonal,qstart,qend);
#endif
  trdiag->trdiagonal = trdiagonal;
  trdiag->qstart = qstart;
  trdiag->qend = qend;
  trdiag->nmismatches = nmismatches;

  assert(qstart < qend);

  return trdiag;
}


/* Guarantees that a chunk is available */
void
Trdiagpool_init (T this) {
#ifdef TRDIAGPOOL_REUSE
  add_new_trdiag_chunk(this);
  add_new_list_chunk(this);
#else
  this->trdiag_cellptr = add_new_trdiag_chunk(this);
  this->trdiag_cellctr = 0;

  this->list_cellptr = add_new_list_chunk(this);
  this->list_cellctr = 0;
#endif

  return;
}


  
