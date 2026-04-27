static char rcsid[] = "$Id: f2a890d54dc96875af5b8d2a314946d9e545718e $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "univdiagpool.h"
#include "univdiagdef.h"

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


#define T Univdiagpool_T
struct T {
  List_T univdiag_chunks;
#ifdef UNIVDIAGPOOL_REUSE
  List_T univdiag_free_cells;
#else
  struct Univdiag_T *univdiag_cellptr;
  int univdiag_cellctr;
#endif

  List_T list_chunks;
#ifdef UNIVDIAGPOOL_REUSE
  List_T list_free_cells;
#else
  struct List_T *list_cellptr;
  int list_cellctr;
#endif
};


#if defined(UNIVDIAGPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_univdiag_chunk (T this, struct Univdiag_T *free_cell) {
  int chunki;
  List_T p;
  struct Univdiag_T *chunk;

  for (p = this->univdiag_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Univdiag_T *) List_head(p);
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
check_univdiag_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->univdiag_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->univdiag_free_cells) == CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->univdiag_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_univdiag_chunk(this,(struct Univdiag_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Univdiag_T cells leaked in univdiagpool chunk %d\n",
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
	fprintf(stderr,"%d out of %d List_T cells leaked in univdiagpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Univdiagpool_reset_memory (T this) {
  struct Univdiag_T *univdiag_chunk;
  struct List_T *list_chunk;

#if defined(UNIVDIAGPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_univdiag_memory(this);
  check_list_memory(this);
#endif

  while (List_next(this->univdiag_chunks) != NULL) {
    this->univdiag_chunks = List_pop_keep(this->univdiag_chunks,(void **) &univdiag_chunk);
    FREE_KEEP(univdiag_chunk);
  }
#ifdef UNIVDIAGPOOL_REUSE
  int celli;
  univdiag_chunk = (struct Univdiag_T *) List_head(this->univdiag_chunks);
  List_free_keep(&this->univdiag_free_cells);
  this->univdiag_free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->univdiag_free_cells = List_push_keep(this->univdiag_free_cells,(void *) &(univdiag_chunk[celli]));
  }
#else
  this->univdiag_cellptr = (struct Univdiag_T *) List_head(this->univdiag_chunks);
  this->univdiag_cellctr = 0;
#endif


  while (List_next(this->list_chunks) != NULL) {
    this->list_chunks = List_pop_keep(this->list_chunks,(void **) &list_chunk);
    FREE_KEEP(list_chunk);
  }
#ifdef UNIVDIAGPOOL_REUSE
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
Univdiagpool_free (T *old) {
  struct Univdiag_T *univdiag_chunk;
  struct List_T *list_chunk;

  while ((*old)->univdiag_chunks != NULL) {
    (*old)->univdiag_chunks = List_pop_keep((*old)->univdiag_chunks,(void **) &univdiag_chunk);
    FREE_KEEP(univdiag_chunk);
  }
#ifdef UNIVDIAGPOOL_REUSE
  List_free_keep(&(*old)->univdiag_free_cells);
#endif


  while ((*old)->list_chunks != NULL) {
    (*old)->list_chunks = List_pop_keep((*old)->list_chunks,(void **) &list_chunk);
    FREE_KEEP(list_chunk);
  }
#ifdef UNIVDIAGPOOL_REUSE
  List_free_keep(&(*old)->list_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Univdiagpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->univdiag_chunks = (List_T) NULL;
#ifdef UNIVDIAGPOOL_REUSE
  new->univdiag_free_cells = (List_T) NULL;
#else
  new->univdiag_cellctr = 0;
#endif

  new->list_chunks = (List_T) NULL;
#ifdef UNIVDIAGPOOL_REUSE
  new->list_free_cells = (List_T) NULL;
#else
  new->list_cellctr = 0;
#endif

  return new;
}


#ifdef UNIVDIAGPOOL_REUSE
void
Univdiagpool_free_univdiag (Univdiag_T *old, T this
#ifdef UNIVDIAGPOOL_TRACE
			    , const char *file, int line
#endif
			    ) {
  this->univdiag_free_cells = List_push_keep(this->univdiag_free_cells,(void *) *old);
#ifdef UNIVDIAGPOOL_TRACE
    printf("Univdiagpool/univdiag: Freed %p -- Univdiagpool_free_univdiag called by %s:%d\n",*old,file,line);
#endif

  *old = (Univdiag_T) NULL;
  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
void
Univdiagpool_free_list (List_T *list, T this
#ifdef UNIVDIAGPOOL_TRACE
			, const char *file, int line
#endif
			) {
  List_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) prev);
#ifdef UNIVDIAGPOOL_TRACE
    printf("Univdiagpool/list: Freed %p -- Univdiagpool_free_list called by %s:%d\n",prev,file,line);
#endif
  }

  *list = (List_T) NULL;

  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
void
Univdiagpool_gc (List_T *list, T this
#ifdef UNIVDIAGPOOL_TRACE
		 , const char *file, int line
#endif
		 ) {
  List_T prev;
  Univdiag_T univdiag;

  while ((prev = *list) != NULL) {
    univdiag = (Univdiag_T) prev->first;
    this->univdiag_free_cells = List_push_keep(this->univdiag_free_cells,(void *) univdiag);
#ifdef UNIVDIAGPOOL_TRACE
    printf("Univdiagpool/univdiag: Freed %p -- Univdiagpool_gc called by %s:%d\n",univdiag,file,line);
#endif

    *list = prev->rest;
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) prev);
#ifdef UNIVDIAGPOOL_TRACE
    printf("Univdiagpool/list: Freed %p -- Univdiagpool_gc called by %s:%d\n",prev,file,line);
#endif
  }

  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
List_T
Univdiagpool_pop (List_T list, T this, Univdiag_T *x
#ifdef UNIVDIAGPOOL_TRACE
		  , const char *file, int line
#endif
		  ) {
  List_T head;

  if (list != NULL) {
    this->list_free_cells = List_push_keep(this->list_free_cells,(void *) list);
#ifdef UNIVDIAGPOOL_TRACE
    printf("Univdiagpool/list: Freed %p -- Univdiagpool_pop called by %s:%d\n",list,file,line);
#endif
    head = list->rest;
    *x = (Univdiag_T) list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct Univdiag_T *
add_new_univdiag_chunk (T this) {
  struct Univdiag_T *chunk;

  chunk = (struct Univdiag_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct Univdiag_T));
  this->univdiag_chunks = List_push_keep(this->univdiag_chunks,(void *) chunk);
#ifdef UNIVDIAGPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->univdiag_free_cells = List_push_keep(this->univdiag_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of univdiags.  Ptr for univdiag %d is %p\n",
	       this->nunivdiags,chunk));

  return chunk;
}

static struct List_T *
add_new_list_chunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->list_chunks = List_push_keep(this->list_chunks,(void *) chunk);
#ifdef UNIVDIAGPOOL_REUSE
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
Univdiagpool_push (List_T list, T this, int qstart, int qend, int nmismatches, Univcoord_T univdiagonal
#ifdef UNIVDIAGPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  List_T listcell;
  Univdiag_T univdiag;

#ifdef UNIVDIAGPOOL_REUSE
  if (this->univdiag_free_cells == (List_T) NULL) {
    add_new_univdiag_chunk(this);
  }
  this->univdiag_free_cells = List_pop_keep(this->univdiag_free_cells,(void **) &univdiag);
#else
  if (this->univdiag_cellctr >= CHUNKSIZE) {
    this->univdiag_cellptr = add_new_univdiag_chunk(this);
    this->univdiag_cellctr = 0;
  }
  this->univdiag_cellctr += 1;
  univdiag = this->univdiag_cellptr++;
#endif

#ifdef UNIVDIAGPOOL_TRACE
  printf("Univdiagpool/univdiag: Allocated %p -- Univdiagpool_push called by %s:%d, contents %u %d..%d\n",
	 univdiag,file,line,univdiagonal,qstart,qend);
#endif
  univdiag->univdiagonal = univdiagonal;
  univdiag->qstart = qstart;
  univdiag->qend = qend;
  univdiag->nmismatches = nmismatches;

  assert(qstart < qend);


#ifdef UNIVDIAGPOOL_REUSE
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

#ifdef UNIVDIAGPOOL_TRACE
  printf("Univdiagpool/list: Allocated %p -- Univdiagpool_push called by %s:%d\n",listcell,file,line);
#endif

  listcell->first = (void *) univdiag;
  listcell->rest = list;

  return listcell;
}


List_T
Univdiagpool_push_existing (List_T list, T this, Univdiag_T univdiag
#ifdef UNIVDIAGPOOL_TRACE
			    , const char *file, int line
#endif
			    ) {
  List_T listcell;

#ifdef UNIVDIAGPOOL_REUSE
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

#ifdef UNIVDIAGPOOL_TRACE
  printf("Univdiagpool/list: Allocated %p -- Univdiagpool_push_existing called by %s:%d, contents %u %d..%d\n",
	 listcell,file,line,univdiag->univdiagonal,univdiag->qstart,univdiag->qend);
#endif

  listcell->first = (void *) univdiag;
  listcell->rest = list;

  return listcell;
}


Univdiag_T
Univdiagpool_new_univdiag (T this, int qstart, int qend, int nmismatches, Univcoord_T univdiagonal
#ifdef UNIVDIAGPOOL_TRACE
			   , const char *file, int line
#endif	
		   ) {
  Univdiag_T univdiag;

#ifdef UNIVDIAGPOOL_REUSE
  if (this->univdiag_free_cells == (List_T) NULL) {
    add_new_univdiag_chunk(this);
  }
  this->univdiag_free_cells = List_pop_keep(this->univdiag_free_cells,(void **) &univdiag);
#else
  if (this->univdiag_cellctr >= CHUNKSIZE) {
    this->univdiag_cellptr = add_new_univdiag_chunk(this);
    this->univdiag_cellctr = 0;
  }
  this->univdiag_cellctr += 1;
  univdiag = this->univdiag_cellptr++;
#endif

#ifdef UNIVDIAGPOOL_TRACE
  printf("Univdiagpool/univdiag: Allocated %p -- Univdiagpool_new_univdiag called by %s:%d, contents %u %d..%d\n",
	 univdiag,file,line,univdiagonal,qstart,qend);
#endif
  univdiag->univdiagonal = univdiagonal;
  univdiag->qstart = qstart;
  univdiag->qend = qend;
  univdiag->nmismatches = nmismatches;

  assert(qstart < qend);

  return univdiag;
}


/* Guarantees that a chunk is available */
void
Univdiagpool_init (T this) {
#ifdef UNIVDIAGPOOL_REUSE
  add_new_univdiag_chunk(this);
  add_new_list_chunk(this);
#else
  this->univdiag_cellptr = add_new_univdiag_chunk(this);
  this->univdiag_cellctr = 0;

  this->list_cellptr = add_new_list_chunk(this);
  this->list_cellctr = 0;
#endif

  return;
}


  
