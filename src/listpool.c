static char rcsid[] = "$Id: b37459c11cd0c6b2abad313d4c882c4dd67956be $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "listpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"


/* Used for lists of Substring_T and Junction_T objects in Stage3end_T objects */
/* Used for lists of Path_T, Uinttable_T, and Univdiag_T objects in path-solve.c */
/* Used for lists of Elt_T objects in extension-search.c */


#define CHUNKSIZE 16384

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Listpool_T
struct T {
  List_T chunks;
#ifdef LISTPOOL_REUSE
  List_T free_cells;
#else
  struct List_T *cellptr;
  int cellctr;
#endif
};


#if defined(LISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_chunk (T this, struct List_T *free_cell) {
  int chunki;
  List_T p;
  struct List_T *chunk;

  for (p = this->chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
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
      chunki = find_chunk(this,(struct List_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < CHUNKSIZE) {
	fprintf(stderr,"%d out of %d List_T cells leaked in listpool chunk %d\n",
		CHUNKSIZE - nfreecells[chunki],CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Listpool_reset_memory (T this) {
  struct List_T *chunk;

#if defined(LISTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_memory(this);
#endif

  while (List_next(this->chunks) != NULL) {
    this->chunks = List_pop_keep(this->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef LISTPOOL_REUSE
  int celli;

  chunk = (struct List_T *) List_head(this->chunks);
  List_free_keep(&this->free_cells);
  this->free_cells = (List_T) NULL;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#else
  this->cellptr = (struct List_T *) List_head(this->chunks);
  this->cellctr = 0;
#endif

  return;
}

void
Listpool_free (T *old) {
  struct List_T *chunk;

  while ((*old)->chunks != NULL) {
    (*old)->chunks = List_pop_keep((*old)->chunks,(void **) &chunk);
    FREE_KEEP(chunk);
  }
#ifdef LISTPOOL_REUSE
  List_free_keep(&(*old)->free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Listpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->chunks = (List_T) NULL;
#ifdef LISTPOOL_REUSE
  new->free_cells = (List_T) NULL;
#else
  new->cellctr = 0;
#endif

  return new;
}


#ifdef LISTPOOL_REUSE
void
Listpool_free_list (List_T *list, T this
#ifdef LISTPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  List_T prev;

  while ((prev = *list) != NULL) {
    *list = prev->rest;
    this->free_cells = List_push_keep(this->free_cells,(void *) prev);
#ifdef LISTPOOL_TRACE
    printf("Listpool: Freed %p -- Listpool_free_list called by %s:%d\n",prev,file,line);
#endif
  }

  *list = (List_T) NULL;

  return;
}
#endif


#ifdef LISTPOOL_REUSE
List_T
Listpool_pop (List_T list, T this, void **contents
#ifdef LISTPOOL_TRACE
	      , const char *file, int line
#endif
	      ) {
  List_T head;

  if (list != NULL) {
    this->free_cells = List_push_keep(this->free_cells,(void *) list);
#ifdef LISTPOOL_TRACE
    printf("Listpool: Freed %p -- Listpool_pop called by %s:%d\n",list,file,line);
#endif
    head = list->rest;
    *contents = list->first;
    return head;
  } else {
    return list;
  }
}
#endif


static struct List_T *
add_new_chunk (T this) {
  struct List_T *chunk;

  chunk = (struct List_T *) MALLOC_KEEP(CHUNKSIZE*sizeof(struct List_T));
  this->chunks = List_push_keep(this->chunks,(void *) chunk);
#ifdef LISTPOOL_REUSE
  int celli;
  for (celli = CHUNKSIZE - 1; celli >= 0; celli--) {
    this->free_cells = List_push_keep(this->free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of listcells.  Ptr for chunk %d is %p\n",
	       List_length(this->chunks),chunk));

  return chunk;
}


List_T
Listpool_push (List_T list, T this, void *contents
#ifdef LISTPOOL_TRACE
	       , const char *file, int line
#endif
	       ) {
  List_T listcell;

#ifdef LISTPOOL_REUSE
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

#ifdef LISTPOOL_TRACE
  printf("Listpool: Allocated %p -- Listpool_push called by %s:%d\n",listcell,file,line);
#endif

  listcell->first = contents;
  listcell->rest = list;

  return listcell;
}

List_T
Listpool_copy (List_T source, T this) {
  List_T dest = NULL;

  while (source != NULL) {
    dest = Listpool_push(dest,this,/*orig*/source->first
			 listpool_trace(__FILE__,__LINE__));
    source = source->rest;
  }
  return List_reverse(dest);
}


/* Guarantees that a chunk is available */
void
Listpool_init (T this) {
#ifdef LISTPOOL_REUSE
  add_new_chunk(this);
#else
  this->cellptr = add_new_chunk(this);;
  this->cellctr = 0;
#endif

  return;
}

