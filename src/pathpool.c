static char rcsid[] = "$Id: cdc2135617406448003a1647023edeacf888c5f4 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pathpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define PATH_CHUNKSIZE 1024
#define JUNCTION_CHUNKSIZE 2048
#define ALTSPLICE_CHUNKSIZE 512
#define STRING_CHUNKSIZE 65536


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Pathpool_T
struct T {
  List_T path_chunks;
#ifdef PATHPOOL_REUSE
  List_T path_free_cells;
#else
  struct Path_T *path_cellptr;
  int path_cellctr;
#endif

  List_T junction_chunks;
#ifdef PATHPOOL_REUSE
  List_T junction_free_cells;
#else
  struct Junction_T *junction_cellptr;
  int junction_cellctr;
#endif

  List_T altsplice_chunks;
#ifdef PATHPOOL_REUSE
  List_T altsplice_free_cells;
#else
  struct Altsplice_T *altsplice_cellptr;
  int altsplice_cellctr;
#endif

  /* Handles variable sizes of strings, so cannot reuse */
  List_T string_chunks;
  char *string_cellptr;
  int string_cellctr;
  int string_chunksize;
};


#if defined(PATHPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_path_chunk (T this, struct Path_T *free_cell) {
  int chunki;
  List_T p;
  struct Path_T *chunk;

  for (p = this->path_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Path_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[PATH_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}

static int
find_junction_chunk (T this, struct Junction_T *free_cell) {
  int chunki;
  List_T p;
  struct Junction_T *chunk;

  for (p = this->junction_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Junction_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[JUNCTION_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}

static int
find_altsplice_chunk (T this, struct Altsplice_T *free_cell) {
  int chunki;
  List_T p;
  struct Altsplice_T *chunk;

  for (p = this->altsplice_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Altsplice_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[ALTSPLICE_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}


/* Checks to see if all memory has been returned to free_cells */
static void
check_path_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->path_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->path_free_cells) == PATH_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->path_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_path_chunk(this,(struct Path_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < PATH_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Path_T cells leaked in pathpool chunk %d\n",
		PATH_CHUNKSIZE - nfreecells[chunki],PATH_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}

static void
check_junction_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->junction_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->junction_free_cells) == JUNCTION_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->junction_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_junction_chunk(this,(struct Junction_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < JUNCTION_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Junction_T cells leaked in pathpool chunk %d\n",
		JUNCTION_CHUNKSIZE - nfreecells[chunki],JUNCTION_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}

static void
check_altsplice_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->altsplice_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->altsplice_free_cells) == ALTSPLICE_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->altsplice_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_altsplice_chunk(this,(struct Altsplice_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < ALTSPLICE_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Altsplice_T cells leaked in pathpool chunk %d\n",
		ALTSPLICE_CHUNKSIZE - nfreecells[chunki],ALTSPLICE_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Pathpool_reset_memory (T this) {
  struct Path_T *path_chunk;
  struct Junction_T *junction_chunk;
  struct Altsplice_T *altsplice_chunk;
  char *string_chunk;

#if defined(PATHPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_path_memory(this);
  check_junction_memory(this);
  check_altsplice_memory(this);
#endif

  while (List_next(this->path_chunks) != NULL) {
    this->path_chunks = List_pop_keep(this->path_chunks,(void **) &path_chunk);
    FREE_KEEP(path_chunk);
  }
#ifdef PATHPOOL_REUSE
  int celli;
  path_chunk = (struct Path_T *) List_head(this->path_chunks);
  List_free_keep(&this->path_free_cells);
  this->path_free_cells = (List_T) NULL;
  for (celli = PATH_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->path_free_cells = List_push_keep(this->path_free_cells,(void *) &(path_chunk[celli]));
  }
#else
  this->path_cellptr = (struct Path_T *) List_head(this->path_chunks);
  this->path_cellctr = 0;
#endif


  while (List_next(this->junction_chunks) != NULL) {
    this->junction_chunks = List_pop_keep(this->junction_chunks,(void **) &junction_chunk);
    FREE_KEEP(junction_chunk);
  }
#ifdef PATHPOOL_REUSE
  junction_chunk = (struct Junction_T *) List_head(this->junction_chunks);
  List_free_keep(&this->junction_free_cells);
  this->junction_free_cells = (List_T) NULL;
  for (celli = JUNCTION_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->junction_free_cells = List_push_keep(this->junction_free_cells,(void *) &(junction_chunk[celli]));
  }
#else
  this->junction_cellptr = (struct Junction_T *) List_head(this->junction_chunks);
  this->junction_cellctr = 0;
#endif  
  

  while (List_next(this->altsplice_chunks) != NULL) {
    this->altsplice_chunks = List_pop_keep(this->altsplice_chunks,(void **) &altsplice_chunk);
    FREE_KEEP(altsplice_chunk);
  }
#ifdef PATHPOOL_REUSE
  altsplice_chunk = (struct Altsplice_T *) List_head(this->altsplice_chunks);
  List_free_keep(&this->altsplice_free_cells);
  this->altsplice_free_cells = (List_T) NULL;
  for (celli = ALTSPLICE_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->altsplice_free_cells = List_push_keep(this->altsplice_free_cells,(void *) &(altsplice_chunk[celli]));
  }
#else
  this->altsplice_cellptr = (struct Altsplice_T *) List_head(this->altsplice_chunks);
  this->altsplice_cellctr = 0;
#endif


  while (List_next(this->string_chunks) != NULL) {
    this->string_chunks = List_pop_keep(this->string_chunks,(void **) &string_chunk);
    FREE_KEEP(string_chunk);
  }
  this->string_cellptr = (char *) List_head(this->string_chunks);
  this->string_chunksize = STRING_CHUNKSIZE;
  this->string_cellctr = 0;

  return;
}

void
Pathpool_free (T *old) {
  struct Path_T *path_chunk;
  struct Junction_T *junction_chunk;
  struct Altsplice_T *altsplice_chunk;
  char *string_chunk;

  while ((*old)->path_chunks != NULL) {
    (*old)->path_chunks = List_pop_keep((*old)->path_chunks,(void **) &path_chunk);
    FREE_KEEP(path_chunk);
  }
#ifdef PATHPOOL_REUSE
  List_free_keep(&(*old)->path_free_cells);
#endif

  while ((*old)->junction_chunks != NULL) {
    (*old)->junction_chunks = List_pop_keep((*old)->junction_chunks,(void **) &junction_chunk);
    FREE_KEEP(junction_chunk);
  }
#ifdef PATHPOOL_REUSE
  List_free_keep(&(*old)->junction_free_cells);
#endif

  while ((*old)->altsplice_chunks != NULL) {
    (*old)->altsplice_chunks = List_pop_keep((*old)->altsplice_chunks,(void **) &altsplice_chunk);
    FREE_KEEP(altsplice_chunk);
  }
#ifdef PATHPOOL_REUSE
  List_free_keep(&(*old)->altsplice_free_cells);
#endif  

  while ((*old)->string_chunks != NULL) {
    (*old)->string_chunks = List_pop_keep((*old)->string_chunks,(void **) &string_chunk);
    FREE_KEEP(string_chunk);
  }

  FREE_KEEP(*old);

  return;
}


T
Pathpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->path_chunks = (List_T) NULL;
#ifdef PATHPOOL_REUSE
  new->path_free_cells = (List_T) NULL;
#else
  new->path_cellctr = 0;
#endif

  new->junction_chunks = (List_T) NULL;
#ifdef PATHPOOL_REUSE
  new->junction_free_cells = (List_T) NULL;
#else
  new->junction_cellctr = 0;
#endif

  new->altsplice_chunks = (List_T) NULL;
#ifdef PATHPOOL_REUSE
  new->altsplice_free_cells = (List_T) NULL;
#else
  new->altsplice_cellctr = 0;
#endif

  new->string_chunks = (List_T) NULL;
  new->string_cellctr = 0;

  return new;
}


#ifdef PATHPOOL_REUSE
void
Pathpool_free_path (Path_T *old, T this
#ifdef PATHPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  this->path_free_cells = List_push_keep(this->path_free_cells,(void *) *old);
#ifdef PATHPOOL_TRACE
  printf("Pathpool/path: Freed %p -- Pathpool_free_path called by %s:%d\n",*old,file,line);
#endif

  *old = (Path_T) NULL;

  return;
}  
#endif


#ifdef PATHPOOL_REUSE
void
Pathpool_free_junction (Junction_T *old, T this
#ifdef PATHPOOL_TRACE
			, const char *file, int line
#endif
			) {
  this->junction_free_cells = List_push_keep(this->junction_free_cells,(void *) *old);
#ifdef PATHPOOL_TRACE
  printf("Pathpool/junction: Freed %p -- Pathpool_free_junction called by %s:%d\n",*old,file,line);
#endif

  *old = (Junction_T) NULL;
  return;
}  
#endif


#ifdef PATHPOOL_REUSE
void
Pathpool_free_altsplice (Altsplice_T *old, T this
#ifdef PATHPOOL_TRACE
			 , const char *file, int line
#endif
			 ) {
  this->altsplice_free_cells = List_push_keep(this->altsplice_free_cells,(void *) *old);
#ifdef PATHPOOL_TRACE
  printf("Pathpool/altsplice: Freed %p -- Pathpool_free_altsplice called by %s:%d\n",*old,file,line);
#endif

  *old = (Altsplice_T) NULL;
  return;
}  
#endif


static struct Path_T *
add_new_path_chunk (T this) {
  struct Path_T *chunk;

  chunk = (struct Path_T *) MALLOC_KEEP(PATH_CHUNKSIZE*sizeof(struct Path_T));
  this->path_chunks = List_push_keep(this->path_chunks,(void *) chunk);
#ifdef PATHPOOL_REUSE
  int celli;
  for (celli = PATH_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->path_free_cells = List_push_keep(this->path_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of path_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->path_chunks),chunk));

  return chunk;
}

static struct Junction_T *
add_new_junction_chunk (T this) {
  struct Junction_T *chunk;

  chunk = (struct Junction_T *) MALLOC_KEEP(JUNCTION_CHUNKSIZE*sizeof(struct Junction_T));
  this->junction_chunks = List_push_keep(this->junction_chunks,(void *) chunk);
#ifdef PATHPOOL_REUSE
  int celli;
  for (celli = JUNCTION_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->junction_free_cells = List_push_keep(this->junction_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of junction_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->junction_chunks),chunk));

  return chunk;
}

static struct Altsplice_T *
add_new_altsplice_chunk (T this) {
  struct Altsplice_T *chunk;

  chunk = (struct Altsplice_T *) MALLOC_KEEP(ALTSPLICE_CHUNKSIZE*sizeof(struct Altsplice_T));
  this->altsplice_chunks = List_push_keep(this->altsplice_chunks,(void *) chunk);
#ifdef PATHPOOL_REUSE
  int celli;
  for (celli = ALTSPLICE_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->altsplice_free_cells = List_push_keep(this->altsplice_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of altsplice_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->altsplice_chunks),chunk));

  return chunk;
}


static char *
add_new_string_chunk (T this, int nchars) {
  char *chunk;

  if (nchars > STRING_CHUNKSIZE) {
    chunk = (char *) MALLOC_KEEP(nchars*sizeof(char));
    this->string_chunksize = nchars;
  } else {
    chunk = (char *) MALLOC_KEEP(STRING_CHUNKSIZE*sizeof(char));
    this->string_chunksize = STRING_CHUNKSIZE;
  }
  this->string_chunks = List_push_keep(this->string_chunks,(void *) chunk);
  debug(printf("Adding a new chunk of string_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->string_chunks),chunk));

  return chunk;
}


Path_T
Pathpool_new_path (T this
#ifdef PATHPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  Path_T new;

#ifdef PATHPOOL_REUSE
  if (this->path_free_cells == (List_T) NULL) {
    add_new_path_chunk(this);
  }
  this->path_free_cells = List_pop_keep(this->path_free_cells,(void **) &new);
#else
  if (this->path_cellctr >= PATH_CHUNKSIZE) {
    this->path_cellptr = add_new_path_chunk(this);
    this->path_cellctr = 0;
  }
  this->path_cellctr += 1;
  new = this->path_cellptr++;
#endif

#ifdef PATHPOOL_TRACE
  printf("Pathpool/path: Allocated %p -- Pathpool_new_path called by %s:%d\n",new,file,line);
#endif

  return new;
}  


Junction_T
Pathpool_new_junction (T this
#ifdef PATHPOOL_TRACE
		       , const char *file, int line
#endif
		       ) {
  Junction_T new;

#ifdef PATHPOOL_REUSE
  if (this->junction_free_cells == (List_T) NULL) {
    add_new_junction_chunk(this);
  }
  this->junction_free_cells = List_pop_keep(this->junction_free_cells,(void **) &new);
#else
  if (this->junction_cellctr >= JUNCTION_CHUNKSIZE) {
    this->junction_cellptr = add_new_junction_chunk(this);
    this->junction_cellctr = 0;
  }
  this->junction_cellctr += 1;
  new = this->junction_cellptr++;
#endif

#ifdef PATHPOOL_TRACE
  printf("Pathpool/junction: Allocated %p -- Pathpool_new_junction called by %s:%d\n",new,file,line);
#endif

  return new;
}  


Altsplice_T
Pathpool_new_altsplice (T this
#ifdef PATHPOOL_TRACE
			, const char *file, int line
#endif
			) {
  Altsplice_T new;

#ifdef PATHPOOL_REUSE
  if (this->altsplice_free_cells == (List_T) NULL) {
    add_new_altsplice_chunk(this);
  }
  this->altsplice_free_cells = List_pop_keep(this->altsplice_free_cells,(void **) &new);
#else
  if (this->altsplice_cellctr >= ALTSPLICE_CHUNKSIZE) {
    this->altsplice_cellptr = add_new_altsplice_chunk(this);
    this->altsplice_cellctr = 0;
  }
  this->altsplice_cellctr += 1;
  new = this->altsplice_cellptr++;
#endif

#ifdef PATHPOOL_TRACE
  printf("Pathpool/altsplice: Allocated %p -- Pathpool_new_altsplice called by %s:%d\n",new,file,line);
#endif

  return new;
}  


char *
Pathpool_new_string (T this, int nchars) {
  char *string;

  if (this->string_cellctr + nchars > this->string_chunksize) {
    /* this->string_cellptr = add_new_string_chunk(this,nchars); */
    /* inlined add_new_string_chunk */
    if (nchars > STRING_CHUNKSIZE) {
      this->string_cellptr = (char *) MALLOC_KEEP(nchars*sizeof(char));
      this->string_chunksize = nchars;
    } else {
      this->string_cellptr = (char *) MALLOC_KEEP(STRING_CHUNKSIZE*sizeof(char));
      this->string_chunksize = STRING_CHUNKSIZE;
    }
    this->string_chunks = List_push_keep(this->string_chunks,(void *) this->string_cellptr);
    
    this->string_cellctr = 0;
  }

  this->string_cellctr += nchars;
  string = this->string_cellptr;
  this->string_cellptr += nchars;

  return string;
}  



/* Guarantees that a chunk is available */
void
Pathpool_init (T this) {
#ifdef PATHPOOL_REUSE
  add_new_path_chunk(this);
  add_new_junction_chunk(this);
  add_new_altsplice_chunk(this);
#else
this->path_cellptr = add_new_path_chunk(this);
  this->path_cellctr = 0;

  this->junction_cellptr = add_new_junction_chunk(this);
  this->junction_cellctr = 0;

  this->altsplice_cellptr = add_new_altsplice_chunk(this);
  this->altsplice_cellctr = 0;
#endif

  this->string_cellptr = add_new_string_chunk(this,STRING_CHUNKSIZE);
  this->string_cellctr = 0;
  this->string_chunksize = STRING_CHUNKSIZE;

  return;
}

