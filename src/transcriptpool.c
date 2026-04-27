static char rcsid[] = "$Id: ed87d8eee6034d59755471bcbd9cfb38180e37e2 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "transcriptpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"


#define TRANSCRIPT_CHUNKSIZE 1024
#define EXON_CHUNKSIZE 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Transcriptpool_T
struct T {
  List_T transcript_chunks;
#ifdef TRANSCRIPTPOOL_REUSE
  List_T transcript_free_cells;
#else
  struct Transcript_T *transcript_cellptr;
  int transcript_cellctr;
#endif

  List_T exon_chunks;
#ifdef TRANSCRIPTPOOL_REUSE
  List_T exon_free_cells;
#else
  struct Exon_T *exon_cellptr;
  int exon_cellctr;
#endif
};


#if defined(TRANSCRIPTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
static int
find_transcript_chunk (T this, struct Transcript_T *free_cell) {
  int chunki;
  List_T p;
  struct Transcript_T *chunk;

  for (p = this->transcript_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Transcript_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[TRANSCRIPT_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}

static int
find_exon_chunk (T this, struct Exon_T *free_cell) {
  int chunki;
  List_T p;
  struct Exon_T *chunk;

  for (p = this->exon_chunks, chunki = 0; p != NULL; p = List_next(p), chunki++) {
    chunk = (struct Exon_T *) List_head(p);
    if (free_cell >= &(chunk[0]) && free_cell < &(chunk[EXON_CHUNKSIZE])) {
      return chunki;
    }
  }

  fprintf(stderr,"Could not find chunk for %p\n",free_cell);
  abort();
  return -1;
}


/* Checks to see if all memory has been returned to free_cells */
static void
check_transcript_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->transcript_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->transcript_free_cells) == TRANSCRIPT_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->transcript_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_transcript_chunk(this,(struct Transcript_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < TRANSCRIPT_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Transcript_T cells leaked in transcriptpool chunk %d\n",
		TRANSCRIPT_CHUNKSIZE - nfreecells[chunki],TRANSCRIPT_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}

static void
check_exon_memory (T this) {
  List_T p;
  int nchunks, chunki;
  int *nfreecells;

  if ((nchunks = List_length(this->exon_chunks)) == 0) {
    /* Skip */
  } else if (List_length(this->exon_free_cells) == EXON_CHUNKSIZE*nchunks) {
    /* Looks okay */
  } else {
    nfreecells = (int *) CALLOC(nchunks,sizeof(int));
  
    for (p = this->exon_free_cells; p != NULL; p = List_next(p)) {
      chunki = find_exon_chunk(this,(struct Exon_T *) List_head(p));
      nfreecells[chunki] += 1;
    }

    for (chunki = 0; chunki < nchunks; chunki++) {
      if (nfreecells[chunki] < EXON_CHUNKSIZE) {
	fprintf(stderr,"%d out of %d Exon_T cells leaked in transcriptpool chunk %d\n",
		EXON_CHUNKSIZE - nfreecells[chunki],EXON_CHUNKSIZE,chunki);
      }
    }

    FREE(nfreecells);
  }

  return;
}
#endif


void
Transcriptpool_reset_memory (T this) {
  struct Transcript_T *transcript_chunk;
  struct Exon_T *exon_chunk;

#if defined(TRANSCRIPTPOOL_REUSE) && defined(CHECK_ASSERTIONS)
  check_transcript_memory(this);
  check_exon_memory(this);
#endif

  while (List_next(this->transcript_chunks) != NULL) {
    this->transcript_chunks = List_pop_keep(this->transcript_chunks,(void **) &transcript_chunk);
    FREE_KEEP(transcript_chunk);
  }
#ifdef TRANSCRIPTPOOL_REUSE
  int celli;
  transcript_chunk = (struct Transcript_T *) List_head(this->transcript_chunks);
  List_free_keep(&this->transcript_free_cells);
  this->transcript_free_cells = (List_T) NULL;
  for (celli = TRANSCRIPT_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->transcript_free_cells = List_push_keep(this->transcript_free_cells,(void *) &(transcript_chunk[celli]));
  }
#else
  this->transcript_cellptr = (struct Transcript_T *) List_head(this->transcript_chunks);
  this->transcript_cellctr = 0;
#endif

  while (List_next(this->exon_chunks) != NULL) {
    this->exon_chunks = List_pop_keep(this->exon_chunks,(void **) &exon_chunk);
    FREE_KEEP(exon_chunk);
  }
#ifdef TRANSCRIPTPOOL_REUSE
  exon_chunk = (struct Exon_T *) List_head(this->exon_chunks);
  List_free_keep(&this->exon_free_cells);
  this->exon_free_cells = (List_T) NULL;
  for (celli = EXON_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->exon_free_cells = List_push_keep(this->exon_free_cells,(void *) &(exon_chunk[celli]));
  }
#else
  this->exon_cellptr = (struct Exon_T *) List_head(this->exon_chunks);
  this->exon_cellctr = 0;
#endif
  
  return;
}


void
Transcriptpool_free (T *old) {
  struct Transcript_T *transcript_chunk;
  struct Exon_T *exon_chunk;

  while ((*old)->transcript_chunks != NULL) {
    (*old)->transcript_chunks = List_pop_keep((*old)->transcript_chunks,(void **) &transcript_chunk);
    FREE_KEEP(transcript_chunk);
  }
#ifdef TRANSCRIPTPOOL_REUSE
  List_free_keep(&(*old)->transcript_free_cells);
#endif

  while ((*old)->exon_chunks != NULL) {
    (*old)->exon_chunks = List_pop_keep((*old)->exon_chunks,(void **) &exon_chunk);
    FREE_KEEP(exon_chunk);
  }
#ifdef TRANSCRIPTPOOL_REUSE
  List_free_keep(&(*old)->exon_free_cells);
#endif

  FREE_KEEP(*old);

  return;
}


T
Transcriptpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->transcript_chunks = (List_T) NULL;
#ifdef TRANSCRIPTPOOL_REUSE
  new->transcript_free_cells = (List_T) NULL;
#else
  new->transcript_cellctr = 0;
#endif

  new->exon_chunks = (List_T) NULL;
#ifdef TRANSCRIPTPOOL_REUSE
  new->exon_free_cells = (List_T) NULL;
#else
  new->exon_cellctr = 0;
#endif

  return new;
}


#ifdef TRANSCRIPTPOOL_REUSE
void
Transcriptpool_free_transcript (Transcript_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
				, const char *file, int line
#endif
				) {
  this->transcript_free_cells = List_push_keep(this->transcript_free_cells,(void *) *old);
#ifdef TRANSCRIPTPOOL_TRACE
  printf("Transcriptpool/transcript: Freed %p -- Transcriptpool_free_transcript called by %s:%d\n",
	 *old,file,line);
#endif

  *old = (Transcript_T) NULL;
  return;
}  
#endif


#ifdef TRANSCRIPTPOOL_REUSE
void
Transcriptpool_free_exon (Exon_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
			  , const char *file, int line
#endif
			  ) {
  this->exon_free_cells = List_push_keep(this->exon_free_cells,(void *) *old);
#ifdef TRANSCRIPTPOOL_TRACE
  printf("Transcriptpool/exon: Freed %p -- Transcriptpool_free_exon called by %s:%d\n",*old,file,line);
#endif

  *old = (Exon_T) NULL;
  return;
}  
#endif


static struct Transcript_T *
add_new_transcript_chunk (T this) {
  struct Transcript_T *chunk;

  chunk = (struct Transcript_T *) MALLOC_KEEP(TRANSCRIPT_CHUNKSIZE*sizeof(struct Transcript_T));
  this->transcript_chunks = List_push_keep(this->transcript_chunks,(void *) chunk);
#ifdef TRANSCRIPTPOOL_REUSE
  int celli;
  for (celli = TRANSCRIPT_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->transcript_free_cells = List_push_keep(this->transcript_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of transcript_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->transcript_chunks),chunk));

  return chunk;
}

static struct Exon_T *
add_new_exon_chunk (T this) {
  struct Exon_T *chunk;

  chunk = (struct Exon_T *) MALLOC_KEEP(EXON_CHUNKSIZE*sizeof(struct Exon_T));
  this->exon_chunks = List_push_keep(this->exon_chunks,(void *) chunk);
#ifdef TRANSCRIPTPOOL_REUSE
  int celli;
  for (celli = EXON_CHUNKSIZE - 1; celli >= 0; celli--) {
    this->exon_free_cells = List_push_keep(this->exon_free_cells,(void *) &(chunk[celli]));
  }
#endif

  debug(printf("Adding a new chunk of exon_cells.  Ptr for chunk %d is %p\n",
	       List_length(this->exon_chunks),chunk));

  return chunk;
}

Transcript_T
Transcriptpool_new_transcript (T this
#ifdef TRANSCRIPTPOOL_TRACE
			       , const char *file, int line
#endif
			       ) {
  Transcript_T new;

#ifdef TRANSCRIPTPOOL_REUSE
  if (this->transcript_free_cells == (List_T) NULL) {
    add_new_transcript_chunk(this);
  }
  this->transcript_free_cells = List_pop_keep(this->transcript_free_cells,(void **) &new);
#else
  if (this->transcript_cellctr >= TRANSCRIPT_CHUNKSIZE) {
    this->transcript_cellptr = add_new_transcript_chunk(this);
    this->transcript_cellctr = 0;
  }
  this->transcript_cellctr += 1;
  new = this->transcript_cellptr++;
#endif

#ifdef TRANSCRIPTPOOL_TRACE
  printf("Transcriptpool/transcript: Allocated %p -- Transcriptpool_new_transcript called by %s:%d\n",new,file,line);
#endif

  return new;
}  


Exon_T
Transcriptpool_new_exon (T this
#ifdef TRANSCRIPTPOOL_TRACE
			 , const char *file, int line
#endif
			 ) {
  Exon_T new;

#ifdef TRANSCRIPTPOOL_REUSE
  if (this->exon_free_cells == (List_T) NULL) {
    add_new_exon_chunk(this);
  }
  this->exon_free_cells = List_pop_keep(this->exon_free_cells,(void **) &new);
#else
  if (this->exon_cellctr >= EXON_CHUNKSIZE) {
    this->exon_cellptr = add_new_exon_chunk(this);
    this->exon_cellctr = 0;
  }
  this->exon_cellctr += 1;
  new = this->exon_cellptr++;
#endif

#ifdef TRANSCRIPTPOOL_TRACE
  printf("Transcriptpool/exon: Allocated %p -- Transcriptpool_new_exon called by %s:%d\n",new,file,line);
#endif

  return new;
}  


/* Guarantees that a chunk is available */
void
Transcriptpool_init (T this) {
#ifdef TRANSCRIPTPOOL_REUSE
  add_new_transcript_chunk(this);
  add_new_exon_chunk(this);
#else
  this->transcript_cellptr = add_new_transcript_chunk(this);
  this->transcript_cellctr = 0;
  
  this->exon_cellptr = add_new_exon_chunk(this);
  this->exon_cellctr = 0;
#endif
  return;
}

