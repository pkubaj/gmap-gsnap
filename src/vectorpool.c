#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "vectorpool.h"
#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "comp.h"
#include "list.h"

#define INT_CHUNKSIZE 2048	
#define UINT_CHUNKSIZE 2048	
#define UNIVCOORD_CHUNKSIZE 2048
#define DOUBLE_CHUNKSIZE 1024


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* For mechanics of memory allocation and deallocation */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Vectorpool_T
struct T {
  int int_chunksize;
  int int_cellctr;
  int *int_cellptr;
  List_T int_chunks;

  int uint_chunksize;
  int uint_cellctr;
  unsigned int *uint_cellptr;
  List_T uint_chunks;

  int univcoord_chunksize;
  int univcoord_cellctr;
  Univcoord_T *univcoord_cellptr;
  List_T univcoord_chunks;

  int double_chunksize;
  int double_cellctr;
  double *double_cellptr;
  List_T double_chunks;
};


void
Vectorpool_reset_memory (T this) {
  int *int_chunk;
  unsigned int *uint_chunk;
  Univcoord_T *univcoord_chunk;
  double *double_chunk;

  /* Previously List_next(this->int_chunks) != NULL, and so on */
  while (this->int_chunks != NULL) {
    this->int_chunks = List_pop(this->int_chunks,(void **) &int_chunk);
    FREE(int_chunk);
  }
  /* this->int_cellptr = (int *) List_head(this->int_chunks); */
  this->int_cellptr = (int *) NULL;
  this->int_chunksize = INT_CHUNKSIZE;
  this->int_cellctr = 0;


  while (this->uint_chunks != NULL) {
    this->uint_chunks = List_pop(this->uint_chunks,(void **) &uint_chunk);
    FREE(uint_chunk);
  }
  /* this->uint_cellptr = (unsigned int *) List_head(this->uint_chunks); */
  this->uint_cellptr = (unsigned int *) NULL;
  this->uint_chunksize = UINT_CHUNKSIZE;
  this->uint_cellctr = 0;


  while (this->univcoord_chunks != NULL) {
    this->univcoord_chunks = List_pop(this->univcoord_chunks,(void **) &univcoord_chunk);
    FREE(univcoord_chunk);
  }
  /* this->univcoord_cellptr = (Univcoord_T *) List_head(this->univcoord_chunks); */
  this->univcoord_cellptr = (Univcoord_T *) NULL;
  this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  this->univcoord_cellctr = 0;
  

  while (this->double_chunks != NULL) {
    this->double_chunks = List_pop(this->double_chunks,(void **) &double_chunk);
    FREE(double_chunk);
  }
  /* this->double_cellptr = (double *) List_head(this->double_chunks); */
  this->double_cellptr = (double *) NULL;
  this->double_chunksize = DOUBLE_CHUNKSIZE;
  this->double_cellctr = 0;

  return;
}

void
Vectorpool_free (T *old) {
  int *int_chunk;
  unsigned int *uint_chunk;
  Univcoord_T *univcoord_chunk;
  double *double_chunk;

  while ((*old)->int_chunks != NULL) {
    (*old)->int_chunks = List_pop((*old)->int_chunks,(void **) &int_chunk);
    FREE(int_chunk);
  }

  while ((*old)->uint_chunks != NULL) {
    (*old)->uint_chunks = List_pop((*old)->uint_chunks,(void **) &uint_chunk);
    FREE(uint_chunk);
  }

  while ((*old)->univcoord_chunks != NULL) {
    (*old)->univcoord_chunks = List_pop((*old)->univcoord_chunks,(void **) &univcoord_chunk);
    FREE(univcoord_chunk);
  }

  while ((*old)->double_chunks != NULL) {
    (*old)->double_chunks = List_pop((*old)->double_chunks,(void **) &double_chunk);
    FREE(double_chunk);
  }

  FREE_KEEP(*old);

  return;
}


#if 0
/* Code is inlined below */
static int *
add_new_int_chunk (T this, int nints) {
  int *int_chunk;

  if (nints > INT_CHUNKSIZE) {
    int_chunk = (int *) MALLOC(nints*sizeof(int));
    this->int_chunksize = nints;
  } else {
    int_chunk = (int *) MALLOC(INT_CHUNKSIZE*sizeof(int));
    this->int_chunksize = INT_CHUNKSIZE;
  }
  this->int_chunks = List_push(this->int_chunks,(void *) int_chunk);
  debug1(printf("Adding a new chunk of int_cells.  Ptr for chunk %d is %p\n",
		List_length(this->int_chunks),chunk));

  return int_chunk;
}
#endif


#if 0
/* Code is inlined below */
static unsigned int *
add_new_uint_chunk (T this, int nuints) {
  unsigned int *uint_chunk;

  if (nuints > UINT_CHUNKSIZE) {
    uint_chunk = (unsigned int *) MALLOC(nuints*sizeof(unsigned int));
    this->uint_chunksize = nuints;
  } else {
    uint_chunk = (unsigned int *) MALLOC(UINT_CHUNKSIZE*sizeof(unsigned int));
    this->uint_chunksize = UINT_CHUNKSIZE;
  }
  this->uint_chunks = List_push(this->uint_chunks,(void *) uint_chunk);
  debug1(printf("Adding a new chunk of uint_cells.  Ptr for chunk %d is %p\n",
		List_length(this->uint_chunks),chunk));

  return uint_chunk;
}
#endif


#if 0
/* Code is inlined below */
static Univcoord_T *
add_new_univcoord_chunk (T this, int nunivcoords) {
  Univcoord_T *univcoord_chunk;

  if (nunivcoords > UNIVCOORD_CHUNKSIZE) {
    univcoord_chunk = (Univcoord_T *) MALLOC(nunivcoords*sizeof(Univcoord_T));
    this->univcoord_chunksize = nunivcoords;
  } else {
    univcoord_chunk = (Univcoord_T *) MALLOC(UNIVCOORD_CHUNKSIZE*sizeof(Univcoord_T));
    this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  }
  this->univcoord_chunks = List_push(this->univcoord_chunks,(void *) univcoord_chunk);
  debug1(printf("Adding a new chunk of univcoord_cells.  Ptr for chunk %d is %p\n",
		List_length(this->univcoord_chunks),chunk));

  return univcoord_chunk;
}
#endif


#if 0
/* Code is inlined below */
static double *
add_new_double_chunk (T this, int ndoubles) {
  double *double_chunk;

  if (ndoubles > DOUBLE_CHUNKSIZE) {
    double_chunk = (double *) MALLOC(ndoubles*sizeof(double));
    this->double_chunksize = ndoubles;
  } else {
    double_chunk = (double *) MALLOC(DOUBLE_CHUNKSIZE*sizeof(double));
    this->double_chunksize = DOUBLE_CHUNKSIZE;
  }
  this->double_chunks = List_push(this->double_chunks,(void *) double_chunk);
  debug1(printf("Adding a new chunk of double_cells.  Ptr for chunk %d is %p\n",
		List_length(this->double_chunks),chunk));

  return double_chunk;
}
#endif


T
Vectorpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->int_cellptr = (int *) NULL;
  new->int_cellctr = 0;
  new->int_chunks = (List_T) NULL;

  new->uint_cellptr = (unsigned int *) NULL;
  new->uint_cellctr = 0;
  new->uint_chunks = (List_T) NULL;

  new->univcoord_cellptr = (Univcoord_T *) NULL;
  new->univcoord_cellctr = 0;
  new->univcoord_chunks = (List_T) NULL;

  new->double_cellptr = (double *) NULL;
  new->double_cellctr = 0;
  new->double_chunks = (List_T) NULL;

  return new;
}


int *
Vectorpool_new_intvector (T this, int nints) {
  int *intvector;

  if (this->int_cellptr == (int *) NULL ||
      this->int_cellctr + nints > this->int_chunksize) {
    /* this->int_cellptr = add_new_int_chunk(this,nints); */
    /* inlined add_new_int_chunk */
    if (nints > this->int_chunksize) {
      this->int_cellptr = (int *) MALLOC(nints*sizeof(int));
      this->int_chunksize = nints;
    } else {
      this->int_cellptr = (int *) MALLOC(this->int_chunksize*sizeof(int));
    }
    this->int_chunks = List_push(this->int_chunks,(void *) this->int_cellptr);
    
    this->int_cellctr = 0;
  }

  this->int_cellctr += nints;

  intvector = this->int_cellptr;
  this->int_cellptr += nints;
  return intvector;
}  


unsigned int *
Vectorpool_new_uintvector (T this, int nuints) {
  unsigned int *uintvector;

  if (this->uint_cellptr == (unsigned int *) NULL ||
      this->uint_cellctr + nuints > this->uint_chunksize) {
    /* this->uint_cellptr = add_new_uint_chunk(this,nuints); */
    /* inlined add_new_uint_chunk */
    if (nuints > this->uint_chunksize) {
      this->uint_cellptr = (unsigned int *) MALLOC(nuints*sizeof(unsigned int));
      this->uint_chunksize = nuints;
    } else {
      this->uint_cellptr = (unsigned int *) MALLOC(this->uint_chunksize*sizeof(unsigned int));
    }
    this->uint_chunks = List_push(this->uint_chunks,(void *) this->uint_cellptr);
    
    this->uint_cellctr = 0;
  }

  this->uint_cellctr += nuints;

  uintvector = this->uint_cellptr;
  this->uint_cellptr += nuints;
  return uintvector;
}  


Univcoord_T *
Vectorpool_new_univcoordvector (T this, int nunivcoords) {
  Univcoord_T *univcoordvector;

  if (this->univcoord_cellptr == (Univcoord_T *) NULL ||
      this->univcoord_cellctr + nunivcoords > this->univcoord_chunksize) {
    /* this->univcoord_cellptr = add_new_univcoord_chunk(this,nunivcoords); */
    /* inlined add_new_univcoord_chunk */
    if (nunivcoords > this->univcoord_chunksize) {
      this->univcoord_cellptr = (Univcoord_T *) MALLOC(nunivcoords*sizeof(Univcoord_T));
      this->univcoord_chunksize = nunivcoords;
    } else {
      this->univcoord_cellptr = (Univcoord_T *) MALLOC(this->univcoord_chunksize*sizeof(Univcoord_T));
    }
    this->univcoord_chunks = List_push(this->univcoord_chunks,(void *) this->univcoord_cellptr);
    
    this->univcoord_cellctr = 0;
  }

  this->univcoord_cellctr += nunivcoords;

  univcoordvector = this->univcoord_cellptr;
  this->univcoord_cellptr += nunivcoords;
  return univcoordvector;
}  


double *
Vectorpool_new_doublevector (T this, int ndoubles) {
  double *doublevector;

  if (this->double_cellptr == (double *) NULL ||
      this->double_cellctr + ndoubles > this->double_chunksize) {
    /* this->double_cellptr = add_new_double_chunk(this,ndoubles); */
    /* inlined add_new_double_chunk */
    if (ndoubles > this->double_chunksize) {
      this->double_cellptr = (double *) MALLOC(ndoubles*sizeof(double));
      this->double_chunksize = ndoubles;
    } else {
      this->double_cellptr = (double *) MALLOC(this->double_chunksize*sizeof(double));
    }
    this->double_chunks = List_push(this->double_chunks,(void *) this->double_cellptr);
    
    this->double_cellctr = 0;
  }

  this->double_cellctr += ndoubles;

  doublevector = this->double_cellptr;
  this->double_cellptr += ndoubles;
  return doublevector;
}  



/* Previously started with an initial chunk */
void
Vectorpool_init (T this) {
  /* this->int_cellptr = add_new_int_chunk(this,INT_CHUNKSIZE); */
  this->int_cellptr = (int *) NULL;
  this->int_chunksize = INT_CHUNKSIZE;
  this->int_cellctr = 0;

  /* this->uint_cellptr = add_new_uint_chunk(this,UINT_CHUNKSIZE); */
  this->uint_cellptr = (unsigned int *) NULL;
  this->uint_chunksize = UINT_CHUNKSIZE;
  this->uint_cellctr = 0;

  /* this->univcoord_cellptr = add_new_univcoord_chunk(this,UNIVCOORD_CHUNKSIZE); */
  this->univcoord_cellptr = (Univcoord_T *) NULL;
  this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  this->univcoord_cellctr = 0;

  /* this->double_cellptr = add_new_double_chunk(this,DOUBLE_CHUNKSIZE); */
  this->double_cellptr = (double *) NULL;
  this->double_chunksize = DOUBLE_CHUNKSIZE;
  this->double_cellctr = 0;

  return;
}

