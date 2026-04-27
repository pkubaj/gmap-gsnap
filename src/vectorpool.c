static char rcsid[] = "$Id: 5cf8de22f27c5a4f69227e4b8eed112256a68c2d $";
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

  while (List_next(this->int_chunks) != NULL) {
    this->int_chunks = List_pop_keep(this->int_chunks,(void **) &int_chunk);
    FREE_KEEP(int_chunk);
  }
  this->int_cellptr = (int *) List_head(this->int_chunks);
  this->int_chunksize = INT_CHUNKSIZE;
  this->int_cellctr = 0;


  while (List_next(this->uint_chunks) != NULL) {
    this->uint_chunks = List_pop_keep(this->uint_chunks,(void **) &uint_chunk);
    FREE_KEEP(uint_chunk);
  }
  this->uint_cellptr = (unsigned int *) List_head(this->uint_chunks);
  this->uint_chunksize = UINT_CHUNKSIZE;
  this->uint_cellctr = 0;


  while (List_next(this->univcoord_chunks) != NULL) {
    this->univcoord_chunks = List_pop_keep(this->univcoord_chunks,(void **) &univcoord_chunk);
    FREE_KEEP(univcoord_chunk);
  }
  this->univcoord_cellptr = (Univcoord_T *) List_head(this->univcoord_chunks);
  this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  this->univcoord_cellctr = 0;
  

  while (List_next(this->double_chunks) != NULL) {
    this->double_chunks = List_pop_keep(this->double_chunks,(void **) &double_chunk);
    FREE_KEEP(double_chunk);
  }
  this->double_cellptr = (double *) List_head(this->double_chunks);
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
    (*old)->int_chunks = List_pop_keep((*old)->int_chunks,(void **) &int_chunk);
    FREE_KEEP(int_chunk);
  }

  while ((*old)->uint_chunks != NULL) {
    (*old)->uint_chunks = List_pop_keep((*old)->uint_chunks,(void **) &uint_chunk);
    FREE_KEEP(uint_chunk);
  }

  while ((*old)->univcoord_chunks != NULL) {
    (*old)->univcoord_chunks = List_pop_keep((*old)->univcoord_chunks,(void **) &univcoord_chunk);
    FREE_KEEP(univcoord_chunk);
  }

  while ((*old)->double_chunks != NULL) {
    (*old)->double_chunks = List_pop_keep((*old)->double_chunks,(void **) &double_chunk);
    FREE_KEEP(double_chunk);
  }

  FREE_KEEP(*old);

  return;
}


static int *
add_new_int_chunk (T this, int nints) {
  int *int_chunk;

  if (nints > INT_CHUNKSIZE) {
    int_chunk = (int *) MALLOC_KEEP(nints*sizeof(int));
    this->int_chunksize = nints;
  } else {
    int_chunk = (int *) MALLOC_KEEP(INT_CHUNKSIZE*sizeof(int));
    this->int_chunksize = INT_CHUNKSIZE;
  }
  this->int_chunks = List_push_keep(this->int_chunks,(void *) int_chunk);
  debug1(printf("Adding a new chunk of int_cells.  Ptr for chunk %d is %p\n",
		List_length(this->int_chunks),chunk));

  return int_chunk;
}


static unsigned int *
add_new_uint_chunk (T this, int nuints) {
  unsigned int *uint_chunk;

  if (nuints > UINT_CHUNKSIZE) {
    uint_chunk = (unsigned int *) MALLOC_KEEP(nuints*sizeof(unsigned int));
    this->uint_chunksize = nuints;
  } else {
    uint_chunk = (unsigned int *) MALLOC_KEEP(UINT_CHUNKSIZE*sizeof(unsigned int));
    this->uint_chunksize = UINT_CHUNKSIZE;
  }
  this->uint_chunks = List_push_keep(this->uint_chunks,(void *) uint_chunk);
  debug1(printf("Adding a new chunk of uint_cells.  Ptr for chunk %d is %p\n",
		List_length(this->uint_chunks),chunk));

  return uint_chunk;
}


static Univcoord_T *
add_new_univcoord_chunk (T this, int nunivcoords) {
  Univcoord_T *univcoord_chunk;

  if (nunivcoords > UNIVCOORD_CHUNKSIZE) {
    univcoord_chunk = (Univcoord_T *) MALLOC_KEEP(nunivcoords*sizeof(Univcoord_T));
    this->univcoord_chunksize = nunivcoords;
  } else {
    univcoord_chunk = (Univcoord_T *) MALLOC_KEEP(UNIVCOORD_CHUNKSIZE*sizeof(Univcoord_T));
    this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  }
  this->univcoord_chunks = List_push_keep(this->univcoord_chunks,(void *) univcoord_chunk);
  debug1(printf("Adding a new chunk of univcoord_cells.  Ptr for chunk %d is %p\n",
		List_length(this->univcoord_chunks),chunk));

  return univcoord_chunk;
}


static double *
add_new_double_chunk (T this, int ndoubles) {
  double *double_chunk;

  if (ndoubles > DOUBLE_CHUNKSIZE) {
    double_chunk = (double *) MALLOC_KEEP(ndoubles*sizeof(double));
    this->double_chunksize = ndoubles;
  } else {
    double_chunk = (double *) MALLOC_KEEP(DOUBLE_CHUNKSIZE*sizeof(double));
    this->double_chunksize = DOUBLE_CHUNKSIZE;
  }
  this->double_chunks = List_push_keep(this->double_chunks,(void *) double_chunk);
  debug1(printf("Adding a new chunk of double_cells.  Ptr for chunk %d is %p\n",
		List_length(this->double_chunks),chunk));

  return double_chunk;
}


T
Vectorpool_new (void) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->int_cellctr = 0;
  new->int_chunks = (List_T) NULL;

  new->uint_cellctr = 0;
  new->uint_chunks = (List_T) NULL;

  new->univcoord_cellctr = 0;
  new->univcoord_chunks = (List_T) NULL;

  new->double_cellctr = 0;
  new->double_chunks = (List_T) NULL;

  return new;
}


int *
Vectorpool_new_intvector (T this, int nints) {
  int *intvector;

  if (this->int_cellctr + nints > this->int_chunksize) {
    /* this->int_cellptr = add_new_int_chunk(this,nints); */
    /* inlined add_new_int_chunk */
    if (nints > INT_CHUNKSIZE) {
      this->int_cellptr = (int *) MALLOC_KEEP(nints*sizeof(int));
      this->int_chunksize = nints;
    } else {
      this->int_cellptr = (int *) MALLOC_KEEP(INT_CHUNKSIZE*sizeof(int));
      this->int_chunksize = INT_CHUNKSIZE;
    }
    this->int_chunks = List_push_keep(this->int_chunks,(void *) this->int_cellptr);
    
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

  if (this->uint_cellctr + nuints > this->uint_chunksize) {
    /* this->uint_cellptr = add_new_uint_chunk(this,nuints); */
    /* inlined add_new_uint_chunk */
    if (nuints > UINT_CHUNKSIZE) {
      this->uint_cellptr = (unsigned int *) MALLOC_KEEP(nuints*sizeof(unsigned int));
      this->uint_chunksize = nuints;
    } else {
      this->uint_cellptr = (unsigned int *) MALLOC_KEEP(UINT_CHUNKSIZE*sizeof(unsigned int));
      this->uint_chunksize = UINT_CHUNKSIZE;
    }
    this->uint_chunks = List_push_keep(this->uint_chunks,(void *) this->uint_cellptr);
    
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

  if (this->univcoord_cellctr + nunivcoords > this->univcoord_chunksize) {
    /* this->univcoord_cellptr = add_new_univcoord_chunk(this,nunivcoords); */
    /* inlined add_new_univcoord_chunk */
    if (nunivcoords > UNIVCOORD_CHUNKSIZE) {
      this->univcoord_cellptr = (Univcoord_T *) MALLOC_KEEP(nunivcoords*sizeof(Univcoord_T));
      this->univcoord_chunksize = nunivcoords;
    } else {
      this->univcoord_cellptr = (Univcoord_T *) MALLOC_KEEP(UNIVCOORD_CHUNKSIZE*sizeof(Univcoord_T));
      this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
    }
    this->univcoord_chunks = List_push_keep(this->univcoord_chunks,(void *) this->univcoord_cellptr);
    
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

  if (this->double_cellctr + ndoubles > this->double_chunksize) {
    /* this->double_cellptr = add_new_double_chunk(this,ndoubles); */
    /* inlined add_new_double_chunk */
    if (ndoubles > DOUBLE_CHUNKSIZE) {
      this->double_cellptr = (double *) MALLOC_KEEP(ndoubles*sizeof(double));
      this->double_chunksize = ndoubles;
    } else {
      this->double_cellptr = (double *) MALLOC_KEEP(DOUBLE_CHUNKSIZE*sizeof(double));
      this->double_chunksize = DOUBLE_CHUNKSIZE;
    }
    this->double_chunks = List_push_keep(this->double_chunks,(void *) this->double_cellptr);
    
    this->double_cellctr = 0;
  }

  this->double_cellctr += ndoubles;
  doublevector = this->double_cellptr;
  this->double_cellptr += ndoubles;

  return doublevector;
}  




/* Guarantees that a chunk is available */
void
Vectorpool_init (T this) {
  this->int_cellptr = add_new_int_chunk(this,INT_CHUNKSIZE);
  this->int_chunksize = INT_CHUNKSIZE;
  this->int_cellctr = 0;

  this->uint_cellptr = add_new_uint_chunk(this,UINT_CHUNKSIZE);
  this->uint_chunksize = UINT_CHUNKSIZE;
  this->uint_cellctr = 0;

  this->univcoord_cellptr = add_new_univcoord_chunk(this,UNIVCOORD_CHUNKSIZE);
  this->univcoord_chunksize = UNIVCOORD_CHUNKSIZE;
  this->univcoord_cellctr = 0;

  this->double_cellptr = add_new_double_chunk(this,DOUBLE_CHUNKSIZE);
  this->double_chunksize = DOUBLE_CHUNKSIZE;
  this->double_cellctr = 0;

  return;
}

