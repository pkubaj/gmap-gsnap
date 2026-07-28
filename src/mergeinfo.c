#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mergeinfo.h"

#include "mem.h"

#define LOCALDB_REGION_SIZE 65536

/* Use keep memory pool, because this is retained if querylength stays the same */
void
Mergeinfo_uint4_free (Mergeinfo_uint4_T *old) {
  FREE_KEEP((*old)->streamsize_copy);
  FREE_KEEP((*old)->combined);
  FREE_KEEP((*old)->totals);
  FREE_KEEP((*old)->heap);
  FREE_KEEP((*old)->nelts);
  FREE_KEEP(*old);
  return;
}


/* Use keep memory pool, because this is retained if querylength stays the same */
Mergeinfo_uint4_T
Mergeinfo_uint4_new (int querylength, Chrpos_T max_localdb_distance) {
  Mergeinfo_uint4_T new = (Mergeinfo_uint4_T) MALLOC_KEEP(sizeof(*new));

  int max_localdb_nregions = (max_localdb_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  int max_nstreams = max_localdb_nregions * querylength;
  int max_heapsize = 2*max_nstreams + 1;

  if (querylength == 0) {
    max_nstreams = 1;
  }

  new->querylength = querylength;
  new->max_nstreams = max_nstreams;
  new->max_heapsize = max_heapsize;

  new->streamsize_copy = (int *) MALLOC_KEEP(max_nstreams * sizeof(int));
  new->combined = (UINT4 **) MALLOC_KEEP(max_nstreams * sizeof(UINT4 *));
  new->totals = (int *) MALLOC_KEEP(max_nstreams * sizeof(int));
  new->heap = (UINT4 **) MALLOC_KEEP((max_heapsize + 1)*sizeof(UINT4 *));
  new->nelts = (int *) MALLOC_KEEP((max_heapsize + 1)*sizeof(int));

  return new;
}

void
Mergeinfo_uint8_free (Mergeinfo_uint8_T *old) {
  FREE_KEEP((*old)->streamsize_copy);
  FREE_KEEP((*old)->combined);
  FREE_KEEP((*old)->totals);
  FREE_KEEP((*old)->heap);
  FREE_KEEP((*old)->nelts);
  FREE_KEEP(*old);
  return;
}


Mergeinfo_uint8_T
Mergeinfo_uint8_new (int querylength, Chrpos_T max_localdb_distance) {
  Mergeinfo_uint8_T new = (Mergeinfo_uint8_T) MALLOC_KEEP(sizeof(*new));

  int max_localdb_nregions = (max_localdb_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  int max_nstreams = max_localdb_nregions * querylength;
  int max_heapsize = 2*max_nstreams + 1;

  if (querylength == 0) {
    max_nstreams = 1;
  }

  new->querylength = querylength;
  new->max_nstreams = max_nstreams;
  new->max_heapsize = max_heapsize;

  new->streamsize_copy = (int *) MALLOC_KEEP(max_nstreams * sizeof(int));
  new->combined = (UINT8 **) MALLOC_KEEP(max_nstreams * sizeof(UINT8 *));
  new->totals = (int *) MALLOC_KEEP(max_nstreams * sizeof(int));
  new->heap = (UINT8 **) MALLOC_KEEP((max_heapsize + 1)*sizeof(UINT8 *));
  new->nelts = (int *) MALLOC_KEEP((max_heapsize + 1)*sizeof(int));

  return new;
}



