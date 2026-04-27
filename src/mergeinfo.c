static char rcsid[] = "$Id: f21debc5a1150ebedc73d1fd9c84450f096f1239 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "mergeinfo.h"

#include "mem.h"

#define LOCALDB_REGION_SIZE 65536

void
Mergeinfo_uint4_free (Mergeinfo_uint4_T *old) {
  FREE((*old)->streamsize_copy);
  FREE((*old)->combined);
  FREE((*old)->totals);
  FREE((*old)->heap);
  FREE((*old)->nelts);
  FREE(*old);
  return;
}


Mergeinfo_uint4_T
Mergeinfo_uint4_new (int querylength, Chrpos_T max_localdb_distance) {
  Mergeinfo_uint4_T new = (Mergeinfo_uint4_T) MALLOC(sizeof(*new));

  int max_localdb_nregions = (max_localdb_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  int max_nstreams = max_localdb_nregions * querylength;
  int max_heapsize = 2*max_nstreams + 1;

  new->querylength = querylength;
  new->max_nstreams = max_nstreams;
  new->max_heapsize = max_heapsize;

  new->streamsize_copy = (int *) MALLOC(max_nstreams * sizeof(int));
  new->combined = (UINT4 **) MALLOC(max_nstreams * sizeof(UINT4 *));
  new->totals = (int *) MALLOC(max_nstreams * sizeof(int));
  new->heap = (UINT4 **) MALLOC((max_heapsize + 1)*sizeof(UINT4 *));
  new->nelts = (int *) MALLOC((max_heapsize + 1)*sizeof(int));

  return new;
}

void
Mergeinfo_uint8_free (Mergeinfo_uint8_T *old) {
  FREE((*old)->streamsize_copy);
  FREE((*old)->combined);
  FREE((*old)->totals);
  FREE((*old)->heap);
  FREE((*old)->nelts);
  FREE(*old);
  return;
}


Mergeinfo_uint8_T
Mergeinfo_uint8_new (int querylength, Chrpos_T max_localdb_distance) {
  Mergeinfo_uint8_T new = (Mergeinfo_uint8_T) MALLOC(sizeof(*new));

  int max_localdb_nregions = (max_localdb_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
  int max_nstreams = max_localdb_nregions * querylength;
  int max_heapsize = 2*max_nstreams + 1;

  new->querylength = querylength;
  new->max_nstreams = max_nstreams;
  new->max_heapsize = max_heapsize;

  new->streamsize_copy = (int *) MALLOC(max_nstreams * sizeof(int));
  new->combined = (UINT8 **) MALLOC(max_nstreams * sizeof(UINT8 *));
  new->totals = (int *) MALLOC(max_nstreams * sizeof(int));
  new->heap = (UINT8 **) MALLOC((max_heapsize + 1)*sizeof(UINT8 *));
  new->nelts = (int *) MALLOC((max_heapsize + 1)*sizeof(int));

  return new;
}



