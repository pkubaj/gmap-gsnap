#ifndef MERGEINFO_INCLUDED
#define MERGEINFO_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct Mergeinfo_uint4_T *Mergeinfo_uint4_T;
typedef struct Mergeinfo_uint8_T *Mergeinfo_uint8_T;
#ifdef LARGE_GENOMES
typedef struct Mergeinfo_uint8_T *Mergeinfo_T;
#else
typedef struct Mergeinfo_uint4_T *Mergeinfo_T;
#endif

#include "types.h"
#include "genomicpos.h"
#include "univcoord.h"


struct Mergeinfo_uint4_T {
  int querylength;

  int max_nstreams;
  int max_heapsize;

  int *streamsize_copy;
  UINT4 **combined;
  int *totals;
  UINT4 **heap;
  int *nelts;
};

struct Mergeinfo_uint8_T {
  int querylength;

  int max_nstreams;
  int max_heapsize;

  int *streamsize_copy;
  UINT8 **combined;
  int *totals;
  UINT8 **heap;
  int *nelts;
};

extern void
Mergeinfo_uint4_free (Mergeinfo_uint4_T *old);
extern Mergeinfo_uint4_T
Mergeinfo_uint4_new (int querylength, Chrpos_T max_localdb_distance);

extern void
Mergeinfo_uint8_free (Mergeinfo_uint8_T *old);
extern Mergeinfo_uint8_T
Mergeinfo_uint8_new (int querylength, Chrpos_T max_localdb_distance);

#endif

