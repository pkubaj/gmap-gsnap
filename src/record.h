/* $Id: 5eccdb765d5f422afffce14e0641423c936b46cd $ */
#ifndef RECORD_INCLUDED
#define RECORD_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "chrnum.h"

typedef struct Record_T *Record_T;
struct Record_T {
  Univcoord_T univdiagonal;		/* Primary sort */
  int qstart;			/* Secondary sort */
  int qend;
  bool anchorp;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;

  Univcoord_T lowpos;
  Univcoord_T highpos;
};

#endif

