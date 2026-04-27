/* $Id: univdiagdef.h 225773 2023-01-27 19:12:35Z twu $ */
#ifndef UNIVDIAGDEF_INCLUDED
#define UNIVDIAGDEF_INCLUDED

#include "bool.h"
#include "univcoord.h"


/* qstart and qend are the genome-normalized coordinates, so qstart
   marks the left coordinate and qend marks the right coordinate.  For
   a plus-strand alignment, qstart = querystart and qend = queryend.
   For a minus-strand alignment qstart = querylength - querystart and
   qend = querylength - queryend. */

#define T Univdiag_T
struct T {
  Univcoord_T univdiagonal;
  int qstart;
  int qend;
  int nmismatches;
};

#undef T
#endif

