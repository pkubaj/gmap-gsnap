/* $Id: f4c6fc6d860ec2a93c99c1b6aed7087313b952bf $ */
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

