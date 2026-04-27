/* $Id: 2cebb8a36675831bde2630c4428ee3ad0bd33930 $ */
#ifndef TRDIAGDEF_INCLUDED
#define TRDIAGDEF_INCLUDED


/* qstart and qend are the genome-normalized coordinates, so qstart
   marks the left coordinate and qend marks the right coordinate.  For
   a plus-strand alignment, qstart = querystart and qend = queryend.
   For a minus-strand alignment qstart = querylength - querystart and
   qend = querylength - queryend. */

#define T Trdiag_T
struct T {
  Trcoord_T trdiagonal;
  int qstart;
  int qend;
  int nmismatches;

  /* int nconsecutive; -- Previously used by extension-search.c */
  /* int intscore; -- Previously used by extension-search */
};

#undef T
#endif

