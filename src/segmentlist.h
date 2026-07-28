#ifndef SEGMENTLIST_INCLUDED
#define SEGMENTLIST_INCLUDED


typedef struct Segmentlist_T *Segmentlist_T;

#include "univcoord.h"
struct Segmentlist_T {
  Segmentlist_T llink;
  Segmentlist_T rlink;

  int index;			/* in the list of exhaustive univdiagonals */
  int qstart;
  int qend;
  Univcoord_T univdiagonal;
  int count;

  /* For dynamic programming */
  int score;
  Segmentlist_T ptr;
};


#define T Segmentlist_T

extern void
Segmentlist_print (T this);

extern T
Segmentlist_prev (T this);

extern T
Segmentlist_next (T this);

extern int
Segmentlist_qstart (T this);

extern int
Segmentlist_qend (T this);

extern Univcoord_T
Segmentlist_univdiagonal (T this);

extern int
Segmentlist_count (T this);


#undef T
#endif


