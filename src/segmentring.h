#ifndef SEGMENTRING_INCLUDED
#define SEGMENTRING_INCLUDED

typedef struct Segmentring_T *Segmentring_T;


#include "segmentlist.h"
#include "bool.h"
#include "univcoord.h"
#include "list.h"
#include "univdiagpool.h"

#define T Segmentring_T

extern T
Segmentring_new (Univcoord_T middle_univdiagonal);

extern void
Segmentring_free (T *old);

extern Segmentlist_T
Segmentring_head (T this);

extern Segmentlist_T
Segmentring_tail (T this);

extern Univcoord_T
Segmentring_middle_univdiagonal (T this);

extern int
Segmentring_length (T this);

extern void
Segmentring_print (T this);

extern void
Segmentring_get (int *index, int *qstart, int *qend, Univcoord_T *univdiagonal,
		 int *count, T this, int i);

extern void
Segmentring_put (T this, int i, int index, int qstart, int qend, Univcoord_T univdiagonal,
		 int count);

extern void
Segmentring_addhi (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
		   int count);

extern void
Segmentring_addlo (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
		   int count);

extern void
Segmentring_add (T this, int pos, int index, int qstart, int qend, Univcoord_T univdiagonal,
		 int count);

extern bool
Segmentring_insert_qstart (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
			   int count, int max_overlap);

extern bool
Segmentring_insert_qend (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
			 int count, int max_overlap);

extern List_T
Segmentring_smooth (bool *usedp, T this, Univdiagpool_T univdiagpool);

#undef T
#endif










