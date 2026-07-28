#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "segmentlist.h"
#include <stdio.h>
#include <stdlib.h>


#define T Segmentlist_T

void
Segmentlist_print (T this) {
  printf("%d..%d at %u with count %d\n",
	 this->qstart,this->qend,this->univdiagonal,this->count);
  return;
}

T
Segmentlist_prev (T this) {
  return this->llink;
}

T
Segmentlist_next (T this) {
  return this->rlink;
}

int
Segmentlist_qstart (T this) {
  return this->qstart;
}

int
Segmentlist_qend (T this) {
  return this->qend;
}

Univcoord_T
Segmentlist_univdiagonal (T this) {
  return this->univdiagonal;
}

int
Segmentlist_count (T this) {
  return this->count;
}

