/* $Id: 0386214e6a8d2747cecd26cd8858eccea5721df6 $ */
#ifndef INTLISTPOOL_INCLUDED
#define INTLISTPOOL_INCLUDED

/* #define INTLISTPOOL_REUSE 1 */
/* #define INTLISTPOOL_TRACE 1 */

#ifdef INTLISTPOOL_TRACE
#define intlistpool_trace(a,b) ,a,b
#else
#define intlistpool_trace(a,b)
#endif

typedef struct Intlistpool_T *Intlistpool_T;

#include "intlist.h"

#define T Intlistpool_T

extern void
Intlistpool_reset_memory (T this);
extern void
Intlistpool_free (T *old);
extern T
Intlistpool_new (void);

#ifdef INTLISTPOOL_REUSE
extern void
Intlistpool_free_list (Intlist_T *old, T this
#ifdef INTLISTPOOL_TRACE
		       , const char *file, int line
#endif
		       );

#else
static inline void
Intlistpool_free_list (Intlist_T *old, T this
#ifdef INTLISTPOOL_TRACE
		       , const char *file, int line
#endif
		       ) {
  (void)(this);
  *old = (Intlist_T) NULL;
  return;
}
#endif


#ifdef INTLISTPOOL_REUSE
extern Intlist_T
Intlistpool_pop (Intlist_T intlist, T this, int *integer
#ifdef INTLISTPOOL_TRACE
		 , const char *file, int line
#endif
		 );

#else
static inline Intlist_T
Intlistpool_pop (Intlist_T intlist, T this, int *integer
#ifdef INTLISTPOOL_TRACE
		 , const char *file, int line
#endif
		 ) {
  (void)(this);
  *integer = intlist->first;
  return intlist->rest;
}
#endif

extern Intlist_T
Intlistpool_push (Intlist_T intlist, T this, int integer
#ifdef INTLISTPOOL_TRACE
		  , const char *file, int line
#endif
		  );

extern Intlist_T
Intlistpool_copy (Intlist_T source, T this);
extern void
Intlistpool_init (T this);

#undef T
#endif
