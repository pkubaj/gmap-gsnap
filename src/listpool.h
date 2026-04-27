/* $Id: e6fba4ed07f91e7a140ff579de8171c4718c7bf7 $ */
#ifndef LISTPOOL_INCLUDED
#define LISTPOOL_INCLUDED

/* #define LISTPOOL_REUSE 1 */
/* #define LISTPOOL_TRACE 1 */

#ifdef LISTPOOL_TRACE
#define listpool_trace(a,b) ,a,b
#else
#define listpool_trace(a,b)
#endif

typedef struct Listpool_T *Listpool_T;

#include "list.h"

#define T Listpool_T

extern void
Listpool_reset_memory (T this);
extern void
Listpool_free (T *old);
extern T
Listpool_new (void);

#ifdef LISTPOOL_REUSE
extern void
Listpool_free_list (List_T *old, T this
#ifdef LISTPOOL_TRACE
		    , const char *file, int line
#endif
		    );

#else
static inline void
Listpool_free_list (List_T *old, T this
#ifdef LISTPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  (void)(this);
  *old = (List_T) NULL;
  return;
}
#endif


#ifdef LISTPOOL_REUSE
extern List_T
Listpool_pop (List_T list, T this, void **contents
#ifdef LISTPOOL_TRACE
	      , const char *file, int line
#endif
	      );

#else
static inline List_T
Listpool_pop (List_T list, T this, void **contents
#ifdef LISTPOOL_TRACE
	      , const char *file, int line
#endif
	      ) {
  (void)(this);
  *contents = list->first;
  return list->rest;
}
#endif


extern List_T
Listpool_push (List_T list, T this, void *contents
#ifdef LISTPOOL_TRACE
	       , const char *file, int line
#endif
	       );

extern List_T
Listpool_copy (List_T source, T this);

extern void
Listpool_init (T this);

#undef T
#endif
