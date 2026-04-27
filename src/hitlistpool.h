/* $Id: 43cba175a8091a3b4b22459b44913e73afcb0d37 $ */
#ifndef HITLISTPOOL_INCLUDED
#define HITLISTPOOL_INCLUDED

/* #define HISTLISTPOOL_REUSE 1 */
/* #define HITLISTPOOL_TRACE 1 */

#ifdef HITLISTPOOL_TRACE
#define hitlistpool_trace(a,b) ,a,b
#else
#define hitlistpool_trace(a,b)
#endif

typedef struct Hitlistpool_T *Hitlistpool_T;

#include "list.h"

#define T Hitlistpool_T

extern void
Hitlistpool_reset_memory (T this);
extern void
Hitlistpool_free (T *old);
extern T
Hitlistpool_new (void);

#ifdef HITLISTPOOL_REUSE
extern void
Hitlistpool_free_list (List_T *old, T this
#ifdef HITLISTPOOL_TRACE
		       , const char *file, int line
#endif
		       );

#else
static inline void
Hitlistpool_free_list (List_T *old, T this
#ifdef HITLISTPOOL_TRACE
		       , const char *file, int line
#endif
		       ) {
  (void)(this);
  *old = (List_T) NULL;
  return;
}
#endif


#ifdef HITLISTPOOL_REUSE
extern List_T
Hitlist_pop (List_T list, T this, void **contents
#ifdef HITLISTPOOL_TRACE
	     , const char *file, int line
#endif
	     );

#else
static inline List_T
Hitlist_pop (List_T list, T this, void **contents
#ifdef HITLISTPOOL_TRACE
	     , const char *file, int line
#endif
	     ) {
  (void)(this);
  *contents = list->first;
  return list->rest;
}
#endif


extern List_T
Hitlist_push (List_T list, T this, void *contents
#ifdef HITLISTPOOL_TRACE
	      , const char *file, int line
#endif
	      );

extern List_T
Hitlist_copy (List_T source, T this);

extern void
Hitlistpool_init (T this);

#undef T
#endif
