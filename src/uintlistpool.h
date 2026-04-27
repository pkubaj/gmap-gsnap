/* $Id: 400afafd1a18356d55df2b58542abe8d36539fa9 $ */
#ifndef UINTLISTPOOL_INCLUDED
#define UINTLISTPOOL_INCLUDED

/* #define UINTLISTPOOL_REUSE 1 */
/* #define UINTLISTPOOL_TRACE 1 */

#ifdef UINTLISTPOOL_TRACE
#define uintlistpool_trace(a,b) ,a,b
#else
#define uintlistpool_trace(a,b)
#endif

typedef struct Uintlistpool_T *Uintlistpool_T;

#include "uintlist.h"

#define T Uintlistpool_T

extern void
Uintlistpool_reset_memory (T this);
extern void
Uintlistpool_free (T *old);
extern T
Uintlistpool_new (void);

#ifdef UINTLISTPOOL_REUSE
extern void
Uintlistpool_free_list (Uintlist_T *old, T this
#ifdef UINTLISTPOOL_TRACE
			, const char *file, int line
#endif
			);

#else
static inline void
Uintlistpool_free_list (Uintlist_T *old, T this
#ifdef UINTLISTPOOL_TRACE
			, const char *file, int line
#endif
			) {
  (void)(this);
  *old = (Uintlist_T) NULL;
  return;
}
#endif


#ifdef UINTLISTPOOL_REUSE
extern Uintlist_T
Uintlistpool_pop (Uintlist_T list, T this, unsigned int *integer
#ifdef UINTLISTPOOL_TRACE
		  , const char *file, int line
#endif
		  );

#else
static inline Uintlist_T
Uintlistpool_pop (Uintlist_T list, T this, unsigned int *integer
#ifdef UINTLISTPOOL_TRACE
		  , const char *file, int line
#endif
		  ) {
  (void)(this);
  *integer = list->first;
  return list->rest;
}
#endif


extern Uintlist_T
Uintlistpool_push (Uintlist_T list, T this, unsigned int integer
#ifdef UINTLISTPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern Uintlist_T
Uintlistpool_copy (Uintlist_T source, T this);

extern void
Uintlistpool_init (T this);

#undef T
#endif
