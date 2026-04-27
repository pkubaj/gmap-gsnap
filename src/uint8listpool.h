/* $Id: 2057fbcc5a9dcc076de8e040d8475e05bc7a2d7e $ */
#ifndef UINT8LISTPOOL_INCLUDED
#define UINT8LISTPOOL_INCLUDED

/* #define UINT8LISTPOOL_REUSE 1 */
/* #define UINT8LISTPOOL_TRACE 1 */

#ifdef UINT8LISTPOOL_TRACE
#define uint8listpool_trace(a,b) ,a,b
#else
#define uint8listpool_trace(a,b)
#endif

typedef struct Uint8listpool_T *Uint8listpool_T;

#include "uint8list.h"

#define T Uint8listpool_T

extern void
Uint8listpool_reset_memory (T this);
extern void
Uint8listpool_free (T *old);
extern T
Uint8listpool_new (void);

#ifdef UINT8LISTPOOL_REUSE
extern void
Uint8listpool_free_list (Uint8list_T *old, T this
#ifdef UINT8LISTPOOL_TRACE
			 , const char *file, int line
#endif
			 );

#else
static inline void
Uint8listpool_free_list (Uint8list_T *old, T this
#ifdef UINT8LISTPOOL_TRACE
			 , const char *file, int line
#endif
			 ) {
  (void)(this);
  *old = (Uint8list_T) NULL;
  return;
}
#endif

#ifdef UINT8LISTPOOL_REUSE
extern Uint8list_T
Uint8listpool_pop (Uint8list_T list, T this, UINT8 *integer
#ifdef UINT8LISTPOOL_TRACE
		   , const char *file, int line
#endif
		   );

#else
static inline Uint8list_T
Uint8listpool_pop (Uint8list_T list, T this, UINT8 *integer
#ifdef UINT8LISTPOOL_TRACE
		   , const char *file, int line
#endif
		   ) {
  (void)(this);
  *integer = list->first;
  return list->rest;
}
#endif

extern Uint8list_T
Uint8listpool_push (Uint8list_T list, T this, UINT8 integer
#ifdef UINT8LISTPOOL_TRACE
		    , const char *file, int line
#endif
		    );

extern Uint8list_T
Uint8listpool_copy (Uint8list_T source, T this);

extern void
Uint8listpool_init (T this);

#undef T
#endif
