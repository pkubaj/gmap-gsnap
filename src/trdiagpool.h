/* $Id: 201ae6670887ea966d5617138477602727b29c40 $ */
#ifndef TRDIAGPOOL_INCLUDED
#define TRDIAGPOOL_INCLUDED

/* #define TRDIAGPOOL_REUSE 1 */
/* #define TRDIAGPOOL_TRACE 1 */

#ifdef TRDIAGPOOL_TRACE
#define trdiagpool_trace(a,b) ,a,b
#else
#define trdiagpool_trace(a,b)
#endif

typedef struct Trdiagpool_T *Trdiagpool_T;

#include "types.h"
#include "trdiag.h"
#include "list.h"

#define T Trdiagpool_T

extern void
Trdiagpool_reset_memory (T this);
extern void
Trdiagpool_free (T *old);
extern T
Trdiagpool_new (void);
extern Trdiag_T
Trdiag_new (T this, int qstart, int qend, int nmismatches, Trcoord_T trdiagonal);

#ifdef TRDIAGPOOL_REUSE
extern void
Trdiagpool_free_trdiag (Trdiag_T *old, T this
#ifdef TRDIAGPOOL_TRACE
			, const char *file, int line
#endif
			);

#else
static inline void
Trdiagpool_free_trdiag (Trdiag_T *old, T this
#ifdef TRDIAGPOOL_TRACE
			, const char *file, int line
#endif
			) {
  (void)(old);
  (void)(this);
  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
extern void
Trdiagpool_free_list (List_T *old, T this
#ifdef TRDIAGPOOL_TRACE
			, const char *file, int line
#endif
			);

#else
static inline void
Trdiagpool_free_list (List_T *old, T this
#ifdef TRDIAGPOOL_TRACE
		      , const char *file, int line
#endif
		      ) {
  (void)(this);
  *old = (List_T) NULL;
  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
extern void
Trdiagpool_gc (List_T *list, T this
#ifdef TRDIAGPOOL_TRACE
	       , const char *file, int line
#endif
	       );

#else
static inline void
Trdiagpool_gc (List_T *list, T this
#ifdef TRDIAGPOOL_TRACE
	       , const char *file, int line
#endif
	       ) {
  (void)(this);
  *list = (List_T) NULL;
  return;
}
#endif


#ifdef TRDIAGPOOL_REUSE
extern List_T
Trdiagpool_pop (List_T list, T this, Trdiag_T *x
#ifdef TRDIAGPOOL_TRACE
		, const char *file, int line
#endif
		);

#else
static inline List_T
Trdiagpool_pop (List_T list, T this, Trdiag_T *x
#ifdef TRDIAGPOOL_TRACE
		, const char *file, int line
#endif
		) {
  (void)(this);
  *x = list->first;
  return list->rest;
}
#endif


extern List_T
Trdiagpool_push (List_T list, T this, int qstart, int qend, int nmismatches, Trcoord_T trdiagonal
#ifdef TRDIAGPOOL_TRACE
		 , const char *file, int line
#endif
		 );

extern List_T
Trdiagpool_push_existing (List_T list, T this, Trdiag_T trdiag
#ifdef TRDIAGPOOL_TRACE
			  , const char *file, int line
#endif
			  );

extern Trdiag_T
Trdiagpool_new_trdiag (T this, int qstart, int qend, int nmismatches, Trcoord_T trdiagonal
#ifdef TRDIAGPOOL_TRACE
		       , const char *file, int line
#endif
		       );

extern void
Trdiagpool_init (T this);

#undef T
#endif
