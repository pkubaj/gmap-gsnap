/* $Id: 901ba20db2ab43bfcf03583a347cbbde9fc4602e $ */
#ifndef UNIVDIAGPOOL_INCLUDED
#define UNIVDIAGPOOL_INCLUDED

/* #define UNIVDIAGPOOL_REUSE 1 */
/* #define UNIVDIAGPOOL_TRACE 1 */

#ifdef UNIVDIAGPOOL_TRACE
#define univdiagpool_trace(a,b) ,a,b
#else
#define univdiagpool_trace(a,b)
#endif

typedef struct Univdiagpool_T *Univdiagpool_T;

#include "univcoord.h"
#include "univdiag.h"
#include "list.h"

#define T Univdiagpool_T

extern void
Univdiagpool_reset_memory (T this);
extern void
Univdiagpool_free (T *old);
extern T
Univdiagpool_new (void);
extern Univdiag_T
Univdiag_new (T this, int qstart, int qend, int nmismatches, Univcoord_T univdiagonal);

#ifdef UNIVDIAGPOOL_REUSE
extern void
Univdiagpool_free_univdiag (Univdiag_T *old, T this
#ifdef UNIVDIAGPOOL_TRACE
			    , const char *file, int line
#endif
			    );

#else
static inline void
Univdiagpool_free_univdiag (Univdiag_T *old, T this
#ifdef UNIVDIAGPOOL_TRACE
			    , const char *file, int line
#endif
			    ) {
  (void)(old);
  (void)(this);
  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
extern void
Univdiagpool_free_list (List_T *old, T this
#ifdef UNIVDIAGPOOL_TRACE
			, const char *file, int line
#endif
			);

#else
static inline void
Univdiagpool_free_list (List_T *old, T this
#ifdef UNIVDIAGPOOL_TRACE
			, const char *file, int line
#endif
			) {
  (void)(this);
  *old = (List_T) NULL;
  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
extern void
Univdiagpool_gc (List_T *list, T this
#ifdef UNIVDIAGPOOL_TRACE
		 , const char *file, int line
#endif
		 );

#else
static inline void
Univdiagpool_gc (List_T *list, T this
#ifdef UNIVDIAGPOOL_TRACE
		 , const char *file, int line
#endif
		 ) {
  (void)(this);
  *list = (List_T) NULL;
  return;
}
#endif


#ifdef UNIVDIAGPOOL_REUSE
extern List_T
Univdiagpool_pop (List_T list, T this, Univdiag_T *x
#ifdef UNIVDIAGPOOL_TRACE
		  , const char *file, int line
#endif
		  );

#else
static inline List_T
Univdiagpool_pop (List_T list, T this, Univdiag_T *x
#ifdef UNIVDIAGPOOL_TRACE
		  , const char *file, int line
#endif
		  ) {
  (void)(this);
  *x = list->first;
  return list->rest;
}
#endif


extern List_T
Univdiagpool_push (List_T list, T this, int qstart, int qend, int nmismatches, Univcoord_T univdiagonal
#ifdef UNIVDIAGPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern List_T
Univdiagpool_push_existing (List_T list, T this, Univdiag_T univdiag
#ifdef UNIVDIAGPOOL_TRACE
			    , const char *file, int line
#endif
			    );

extern Univdiag_T
Univdiagpool_new_univdiag (T this, int qstart, int qend, int nmismatches, Univcoord_T univdiagonal
#ifdef UNIVDIAGPOOL_TRACE
			   , const char *file, int line
#endif
			   );

extern void
Univdiagpool_init (T this);

#undef T
#endif
