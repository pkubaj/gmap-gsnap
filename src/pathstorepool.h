/* $Id: b77da64b3fdd3d2e6733712433b2e8332da5ac62 $ */
#ifndef PATHSTOREPOOL_INCLUDED
#define PATHSTOREPOOL_INCLUDED

/* #define PATHSTOREPOOL_REUSE 1 */
/* #define PATHSTOREPOOL_TRACE 1 */

#ifdef PATHSTOREPOOL_TRACE
#define pathstorepool_trace(a,b) ,a,b
#else
#define pathstorepool_trace(a,b)
#endif

typedef struct Pathstorepool_T *Pathstorepool_T;

#include "pathstore.h"
#include "univdiag.h"

#define T Pathstorepool_T

extern void
Pathstorepool_reset_memory (T this);
extern void
Pathstorepool_free (T *old);
extern T
Pathstorepool_new (void);

#ifdef PATHSTOREPOOL_REUSE
extern void
Pathstorepool_free_pathstore (Pathstore_T *old, T this
#ifdef PATHSTOREPOOL_TRACE
			  , const char *file, int line
#endif
);

#else
static inline void
Pathstorepool_free_pathstore (Pathstore_T *old, T this
#ifdef PATHSTOREPOOL_TRACE
			  , const char *file, int line
#endif
		    ) {
  (void)(this);
  *old = (Pathstore_T) NULL;
  return;
}
#endif


extern Pathstore_T
Pathstorepool_new_pathstore (T this
#ifdef PATHSTOREPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern void
Pathstorepool_init (T this);

#undef T
#endif
