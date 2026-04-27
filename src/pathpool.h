/* $Id: f8e0abb9052051c5aba0590c0eef27f912933471 $ */
#ifndef PATHPOOL_INCLUDED
#define PATHPOOL_INCLUDED

/* #define PATHPOOL_REUSE 1 */
/* #define PATHPOOL_TRACE 1 */

#ifdef PATHPOOL_TRACE
#define pathpool_trace(a,b) ,a,b
#else
#define pathpool_trace(a,b)
#endif

typedef struct Pathpool_T *Pathpool_T;

#include "path.h"
#include "junction.h"
#include "altsplice.h"

#define T Pathpool_T

extern void
Pathpool_reset_memory (T this);
extern void
Pathpool_free (T *old);
extern T
Pathpool_new (void);

#ifdef PATHPOOL_REUSE
extern void
Pathpool_free_path (Path_T *old, T this
#ifdef PATHPOOL_TRACE
		    , const char *file, int line
#endif
);

#else
static inline void
Pathpool_free_path (Path_T *old, T this
#ifdef PATHPOOL_TRACE
		    , const char *file, int line
#endif
		    ) {
  (void)(this);
  *old = (Path_T) NULL;
  return;
}
#endif

#ifdef PATHPOOL_REUSE
extern void
Pathpool_free_junction (Junction_T *old, T this
#ifdef PATHPOOL_TRACE
			, const char *file, int line
#endif
			);

#else

static inline void
Pathpool_free_junction (Junction_T *old, T this
#ifdef PATHPOOL_TRACE
			, const char *file, int line
#endif
			) {
  (void)(old);
  (void)(this);
  return;
}
#endif

#ifdef PATHPOOL_REUSE
extern void
Pathpool_free_altsplice (Altsplice_T *old, T this
#ifdef PATHPOOL_TRACE
			 , const char *file, int line
#endif
			 );

#else

static inline void
Pathpool_free_altsplice (Altsplice_T *old, T this
#ifdef PATHPOOL_TRACE
			 , const char *file, int line
#endif
			 ) {
  (void)(old);
  (void)(this);
  return;
}
#endif


extern Path_T
Pathpool_new_path (T this
#ifdef PATHPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern Junction_T
Pathpool_new_junction (T this
#ifdef PATHPOOL_TRACE
		       , const char *file, int line
#endif
		       );

extern Altsplice_T
Pathpool_new_altsplice (T this
#ifdef PATHPOOL_TRACE
			, const char *file, int line
#endif
			);

extern char *
Pathpool_new_string (T this, int nchars);

extern void
Pathpool_init (T this);

#undef T
#endif
