/* $Id: f8e0abb9052051c5aba0590c0eef27f912933471 $ */
#ifndef AUXINFOPOOL_INCLUDED
#define AUXINFOPOOL_INCLUDED

/* #define AUXINFOPOOL_REUSE 1 */
/* #define AUXINFOPOOL_TRACE 1 */

#ifdef AUXINFOPOOL_TRACE
#define auxinfopool_trace(a,b) ,a,b
#else
#define auxinfopool_trace(a,b)
#endif

typedef struct Auxinfopool_T *Auxinfopool_T;

#include "auxinfo.h"
#include "univdiag.h"

#define T Auxinfopool_T

extern void
Auxinfopool_reset_memory (T this);
extern void
Auxinfopool_free (T *old);
extern T
Auxinfopool_new (void);

#ifdef AUXINFOPOOL_REUSE
extern void
Auxinfopool_free_auxinfo (Auxinfo_T *old, T this
#ifdef AUXINFOPOOL_TRACE
			  , const char *file, int line
#endif
);

#else
static inline void
Auxinfopool_free_auxinfo (Auxinfo_T *old, T this
#ifdef AUXINFOPOOL_TRACE
			  , const char *file, int line
#endif
		    ) {
  (void)(this);
  *old = (Auxinfo_T) NULL;
  return;
}
#endif


extern Auxinfo_T
Auxinfopool_new_auxinfo (T this
#ifdef AUXINFOPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern void
Auxinfopool_init (T this);

#undef T
#endif
