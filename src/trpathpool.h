/* $Id: 29de6bafdddb6baffd9707ed9ef77fb77b85d160 $ */
#ifndef TRPATHPOOL_INCLUDED
#define TRPATHPOOL_INCLUDED

/* #define TRPATHPOOL_REUSE 1 */
/* #define TRPATHPOOL_TRACE 1 */

#ifdef TRPATHPOOL_TRACE
#define trpathpool_trace(a,b) ,a,b
#else
#define trpathpool_trace(a,b)
#endif

typedef struct Trpathpool_T *Trpathpool_T;

#include "trpath.h"

#define T Trpathpool_T

extern void
Trpathpool_reset_memory (T this);
extern void
Trpathpool_free (T *old);
extern T
Trpathpool_new (void);

#ifdef TRPATHPOOL_REUSE
extern void
Trpathpool_free_trpath (Trpath_T *old, T this
#ifdef TRPATHPOOL_TRACE
			, const char *file, int line
#endif
			);

#else
static inline void
Trpathpool_free_trpath (Trpath_T *old, T this
#ifdef TRPATHPOOL_TRACE
			, const char *file, int line
#endif
			) {
  (void)(old);
  (void)(this);
  return;
}
#endif


extern Trpath_T
Trpathpool_new_trpath (T this
#ifdef TRPATHPOOL_TRACE
		       , const char *file, int line
#endif
		       );

extern void
Trpathpool_init (T this);

#undef T
#endif
