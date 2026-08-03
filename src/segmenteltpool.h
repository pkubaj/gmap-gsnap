#ifndef SEGMENTELTPOOL_INCLUDED
#define SEGMENTELTPOOL_INCLUDED

/* #define SEGMENTELTPOOL_REUSE 1 */
/* #define SEGMENTELTPOOL_TRACE 1 */

#ifdef SEGMENTELTPOOL_TRACE
#define segmenteltpool_trace(a,b) ,a,b
#else
#define segmenteltpool_trace(a,b)
#endif

typedef struct Segmenteltpool_T *Segmenteltpool_T;

#include "segmentelt.h"

#define T Segmenteltpool_T

extern void
Segmenteltpool_reset_memory (T this);
extern void
Segmenteltpool_free (T *old);
extern T
Segmenteltpool_new (void);

#ifdef SEGMENTELTPOOL_REUSE
extern void
Segmenteltpool_free_segmentelt (Segmentelt_T *old, T this
#ifdef SEGMENTELTPOOL_TRACE
			  , const char *file, int line
#endif
);

#else
static inline void
Segmenteltpool_free_segmentelt (Segmentelt_T *old, T this
#ifdef SEGMENTELTPOOL_TRACE
			  , const char *file, int line
#endif
		    ) {
  (void)(this);
  *old = (Segmentelt_T) NULL;
  return;
}
#endif


extern Segmentelt_T
Segmenteltpool_new_segmentelt (T this
#ifdef SEGMENTELTPOOL_TRACE
		   , const char *file, int line
#endif
		   );

extern void
Segmenteltpool_init (T this);

#undef T
#endif
