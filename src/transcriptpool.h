/* $Id: be835bf757de46fbbfcf89d23a112397b2b173b1 $ */
#ifndef TRANSCRIPTPOOL_INCLUDED
#define TRANSCRIPTPOOL_INCLUDED

/* #define TRANSCRIPTPOOL_REUSE 1 */
/* #define TRANSCRIPTPOOL_TRACE 1 */

#ifdef TRANSCRIPTPOOL_TRACE
#define transcriptpool_trace(a,b) ,a,b
#else
#define transcriptpool_trace(a,b)
#endif

typedef struct Transcriptpool_T *Transcriptpool_T;

#include "transcript.h"
#include "exon.h"

#define T Transcriptpool_T

extern void
Transcriptpool_reset_memory (T this);
extern void
Transcriptpool_free (T *old);
extern T
Transcriptpool_new (void);

#ifdef TRANSCRIPTPOOL_REUSE
extern void
Transcriptpool_free_transcript (Transcript_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
				, const char *file, int line
#endif
				);

#else
static inline void
Transcriptpool_free_transcript (Transcript_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
				, const char *file, int line
#endif
				) {
  (void)(this);
  *old = (Transcript_T) NULL;
  return;
}
#endif

#ifdef TRANSCRIPTPOOL_REUSE
extern void
Transcriptpool_free_exon (Exon_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
			  , const char *file, int line
#endif
			  );

#else
static inline void
Transcriptpool_free_exon (Exon_T *old, T this
#ifdef TRANSCRIPTPOOL_TRACE
			 , const char *file, int line
#endif
			  ) {
  (void)(this);
  *old = (Exon_T) NULL;
  return;
}
#endif


extern Transcript_T
Transcriptpool_new_transcript (T this
#ifdef TRANSCRIPTPOOL_TRACE
			       , const char *file, int line
#endif
			       );

extern Exon_T
Transcriptpool_new_exon (T this
#ifdef TRANSCRIPTPOOL_TRACE
			 , const char *file, int line
#endif
			 );

extern void
Transcriptpool_init (T this);

#undef T
#endif
