/* $Id: request.h 224936 2022-10-08 19:06:28Z twu $ */
#ifndef REQUEST_INCLUDED
#define REQUEST_INCLUDED

#if defined(GFILTER)
#include "shortread.h"
#elif defined(GEXACT)
#include "shortread.h"
#elif defined(GSNAP)
#include "shortread.h"
#else
#include "bool.h"
#include "sequence.h"
#include "genome.h"
#endif

#define T Request_T
typedef struct T *T;

extern unsigned int
Request_id (T this);

#if defined(GEXACT)

extern Shortread_T
Request_queryseq1 (T this);
extern T
Request_new (unsigned int id, Shortread_T queryseq1);

#elif defined(GSNAP) || defined(GFILTER)

extern Shortread_T
Request_queryseq1 (T this);
extern Shortread_T
Request_queryseq2 (T this);
extern T
Request_new (unsigned int id, Shortread_T queryseq1, Shortread_T queryseq2);

#else

extern Genome_T
Request_genome (T this);
extern Genome_T
Request_genomealt (T this);
extern Sequence_T
Request_queryseq (T this);
extern T
Request_new (unsigned int id, Genome_T genome, Genome_T genomealt, Sequence_T queryseq,
	     bool free_genome_p);
extern void
Request_setup (Genome_T global_genome_in, Genome_T global_genomealt_in);

#endif


extern void
Request_free (T *old);

#undef T
#endif
