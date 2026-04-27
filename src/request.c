static char rcsid[] = "$Id: request.c 224936 2022-10-08 19:06:28Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "request.h"
#include "assert.h"
#include "mem.h"

#if defined(GFILTER)
#elif defined(GEXACT)
#elif defined(GSNAP)
#else
static Genome_T global_genome;
static Genome_T global_genomealt;
#endif


#define T Request_T
struct T {
  unsigned int id;
#if defined(GFILTER)
  Shortread_T queryseq1;
  Shortread_T queryseq2;
#elif defined(GEXACT)
  Shortread_T queryseq1;
#elif defined(GSNAP)
  Shortread_T queryseq1;
  Shortread_T queryseq2;
#else
  Genome_T genome;
  Genome_T genomealt;
  Sequence_T queryseq;

  /* In --crossalign mode, we cannot free genome or sequence until all
     requests for that genome and sequence are completed */
  bool free_genome_p;
#endif  
};


unsigned int
Request_id (T this) {
  return this->id;
}

#if defined(GEXACT)

Shortread_T
Request_queryseq1 (T this) {
  return this->queryseq1;
}

T
Request_new (unsigned int id, Shortread_T queryseq1) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->id = id;
  new->queryseq1 = queryseq1;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Shortread_free(&(*old)->queryseq1);
    FREE_IN(*old);
  }
  return;
}

#elif defined(GSNAP) || defined(GFILTER)

Shortread_T
Request_queryseq1 (T this) {
  return this->queryseq1;
}

Shortread_T
Request_queryseq2 (T this) {
  return this->queryseq2;
}

T
Request_new (unsigned int id, Shortread_T queryseq1, Shortread_T queryseq2) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->id = id;
  new->queryseq1 = queryseq1;
  new->queryseq2 = queryseq2;
  return new;
}

void
Request_free (T *old) {
  if (*old) {
    Shortread_free(&(*old)->queryseq1);
    if ((*old)->queryseq2) {
      Shortread_free(&(*old)->queryseq2);
    }
    FREE_IN(*old);
  }
  return;
}


#else  /* GMAP */

Genome_T
Request_genome (T this) {
  if (this->genome == NULL) {
    return global_genome;
  } else {
    return this->genome;
  }
}

Genome_T
Request_genomealt (T this) {
  if (this->genomealt == NULL) {
    return global_genomealt;
  } else {
    return this->genomealt;
  }
}

Sequence_T
Request_queryseq (T this) {
  return this->queryseq;
}

T
Request_new (unsigned int id, Genome_T genome, Genome_T genomealt, Sequence_T queryseq,
	     bool free_genome_p) {
  T new = (T) MALLOC_IN(sizeof(*new));

  new->id = id;
  new->genome = genome;
  new->genomealt = genomealt;
  new->queryseq = queryseq;
  new->free_genome_p = free_genome_p;

  Sequence_incr_nrequests(queryseq);
  
  return new;
}

void
Request_free (T *old) {
  Sequence_T queryseq;

  if (*old) {
    if ((*old)->free_genome_p == true) {
      if ((*old)->genomealt != (*old)->genome) {
	Genome_free(&(*old)->genomealt);
      }
      if ((*old)->genome != NULL) {
	Genome_free(&(*old)->genome);
      }
    }

    queryseq = (*old)->queryseq;
    if (Sequence_decr_nrequests(queryseq) == 0) {
      Sequence_free(&queryseq);
    }
    FREE_IN(*old);
  }
  return;
}

void
Request_setup (Genome_T global_genome_in, Genome_T global_genomealt_in) {
  global_genome = global_genome_in;
  global_genomealt = global_genomealt_in;
  return;
}


#endif	/* not GSNAP */
