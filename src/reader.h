/* $Id: 3490efc9467595dd1c93ccef1612ebf2174a2d8a $ */
#ifndef READER_INCLUDED
#define READER_INCLUDED

typedef struct Reader_T *Reader_T;

/* Still used by pair.c and chimera.c */
typedef enum {FIVE, THREE, MIDDLE} cDNAEnd_T;

#include <stdio.h>
#include "bool.h"
#include "types.h"

struct Reader_T {
  int oligosize;

  int querystart;
  int queryend;
  int querystart_save;
  int queryend_save;

  char *startinit;
  char *startptr;
  char *endptr;

  char *startbound;		/* Saved for reset */
  char *endbound;
};


#define T Reader_T
typedef struct T *T;

static inline int
Reader_oligosize (T this) {
  return this->oligosize;
}

static inline int
Reader_querystart (T this) {
  return this->querystart;
}

static inline int
Reader_queryend (T this) {
  return this->queryend;
}

/* Same as current querypos + oligosize */
static inline int
Reader_startpos (T this) {
  return (this->startptr - this->startinit);
}

static inline int
Reader_endpos (T this) {
  return (this->endptr - this->startinit);
}

extern void
Reader_reset_start (T this, int querypos);
extern void
Reader_reset_end (T this, int querypos);
extern void
Reader_reset_ends (T this);

extern T
Reader_new (char *sequence, int querystart, int queryend, int oligosize);
extern void
Reader_free (T *old);

#ifndef GSNAP
extern char
Reader_getc (T this, cDNAEnd_T cdnaend, int blocksize);
#endif

extern char
Reader_getc_5 (T this);
extern char
Reader_getc_3 (T this);

extern Oligospace_T
Reader_check (T this, int querypos, int indexsize);


#undef T
#endif
