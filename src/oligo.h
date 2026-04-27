/* $Id: f0b9fbe1ca4cd6b60b02bf52c788ff2c9c1f50e8 $ */
#ifndef OLIGO_INCLUDED
#define OLIGO_INCLUDED

#include "bool.h"
#include "genomicpos.h"
#include "indexdb.h"
#include "reader.h"

typedef enum {INIT, DONE, INVALID, VALID} Oligostate_T;

#define WATSON +1
#define CRICK +2

/* Used for debugging only */
extern char *
Oligo_one_nt (Oligospace_T oligo, int oligosize);

extern Oligostate_T
Oligo_next_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward, 
	      Oligospace_T *revcomp, Reader_T reader, int genestrand);
extern Oligostate_T
Oligo_next_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward, 
	      Oligospace_T *revcomp, Reader_T reader, int genestrand);

extern Oligostate_T
Oligo_skip_5 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip);
extern Oligostate_T
Oligo_skip_3 (Oligostate_T last_state, int *querypos, Oligospace_T *forward,
	      Oligospace_T *revcomp, Reader_T reader, int genestrand, int nskip);

extern void
Oligo_setup (int mode);

#endif
