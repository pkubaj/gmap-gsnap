#ifndef SIGNATURE_INCLUDED
#define SIGNATURE_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "types.h"

extern void
Signature_setup (Width_T ksize_in, Width_T msize_in);

extern Signature_T
Signature_get (Oligospace_T kmer);

extern Signature_T
Signature_revise_right (Oligospace_T kmer, Signature_T prev_signature);

extern Signature_T
Signature_revise_left (Oligospace_T kmer, Signature_T prev_signature);

#endif
