#ifndef EF64_INCLUDED
#define EF64_INCLUDED

#include "bool.h"
#include "types.h"
#include "uintlist.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "iit-read-univ.h"

#define T EF64_T
typedef struct T *T;

extern void
EF64_free (T *old);
extern T
EF64_new_from_oligos (Oligospace_T *repetitive_oligos, uint64_t nrepetitive, int oligosize);
extern T
EF64_new_from_univcoordlist (Univcoordlist_T positions, const uint64_t nbits);
extern T
EF64_new_from_chromosome_iit (Univ_IIT_T chromosome_iit);
extern T
EF64_new_from_interleaved_univcoords (Univcoord_T *endpoints, uint64_t nintervals, const uint64_t nbits);
extern uint64_t
EF64_rank (T this, const size_t k);
extern void
EF64_two_ranks (uint64_t *rank1, uint64_t *rank2, T this, const size_t k1, const size_t k2);
extern size_t
EF64_select (T this, const uint64_t rank);
extern uint64_t
EF64_next (T this, const uint64_t rank, size_t *next);
extern Trnum_T
EF64_trnum (Trcoord_T *troffset, Trcoord_T *trhigh, T this,
	    Trcoord_T low_trcoord, Trcoord_T high_trcoord);
extern Chrnum_T
EF64_chrnum (Univcoord_T *chroffset, Univcoord_T *chrhigh, T this,
	     Univcoord_T low_univcoord, Univcoord_T high_univcoord);
extern void
EF64_chrbounds (Univcoord_T *chroffset, Univcoord_T *chrhigh, T this, Chrnum_T chrnum);
extern bool
EF64_presentp (size_t value, T this);

#undef T
#endif

