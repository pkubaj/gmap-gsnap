#ifndef SELECT64_ZEROES_INCLUDED
#define SELECT64_ZEROES_INCLUDED

#include "types.h"

#define T Select64_zeroes_T
typedef struct T *T;

extern void
Select64_zeroes_free (T *old);
extern T
Select64_zeroes_new (uint64_t *const bits, const uint64_t nbits);
extern uint64_t
Select64_zeroes_select (T this, const uint64_t rank);
extern uint64_t
Select64_zeroes_next (T this, const uint64_t rank, uint64_t *const next);

#undef T
#endif
