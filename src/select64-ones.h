#ifndef SELECT64_ONES_INCLUDED
#define SELECT64_ONES_INCLUDED

#include "types.h"

#define T Select64_ones_T
typedef struct T *T;

extern void
Select64_ones_free (T *old);
extern T
Select64_ones_new (uint64_t *const bits, const uint64_t nbits);
extern uint64_t
Select64_ones_select (T this, const uint64_t rank);
extern uint64_t
Select64_ones_next (T this, const uint64_t rank, uint64_t *const next);

#undef T
#endif
