#ifndef BITPACK64_INCR_INCLUDED
#define BITPACK64_INCR_INCLUDED
#include "types.h"

/* For reading direct-coded bitstreams.  Horizontal format. */
extern UINT4
Bitpack64_incr (Signature_T signature, UINT4 *ptrs, UINT4 *comp);

extern UINT4
Bitpack64_incr_bitpack (Signature_T signature, int packsize, UINT4 *bitpack);

extern void
Bitpack64_add (Signature_T signature, UINT4 *ptrs, UINT4 *comp, UINT4 increment);

extern void
Bitpack64_add_bitpack (Signature_T signature, int packsize, UINT4 *bitpack, UINT4 increment);

extern UINT4 *
Bitpack64_realloc_one (int packsize, UINT4 *bitpack);

extern UINT4 *
Bitpack64_realloc_multiple (int old_packsize, int new_packsize, UINT4 *bitpack);

#endif
