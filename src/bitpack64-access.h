#ifndef BITPACK64_ACCESS_INCLUDED
#define BITPACK64_ACCESS_INCLUDED
#include "types.h"
#include "bool.h"

/* For reading direct-coded bitstreams */
extern UINT4
Bitpack64_access (Signature_T signature, UINT4 *ptrs, UINT4 *comp);

extern UINT4
Bitpack64_access_bitpack (Signature_T signature, int packsize, UINT4 *bitpack);

extern bool
Bitpack64_access_filledp (Signature_T signature, int packsize, UINT4 *bitpack);

/* Called only by atoiindex and cmetindex */
extern int
Bitpack64_access_new_packsize (Signature_T signature, int old_packsize, UINT4 *bitpack, int increment);

extern void
Bitpack64_extract_bitpack (UINT4 *out, int packsize, UINT4 *bitpack);

#endif
