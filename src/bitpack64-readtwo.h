#ifndef BITPACK64_READTWO_INCLUDED
#define BITPACK64_READTWO_INCLUDED
#include "types.h"

/* For reading differential-coded bitstreams */

extern UINT4
Bitpack64_read_two (UINT4 *end0, Signature_T signature,
		    UINT4 *bitpackptrs, UINT4 *bitpackcomp);

extern UINT8
Bitpack64_read_two_huge (UINT8 *end0, Signature_T signature,
			 UINT4 *bitpackpages, UINT4 *bitpackptrs, UINT4 *bitpackcomp);

#endif
