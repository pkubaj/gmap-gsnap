#ifndef MD5_INCLUDED
#define MD5_INCLUDED

#include <stdio.h>
#include "filestring.h"

extern unsigned char *
MD5_compute (unsigned char *input, int input_len);
extern void
MD5_print (Filestring_T fp, unsigned char *digest);

#endif

