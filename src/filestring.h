/* $Id: filestring.h 224951 2022-10-08 20:03:19Z twu $ */
#ifndef FILESTRING_INCLUDED
#define FILESTRING_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "samflags.h"

#define FPRINTF Filestring_put
#define PUTC Filestring_putc


#define T Filestring_T
typedef struct T *T;

extern void
Filestring_setup (bool split_simple_p_in);
extern unsigned int
Filestring_id (T this);
extern void
Filestring_set_split_output (T this, bool concordant_softclipped_p, int split_output);
extern SAM_split_output_type
Filestring_split_output (T this);
extern T
Filestring_new ();
extern void
Filestring_free (T *old, bool free_string_p);
extern void
Filestring_stringify (T this);
extern char *
Filestring_string (T this);
extern void
Filestring_print (FILE *fp, T this);
extern char *
Filestring_get (int *strlength, T this);
extern void
Filestring_put (T this, const char *format, ...);
extern void
Filestring_putc (char c, T this);
extern void
Filestring_puts (T this, char *string, int strlength);
extern void
Filestring_merge (T dest, T source);


#undef T
#endif


