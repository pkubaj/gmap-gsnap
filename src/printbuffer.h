/* $Id: fcc41c4373bc0cf8ef976a0a5a636b817cb59e9f $ */
#ifndef PRINTBUFFER_INCLUDED
#define PRINTBUFFER_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "samflags.h"

#define T Printbuffer_T
typedef struct T *T;

extern unsigned int
Printbuffer_nlines (T this);

#if defined(GFILTER)
extern T
Printbuffer_new (char *output_root, bool paired_end_p, bool appendp);

#else
extern T
Printbuffer_new (bool split_simple_p, char *split_output_root, char *output_file);
#endif

extern void
Printbuffer_free (T *old);


#if defined(GFILTER)
void
Printbuffer_store (T this, char *string1, char *string2);
#else
extern void
Printbuffer_store (T this, SAM_split_output_type split_output, char *string,
		   char *string_failedinput
#ifdef GSNAP
		   , char *string_failedinput_1, char *string_failedinput_2
#endif
		   );
#endif

#if defined(GFILTER)
extern void
Printbuffer_print (T this);

#else
extern void
Printbuffer_print (T this, FILE **outputs, FILE *output_failedinput,
#ifdef GSNAP
		   FILE *output_failedinput_1, FILE *output_failedinput_2,
#endif
		   bool paired_end_p, bool appendp);
#endif

#undef T
#endif

