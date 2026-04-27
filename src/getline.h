#ifndef GETLINE_INCLUDED
#define GETLINE_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_ZLIB, HAVE_BZLIB  */
#endif

#include <stdio.h>
#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif

extern char *
Getline (FILE *fp);

#ifdef HAVE_ZLIB
extern char *
Getline_gzip (gzFile fp);
#endif

#ifdef HAVE_BZLIB
extern char *
Getline_bzip2 (Bzip2_T fp);
#endif

extern char *
Getline_wlength (int *string_length, FILE *fp);
extern char *
Getline_wlinefeed (int *string_length, FILE *fp);
extern char *
Getline_wendtab (FILE *fp);

#endif


