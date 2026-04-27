/* $Id: 2ef7587de59aa67bed1dc86675c02ee1cb8fad0e $ */
#ifndef INBUFFER_INCLUDED
#define INBUFFER_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_ZLIB, HAVE_BZLIB */
#endif

#include <stdio.h>
#include "bool.h"
#include "outbuffer.h"
#include "request.h"

#if defined(GFILTER) || defined(GEXACT) || defined(GSNAP)
#include "shortread.h"
#else
#include "sequence.h"
#include "genome.h"
#include "list.h"
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#define T Inbuffer_T
typedef struct T *T;

extern void
Inbuffer_setup (bool single_cell_p_in, bool filter_if_both_p_in);

#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
extern T
Inbuffer_cmdline (char *queryseq, int querylength, Genome_T genome, Genome_T genomealt);
#endif

extern T
Inbuffer_new (int nextchar, FILE *input,

#if defined(GEXACT)
#ifdef HAVE_ZLIB
	      gzFile gzipped,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped,
#endif

#elif defined(GSNAP) || defined(GFILTER)
	      FILE *input2,
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
	      bool interleavedp,
#endif
	      char *read_files_command, char **files, int nfiles, unsigned int nspaces,
#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
	      bool user_pairalign_p, List_T user_genomes,
#endif
	      unsigned int part_modulus, unsigned int part_interval);

extern void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer);

extern void
Inbuffer_free (T *old);

extern unsigned int
Inbuffer_fill_init (T this);

extern Request_T
Inbuffer_get_request (T this);

#undef T
#endif

