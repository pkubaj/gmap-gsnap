/* $Id: compress-write.h 225249 2022-11-06 05:11:32Z twu $ */
#ifndef COMPRESS_WRITE_INCLUDED
#define COMPRESS_WRITE_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "types.h"		/* Needed also for HAVE_64_BIT */
#include "univcoord.h"
#include "genomicpos.h"

extern int
Compress_get_char (FILE *sequence_fp, Univcoord_T position, bool uncompressedp);
extern int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Univcoord_T startpos,
		      Univcoord_T endpos, int index1part);
extern int
Compress_update_memory (int nbadchars, Genomecomp_T *genomecomp, char *gbuffer, Univcoord_T startpos,
			Univcoord_T endpos);
extern void
Compress_compress (char **files, int nfiles, bool stdin_p);
extern void
Compress_unshuffle (FILE *out, FILE *in);
extern void
Compress_unshuffle_bits (FILE *high_out, FILE *low_out, FILE *flags_out, FILE *in);
extern Genomecomp_T *
Compress_create_blocks_comp (char *genomicseg, Univcoord_T genomelength);
extern Genomecomp_T *
Compress_create_blocks_bits (Genomecomp_T *genomecomp, Univcoord_T genomelength);

extern void
Compress_cat (FILE *out, char **files, Univcoord_T *genomelengths, int nfiles);

#undef T
#endif

