/* $Id: samheader.h 224935 2022-10-08 19:05:51Z twu $ */
#ifndef SAMHEADER_INCLUDED
#define SAMHEADER_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "iit-read-univ.h"
#include "outputtype.h"
#include "genome.h"
#include "samflags.h"
#include "filestring.h"

extern void
SAM_header_setup (int argc_in, char **argv_in, int optind_in,
		  int nworkers_in, bool orderedp_in,
		  Univ_IIT_T chromosome_iit_in,
#if defined(GEXACT)
#elif defined(GSNAP)
		  Outputtype_T output_type_in,
#else
		  Printtype_T printtype_in, Genome_T global_genome_in,
#endif
		  bool sam_headers_p_in, char *sam_read_group_id_in, char *sam_read_group_name_in,
		  char *sam_read_group_library_in, char *sam_read_group_platform_in);

extern FILE *
SAM_header_fopen (SAM_split_output_type split_output, bool split_simple_p,
		  char *fileroot, bool paired_end_p, bool appendp);

extern Filestring_T
SAM_header_change_HD_tosorted (FILE *input, int headerlen);

extern void
SAM_header_print_all (FILE *output);

extern int
SAM_header_length (int *lastchar, FILE *fp);

#endif

