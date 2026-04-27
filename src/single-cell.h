/* $Id: 16093b05f3c2d32b4110a08f518e0b40b75a8ffb $ */
#ifndef SINGLE_CELL_INCLUDED
#define SINGLE_CELL_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bool.h"
#include "filestring.h"
#include "shortread.h"


extern void
Single_cell_compute_priors (char *whitelist_file,
			    char *read_files_command, bool gunzip_p, bool bunzip2_p,
			    char **files, int nfiles);

extern void
Single_cell_print_fields (Filestring_T fp, Shortread_T infoseq);

extern void
Single_cell_setup (int wellpos_in);

extern void
Single_cell_cleanup ();

#endif

