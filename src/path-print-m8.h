/* $Id: 36a2548040e875db9f06a7d8800599205010b237 $ */
#ifndef PATH_PRINT_M8_INCLUDED
#define PATH_PRINT_M8_INCLUDED

#include "path.h"

#include "filestring.h"
#include "listpool.h"

#include "iit-read-univ.h"


extern void
Path_print_m8 (Filestring_T fp, Path_T path, char *accession, char *acc_suffix,
	       bool invertp, Listpool_T listpool);

extern void
Path_print_m8_setup (Univ_IIT_T chromosome_iit_in);

#undef T
#endif

