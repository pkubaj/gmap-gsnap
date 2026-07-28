/* $Id: ce74bcdf03604b5d622a9941ae86c98ba5387882 $ */
#ifndef PATH_PRINT_ALIGNMENT_INCLUDED
#define PATH_PRINT_ALIGNMENT_INCLUDED

#include "path.h"
#include "pathpair.h"

#include "shortread.h"
#include "filestring.h"
#include "listpool.h"

#include "iit-read-univ.h"


extern void
Path_print_alignment (Filestring_T fp, Path_T path, Pathpair_T pathpair,
		      Shortread_T shortread, bool invertp, Listpool_T listpool);

extern void
Path_print_alignment_setup (Univ_IIT_T chromosome_iit_in, bool method_print_p_in,
			    bool print_univdiagonal_p_in, bool print_nsnpdiffs_p_in);

#undef T
#endif

