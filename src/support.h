/* $Id$ */
#ifndef SUPPORT_INCLUDED
#define SUPPORT_INCLUDED

typedef struct Support_T *Support_T;

#include <stdio.h>

#include "bool.h"
#include "path.h"
#include "pathpair.h"

#include "ef64.h"
#include "iit-read-univ.h"
#include "univcoord.h"
#include "knownsplicing.h"
#include "genomebits.h"

#define T Support_T

extern T
Support_new (Knownsplicing_T knownsplicing);

extern void
Support_free (T *old);

extern void
Support_revise_single (T this, Path_T *patharray, int npaths);

extern void
Support_revise_paired (T this, Pathpair_T *pathpairarray, int npaths);

extern bool
Support_check (T this, Path_T path);

extern void
Support_reset_counts (T this);

extern void
Support_filter (T this, EF64_T chromosome_ef64, Univcoord_T genomelength, int threshold);

extern void
Support_convert (T this, Univcoord_T genomelength);

extern T
Support_read (FILE *fp, EF64_T chromosome_ef64, Univ_IIT_T chromosome_iit);

extern bool
Support_presentp (T this,
		  Univcoord_T splice_position_i, Univcoord_T splice_position_j,
		  bool plusp, bool sense_forward_p);

extern void
Support_dump (FILE *fp, T this, EF64_T chromosome_ef64, Univ_IIT_T chromosome_iit,
	      Genomebits_T ref, bool include_known_p);

#undef T
#endif

