/* $Id$ */
#ifndef GENOMEBITS_TRIM_INCLUDED
#define GENOMEBITS_TRIM_INCLUDED

#include "genomebits.h"
#include "bool.h"
#include "mode.h"
#include "compress.h"
#include "univcoord.h"


#define T Genomebits_T

extern int
Genomebits_trim_qend (int *nmismatches_to_trimpos,
		      Compress_T query_compress, T ref,
		      Univcoord_T univdiagonal, int querylength,
		      int pos5, int pos3, bool plusp, int genestrand);

extern int
Genomebits_trim_qstart (int *nmismatches_to_trimpos,
			Compress_T query_compress, T ref,
			Univcoord_T univdiagonal, int querylength,
			int pos5, int pos3, bool plusp, int genestrand);

extern void
Genomebits_trim_setup (Mode_T mode, bool maskedp);

#undef T
#endif



