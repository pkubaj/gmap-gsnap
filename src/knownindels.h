/* $Id: knownindels.h 225283 2022-11-07 05:51:06Z twu $ */
#ifndef KNOWNINDELS_INCLUDED
#define KNOWNINDELS_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "univcoord.h"
#include "iit-read.h"
#include "iit-read-univ.h"
#include "ef64.h"


#define T Knownindels_T
typedef struct T *T;

extern void
Knownindels_free (T *old);

extern T
Knownindels_new (Univcoordtable_T indel_table, Univcoord_T genomelength,
		 FILE *dump_indels_fp, Univ_IIT_T chromosome_iit, EF64_T chromosome_ef64);

extern T
Knownindels_new_from_dump (FILE *fp, Univcoord_T genomelength);

extern int
Knownindels_find_lowest (int *indel_pos, T this,
			 Univcoord_T univdiagonal, int querylength,
			 int pos5, int pos3);

extern int
Knownindels_find_highest (int *indel_pos, T this,
			  Univcoord_T univdiagonal, int querylength,
			  int pos5, int pos3);

#undef T
#endif
