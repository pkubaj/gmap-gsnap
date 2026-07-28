/* $Id: 40706c502de4dd770bd6d60d54f0c09bf4c992ba $ */
#ifndef CHRNUM_INCLUDED
#define CHRNUM_INCLUDED

typedef int Chrnum_T;
typedef unsigned int Trnum_T;

#include "bool.h"
#include "iit-read-univ.h"
#include "types.h"
#include "genomicpos.h"
#include "univcoord.h"

extern char *
Chrnum_to_string (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern char *
Chrnum_to_string_signed (Chrnum_T chrnum, Univ_IIT_T chromosome_iit, bool watsonp);
extern Chrpos_T
Chrnum_length (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern Univcoord_T
Chrnum_offset (Chrnum_T chrnum, Univ_IIT_T chromosome_iit);
extern void
Chrnum_print_position (Univcoord_T position, Univ_IIT_T chromosome_iit);

#endif
