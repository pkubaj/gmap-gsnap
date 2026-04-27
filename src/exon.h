/* $Id: 77da5faba4296e02d86b358a02240cfd819bc2e7 $ */
#ifndef EXON_INCLUDED
#define EXON_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

typedef struct Exon_T *Exon_T;

#include "list.h"
#include "transcriptpool.h"
#include "listpool.h"
#include "filestring.h"

#define T Exon_T
struct T {
  int exoni;
  char firstchar, lastchar;
};

extern void
Exon_free (T *old, Transcriptpool_T transcriptpool);

extern void
Exon_list_gc (List_T *exons, Listpool_T listpool, Transcriptpool_T transcriptpool);

extern T
Exon_new (char firstchar, int exoni, char lastchar,
	  Transcriptpool_T transcriptpool);

extern List_T
Exon_list_copy (List_T old, Transcriptpool_T transcriptpool, Listpool_T listpool);

extern bool
Exon_list_consecutivep (List_T exons);
extern bool
Exon_list_validp (bool *repairablep, List_T exons);
extern void
Exon_print_list (Filestring_T fp, List_T exons);

extern void
Exon_print_list_stdout (List_T exons);

#undef T
#endif

