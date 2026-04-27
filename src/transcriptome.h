/* $Id: 5b895138c91ecfce3c9aecff3ef590178cbd44da $ */
#ifndef TRANSCRIPTOME_INCLUDED
#define TRANSCRIPTOME_INCLUDED
#include "bool.h"
#include "types.h"
#include "chrnum.h"

#define T Transcriptome_T
typedef struct T *T;

extern void
Transcriptome_free (T *old);
extern T
Transcriptome_new (char *genomesubdir, char *genome_fileroot, char *transcriptome_fileroot,
		   bool sharedp);
extern Chrnum_T
Transcriptome_chrnum (int *transcript_genestrand, T this, int trnum);
extern int
Transcriptome_nexons (T this, int trnum);
extern int
Transcriptome_exons (int **exonbounds, Chrpos_T **exonstarts, T this, Trnum_T trnum);
extern Trnum_T
Transcriptome_trnum (int *nexons, int **exonbounds, Chrpos_T **exonstarts, T this, int map_index);

#undef T
#endif


