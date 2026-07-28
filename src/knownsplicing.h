/* $Id: a7bae684e735a2adca0f7bceb421dbba8111fb59 $ */
#ifndef KNOWNSPLICING_INCLUDED
#define KNOWNSPLICING_INCLUDED

#include <stdio.h>
#include "bool.h"
#include "genomicpos.h"
#include "univcoord.h"
#include "uintlist.h"
#include "compress.h"
#include "genomebits.h"
#include "iit-read.h"
#include "iit-read-univ.h"
#include "transcriptome.h"
#include "ef64.h"


#define T Knownsplicing_T
typedef struct T *T;
struct T {
  EF64_T donor_ef;
  Univcoord_T *donor_endpoints;
  int donor_nintervals;

  EF64_T acceptor_ef;
  Univcoord_T *acceptor_endpoints;
  int acceptor_nintervals;

  EF64_T antidonor_ef;
  Univcoord_T *antidonor_endpoints;
  int antidonor_nintervals;

  EF64_T antiacceptor_ef;
  Univcoord_T *antiacceptor_endpoints;
  int antiacceptor_nintervals;
};

extern void
Knownsplicing_free (T *old);

extern T
Knownsplicing_new (Univcoordlist_T donor_startpoints, Univcoordlist_T donor_partners,
		   Univcoordlist_T acceptor_startpoints, Univcoordlist_T acceptor_partners,
		   Univcoordlist_T antidonor_startpoints, Univcoordlist_T antidonor_partners,
		   Univcoordlist_T antiacceptor_startpoints, Univcoordlist_T antiacceptor_partners,
		   Univcoordtableuint_T donor_table, Univcoordtableuint_T acceptor_table, 
		   Univcoordtableuint_T antidonor_table, Univcoordtableuint_T antiacceptor_table, 
		   Univcoord_T genomelength, Univ_IIT_T chromosome_iit,
		   EF64_T chromosome_ef64, bool intron_level_p);

extern T
Knownsplicing_new_from_dump (FILE *fp, Univcoord_T genomelength);

extern int
Knownsplicing_nintervals (T this);

extern Univcoord_T *
Knownsplicing_donors (uint64_t *low_rank, uint64_t *high_rank, T this,
		      Univcoord_T univdiagonal, int querylength, int pos5, int pos3);
extern Univcoord_T *
Knownsplicing_acceptors (uint64_t *low_rank, uint64_t *high_rank, T this,
			 Univcoord_T univdiagonal, int querylength, int pos5, int pos3);
extern Univcoord_T *
Knownsplicing_antidonors (uint64_t *low_rank, uint64_t *high_rank, T this,
			  Univcoord_T univdiagonal, int querylength, int pos5, int pos3);
extern Univcoord_T *
Knownsplicing_antiacceptors (uint64_t *low_rank, uint64_t *high_rank, T this,
			     Univcoord_T univdiagonal, int querylength, int pos5, int pos3);

extern bool
Knownsplicing_intron_p (T this, Univcoord_T donor_position, Univcoord_T acceptor_position);

extern bool
Knownsplicing_antiintron_p (T this, Univcoord_T acceptor_position, Univcoord_T donor_position);

extern T
Knownsplicing_from_splicing_iit (IIT_T splicing_iit, int *splicing_divint_crosstable,
				 int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
				 bool intron_level_p);
extern T
Knownsplicing_from_transcriptome (Transcriptome_T transcriptome, int nalignments,
				  EF64_T chromosome_ef64, Univcoord_T genomelength, bool intron_level_p);

#undef T
#endif
