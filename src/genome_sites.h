/* $Id: 65d0e2add68a59988f4d2e8f599aa3b7ab7b167d $ */
#ifndef GENOME_SITES_INCLUDED
#define GENOME_SITES_INCLUDED

#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "genomicpos.h"
#include "genome.h"
#include "intlistpool.h"

#define T Genome_T

extern void
Genome_sites_setup (T genome_in, T genomealt_in);

extern int
Genome_donor_sites (int *sites, int *types, int *old_knownpos, int *old_knowni,
		    Univcoord_T left, int pos5, int pos3, Intlistpool_T intlistpool);
extern int
Genome_acceptor_sites (int *sites, int *types, int *old_knownpos, int *old_knowni,
		       Univcoord_T left, int pos5, int pos3, Intlistpool_T intlistpool);

extern int
Genome_antidonor_sites (int *sites, int *types, int *old_knownpos, int *old_knowni,
			Univcoord_T left, int pos5, int pos3, Intlistpool_T intlistpool);

extern int
Genome_antiacceptor_sites (int *sites, int *types, int *old_knownpos, int *old_knowni,
			   Univcoord_T left, int pos5, int pos3, Intlistpool_T intlistpool);


extern void
Genome_fill_donor_gtgc_sites (bool *sitep, Univcoord_T left, int pos5, int pos3);
extern void
Genome_fill_acceptor_ag_sites (bool *sitep, Univcoord_T left, int pos5, int pos3);
extern void
Genome_fill_antidonor_gtgc_sites (bool *sitep, Univcoord_T left, int pos5, int pos3);
extern void
Genome_fill_antiacceptor_ag_sites (bool *sitep, Univcoord_T left, int pos5, int pos3);


#undef T
#endif

