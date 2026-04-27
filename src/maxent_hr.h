/* $Id: maxent_hr.h 223511 2020-11-14 15:50:08Z twu $ */
#ifndef MAXENT_HR_INCLUDED
#define MAXENT_HR_INCLUDED

#include "genomicpos.h"
#include "types.h"
#include "univcoord.h"
#include "genome.h"


#ifdef GSNAP
extern void
Maxent_hr_setup (Genome_T genome_in, Genome_T genomealt_in);
#endif

extern double
Maxent_hr_donor_prob (
#ifndef GSNAP
		      Genome_T genome, Genome_T genomealt,
#endif
		      Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_acceptor_prob (
#ifndef GSNAP
			 Genome_T genome, Genome_T genomealt,
#endif
			 Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antidonor_prob (
#ifndef GSNAP
			  Genome_T genome, Genome_T genomealt,
#endif
			  Univcoord_T splice_pos, Univcoord_T chroffset);

extern double
Maxent_hr_antiacceptor_prob (
#ifndef GSNAP
			     Genome_T genome, Genome_T genomealt,
#endif
			     Univcoord_T splice_pos, Univcoord_T chroffset);

#endif

