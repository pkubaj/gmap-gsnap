/* $Id: 0328df3dcc63b39fa649824ba8d0873cd57ff87a $ */
#ifndef GENOME_CANONICAL_INCLUDED
#define GENOME_CANONICAL_INCLUDED

#include "bool.h"
#include "univcoord.h"
#include "genome.h"

#define T Genome_T

extern bool
Genome_sense_canonicalp (T genome, T genomealt,
			 Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			 Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			 Univcoord_T chroffset);

extern bool
Genome_antisense_canonicalp (T genome, T genomealt,
			     Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			     Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			     Univcoord_T chroffset);

#undef T
#endif

