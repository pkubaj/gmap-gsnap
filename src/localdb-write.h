#ifndef LOCALDB_WRITE_INCLUDED
#define LOCALDB_WRITE_INCLUDED

#include "genome.h"
#include "univcoord.h"

extern void
Localdb_write (char *saindex16file, char *sarray16file,
	       char *sarray8file, char *sasort16file,
	       Genome_T genomecomp, Univcoord_T genomelength);

#endif


