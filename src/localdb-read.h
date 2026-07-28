#ifndef LOCALDB_READ_INCLUDED
#define LOCALDB_READ_INCLUDED

typedef struct Localdb_T *Localdb_T;

#include "access.h"
#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "univdiag.h"
#include "compress.h"
#include "genome.h"
#include "genomebits.h"
#include "stage1hr.h"

#define T Localdb_T

extern void
Localdb_setup (Genome_T genomecomp_in, Univcoord_T genomelength_in,
	       int index1part_in, int index1interval_in);

extern void
Localdb_free (T *old);

extern T
Localdb_new (char *genomesubdir, char *fileroot,
	     Access_mode_T localdb_access, bool sharedp,
	     bool multiple_sequences_p, bool preload_shared_memory_p, bool unload_shared_memory_p);

extern int
Localdb_get (bool *sortedp, bool *trimmedp,
	     int *matchlength, int *nmismatches, Univdiag_T *univdiags_alloc,
	     T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
	     Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
	     char *queryptr, int pos5, int pos3, int querylength,
	     Univcoord_T anchor_univdiagonal, Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, int genestrand,
	     Genomebits_T genomebits, int nmismatches_allowed,
	     bool extend5p, bool trim5p, bool trim3p);

extern Univcoord_T
Localdb_get_one_low (Univdiag_T *univdiags_alloc,
		     T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		     Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		     char *queryptr, int querylength,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		     Compress_T query_compress, bool plusp, int genestrand,
		     Genomebits_T genomebits, int nmismatches_allowed);

extern Univcoord_T
Localdb_get_one_high (Univdiag_T *univdiags_alloc,
		      T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		      Stage1_T stage1, int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		      char *queryptr, int querylength,
		      Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		      Compress_T query_compress, bool plusp, int genestrand,
		      Genomebits_T genomebits, int nmismatches_allowed);

#undef T
#endif


