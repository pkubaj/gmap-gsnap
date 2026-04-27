/* $Id: ae1e3813d723eb653b894a2076ba1c5b997a0be7 $ */
#ifndef INDEL_INCLUDED
#define INDEL_INCLUDED

typedef struct Indelinfo_T *Indelinfo_T;

#include "bool.h"
#include "list.h"
#include "chrnum.h"
#include "genomicpos.h"
#include "compress.h"
#include "genomebits.h"

struct Indelinfo_T {
  int *int_memory;

  int *mismatch_positions_left;
  int *mismatch_positions_right;
};

extern void
Indel_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
	     int min_indel_end_matches_in, bool maskedp_in);

extern void
Indelinfo_free (Indelinfo_T *old);
extern Indelinfo_T
Indelinfo_new (int querylength);

extern int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
				Univcoord_T univdiagonal_i, int indels, Univcoord_T chrhigh,
				int *mismatch_positions_left, int nmismatches_left,
				int *mismatch_positions_right, int nmismatches_right,
				Genomebits_T ome, Genomebits_T ome_alt, Compress_T query_compress,
				int pos5, int pos3, int querylength, Indelinfo_T indelinfo,
				bool plusp, int genestrand, bool want_lowest_coordinate_p);

int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
			       Univcoord_T univdiagonal_i, int indels, Univcoord_T chrhigh,
			       int *mismatch_positions_left, int nmismatches_left,
			       int *mismatch_positions_right, int nmismatches_right,
			       Genomebits_T ome, Genomebits_T ome_alt, Compress_T query_compress,
			       int pos5, int pos3, int querylength, Indelinfo_T indelinfo,
			       bool plusp, int genestrand, bool want_lowest_coordinate_p);

#endif

