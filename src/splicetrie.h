/* $Id: splicetrie.h 223511 2020-11-14 15:50:08Z twu $ */
#ifndef SPLICETRIE_INCLUDED
#define SPLICETRIE_INCLUDED

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "list.h"
#include "intlist.h"
#include "splicetrie_build.h"	/* For Splicetype_T */
#include "compress.h"

#include "dynprog.h"
#include "pairpool.h"


extern void
Splicetrie_setup (Univcoord_T *splicesites_in, Genomecomp_T *splicefrags_ref_in, Genomecomp_T *splicefrags_alt_in,
		  Trieoffset_T *trieoffsets_obs_in, Triecontent_T *triecontents_obs_in,
		  Trieoffset_T *trieoffsets_max_in, Triecontent_T *triecontents_max_in,
		  bool snpp_in, bool amb_closest_p_in, bool amb_clip_p_in, int min_shortend_in);

extern List_T
Splicetrie_solve_end5 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *revsequence1, char *revsequenceuc1,
		       int length1, int length2, int revoffset1, int revoffset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		       int extraband_end, double defect_rate);

extern List_T
Splicetrie_solve_end3 (List_T best_pairs, Triecontent_T *triecontents, Trieoffset_T *trieoffsets, int j,
		       Univcoord_T knownsplice_limit_low, Univcoord_T knownsplice_limit_high,

		       int *finalscore, int *nmatches, int *nmismatches,
		       int *nopens, int *nindels, bool *knownsplicep, int *ambig_end_length,
		       int *threshold_miss_score, int obsmax_penalty, int perfect_score,

		       Univcoord_T anchor_splicesite, char *splicejunction, char *splicejunction_alt,
		       int splicelength, int contlength, Splicetype_T far_splicetype,
		       Univcoord_T chroffset, Univcoord_T chrhigh,
		       int *dynprogindex, Dynprog_T dynprog, 
		       char *sequence1, char *sequenceuc1,
		       int length1, int length2, int offset1, int offset2,
		       int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		       Genome_T genome, Genome_T genomealt, Pairpool_T pairpool,
		       int extraband_end, double defect_rate);

#endif
