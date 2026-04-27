/* $Id: dynprog_genome.h 224700 2021-11-08 16:46:32Z twu $ */
#ifndef DYNPROG_GENOME_INCLUDED
#define DYNPROG_GENOME_INCLUDED

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "pairpool.h"
#include "chrnum.h"
#include "iit-read.h"
#include "types.h"
#include "dynprog.h"

#define T Dynprog_T

extern void
Dynprog_genome_setup (bool novelsplicingp_in,
		      IIT_T splicing_iit_in, int *splicing_divint_crosstable_in,
		      int donor_typeint_in, int acceptor_typeint_in,
		      int user_open_in, int user_extend_in, bool user_dynprog_p_in);


extern List_T
Dynprog_genome_gap (int *dynprogindex, int *new_leftgenomepos, int *new_rightgenomepos,
		    double *left_prob, double *right_prob,

		    int *traceback_score, int *nmatches, int *nmismatches, int *nopens,
		    int *nindels, int *exonhead, int *introntype,

		    T dynprogL, T dynprogR,
		    char *rsequence, char *rsequenceuc, int rlength, int glengthL, int glengthR, 
		    int roffset, int goffsetL, int rev_goffsetR, 
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    int cdna_direction, bool watsonp, int genestrand, bool jump_late_p,
		    Genome_T genome, Genome_T genomealt, Pairpool_T pairpool, int extraband_paired,
		    double defect_rate, int maxpeelback, bool halfp, bool finalp);


#undef T
#endif

