#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memset() */
#include <math.h>
#include <ctype.h>		/* for tolower() */
#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"

#include "genomebits_mismatches.h" /* For MISMATCH_EXTRA */

#include "reader.h"
#include "oligo.h"

#include "list.h"
#include "complement.h"
#include "compress.h"
#include "extension-search.h"	/* For Elt_gc */
#include "tr-extension-search.h"	/* For Tr_elt_gc */

#include "path.h"
#include "path-solve.h"
#include "path-eval.h"


/* Note FORMULA: formulas for querypos <-> diagonal (diagterm in call to Indexdb_read) are:

plus: diagonal = position + querylength - querypos
minus: diagonal = position + querypos + index1part

For minus, the index1part is needed in call to Indexdb_read because
position is stored at beginning of plus oligomer, which corresponds to
end of minus oligomer.  As a result, we have the following formulas:

high genomic position = diagonal (corresponds to querypos =
querylength for plus, and querypos = 0 for minus)

low genomic position = diagonal - querylength (corresponds to querypos
= 0 for plus, and querypos = querylength for minus)

Holds when we use Reader_T to read from 5' end of forward query and 3'
end of revcomp query simultaneously.  If we create a queryrc sequence,
then we can use just the plus formula, and convert the query
coordinates later.

*/


/* Affects only transcriptome-guided genomic alignment on paired-end
   reads.  If one end stops at transcriptome results and the other end
   requires genomic results, continues the first end to find genomic
   results.  Eliminates the greedy and potentially false alignment of
   one end to a transcript.  However, if transcriptome procedures miss
   alignments, then this adds time and memory significantly. */
/* #define AVOID_UNEVEN_LEVELS 1 */

#define LOCALDB_REGION_SIZE 65536

#define NO_EXTENSIONS_BEFORE_ZERO 1

#define ALLOW_MIDDLE_ALIGNMENTS 1

/* #define EXTRACT_GENOMICSEG 1 */
#ifdef EXTRACT_GENOMICSEG
#define MAX_INDEXSIZE 8
#endif


/* MAX_NALIGNMENTS of 2 vs 1 gets 1600 improvements in 275,000 reads */
/* MAX_NALIGNMENTS of 3 vs 2 gets 96 improvements in 275,000 reads */
#define MAX_NALIGNMENTS 3

#define MAX_ALLOCATION 200

#define PAIRMAX_ADDITIONAL 10000 /* Allows for finding of unpaired GMAP alignments beyond pairmax */

/* static int kmer_search_sizelimit = 100; */
/* static int stage1hr_sizelimit = 3000; */
/* static int extension_search_sizelimit = 3000; */

static Indexdb_T indexdb_fwd;
static Indexdb_T indexdb_rev;
static Indexdb_T indexdb_tr;

static EF64_T repetitive_ef64;

static int index1part;
static int index1interval;
static int index1part_tr;
static int index1interval_tr;

static int leftreadshift;
static Oligospace_T oligobase_mask; /* same as kmer_mask */
static int leftreadshift_tr;
static Oligospace_T oligobase_mask_tr; /* same as kmer_mask */

static Chrpos_T positive_gap_distance;

static Transcriptome_T transcriptome;


#define A_CHAR 0x0
#define C_CHAR 0x1
#define G_CHAR 0x2
#define T_CHAR 0x3


/* Originally allowed only 1, to print only unique translocations.
   But need to allow enough to avoid missing some translocations. */
/* For transcript splicing, need to increase MAXCHIMERAPATHS */
/* #define MAXCHIMERAPATHS 100 */
#define MAXCHIMERAPATHS 10000

#define NREQUIRED_FAST 2	/* For candidate generation using
				   multimiss.  A value of 2 implies 
				   specificity of a 24-mer, which
				   should be low for a human-sized
				   genome */

#define MAX_INDEX1INTERVAL 3



/* Overall flow */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Stage1_list_paths */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Stage1_init_end_gen */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Stage1_init_end_oligos_tr and Stage1_init_end_positions_tr */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* fill all oligos */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* fill all positions */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* consolidate_paired_results and choose_among_paired */ 
#ifdef DEBUG16
#define debug16(x) x
#else
#define debug16(x)
#endif


#define T Stage1_T


static void
Stage1_init_computation (T this) {

  this->sense_trnums = (Trnum_T *) NULL;
  this->sense_troffsets = (Trcoord_T *) NULL;
  this->sense_trhighs = (Trcoord_T *) NULL;
  this->_sense_trdiagonals = (Trcoord_T *) NULL;
  this->sense_tstarts = (int *) NULL;
  this->sense_tends = (int *) NULL;
  this->n_sense_trdiagonals = 0;

  this->antisense_trnums = (Trnum_T *) NULL;
  this->antisense_troffsets = (Trcoord_T *) NULL;
  this->antisense_trhighs = (Trcoord_T *) NULL;
  this->_antisense_trdiagonals = (Trcoord_T *) NULL;
  this->antisense_tstarts = (int *) NULL;
  this->antisense_tends = (int *) NULL;
  this->n_antisense_trdiagonals = 0;

  this->partial_sense_trpaths = (List_T) NULL;
  this->partial_antisense_trpaths = (List_T) NULL;
  this->sense_trpaths = (List_T) NULL;
  this->antisense_trpaths = (List_T) NULL;
  this->exact_paths_computed_p = false;

  this->all_univdiagonals_gplus = (Univcoord_T *) NULL;
  this->all_univdiagonals_gminus = (Univcoord_T *) NULL;
  this->all_auxinfo_gplus = (Auxinfo_T *) NULL;
  this->all_auxinfo_gminus = (Auxinfo_T *) NULL;
  this->all_pathstores_gplus = (Pathstore_T *) NULL;
  this->all_pathstores_gminus = (Pathstore_T *) NULL;
  this->all_nunivdiagonals_gplus = 0;
  this->all_nunivdiagonals_gminus = 0;

  this->extension_gplus = (Univcoord_T *) NULL;
  this->extension_qstart_gplus = (int *) NULL;
  this->extension_qend_gplus = (int *) NULL;
  this->nextension_gplus = 0;

  this->extension_gminus = (Univcoord_T *) NULL;
  this->extension_qstart_gminus = (int *) NULL;
  this->extension_qend_gminus = (int *) NULL;
  this->nextension_gminus = 0;

  this->exhaustive_gplus = (Univcoord_T *) NULL; /* aligned */
  this->exhaustive_qstart_gplus = (int *) NULL;
  this->exhaustive_qend_gplus = (int *) NULL;
  this->exhaustive_counts_gplus = (int *) NULL;
  this->nexhaustive_gplus = 0;

  this->exhaustive_gminus = (Univcoord_T *) NULL; /* aligned */
  this->exhaustive_qstart_gminus = (int *) NULL;
  this->exhaustive_qend_gminus = (int *) NULL;
  this->exhaustive_counts_gminus = (int *) NULL;
  this->nexhaustive_gminus = 0;

  /* this->unsolved_sense_paths_gplus = (List_T) NULL; */
  /* this->unsolved_sense_paths_gminus = (List_T) NULL; */
  /* this->unsolved_antisense_paths_gplus = (List_T) NULL; */
  /* this->unsolved_antisense_paths_gminus = (List_T) NULL; */

  /* this->unextended_sense_paths_gplus = (List_T) NULL; */
  /* this->unextended_sense_paths_gminus = (List_T) NULL; */
  /* this->unextended_antisense_paths_gplus = (List_T) NULL; */
  /* this->unextended_antisense_paths_gminus = (List_T) NULL; */


  /* Completed paths are stored as arrays */
  this->sense_paths_gplus = (Path_T *) NULL;
  this->sense_paths_gminus = (Path_T *) NULL;
  this->antisense_paths_gplus = (Path_T *) NULL;
  this->antisense_paths_gminus = (Path_T *) NULL;


  /* Uses Listpool_T procedures */
  this->queryfwd_plus_set = (List_T) NULL;
  this->queryfwd_minus_set = (List_T) NULL;
  this->queryrev_plus_set = (List_T) NULL;
  this->queryrev_minus_set = (List_T) NULL;

  this->tr_queryfwd_plus_set = (List_T) NULL;
  this->tr_queryfwd_minus_set = (List_T) NULL;
  this->tr_queryrev_plus_set = (List_T) NULL;
  this->tr_queryrev_minus_set = (List_T) NULL;

  this->exact_paths_computed_p = false;
  this->all_oligos_gen_filledp = false;
  this->all_positions_gen_filledp = false;

  return;
}


T
Stage1_new (bool first_read_p) {
  T new = (T) MALLOC(sizeof(*new));

  new->first_read_p = first_read_p;

  new->reader = Reader_new(/*oligosize*/index1part);

  if (transcriptome == NULL) {
    new->tr_reader = (Reader_T) NULL;
  } else if (index1part_tr == index1part) {
    new->tr_reader = new->reader;
  } else {
    new->tr_reader = Reader_new(/*oligosize*/index1part_tr);
  }

  Stage1_init_computation(new);

  return new;
}


#ifdef DEBUG1
void
Stage1_dump (T this, int querylength) {
  int query_lastpos = querylength - index1part, querypos;
  int i;

  printf("Stage1_dump\n");
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->validp[querypos] == false) {
      printf("plus  %d invalid:\n",querypos);
    } else if (this->plus_retrievedp[querypos] == true) {
      printf("plus  %d +%d %s (%d):",querypos,
	     this->plus_diagterms[querypos],
	     Oligo_one_nt(this->forward_oligos[querypos],index1part),
	     this->plus_npositions[querypos]);
#ifdef RAW
      for (i = 0; i < this->plus_npositions[querypos]; i++) {
	printf(" %u",this->plus_positions[querypos][i]);
      }
#else
      for (i = 0; i < this->plus_npositions[querypos]; i++) {
	printf(" %u",this->plus_positions[querypos][i] + /*diagterm*/(querylength - querypos));
      }
#endif
      printf("\n");
    }
    if (this->validp[querypos] == false) {
      printf("minus %d invalid:\n",querypos);
    } else if (this->minus_retrievedp[querypos] == true) {
      printf("minus %d +%d %s (%d):",
	     querypos,this->minus_diagterms[querypos],
	     Oligo_one_nt(this->revcomp_oligos[querypos],index1part),
	     this->minus_npositions[querypos]);
#ifdef RAW
      for (i = 0; i < this->minus_npositions[querypos]; i++) {
	printf(" %u",this->minus_positions[querypos][i]);
      }
#else
      for (i = 0; i < this->minus_npositions[querypos]; i++) {
	printf(" %u",this->minus_positions[querypos][i] + /*diagterm*/(querypos + index1part));
      }
#endif
      printf("\n");
    }
  }
  printf("\n");
  return;
}
#endif


#ifdef DEBUG1
void
Stage1_dump_tr (T this, int querylength) {
  int query_lastpos = querylength - index1part_tr, querypos;
  int i;

  printf("Stage1_dump_tr\n");

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->tr_validp[querypos] == false) {
      printf("plus  %d invalid:\n",querypos);
    } else if (this->tr_plus_retrievedp[querypos] == true) {
      printf("plus  %d +%d %s (%d):",
	     querypos,this->tr_plus_diagterms[querypos],
	     Oligo_one_nt(this->tr_forward_oligos[querypos],index1part_tr),
	     this->tr_plus_npositions[querypos]);
#ifdef RAW
      for (i = 0; i < this->tr_plus_npositions[querypos]; i++) {
	printf(" %u",this->tr_plus_positions[querypos][i]);
      }
#else
      for (i = 0; i < this->tr_plus_npositions[querypos]; i++) {
	printf(" %u",this->tr_plus_positions[querypos][i] + /*diagterm*/(querylength - querypos));
      }
#endif
      printf("\n");
    }
  }
  printf("\n");

  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->tr_validp[querypos] == false) {
      printf("minus %d invalid:\n",querypos);
    } else if (this->tr_minus_retrievedp[querypos] == true) {
      printf("minus %d +%d %s (%d):",
	     querypos,this->tr_minus_diagterms[querypos],
	     Oligo_one_nt(this->tr_revcomp_oligos[querypos],index1part_tr),
	     this->tr_minus_npositions[querypos]);
#ifdef RAW
      for (i = 0; i < this->tr_minus_npositions[querypos]; i++) {
	printf(" %u",this->tr_minus_positions[querypos][i]);
      }
#else
      for (i = 0; i < this->tr_minus_npositions[querypos]; i++) {
	printf(" %u",this->tr_minus_positions[querypos][i] + /*diagterm*/(querypos + index1part_tr));
      }
#endif
      printf("\n");
    }
  }
  printf("\n");

  return;
}
#endif


/* Fills in ends for genomic exact search */
/* Skips repetitive oligos to find querystart and queryend */
void
Stage1_init_end_gen (int *querystart, int *queryend, T this,
		     int querylength, int genestrand) {
  Oligostate_T last_state;
  Oligospace_T forward, revcomp, forward_oligo, revcomp_oligo;
  int querypos;
  int mod;
  bool donep;

  assert(querylength >= index1part + index1interval - 1);
  /* query_lastpos = querylength - index1part; */

  Reader_reset_ends(this->reader);
  last_state = INIT;
  forward = revcomp = 0;
  querypos = 0;

  donep = false;
  while (donep == false &&
	 (last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,
				    this->reader,genestrand)) != DONE) {
    if (last_state != VALID) {
      /* querypos is not defined when last_state != VALID */
      debug8(printf("oligo at querypos %d is not valid\n",querypos));

    } else if (EF64_presentp(forward & oligobase_mask,repetitive_ef64) == true) {
      debug8(printf("oligo at querypos %d is repetitive\n",querypos));

    } else {
      donep = true;
    }
  }

  if (donep == false) {
    *querystart = querylength;
    *queryend = 0;
    return;
  } else {
    *querystart = querypos;
    mod = 0;

    this->validp[querypos] = true;
    forward_oligo = this->forward_oligos[querypos] = forward & oligobase_mask;
#ifdef LARGE_GENOMES      
    this->plus_npositions[querypos] =
      Indexdb_largeptr(&(this->plus_positions_high[querypos]),&(this->plus_positions[querypos]),
		       /*plus_indexdb*/indexdb_fwd,forward_oligo);
#else
    this->plus_npositions[querypos] =
      Indexdb_ptr(&(this->plus_positions[querypos]),/*plus_indexdb*/indexdb_fwd,forward_oligo);
#endif
    this->plus_diagterms[querypos] = querylength - querypos; /* plus */
    this->plus_retrievedp[querypos] = true;
    debug2(printf("(1) plus_npositions_end5[%d] = %d, oligo %016lX\n",
		  querypos,this->plus_npositions[querypos],forward_oligo));
      
    revcomp_oligo = this->revcomp_oligos[querypos] = (revcomp >> leftreadshift) & oligobase_mask;
#ifdef LARGE_GENOMES
    this->minus_npositions[querypos] =
      Indexdb_largeptr(&(this->minus_positions_high[querypos]),&(this->minus_positions[querypos]),
		       /*minus_indexdb*/indexdb_rev,revcomp_oligo);
#else
    this->minus_npositions[querypos] =
      Indexdb_ptr(&(this->minus_positions[querypos]),/*minus_indexdb*/indexdb_rev,revcomp_oligo);
#endif
    this->minus_diagterms[querypos] = querypos + index1part; /* minus */
    this->minus_retrievedp[querypos] = true;
    mod++;

    while (mod < index1interval &&
	   (last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,
				      this->reader,genestrand)) != DONE) {
      this->validp[querypos] = true;

      forward_oligo = this->forward_oligos[querypos] = forward & oligobase_mask;
#ifdef LARGE_GENOMES      
      this->plus_npositions[querypos] =
	Indexdb_largeptr(&(this->plus_positions_high[querypos]),&(this->plus_positions[querypos]),
			 /*plus_indexdb*/indexdb_fwd,forward_oligo);
#else
      this->plus_npositions[querypos] =
	Indexdb_ptr(&(this->plus_positions[querypos]),/*plus_indexdb*/indexdb_fwd,forward_oligo);
#endif
      this->plus_diagterms[querypos] = querylength - querypos; /* plus */
      this->plus_retrievedp[querypos] = true;
      debug2(printf("(1) plus_npositions_end5[%d] = %d, oligo %016lX\n",
		    querypos,this->plus_npositions[querypos],forward_oligo));
      
      revcomp_oligo = this->revcomp_oligos[querypos] = (revcomp >> leftreadshift) & oligobase_mask;
#ifdef LARGE_GENOMES
      this->minus_npositions[querypos] =
	Indexdb_largeptr(&(this->minus_positions_high[querypos]),&(this->minus_positions[querypos]),
			 /*minus_indexdb*/indexdb_rev,revcomp_oligo);
#else
      this->minus_npositions[querypos] =
	Indexdb_ptr(&(this->minus_positions[querypos]),/*minus_indexdb*/indexdb_rev,revcomp_oligo);
#endif
      this->minus_diagterms[querypos] = querypos + index1part; /* minus */
      this->minus_retrievedp[querypos] = true;
      debug2(printf("(2) minus_npositions_end5[%d] = %d, oligo %016lX\n",
		    querypos,this->minus_npositions[querypos],revcomp_oligo));
      
      debug2(printf("5' end: %s %s: %d plus positions, %d minus positions, genestrand %d\n",
		    Oligo_one_nt(forward_oligo,index1part),Oligo_one_nt(revcomp_oligo,index1part),
		    this->plus_npositions[querypos],this->minus_npositions[querypos],genestrand));
      mod++;
    }
  }


  /* mod here is relative to query_lastpos.  Kmer_search procedures
     need to find correspondences between mod5 and mod3 */
  /* query_lastpos = querylength - index1part; */
  Reader_reset_ends(this->reader);
  last_state = INIT;
  forward = revcomp = 0;

  donep = false;
  while (donep == false &&
	 (last_state = Oligo_next_3(last_state,&querypos,&forward,&revcomp,
				    this->reader,genestrand)) != DONE) {
    if (last_state != VALID) {
      /* querypos is not defined when last_state != VALID */
      debug8(printf("oligo at querypos %d is not valid\n",querypos));

    } else if (EF64_presentp(revcomp & oligobase_mask,repetitive_ef64) == true) {
      debug8(printf("oligo at querypos %d is repetitive\n",querypos));

    } else {
      donep = true;
    }
  }

  if (donep == false) {
    *querystart = querylength;
    *queryend = 0;
    return;

  } else {
    *queryend = querypos;
    mod = 0;

    this->validp[querypos] = true;

    forward_oligo = this->forward_oligos[querypos] = (forward >> leftreadshift) & oligobase_mask;
#ifdef LARGE_GENOMES
    this->plus_npositions[querypos] =
      Indexdb_largeptr(&(this->plus_positions_high[querypos]),&(this->plus_positions[querypos]),
		       /*plus_indexdb*/indexdb_fwd,forward_oligo);
#else
    this->plus_npositions[querypos] =
      Indexdb_ptr(&(this->plus_positions[querypos]),/*plus_indexdb*/indexdb_fwd,forward_oligo);
#endif
    this->plus_diagterms[querypos] = querylength - querypos; /* plus */
    this->plus_retrievedp[querypos] = true;
    debug2(printf("(3) plus_npositions[%d] = %d, oligo %016lX\n",
		  querypos,this->plus_npositions[querypos],forward_oligo));

    revcomp_oligo = this->revcomp_oligos[querypos] = revcomp & oligobase_mask;
#ifdef LARGE_GENOMES
    this->minus_npositions[querypos] =
      Indexdb_largeptr(&(this->minus_positions_high[querypos]),&(this->minus_positions[querypos]),
		       /*minus_indexdb*/indexdb_rev,revcomp_oligo);
#else
    this->minus_npositions[querypos] =
      Indexdb_ptr(&(this->minus_positions[querypos]),/*minus_indexdb*/indexdb_rev,revcomp_oligo);
#endif
    this->minus_diagterms[querypos] = querypos + index1part; /* minus */
    this->minus_retrievedp[querypos] = true;

    mod++;

    while (mod < index1interval &&
	   (last_state = Oligo_next_3(last_state,&querypos,&forward,&revcomp,
				      this->reader,genestrand)) != DONE) {
      this->validp[querypos] = true;

      forward_oligo = this->forward_oligos[querypos] = (forward >> leftreadshift) & oligobase_mask;
#ifdef LARGE_GENOMES
      this->plus_npositions[querypos] =
	Indexdb_largeptr(&(this->plus_positions_high[querypos]),&(this->plus_positions[querypos]),
			 /*plus_indexdb*/indexdb_fwd,forward_oligo);
#else
      this->plus_npositions[querypos] =
	Indexdb_ptr(&(this->plus_positions[querypos]),/*plus_indexdb*/indexdb_fwd,forward_oligo);
#endif
      this->plus_diagterms[querypos] = querylength - querypos; /* plus */
      this->plus_retrievedp[querypos] = true;
      debug2(printf("(3) plus_npositions[%d] = %d, oligo %016lX\n",
		    querypos,this->plus_npositions[querypos],forward_oligo));
      
      revcomp_oligo = this->revcomp_oligos[querypos] = revcomp & oligobase_mask;
#ifdef LARGE_GENOMES
      this->minus_npositions[querypos] =
	Indexdb_largeptr(&(this->minus_positions_high[querypos]),&(this->minus_positions[querypos]),
			 /*minus_indexdb*/indexdb_rev,revcomp_oligo);
#else
      this->minus_npositions[querypos] =
	Indexdb_ptr(&(this->minus_positions[querypos]),/*minus_indexdb*/indexdb_rev,revcomp_oligo);
#endif
      this->minus_diagterms[querypos] = querypos + index1part; /* minus */
      this->minus_retrievedp[querypos] = true;
      debug2(printf("(4) minus_npositions[%d] = %d, oligo %016lX\n",
		    querypos,this->minus_npositions[querypos],revcomp_oligo));

      debug2(printf("3' end: %s %s: %d plus positions, %d minus positions, genestrand %d\n",
		    Oligo_one_nt(forward_oligo,index1part),Oligo_one_nt(revcomp_oligo,index1part),
		    this->plus_npositions[querypos],this->minus_npositions[querypos],genestrand));
      mod++;
    }
  }

  return;
}


void
Stage1_init_end_tr (T this, int querylength) {
  Oligospace_T forward, revcomp, forward_oligo, revcomp_oligo;
  int querypos, query_lastpos;

  Reader_reset_ends(this->tr_reader);
  forward = revcomp = 0;

  if (Oligo_next_5(/*last_state*/INIT,&querypos,&forward,&revcomp,
		   this->tr_reader,/*genestrand*/0) == DONE) {
    /* Skip */

  } else if (querypos == 0) {
    this->tr_validp[querypos] = true;

    forward_oligo = this->tr_forward_oligos[querypos] = forward & oligobase_mask_tr;
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[querypos] = querylength - querypos; /* plus */
    this->tr_plus_retrievedp[querypos] = true;

    revcomp_oligo = this->tr_revcomp_oligos[querypos] = (revcomp >> leftreadshift_tr) & oligobase_mask_tr;
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
    this->tr_minus_diagterms[querypos] = querypos + index1part_tr; /* minus */
    this->tr_minus_retrievedp[querypos] = true;

    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }


  Reader_reset_ends(this->tr_reader);
  forward = revcomp = 0;
  if (Oligo_next_3(/*last_state*/INIT,&querypos,&forward,&revcomp,
		   this->tr_reader,/*genestrand*/0) == DONE) {
    /* Skip */

  } else if (querypos == (query_lastpos = querylength - index1part_tr)) {
    this->tr_validp[querypos] = true;
    forward_oligo = this->tr_forward_oligos[querypos] = (forward >> leftreadshift_tr) & oligobase_mask_tr;
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[querypos] = querylength - querypos; /* plus */
    this->tr_plus_retrievedp[querypos] = true;
      
    revcomp_oligo = this->tr_revcomp_oligos[querypos] = revcomp & oligobase_mask_tr;
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
    this->tr_minus_diagterms[querypos] = querypos + index1part_tr; /* minus */
    this->tr_minus_retrievedp[querypos] = true;
    
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }
	 
  return;
}


#if 0
void
Stage1_init_end_positions_tr (T this, int querylength) {
  Oligospace_T forward_oligo, revcomp_oligo;
  int query_lastpos;

  if (this->tr_validp[0] == true) {
    forward_oligo = this->tr_forward_oligos[0];
    this->tr_plus_npositions[0] =
      Indexdb_ptr(&(this->tr_plus_positions[0]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[0] = querylength; /* querylength - querypos */
    this->tr_plus_retrievedp[0] = true;

    revcomp_oligo = this->tr_revcomp_oligos[0];
    this->tr_minus_npositions[0] =
      Indexdb_ptr(&(this->tr_minus_positions[0]),indexdb_tr,revcomp_oligo);
    this->tr_minus_diagterms[0] = index1part_tr; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[0] = true;

    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[0],this->tr_minus_npositions[0]));
  }

  query_lastpos = querylength - index1part_tr;
  if (this->tr_validp[query_lastpos] == true) {
    forward_oligo = this->tr_forward_oligos[query_lastpos];
    this->tr_plus_npositions[query_lastpos] =
      Indexdb_ptr(&(this->tr_plus_positions[query_lastpos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[query_lastpos] = index1part_tr; /* querylength - querypos */
    this->tr_plus_retrievedp[query_lastpos] = true;
      
    revcomp_oligo = this->tr_revcomp_oligos[query_lastpos];
    this->tr_minus_npositions[query_lastpos] =
      Indexdb_ptr(&(this->tr_minus_positions[query_lastpos]),indexdb_tr,revcomp_oligo);
     this->tr_minus_diagterms[query_lastpos] = querylength; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[query_lastpos] = true;
    
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[query_lastpos],this->tr_minus_npositions[query_lastpos]));
  }

  return;
}
#endif


#if 0
void
Stage1_init_end2_positions_tr (T this, int querylength) {
  Oligospace_T forward_oligo, revcomp_oligo;
  int query_lastpos, querypos;

  /* querypos = 0; */
  if (this->tr_validp[0] == true) {
    forward_oligo = this->tr_forward_oligos[0];
    this->tr_plus_npositions[0] =
      Indexdb_ptr(&(this->tr_plus_positions[0]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[0] = querylength; /* querylength - querypos */
    this->tr_plus_retrievedp[0] = true;

    revcomp_oligo = this->tr_revcomp_oligos[0];
    this->tr_minus_npositions[0] =
      Indexdb_ptr(&(this->tr_minus_positions[0]),indexdb_tr,revcomp_oligo);
    this->tr_minus_diagterms[0] = index1part_tr; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[0] = true;

    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[0],this->tr_minus_npositions[0]));
  }

  querypos = index1part_tr;
  if (this->tr_validp[querypos] == true) {
    forward_oligo = this->tr_forward_oligos[querypos];
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[querypos] = querylength - index1part_tr; /* querylength - querypos */
    this->tr_plus_retrievedp[querypos] = true;

    revcomp_oligo = this->tr_revcomp_oligos[querypos];
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
    this->tr_minus_diagterms[querypos] = index1part_tr + index1part_tr; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[querypos] = true;

    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }


  query_lastpos = querylength - index1part_tr;
  if (this->tr_validp[query_lastpos] == true) {
    forward_oligo = this->tr_forward_oligos[query_lastpos];
    this->tr_plus_npositions[query_lastpos] =
      Indexdb_ptr(&(this->tr_plus_positions[query_lastpos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[query_lastpos] = index1part_tr; /* querylength - querypos */
    this->tr_plus_retrievedp[query_lastpos] = true;
      
    revcomp_oligo = this->tr_revcomp_oligos[query_lastpos];
    this->tr_minus_npositions[query_lastpos] =
      Indexdb_ptr(&(this->tr_minus_positions[query_lastpos]),indexdb_tr,revcomp_oligo);
     this->tr_minus_diagterms[query_lastpos] = querylength; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[query_lastpos] = true;
    
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[query_lastpos],this->tr_minus_npositions[query_lastpos]));
  }

  querypos = query_lastpos - index1part_tr;
  if (this->tr_validp[querypos] == true) {
    forward_oligo = this->tr_forward_oligos[querypos];
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_diagterms[querypos] = index1part_tr + index1part_tr; /* querylength - querypos */
    this->tr_plus_retrievedp[querypos] = true;
      
    revcomp_oligo = this->tr_revcomp_oligos[querypos];
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
     this->tr_minus_diagterms[querypos] = querylength - index1part_tr; /* querypos + index1part_tr */
    this->tr_minus_retrievedp[querypos] = true;
    
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }

  return;
}
#endif


#if 0
/* For Transcriptome_exact2.  Works if we called init_end_positions_tr for TR_EXACT1 before */
void
Stage1_init_end2_positions_tr (T this, int querylength) {
  Oligospace_T forward_oligo, revcomp_oligo;
  int querypos, query_lastpos;

  query_lastpos = querylength - index1part_tr;
  if ((querypos = index1part_tr) > query_lastpos) {
    /* Skip */
  } else if (this->tr_validp[querypos] == true) {
    forward_oligo = this->tr_forward_oligos[querypos];
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_retrievedp[querypos] = true;

    revcomp_oligo = this->tr_revcomp_oligos[querypos];
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
    this->tr_minus_retrievedp[querypos] = true;

    debug3(printf("5' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }

  if ((querypos = query_lastpos - index1part_tr) < 0) {
    /* Skip */
  } else if (this->tr_validp[querypos] == true) {
    forward_oligo = this->tr_forward_oligos[querypos];
    this->tr_plus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_plus_positions[querypos]),indexdb_tr,forward_oligo);
    this->tr_plus_retrievedp[querypos] = true;
      
    revcomp_oligo = this->tr_revcomp_oligos[querypos];
    this->tr_minus_npositions[querypos] =
      Indexdb_ptr(&(this->tr_minus_positions[querypos]),indexdb_tr,revcomp_oligo);
    this->tr_minus_retrievedp[querypos] = true;
    
    debug3(printf("3' end: %s %s: %d plus positions, %d minus positions\n",
		  Oligo_one_nt(forward_oligo,index1part_tr),Oligo_one_nt(revcomp_oligo,index1part_tr),
		  this->tr_plus_npositions[querypos],this->tr_minus_npositions[querypos]));
  }

  return;
}
#endif



#if 0
/* Puts values of positions_end5 and positions_end3 from
   Kmer_search_exact into positions for other methods */
/* Now performed within Stage1_init */
void
Stage1_integrate_end_positions (T this, int querylength) {
  int query_lastpos = querylength - index1part;
  int mod;

  for (mod = 0; mod < index1interval; mod++) {
    /* querypos = mod; */
    this->plus_validp[mod] = true;
    this->plus_retrievedp[mod] = true;
#ifdef LARGE_GENOMES
    this->plus_positions_high[mod] = this->plus_positions_high_end5[mod];
#endif
    this->plus_positions[mod] = this->plus_positions_end5[mod];
    this->plus_npositions[mod] = this->plus_npositions_end5[mod];

    this->plus_validp[query_lastpos-mod] = true;
    this->plus_retrievedp[query_lastpos-mod] = true;
#ifdef LARGE_GENOMES
    this->plus_positions_high[query_lastpos-mod] = this->plus_positions_high_end3[mod];
#endif
    this->plus_positions[query_lastpos-mod] = this->plus_positions_end3[mod];
    this->plus_npositions[query_lastpos-mod] = this->plus_npositions_end3[mod];

    /* Using new sarray and segment-based conventions */
    this->minus_validp[query_lastpos-mod] = true;
    this->minus_retrievedp[query_lastpos-mod] = true;
#ifdef LARGE_GENOMES
    this->minus_positions_high[query_lastpos-mod] = this->minus_positions_high_end5[mod];
#endif
    this->minus_positions[query_lastpos-mod] = this->minus_positions_end5[mod];
    this->minus_npositions[query_lastpos-mod] = this->minus_npositions_end5[mod];

    this->minus_validp[mod] = true;
    this->minus_retrievedp[mod] = true;
#ifdef LARGE_GENOMES
    this->minus_positions_high[mod] = this->minus_positions_high_end3[mod];
#endif
    this->minus_positions[mod] = this->minus_positions_end3[mod];
    this->minus_npositions[mod] = this->minus_npositions_end3[mod];

#if 0
    printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	   mod,this->plus_positions[mod],this->plus_npositions[mod]);
    printf("Initializing plus_positions[%d] to be %p, with count of %d\n",
	 query_lastpos-mod,this->plus_positions[query_lastpos-mod],this->plus_npositions[query_lastpos-mod]);

    printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	   mod,this->minus_positions[mod],this->minus_npositions[mod]);
    printf("Initializing minus_positions[%d] to be %p, with count of %d\n",
	 query_lastpos-mod,this->minus_positions[query_lastpos-mod],this->minus_npositions[query_lastpos-mod]);
#endif
  }

  /* Stage1_dump(this,querylength); */

  return;
}
#endif


/* Instead of Univcoordtable, we need something big enough to store an Oligospace_T object */
/* Worrying about repetitive oligos slows down program */
void
Stage1_fill_all_oligos_gen (T this, int querylength, int genestrand) {
  int querypos;
  Oligostate_T last_state = INIT;
  Oligospace_T forward = 0, revcomp = 0;

#if 0
  Univcoordtable_T oligo_seenp;
  Oligospace_T oligo;
  oligo_seenp = Univcoordtable_new(/*hint*/querylength);
#endif

  /* query_lastpos = querylength - index1part; */
  Reader_reset_ends(this->reader);

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  querypos = 0;
  while ((last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,
				    this->reader,genestrand)) != DONE) {
    if (last_state != VALID) {
      /* querypos is not defined when last_state != VALID */
      debug8(printf("oligo at querypos %d is not valid\n",querypos));

    } else if (EF64_presentp(forward & oligobase_mask,repetitive_ef64) == true) {
      debug8(printf("oligo at querypos %d is repetitive\n",querypos));

    } else {
      /* querypos_rc = query_lastpos - querypos; */
      /* Previously assigned revcomp oligo to minus_oligos[querypos_rc] */
      /* oligo = */ this->forward_oligos[querypos] = forward & oligobase_mask;
      this->revcomp_oligos[querypos] = (revcomp >> leftreadshift) & oligobase_mask;
      debug8(printf("Putting forward oligo %016lX and revcomp oligo %016lX at querypos %d\n",
		    this->forward_oligos[querypos],this->revcomp_oligos[querypos],querypos));

#if 0
      if (Univcoordtable_get(oligo_seenp,oligo) != NULL) {
	/* Handling repetitive sequences */
	debug8(printf("oligo at plus %d already seen, so marking as invalid\n",querypos));
	this->validp[querypos] = false;
      } else {
	this->validp[querypos] = true;
	Univcoordtable_put(oligo_seenp,oligo,(void *) true);
      }
#else
      this->validp[querypos] = true;
#endif
    }
  }

#if 0
  Univcoordtable_free(&oligo_seenp);
#endif

  this->all_oligos_gen_filledp = true;

  return;
}


/* Don't need to worry about repetitive oligos for transcriptome
   methods, and it slows down program */
void
Stage1_fill_all_oligos_tr (T this, int
			   querylength) {
  int querypos;
  Oligostate_T last_state = INIT;
  Oligospace_T forward = 0, revcomp = 0;

#if 0
  Univcoordtable_T oligo_seenp;
  Oligospace_T oligo;
  oligo_seenp = Univcoordtable_new(/*hint*/querylength);
#endif

  /* query_lastpos = querylength - index1part; */
  Reader_reset_ends(this->tr_reader);

  /* Note: leftshifting is done here, rather than in Oligo_lookup */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  querypos = 0;
  while ((last_state = Oligo_next_5(last_state,&querypos,&forward,&revcomp,
				    this->tr_reader,/*genestrand*/0)) != DONE) {
    if (last_state != VALID) {
      /* querypos is not defined when last_state != VALID */
      debug8(printf("oligo at querypos %d is not valid\n",querypos));
    } else {
      /* querypos_rc = query_lastpos - querypos; */
      /* Previously assigned revcomp oligo to minus_oligos[querypos_rc] */
      /* oligo = */ this->tr_forward_oligos[querypos] = forward & oligobase_mask_tr;
      this->tr_revcomp_oligos[querypos] = (revcomp >> leftreadshift_tr) & oligobase_mask_tr;
      debug8(printf("Putting forward oligo %016lX and revcomp oligo %016lX at querypos %d\n",
		    this->tr_forward_oligos[querypos],this->tr_revcomp_oligos[querypos],querypos));

#if 0
      if (Univcoordtable_get(oligo_seenp,oligo) != NULL) {
	/* Handling repetitive sequences */
	debug8(printf("oligo at plus %d already seen, so marking as invalid\n",querypos));
	this->tr_validp[querypos] = false;
      } else {
	this->tr_validp[querypos] = true;
	Univcoordtable_put(oligo_seenp,oligo,(void *) true);
      }
#else
      this->tr_validp[querypos] = true;
#endif
    }
  }

#if 0
  Univcoordtable_free(&oligo_seenp);
#endif

  if (index1part_tr == index1part) {
    this->all_oligos_gen_filledp = true;
  }

  return;
}


void
Stage1_fill_all_positions_gen (int *total_npositions_plus, int *total_npositions_minus,
			       T this, int querylength, int genestrand) {
  int query_lastpos, querypos;
  Indexdb_T plus_indexdb, minus_indexdb;
  int npositions;

  debug9(printf("Filling all positions for %p\n",this));

  if (genestrand == +2) {
    plus_indexdb = indexdb_rev;
    minus_indexdb = indexdb_fwd;
  } else {
    plus_indexdb = indexdb_fwd;
    minus_indexdb = indexdb_rev;
  }

  query_lastpos = querylength - index1part;

  /* *max_npositions_plus = 0; */
  /* *max_npositions_minus = 0; */
  *total_npositions_plus = 0;
  *total_npositions_minus = 0;

  /* Assumes that forward_oligos and revcomp_oligos have been filled
     in (by Stage1_fill_all_oligos */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    debug9(printf("querypos %d, validp %d, plus retrievedp %d, minus retrievedp %d, ",
		  querypos,this->validp[querypos],this->plus_retrievedp[querypos],this->minus_retrievedp[querypos]));

    if (this->validp[querypos] == false) {
      /* Forward and revcomp oligos not valid */
      this->plus_npositions[querypos] = 0;
      this->minus_npositions[querypos] = 0;

    } else {
      this->plus_diagterms[querypos] = querylength - querypos;
      this->minus_diagterms[querypos] = querypos + index1part;

      if (this->plus_retrievedp[querypos] == true) {
	/* No need to do anything */
#if 0
	if ((npositions = this->plus_npositions[querypos]) > *max_npositions_plus) {
	  *max_npositions_plus = npositions;
	}
	*total_npositions_plus += npositions;
#endif
	*total_npositions_plus += this->plus_npositions[querypos];

      } else {
#ifdef LARGE_GENOMES
	npositions = this->plus_npositions[querypos] = 
	  Indexdb_largeptr(&this->plus_positions_high[querypos],&this->plus_positions[querypos],
			   plus_indexdb,this->forward_oligos[querypos]);
#else
	npositions = this->plus_npositions[querypos] = 
	  Indexdb_ptr(&this->plus_positions[querypos],plus_indexdb,this->forward_oligos[querypos]);
#endif

#if 0
	if (npositions > *max_npositions_plus) {
	  *max_npositions_plus = npositions;
	}
#endif
	*total_npositions_plus += npositions;

	this->plus_retrievedp[querypos] = true;
      }

      if (this->minus_retrievedp[querypos] == true) {
	/* No need to do anything */
#if 0
	if ((npositions = this->minus_npositions[querypos]) > *max_npositions_minus) {
	  *max_npositions_minus = npositions;
	}
	*total_npositions_minus += npositions;
#endif
	*total_npositions_minus += this->minus_npositions[querypos];

      } else {
#ifdef LARGE_GENOMES
	npositions = this->minus_npositions[querypos] =
	  Indexdb_largeptr(&this->minus_positions_high[querypos],&this->minus_positions[querypos],
			   minus_indexdb,this->revcomp_oligos[querypos]);
#else
	npositions = this->minus_npositions[querypos] =
	  Indexdb_ptr(&this->minus_positions[querypos],minus_indexdb,this->revcomp_oligos[querypos]);
#endif

#if 0
	if (npositions > *max_npositions_minus) {
	  *max_npositions_minus = npositions;
	}
#endif
	*total_npositions_minus += npositions;

	this->minus_retrievedp[querypos] = true;
      }
    }
    debug9(printf("plus npositions %d, minus npositions %d\n",
		  this->plus_npositions[querypos],this->minus_npositions[querypos]));
  }

  this->all_positions_gen_filledp = true;

  return;
}


void
Stage1_fill_all_positions_tr (T this, int querylength) {
  int query_lastpos, querypos;

  query_lastpos = querylength - index1part_tr;

  /* Assumes that forward_oligos and revcomp_oligos have been filled
     in (by Stage1_fill_all_oligos */
  /* Format is 010llX because 19-mer is maximum k-mer size, which would require 10 chars */
  /* debug(printf("oligobase_mask: %010llX\n",oligobase_mask)); */
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->tr_validp[querypos] == false) {
      /* Forward and revcomp oligos not valid */
    } else {
      if (this->tr_plus_retrievedp[querypos] == true) {
	/* No need to do anything */
      } else {
	this->tr_plus_npositions[querypos] = 
	  Indexdb_ptr(&this->tr_plus_positions[querypos],indexdb_tr,this->tr_forward_oligos[querypos]);
	this->tr_plus_diagterms[querypos] = querylength - querypos;
	this->tr_plus_retrievedp[querypos] = true;
      }

      if (this->tr_minus_retrievedp[querypos] == true) {
	/* No need to do anything */
      } else {
	this->tr_minus_npositions[querypos] =
	  Indexdb_ptr(&this->tr_minus_positions[querypos],indexdb_tr,this->tr_revcomp_oligos[querypos]);
	this->tr_minus_diagterms[querypos] = querypos + index1part_tr;
	this->tr_minus_retrievedp[querypos] = true;
      }
    }
  }

  return;
}


#ifdef DEBUG1
static void
list_trpaths (List_T list, char *destination) {
  List_T p;
  Trpath_T trpath;

  for (p = list; p != NULL; p = List_next(p)) {
    trpath = (Trpath_T) List_head(p);
    printf("Destination %s: ",destination);
    Trpath_print(trpath);
  }

  return;
}
#endif


#ifdef DEBUG1
void
Stage1_list_trpaths (T this) {

  printf("Dump of trpaths\n");
  list_trpaths(this->sense_trpaths,"sense_trpaths");
  list_trpaths(this->antisense_trpaths,"antisense_trpaths");
  printf("End dump\n");

  return;
}
#endif


#if 0
static void
list_paths (List_T list, char *destination, bool expected_sensedir) {
  List_T p;
  Path_T path;

  for (p = list; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    printf("Destination %s: ",destination);
    Path_print(path);
    assert(path->sensedir == expected_sensedir);
  }

  return;
}

static void
list_paths_array (Path_T *paths, int n, char *destination, bool expected_sensedir) {
  int i;
  Path_T path;

  for (i = 0; i < n; i++) {
    path = paths[i];
    printf("Destination %s: ",destination);
    Path_print(path);
    assert(path->sensedir == expected_sensedir);
  }

  return;
}
#endif


#if 0
void
Stage1_list_coords (T this) {
  int i;
  
  printf("=====%s Dump of coords=====\n",this->first_read_p ? "5'" : "3'");
  for (i = 0; i < this->nunique_sense_coords_gplus; i++) {
    printf("sense plus %u\n",this->sense_coords_gplus[i]);
  }
  for (i = 0; i < this->nunique_antisense_coords_gplus; i++) {
    printf("antisense plus %u\n",this->antisense_coords_gplus[i]);
  }
  for (i = 0; i < this->nunique_sense_coords_gminus; i++) {
    printf("sense minus %u\n",this->sense_coords_gminus[i]);
  }
  for (i = 0; i < this->nunique_antisense_coords_gminus; i++) {
    printf("antisense minus %u\n",this->antisense_coords_gminus[i]);
  }
  printf("=====%s End of coords=====\n",this->first_read_p ? "5'" : "3'");

  return;
}
#endif



void
Stage1_list_extension (T this) {
  int i;

  printf("=====%s Dump of extension=====\n",this->first_read_p ? "5'" : "3'");
  printf(">plus extension\n");
  for (i = 0; i < this->nextension_gplus; i++) {
    printf("%u %d..%d\n",this->extension_gplus[i],this->extension_qstart_gplus[i],this->extension_qend_gplus[i]);
  }
  printf(">minus extension\n");
  for (i = 0; i < this->nextension_gminus; i++) {
    printf("%u %d..%d\n",this->extension_gminus[i],this->extension_qstart_gminus[i],this->extension_qend_gminus[i]);
  }
  printf("=====%s End of extension=====\n",this->first_read_p ? "5'" : "3'");

  return;
}


void
Stage1_list_exhaustive (T this) {
  int i;

  printf("=====%s Dump of exhaustive=====\n",this->first_read_p ? "5'" : "3'");
  printf(">plus exhaustive\n");
  for (i = 0; i < this->nexhaustive_gplus; i++) {
    printf("%u %d..%d %d\n",
	   this->exhaustive_gplus[i],this->exhaustive_qstart_gplus[i],this->exhaustive_qend_gplus[i],
	   this->exhaustive_counts_gplus[i]);
  }
  printf(">minus exhaustive\n");
  for (i = 0; i < this->nexhaustive_gminus; i++) {
    printf("%u %d..%d %d\n",
	   this->exhaustive_gminus[i],this->exhaustive_qstart_gminus[i],this->exhaustive_qend_gminus[i],
	   this->exhaustive_counts_gminus[i]);
  }
  printf("=====%s End of exhaustive=====\n",this->first_read_p ? "5'" : "3'");

  return;
}

void
Stage1_list_all_univdiagonals (T this) {
  int i;
  Auxinfo_T auxinfo;
  Pathstore_T pathstore;
  List_T p;
  Path_T path;

  printf("=====%s Dump of all univdiagonals=====\n",this->first_read_p ? "5'" : "3'");
  printf(">plus all univdiagonals (%d)\n",this->all_nunivdiagonals_gplus);
  for (i = 0; i < this->all_nunivdiagonals_gplus; i++) {
    printf("%u.  Auxinfo methods:",this->all_univdiagonals_gplus[i]);
    for (auxinfo = this->all_auxinfo_gplus[i]; auxinfo != NULL; auxinfo = auxinfo->rest) {
      printf(" %s",Method_string(auxinfo->method));
    }
    printf("\n");

    pathstore = this->all_pathstores_gplus[i];
    for (p = pathstore->unextended_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Unextended sense: "); Path_print(path);
    }
    for (p = pathstore->unextended_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Unextended antisense: "); Path_print(path);
    }
    for (p = pathstore->complete_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Sense: "); Path_print(path);
    }
    for (p = pathstore->complete_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Antisense: "); Path_print(path);
    }
    for (p = pathstore->best_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Best sense: "); Path_print(path);
    }
    for (p = pathstore->best_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Best antisense: "); Path_print(path);
    }
  }

  printf(">minus all univdiagonals (%d)\n",this->all_nunivdiagonals_gminus);
  for (i = 0; i < this->all_nunivdiagonals_gminus; i++) {
    printf("%u.  Auxinfo methods:",this->all_univdiagonals_gminus[i]);
    for (auxinfo = this->all_auxinfo_gminus[i]; auxinfo != NULL; auxinfo = auxinfo->rest) {
      printf(" %s",Method_string(auxinfo->method));
    }
    printf("\n");

    pathstore = this->all_pathstores_gminus[i];
    for (p = pathstore->unextended_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Unextended sense: "); Path_print(path);
    }
    for (p = pathstore->unextended_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Unextended antisense: "); Path_print(path);
    }
    for (p = pathstore->complete_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Sense: "); Path_print(path);
    }
    for (p = pathstore->complete_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Antisense: "); Path_print(path);
    }
    for (p = pathstore->best_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Best sense: "); Path_print(path);
    }
    for (p = pathstore->best_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      printf("  Best antisense: "); Path_print(path);
    }
    printf("\n");
  }
  printf("=====%s End of all univdiagonals=====\n",this->first_read_p ? "5'" : "3'");

  return;
}


void
Stage1_trdiagonals_gc (T this) {

  FREE(this->sense_trnums);
  FREE(this->sense_troffsets);
  FREE(this->sense_trhighs);
  this->sense_trnums = (Trnum_T *) NULL;
  this->sense_troffsets = (Trcoord_T *) NULL;
  this->sense_trhighs = (Trcoord_T *) NULL;

  FREE_ALIGN(this->_sense_trdiagonals);
  this->_sense_trdiagonals = (Trcoord_T *) NULL;
  this->n_sense_trdiagonals = 0;

  FREE(this->sense_tstarts);
  FREE(this->sense_tends);
  this->sense_tstarts = (int *) NULL;
  this->sense_tends = (int *) NULL;

  FREE(this->antisense_trnums);
  FREE(this->antisense_troffsets);
  FREE(this->antisense_trhighs);
  this->antisense_trnums = (Trnum_T *) NULL;
  this->antisense_troffsets = (Trcoord_T *) NULL;
  this->antisense_trhighs = (Trcoord_T *) NULL;

  FREE_ALIGN(this->_antisense_trdiagonals);
  this->_antisense_trdiagonals = (Trcoord_T *) NULL;
  this->n_antisense_trdiagonals = 0;

  FREE(this->antisense_tstarts);
  FREE(this->antisense_tends);
  this->antisense_tstarts = (int *) NULL;
  this->antisense_tends = (int *) NULL;

  return;
}


/* Takes care of items that depend on pools or allocated through
   computations, so need to call this before resetting memory pools.
   Other fields in Stage3_free_nonpool can be handled after resetting
   memory pools */
static void
Stage1_free_computation (T this, Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
			 Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool, Intlistpool_T intlistpool,
			 Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool, Trpathpool_T trpathpool,
			 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Stage1_trdiagonals_gc(this);

  FREE_ALIGN(this->all_univdiagonals_gplus);
  FREE_ALIGN(this->all_univdiagonals_gminus);
  this->all_univdiagonals_gplus = (Univcoord_T *) NULL;
  this->all_univdiagonals_gminus = (Univcoord_T *) NULL;

  FREE(this->extension_gplus); /* Not aligned */
  FREE(this->extension_qstart_gplus);
  FREE(this->extension_qend_gplus);
  this->extension_gplus = (Univcoord_T *) NULL;
  this->extension_qstart_gplus = (int *) NULL;
  this->extension_qend_gplus = (int *) NULL;
  this->nextension_gplus = 0;

  FREE(this->extension_gminus);
  FREE(this->extension_qstart_gminus);
  FREE(this->extension_qend_gminus);
  this->extension_gminus = (Univcoord_T *) NULL;
  this->extension_qstart_gminus = (int *) NULL;
  this->extension_qend_gminus = (int *) NULL;
  this->nextension_gminus = 0;

  FREE_ALIGN(this->exhaustive_gplus);
  FREE(this->exhaustive_qstart_gplus);
  FREE(this->exhaustive_qend_gplus);
  FREE(this->exhaustive_counts_gplus);
  this->exhaustive_gplus = (Univcoord_T *) NULL; /* aligned */
  this->exhaustive_qstart_gplus = (int *) NULL;
  this->exhaustive_qend_gplus = (int *) NULL;
  this->exhaustive_counts_gplus = (int *) NULL;
  this->nexhaustive_gplus = 0;

  FREE_ALIGN(this->exhaustive_gminus);
  FREE(this->exhaustive_qstart_gminus);
  FREE(this->exhaustive_qend_gminus);
  FREE(this->exhaustive_counts_gminus);
  this->exhaustive_gminus = (Univcoord_T *) NULL; /* aligned */
  this->exhaustive_qstart_gminus = (int *) NULL;
  this->exhaustive_qend_gminus = (int *) NULL;
  this->exhaustive_counts_gminus = (int *) NULL;
  this->nexhaustive_gminus = 0;

  Pathstore_gc(this->all_pathstores_gplus,this->all_nunivdiagonals_gplus,
	       pathstorepool,intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(this->all_pathstores_gminus,this->all_nunivdiagonals_gminus,
	       pathstorepool,intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  this->all_pathstores_gplus = (Pathstore_T *) NULL;
  this->all_pathstores_gminus = (Pathstore_T *) NULL;

  Auxinfo_gc(this->all_auxinfo_gplus,this->all_nunivdiagonals_gplus,
	     univdiagpool,auxinfopool);
  Auxinfo_gc(this->all_auxinfo_gminus,this->all_nunivdiagonals_gminus,
	     univdiagpool,auxinfopool);
  this->all_auxinfo_gplus = (Auxinfo_T *) NULL;
  this->all_auxinfo_gminus = (Auxinfo_T *) NULL;

  Trpath_gc(&this->partial_sense_trpaths,
	    intlistpool,uintlistpool,listpool,pathpool,trpathpool,hitlistpool);
  Trpath_gc(&this->partial_antisense_trpaths,
	    intlistpool,uintlistpool,listpool,pathpool,trpathpool,hitlistpool);
  this->partial_sense_trpaths = (List_T) NULL;
  this->partial_antisense_trpaths = (List_T) NULL;

  Trpath_gc(&this->sense_trpaths,
	    intlistpool,uintlistpool,listpool,pathpool,trpathpool,hitlistpool);
  Trpath_gc(&this->antisense_trpaths,
	    intlistpool,uintlistpool,listpool,pathpool,trpathpool,hitlistpool);
  this->sense_trpaths = (List_T) NULL;
  this->antisense_trpaths = (List_T) NULL;

  Elt_gc(&this->queryfwd_plus_set,listpool);
  Elt_gc(&this->queryfwd_minus_set,listpool);
  Elt_gc(&this->queryrev_plus_set,listpool);
  Elt_gc(&this->queryrev_minus_set,listpool);

  Tr_elt_gc(&this->tr_queryfwd_plus_set,listpool,trdiagpool);
  Tr_elt_gc(&this->tr_queryfwd_minus_set,listpool,trdiagpool);
  Tr_elt_gc(&this->tr_queryrev_plus_set,listpool,trdiagpool);
  Tr_elt_gc(&this->tr_queryrev_minus_set,listpool,trdiagpool);

  return;
}



static void
Stage1_free_sequence (T this) {
  FREE(this->mismatch_positions_alloc);
  FREE(this->positions_alloc);
  Indelinfo_free(&this->indelinfo);
  Spliceinfo_free(&this->spliceinfo);

  Mergeinfo_uint4_free(&this->mergeinfo_tr);
#ifdef LARGE_GENOMES
  Mergeinfo_uint8_free(&this->mergeinfo);
#else
  Mergeinfo_uint4_free(&this->mergeinfo);
#endif

#ifdef LARGE_GENOMES
  FREE(this->stream_high_alloc);
  FREE(this->stream_low_alloc);
#endif

  FREE(this->streamptr_alloc);
  FREE(this->streamsize_alloc);

  FREE(this->querypos_diagterm_alloc);
  FREE(this->retrievedp_allocated);

#ifdef LARGE_GENOMES
  FREE(this->positions_high_allocated);
#endif
  FREE(this->positions_allocated);
  FREE(this->npositions_allocated);

  FREE(this->plus_diagterms);
  FREE(this->minus_diagterms);

  if (transcriptome == NULL) {
    /* Skip */

  } else {
    FREE(this->tr_retrievedp_allocated);
    FREE(this->tr_positions_allocated);
    FREE(this->tr_npositions_allocated);

    FREE(this->tr_plus_diagterms);
    FREE(this->tr_minus_diagterms);
  }

  FREE(this->validp);
  FREE(this->forward_oligos);
  FREE(this->revcomp_oligos);

  if (transcriptome == NULL) {
    /* Skip */
  } else if (index1part_tr == index1part) {
    /* Skip */
  } else {
    FREE(this->tr_validp);
    FREE(this->tr_forward_oligos);
    FREE(this->tr_revcomp_oligos);
 }

  return;
}


void
Stage1_free (T *old, Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
	     Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool, Intlistpool_T intlistpool,
	     Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Trpathpool_T trpathpool,
	     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Stage1_free_sequence(*old);
  Stage1_free_computation(*old,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
			  intlistpool,uintlistpool,univcoordlistpool,
			  listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

  if (transcriptome == NULL) {
    /* Skip */
  } else if (index1part_tr == index1part) {
    /* Skip */
  } else {
    Reader_free(&(*old)->tr_reader);
  }

  Reader_free(&(*old)->reader);
  FREE(*old);

  return;
}


/* Previously, tried to re-use old memory if querylength has not
   changed, but this seems to result in growth in virtual memory usage */
void
Stage1_set_sequence (T this, Shortread_T queryseq, Localdb_T localdb) {

  char *queryuc_ptr = Shortread_queryuc_ptr(queryseq);
  int querylength = Shortread_fulllength(queryseq);

  int overhang = index1interval - 1;
  int overhang_tr = index1interval_tr - 1;

  /* Maximum number of localdb regions possible, based on region size of 65536 bp */
  int max_nstreams, max_localdb_nregions;

  Reader_set_sequence(this->reader,queryuc_ptr,/*querystart*/0,/*queryend*/querylength);

  /* Need to allocate (max_mismatches+MISMATCH_EXTRA), where
     max_mismatches is provided to Genomebits_mismatches_left or
     Genome_mismatches_right */
  this->mismatch_positions_alloc = (int *) MALLOC((querylength+MISMATCH_EXTRA)*sizeof(int));
  this->positions_alloc = (int *) MALLOC((querylength+1)*sizeof(int));
  this->indelinfo = Indelinfo_new(querylength);
  this->spliceinfo = Spliceinfo_new(querylength);

  this->mergeinfo_tr = Mergeinfo_uint4_new(querylength,/*max_localdb_distance*/0);
#ifdef LARGE_GENOMES
  this->mergeinfo = Mergeinfo_uint8_new(querylength,positive_gap_distance);
#else
  this->mergeinfo = Mergeinfo_uint4_new(querylength,positive_gap_distance);
#endif

  /* Memory previously allocated for Segment_identify in segment-search.c, and
     Merge_diagonals in kmer-search.c (which needs four sets of
     arrays) */
#ifdef LARGE_GENOMES
  this->stream_high_alloc = (unsigned char **) MALLOC(4*querylength*sizeof(unsigned char *));
  this->gplus_stream_high_array_5 = &(this->stream_high_alloc[0]);
  this->gminus_stream_high_array_5 = &(this->stream_high_alloc[querylength]);
  this->gplus_stream_high_array_3 = &(this->stream_high_alloc[2*querylength]);
  this->gminus_stream_high_array_3 = &(this->stream_high_alloc[3*querylength]);

  this->stream_low_alloc = (UINT4 **) MALLOC(4*querylength*sizeof(UINT4 *));
  this->gplus_stream_low_array_5 = &(this->stream_low_alloc[0]);
  this->gminus_stream_low_array_5 = &(this->stream_low_alloc[querylength]);
  this->gplus_stream_low_array_3 = &(this->stream_low_alloc[2*querylength]);
  this->gminus_stream_low_array_3 = &(this->stream_low_alloc[3*querylength]);
#endif


  if (querylength == 0) {
    max_nstreams = 1;
  } else if (localdb == NULL) {
    max_nstreams = 1;		/* Should be 0, but don't want to allocate 0 bytes */
  } else {
    max_localdb_nregions = (positive_gap_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
    max_nstreams = max_localdb_nregions * querylength;
  }

  if (transcriptome != NULL) {
    /* For tplus and tminus */
    if (max_nstreams < 2 * querylength) {
      max_nstreams = 2 * querylength;
    }
  } else {
    if (max_nstreams < querylength) {
      max_nstreams = querylength;
    }
  }

  this->max_nstreams = max_nstreams;
  this->streamptr_alloc = (Univcoord_T **) MALLOC(this->max_nstreams*sizeof(Univcoord_T *));
  this->streamsize_alloc = (int *) MALLOC(this->max_nstreams * sizeof(int));

  /* Can assign pointers, even if transcriptome is NULL */
  this->tplus_streamsize_array = &(this->streamsize_alloc[0]);
  this->tminus_streamsize_array = &(this->streamsize_alloc[querylength]);

  if (querylength == 0) {
    this->querypos_diagterm_alloc = (int *) NULL;
    this->tplus_diagterm_array = (int *) NULL;
    this->tminus_diagterm_array = (int *) NULL;
  } else {
    this->querypos_diagterm_alloc = (int *) MALLOC(2*querylength*sizeof(int));
    this->tplus_diagterm_array = &(this->querypos_diagterm_alloc[0]);
    this->tminus_diagterm_array = &(this->querypos_diagterm_alloc[querylength]);
  }

  this->retrievedp_allocated = (bool *) CALLOC(2 * (querylength+overhang),sizeof(bool));
  this->plus_retrievedp = &(this->retrievedp_allocated[overhang]);
  this->minus_retrievedp = &(this->retrievedp_allocated[(querylength+overhang)+overhang]);

#ifdef LARGE_GENOMES
  this->positions_high_allocated = (unsigned char **) CALLOC(2 * (querylength+overhang),sizeof(unsigned char *));
  this->plus_positions_high = &(this->positions_high_allocated[overhang]);
  this->minus_positions_high = &(this->positions_high_allocated[(querylength+overhang)+overhang]);
#endif
  this->positions_allocated = (UINT4 **) CALLOC(2 * (querylength+overhang),sizeof(UINT4 *));
  this->plus_positions = &(this->positions_allocated[overhang]);
  this->minus_positions = &(this->positions_allocated[(querylength+overhang)+overhang]);

  this->npositions_allocated = (int *) CALLOC(2 * (querylength+overhang),sizeof(int));
  this->plus_npositions = &(this->npositions_allocated[overhang]);
  this->minus_npositions = &(this->npositions_allocated[(querylength+overhang)+overhang]);

  if (querylength == 0) {
    this->plus_diagterms = (int *) NULL;
    this->minus_diagterms = (int *) NULL;
  } else {
    this->plus_diagterms = (int *) MALLOC(querylength*sizeof(int));
    this->minus_diagterms = (int *) MALLOC(querylength*sizeof(int));
  }

  if (transcriptome == NULL) {
    this->tr_retrievedp_allocated = (bool *) NULL;
    this->tr_plus_retrievedp = this->tr_minus_retrievedp = (bool *) NULL;

    this->tr_positions_allocated = (UINT4 **) NULL;
    this->tr_plus_positions = this->tr_minus_positions = (UINT4 **) NULL;

    this->tr_npositions_allocated = (int *) NULL;
    this->tr_plus_npositions = this->tr_minus_npositions = (int *) NULL;

    this->tr_plus_diagterms = (int *) NULL;
    this->tr_minus_diagterms = (int *) NULL;

  } else {
    this->tr_retrievedp_allocated = (bool *) CALLOC(2 * (querylength+overhang_tr),sizeof(bool));
    this->tr_plus_retrievedp = &(this->tr_retrievedp_allocated[overhang_tr]);
    this->tr_minus_retrievedp = &(this->tr_retrievedp_allocated[(querylength+overhang_tr)+overhang_tr]);

    this->tr_positions_allocated = (UINT4 **) CALLOC(2 * (querylength+overhang_tr),sizeof(UINT4 *));
    this->tr_plus_positions = &(this->tr_positions_allocated[overhang_tr]);
    this->tr_minus_positions = &(this->tr_positions_allocated[(querylength+overhang_tr)+overhang_tr]);

    this->tr_npositions_allocated = (int *) CALLOC(2 * (querylength+overhang_tr),sizeof(int));
    this->tr_plus_npositions = &(this->tr_npositions_allocated[overhang_tr]);
    this->tr_minus_npositions = &(this->tr_npositions_allocated[(querylength+overhang_tr)+overhang_tr]);

    this->tr_plus_diagterms = (int *) MALLOC(querylength*sizeof(int));
    this->tr_minus_diagterms = (int *) MALLOC(querylength*sizeof(int));
  }

  Reader_set_sequence(this->reader,queryuc_ptr,/*querystart*/0,/*queryend*/querylength);
  
  this->validp = (bool *) CALLOC(querylength+overhang,sizeof(bool));
  this->forward_oligos = (Oligospace_T *) MALLOC((querylength+overhang)*sizeof(Oligospace_T));
  this->revcomp_oligos = (Oligospace_T *) MALLOC((querylength+overhang)*sizeof(Oligospace_T));

  if (transcriptome == NULL) {
    this->tr_reader = (Reader_T) NULL;
    this->tr_validp = (bool *) NULL;
    this->tr_forward_oligos = this->tr_revcomp_oligos = (Oligospace_T *) NULL;

  } else if (index1part_tr == index1part) {
    this->tr_validp = this->validp;
    this->tr_forward_oligos = this->forward_oligos;
    this->tr_revcomp_oligos = this->revcomp_oligos;

  } else {
    /* tr_reader is distinct from reader */
    Reader_set_sequence(this->tr_reader,queryuc_ptr,/*querystart*/0,/*queryend*/querylength);

    this->tr_validp = (bool *) CALLOC(querylength+overhang_tr,sizeof(bool));
    this->tr_forward_oligos = (Oligospace_T *) MALLOC((querylength+overhang_tr)*sizeof(Oligospace_T));
    this->tr_revcomp_oligos = (Oligospace_T *) MALLOC((querylength+overhang_tr)*sizeof(Oligospace_T));
 }

  return;
}



bool
Stage1_collect_paths (List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,
		      T this, Hitlistpool_T hitlistpool) {
  bool foundp = false;

  Pathstore_collect_paths(&foundp,&(*sense_paths_gplus),&(*antisense_paths_gplus),
			  this->all_pathstores_gplus,this->all_nunivdiagonals_gplus,hitlistpool);
  Pathstore_collect_paths(&foundp,&(*sense_paths_gminus),&(*antisense_paths_gminus),
			  this->all_pathstores_gminus,this->all_nunivdiagonals_gminus,hitlistpool);

  return foundp;
}


void
Stage1hr_setup (Indexdb_T indexdb_fwd_in, Indexdb_T indexdb_rev_in, Indexdb_T indexdb_tr_in,
		EF64_T repetitive_ef64_in, int index1part_in, int index1interval_in,
		int index1part_tr_in, int index1interval_tr_in, 
		int max_deletionlen, Chrpos_T shortsplicedist,
		Transcriptome_T transcriptome_in) {

  indexdb_fwd = indexdb_fwd_in;
  indexdb_rev = indexdb_rev_in;
  indexdb_tr = indexdb_tr_in;

  repetitive_ef64 = repetitive_ef64_in;

  index1part = index1part_in;
  index1interval = index1interval_in;
  index1part_tr = index1part_tr_in;
  index1interval_tr = index1interval_tr_in;

#ifdef HAVE_64_BIT
  leftreadshift = 64 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
  leftreadshift_tr = 64 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#else
  leftreadshift = 32 - index1part - index1part;
  oligobase_mask = ~(~ (Oligospace_T) 0 << 2*index1part);
  leftreadshift_tr = 32 - index1part_tr - index1part_tr;
  oligobase_mask_tr = ~(~ (Oligospace_T) 0 << 2*index1part_tr);
#endif

  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;
  transcriptome = transcriptome_in;

  return;
}
