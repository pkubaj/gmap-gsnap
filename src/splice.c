#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "splice.h"

#include <stdio.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "sense.h"

#include "genomebits_count.h"
#include "genomebits_mismatches.h"
#include "genomebits_trim.h"
#include "spliceends.h"

#include "genome.h"
#include "genome_sites.h"
#include "maxent.h"
#include "maxent_hr.h"
#include "univcoord.h"
#include "complement.h"


/* Causes problems with counting mismatches */
/* #define TRIM_AT_CHROMOSOME_BOUNDS 1 */

#define MIN_INTRONLEN 30	/* Same as in junction.c */

#define EXTRA_MISMATCHES 8
#define EXTRA_MATCHES 8 /* Using a special calculation for cases with a potential middle exon or spliceindel */
#define MIN_EXONLEN 3		/* Was 9, but allowing for microexons of 3 bp */
#define MISMATCHES_AT_SITE 2

/* #define ALLOW_ATYPICAL_MIDDLE 1 */

#define LOWPROB_SUPPORT 20
#define MIN_SUPPORT_SPLICE 6	/* First level in sufficient_support_p */
#define MIN_SUPPORT_SPLICE_PLUS_INDEL 12

#define MIN_PROB 0.85		/* For non-salvage */
#define PROB_SLOP 0.2
#define MISMATCHES_SLOP 1

/* #define MIN_SPLICE_PROB 0.4 */	/* Skip if both probs are less than this */
#define MIN_SPLICE_PLUS_INDEL_PROB 0.5 /* Skip if either prob is less than this */


#if 0
/* Creates issues with ambiguous substrings */
#define LOCALSPLICING_NMATCHES_SLOP 1
#else
#define LOCALSPLICING_NMATCHES_SLOP 0
#endif
#define LOCALSPLICING_PROB_SLOP 0.05

#define SLOP 1

#define add_bounded(x,plusterm,highbound) ((x + (plusterm) >= highbound) ? (highbound - 1) : x + (plusterm))
#define subtract_bounded(x,minusterm,lowbound) ((x < lowbound + (minusterm)) ? lowbound : x - (minusterm))

/* Splice_resolve */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Splice_resolve_fusion */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Splice_resolve_fusion, counting mismatches */
#ifdef DEBUG2A
#define debug2a(x) x
#else
#define debug2a(x)
#endif

/* known splicing */
#ifdef DEBUG4S
#define debug4s(x) x
#else
#define debug4s(x)
#endif

/* Splice_accept_p */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



#define MAX_INSERTIONLEN 3
#define MAX_DELETIONLEN 6

static bool novelsplicingp;
static bool splice_evalp;

static Support_T support_tables = NULL;

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;

static char complCode[128] = COMPLEMENT_LC;


#ifdef DEBUG2
static void
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return;
}
#endif


void
Splice_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
	      int max_insertionlen_in, int max_deletionlen_in,
	      bool novelsplicingp_in, bool splice_evalp_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  /* max_insertionlen = max_insertionlen_in; */
  /* max_deletionlen = max_deletionlen_in; */
  novelsplicingp = novelsplicingp_in;
  splice_evalp = splice_evalp_in;

  return;
}


void
Splice_pass2_setup (Support_T support_tables_in) {

  support_tables = support_tables_in;

  return;
}


void
Spliceinfo_free (Spliceinfo_T *old) {
  FREE((*old)->mismatch_positions_left1);
  FREE((*old)->mismatch_positions_right1);
  FREE((*old)->mismatch_positions_left2);
  FREE((*old)->mismatch_positions_right2);

  FREE((*old)->segmenti_sites_alloc1);
  FREE((*old)->segmenti_types_alloc1);
  FREE((*old)->segmenti_sites_alloc2);
  FREE((*old)->segmenti_types_alloc2);
  FREE((*old)->segmenti_probs_alloc);

  FREE((*old)->segmentk1_sites_alloc1);
  FREE((*old)->segmentk1_types_alloc1);
  FREE((*old)->segmentk1_sites_alloc2);
  FREE((*old)->segmentk1_types_alloc2);
  FREE((*old)->segmentk1_probs_alloc);

  FREE((*old)->segmentk2_sites_alloc1);
  FREE((*old)->segmentk2_types_alloc1);
  FREE((*old)->segmentk2_sites_alloc2);
  FREE((*old)->segmentk2_types_alloc2);
  FREE((*old)->segmentk2_probs_alloc);

  FREE((*old)->segmentj_sites_alloc1);
  FREE((*old)->segmentj_types_alloc1);
  FREE((*old)->segmentj_sites_alloc2);
  FREE((*old)->segmentj_types_alloc2);
  FREE((*old)->segmentj_probs_alloc);

  FREE(*old);

  return;
}


/* The only information accessed externally is ambig_qstarts and ambig_qends */
Spliceinfo_T
Spliceinfo_new (int querylength) {
  Spliceinfo_T new = (Spliceinfo_T) MALLOC(sizeof(*new));

  /* MISMATCH_EXTRA defined in genomebits_mismatches.h */
  new->mismatch_positions_left1 = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));
  new->mismatch_positions_right1 = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));
  new->mismatch_positions_left2 = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));
  new->mismatch_positions_right2 = (int *) MALLOC((querylength + MISMATCH_EXTRA)*sizeof(int));

  new->segmenti_sites_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmenti_types_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmenti_sites_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmenti_types_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmenti_probs_alloc = (double *) MALLOC((querylength + 1)*sizeof(double));

  new->segmentk1_sites_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk1_types_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk1_sites_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk1_types_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk1_probs_alloc = (double *) MALLOC((querylength + 1)*sizeof(double));

  new->segmentk2_sites_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk2_types_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk2_sites_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk2_types_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentk2_probs_alloc = (double *) MALLOC((querylength + 1)*sizeof(double));

  new->segmentj_sites_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentj_types_alloc1 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentj_sites_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentj_types_alloc2 = (int *) MALLOC((querylength + 1)*sizeof(int));
  new->segmentj_probs_alloc = (double *) MALLOC((querylength + 1)*sizeof(double));

  return new;
}

/* Note: contents of spliceinfo are filled in by kmer-search and path-solve procedures */


#if 0
/* Same as in spliceindel.c */
static bool
sufficient_support_p (int adj_support, double splice_prob) {
  debug1(printf("Checking for sufficient splice prob, based on adj_support %d and splice prob %.2f\n",
		adj_support,splice_prob));

  if (splice_prob > 0.95) {
    return (adj_support >= 6) ? true : false; /* threshold set to MIN_SUPPORT_SPLICE */

  } else if (splice_prob > 0.90) {
    return (adj_support >= 8) ? true : false;

  } else if (splice_prob > 0.80) {
    return (adj_support >= 12) ? true : false;

  } else if (splice_prob > 0.50) {
    /* Previously was 15 */
    return (adj_support >= 20) ? true : false;

  } else {
    /* Previously was 20 */
    return (adj_support >= 25) ? true : false;
  }
}
#endif


/* Can have lower splice prob on one side when the distances are small */
/* Note: supporti and supportj may not correspond to splice_prob_1 and splice_prob_2 */
bool
Splice_accept_p (int splice_qpos, int querylength,
		 double splice_prob_i, double splice_prob_j,
		 Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		 int supporti, int nmismatchesi, int supportj, int nmismatchesj,
		 bool plusp, bool sense_forward_p, bool innerp) {

  Univcoord_T splice_position_i, splice_position_j;
  Chrpos_T splice_distance;
  bool acceptp;

  splice_distance = univdiagonal_j - univdiagonal_i;
  debug9(printf("Checking for sufficient splice probs %f and %f, based on splice distance %u and support %d and %d, innerp %d\n",
		splice_prob_i,splice_prob_j,splice_distance,supporti,supportj,innerp));

  if (support_tables != NULL) {

    /* Use best support from the entire dataset, rather than from this read, and rely upon that */
    splice_position_i = univdiagonal_i - (Univcoord_T) querylength + (Univcoord_T) splice_qpos;
    splice_position_j = univdiagonal_j - (Univcoord_T) querylength + (Univcoord_T) splice_qpos;

    debug9(printf(" => Support_presentp returning %d\n",
		  Support_presentp(support_tables,splice_position_i,splice_position_j,
				   plusp,sense_forward_p)));
    return Support_presentp(support_tables,splice_position_i,splice_position_j,
			    plusp,sense_forward_p);

  } else if (splice_evalp == false) {
    return true;

  } else if (splice_prob_i > 1.0 && splice_prob_j > 1.0) {
    /* Known intron (Note that we are allowing all of them in Trpath_convert procedures.) */
    debug9(printf("Known intron => accept\n"));
    return true;

  } else if (splice_prob_i < 0.5 && splice_prob_j < 0.5) {
    debug9(printf("Bad intron => reject\n"));
    return false;

#if 0
  } else if (strictp == false && (supporti + supportj) < 4*(nmismatchesi + nmismatchesj)) {
    debug9(printf("Total nmismatches poor\n"));
    return false;
    
  } else if (strictp == true && (supporti < 4*nmismatchesi || supportj < 4*nmismatchesj)) {
    debug9(printf("Individual nmismatches poor\n"));
    return false;
#endif

  } else if (splice_prob_i < 0.2 || splice_prob_j < 0.2) {
    /* One splice site is very poor */
    debug9(printf("Splice with one very poor probability => reject\n"));
    return false;

  } else if (supporti < 4*nmismatchesi || supportj < 4*nmismatchesj) {
    debug9(printf("Individual nmismatches poor\n"));
    return false;

  } else if (/*innerp == false &&*/splice_distance >= HUGE_SPLICE) {
    /* Special requirements for large splices, inner or outer */
    if (splice_prob_i >= 0.8 && splice_prob_j >= 0.8 &&
	supporti >= 20 && supportj >= 20) {
      /* Sufficient probabilities and support for a huge splice */
      debug9(printf("Probabilities and support are sufficient for a huge splice => accept\n"));
      return true;
    } else {
      debug9(printf("Probabilities and support are insufficient for a huge splice => reject\n"));
      return false;
    }

  } else if (innerp == false && splice_distance >= MEDIUM_SPLICE) {
    /* Special requirements for large outer splices */
    if (splice_prob_i >= 0.8 && splice_prob_j >= 0.8 &&
	supporti >= 12 && supportj >= 12) {
      /* Sufficient probabilities and support for a medium splice */
      debug9(printf("Probabilities and support are sufficient for a medium splice => accept\n"));
      return true;
    } else {
      debug9(printf("Probabilities and support are insufficient for a medium splice => reject\n"));
      return false;
    }

  } else if (splice_prob_i >= 0.8 && splice_prob_j >= 0.8) {
    /* Both are good */
    debug9(printf("Both probabilities are good => accept\n"));
    return true;

  } else if (supporti >= 20 && supportj >= 20) {
    /* Sufficient support */
    debug9(printf("Both supports are long => accept\n"));
    return true;

  } else if (splice_distance > SMALL_SPLICE) {
    debug9(printf("Long splice needs two good probabilities or support => reject\n"));
    return false;

  } else if (splice_prob_i < 0.8 && splice_prob_j < 0.8) {
    /* Both are poor */
    debug9(printf("Short splice needs one good probability or support => reject\n"));
    return false;

  } else if (supporti < 8 || supportj < 8) {
    /* Allow one splice prob to be low if there is sufficient support */
    /* Lack support on at least one end */
    debug9(printf("Short splice but one poor probability needs support of 8 => reject\n"));
    return false;

  } else {
    /* Has support on both ends */
    debug9(printf("Short splice but one poor probability has support of 8 => accept\n"));
    return true;
  }
}


#if 0
/* Do not compare against true or false */
/* Want loose criterion, otherwise, we incur slowdown from having to
   run GSNAP algorithm */
static int
sufficient_splice_prob_local (int support, int nmismatches, double spliceprob) {
  support -= 3*nmismatches;
  if (support <= 9) {
    return (spliceprob > 0.80);
  } else if (support <= 12) {
    return (spliceprob > 0.70);
  } else if (support <= 15) {
    return (spliceprob > 0.60);
  } else if (support <= 25) {
    return (spliceprob > 0.50);
  } else {
    return (spliceprob > 0.40);
  }
}
#endif


/* segmenti_types is positive for a known site, -1 for a GTGC-AG site, and -2 for an AT-AC site */
static int
compute_donor_sites (double *optimal_prob, int **segmenti_sites, int **segmenti_types, double **segmenti_probs,
		     int *sites_alloc1, int *types_alloc1,
		     int *sites_alloc2, int *types_alloc2,
		     double *probs_alloc, int splice_qpos_low, int splice_qpos_high, int querylength,
		     Univcoord_T univdiagonal_i, Univcoord_T chroffset,
		     Knownsplicing_T knownsplicing, Intlistpool_T intlistpool) {

  int segmenti_nsites = 0, i;
  Univcoord_T *endpoints;
  uint64_t low_rank, high_rank, next_rank, rank;
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  double prob;

  *segmenti_sites = sites_alloc1;
  *segmenti_types = types_alloc1;

  if (splice_qpos_low >= splice_qpos_high) {
    return 0;
  }

  if (knownsplicing != NULL) {
    endpoints = Knownsplicing_donors(&low_rank,&high_rank,knownsplicing,
				     univdiagonal_i,querylength,splice_qpos_low,splice_qpos_high);
    debug1(printf("Knownsplicing donors at %u+%d..%d yields low_rank %lu to high_rank %lu\n",
		  segmenti_left,splice_qpos_low,splice_qpos_high,low_rank,high_rank));
    /* Collapse duplicates because we don't care about partners */
    rank = low_rank;
    while (rank < high_rank) {
      debug4s(printf("Setting known donor %d for segmenti at %u..%u\n",(int) rank,endpoints[2*rank],endpoints[2*rank+1]));
      (*segmenti_types)[segmenti_nsites] = (int) rank;
      (*segmenti_sites)[segmenti_nsites++] = endpoints[2*rank] - segmenti_left;

      next_rank = rank + 1;
      while (next_rank < high_rank && endpoints[2*next_rank] == endpoints[2*rank]) {
	next_rank += 1;
      }
      rank = next_rank;
    }

#ifdef DEBUG1
    printf("Found %d known donori sites:",segmenti_nsites);
    for (i = 0; i < segmenti_nsites; i++) {
      printf(" %d",(*segmenti_sites)[i]);
      if ((*segmenti_types)[i] >= 0) {
	printf(" [#%d]",(*segmenti_types)[i]);
      }
    }
    printf("\n");
#endif
  }

  if (novelsplicingp == true) {
    sites_alloc1[segmenti_nsites] = querylength; /* Needed to terminate search in Genome_sites */

    *segmenti_sites = sites_alloc2;
    *segmenti_types = types_alloc2;
    segmenti_nsites = Genome_donor_sites(*segmenti_sites,*segmenti_types,
					 /*old_knownpos*/sites_alloc1,/*old_knowni*/types_alloc1,
					 segmenti_left,/*pos5*/splice_qpos_low,/*pos3*/splice_qpos_high,
					 intlistpool);
  }

  *optimal_prob = 0.0;
  *segmenti_probs = probs_alloc;
  for (i = 0; i < segmenti_nsites; i++) {
    if ((*segmenti_types)[i] >= 0) {
      (*segmenti_probs)[i] = 1.0;
      *optimal_prob = 1.0;
    } else {
      prob = (*segmenti_probs)[i] = Maxent_hr_donor_prob(segmenti_left + (*segmenti_sites)[i],chroffset);
      if (prob > *optimal_prob) {
	*optimal_prob = prob;
      }
    }
  }

#ifdef DEBUG1
  printf("Found %d donori sites for univdiagonal %u, %d..%d\n",
	 segmenti_nsites,univdiagonal_i,splice_qpos_low,splice_qpos_high);
  for (i = 0; i < segmenti_nsites; i++) {
    printf(" %d",(*segmenti_sites)[i]);
    if ((*segmenti_types)[i] >= 0) {
      printf(" [#%d]",(*segmenti_types)[i]);
    } else if ((*segmenti_types)[i] == -1) {
      printf(" [GTGC]");
    } else if ((*segmenti_types)[i] == -2) {
      printf(" [AT]");
    }
    printf(" (%.6f)",Maxent_hr_donor_prob(segmenti_left + (*segmenti_sites)[i],chroffset));
  }
  printf("\n");
#endif

  return segmenti_nsites;
}


static int
compute_acceptor_sites (double *optimal_prob, int **segmentj_sites, int **segmentj_types, double **segmentj_probs,
			int *sites_alloc1, int *types_alloc1,
			int *sites_alloc2, int *types_alloc2, 
			double *probs_alloc,
			int splice_qpos_low, int splice_qpos_high, int querylength,
			Univcoord_T univdiagonal_j, Univcoord_T chroffset,
			Knownsplicing_T knownsplicing, Intlistpool_T intlistpool) {

  int segmentj_nsites = 0, i;
  Univcoord_T *endpoints;
  uint64_t low_rank, high_rank, next_rank, rank;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
  double prob;

  *segmentj_sites = sites_alloc1;
  *segmentj_types = types_alloc1;

  if (splice_qpos_low >= splice_qpos_high) {
    return 0;
  }

  if (knownsplicing != NULL) {
    endpoints = Knownsplicing_acceptors(&low_rank,&high_rank,knownsplicing,
					univdiagonal_j,querylength,splice_qpos_low,splice_qpos_high);
    debug1(printf("Knownsplicing acceptors at %u+%d..%d yields low_rank %lu to high_rank %lu\n",
		  segmentj_left,splice_qpos_low,splice_qpos_high,low_rank,high_rank));
    /* Collapse duplicates because we don't care about partners */
    rank = low_rank;
    while (rank < high_rank) {
      debug4s(printf("Setting known acceptor %d for segmentj at %u..%u\n",(int) rank,endpoints[2*rank],endpoints[2*rank+1]));
      (*segmentj_types)[segmentj_nsites] = (int) rank;
      (*segmentj_sites)[segmentj_nsites++] = endpoints[2*rank] - segmentj_left;

      next_rank = rank + 1;
      while (next_rank < high_rank && endpoints[2*next_rank] == endpoints[2*rank]) {
	next_rank += 1;
      }
      rank = next_rank;
    }
#ifdef DEBUG1
    printf("Found %d known acceptorj sites:",segmentj_nsites);
    for (i = 0; i < segmentj_nsites; i++) {
      printf(" %d",(*segmentj_sites)[i]);
      if ((*segmentj_types)[i] >= 0) {
	printf(" [#%d]",(*segmentj_types)[i]);
      }
    }
    printf("\n");
#endif
  }

  if (novelsplicingp == true) {
    sites_alloc1[segmentj_nsites] = querylength; /* Needed to terminate search in Genome_sites */

    *segmentj_sites = sites_alloc2;
    *segmentj_types = types_alloc2;
    segmentj_nsites = Genome_acceptor_sites(*segmentj_sites,*segmentj_types,
					    /*old_knownpos*/sites_alloc1,/*old_knowni*/types_alloc1,
					    segmentj_left,/*pos5*/splice_qpos_low,/*pos3*/splice_qpos_high,
					    intlistpool);
  }

  *optimal_prob = 0.0;
  *segmentj_probs = probs_alloc;
  for (i = 0; i < segmentj_nsites; i++) {
    if ((*segmentj_types)[i] >= 0) {
      (*segmentj_probs)[i] = 1.0;
      *optimal_prob = 1.0;
    } else {
      prob = (*segmentj_probs)[i] = Maxent_hr_acceptor_prob(segmentj_left + (*segmentj_sites)[i],chroffset);
      if (prob > *optimal_prob) {
	*optimal_prob = prob;
      }
    }
  }

#ifdef DEBUG1
  printf("Found %d acceptorj sites for univdiagonal %u, %d..%d\n",
	 segmentj_nsites,univdiagonal_j,splice_qpos_low,splice_qpos_high);
  for (i = 0; i < segmentj_nsites; i++) {
    printf(" %d",(*segmentj_sites)[i]);
    if ((*segmentj_types)[i] >= 0) {
      printf(" [#%d]",(*segmentj_types)[i]);
    } else if ((*segmentj_types)[i] == -1) {
      printf(" [AG]");
    } else if ((*segmentj_types)[i] == -2) {
      printf(" [AC]");
    }
    /* printf(" (%.6f)",Maxent_hr_acceptor_prob(segmentj_left + (*segmentj_sites)[i],chroffset)); */
    printf(" (%.6f)",(*segmentj_probs)[i]);
  }
  printf("\n");
#endif

  return segmentj_nsites;
}


static int
compute_antiacceptor_sites (double *optimal_prob, int **segmenti_sites, int **segmenti_types, double **segmenti_probs,
			    int *sites_alloc1, int *types_alloc1,
			    int *sites_alloc2, int *types_alloc2,
			    double *probs_alloc,
			    int splice_qpos_low, int splice_qpos_high, int querylength,
			    Univcoord_T univdiagonal_i, Univcoord_T chroffset,
			    Knownsplicing_T knownsplicing, Intlistpool_T intlistpool) {

  int segmenti_nsites = 0, i;
  Univcoord_T *endpoints;
  uint64_t low_rank, high_rank, next_rank, rank;
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  double prob;

  *segmenti_sites = sites_alloc1;
  *segmenti_types = types_alloc1;

  if (splice_qpos_low >= splice_qpos_high) {
    return 0;
  }

  if (knownsplicing != NULL) {
    endpoints = Knownsplicing_antiacceptors(&low_rank,&high_rank,knownsplicing,
					    univdiagonal_i,querylength,splice_qpos_low,splice_qpos_high);
    debug1(printf("Knownsplicing antiacceptors at %u+%d..%d yields low_rank %lu to high_rank %lu\n",
		  segmenti_left,splice_qpos_low,splice_qpos_high,low_rank,high_rank));
    /* Collapse duplicates because we don't care about partners */
    rank = low_rank;
    while (rank < high_rank) {
      debug4s(printf("Setting known antiacceptor %d for segmenti at %u..%u\n",(int) rank,endpoints[2*rank],endpoints[2*rank+1]));
      (*segmenti_types)[segmenti_nsites] = (int) rank;
      (*segmenti_sites)[segmenti_nsites++] = endpoints[2*rank] - segmenti_left;
      
      next_rank = rank + 1;
      while (next_rank < high_rank && endpoints[2*next_rank] == endpoints[2*rank]) {
	next_rank += 1;
      }
      rank = next_rank;
    }
#ifdef DEBUG1
    printf("Found %d known antiacceptori sites:",segmenti_nsites);
    for (i = 0; i < segmenti_nsites; i++) {
      printf(" %d",(*segmenti_sites)[i]);
      if ((*segmenti_types)[i] >= 0) {
	printf(" [#%d]",(*segmenti_types)[i]);
      }
      printf("\n");
    }
#endif
  }

  if (novelsplicingp == true) {
    sites_alloc1[segmenti_nsites] = querylength; /* Needed to terminate search in Genome_sites */

    *segmenti_sites = sites_alloc2;
    *segmenti_types = types_alloc2;
    segmenti_nsites = Genome_antiacceptor_sites(*segmenti_sites,*segmenti_types,
						/*old_knownpos*/sites_alloc1,/*old_knowni*/types_alloc1,
						segmenti_left,splice_qpos_low,splice_qpos_high,intlistpool);
  }

  *optimal_prob = 0.0;
  *segmenti_probs = probs_alloc;
  for (i = 0; i < segmenti_nsites; i++) {
    if ((*segmenti_types)[i] >= 0) {
      (*segmenti_probs)[i] = 1.0;
      *optimal_prob = 1.0;
    } else {
      prob = (*segmenti_probs)[i] = Maxent_hr_antiacceptor_prob(segmenti_left + (*segmenti_sites)[i],chroffset);
      if (prob > *optimal_prob) {
	*optimal_prob = prob;
      }
    }
  }

#ifdef DEBUG1
  printf("Found %d antiacceptori sites for univdiagonal %u, %d..%d\n",
	 segmenti_nsites,univdiagonal_i,splice_qpos_low,splice_qpos_high);
  for (i = 0; i < segmenti_nsites; i++) {
    printf(" %d",(*segmenti_sites)[i]);
    if ((*segmenti_types)[i] >= 0) {
      printf(" [#%d]",(*segmenti_types)[i]);
      } else if ((*segmenti_types)[i] == -1) {
      printf(" [AG]");
    } else if ((*segmenti_types)[i] == -2) {
      printf(" [AC]");
    }
    printf(" (%.6f)",Maxent_hr_antiacceptor_prob(segmenti_left + (*segmenti_sites)[i],chroffset));
  }
  printf("\n");
#endif

  return segmenti_nsites;
}


static int
compute_antidonor_sites (double *optimal_prob, int **segmentj_sites, int **segmentj_types, double **segmentj_probs,
			 int *sites_alloc1, int *types_alloc1,
			 int *sites_alloc2, int *types_alloc2,
			 double *probs_alloc,
			 int splice_qpos_low, int splice_qpos_high, int querylength,
			 Univcoord_T univdiagonal_j, Univcoord_T chroffset,
			 Knownsplicing_T knownsplicing, Intlistpool_T intlistpool) {

  int segmentj_nsites = 0, i;
  Univcoord_T *endpoints;
  uint64_t low_rank, high_rank, next_rank, rank;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
  double prob;

  *segmentj_sites = sites_alloc1;
  *segmentj_types = types_alloc1;

  if (splice_qpos_low >= splice_qpos_high) {
    return 0;
  }

  if (knownsplicing != NULL) {
    endpoints = Knownsplicing_antidonors(&low_rank,&high_rank,knownsplicing,
					 univdiagonal_j,querylength,splice_qpos_low,splice_qpos_high);
    debug1(printf("Knownsplicing antidonors at %u+%d..%d yields low_rank %lu to high_rank %lu\n",
		  segmentj_left,splice_qpos_low,splice_qpos_high,low_rank,high_rank));
    /* Collapse duplicates because we don't care about partners */
    rank = low_rank;
    while (rank < high_rank) {
      debug4s(printf("Setting known antidonor %d for segmentj at %u..%u\n",(int) rank,endpoints[2*rank],endpoints[2*rank+1]));
      (*segmentj_types)[segmentj_nsites] = (int) rank;
      (*segmentj_sites)[segmentj_nsites++] = endpoints[2*rank] - segmentj_left;

      next_rank = rank + 1;
      while (next_rank < high_rank && endpoints[2*next_rank] == endpoints[2*rank]) {
	next_rank += 1;
      }
      rank = next_rank;
    }

#ifdef DEBUG1
    printf("Found %d known antidonorj sites:",segmentj_nsites);
    for (i = 0; i < segmentj_nsites; i++) {
      printf(" %d",(*segmentj_sites)[i]);
      if ((*segmentj_types)[i] >= 0) {
	printf(" [#%d]",(*segmentj_types)[i]);
      }
    }
    printf("\n");
#endif
  }

  if (novelsplicingp == true) {
    sites_alloc1[segmentj_nsites] = querylength; /* Needed to terminate search in Genome_sites */

    *segmentj_sites = sites_alloc2;
    *segmentj_types = types_alloc2;
    segmentj_nsites = Genome_antidonor_sites(*segmentj_sites,*segmentj_types,
					     /*old_knownpos*/sites_alloc1,/*old_knowni*/types_alloc1,
					     segmentj_left,splice_qpos_low,splice_qpos_high,intlistpool);

  }

  *optimal_prob = 0.0;
  *segmentj_probs = probs_alloc;
  for (i = 0; i < segmentj_nsites; i++) {
    if ((*segmentj_types)[i] >= 0) {
      (*segmentj_probs)[i] = 1.0;
      *optimal_prob = 1.0;
    } else {
      prob = (*segmentj_probs)[i] = Maxent_hr_antidonor_prob(segmentj_left + (*segmentj_sites)[i],chroffset);
      if (prob > *optimal_prob) {
	*optimal_prob = prob;
      }
    }
  }

#ifdef DEBUG1
  printf("Found %d antidonorj sites for univdiagonal %u, %d..%d\n",
	 segmentj_nsites,univdiagonal_j,splice_qpos_low,splice_qpos_high);
  for (i = 0; i < segmentj_nsites; i++) {
    printf(" %d",(*segmentj_sites)[i]);
    if ((*segmentj_types)[i] >= 0) {
      printf(" [#%d]",(*segmentj_types)[i]);
    } else if ((*segmentj_types)[i] == -1) {
      printf(" [GTGC]");
    } else if ((*segmentj_types)[i] == -2) {
      printf(" [AT]");
    }
    printf(" (%.6f)",Maxent_hr_antidonor_prob(segmentj_left + (*segmentj_sites)[i],chroffset));
  }
  printf("\n");
#endif

  return segmentj_nsites;
}


void
donor_dinucleotide (char *donor1, char *donor2, Univcoord_T left, int splice_querypos,
		    int querylength, bool plusp, bool sense_forward_p) {
  char donor1_alt, donor2_alt;

  if (plusp == sense_forward_p) {
    if (plusp == true) {
      *donor1 = Genome_get_char(&donor1_alt,left + splice_querypos);
      *donor2 = Genome_get_char(&donor2_alt,left + splice_querypos + 1);
    } else {
      *donor1 = Genome_get_char(&donor1_alt,left + (querylength - splice_querypos));
      *donor2 = Genome_get_char(&donor2_alt,left + (querylength - splice_querypos) + 1);
    }

  } else {
    if (plusp == true) {
      *donor1 = complCode[(int) Genome_get_char(&donor1_alt,left + splice_querypos - 1)];
      *donor2 = complCode[(int) Genome_get_char(&donor2_alt,left + splice_querypos - 2)];
    } else {
      *donor1 = complCode[(int) Genome_get_char(&donor1_alt,left + (querylength - splice_querypos) - 1)];
      *donor2 = complCode[(int) Genome_get_char(&donor2_alt,left + (querylength - splice_querypos) - 2)];
    }
  }

  return;
}

void
acceptor_dinucleotide (char *acceptor1, char *acceptor2, Univcoord_T left, int splice_querypos,
		       int querylength, bool plusp, bool sense_forward_p) {
  char acceptor1_alt, acceptor2_alt;

  if (plusp == sense_forward_p) {
    if (plusp == true) {
      *acceptor1 = Genome_get_char(&acceptor1_alt,left + splice_querypos - 1);
      *acceptor2 = Genome_get_char(&acceptor2_alt,left + splice_querypos - 2);
    } else {
      *acceptor1 = Genome_get_char(&acceptor1_alt,left + (querylength - splice_querypos) - 1);
      *acceptor2 = Genome_get_char(&acceptor2_alt,left + (querylength - splice_querypos) - 2);
    }

  } else {
    if (plusp == true) {
      *acceptor1 = complCode[(int) Genome_get_char(&acceptor1_alt,left + splice_querypos)];
      *acceptor2 = complCode[(int) Genome_get_char(&acceptor2_alt,left + splice_querypos + 1)];
    } else {
      *acceptor1 = complCode[(int) Genome_get_char(&acceptor1_alt,left + (querylength - splice_querypos))];
      *acceptor2 = complCode[(int) Genome_get_char(&acceptor2_alt,left + (querylength - splice_querypos) + 1)];
    }
  }

  return;
}


#ifndef CHECK_ASSERTIONS
static inline void
check_ascending (int *positions, int n) {
  return;
}

static inline void
check_descending (int *positions, int n) {
  return;
}

#else

static void
check_ascending (int *positions, int n) {
  int prevpos;
  int i;

  prevpos = positions[0];
  for (i = 1; i < n; i++) {
    if (positions[i] <= prevpos) {
      printf("Expecting ascending, but at %d, got %d <= %d\n",
	     i,positions[i],prevpos);
      abort();
    }
    prevpos = positions[i];
  }
 
  return;
}

static void
check_descending (int *positions, int n) {
  int prevpos;
  int i;

  prevpos = positions[0];
  for (i = 1; i < n; i++) {
    if (positions[i] >= prevpos) {
      printf("Expecting descending, but at %d, got %d >= %d\n",
	     i,positions[i],prevpos);
      abort();
    }
    prevpos = positions[i];
  }
 
  return;
}

#endif



/* All sites and positions are on querypos coordinates, not qpos */
/* segmentD_sites and segmentA_sites are ascending */
/* mismatch_positions_donor are ascending, but mismatch_positions_acceptor are descending */
static int
splice_sense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
	      double *best_donor_prob, double *best_acceptor_prob,

	      int *best_nmismatches_D, int *best_nmismatches_A,
	      int *best_ref_nmismatches_D, int *best_ref_nmismatches_A,

	      Univcoord_T univdiagonal_D, Univcoord_T univdiagonal_A,
	      Univcoord_T chroffset_D, Univcoord_T chroffset_A,
	      int plusDp, int plusAp, int querylength,
	      int *mismatch_positions_donor, int nmismatches_donor,
	      int *mismatch_positions_acceptor, int nmismatches_acceptor,
	      
	      int *segmentD_sites, int *segmentA_sites,
	      int *segmentD_types, int *segmentA_types,
	      double *segmentD_probs, double *segmentA_probs,
	      int segmentD_nsites, int segmentA_nsites) {

  int best_splice_querypos = -1, splice_querypos, i, j;
  int type;

  int best_nmismatches = querylength, nmismatches,
    segmentD_nmismatches, segmentA_nmismatches;
  double best_prob = 0.0, donor_prob, acceptor_prob;
  /* double best_probi, best_probj; */

  /* int adj_supporti, adj_supportj, supporti, supportj; */
  Univcoord_T segmentD_left, segmentA_left;


  check_ascending(segmentD_sites,segmentD_nsites);
  check_ascending(segmentA_sites,segmentA_nsites);
  check_ascending(mismatch_positions_donor,nmismatches_donor);
  check_descending(mismatch_positions_acceptor,nmismatches_acceptor);


  debug1(printf("Entered splice_sense with %d donor sites and %d acceptor sites\n",
		segmentD_nsites,segmentA_nsites));
  segmentD_left = univdiagonal_D - querylength;
  segmentA_left = univdiagonal_A - querylength;

  *best_donor_prob = *best_acceptor_prob = 0.0;

  segmentD_nmismatches = 0;
  segmentA_nmismatches = nmismatches_acceptor;
  i = j = 0;

  while (i < segmentD_nsites && j < segmentA_nsites) {
    debug1(printf("i [%d/%d]: %d (%f) and j [%d/%d]: %d (%f)\n",
		  i,segmentD_nsites,segmentD_sites[i],segmentD_probs[i],
		  j,segmentA_nsites,segmentA_sites[j],segmentA_probs[j]));
    if ((splice_querypos = segmentD_sites[i]) < segmentA_sites[j]) {
      i++;
      
    } else if (splice_querypos > segmentA_sites[j]) {
      j++;
      
    } else if ((segmentD_types[i] < TYPE1_SITE || segmentA_types[j] < TYPE1_SITE) &&
	       segmentD_types[i] != segmentA_types[j]) {
      debug1(printf("types %d and %d not compatible\n",segmentD_types[i],segmentA_types[j]));
      i++; j++;

    } else {
      debug1(printf("splice matches at %d\n",splice_querypos));
      while (segmentD_nmismatches <= nmismatches_donor && mismatch_positions_donor[segmentD_nmismatches] < splice_querypos) {
	segmentD_nmismatches++;
      }
      while (segmentA_nmismatches - 1 >= 0 && mismatch_positions_acceptor[segmentA_nmismatches - 1] < splice_querypos) {
	segmentA_nmismatches--;
      }
      
      debug1(printf("Evaluating plus splice querypos %d\n",splice_querypos));
      debug1(printf("%d mismatches on segmentD (..%d)\n",segmentD_nmismatches,splice_querypos));
      debug1(printf("%d mismatches on segmentA (%d..)\n",segmentA_nmismatches,splice_querypos));
      
      if (segmentD_types[i] >= 0) {
	type = KNOWN_SITE;
      } else if (segmentD_types[i] == -1) {
	type = TYPE1_SITE;
      } else if (segmentD_types[i] == -2) {
	type = TYPE2_SITE;
      }
      if (segmentA_types[j] >= 0) {
	type = KNOWN_SITE;
      } else if (segmentA_types[j] == -1) {
	type = TYPE1_SITE;
      } else if (segmentA_types[j] == -2) {
	type = TYPE2_SITE;
      }

      donor_prob = segmentD_probs[i];
      acceptor_prob = segmentA_probs[j];

      
      if ((nmismatches = segmentD_nmismatches + segmentA_nmismatches) < best_nmismatches) {
	best_nmismatches = nmismatches;
	best_splice_querypos = splice_querypos;
	*best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	*best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	*best_donor_prob = donor_prob;
	*best_acceptor_prob = acceptor_prob;
	best_prob = donor_prob + acceptor_prob;
	*best_type = type;
	
      } else if (nmismatches == best_nmismatches) {
	if (donor_prob + acceptor_prob > best_prob) {
	  best_splice_querypos = splice_querypos;
	  *best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	  *best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	  *best_donor_prob = donor_prob;
	  *best_acceptor_prob = acceptor_prob;
	  best_prob = donor_prob + acceptor_prob;
	*best_type = type;
	}
      }
      i++; j++;
    }
  }

  if (best_splice_querypos >= 0) {
    donor_dinucleotide(&(*donor1),&(*donor2),segmentD_left,best_splice_querypos,
		       querylength,plusDp,/*sense_forward_p*/true);
    acceptor_dinucleotide(&(*acceptor1),&(*acceptor2),segmentA_left,best_splice_querypos,
			  querylength,plusAp,/*sense_forward_p*/true);
  }

  debug1(printf("splice_sense returning %d (%f..%f)\n",
		best_splice_querypos,*best_donor_prob,*best_acceptor_prob));
  return best_splice_querypos;
}


/* All sites and positions are on querypos coordinates, not qpos */
/* segmentA_sites and segmentD_sites are ascending */
/* mismatch_positions_acceptor are ascending, but mismatch_positions_donor are descending */
static int
splice_antisense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		  double *best_donor_prob, double *best_acceptor_prob,		      

		  int *best_nmismatches_A, int *best_nmismatches_D,
		  int *best_ref_nmismatches_A, int *best_ref_nmismatches_D,
		  Univcoord_T univdiagonal_A, Univcoord_T univdiagonal_D,
		  Univcoord_T chroffset_A, Univcoord_T chroffset_D,
		  bool plusAp, bool plusDp, int querylength,

		  int *mismatch_positions_acceptor, int nmismatches_acceptor,
		  int *mismatch_positions_donor, int nmismatches_donor,

		  int *segmentA_sites, int *segmentD_sites,
		  int *segmentA_types, int *segmentD_types,
		  double *segmentA_probs, double *segmentD_probs,
		  int segmentA_nsites, int segmentD_nsites) {

  int best_splice_querypos = -1, splice_querypos, i, j;
  int type;

  int best_nmismatches = querylength, nmismatches,
    segmentA_nmismatches, segmentD_nmismatches;
  double best_prob = 0.0, donor_prob, acceptor_prob;
  /* double best_probi, best_probj; */

  /* int adj_supporti, adj_supportj, supporti, supportj; */
  Univcoord_T segmentA_left, segmentD_left;
  

  check_ascending(segmentA_sites,segmentA_nsites);
  check_ascending(segmentD_sites,segmentD_nsites);
  check_ascending(mismatch_positions_acceptor,nmismatches_acceptor);
  check_descending(mismatch_positions_donor,nmismatches_donor);


  debug1(printf("Entered splice_antisense\n"));
  segmentA_left = univdiagonal_A - querylength;
  segmentD_left = univdiagonal_D - querylength;

  *best_donor_prob = *best_acceptor_prob = 0.0;

  /* Match up sites */
  segmentA_nmismatches = 0;
  segmentD_nmismatches = nmismatches_donor;
  i = j = 0;

  while (i < segmentA_nsites && j < segmentD_nsites) {
    debug1(printf("i [%d/%d]: %d (%f) and j [%d/%d]: %d (%f)\n",
		  i,segmentA_nsites,segmentA_sites[i],segmentA_probs[i],
		  j,segmentD_nsites,segmentD_sites[j],segmentD_probs[j]));
    if ((splice_querypos = segmentA_sites[i]) < segmentD_sites[j]) {
      i++;
      
    } else if (splice_querypos > segmentD_sites[j]) {
      j++;
      
    } else if ((segmentA_types[i] < TYPE1_SITE || segmentD_types[j] < TYPE1_SITE) &&
	       segmentA_types[i] != segmentD_types[j]) {
      debug1(printf("types %d and %d not compatible\n",segmentA_types[i],segmentD_types[j]));
      i++; j++;

    } else {
      debug1(printf("splice matches at %d\n",splice_querypos));
      while (segmentA_nmismatches <= nmismatches_acceptor && mismatch_positions_acceptor[segmentA_nmismatches] < splice_querypos) {
	segmentA_nmismatches++;
      }
      while (segmentD_nmismatches - 1 >= 0 && mismatch_positions_donor[segmentD_nmismatches - 1] < splice_querypos) {
	segmentD_nmismatches--;
      }
      
      debug1(printf("Evaluating splice querypos %d\n",splice_querypos));
      debug1(printf("%d mismatches on segmentA (..%d)\n",segmentA_nmismatches,splice_querypos));
      debug1(printf("%d mismatches on segmentD (%d..)\n",segmentD_nmismatches,splice_querypos));
      
      if (segmentA_types[i] >= 0) {
	type = KNOWN_SITE;
      } else if (segmentA_types[i] == -1) {
	type = TYPE1_SITE;
      } else if (segmentA_types[i] == -2) {
	type = TYPE2_SITE;
      }
      if (segmentD_types[j] >= 0) {
	type = KNOWN_SITE;
      } else if (segmentD_types[j] == -1) {
	type = TYPE1_SITE;
      } else if (segmentD_types[j] == -2) {
	type = TYPE2_SITE;
      }

      acceptor_prob = segmentA_probs[i];
      donor_prob = segmentD_probs[j];

      if ((nmismatches = segmentA_nmismatches + segmentD_nmismatches) < best_nmismatches) {
	best_nmismatches = nmismatches;
	best_splice_querypos = splice_querypos;
	*best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	*best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	*best_donor_prob = donor_prob;
	*best_acceptor_prob = acceptor_prob;
	best_prob = donor_prob + acceptor_prob;
	*best_type = type;
	
      } else if (nmismatches == best_nmismatches) {
	if (donor_prob + acceptor_prob > best_prob) {
	  best_splice_querypos = splice_querypos;
	  *best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	  *best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	  *best_donor_prob = donor_prob;
	  *best_acceptor_prob = acceptor_prob;
	  best_prob = donor_prob + acceptor_prob;
	  *best_type = type;
	}
      }
      i++; j++;
    }
  }

  if (best_splice_querypos >= 0) {
    acceptor_dinucleotide(&(*acceptor1),&(*acceptor2),segmentA_left,best_splice_querypos,
			  querylength,plusAp,/*sense_forward_p*/false);
    donor_dinucleotide(&(*donor1),&(*donor2),segmentD_left,best_splice_querypos,
		       querylength,plusDp,/*sense_forward_p*/false);
  }

  debug1(printf("splice_antisense returning %d (%f..%f)\n",
		best_splice_querypos,*best_donor_prob,*best_acceptor_prob));
  return best_splice_querypos;
}


/* Looks primarily at support, then splice prob */
/* mismatch_positions_donor are ascending, but mismatch_positions_acceptor are descending */

static int
atypical_sense (char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		double *best_donor_prob, double *best_acceptor_prob,

		int *best_nmismatches_D, int *best_nmismatches_A,
		int *best_ref_nmismatches_D, int *best_ref_nmismatches_A,
		 
		Univcoord_T univdiagonal_D, Univcoord_T univdiagonal_A,
		Univcoord_T chroffset_D, Univcoord_T chroffset_A,
		bool plusDp, bool plusAp, int querylength,
		int *mismatch_positions_donor, int nmismatches_donor,
		int *mismatch_positions_acceptor, int nmismatches_acceptor,

		int splice_querypos_low, int splice_querypos_high) {

  int donori_save, acceptori_save;
  int best_nmismatches, nmismatches, segmentD_nmismatches, segmentA_nmismatches;
  
  int best_splice_querypos, splice_querypos;

  double best_prob, prob, donor_prob, acceptor_prob;
  Univcoord_T segmentD_left, segmentA_left;


  debug1(printf("Entered atypical_sense\n"));
  /* return -1; */

  check_ascending(mismatch_positions_donor,nmismatches_donor);
  check_descending(mismatch_positions_acceptor,nmismatches_acceptor);

  segmentD_left = univdiagonal_D - querylength;
  segmentA_left = univdiagonal_A - querylength;

  /* Could have splice_querypos_low == splice_querypos_high in a short region */
  /* assert(splice_querypos_low < splice_querypos_high); */
  if (splice_querypos_low >= splice_querypos_high) {
    return -1;
  }


  /* Find best_nmismatches */
  segmentD_nmismatches = 0;
  while (segmentD_nmismatches < nmismatches_donor &&
	 mismatch_positions_donor[segmentD_nmismatches] < splice_querypos_low) {
    segmentD_nmismatches++;
  }
  donori_save = segmentD_nmismatches;

  segmentA_nmismatches = 0;
  while (segmentA_nmismatches < nmismatches_acceptor &&
	 mismatch_positions_acceptor[segmentA_nmismatches] >= splice_querypos_low) {
    segmentA_nmismatches++;
  }
  acceptori_save = segmentA_nmismatches;

  best_nmismatches = nmismatches_donor + nmismatches_acceptor;

  for (splice_querypos = splice_querypos_low; splice_querypos < splice_querypos_high; splice_querypos++) {
    if (segmentD_nmismatches < nmismatches_donor &&
	mismatch_positions_donor[segmentD_nmismatches] < splice_querypos) {
      segmentD_nmismatches++;
    }
    if (segmentA_nmismatches - 1 >= 0 &&
	mismatch_positions_acceptor[segmentA_nmismatches - 1] < splice_querypos) {
      segmentA_nmismatches--;
    }
#if 1
    debug1(printf("At splice_querypos %d, have %d + %d nmismatches\n",
		  splice_querypos,segmentD_nmismatches,segmentA_nmismatches));
#endif
    if ((nmismatches = segmentD_nmismatches + segmentA_nmismatches) < best_nmismatches) {
      best_nmismatches = nmismatches;
    }
  }
  best_nmismatches += 1;	/* Allow for mismatches */


  /* Find best splice probs within positions with best_nmismatches or fewer */
  segmentD_nmismatches = donori_save;
  segmentA_nmismatches = acceptori_save;
  best_prob = 0.0;

  for (splice_querypos = splice_querypos_low; splice_querypos < splice_querypos_high; splice_querypos++) {
    if (segmentD_nmismatches < nmismatches_donor &&
	mismatch_positions_donor[segmentD_nmismatches] < splice_querypos) {
      segmentD_nmismatches++;
    }
    if (segmentA_nmismatches - 1 >= 0 &&
	mismatch_positions_acceptor[segmentA_nmismatches - 1] < splice_querypos) {
      segmentA_nmismatches--;
    }
    if (segmentD_nmismatches + segmentA_nmismatches <= best_nmismatches) {
      if (plusDp == /*sense_forward_p*/true) {
	donor_prob = Maxent_hr_donor_prob(segmentD_left + splice_querypos,chroffset_D);
      } else {
	/* Appears to be correct */
	donor_prob = Maxent_hr_antidonor_prob(segmentD_left + (querylength - splice_querypos),chroffset_D);
      }
      if (plusAp == /*sense_forward_p*/true) {
	acceptor_prob = Maxent_hr_acceptor_prob(segmentA_left + splice_querypos,chroffset_A);
      } else {
	/* Appears to be correct */
	acceptor_prob = Maxent_hr_antiacceptor_prob(segmentA_left + (querylength - splice_querypos),chroffset_A);
      }
      debug1(printf("  At splice_querypos %d, probs are %f donor and %f acceptor, plusp %d/%d, sense\n",
		    splice_querypos,donor_prob,acceptor_prob,plusDp,plusAp));
      if ((prob = donor_prob + acceptor_prob) > best_prob) {
	best_splice_querypos = splice_querypos;
	*best_donor_prob = donor_prob;
	*best_acceptor_prob = acceptor_prob;
	*best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	*best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	best_prob = prob;
      }
    }
  }

#if 0
  if (trim5p == false) {
    *trimpos5 = pos5;
  } else {
    *trimpos5 = Genomebits_trim_qstart(&(*best_nmismatches_i),query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/best_splice_querypos,plusp,/*genestrand*/0);
  }

  if (trim3p == false) {
    *trimpos3 = pos3;
  } else {
    *trimpos3 = Genomebits_trim_qend(&(*best_nmismatches_j),query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/best_splice_querypos,pos3,plusp,/*genestrand*/0);
  }
#endif


  if (best_splice_querypos >= 0) {
    donor_dinucleotide(&(*donor1),&(*donor2),segmentD_left,best_splice_querypos,
		       querylength,plusDp,/*sense_forward_p*/true);
    acceptor_dinucleotide(&(*acceptor1),&(*acceptor2),segmentA_left,best_splice_querypos,
			  querylength,plusAp,/*sense_forward_p*/true);
  }

  debug1(printf("atypical_sense returning splice_querypos %d, with %d + %d nmismatches\n",
		best_splice_querypos,*best_nmismatches_D,*best_nmismatches_A));
  return best_splice_querypos;
}


/* segmentA_positions and segmentD_positions are ascending */
/* mismatch_positions_acceptor are ascending, but mismatch_positions_donor are descending */
static int
atypical_antisense (char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		    double *best_donor_prob, double *best_acceptor_prob,

		    int *best_nmismatches_A, int *best_nmismatches_D,
		    int *best_ref_nmismatches_A, int *best_ref_nmismatches_D,
		 
		    Univcoord_T univdiagonal_A, Univcoord_T univdiagonal_D,
		    Univcoord_T chroffset_A, Univcoord_T chroffset_D,
		    bool plusAp, bool plusDp, int querylength,
		    
		    int *mismatch_positions_acceptor, int nmismatches_acceptor,
		    int *mismatch_positions_donor, int nmismatches_donor,
		    
		    int splice_querypos_low, int splice_querypos_high) {
  
  int acceptori_save, donori_save;
  int best_nmismatches, nmismatches, segmentA_nmismatches, segmentD_nmismatches;
  
  int best_splice_querypos, splice_querypos;

  double best_prob, prob, donor_prob, acceptor_prob;
  Univcoord_T segmentA_left, segmentD_left;


  debug1(printf("Entered atypical_antisense\n"));
  /* return -1; */

  check_ascending(mismatch_positions_acceptor,nmismatches_acceptor);
  check_descending(mismatch_positions_donor,nmismatches_donor);

  segmentA_left = univdiagonal_A - querylength;
  segmentD_left = univdiagonal_D - querylength;

  /* Could have splice_querypos_low == splice_querypos_high in a short region */
  /* assert(splice_querypos_low < splice_querypos_high); */
  if (splice_querypos_low >= splice_querypos_high) {
    return -1;
  }


  /* Find best_nmismatches */
  segmentA_nmismatches = 0;
  while (segmentA_nmismatches < nmismatches_acceptor &&
	 mismatch_positions_acceptor[segmentA_nmismatches] < splice_querypos_low) {
    segmentA_nmismatches++;
  }
  acceptori_save = segmentA_nmismatches;

  segmentD_nmismatches = 0;
  while (segmentD_nmismatches < nmismatches_donor &&
	 mismatch_positions_donor[segmentD_nmismatches] >= splice_querypos_low) {
    segmentD_nmismatches++;
  }
  donori_save = segmentD_nmismatches;

  best_nmismatches = nmismatches_acceptor + nmismatches_donor;

  for (splice_querypos = splice_querypos_low; splice_querypos < splice_querypos_high; splice_querypos++) {
    if (segmentA_nmismatches < nmismatches_acceptor &&
	mismatch_positions_acceptor[segmentA_nmismatches] < splice_querypos) {
      segmentA_nmismatches++;
    }
    if (segmentD_nmismatches - 1 >= 0 &&
	mismatch_positions_donor[segmentD_nmismatches - 1] < splice_querypos) {
      segmentD_nmismatches--;
    }
#if 1
    debug1(printf("At splice_querypos %d, have %d + %d nmismatches\n",
		  splice_querypos,segmentA_nmismatches,segmentD_nmismatches));
#endif
    if ((nmismatches = segmentA_nmismatches + segmentD_nmismatches) < best_nmismatches) {
      best_nmismatches = nmismatches;
    }
  }
  best_nmismatches += 1;	/* Allow for mismatches */


  /* Find best splice probs within positions with best_nmismatches or fewer */
  segmentA_nmismatches = acceptori_save;
  segmentD_nmismatches = donori_save;
  best_prob = 0.0;

  for (splice_querypos = splice_querypos_low; splice_querypos < splice_querypos_high; splice_querypos++) {
    if (segmentA_nmismatches < nmismatches_acceptor &&
	mismatch_positions_acceptor[segmentA_nmismatches] < splice_querypos) {
	segmentA_nmismatches++;
    }
    if (segmentD_nmismatches - 1 >= 0 &&
	mismatch_positions_donor[segmentD_nmismatches - 1] < splice_querypos) {
      segmentD_nmismatches--;
    }
    if (segmentA_nmismatches + segmentD_nmismatches <= best_nmismatches) {
      if (plusAp == /*sense_forward_p*/false) {
	/* Appears to be correct */
	acceptor_prob = Maxent_hr_acceptor_prob(segmentA_left + (querylength - splice_querypos),chroffset_A);
      } else {
	acceptor_prob = Maxent_hr_antiacceptor_prob(segmentA_left + splice_querypos,chroffset_A);
      }
      if (plusDp == /*sense_forward_p*/false) {
	/* Appears to be correct */
	donor_prob = Maxent_hr_donor_prob(segmentD_left + (querylength - splice_querypos),chroffset_D);
      } else {
	donor_prob = Maxent_hr_antidonor_prob(segmentD_left + splice_querypos,chroffset_D);
      }
      debug1(printf("  At splice_querypos %d, probs are %f acceptor and %f donor, plusp %d/%d, antisense\n",
		    splice_querypos,acceptor_prob,donor_prob,plusAp,plusDp));
      if ((prob = acceptor_prob + donor_prob) > best_prob) {
	best_splice_querypos = splice_querypos;
	*best_acceptor_prob = acceptor_prob;
	*best_donor_prob = donor_prob;
	*best_nmismatches_A = *best_ref_nmismatches_A = segmentA_nmismatches;
	*best_nmismatches_D = *best_ref_nmismatches_D = segmentD_nmismatches;
	best_prob = prob;
      }
    }
  }

#if 0
  if (trim5p == false) {
    *trimpos5 = pos5;
  } else {
    *trimpos5 = Genomebits_trim_qstart(&(*best_nmismatches_i),query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/best_splice_querypos,plusp,/*genestrand*/0);
  }

  if (trim3p == false) {
    *trimpos3 = pos3;
  } else {
    *trimpos3 = Genomebits_trim_qend(&(*best_nmismatches_j),query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/best_splice_querypos,pos3,plusp,/*genestrand*/0);
  }
#endif

  if (best_splice_querypos >= 0) {
    acceptor_dinucleotide(&(*acceptor1),&(*acceptor2),segmentA_left,best_splice_querypos,
			  querylength,plusAp,/*sense_forward_p*/false);
    donor_dinucleotide(&(*donor1),&(*donor2),segmentD_left,best_splice_querypos,
		       querylength,plusDp,/*sense_forward_p*/false);
  }

  debug1(printf("atypical_antisense returning splice_qpos %d, with %d + %d nmismatches\n",
		best_splice_querypos,*best_nmismatches_A,*best_nmismatches_D));
  return best_splice_querypos;
}



/* Use querylength to invert sites, so splice occurs after splice_querypos */
static void
invert_sites (int *sites, int *types, double *probs, int n, int querylength) {
  int i, j;
  int temp;
  double temp_prob;

  for (i = 0, j = n - 1; i < n/2; i++, j--) {
    temp = querylength - sites[i];
    sites[i] = querylength - sites[j];
    sites[j] = temp;

    temp = types[i];
    types[i] = types[j];
    types[j] = temp;

    temp_prob = probs[i];
    probs[i] = probs[j];
    probs[j] = temp_prob;
  }
  if (i == j) {
    sites[i] = querylength - sites[i];
  }

  return;
}


/* Use (querylength - 1) to invert positions so they mirror the original positions */
/* mismatch_positions have entries from 0 through n */
static void
invert_mismatch_positions (int *positions, int n, int querylength) {
  int i;

  for (i = 0; i <= n; i++) {
    positions[i] = (querylength - 1) - positions[i];
  }

  return;
}

#if defined(DEBUG1) || defined(DEBUG2)
static void
print_diffs (char *string1, char *string2, int querylength) {
  int i;

  for (i = 0; i < querylength; i++) {
    if (string1[i] == string2[i]) {
      printf("|");
    } else {
      printf(" ");
    }
  }
}
#endif


static Univcoord_T
find_middle_exon (int *best_type_i, int *best_type_j, int *best_splice_qpos_i, int *best_splice_qpos_j,
		  int *best_nmismatches_i, int *best_nmismatches_middle, int *best_nmismatches_j,
		  int *best_ref_nmismatches_i, int *best_ref_nmismatches_middle, int *best_ref_nmismatches_j,
		  double *best_donor1_prob, double *best_acceptor1_prob,
		  double *best_donor2_prob, double *best_acceptor2_prob,

		  Univcoord_T *middle_univdiagonals, int n_middle_univdiagonals,

		  Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		  int pos5, int pos3, int splice_qpos_5, int splice_qpos_3, int querylength,

		  int *mismatch_positions_left1, int nmismatches_left1,
		  int *mismatch_positions_right2, int nmismatches_right2,

		  int *segmenti_sites, int *segmentj_sites,
		  int *segmenti_types, int *segmentj_types,
		  double *segmenti_probs, double *segmentj_probs,
		  int segmenti_nsites, int segmentj_nsites,

		  Compress_T query_compress, Univcoord_T chroffset,

		  Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing, Intlistpool_T intlistpool,
		  bool plusp, bool sense_forward_p, int genestrand) {

  Univcoord_T best_middle_univdiagonal = 0, middle_univdiagonal;
  int typei, typej, type1, type2;
  int splice_qpos_i, splice_qpos_j, splice_querypos_1, splice_querypos_2;
  int best_nmismatches = querylength, nmismatches,
    nmismatches_i, nmismatches_middle, nmismatches_j,
    ref_nmismatches_i, ref_nmismatches_middle, ref_nmismatches_j;
  int ignore_nmismatches, ignore_ref_nmismatches;

  int *mismatch_positions_right1, *mismatch_positions_left2;
  int nmismatches_right1, nmismatches_left2;

  int *segmentk1_sites, *segmentk2_sites;
  int *segmentk1_types, *segmentk2_types;
  double *segmentk1_probs, *segmentk2_probs;
  int segmentk1_nsites, segmentk2_nsites;

  char donor1, donor2, acceptor1, acceptor2;
  double donor1_prob, acceptor1_prob, donor2_prob, acceptor2_prob,
    optimal_acceptor1_prob, optimal_donor2_prob;
  int i;
  

  return 0;

  debug1(printf("Entering find_middle_exon with pos5 %d, splice_qpos_5 %d, splice_qpos_3 %d, pos3 %d, plusp %d, sense_forward_p %d\n",
		pos5,splice_qpos_5,splice_qpos_3,pos3,plusp,sense_forward_p));

  mismatch_positions_right1 = spliceinfo->mismatch_positions_right1; /* Use allocated memory */
  mismatch_positions_left2 = spliceinfo->mismatch_positions_left2; /* Use allocated memory */

  for (i = 0; i < n_middle_univdiagonals; i++) {
    middle_univdiagonal = middle_univdiagonals[i];
    debug1(printf("Testing middle_univdiagonal %u\n",middle_univdiagonal));

#ifdef DEBUG1
    char *gbuffer = (char *) CALLOC(querylength+1,sizeof(char));
    Genome_fill_buffer(middle_univdiagonal - querylength,querylength,gbuffer);
    char *queryseq = Compress_queryseq(query_compress,querylength);
    
    printf("gm: %s\n",gbuffer);
    printf("    ");
    print_diffs(gbuffer,queryseq,querylength);
    printf("\n");
    printf("q:  %s\n",queryseq);
    FREE(queryseq);
    FREE(gbuffer);
#endif

    nmismatches_right1 =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_right1,/*max_mismatches*/pos3 - pos5,
					       genomebits,genomebits_alt,query_compress,
					       middle_univdiagonal,querylength,
					       pos5,/*pos3*/splice_qpos_3,plusp,genestrand);
    
    nmismatches_left2 =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left2,/*max_mismatches*/pos3 - pos5,
					      genomebits,genomebits_alt,query_compress,
					      middle_univdiagonal,querylength,
					      /*pos5*/splice_qpos_5,pos3,plusp,genestrand);

    if (plusp == sense_forward_p) {
      /* geneplus */
      segmentk1_nsites =
	compute_acceptor_sites(&optimal_acceptor1_prob,&segmentk1_sites,&segmentk1_types,&segmentk1_probs,
			       spliceinfo->segmentk1_sites_alloc1,spliceinfo->segmentk1_types_alloc1,
			       spliceinfo->segmentk1_sites_alloc2,spliceinfo->segmentk1_types_alloc2,
			       spliceinfo->segmentk1_probs_alloc,
			       /*splice_qpos_low*/subtract_bounded(splice_qpos_5,EXTRA_MATCHES,pos5 + 1),
			       /*splice_qpos_high*/add_bounded(splice_qpos_5,EXTRA_MISMATCHES,splice_qpos_3 - 1),
			       querylength,middle_univdiagonal,chroffset,knownsplicing,intlistpool);
      segmentk2_nsites =
	compute_donor_sites(&optimal_donor2_prob,&segmentk2_sites,&segmentk2_types,&segmentk2_probs,
			    spliceinfo->segmentk2_sites_alloc1,spliceinfo->segmentk2_types_alloc1,
			    spliceinfo->segmentk2_sites_alloc2,spliceinfo->segmentk2_types_alloc2,
			    spliceinfo->segmentk2_probs_alloc,
			    /*splice_qpos_low*/subtract_bounded(splice_qpos_3,EXTRA_MISMATCHES,splice_qpos_5 + 1),
			    /*splice_qpos_high*/add_bounded(splice_qpos_3,EXTRA_MATCHES,pos3 - 1),
			    querylength,middle_univdiagonal,chroffset,knownsplicing,intlistpool);
    } else {
      /* geneminus */
      segmentk1_nsites =
	compute_antidonor_sites(&optimal_donor2_prob,&segmentk1_sites,&segmentk1_types,&segmentk1_probs,
				spliceinfo->segmentk1_sites_alloc1,spliceinfo->segmentk1_types_alloc1,
				spliceinfo->segmentk1_sites_alloc2,spliceinfo->segmentk1_types_alloc2,
				spliceinfo->segmentk1_probs_alloc,
				/*splice_qpos_low*/subtract_bounded(splice_qpos_5,EXTRA_MATCHES,pos5 + 1),
				/*splice_qpos_high*/add_bounded(splice_qpos_5,EXTRA_MISMATCHES,splice_qpos_3 - 1),
				querylength,middle_univdiagonal,chroffset,knownsplicing,intlistpool);
      segmentk2_nsites =
	compute_antiacceptor_sites(&optimal_acceptor1_prob,&segmentk2_sites,&segmentk2_types,&segmentk2_probs,
				   spliceinfo->segmentk2_sites_alloc1,spliceinfo->segmentk2_types_alloc1,
				   spliceinfo->segmentk2_sites_alloc2,spliceinfo->segmentk2_types_alloc2,
				   spliceinfo->segmentk2_probs_alloc,
				   /*splice_qpos_low*/subtract_bounded(splice_qpos_3,EXTRA_MISMATCHES,splice_qpos_5 + 1),
				   /*splice_qpos_high*/add_bounded(splice_qpos_3,EXTRA_MATCHES,pos3 - 1),
				   querylength,middle_univdiagonal,chroffset,knownsplicing,intlistpool);
    }
    if (plusp == false) {
      /* Already should have inverted mismatch_positions_left1 and
	 mismatch_positions_right2, as well as segmenti_sites and
	 segmentj_sites.  Here we are inverting for the middle exon. */
      invert_mismatch_positions(mismatch_positions_right1,nmismatches_right1,querylength);
      invert_mismatch_positions(mismatch_positions_left2,nmismatches_left2,querylength);
      invert_sites(segmentk1_sites,segmentk1_types,segmentk1_probs,segmentk1_nsites,querylength);
      invert_sites(segmentk2_sites,segmentk2_types,segmentk2_probs,segmentk2_nsites,querylength);
    }

    if (sense_forward_p == true) {
      if (plusp == true) {
	/* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
	if ((splice_querypos_1 =
	     splice_sense(&type1,&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			  /*D*/&nmismatches_i,/*A*/&ignore_nmismatches,
			  /*D*/&ref_nmismatches_i,/*A*/&ignore_ref_nmismatches,
			  /*D*/univdiagonal_i,/*A*/middle_univdiagonal,/*D*/chroffset,/*A*/chroffset,
			  /*plusDp*/true,/*plusAp*/true,querylength,
			  /*donor*/mismatch_positions_left1,/*donor*/nmismatches_left1,
			  /*acceptor*/mismatch_positions_right1,/*acceptor*/nmismatches_right1,
			  /*segmentD_sites*/segmenti_sites,/*segmentA_sites*/segmentk1_sites,
			  /*segmentD_types*/segmenti_types,/*segmentA_types*/segmentk1_types,
			  /*segmentD_probs*/segmenti_probs,/*segmentA_probs*/segmentk1_probs,
			  /*segmentD_nsites*/segmenti_nsites,/*segmentA_nsites*/segmentk1_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_1 =
	    atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			   /*D*/&nmismatches_i,/*A*/&ignore_nmismatches,
			   /*D*/&ref_nmismatches_i,/*A*/&ignore_ref_nmismatches,
			   /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			   /*plusDp*/true,/*plusAp*/true,querylength,
			   /*donor*/mismatch_positions_left1,/*donor*/nmismatches_left1,
			   /*acceptor*/mismatch_positions_right1,/*acceptor*/nmismatches_right1,
			   /*splice_querypos_low*/pos5,/*splice_querypos_high*/splice_qpos_5);
#endif
	}
	
	if ((splice_querypos_2 =
	     splice_sense(&type2,&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			  /*D*/&ignore_nmismatches,/*A*/&nmismatches_j,
			  /*D*/&ignore_ref_nmismatches,/*A*/&ref_nmismatches_j,
			  /*D*/middle_univdiagonal,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			  /*plusDp*/true,/*plusAp*/true,querylength,
			  /*donor*/mismatch_positions_left2,/*donor*/nmismatches_left2,
			  /*acceptor*/mismatch_positions_right2,/*acceptor*/nmismatches_right2,
			  /*segmentD_sites*/segmentk2_sites,/*segmentA_sites*/segmentj_sites,
			  /*segmentD_types*/segmentk2_types,/*segmentA_types*/segmentj_types,
			  /*segmentD_probs*/segmentk2_probs,/*segmentA_probs*/segmentj_probs,
			  /*segmentD_nsites*/segmentk2_nsites,/*segmentA_nsites*/segmentj_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_2 = 
	    atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			   /*D*/&ignore_nmismatches,/*A*/&nmismatches_j,
			   /*D*/&ignore_ref_nmismatches,/*A*/&ref_nmismatches_j,
			   /*D*/middle_univdiagonal,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			   /*plusDp*/true,/*plusAp*/true,querylength,
			   /*donor*/mismatch_positions_left2,/*donor*/nmismatches_left2,
			   /*acceptor*/mismatch_positions_right2,/*acceptor*/nmismatches_right2,
			   /*splice_querypos_low*/splice_qpos_3,/*splice_querypos_high*/pos3);
#endif
	}

      } else {
	/* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	if ((splice_querypos_1 =
	     splice_sense(&type1,&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			  /*D*/&nmismatches_j,/*A*/&ignore_nmismatches,
			  /*D*/&ref_nmismatches_j,/*A*/&ignore_ref_nmismatches,
			  /*D*/univdiagonal_j,/*A*/middle_univdiagonal,/*D*/chroffset,/*A*/chroffset,
			  /*plusDp*/false,/*plusAp*/false,querylength,
			  /*donor*/mismatch_positions_right2,nmismatches_right2,
			  /*acceptor*/mismatch_positions_left2,/*acceptor*/nmismatches_left2,
			  /*segmentD_sites*/segmentj_sites,/*segmentA_sites*/segmentk2_sites,
			  /*segmentD_types*/segmentj_types,/*segmentA_types*/segmentk2_types,
			  /*segmentD_probs*/segmentj_probs,/*segmentA_probs*/segmentk2_probs,
			  /*segmentD_nsites*/segmentj_nsites,/*segmentA_nsites*/segmentk2_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_1 =
	    atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			   /*D*/&nmismatches_j,/*A*/&ignore_nmismatches,
			   /*D*/&ref_nmismatches_j,/*A*/&ignore_ref_nmismatches,
			   /*D*/univdiagonal_j,/*A*/middle_univdiagonal,/*D*/chroffset,/*A*/chroffset,
			   /*plusDp*/false,/*plusAp*/false,querylength,
			   /*donor*/mismatch_positions_right2,nmismatches_right2,
			   /*acceptor*/mismatch_positions_left2,/*acceptor*/nmismatches_left2,
			   /*splice_querypos_low*/querylength - pos3,
			   /*splice_querypos_high*/querylength - splice_qpos_3);
#endif
	}

	if ((splice_querypos_2 =
	     splice_sense(&type2,&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			  /*D*/&ignore_nmismatches,/*A*/&nmismatches_i,
			  /*D*/&ignore_ref_nmismatches,/*A*/&ref_nmismatches_i,
			  /*D*/middle_univdiagonal,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			  /*plusDp*/false,/*plusAp*/false,querylength,
			  /*donor*/mismatch_positions_right1,/*donor*/nmismatches_right1,
			  /*acceptor*/mismatch_positions_left1,/*acceptor*/nmismatches_left1,
			  /*segmentD_sites*/segmentk1_sites,/*segmentA_sites*/segmenti_sites,
			  /*segmentD_types*/segmentk1_types,/*segmentA_types*/segmenti_types,
			  /*segmentD_probs*/segmentk1_probs,/*segmentA_probs*/segmenti_probs,
			  /*segmentD_nsites*/segmentk1_nsites,/*segmentA_nsites*/segmenti_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_2 =
	    atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			   /*D*/&ignore_nmismatches,/*A*/&nmismatches_i,
			   /*D*/&ignore_ref_nmismatches,/*A*/&ref_nmismatches_i,
			   /*D*/middle_univdiagonal,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			   /*plusDp*/false,/*plusAp*/false,querylength,
			   /*donor*/mismatch_positions_right1,/*donor*/nmismatches_right1,
			   /*acceptor*/mismatch_positions_left1,/*acceptor*/nmismatches_left1,
			   /*splice_querypos_low*/querylength - pos5,
			   /*splice_querypos_high*/querylength - splice_qpos_5);
#endif
	}
      }

    } else {
      if (plusp == true) {
	/* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	if ((splice_querypos_1 =
	     splice_antisense(&type1,&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			      /*A*/&nmismatches_i,/*D*/&ignore_nmismatches,
			      /*A*/&ref_nmismatches_i,/*D*/&ignore_ref_nmismatches,
			      /*A*/univdiagonal_i,/*D*/middle_univdiagonal,/*A*/chroffset,/*D*/chroffset,
			      /*plusAp*/true,/*plusDp*/true,querylength,
			      /*acceptor*/mismatch_positions_left1,/*acceptor*/nmismatches_left1,
			      /*donor*/mismatch_positions_right1,/*donor*/nmismatches_right1,
			      /*segmentA_sites*/segmenti_sites,/*segmentD_sites*/segmentk1_sites,
			      /*segmentA_types*/segmenti_types,/*segmentD_types*/segmentk1_types,
			      /*segmentA_probs*/segmenti_probs,/*segmentD_probs*/segmentk1_probs,
			      /*segmentA_nsites*/segmenti_nsites,/*segmentD_nsites*/segmentk1_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_1 =
	    atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			       /*A*/&nmismatches_i,/*D*/&ignore_nmismatches,
			       /*A*/&ref_nmismatches_i,/*D*/&ignore_ref_nmismatches,
			       /*A*/univdiagonal_i,/*D*/middle_univdiagonal,/*A*/chroffset,/*D*/chroffset,
			       /*plusAp*/true,/*plusDp*/true,querylength,
			       /*donor*/mismatch_positions_left1,/*donor*/nmismatches_left1,
			       /*acceptor*/mismatch_positions_right1,/*acceptor*/nmismatches_right1,
			       /*splice_querypos_low*/pos5,/*splice_querypos_high*/splice_qpos_5);
#endif
	}
	
	if ((splice_querypos_2 =
	     splice_antisense(&type2,&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			      /*A*/&ignore_nmismatches,/*D*/&nmismatches_j,
			      /*A*/&ignore_ref_nmismatches,/*D*/&ref_nmismatches_j,
			      /*A*/middle_univdiagonal,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			      /*plusAp*/true,/*plusDp*/true,querylength,
			      /*acceptor*/mismatch_positions_left2,/*acceptor*/nmismatches_left2,
			      /*donor*/mismatch_positions_right2,/*donor*/nmismatches_right2,
			      /*segmentA_sites*/segmentk2_sites,/*segmentD_sites*/segmentj_sites,
			      /*segmentA_types*/segmentk2_types,/*segmentD_types*/segmentj_types,
			      /*segmentA_probs*/segmentk2_probs,/*segmentD_probs*/segmentj_probs,
			      /*segmentA_nsites*/segmentk2_nsites,/*segmentD_nsites*/segmentj_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_2 =
	    atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			       /*A*/&ignore_nmismatches,/*D*/&nmismatches_j,
			       /*A*/&ignore_ref_nmismatches,/*D*/&ref_nmismatches_j,
			       /*A*/middle_univdiagonal,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			       /*plusAp*/true,/*plusDp*/true,querylength,
			       /*donor*/mismatch_positions_left2,/*donor*/nmismatches_left2,
			       /*acceptor*/mismatch_positions_right2,/*acceptor*/nmismatches_right2,
			       /*splice_querypos_low*/splice_qpos_3,/*splice_querypos_high*/pos3);
#endif
	}

      } else {
	/* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
	if ((splice_querypos_1 =
	     splice_antisense(&type1,&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			      /*A*/&nmismatches_j,/*D*/&ignore_nmismatches,
			      /*A*/&ref_nmismatches_j,/*D*/&ignore_ref_nmismatches,
			      /*A*/univdiagonal_j,/*D*/middle_univdiagonal,/*A*/chroffset,/*D*/chroffset,
			      /*plusAp*/false,/*plusDp*/false,querylength,
			      /*acceptor*/mismatch_positions_right2,/*acceptor*/nmismatches_right2,
			      /*donor*/mismatch_positions_left2,/*donor*/nmismatches_left2,
			      /*segmentA_sites*/segmentj_sites,/*segmentD_sites*/segmentk2_sites,
			      /*segmentA_types*/segmentj_types,/*segmentD_types*/segmentk2_types,
			      /*segmentA_probs*/segmentj_probs,/*segmentD_probs*/segmentk2_probs,
			      /*segmentA_nsites*/segmentj_nsites,/*segmentD_nsites*/segmentk2_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_1 =
	    atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&donor1_prob,&acceptor1_prob,
			       /*A*/&nmismatches_j,/*D*/&ignore_nmismatches,
			       /*A*/&ref_nmismatches_j,/*D*/&ignore_ref_nmismatches,
			       /*A*/univdiagonal_j,/*D*/middle_univdiagonal,/*A*/chroffset,/*D*/chroffset,
			       /*plusAp*/false,/*plusDp*/false,querylength,
			       /*donor*/mismatch_positions_right2,nmismatches_right2,
			       /*acceptor*/mismatch_positions_left2,/*acceptor*/nmismatches_left2,
			       /*splice_querypos_low*/querylength - pos3,
			       /*splice_querypos_high*/querylength - splice_qpos_3);
#endif
	}
	  
	if ((splice_querypos_2 =
	     splice_antisense(&type2,&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			      /*A*/&ignore_nmismatches,/*D*/&nmismatches_i,
			      /*A*/&ignore_ref_nmismatches,/*D*/&ref_nmismatches_i,
			      /*A*/middle_univdiagonal,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			      /*plusAp*/false,/*plusDp*/false,querylength,
			      /*acceptor*/mismatch_positions_right1,/*acceptor*/nmismatches_right1,
			      /*donor*/mismatch_positions_left1,/*donor*/nmismatches_left1,
			      /*segmentA_sites*/segmentk1_sites,/*segmentD_sites*/segmenti_sites,
			      /*segmentA_types*/segmentk1_types,/*segmentD_types*/segmenti_types,
			      /*segmentA_probs*/segmentk1_probs,/*segmentD_probs*/segmenti_probs,
			      /*segmentA_nsites*/segmentk1_nsites,/*segmentD_nsites*/segmenti_nsites)) <= 0) {
#ifdef ALLOW_ATYPICAL_MIDDLE
	  splice_querypos_2 =
	    atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&donor2_prob,&acceptor2_prob,
			       /*A*/&ignore_nmismatches,/*D*/&nmismatches_i,
			       /*A*/&ignore_ref_nmismatches,/*D*/&ref_nmismatches_i,
			       /*A*/middle_univdiagonal,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			       /*plusAp*/false,/*plusDp*/false,querylength,
			       /*donor*/mismatch_positions_right1,/*donor*/nmismatches_right1,
			       /*acceptor*/mismatch_positions_left1,/*acceptor*/nmismatches_left1,
			       /*splice_querypos_low*/querylength - splice_qpos_5,
			       /*splice_querypos_high*/querylength - pos5);
#endif
	}
      }
    }
    
    if (splice_querypos_1 <= 0 || splice_querypos_2 <= 0) {
      splice_qpos_i = -1;
      splice_qpos_j = -1;
    } else if (plusp == true) {
      splice_qpos_i = splice_querypos_1;
      splice_qpos_j = splice_querypos_2;
      typei = type1; typej = type2;
    } else {
      splice_qpos_i = querylength - splice_querypos_2;
      splice_qpos_j = querylength - splice_querypos_1;
      typei = type2; typej = type1;
    }

    debug1(printf("splice_querypos_1 %d, splice_querypos_2 %d => splice_qpos_i %d, splice_qpos_j %d\n",
		  splice_querypos_1,splice_querypos_2,splice_qpos_i,splice_qpos_j));

    if (splice_qpos_i >= 0 && splice_qpos_j >= 0 && splice_qpos_i < splice_qpos_j) {
      nmismatches_middle =
	Genomebits_count_mismatches_substring(&ref_nmismatches_middle,genomebits,genomebits_alt,
					      query_compress,/*univdiagonal*/middle_univdiagonal,querylength,
					      /*pos5*/splice_qpos_i,/*pos3*/splice_qpos_j,
					      plusp,genestrand);

      debug1(printf("For splices %d and %d, total_nmismatches: %d + %d + %d\n",
		    splice_qpos_i,splice_qpos_j,nmismatches_i,nmismatches_middle,nmismatches_j));
      if ((nmismatches = nmismatches_i + nmismatches_middle + nmismatches_j) < best_nmismatches) {
	best_middle_univdiagonal = middle_univdiagonal;
	*best_type_i = typei;
	*best_type_j = typej;

	*best_splice_qpos_i = splice_qpos_i;
	*best_splice_qpos_j = splice_qpos_j;
	*best_nmismatches_i = *best_ref_nmismatches_i = nmismatches_i;
	*best_nmismatches_middle = nmismatches_middle;
	*best_ref_nmismatches_middle = ref_nmismatches_middle;
	*best_nmismatches_j = *best_ref_nmismatches_j = nmismatches_j;
	*best_donor1_prob = donor1_prob;
	*best_acceptor1_prob = acceptor1_prob;
	*best_donor2_prob = donor2_prob;
	*best_acceptor2_prob = acceptor2_prob;
	best_nmismatches = nmismatches;
      }
    }
  }

#ifdef DEBUG1
  if (best_middle_univdiagonal == 0) {
    printf("No middle exon found\n");
  } else {
    printf("Best candidate middle exon is %d..%d for univdiagonal %u with %d segmenti, %d middle, and %d segmentj mismatches\n",
	   *best_splice_qpos_i,*best_splice_qpos_j,best_middle_univdiagonal,*best_nmismatches_i,*best_nmismatches_middle,*best_nmismatches_j);
  }
#endif

  return best_middle_univdiagonal;
}


/* Get better results if we recompute trimpos5 and trimpos3 based on the indel_pos and splice_qpos from spliceindel */
/* Insertions have negative values for nindels, and deletions positive, because
   we flip signs.  Path_solve procedures expect this convention. */
static int
spliceindel_resolve (int *best_type, int *trimpos5, int *trimpos3, int *best_nindels, int *best_indel_pos,
		     int *best_nmismatches_i, int *best_nmismatches_j, int *best_nmismatches_indel,
		     int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
		     int *best_ref_nmismatches_indel, double *best_donor_prob, double *best_acceptor_prob,
		     
		     Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		     Compress_T query_compress, bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		     int *mismatch_positions_left, int nmismatches_left,
		     int *mismatch_positions_right, int nmismatches_right,

		     int *segmenti_sites, int *segmentj_sites, int *segmenti_types, int *segmentj_types,
		     double *segmenti_probs, double *segmentj_probs, int segmenti_nsites, int segmentj_nsites,
		     double optimal_segmenti_prob, double optimal_segmentj_prob,
		     
		     int pos5, int pos3, int querylength,
		     Indelinfo_T indelinfo, bool sense_forward_p, int genestrand,
		     bool trim5p, bool trim3p) {

  int best_splice_qpos, splice_qpos, i, j;
  int best_nmismatches, nmismatches, nmismatches1, nmismatches2, ref_nmismatches1, ref_nmismatches2;

  int nindels, indel_pos, ninserts, max_deletionlen;
  int nmismatches_indel, nmismatches_i, nmismatches_j;
  double best_prob, prob, best_probi, best_probj, probi, probj;


  debug1(printf("Entered spliceindel_resolve with plusp %d, sense_forward_p %d, range %d..%d, optimal probs %f and %f\n",
		plusp,sense_forward_p,pos5,pos3,optimal_segmenti_prob,optimal_segmentj_prob));

  assert(segmenti_nsites > 0);
  assert(segmentj_nsites > 0);
  check_ascending(mismatch_positions_left,nmismatches_left);
  check_descending(mismatch_positions_right,nmismatches_right);

  /* segmenti_left = univdiagonal_i - querylength; */
  /* segmentj_left = univdiagonal_j - querylength; */


  *best_nindels = 0;
  *best_indel_pos = -1;
  best_splice_qpos = -1;
  best_prob = 0.0;
  best_nmismatches = querylength;  /* was nmismatches_allowed + 1 */


  debug1(printf("Comparing max_deletionlen %d with %d\n",MAX_DELETIONLEN,(int) (univdiagonal_j - univdiagonal_i)));
  if ((max_deletionlen = MAX_DELETIONLEN) > (int) (univdiagonal_j - univdiagonal_i)) {
    /* The sum of the deletion length and splice distance must be less than the distance in the genome */
    max_deletionlen = (int) (univdiagonal_j - univdiagonal_i);
  }

  /* Left anchor (splice) */
  /* Advancing segmenti splice site from low qpos to high qpos.  Same direction for segmentj */
  /* All splicesites run from low qpos to high qpos */
  debug1(printf(">Left anchor -- splice on segmenti and indel on segmentj:\n"));
  nmismatches_i = 0;
  i = 0;
  j = 0;

  /* Count mismatches, which are also from low qpos to high qpos */
  splice_qpos = segmenti_sites[i];
  while (nmismatches_i <= nmismatches_left && mismatch_positions_left[nmismatches_i] < splice_qpos) {
    nmismatches_i++;
  }

#if 0
  if (nmismatches_i > nmismatches_allowed) {
    printf("(1) nmismatches_i for splice_qpos %d is >= %d (exceeds nmismatches_allowed)\n",
	   splice_qpos,nmismatches_i);
  } else {
    printf("(1) nmismatches_i for splice_qpos %d is %d\n",splice_qpos,nmismatches_i);
  }
#endif
  
  /* mismatches_i does not hold for indels */
  while (i < segmenti_nsites /*&& nmismatches_i <= nmismatches_allowed*/) {
    probi = segmenti_probs[i];

    debug1(printf("Probi at %d is %f\n",splice_qpos,probi));

    if (probi >= SPLICEINDEL_PROB || probi == optimal_segmenti_prob) {
      /* Backup to low qpos, and then advance forward */
      while (j - 1 >= 0 && splice_qpos - segmentj_sites[j-1] <= max_deletionlen) {
	debug1(printf("Backing up j to %d because splice_qpos %d - %d < max_deletionlen %d\n",
		      j - 1,splice_qpos,segmentj_sites[j-1],max_deletionlen));
	j--;
      }

      /* Advance */
      while (j < segmentj_nsites && splice_qpos - segmentj_sites[j] > max_deletionlen) {
	j++;
      }

      /* Deletions on segmentj */
      while (j < segmentj_nsites && segmentj_sites[j] < splice_qpos) {
	assert(splice_qpos - segmentj_sites[j] <= max_deletionlen);

	if ((segmenti_types[i] < TYPE1_SITE || segmentj_types[j] < TYPE1_SITE) &&
	    segmenti_types[i] != segmentj_types[j]) {
	  /* Not compatible */
	} else {
	  probj = segmentj_probs[j];

	  debug1(printf("Deletion: probj at %d is %f\n",segmentj_sites[j],probj));
	  
	  if (probj >= SPLICEINDEL_PROB || probj == optimal_segmentj_prob) {
	    nindels = segmentj_sites[j] - splice_qpos; /* Deletions should be negative */
	    debug1(printf("univdiagonal_i %u, univdiagonal_j %u, nindels %d\n",
			  univdiagonal_i,univdiagonal_j,nindels));
	    /* debug1(printf("(1) Trying deletion on segmentj of %d from %d to %d\n",nindels,splice_qpos,segmentj_sites[j])); */
	    debug1(printf("(1) Trying deletion on segmentj of %d from %d to %d\n",nindels,splice_qpos,querylength));
	    /* ? Can re-use mismatch_positions_right because it is based on (segmentj_left + nindels) - nindels */
	    
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,&ref_nmismatches1,&ref_nmismatches2,
							   univdiagonal_j + nindels,nindels,chrhigh,
							   /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right, previously re-used*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,
							   query_compress,/*pos5*/splice_qpos,pos3,
							   querylength,indelinfo,plusp,genestrand,
							   /*want_lowest_coordinate_p*/true)) < 0) {
	      debug1(printf(" => could not find deletion on segmentj\n"));
	    } else {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos < pos3);
	      /* support_indel = indel_pos - splice_qpos; */
	      nmismatches_indel = nmismatches1;	  /* From splice to indel */
	      nmismatches_j = nmismatches2;
	      debug1(printf(" => splice_qpos %d, indel_pos %d, nmismatches_i %d, mismatches_indel %d, mismatches_j %d, prob %f",
			    splice_qpos,indel_pos,nmismatches_i,nmismatches_indel,nmismatches_j,probi+probj));
	      prob = probi + probj;
	      /* Previously looked for best_prob, but as long as the
		 probs are sufficient, we should minimize nmismatches */
	      
	      if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf(" **"));
		if (segmentj_types[j] >= 0) {
		  *best_type = KNOWN_SITE;
		} else if (segmentj_types[j] == -1) {
		  *best_type = TYPE1_SITE;
		} else if (segmentj_types[j] == -2) {
		  *best_type = TYPE2_SITE;
		}

		*best_nindels = -nindels; /* Flip sign */
		*best_indel_pos = indel_pos;
		*best_nmismatches_i = *best_ref_nmismatches_i = nmismatches_i;
		*best_nmismatches_indel = *best_ref_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_j = *best_ref_nmismatches_j = nmismatches_j;
		/* *best_trimpos3 = trimpos3; */
		best_splice_qpos = splice_qpos;
		/* best_types_i = segmenti_types[i]; */
		/* best_types_j = segmentj_types[j]; */
		best_probi = probi;
		best_probj = probj;
		best_prob = prob;
		best_nmismatches = nmismatches;
	      }
	      debug1(printf("\n"));
	    }
	  }
	}
	  
	j++;
      }

      if (j < segmentj_nsites && segmentj_sites[j] == splice_qpos) {
	/* Splice without indel */
	j++;
      }

      /* Insertions on segmentj */
      while (j < segmentj_nsites && segmentj_sites[j] - splice_qpos <= MAX_INSERTIONLEN) {
	if ((segmenti_types[i] < TYPE1_SITE || segmentj_types[j] < TYPE1_SITE) &&
	    segmenti_types[i] != segmentj_types[j]) {
	  /* Not compatible */
	} else {
	  probj = segmentj_probs[j];

	  debug1(printf("Insertion: probj at %d is %f\n",segmentj_sites[j],probj));
	  if (probj >= SPLICEINDEL_PROB || probj == optimal_segmentj_prob) {
	    nindels = segmentj_sites[j] - splice_qpos; /* Insertions should be positive */
	    debug1(printf("univdiagonal_i %u, univdiagonal_j %u, nindels %d\n",
			  univdiagonal_i,univdiagonal_j,nindels));
	    /* debug1(printf("(1) Trying insertion on segmentj of %d from %d to %d\n",nindels,splice_qpos,segmentj_sites[j])); */
	    debug1(printf("(1) Trying insertion on segmentj of %d from %d to %d\n",nindels,splice_qpos,querylength));
	    /* ? Can re-use mismatch_positions_right because it is based on (segmentj_left + nindels) - nindels */
	    
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,&ref_nmismatches1,&ref_nmismatches2,
							    univdiagonal_j + nindels,nindels,chrhigh,
							    /*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right, previously re-used*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,
							    query_compress,/*pos5*/splice_qpos,pos3,
							    querylength,indelinfo,plusp,genestrand,
							    /*want_lowest_coordinate_p*/true)) < 0) {
	      debug1(printf(" => Could not find insertion on segmentj\n"));
	    } else {
	      assert(indel_pos > splice_qpos);
	      assert(indel_pos + nindels < pos3);
	      /* support_indel = indel_pos - splice_qpos; */
	      nmismatches_indel = nmismatches1;	  /* From splice to indel */
	      nmismatches_j = nmismatches2;
	      debug1(printf(" => splice_qpos %d, indel_pos %d, mismatches_i %d, mismatches_indel %d, mismatches_j %d, prob %f",
			    splice_qpos,indel_pos,nmismatches_i,nmismatches_indel,nmismatches_j,probi+probj));
	      prob = probi + probj;
	      /* Previously looked for best_prob, but as long as the
		 probs are sufficient, we should minimize nmismatches */
	      if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {		
		debug1(printf(" **"));
		if (segmentj_types[j] >= 0) {
		  *best_type = KNOWN_SITE;
		} else if (segmentj_types[j] == -1) {
		  *best_type = TYPE1_SITE;
		} else if (segmentj_types[j] == -2) {
		  *best_type = TYPE2_SITE;
		}

		*best_nindels = -nindels; /* Flip sign */
		*best_indel_pos = indel_pos;
		*best_nmismatches_i = *best_ref_nmismatches_i = nmismatches_i;
		*best_nmismatches_indel = *best_ref_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_j = *best_ref_nmismatches_j = nmismatches_j;
		/* *best_trimpos3 = trimpos3; */
		best_splice_qpos = splice_qpos;
		/* best_types_i = segmenti_types[i]; */
		/* best_types_j = segmentj_types[j]; */
		best_probi = probi;
		best_probj = probj;
		best_prob = prob;
		best_nmismatches = nmismatches;
	      }
	      debug1(printf("\n"));
	    }
	  }
	}

	j++;
      }
    }

    if (++i < segmenti_nsites) {
      /* Count mismatches, which are also from low qpos to high qpos */
      splice_qpos = segmenti_sites[i];
      while (nmismatches_i <= nmismatches_left && mismatch_positions_left[nmismatches_i] < splice_qpos) {
	nmismatches_i++;
      }

#if 0
      if (nmismatches_i > nmismatches_allowed) {
	printf("(2) nmismatches_i for splice_qpos %d is >= %d (exceeds nmismatches_allowed)\n",
	       splice_qpos,nmismatches_i);
      } else {
	printf("(2) nmismatches_i for splice_qpos %d is %d\n",splice_qpos,nmismatches_i);
      }
#endif
    }
  }
    

  /* Right anchor (splice) */
  /* Advancing segmentj splice site from high qpos to low qpos.  Same direction for segmenti */
  /* All splicesites run from low qpos to high qpos */
  debug1(printf(">Right anchor -- splice on segmentj and indel on segmenti:\n"));
  nmismatches_j = 0;
  i = segmenti_nsites - 1;
  j = segmentj_nsites - 1;

  /* Count mismatches, which are also from high qpos to low qpos */
  splice_qpos = segmentj_sites[j];
  while (nmismatches_j <= nmismatches_right && mismatch_positions_right[nmismatches_j] >= splice_qpos) {
    nmismatches_j++;
  }

#if 0
  if (nmismatches_j > nmismatches_allowed) {
    printf("(1) nmismatches_j for splice_qpos %d is >= %d (exceeds nmismatches_allowed)\n",
	   splice_qpos,nmismatches_j);
  } else {
    printf("(1) nmismatches_j for splice_qpos %d is %d\n",splice_qpos,nmismatches_j);
  }
#endif

  /* nmismatches_j does not hold for indels */
  while (j >= 0 /*&& nmismatches_j <= nmismatches_allowed*/) {
    probj = segmentj_probs[j];

    debug1(printf("Probj at %d is %f\n",splice_qpos,probj));

    if (probj >= SPLICEINDEL_PROB || probj == optimal_segmentj_prob) {
      /* Backup to high qpos, and then advance to low qpos */
      while (i + 1 < segmenti_nsites && segmenti_sites[i+1] - splice_qpos <= max_deletionlen) {
	debug1(printf("Backing up i to %d because %d - splice_qpos %d < max_deletionlen %d\n",
		      i + 1,segmenti_sites[i+1],splice_qpos,MAX_INSERTIONLEN));
	i++;
      }

      /* Advance */
      while (i >= 0 && segmenti_nsites && segmenti_sites[i] - splice_qpos > max_deletionlen) {
	i--;
      }

      /* Deletions on segmenti */
      while (i >= 0 && segmenti_sites[i] > splice_qpos) {
	assert(segmenti_sites[i] - splice_qpos <= max_deletionlen);

	if ((segmenti_types[i] < TYPE1_SITE || segmentj_types[j] < TYPE1_SITE) &&
	    segmenti_types[i] != segmentj_types[j]) {
	  /* Not compatible */
	} else {
	  probi = segmenti_probs[i];

	  debug1(printf("Deletion: probi at %d is %f\n",segmenti_sites[i],probi));

	  if (probi >= SPLICEINDEL_PROB || probi == optimal_segmenti_prob) {
	    nindels = splice_qpos - segmenti_sites[i]; /* Deletions should be negative */
	    debug1(printf("univdiagonal_i %u, univdiagonal_j %u, nindels %d\n",
			  univdiagonal_i,univdiagonal_j,nindels));
	    /* debug1(printf("(2) Trying deletion on segmenti of %d from %d to %d\n",nindels,segmenti_sites[i],splice_qpos)); */
	    debug1(printf("(2) Trying deletion on segmenti of %d from %d to %d\n",nindels,0,splice_qpos));
	    /* ? Can re-use mismatch_positions_left because it is based on segmenti_left */
	    
	    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches1,&nmismatches2,&ref_nmismatches1,&ref_nmismatches2,
							   univdiagonal_i,nindels,chrhigh,
							   /*mismatch_positions_left, previously re-used*/NULL,/*nmismatches_left*/0,
							   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							   /*ome*/genomebits,/*ome_alt*/genomebits_alt,
							   query_compress,pos5,/*pos3*/splice_qpos,
							   querylength,indelinfo,plusp,genestrand,
							   /*want_lowest_coordinate_p*/true)) < 0) {
	      debug1(printf(" => could not find deletion on segmenti\n"));
	    } else {
	      assert(indel_pos > pos5);
	      assert(indel_pos < splice_qpos);
	      /* support_indel = splice_qpos - indel_pos; */
	      nmismatches_i = nmismatches1;
	      nmismatches_indel = nmismatches2; /* From indel to splice */
	      debug1(printf(" => indel_pos %d, splice_qpos %d, mismatches_i %d, mismatches_indel %d, nmismatches_j %d, prob %f",
			    indel_pos,splice_qpos,nmismatches_i,nmismatches_indel,nmismatches_j,probi+probj));
	      prob = probi + probj;
	      /* Previously looked for best_prob, but as long as the
		 probs are sufficient, we should minimize nmismatches */
	      if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf(" **"));
		if (segmenti_types[i] >= 0) {
		  *best_type = KNOWN_SITE;
		} else if (segmenti_types[i] == -1) {
		  *best_type = TYPE1_SITE;
		} else if (segmenti_types[i] == -2) {
		  *best_type = TYPE2_SITE;
		}

		*best_nindels = -nindels; /* Flip sign */
		*best_indel_pos = indel_pos;
		*best_nmismatches_i = *best_ref_nmismatches_i = nmismatches_i;
		*best_nmismatches_indel = *best_ref_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_j = *best_ref_nmismatches_j = nmismatches_j;
		/* *best_trimpos5 = trimpos5; */
		best_splice_qpos = splice_qpos;
		/* best_types_i = segmenti_types[i]; */
		/* best_types_j = segmentj_types[j]; */
		best_probi = probi;
		best_probj = probj;
		best_prob = prob;
		best_nmismatches = nmismatches;
	      }
	      debug1(printf("\n"));
	    }
	  }
	}

	i--;
      }

      if (i >= 0 && segmenti_sites[i] == splice_qpos) {
	/* Splice without indel */
	i--;
      }

      /* Insertions on segmenti */
      while (i >= 0 && splice_qpos - segmenti_sites[i] <= MAX_INSERTIONLEN) {
	if ((segmenti_types[i] < TYPE1_SITE || segmentj_types[j] < TYPE1_SITE) &&
	    segmenti_types[i] != segmentj_types[j]) {
	  /* Not compatible */
	} else {
	  probi = segmenti_probs[i];

	  debug1(printf("Insertion: probi at %d is %f\n",segmenti_sites[i],probi));

	  if (probi >= SPLICEINDEL_PROB || probi == optimal_segmenti_prob) {
	    nindels = splice_qpos - segmenti_sites[i]; /* Insertions should be positive */
	    debug1(printf("univdiagonal_i %u, univdiagonal_j %u, nindels %d\n",
			  univdiagonal_i,univdiagonal_j,nindels));
	    /* debug1(printf("(2) Trying insertion on segmenti of %d from %d to %d\n",nindels,segmenti_sites[i],splice_qpos)); */
	    debug1(printf("(2) Trying insertion on segmenti of %d from %d to %d\n",nindels,0,splice_qpos));
	    /* ? Can re-use mismatch_positions_left because it is based on segmenti_left */
	    
	    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches1,&nmismatches2,&ref_nmismatches1,&ref_nmismatches2,
							    univdiagonal_i,nindels,chrhigh,
							    /*mismatch_positions_left, previously re-used*/NULL,/*nmismatches_left*/0,
							    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
							    /*ome*/genomebits,/*ome_alt*/genomebits_alt,
							    query_compress,pos5,/*pos3*/splice_qpos,
							    querylength,indelinfo,plusp,genestrand,
							    /*want_lowest_coordinate_p*/true)) < 0) {
	      debug1(printf(" => could not find insertion on segmenti\n"));
	    } else {
	      assert(indel_pos > pos5);
	      assert(indel_pos + nindels < splice_qpos);
	      /* support_indel = splice_qpos - (indel_pos + nindels); */
	      nmismatches_i = nmismatches1;
	      nmismatches_indel = nmismatches2; /* From indel to splice */
	      debug1(printf(" => indel_pos %d, splice_qpos %d, mismatches_i %d, mismatches_indel %d, mismatches_j %d, prob %f",
			    indel_pos,splice_qpos,nmismatches_i,nmismatches_indel,nmismatches_j,probi+probj));
	      prob = probi + probj;
	      /* Previously looked for best_prob, but as long as the
		 probs are sufficient, we should minimize nmismatches */
	      if ((nmismatches = nmismatches_indel + nmismatches_i + nmismatches_j) < best_nmismatches) {
		debug1(printf(" **"));
		if (segmenti_types[i] >= 0) {
		  *best_type = KNOWN_SITE;
		} else if (segmenti_types[i] == -1) {
		  *best_type = TYPE1_SITE;
		} else if (segmenti_types[i] == -2) {
		  *best_type = TYPE2_SITE;
		}

		*best_nindels = -nindels; /* Flip sign */
		*best_indel_pos = indel_pos;
		*best_nmismatches_i = *best_ref_nmismatches_i = nmismatches_i;
		*best_nmismatches_indel = *best_ref_nmismatches_indel = nmismatches_indel;
		*best_nmismatches_j = *best_ref_nmismatches_j = nmismatches_j;
		/* *best_trimpos5 = trimpos5; */
		best_splice_qpos = splice_qpos;
		/* best_types_i = segmenti_types[i]; */
		/* best_types_j = segmentj_types[j]; */
		best_probi = probi;
		best_probj = probj;
		best_prob = prob;
		best_nmismatches = nmismatches;
	      }
	      debug1(printf("\n"));
	    }
	  }
	}

	i--;
      }

    }
      
    if (--j >= 0) {
      /* Count mismatches, which are also from high qpos to low qpos */
      splice_qpos = segmentj_sites[j];
      while (nmismatches_j <= nmismatches_right && mismatch_positions_right[nmismatches_j] >= splice_qpos) {
	nmismatches_j++;
      }

#if 0
      if (nmismatches_j > nmismatches_allowed) {
	printf("(2) nmismatches_j for splice_qpos %d is >= %d (exceeds nmismatches_allowed)\n",
	       splice_qpos,nmismatches_j);
      } else {
	printf("(2) nmismatches_j for splice_qpos %d is %d\n",splice_qpos,nmismatches_j);
      }
#endif
    }
  }
  
  if (plusp == sense_forward_p) {
    *best_donor_prob = best_probi;
    *best_acceptor_prob = best_probj;
  } else {
    *best_donor_prob = best_probj;
    *best_acceptor_prob = best_probi;
  }

  if (best_splice_qpos <= 0) {
    /* Do not trim */

  } else if (*best_indel_pos < best_splice_qpos) {
    if (trim5p == false) {
      *trimpos5 = pos5;
    } else if ((*best_indel_pos) - pos5 < 8) {
      *trimpos5 = pos5;
    } else {
      debug1(printf("(1) Computing trimpos5 from 0 to %d\n",*best_indel_pos));
      *trimpos5 = Genomebits_trim_qstart(&(*best_nmismatches_i),query_compress,genomebits,
					 /*univdiagonal*/univdiagonal_i,querylength,
					 pos5,/*pos3*/*best_indel_pos,plusp,genestrand);
      *best_ref_nmismatches_i = *best_nmismatches_i;
    }

    if (trim3p == false) {
      *trimpos3 = pos3;
    } else if (pos3 - (*best_indel_pos) < 8) {
      *trimpos3 = pos3;
    } else {
      debug1(printf("(1) Computing trimpos3 from %d to %d\n",best_splice_qpos,querylength));
      *trimpos3 = Genomebits_trim_qend(&(*best_nmismatches_j),query_compress,genomebits,
				       /*univdiagonal*/univdiagonal_j,querylength,
				       /*pos5*/best_splice_qpos,pos3,plusp,genestrand);
      *best_ref_nmismatches_j = *best_nmismatches_j;
    }

  } else {
    if (trim5p == false) {
      *trimpos5 = pos5;
    } else if (best_splice_qpos - pos5 < 8) {
      *trimpos5 = pos5;
    } else {
      debug1(printf("(2) Computing trimpos5 from 0 to %d\n",best_splice_qpos));
      *trimpos5 = Genomebits_trim_qstart(&(*best_nmismatches_i),query_compress,genomebits,
					 /*univdiagonal*/univdiagonal_i,querylength,
					 pos5,/*pos3*/best_splice_qpos,plusp,genestrand);
      *best_ref_nmismatches_i = *best_nmismatches_i;
    }

    ninserts = (*best_nindels < 0) ? -(*best_nindels) : 0;
    if (trim3p == false) {
      *trimpos3 = pos3;
    } else if (pos3 - ((*best_indel_pos) + ninserts) < 8) {
      *trimpos3 = pos3;
    } else {
      /* Insertions are now negative, and deletions positive, because
	 we have flipped signs */
      debug1(printf("(2) Computing trimpos3 from %d + %d to %d\n",*best_indel_pos,ninserts,querylength));
      *trimpos3 = Genomebits_trim_qend(&(*best_nmismatches_j),query_compress,genomebits,
				       /*univdiagonal*/univdiagonal_j,querylength,
				       /*pos5*/*best_indel_pos + ninserts,pos3,plusp,genestrand);
      *best_ref_nmismatches_j = *best_nmismatches_j;
    }
  }

#ifdef DEBUG1
  if (best_splice_qpos <= 0) {
    printf("=> No spliceindel\n");
    assert(*best_nindels == 0);
    assert(*best_indel_pos == -1);
  } else {
    printf("spliceindel_resolve returning splice_qpos %d, nindels %d indel_pos %d, mismatches %d + %d + %d, trimpos5 %d, trimpos3 %d\n",
	   best_splice_qpos,*best_nindels,*best_indel_pos,
	   *best_nmismatches_i,*best_nmismatches_indel,*best_nmismatches_j,
	   *trimpos5,*trimpos3);
  }
#endif

  return best_splice_qpos;
}


#if 0
/* ? Needs to be at least 3 to handle indels near splice site.  Now we
   are using IGNORE_MATCHES for spliceindel */

/* Note: types holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

/* We set check_support_p to be false for inner resolve, where we can
   expect short pieces, and have some certainty in the short
   concordant region that they are correct */
int
Splice_resolve (int *best_type, int *typei, int *typej, int *trimpos5, int *trimpos3,
		Univcoord_T *middle_univdiagonal, int *splice_qpos_i, int *splice_qpos_j,
		int *best_nindels, int *best_indel_pos,
		int *best_nmismatches_i, int *nmismatches_middle, int *best_nmismatches_j,
		int *best_nmismatches_indel,

		int *best_ref_nmismatches_i, int *ref_nmismatches_middle, int *best_ref_nmismatches_j,
		int *best_ref_nmismatches_indel,
		double *best_donor1_prob, double *best_acceptor1_prob,
		double *best_donor2_prob, double *best_acceptor2_prob,

		Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		Stage1_T stage1, Compress_T query_compress, char *queryptr,
		bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		Localdb_T localdb, int localdb_nmismatches_allowed,
		
		int pos5, int pos3, int querylength,
		Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		Intlistpool_T intlistpool,
		bool sense_forward_p, int genestrand, bool trim5p, bool trim3p, bool innerp,
		bool allow_indels_p) {

  int type;
  int best_splice_qpos, bound5, bound3, best_splice_querypos = -1;
  int splice_region_low, splice_region_high, splice_qpos_qstart, splice_qpos_qend;
  int splice_qpos_i_low, splice_qpos_i_high, splice_qpos_j_low, splice_qpos_j_high;
  int *mismatch_positions_left, *mismatch_positions_right;
  int nmismatches_left, nmismatches_right;
  char donor1, donor2, acceptor1, acceptor2;
  int ignore_nmismatches;

  int splice_qpos, indel_pos, spliceindel_trimpos5, spliceindel_trimpos3;
  int nindels, spliceindel_nmismatches_i, spliceindel_nmismatches_j, nmismatches_indel,
    spliceindel_ref_nmismatches_i, spliceindel_ref_nmismatches_j, ref_nmismatches_indel;
  double spliceindel_donor1_prob, spliceindel_acceptor2_prob;

  /* int best_nmismatches, nmismatches, segmenti_nmismatches, segmentj_nmismatches; */
  /* double best_probi, best_probj; */
  double optimal_donor1_prob, optimal_acceptor2_prob, optimal_segmenti_prob, optimal_segmentj_prob;
  /* double best_prob, donor_prob, acceptor_prob; */

  int *segmenti_sites, *segmentj_sites;
  int *segmenti_types, *segmentj_types;
  double *segmenti_probs, *segmentj_probs;
  int segmenti_nsites, segmentj_nsites;

  int supporti, supportj;
  
  int nmismatches_to_trimpos;

  Univcoord_T *middle_univdiagonals;
  int n_middle_univdiagonals;

  bool invertedp = false;

#ifdef DEBUG1
  char *gbuffer1, *gbuffer2;
  static int ncalls = 0;
#endif


  *middle_univdiagonal = (Univcoord_T) 0;
  *best_nmismatches_i = *best_nmismatches_j = querylength;
  *best_nindels = *best_nmismatches_indel = 0;
  *best_indel_pos = -1;
  *best_donor1_prob = *best_acceptor1_prob = 0.0;
  *best_donor2_prob = *best_acceptor2_prob = 0.0;

#if defined(DEBUG1) || defined(TRIM_AT_CHROMOSOME_BOUNDS)
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
#endif

#ifdef DEBUG1
  printf("Splice_resolve, plusp %d, compress_fwdp %d, sense_forward_p %d, call %d: Getting genome at univdiagonal_i %u [%u] and univdiagonal_j %u [%u] (diff: %d), range %d..%d\n",
	 plusp,Compress_fwdp(query_compress),sense_forward_p,++ncalls,
	 univdiagonal_i,univdiagonal_i - chroffset,univdiagonal_j,univdiagonal_j - chroffset,
	 univdiagonal_j-univdiagonal_i,pos5,pos3);
#endif
  assert(Compress_fwdp(query_compress) == plusp);
  assert(univdiagonal_j > univdiagonal_i);


#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonal_i + pos5 >= chroffset + querylength) ? pos5 : (int) (chroffset - segmenti_left);
  pos3 = (univdiagonal_j + pos3 <= chrhigh + querylength) ? pos3 : (int) (chrhigh - segmentj_left);
#endif


  /* Require separation from endpoints */
#if 0
  splice_qpos_start = pos5 + 1;
  splice_qpos_end = pos3 - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > (querylength - 1) - min_shortend) {
    splice_qpos_end = (querylength - 1) - min_shortend;
  }
#endif

  /* Determine feasibility of a splice */

  /* New method, different from relying on nmismatches_allowed */
  mismatch_positions_left = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_right = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  /* ascending */
  nmismatches_left =
    Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					    genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,pos5,pos3,plusp,genestrand);

  /* descending */
  nmismatches_right =
    Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					     genomebits,genomebits_alt,query_compress,
					     univdiagonal_j,querylength,pos5,pos3,plusp,genestrand);

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  /* was splice_qpos_high */
  splice_qpos_qend = Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_left,nmismatches_left,
						   pos5,pos3,querylength);
  /* was splice_qpos_low */
  splice_qpos_qstart = Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_right,nmismatches_right,
						    pos5,pos3);
  debug1(printf("trim_nosplice for pos5 %d, pos3 %d yields low %d and high %d\n",
		pos5,pos3,splice_qpos_qstart,splice_qpos_qend));

  /* Trimming */
  if (trim5p == false) {
    *trimpos5 = pos5;
  } else if (splice_qpos_qend - pos5 < 8) {
    *trimpos5 = pos5;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qstart with %d..%d",pos5,splice_qpos_qend));
    *trimpos5 = Genomebits_trim_qstart(&ignore_nmismatches,query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/splice_qpos_qend,plusp,genestrand);
    debug1(printf(" to yield %d..\n",*trimpos5));
  }

  if (trim3p == false) {
    *trimpos3 = pos3;
  } else if (pos3 - splice_qpos_qstart < 8) {
    *trimpos3 = pos3;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qend with %d..%d",splice_qpos_qstart,pos3));
    *trimpos3 = Genomebits_trim_qend(&ignore_nmismatches,query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/splice_qpos_qstart,pos3,plusp,genestrand);
    debug1(printf(" to yield ..%d\n",*trimpos3));
  }

  if (*trimpos5 != pos5) {
    /* ascending */
    nmismatches_left =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					      genomebits,genomebits_alt,query_compress,
					      univdiagonal_i,querylength,*trimpos5,pos3,plusp,genestrand);
  }

  if (*trimpos3 != pos3) {
    /* descending */
    nmismatches_right =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					       genomebits,genomebits_alt,query_compress,
					       univdiagonal_j,querylength,pos5,*trimpos3,plusp,genestrand);
  }

  pos5 = *trimpos5;
  pos3 = *trimpos3;

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );


  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_qstart));
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_qend));


#ifdef DEBUG1
  gbuffer1 = (char *) CALLOC(querylength+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonal_i - querylength,querylength,gbuffer1);
  Genome_fill_buffer(univdiagonal_j - querylength,querylength,gbuffer2);
  char *queryseq = Compress_queryseq(query_compress,querylength);

  printf("Splice_resolve, plusp %d: Getting genome at left %llu and %llu\n",
	 plusp,(unsigned long long) segmenti_left,(unsigned long long) segmentj_left);
  printf("g1: %s\n",gbuffer1);
  printf("    ");
  print_diffs(gbuffer1,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbuffer2,queryseq,querylength);
  printf("\n");
  printf("g2: %s\n",gbuffer2);
  FREE(queryseq);
  FREE(gbuffer2);
  FREE(gbuffer1);
#endif



  if (splice_qpos_qstart <= splice_qpos_qend + MISMATCHES_AT_SITE) {
    /* Appears to have matches to both splice sites.  Span just the crossing parts */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

  } else {
    /* There is a gap in matches.  Span the entire region */
    splice_qpos_i_low = splice_qpos_j_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);
  }

  debug1(printf("  Splice i bounds: %d..%d\n",splice_qpos_i_low,splice_qpos_i_high));
  debug1(printf("  Splice j bounds: %d..%d\n",splice_qpos_j_low,splice_qpos_j_high));

  splice_region_low = /*min*/ (splice_qpos_i_low < splice_qpos_j_low) ? splice_qpos_i_low : splice_qpos_j_low;
  splice_region_high = /*max*/ (splice_qpos_i_high > splice_qpos_j_high) ? splice_qpos_i_high : splice_qpos_j_high;

  if (plusp == sense_forward_p) {
    /* geneplus */
    segmenti_nsites = compute_donor_sites(&optimal_donor1_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
					  spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					  spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					  spliceinfo->segmenti_probs_alloc,
					  splice_region_low,splice_region_high,
					  querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_acceptor_sites(&optimal_acceptor2_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					     spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					     spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					     spliceinfo->segmentj_probs_alloc,
					     splice_region_low,splice_region_high,
					     querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_donor1_prob;
    optimal_segmentj_prob = optimal_acceptor2_prob;

  } else {
    /* geneminus */
    segmenti_nsites = compute_antiacceptor_sites(&optimal_acceptor2_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
						 spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						 spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						 spliceinfo->segmenti_probs_alloc,
						 splice_region_low,splice_region_high,
						 querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_antidonor_sites(&optimal_donor1_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					      spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					      spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					      spliceinfo->segmentj_probs_alloc,
					      splice_region_low,splice_region_high,
					      querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_acceptor2_prob;
    optimal_segmentj_prob = optimal_donor1_prob;
  }


  best_splice_qpos = -1;

  /* Try standard splicing */
  debug1(printf("Computing splice sites\n"));
    
  /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
  if (plusp == true) {
    /* Skip */
    check_ascending(mismatch_positions_left,nmismatches_left);
    check_descending(mismatch_positions_right,nmismatches_right);

#if 0
  } else if (invertedp == true) {
    /* Already inverted */
    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);
#endif

  } else {
    invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
    invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
    invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);

    invertedp = true;
  }
    
  debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

  if (sense_forward_p == true) {
    if (plusp == true) {
      /* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
		     /*D*/&(*best_nmismatches_i),/*A*/&(*best_nmismatches_j),
		     /*D*/&(*best_ref_nmismatches_i),/*A*/&(*best_ref_nmismatches_j),
		     /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/true,/*plusAp*/true,querylength,
		     /*donor*/mismatch_positions_left,/*donor*/nmismatches_left,
		     /*acceptor*/mismatch_positions_right,/*acceptor*/nmismatches_right,
		     /*D*/segmenti_sites,/*A*/segmentj_sites,
		     /*D*/segmenti_types,/*A*/segmentj_types,
		     /*D*/segmenti_probs,/*A*/segmentj_probs,
		     /*D*/segmenti_nsites,/*A*/segmentj_nsites);
	
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	*best_type = type;
      }
      
    } else {
      /* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
		     /*D*/&(*best_nmismatches_j),/*A*/&(*best_nmismatches_i),
		     /*D*/&(*best_ref_nmismatches_j),/*A*/&(*best_ref_nmismatches_i),
		     /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/false,/*plusAp*/false,querylength,
		     /*donor*/mismatch_positions_right,/*donor*/nmismatches_right,
		     /*acceptor*/mismatch_positions_left,/*acceptor*/nmismatches_left,
		     /*D*/segmentj_sites,/*A*/segmenti_sites,
		     /*D*/segmentj_types,/*A*/segmenti_types,
		     /*D*/segmentj_probs,/*A*/segmenti_probs,
		     /*D*/segmentj_nsites,/*A*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	*best_type = type;
      }
    }
    
  } else {
    if (plusp == true) {
      /* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			 /*A*/&(*best_nmismatches_i),/*D*/&(*best_nmismatches_j),
			 /*A*/&(*best_ref_nmismatches_i),/*D*/&(*best_ref_nmismatches_j),
			 /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/plusp,/*plusDp*/plusp,querylength,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*A*/segmenti_sites,/*D*/segmentj_sites,
			 /*A*/segmenti_types,/*D*/segmentj_types,
			 /*A*/segmenti_probs,/*D*/segmentj_probs,
			 /*A*/segmenti_nsites,/*D*/segmentj_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	*best_type = type;
      }
      
    } else {
      /* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			 /*A*/&(*best_nmismatches_j),/*D*/&(*best_nmismatches_i),
			 /*A*/&(*best_ref_nmismatches_j),/*D*/&(*best_ref_nmismatches_i),
			 /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/false,/*plusDp*/false,querylength,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*A*/segmentj_sites,/*D*/segmenti_sites,
			 /*A*/segmentj_types,/*D*/segmenti_types,
			 /*A*/segmentj_probs,/*D*/segmenti_probs,
			 /*A*/segmentj_nsites,/*D*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	*best_type = type;
      }
    }
  }


  if (best_splice_qpos <= 0 &&
      splice_qpos_qstart >= splice_qpos_qend + MIN_EXONLEN) {
    /* Possible middle exon or multiple sequencing errors */
    debug1(printf("splice_qpos_qstart %d >= splice_qpos_qend %d + %d, so trying to find a middle exon\n",
		  splice_qpos_qstart,splice_qpos_qend,MIN_EXONLEN));

    if (plusp == true) {
      n_middle_univdiagonals =
	Spliceends_middle_plus(&middle_univdiagonals,
			       stage1,streamspace_max_alloc,streamspace_alloc,
			       /*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
			       /*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
			       query_compress,queryptr,
			       univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
			       localdb_nmismatches_allowed);
      /* mismatch_positions_left is ascending */
      /* mismatch_positions_right is descending */

    } else {
      n_middle_univdiagonals =
	Spliceends_middle_minus(&middle_univdiagonals,
				stage1,streamspace_max_alloc,streamspace_alloc,
				/*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
				/*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
				query_compress,queryptr,
				univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
				localdb_nmismatches_allowed);
    }


    debug1(printf("The two ends cannot meet, and with knownsplicing trying to find a middle exon\n"));
    if (n_middle_univdiagonals == 0) {
      /* Skip: middle_univdiagonals is unassigned */

    } else {
      if (plusp == true) {
	/* Skip */
	check_ascending(mismatch_positions_left,nmismatches_left);
	check_descending(mismatch_positions_right,nmismatches_right);

      } else if (invertedp == true) {
	/* Already inverted */
	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

      } else {
	invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
	invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
	invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
	invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

	invertedp = true;
      }
      
      if ((*middle_univdiagonal =
	   find_middle_exon(&(*typei),&(*typej),&(*splice_qpos_i),&(*splice_qpos_j),
			    &(*best_nmismatches_i),&(*nmismatches_middle),&(*best_nmismatches_j),
			    &(*best_ref_nmismatches_i),&(*ref_nmismatches_middle),&(*best_ref_nmismatches_j),
			    &(*best_donor1_prob),&(*best_acceptor1_prob),
			    &(*best_donor2_prob),&(*best_acceptor2_prob),
			    middle_univdiagonals,n_middle_univdiagonals,univdiagonal_i,univdiagonal_j,
			    pos5,pos3,/*splice_qpos_5*/splice_qpos_qend,/*splice_qpos_3*/splice_qpos_qstart,querylength,
			    
			    mismatch_positions_left,nmismatches_left,
			    mismatch_positions_right,nmismatches_right,
			    
			    segmenti_sites,segmentj_sites,
			    segmenti_types,segmentj_types,
			    segmenti_probs,segmentj_probs,
			    segmenti_nsites,segmentj_nsites,
			    
			    query_compress,chroffset,spliceinfo,knownsplicing,intlistpool,
			    plusp,sense_forward_p,genestrand)) == 0) {
	FREE(middle_univdiagonals);

      } else {
	/* find_middle_exon returns values as qpos, not querypos */
	debug1(printf("Found a middle exon\n"));
	FREE(middle_univdiagonals);
	
	/* Do not trim */
	/* *trimpos5 = pos5; */
	/* *trimpos3 = pos3; */
	return -1;		/* This plus *middle_univdiagonal indicates a success */
      }
    }
  }

  debug1(printf("(1) Before spliceindel_resolve, have splice_qpos %d, nsites %d and %d, probs %f vs %f and %f vs %f\n",
		best_splice_qpos,segmenti_nsites,segmentj_nsites,
		*best_donor1_prob,optimal_donor1_prob,*best_acceptor2_prob,optimal_acceptor2_prob));
  if (allow_indels_p == false) {
    /* Skip */
    debug1(printf("(1) Indels not allowed\n"));

  } else if (segmenti_nsites == 0 || segmentj_nsites == 0) {
    debug1(printf("Cannot call spliceindel_resolve because nsites %d and %d\n",
		  segmenti_nsites,segmentj_nsites));

  } else if (*best_donor1_prob == optimal_donor1_prob && *best_acceptor2_prob == optimal_acceptor2_prob) {
    debug1(printf("No need to call spliceindel_resolve because exact probs are optimal\n"));
    
  } else {
    /* Try splice plus indel */

    if (invertedp == true) {
      /* Computes in qpos, not querypos, so undo inversions */
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = false;
    }

    debug1(printf("For comparison, exact nmismatches are %d + %d\n",
		  *best_nmismatches_i,*best_nmismatches_j));
    if ((splice_qpos =
	 spliceindel_resolve(&type,&spliceindel_trimpos5,&spliceindel_trimpos3,&nindels,&indel_pos,
			     &spliceindel_nmismatches_i,&spliceindel_nmismatches_j,&nmismatches_indel,
			     &spliceindel_ref_nmismatches_i,&spliceindel_ref_nmismatches_j,
			     &ref_nmismatches_indel,
			     &spliceindel_donor1_prob,&spliceindel_acceptor2_prob,
			     
			     univdiagonal_i,univdiagonal_j,
			     query_compress,plusp,chroffset,chrhigh,
			     
			     mismatch_positions_left,nmismatches_left,
			     mismatch_positions_right,nmismatches_right,
			     
			     segmenti_sites,segmentj_sites,segmenti_types,segmentj_types,
			     segmenti_probs,segmentj_probs,segmenti_nsites,segmentj_nsites,
			     optimal_segmenti_prob,optimal_segmentj_prob,
			     
			     pos5,pos3,querylength,indelinfo,sense_forward_p,genestrand,
			     trim5p,trim3p)) <= 0) {
      /* No spliceindel found */

    } else if (best_splice_qpos > 0 &&
	       spliceindel_nmismatches_i + nmismatches_indel + spliceindel_nmismatches_j >=
	       *best_nmismatches_i + *best_nmismatches_j) {
      /* Spliceindel not as good as exact splice */

    } else if (abs(splice_qpos - indel_pos) - 4*nmismatches_indel < 0) {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    *best_nmismatches_i,*best_nmismatches_j));
      debug1(printf("Indel segment has too many mismatches\n"));
      /* Fall through to atypical splice */

    } else {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    *best_nmismatches_i,*best_nmismatches_j));

      *best_type = type;
      *trimpos5 = spliceindel_trimpos5;
      *trimpos3 = spliceindel_trimpos3;
      *best_nindels = nindels;
      *best_indel_pos = indel_pos;
      *best_nmismatches_i = spliceindel_nmismatches_i;
      *best_nmismatches_indel = nmismatches_indel;
      *best_nmismatches_j = spliceindel_nmismatches_j;
      *best_ref_nmismatches_i = spliceindel_ref_nmismatches_i;
      *best_ref_nmismatches_indel = ref_nmismatches_indel;
      *best_ref_nmismatches_j = spliceindel_ref_nmismatches_j;
      *best_donor1_prob = spliceindel_donor1_prob;
      *best_acceptor2_prob = spliceindel_acceptor2_prob;
      best_splice_qpos = splice_qpos;
    }
  }

#if 0
  if (best_splice_qpos <= 0) {
    /* Try atypical splicing */
    debug1(printf("Computing splice sites\n"));
    
    /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
    if (plusp == true) {
      /* Skip */
    } else if (invertedp == true) {
      /* Skip */
    } else {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = true;
    }
    
    debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

    if (splice_qpos_i_low < splice_qpos_j_low) {
      splice_qpos_qstart = splice_qpos_i_low;
    } else {
      splice_qpos_qstart = splice_qpos_j_low;
    }

    if (splice_qpos_i_high > splice_qpos_j_high) {
      splice_qpos_qend = splice_qpos_i_high;
    } else {
      splice_qpos_qend = splice_qpos_j_high;
    }

    if (sense_forward_p == true) {
      if (plusp == true) {
	/* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			 /*D*/&(*best_nmismatches_i),/*A*/&(*best_nmismatches_j),
			 /*D*/&(*best_ref_nmismatches_i),/*A*/&(*best_ref_nmismatches_j),
			 /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/true,/*plusAp*/true,querylength,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);
	  
	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			 /*D*/&(*best_nmismatches_j),/*A*/&(*best_nmismatches_i),
			 /*D*/&(*best_ref_nmismatches_j),/*A*/&(*best_ref_nmismatches_i),
			 /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/false,/*plusAp*/false,querylength,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*splice_querypos_low*/querylength - splice_qpos_qend,
			 /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }

    } else {
      if (plusp == true) {
	/* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			     /*A*/&(*best_nmismatches_i),/*D*/&(*best_nmismatches_j),
			     /*A*/&(*best_ref_nmismatches_i),/*D*/&(*best_ref_nmismatches_j),
			     /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/true,/*plusDp*/true,querylength,
			     /*acceptor*/mismatch_positions_left,nmismatches_left,
			     /*donor*/mismatch_positions_right,nmismatches_right,
			     /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor1_prob),&(*best_acceptor2_prob),
			     /*A*/&(*best_nmismatches_j),/*D*/&(*best_nmismatches_i),
			     /*A*/&(*best_ref_nmismatches_j),/*D*/&(*best_ref_nmismatches_i),
			     /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/false,/*plusDp*/false,querylength,
			     /*acceptor*/mismatch_positions_right,nmismatches_right,
			     /*donor*/mismatch_positions_left,nmismatches_left,
			     /*splice_querypos_low*/querylength - splice_qpos_qend,
			     /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }
    }
  }    
#endif


#if 0
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < best_splice_qpos) {
    lefti++;
  }
  *best_nmismatches_i = *best_ref_nmismatches_i = lefti;
  
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= best_splice_qpos) {
    righti++;
  }
  *best_nmismatches_j = *best_ref_nmismatches_j = righti;
#endif


  if (best_splice_qpos <= 0) {
    return -1;			/* Fail */

  } else {
#if 0
    if (plusp == sense_forward_p) {
      best_probi = *best_donor1_prob;
      best_probj = *best_acceptor2_prob;
    } else {
      best_probi = *best_acceptor2_prob;
      best_probj = *best_donor1_prob;
    }
#endif
    debug1(printf("donor1 %f, acceptor2 %f\n",*best_donor1_prob,*best_acceptor2_prob));

    /* Check that the exon-intron boundaries are okay */
    if (invertedp == true) {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    }

#ifdef DEBUG1
    printf("Checking bounds\n");
    printf("%d mismatches on left:",nmismatches_left);
    for (int i = 0; i <= nmismatches_left; i++) {
      printf(" %d",mismatch_positions_left[i]);
    }
    printf("\n");

    printf("%d mismatches on right:",nmismatches_right);
    for (int i = 0; i <= nmismatches_right; i++) {
      printf(" %d",mismatch_positions_right[i]);
    }
    printf("\n");
#endif

#if 0
    bound5 = Spliceends_trim_qend_nosplice(&ignore_nmismatches,mismatch_positions_left,nmismatches_left,
					   *trimpos5,best_splice_qpos,querylength);
    bound3 = Spliceends_trim_qstart_nosplice(&ignore_nmismatches,mismatch_positions_right,nmismatches_right,
					     best_splice_qpos,*trimpos3);

    debug1(printf("Comparing bound5 %d with splice_qpos %d\n",bound5,best_splice_qpos));
    debug1(printf("Comparing bound3 %d with splice_qpos %d\n",bound3,best_splice_qpos));
#endif

    supporti = best_splice_qpos - (*trimpos5);
    supportj = (*trimpos3) - best_splice_qpos;
#if 0
    if (plusp == sense_forward_p) {
      *adj_donor_support = supporti - 4*(*best_nmismatches_i);
      *adj_acceptor_support = supportj - 4*(*best_nmismatches_j);
    } else {
      *adj_donor_support = supportj - 4*(*best_nmismatches_j);
      *adj_acceptor_support = supporti - 4*(*best_nmismatches_i);
    }
#endif

    /* For local splicing, we check support before probability */
#if 0
    if (bound5 < best_splice_qpos - 2 || bound3 > best_splice_qpos + 2) {
      debug1(printf("Rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, mismatches_indel %d, based on bounds\n\n",
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj,*best_nmismatches_indel));
      return -1;
    }
#endif
      
#if 0
    if (supporti >= 25 && supportj >= 25 && (*best_nmismatches_i) + (*best_nmismatches_j) <= 1 &&
	adj_supporti > 0 && adj_supportj > 0) {
      /* Can give rise to very bad splices */
      debug1(printf("For %u and %u, accepting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, mismatches_indel %d, based on support\n\n",
		    univdiagonal_i,univdiagonal_j,
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj,*best_nmismatches_indel));
      return best_splice_qpos;
    }
#endif

#if 0
    if (/*check_support_p == true -- always true &&*/ sufficient_support_p(adj_supporti,best_probi) == false) {
      /* Since we are allowing atypical splices, should always check support */
      debug1(printf("Rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on adj support and prob\n\n",
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj));
      return -1;		/* Fail */

    } else if (/*check_support_p == true -- always true &&*/ sufficient_support_p(adj_supportj,best_probj) == false) {
      /* Since we are allowing atypical splices, should always check support */
      debug1(printf("Rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on adj support and prob\n\n",
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj));
      return -1;		/* Fail */
    }
#endif

#if 0
    if (/*salvagep == false && */ best_probi < 0.5 || best_probj < 0.5) {
      /* Previously used conjunction, but that can lead to the wrong sense reported */
      /* Previously used SPLICE_PROB_LOW, or 0.2 */
      debug1(printf("Rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on probs\n\n",
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj));
      return -1;		/* Fail */

    } else {
      debug1(printf("For %u and %u, accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    univdiagonal_i,univdiagonal_j,best_splice_qpos,*best_nindels,
		    (*best_nindels < 0) ? "insertion" : "deletion",
		    *best_indel_pos,supporti,*best_nmismatches_i,best_probi,
		    supportj,*best_nmismatches_j,best_probj));
      return best_splice_qpos;
    }
#endif

    if (Splice_accept_p(best_splice_qpos,querylength,
			*best_donor1_prob,*best_acceptor2_prob,
			univdiagonal_i,univdiagonal_j,
			supporti,*best_nmismatches_i,supportj,*best_nmismatches_j,
			plusp,sense_forward_p,innerp) == false) {
      debug1(printf("Splice_resolve_qstart on dist %u rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on probs\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,*best_nindels,*best_indel_pos,supporti,*best_nmismatches_i,*best_donor1_prob,
		    supportj,*best_nmismatches_j,*best_acceptor2_prob));
      return -1;		/* Fail */

    } else {
      debug1(printf("Splice_resolve_qstart on dist %u accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,*best_nindels,(*best_nindels < 0) ? "insertion" : (*best_nindels > 0) ? "deletion" : "",
		    *best_indel_pos,supporti,*best_nmismatches_i,*best_donor1_prob,
		    supportj,*best_nmismatches_j,*best_acceptor2_prob));
      return best_splice_qpos;
    }
  }
}
#endif



/* ? Needs to be at least 3 to handle indels near splice site.  Now we
   are using IGNORE_MATCHES for spliceindel */

/* Note: types holds joffset + j + 1, so 0 represents no known site
   and values greater than 0 represent a known site.  Need to subtract
   1 to obtain joffset + j. */

/* We set check_support_p to be false for inner resolve, where we can
   expect short pieces, and have some certainty in the short
   concordant region that they are correct */
int
Splice_resolve_qstart (Univcoordlist_T *univdiagonals, Intlist_T *nmismatches, Intlist_T *ref_nmismatches,
		       Intlist_T *endpoints, List_T *junctions,

		       Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		       Stage1_T stage1, Compress_T query_compress, char *queryptr,
		       bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		       int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       
		       Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		       Localdb_T localdb, int localdb_nmismatches_allowed,

		       int pos5, int pos3, int querylength, int cont_trimpos, int cont_nmismatches,

		       Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		       Univcoordlistpool_T univcoordlistpool, Intlistpool_T intlistpool,
		       Listpool_T listpool, Pathpool_T pathpool,

		       int sensedir, int genestrand, bool trim5p, bool trim3p, bool innerp,
		       bool allow_indels_p, bool find_middle_p) {

  int type;
  Univcoord_T middle_univdiagonal = (Univcoord_T) 0;
  int best_nmismatches_i = querylength, best_nmismatches_j = querylength;
  int best_ref_nmismatches_i = querylength, best_ref_nmismatches_j = querylength;

  int best_nindels = 0, best_nmismatches_indel = 0, best_ref_nmismatches_indel = 0, ninserts;
  int best_indel_pos = -1;
  double best_donor1_prob = 0.0, best_acceptor2_prob = 0.0;
  int trimpos5, trimpos3;

  Chrpos_T splice_distance;

#if 1
  /* For middle exon */
  Chrpos_T splice_distance_i, splice_distance_j;
  double best_donor2_prob = 0.0, best_acceptor1_prob = 0.0;
  int splice_qpos_i, splice_qpos_j;
  int nmismatches_middle = querylength, ref_nmismatches_middle = querylength;
  int typei, typej;
  int n_middle_univdiagonals;
  Univcoord_T *middle_univdiagonals;
#endif  

  int best_splice_qpos, best_splice_querypos = -1;
  int splice_qpos_qstart, splice_qpos_qend;
  int splice_qpos_i_low, splice_qpos_i_high, splice_qpos_j_low, splice_qpos_j_high;
  int *mismatch_positions_left, *mismatch_positions_right;
  int nmismatches_left, nmismatches_right;
  char donor1, donor2, acceptor1, acceptor2;
  int ignore_nmismatches;

  int splice_qpos, indel_pos, spliceindel_trimpos5, spliceindel_trimpos3;
  int nindels, spliceindel_nmismatches_i, spliceindel_nmismatches_j, nmismatches_indel,
    spliceindel_ref_nmismatches_i, spliceindel_ref_nmismatches_j, ref_nmismatches_indel;
  double spliceindel_donor1_prob, spliceindel_acceptor2_prob;

  /* int best_nmismatches, nmismatches, segmenti_nmismatches, segmentj_nmismatches; */
  /* double best_probi, best_probj; */
  double optimal_donor1_prob, optimal_acceptor2_prob, optimal_segmenti_prob, optimal_segmentj_prob;
  /* double best_prob, donor_prob, acceptor_prob; */

  int *segmenti_sites, *segmentj_sites;
  int *segmenti_types, *segmentj_types;
  double *segmenti_probs, *segmentj_probs;
  int segmenti_nsites, segmentj_nsites;

  int supporti, supportj;
  
  int nmismatches_to_trimpos;
  int nmatches_trim_i, nmatches_trim_j,
    nmismatches_trim_i, nmismatches_trim_j, ignore;
  int endpoint_i, endpoint_j;

  int nmatches_stay, nmismatches_stay;
  int score_stay, score_splice;

  bool invertedp = false;
  bool sense_forward_p = (sensedir == SENSE_FORWARD);

#ifdef DEBUG1
  char *gbuffer1, *gbuffer2;
  static int ncalls = 0;
#endif

#if defined(DEBUG1) || defined(TRIM_AT_CHROMOSOME_BOUNDS)
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
#endif

#ifdef DEBUG1
  printf("Splice_resolve_qstart (%s), plusp %d, compress_fwdp %d, sense_forward_p %d, call %d: Getting genome at univdiagonal_i %u [%u] and univdiagonal_j %u [%u] (diff: %d), range %d..%d\n",
	 sensedir == SENSE_FORWARD ? "sense" : "antisense",
	 plusp,Compress_fwdp(query_compress),sense_forward_p,++ncalls,
	 univdiagonal_i,univdiagonal_i - chroffset,univdiagonal_j,univdiagonal_j - chroffset,
	 univdiagonal_j-univdiagonal_i,pos5,pos3);
#endif
  assert(Compress_fwdp(query_compress) == plusp);
  assert(univdiagonal_j > univdiagonal_i);


#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonal_i + pos5 >= chroffset + querylength) ? pos5 : (int) (chroffset - segmenti_left);
  pos3 = (univdiagonal_j + pos3 <= chrhigh + querylength) ? pos3 : (int) (chrhigh - segmentj_left);
#endif


  /* Require separation from endpoints */
#if 0
  splice_qpos_start = pos5 + 1;
  splice_qpos_end = pos3 - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > (querylength - 1) - min_shortend) {
    splice_qpos_end = (querylength - 1) - min_shortend;
  }
#endif

  /* Determine feasibility of a splice */

  /* New method, different from relying on nmismatches_allowed */
  mismatch_positions_left = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_right = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  /* ascending */
  nmismatches_left =
    Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					    genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,pos5,pos3,plusp,genestrand);

  /* descending */
  nmismatches_right =
    Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					     genomebits,genomebits_alt,query_compress,
					     univdiagonal_j,querylength,pos5,pos3,plusp,genestrand);

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );


  splice_qpos_qend = Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_left,nmismatches_left,
						   pos5,pos3,querylength);
  splice_qpos_qstart = Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_right,nmismatches_right,
						    pos5,pos3);
  debug1(printf("trim_nosplice for pos5 %d, pos3 %d yields low %d and high %d\n",
		pos5,pos3,splice_qpos_qstart,splice_qpos_qend));

  /* Trimming */
  if (trim5p == false) {
    trimpos5 = pos5;
  } else if (splice_qpos_qend - pos5 < 8) {
    trimpos5 = pos5;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qstart with %d..%d",pos5,splice_qpos_qend));
    trimpos5 = Genomebits_trim_qstart(&ignore_nmismatches,query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/splice_qpos_qend,plusp,genestrand);
    debug1(printf(" to yield %d..\n",trimpos5));
  }

  if (trim3p == false) {
    trimpos3 = pos3;
  } else if (pos3 - splice_qpos_qstart < 8) {
    trimpos3 = pos3;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qend with %d..%d",splice_qpos_qstart,pos3));
    trimpos3 = Genomebits_trim_qend(&ignore_nmismatches,query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/splice_qpos_qstart,pos3,plusp,genestrand);
    debug1(printf(" to yield ..%d\n",trimpos3));
  }

  if (trimpos5 != pos5) {
    /* ascending */
    nmismatches_left =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					      genomebits,genomebits_alt,query_compress,
					      univdiagonal_i,querylength,trimpos5,pos3,plusp,genestrand);
  }

  if (trimpos3 != pos3) {
    /* descending */
    nmismatches_right =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					       genomebits,genomebits_alt,query_compress,
					       univdiagonal_j,querylength,pos5,trimpos3,plusp,genestrand);
  }

  pos5 = trimpos5;
  pos3 = trimpos3;

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

#if 0
  /* From Splice_nomiddle */
  if (splice_qpos_qstart < pos5 + 1) {
    splice_qpos_qstart = pos5 + 1;
  }
  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_qstart));

  if (splice_qpos_qend > pos3 - 1) {
    splice_qpos_qend = pos3 - 1;
  }
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_qend));
#endif

  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_qstart));
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_qend));


#ifdef DEBUG1
  gbuffer1 = (char *) CALLOC(querylength+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonal_i - querylength,querylength,gbuffer1);
  Genome_fill_buffer(univdiagonal_j - querylength,querylength,gbuffer2);
  char *queryseq = Compress_queryseq(query_compress,querylength);

  printf("Splice_resolve, plusp %d: Getting genome at left %llu and %llu\n",
	 plusp,(unsigned long long) segmenti_left,(unsigned long long) segmentj_left);
  printf("g1: %s\n",gbuffer1);
  printf("    ");
  print_diffs(gbuffer1,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbuffer2,queryseq,querylength);
  printf("\n");
  printf("g2: %s\n",gbuffer2);
  FREE(queryseq);
  FREE(gbuffer2);
  FREE(gbuffer1);
#endif


  debug1(printf("find_middle_p: %d\n",find_middle_p));

#if 1
  if (0 && find_middle_p == false) {
    /* From Splice_nomiddle */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

  } else if (splice_qpos_qstart <= splice_qpos_qend + MISMATCHES_AT_SITE) {
    /* Appears to have matches to both splice sites.  Span just the crossing parts */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

    /* Set both to min */
    if (splice_qpos_i_low < splice_qpos_j_low) {
      splice_qpos_j_low = splice_qpos_i_low;
    } else {
      splice_qpos_i_low = splice_qpos_j_low;
    }

    /* Set both to max */
    if (splice_qpos_i_high > splice_qpos_j_high) {
      splice_qpos_j_high = splice_qpos_i_high;
    } else {
      splice_qpos_i_high = splice_qpos_j_high;
    }

  } else {
    /* There is a gap in matches.  Span the entire region */
    splice_qpos_i_low = splice_qpos_j_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);
  }

#else
  /* This code is not in Splice_nomiddle */
  if (splice_qpos_qstart <= splice_qpos_qend + MISMATCHES_AT_SITE) {
    /* Appears to have matches to both splice sites.  Span just the crossing parts */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

  } else {
    /* There is a gap in matches.  Span the entire region */
    splice_qpos_i_low = splice_qpos_j_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);
  }

  splice_region_low = /*min*/ (splice_qpos_i_low < splice_qpos_j_low) ? splice_qpos_i_low : splice_qpos_j_low;
  splice_region_high = /*max*/ (splice_qpos_i_high > splice_qpos_j_high) ? splice_qpos_i_high : splice_qpos_j_high;

  /* Good */
  debug1(printf("  Splice region: %d..%d\n",splice_region_low,splice_region_high));
#endif


  debug1(printf("  Splice i bounds: %d..%d\n",splice_qpos_i_low,splice_qpos_i_high));
  debug1(printf("  Splice j bounds: %d..%d\n",splice_qpos_j_low,splice_qpos_j_high));


  if (plusp == sense_forward_p) {
    /* geneplus */
    segmenti_nsites = compute_donor_sites(&optimal_donor1_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
					  spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					  spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					  spliceinfo->segmenti_probs_alloc,
#if 0
					  splice_region_low,splice_region_high,
#else
					  splice_qpos_i_low,splice_qpos_i_high,
#endif
					  querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_acceptor_sites(&optimal_acceptor2_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					     spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					     spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					     spliceinfo->segmentj_probs_alloc,
#if 0
					     splice_region_low,splice_region_high,
#else
					     splice_qpos_j_low,splice_qpos_j_high,
#endif
					     querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_donor1_prob;
    optimal_segmentj_prob = optimal_acceptor2_prob;

  } else {
    /* geneminus */
    segmenti_nsites = compute_antiacceptor_sites(&optimal_acceptor2_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
						 spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						 spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						 spliceinfo->segmenti_probs_alloc,
#if 0
						 splice_region_low,splice_region_high,
#else
						 splice_qpos_i_low,splice_qpos_i_high,
#endif
						 querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_antidonor_sites(&optimal_donor1_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					      spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					      spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					      spliceinfo->segmentj_probs_alloc,
#if 0
					      splice_region_low,splice_region_high,
#else
					      splice_qpos_j_low,splice_qpos_j_high,
#endif
					      querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_acceptor2_prob;
    optimal_segmentj_prob = optimal_donor1_prob;
  }


  best_splice_qpos = -1;

  /* Try standard splicing */
  debug1(printf("Computing splice sites\n"));
    
  /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
  if (plusp == true) {
    /* Skip */
    check_ascending(mismatch_positions_left,nmismatches_left);
    check_descending(mismatch_positions_right,nmismatches_right);

#if 0
  } else if (invertedp == true) {
    /* Already inverted */
    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);
#endif

  } else {
    invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
    invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
    invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);

    invertedp = true;
  }
    
  debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

  if (sense_forward_p == true) {
    if (plusp == true) {
      /* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
		     /*D*/&best_nmismatches_i,/*A*/&best_nmismatches_j,
		     /*D*/&best_ref_nmismatches_i,/*A*/&best_ref_nmismatches_j,
		     /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/true,/*plusAp*/true,querylength,
		     /*donor*/mismatch_positions_left,/*donor*/nmismatches_left,
		     /*acceptor*/mismatch_positions_right,/*acceptor*/nmismatches_right,
		     /*D*/segmenti_sites,/*A*/segmentj_sites,
		     /*D*/segmenti_types,/*A*/segmentj_types,
		     /*D*/segmenti_probs,/*A*/segmentj_probs,
		     /*D*/segmenti_nsites,/*A*/segmentj_nsites);
	
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	/* best_type = type; */
      }
      
    } else {
      /* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
		     /*D*/&best_nmismatches_j,/*A*/&best_nmismatches_i,
		     /*D*/&best_ref_nmismatches_j,/*A*/&best_ref_nmismatches_i,
		     /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/false,/*plusAp*/false,querylength,
		     /*donor*/mismatch_positions_right,/*donor*/nmismatches_right,
		     /*acceptor*/mismatch_positions_left,/*acceptor*/nmismatches_left,
		     /*D*/segmentj_sites,/*A*/segmenti_sites,
		     /*D*/segmentj_types,/*A*/segmenti_types,
		     /*D*/segmentj_probs,/*A*/segmenti_probs,
		     /*D*/segmentj_nsites,/*A*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	/* best_type = type; */
      }
    }
    
  } else {
    if (plusp == true) {
      /* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*A*/&best_nmismatches_i,/*D*/&best_nmismatches_j,
			 /*A*/&best_ref_nmismatches_i,/*D*/&best_ref_nmismatches_j,
			 /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/plusp,/*plusDp*/plusp,querylength,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*A*/segmenti_sites,/*D*/segmentj_sites,
			 /*A*/segmenti_types,/*D*/segmentj_types,
			 /*A*/segmenti_probs,/*D*/segmentj_probs,
			 /*A*/segmenti_nsites,/*D*/segmentj_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	/* best_type = type; */
      }
      
    } else {
      /* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*A*/&best_nmismatches_j,/*D*/&best_nmismatches_i,
			 /*A*/&best_ref_nmismatches_j,/*D*/&best_ref_nmismatches_i,
			 /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/false,/*plusDp*/false,querylength,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*A*/segmentj_sites,/*D*/segmenti_sites,
			 /*A*/segmentj_types,/*D*/segmenti_types,
			 /*A*/segmentj_probs,/*D*/segmenti_probs,
			 /*A*/segmentj_nsites,/*D*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	/* best_type = type; */
      }
    }
  }


  if (find_middle_p == false) {
    /* Skip */

  } else if (best_splice_qpos <= 0 &&
	     splice_qpos_qstart >= splice_qpos_qend + MIN_EXONLEN) {
    /* Possible middle exon or multiple sequencing errors */
    debug1(printf("splice_qpos_qstart %d >= splice_qpos_qend %d + %d, so trying to find a middle exon\n",
		  splice_qpos_qstart,splice_qpos_qend,MIN_EXONLEN));

    if (plusp == true) {
      n_middle_univdiagonals =
	Spliceends_middle_plus(&middle_univdiagonals,
			       stage1,streamspace_max_alloc,streamspace_alloc,
			       /*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
			       /*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
			       query_compress,queryptr,
			       univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
			       localdb_nmismatches_allowed);
      /* mismatch_positions_left is ascending */
      /* mismatch_positions_right is descending */

    } else {
      n_middle_univdiagonals =
	Spliceends_middle_minus(&middle_univdiagonals,
				stage1,streamspace_max_alloc,streamspace_alloc,
				/*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
				/*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
				query_compress,queryptr,
				univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
				localdb_nmismatches_allowed);
    }


    debug1(printf("The two ends cannot meet, and with knownsplicing trying to find a middle exon\n"));
    if (n_middle_univdiagonals == 0) {
      /* Skip: middle_univdiagonals is unassigned */

    } else {
      if (plusp == true) {
	/* Skip */
	check_ascending(mismatch_positions_left,nmismatches_left);
	check_descending(mismatch_positions_right,nmismatches_right);

      } else if (invertedp == true) {
	/* Already inverted */
	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

      } else {
	invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
	invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
	invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
	invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

	invertedp = true;
      }
      
      if ((middle_univdiagonal =
	   find_middle_exon(&typei,&typej,&splice_qpos_i,&splice_qpos_j,
			    &best_nmismatches_i,&nmismatches_middle,&best_nmismatches_j,
			    &best_ref_nmismatches_i,&ref_nmismatches_middle,&best_ref_nmismatches_j,
			    &best_donor1_prob,&best_acceptor1_prob,
			    &best_donor2_prob,&best_acceptor2_prob,
			    middle_univdiagonals,n_middle_univdiagonals,univdiagonal_i,univdiagonal_j,
			    pos5,pos3,/*splice_qpos_5*/splice_qpos_qend,/*splice_qpos_3*/splice_qpos_qstart,querylength,
			    
			    mismatch_positions_left,nmismatches_left,
			    mismatch_positions_right,nmismatches_right,
			    
			    segmenti_sites,segmentj_sites,
			    segmenti_types,segmentj_types,
			    segmenti_probs,segmentj_probs,
			    segmenti_nsites,segmentj_nsites,
			    
			    query_compress,chroffset,spliceinfo,knownsplicing,intlistpool,
			    plusp,sense_forward_p,genestrand)) != 0) {
	/* find_middle_exon returns values as qpos, not querypos */
	debug1(printf("(1) Splice_resolve_qstart (%s): found a middle exon %u at splice_qpos %d and %d\n",
		      sensedir == SENSE_FORWARD ? "sense" : "antisense",
		      middle_univdiagonal,splice_qpos_i,splice_qpos_j));

	FREE(middle_univdiagonals);
	
#if 0
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
#else
	*univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
						univcoordlistpool_trace(__FILE__,__LINE__));
	*univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_i
						univcoordlistpool_trace(__FILE__,__LINE__));
#endif
					      
	/* For qstart, push j first, then push i */
#if 0
	Intlist_head_set(newpath->nmismatches,nmismatches_j);
	Intlist_head_set(newpath->ref_nmismatches,ref_nmismatches_j);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_middle
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_middle
						    intlistpool_trace(__FILE__,__LINE__));
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_i
						    intlistpool_trace(__FILE__,__LINE__));
#else
	*nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_j
					intlistpool_trace(__FILE__,__LINE__));
	*nmismatches = Intlistpool_push(*nmismatches,intlistpool,nmismatches_middle
					intlistpool_trace(__FILE__,__LINE__));
	*nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_i
					intlistpool_trace(__FILE__,__LINE__));
				      
	*ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_j
					    intlistpool_trace(__FILE__,__LINE__));
	*ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,ref_nmismatches_middle
					    intlistpool_trace(__FILE__,__LINE__));
	*ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_i
					    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
	Intlist_head_set(newpath->endpoints,splice_qpos_j);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos_i
					      intlistpool_trace(__FILE__,__LINE__));
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos5
					      intlistpool_trace(__FILE__,__LINE__));
#else
	*endpoints = Intlistpool_push(NULL,intlistpool,splice_qpos_j
				      intlistpool_trace(__FILE__,__LINE__));
	*endpoints = Intlistpool_push(*endpoints,intlistpool,splice_qpos_i
				      intlistpool_trace(__FILE__,__LINE__));
	*endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos5
				      intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
	splice_distance_j = univdiagonal - middle_univdiagonal;
	splice_distance_i = middle_univdiagonal - low_univdiagonal;

	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance_j,sensedir,
									donor2_prob,acceptor2_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance_i,sensedir,
									donor1_prob,acceptor1_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
#else
	splice_distance_j = univdiagonal_j - middle_univdiagonal;
	splice_distance_i = middle_univdiagonal - univdiagonal_i;
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance_j,sensedir,
								best_donor2_prob,best_acceptor2_prob,
								pathpool)
				   listpool_trace(__FILE__,__LINE__));
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_splice(splice_distance_i,sensedir,
								best_donor1_prob,best_acceptor1_prob,
								pathpool)
				   listpool_trace(__FILE__,__LINE__));
#endif

	return 2;			/* 2 segments */

      } else {
	FREE(middle_univdiagonals);
	/* Proceed below */
      }
    }
  }

  debug1(printf("(1) Before spliceindel_resolve, have splice_qpos %d, nsites %d and %d, probs %f vs %f and %f vs %f\n",
		best_splice_qpos,segmenti_nsites,segmentj_nsites,
		best_donor1_prob,optimal_donor1_prob,best_acceptor2_prob,optimal_acceptor2_prob));
  if (allow_indels_p == false) {
    /* Skip */
    debug1(printf("(1) Indels not allowed\n"));

  } else if (segmenti_nsites == 0 || segmentj_nsites == 0) {
    debug1(printf("Cannot call spliceindel_resolve because nsites %d and %d\n",
		  segmenti_nsites,segmentj_nsites));

  } else if (best_donor1_prob == optimal_donor1_prob && best_acceptor2_prob == optimal_acceptor2_prob) {
    debug1(printf("No need to call spliceindel_resolve because exact probs are optimal\n"));
    
  } else {
    /* Try splice plus indel */

    if (invertedp == true) {
      /* Computes in qpos, not querypos, so undo inversions */
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = false;
    }

    debug1(printf("For comparison, exact nmismatches are %d + %d\n",
		  best_nmismatches_i,best_nmismatches_j));
    if ((splice_qpos =
	 spliceindel_resolve(&type,&spliceindel_trimpos5,&spliceindel_trimpos3,&nindels,&indel_pos,
			     &spliceindel_nmismatches_i,&spliceindel_nmismatches_j,&nmismatches_indel,
			     &spliceindel_ref_nmismatches_i,&spliceindel_ref_nmismatches_j,
			     &ref_nmismatches_indel,
			     &spliceindel_donor1_prob,&spliceindel_acceptor2_prob,
			     
			     univdiagonal_i,univdiagonal_j,
			     query_compress,plusp,chroffset,chrhigh,
			     
			     mismatch_positions_left,nmismatches_left,
			     mismatch_positions_right,nmismatches_right,
			     
			     segmenti_sites,segmentj_sites,segmenti_types,segmentj_types,
			     segmenti_probs,segmentj_probs,segmenti_nsites,segmentj_nsites,
			     optimal_segmenti_prob,optimal_segmentj_prob,
			     
			     pos5,pos3,querylength,indelinfo,sense_forward_p,genestrand,
			     trim5p,trim3p)) <= 0) {
      /* No spliceindel found */
      debug1(printf("Splice_resolve_qstart: No spliceindel found\n"));

    } else if (best_splice_qpos > 0 &&
	       spliceindel_nmismatches_i + nmismatches_indel + spliceindel_nmismatches_j >=
	       best_nmismatches_i + best_nmismatches_j) {
      /* Spliceindel not as good as exact splice */
      debug1(printf("Splice_resolve_qstart: Spliceindel not as good as exact splice\n"));

    } else if (abs(splice_qpos - indel_pos) - 4*nmismatches_indel < 0) {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    best_nmismatches_i,best_nmismatches_j));
      debug1(printf("Splice_resolve_qstart: Indel segment has too many mismatches\n"));
      /* Fall through to atypical splice */

    } else {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    best_nmismatches_i,best_nmismatches_j));
      debug1(printf("Splice_resolve_qstart: Spliceindel acceptable\n"));

      /* best_type = type; */
      trimpos5 = spliceindel_trimpos5;
      trimpos3 = spliceindel_trimpos3;
      best_nindels = nindels;
      best_indel_pos = indel_pos;
      best_nmismatches_i = spliceindel_nmismatches_i;
      best_nmismatches_indel = nmismatches_indel;
      best_nmismatches_j = spliceindel_nmismatches_j;
      best_ref_nmismatches_i = spliceindel_ref_nmismatches_i;
      best_ref_nmismatches_indel = ref_nmismatches_indel;
      best_ref_nmismatches_j = spliceindel_ref_nmismatches_j;
      best_donor1_prob = spliceindel_donor1_prob;
      best_acceptor2_prob = spliceindel_acceptor2_prob;
      best_splice_qpos = splice_qpos;
    }
  }

#if 0
  if (best_splice_qpos <= 0) {
    /* Try atypical splicing */
    debug1(printf("Computing splice sites\n"));
    
    /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
    if (plusp == true) {
      /* Skip */
    } else if (invertedp == true) {
      /* Skip */
    } else {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = true;
    }
    
    debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

    if (splice_qpos_i_low < splice_qpos_j_low) {
      splice_qpos_qstart = splice_qpos_i_low;
    } else {
      splice_qpos_qstart = splice_qpos_j_low;
    }

    if (splice_qpos_i_high > splice_qpos_j_high) {
      splice_qpos_qend = splice_qpos_i_high;
    } else {
      splice_qpos_qend = splice_qpos_j_high;
    }

    if (sense_forward_p == true) {
      if (plusp == true) {
	/* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*D*/&best_nmismatches_i,/*A*/&best_nmismatches_j,
			 /*D*/&best_ref_nmismatches_i,/*A*/&best_ref_nmismatches_j,
			 /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/true,/*plusAp*/true,querylength,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);
	  
	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*D*/&best_nmismatches_j,/*A*/&best_nmismatches_i,
			 /*D*/&best_ref_nmismatches_j,/*A*/&best_ref_nmismatches_i,
			 /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/false,/*plusAp*/false,querylength,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*splice_querypos_low*/querylength - splice_qpos_qend,
			 /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }

    } else {
      if (plusp == true) {
	/* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			     /*A*/&best_nmismatches_i,/*D*/&best_nmismatches_j,
			     /*A*/&best_ref_nmismatches_i,/*D*/&best_ref_nmismatches_j,
			     /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/true,/*plusDp*/true,querylength,
			     /*acceptor*/mismatch_positions_left,nmismatches_left,
			     /*donor*/mismatch_positions_right,nmismatches_right,
			     /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			     /*A*/&best_nmismatches_j,/*D*/&best_nmismatches_i,
			     /*A*/&best_ref_nmismatches_j,/*D*/&best_ref_nmismatches_i,
			     /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/false,/*plusDp*/false,querylength,
			     /*acceptor*/mismatch_positions_right,nmismatches_right,
			     /*donor*/mismatch_positions_left,nmismatches_left,
			     /*splice_querypos_low*/querylength - splice_qpos_qend,
			     /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }
    }
  }    
#endif


#if 0
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < best_splice_qpos) {
    lefti++;
  }
  best_nmismatches_i = best_ref_nmismatches_i = lefti;
  
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= best_splice_qpos) {
    righti++;
  }
  best_nmismatches_j = best_ref_nmismatches_j = righti;
#endif


  if (best_splice_qpos <= 0) {
    assert(middle_univdiagonal == 0);
    return 0;

  } else {
#if 0
    if (plusp == sense_forward_p) {
      best_probi = best_donor1_prob;
      best_probj = best_acceptor2_prob;
    } else {
      best_probi = best_acceptor2_prob;
      best_probj = best_donor1_prob;
    }
#endif
    debug1(printf("donor1 %f, acceptor2 %f\n",best_donor1_prob,best_acceptor2_prob));

    /* Check that the exon-intron boundaries are okay */
    if (invertedp == true) {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    }

#ifdef DEBUG1
    printf("Checking bounds\n");
    printf("%d mismatches on left:",nmismatches_left);
    for (int i = 0; i <= nmismatches_left; i++) {
      printf(" %d",mismatch_positions_left[i]);
    }
    printf("\n");

    printf("%d mismatches on right:",nmismatches_right);
    for (int i = 0; i <= nmismatches_right; i++) {
      printf(" %d",mismatch_positions_right[i]);
    }
    printf("\n");
#endif

    supporti = best_splice_qpos - trimpos5;
    supportj = trimpos3 - best_splice_qpos;

#ifdef USE_TEST1
    /* Benchmarking indicates this is generally worse, except cases of high sequence error */
    if (best_nindels == 0) {
      endpoint_i = endpoint_j = best_splice_qpos;
    } else if (best_splice_qpos < best_indel_pos) {
      endpoint_i = best_splice_qpos;
      endpoint_j = best_indel_pos;
    } else {
      endpoint_i = best_indel_pos;
      endpoint_j = best_splice_qpos;
    }

    /* Test 1 */
    if (endpoint_i <= splice_qpos_qend) {
      nmatches_trim_i = 1;
      nmismatches_trim_i = 0;
    } else {
      nmismatches_trim_i =
	Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					      query_compress,/*univdiagonal*/univdiagonal_i,querylength,
					      /*pos5*/splice_qpos_qend,/*pos3*/endpoint_i,
					      plusp,genestrand);
      nmatches_trim_i = (endpoint_i - splice_qpos_qend) - nmismatches_trim_i;
    }

    if (endpoint_j >= splice_qpos_qstart) {
      nmatches_trim_j = 1;
      nmismatches_trim_j = 0;
    } else {
      nmismatches_trim_j =
	Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					      query_compress,/*univdiagonal*/univdiagonal_j,querylength,
					      /*pos5*/endpoint_j,/*pos3*/splice_qpos_qstart,
					      plusp,genestrand);
      nmatches_trim_j = (splice_qpos_qstart - endpoint_j) - nmismatches_trim_j;
    }

    debug1(printf("Solution has %d matches and %d mismatches in trim_i\n",
		  nmatches_trim_i,nmismatches_trim_i));
    debug1(printf("Solution has %d matches and %d mismatches in trim_j\n",
		  nmatches_trim_j,nmismatches_trim_j));

    if (nmismatches_trim_i >= nmatches_trim_i) {
      return 0;

    } else if (nmismatches_trim_j >= nmatches_trim_j) {
      return 0;
    }
#endif

    /* Test 2 */
    /* Eventually want this */
    /* See what happens we continue on univdiagonal_j */
    score_splice = (supporti - best_nmismatches_i) - best_nmismatches_i;

    nmismatches_stay =
      Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					    query_compress,/*univdiagonal*/univdiagonal_j,querylength,
					    /*pos5*/trimpos5,/*pos3*/best_splice_qpos,
					    plusp,genestrand);
    nmatches_stay = (best_splice_qpos - trimpos5) - nmismatches_stay;
    score_stay = nmatches_stay - nmismatches_stay;
    debug1(printf("splice_qpos %d is outside transition %d..%d.  Staying on univdiagonal_j yields %d matches and %d mismatches => score %d vs splice score %d\n",
		  best_splice_qpos,splice_qpos_qend,splice_qpos_qstart,nmatches_stay,nmismatches_stay,score_stay,score_splice));
    if (score_stay >= score_splice - 6) {
      debug1(printf("Splice_resolve_qstart: Okay to stay on univdiagonal_j\n\n"));
      return 0;
    }

    if (Splice_accept_p(best_splice_qpos,querylength,
			best_donor1_prob,best_acceptor2_prob,
			univdiagonal_i,univdiagonal_j,
			supporti,best_nmismatches_i,supportj,best_nmismatches_j,
			plusp,sense_forward_p,innerp) == false) {
      debug1(printf("Splice_resolve_qstart on dist %u rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on probs\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,
		    supportj,best_nmismatches_j,best_acceptor2_prob));
      assert(middle_univdiagonal == 0);
      return 0;		/* Fail */

#if 1
    } else if (cont_trimpos >= 0 &&
	       (pos3 - trimpos5) - best_nmismatches_i - best_nmismatches_j - best_nmismatches_indel <= (pos3 - cont_trimpos) - cont_nmismatches) {
      debug1(printf("Solution yields fewer matches than continuation: (%d - %d) - %d - %d - %d vs (%d - %d) - %d\n\n",
		    pos3,trimpos5,best_nmismatches_i,best_nmismatches_j,best_nmismatches_indel,pos3,cont_trimpos,cont_nmismatches));
      return 0;
#endif

    } else if (best_nindels == 0) {
      /* Splice only */
      debug1(printf("(2) Splice_resolve_qstart (%s) on dist %u accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",/*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));

      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal_i
					      univcoordlistpool_trace(__FILE__,__LINE__));

      /* For qstart, push j first, then push i */
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));
      
      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
      
      *endpoints = Intlistpool_push(NULL,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos5
				    intlistpool_trace(__FILE__,__LINE__));

      splice_distance = univdiagonal_j - univdiagonal_i;
      *junctions = Listpool_push(NULL,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));
      
      return 1;			/* 1 segment */
      
    } else if (best_splice_qpos < best_indel_pos) {
      /* Push indel (based on left) then splice.  splice is distal, indel is medial */
      debug1(printf("(3) Splice_resolve_qstart (%s) on dist %u accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",/*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));

#if 0
      middle_univdiagonal = univdiagonal - nindels; /* nindels = univdiagonal - middle_univdiagonal */
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
#else
      middle_univdiagonal = univdiagonal_j - nindels; /* nindels = univdiagonal - middle_univdiagonal */
      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
      *univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_i
					      univcoordlistpool_trace(__FILE__,__LINE__));
#endif

      /* For qstart, push j first, then push i */
#if 0
      Intlist_head_set(newpath->nmismatches,nmismatches_j);
      Intlist_head_set(newpath->ref_nmismatches,nmismatches_j);
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_indel
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_indel
						  intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_i
						  intlistpool_trace(__FILE__,__LINE__));
#else
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_indel
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));

      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_indel
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      Intlist_head_set(newpath->endpoints,indel_pos);
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos
					    intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos5
					    intlistpool_trace(__FILE__,__LINE__));
#else
      *endpoints = Intlistpool_push(NULL,intlistpool,best_indel_pos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos5
				    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      /* Indel first */
      if (nindels < 0) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(-nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	/* deletionpos = low_left + splice_distance + indel_pos; -- qstart; add splice_distance if indel first */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      }
      
      /* Splice second */
      splice_distance = middle_univdiagonal - low_univdiagonal;
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
								      donor1_prob,acceptor2_prob,pathpool)
					 listpool_trace(__FILE__,__LINE__));
#else
      if (best_nindels < 0) {
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_insertion(-best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      } else {
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_deletion(best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      }
	  
      splice_distance = middle_univdiagonal - univdiagonal_i;
      *junctions = Listpool_push(*junctions,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));
#endif

      return 2;
	  
    } else {
      /* Push splice then indel (based on low_univdiagonal).  indel is distal, splice is medial */
      debug1(printf("(4) Splice_resolve_qstart (%s) on dist %u accepting splice pos %d, distance %u, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",
		    /*dist*/univdiagonal_j - univdiagonal_i,best_splice_qpos,
		    /*splice_distance*/univdiagonal_j - middle_univdiagonal,
		    best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));

#if 0
      middle_univdiagonal = low_univdiagonal + nindels; /* nindels = middle_univdiagonal - low_univdiagonal */
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,low_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
#else
      middle_univdiagonal = univdiagonal_i + nindels;
      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
      *univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_i
					      univcoordlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      Intlist_head_set(newpath->nmismatches,nmismatches_j);
      Intlist_head_set(newpath->ref_nmismatches,ref_nmismatches_j);
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_indel
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_indel
						  intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_i
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_i
						  intlistpool_trace(__FILE__,__LINE__));
#else
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_indel
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));
      
      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_indel
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      Intlist_head_set(newpath->endpoints,splice_qpos);
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos
					    intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos5   /* was low_qstart */
					    intlistpool_trace(__FILE__,__LINE__));
#else
      *endpoints = Intlistpool_push(NULL,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,best_indel_pos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos5
				    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      /* Splice first */
      splice_distance = univdiagonal - middle_univdiagonal;
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
								      donor1_prob,acceptor2_prob,pathpool)
					 listpool_trace(__FILE__,__LINE__));
	  
      /* Indel second */
      if (nindels < 0) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(-nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	/* deletionpos = low_left + indel_pos; -- qstart; do not add splice_distance if indel second */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      }
#else
      splice_distance = univdiagonal_j - middle_univdiagonal;
      *junctions = Listpool_push(NULL,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));
      
      if (best_nindels < 0) {
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_insertion(-best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      } else {
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_deletion(best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      }
#endif

      return 2;		/* 2 segments */
    }
  }
}


int
Splice_resolve_qend (Univcoordlist_T *univdiagonals, Intlist_T *nmismatches, Intlist_T *ref_nmismatches,
		     Intlist_T *endpoints, List_T *junctions,

		     Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		     Stage1_T stage1, Compress_T query_compress, char *queryptr,
		     bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,
		     
		     int streamspace_max_alloc, Univcoord_T *streamspace_alloc,

		     Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		     Localdb_T localdb, int localdb_nmismatches_allowed,

		     int pos5, int pos3, int querylength, int cont_trimpos, int cont_nmismatches,
		     Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     Univcoordlistpool_T univcoordlistpool, Intlistpool_T intlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,

		     int sensedir, int genestrand, bool trim5p, bool trim3p, bool innerp,
		     bool allow_indels_p, bool find_middle_p) {

  int type;
  Univcoord_T middle_univdiagonal = (Univcoord_T) 0;
  int best_nmismatches_i = querylength, best_nmismatches_j = querylength;
  int best_ref_nmismatches_i = querylength, best_ref_nmismatches_j = querylength;

  int best_nindels = 0, best_nmismatches_indel = 0, best_ref_nmismatches_indel = 0, ninserts;
  int best_indel_pos = -1;
  double best_donor1_prob = 0.0, best_acceptor2_prob = 0.0;
  int trimpos5, trimpos3;

  Chrpos_T splice_distance;

#if 1
  /* For middle exon */
  Chrpos_T splice_distance_i, splice_distance_j;
  double best_donor2_prob = 0.0, best_acceptor1_prob = 0.0;
  int splice_region_low, splice_region_high, splice_qpos_i, splice_qpos_j;
  int nmismatches_middle = querylength, ref_nmismatches_middle = querylength;
  int typei, typej;
  int n_middle_univdiagonals;
  Univcoord_T *middle_univdiagonals;
#endif  

  int best_splice_qpos, best_splice_querypos = -1;
  int splice_qpos_qstart, splice_qpos_qend;
  int splice_qpos_i_low, splice_qpos_i_high, splice_qpos_j_low, splice_qpos_j_high;
  int *mismatch_positions_left, *mismatch_positions_right;
  int nmismatches_left, nmismatches_right;
  char donor1, donor2, acceptor1, acceptor2;
  int ignore_nmismatches;

  int splice_qpos, indel_pos, spliceindel_trimpos5, spliceindel_trimpos3;
  int nindels, spliceindel_nmismatches_i, spliceindel_nmismatches_j, nmismatches_indel,
    spliceindel_ref_nmismatches_i, spliceindel_ref_nmismatches_j, ref_nmismatches_indel;
  double spliceindel_donor1_prob, spliceindel_acceptor2_prob;

  /* int best_nmismatches, nmismatches, segmenti_nmismatches, segmentj_nmismatches; */
  /* double best_probi, best_probj; */
  double optimal_donor1_prob, optimal_acceptor2_prob, optimal_segmenti_prob, optimal_segmentj_prob;
  /* double best_prob, donor_prob, acceptor_prob; */

  int *segmenti_sites, *segmentj_sites;
  int *segmenti_types, *segmentj_types;
  double *segmenti_probs, *segmentj_probs;
  int segmenti_nsites, segmentj_nsites;

  int supporti, supportj;
  
  int nmismatches_to_trimpos;
  int nmatches_trim_i, nmatches_trim_j,
    nmismatches_trim_i, nmismatches_trim_j, ignore;
  int endpoint_i, endpoint_j;

  int nmatches_stay, nmismatches_stay;
  int score_stay, score_splice;

  bool invertedp = false;
  bool sense_forward_p = (sensedir == SENSE_FORWARD);

#ifdef DEBUG1
  char *gbuffer1, *gbuffer2;
  static int ncalls = 0;
#endif

#if defined(DEBUG1) || defined(TRIM_AT_CHROMOSOME_BOUNDS)
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
#endif

#ifdef DEBUG1
  printf("Splice_resolve_qend (%s), plusp %d, compress_fwdp %d, sense_forward_p %d, call %d: Getting genome at univdiagonal_i %u [%u] and univdiagonal_j %u [%u] (diff: %d), range %d..%d\n",
	 sensedir == SENSE_FORWARD ? "sense" : "antisense",
	 plusp,Compress_fwdp(query_compress),sense_forward_p,++ncalls,
	 univdiagonal_i,univdiagonal_i - chroffset,univdiagonal_j,univdiagonal_j - chroffset,
	 univdiagonal_j-univdiagonal_i,pos5,pos3);
#endif
  assert(Compress_fwdp(query_compress) == plusp);
  assert(univdiagonal_j > univdiagonal_i);


#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonal_i + pos5 >= chroffset + querylength) ? pos5 : (int) (chroffset - segmenti_left);
  pos3 = (univdiagonal_j + pos3 <= chrhigh + querylength) ? pos3 : (int) (chrhigh - segmentj_left);
#endif


  /* Require separation from endpoints */
#if 0
  splice_qpos_start = pos5 + 1;
  splice_qpos_end = pos3 - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > (querylength - 1) - min_shortend) {
    splice_qpos_end = (querylength - 1) - min_shortend;
  }
#endif

  /* Determine feasibility of a splice */

  /* New method, different from relying on nmismatches_allowed */
  mismatch_positions_left = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_right = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  /* ascending */
  nmismatches_left =
    Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					    genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,pos5,pos3,plusp,genestrand);

  /* descending */
  nmismatches_right =
    Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					     genomebits,genomebits_alt,query_compress,
					     univdiagonal_j,querylength,pos5,pos3,plusp,genestrand);

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );


  /* was splice_qpos_high */
  splice_qpos_qend = Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_left,nmismatches_left,
						   pos5,pos3,querylength);
  /* was splice_qpos_low */
  splice_qpos_qstart = Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_right,nmismatches_right,
						    pos5,pos3);
  debug1(printf("trim_nosplice for pos5 %d, pos3 %d yields low %d and high %d\n",
		pos5,pos3,splice_qpos_qstart,splice_qpos_qend));

  /* Trimming */
  if (trim5p == false) {
    trimpos5 = pos5;
  } else if (splice_qpos_qend - pos5 < 8) {
    trimpos5 = pos5;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qstart with %d..%d",pos5,splice_qpos_qend));
    trimpos5 = Genomebits_trim_qstart(&ignore_nmismatches,query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/splice_qpos_qend,plusp,genestrand);
    debug1(printf(" to yield %d..\n",trimpos5));
  }

  if (trim3p == false) {
    trimpos3 = pos3;
  } else if (pos3 - splice_qpos_qstart < 8) {
    trimpos3 = pos3;
  } else {
    debug1(printf("(1) Calling Genomebits_trim_qend with %d..%d",splice_qpos_qstart,pos3));
    trimpos3 = Genomebits_trim_qend(&ignore_nmismatches,query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/splice_qpos_qstart,pos3,plusp,genestrand);
    debug1(printf(" to yield ..%d\n",trimpos3));
  }

  if (trimpos5 != pos5) {
    /* ascending */
    nmismatches_left =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					      genomebits,genomebits_alt,query_compress,
					      univdiagonal_i,querylength,trimpos5,pos3,plusp,genestrand);
  }

  if (trimpos3 != pos3) {
    /* descending */
    nmismatches_right =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					       genomebits,genomebits_alt,query_compress,
					       univdiagonal_j,querylength,pos5,trimpos3,plusp,genestrand);
  }

  pos5 = trimpos5;
  pos3 = trimpos3;

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

#if 0
  /* From Splice_nomiddle */
  if (splice_qpos_qstart < pos5 + 1) {
    splice_qpos_qstart = pos5 + 1;
  }
  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_qstart));

  if (splice_qpos_qend > pos3 - 1) {
    splice_qpos_qend = pos3 - 1;
  }
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_qend));
#endif

  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_qstart));
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_qend));


#ifdef DEBUG1
  gbuffer1 = (char *) CALLOC(querylength+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonal_i - querylength,querylength,gbuffer1);
  Genome_fill_buffer(univdiagonal_j - querylength,querylength,gbuffer2);
  char *queryseq = Compress_queryseq(query_compress,querylength);

  printf("Splice_resolve_qend, plusp %d: Getting genome at left %llu and %llu\n",
	 plusp,(unsigned long long) segmenti_left,(unsigned long long) segmentj_left);
  printf("g1: %s\n",gbuffer1);
  printf("    ");
  print_diffs(gbuffer1,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbuffer2,queryseq,querylength);
  printf("\n");
  printf("g2: %s\n",gbuffer2);
  FREE(queryseq);
  FREE(gbuffer2);
  FREE(gbuffer1);
#endif


  debug1(printf("find_middle_p: %d\n",find_middle_p));

#if 1
  if (0 && find_middle_p == false) {
    /* From Splice_nomiddle */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

  } else if (splice_qpos_qstart <= splice_qpos_qend + MISMATCHES_AT_SITE) {
    /* Appears to have matches to both splice sites.  Span just the crossing parts */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

    /* Set both to min */
    if (splice_qpos_i_low < splice_qpos_j_low) {
      splice_qpos_j_low = splice_qpos_i_low;
    } else {
      splice_qpos_i_low = splice_qpos_j_low;
    }

    /* Set both to max */
    if (splice_qpos_i_high > splice_qpos_j_high) {
      splice_qpos_j_high = splice_qpos_i_high;
    } else {
      splice_qpos_i_high = splice_qpos_j_high;
    }

  } else {
    /* There is a gap in matches.  Span the entire region */
    splice_qpos_i_low = splice_qpos_j_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);
  }
#else
  if (splice_qpos_qstart <= splice_qpos_qend + MISMATCHES_AT_SITE) {
    /* Appears to have matches to both splice sites.  Span just the crossing parts */
    splice_qpos_i_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = add_bounded(splice_qpos_qend,EXTRA_MISMATCHES,pos3 - 1);
    splice_qpos_j_low = subtract_bounded(splice_qpos_qstart,EXTRA_MISMATCHES,pos5 + 1);
    splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);

  } else {
    /* There is a gap in matches.  Span the entire region */
    splice_qpos_i_low = splice_qpos_j_low = subtract_bounded(splice_qpos_qend,EXTRA_MATCHES,pos5 + 1);
    splice_qpos_i_high = splice_qpos_j_high = add_bounded(splice_qpos_qstart,EXTRA_MATCHES,pos3 - 1);
  }

  splice_region_low = /*min*/ (splice_qpos_i_low < splice_qpos_j_low) ? splice_qpos_i_low : splice_qpos_j_low;
  splice_region_high = /*max*/ (splice_qpos_i_high > splice_qpos_j_high) ? splice_qpos_i_high : splice_qpos_j_high;

  debug1(printf("  Splice region: %d..%d\n",splice_region_low,splice_region_high));
#endif

  debug1(printf("  Splice i bounds: %d..%d\n",splice_qpos_i_low,splice_qpos_i_high));
  debug1(printf("  Splice j bounds: %d..%d\n",splice_qpos_j_low,splice_qpos_j_high));


  if (plusp == sense_forward_p) {
    /* geneplus */
    segmenti_nsites = compute_donor_sites(&optimal_donor1_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
					  spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					  spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					  spliceinfo->segmenti_probs_alloc,
#if 0
					  splice_region_low,splice_region_high,
#else
					  splice_qpos_i_low,splice_qpos_i_high,
#endif
					  querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_acceptor_sites(&optimal_acceptor2_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					     spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					     spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					     spliceinfo->segmentj_probs_alloc,
#if 0
					     splice_region_low,splice_region_high,
#else
					     splice_qpos_j_low,splice_qpos_j_high,
#endif
					     querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_donor1_prob;
    optimal_segmentj_prob = optimal_acceptor2_prob;

  } else {
    /* geneminus */
    segmenti_nsites = compute_antiacceptor_sites(&optimal_acceptor2_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
						 spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						 spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						 spliceinfo->segmenti_probs_alloc,
#if 0
						 splice_region_low,splice_region_high,
#else
						 splice_qpos_i_low,splice_qpos_i_high,
#endif
						 querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_antidonor_sites(&optimal_donor1_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					      spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					      spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					      spliceinfo->segmentj_probs_alloc,
#if 0
					      splice_region_low,splice_region_high,
#else
					      splice_qpos_j_low,splice_qpos_j_high,
#endif
					      querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_acceptor2_prob;
    optimal_segmentj_prob = optimal_donor1_prob;
  }


  best_splice_qpos = -1;

  /* Try standard splicing */
  debug1(printf("Computing splice sites\n"));
    
  /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
  if (plusp == true) {
    /* Skip */
    check_ascending(mismatch_positions_left,nmismatches_left);
    check_descending(mismatch_positions_right,nmismatches_right);

#if 0
  } else if (invertedp == true) {
    /* Already inverted */
    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);
#endif

  } else {
    invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
    invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
    invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

    check_descending(mismatch_positions_left,nmismatches_left);
    check_ascending(mismatch_positions_right,nmismatches_right);

    invertedp = true;
  }
    
  debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

  if (sense_forward_p == true) {
    if (plusp == true) {
      /* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
		     /*D*/&best_nmismatches_i,/*A*/&best_nmismatches_j,
		     /*D*/&best_ref_nmismatches_i,/*A*/&best_ref_nmismatches_j,
		     /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/true,/*plusAp*/true,querylength,
		     /*donor*/mismatch_positions_left,/*donor*/nmismatches_left,
		     /*acceptor*/mismatch_positions_right,/*acceptor*/nmismatches_right,
		     /*D*/segmenti_sites,/*A*/segmentj_sites,
		     /*D*/segmenti_types,/*A*/segmentj_types,
		     /*D*/segmenti_probs,/*A*/segmentj_probs,
		     /*D*/segmenti_nsites,/*A*/segmentj_nsites);
	
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	/* best_type = type; */
      }
      
    } else {
      /* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
		     /*D*/&best_nmismatches_j,/*A*/&best_nmismatches_i,
		     /*D*/&best_ref_nmismatches_j,/*A*/&best_ref_nmismatches_i,
		     /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/false,/*plusAp*/false,querylength,
		     /*donor*/mismatch_positions_right,/*donor*/nmismatches_right,
		     /*acceptor*/mismatch_positions_left,/*acceptor*/nmismatches_left,
		     /*D*/segmentj_sites,/*A*/segmenti_sites,
		     /*D*/segmentj_types,/*A*/segmenti_types,
		     /*D*/segmentj_probs,/*A*/segmenti_probs,
		     /*D*/segmentj_nsites,/*A*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	/* best_type = type; */
      }
    }
    
  } else {
    if (plusp == true) {
      /* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*A*/&best_nmismatches_i,/*D*/&best_nmismatches_j,
			 /*A*/&best_ref_nmismatches_i,/*D*/&best_ref_nmismatches_j,
			 /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/plusp,/*plusDp*/plusp,querylength,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*A*/segmenti_sites,/*D*/segmentj_sites,
			 /*A*/segmenti_types,/*D*/segmentj_types,
			 /*A*/segmenti_probs,/*D*/segmentj_probs,
			 /*A*/segmenti_nsites,/*D*/segmentj_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	/* best_type = type; */
      }
      
    } else {
      /* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*A*/&best_nmismatches_j,/*D*/&best_nmismatches_i,
			 /*A*/&best_ref_nmismatches_j,/*D*/&best_ref_nmismatches_i,
			 /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/false,/*plusDp*/false,querylength,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*A*/segmentj_sites,/*D*/segmenti_sites,
			 /*A*/segmentj_types,/*D*/segmenti_types,
			 /*A*/segmentj_probs,/*D*/segmenti_probs,
			 /*A*/segmentj_nsites,/*D*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	/* best_type = type; */
      }
    }
  }


  if (find_middle_p == false) {
    /* Skip */

  } else if (best_splice_qpos <= 0 &&
	     splice_qpos_qstart >= splice_qpos_qend + MIN_EXONLEN) {
    /* Possible middle exon or multiple sequencing errors */
    debug1(printf("splice_qpos_qstart %d >= splice_qpos_qend %d + %d, so trying to find a middle exon\n",
		  splice_qpos_qstart,splice_qpos_qend,MIN_EXONLEN));

    if (plusp == true) {
      n_middle_univdiagonals =
	Spliceends_middle_plus(&middle_univdiagonals,
			       stage1,streamspace_max_alloc,streamspace_alloc,
			       /*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
			       /*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
			       query_compress,queryptr,
			       univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
			       localdb_nmismatches_allowed);
      /* mismatch_positions_left is ascending */
      /* mismatch_positions_right is descending */

    } else {
      n_middle_univdiagonals =
	Spliceends_middle_minus(&middle_univdiagonals,
				stage1,streamspace_max_alloc,streamspace_alloc,
				/*qstart*/splice_qpos_qend,/*qend*/splice_qpos_qstart,querylength,
				/*low_univdiagonal*/univdiagonal_i + 1,/*high_univdiagonal*/univdiagonal_j - 1,
				query_compress,queryptr,
				univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb,
				localdb_nmismatches_allowed);
    }


    debug1(printf("The two ends cannot meet, and with knownsplicing trying to find a middle exon\n"));
    if (n_middle_univdiagonals == 0) {
      /* Skip: middle_univdiagonals is unassigned */

    } else {
      if (plusp == true) {
	/* Skip */
	check_ascending(mismatch_positions_left,nmismatches_left);
	check_descending(mismatch_positions_right,nmismatches_right);

      } else if (invertedp == true) {
	/* Already inverted */
	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

      } else {
	invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
	invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
	invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
	invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);

	check_descending(mismatch_positions_left,nmismatches_left);
	check_ascending(mismatch_positions_right,nmismatches_right);

	invertedp = true;
      }
      
      if ((middle_univdiagonal =
	   find_middle_exon(&typei,&typej,&splice_qpos_i,&splice_qpos_j,
			    &best_nmismatches_i,&nmismatches_middle,&best_nmismatches_j,
			    &best_ref_nmismatches_i,&ref_nmismatches_middle,&best_ref_nmismatches_j,
			    &best_donor1_prob,&best_acceptor1_prob,
			    &best_donor2_prob,&best_acceptor2_prob,
			    middle_univdiagonals,n_middle_univdiagonals,univdiagonal_i,univdiagonal_j,
			    pos5,pos3,/*splice_qpos_5*/splice_qpos_qend,/*splice_qpos_3*/splice_qpos_qstart,querylength,
			    
			    mismatch_positions_left,nmismatches_left,
			    mismatch_positions_right,nmismatches_right,
			    
			    segmenti_sites,segmentj_sites,
			    segmenti_types,segmentj_types,
			    segmenti_probs,segmentj_probs,
			    segmenti_nsites,segmentj_nsites,
			    
			    query_compress,chroffset,spliceinfo,knownsplicing,intlistpool,
			    plusp,sense_forward_p,genestrand)) != 0) {
	/* find_middle_exon returns values as qpos, not querypos */
	debug1(printf("(1) Splice_resolve_qend (%s): found a middle exon %u at splice qpos %d and %d\n",
		      sensedir == SENSE_FORWARD ? "sense" : "antisense",
		      middle_univdiagonal,splice_qpos_i,splice_qpos_j));

	FREE(middle_univdiagonals);
	
#if 0
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
	newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
							univcoordlistpool_trace(__FILE__,__LINE__));
#else
	*univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
						univcoordlistpool_trace(__FILE__,__LINE__));
	*univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_j
						univcoordlistpool_trace(__FILE__,__LINE__));
#endif

	/* For qend, push i first, then push j */
#if 0
	Intlist_head_set(newpath->nmismatches,nmismatches_i);
	Intlist_head_set(newpath->ref_nmismatches,ref_nmismatches_i);
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_middle
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_middle
						    intlistpool_trace(__FILE__,__LINE__));
	newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
						intlistpool_trace(__FILE__,__LINE__));
	newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_j
						    intlistpool_trace(__FILE__,__LINE__));
#else    
	*nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_i
					intlistpool_trace(__FILE__,__LINE__));
	*nmismatches = Intlistpool_push(*nmismatches,intlistpool,nmismatches_middle
					intlistpool_trace(__FILE__,__LINE__));
	*nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_j
					intlistpool_trace(__FILE__,__LINE__));
  
	*ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_i
					    intlistpool_trace(__FILE__,__LINE__));
	*ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,ref_nmismatches_middle
					    intlistpool_trace(__FILE__,__LINE__));
	*ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_j
					    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
	Intlist_head_set(newpath->endpoints,splice_qpos_i);
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos_j
					      intlistpool_trace(__FILE__,__LINE__));
	newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos3
					      intlistpool_trace(__FILE__,__LINE__));
#else
	*endpoints = Intlistpool_push(NULL,intlistpool,splice_qpos_i
				      intlistpool_trace(__FILE__,__LINE__));
	*endpoints = Intlistpool_push(*endpoints,intlistpool,splice_qpos_j
				      intlistpool_trace(__FILE__,__LINE__));
	*endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos3
				      intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
	splice_distance_j = high_univdiagonal - middle_univdiagonal;
	splice_distance_i = middle_univdiagonal - univdiagonal;
    
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance_i,sensedir,
									donor1_prob,acceptor1_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_splice(splice_distance_j,sensedir,
									donor2_prob,acceptor2_prob,
									pathpool)
					   listpool_trace(__FILE__,__LINE__));
#else
	splice_distance_j = univdiagonal_j - middle_univdiagonal;
	splice_distance_i = middle_univdiagonal - univdiagonal_i;
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_splice(splice_distance_i,sensedir,
								best_donor2_prob,best_acceptor2_prob,
								pathpool)
				   listpool_trace(__FILE__,__LINE__));
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_splice(splice_distance_j,sensedir,
								best_donor1_prob,best_acceptor1_prob,
								pathpool)
				   listpool_trace(__FILE__,__LINE__));
#endif
    
	return 2;

      } else {
	FREE(middle_univdiagonals);
	/* Proceed below */
      }
    }
  }

  debug1(printf("(1) Before spliceindel_resolve, have splice_qpos %d, nsites %d and %d, probs %f vs %f and %f vs %f\n",
		best_splice_qpos,segmenti_nsites,segmentj_nsites,
		best_donor1_prob,optimal_donor1_prob,best_acceptor2_prob,optimal_acceptor2_prob));
  if (allow_indels_p == false) {
    /* Skip */
    debug1(printf("(1) Indels not allowed\n"));

  } else if (segmenti_nsites == 0 || segmentj_nsites == 0) {
    debug1(printf("Cannot call spliceindel_resolve because nsites %d and %d\n",
		  segmenti_nsites,segmentj_nsites));

  } else if (best_donor1_prob == optimal_donor1_prob && best_acceptor2_prob == optimal_acceptor2_prob) {
    debug1(printf("No need to call spliceindel_resolve because exact probs are optimal\n"));
    
  } else {
    /* Try splice plus indel */

    if (invertedp == true) {
      /* Computes in qpos, not querypos, so undo inversions */
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = false;
    }

    debug1(printf("For comparison, exact nmismatches are %d + %d\n",
		  best_nmismatches_i,best_nmismatches_j));
    if ((splice_qpos =
	 spliceindel_resolve(&type,&spliceindel_trimpos5,&spliceindel_trimpos3,&nindels,&indel_pos,
			     &spliceindel_nmismatches_i,&spliceindel_nmismatches_j,&nmismatches_indel,
			     &spliceindel_ref_nmismatches_i,&spliceindel_ref_nmismatches_j,
			     &ref_nmismatches_indel,
			     &spliceindel_donor1_prob,&spliceindel_acceptor2_prob,
			     
			     univdiagonal_i,univdiagonal_j,
			     query_compress,plusp,chroffset,chrhigh,
			     
			     mismatch_positions_left,nmismatches_left,
			     mismatch_positions_right,nmismatches_right,
			     
			     segmenti_sites,segmentj_sites,segmenti_types,segmentj_types,
			     segmenti_probs,segmentj_probs,segmenti_nsites,segmentj_nsites,
			     optimal_segmenti_prob,optimal_segmentj_prob,
			     
			     pos5,pos3,querylength,indelinfo,sense_forward_p,genestrand,
			     trim5p,trim3p)) <= 0) {
      /* No spliceindel found */
      debug1(printf("Splice_resolve_qend: No spliceindel found\n"));

    } else if (best_splice_qpos > 0 &&
	       spliceindel_nmismatches_i + nmismatches_indel + spliceindel_nmismatches_j >=
	       best_nmismatches_i + best_nmismatches_j) {
      /* Spliceindel not as good as exact splice */
      debug1(printf("Splice_resolve_qend: Spliceindel not as good as exact splice\n"));

    } else if (abs(splice_qpos - indel_pos) - 4*nmismatches_indel < 0) {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    best_nmismatches_i,best_nmismatches_j));
      debug1(printf("Splice_resolve_qend: Indel segment has too many mismatches\n"));
      /* Fall through to atypical splice */

    } else {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    best_nmismatches_i,best_nmismatches_j));
      debug1(printf("Splice_resolve_qend: Spliceindel acceptable\n"));

      /* best_type = type; */
      trimpos5 = spliceindel_trimpos5;
      trimpos3 = spliceindel_trimpos3;
      best_nindels = nindels;
      best_indel_pos = indel_pos;
      best_nmismatches_i = spliceindel_nmismatches_i;
      best_nmismatches_indel = nmismatches_indel;
      best_nmismatches_j = spliceindel_nmismatches_j;
      best_ref_nmismatches_i = spliceindel_ref_nmismatches_i;
      best_ref_nmismatches_indel = ref_nmismatches_indel;
      best_ref_nmismatches_j = spliceindel_ref_nmismatches_j;
      best_donor1_prob = spliceindel_donor1_prob;
      best_acceptor2_prob = spliceindel_acceptor2_prob;
      best_splice_qpos = splice_qpos;
    }
  }

#if 0
  if (best_splice_qpos <= 0) {
    /* Try atypical splicing */
    debug1(printf("Computing splice sites\n"));
    
    /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
    if (plusp == true) {
      /* Skip */
    } else if (invertedp == true) {
      /* Skip */
    } else {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = true;
    }
    
    debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

    if (splice_qpos_i_low < splice_qpos_j_low) {
      splice_qpos_qstart = splice_qpos_i_low;
    } else {
      splice_qpos_qstart = splice_qpos_j_low;
    }

    if (splice_qpos_i_high > splice_qpos_j_high) {
      splice_qpos_qend = splice_qpos_i_high;
    } else {
      splice_qpos_qend = splice_qpos_j_high;
    }

    if (sense_forward_p == true) {
      if (plusp == true) {
	/* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*D*/&best_nmismatches_i,/*A*/&best_nmismatches_j,
			 /*D*/&best_ref_nmismatches_i,/*A*/&best_ref_nmismatches_j,
			 /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/true,/*plusAp*/true,querylength,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);
	  
	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			 /*D*/&best_nmismatches_j,/*A*/&best_nmismatches_i,
			 /*D*/&best_ref_nmismatches_j,/*A*/&best_ref_nmismatches_i,
			 /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/false,/*plusAp*/false,querylength,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*splice_querypos_low*/querylength - splice_qpos_qend,
			 /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }

    } else {
      if (plusp == true) {
	/* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			     /*A*/&best_nmismatches_i,/*D*/&best_nmismatches_j,
			     /*A*/&best_ref_nmismatches_i,/*D*/&best_ref_nmismatches_j,
			     /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/true,/*plusDp*/true,querylength,
			     /*acceptor*/mismatch_positions_left,nmismatches_left,
			     /*donor*/mismatch_positions_right,nmismatches_right,
			     /*splice_querypos_low*/splice_qpos_qstart,/*splice_querypos_high*/splice_qpos_qend);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&best_donor1_prob,&best_acceptor2_prob,
			     /*A*/&best_nmismatches_j,/*D*/&best_nmismatches_i,
			     /*A*/&best_ref_nmismatches_j,/*D*/&best_ref_nmismatches_i,
			     /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/false,/*plusDp*/false,querylength,
			     /*acceptor*/mismatch_positions_right,nmismatches_right,
			     /*donor*/mismatch_positions_left,nmismatches_left,
			     /*splice_querypos_low*/querylength - splice_qpos_qend,
			     /*splice_querypos_high*/querylength - splice_qpos_qstart);

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }
    }
  }    
#endif


#if 0
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < best_splice_qpos) {
    lefti++;
  }
  best_nmismatches_i = best_ref_nmismatches_i = lefti;
  
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= best_splice_qpos) {
    righti++;
  }
  best_nmismatches_j = best_ref_nmismatches_j = righti;
#endif


  if (best_splice_qpos <= 0) {
    assert(middle_univdiagonal == 0);
    return 0;			/* Fail */

  } else {
#if 0
    if (plusp == sense_forward_p) {
      best_probi = best_donor1_prob;
      best_probj = best_acceptor2_prob;
    } else {
      best_probi = best_acceptor2_prob;
      best_probj = best_donor1_prob;
    }
#endif
    debug1(printf("donor1 %f, acceptor2 %f\n",best_donor1_prob,best_acceptor2_prob));

    /* Check that the exon-intron boundaries are okay */
    if (invertedp == true) {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    }

#ifdef DEBUG1
    printf("Checking bounds\n");
    printf("%d mismatches on left:",nmismatches_left);
    for (int i = 0; i <= nmismatches_left; i++) {
      printf(" %d",mismatch_positions_left[i]);
    }
    printf("\n");

    printf("%d mismatches on right:",nmismatches_right);
    for (int i = 0; i <= nmismatches_right; i++) {
      printf(" %d",mismatch_positions_right[i]);
    }
    printf("\n");
#endif

    supporti = best_splice_qpos - trimpos5;
    supportj = trimpos3 - best_splice_qpos;

#ifdef USE_TEST1
    /* Benchmarking indicates this is generally worse, except cases of high sequence error */
    if (best_nindels == 0) {
      endpoint_i = endpoint_j = best_splice_qpos;
    } else if (best_splice_qpos < best_indel_pos) {
      endpoint_i = best_splice_qpos;
      endpoint_j = best_indel_pos;
    } else {
      endpoint_i = best_indel_pos;
      endpoint_j = best_splice_qpos;
    }

    /* Test 1 */
    if (endpoint_i <= splice_qpos_qend) {
      nmatches_trim_i = 1;
      nmismatches_trim_i = 0;
    } else {
      nmismatches_trim_i =
	Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					      query_compress,/*univdiagonal*/univdiagonal_i,querylength,
					      /*pos5*/splice_qpos_qend,/*pos3*/endpoint_i,
					      plusp,genestrand);
      nmatches_trim_i = (endpoint_i - splice_qpos_qend) - nmismatches_trim_i;
    }

    if (endpoint_j >= splice_qpos_qstart) {
      nmatches_trim_j = 1;
      nmismatches_trim_j = 0;
    } else {
      nmismatches_trim_j =
	Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					      query_compress,/*univdiagonal*/univdiagonal_j,querylength,
					      /*pos5*/endpoint_j,/*pos3*/splice_qpos_qstart,
					      plusp,genestrand);
      nmatches_trim_j = (splice_qpos_qstart - endpoint_j) - nmismatches_trim_j;
    }

    debug1(printf("Solution has %d matches and %d mismatches in trim_i\n",
		  nmatches_trim_i,nmismatches_trim_i));
    debug1(printf("Solution has %d matches and %d mismatches in trim_j\n",
		  nmatches_trim_j,nmismatches_trim_j));

    if (nmismatches_trim_i >= nmatches_trim_i) {
      return 0;

    } else if (nmismatches_trim_j >= nmatches_trim_j) {
      return 0;
    }
#endif

    /* Test 2 */
    /* Eventually want this */
    /* See what happens we continue on univdiagonal_i */
    score_splice = (supportj - best_nmismatches_j) - best_nmismatches_j;

    nmismatches_stay =
      Genomebits_count_mismatches_substring(&ignore,genomebits,genomebits_alt,
					    query_compress,/*univdiagonal*/univdiagonal_i,querylength,
					    /*pos5*/best_splice_qpos,/*pos3*/trimpos3,
					    plusp,genestrand);
    nmatches_stay = (trimpos3 - best_splice_qpos) - nmismatches_stay;
    score_stay = nmatches_stay - nmismatches_stay;
    debug1(printf("splice_qpos %d is outside transition %d..%d.  Staying on univdiagonal_j yields %d matches and %d mismatches => score %d vs splice score %d\n",
		  best_splice_qpos,splice_qpos_qend,splice_qpos_qstart,nmatches_stay,nmismatches_stay,score_stay,score_splice));
    if (score_stay >= score_splice - 6) {
      debug1(printf("Splice_resolve_qend: Okay to stay on univdiagonal_i\n\n"));
      return 0;
    }

    if (Splice_accept_p(best_splice_qpos,querylength,
			       best_donor1_prob,best_acceptor2_prob,
			       univdiagonal_i,univdiagonal_j,
			       supporti,best_nmismatches_i,supportj,best_nmismatches_j,
			       plusp,sense_forward_p,innerp) == false) {
      debug1(printf("Splice_resolve_qend on dist %u rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f, based on probs\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,
		    supportj,best_nmismatches_j,best_acceptor2_prob));
      assert(middle_univdiagonal == 0);
      return 0;		/* Fail */

#if 1
    } else if (cont_trimpos >= 0 &&
	       trimpos3 - best_nmismatches_i - best_nmismatches_j - best_nmismatches_indel <= cont_trimpos - cont_nmismatches) {
      debug1(printf("Solution yields fewer matches than continuation: %d - %d - %d - %d vs %d - %d\n\n",
		    trimpos3,best_nmismatches_i,best_nmismatches_j,best_nmismatches_indel,cont_trimpos,cont_nmismatches));
      return 0;
#endif

    } else if (best_nindels == 0) {
      /* Splice only */
      debug1(printf("(2) Splice_resolve_qend (%s) on dist %u accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",/*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));
      
      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal_j
					      univcoordlistpool_trace(__FILE__,__LINE__));

      /* For qend, push i first, then push j */
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));

      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));

      *endpoints = Intlistpool_push(NULL,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos3
				    intlistpool_trace(__FILE__,__LINE__));
      
      splice_distance = univdiagonal_j - univdiagonal_i;
      *junctions = Listpool_push(NULL,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));

      return 1;			/* 1 segment */
      
    } else if (best_indel_pos < best_splice_qpos) {
      /* Push indel (based on left) then splice.  indel is medial, splice is distal. */
      debug1(printf("(3) Splice_resolve_qend (%s) on dist %u accepting splice pos %d, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",/*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));

#if 0
      middle_univdiagonal = univdiagonal + nindels; /* nindels = middle_univdiagonal - univdiagonal */
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
#else
      middle_univdiagonal = univdiagonal_i + nindels;
      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
      *univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_j
					      univcoordlistpool_trace(__FILE__,__LINE__));
#endif      

      /* For qend, push i first, then push j */
#if 0
      Intlist_head_set(newpath->nmismatches,nmismatches_i);
      Intlist_head_set(newpath->ref_nmismatches,nmismatches_i);
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_indel
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_indel
						  intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_j
						  intlistpool_trace(__FILE__,__LINE__));
#else
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_indel
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));
				      
      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_indel
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));
#endif
      
#if 0
      Intlist_head_set(newpath->endpoints,indel_pos);
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,splice_qpos
					    intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos3   /* was high_qend */
					    intlistpool_trace(__FILE__,__LINE__));
#else      
      *endpoints = Intlistpool_push(NULL,intlistpool,best_indel_pos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos3
				    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      /* Indel first */
      if (nindels < 0) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(-nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	/* deletionpos = left + indel_pos; -- qend; do not add splice_distance if indel first */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      }
      
      /* Splice second */
      splice_distance = high_univdiagonal - middle_univdiagonal;
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
								      donor1_prob,acceptor2_prob,pathpool)
					 listpool_trace(__FILE__,__LINE__));
#else
      if (best_nindels < 0) {
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_insertion(-best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      } else {
	*junctions = Listpool_push(NULL,listpool,
				   (void *) Junction_new_deletion(best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      }
	  
      splice_distance = univdiagonal_j - middle_univdiagonal;
      *junctions = Listpool_push(*junctions,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));
#endif

      return 2;			/* 2 segments */
      
    } else {
      /* Push splice then indel (based on high_univdiagonal).  splice is medial, indel is distal. */
      debug1(printf("(4) Splice_resolve_qend (%s) on dist %u accepting splice pos %d, distance %u, nindels %d (%s), indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    sensedir == SENSE_FORWARD ? "sense" : "antisense",
		    /*dist*/univdiagonal_j - univdiagonal_i,best_splice_qpos,
		    /*splice_distance*/ middle_univdiagonal - univdiagonal_i,
		    best_nindels,(best_nindels < 0) ? "insertion" : (best_nindels > 0) ? "deletion" : "",
		    best_indel_pos,supporti,best_nmismatches_i,best_donor1_prob,supportj,best_nmismatches_j,best_acceptor2_prob));

#if 0
      middle_univdiagonal = high_univdiagonal - nindels; /* nindels = high_univdiagonal - middle_univdiagonal */
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,middle_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
      newpath->univdiagonals = Univcoordlistpool_push(newpath->univdiagonals,univcoordlistpool,high_univdiagonal
						      univcoordlistpool_trace(__FILE__,__LINE__));
#else
      middle_univdiagonal = univdiagonal_j - nindels;  /* nindels = middle_univdiagonal - low_univdiagonal */
      *univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,middle_univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
      *univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,univdiagonal_j
					      univcoordlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      Intlist_head_set(newpath->nmismatches,nmismatches_i);
      Intlist_head_set(newpath->ref_nmismatches,ref_nmismatches_i);
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_indel
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,nmismatches_indel
						  intlistpool_trace(__FILE__,__LINE__));
      newpath->nmismatches = Intlistpool_push(newpath->nmismatches,intlistpool,nmismatches_j
					      intlistpool_trace(__FILE__,__LINE__));
      newpath->ref_nmismatches = Intlistpool_push(newpath->ref_nmismatches,intlistpool,ref_nmismatches_j
						  intlistpool_trace(__FILE__,__LINE__));
#else      
      *nmismatches = Intlistpool_push(NULL,intlistpool,best_nmismatches_i
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_indel
				      intlistpool_trace(__FILE__,__LINE__));
      *nmismatches = Intlistpool_push(*nmismatches,intlistpool,best_nmismatches_j
				      intlistpool_trace(__FILE__,__LINE__));
      
      *ref_nmismatches = Intlistpool_push(NULL,intlistpool,best_ref_nmismatches_i
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_indel
					  intlistpool_trace(__FILE__,__LINE__));
      *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,best_ref_nmismatches_j
					  intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0
      Intlist_head_set(newpath->endpoints,splice_qpos);
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,indel_pos
					    intlistpool_trace(__FILE__,__LINE__));
      newpath->endpoints = Intlistpool_push(newpath->endpoints,intlistpool,trimpos3   /* was high_qend */
					    intlistpool_trace(__FILE__,__LINE__));
#else
      *endpoints = Intlistpool_push(NULL,intlistpool,best_splice_qpos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,best_indel_pos
				    intlistpool_trace(__FILE__,__LINE__));
      *endpoints = Intlistpool_push(*endpoints,intlistpool,trimpos3
				    intlistpool_trace(__FILE__,__LINE__));
#endif

#if 0      
      /* Splice first */
      splice_distance = middle_univdiagonal - univdiagonal;
      newpath->junctions = Listpool_push(newpath->junctions,listpool,
					 (void *) Junction_new_splice(splice_distance,SENSE_FORWARD,
								      donor1_prob,acceptor2_prob,pathpool)
					 listpool_trace(__FILE__,__LINE__));
      
      /* Indel second */
      if (nindels < 0) {
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_insertion(-nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      } else {
	/* deletionpos = left + splice_distance + indel_pos; -- qend; add splice_distance if indel second */
	newpath->junctions = Listpool_push(newpath->junctions,listpool,
					   (void *) Junction_new_deletion(nindels,pathpool)
					   listpool_trace(__FILE__,__LINE__));
      }
#else
      splice_distance = middle_univdiagonal - univdiagonal_i;
      *junctions = Listpool_push(NULL,listpool,
				 (void *) Junction_new_splice(splice_distance,sensedir,
							      best_donor1_prob,best_acceptor2_prob,pathpool)
				 listpool_trace(__FILE__,__LINE__));
      
      if (best_nindels < 0) {
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_insertion(-best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      } else {
	*junctions = Listpool_push(*junctions,listpool,
				   (void *) Junction_new_deletion(best_nindels,pathpool)
				   listpool_trace(__FILE__,__LINE__));
      }
#endif
      
      return 2;		/* 2 segments */
    }
  }
}


#if 1
/* Called by Path_qstart_resolve and Path_qend_resolve */
int
Splice_nomiddle (int *best_type, int *trimpos5, int *trimpos3, int *supporti, int *supportj,
		 int *best_nindels, int *best_indel_pos,
		 int *best_nmismatches_i, int *best_nmismatches_j,
		 int *best_nmismatches_indel,
		 
		 int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
		 int *best_ref_nmismatches_indel,
		 double *best_donor_prob, double *best_acceptor_prob,

		 Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		 Compress_T query_compress,
		 bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		 int pos5, int pos3, int querylength,
		 Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		 Intlistpool_T intlistpool,
		 bool sense_forward_p, int genestrand, bool trim5p, bool trim3p, bool innerp,
		 bool allow_indels_p) {

  int type;
  int best_splice_qpos, best_splice_querypos = -1;
  int splice_qpos_low, splice_qpos_high;
  int *mismatch_positions_left, *mismatch_positions_right;
  int nmismatches_left, nmismatches_right;
  char donor1, donor2, acceptor1, acceptor2;
  int ignore_nmismatches;

  int splice_qpos, indel_pos, spliceindel_trimpos5, spliceindel_trimpos3;
  int nindels, spliceindel_nmismatches_i, spliceindel_nmismatches_j, nmismatches_indel,
    spliceindel_ref_nmismatches_i, spliceindel_ref_nmismatches_j, ref_nmismatches_indel;
  double spliceindel_donor_prob, spliceindel_acceptor_prob;

  /* int best_nmismatches, nmismatches, segmenti_nmismatches, segmentj_nmismatches; */
  double best_probi, best_probj;
  double optimal_donor_prob, optimal_acceptor_prob, optimal_segmenti_prob, optimal_segmentj_prob;
  /* double best_prob, donor_prob, acceptor_prob; */

  int *segmenti_sites, *segmentj_sites;
  int *segmenti_types, *segmentj_types;
  double *segmenti_probs, *segmentj_probs;
  int segmenti_nsites, segmentj_nsites;

  int nmismatches_to_trimpos;

  bool invertedp = false;

#ifdef DEBUG1
  char *gbuffer1, *gbuffer2;
  static int ncalls = 0;
#endif


  *best_nindels = 0;
  *best_indel_pos = -1;
  *best_nmismatches_i = *best_nmismatches_j = querylength;
  *best_donor_prob = *best_acceptor_prob = 0.0;

#if defined(DEBUG1) || defined(TRIM_AT_CHROMOSOME_BOUNDS)
  Univcoord_T segmenti_left = univdiagonal_i - querylength;
  Univcoord_T segmentj_left = univdiagonal_j - querylength;
#endif

#ifdef DEBUG1
  printf("Splice_nomiddle plusp %d, compress_fwdp %d, sense_forward_p %d, call %d: Getting genome at univdiagonal_i %u [%u] and univdiagonal_j %u [%u] (diff: %d), range %d..%d\n",
	 plusp,Compress_fwdp(query_compress),sense_forward_p,++ncalls,
	 univdiagonal_i,univdiagonal_i - chroffset,univdiagonal_j,univdiagonal_j - chroffset,
	 univdiagonal_j-univdiagonal_i,pos5,pos3);
#endif
  assert(Compress_fwdp(query_compress) == plusp);
  assert(univdiagonal_j > univdiagonal_i);


#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonal_i + pos5 >= chroffset + querylength) ? pos5 : (int) (chroffset - segmenti_left);
  pos3 = (univdiagonal_j + pos3 <= chrhigh + querylength) ? pos3 : (int) (chrhigh - segmentj_left);
#endif


  /* Require separation from endpoints */
#if 0
  splice_qpos_start = pos5 + 1;
  splice_qpos_end = pos3 - 1;
  if (splice_qpos_start < min_shortend) {
    splice_qpos_start = min_shortend;
  }
  if (splice_qpos_end > (querylength - 1) - min_shortend) {
    splice_qpos_end = (querylength - 1) - min_shortend;
  }
#endif

  /* Determine feasibility of a splice */

  /* New method, different from relying on nmismatches_allowed */
  mismatch_positions_left = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_right = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  nmismatches_left =
    Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					    genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,pos5,pos3,plusp,genestrand);

  nmismatches_right =
    Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					     genomebits,genomebits_alt,query_compress,
					     univdiagonal_j,querylength,pos5,pos3,plusp,genestrand);

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );


  splice_qpos_high = Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_left,nmismatches_left,
						   pos5,pos3,querylength);
  splice_qpos_low = Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_right,nmismatches_right,
						    pos5,pos3);

  /* Trimming */
  if (trim5p == false) {
    *trimpos5 = pos5;
  } else if (splice_qpos_high - pos5 < 8) {
    *trimpos5 = pos5;
  } else {
    debug1(printf("(2) Calling Genomebits_trim_qstart with %d..%d",pos5,splice_qpos_high));
    *trimpos5 = Genomebits_trim_qstart(&ignore_nmismatches,query_compress,genomebits,
				       univdiagonal_i,querylength,
				       pos5,/*pos3*/splice_qpos_high,plusp,genestrand);
    debug1(printf(" to yield %d..\n",*trimpos5));
  }

  if (trim3p == false) {
    *trimpos3 = pos3;
  } else if (pos3 - splice_qpos_low < 8) {
    *trimpos3 = pos3;
  } else {
    debug1(printf("(2) Calling Genomebits_trim_qend with %d..%d",splice_qpos_low,pos3));
    *trimpos3 = Genomebits_trim_qend(&ignore_nmismatches,query_compress,genomebits,
				     univdiagonal_j,querylength,
				     /*pos5*/splice_qpos_low,pos3,plusp,genestrand);
    debug1(printf(" to yield ..%d\n",*trimpos3));
  }

  if (*trimpos5 != pos5) {
    nmismatches_left =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_left,/*max_mismatches*/pos3 - pos5,
					      genomebits,genomebits_alt,query_compress,
					      univdiagonal_i,querylength,*trimpos5,pos3,plusp,genestrand);
  }

  if (*trimpos3 != pos3) {
    nmismatches_right =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_right,/*max_mismatches*/pos3 - pos5,
					       genomebits,genomebits_alt,query_compress,
					       univdiagonal_j,querylength,pos5,*trimpos3,plusp,genestrand);
  }

  pos5 = *trimpos5;
  pos3 = *trimpos3;

  debug1(
	 printf("(3) %d mismatches on left from %d to %d at:",nmismatches_left,pos5,pos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  debug1(
	 printf("(3) %d mismatches on right from %d to %d at:",nmismatches_right,pos3,pos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );


  if (splice_qpos_low < pos5 + 1) {
    splice_qpos_low = pos5 + 1;
  }
  debug1(printf("  Splice low qpos bound for exact splice: %d\n",splice_qpos_low));

  if (splice_qpos_high > pos3 - 1) {
    splice_qpos_high = pos3 - 1;
  }
  debug1(printf("  Splice high qpos bound for exact splice: %d\n",splice_qpos_high));


#ifdef DEBUG1
  gbuffer1 = (char *) CALLOC(querylength+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonal_i - querylength,querylength,gbuffer1);
  Genome_fill_buffer(univdiagonal_j - querylength,querylength,gbuffer2);
  char *queryseq = Compress_queryseq(query_compress,querylength);

  printf("Splice_nomiddle, plusp %d: Getting genome at left %llu and %llu\n",
	 plusp,(unsigned long long) segmenti_left,(unsigned long long) segmentj_left);
  printf("g1: %s\n",gbuffer1);
  printf("    ");
  print_diffs(gbuffer1,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbuffer2,queryseq,querylength);
  printf("\n");
  printf("g2: %s\n",gbuffer2);
  FREE(queryseq);
  FREE(gbuffer2);
  FREE(gbuffer1);
#endif


#if 0
  printf("Splice_nomiddle computing %d..%d and %d..%d\n",
	 /*splice_qpos_low*/subtract_bounded(splice_qpos_high,EXTRA_MATCHES,pos5 + 1),
	 /*splice_qpos_high*/add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1),
	 /*splice_qpos_low*/subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1),
	 /*splice_qpos_high*/add_bounded(splice_qpos_low,EXTRA_MATCHES,pos3 - 1));
#endif

  if (plusp == sense_forward_p) {
    /* geneplus */
    segmenti_nsites = compute_donor_sites(&optimal_donor_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
					  spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					  spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					  spliceinfo->segmenti_probs_alloc,
					  /*splice_qpos_low*/subtract_bounded(splice_qpos_high,EXTRA_MATCHES,pos5 + 1),
					  /*splice_qpos_high*/add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1),
					  querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_acceptor_sites(&optimal_acceptor_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					     spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					     spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					     spliceinfo->segmentj_probs_alloc,
					     /*splice_qpos_low*/subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1),
					     /*splice_qpos_high*/add_bounded(splice_qpos_low,EXTRA_MATCHES,pos3 - 1),
					     querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_donor_prob;
    optimal_segmentj_prob = optimal_acceptor_prob;

  } else {
    /* geneminus */
    segmenti_nsites = compute_antiacceptor_sites(&optimal_acceptor_prob,&segmenti_sites,&segmenti_types,&segmenti_probs,
						 spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						 spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						 spliceinfo->segmenti_probs_alloc,
						 /*splice_qpos_low*/subtract_bounded(splice_qpos_high,EXTRA_MATCHES,pos5 + 1),
						 /*splice_qpos_high*/add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1),
						 querylength,univdiagonal_i,chroffset,knownsplicing,intlistpool);
    
    segmentj_nsites = compute_antidonor_sites(&optimal_donor_prob,&segmentj_sites,&segmentj_types,&segmentj_probs,
					      spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					      spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					      spliceinfo->segmentj_probs_alloc,
					      /*splice_qpos_low*/subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1),
					      /*splice_qpos_high*/add_bounded(splice_qpos_low,EXTRA_MATCHES,pos3 - 1),
					      querylength,univdiagonal_j,chroffset,knownsplicing,intlistpool);
    optimal_segmenti_prob = optimal_acceptor_prob;
    optimal_segmentj_prob = optimal_donor_prob;
  }


#if 0
  /* This code is specific to Splice_nomiddle, and not in Splice_resolve_qstart or Splice_resolve_qend */
  if (splice_qpos_low >= splice_qpos_high + MIN_EXONLEN) {
    /* Not allowing for middle exons */
    debug1(printf("Since splice_qpos_low %d >= splice_qpos_high %d + MIN_EXONLEN %d, returning -1\n",
		  splice_qpos_low,splice_qpos_high,MIN_EXONLEN));
    return -1;
  } else {
    best_splice_qpos = -1;
  }
#else
  best_splice_qpos = -1;
#endif

  /* Splice_nomiddle required splice_qpos_low <= splice_qpos_high + MISMATCHES_AT_SITE */
  
  /* Try standard splicing */
  debug1(printf("Computing splice sites\n"));
    
  /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
  if (plusp == false) {
    invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
    invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
    invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
    invertedp = true;
  }
    
  debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

  if (sense_forward_p == true) {
    if (plusp == true) {
      /* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
		     /*D*/&(*best_nmismatches_i),/*A*/&(*best_nmismatches_j),
		     /*D*/&(*best_ref_nmismatches_i),/*A*/&(*best_ref_nmismatches_j),
		     /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/true,/*plusAp*/true,querylength,
		     /*donor*/mismatch_positions_left,/*donor*/nmismatches_left,
		     /*acceptor*/mismatch_positions_right,/*acceptor*/nmismatches_right,
		     /*D*/segmenti_sites,/*A*/segmentj_sites,
		     /*D*/segmenti_types,/*A*/segmentj_types,
		     /*D*/segmenti_probs,/*A*/segmentj_probs,
		     /*D*/segmenti_nsites,/*A*/segmentj_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	*best_type = type;
      }

    } else {
      /* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_sense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
		     /*D*/&(*best_nmismatches_j),/*A*/&(*best_nmismatches_i),
		     /*D*/&(*best_ref_nmismatches_j),/*A*/&(*best_ref_nmismatches_i),
		     /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
		     /*plusDp*/false,/*plusAp*/false,querylength,
		     /*donor*/mismatch_positions_right,/*donor*/nmismatches_right,
		     /*acceptor*/mismatch_positions_left,/*acceptor*/nmismatches_left,
		     /*D*/segmentj_sites,/*A*/segmenti_sites,
		     /*D*/segmentj_types,/*A*/segmenti_types,
		     /*D*/segmentj_probs,/*A*/segmenti_probs,
		     /*D*/segmentj_nsites,/*A*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	*best_type = type;
      }
    }
    
  } else {
    if (plusp == true) {
      /* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			 /*A*/&(*best_nmismatches_i),/*D*/&(*best_nmismatches_j),
			 /*A*/&(*best_ref_nmismatches_i),/*D*/&(*best_ref_nmismatches_j),
			 /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/plusp,/*plusDp*/plusp,querylength,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*A*/segmenti_sites,/*D*/segmentj_sites,
			 /*A*/segmenti_types,/*D*/segmentj_types,
			 /*A*/segmenti_probs,/*D*/segmentj_probs,
			 /*A*/segmenti_nsites,/*D*/segmentj_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = best_splice_querypos;
	*best_type = type;
      }
      
    } else {
      /* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
      best_splice_querypos =
	splice_antisense(&type,&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			 /*A*/&(*best_nmismatches_j),/*D*/&(*best_nmismatches_i),
			 /*A*/&(*best_ref_nmismatches_j),/*D*/&(*best_ref_nmismatches_i),
			 /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			 /*plusAp*/false,/*plusDp*/false,querylength,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*A*/segmentj_sites,/*D*/segmenti_sites,
			 /*A*/segmentj_types,/*D*/segmenti_types,
			 /*A*/segmentj_probs,/*D*/segmenti_probs,
			 /*A*/segmentj_nsites,/*D*/segmenti_nsites);
      
      if (best_splice_querypos <= 0) {
	best_splice_qpos = -1;
      } else {
	best_splice_qpos = querylength - best_splice_querypos;
	*best_type = type;
      }
    }
  }


  debug1(printf("(2) Before spliceindel_resolve, have splice_qpos %d, nsites %d and %d, probs %f vs %f and %f vs %f\n",
		best_splice_qpos,segmenti_nsites,segmentj_nsites,
		*best_donor_prob,optimal_donor_prob,*best_acceptor_prob,optimal_acceptor_prob));

  if (allow_indels_p == false) {
    /* Skip */
    debug1(printf("(2) Indels not allowed\n"));

  } else if (segmenti_nsites == 0 || segmentj_nsites == 0) {
    debug1(printf("Cannot call spliceindel_resolve because nsites %d and %d\n",
		  segmenti_nsites,segmentj_nsites));

  } else if (*best_donor_prob == optimal_donor_prob && *best_acceptor_prob == optimal_acceptor_prob) {
    debug1(printf("No need to call spliceindel_resolve because exact probs are optimal\n"));

  } else {
    /* Try splice plus indel */
    /* Computes in qpos, not querypos, so no need for inversions */

    if (invertedp == true) {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = false;
    }

    debug1(printf("For comparison, exact nmismatches are %d + %d\n",
		  *best_nmismatches_i,*best_nmismatches_j));
    if ((splice_qpos =
	 spliceindel_resolve(&type,&spliceindel_trimpos5,&spliceindel_trimpos3,&nindels,&indel_pos,
			     &spliceindel_nmismatches_i,&spliceindel_nmismatches_j,&nmismatches_indel,
			     &spliceindel_ref_nmismatches_i,&spliceindel_ref_nmismatches_j,
			     &ref_nmismatches_indel,
			     &spliceindel_donor_prob,&spliceindel_acceptor_prob,
			  
			     univdiagonal_i,univdiagonal_j,
			     query_compress,plusp,chroffset,chrhigh,
			     
			     mismatch_positions_left,nmismatches_left,
			     mismatch_positions_right,nmismatches_right,
			     
			     segmenti_sites,segmentj_sites,segmenti_types,segmentj_types,
			     segmenti_probs,segmentj_probs,segmenti_nsites,segmentj_nsites,
			     optimal_segmenti_prob,optimal_segmentj_prob,
			     
			     pos5,pos3,querylength,indelinfo,sense_forward_p,genestrand,
			     trim5p,trim3p)) <= 0) {
      /* No spliceindel found */
      debug1(printf("Splice_nomiddle: No spliceindel found\n"));

    } else if (best_splice_qpos > 0 &&
	       spliceindel_nmismatches_i + nmismatches_indel + spliceindel_nmismatches_j >=
	       *best_nmismatches_i + *best_nmismatches_j) {
      /* Spliceindel not as good as exact splice */
      debug1(printf("Splice_nomiddle: Spliceindel not as good as exact splice\n"));

    } else if (abs(splice_qpos - indel_pos) - 4*nmismatches_indel < 0) {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    *best_nmismatches_i,*best_nmismatches_j));
      debug1(printf("Splice_nomiddle: Indel segment has too many mismatches\n"));

    } else {
      debug1(printf("Have nmismatches %d + %d + %d vs %d vs %d\n",
		    spliceindel_nmismatches_i,nmismatches_indel,spliceindel_nmismatches_j,
		    *best_nmismatches_i,*best_nmismatches_j));
      debug1(printf("Splice_nomiddle: Spliceindel acceptable\n"));

      *best_type = type;
      *trimpos5 = spliceindel_trimpos5;
      *trimpos3 = spliceindel_trimpos3;
      *best_nindels = nindels;
      *best_indel_pos = indel_pos;
      *best_nmismatches_i = spliceindel_nmismatches_i;
      *best_nmismatches_indel = nmismatches_indel;
      *best_nmismatches_j = spliceindel_nmismatches_j;
      *best_ref_nmismatches_i = spliceindel_ref_nmismatches_i;
      *best_ref_nmismatches_indel = ref_nmismatches_indel;
      *best_ref_nmismatches_j = spliceindel_ref_nmismatches_j;
      *best_donor_prob = spliceindel_donor_prob;
      *best_acceptor_prob = spliceindel_acceptor_prob;
      best_splice_qpos = splice_qpos;
    }
  }

#if 0
  /* Don't consider atypical splice for the ends of an alignment */
  if (best_splice_qpos <= 0 && splice_qpos_low <= splice_qpos_high + MISMATCHES_AT_SITE) {
    /* Try atypical splicing */
    debug1(printf("Computing splice sites\n"));
    
    /* Convert from qpos to querypos before calling splice_sense or splice_antisense */
    if (plusp == true) {
      /* Skip */
    } else if (invertedp == true) {
      /* Skip */
    } else {
      invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
      invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
      invert_sites(segmenti_sites,segmenti_types,segmenti_probs,segmenti_nsites,querylength);
      invert_sites(segmentj_sites,segmentj_types,segmentj_probs,segmentj_nsites,querylength);
      invertedp = true;
    }
    
    debug1(printf("plusp %d, sense_forward_p %d\n",plusp,sense_forward_p));

    if (sense_forward_p == true) {
      if (plusp == true) {
	/* sense, gplus: geneplus.  donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			 /*D*/&(*best_nmismatches_i),/*A*/&(*best_nmismatches_j),
			 /*D*/&(*best_ref_nmismatches_i),/*A*/&(*best_ref_nmismatches_j),
			 /*D*/univdiagonal_i,/*A*/univdiagonal_j,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/true,/*plusAp*/true,querylength,
			 /*donor*/mismatch_positions_left,nmismatches_left,
			 /*acceptor*/mismatch_positions_right,nmismatches_right,
			 subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1),
			 add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1));
	  
	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* sense, gminus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_sense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			 /*D*/&(*best_nmismatches_j),/*A*/&(*best_nmismatches_i),
			 /*D*/&(*best_ref_nmismatches_j),/*A*/&(*best_ref_nmismatches_i),
			 /*D*/univdiagonal_j,/*A*/univdiagonal_i,/*D*/chroffset,/*A*/chroffset,
			 /*plusDp*/false,/*plusAp*/false,querylength,
			 /*donor*/mismatch_positions_right,nmismatches_right,
			 /*acceptor*/mismatch_positions_left,nmismatches_left,
			 querylength - add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1),
			 querylength - subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1));

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }

    } else {
      if (plusp == true) {
	/* antisense, gplus: geneminus.  donor is univdiagonal_j and acceptor is univdiagonal_i */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			     /*A*/&(*best_nmismatches_i),/*D*/&(*best_nmismatches_j),
			     /*A*/&(*best_ref_nmismatches_i),/*D*/&(*best_ref_nmismatches_j),
			     /*A*/univdiagonal_i,/*D*/univdiagonal_j,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/true,/*plusDp*/true,querylength,
			     /*acceptor*/mismatch_positions_left,nmismatches_left,
			     /*donor*/mismatch_positions_right,nmismatches_right,
			     subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1),
			     add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1));

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = best_splice_querypos;
	}

      } else {
	/* antisense, gminus: geneplus, donor is univdiagonal_i and acceptor is univdiagonal_j */
	best_splice_querypos =
	  atypical_antisense(&donor1,&donor2,&acceptor1,&acceptor2,&(*best_donor_prob),&(*best_acceptor_prob),
			     /*A*/&(*best_nmismatches_j),/*D*/&(*best_nmismatches_i),
			     /*A*/&(*best_ref_nmismatches_j),/*D*/&(*best_ref_nmismatches_i),
			     /*A*/univdiagonal_j,/*D*/univdiagonal_i,/*A*/chroffset,/*D*/chroffset,
			     /*plusAp*/false,/*plusDp*/false,querylength,
			     /*acceptor*/mismatch_positions_right,nmismatches_right,
			     /*donor*/mismatch_positions_left,nmismatches_left,
			     querylength - add_bounded(splice_qpos_high,EXTRA_MISMATCHES,pos3 - 1),
			     querylength - subtract_bounded(splice_qpos_low,EXTRA_MISMATCHES,pos5 + 1));

	if (best_splice_querypos <= 0) {
	  best_splice_qpos = -1;
	} else {
	  best_splice_qpos = querylength - best_splice_querypos;
	}
      }
    }
  }    
#endif


#if 0
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < best_splice_qpos) {
    lefti++;
  }
  *best_nmismatches_i = *best_ref_nmismatches_i = lefti;
  
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= best_splice_qpos) {
    righti++;
  }
  *best_nmismatches_j = *best_ref_nmismatches_j = righti;
#endif


  if (best_splice_qpos <= 0) {
    return -1;			/* Fail */

  } else {
    if (plusp == sense_forward_p) {
      best_probi = *best_donor_prob;
      best_probj = *best_acceptor_prob;
    } else {
      best_probi = *best_acceptor_prob;
      best_probj = *best_donor_prob;
    }

    *supporti = best_splice_qpos - (*trimpos5);
    *supportj = (*trimpos3) - best_splice_qpos;
#if 0
    if (plusp == sense_forward_p) {
      *adj_donor_support = supporti - 4*(*best_nmismatches_i);
      *adj_acceptor_support = supportj - 4*(*best_nmismatches_j);
    } else {
      *adj_donor_support = supportj - 4*(*best_nmismatches_j);
      *adj_acceptor_support = supporti - 4*(*best_nmismatches_i);
    }
#endif

    /* This procedure is called for transcriptome-guided paths, but we still want to check for correctness */
    if (Splice_accept_p(best_splice_qpos,querylength,
			best_probi,best_probj,univdiagonal_i,univdiagonal_j,
			*supporti,*best_nmismatches_i,*supportj,*best_nmismatches_j,
			plusp,sense_forward_p,innerp) == false) {
      debug1(printf("Splice_nomiddle on dist %u rejecting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,*best_nindels,*best_indel_pos,*supporti,*best_nmismatches_i,best_probi,
		    *supportj,*best_nmismatches_j,best_probj));
      return -1;

    } else {
      debug1(printf("Splice_nomiddle on dist %u accepting splice pos %d, nindels %d, indel pos %d with supporti %d, mismatches %d, probi %f and supportj %d, mismatches %d, probj %f\n\n",
		    /*splice_distance*/univdiagonal_j - univdiagonal_i,
		    best_splice_qpos,*best_nindels,*best_indel_pos,*supporti,*best_nmismatches_i,best_probi,
		    *supportj,*best_nmismatches_j,best_probj));
      return best_splice_qpos;
    }
  }
}
#endif



#if 0
/* Replaced by Splice_fusion_sense and Splice_fusion_antisense */
int
Splice_resolve_fusion (char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		       int *best_nindels, int *best_indel_pos,
		       int *best_nmismatches_5, int *best_nmismatches_3,

		       int *best_ref_nmismatches_5, int *best_ref_nmismatches_3,
		       double *best_donor_prob, double *best_acceptor_prob,
		       
		       Univcoord_T univdiagonal5, Univcoord_T univdiagonal3,
		       Compress_T query_compress_5, bool plus5p, Univcoord_T chroffset5,
		       Compress_T query_compress_3, bool plus3p, Univcoord_T chroffset3,
		     
		       int querypos5, int querypos3, int querylength,
		       Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		       bool sense_forward_p, int genestrand,
		       int nmismatches_allowed, int max_insertionlen, int max_deletionlen) {

  int splice_querypos_low, splice_querypos_high;
  int *mismatch_positions_left, *mismatch_positions_right;
  int nmismatches_left, nmismatches_right;
  int mismatch5_adj, mismatch3_adj;
  
  int best_splice_querypos;

  /* int best_nmismatches, nmismatches, segment5_nmismatches, segment3_nmismatches; */
  double best_prob5, best_prob3;
  /* double best_prob, donor_prob, acceptor_prob; */

  int segment5_nsites, segment3_nsites;
  int *segment5_sites, *segment3_sites;
  int *segment5_types, *segment3_types;

  int adj_support5, adj_support3, support5, support3;
  Univcoord_T segment5_left, segment3_left;
 
 
  /* TODO: FIX.  TURNING OFF FOR NOW */
  return -1;

#ifdef DEBUG2
  char *gbuffer1, *gbuffer2;
  static int ncalls = 0;
#endif

#ifdef DEBUG2
  segment5_left = univdiagonal5 - querylength;
  segment3_left = univdiagonal3 - querylength;

  printf("Splice_resolve_fusion, plusp %d and %d, sense_forward_p %d, call %d: Getting genome at segment5_left %u [%u] and segment3_left %u [%u], range %d..%d, querylength %d\n",
	 plus5p,plus3p,sense_forward_p,++ncalls,
	 segment5_left,segment5_left - chroffset5,segment3_left,segment3_left - chroffset3,
	 querypos5,querypos3,querylength);

#endif

  assert(Compress_fwdp(query_compress_5) == plus5p);
  assert(Compress_fwdp(query_compress_3) == plus3p);


#if 0
  /* TRIM QUERY AT CHROMOSOME BOUNDS */
  univdiagonal5 = segment5_left + (Univcoord_T) querylength;
  univdiagonal3 = segment3_left + (Univcoord_T) querylength;
  if (plus5p == true) {
    pos5 = querypos5;
  } else {
    pos5 = (querylength - 1) - querypos5;
  }
  if (plus3p == true) {
    pos3 = querypos3;
  } else {
    pos3 = (querylength - 1) - querypos3;
  }
  pos5 = (univdiagonal5 + pos5 >= chroffset5 + querylength) ? pos5 : (int) (chroffset5 - segment5_left);
  pos3 = (univdiagonal3 + pos3 <= chrhigh3 + querylength) ? pos3 : (int) (chrhigh3 - segment3_left);
#endif

  /* Determine feasibility of a splice */
  mismatch_positions_left = spliceinfo->mismatch_positions_left1; /* Use allocated memory */

  /* Don't use nmismatches_allowed, because we want to count all mismatches for proper accounting and decision-making */
  if (plus5p == true) {
    debug2(printf("nmismatches_allowed is %d.  Calling Genomebits_mismatches_fromleft over %d..%d\n",
		  nmismatches_allowed,querypos5,querypos3));
    nmismatches_left = Genomebits_mismatches_fromleft(mismatch_positions_left,/*nmismatches_allowed*/querylength,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
						      query_compress_5,univdiagonal5,querylength,
						      querypos5,querypos3,/*plusp*/true,genestrand);
    mismatch5_adj = 0;				  /* splice occurs to right of querypos */
  } else {
    debug2(printf("nmismatches_allowed is %d.  Calling Genomebits_mismatches_fromleft over %d..%d\n",
		  nmismatches_allowed,querylength - querypos3,querylength - querypos5));
    nmismatches_left = Genomebits_mismatches_fromleft(mismatch_positions_left,/*nmismatches_allowed*/querylength,
						      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
						      query_compress_5,univdiagonal5,querylength,
						      querylength - querypos3,querylength - querypos5,
						      /*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_left,nmismatches_left,querylength);
    mismatch5_adj = 1;		/* splice occurs to left of querypos */
  }

  debug2(
	 printf("(4) %d mismatches on left from %d to %d at:",nmismatches_left,querypos5,querypos3);
	 for (int i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  mismatch_positions_right = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  /* Don't use nmismatches_allowed, because we want to count all mismatches for proper accounting and decision-making */
  if (plus3p == true) {
    debug2(printf("nmismatches_allowed is %d.  Calling Genomebits_mismatches_fromright over %d..%d\n",
		  nmismatches_allowed,querypos5,querypos3));
    nmismatches_right = Genomebits_mismatches_fromright(mismatch_positions_right,/*nmismatches_allowed*/querylength,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,
							query_compress_3,univdiagonal3,querylength,
							querypos5,querypos3,/*plusp*/true,genestrand);
    mismatch3_adj = 0;				  /* splice occurs to right of querypos */
  } else {
    debug2(printf("nmismatches_allowed is %d.  Calling Genomebits_mismatches_fromright over %d..%d\n",
		  nmismatches_allowed,querylength - querypos3,querylength - querypos5));
    nmismatches_right = Genomebits_mismatches_fromright(mismatch_positions_right,/*nmismatches_allowed*/querylength,
							/*ome*/genomebits,/*ome_alt*/genomebits_alt,
							query_compress_3,univdiagonal3,querylength,
							querylength - querypos3,querylength - querypos5,
							/*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_right,nmismatches_right,querylength);
    mismatch3_adj = 1;		/* splice occurs to left of querypos */
  }

  debug2(
	 printf("(4) %d mismatches on right from %d to %d at:",nmismatches_right,querypos3,querypos5);
	 for (int i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

#ifdef DEBUG2
  gbuffer1 = (char *) CALLOC(querylength+1,sizeof(char));
  gbuffer2 = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonal5 - querylength,querylength,gbuffer1);
  Genome_fill_buffer(univdiagonal3 - querylength,querylength,gbuffer2);
  if (plus3p != plus5p) {
    make_complement_inplace(gbuffer2,querylength);
  }
    
  char *queryseq = Compress_queryseq(query_compress_5,querylength);

  printf("Splice_resolve_fusion, plusp %d and %d, nmismatches_allowed %d: Getting genome at left %llu and %llu\n",
	 plus5p,plus3p,nmismatches_allowed,(unsigned long long) segment5_left,(unsigned long long) segment3_left);
  printf("g1: %s\n",gbuffer1);
  printf("    ");
  print_diffs(gbuffer1,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbuffer2,queryseq,querylength);
  printf("\n");
  printf("g2: %s\n",gbuffer2);
  FREE(queryseq);
  FREE(gbuffer2);
  FREE(gbuffer1);
#endif


  /* Find low bound for splice pos based on nmismatches_right */
  /* Since we invert the positions, need to set the value above at [nmismatches_allowed] */
  if (nmismatches_right < nmismatches_allowed) {
    splice_querypos_low = mismatch_positions_right[nmismatches_right] - 1; /* ? need - 1 */
  } else {
    splice_querypos_low = mismatch_positions_right[nmismatches_allowed] - 1; /* ? need - 1 */
  }
  debug2(printf("  Splice low bound: %d",splice_querypos_low));
  if (splice_querypos_low < querypos5 + 1) {
    splice_querypos_low = querypos5 + 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_low));
  }
  debug2(printf("\n"));

  /* Find high bound for splice pos based on nmismatches_left */
  if (nmismatches_left < nmismatches_allowed) {
    splice_querypos_high = mismatch_positions_left[nmismatches_left] + 1; /* ? need + 1 */;
  } else {
    splice_querypos_high = mismatch_positions_left[nmismatches_allowed] + 1; /* ? need + 1 */;
  }
  debug2(printf("  Splice high bound: %d",splice_querypos_high));
  if (splice_querypos_high > querypos3 - 1) {
    splice_querypos_high = querypos3 - 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_high));
  }
  debug2(printf("\n"));

  if (splice_querypos_low >= splice_querypos_high) {
    /* The two ends cannot meet */
    *best_nmismatches_5 = *best_nmismatches_3 = -1; /* Indicates that calling procedure needs to compute numbers of mismatches  */
    return -1;
  }


  segment5_left = univdiagonal5 - querylength;
  segment3_left = univdiagonal3 - querylength;

  best_splice_querypos = -1;
  if (sense_forward_p == true) {
    if (plus5p == /*sense_forward_p*/true) {
      /* (1) */
      segment5_nsites = compute_donor_sites(&optimal_donor_prob,&segment5_sites,&segment5_types,&segment5_probs,
					    spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					    spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					    spliceinfo->segmenti_probs_alloc,
					    /*splice_qpos_low*/splice_querypos_low,/*splice_qpos_high*/splice_querypos_high,querylength,
					    univdiagonal5,chroffset5,knownsplicing,intlistpool);
    } else {
      /* (2) */
      segment5_nsites = compute_antidonor_sites(&optimal_donor_prob,&segment5_sites,&segment5_types,&segment5_probs,
						spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						spliceinfo->segmenti_probs_alloc,
						/*splice_qpos_low*/querylength - splice_querypos_high,
						/*splice_qpos_high*/querylength - splice_querypos_low,querylength,
						univdiagonal5,chroffset5,knownsplicing,intlistpool);
      invert_sites(segment5_sites,segment5_types,segment5_probs,segment5_nsites,querylength);
    }

    if (plus3p == /*sense_forward_p*/true) {
      /* (3) */
      segment3_nsites = compute_acceptor_sites(&optimal_acceptor_prob,&segment3_sites,&segment3_types,&segment3_probs,
					       spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					       spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					       spliceinfo->segmentj_probs_alloc,
					       /*splice_qpos_low*/splice_querypos_low,/*splice_qpos_high*/splice_querypos_high,querylength,
					       univdiagonal3,chroffset3,knownsplicing,intlistpool);

    } else {
      /* (4) */
      segment3_nsites = compute_antiacceptor_sites(&optimal_acceptor_prob,&segment3_sites,&segment3_types,&segment3_probs,
						   spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
						   spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
						   spliceinfo->segmentj_sites_alloc,
						   /*splice_qpos_low*/querylength - splice_querypos_high,
						   /*splice_qpos_high*/querylength - splice_querypos_low,querylength,
						   univdiagonal3,chroffset3,knownsplicing,intlistpool);
      invert_sites(segment3_sites,segment3_types,segment3_probs,segment3_nsites,querylength);
    }

    best_splice_querypos =
      splice_sense(&(*donor1),&(*donor2),&(*acceptor1),&(*acceptor2),&(*best_donor_prob),&(*best_acceptor_prob),
		   /*D*/&(*best_nmismatches_5),/*A*/&(*best_nmismatches_3),
		   /*D*/&(*best_ref_nmismatches_5),/*A*/&(*best_ref_nmismatches_3),
		   /*D*/univdiagonal5,/*A*/univdiagonal3,/*D*/chroffset5,/*A*/chroffset3,
		   plus5p,plus3p,querylength,
		   /*donor*/mismatch_positions_left,/*donor*/nmismatches_left,
		   /*acceptor*/mismatch_positions_right,/*acceptor*/nmismatches_right,
		   /*D*/segment5_sites,/*A*/segment3_sites,
		   /*D*/segment5_types,/*A*/segment3_types,
		   /*D*/segment5_probs,/*A*/segment3_probs,
		   /*D*/segment5_nsites,/*A*/segment3_nsites);

  } else {
    if (plus5p == /*sense_forward_p*/false) {
      /* (5) */
      segment5_nsites = compute_acceptor_sites(&optimal_acceptor_prob,&segment5_sites,&segment5_types,&segment5_probs,
					       spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					       spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					       spliceinfo->segmenti_probs_alloc,
					       /*splice_qpos_low*/querylength - splice_querypos_high,
					       /*splice_qpos_high*/querylength - splice_querypos_low,querylength,
					       univdiagonal5,chroffset5,knownsplicing,intlistpool);
      invert_sites(segment5_sites,segment5_types,segment5_probs,segment5_nsites,querylength);

    } else {
      /* (6) */
      segment5_nsites = compute_antiacceptor_sites(&optimal_acceptor_prob,&segment5_sites,&segment5_types,&segment5_probs,
						   spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						   spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						   spliceinfo->segmenti_probs_alloc,
						   /*splice_qpos_low*/splice_querypos_low,/*splice_qpos_high*/splice_querypos_high,querylength,
						   univdiagonal5,chroffset5,knownsplicing,intlistpool);
    }

    if (plus3p == /*sense_forward_p*/false) {
      /* (7) */
      segment3_nsites = compute_donor_sites(&optimal_donor_prob,&segment3_sites,&segment3_types,&segment3_probs,
					    spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					    spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					    spliceinfo->segmentj_probs_alloc,
					    /*splice_qpos_low*/querylength - splice_querypos_high,
					    /*splice_qpos_high*/querylength - splice_querypos_low,querylength,
					    univdiagonal3,chroffset3,knownsplicing,intlistpool);
      invert_sites(segment3_sites,segment3_types,segment3_probs,segment3_nsites,querylength);

    } else {
      /* (8) */
      segment3_nsites = compute_antidonor_sites(&optimal_donor_prob,&segment3_sites,&segment3_types,&segment3_probs,
						spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
						spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
						spliceinfo->segmentj_probs_alloc,
						/*splice_qpos_low*/splice_querypos_low,/*splice_qpos_high*/splice_querypos_high,querylength,
						univdiagonal3,chroffset3,knownsplicing,intlistpool);
    }

    best_splice_querypos =
      splice_antisense(&(*donor1),&(*donor2),&(*acceptor1),&(*acceptor2),&(*best_donor_prob),&(*best_acceptor_prob),
		       /*A*/&(*best_nmismatches_5),/*D*/&(*best_nmismatches_3),
		       /*A*/&(*best_ref_nmismatches_5),/*D*/&(*best_ref_nmismatches_3),
		       /*A*/univdiagonal5,/*D*/univdiagonal3,/*A*/chroffset5,/*D*/chroffset3,
		       /*plusAp*/plus5p,/*plusDp*/plus3p,querylength,
		       /*acceptor*/mismatch_positions_left,/*acceptor*/nmismatches_left,
		       /*donor*/mismatch_positions_right,/*donor*/nmismatches_right,
		       /*A*/segment5_sites,/*D*/segment3_sites,
		       /*A*/segment5_types,/*D*/segment3_types,
		       /*A*/segment5_probs,/*D*/segment3_probs,
		       /*A*/segment5_nsites,/*D*/segment3_nsites);
  }


  if (best_splice_querypos >= 0) {
#if 0
    lefti = 0;
    while (lefti < nmismatches_left && mismatch_positions_left[lefti] < best_splice_qpos) {
      lefti++;
    }
    *best_nmismatches_i = *best_ref_nmismatches_i = lefti;
    
    righti = 0;
    while (righti < nmismatches_right && mismatch_positions_right[righti] >= best_splice_qpos) {
      righti++;
    }
    *best_nmismatches_j = *best_ref_nmismatches_j = righti;
#endif

    if (sense_forward_p == true) {
      best_prob5 = *best_donor_prob;
      best_prob3 = *best_acceptor_prob;
    } else {
      best_prob5 = *best_acceptor_prob;
      best_prob3 = *best_donor_prob;
    }

    support5 = best_splice_querypos - querypos5;
    support3 = querypos3 - best_splice_querypos;

    adj_support5 = support5 - 4*(*best_nmismatches_5);
    adj_support3 = support3 - 4*(*best_nmismatches_3);

    /* For fusions, we check probability before support */
    if (best_prob5 < SPLICE_PROB_LOW || best_prob3 < SPLICE_PROB_LOW) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on probs\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      /* Fall through */

    } else if (support5 >= 25 && support3 >= 25 && (*best_nmismatches_5) + (*best_nmismatches_3) <= 1) {
      debug2(printf("Accepting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on support\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return best_splice_querypos;

    } else if (sufficient_support_p(adj_support5,best_prob5) == false) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      /* Fall through */

    } else if (sufficient_support_p(adj_support3,best_prob3) == false) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      /* Fall through */

    } else {
      debug2(printf("Accepting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return best_splice_querypos;
    }
  }

  /* If we resort to atypical splice, then need to clear indel information */
  /* Not allowing atypical splices for gene fusions */

  *best_nindels = 0;
  *best_indel_pos = -1;
  return -1;
}
#endif


int
Splice_fusion_sense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		     double *best_donor_prob, double *best_acceptor_prob,

		     int *best_nmismatches_D, int *best_nmismatches_A,
		     int *best_ref_nmismatches_D, int *best_ref_nmismatches_A,
		       
		     Univcoord_T univdiagonalD, Univcoord_T univdiagonalA,
		     Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
		     Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
		     
		     int queryposD, int queryposA, int querylength,
		     Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     Intlistpool_T intlistpool, int genestrand) {

  int best_splice_querypos = -1;
  int splice_querypos_low, splice_querypos_high;
  int splice_queryposD_low, splice_queryposD_high,
    splice_queryposA_low, splice_queryposA_high;
  int splice_qposD_low, splice_qposD_high,
    splice_qposA_low, splice_qposA_high;

  int *mismatch_positions_donor, *mismatch_positions_acceptor;
  int nmismatches_donor, nmismatches_acceptor, nmismatches_to_trimpos;
  
  int *segmentD_sites, *segmentA_sites;
  int *segmentD_types, *segmentA_types;
  double *segmentD_probs, *segmentA_probs;
  int segmentD_nsites, segmentA_nsites;

  double optimal_donor_prob, optimal_acceptor_prob;
  int supportD, supportA;
 
#ifdef DEBUG2
  char *gbufferD, *gbufferA;
  static int ncalls = 0;
#endif

#ifdef DEBUG2
  printf("Splice_fusion_sense, plusp %d and %d, call %d: Getting genome at univdiagonalD %u and univdiagonalA %u, range %d..%d, querylength %d\n",
	 plusDp,plusAp,++ncalls,univdiagonalD,univdiagonalA,
	 queryposD,queryposA,querylength);
#endif

  assert(Compress_fwdp(query_compress_D) == plusDp);
  assert(Compress_fwdp(query_compress_A) == plusAp);

  mismatch_positions_donor = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_acceptor = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  if (plusDp == true) {
    nmismatches_donor =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_donor,/*max_mismatches*/queryposA - queryposD,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					      query_compress_D,univdiagonalD,querylength,
					      /*pos5*/queryposD,/*pos3*/queryposA,
					      /*plusp*/true,genestrand);
  } else {
    nmismatches_donor =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_donor,/*max_mismatches*/queryposA - queryposD,
					       /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					       query_compress_D,univdiagonalD,querylength,
					       /*pos5*/querylength - queryposA,/*pos3*/querylength - queryposD,
					       /*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_donor,nmismatches_donor,querylength);
  }

  if (plusAp == true) {
    nmismatches_acceptor =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_acceptor,/*max_mismatches*/queryposA - queryposD,
					       /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					       query_compress_A,univdiagonalA,querylength,
					       /*pos5*/queryposD,/*pos3*/queryposA,
					       /*plusp*/true,genestrand);
  } else {
    nmismatches_acceptor =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_acceptor,/*max_mismatches*/queryposA - queryposD,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					      query_compress_A,univdiagonalA,querylength,
					      /*pos5*/querylength - queryposA,/*pos3*/querylength - queryposD,
					      /*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_acceptor,nmismatches_acceptor,querylength);
  }

  debug2(
	 printf("(4) %d mismatches on donor from %d to %d at:",nmismatches_donor,queryposD,queryposA);
	 for (int i = 0; i <= nmismatches_donor; i++) {
	   printf(" %d",mismatch_positions_donor[i]);
	 }
	 printf("\n");
	 );

  debug2(
	 printf("(4) %d mismatches on acceptor from %d to %d at:",nmismatches_acceptor,queryposA,queryposD);
	 for (int i = 0; i <= nmismatches_acceptor; i++) {
	   printf(" %d",mismatch_positions_acceptor[i]);
	 }
	 printf("\n");
	 );

#ifdef DEBUG2
  gbufferD = (char *) CALLOC(querylength+1,sizeof(char));
  gbufferA = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonalD - querylength,querylength,gbufferD);
  Genome_fill_buffer(univdiagonalA - querylength,querylength,gbufferA);
  if (plusAp != plusDp) {
    make_complement_inplace(gbufferA,querylength);
  }
    
  char *queryseq = Compress_queryseq(query_compress_D,querylength);

  printf("Splice_fusion_sense, plusp %d and %d: Getting genome at %llu and %llu\n",
	 plusDp,plusAp,(unsigned long long) univdiagonalD,(unsigned long long) univdiagonalA);
  printf("gD: %s\n",gbufferD);
  printf("    ");
  print_diffs(gbufferD,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbufferA,queryseq,querylength);
  printf("\n");
  printf("gA: %s\n",gbufferA);
  FREE(queryseq);
  FREE(gbufferA);
  FREE(gbufferD);
#endif

  /* Compute as querypos, not qpos.  Using positions after inverting */
  splice_querypos_low =
    Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_acceptor,nmismatches_acceptor,
				      queryposD,queryposA);
  splice_querypos_high =
    Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_donor,nmismatches_donor,
				  queryposD,queryposA,querylength);

  debug2(printf("  Splice low bound: %d",splice_querypos_low));
  if (splice_querypos_low < queryposD + 1) {
    splice_querypos_low = queryposD + 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_low));
  }
  debug2(printf("\n"));

  debug2(printf("  Splice high bound: %d",splice_querypos_high));
  if (splice_querypos_high > queryposA - 1) {
    splice_querypos_high = queryposA - 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_high));
  }
  debug2(printf("\n"));

  if (splice_querypos_low >= splice_querypos_high + MISMATCHES_AT_SITE) {
    /* The two ends cannot meet */
    return -1;

  } else {
    splice_queryposD_low = subtract_bounded(splice_querypos_high,EXTRA_MATCHES,queryposD + 1);
    splice_queryposD_high = add_bounded(splice_querypos_high,EXTRA_MISMATCHES,queryposA - 1);
    splice_queryposA_low = subtract_bounded(splice_querypos_low,EXTRA_MISMATCHES,queryposD + 1);
    splice_queryposA_high = add_bounded(splice_querypos_low,EXTRA_MATCHES,queryposA - 1);
  }

  if (plusDp == true) {
    splice_qposD_low = splice_queryposD_low;
    splice_qposD_high = splice_queryposD_high;
  } else {
    splice_qposD_low = querylength - splice_queryposD_high;
    splice_qposD_high = querylength - splice_queryposD_low;
  }

  if (plusAp == true) {
    splice_qposA_low = splice_queryposA_low;
    splice_qposA_high = splice_queryposA_high;
  } else {
    splice_qposA_low = querylength - splice_queryposA_high;
    splice_qposA_high = querylength - splice_queryposA_low;
  }

  if (plusDp == /*sense_forward_p*/true) {
    debug2(printf("Computing donor sites from %d to %d\n",splice_qposD_low,splice_qposD_high));
    segmentD_nsites = compute_donor_sites(&optimal_donor_prob,&segmentD_sites,&segmentD_types,&segmentD_probs,
					  spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					  spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					  spliceinfo->segmenti_probs_alloc,splice_qposD_low,splice_qposD_high,
					  querylength,univdiagonalD,chroffset_D,knownsplicing,intlistpool);

  } else {
    debug2(printf("Computing antidonor sites from %d to %d\n",splice_qposD_low,splice_qposD_high));
    segmentD_nsites = compute_antidonor_sites(&optimal_donor_prob,&segmentD_sites,&segmentD_types,&segmentD_probs,
					      spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					      spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					      spliceinfo->segmenti_probs_alloc,splice_qposD_low,splice_qposD_high,
					      querylength,univdiagonalD,chroffset_D,knownsplicing,intlistpool);
  }
  
  if (plusAp == /*sense_forward_p*/true) {
    debug2(printf("Computing acceptor sites from %d to %d\n",splice_qposA_low,splice_qposA_high));
    segmentA_nsites = compute_acceptor_sites(&optimal_acceptor_prob,&segmentA_sites,&segmentA_types,&segmentA_probs,
					     spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					     spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					     spliceinfo->segmentj_probs_alloc,splice_qposA_low,splice_qposA_high,
					     querylength,univdiagonalA,chroffset_A,knownsplicing,intlistpool);
  } else {
    debug2(printf("Computing antiacceptor sites from %d to %d\n",splice_qposA_low,splice_qposA_high));
    segmentA_nsites = compute_antiacceptor_sites(&optimal_acceptor_prob,&segmentA_sites,&segmentA_types,&segmentA_probs,
						 spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
						 spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
						 spliceinfo->segmentj_probs_alloc,splice_qposA_low,splice_qposA_high,
						 querylength,univdiagonalA,chroffset_A,knownsplicing,intlistpool);
  }

  if (plusDp == false) {
    invert_sites(segmentD_sites,segmentD_types,segmentD_probs,segmentD_nsites,querylength);
  }
  if (plusAp == false) {
    invert_sites(segmentA_sites,segmentA_types,segmentA_probs,segmentA_nsites,querylength);
  }
  
  best_splice_querypos =
    splice_sense(&(*best_type),&(*donor1),&(*donor2),&(*acceptor1),&(*acceptor2),&(*best_donor_prob),&(*best_acceptor_prob),
		 /*D*/&(*best_nmismatches_D),/*A*/&(*best_nmismatches_A),
		 /*D*/&(*best_ref_nmismatches_D),/*A*/&(*best_ref_nmismatches_A),
		 /*D*/univdiagonalD,/*A*/univdiagonalA,/*D*/chroffset_D,/*A*/chroffset_A,
		 plusDp,plusAp,querylength,
		 /*donor*/mismatch_positions_donor,/*donor*/nmismatches_donor,
		 /*acceptor*/mismatch_positions_acceptor,/*acceptor*/nmismatches_acceptor,
		 /*D*/segmentD_sites,/*A*/segmentA_sites,
		 /*D*/segmentD_types,/*A*/segmentA_types,
		 /*D*/segmentD_probs,/*A*/segmentA_probs,
		 /*D*/segmentD_nsites,/*A*/segmentA_nsites);

  if (best_splice_querypos < 0) {
    return -1;
  } else {
    supportD = best_splice_querypos - queryposD;
    supportA = queryposA - best_splice_querypos;

    /* adj_supportD = supportD - 4*(*best_nmismatches_D); */
    /* adj_supportA = supportA - 4*(*best_nmismatches_A); */

    /* For fusions, we check probability before support */
    if ((*best_donor_prob) < SPLICE_PROB_LOW || (*best_acceptor_prob) < SPLICE_PROB_LOW) {
      debug2(printf("Rejecting %c%c-%c%c splice pos %d, with supportD %d, mismatches %d and supportA %d, mismatches %d, probs %f and %f, based on probs\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportD,*best_nmismatches_D,
		    supportA,*best_nmismatches_A,*best_donor_prob,*best_acceptor_prob));
      return -1;

    } else if (supportD >= 25 && supportA >= 25 && (*best_nmismatches_D) + (*best_nmismatches_A) <= 1) {
      debug2(printf("Accepting %c%c-%c%c splice pos %d, with supportD %d, mismatches %d and supportA %d, mismatches %d, probs %f and %f, based on support\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportD,*best_nmismatches_D,
		    supportA,*best_nmismatches_A,*best_donor_prob,*best_acceptor_prob));
      return best_splice_querypos;

#if 0
    } else if (sufficient_support_p(adj_support5,best_prob5) == false) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return -1;

    } else if (sufficient_support_p(adj_support3,best_prob3) == false) {
      debug2(printf("Rejecting %c%c-%c%c splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return -1;
#endif

    } else {
      debug2(printf("Accepting %c%c-%c%c splice pos %d, with supportD %d, mismatches %d and supportA %d, mismatches %d, probs %f and %f\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportD,*best_nmismatches_D,
		    supportA,*best_nmismatches_A,*best_donor_prob,*best_acceptor_prob));
      return best_splice_querypos;
    }
  }
}


int
Splice_fusion_antisense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
			 double *best_donor_prob, double *best_acceptor_prob,

			 int *best_nmismatches_A, int *best_nmismatches_D,
			 int *best_ref_nmismatches_A, int *best_ref_nmismatches_D,
		       
			 Univcoord_T univdiagonalA, Univcoord_T univdiagonalD,
			 Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
			 Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
			 
			 int queryposA, int queryposD, int querylength,
			 Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
			 Intlistpool_T intlistpool, int genestrand) {

  int best_splice_querypos = -1;
  int splice_querypos_low, splice_querypos_high;
  int splice_queryposA_low, splice_queryposA_high,
    splice_queryposD_low, splice_queryposD_high;
  int splice_qposA_low, splice_qposA_high,
    splice_qposD_low, splice_qposD_high;

  int *mismatch_positions_acceptor, *mismatch_positions_donor;
  int nmismatches_acceptor, nmismatches_donor, nmismatches_to_trimpos;
  
  int *segmentA_sites, *segmentD_sites;
  int *segmentA_types, *segmentD_types;
  double *segmentA_probs, *segmentD_probs;
  int segmentA_nsites, segmentD_nsites;

  double optimal_donor_prob, optimal_acceptor_prob;
  int supportA, supportD;
 
#ifdef DEBUG2
  char *gbufferA, *gbufferD;
  static int ncalls = 0;
#endif

#ifdef DEBUG2
  printf("Splice_fusion_antisense, plusp %d and %d, call %d: Getting genome at univdiagonalA %u and univdiagonalD %u, range %d..%d, querylength %d\n",
	 plusAp,plusDp,++ncalls,univdiagonalA,univdiagonalD,
	 queryposA,queryposD,querylength);
#endif

  assert(Compress_fwdp(query_compress_A) == plusAp);
  assert(Compress_fwdp(query_compress_D) == plusDp);

  mismatch_positions_acceptor = spliceinfo->mismatch_positions_left1; /* Use allocated memory */
  mismatch_positions_donor = spliceinfo->mismatch_positions_right2; /* Use allocated memory */

  if (plusAp == true) {
    nmismatches_acceptor =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_acceptor,/*max_mismatches*/queryposD - queryposA,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					      query_compress_A,univdiagonalA,querylength,
					      /*pos5*/queryposA,/*pos3*/queryposD,
					      /*plusp*/true,genestrand);
  } else {
    nmismatches_acceptor =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_acceptor,/*max_mismatches*/queryposD - queryposA,
					       /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					       query_compress_A,univdiagonalA,querylength,
					       /*pos5*/querylength - queryposD,/*pos3*/querylength - queryposA,
					       /*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_acceptor,nmismatches_acceptor,querylength);
  }

  if (plusDp == true) {
    nmismatches_donor =
      Genomebits_mismatches_fromright_for_trim(mismatch_positions_donor,/*max_mismatches*/queryposD - queryposA,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					      query_compress_D,univdiagonalD,querylength,
					      /*pos5*/queryposA,/*pos3*/queryposD,
					      /*plusp*/true,genestrand);
  } else {
    nmismatches_donor =
      Genomebits_mismatches_fromleft_for_trim(mismatch_positions_donor,/*max_mismatches*/queryposD - queryposA,
					      /*ome*/genomebits,/*ome_alt*/genomebits_alt,
					      query_compress_D,univdiagonalD,querylength,
					      /*pos5*/querylength - queryposD,/*pos3*/querylength - queryposA,
					      /*plusp*/false,genestrand);
    invert_mismatch_positions(mismatch_positions_donor,nmismatches_donor,querylength);
  }

  debug2(
	 printf("(4) %d mismatches on acceptor from %d to %d at:",nmismatches_acceptor,queryposD,queryposA);
	 for (int i = 0; i <= nmismatches_acceptor; i++) {
	   printf(" %d",mismatch_positions_acceptor[i]);
	 }
	 printf("\n");
	 );

  debug2(
	 printf("(4) %d mismatches on donor from %d to %d at:",nmismatches_donor,queryposA,queryposD);
	 for (int i = 0; i <= nmismatches_donor; i++) {
	   printf(" %d",mismatch_positions_donor[i]);
	 }
	 printf("\n");
	 );

#ifdef DEBUG2
  gbufferA = (char *) CALLOC(querylength+1,sizeof(char));
  gbufferD = (char *) CALLOC(querylength+1,sizeof(char));
  Genome_fill_buffer(univdiagonalA - querylength,querylength,gbufferA);
  Genome_fill_buffer(univdiagonalD - querylength,querylength,gbufferD);
  if (plusDp != plusAp) {
    make_complement_inplace(gbufferD,querylength);
  }
    
  char *queryseq = Compress_queryseq(query_compress_A,querylength);

  printf("Splice_fusion_antisense, plusp %d and %d: Getting genome at left %llu and %llu\n",
	 plusAp,plusDp,(unsigned long long) univdiagonalA,(unsigned long long) univdiagonalD);
  printf("gA: %s\n",gbufferA);
  printf("    ");
  print_diffs(gbufferA,queryseq,querylength);
  printf("\n");
  printf("q:  %s\n",queryseq);
  printf("    ");
  print_diffs(gbufferD,queryseq,querylength);
  printf("\n");
  printf("gD: %s\n",gbufferD);
  FREE(queryseq);
  FREE(gbufferD);
  FREE(gbufferA);
#endif

  /* Compute as querypos, not qpos.  Using positions after inverting */
  splice_querypos_low =
    Spliceends_trim_qstart_nosplice(&nmismatches_to_trimpos,mismatch_positions_donor,nmismatches_donor,
				      queryposA,queryposD);
  splice_querypos_high =
    Spliceends_trim_qend_nosplice(&nmismatches_to_trimpos,mismatch_positions_acceptor,nmismatches_acceptor,
				  queryposA,queryposD,querylength);

  debug2(printf("  Splice low bound: %d",splice_querypos_low));
  if (splice_querypos_low < queryposA + 1) {
    splice_querypos_low = queryposA + 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_low));
  }
  debug2(printf("\n"));

  debug2(printf("  Splice high bound: %d",splice_querypos_high));
  if (splice_querypos_high > queryposD - 1) {
    splice_querypos_high = queryposD - 1;
    debug2(printf(" -- changed to %d\n",splice_querypos_high));
  }
  debug2(printf("\n"));

  if (splice_querypos_low >= splice_querypos_high + MISMATCHES_AT_SITE) {
    /* The two ends cannot meet */
    return -1;

  } else {
    splice_queryposA_low = subtract_bounded(splice_querypos_low,EXTRA_MISMATCHES,queryposA + 1);
    splice_queryposA_high = add_bounded(splice_querypos_low,EXTRA_MATCHES,queryposD - 1);
    splice_queryposD_low = subtract_bounded(splice_querypos_high,EXTRA_MATCHES,queryposA + 1);
    splice_queryposD_high = add_bounded(splice_querypos_high,EXTRA_MISMATCHES,queryposD - 1);
  }

  if (plusAp == true) {
    splice_qposA_low = splice_queryposA_low;
    splice_qposA_high = splice_queryposA_high;
  } else {
    splice_qposA_low = querylength - splice_queryposA_high;
    splice_qposA_high = querylength - splice_queryposA_low;
  }

  if (plusDp == true) {
    splice_qposD_low = splice_queryposD_low;
    splice_qposD_high = splice_queryposD_high;
  } else {
    splice_qposD_low = querylength - splice_queryposD_high;
    splice_qposD_high = querylength - splice_queryposD_low;
  }

  if (plusAp == /*sense_forward_p*/false) {
    debug2(printf("Computing acceptor sites from %d to %d\n",splice_qposA_low,splice_qposA_high));
    segmentA_nsites = compute_acceptor_sites(&optimal_acceptor_prob,&segmentA_sites,&segmentA_types,&segmentA_probs,
					     spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
					     spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
					     spliceinfo->segmenti_probs_alloc,splice_qposA_low,splice_qposA_high,
					     querylength,univdiagonalA,chroffset_A,knownsplicing,intlistpool);
  } else {
    debug2(printf("Computing antiacceptor sites from %d to %d\n",splice_qposA_low,splice_qposA_high));
    segmentA_nsites = compute_antiacceptor_sites(&optimal_acceptor_prob,&segmentA_sites,&segmentA_types,&segmentA_probs,
						 spliceinfo->segmenti_sites_alloc1,spliceinfo->segmenti_types_alloc1,
						 spliceinfo->segmenti_sites_alloc2,spliceinfo->segmenti_types_alloc2,
						 spliceinfo->segmenti_probs_alloc,splice_qposA_low,splice_qposA_high,
						 querylength,univdiagonalA,chroffset_A,knownsplicing,intlistpool);
  }

  if (plusDp == /*sense_forward_p*/false) {
    debug2(printf("Computing donor sites from %d to %d\n",splice_qposD_low,splice_qposD_high));
    segmentD_nsites = compute_donor_sites(&optimal_donor_prob,&segmentD_sites,&segmentD_types,&segmentD_probs,
					  spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					  spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					  spliceinfo->segmentj_probs_alloc,splice_qposD_low,splice_qposD_high,
					  querylength,univdiagonalD,chroffset_D,knownsplicing,intlistpool);

  } else {
    debug2(printf("Computing antidonor sites from %d to %d\n",splice_qposD_low,splice_qposD_high));
    segmentD_nsites = compute_antidonor_sites(&optimal_donor_prob,&segmentD_sites,&segmentD_types,&segmentD_probs,
					      spliceinfo->segmentj_sites_alloc1,spliceinfo->segmentj_types_alloc1,
					      spliceinfo->segmentj_sites_alloc2,spliceinfo->segmentj_types_alloc2,
					      spliceinfo->segmentj_probs_alloc,splice_qposD_low,splice_qposD_high,
					      querylength,univdiagonalD,chroffset_D,knownsplicing,intlistpool);
  }
  
  if (plusAp == false) {
    invert_sites(segmentA_sites,segmentA_types,segmentA_probs,segmentA_nsites,querylength);
  }
  if (plusDp == false) {
    invert_sites(segmentD_sites,segmentD_types,segmentD_probs,segmentD_nsites,querylength);
  }
  
  best_splice_querypos =
    splice_antisense(&(*best_type),&(*donor1),&(*donor2),&(*acceptor1),&(*acceptor2),&(*best_donor_prob),&(*best_acceptor_prob),
		     /*A*/&(*best_nmismatches_A),/*D*/&(*best_nmismatches_D),
		     /*A*/&(*best_ref_nmismatches_A),/*D*/&(*best_ref_nmismatches_D),
		     /*A*/univdiagonalA,/*D*/univdiagonalD,/*A*/chroffset_A,/*D*/chroffset_D,
		     plusAp,plusDp,querylength,
		     /*acceptor*/mismatch_positions_acceptor,/*acceptor*/nmismatches_acceptor,
		     /*donor*/mismatch_positions_donor,/*donor*/nmismatches_donor,
		     /*A*/segmentA_sites,/*D*/segmentD_sites,
		     /*A*/segmentA_types,/*D*/segmentD_types,
		     /*A*/segmentA_probs,/*D*/segmentD_probs,
		     /*A*/segmentA_nsites,/*D*/segmentD_nsites);

  if (best_splice_querypos < 0) {
    return -1;
  } else {
    supportA = best_splice_querypos - queryposA;
    supportD = queryposD - best_splice_querypos;

    /* adj_supportA = supportA - 4*(*best_nmismatches_A); */
    /* adj_supportD = supportD - 4*(*best_nmismatches_D); */

    /* For fusions, we check probability before support */
    if ((*best_acceptor_prob) < SPLICE_PROB_LOW || (*best_donor_prob) < SPLICE_PROB_LOW) {
      debug2(printf("Rejecting %c%c-%c%c splice pos %d, with supportA %d, mismatches %d and supportD %d, mismatches %d, probs %f and %f, based on probs\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportA,*best_nmismatches_A,
		    supportD,*best_nmismatches_D,*best_donor_prob,*best_acceptor_prob));
      return -1;

    } else if (supportA >= 25 && supportD >= 25 && (*best_nmismatches_A) + (*best_nmismatches_D) <= 1) {
      debug2(printf("Accepting %c%c-%c%c splice pos %d, with supportA %d, mismatches %d and supportD %d, mismatches %d, probs %f and %f, based on support\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportA,*best_nmismatches_A,
		    supportD,*best_nmismatches_D,*best_donor_prob,*best_acceptor_prob));
      return best_splice_querypos;

#if 0
    } else if (sufficient_support_p(adj_support5,best_prob5) == false) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return -1;

    } else if (sufficient_support_p(adj_support3,best_prob3) == false) {
      debug2(printf("Rejecting splice pos %d, nindels %d, indel pos %d with support5 %d, mismatches %d, prob5 %f and support3 %d, mismatches %d, prob3 %f, based on adj support and prob\n",
		    best_splice_querypos,*best_nindels,*best_indel_pos,support5,*best_nmismatches_5,best_prob5,
		    support3,*best_nmismatches_3,best_prob3));
      return -1;
#endif

    } else {
      debug2(printf("Accepting %c%c-%c%c splice pos %d, with supportA %d, mismatches %d and supportD %d, mismatches %d, probs %f and %f\n",
		    *donor1,*donor2,*acceptor2,*acceptor1,
		    best_splice_querypos,supportA,*best_nmismatches_A,
		    supportD,*best_nmismatches_D,*best_donor_prob,*best_acceptor_prob));
      return best_splice_querypos;
    }
  }
}


