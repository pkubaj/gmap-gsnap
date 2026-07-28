#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "indel.h"

#include "assert.h"
#include "mem.h"
#include "genome.h"		/* For Genome_fill_buffer */
#include "genomebits_count.h"
#include "genomebits_mismatches.h"
#include "intron.h"


/* Causes problems with counting mismatches */
/* #define TRIM_AT_CHROMOSOME_BOUNDS 1 */


/* Resolve indels */ 
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Solve end indels */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif


static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;

static int min_indel_end_matches;
static bool maskedp = false;

void
Indel_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
	     int min_indel_end_matches_in, bool maskedp_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  min_indel_end_matches = min_indel_end_matches_in;
  maskedp = maskedp_in;

  return;
}


/* Use keep memory pool because we retain this memory if the querylength remains the same */
void
Indelinfo_free (Indelinfo_T *old) {
  FREE((*old)->int_memory);
  FREE(*old);
  return;
}

/* Use keep memory pool because we retain this memory if the querylength remains the same */
Indelinfo_T
Indelinfo_new (int querylength) {
  Indelinfo_T new = (Indelinfo_T) MALLOC(sizeof(*new));
  int a, b = (querylength+MISMATCH_EXTRA);

  new->int_memory = (int *) MALLOC((2*b) * sizeof(int));

  /* 2 arrays of (querylength+MISMATCH_EXTRA) ints */
  new->mismatch_positions_left = &(new->int_memory[0]); a = b;
  new->mismatch_positions_right = &(new->int_memory[a]); /*a += b;*/

  return new;
}


/* For transcriptome alignments, plusp (for tplusp) may be true or false; want_lowest_coordinate_p is based on gplus */
/* For alignments via middle_path with ascending univdiags, plusp should be true */
/* indels is positive here */
/* Ideally, we should optimize based on ref_nmismatches rather than nmismatches, but this seems difficult */

/* TODO: Make use of knownindels */
int
Indel_resolve_middle_insertion (int *best_nmismatches_i, int *best_nmismatches_j,
				int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
				Univcoord_T univdiagonal_i, int indels, Univcoord_T chrhigh,
				int *mismatch_positions_left, int nmismatches_left,
				int *mismatch_positions_right, int nmismatches_right,
				Genomebits_T ome, Genomebits_T ome_alt, Compress_T query_compress,
				int pos5, int pos3, int querylength, Indelinfo_T indelinfo,
				bool plusp, int genestrand, bool want_lowest_coordinate_p) {

  int best_indel_pos, best_indel_pos_case1, best_indel_pos_case2,
    indel_qpos_low, indel_qpos_high, indel_pos;
#ifdef DEBUG2
  Univcoord_T left;
  char *gbuffer, *queryseq;
  int i;
#endif
  int nmismatches_left_case1, nmismatches_right_case1,
    nmismatches_left_case2, nmismatches_right_case2;
  int best_sum, sum, lefti, righti;
  Univcoord_T univdiagonal_j;
  int max_mismatches_allowed = querylength;

  debug2(printf("Entered Indel_resolve_middle_insertion at univdiagonal_i %u, indels %d, with pos5 %u and pos3 %u.  Want lowest: %d\n",
		univdiagonal_i,indels,pos5,pos3,want_lowest_coordinate_p));

  assert(indels > 0);		/* Insertions are positive */

  univdiagonal_j = univdiagonal_i - indels;

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonali + pos5 >= chroffset + (Univcoord_T) querylength) ? pos5 : (int) (chroffset - left);
  pos3 = (univdiagonalj + pos3 <= chrhigh + (Univcoord_T) querylength) ? pos3 : (int) (chrhigh - (left - indels));
#endif

  if (pos5 >= pos3) {
    return -1;
  } else if (univdiagonal_i >= chrhigh) {
    /* Don't attempt to align past chromosome bounds */
    return -1;
  } else if (mismatch_positions_left != NULL) {
    /* ? Does not work for spliceindels */
    /* Re-use the results from the caller (end indels or spliceindel).  Have checked, and this is okay. */
    debug2(printf("Re-using mismatch_positions_left from the caller\n"));
  } else {
    mismatch_positions_left = indelinfo->mismatch_positions_left; /* Use allocated memory */
    debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromleft over %d..%d\n",
		  max_mismatches_allowed,pos5,pos3));

    nmismatches_left = Genomebits_mismatches_fromleft(mismatch_positions_left,max_mismatches_allowed,
						  ome,ome_alt,query_compress,univdiagonal_i,querylength,
						  pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_left <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
  }
  debug2(
	 printf("%d mismatches on left in %d..%d at:",nmismatches_left,pos5,pos3);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  if (mismatch_positions_right != NULL) {
    /* ? Does not work for spliceindels */
    /* Re-use the results from the caller (end indels or spliceindel).  Have checked, and this is okay. */
    debug2(printf("Re-using mismatch_positions_right from the caller\n"));
  } else {
    mismatch_positions_right = indelinfo->mismatch_positions_right; /* Use allocated memory */
    debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromright over %d..%d\n",
		  max_mismatches_allowed,pos5,pos3));

    nmismatches_right = Genomebits_mismatches_fromright(mismatch_positions_right,max_mismatches_allowed,
						    ome,ome_alt,query_compress,univdiagonal_j,querylength,
						    pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_right <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
  }
  debug2(
	 printf("%d mismatches on right in %d..%d at:",nmismatches_right,pos5,pos3);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  /* 1. Find low bound for indel pos based on nmismatches_right */
  if (nmismatches_right < max_mismatches_allowed) {
    indel_qpos_low = mismatch_positions_right[nmismatches_right] - indels + 1;
  } else {
    indel_qpos_low = mismatch_positions_right[max_mismatches_allowed] - indels + 1;
    debug2(printf("indel_qpos_low is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_low));
  }

  debug2(printf("  Comparing indel low bound %d with pos5 %d + 1\n",indel_qpos_low,pos5));
  if (indel_qpos_low < pos5 + 1) {
    /* Require separation from previous segment */
    indel_qpos_low = pos5 + 1;
  }
  if (indel_qpos_low < min_indel_end_matches) {
    /* Require separation from start of query */
    indel_qpos_low = min_indel_end_matches;
  }
  debug2(printf("  Indel low bound: %d\n",indel_qpos_low));


  /* 2. Find high bound for indel pos based on nmismatches_left */
  if (nmismatches_left < max_mismatches_allowed) {
    indel_qpos_high = mismatch_positions_left[nmismatches_left];
  } else {
    indel_qpos_high = mismatch_positions_left[max_mismatches_allowed];
    debug2(printf("indel_qpos_high is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_high));
  }

  debug2(printf("  Comparing indel high bound %d with pos3 %d - 1 - nindels %d\n",
		indel_qpos_high,pos3,indels));
  if (indel_qpos_high > (pos3 - 1) - indels) {
    /* Require separation from next segment */
    indel_qpos_high = (pos3 - 1) - indels;
  }
  if (indel_qpos_high > (querylength - min_indel_end_matches) - indels) {
    /* Require separation from end of query */
    indel_qpos_high = (querylength - min_indel_end_matches) - indels;
  }
  debug2(printf("  Indel high bound: %d\n",indel_qpos_high));

  if (indel_qpos_low > indel_qpos_high) {
    debug2(printf("Indel_resolve_middle_insertion returning -1 before search\n"));
    return -1;
  }


  /* query has insertion.  Get |indels| less from genome; trim from left. */
  /* left = ptr->diagonal - querylength; */

#ifdef DEBUG2
  gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
  left = univdiagonal_i - querylength;
  Genome_fill_buffer(left+indels,querylength-indels,gbuffer);
  queryseq = Compress_queryseq(query_compress,querylength);

  printf("solve_middle_indel, plusp %d, want low %d, insertion: Getting genome at left %llu + indels %d = %llu\n",
	 plusp,want_lowest_coordinate_p,(unsigned long long) left,indels,(unsigned long long) left+indels);
  printf("g1: %s\n",gbuffer);
  printf("q:  %s\n",queryseq);
  printf("g2: %s\n",&(gbuffer[indels]));
  FREE(queryseq);
  FREE(gbuffer);
#endif


  if (want_lowest_coordinate_p == true) {
    best_indel_pos_case1 = best_indel_pos_case2 = querylength;
  } else {
    best_indel_pos_case1 = best_indel_pos_case2 = -1;
  }

  /* Case 1: Use nmismatches_right as decision points */
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] > indel_qpos_high + indels - 1) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    /* Adjust because indel_pos is mismatch_positions_right - indels + 1 */
    righti++;
  }

  lefti = nmismatches_left - 1;
  while (lefti >= 0 && mismatch_positions_left[lefti] > indel_qpos_high) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    lefti--;
  }

  /* Search in descending qpos */
  best_sum = querylength;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= indel_qpos_low + indels - 1) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti] - indels) {
      lefti--;
    }
    sum = righti + lefti + 1;
  
    debug2(printf("  (Case 1) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]-indels+1));
    if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == true)) {
      indel_pos = mismatch_positions_right[righti] - indels + 1;
      assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);
      best_indel_pos_case1 = indel_pos;
      nmismatches_right_case1 = righti;
      nmismatches_left_case1 = lefti + 1;
      debug2(printf("**"));
      best_sum = sum;
    }
    righti++;
  }
  debug2(printf("\n"));


  /* Case 2: Use nmismatches_left as decision points */
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < indel_qpos_low) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    lefti++;
  }

  righti = nmismatches_right - 1;
  while (righti >= 0 && mismatch_positions_right[righti] < indel_qpos_low + indels - 1) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    /* Adjust because indel_pos is mismatch_positions_right - indels + 1 */
    righti--;
  }

  /* Search in ascending qpos */
  best_sum = querylength;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] <= indel_qpos_high) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti] + indels) {
      righti--;
    }
    sum = lefti + righti + 1;

    debug2(printf("  (Case 2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == false)) {
      indel_pos = mismatch_positions_left[lefti];
      best_indel_pos_case2 = indel_pos;
      nmismatches_right_case2 = righti + 1;
      nmismatches_left_case2 = lefti;
      debug2(printf("**"));
      best_sum = sum;
    }
    lefti++;
  }
  debug2(printf("\n"));

  if (want_lowest_coordinate_p == true) {
    if (best_indel_pos_case1 <= best_indel_pos_case2) {
      debug2(printf("%d <= %d, so case1 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
    } else {
      debug2(printf("%d > %d, so case2 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
    }
  } else {
    if (best_indel_pos_case1 >= best_indel_pos_case2) {
      debug2(printf("%d >= %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
    } else {
      debug2(printf("%d < %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
    }
  }

  if (best_indel_pos == -1 || best_indel_pos == querylength) {
    debug2(printf("Indel_resolve_middle_insertion returning -1 after search\n"));
    return -1;
  }

  if (maskedp == false) {
    *best_nmismatches_i = *best_ref_nmismatches_i;
    *best_nmismatches_j = *best_ref_nmismatches_j;

  } else {
    /* Need to get masked information */
    *best_nmismatches_i =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_i),genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,
					    pos5,/*pos3*/best_indel_pos,plusp,genestrand);
    *best_nmismatches_j =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_j),genomebits,genomebits_alt,query_compress,
					    univdiagonal_j,querylength,
					    /*pos5*/best_indel_pos+indels,pos3,plusp,genestrand);
  }

  debug2(printf("Indel_resolve_middle_insertion returning %d with mismatches %d+%d\n\n",
		best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
  return best_indel_pos;
}


/* For transcriptome alignments, plusp may be true or false */
/* For alignments via middle_path with ascending univdiags, plusp should be true */
/* indels is negative here */
/* Ideally, we should optimize based on ref_nmismatches rather than nmismatches, but this seems difficult */

/* Caller can provide either mismatch_positions_left or
   mismatch_positions_right to save on computation */

/* TODO: Make use of knownindels */
int
Indel_resolve_middle_deletion (int *best_nmismatches_i, int *best_nmismatches_j,
			       int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
			       Univcoord_T univdiagonal_i, int indels, Univcoord_T chrhigh,
			       int *mismatch_positions_left, int nmismatches_left,
			       int *mismatch_positions_right, int nmismatches_right,
			       Genomebits_T ome, Genomebits_T ome_alt, Compress_T query_compress,
			       int pos5, int pos3, int querylength, Indelinfo_T indelinfo,
			       bool plusp, int genestrand, bool want_lowest_coordinate_p) {
  int best_indel_pos, best_indel_pos_case1, best_indel_pos_case2,
    indel_qpos_low, indel_qpos_high, indel_pos;
#ifdef DEBUG2
  Univcoord_T left;
  char *gbuffer, *queryseq;
  int i;
#endif
  int nmismatches_left_case1, nmismatches_right_case1,
    nmismatches_left_case2, nmismatches_right_case2;
  int best_sum, sum, lefti, righti;
  Univcoord_T univdiagonal_j;
  int max_mismatches_allowed = querylength;


  univdiagonal_j = univdiagonal_i - indels;

  debug2(printf("Entered Indel_resolve_middle_deletion at univdiagonal %u, indels %d, with pos5 %u and pos3 %u.  Want lowest: %d\n",
		univdiagonal_i,indels,pos5,pos3,want_lowest_coordinate_p));

  assert(indels < 0);		/* Deletions are negative */

#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  pos5 = (univdiagonali + pos5 >= chroffset + (Univcoord_T) querylength) ? pos5 : (int) (chroffset - left);
  pos3 = (univdiagonalj + pos3 <= chrhigh + (Univcoord_T) querylength) ? pos3 : (int) (chrhigh - (left - indels));
#endif
    
  if (pos5 >= pos3) {
    return -1;
  } else if (univdiagonal_j >= chrhigh) {
    /* Don't attempt to align past chromosome bounds */
    return -1;
  } else if (mismatch_positions_left != NULL) {
    /* ? Does not work for spliceindels */
    /* Re-use the results from the caller (end indels or spliceindel).  Have checked, and this is okay. */
    debug2(printf("Reusing mismatch_positions_left\n"));
  } else {
    mismatch_positions_left = indelinfo->mismatch_positions_left; /* Use allocated memory */
    debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromleft over %d..%d\n",
		  max_mismatches_allowed,pos5,pos3));

    nmismatches_left = Genomebits_mismatches_fromleft(mismatch_positions_left,max_mismatches_allowed,
						  ome,ome_alt,query_compress,univdiagonal_i,querylength,
						  pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_left <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
  }
  debug2(
	 printf("%d mismatches on left in %d..%d at:",nmismatches_left,pos5,pos3);
	 for (i = 0; i <= nmismatches_left; i++) {
	   printf(" %d",mismatch_positions_left[i]);
	 }
	 printf("\n");
	 );

  if (mismatch_positions_right != NULL) {
    /* ? Does not work for spliceindels */
    /* Re-use the results from the caller (end indels or spliceindel).  Have checked, and this is okay. */
    debug2(printf("Reusing mismatch_positions_right\n"));
  } else {
    mismatch_positions_right = indelinfo->mismatch_positions_right; /* Use allocated memory */
    debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromright over %d..%d\n",
		  max_mismatches_allowed,pos5,pos3));

    nmismatches_right = Genomebits_mismatches_fromright(mismatch_positions_right,max_mismatches_allowed,
						    ome,ome_alt,query_compress,univdiagonal_j,querylength,
						    pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_right <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
  }
  debug2(
	 printf("%d mismatches on right in %d..%d at:",nmismatches_right,pos5,pos3);
	 for (i = 0; i <= nmismatches_right; i++) {
	   printf(" %d",mismatch_positions_right[i]);
	 }
	 printf("\n");
	 );

  /* 1. Find low bound for indel pos based on nmismatches_right */
  if (nmismatches_right < max_mismatches_allowed) {
    indel_qpos_low = mismatch_positions_right[nmismatches_right];
  } else {
    indel_qpos_low = mismatch_positions_right[max_mismatches_allowed];
    debug2(printf("indel_qpos_low is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_low));
  }

  debug2(printf("  Comparing indel low bound %d with pos5 %d + 1\n",indel_qpos_low,pos5));
  if (indel_qpos_low < pos5 + 1) {
    /* Require separation from previous segment */
    indel_qpos_low = pos5 + 1;
  }
  if (indel_qpos_low < min_indel_end_matches) {
    /* Require separation from start of query */
    indel_qpos_low = min_indel_end_matches;
  }
  debug2(printf("  Indel low bound: %d\n",indel_qpos_low));


  /* 2. Find high bound for indel pos based on nmismatches_left */
  if (nmismatches_left < max_mismatches_allowed) {
    indel_qpos_high = mismatch_positions_left[nmismatches_left];
  } else {
    indel_qpos_high = mismatch_positions_left[max_mismatches_allowed];
    debug2(printf("indel_qpos_high is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_high));
  }

  debug2(printf("  Comparing indel high bound %d with pos3 %d - 1\n",indel_qpos_high,pos3));
  if (indel_qpos_high > pos3 - 1) {
    /* Require separation from next segment */
    indel_qpos_high = pos3 - 1;
  }
  if (indel_qpos_high > querylength - min_indel_end_matches) {
    /* Require separation from end of query */
    indel_qpos_high = querylength - min_indel_end_matches;
  }
  debug2(printf("  Indel high bound: %d\n",indel_qpos_high));

  if (indel_qpos_low > indel_qpos_high) {
    debug2(printf("Indel_resolve_middle_deletion returning -1 before search\n"));
    return -1;
  }
 

  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

#ifdef DEBUG2
  gbuffer = (char *) CALLOC(querylength-indels+1,sizeof(char));
  left = univdiagonal_i - querylength;
  Genome_fill_buffer(left,querylength-indels,gbuffer);
  printf("solve_middle_indel, plusp %d, want low %d, deletion (indels %d), max_mismatches_allowed %d: Getting genome at left %llu\n",
	 plusp,want_lowest_coordinate_p,indels,max_mismatches_allowed,(unsigned long long) left);
  printf("g1: %s\n",gbuffer);
  printf("q:  ");
  Compress_print_queryseq(query_compress,/*pos5*/0,/*pos3*/querylength);
  printf("\n");
  printf("g2: %s\n",&(gbuffer[-indels]));
  FREE(gbuffer);
#endif


  if (want_lowest_coordinate_p == true) {
    best_indel_pos_case1 = best_indel_pos_case2 = querylength;
  } else {
    best_indel_pos_case1 = best_indel_pos_case2 = -1;
  }

  /* Case 1: Use nmismatches_right as decision points */
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] > indel_qpos_high - 1) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    /* Subtract 1 because indel_pos is mismatch_positions_right + 1 */
    righti++;
  }

  lefti = nmismatches_left - 1;
  while (lefti >= 0 && mismatch_positions_left[lefti] > indel_qpos_high) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    lefti--;
  }

  /* Search in descending qpos */
  best_sum = querylength;
  while (righti < nmismatches_right && mismatch_positions_right[righti] >= indel_qpos_low - 1) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
      lefti--;
    }
    sum = righti + lefti + 1;

    debug2(printf("  (Case 1) sum %d=%d+%d at indel_pos %d.",
		  sum,righti,lefti+1,mismatch_positions_right[righti]+1));
    if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == true)) {
      indel_pos = mismatch_positions_right[righti] + 1;
      assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);
      best_indel_pos_case1 = indel_pos;
      nmismatches_right_case1 = righti;
      nmismatches_left_case1 = lefti + 1;
      debug2(printf("**"));
      best_sum = sum;
    }
    righti++;
  }
  debug2(printf("\n"));


  /* Case 2: Use nmismatches_left as decision points */
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < indel_qpos_low) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    lefti++;
  }
  
  righti = nmismatches_right - 1;
  while (righti >= 0 && mismatch_positions_right[righti] < indel_qpos_low - 1) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    /* Subtract 1 because indel_pos is mismatch_positions_right + 1 */
    righti--;
  }
  
  /* Search in ascending qpos */
  best_sum = querylength;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] <= indel_qpos_high) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
      righti--;
    }
    sum = lefti + righti + 1;
    
    debug2(printf("  (Case 2) sum %d=%d+%d at indel_pos %d.",
		  sum,lefti,righti+1,mismatch_positions_left[lefti]));
    indel_pos = mismatch_positions_left[lefti];
    assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);
    if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == false)) {
      best_indel_pos_case2 = indel_pos;
      nmismatches_left_case2 = lefti;
      nmismatches_right_case2 = righti + 1;
      debug2(printf("**"));
      best_sum = sum;
    }
    lefti++;
  }
  debug2(printf("\n"));

  if (want_lowest_coordinate_p == true) {
    if (best_indel_pos_case1 <= best_indel_pos_case2) {
      debug2(printf("%d <= %d, so case1 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
    } else {
      debug2(printf("%d > %d, so case2 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
    }
  } else {
    if (best_indel_pos_case1 >= best_indel_pos_case2) {
      debug2(printf("%d >= %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
    } else {
      debug2(printf("%d < %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
    }
  }

  if (best_indel_pos == -1 || best_indel_pos == querylength) {
    debug2(printf("Indel_resolve_middle_deletion returning -1 after search\n"));
    return -1;
  }

  if (maskedp == false) {
    *best_nmismatches_i = *best_ref_nmismatches_i;
    *best_nmismatches_j = *best_ref_nmismatches_j;

  } else {
    /* Need to get masked information */
    *best_nmismatches_i =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_i),genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,
					    pos5,/*pos3*/best_indel_pos,plusp,genestrand);
    *best_nmismatches_j =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_j),genomebits,genomebits_alt,query_compress,
					    univdiagonal_j,querylength,
					    /*pos5*/best_indel_pos,pos3,plusp,genestrand);
  }

  debug2(printf("Indel_resolve_middle_deletion returning %d with nmismatches %d+%d\n",
		best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
  return best_indel_pos;
}


#if 0
/* TODO: Make use of knownindels */
int
Indel_resolve_middle_deletion_or_splice (int *best_introntype, int *best_nmismatches_i, int *best_nmismatches_j,
					 int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
					 Univcoord_T univdiagonal_i, int indels,
					 Genomebits_T ome, Genomebits_T ome_alt, Compress_T query_compress,
					 Univcoord_T chroffset, Univcoord_T chrhigh,
					 int pos5, int pos3, int querylength, Indelinfo_T indelinfo,
					 bool plusp, int genestrand, int min_intronlength,
					 int max_mismatches_allowed, bool want_lowest_coordinate_p) {
  int best_indel_pos, best_indel_pos_case1, best_indel_pos_case2,
    indel_qpos_low, indel_qpos_high, indel_pos;
  int nmismatches_left, nmismatches_right;
#ifdef DEBUG2
  char *gbuffer2;
  int querystart, queryend;
  int i;
#endif
  int nmismatches_left_case1, nmismatches_right_case1,
    nmismatches_left_case2, nmismatches_right_case2;
  int introntype_case1, introntype_case2;
  int best_sum, sum, lefti, righti;
  int best_intron_level, intron_level;
  int introntype;

  int *mismatch_positions_left, *mismatch_positions_right;
  char *gbuffer, left1, left2, right2, right1;
  int length;
  Univcoord_T univdiagonal_j;


  univdiagonal_j = univdiagonal_i - indels;

  debug2(printf("Entered Indel_resolve_middle_deletion_or_splice at %u with pos5 %u and pos3 %u\n",
		left,pos5,pos3));

  mismatch_positions_left = indelinfo->mismatch_positions_left; /* Use allocated memory */
  mismatch_positions_right = indelinfo->mismatch_positions_right; /* Use allocated memory */


  /* query has deletion.  Get |indels| more from genome; add to right. */
  /* left = ptr->diagonal - querylength; */

  assert(indels < 0);		/* Deletions and splices are negative */

#ifdef DEBUG2
  queryend = querylength - indels; /* Extends past querylength */
  gbuffer2 = (char *) CALLOC(queryend+1,sizeof(char));

  /* univdiagonal = left + (Univcoord_T) querylength; */
  /* queryend = (univdiagonal + queryend <= chrhigh + querylength) ? queryend : (int) (chrhigh - left); */
  Genome_fill_buffer(left + pos5,/*length*/queryend - pos5,gbuffer2);
  queryseq = Compress_queryseq(query_compress,querylength);

  debug2(printf("solve_middle_indel, plusp %d, want low %d, deletion (indels %d), max_mismatches_allowed %d: Getting genome at diagonal - querylength %d = %llu\n",
		plusp,want_lowest_coordinate_p,indels,max_mismatches_allowed,querylength,(unsigned long long) left));
  printf("g1: %s\n",gbuffer2);
  printf("q:  %s\n",queryseq);
  printf("g2: %s\n",&(gbuffer2[-indels]));
  FREE(queryseq);
  FREE(gbuffer2);
#endif

    
#ifdef TRIM_AT_CHROMOSOME_BOUNDS
  Univcoord_T univdiagonali, univdiagonalj;
  univdiagonali = left + (Univcoord_T) querylength;
  univdiagonalj = (left - indels) + (Univcoord_T) querylength;
  pos5 = (univdiagonali + pos5 >= chroffset + (Univcoord_T) querylength) ? pos5 : (int) (chroffset - left);
  pos3 = (univdiagonalj + pos3 <= chrhigh + (Univcoord_T) querylength) ? pos3 : (int) (chrhigh - (left - indels));
  /* Genome_fill_buffer(left + querystart,queryend - querystart,gbuffer); */
#endif

  debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromleft over %d..%d\n",
		max_mismatches_allowed,pos5,pos3));

  if (pos5 >= pos3) {
    return -1;
  } else {
    nmismatches_left = Genomebits_mismatches_fromleft(mismatch_positions_left,max_mismatches_allowed,
						  ome,ome_alt,query_compress,univdiagonal_i,querylength,
						  pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_left <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
    debug2(
	   printf("%d mismatches on left in %d..%d at:",nmismatches_left,pos5,pos3);
	   for (i = 0; i <= nmismatches_left; i++) {
	     printf(" %d",mismatch_positions_left[i]);
	   }
	   printf("\n");
	 );

    /* No need to check chromosome bounds */
    debug2(printf("max_mismatches_allowed is %d.  Calling Genomebits_mismatches_fromright over %d..%d\n",
		  max_mismatches_allowed,pos5,pos3));
    
    nmismatches_right = Genomebits_mismatches_fromright(mismatch_positions_right,max_mismatches_allowed,
						    ome,ome_alt,query_compress,univdiagonal_j,querylength,
						    pos5,pos3,plusp,genestrand);
    /* assert(nmismatches_right <= max_mismatches_allowed + 1); -- Doesn't work with SIMD decode procedures */
    debug2(
	   printf("%d mismatches on right in %d..%d at:",nmismatches_right,pos5,pos3);
	   for (i = 0; i <= nmismatches_right; i++) {
	     printf(" %d",mismatch_positions_right[i]);
	   }
	   printf("\n");
	   );
  }


  /* Find low bound for indel pos based on nmismatches_right */
  if (nmismatches_right < max_mismatches_allowed) {
    indel_qpos_low = mismatch_positions_right[nmismatches_right];
  } else {
    indel_qpos_low = mismatch_positions_right[max_mismatches_allowed];
    debug2(printf("indel_qpos_low is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_low));
  }
  if (indel_qpos_low < pos5 + 1) {
    /* Require separation from previous segment */
    indel_qpos_low = pos5 + 1;
  }
  if (indel_qpos_low < min_indel_end_matches) {
    /* Require separation from start of query */
    indel_qpos_low = min_indel_end_matches;
  }
  debug2(printf("  Indel low bound: %d\n",indel_qpos_low));

  /* Find high bound for indel pos based on nmismatches_left */
  if (nmismatches_left < max_mismatches_allowed) {
    indel_qpos_high = mismatch_positions_left[nmismatches_left];
  } else {
    indel_qpos_high = mismatch_positions_left[max_mismatches_allowed];
    debug2(printf("indel_qpos_high is determined by max_mismatches_allowed %d => %d\n",
		  max_mismatches_allowed,indel_qpos_high));
  }
  if (indel_qpos_high > pos3 - 1) {
    /* Require separation from next segment */
    indel_qpos_high = pos3 - 1;
  }
  if (indel_qpos_high > querylength - min_indel_end_matches) {
    /* Require separation from end of query */
    indel_qpos_high = querylength - min_indel_end_matches;
  }
  debug2(printf("  Indel high bound: %d\n",indel_qpos_high));

  if (indel_qpos_low > indel_qpos_high) {
    debug2(printf("Indel_resolve_middle_deletion_or_splice returning -1 before search\n"));
    return -1;
  } else if (-indels >= min_intronlength) {
    /* Obtaining gbuffer is less expensive than repeated calls to Genome_bits_get_char */
    /* lowest nucleotide position would be indel_qpos_low */
    /* highest nucleotide position would be indel_qpos_high - indels - 1 */
    length = (indel_qpos_high - indels - 1) - indel_qpos_low + 1;
    gbuffer = (char *) MALLOC((length + 1) * sizeof(char));
    Genome_fill_buffer(left + indel_qpos_low,length,gbuffer);
  }


  if (want_lowest_coordinate_p == true) {
    best_indel_pos_case1 = best_indel_pos_case2 = querylength;
  } else {
    best_indel_pos_case1 = best_indel_pos_case2 = -1;
  }

  /* Case 1: Use nmismatches_right as decision points */
  righti = 0;
  while (righti < nmismatches_right && mismatch_positions_right[righti] > indel_qpos_high - 1) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    /* Subtract 1 because indel_pos is mismatch_positions_right + 1 */
    righti++;
  }

  lefti = nmismatches_left - 1;
  while (lefti >= 0 && mismatch_positions_left[lefti] > indel_qpos_high) {
    /* Skip leftward to the feasible region, where the number of mismatches is accurate */
    lefti--;
  }

  /* Search in descending qpos */
  best_sum = querylength;
  best_intron_level = Intron_level(NONINTRON);

  while (righti < nmismatches_right && mismatch_positions_right[righti] >= indel_qpos_low - 1) {
    while (lefti >= 0 && mismatch_positions_left[lefti] > mismatch_positions_right[righti]) {
      lefti--;
    }
    sum = righti + lefti + 1;

    indel_pos = mismatch_positions_right[righti] + 1;
    assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);

    if (-indels < min_intronlength) {
      debug2(printf("  (Case 1 del) sum %d=%d+%d at indel_pos %d.",
		    sum,righti,lefti+1,mismatch_positions_right[righti]+1));
      if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == true)) {
	indel_pos = mismatch_positions_right[righti] + 1;
	assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);
	best_indel_pos_case1 = indel_pos;
	nmismatches_right_case1 = righti;
	nmismatches_left_case1 = lefti + 1;
	introntype_case1 = NONINTRON;
	debug2(printf("**"));
	best_sum = sum;
	best_intron_level = Intron_level(NONINTRON);
      }

    } else if (sum > best_sum) {
      /* Skip getting intron info */
      debug2(printf("  (Case 1 intron) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
    } else {
      /* Account for introntype in cases of ties */
      left1 = gbuffer[indel_pos - indel_qpos_low];
      left2 = gbuffer[indel_pos+1 - indel_qpos_low];
      right2 = gbuffer[indel_pos-indels-2 - indel_qpos_low];
      right1 = gbuffer[indel_pos-indels-1 - indel_qpos_low];
      /* assert(left1 == Genome_bits_get_char(left+indel_pos)); */
      /* assert(left2 == Genome_bits_get_char(left+indel_pos+1)); */
      /* assert(right2 == Genome_bits_get_char(left+indel_pos-indels-2)); */
      /* assert(right1 == Genome_bits_get_char(left+indel_pos-indels-1)); */

      introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
      intron_level = Intron_level(introntype);
      debug2(printf("  (Case 1 intron) sum %d=%d+%d at indel_pos %d (%c%c-%c%c, type %s).",
		    sum,righti,lefti+1,mismatch_positions_right[righti]+1,
		    left1,left2,right2,right1,Intron_type_string(introntype)));
      if (sum < best_sum ||
	  (sum == best_sum && intron_level > best_intron_level) ||
	  (sum == best_sum && intron_level == best_intron_level &&
	   want_lowest_coordinate_p == true)) {
	indel_pos = mismatch_positions_right[righti] + 1;

	best_indel_pos_case1 = indel_pos;
	nmismatches_right_case1 = righti;
	nmismatches_left_case1 = lefti + 1;
	introntype_case1 = introntype;
	debug2(printf("**"));
	best_sum = sum;
	best_intron_level = intron_level;
      }
    }
    righti++;
  }
  debug2(printf("\n"));


  /* Case 2: Use nmismatches_left as decision points */
  lefti = 0;
  while (lefti < nmismatches_left && mismatch_positions_left[lefti] < indel_qpos_low) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    lefti++;
  }

  righti = nmismatches_right - 1;
  while (righti >= 0 && mismatch_positions_right[righti] < indel_qpos_low - 1) {
    /* Skip rightward to the feasible region, where the number of mismatches is accurate */
    /* Subtract 1 because indel_pos is mismatch_positions_right + 1 */
    righti--;
  }

  /* Search in ascending qpos */
  best_sum = querylength;
  best_intron_level = Intron_level(NONINTRON);

  while (lefti < nmismatches_left && mismatch_positions_left[lefti] <= indel_qpos_high) {
    while (righti >= 0 && mismatch_positions_right[righti] < mismatch_positions_left[lefti]) {
      righti--;
    }
    sum = lefti + righti + 1;

    indel_pos = mismatch_positions_left[lefti];
    assert(indel_pos >= indel_qpos_low && indel_pos <= indel_qpos_high);

    if (-indels < min_intronlength) {
      debug2(printf("  (Case 2 del) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
      if (sum < best_sum || (sum == best_sum && want_lowest_coordinate_p == false)) {
	best_indel_pos_case2 = indel_pos;
	nmismatches_left_case2 = lefti;
	nmismatches_right_case2 = righti + 1;
	introntype_case2 = NONINTRON;
	debug2(printf("**"));
	best_sum = sum;
	best_intron_level = Intron_level(NONINTRON);
      }

    } else if (sum > best_sum) {
      /* Skip getting intron info */
      debug2(printf("  (Case 2 intron) sum %d=%d+%d at indel_pos %d.",
		    sum,lefti,righti+1,mismatch_positions_left[lefti]));
    } else {
      /* Account for introntype in cases of ties */
      left1 = gbuffer[indel_pos - indel_qpos_low];
      left2 = gbuffer[indel_pos+1 - indel_qpos_low];
      right2 = gbuffer[indel_pos-indels-2 - indel_qpos_low];
      right1 = gbuffer[indel_pos-indels-1 - indel_qpos_low];
      /* assert(left1 == Genome_bits_get_char(left+indel_pos)); */
      /* assert(left2 == Genome_bits_get_char(left+indel_pos+1)); */
      /* assert(right2 == Genome_bits_get_char(left+indel_pos-indels-2)); */
      /* assert(right1 == Genome_bits_get_char(left+indel_pos-indels-1)); */

      introntype = Intron_type(left1,left2,right2,right1,left1,left2,right2,right1,/*cdna_direction*/0);
      intron_level = Intron_level(introntype);
      debug2(printf("  (Case 2 intron) sum %d=%d+%d at indel_pos %d (%c%c-%c%c, type %s).",
		    sum,lefti,righti+1,mismatch_positions_left[lefti],
		    left1,left2,right2,right1,Intron_type_string(introntype)));
      if (sum < best_sum ||
	  (sum == best_sum && intron_level > best_intron_level) ||
	  (sum == best_sum && intron_level == best_intron_level &&
	   want_lowest_coordinate_p == false)) {
	best_indel_pos_case2 = indel_pos;
	nmismatches_left_case2 = lefti;
	nmismatches_right_case2 = righti + 1;
	introntype_case2 = introntype;
	debug2(printf("**"));
	best_sum = sum;
	best_intron_level = intron_level;
      }
    }
    lefti++;
  }
  debug2(printf("\n"));

  if (want_lowest_coordinate_p == true) {
    if (best_indel_pos_case1 <= best_indel_pos_case2) {
      debug2(printf("%d <= %d, so case1 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
      *best_introntype = introntype_case1;
    } else {
      debug2(printf("%d > %d, so case2 is lowest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
      *best_introntype = introntype_case2;
    }
  } else {
    if (best_indel_pos_case1 >= best_indel_pos_case2) {
      debug2(printf("%d >= %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case1;
      *best_ref_nmismatches_i = nmismatches_left_case1;
      *best_ref_nmismatches_j = nmismatches_right_case1;
      *best_introntype = introntype_case1;
    } else {
      debug2(printf("%d < %d, so case1 is highest\n",best_indel_pos_case1,best_indel_pos_case2));
      best_indel_pos = best_indel_pos_case2;
      *best_ref_nmismatches_i = nmismatches_left_case2;
      *best_ref_nmismatches_j = nmismatches_right_case2;
      *best_introntype = introntype_case2;
    }
  }

  if (-indels >= min_intronlength) {
    FREE(gbuffer);
  }


  if (best_indel_pos == -1 || best_indel_pos == querylength) {
    debug2(printf("Indel_resolve_middle_deletion_or_splice returning -1 after search\n"));
    return -1;
  }

  if (maskedp == false) {
    *best_nmismatches_i = *best_ref_nmismatches_i;
    *best_nmismatches_j = *best_ref_nmismatches_j;
  } else if (best_indel_pos < 0) {
    *best_nmismatches_i = *best_ref_nmismatches_i;
    *best_nmismatches_j = *best_ref_nmismatches_j;
  } else {
    *best_nmismatches_i =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_i),genomebits,genomebits_alt,query_compress,
					    univdiagonal_i,querylength,
					    pos5,/*pos3*/best_indel_pos,plusp,genestrand);
    *best_nmismatches_j =
      Genomebits_count_mismatches_substring(&(*best_ref_nmismatches_j),genomebits,genomebits_alt,query_compress,
					    univdiagonal_j,querylength,
					    /*pos5*/best_indel_pos,pos3,plusp,genestrand);
  }

  debug2(printf("Indel_resolve_middle_deletion_or_splice returning %d with nmismatches %d+%d\n\n",
		best_indel_pos,*best_nmismatches_i,*best_nmismatches_j));
  return best_indel_pos;
}
#endif


