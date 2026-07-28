static char rcsid[] = "$Id: 98f237f3ec23f7d73036cd3c52ecad334f9f7f8f $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_trim.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"

#include "simd.h"


#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Previously set mismatch score to be -4, but poor sequencing quality does better with -3 */
#define TRIM_MATCH_SCORE 1
#define TRIM_MISMATCH_SCORE_MULT -3 /* Requires 3 matches to compensate */


static Diffproc_32_T block_diff_32;
static Diffproc_snp_32_T block_diff_snp_32;

static Diffproc_64_T block_diff_64; 
static Diffproc_snp_64_T block_diff_snp_64; 

#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_128; 
static Diffproc_snp_128_T block_diff_snp_128; 
#endif


#define T Genomebits_T



/* Integrates former Spliceends_trim_qend_nosplice.  Returns trimpos */
/* We set query_unk_mismatch_p to false because trimming query N's can
   affect the --clip-overlap feature.  Also, some sequences are poor
   with lots of N's, indicating that the sequence is not reliable */
/* genome_unk_mismatch_p needs to be true so circular alignments around the origin are favored */
/* Follows nmismatches_fromleft from genomebits_mismatches.c */
int
Genomebits_trim_qend (int *nmismatches_to_trimpos,
		      Compress_T query_compress, T ref,
		      Univcoord_T univdiagonal, int querylength,
		      int pos5, int pos3, bool plusp, int genestrand) {

  int trimpos = pos3, prevpos = pos5 - 1, pos;
  int max_score = (pos3 - pos5)*TRIM_MISMATCH_SCORE_MULT, score = 0;
  bool query_unk_mismatch_p = false;
  bool genome_unk_mismatch_p = true;

  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
  int relpos;

  assert(pos5 >= 0);
  assert(pos5 <= pos3);		/* Can be equal in get_exhaustive */
  assert(pos3 <= querylength);

#if 0
  if (allow_soft_clips_p == false) {
  }
#endif

  *nmismatches_to_trimpos = 0;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug1(
	printf("\n\n");
	printf("Entered Genomebits_trim_qend\n");
	printf("Genome (in mismatches_fromleft) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug1(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug1(printf("Query shifted %d:\n",nshift));
  debug1(Compress_print(query_compress,nshift,pos5,pos3));

  /* For rightward scanning */
  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5;
  debug1(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug1(printf("bothdisc:   %08X\n",diff_32));

    while (nonzero_p_32(diff_32)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_32(diff_32));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos3 - pos) < 0) {
	debug1(printf(" redemption: %d => terminate",pos3 - pos));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos,*nmismatches_to_trimpos));
	return trimpos;		/* At the mismatch for qend */
      }
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_32 = clear_lowbit_32(diff_32,relpos);
    }

    if (*nmismatches_to_trimpos == nmismatches - 1) {
      /* If last mismatch compensated for previous, then take the last
	 segment, regardless of whether it compensates for the last
	 mismatch */
      debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,nmismatches));
      trimpos = pos3;
      *nmismatches_to_trimpos += 1;

    } else {
      /* Final segment */
      pos = pos3;
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	/* max_score = score; */
#if 0
      } else if (score + /*redemption*/(pos3 - pos) < 0) {
	debug1(printf(" redemption: %d => terminate",pos3 - pos));
	donep = true;
#endif
      }
      debug1(printf("\n"));
    }

    debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		 trimpos,*nmismatches_to_trimpos));
    return trimpos;		/* At the mismatch for qend */

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug1(printf("bothdisc:   %016lX\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_64(diff_64));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos3 - pos) < 0) {
	debug1(printf(" redemption: %d => terminate",pos3 - pos));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos,*nmismatches_to_trimpos));
	return trimpos;		/* At the mismatch for qend */
      }
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_64 = clear_lowbit_64(diff_64,relpos);
    }

    if (*nmismatches_to_trimpos == nmismatches - 1) {
      /* If last mismatch compensated for previous, then take the last
	 segment, regardless of whether it compensates for the last
	 mismatch */
      debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,nmismatches));
      trimpos = pos3;
      *nmismatches_to_trimpos += 1;

    } else {
      /* Final segment */
      pos = pos3;
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	/* max_score = score; */
#if 0
      } else if (score + /*redemption*/(pos3 - pos) < 0) {
	debug1(printf(" redemption: %d => terminate",pos3 - pos));
	donep = true;
#endif
      }
      debug1(printf("\n"));
    }

    debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		 trimpos,*nmismatches_to_trimpos));
    return trimpos;		/* At the mismatch for qend */

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug1(printf("startdisc:  %016lX\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_64(diff_64));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos3 - pos) < 0) {
	debug1(printf(" redemption: %d => terminate",pos3 - pos));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos,*nmismatches_to_trimpos));
	return trimpos;		/* At the mismatch for qend */
      }
	
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_64 = clear_lowbit_64(diff_64,relpos);
    }

    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    offset += 64;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug1(printf("nodisc:     %016lX\n",diff_64));
      
      while (nonzero_p_64(diff_64)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_64(diff_64));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos3 - pos) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos3 - pos));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos,*nmismatches_to_trimpos));
	  return trimpos;		/* At the mismatch for qend */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_64 = clear_lowbit_64(diff_64,relpos);
      }

      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug1(printf("enddisc:    %016lX\n",diff_64));

      while (nonzero_p_64(diff_64)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_64(diff_64));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos3 - pos) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos3 - pos));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos,*nmismatches_to_trimpos));
	  return trimpos;		/* At the mismatch for qend */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_64 = clear_lowbit_64(diff_64,relpos);
      }

      if (*nmismatches_to_trimpos == nmismatches - 1) {
	/* If last mismatch compensated for previous, then take the last
	   segment, regardless of whether it compensates for the last
	   mismatch */
	debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		      *nmismatches_to_trimpos,nmismatches));
	trimpos = pos3;
	*nmismatches_to_trimpos += 1;
      } else {
	/* Final segment */
	pos = pos3;
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  /* max_score = score; */
#if 0
	} else if (score + /*redemption*/(pos3 - pos) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos3 - pos));
	  donep = true;
#endif
	}
	debug1(printf("\n"));
      }

      debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		   trimpos,*nmismatches_to_trimpos));
      return trimpos;		/* At the mismatch for qend */

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug1(printf("enddisc:    %08X\n",diff_32));

      while (nonzero_p_32(diff_32)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset + (relpos = count_trailing_zeroes_32(diff_32));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos3 - pos) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos3 - pos));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos,*nmismatches_to_trimpos));
	  return trimpos;		/* At the mismatch for qend */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_32 = clear_lowbit_32(diff_32,relpos);
      }

      if (*nmismatches_to_trimpos == nmismatches - 1) {
	/* If last mismatch compensated for previous, then take the last
	   segment, regardless of whether it compensates for the last
	   mismatch */
	debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		      *nmismatches_to_trimpos,nmismatches));
	trimpos = pos3;
	*nmismatches_to_trimpos += 1;
      } else {
	/* Final segment */
	pos = pos3;
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  /* max_score = score; */
#if 0
	} else if (score + /*redemption*/(pos3 - pos) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos3 - pos));
	  donep = true;
#endif
	}
	debug1(printf("\n"));
      }

      debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		   trimpos,*nmismatches_to_trimpos));
      return trimpos;		/* At the mismatch for qend */
    }
  }
}


/* Integrates former Spliceends_trim_qstart_nosplice.  Returns trimpos */
/* We set query_unk_mismatch_p to false because trimming query N's can
   affect the --clip-overlap feature.  Also, some sequences are poor
   with lots of N's, indicating that the sequence is not reliable */
/* genome_unk_mismatch_p needs to be true so circular alignments around the origin are favored */
/* Follows nmismatches_fromright from genomebits_mismatches.c */
int
Genomebits_trim_qstart (int *nmismatches_to_trimpos,
			Compress_T query_compress, T ref,
			Univcoord_T univdiagonal, int querylength,
			int pos5, int pos3, bool plusp, int genestrand) {

  int trimpos = pos5, prevpos = pos3, pos;
  int max_score = (pos3 - pos5)*TRIM_MISMATCH_SCORE_MULT, score = 0;
  bool query_unk_mismatch_p = false;
  bool genome_unk_mismatch_p = true;

  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *start_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
  int relpos;

  assert(pos5 <= pos3);		/* Can be equal in get_exhaustive */
  assert(pos3 <= querylength);


#if 0
  if (allow_soft_clips_p == false) {
  }
#endif

  *nmismatches_to_trimpos = 0;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug1(
	printf("\n\n");
	printf("Entered Genomebits_trim_qstart\n");
	printf("Genome (in mismatches_fromright) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  endblocki = (left+pos3)/32U;
  startblocki = (left+pos5)/32U;

  debug1(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos3/*or pos3 - 1?*/);
  debug1(printf("Query shifted %d:\n",nshift));
  debug1(Compress_print(query_compress,nshift,pos5,pos3));

  /* For leftward scanning */
  ref_high_ptr = &(ref->high_blocks[endblocki]);
  ref_low_ptr = &(ref->low_blocks[endblocki]);
  ref_flags_ptr = &(ref->flags_blocks[endblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = (32 - enddiscard) + pos3 - 1;
  debug1(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = clear_start_32(diff_32,startdiscard);
    debug1(printf("bothdisc:   %08X\n",diff_32));

    while (nonzero_p_32(diff_32)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_32(diff_32));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos + 1,*nmismatches_to_trimpos));
	return trimpos + 1;		/* One position after the mismatch for qstart */
      }
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_32 = clear_highbit_32(diff_32,relpos);
    }

    if (pos3 - pos5 > nmismatches && *nmismatches_to_trimpos == nmismatches - 1) {
      /* If there was a match, and last mismatch compensated for
	 previous, then take the last segment, regardless of whether
	 it compensates for the last mismatch */
      debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,nmismatches));
      trimpos = pos5 - 1;
      *nmismatches_to_trimpos += 1;

    } else {
      /* Final segment */
      pos = pos5 - 1;
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	/* max_score = score; */
#if 0
      } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	donep = true;
#endif
      }
      debug1(printf("\n"));
    }

    debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		 trimpos + 1,*nmismatches_to_trimpos));
    return trimpos + 1;		/* One position after the mismatch for qstart */

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug1(printf("bothdisc:   %016lX\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_64(diff_64));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos + 1,*nmismatches_to_trimpos));
	return trimpos + 1;		/* One position after the mismatch for qstart */
      }
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_64 = clear_highbit_64(diff_64,relpos);
    }

    if (pos3 - pos5 > nmismatches && *nmismatches_to_trimpos == nmismatches - 1) {
      /* If there was a match, and last mismatch compensated for
	 previous, then take the last segment, regardless of whether
	 it compensates for the last mismatch */
      debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		  *nmismatches_to_trimpos,nmismatches));
      trimpos = pos5 - 1;
      *nmismatches_to_trimpos += 1;

    } else {
      /* Final segment */
      pos = pos5 - 1;
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	/* max_score = score; */
#if 0
      } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	donep = true;
#endif
      }
      debug1(printf("\n"));
    }

    debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		 trimpos + 1,*nmismatches_to_trimpos));
    return trimpos + 1;		/* One position after the mismatch for qstart */

  } else {
    /* Multiple words */
    start_ptr = &(ref->high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug1(printf("enddisc:    %016lX\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_64(diff_64));
      score += TRIM_MISMATCH_SCORE_MULT;
      score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
      debug1(printf("pos %d, score %d",pos,score));
      if (score >= max_score) {
	debug1(printf(" **"));
	trimpos = pos;
	*nmismatches_to_trimpos = nmismatches;
	max_score = score;
      } else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		      trimpos + 1,*nmismatches_to_trimpos));
	return trimpos + 1;		/* One position after the mismatch for qstart */
      }
      debug1(printf("\n"));
      prevpos = pos;
      nmismatches++;

      diff_64 = clear_highbit_64(diff_64,relpos);
    }

    query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
    ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
    offset -= 64;

    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug1(printf("nodisc:     %016lX\n",diff_64));

      while (nonzero_p_64(diff_64)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_64(diff_64));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos + 1,*nmismatches_to_trimpos));
	  return trimpos + 1;		/* One position after the mismatch for qstart */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_64 = clear_highbit_64(diff_64,relpos);
      }

      query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
      ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
      offset -= 64;
    }

    if (ref_high_ptr == start_ptr + 1) {
      /* Start 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_start_64(diff_64,startdiscard);
      debug1(printf("startdisc:  %016lX\n",diff_64));

      while (nonzero_p_64(diff_64)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_64(diff_64));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos + 1,*nmismatches_to_trimpos));
	  return trimpos + 1;		/* One position after the mismatch for qstart */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_64 = clear_highbit_64(diff_64,relpos);
      }

      if (*nmismatches_to_trimpos == nmismatches - 1) {
	/* If last mismatch compensated for previous, then take the last
	   segment, regardless of whether it compensates for the last
	   mismatch */
	debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		      *nmismatches_to_trimpos,nmismatches));
	trimpos = pos5 - 1;
	*nmismatches_to_trimpos += 1;
	
      } else {
	/* Final segment */
	pos = pos5 - 1;
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  /* max_score = score; */
#if 0
	} else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	  donep = true;
#endif
	}
	debug1(printf("\n"));
      }

      debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		   trimpos + 1,*nmismatches_to_trimpos));
      return trimpos + 1;		/* One position after the mismatch for qstart */

    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_start_32(diff_32,startdiscard);
      debug1(printf("startdisc:  %08X\n",diff_32));

      while (nonzero_p_32(diff_32)) {
	pos = /*mismatch_positions[nmismatches++] =*/ offset - (relpos = count_leading_zeroes_32(diff_32));
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  max_score = score;
	} else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	  debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
			trimpos + 1,*nmismatches_to_trimpos));
	  return trimpos + 1;		/* One position after the mismatch for qstart */
	}
	debug1(printf("\n"));
	prevpos = pos;
	nmismatches++;

	diff_32 = clear_highbit_32(diff_32,relpos);
      }

      if (*nmismatches_to_trimpos == nmismatches - 1) {
	/* If last mismatch compensated for previous, then take the last
	   segment, regardless of whether it compensates for the last
	   mismatch */
	debug1(printf("Last mismatch compensates because %d nmismatches == total %d - 1\n",
		      *nmismatches_to_trimpos,nmismatches));
	trimpos = pos5 - 1;
	*nmismatches_to_trimpos += 1;
	
      } else {
	/* Final segment */
	pos = pos5 - 1;
	score += TRIM_MISMATCH_SCORE_MULT;
	score += (prevpos - pos - 1)*TRIM_MATCH_SCORE;
	debug1(printf("pos %d, score %d",pos,score));
	if (score >= max_score) {
	  debug1(printf(" **"));
	  trimpos = pos;
	  *nmismatches_to_trimpos = nmismatches;
	  /* max_score = score; */
#if 0
	} else if (score + /*redemption*/(pos + 1 - pos5) < 0) {
	  debug1(printf(" redemption: %d => terminate",pos + 1 - pos5));
	  donep = true;
#endif
	}
	debug1(printf("\n"));
      }

      debug1(printf("Returning trimpos %d, %d nmismatches to trimpos\n",
		   trimpos + 1,*nmismatches_to_trimpos));
      return trimpos + 1;		/* One position after the mismatch for qstart */
    }
  }
}


void
Genomebits_trim_setup (Mode_T mode, bool maskedp) {

  switch (mode) {
  case STANDARD:
    block_diff_32 = block_diff_standard_32;
    block_diff_64 = block_diff_standard_64;
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_standard_128;
#endif
    break;

  case CMET_STRANDED: case CMET_NONSTRANDED:
    block_diff_32 = block_diff_cmet_32;
    block_diff_64 = block_diff_cmet_64;
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_cmet_128;
#endif
    break;

  case ATOI_STRANDED: case ATOI_NONSTRANDED:
    block_diff_32 = block_diff_atoi_32;
    block_diff_64 = block_diff_atoi_64;
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_atoi_128;
#endif
    break;

  case TTOC_STRANDED: case TTOC_NONSTRANDED:
    block_diff_32 = block_diff_ttoc_32;
    block_diff_64 = block_diff_ttoc_64;
#ifdef HAVE_SSE2
    block_diff_128 = block_diff_ttoc_128;
#endif
    break;
  default: fprintf(stderr,"Mode %d not recognized\n",mode); abort();
  }

#ifndef GSNAP
  if (maskedp == true) {
    block_diff_snp_32 = block_diff_standard_masked_32;
    block_diff_snp_64 = block_diff_standard_masked_64;
  } else {
    block_diff_snp_32 = block_diff_standard_wildcard_32;
    block_diff_snp_64 = block_diff_standard_wildcard_64;
  }
#ifdef HAVE_SSE2
  if (maskedp == true) {
    block_diff_snp_128 = block_diff_standard_masked_128;
  } else {
    block_diff_snp_128 = block_diff_standard_wildcard_128;
  }
#endif

#else
  switch (mode) {
  case STANDARD:
    if (maskedp == true) {
      block_diff_snp_32 = block_diff_standard_masked_32;
      block_diff_snp_64 = block_diff_standard_masked_64;
    } else {
      block_diff_snp_32 = block_diff_standard_wildcard_32;
      block_diff_snp_64 = block_diff_standard_wildcard_64;
    }
#ifdef HAVE_SSE2
    if (maskedp == true) {
      block_diff_snp_128 = block_diff_standard_masked_128;
    } else {
      block_diff_snp_128 = block_diff_standard_wildcard_128;
    }
#endif
    break;

  case CMET_STRANDED: case CMET_NONSTRANDED:
    block_diff_snp_32 = block_diff_cmet_snp_32;
    block_diff_snp_64 = block_diff_cmet_snp_64;
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_cmet_snp_128;
#endif
    break;

  case ATOI_STRANDED: case ATOI_NONSTRANDED:
    block_diff_snp_32 = block_diff_atoi_snp_32;
    block_diff_snp_64 = block_diff_atoi_snp_64;
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_atoi_snp_128;
#endif
    break;

  case TTOC_STRANDED: case TTOC_NONSTRANDED:
    block_diff_snp_32 = block_diff_ttoc_snp_32;
    block_diff_snp_64 = block_diff_ttoc_snp_64;
#ifdef HAVE_SSE2
    block_diff_snp_128 = block_diff_ttoc_snp_128;
#endif
    break;
  default: fprintf(stderr,"Mode %d not recognized\n",mode); abort();
  }
#endif

  return;
}



