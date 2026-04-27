static char rcsid[] = "$Id: 3b377f490b8a27f201fb233bb759cc1ae2b3a027 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_consec.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* These are global values, used for alignment.  Previously for
   trimming, treated query N's as mismatches, but this is not correct
   for query N's in the middle of the read.  Also, trimming query N's
   can affect the --clip-overlap feature. */
static bool query_unk_mismatch_p = false; /* Needs to be false to work with suffix arrays */
static bool genome_unk_mismatch_p = true;

static Diffproc_32_T block_diff_32;
static Diffproc_64_T block_diff_64; 
#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_128; 
#endif


#define T Genomebits_T

/* Counts matches from pos5 to pos3 up to first mismatch.  Modified from mismatches_left */
/* Also returns the genomic char at the mismatch position, unless no
   mismatches are found from pos5 to pos3.  Caller should check if the
   return value equals (pos3 - pos5). The mismatch char is intended
   for use when a suffix array search calls this procedure. */

static char CHARTABLE[] = "ACGTXXXX";

/* Used by suffix array, which needs the mismatch char.  query unk and
   genome unk must both be marked as a mismatch for the suffix array
   to work */
int
Genomebits_consecutive_matches_wmm (char *mismatch_char, T ref, Compress_T query_compress,
				    Univcoord_T univdiagonal, int querylength,
				    int pos5, int pos3, bool plusp, int genestrand) {
  int mismatch_position_from_pos5, offset, nshift, relpos;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr, high, low, flags;
  int columni;
  int idx;
  UINT4 diff_32;
  UINT8 diff_64;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Genome (in consecutive_matches_wmm):\n");
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == plusp);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query (%s) shifted %d:\n",Compress_fwdp(query_compress) == true ? "fwd" : "rev",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard;	/* For mismatch_positions_from_pos5, we do not need to add pos5 */
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d, offset %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32,offset));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X",diff_32));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_32(diff_32)));

    if (nonzero_p_32(diff_32)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_32(diff_32));

      high = ref_high_ptr[0]; low = ref_low_ptr[0]; flags = ref_flags_ptr[0];
#ifndef HAVE_SSE2
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
      idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
      *mismatch_char = CHARTABLE[idx];

      debug(printf("returning %d - %d consecutive matches and %c\n",
		   pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
      return mismatch_position_from_pos5;

    } else {
      /* No need to return mismatch_char */
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:   %016lX",diff_64));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));

      columni = relpos / 32U;
      high = ref_high_ptr[columni]; low = ref_low_ptr[columni]; flags = ref_flags_ptr[columni];
      relpos %= 32;
#ifndef HAVE_SSE2
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
      idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
      *mismatch_char = CHARTABLE[idx];

      debug(printf("returning %d - %d consecutive matches and %c\n",
		   pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
      return mismatch_position_from_pos5;

    } else {
      /* No need to return mismatch_char */
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX",diff_64));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));

      columni = relpos / 32U;
      high = ref_high_ptr[columni]; low = ref_low_ptr[columni]; flags = ref_flags_ptr[columni];
      relpos %= 32;
#ifndef HAVE_SSE2
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
      idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
      idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
      *mismatch_char = CHARTABLE[idx];

      debug(printf("returning %d - %d consecutive matches and %c\n",
		   pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
      return mismatch_position_from_pos5;
    } else {
      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
      /* No discards */
      debug(printf("nodisc:     %016lX",diff_64));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));

	columni = relpos / 32U;
	high = ref_high_ptr[columni]; low = ref_low_ptr[columni]; flags = ref_flags_ptr[columni];
	relpos %= 32;
#ifndef HAVE_SSE2
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
	idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
	*mismatch_char = CHARTABLE[idx];

	debug(printf("returning %d - %d consecutive matches and %c\n",
		     pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
	return mismatch_position_from_pos5;
      } else {
	query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
	ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
	offset += 64;
      }
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX",diff_64));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));

	columni = relpos / 32U;
	high = ref_high_ptr[columni]; low = ref_low_ptr[columni]; flags = ref_flags_ptr[columni];
	relpos %= 32;
#ifndef HAVE_SSE2
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
	idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
	*mismatch_char = CHARTABLE[idx];

	debug(printf("returning %d - %d consecutive matches and %c\n",
		     pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
	return mismatch_position_from_pos5;

      } else {
	/* No need to return mismatch_char */
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,/*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X",diff_32));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_32(diff_32)));

      if (nonzero_p_32(diff_32)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_32(diff_32));

	high = ref_high_ptr[0]; low = ref_low_ptr[0]; flags = ref_flags_ptr[0];
#ifndef HAVE_SSE2
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
	idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
	idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
	*mismatch_char = CHARTABLE[idx];

	debug(printf("returning %d - %d consecutive matches and %c\n",
		     pos5 + mismatch_position_from_pos5,pos5,*mismatch_char));
	return mismatch_position_from_pos5;

      } else {
	/* No need to return mismatch_char */
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }
    }
  }
}


/* Used by extension search */
/* genome_unk_mismatch_p needs to be true so circular alignments across origin are favored */
int
Genomebits_consecutive_matches_rightward (T ref, Compress_T query_compress,
					  Univcoord_T univdiagonal, int querylength,
					  int pos5, int pos3, bool plusp, int genestrand) {
  int mismatch_position_from_pos5, offset, nshift, relpos;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Genome (in consecutive_matches_rightward):\n");
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == plusp);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query (%s) shifted %d:\n",Compress_fwdp(query_compress) == true ? "fwd" : "rev",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d, offset %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32,offset));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X",diff_32));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_32(diff_32)));

    if (nonzero_p_32(diff_32)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_32(diff_32));
      debug(printf("1 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
      return mismatch_position_from_pos5;

    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:   %016lX",diff_64));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));
      debug(printf("1 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
      return mismatch_position_from_pos5;

    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX",diff_64));
    debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));
      debug(printf("2 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
      return mismatch_position_from_pos5;
    } else {
      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      /* No discards */
      debug(printf("nodisc:     %016lX",diff_64));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));
	debug(printf("3 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
	return mismatch_position_from_pos5;
      } else {
	query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
	ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
	offset += 64;
      }
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX",diff_64));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_64(diff_64));
	debug(printf("4 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
	return mismatch_position_from_pos5;

      } else {
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X",diff_32));
      debug(printf(" => %d trailing zeroes\n",count_trailing_zeroes_32(diff_32)));

      if (nonzero_p_32(diff_32)) {
	mismatch_position_from_pos5 = offset + (relpos = count_trailing_zeroes_32(diff_32));
	debug(printf("4 returning %d - %d consecutive matches\n",pos5 + mismatch_position_from_pos5,pos5));
	return mismatch_position_from_pos5;

      } else {
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }
    }
  }
}


/* Used by extension search */
/* Counts matches from pos3 to pos5 up to first mismatch.  Modified from mismatches_right */
/* genome_unk_mismatch_p needs to be true so circular alignments across origin are favored */
int
Genomebits_consecutive_matches_leftward (T ref, Compress_T query_compress,
					 Univcoord_T univdiagonal, int querylength,
					 int pos5, int pos3, bool plusp, int genestrand) {
  int mismatch_position_from_pos3, offset, nshift, relpos;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *start_ptr;
  UINT4 diff_32;
  UINT8 diff_64;


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Genome (in consecutive_matches_leftward):\n");
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == plusp);

  endblocki = (left+pos3)/32U;
  startblocki = (left+pos5)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos3/*or pos3 - 1?*/);
  debug(printf("Query (%s) shifted %d:\n",Compress_fwdp(query_compress) == true ? "fwd" : "rev",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  /* For leftward scanning */
  ref_high_ptr = &(ref->high_blocks[endblocki]);
  ref_low_ptr = &(ref->low_blocks[endblocki]);
  ref_flags_ptr = &(ref->flags_blocks[endblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = (32 - enddiscard); /* For mismatch_position_from_pos3, do not need to add pos3 - 1 */
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d, offset %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32,offset));
  
  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = clear_start_32(diff_32,startdiscard);
    debug(printf("bothdisc:   %08X",diff_32));
    debug(printf(" => %d leading zeroes\n",count_leading_zeroes_32(diff_32)));

    if (nonzero_p_32(diff_32)) {
      mismatch_position_from_pos3 = (relpos = count_leading_zeroes_32(diff_32)) - offset;
      debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
      return mismatch_position_from_pos3;
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("bothdisc:   %016lX",diff_64));
    debug(printf(" => %d leading zeroes\n",count_leading_zeroes_64(diff_64)));
    
    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos3 = (relpos = count_leading_zeroes_64(diff_64)) - offset;
      debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
      return mismatch_position_from_pos3;
    } else {
      debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
      return (pos3 - pos5);
    }

  } else {
    /* Multiple words */
    start_ptr = &(ref->high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("enddisc:    %016lX",diff_64));
    debug(printf(" => %d leading zeroes\n",count_leading_zeroes_64(diff_64)));
    
    if (nonzero_p_64(diff_64)) {
      mismatch_position_from_pos3 = (relpos = count_leading_zeroes_64(diff_64)) - offset;
      debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
      return mismatch_position_from_pos3;
    } else {
      query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
      ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
      offset -= 64;
    }

    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      /* No discards */
      debug(printf("nodisc:     %016lX",diff_64));
      debug(printf(" => %d leading zeroes\n",count_leading_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos3 = (relpos = count_leading_zeroes_64(diff_64)) - offset;
	debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
	return mismatch_position_from_pos3;
      } else {
	query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
	ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
	offset -= 64;
      }
    }

    if (ref_high_ptr == start_ptr + 1) {
      /* Start 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      diff_64 = clear_start_64(diff_64,startdiscard);
      debug(printf("startdisc:  %016lX",diff_64));
      debug(printf(" => %d leading zeroes\n",count_leading_zeroes_64(diff_64)));

      if (nonzero_p_64(diff_64)) {
	mismatch_position_from_pos3 = (relpos = count_leading_zeroes_64(diff_64)) - offset;
	debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
	return mismatch_position_from_pos3;
      } else {
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }
      
    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,/*genome_unk_mismatch_p*/true);
      diff_32 = clear_start_32(diff_32,startdiscard);
      debug(printf("startdisc:  %08X",diff_32));
      debug(printf(" => %d leading zeroes\n",count_leading_zeroes_32(diff_32)));

      if (nonzero_p_32(diff_32)) {
	mismatch_position_from_pos3 = (relpos = count_leading_zeroes_32(diff_32)) - offset;
	debug(printf("returning %d - %d consecutive matches\n",pos3,pos3 - mismatch_position_from_pos3));
	return mismatch_position_from_pos3;
      } else {
	debug(printf("Would return %d - %d consecutive matches\n",pos3,pos5));
	return (pos3 - pos5);
      }
    }
  }
}


void
Genomebits_consec_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			 Mode_T mode) {

  query_unk_mismatch_p = query_unk_mismatch_p_in;
  genome_unk_mismatch_p = genome_unk_mismatch_p_in;

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

  return;
}


