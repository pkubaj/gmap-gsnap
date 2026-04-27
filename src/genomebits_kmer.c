static char rcsid[] = "$Id: 3b4a4ebf27e89a015546ceb367d23b5a51278f2d $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_kmer.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

static bool genome_unk_mismatch_p = true;

static Diffproc_32_T block_diff_32;
static Diffproc_64_T block_diff_64; 
#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_128; 
#endif


#define T Genomebits_T

/* from left, or a rightward search */
int
Genomebits_first_kmer_left (int *nmismatches5, T ref, Compress_T query_compress,
			    Univcoord_T univdiagonal, int querylength,
			    int pos5, int pos3, bool plusp, int genestrand,
			    bool query_unk_mismatch_p, int kmer) {
  int mismatch_position, last_mismatch, offset, nshift, relpos;
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
	printf("Entered Genome_first_kmer_left with plusp %d, pos5 %d and pos3 %d\n",plusp,pos5,pos3);
	printf("Genome (in first_kmer_left) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(pos5 < pos3);

  *nmismatches5 = -1;

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query (fwdp %d) shifted %d:\n",Compress_fwdp(query_compress),nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  last_mismatch = mismatch_position = pos5 - 1;

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      *nmismatches5 += 1;
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);

    while (nonzero_p_64(diff_64) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      *nmismatches5 += 1;
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      *nmismatches5 += 1;
      return (mismatch_position + 1);
    }

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);

    while (nonzero_p_64(diff_64) && mismatch_position - last_mismatch <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      *nmismatches5 += 1;
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    debug(printf("2 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (mismatch_position - last_mismatch > kmer) {
      return (last_mismatch + 1);
    } else {
      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */

      while (nonzero_p_64(diff_64) && mismatch_position - last_mismatch <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
	*nmismatches5 += 1;
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
      debug(printf("3 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
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
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);

      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      } else {
	*nmismatches5 += 1;
	return (mismatch_position + 1);
      }
      
    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);

      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (mismatch_position - last_mismatch > kmer) {
	return (last_mismatch + 1);
      } else {
	*nmismatches5 += 1;
	return (mismatch_position + 1);
      }
    }
  }
}


/* from right, or a leftward search */
int
Genomebits_first_kmer_right (int *nmismatches3, T ref, Compress_T query_compress,
			     Univcoord_T univdiagonal, int querylength,
			     int pos5, int pos3, bool plusp, int genestrand,
			     bool query_unk_mismatch_p, int kmer) {
  int mismatch_position, last_mismatch, offset, nshift, relpos;
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

  assert(pos5 < pos3);
  debug(
	printf("\n\n");
	printf("Entered Genome_first_kmer_right with plusp %d, pos5 %d and pos3 %d\n",plusp,pos5,pos3);
	printf("Genome (in first_kmer_right) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  *nmismatches3 = -1;

  endblocki = (left+pos3)/32U;
  startblocki = (left+pos5)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos3/*or pos3 - 1?*/);
  debug(printf("Query (fwdp %d) shifted %d:\n",Compress_fwdp(query_compress),nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  /* For leftward scanning */
  ref_high_ptr = &(ref->high_blocks[endblocki]);
  ref_low_ptr = &(ref->low_blocks[endblocki]);
  ref_flags_ptr = &(ref->flags_blocks[endblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = (32 - enddiscard) + pos3 - 1;
  debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

  last_mismatch = mismatch_position = pos3;

  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = clear_start_32(diff_32,startdiscard);

    while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
      *nmismatches3 += 1;
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    startdiscard += 32;
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);

    while (nonzero_p_64(diff_64) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_64(diff_64));
      *nmismatches3 += 1;
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    debug(printf("1 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      *nmismatches3 += 1;
      return mismatch_position;
    }

  } else {
    /* Multiple words */
    start_ptr = &(ref->high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);

    while (nonzero_p_64(diff_64) && last_mismatch - mismatch_position <= kmer) {
      last_mismatch = mismatch_position;
      mismatch_position = offset - (relpos = count_leading_zeroes_64(diff_64));
      *nmismatches3 += 1;
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    debug(printf("2 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
    if (last_mismatch - mismatch_position > kmer) {
      return last_mismatch;
    } else {
      query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
      ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
      offset -= 64;
    }

    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */

      while (nonzero_p_64(diff_64) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_64(diff_64));
	*nmismatches3 += 1;
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      debug(printf("3 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
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
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_start_64(diff_64,startdiscard);

      while (nonzero_p_64(diff_64) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_64(diff_64));
	*nmismatches3 += 1;
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      } else {
	*nmismatches3 += 1;
	return mismatch_position;
      }

    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_start_32(diff_32,startdiscard);

      while (nonzero_p_32(diff_32) && last_mismatch - mismatch_position <= kmer) {
	last_mismatch = mismatch_position;
	mismatch_position = offset - (relpos = count_leading_zeroes_32(diff_32));
	*nmismatches3 += 1;
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      debug(printf("4 mismatch_position %d, last_mismatch %d\n",mismatch_position,last_mismatch));
      if (last_mismatch - mismatch_position > kmer) {
	return last_mismatch;
      } else {
	*nmismatches3 += 1;
	return mismatch_position;
      }
    }
  }
}


void
Genomebits_kmer_setup (bool genome_unk_mismatch_p_in, Mode_T mode) {

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


