static char rcsid[] = "$Id: ee8cb5519252acf26e38798f6d054e210a15df01 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_mismatches.h"
#include "genomebits_decode.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_AVX2
#include <immintrin.h>
#endif
#endif


/* Results */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

static bool query_unk_mismatch_p = false; /* Needs to be false to work with suffix arrays */
static bool genome_unk_mismatch_p = true;

static Diffproc_32_T block_diff_32;
static Diffproc_snp_32_T block_diff_snp_32;

static Diffproc_64_T block_diff_64; 
static Diffproc_snp_64_T block_diff_snp_64; 

#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_128; 
static Diffproc_snp_128_T block_diff_snp_128; 
#endif


#define T Genomebits_T

/* Can return a value in 0..(max_mismatches+1) */
/* mismatch_positions must have (max_mismatches+1) slots available */
static int
mismatches_fromleft (int *mismatch_positions, int max_mismatches, T ref, Compress_T query_compress,
		     Univcoord_T univdiagonal, int querylength,
		     int pos5, int pos3, bool plusp, int genestrand,
		     bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Entered mismatches_fromleft with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_fromleft) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);


  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    return Genomebits_decode_trailing_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
#endif
    return nmismatches;

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    return Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
    if ((nmismatches = Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches)) > max_mismatches) {
      return nmismatches;
    }
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
#endif

    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    offset += 64;


    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));
      
#ifdef HAVE_SSE2
      if ((nmismatches = Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches)) > max_mismatches) {
	return nmismatches;
      }
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
#endif

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
      debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
      return Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
      return nmismatches;
#endif

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X\n",diff_32));

#ifdef HAVE_SSE2
      return Genomebits_decode_trailing_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      return nmismatches;
#endif
    }
  }
}


/* Returns nmismatches_fromboth */
static int
mismatches_fromleft_alt (int *mismatch_positions, int max_mismatches, T ref, T alt,
			 Compress_T query_compress, Univcoord_T univdiagonal, int querylength,
			 int pos5, int pos3, bool plusp, int genestrand,
			 bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  Genomecomp_T *alt_high_ptr, *alt_low_ptr, *alt_flags_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Entered mismatches_fromleft_alt with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_fromleft_alt) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);
  alt_high_ptr = &(alt->high_blocks[startblocki]);
  alt_low_ptr = &(alt->low_blocks[startblocki]);
  alt_flags_ptr = &(alt->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    return Genomebits_decode_trailing_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
#endif
    return nmismatches;

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    return Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
    if ((nmismatches = Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches)) > max_mismatches) {
      return nmismatches;
    }
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
#endif

    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    alt_high_ptr += 2; alt_low_ptr += 2; alt_flags_ptr += 2;
    offset += 64;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));

#ifdef HAVE_SSE2
      if ((nmismatches = Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches)) > max_mismatches) {
	return nmismatches;
      }
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
#endif

      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      alt_high_ptr += 2; alt_low_ptr += 2; alt_flags_ptr += 2;
      offset += 64;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
      return Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
      return nmismatches;
#endif

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X\n",diff_32));

#ifdef HAVE_SSE2
      return Genomebits_decode_trailing_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
      return nmismatches;
#endif
    }
  }
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches */
/* If request max_mismatches 3, could return 0, 1, 2, 3, or 4.  Array
   contains one more element than the return value.  Therefore, need
   to supply array with (max_mismatches + 2) spaces */
int
Genomebits_mismatches_fromleft (int *mismatch_positions, int max_mismatches, T ref, T alt,
				Compress_T query_compress,
				Univcoord_T univdiagonal, int querylength,
				int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG0
  int i;
#endif

  assert(pos5 < pos3);

  if (alt == NULL) {
    nmismatches = mismatches_fromleft(&(*mismatch_positions),max_mismatches,ref,query_compress,
				      univdiagonal,querylength,pos5,pos3,plusp,genestrand,
				      query_unk_mismatch_p,genome_unk_mismatch_p);
  } else {

    nmismatches = mismatches_fromleft_alt(&(*mismatch_positions),max_mismatches,ref,alt,query_compress,
					  univdiagonal,querylength,pos5,pos3,plusp,genestrand,
					  query_unk_mismatch_p,genome_unk_mismatch_p);
  }
  mismatch_positions[nmismatches] = pos3;
  debug0(printf("Genomebits_mismatches_fromleft with left %u, pos5 %d, pos3 %d\n",left,pos5,pos3));
  debug0(
	 printf("%d mismatches on left: ",nmismatches);
	 for (i = 0; i <= nmismatches; i++) {
	   printf("%d ",mismatch_positions[i]);
	 }
	 printf("\n");
	 );
  
  return nmismatches;
}


/* We set query_unk_mismatch_p to false because trimming query N's
   can affect the --clip-overlap feature. */
/* genome_unk_mismatch_p needs to be true so circular alignments around the origin are favored */
int
Genomebits_mismatches_fromleft_for_trim (int *mismatch_positions, int max_mismatches, T ref, T alt,
					 Compress_T query_compress,
					 Univcoord_T univdiagonal, int querylength,
					 int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

  assert(pos5 < pos3);

  if (alt == NULL) {
    nmismatches = mismatches_fromleft(&(*mismatch_positions),max_mismatches,ref,query_compress,
				      univdiagonal,querylength,pos5,pos3,
				      plusp,genestrand,/*query_unk_mismatch_p*/false,
				      /*genome_unk_mismatch_p*/true);
  } else {
    nmismatches = mismatches_fromleft_alt(&(*mismatch_positions),max_mismatches,ref,alt,query_compress,
					  univdiagonal,querylength,pos5,pos3,
					  plusp,genestrand,/*query_unk_mismatch_p*/false,
					  /*genome_unk_mismatch_p*/true);
  }
  mismatch_positions[nmismatches] = pos3;
  debug(
	printf("%d mismatches on left: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  
  return nmismatches;
}


static int
mismatches_fromright (int *mismatch_positions, int max_mismatches, T ref, Compress_T query_compress,
		      Univcoord_T univdiagonal, int querylength,
		      int pos5, int pos3, bool plusp, int genestrand,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *start_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Entered mismatches_fromright with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_fromright) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  endblocki = (left+pos3)/32U;
  startblocki = (left+pos5)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos3/*or pos3 - 1?*/);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  /* For leftward scanning */
  ref_high_ptr = &(ref->high_blocks[endblocki]);
  ref_low_ptr = &(ref->low_blocks[endblocki]);
  ref_flags_ptr = &(ref->flags_blocks[endblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = (32 - enddiscard) + pos3 - 1;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = clear_start_32(diff_32,startdiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    return Genomebits_decode_leading_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    return Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    return nmismatches;
#endif

  } else {
    /* Multiple words */
    start_ptr = &(ref->high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
			      query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
    if ((nmismatches = Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches)) > max_mismatches) {
      return nmismatches;
    }
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
#endif

    query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
    ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
    offset -= 64;


    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_64)(query_high_shifted - 1,query_low_shifted - 1,
				query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));

#ifdef HAVE_SSE2
      if ((nmismatches = Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches)) > max_mismatches) {
	return nmismatches;
      }
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
#endif

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
      debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
      return Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      return nmismatches;
#endif

    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_start_32(diff_32,startdiscard);
      debug(printf("startdisc:  %08X\n",diff_32));

#ifdef HAVE_SSE2
      return Genomebits_decode_leading_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      return nmismatches;
#endif
    }
  }
}


/* Returns nmismatches_fromboth */
static int
mismatches_fromright_alt (int *mismatch_positions, int max_mismatches, T ref, T alt,
			  Compress_T query_compress, Univcoord_T univdiagonal, int querylength,
			  int pos5, int pos3, bool plusp, int genestrand,
			  bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, offset, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *start_ptr;
  Genomecomp_T *alt_high_ptr, *alt_low_ptr, *alt_flags_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Entered mismatches_fromright_alt with %d max_mismatches\n",max_mismatches);
	printf("Genome (in mismatches_fromright_alt) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  endblocki = (left+pos3)/32U;
  startblocki = (left+pos5)/32U;

  debug(printf("left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos3/*or pos3 - 1?*/);
  debug(printf("Query shifted %d:\n",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  /* For leftward scanning */
  ref_high_ptr = &(ref->high_blocks[endblocki]);
  ref_low_ptr = &(ref->low_blocks[endblocki]);
  ref_flags_ptr = &(ref->flags_blocks[endblocki]);
  alt_high_ptr = &(alt->high_blocks[endblocki]);
  alt_low_ptr = &(alt->low_blocks[endblocki]);
  alt_flags_ptr = &(alt->flags_blocks[endblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = (32 - enddiscard) + pos3 - 1;
  debug(printf("nshift = %d, startdiscard = %u, enddiscard = %u\n",nshift,startdiscard,enddiscard));

  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    return Genomebits_decode_leading_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
    return nmismatches;
#endif

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_snp_64)(query_high_shifted - 1,query_low_shifted - 1,
				  query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				  alt_high_ptr - 1,alt_low_ptr - 1,alt_flags_ptr - 1,
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    return Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    return nmismatches;
#endif
    
  } else {
    /* Multiple words */
    start_ptr = &(ref->high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_snp_64)(query_high_shifted - 1,query_low_shifted - 1,
				  query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				  alt_high_ptr - 1,alt_low_ptr - 1,alt_flags_ptr - 1,
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
    if ((nmismatches = Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches)) > max_mismatches) {
      return nmismatches;
    }
#else
    while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
    if (nmismatches > max_mismatches) {
      return nmismatches;
    }
#endif

    query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
    ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
    alt_high_ptr -= 2; alt_low_ptr -= 2; alt_flags_ptr -= 2;
    offset -= 64;


    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_snp_64)(query_high_shifted - 1,query_low_shifted - 1,
				    query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				    alt_high_ptr - 1,alt_low_ptr - 1,alt_flags_ptr - 1,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));

#ifdef HAVE_SSE2
      if ((nmismatches = Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches)) > max_mismatches) {
	return nmismatches;
      }
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      if (nmismatches > max_mismatches) {
	return nmismatches;
      }
#endif

      query_high_shifted -= 2; query_low_shifted -= 2; query_flags_shifted -= 2;
      ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
      alt_high_ptr -= 2; alt_low_ptr -= 2; alt_flags_ptr -= 2;
      offset -= 64;
    }

    if (ref_high_ptr == start_ptr + 1) {
      /* Start 64-bit */
      diff_64 = (block_diff_snp_64)(query_high_shifted - 1,query_low_shifted - 1,
				    query_flags_shifted - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				    alt_high_ptr - 1,alt_low_ptr - 1,alt_flags_ptr - 1,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_start_64(diff_64,startdiscard);
      debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
      return Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
      return nmismatches;
#endif

    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_start_32(diff_32,startdiscard);
      debug(printf("startdisc:  %08X\n",diff_32));

#ifdef HAVE_SSE2
      return Genomebits_decode_leading_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches <= max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
      return nmismatches;
#endif
    }
  }
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches */
int
Genomebits_mismatches_fromright (int *mismatch_positions, int max_mismatches, T ref, T alt,
				 Compress_T query_compress,
				 Univcoord_T univdiagonal, int querylength,
				 int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG0
  int i;
#endif

  assert(pos5 < pos3);

  if (alt == NULL) {
    nmismatches = mismatches_fromright(&(*mismatch_positions),max_mismatches,ref,query_compress,
				       univdiagonal,querylength,pos5,pos3,plusp,genestrand,
				       query_unk_mismatch_p,genome_unk_mismatch_p);
  } else {
    nmismatches = mismatches_fromright_alt(&(*mismatch_positions),max_mismatches,ref,alt,query_compress,
					   univdiagonal,querylength,pos5,pos3,plusp,genestrand,
					   query_unk_mismatch_p,genome_unk_mismatch_p);
  }
  mismatch_positions[nmismatches] = pos5 - 1;
  debug0(printf("Genomebits_mismatches_fromleft with left %u, pos5 %d, pos3 %d\n",left,pos5,pos3));
  debug0(
	 printf("%d mismatches on right: ",nmismatches);
	 for (i = 0; i <= nmismatches; i++) {
	   printf("%d ",mismatch_positions[i]);
	 }
	 printf("\n");
	 );
  return nmismatches;
}


/* Returns mismatch_positions[0..nmismatches], where nmismatches <= max_mismatches */
/* We set query_unk_mismatch_p to false because trimming query N's
   can affect the --clip-overlap feature. */
/* genome_unk_mismatch_p needs to be true so circular alignments around the origin are favored */
int
Genomebits_mismatches_fromright_for_trim (int *mismatch_positions, int max_mismatches, T ref, T alt,
					  Compress_T query_compress,
					  Univcoord_T univdiagonal, int querylength,
					  int pos5, int pos3, bool plusp, int genestrand) {
  int nmismatches;
#ifdef DEBUG
  int i;
#endif

  assert(pos5 < pos3);

  if (alt == NULL) {
    nmismatches = mismatches_fromright(&(*mismatch_positions),max_mismatches,ref,query_compress,
				       univdiagonal,querylength,pos5,pos3,
				       plusp,genestrand,/*query_unk_mismatch_p*/false,
				       /*genome_unk_mismatch_p*/true);
  } else {
    nmismatches = mismatches_fromright_alt(&(*mismatch_positions),max_mismatches,ref,alt,query_compress,
					   univdiagonal,querylength,pos5,pos3,
					   plusp,genestrand,/*query_unk_mismatch_p*/false,
					   /*genome_unk_mismatch_p*/true);
  }
  mismatch_positions[nmismatches] = pos5 - 1;
  debug(
	printf("%d mismatches on right: ",nmismatches);
	for (i = 0; i <= nmismatches; i++) {
	  printf("%d ",mismatch_positions[i]);
	}
	printf("\n");
	);
  return nmismatches;
}

void
Genomebits_mismatches_setup (bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			     Mode_T mode, bool maskedp) {

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



