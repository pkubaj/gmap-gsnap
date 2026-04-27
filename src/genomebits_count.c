static char rcsid[] = "$Id: 60c5707090c1aa0c79a99267135a1d4436bd0ab9 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_count.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* mark_mismatches */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Genomebits_T

static T ref;
static T alt;

static bool md_report_snps_p;
static bool maskedp;
static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = false; /* Needs to be false for path-eval assertions, but needs to be true for circular alignments and mark_mismatches */

static Diffproc_32_T block_diff_32;
static Diffproc_snp_32_T block_diff_snp_32;

static Diffproc_64_T block_diff_64;
static Diffproc_snp_64_T block_diff_snp_64;

#ifdef HAVE_SSE2
static Diffproc_128_T block_diff_128;
static Diffproc_snp_128_T block_diff_snp_128;
#endif


#define T Genomebits_T

static int
count_mismatches_substring (T ref, Compress_T query_compress,
			    Univcoord_T univdiagonal, int querylength,
			    int pos5, int pos3, bool plusp, int genestrand,
			    bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, nshift;
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
	printf("Genome from %u+%d to %u+%d, query_unk %d, genome_unk %d:\n",
	       left,pos5,left,pos3,query_unk_mismatch_p,genome_unk_mismatch_p);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == plusp);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("univdiagonal %u, left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       univdiagonal,left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query (%s) shifted %d:\n",Compress_fwdp(query_compress) == true ? "fwd" : "rev",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32;
  enddiscard = (left+pos3) % 32;
  debug(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X",diff_32));
    debug(printf(" => %d mismatches\n",popcount_ones_32(diff_32)));

    debug(printf("returning %d mismatches\n",popcount_ones_32(diff_32)));
    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:  %016lX",diff_64));
    debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));

    debug(printf("returning %lld mismatches\n",popcount_ones_64(diff_64)));
    return popcount_ones_64(diff_64);

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX",diff_64));

    nmismatches = popcount_ones_64(diff_64);
    debug(printf(" => %d mismatches\n",nmismatches));

    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX",diff_64));

      nmismatches += popcount_ones_64(diff_64);
      debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));
      
      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX",diff_64));
      debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));

      debug(printf("returning %lld mismatches\n",nmismatches + popcount_ones_64(diff_64)));
      return nmismatches + popcount_ones_64(diff_64);

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X",diff_32));
      debug(printf(" => %d mismatches\n",popcount_ones_32(diff_32)));

      debug(printf("returning %d mismatches\n",nmismatches + popcount_ones_32(diff_32)));
      return nmismatches + popcount_ones_32(diff_32);
    }
  }
}

static int
count_mismatches_substring_snps (T ref, T alt, Compress_T query_compress,
				 Univcoord_T univdiagonal, int querylength,
				 int pos5, int pos3, bool plusp, int genestrand,
				 bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, nshift;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  Genomecomp_T *alt_high_ptr, *alt_low_ptr, *alt_flags_ptr;
  UINT4 diff_32;
  UINT8 diff_64;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug(
	printf("\n\n");
	printf("Genome (in count_mismatches_substring_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == plusp);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug(printf("univdiagonal %u, left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
	       univdiagonal,left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug(printf("Query (%s) shifted %d:\n",Compress_fwdp(query_compress) == true ? "fwd" : "rev",nshift));
  debug(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);
  alt_high_ptr = &(alt->high_blocks[startblocki]);
  alt_low_ptr = &(alt->low_blocks[startblocki]);
  alt_flags_ptr = &(alt->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
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

    debug(printf(" => returning %d mismatches\n",popcount_ones_32(diff_32)));
    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);

    debug(printf(" => returning %llu mismatches\n",popcount_ones_64(diff_64)));
    return popcount_ones_64(diff_64);

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    nmismatches = popcount_ones_64(diff_64);
    debug(printf("startdisc: %016lX (%d mismatches)",diff_64,nmismatches));

    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    alt_high_ptr += 2; alt_low_ptr += 2; alt_flags_ptr += 2;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      nmismatches += popcount_ones_64(diff_64);
      debug(printf("  nodisc: %016lX (%lld mismatches)",diff_64,popcount_ones_64(diff_64)));

      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      alt_high_ptr += 2; alt_low_ptr += 2; alt_flags_ptr += 2;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("  enddisc: %016lX (%lld mismatches)",diff_64,popcount_ones_64(diff_64)));

      debug(printf(" => returning %llu mismatches\n",nmismatches + popcount_ones_64(diff_64)));
      return nmismatches + popcount_ones_64(diff_64);
      
    } else {
      /* End 32-bit */
      diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,
				    plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);

      debug(printf(" => returning %d mismatches\n",nmismatches + popcount_ones_32(diff_32)));
      return nmismatches + popcount_ones_32(diff_32);
    }
  }
}


/* left is where the start of the query matches.  pos5 is where we
   want to start comparing in the query.  pos3 is just after where we
   want to stop comparing in the query, i.e., stop at (pos3-1)
   inclusive */
int
Genomebits_count_mismatches_substring (int *ref_mismatches, T ref, T alt, Compress_T query_compress,
				       Univcoord_T univdiagonal, int querylength,
				       int pos5, int pos3, bool plusp, int genestrand) {
  assert(pos5 <= pos3);

  if (alt == NULL) {
    *ref_mismatches = count_mismatches_substring(ref,query_compress,univdiagonal,querylength,
						 pos5,pos3,plusp,genestrand,
						 query_unk_mismatch_p,genome_unk_mismatch_p);
    return *ref_mismatches;

  } else if (maskedp == false) {
    /* Rely solely on SNP-tolerant alignment */
    *ref_mismatches = count_mismatches_substring_snps(ref,alt,query_compress,
						      univdiagonal,querylength,
						      pos5,pos3,plusp,genestrand,
						      query_unk_mismatch_p,genome_unk_mismatch_p);
    return *ref_mismatches;

  } else {
    /* Mask genomic N's, and do not count as mismatches */
    /* Return only mismatches to exons */
    *ref_mismatches = count_mismatches_substring_snps(ref,alt,query_compress,
						      univdiagonal,querylength,
						      pos5,pos3,plusp,genestrand,
						      query_unk_mismatch_p,/*genome_unk_mismatch_p*/false);
    return *ref_mismatches;
  }
}


/************************************************************************
 *  Marking
 ************************************************************************/

static char lowercase_chartable[] = "acgtxxxx";
static char lowercase_revcomp_chartable[] = "tgcaxxxx";

/* Derived from mismatches_left() */
/* blocks can be ref_blocks or alt_blocks (used when --use-mask is specified) */
/* Input genomic sequence originally starts off as query, so procedure
   needs to replace mismatches with ref sequence */
static int
mark_mismatches (char *genomic, T ref, Compress_T query_compress,
		 Univcoord_T univdiagonal, int querylength,
		 int pos5, int pos3, bool segment_plusp, bool query_plusp, int genestrand,
		 bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, mismatch_position, offset, nshift, relpos;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
  int idx;


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug1(
	 printf("\n\n");
	 printf("genomic before mismatches (same as query) = %s\n",genomic);
	 printf("Genome (in mark_mismatches) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	 Genomebits_print(ref,left+pos5,left+pos3);
	 printf("\n");
	 );
  
  assert(Compress_fwdp(query_compress) == segment_plusp);

  startblocki =(left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug1(printf("unidiagonal %u, left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
		univdiagonal,left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug1(printf("Query shifted %d:\n",nshift));
  debug1(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5; /* for mismatch_position from qpos 0 */
  debug1(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug1(printf("diff_32 %08X\n",diff_32));

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      idx = ((*ref_flags_ptr >> relpos) & 0x1) << 2 |
	((*ref_high_ptr >> relpos) & 0x1) << 1 |
	((*ref_low_ptr >> relpos) & 0x1);
      debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    debug1(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug1(printf("diff_64 %016X\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
      idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	((cast64(ref_low_ptr) >> relpos) & 0x1);
      debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    debug1(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
			      query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug1(printf("diff_64 %016X\n",diff_64));

    while (nonzero_p_64(diff_64)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
      idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	((cast64(ref_low_ptr) >> relpos) & 0x1);
      debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    offset += 64;


    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug1(printf("diff_64 %016X\n",diff_64));

      while (nonzero_p_64(diff_64)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
	idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	  ((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	  ((cast64(ref_low_ptr) >> relpos) & 0x1);
	debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
      query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug1(printf("diff_64 %016X\n",diff_64));

      while (nonzero_p_64(diff_64)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
	idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	  ((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	  ((cast64(ref_low_ptr) >> relpos) & 0x1);
	debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
      debug1(printf("genomic = %s\n",genomic));
      return nmismatches;

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_shifted,query_low_shifted,
				query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug1(printf("diff_32 %08X\n",diff_32));

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	idx = ((*ref_flags_ptr >> relpos) & 0x1) << 2 |
	  ((*ref_high_ptr >> relpos) & 0x1) << 1 |
	  ((*ref_low_ptr >> relpos) & 0x1);
	debug1(printf("mismatch position %d, idx %d\n",mismatch_position,idx));
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
      debug1(printf("genomic = %s\n",genomic));
      return nmismatches;
    }
  }
}


/* Derived from mismatches_left_alt() */
/* Returns nmismatches_both */

/* Should be the same as mark_mismatches, except it calls
   block_diff_snp procedures.  Since genomic originally starts off as
   query, it needs to repl ace mismatches with ref sequence */
static int
mark_mismatches_snps (char *genomic, T ref, T alt, Compress_T query_compress,
		      Univcoord_T univdiagonal, int querylength,
		      int pos5, int pos3, bool segment_plusp, bool query_plusp, int genestrand,
		      bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0, mismatch_position, offset, nshift, relpos;
  int startdiscard, enddiscard;
  Univcoord_T left, startblocki, endblocki;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  Genomecomp_T *alt_high_ptr, *alt_low_ptr, *alt_flags_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
  int idx;

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }

  debug1(
	printf("\n\n");
	printf("genomic before mismatches (same as query) = %s\n",genomic);
	printf("Genome (in mark_mismatches_snps) from %u+%d to %u+%d:\n",left,pos5,left,pos3);
	Genomebits_print(ref,left+pos5,left+pos3);
	printf("\n");
	);

  assert(Compress_fwdp(query_compress) == segment_plusp);

  startblocki = (left+pos5)/32U;
  endblocki = (left+pos3)/32U;

  debug1(printf("univdiagonal %u, left = %u, pos5 = %d, pos3 = %d, startblocki = %u, endblocki = %u\n",
		univdiagonal,left,pos5,pos3,startblocki,endblocki));

  nshift = left % 32U;
  Compress_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,query_compress,
		 nshift,/*initpos*/pos5);
  debug1(printf("Query shifted %d:\n",nshift));
  debug1(Compress_print(query_compress,nshift,pos5,pos3));

  ref_high_ptr = &(ref->high_blocks[startblocki]);
  ref_low_ptr = &(ref->low_blocks[startblocki]);
  ref_flags_ptr = &(ref->flags_blocks[startblocki]);
  alt_high_ptr = &(alt->high_blocks[startblocki]);
  alt_low_ptr = &(alt->low_blocks[startblocki]);
  alt_flags_ptr = &(alt->flags_blocks[startblocki]);

  startdiscard = (left+pos5) % 32U;
  enddiscard = (left+pos3) % 32U;
  offset = -startdiscard + pos5; /* For mismatch_position from qpos 0 */
  debug1(printf("nshift = %d, startdiscard = %d, enddiscard = %d or %d\n",
	       nshift,startdiscard,enddiscard,enddiscard + 32));

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);

    while (nonzero_p_32(diff_32)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
      idx = ((*ref_flags_ptr >> relpos) & 0x1) << 2 |
	((*ref_high_ptr >> relpos) & 0x1) << 1 |
	((*ref_low_ptr >> relpos) & 0x1);
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    debug1(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);

    while (nonzero_p_64(diff_64)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
      idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	((cast64(ref_low_ptr) >> relpos) & 0x1);
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    debug1(printf("genomic = %s\n",genomic));
    return nmismatches;

  } else {
    /* Multiple words */
    end_ptr = &(ref->high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				  query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				  alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				  segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);

    while (nonzero_p_64(diff_64)) {
      mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
      idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	((cast64(ref_low_ptr) >> relpos) & 0x1);
      if (segment_plusp != query_plusp) {
	mismatch_position = (querylength - 1) - mismatch_position;
	genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
      } else {
	genomic[mismatch_position] = lowercase_chartable[idx];
      }
      nmismatches++;
    }
    query_high_shifted += 2; query_low_shifted += 2; query_flags_shifted += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    alt_high_ptr += 2; alt_low_ptr += 2; alt_flags_ptr += 2;
    offset += 64;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_snp_64)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				    segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */

      while (nonzero_p_64(diff_64)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
	idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	  ((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	  ((cast64(ref_low_ptr) >> relpos) & 0x1);
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
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
				    segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);

      while (nonzero_p_64(diff_64)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
	idx = ((cast64(ref_flags_ptr) >> relpos) & 0x1) << 2 |
	  ((cast64(ref_high_ptr) >> relpos) & 0x1) << 1 |
	  ((cast64(ref_low_ptr) >> relpos) & 0x1);
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
      return nmismatches;

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_snp_32)(query_high_shifted,query_low_shifted,
				    query_flags_shifted,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				    alt_high_ptr,alt_low_ptr,alt_flags_ptr,			      
				    segment_plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);

      while (nonzero_p_32(diff_32)) {
	mismatch_position = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
	idx = ((*ref_flags_ptr >> relpos) & 0x1) << 2 |
	  ((*ref_high_ptr >> relpos) & 0x1) << 1 |
	  ((*ref_low_ptr >> relpos) & 0x1);
	if (segment_plusp != query_plusp) {
	  mismatch_position = (querylength - 1) - mismatch_position;
	  genomic[mismatch_position] = lowercase_revcomp_chartable[idx];
	} else {
	  genomic[mismatch_position] = lowercase_chartable[idx];
	}
	nmismatches++;
      }
      return nmismatches;
    }
  }
}


int
Genomebits_mark_mismatches (int *nmatches_exonic, char *genomic,
			    Compress_T query_compress,
			    Univcoord_T univdiagonal, int querylength,
			    int pos5, int pos3, bool segment_plusp, bool query_plusp,
			    int genestrand) {
  assert(pos5 < pos3);

  if (alt == NULL) {
    *nmatches_exonic = 0;
    return mark_mismatches(&(*genomic),ref,query_compress,
			   univdiagonal,querylength,pos5,pos3,
			   segment_plusp,query_plusp,genestrand,
			   /*query_unk_mismatch_p*/true,/*genome_unk_mismatch_p*/true);

  } else if (maskedp == true) {
    *nmatches_exonic = (pos3 - pos5) -
      count_mismatches_substring(alt,query_compress,
				 univdiagonal,querylength,
				 pos5,pos3,segment_plusp,genestrand,
				 query_unk_mismatch_p,
				 /*genome_unk_mismatch_p*/true);

    /* Returning and marking only mismatches in exon region */
    return mark_mismatches(&(*genomic),alt,query_compress,
			   univdiagonal,querylength,pos5,pos3,
			   segment_plusp,query_plusp,genestrand,
			   query_unk_mismatch_p,/*genome_unk_mismatch_p*/false);

  } else if (md_report_snps_p == true) {
    /* Mark relative to ref, not using alt */
    *nmatches_exonic = 0;
    return mark_mismatches(&(*genomic),ref,query_compress,
			   univdiagonal,querylength,pos5,pos3,
			   segment_plusp,query_plusp,genestrand,
			   query_unk_mismatch_p,genome_unk_mismatch_p);


  } else {
    /* Mark in with snp-tolerance (ignoring known SNPs) */
    *nmatches_exonic = 0;
    return mark_mismatches_snps(&(*genomic),ref,alt,query_compress,
				univdiagonal,querylength,pos5,pos3,
				segment_plusp,query_plusp,genestrand,
				query_unk_mismatch_p,genome_unk_mismatch_p);
  }
}




void
Genomebits_count_setup (T ref_in, T alt_in,
			bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			Mode_T mode, bool md_report_snps_p_in, bool maskedp_in) {

  ref = ref_in;
  alt = alt_in;

  query_unk_mismatch_p = query_unk_mismatch_p_in;
  genome_unk_mismatch_p = genome_unk_mismatch_p_in;
  md_report_snps_p = md_report_snps_p_in;
  maskedp = maskedp_in;


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



