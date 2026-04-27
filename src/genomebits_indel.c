static char rcsid[] = "$Id: e7c313e61f1e8e849e04c050671254c6f1981c3c $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "genomebits_indel.h"
#include "genomebits_decode.h"

#include <stdio.h>

#include "assert.h"
#include "except.h"


/* Taken from substring.c */
/* Designed to allow 1 match to offset 1 mismatch.  To handle 2 matches vs 2 mismatches, penalize for multiple mismatches */
#define TRIM_MATCH_SCORE 1
#define TRIM_MISMATCH_SCORE_LAST -1 /* Requires 1 match to compensate */
#define TRIM_MISMATCH_SCORE_MULT -4 /* Requires 4 matches to compensate */
#define MIN_INITIAL_SCORE 5
#define MIN_OVERALL_SCORE 2


/* Requires DEBUG1 in compress.c */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Trimming at ends */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


#define T Genomebits_T

static int max_insertionlen;
static int max_deletionlen;


static bool query_unk_mismatch_p = false;
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


#ifdef DEBUG
static void
write_chars (Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags) {
  char Buffer[33];
  int i;

  Buffer[32] = '\0';
  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < 32; i++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }

  if (flags != 0U) {
    for (i = 0; i < 32; i++) {
      if (flags & 0x01) {
	Buffer[i] = 'N';
      }
      flags >>= 1;
    }
  }

  printf("%s",Buffer);
  return;
}
#endif


#ifdef DEBUG
static void
shifted_print (Genomecomp_T *high_ptr, Genomecomp_T *low_ptr, Genomecomp_T *flags_ptr,
	       int startblocki, int endblocki) {
  Genomecomp_T high, low, flags;
  int blocki;

  for (blocki = startblocki; blocki <= endblocki; blocki++) {
    high = high_ptr[blocki]; low = low_ptr[blocki]; flags = flags_ptr[blocki];
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");
  }

  printf("\n");
  return;
}
#endif


/* For deletions, let L = pos3 - pos5, be the length of the fragment
   being tested.  Then the probability of matching L nucleotides
   exactly is $p = (1/4)^L$.  If we test N times, then the probability
   of seeing a hit is $(1 - (1-p)^N$.  For a threshold of 0.01, we
   compute the allowable deletion lengths as follows:

   for (L in 1:10) {
     p <- (0.25)^L
     N <- 0
     while (N + 1 <= 100 && 1 - (1 - p)^(N + 1) < 0.01) {
       N <- N + 1
     }
     cat(L,N,"\n")
   }
     
   1 0
   2 0
   3 0
   4 2
   5 10
   6 41
   7 100

   If we allow one mismatch, the probability of seeing (L - 1) matches
   and one mismatch is (L choose 1)*4^(1)*1^(L - 1)/4^L =
   L*(0.25)^(L-1).  Performing the a similar computation, we get:

   1 0
   2 0
   3 0
   4 0
   5 0
   6 1
   7 5
   8 20
   9 73
   10 100
*/


static int
compute_max_deletionlen_exact (int fragment_length) {
  int deletionlen;

  if (fragment_length <= 3) {
    deletionlen =  0;
  } else if (fragment_length == 4) {
    deletionlen =  2;
  } else if (fragment_length == 5) {
    deletionlen =  10;
  } else {
    deletionlen =  30;
  }

  if (deletionlen > max_deletionlen) {
    return max_deletionlen;
  } else {
    return deletionlen;
  }
}

static int
compute_max_deletionlen_one_mm (int fragment_length) {
  int deletionlen;

  if (fragment_length <= 5) {
    deletionlen = 0;
  } else if (fragment_length == 6) {
    deletionlen = 1;
  } else if (fragment_length == 7) {
    deletionlen = 5;
  } else if (fragment_length == 8) {
    deletionlen = 20;
  } else {
    deletionlen = 30;
  }

  if (deletionlen > max_deletionlen) {
    return max_deletionlen;
  } else {
    return deletionlen;
  }
}


/* For insertions, let L = pos3 - pos5, be the length of the fragment
   being tested, and consider an insertion of length I, which means
   that I nucleotides are ignored and (L-I) are allowed to match.  The
   probability of those (L-I) nucleotides matching exactly is $p(I) =
   (1/4)^(L-I)$.  If we test N times, then the probability of seeing a
   hit is $1 - \PI_{I=1}^N (1 - p(I))$.  For a threshold of 0.01, we
   compute the allowable insertion lengths as follows:

   prob.insertion <- function (L, N) {
       p <- 1
       for (I in 1:N) {
           p <- p * (1 - (0.25)^(L - I))
       }
       1 - p
   }

   for (L in 1:20) {
       N <- 0
       while (N + 1 < L && prob.insertion(L,N+1) < 0.01) {
          N <- N + 1
       }
       cat(L,N,"\n")
   }

   1 0 
   2 0 
   3 0 
   4 0 
   5 1 
   6 2 
   7 3 
   8 4 
   9 5 
   10 6 
   11 7 
   12 8 
   13 9 
   14 10 
   15 11 
   16 12 
   17 13 
   18 14 
   19 15 
   20 16 

   which shows that we can allow up to an insertion of (L - 4).  If we
   allow 1 mismatch, then we get

   prob.insertion <- function (L, N) {
       p <- 1
       for (I in 1:N) {
           p <- p * (1 - (L - I)*(0.25)^(L - I - 1))
       }
       1 - p
   }

   1 0 
   2 0 
   3 0 
   4 0 
   5 0 
   6 0 
   7 1 
   8 2 
   9 3 
   10 4 
   11 5 
   12 6 
   13 7 
   14 8 
   15 9 
   16 10 
   17 11 
   18 12 
   19 13 
   20 14 

   or essentially an insertion of (L - 6).

*/


static inline int
compute_max_insertionlen_exact (int fragment_length) {
  int insertionlen = fragment_length - 4;

  if (insertionlen > max_insertionlen) {
    return max_insertionlen;
  } else {
    return insertionlen;
  }
}

static inline int
compute_max_insertionlen_one_mm (int fragment_length) {
  int insertionlen = fragment_length - 6;

  if (insertionlen > max_insertionlen) {
    return max_insertionlen;
  } else {
    return insertionlen;
  }
}



/* Modified from Compress_shift in compress.c */
/* Shifts region startpos..endpos to be at bit 0 */
static int
genomebits_shift (Genomecomp_T **genome_high_shifted, Genomecomp_T **genome_low_shifted,
		  Genomecomp_T **genome_flags_shifted, T genome,
		  Univcoord_T startpos, Univcoord_T endpos) {
  int nwords, wordi;
  int leftshift, rightshift;
  Genomecomp_T *high_ptr, *low_ptr, *flags_ptr;
  Univcoord_T startblocki, endblocki, blocki;

  startblocki = startpos/32U;
  endblocki = endpos/32U;
  nwords = endblocki - startblocki + 1;

  debug(printf("genomebits_shift called with startpos %u and endpos %u => %d words\n",
	       startpos,endpos,nwords));

  *genome_high_shifted = (Genomecomp_T *) MALLOC(nwords*sizeof(Genomecomp_T));
  *genome_low_shifted = (Genomecomp_T *) MALLOC(nwords*sizeof(Genomecomp_T));
  *genome_flags_shifted = (Genomecomp_T *) MALLOC(nwords*sizeof(Genomecomp_T));

  high_ptr = &((*genome_high_shifted)[0]);
  low_ptr = &((*genome_low_shifted)[0]);
  flags_ptr = &((*genome_flags_shifted)[0]);
  
  rightshift = startpos % 32;
  leftshift = 32 - rightshift;

  debug(printf("startblocki %d, endblocki %d\n",startblocki,endblocki));

  /* Shift genome leftward (by rightshift) */
  blocki = startblocki;
  if (rightshift == 0) {
    for (wordi = 0; wordi < nwords; wordi++) {
      *high_ptr++ = genome->high_blocks[blocki];
      *low_ptr++ = genome->low_blocks[blocki];
      *flags_ptr++ = genome->flags_blocks[blocki];
      blocki++;
    }
  } else {
    for (wordi = 0; wordi < nwords; wordi++) {
      *high_ptr++ = (genome->high_blocks[blocki] >> rightshift) | (genome->high_blocks[blocki+1] << leftshift);
      *low_ptr++ = (genome->low_blocks[blocki] >> rightshift) | (genome->low_blocks[blocki+1] << leftshift);
      *flags_ptr++ = (genome->flags_blocks[blocki] >> rightshift) | (genome->flags_blocks[blocki+1] << leftshift);
      blocki++;
    }
  }

  return nwords;
}


static int
querybits_shift (Genomecomp_T **query_high_shifted, Genomecomp_T **query_low_shifted,
		 Genomecomp_T **query_flags_shifted, Genomecomp_T *query_high_blocks,
		 Genomecomp_T *query_low_blocks, Genomecomp_T *query_flags_blocks,
		 int pos5, int pos3) {
  int nwords, wordi, nbits;
  int leftshift, rightshift;
  Genomecomp_T *high_ptr, *low_ptr, *flags_ptr;
  int startblocki, blocki;

  nbits = pos3 - pos5;
  nwords = (nbits + 31)/32U;

  debug(printf("querybits_shift called with pos5 %d and pos3 %d => %d words\n",
	       pos5,pos3,nwords));

  /* Include an extra word at the end so we can use it in leftshift/rightshift */
  *query_high_shifted = (Genomecomp_T *) MALLOC((nwords + 1)*sizeof(Genomecomp_T));
  *query_low_shifted = (Genomecomp_T *) MALLOC((nwords + 1)*sizeof(Genomecomp_T));
  *query_flags_shifted = (Genomecomp_T *) MALLOC((nwords + 1)*sizeof(Genomecomp_T));

  high_ptr = &((*query_high_shifted)[0]);
  low_ptr = &((*query_low_shifted)[0]);
  flags_ptr = &((*query_flags_shifted)[0]);
  
  rightshift = pos5 % 32;
  leftshift = 32 - rightshift;

  startblocki = pos5/32;
  debug(int endblocki = pos3/32);

  debug(printf("startblocki %d, endblocki %d\n",startblocki,endblocki));

  /* Shift genome leftward (by rightshift) */
  blocki = startblocki;
  if (rightshift == 0) {
    for (wordi = 0; wordi < nwords; wordi++) {
      *high_ptr++ = query_high_blocks[blocki];
      *low_ptr++ = query_low_blocks[blocki];
      *flags_ptr++ = query_flags_blocks[blocki];
      blocki++;
    }

  } else {
    for (wordi = 0; wordi < nwords; wordi++) {
      *high_ptr++ = (query_high_blocks[blocki] >> rightshift) | (query_high_blocks[blocki+1] << leftshift);
      *low_ptr++ = (query_low_blocks[blocki] >> rightshift) | (query_low_blocks[blocki+1] << leftshift);
      *flags_ptr++ = (query_flags_blocks[blocki] >> rightshift) | (query_flags_blocks[blocki+1] << leftshift);
      blocki++;
    }
  }

  /* Fill last word with Ns */
  *high_ptr = 0U;
  *low_ptr = 0U;
  *flags_ptr = ~0U;

  return nwords;
}


/* Modified from count_mismatches_substring in genomebits_count.c */
static int
count_mismatches_substring (Genomecomp_T *query_high_adj, Genomecomp_T *query_low_adj,
			    Genomecomp_T *query_flags_adj, Genomecomp_T *ref_high_blocks,
			    Genomecomp_T *ref_low_blocks, Genomecomp_T *ref_flags_blocks,
			    int startblocki, int endblocki, int startdiscard, int enddiscard,
			    bool plusp, int genestrand,
			    bool query_unk_mismatch_p, bool genome_unk_mismatch_p) {
  int nmismatches = 0;
  Genomecomp_T *query_high_ptr, *query_low_ptr, *query_flags_ptr;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr;
  Genomecomp_T *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;


  debug(printf("Entering count_mismatches_substring with startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
	       startblocki,endblocki,startdiscard,enddiscard));

  query_high_ptr = &(query_high_adj[startblocki]);
  query_low_ptr = &(query_low_adj[startblocki]);
  query_flags_ptr = &(query_flags_adj[startblocki]);

  ref_high_ptr = &(ref_high_blocks[startblocki]);
  ref_low_ptr = &(ref_low_blocks[startblocki]);
  ref_flags_ptr = &(ref_flags_blocks[startblocki]);

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X",diff_32));
    debug(printf(" => %d mismatches\n",popcount_ones_32(diff_32)));

    debug(printf("returning %d mismatches\n",popcount_ones_32(diff_32)));
    return popcount_ones_32(diff_32);

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:  %016lX",diff_64));
    debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));

    debug(printf("returning %lld mismatches\n",popcount_ones_64(diff_64)));
    return popcount_ones_64(diff_64);

  } else {
    /* Multiple words */
    end_ptr = &(ref_high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX",diff_64));

    nmismatches = popcount_ones_64(diff_64);
    debug(printf(" => %d mismatches\n",nmismatches));

    query_high_ptr += 2; query_low_ptr += 2; query_flags_ptr += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;

    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX",diff_64));

      nmismatches += popcount_ones_64(diff_64);
      debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));
      
      query_high_ptr += 2; query_low_ptr += 2; query_flags_ptr += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX",diff_64));
      debug(printf(" => %lld mismatches\n",popcount_ones_64(diff_64)));

      debug(printf("returning %lld mismatches\n",nmismatches + popcount_ones_64(diff_64)));
      return nmismatches + popcount_ones_64(diff_64);

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
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
solve_high_short (int *nmismatches_to_trimpos, Univcoord_T univdiagonal, int querylength,
		  int pos5, int pos3, Compress_T query_compress,
		  Genomebits_T omebits, bool plusp, int genestrand) {
  int best_adj = 0;
  int base_nmismatches, nmismatches;
  int adj;
  int max_insertionlen, max_deletionlen;
  int genome_nwords, wordi;
  int leftshift, rightshift;
  Univcoord_T left, startpos, endpos;
  Genomecomp_T *ref_high_shifted, *ref_low_shifted, *ref_flags_shifted;
  Genomecomp_T *query_high_blocks, *query_low_blocks, *query_flags_blocks;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *query_high_adj, *query_low_adj, *query_flags_adj;
  int startblocki, endblocki, blocki;
  int startdiscard, enddiscard;

  debug(printf("solve_high_short called with univdiagonal %u, pos %d..%d\n",
	       univdiagonal,pos5,pos3));

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }
  max_deletionlen = compute_max_deletionlen_exact(pos3 - pos5);

  max_insertionlen = compute_max_insertionlen_exact(pos3 - pos5);
  
  if (max_insertionlen <= 0 && max_deletionlen <= 0) {
    return 0;
  } else {
    startpos = left + pos5;
    endpos = left + pos3 + max_deletionlen;
    genome_nwords = genomebits_shift(&ref_high_shifted,&ref_low_shifted,&ref_flags_shifted,
				     /*genome*/omebits,startpos,endpos);
    debug(printf("Genome:\n"));
    debug(shifted_print(ref_high_shifted,ref_low_shifted,ref_flags_shifted,
			/*startblocki*/0,/*endblocki*/genome_nwords - 1));
  }

  
  /* Normally we shift query right, but here we shift it leftward */
  /* Get baseline */
  Compress_shift(&query_high_blocks,&query_low_blocks,&query_flags_blocks,query_compress,
		 /*nshift*/0,/*initpos*/0);
  debug(printf("Query baseline:\n"));
  debug(Compress_print(query_compress,/*nshift*/0,pos5,pos3));
#ifdef DEBUG
  int query_nwords =
#endif
    querybits_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,
		    query_high_blocks,query_low_blocks,query_flags_blocks,
		    pos5,pos3);
  debug(printf("Query (%d words):\n",query_nwords));
  debug(shifted_print(query_high_shifted,query_low_shifted,query_flags_shifted,
		      /*startblocki*/0,/*endblocki*/query_nwords - 1));
  
  query_high_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_low_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_flags_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));


  /* No indel.  No shift */
  startblocki = 0;
  startdiscard = 0;
  endblocki = ((pos3 - pos5) - /*adj*/0)/32;
  enddiscard = ((pos3 - pos5) - /*adj*/0) % 32;

  /* base_nmismatches should be from pos5 to pos3, essentially continuation_nmismatches */
  *nmismatches_to_trimpos = base_nmismatches =
    count_mismatches_substring(query_high_shifted,query_low_shifted,
			       query_flags_shifted,ref_high_shifted,
			       ref_low_shifted,ref_flags_shifted,
			       startblocki,endblocki,startdiscard,enddiscard,
			       plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
  debug(printf("No indel, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
	       startblocki,endblocki,startdiscard,enddiscard,base_nmismatches));
  debug(shifted_print(query_high_shifted,query_low_shifted,query_flags_shifted,startblocki,endblocki));


  /* Deletions */
  for (adj = 1; adj <= max_deletionlen; adj++) {
    startblocki = adj/32;
    startdiscard = adj % 32;
    endblocki = (adj + (pos3 - pos5))/32;
    enddiscard = (adj + (pos3 - pos5)) % 32;
    debug(printf(">adj %d/%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 adj,max_deletionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = adj % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }

    if ((nmismatches = count_mismatches_substring(query_high_adj,query_low_adj,
						  query_flags_adj,ref_high_shifted,
						  ref_low_shifted,ref_flags_shifted,
						  startblocki,endblocki,startdiscard,enddiscard,
						  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p)) < *nmismatches_to_trimpos) {
      *nmismatches_to_trimpos = nmismatches;
      best_adj = +adj;
    }
    debug(printf("adj %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
		 adj,startblocki,endblocki,startdiscard,enddiscard,nmismatches));
    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  if (best_adj == 0) {
    /* Nothing found */
  } else if (*nmismatches_to_trimpos > 1) {
    /* Too many mismatches, so ignore */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else if (*nmismatches_to_trimpos == 0) {
    /* Accept */
  } else if (best_adj > compute_max_deletionlen_one_mm(pos3 - pos5)) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else {
    /* Accept */
  }


  /* Insertions */
  startblocki = 0;
  startdiscard = 0;
  for (adj = 1; adj <= max_insertionlen; adj++) {
    endblocki = ((pos3 - pos5) - adj)/32;
    enddiscard = ((pos3 - pos5) - adj) % 32;
    debug(printf(">adj %d/%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 -adj,-max_insertionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* This is the only case where we shift leftward */
    /* Shift query fragment leftward (by rightshift) */
    wordi = 0;
    if ((rightshift = adj % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      leftshift = 32 - rightshift;
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing rightshift of %d on wordi %d and leftshift of %d on wordi %d\n",
		     blocki,rightshift,wordi,leftshift,wordi+1));
	query_high_adj[blocki] = (query_high_shifted[wordi] >> rightshift) | (query_high_shifted[wordi+1] << leftshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] >> rightshift) | (query_low_shifted[wordi+1] << leftshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] >> rightshift) | (query_flags_shifted[wordi+1] << leftshift);
	wordi++;
      }
    }
    
    if ((nmismatches = count_mismatches_substring(query_high_adj,query_low_adj,
						  query_flags_adj,ref_high_shifted,
						  ref_low_shifted,ref_flags_shifted,
						  startblocki,endblocki,startdiscard,enddiscard,
						  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p)) < *nmismatches_to_trimpos) {
      *nmismatches_to_trimpos = nmismatches;
      best_adj = -adj;
    }
    debug(printf("adj %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
		 -adj,startblocki,endblocki,startdiscard,enddiscard,nmismatches));
    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }


  if (best_adj >= 0) {
    /* No insertion found */
  } else if (*nmismatches_to_trimpos > 1) {
    /* Too many mismatches, so ignore */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else if (*nmismatches_to_trimpos == 0) {
    /* Accept */
  } else if (-best_adj > compute_max_insertionlen_one_mm(pos3 - pos5)) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else {
    /* Accept */
  }


  FREE(query_flags_adj);
  FREE(query_low_adj);
  FREE(query_high_adj);

  FREE(query_flags_shifted);
  FREE(query_low_shifted);
  FREE(query_high_shifted);

  FREE(ref_flags_shifted);
  FREE(ref_low_shifted);
  FREE(ref_high_shifted);

  return best_adj;
}


static int
solve_low_short (int *nmismatches_to_trimpos, Univcoord_T univdiagonal, int querylength,
		 int pos5, int pos3, Compress_T query_compress,
		 Genomebits_T omebits, bool plusp, int genestrand) {
  int best_adj = 0;
  int base_nmismatches, nmismatches;
  int adj;
  int max_insertionlen, max_deletionlen;
  int genome_nwords, wordi;
  int leftshift, rightshift;
  Univcoord_T left, startpos, endpos;
  Genomecomp_T *ref_high_shifted, *ref_low_shifted, *ref_flags_shifted;
  Genomecomp_T *query_high_blocks, *query_low_blocks, *query_flags_blocks;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *query_high_adj, *query_low_adj, *query_flags_adj;
  int startblocki, endblocki, blocki;
  int startdiscard, enddiscard;


  debug(printf("solve_low_short called with univdiagonal %u, pos %d..%d\n",
	       univdiagonal,pos5,pos3));


  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
    max_deletionlen = 0;
  } else {
    max_deletionlen = compute_max_deletionlen_exact(pos3 - pos5);
    if (univdiagonal + pos5 < (Univcoord_T) (querylength + max_deletionlen)) {
      max_deletionlen = 0;
    }
  }

  max_insertionlen = compute_max_insertionlen_exact(pos3 - pos5);
  
  if (max_insertionlen <= 0 && max_deletionlen <= 0) {
    return 0;
  } else {
    startpos = left + pos5 - max_deletionlen;
    endpos = left + pos3;
    genome_nwords = genomebits_shift(&ref_high_shifted,&ref_low_shifted,&ref_flags_shifted,
				     /*genome*/omebits,startpos,endpos);
    debug(printf("Genome:\n"));
    debug(shifted_print(ref_high_shifted,ref_low_shifted,ref_flags_shifted,
			/*startblocki*/0,/*endblocki*/genome_nwords - 1));
  }

  
  /* Normally we shift query right, but here we shift it leftward */
  /* Get baseline */
  Compress_shift(&query_high_blocks,&query_low_blocks,&query_flags_blocks,query_compress,
		 /*nshift*/0,/*initpos*/0);
  debug(printf("Query baseline:\n"));
  debug(Compress_print(query_compress,/*nshift*/0,pos5,pos3));
#ifdef DEBUG
  int query_nwords =
#endif
    querybits_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,
		    query_high_blocks,query_low_blocks,query_flags_blocks,
		    pos5,pos3);
  debug(printf("Query (%d words):\n",query_nwords));
  debug(shifted_print(query_high_shifted,query_low_shifted,query_flags_shifted,
		      /*startblocki*/0,/*endblocki*/query_nwords - 1));

  query_high_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_low_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_flags_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));


  /* No indel.  No shift. */
  startblocki = (max_deletionlen + /*adj*/0) / 32;
  startdiscard = (max_deletionlen + /*adj*/0) % 32;
  endblocki = (max_deletionlen + (pos3 - pos5))/32;
  enddiscard = (max_deletionlen + (pos3 - pos5)) % 32;

  /* Shift query fragment rightward (by leftshift) */
  wordi = 0;
  if ((leftshift = (max_deletionlen + /*adj*/0) % 32) == 0) {
    for (blocki = startblocki; blocki <= endblocki; blocki++) {
      debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi];
      query_low_adj[blocki] = query_low_shifted[wordi];
      query_flags_adj[blocki] = query_flags_shifted[wordi];
      wordi++;
    }

  } else {
    rightshift = 32 - leftshift;
    blocki = startblocki;
    debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		 blocki,leftshift,wordi));
    query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
    query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
    query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
    wordi++;

    for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
      debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		   blocki,leftshift,wordi,rightshift,wordi-1));
      query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
      query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
      query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
      wordi++;
    }
  }

  *nmismatches_to_trimpos = base_nmismatches =
    count_mismatches_substring(query_high_adj,query_low_adj,
			       query_flags_adj,ref_high_shifted,
			       ref_low_shifted,ref_flags_shifted,
			       startblocki,endblocki,startdiscard,enddiscard,
			       plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
  debug(printf("No indel, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
	       startblocki,endblocki,startdiscard,enddiscard,base_nmismatches));
  debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));


  /* Deletions */
  for (adj = 1; adj <= max_deletionlen; adj++) {
    startblocki = (max_deletionlen - adj) / 32;
    startdiscard = (max_deletionlen - adj) % 32;
    endblocki = (max_deletionlen - adj + (pos3 - pos5)) / 32;
    enddiscard = (max_deletionlen - adj + (pos3 - pos5)) % 32;
    debug(printf(">adj %d/%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 adj,max_deletionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = (max_deletionlen - adj) % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }
      
    if ((nmismatches = count_mismatches_substring(query_high_adj,query_low_adj,
						  query_flags_adj,ref_high_shifted,
						  ref_low_shifted,ref_flags_shifted,
						  startblocki,endblocki,startdiscard,enddiscard,
						  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p)) < *nmismatches_to_trimpos) {
      *nmismatches_to_trimpos = nmismatches;
      best_adj = +adj;
    }
    debug(printf("adj %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
		 adj,startblocki,endblocki,startdiscard,enddiscard,nmismatches));
    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  if (best_adj == 0) {
    /* Nothing found */
  } else if (*nmismatches_to_trimpos > 1) {
    /* Too many mismatches, so ignore */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else if (*nmismatches_to_trimpos == 0) {
    /* Accept */
  } else if (best_adj > compute_max_deletionlen_one_mm(pos3 - pos5)) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else {
    /* Accept */
  }


  /* Insertions */
  endblocki = (max_deletionlen + (pos3 - pos5))/32;
  enddiscard = (max_deletionlen + (pos3 - pos5)) % 32;
  for (adj = 1; adj <= max_insertionlen; adj++) {
    startblocki = (max_deletionlen + adj) / 32;
    startdiscard = (max_deletionlen + adj) % 32;
    debug(printf(">adj %d/%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 -adj,-max_insertionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = (max_deletionlen + adj) % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }

    if ((nmismatches = count_mismatches_substring(query_high_adj,query_low_adj,
						  query_flags_adj,ref_high_shifted,
						  ref_low_shifted,ref_flags_shifted,
						  startblocki,endblocki,startdiscard,enddiscard,
						  plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p)) < *nmismatches_to_trimpos) {
      *nmismatches_to_trimpos = nmismatches;
      best_adj = -adj;
    }
    debug(printf("adj %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, nmismatches %d\n",
		 -adj,startblocki,endblocki,startdiscard,enddiscard,nmismatches));
    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  if (best_adj >= 0) {
    /* No insertion found */
  } else if (*nmismatches_to_trimpos > 1) {
    /* Too many mismatches, so ignore */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else if (*nmismatches_to_trimpos == 0) {
    /* Accept */
  } else if (-best_adj > compute_max_insertionlen_one_mm(pos3 - pos5)) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = base_nmismatches;
  } else {
    /* Accept */
  }


  FREE(query_flags_adj);
  FREE(query_low_adj);
  FREE(query_high_adj);

  FREE(query_flags_shifted);
  FREE(query_low_shifted);
  FREE(query_high_shifted);

  FREE(ref_flags_shifted);
  FREE(ref_low_shifted);
  FREE(ref_high_shifted);

  return best_adj;
}


/************************************************************************
 * High or qend indels, with mismatches starting on left
 ************************************************************************/

/* Modified from genomebits_mismatches.c */
static int
mismatches_fromleft (int *mismatch_positions, int max_mismatches, int rightmost_value,
		     Genomecomp_T *query_high_adj, Genomecomp_T *query_low_adj,
		     Genomecomp_T *query_flags_adj, Genomecomp_T *ref_high_blocks,
		     Genomecomp_T *ref_low_blocks, Genomecomp_T *ref_flags_blocks,
		     int startblocki, int endblocki, int startdiscard, int enddiscard,
		     bool plusp, int genestrand) {
  int nmismatches = 0, offset;
  Genomecomp_T *query_high_ptr, *query_low_ptr, *query_flags_ptr;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *end_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif

  /* Matching the parameters for Genomebits_mismatches_left_trim */
  bool query_unk_mismatch_p = false;
  bool genome_unk_mismatch_p = false;

  debug(printf("Entering mismatches_fromleft with startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, max_mismatches %d\n",
	       startblocki,endblocki,startdiscard,enddiscard,max_mismatches));

  query_high_ptr = &(query_high_adj[startblocki]);
  query_low_ptr = &(query_low_adj[startblocki]);
  query_flags_ptr = &(query_flags_adj[startblocki]);

  ref_high_ptr = &(ref_high_blocks[startblocki]);
  ref_low_ptr = &(ref_low_blocks[startblocki]);
  ref_flags_ptr = &(ref_flags_blocks[startblocki]);

  offset = -startdiscard;

  if (endblocki == startblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_32 = clear_start_32(diff_32,startdiscard);
    diff_32 = clear_end_32(diff_32,enddiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_trailing_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
      diff_32 = clear_lowbit_32(diff_32,relpos);
    }
#endif
    mismatch_positions[nmismatches] = rightmost_value;
    return nmismatches;

  } else if (endblocki == startblocki + 1) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
#endif
    mismatch_positions[nmismatches] = rightmost_value;
    return nmismatches;

  } else {
    /* Multiple words */
    end_ptr = &(ref_high_blocks[endblocki]);

    /* Start word */
    diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_trailing_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
      diff_64 = clear_lowbit_64(diff_64,relpos);
    }
#endif
    if (nmismatches >= max_mismatches) {
      mismatch_positions[nmismatches] = rightmost_value;
      return nmismatches;
    }

    query_high_ptr += 2; query_low_ptr += 2; query_flags_ptr += 2;
    ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
    offset += 64;


    /* Middle words */
    while (ref_high_ptr + 2 <= end_ptr) {
      diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));
      
#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
#endif
      if (nmismatches >= max_mismatches) {
	mismatch_positions[nmismatches] = rightmost_value;
	return nmismatches;
      }

      query_high_ptr += 2; query_low_ptr += 2; query_flags_ptr += 2;
      ref_high_ptr += 2; ref_low_ptr += 2; ref_flags_ptr += 2;
      offset += 64;
    }

    if (ref_high_ptr + 1 == end_ptr) {
      /* End 64-bit */
      diff_64 = (block_diff_64)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_end_64(diff_64,enddiscard + 32);
      debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_trailing_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_64(diff_64));
	diff_64 = clear_lowbit_64(diff_64,relpos);
      }
#endif
      mismatch_positions[nmismatches] = rightmost_value;
      return nmismatches;

    } else {
      /* End 32-bit */
      diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_end_32(diff_32,enddiscard);
      debug(printf("enddisc:    %08X\n",diff_32));

#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_trailing_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset + (relpos = count_trailing_zeroes_32(diff_32));
	diff_32 = clear_lowbit_32(diff_32,relpos);
      }
#endif
      mismatch_positions[nmismatches] = rightmost_value;
      return nmismatches;
    }
  }
}


/* Modified from Substring_trim_qend_nosplice */
/* was left_trim */
static int
qend_trim (int *fragment_nmismatches, int *mismatch_positions, int total_nmismatches,
	   int *nmatches_baseline, int fragment_length) {
  /* max_score of 2 sets an expectation of more matches than mismatches */
  int max_score = MIN_OVERALL_SCORE, score;
  int trimpos, prevpos, pos;
  int i;

  debug8(printf("Entering qend_trim with fragment length %d\n",fragment_length));

  *fragment_nmismatches = 0;
  if (total_nmismatches == 0) {
    debug8(printf("No mismatches: Trim qend trimpos %d, fragment_nmismatches 0\n",
		  fragment_length));
    return fragment_length;
  }

  /* Relative to pos5 */
  /* -1 | mismatch positions | +fragment_length */
  trimpos = 0;
  prevpos = -1;
  pos = mismatch_positions[0];
  /* Don't add mismatch initially because we stop before the mismatch */
  score = (pos - prevpos - 1)*TRIM_MATCH_SCORE /*+ TRIM_MISMATCH_SCORE_MULT*/;
  *fragment_nmismatches = 0;
  debug8(printf("initial pos %d, score %d",pos,score));
  if (score - nmatches_baseline[pos] >= MIN_INITIAL_SCORE) {
    debug8(printf(" **"));
    trimpos = pos;
    max_score = score - nmatches_baseline[pos];
  }
  debug8(printf("\n"));
  prevpos = pos;

  for (i = 1; i < total_nmismatches; i++) {
    pos = mismatch_positions[i];
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug8(printf("pos %d, score %d",pos,score));
    if (score - nmatches_baseline[pos] >= max_score) {
      debug8(printf(" **"));
      trimpos = pos;
      *fragment_nmismatches = i;
      max_score = score - nmatches_baseline[pos];
    }
    debug8(printf("\n"));
    prevpos = pos;
  }


  if (*fragment_nmismatches == total_nmismatches - 1) {
    /* If last mismatch compensated for previous, then take the last
       segment, regardless of whether it compensates for the last
       mismatch */
    trimpos = fragment_length;
    *fragment_nmismatches += 1;

  } else {
    /* See if last segment compensates */
    pos = fragment_length;
    debug8(printf("last segment has matches from prevpos %d to pos %d\n",prevpos,pos));
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug8(printf("pos %d, score %d",pos,score));
    if (score - nmatches_baseline[pos] >= max_score) {
      debug8(printf(" **"));
      trimpos = pos;
      *fragment_nmismatches = i;
      /* max_score = score - nmatches_baseline[pos]; */
    }
    debug8(printf("\n"));
    /* prevpos = pos; */
  }

  debug8(printf("Trim qend trimpos %d, fragment_nmismatches %d\n",trimpos,*fragment_nmismatches));
  return trimpos;		/* At the mismatch for qend */
}


/* Returns adj */
static int
solve_high_long (int *best_trimpos, int *nmismatches_to_trimpos,
		 Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
		 Compress_T query_compress, int *mismatch_positions_alloc,
		 Genomebits_T omebits, bool plusp, int genestrand) {
  int best_adj = 0;
  int nmismatches, fragment_nmismatches;
  int ins, del;
  int max_insertionlen, max_deletionlen;
  int genome_nwords, wordi;
  int leftshift, rightshift;
  Univcoord_T left, startpos, endpos;
  Genomecomp_T *ref_high_shifted, *ref_low_shifted, *ref_flags_shifted;
  Genomecomp_T *query_high_blocks, *query_low_blocks, *query_flags_blocks;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *query_high_adj, *query_low_adj, *query_flags_adj;
  int startblocki, endblocki, blocki;
  int startdiscard, enddiscard;
  int fragment_length, trimpos;

  int *nmatches_baseline;
  int nmatches, i, k;
  int pos, prevpos;
  int best_trimpos_adj;
  

  debug(printf("solve_high_long called with univdiagonal %u, pos %d..%d\n",
	       univdiagonal,pos5,pos3));

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
  }
  max_deletionlen = compute_max_deletionlen_exact(pos3 - pos5);

  fragment_length = pos3 - pos5;
  max_insertionlen = compute_max_insertionlen_exact(pos3 - pos5);
  
  if (max_insertionlen <= 0 && max_deletionlen <= 0) {
    return 0;
  } else {
    startpos = left + pos5;
    endpos = left + pos3 + max_deletionlen;
    genome_nwords = genomebits_shift(&ref_high_shifted,&ref_low_shifted,&ref_flags_shifted,
				     /*genome*/omebits,startpos,endpos);
    debug(printf("Genome:\n"));
    debug(shifted_print(ref_high_shifted,ref_low_shifted,ref_flags_shifted,
			/*startblocki*/0,/*endblocki*/genome_nwords - 1));
  }

  
  /* Normally we shift query right, but here we shift it leftward */

  /* Get baseline */
  Compress_shift(&query_high_blocks,&query_low_blocks,&query_flags_blocks,query_compress,
		 /*nshift*/0,/*initpos*/0);
  debug(printf("Query baseline:\n"));
  debug(Compress_print(query_compress,/*nshift*/0,pos5,pos3));
#ifdef DEBUG
  int query_nwords =
#endif
    querybits_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,
		    query_high_blocks,query_low_blocks,query_flags_blocks,
		    pos5,pos3);
  debug(printf("Query (%d words):\n",query_nwords));
  debug(shifted_print(query_high_shifted,query_low_shifted,query_flags_shifted,
		      /*startblocki*/0,/*endblocki*/query_nwords - 1));
  
  query_high_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_low_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_flags_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));


  /* Compute nmatches_baseline, the number of matches we would have without any indels.
     Used to provide a comparison for the indel scores */
  startblocki = 0;
  startdiscard = 0;
  endblocki = fragment_length/32;
  enddiscard = fragment_length % 32;
  debug(printf(">baseline 0 (high), startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
	       startblocki,endblocki,startdiscard,enddiscard));

  wordi = 0;
  for (blocki = startblocki; blocki <= endblocki; blocki++) {
    debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
    query_high_adj[blocki] = query_high_shifted[wordi];
    query_low_adj[blocki] = query_low_shifted[wordi];
    query_flags_adj[blocki] = query_flags_shifted[wordi];
    wordi++;
  }

  nmismatches = mismatches_fromleft(mismatch_positions_alloc,/*max_mismatches*/fragment_length,
				    /*rightmost_value*/fragment_length,
				    query_high_adj,query_low_adj,
				    query_flags_adj,ref_high_shifted,
				    ref_low_shifted,ref_flags_shifted,
				    startblocki,endblocki,startdiscard,enddiscard,
				    plusp,genestrand);

  nmatches = 0;
  nmatches_baseline = (int *) MALLOC((fragment_length+1)*sizeof(int));
  prevpos = -1;
  for (i = 0; i < nmismatches; i++) {
    pos = mismatch_positions_alloc[i];
    for (k = prevpos + 1; k < pos; k++) {
      nmatches_baseline[k] = ++nmatches;
    }
    nmatches_baseline[k] = nmatches;
    prevpos = pos;
  }
  for (k = prevpos + 1; k < fragment_length; k++) {
    nmatches_baseline[k] = ++nmatches;
  }
  nmatches_baseline[k] = nmatches;
  
#ifdef DEBUG8
  for (i = 0; i <= fragment_length; i++) {
    printf("%d: %d\n",i,nmatches_baseline[i]);
  }
#endif
  /* End of baseline */


  *nmismatches_to_trimpos = fragment_length;
  *best_trimpos = 0;		/* relative to pos5 */
  best_trimpos_adj = pos5;

  /* Deletions */
  for (del = 1; del <= max_deletionlen; del++) {
    startblocki = del/32;
    startdiscard = del % 32;
    endblocki = (del + (pos3 - pos5))/32;
    enddiscard = (del + (pos3 - pos5)) % 32;
    debug(printf(">del +%d/%d (high), startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 del,max_deletionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = del % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }

    nmismatches = mismatches_fromleft(mismatch_positions_alloc,/*max_mismatches*/fragment_length,
				      /*rightmost_value*/fragment_length,
				      query_high_adj,query_low_adj,
				      query_flags_adj,ref_high_shifted,
				      ref_low_shifted,ref_flags_shifted,
				      startblocki,endblocki,startdiscard,enddiscard,
				      plusp,genestrand);
    assert(mismatch_positions_alloc[nmismatches] == fragment_length);

#ifdef DEBUG
    printf("%d mismatches: ",nmismatches);
    for (int i = 0; i < nmismatches; i++) {
      printf(" %d",mismatch_positions_alloc[i]);
    }
    printf("\n");
#endif

    trimpos = qend_trim(&fragment_nmismatches,mismatch_positions_alloc,nmismatches,
			nmatches_baseline,fragment_length);
    debug(printf("del +%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, trimpos %d, fragment_nmismatches %d\n",
		 del,startblocki,endblocki,startdiscard,enddiscard,trimpos,fragment_nmismatches));

    if (trimpos > *best_trimpos) {
      *best_trimpos = trimpos;
      best_trimpos_adj = pos5;
      *nmismatches_to_trimpos = fragment_nmismatches;
      best_adj = +del;
    }

    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  debug(printf("After deletions, best_adj %d, best_trimpos %d, nmismatches_to_trimpos %d\n",
	       best_adj,*best_trimpos,*nmismatches_to_trimpos));
	 
  if (best_adj == 0) {
    /* Nothing found */

  } else if ((*best_trimpos) - (*nmismatches_to_trimpos) < 4) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = fragment_length;

  } else {
    /* Accept */
  }


  /* Insertions */
  startblocki = 0;
  startdiscard = 0;
  for (ins = 1; ins <= max_insertionlen; ins++) {
    endblocki = ((pos3 - pos5) - ins)/32;
    enddiscard = ((pos3 - pos5) - ins) % 32;
    debug(printf(">ins %d/%d (high), startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 -ins,-max_insertionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* This is the only case where we shift leftward */
    /* Shift query fragment leftward (by rightshift) */
    wordi = 0;
    if ((rightshift = ins % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      leftshift = 32 - rightshift;
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing rightshift of %d on wordi %d and leftshift of %d on wordi %d\n",
		     blocki,rightshift,wordi,leftshift,wordi+1));
	query_high_adj[blocki] = (query_high_shifted[wordi] >> rightshift) | (query_high_shifted[wordi+1] << leftshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] >> rightshift) | (query_low_shifted[wordi+1] << leftshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] >> rightshift) | (query_flags_shifted[wordi+1] << leftshift);
	wordi++;
      }
    }
    
    nmismatches = mismatches_fromleft(mismatch_positions_alloc,
				      /*max_mismatches*/fragment_length - ins,
				      /*rightmost_value*/(fragment_length - ins),
				      query_high_adj,query_low_adj,
				      query_flags_adj,ref_high_shifted,
				      ref_low_shifted,ref_flags_shifted,
				      startblocki,endblocki,startdiscard,enddiscard,
				      plusp,genestrand);
    assert(mismatch_positions_alloc[nmismatches] == fragment_length - ins);

#ifdef DEBUG
    printf("%d mismatches: ",nmismatches);
    for (int i = 0; i < nmismatches; i++) {
      printf(" %d",mismatch_positions_alloc[i]);
    }
    printf("\n");
#endif

    trimpos = qend_trim(&fragment_nmismatches,mismatch_positions_alloc,nmismatches,
			nmatches_baseline,/*fragment_length*/(fragment_length - ins));
    debug(printf("ins %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, trimpos %d, fragment_nmismatches %d\n",
		 -ins,startblocki,endblocki,startdiscard,enddiscard,trimpos,fragment_nmismatches));

    /* Previously required ins < fragment_nmismatches */
    if (trimpos > *best_trimpos) {
      *best_trimpos = trimpos;
      best_trimpos_adj = pos5 - ins;
      if ((*nmismatches_to_trimpos = fragment_nmismatches - ins) < 0) {
	*nmismatches_to_trimpos = 0;
      }
      best_adj = -ins;
    }
    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }


  if (best_adj >= 0) {
    /* No insertion found */
    /* Accept */
  } else if ((*best_trimpos) - (*nmismatches_to_trimpos) < 4) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = fragment_length;
  } else {
    /* Accept */
  }


  FREE(nmatches_baseline);

  FREE(query_flags_adj);
  FREE(query_low_adj);
  FREE(query_high_adj);

  FREE(query_flags_shifted);
  FREE(query_low_shifted);
  FREE(query_high_shifted);

  FREE(ref_flags_shifted);
  FREE(ref_low_shifted);
  FREE(ref_high_shifted);

  *best_trimpos += best_trimpos_adj;
  assert(*best_trimpos >= 0 && *best_trimpos <= querylength);
  return best_adj;
}


/* Returns the number of nindels.  A value of 0 indicates failure. */
int
Genomebits_indel_solve_high (int *best_trimpos, int *nmismatches_to_trimpos,
			     Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
			     Compress_T query_compress, int *mismatch_positions_alloc,
			     Genomebits_T omebits, Genomebits_T omebits_alt,
			     bool plusp, int genestrand) {
  int adj;
  int fragment_length = pos3 - pos5;

  if (fragment_length > 6) {
    if ((adj = solve_high_long(&(*best_trimpos),&(*nmismatches_to_trimpos),
			       univdiagonal,querylength,pos5,pos3,
			       query_compress,mismatch_positions_alloc,omebits,
			       plusp,genestrand)) != 0) {
      return adj;


    } else if (omebits_alt != NULL &&
	       (adj = solve_high_long(&(*best_trimpos),&(*nmismatches_to_trimpos),
				      univdiagonal,querylength,pos5,pos3,
				      query_compress,mismatch_positions_alloc,omebits_alt,
				      plusp,genestrand)) != 0) {
      return adj;

    } else {
      return 0;
    }

  } else if ((adj = solve_high_short(&(*nmismatches_to_trimpos),
				     univdiagonal,querylength,pos5,pos3,
				     query_compress,omebits,plusp,genestrand)) != 0) {
    *best_trimpos = pos3;
    return adj;

  } else if (omebits_alt != NULL &&
	     (adj = solve_high_short(&(*nmismatches_to_trimpos),
				     univdiagonal,querylength,pos5,pos3,
				     query_compress,omebits_alt,plusp,genestrand)) != 0) {
    *best_trimpos = pos3;
    return adj;

  } else {
    return 0;
  }
}


/************************************************************************
 * Low or qstart indels, with mismatches starting on right
 ************************************************************************/

/* Modified from genomebits_mismatches.c */
static int
mismatches_fromright (int *mismatch_positions, int max_mismatches, int leftmost_value,
		      Genomecomp_T *query_high_adj, Genomecomp_T *query_low_adj,
		      Genomecomp_T *query_flags_adj, Genomecomp_T *ref_high_blocks,
		      Genomecomp_T *ref_low_blocks, Genomecomp_T *ref_flags_blocks,
		      int startblocki, int endblocki, int startdiscard, int enddiscard,
		      bool plusp, int genestrand) {
  int nmismatches = 0, offset;
  Genomecomp_T *query_high_ptr, *query_low_ptr, *query_flags_ptr;
  Genomecomp_T *ref_high_ptr, *ref_low_ptr, *ref_flags_ptr, *start_ptr;
  UINT4 diff_32;
  UINT8 diff_64;
#if !defined(HAVE_SSE2)
  int relpos;
#endif

  /* Matching the parameters for Genomebits_mismatches_right_trim */
  bool query_unk_mismatch_p = false;
  bool genome_unk_mismatch_p = false;

  debug(printf("Entering mismatches_fromright with startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, max_mismatches %d\n",
	       startblocki,endblocki,startdiscard,enddiscard,max_mismatches));

  /* For leftward scanning */
  query_high_ptr = &(query_high_adj[endblocki]);
  query_low_ptr = &(query_low_adj[endblocki]);
  query_flags_ptr = &(query_flags_adj[endblocki]);

  ref_high_ptr = &(ref_high_blocks[endblocki]);
  ref_low_ptr = &(ref_low_blocks[endblocki]);
  ref_flags_ptr = &(ref_flags_blocks[endblocki]);

  /* offset = (32 - enddiscard) + pos3 - 1; */
  offset = 32 - enddiscard;	/* relative to pos3 */

  if (startblocki == endblocki) {
    /* Single 32-bit */
    diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
			      query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    debug(printf("before32:   %08X\n",diff_32));
    diff_32 = clear_end_32(diff_32,enddiscard);
    diff_32 = clear_start_32(diff_32,startdiscard);
    debug(printf("bothdisc:   %08X\n",diff_32));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_leading_32(mismatch_positions,/*nmismatches*/0,diff_32,offset,max_mismatches);
#else
    while (nonzero_p_32(diff_32) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
      diff_32 = clear_highbit_32(diff_32,relpos);
    }
#endif
    mismatch_positions[nmismatches] = leftmost_value;
    return nmismatches;

  } else if (startblocki + 1 == endblocki) {
    /* Single 64-bit */
    diff_64 = (block_diff_64)(query_high_ptr - 1,query_low_ptr - 1,
			      query_flags_ptr - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    debug(printf("before64:   %016lX\n",diff_64));
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    diff_64 = clear_start_64(diff_64,startdiscard);
    debug(printf("bothdisc:   %016lX\n",diff_64));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
#endif
    mismatch_positions[nmismatches] = leftmost_value;
    return nmismatches;

  } else {
    /* Multiple words */
    start_ptr = &(ref_high_blocks[startblocki]);

    /* End word */
    diff_64 = (block_diff_64)(query_high_ptr - 1,query_low_ptr - 1,
			      query_flags_ptr - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
			      plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
    diff_64 = clear_end_64(diff_64,enddiscard + 32);
    debug(printf("enddisc:    %016lX\n",diff_64));

#ifdef HAVE_SSE2
    nmismatches = Genomebits_decode_leading_64(mismatch_positions,/*nmismatches*/0,diff_64,offset,max_mismatches);
#else
    while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
      mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
      diff_64 = clear_highbit_64(diff_64,relpos);
    }
#endif
    if (nmismatches >= max_mismatches) {
      mismatch_positions[nmismatches] = leftmost_value;
      return nmismatches;
    }

    query_high_ptr -= 2; query_low_ptr -= 2; query_flags_ptr -= 2;
    ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
    offset -= 64;


    /* Middle words */
    while (ref_high_ptr >= start_ptr + 2) {
      diff_64 = (block_diff_64)(query_high_ptr - 1,query_low_ptr - 1,
				query_flags_ptr - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      /* No discards */
      debug(printf("nodisc:     %016lX\n",diff_64));

#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
#endif
      if (nmismatches >= max_mismatches) {
	mismatch_positions[nmismatches] = leftmost_value;
	return nmismatches;
      }

      query_high_ptr -= 2; query_low_ptr -= 2; query_flags_ptr -= 2;
      ref_high_ptr -= 2; ref_low_ptr -= 2; ref_flags_ptr -= 2;
      offset -= 64;
    }

    if (ref_high_ptr == start_ptr + 1) {
      /* Start 64-bit */
      diff_64 = (block_diff_64)(query_high_ptr - 1,query_low_ptr - 1,
				query_flags_ptr - 1,ref_high_ptr - 1,ref_low_ptr - 1,ref_flags_ptr - 1,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_64 = clear_start_64(diff_64,startdiscard);
      debug(printf("startdisc:  %016lX\n",diff_64));

#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_leading_64(&(mismatch_positions[nmismatches]),nmismatches,diff_64,offset,max_mismatches);
#else
      while (nonzero_p_64(diff_64) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_64(diff_64));
	diff_64 = clear_highbit_64(diff_64,relpos);
      }
#endif
      mismatch_positions[nmismatches] = leftmost_value;
      return nmismatches;

    } else {
      /* Start 32-bit */
      diff_32 = (block_diff_32)(query_high_ptr,query_low_ptr,
				query_flags_ptr,ref_high_ptr,ref_low_ptr,ref_flags_ptr,
				plusp,genestrand,query_unk_mismatch_p,genome_unk_mismatch_p);
      diff_32 = clear_start_32(diff_32,startdiscard);
      debug(printf("startdisc:  %08X\n",diff_32));

#ifdef HAVE_SSE2
      nmismatches = Genomebits_decode_leading_32(&(mismatch_positions[nmismatches]),nmismatches,diff_32,offset,max_mismatches);
#else
      while (nonzero_p_32(diff_32) && nmismatches < max_mismatches) {
	mismatch_positions[nmismatches++] = offset - (relpos = count_leading_zeroes_32(diff_32));
	diff_32 = clear_highbit_32(diff_32,relpos);
      }
#endif
      mismatch_positions[nmismatches] = leftmost_value;
      return nmismatches;
    }
  }
}


/* Modified from trim_qstart_nosplice */
/* was right_trim */
static int
qstart_trim (int *fragment_nmismatches, int *mismatch_positions, int total_nmismatches,
	     int *nmatches_baseline, int fragment_length) {
  /* max_score of 2 sets an expectation of more matches than mismatches */
  int max_score = MIN_OVERALL_SCORE, score;
  int trimpos, prevpos, pos;
  int i;

  debug8(printf("Entering qstart_trim with fragment length %d\n",fragment_length));
  assert(fragment_length > 0);

  *fragment_nmismatches = 0;
  if (total_nmismatches == 0) {
    debug8(printf("No mismatches: Trim qstart trimpos %d, fragment_nmismatches %d\n",
		  -fragment_length,0));
    return -fragment_length;
  }

  /* Relative to pos3 */
  /* +1 | mismatch positions (negative values) | -fragment_length */

  trimpos = 0;
  prevpos = -1;
  pos = -(mismatch_positions[0]);
  /* Don't add mismatch initially because we stop before the mismatch */
  score = (pos - prevpos - 1)*TRIM_MATCH_SCORE /*+ TRIM_MISMATCH_SCORE_MULT*/;
  *fragment_nmismatches = 0;
  debug8(printf("initial pos %d, score %d",pos,score));
  if (score - nmatches_baseline[pos] >= MIN_INITIAL_SCORE) {
    debug8(printf(" **"));
    trimpos = pos;
    max_score = score - nmatches_baseline[pos];
  }
  debug8(printf("\n"));
  prevpos = pos;

  for (i = 1; i < total_nmismatches; i++) {
    pos = -(mismatch_positions[i]);
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug8(printf("pos %d, score %d",pos,score));
    if (score - nmatches_baseline[pos] >= max_score) {
      debug8(printf(" **"));
      trimpos = pos;
      *fragment_nmismatches = i;
      max_score = score - nmatches_baseline[pos];
    }
    debug8(printf("\n"));
    prevpos = pos;
  }

  if (*fragment_nmismatches == total_nmismatches - 1) {
    /* If last mismatch compensated for previous, then take the last
       segment, regardless of whether it compensates for the last
       mismatch */
    trimpos = +fragment_length;
    *fragment_nmismatches = total_nmismatches;

  } else {
    /* See if last segment compensates */
    pos = +fragment_length;
    debug8(printf("last segment has matches from prevpos %d to pos %d\n",prevpos,pos));
    score += TRIM_MISMATCH_SCORE_MULT;
    score += (pos - prevpos - 1)*TRIM_MATCH_SCORE;
    debug8(printf("final pos %d, score %d",pos,score));
    if (score - nmatches_baseline[pos] >= max_score) {
      debug8(printf(" **"));
      trimpos = pos;
      *fragment_nmismatches = total_nmismatches;
      /* max_score = score - nmatches_baseline[pos]; */
    }
    debug8(printf("\n"));
    /* prevpos = pos; */
  }

  debug8(printf("Trim qstart trimpos %d, fragment_nmismatches %d\n",
		-trimpos,*fragment_nmismatches));

  return -trimpos;		/* This is one position after the mismatch for qstart */
}


static int
solve_low_long (int *best_trimpos, int *nmismatches_to_trimpos,
		Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
		Compress_T query_compress, int *mismatch_positions_alloc,
		Genomebits_T omebits, bool plusp, int genestrand) {
  int best_adj = 0;
  int nmismatches, fragment_nmismatches;
  int ins, del;
  int max_insertionlen, max_deletionlen;
  int genome_nwords, wordi;
  int leftshift, rightshift;
  Univcoord_T left, startpos, endpos;
  Genomecomp_T *ref_high_shifted, *ref_low_shifted, *ref_flags_shifted;
  Genomecomp_T *query_high_blocks, *query_low_blocks, *query_flags_blocks;
  Genomecomp_T *query_high_shifted, *query_low_shifted, *query_flags_shifted;
  Genomecomp_T *query_high_adj, *query_low_adj, *query_flags_adj;
  int startblocki, endblocki, blocki;
  int startdiscard, enddiscard;
  int fragment_length, trimpos;

  int *nmatches_baseline;
  int nmatches, i, k;
  int pos, prevpos;
  int best_trimpos_adj;


  debug(printf("solve_low_long called with univdiagonal %u, pos %d..%d\n",univdiagonal,pos5,pos3));

  /* Check for alignments before genome start */
  left = univdiagonal - querylength;
  if (univdiagonal + pos5 < (Univcoord_T) querylength) {
    pos5 = querylength - univdiagonal;
    max_deletionlen = 0;
  } else {
    max_deletionlen = compute_max_deletionlen_exact(pos3 - pos5);
    if (univdiagonal + pos5 < (Univcoord_T) (querylength + max_deletionlen)) {
      max_deletionlen = 0;
    }
  }

  fragment_length = pos3 - pos5;
  max_insertionlen = compute_max_insertionlen_exact(pos3 - pos5);
  
  if (max_insertionlen <= 0 && max_deletionlen <= 0) {
    return 0;
  } else {
    startpos = left + pos5 - max_deletionlen;
    endpos = left + pos3;
    genome_nwords = genomebits_shift(&ref_high_shifted,&ref_low_shifted,&ref_flags_shifted,
				     /*genome*/omebits,startpos,endpos);
    debug(printf("Genome:\n"));
    debug(shifted_print(ref_high_shifted,ref_low_shifted,ref_flags_shifted,
			/*startblocki*/0,/*endblocki*/genome_nwords - 1));
  }

  
  /* Normally we shift query right, but here we shift it leftward */

  /* Get baseline */
  Compress_shift(&query_high_blocks,&query_low_blocks,&query_flags_blocks,query_compress,
		 /*nshift*/0,/*initpos*/0);
  debug(printf("Query baseline:\n"));
  debug(Compress_print(query_compress,/*nshift*/0,pos5,pos3));
#ifdef DEBUG
  int query_nwords =
#endif
    querybits_shift(&query_high_shifted,&query_low_shifted,&query_flags_shifted,
		    query_high_blocks,query_low_blocks,query_flags_blocks,
		    pos5,pos3);
  debug(printf("Query (%d words):\n",query_nwords));
  debug(shifted_print(query_high_shifted,query_low_shifted,query_flags_shifted,
		      /*startblocki*/0,/*endblocki*/query_nwords - 1));

  query_high_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_low_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));
  query_flags_adj = (Genomecomp_T *) MALLOC(genome_nwords*sizeof(Genomecomp_T));

  /* Compute nmatches_baseline, the number of matches we would have without any indels.
     Used to provide a comparison for the indel scores */
  startblocki = max_deletionlen / 32;
  startdiscard = max_deletionlen % 32;
  endblocki = (max_deletionlen + fragment_length) / 32;
  enddiscard = (max_deletionlen + fragment_length) % 32;

  wordi = 0;
  if ((leftshift = max_deletionlen % 32) == 0) {
    for (blocki = startblocki; blocki <= endblocki; blocki++) {
      debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi];
      query_low_adj[blocki] = query_low_shifted[wordi];
      query_flags_adj[blocki] = query_flags_shifted[wordi];
      wordi++;
    }

  } else {
    rightshift = 32 - leftshift;
    blocki = startblocki;
    debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		 blocki,leftshift,wordi));
    query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
    query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
    query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
    wordi++;

    for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
      debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		   blocki,leftshift,wordi,rightshift,wordi-1));
      query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
      query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
      query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
      wordi++;
    }
  }

  nmismatches = mismatches_fromright(mismatch_positions_alloc,/*max_mismatches*/fragment_length,
				     /*leftmost_value*/-(fragment_length),
				     query_high_adj,query_low_adj,
				     query_flags_adj,ref_high_shifted,
				     ref_low_shifted,ref_flags_shifted,
				     startblocki,endblocki,startdiscard,enddiscard,
				     plusp,genestrand);

  nmatches = 0;
  nmatches_baseline = (int *) MALLOC((fragment_length+1)*sizeof(int));
  prevpos = -1;
  for (i = 0; i < nmismatches; i++) {
    pos = -(mismatch_positions_alloc[i]);
    for (k = prevpos + 1; k < pos; k++) {
      nmatches_baseline[k] = ++nmatches;
    }
    nmatches_baseline[k] = nmatches;
    prevpos = pos;
  }
  for (k = prevpos + 1; k < fragment_length; k++) {
    nmatches_baseline[k] = nmatches;
  }
  nmatches_baseline[k] = nmatches;
  
#ifdef DEBUG8
  for (i = 0; i <= fragment_length; i++) {
    printf("%d: %d\n",i,nmatches_baseline[i]);
  }
#endif
  /* End of baseline */


  *nmismatches_to_trimpos = fragment_length;
  *best_trimpos = fragment_length; /* relative to pos3 */
  best_trimpos_adj = pos3;

  /* Deletions */
  for (del = 1; del <= max_deletionlen; del++) {
    startblocki = (max_deletionlen - del) / 32;
    startdiscard = (max_deletionlen - del) % 32;
    endblocki = (max_deletionlen - del + (pos3 - pos5)) / 32;
    enddiscard = (max_deletionlen - del + (pos3 - pos5)) % 32;
    debug(printf(">del +%d/%d (low), startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 del,max_deletionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = (max_deletionlen - del) % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }
      
    nmismatches = mismatches_fromright(mismatch_positions_alloc,/*max_mismatches*/fragment_length,
				       /*leftmost_value*/-(fragment_length),
				       query_high_adj,query_low_adj,
				       query_flags_adj,ref_high_shifted,
				       ref_low_shifted,ref_flags_shifted,
				       startblocki,endblocki,startdiscard,enddiscard,
				       plusp,genestrand);
    assert(mismatch_positions_alloc[nmismatches] == -fragment_length);
#ifdef DEBUG
    printf("%d mismatches: ",nmismatches);
    for (int i = 0; i < nmismatches; i++) {
      printf(" %d",mismatch_positions_alloc[i]);
    }
    printf("\n");
#endif

    trimpos = qstart_trim(&fragment_nmismatches,mismatch_positions_alloc,nmismatches,
			  nmatches_baseline,fragment_length);
    debug(printf("del +%d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, trimpos %d, fragment_nmismatches %d\n",
		 del,startblocki,endblocki,startdiscard,enddiscard,trimpos,fragment_nmismatches));

    if (trimpos < *best_trimpos) {
      *best_trimpos = trimpos;
      best_trimpos_adj = pos3;
      *nmismatches_to_trimpos = fragment_nmismatches;
      best_adj = +del;
    }

    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  debug(printf("After deletions, best_adj %d, best_trimpos %d, nmismatches_to_trimpos %d\n",
	       best_adj,*best_trimpos,*nmismatches_to_trimpos));

  /* best_trimpos is negative */
  if (best_adj == 0) {
    /* Nothing found */
  } else if (-(*best_trimpos) - (*nmismatches_to_trimpos) < 4) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = fragment_length;
  } else {
    /* Accept */
  }


  /* Insertions */
  endblocki = (max_deletionlen + (pos3 - pos5))/32;
  enddiscard = (max_deletionlen + (pos3 - pos5)) % 32;
  for (ins = 1; ins <= max_insertionlen; ins++) {
    startblocki = (max_deletionlen + ins) / 32;
    startdiscard = (max_deletionlen + ins) % 32;
    debug(printf(">ins %d/%d (low), startblocki %d, endblocki %d, startdiscard %d, enddiscard %d\n",
		 -ins,-max_insertionlen,startblocki,endblocki,startdiscard,enddiscard));

    /* Shift query fragment rightward (by leftshift) */
    wordi = 0;
    if ((leftshift = (max_deletionlen + ins) % 32) == 0) {
      for (blocki = startblocki; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, copying wordi %d\n",blocki,wordi));
	query_high_adj[blocki] = query_high_shifted[wordi];
	query_low_adj[blocki] = query_low_shifted[wordi];
	query_flags_adj[blocki] = query_flags_shifted[wordi];
	wordi++;
      }

    } else {
      rightshift = 32 - leftshift;
      blocki = startblocki;
      debug(printf("For startblocki %d, doing leftshift of %d on wordi %d\n",
		   blocki,leftshift,wordi));
      query_high_adj[blocki] = query_high_shifted[wordi] << leftshift;
      query_low_adj[blocki] = query_low_shifted[wordi] << leftshift;
      query_flags_adj[blocki] = query_flags_shifted[wordi] << leftshift;
      wordi++;

      for (blocki = startblocki + 1; blocki <= endblocki; blocki++) {
	debug(printf("For blocki %d, doing leftshift of %d on wordi %d and rightshift of %d on wordi %d\n",
		     blocki,leftshift,wordi,rightshift,wordi-1));
	query_high_adj[blocki] = (query_high_shifted[wordi] << leftshift) | (query_high_shifted[wordi-1] >> rightshift);
	query_low_adj[blocki] = (query_low_shifted[wordi] << leftshift) | (query_low_shifted[wordi-1] >> rightshift);
	query_flags_adj[blocki] = (query_flags_shifted[wordi] << leftshift) | (query_flags_shifted[wordi-1] >> rightshift);
	wordi++;
      }
    }

    nmismatches = mismatches_fromright(mismatch_positions_alloc,
				       /*max_mismatches*/fragment_length - ins,
				       /*leftmost_value*/-(fragment_length - ins),
				       query_high_adj,query_low_adj,
				       query_flags_adj,ref_high_shifted,
				       ref_low_shifted,ref_flags_shifted,
				       startblocki,endblocki,startdiscard,enddiscard,
				       plusp,genestrand);
    assert(mismatch_positions_alloc[nmismatches] == -(fragment_length - ins));
#ifdef DEBUG
    printf("%d mismatches: ",nmismatches);
    for (int i = 0; i < nmismatches; i++) {
      printf(" %d",mismatch_positions_alloc[i]);
    }
    printf("\n");
#endif

    trimpos = qstart_trim(&fragment_nmismatches,mismatch_positions_alloc,nmismatches,
			  nmatches_baseline,/*fragment_length*/fragment_length - ins);
    debug(printf("ins %d, startblocki %d, endblocki %d, startdiscard %d, enddiscard %d, trimpos %d, fragment_nmismatches %d\n",
		 -ins,startblocki,endblocki,startdiscard,enddiscard,trimpos,fragment_nmismatches));

    if (trimpos < *best_trimpos) {
      *best_trimpos = trimpos;
      best_trimpos_adj = pos3 - ins;
      *nmismatches_to_trimpos = fragment_nmismatches;
      best_adj = -ins;
    } 

    debug(shifted_print(query_high_adj,query_low_adj,query_flags_adj,startblocki,endblocki));
  }

  debug(printf("After insertions, best_adj %d, best_trimpos %d, nmismatches_to_trimpos %d\n",
	       best_adj,*best_trimpos,*nmismatches_to_trimpos));

  /* best_trimpos is negative */
  if (best_adj >= 0) {
    /* No insertion found */
  } else if (-(*best_trimpos) - (*nmismatches_to_trimpos) < 4) {
    /* Reject */
    best_adj = 0;
    *nmismatches_to_trimpos = fragment_length;
  } else {
    /* Accept */
  }


  FREE(nmatches_baseline);

  FREE(query_flags_adj);
  FREE(query_low_adj);
  FREE(query_high_adj);

  FREE(query_flags_shifted);
  FREE(query_low_shifted);
  FREE(query_high_shifted);

  FREE(ref_flags_shifted);
  FREE(ref_low_shifted);
  FREE(ref_high_shifted);


  *best_trimpos += best_trimpos_adj;
  assert(*best_trimpos >= 0 && *best_trimpos <= querylength);
  return best_adj;
}


/* Returns the number of nindels.  A value of 0 indicates failure. */
int
Genomebits_indel_solve_low (int *best_trimpos, int *nmismatches_to_trimpos,
			    Univcoord_T univdiagonal, int querylength, int pos5, int pos3,
			    Compress_T query_compress, int *mismatch_positions_alloc,
			    Genomebits_T omebits, Genomebits_T omebits_alt,
			    bool plusp, int genestrand) {
  int adj;
  int fragment_length = pos3 - pos5;

  if (fragment_length > 6) {
    if ((adj = solve_low_long(&(*best_trimpos),&(*nmismatches_to_trimpos),
			      univdiagonal,querylength,pos5,pos3,
			      query_compress,mismatch_positions_alloc,omebits,
			      plusp,genestrand)) != 0) {
      return adj;
      
    } else if (omebits_alt != NULL &&
	       (adj = solve_low_long(&(*best_trimpos),&(*nmismatches_to_trimpos),
				     univdiagonal,querylength,pos5,pos3,
				     query_compress,mismatch_positions_alloc,omebits_alt,
				     plusp,genestrand)) != 0) {
      return adj;

    } else {
      return 0;
    }

  } else if ((adj = solve_low_short(&(*nmismatches_to_trimpos),
				    univdiagonal,querylength,pos5,pos3,
				    query_compress,omebits,plusp,genestrand)) != 0) {
    *best_trimpos = pos5;
    return adj;

  } else if (omebits_alt != NULL &&
	     (adj = solve_low_short(&(*nmismatches_to_trimpos),
				    univdiagonal,querylength,pos5,pos3,
				    query_compress,omebits_alt,plusp,genestrand)) != 0) {	     
    *best_trimpos = pos5;
    return adj;

  } else {
    return 0;
  }
}


void
Genomebits_indel_setup (int max_insertionlen_in, int max_deletionlen_in,
			bool query_unk_mismatch_p_in, bool genome_unk_mismatch_p_in,
			Mode_T mode, bool maskedp) {

  max_insertionlen = max_insertionlen_in;
  max_deletionlen = max_deletionlen_in;
  
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



