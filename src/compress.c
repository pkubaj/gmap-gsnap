#define DEBUG1 1
static char rcsid[] = "$Id: cd17908a7856ea455ccae738b96b631530d1523d $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "compress.h"


#include <stdio.h>
#include <stdlib.h>		/* For posix_memalign */
#include <stddef.h>
#include <string.h>

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSSE3
#include <tmmintrin.h>
#endif
#if defined(HAVE_AVX2) && defined(HAVE_MM256_LOADSTORE_U2)
#include <immintrin.h>
#endif
#endif


#include "assert.h"
#include "mem.h"		/* For Compress_new_fwd and Compress_new_rev */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Compress_print */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* SIMD */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define NSHIFTS 32

#define T Compress_T
struct T {
  bool fwdp;			/* For debugging */

  int length;			/* For Compress_print_blocks */
  int nwords;			/* Number of 32-bit words  */

  Genomecomp_T *high_blocks;
  Genomecomp_T *low_blocks;
  Genomecomp_T *flags_blocks;

  Genomecomp_T *high_shift_array[NSHIFTS];
  Genomecomp_T *low_shift_array[NSHIFTS];
  Genomecomp_T *flags_shift_array[NSHIFTS];

  bool flag_set_p;
  bool availp[NSHIFTS];
};


void
Compress_free (T *old) {

#if defined(HAVE_ARM)
  FREE((*old)->flags_shift_array[0]);
  FREE((*old)->low_shift_array[0]);
  FREE((*old)->high_shift_array[0]);

#elif !defined(HAVE_SSE2)
  FREE((*old)->flags_shift_array[0]);
  FREE((*old)->low_shift_array[0]);
  FREE((*old)->high_shift_array[0]);

#else
  _mm_free((*old)->flags_shift_array[0]);
  _mm_free((*old)->low_shift_array[0]);
  _mm_free((*old)->high_shift_array[0]);
#endif

  FREE(*old);

  return;
}


#ifdef DEBUG
static void
print_vector_hex (__m128i x) {
  printf("%08X %08X %08X %08X\n",
	 _mm_extract_epi32(x,3),_mm_extract_epi32(x,2),_mm_extract_epi32(x,1),_mm_extract_epi32(x,0));
  return;
}
#endif


#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
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


#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
void
Compress_print (T this, int nshift, int pos5, int pos3) {
  Genomecomp_T *high_blocks, *low_blocks, *flags_blocks;
  int ptr, ptr5, ptr3, endptr, i;
  int startdiscard, enddiscard;

  high_blocks = this->high_shift_array[nshift];
  low_blocks = this->low_shift_array[nshift];
  flags_blocks = this->flags_shift_array[nshift];

  ptr = (nshift + 0)/32U;
  ptr5 = (nshift + pos5)/32U;
  ptr3 = (nshift + pos3)/32U;
  endptr = (nshift + this->length)/32U;

  startdiscard = (nshift + pos5) % 32;
  enddiscard = (nshift + pos3) % 32;

  while (ptr <= endptr) {
    if (ptr == ptr5) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 -- spacing*/
      printf("                                              \t");
      /* printf("%u\t",wordi); */
      for (i = 0; i < startdiscard; i++) {
	printf("v");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
    }

    printf("high: %08X  low: %08X  flags: %08X\t",
	   high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr]);
    write_chars(high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr]);
    printf("\n");

    if (ptr == ptr3) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 -- spacing */
      printf("                                              \t");
      /* printf("%u\t",wordi); */
      for (i = 0; i < enddiscard; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("^");
      }
      printf("\n");
    }

    ptr += 1;
  }

  printf("\n");
  return;
}
#endif


#ifdef DEBUG1
static void
write_chars_wdiscard (Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags,
		      int startdiscard, int enddiscard) {
  char Buffer[33];
  int k = 0, i;

  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < startdiscard; i++) {
    high >>= 1;
    low >>= 1;
  }
  for (i = startdiscard, k = 0; i < enddiscard; i++, k++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: Buffer[k] = 'A'; break;
    case 1U: Buffer[k] = 'C'; break;
    case 2U: Buffer[k] = 'G'; break;
    case 3U: Buffer[k] = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }
  Buffer[k] = '\0';

  if (flags != 0U) {
    for (i = 0; i < startdiscard; i++) {
      flags >>= 1;
    }
    for (i = startdiscard, k = 0; i < enddiscard; i++, k++) {
      if (flags & 0x01) {
	Buffer[k] = 'N';
      }
      flags >>= 1;
    }
  }

  printf("%s",Buffer);
  return;
}


static int
stream_chars_wdiscard (char *ptr, Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags,
		       int startdiscard, int enddiscard) {
  char *start = ptr;
  int k = 0, i;

  /* printf("%08X %08X %08X => ",high,low,flags); */

  for (i = 0; i < startdiscard; i++) {
    high >>= 1;
    low >>= 1;
  }
  for (i = startdiscard, k = 0; i < enddiscard; i++, k++) {
    switch (((high & 0x01) << 1) | (low & 0x01)) {
    case 0U: *ptr++ = 'A'; break;
    case 1U: *ptr++ = 'C'; break;
    case 2U: *ptr++ = 'G'; break;
    case 3U: *ptr++ = 'T'; break;
    default: abort();
    }
    high >>= 1;
    low >>= 1;
  }

  if (flags != 0U) {
    for (i = 0; i < startdiscard; i++) {
      flags >>= 1;
    }
    for (i = startdiscard, k = 0, ptr = start; i < enddiscard; i++, k++, ptr++) {
      if (flags & 0x01) {
	*ptr = 'N';
      }
      flags >>= 1;
    }
  }

  return (ptr - start);
}
#endif


#ifdef DEBUG1
void
Compress_print_queryseq (T this, int pos5, int pos3) {
  Genomecomp_T *high_blocks, *low_blocks, *flags_blocks;
  int ptr, ptr5, ptr3;
  int startdiscard, enddiscard;
  int nshift = 0;

  high_blocks = this->high_shift_array[nshift];
  low_blocks = this->low_shift_array[nshift];
  flags_blocks = this->flags_shift_array[nshift];

  /* ptr = (nshift + 0)/32U; */
  ptr5 = (nshift + pos5)/32U;
  ptr3 = (nshift + pos3)/32U;
  /* endptr = (nshift + this->length)/32U; */

  startdiscard = (nshift + pos5) % 32;
  enddiscard = (nshift + pos3) % 32;

  ptr = ptr5;
  write_chars_wdiscard(high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
		       startdiscard,/*enddiscard*/32);
  ptr += 1;
  while (ptr < ptr3) {
    write_chars_wdiscard(high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
			 /*startdiscard*/0,/*enddiscard*/32);
    ptr += 1;
  }
  write_chars_wdiscard(high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
		       /*startdiscard*/0,enddiscard);

  return;
}


char *
Compress_queryseq (T this, int querylength) {
  char *queryseq, *p;
  Genomecomp_T *high_blocks, *low_blocks, *flags_blocks;
  int ptr, ptr5, ptr3;
  int startdiscard, enddiscard;
  int nshift = 0;
  int pos5 = 0, pos3 = querylength;

  p = queryseq = (char *) MALLOC((querylength+1)*sizeof(char));

  high_blocks = this->high_shift_array[nshift];
  low_blocks = this->low_shift_array[nshift];
  flags_blocks = this->flags_shift_array[nshift];

  /* ptr = (nshift + 0)/32U; */
  ptr5 = (nshift + pos5)/32U;
  ptr3 = (nshift + pos3)/32U;
  /* endptr = (nshift + this->length)/32U; */

  startdiscard = (nshift + pos5) % 32;
  enddiscard = (nshift + pos3) % 32;

  ptr = ptr5;
  p += stream_chars_wdiscard(p,high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
			     startdiscard,/*enddiscard*/32);
  ptr += 1;
  while (ptr < ptr3) {
    p += stream_chars_wdiscard(p,high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
			       /*startdiscard*/0,/*enddiscard*/32);
    ptr += 1;
  }
  p += stream_chars_wdiscard(p,high_blocks[ptr],low_blocks[ptr],flags_blocks[ptr],
			     /*startdiscard*/0,enddiscard);

  *p = '\0';
  return queryseq;
}
#endif



/*                   87654321 */
#define LEFT_SET   0x80000000
#define LEFT_CLEAR 0x00000000

#define EXTRA_WORDS 1 /* 1 so we can shift to next word */

T
Compress_new_fwd (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T position;
  int wordi;
  int c, i;
  int in_counter = 0;
  int nblocks;			/* Number of 128-bit registers */

  new->fwdp = true;
  new->length = length;
  new->flag_set_p = false;


#if defined(HAVE_ARM)
  /* Align to 128-bit boundaries, so we can use aligned SIMD loads for high_blocks, low_blocks, and flags_blocks */
  nblocks = (length+127)/128U;
  new->nwords = nblocks * 4;

  posix_memalign((void **) &(new->high_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  posix_memalign((void **) &(new->low_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  posix_memalign((void **) &(new->flags_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));

#elif !defined(HAVE_SSE2)
  /* Using 32-bit quantities */
  nblocks = (length+31)/32U;
  new->nwords = nblocks * 1;

  new->high_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  new->low_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  new->flags_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));

#else  
  /* Align to 128-bit boundaries, so we can use aligned SIMD loads for high_blocks, low_blocks, and flags_blocks */
  nblocks = (length+127)/128U;
  new->nwords = nblocks * 4;

  new->high_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
  new->low_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
  new->flags_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
#endif

  /* Fill nshift of 0 with the original sequence; these will be aligned to a 128-bit boundary */
  new->high_blocks = new->high_shift_array[0];
  new->low_blocks = new->low_shift_array[0];
  new->flags_blocks = new->flags_shift_array[0];
  new->availp[0] = true;	/* Same as new->blocks */

  for (i = 1; i < NSHIFTS; i++) {
    new->high_shift_array[i] = &(new->high_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->low_shift_array[i] = &(new->low_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->flags_shift_array[i] = &(new->flags_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->availp[i] = false;
  }


  wordi = 0;
  position = 0U;
  while (position < length) {
#if !defined(HAVE_SSE2)
    high = low = flags = 0U;
    in_counter = 0;
    while (position < length && in_counter < 32) {
      c = gbuffer[position++];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'A': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'C': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'G':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'T':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }
      
    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    new->high_blocks[wordi] = high;
    new->low_blocks[wordi] = low;
    new->flags_blocks[wordi] = flags;
    wordi++;
    
#else
    for (i = 0; i < 4; i++) {
      high = low = flags = 0U;
      in_counter = 0;
      while (position < length && in_counter < 32) {
	c = gbuffer[position++];
	high >>= 1;
	low >>= 1;
	flags >>= 1;
	
	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case 'A': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'C': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	case 'G':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'T':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
	}
	in_counter++;
      }
      
      while (in_counter < 32) {
	high >>= 1;
	low >>= 1;
	flags >>= 1;
	in_counter++;
      }

      new->high_blocks[wordi + i] = high;
      new->low_blocks[wordi + i] = low;
      new->flags_blocks[wordi + i] = flags;
      if (flags) {
	new->flag_set_p = true;
      }
    }

    wordi += 4;
#endif
  }

  /* For the extra word at the end */
  new->high_blocks[wordi] = 0U;
  new->low_blocks[wordi] = 0U;
  new->flags_blocks[wordi] = 0U;

  debug0(printf("Compress_new_fwd\n"));
  debug0(Compress_print_blocks(new->blocks,new->nblocks));
  debug0(printf("\n"));

  return new;
}

T
Compress_new_rev (char *gbuffer, Chrpos_T length) {
  T new = (T) MALLOC(sizeof(*new));
  Genomecomp_T high, low, flags;
  Chrpos_T position;
  int wordi;
  int c, i;
  int in_counter = 0;
  int nblocks;			/* Number of 128-bit registers */

  new->fwdp = false;
  new->length = length;
  new->flag_set_p = false;

#if defined(HAVE_ARM)
  /* Align to 128-bit boundaries, so we can use aligned SIMD loads for high_blocks, low_blocks, and flags_blocks */
  nblocks = (length+127)/128U;
  new->nwords = nblocks * 4;

  posix_memalign((void **) &(new->high_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  posix_memalign((void **) &(new->low_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  posix_memalign((void **) &(new->flags_shift_array[0]),/*alignment*/16,NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));

#elif !defined(HAVE_SSE2)
  /* Using 32-bit quantities */
  nblocks = (length+31)/32U;
  new->nwords = nblocks * 1;

  new->high_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  new->low_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));
  new->flags_shift_array[0] = (Genomecomp_T *) MALLOC(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T));

#else  
  nblocks = (length+127)/128U;
  new->nwords = nblocks * 4;

  /* Align to 128-bit boundaries, so we can use aligned SIMD loads for high_blocks, low_blocks, and flags_blocks */
  new->high_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
  new->low_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
  new->flags_shift_array[0] = (Genomecomp_T *) _mm_malloc(NSHIFTS * (new->nwords + EXTRA_WORDS) * sizeof(Genomecomp_T),16);
#endif

  /* Fill nshift of 0 with the original sequence; these will be aligned to a 128-bit boundary */
  new->high_blocks = new->high_shift_array[0];
  new->low_blocks = new->low_shift_array[0];
  new->flags_blocks = new->flags_shift_array[0];
  new->availp[0] = true;

  for (i = 1; i < NSHIFTS; i++) {
    new->high_shift_array[i] = &(new->high_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->low_shift_array[i] = &(new->low_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->flags_shift_array[i] = &(new->flags_shift_array[i-1][new->nwords + EXTRA_WORDS]);
    new->availp[i] = false;
  }

  wordi = 0;
  position = length;
  while (position > 0) {

#if !defined(HAVE_SSE2)
    high = low = flags = 0U;
    in_counter = 0;
    while (position > 0 && in_counter < 32) {
      c = gbuffer[--position];
      high >>= 1;
      low >>= 1;
      flags >>= 1;

      /* Assume that gbuffer is upper case */
      switch /*(uppercaseCode[c])*/ (c) {
      case 'T': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'G': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      case 'C':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
      case 'A':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
      default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
      }
      in_counter++;
    }

    while (in_counter < 32) {
      high >>= 1;
      low >>= 1;
      flags >>= 1;
      in_counter++;
    }

    new->high_blocks[wordi] = high;
    new->low_blocks[wordi] = low;
    new->flags_blocks[wordi] = flags;
    wordi++;

#else
    for (i = 0; i < 4; i++) {
      high = low = flags = 0U;
      in_counter = 0;
      while (position > 0 && in_counter < 32) {
	c = gbuffer[--position];
	high >>= 1;
	low >>= 1;
	flags >>= 1;

	/* Assume that gbuffer is upper case */
	switch /*(uppercaseCode[c])*/ (c) {
	case 'T': /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'G': /* high |= LEFT_CLEAR; */    low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	case 'C':    high |= LEFT_SET;      /* low |= LEFT_CLEAR; */ /* flags |= LEFT_CLEAR; */ break;
	case 'A':    high |= LEFT_SET;         low |= LEFT_SET;      /* flags |= LEFT_CLEAR; */ break;
	default:  /* high |= LEFT_CLEAR; */ /* low |= LEFT_CLEAR; */    flags |= LEFT_SET;
	}
	in_counter++;
      }

      while (in_counter < 32) {
	high >>= 1;
	low >>= 1;
	flags >>= 1;
	in_counter++;
      }

      new->high_blocks[wordi + i] = high;
      new->low_blocks[wordi + i] = low;
      new->flags_blocks[wordi + i] = flags;
      if (flags) {
	new->flag_set_p = true;
      }
    }
    wordi += 4;
#endif
  }

  /* For the extra word at the end */
  new->high_blocks[wordi] = 0U;
  new->low_blocks[wordi] = 0U;
  new->flags_blocks[wordi] = 0U;

  debug0(printf("Compress_new_rev\n"));
  debug0(Compress_print_blocks(new->blocks,new->nblocks));
  debug0(printf("\n"));

  return new;
}


#if !defined(HAVE_SSE2) || defined(SSE2_SLLI_CONST_IMM8)
void
Compress_shift (Genomecomp_T **query_high_shifted, Genomecomp_T **query_low_shifted,
		Genomecomp_T **query_flags_shifted, T this, int nshift, int initpos) {
  Genomecomp_T *high_ptr, *low_ptr, *flags_ptr;
  int rightshift, initshift;
  int wordi;

  if (this->availp[nshift] == false) {
    wordi = this->nwords;
    rightshift = 32 - nshift;

    high_ptr = &(this->high_shift_array[nshift][wordi]);
    low_ptr = &(this->low_shift_array[nshift][wordi]);
    flags_ptr = &(this->flags_shift_array[nshift][wordi]);

    while (wordi > 0) {
      *high_ptr-- = (this->high_blocks[wordi] << nshift) | (this->high_blocks[wordi-1] >> rightshift);
      *low_ptr-- = (this->low_blocks[wordi] << nshift) | (this->low_blocks[wordi-1] >> rightshift);
      *flags_ptr-- = (this->flags_blocks[wordi] << nshift) | (this->flags_blocks[wordi-1] >> rightshift);
      wordi--;
    }

    *high_ptr = this->high_blocks[0] << nshift;
    *low_ptr = this->low_blocks[0] << nshift;
    *flags_ptr = this->flags_blocks[0] << nshift;

    this->availp[nshift] = true;
  }

#if 0
  printf("nshift %d\n",nshift);
  for (int wordi = 0; wordi <= this->nwords; wordi++) {
    printf("%08X %08X %08X\n",
	   this->high_shift_array[nshift][wordi],
	   this->low_shift_array[nshift][wordi],
	   this->flags_shift_array[nshift][wordi]);
  }
  exit(9);
#endif

  initshift = (nshift + initpos)/32;

  *query_high_shifted = &(this->high_shift_array[nshift][initshift]);
  *query_low_shifted = &(this->low_shift_array[nshift][initshift]);
  *query_flags_shifted = &(this->flags_shift_array[nshift][initshift]);

  return;
}

#else
void
Compress_shift (Genomecomp_T **query_high_shifted, Genomecomp_T **query_low_shifted,
		Genomecomp_T **query_flags_shifted, T this, int nshift, int initpos) {
  Genomecomp_T *high_shifted, *low_shifted, *flags_shifted, *high_ptr, *low_ptr, *flags_ptr;
  int rightshift, initshift;
  __m128i result, current, prev, old_prev;
#if defined(HAVE_AVX2) && defined(HAVE_MM256_LOADSTORE_U2)
  __m256i result256, current256, prev256, old_prev256;
#endif

  if (this->availp[nshift] == false) {
    rightshift = 32 - nshift;

#if defined(HAVE_AVX2) && defined(HAVE_MM256_LOADSTORE_U2)
    /* High and low */
    high_shifted = &(this->high_shift_array[nshift][this->nwords - 3]);
    low_shifted = &(this->low_shift_array[nshift][this->nwords - 3]);
    high_ptr = &(this->high_blocks[this->nwords - 3]);
    low_ptr = &(this->low_blocks[this->nwords - 3]);
    
    prev256 = _mm256_loadu2_m128i((__m128i *) (high_ptr - 1),(__m128i *) (low_ptr - 1));
    current256 = _mm256_bsrli_epi128(prev256,4); /* Works on 128-bit lanes separately */
    result256 = _mm256_or_si256(_mm256_slli_epi32(current256,nshift),_mm256_srli_epi32(prev256,rightshift));
    _mm256_storeu2_si128((__m128i *) high_shifted,(__m128i *) low_shifted,result256);
    high_shifted -= 4; low_shifted -= 4; high_ptr -= 4; low_ptr -= 4;

    while (high_ptr > this->high_blocks) {
      old_prev256 = prev256;
      prev256 = _mm256_loadu2_m128i((__m128i *) (high_ptr - 1),(__m128i *) (low_ptr - 1));
      current256 = _mm256_alignr_epi8(old_prev256,prev256,4);
      result256 = _mm256_or_si256(_mm256_slli_epi32(current256,nshift),_mm256_srli_epi32(prev256,rightshift));
      _mm256_storeu2_si128((__m128i *) high_shifted,(__m128i *) low_shifted,result256);
      high_shifted -= 4; low_shifted -= 4; high_ptr -= 4; low_ptr -= 4;
    }

#else
    /* High */
    high_shifted = &(this->high_shift_array[nshift][this->nwords - 3]);
    high_ptr = &(this->high_blocks[this->nwords - 3]);

    /* First 128-bit register has 0 bytes in the high word of current */
    prev = _mm_load_si128((__m128i *) (high_ptr - 1)); /* (high_ptr - 1) is aligned */
    current = _mm_srli_si128(prev,4);
    result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
    debug(printf("current:   "); print_vector_hex(current);
	  printf("prev:      "); print_vector_hex(prev);
	  printf("result:    "); print_vector_hex(result));

    _mm_storeu_si128((__m128i *) high_shifted,result);
    high_shifted -= 4; high_ptr -= 4;

    while (high_ptr > this->high_blocks) {
#ifdef HAVE_SSSE3
      old_prev = prev;
      prev = _mm_load_si128((__m128i *) (high_ptr - 1));
      current = _mm_alignr_epi8(old_prev,prev,4);
#else 
      prev = _mm_load_si128((__m128i *) (high_ptr - 1));
      current = _mm_loadu_si128((__m128i *) high_ptr);
#endif
      result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
      _mm_storeu_si128((__m128i *) high_shifted,result);
      high_shifted -= 4; high_ptr -= 4;
    }
    this->high_shift_array[nshift][0] = this->high_blocks[0] << nshift;
    

    /* Low */
    low_shifted = &(this->low_shift_array[nshift][this->nwords - 3]);
    low_ptr = &(this->low_blocks[this->nwords - 3]);

    /* First 128-bit register has 0 bytes in the high word of current */
    prev = _mm_load_si128((__m128i *) (low_ptr - 1));
    current = _mm_srli_si128(prev,4);
    result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
    _mm_storeu_si128((__m128i *) low_shifted,result);
    low_shifted -= 4; low_ptr -= 4;

    while (low_ptr > this->low_blocks) {
#ifdef HAVE_SSSE3
      old_prev = prev;
      prev = _mm_load_si128((__m128i *) (low_ptr - 1));
      current = _mm_alignr_epi8(old_prev,prev,4);
#else 
      prev = _mm_load_si128((__m128i *) (low_ptr - 1));
      current = _mm_loadu_si128((__m128i *) low_ptr);
#endif
      result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
      _mm_storeu_si128((__m128i *) low_shifted,result);
      low_shifted -= 4; low_ptr -= 4;
    }
    this->low_shift_array[nshift][0] = this->low_blocks[0] << nshift;
#endif	/* HAVE_AVX2 and HAVE_MM256_LOADSTORE_U2 */
    
    if (this->flag_set_p == false) {
      /* Point the array to the blocks, which are all zero */
      this->flags_shift_array[nshift] = &(this->flags_blocks[0]);
    } else {
      /* Flags */
      flags_shifted = &(this->flags_shift_array[nshift][this->nwords - 3]);
      flags_ptr = &(this->flags_blocks[this->nwords - 3]);

      /* First 128-bit register has 0 bytes in the high word of current */
      prev = _mm_load_si128((__m128i *) (flags_ptr - 1));
      current = _mm_srli_si128(prev,4);
      result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
      _mm_storeu_si128((__m128i *) flags_shifted,result);
      flags_shifted -= 4; flags_ptr -= 4;

      while (flags_ptr > this->flags_blocks) {
#ifdef HAVE_SSSE3
	old_prev = prev;
	prev = _mm_load_si128((__m128i *) (flags_ptr - 1));
	current = _mm_alignr_epi8(old_prev,prev,4);
#else 
	prev = _mm_load_si128((__m128i *) (flags_ptr - 1));
	current = _mm_loadu_si128((__m128i *) flags_ptr);
#endif
	result = _mm_or_si128(_mm_slli_epi32(current,nshift),_mm_srli_epi32(prev,rightshift));
	_mm_storeu_si128((__m128i *) flags_shifted,result);
	flags_shifted -= 4; flags_ptr -= 4;
      }
      this->flags_shift_array[nshift][0] = this->flags_blocks[0] << nshift;
    }
  
    this->availp[nshift] = true;
  }

#if 0
  printf("nshift %d\n",nshift);
  for (int wordi = 0; wordi <= this->nwords; wordi++) {
    printf("%08X %08X %08X\n",
	   this->high_shift_array[nshift][wordi],
	   this->low_shift_array[nshift][wordi],
	   this->flags_shift_array[nshift][wordi]);
  }
  exit(9);
#endif

  initshift = (nshift + initpos)/32;

  *query_high_shifted = &(this->high_shift_array[nshift][initshift]);
  *query_low_shifted = &(this->low_shift_array[nshift][initshift]);
  *query_flags_shifted = &(this->flags_shift_array[nshift][initshift]);

  return;
}
#endif


bool
Compress_non_acgt (T this) {
  return this->flag_set_p;
}


bool
Compress_fwdp (T this) {
  return this->fwdp;
}

