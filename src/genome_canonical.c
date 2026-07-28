static char rcsid[] = "$Id: 53099d062cc97047832ce83fc0634689edfb4e3a $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genome_canonical.h"

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "except.h"
#include "maxent_hr.h"
#include "popcount.h"
#include "dinucl_bits.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif


/* Splice sites */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* last_dinucleotide_positions */
#ifdef DEBUG3
#define debug3(x) x
#else
#define debug3(x)
#endif

/* prev_dinucleotide_position */
#ifdef DEBUG3A
#define debug3a(x) x
#else
#define debug3a(x)
#endif

/* prev_dinucleotide_positions_rev */
#ifdef DEBUG3B
#define debug3b(x) x
#else
#define debug3b(x)
#endif

/* canonicalp */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif


#define T Genome_T

#define clear_start(diff,startdiscard) (diff & (~0U << (startdiscard)))
#define clear_end(diff,enddiscard) (diff & ~(~0U << (enddiscard)))


/************************************************************************
 *   Splice sites
 ************************************************************************/

static Genomecomp_T
block_find (Genomecomp_T *high_halfsite, Genomecomp_T *low_halfsite, Genomecomp_T *ptr, const Genomecomp_T *splicesite_bits) {
  Genomecomp_T found, compare, flags;

  /* High */
  debug2(printf("Evaluating high %08X and low %08X\n",*ptr,ptr[1]));
#ifdef WORDS_BIGENDIAN
  compare = Bigendian_convert_uint(*ptr++);
#else
  compare = *ptr++;
#endif
  /* Get high_halfsite bit */
  found = splicesite_bits[compare >> 16];
  *high_halfsite = (found & 0x100) >> 8;
  found = (found << 24) | 0x00FFFFFF;
  debug2(printf("  splicesite_bits_3: %08X, high_halfsite %d\n",
		splicesite_bits[compare >> 16] << 24,*high_halfsite));

  found &= (splicesite_bits[compare & 0x0000FFFF] << 16) | 0xFE00FFFF; /* Use FE to allow for high bit */
  debug2(printf("  splicesite_bits_2: %08X\n",splicesite_bits[compare & 0x0000FFFF] << 16));

  /* Low */
#ifdef WORDS_BIGENDIAN
  compare = Bigendian_convert_uint(*ptr++);
#else
  compare = *ptr++;
#endif
  found &= (splicesite_bits[compare >> 16] << 8) | 0xFFFE00FF;
  debug2(printf("  splicesite_bits_1: %08X\n",splicesite_bits[compare >> 16] << 8));

  found &= splicesite_bits[compare & 0x0000FFFF] | 0xFFFFFE00;
  debug2(printf("  splicesite_bits_0: %08X, low_halfsite %d\n",
		splicesite_bits[compare & 0x0000FFFF],found & 0x1));


  /* Flags: N is considered a mismatch */
  debug2(printf("Marking flags: genome %08X ",*ptr));
#ifdef WORDS_BIGENDIAN
  flags = Bigendian_convert_uint(*ptr);
#else
  flags = (*ptr);
#endif

  found &= ~flags;
  found &= ~(flags << 1);

  *low_halfsite = found & 0x00000001;

  *high_halfsite &= ~(flags >> 31);

  /* splicesite_bits are 1-based, so need to reduce values by 1 */
  /* No need to clear top bit */
  found >>= 1;
  debug2(printf(" => found %08X, high_halfsite %d, low_halfsite %d\n",found,*high_halfsite,*low_halfsite));

  return found;
}


static Genomecomp_T
block_find_snp (Genomecomp_T *high_halfsite, Genomecomp_T *low_halfsite, Genomecomp_T *ref_ptr, Genomecomp_T *alt_ptr,
		const Genomecomp_T *splicesite_bits) {
  Genomecomp_T found, ref, alt, flags;

  /* High */
  debug2(printf("Evaluating ref high %08X and low %08X, and alt high %08X and low %08X\n",
		*ref_ptr,ref_ptr[1],*alt_ptr,alt_ptr[1]));
#ifdef WORDS_BIGENDIAN
  ref = Bigendian_convert_uint(*ref_ptr++);
  alt = Bigendian_convert_uint(*alt_ptr++);
#else
  ref = *ref_ptr++;
  alt = *alt_ptr++;
#endif
  /* Get high_halfsite bit */
  found = splicesite_bits[ref >> 16] | splicesite_bits[alt >> 16];
  *high_halfsite = (found & 0x100) >> 8;
  found = (found << 24) | 0x00FFFFFF;
  debug2(printf("  splicesite_bits_3: %08X | %08X, high_halfsite %d\n",
		splicesite_bits[ref >> 16] << 24,splicesite_bits[alt >> 16] << 24,*high_halfsite));

 /* Use FE to allow for high bit */
  found &= ((splicesite_bits[ref & 0x0000FFFF] | splicesite_bits[alt & 0x0000FFFF]) << 16) | 0xFE00FFFF;
  debug2(printf("  splicesite_bits_2: %08X | %08X\n",
		splicesite_bits[ref & 0x0000FFFF] << 16,splicesite_bits[alt & 0x0000FFFF] << 16));

  /* Low */
#ifdef WORDS_BIGENDIAN
  ref = Bigendian_convert_uint(*ref_ptr++);
  alt = Bigendian_convert_uint(*alt_ptr++);
#else
  ref = *ref_ptr++;
  alt = *alt_ptr++;
#endif
  found &= ((splicesite_bits[ref >> 16] | splicesite_bits[alt >> 16]) << 8) | 0xFFFE00FF;
  debug2(printf("  splicesite_bits_1: %08X | %08X\n",
		splicesite_bits[ref >> 16] << 8,splicesite_bits[alt >> 16] << 8));

  found &= (splicesite_bits[ref & 0x0000FFFF] | splicesite_bits[alt & 0x0000FFFF]) | 0xFFFFFE00;
  debug2(printf("  splicesite_bits_0: %08X | %08X, low_halfsite %d\n",
		splicesite_bits[ref & 0x0000FFFF],splicesite_bits[ref & 0x0000FFFF],found & 0x1));


  /* Handle only reference flags, which indicate N */
  /* Reference flags: N is considered a mismatch */
  debug2(printf("Marking flags: genome %08X ",*ref_ptr));
#ifdef WORDS_BIGENDIAN
  flags = Bigendian_convert_uint(*ref_ptr);
#else
  flags = (*ref_ptr);
#endif

  found &= ~flags;
  found &= ~(flags << 1);

  *low_halfsite = found & 0x00000001;

  *high_halfsite &= ~(flags >> 31);

  /* splicesite_bits are 1-based, so need to reduce values by 1 */
  /* No need to clear top bit */
  found >>= 1;
  debug2(printf(" => found %08X, high_halfsite %d, low_halfsite %d\n",found,*high_halfsite,*low_halfsite));

  return found;
}



/* prevpos (startblock) corresponds to a lower genomicpos than pos */
static Univcoord_T
prev_dinucleotide_position (T genome, T genomealt, Univcoord_T pos, Univcoord_T prevpos,
#ifdef DEBUG3A
			    Univcoord_T chroffset,
#endif
			    const Genomecomp_T *splicesite_bits, int splicepos_offset) {
  Univcoord_T foundpos;
  Univcoord_T offset;
  int startdiscard, enddiscard;
  Univcoord_T startblocki, endblocki;
  Genomecomp_T *ptr, *altptr, *start;
  Genomecomp_T found;
  Genomecomp_T high_halfsite, low_halfsite, prev_low_halfsite;
  int relpos;
#ifndef HAVE_BUILTIN_CLZ
  Genomecomp_T top;
#endif

  startblocki = prevpos/32U*3;
  endblocki = pos/32U*3;
  enddiscard = pos % 32;

  debug3a(
	printf("\n\n");
	printf("Genome (in prev_dinucleotide_position): chroffset %llu, pos %llu, prevpos %llu\n",
	       (unsigned long long) chroffset,(unsigned long long) (pos-chroffset),(unsigned long long) (prevpos-chroffset));
	Genome_print_blocks_chrpos(genome->blocks,prevpos,pos,chroffset);
	printf("\n");
	);

  offset = (pos - 1) - enddiscard + 32 + splicepos_offset;

  if (endblocki == startblocki) {
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[endblocki]),&(genomealt->blocks[endblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[endblocki]),splicesite_bits);
    }
    found = clear_end(found,enddiscard);
    startdiscard = prevpos % 32;
    found = clear_start(found,startdiscard);

    if (found != 0U) {
#ifdef HAVE_BUILTIN_CLZ
      foundpos = offset - (relpos = __builtin_clz(found));
#else
      foundpos = offset - (relpos = (top = found >> 16) ? clz_table[top] : 16 + clz_table[found]);
#endif
      debug3a(printf("oneblock: found is %08X => offset %llu - relpos %d (%llu) => returning %llu\n",
		     found,(unsigned long long) offset,relpos,(unsigned long long) (foundpos-chroffset),
		     (unsigned long long) foundpos));
      return foundpos;
    } else {
      debug3a(printf("oneblock: not found\n"));
      return -1;
    }

  } else {
    /* Endblock */
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,
			     &(genome->blocks[endblocki]),&(genomealt->blocks[endblocki]),splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,&(genome->blocks[endblocki]),splicesite_bits);
    }
    found = clear_end(found,enddiscard);

    if (found != 0U) {
#ifdef HAVE_BUILTIN_CLZ
      foundpos = offset - (relpos = __builtin_clz(found));
#else
      foundpos = offset - (relpos = (top = found >> 16) ? clz_table[top] : 16 + clz_table[found]);
#endif
      debug3a(printf("endblock: found is %08X => offset %llu - relpos %d (%llu) => returning %llu\n",
		     found,(unsigned long long) offset,relpos,(unsigned long long) (foundpos-chroffset),
		     (unsigned long long) foundpos));
      return foundpos;
    }

    ptr = &(genome->blocks[endblocki-3]);
    altptr = &(genomealt->blocks[endblocki-3]);
    start = &(genome->blocks[startblocki]);
    offset -= 32;

    while (ptr > start) {
      /* Middle blocks */
      prev_low_halfsite = low_halfsite;
      if (genomealt != genome) {
	found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
      } else {
	found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
      }

      if (high_halfsite & prev_low_halfsite) {
	debug3a(printf("high_halfsite & prev_low_halfsite => offset %llu - 1 (%llu)\n",
		       (unsigned long long) offset,(unsigned long long) offset));
	return offset;
      } else if (found != 0U) {
#ifdef HAVE_BUILTIN_CLZ
	foundpos = offset - (relpos = __builtin_clz(found));
#else
	foundpos = offset - (relpos = (top = found >> 16) ? clz_table[top] : 16 + clz_table[found]);
#endif
	debug3a(printf("middleblock: found is %08X => offset %llu - relpos %d (%llu) => returning %llu\n",
		       found,(unsigned long long) offset,relpos,(unsigned long long) (foundpos-chroffset),
		       (unsigned long long) foundpos));
	return foundpos;
      }

      ptr -= 3;
      altptr -= 3;
      offset -= 32;
    }

    /* Startblock */
    prev_low_halfsite = low_halfsite;
    if (genomealt != genome) {
      found = block_find_snp(&high_halfsite,&low_halfsite,ptr,altptr,splicesite_bits);
    } else {
      found = block_find(&high_halfsite,&low_halfsite,ptr,splicesite_bits);
    }

    if (high_halfsite & prev_low_halfsite) {
      debug3a(printf("high_halfsite & prev_low_halfsite => offset %llu - 1 (%llu) => returning %llu\n",
		     (unsigned long long) offset,(unsigned long long) (offset-1-chroffset),(unsigned long long) offset));
      return offset;
    } else {
      startdiscard = prevpos % 32;
      found = clear_start(found,startdiscard);

      if (found != 0U) {
#ifdef HAVE_BUILTIN_CLZ
	foundpos = offset - (relpos = __builtin_clz(found));
#else
	foundpos = offset - (relpos = (top = found >> 16) ? clz_table[top] : 16 + clz_table[found]);
#endif
	debug3a(printf("startblock: found is %08X => offset %llu - relpos %d (%llu) => returning %llu\n",
		       found,(unsigned long long) offset,relpos,(unsigned long long) (foundpos-chroffset),
		       (unsigned long long) foundpos));
	return foundpos;
      } else {
	debug3a(printf("startblock: not found\n"));
	return -1;
      }
    }
  }
}



bool
Genome_sense_canonicalp (T genome, T genomealt,
			 Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			 Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			 Univcoord_T chroffset) {
  Univcoord_T donorpos, acceptorpos;
  Univcoord_T donor_shift, acceptor_shift;

  debug4(printf("Entered Genome_sense_canonicalp with donor %llu..%llu and acceptor %llu..%llu\n",
		(unsigned long long) (donor_leftbound-chroffset),
		(unsigned long long) (donor_rightbound-chroffset),
		(unsigned long long) (acceptor_leftbound-chroffset),
		(unsigned long long) (acceptor_rightbound-chroffset)));

  donorpos = prev_dinucleotide_position(genome,genomealt,donor_rightbound+1,donor_leftbound,
#ifdef DEBUG3A
					chroffset,
#endif
					donor_gtgc_bits,/*splicepos_offset*/-1);
  acceptorpos = prev_dinucleotide_position(genome,genomealt,acceptor_rightbound-1,acceptor_leftbound-2,
#ifdef DEBUG3A
					   chroffset,
#endif
					   acceptor_ag_bits,/*splicepos_offset*/2);

  while (donorpos != (Univcoord_T) -1 && acceptorpos != (Univcoord_T) -1) {
    donor_shift = donor_rightbound - (donorpos + 1U); /* Shift coordinates to match input */
    acceptor_shift = acceptor_rightbound - acceptorpos; 

    if (donor_shift < acceptor_shift) {
      donorpos = prev_dinucleotide_position(genome,genomealt,((donorpos + 1U)+1)-1,donor_leftbound,
#ifdef DEBUG3A
					    chroffset,
#endif
					    donor_gtgc_bits,/*splicepos_offset*/-1);
    } else if (acceptor_shift < donor_shift) {
      acceptorpos = prev_dinucleotide_position(genome,genomealt,(acceptorpos-1)-1,acceptor_leftbound-2,
#ifdef DEBUG3A
					       chroffset,
#endif
					       acceptor_ag_bits,/*splicepos_offset*/2);
    } else {
      debug4(printf("Found GT/GC-AG intron: donor at %llu (shift %llu) and acceptor at %llu (shift %llu)\n",
		    (unsigned long long) ((donorpos + 1U)-chroffset),(unsigned long long) donor_shift,
		    (unsigned long long) (acceptorpos-chroffset),(unsigned long long) acceptor_shift));

      debug4(printf("donor prob %f, acceptor prob %f\n",
		    Maxent_hr_donor_prob(genome,genomealt,donorpos + 1U,chroffset),
		    Maxent_hr_acceptor_prob(genome,genomealt,acceptorpos,chroffset)));
      if (Maxent_hr_donor_prob(genome,genomealt,donorpos + 1U,chroffset) > 0.9 &&
	  Maxent_hr_acceptor_prob(genome,genomealt,acceptorpos,chroffset) > 0.9) {
	return true;
      } else {
	donorpos = prev_dinucleotide_position(genome,genomealt,((donorpos + 1U)+1)-1,donor_leftbound,
#ifdef DEBUG3A
					      chroffset,
#endif
					      donor_gtgc_bits,/*splicepos_offset*/-1);
	acceptorpos = prev_dinucleotide_position(genome,genomealt,(acceptorpos-1)-1,acceptor_leftbound-2,
#ifdef DEBUG3A
						 chroffset,
#endif
						 acceptor_ag_bits,/*splicepos_offset*/2);
      }
    }
  }


  donorpos = prev_dinucleotide_position(genome,genomealt,donor_rightbound+1,donor_leftbound,
#ifdef DEBUG3A
					chroffset,
#endif
					donor_at_bits,/*splicepos_offset*/-1);
  acceptorpos = prev_dinucleotide_position(genome,genomealt,acceptor_rightbound-1,acceptor_leftbound-2,
#ifdef DEBUG3A
					   chroffset,
#endif
					   acceptor_ac_bits,/*splicepos_offset*/2);

  while (donorpos != (Univcoord_T) -1 && acceptorpos != (Univcoord_T) -1) {
    donor_shift = donor_rightbound - (donorpos + 1U); /* Shift coordinates to match input */
    acceptor_shift = acceptor_rightbound - acceptorpos; 

    if (donor_shift < acceptor_shift) {
      donorpos = prev_dinucleotide_position(genome,genomealt,((donorpos + 1U)+1)-1,donor_leftbound,
#ifdef DEBUG3A
					    chroffset,
#endif
					    donor_at_bits,/*splicepos_offset*/-1);
    } else if (acceptor_shift < donor_shift) {
      acceptorpos = prev_dinucleotide_position(genome,genomealt,(acceptorpos-1)-1,acceptor_leftbound-2,
#ifdef DEBUG3A
					       chroffset,
#endif
					       acceptor_ac_bits,/*splicepos_offset*/2);
    } else {
      debug4(printf("Found AT-AC intron: donor at %llu (shift %llu) and acceptor at %llu (shift %llu)\n",
		    (unsigned long long) ((donorpos + 1U)-chroffset),(unsigned long long) donor_shift,
		    (unsigned long long) (acceptorpos-chroffset),(unsigned long long) acceptor_shift));

      debug4(printf("donor prob %f, acceptor prob %f\n",
		    Maxent_hr_donor_prob(genome,genomealt,donorpos + 1U,chroffset),
		    Maxent_hr_acceptor_prob(genome,genomealt,acceptorpos,chroffset)));
      if (Maxent_hr_donor_prob(genome,genomealt,donorpos + 1U,chroffset) > 0.9 &&
	  Maxent_hr_acceptor_prob(genome,genomealt,acceptorpos,chroffset) > 0.9) {
	return true;
      } else {
	donorpos = prev_dinucleotide_position(genome,genomealt,((donorpos + 1U)+1)-1,donor_leftbound,
#ifdef DEBUG3A
					      chroffset,
#endif
					      donor_at_bits,/*splicepos_offset*/-1);
      acceptorpos = prev_dinucleotide_position(genome,genomealt,(acceptorpos-1)-1,acceptor_leftbound-2,
#ifdef DEBUG3A
					       chroffset,
#endif
					       acceptor_ac_bits,/*splicepos_offset*/2);
      }
    }
  }

  return false;
}


bool
Genome_antisense_canonicalp (T genome, T genomealt,
			     Univcoord_T donor_rightbound, Univcoord_T donor_leftbound,
			     Univcoord_T acceptor_rightbound, Univcoord_T acceptor_leftbound,
			     Univcoord_T chroffset) {
  Univcoord_T donorpos, acceptorpos;
  Univcoord_T donor_shift, acceptor_shift;

  debug4(printf("Entered Genome_antisense_canonicalp with donor %llu..%llu and acceptor %llu..%llu\n",
		(unsigned long long) (donor_leftbound-chroffset),
		(unsigned long long) (donor_rightbound-chroffset),
		(unsigned long long) (acceptor_leftbound-chroffset),
		(unsigned long long) (acceptor_rightbound-chroffset)));

  donorpos = prev_dinucleotide_position(genome,genomealt,donor_rightbound-1,donor_leftbound-2,
#ifdef DEBUG3A
					chroffset,
#endif
					antidonor_gtgc_bits,/*splicepos_offset*/2);

  acceptorpos = prev_dinucleotide_position(genome,genomealt,acceptor_rightbound+1,acceptor_leftbound,
#ifdef DEBUG3A
					   chroffset,
#endif
					   antiacceptor_ag_bits,/*splicepos_offset*/-1);
  
  while (donorpos != (Univcoord_T) -1 && acceptorpos != (Univcoord_T) -1) {
    donor_shift = donor_rightbound - donorpos;
    acceptor_shift = acceptor_rightbound - (acceptorpos + 1U); /* Shift coordinates to match input */

    if (donor_shift < acceptor_shift) {
      donorpos = prev_dinucleotide_position(genome,genomealt,(donorpos-1)-1,donor_leftbound-2,
#ifdef DEBUG3A
					    chroffset,
#endif
					    antidonor_gtgc_bits,/*splicepos_offset*/2);
    } else if (acceptor_shift < donor_shift) {
      acceptorpos = prev_dinucleotide_position(genome,genomealt,((acceptorpos + 1U)+1)-1,acceptor_leftbound,
#ifdef DEBUG3A
					       chroffset,
#endif
					       antiacceptor_ag_bits,/*splicepos_offset*/-1);
    } else {
      debug4(printf("antidonor prob %f, antiacceptor prob %f\n",
		    Maxent_hr_antidonor_prob(genome,genomealt,donorpos,chroffset),
		    Maxent_hr_antiacceptor_prob(genome,genomealt,acceptorpos + 1U,chroffset)));
      if (Maxent_hr_antidonor_prob(genome,genomealt,donorpos,chroffset) > 0.9 &&
	  Maxent_hr_antiacceptor_prob(genome,genomealt,acceptorpos + 1U,chroffset) > 0.9) {
	return true;
      } else {
	donorpos = prev_dinucleotide_position(genome,genomealt,(donorpos-1)-1,donor_leftbound-2,
#ifdef DEBUG3A
					      chroffset,
#endif
					      antidonor_gtgc_bits,/*splicepos_offset*/2);
	acceptorpos = prev_dinucleotide_position(genome,genomealt,((acceptorpos + 1U)+1)-1,acceptor_leftbound,
#ifdef DEBUG3A
						 chroffset,
#endif
						 antiacceptor_ag_bits,/*splicepos_offset*/-1);
      }
    }
  }


  donorpos = prev_dinucleotide_position(genome,genomealt,donor_rightbound-1,donor_leftbound-2,
#ifdef DEBUG3A
					chroffset,
#endif
					antidonor_at_bits,/*splicepos_offset*/2);

  acceptorpos = prev_dinucleotide_position(genome,genomealt,acceptor_rightbound+1,acceptor_leftbound,
#ifdef DEBUG3A
					   chroffset,
#endif
					   antiacceptor_ac_bits,/*splicepos_offset*/-1);
  
  while (donorpos != (Univcoord_T) -1 && acceptorpos != (Univcoord_T) -1) {
    donor_shift = donor_rightbound - donorpos;
    acceptor_shift = acceptor_rightbound - (acceptorpos + 1U); /* Shift coordinates to match input */

    if (donor_shift < acceptor_shift) {
      donorpos = prev_dinucleotide_position(genome,genomealt,(donorpos-1)-1,donor_leftbound-2,
#ifdef DEBUG3A
					    chroffset,
#endif
					    antidonor_at_bits,/*splicepos_offset*/2);
    } else if (acceptor_shift < donor_shift) {
      acceptorpos = prev_dinucleotide_position(genome,genomealt,((acceptorpos + 1U)+1)-1,acceptor_leftbound,
#ifdef DEBUG3A
					       chroffset,
#endif
					       antiacceptor_ac_bits,/*splicepos_offset*/-1);
    } else {
      debug4(printf("antidonor prob %f, antiacceptor prob %f\n",
		    Maxent_hr_antidonor_prob(genome,genomealt,donorpos,chroffset),
		    Maxent_hr_antiacceptor_prob(genome,genomealt,acceptorpos + 1U,chroffset)));
      if (Maxent_hr_antidonor_prob(genome,genomealt,donorpos,chroffset) > 0.9 &&
	  Maxent_hr_antiacceptor_prob(genome,genomealt,acceptorpos + 1U,chroffset) > 0.9) {
	return true;
      } else {
	donorpos = prev_dinucleotide_position(genome,genomealt,(donorpos-1)-1,donor_leftbound-2,
#ifdef DEBUG3A
					      chroffset,
#endif
					      antidonor_at_bits,/*splicepos_offset*/2);
	acceptorpos = prev_dinucleotide_position(genome,genomealt,((acceptorpos + 1U)+1)-1,acceptor_leftbound,
#ifdef DEBUG3A
						 chroffset,
#endif
						 antiacceptor_ac_bits,/*splicepos_offset*/-1);
      }
    }
  }

  return false;
}



