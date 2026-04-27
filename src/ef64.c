/* Code modified from Sux package by Sebastiano Vigna,
   file bits/EliasFano.hpp */

#include "ef64.h"

#include <stdio.h>
#include "assert.h"
#include "mem.h"

#include "select64-common.h"
#include "select64-ones.h"
#include "select64-zeroes.h"


/* EF64_chrnum results */
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

/* Storing positions */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Check correctness of EF_two_ranks */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


static inline void
set (uint64_t *bits, const uint64_t pos) {
  bits[pos / 64] |= 1ULL << pos % 64;

  return;
}

static inline uint64_t
get_bits (uint64_t *bits, const uint64_t start, const int width) {
  int start_word = start / 64;
  int start_bit = start % 64;
  int total_offset = start_bit + width;
  uint64_t result = bits[start_word] >> start_bit;
  return (total_offset <= 64 ? result : result | bits[start_word + 1] << (64 - start_bit)) & ((1ULL << width) - 1);
}

static inline void
set_bits (uint64_t *bits, const uint64_t start, const int width, const uint64_t value) {
  uint64_t start_word = start / 64;
  uint64_t end_word = (start + width - 1) / 64;
  uint64_t start_bit = start % 64;

  if (start_word == end_word) {
    bits[start_word] &= ~(((1ULL << width) - 1) << start_bit);
    bits[start_word] |= value << start_bit;
  } else {
    /* Here start_bit > 0 */
    bits[start_word] &= (1ULL << start_bit) - 1;
    bits[start_word] |= value << start_bit;
    bits[end_word] &= -(1ULL << (width - 64 + start_bit));
    bits[end_word] |= value >> (64 - start_bit);
  }

  return;
}



#define T EF64_T
struct T {
  uint64_t *lower_bits;
  uint64_t *upper_bits;
  Select64_ones_T select_upper;
  Select64_zeroes_T selectz_upper;

  uint64_t num_ones;
  uint64_t nbits;
  /* Highest position observed, needed in EF64_rank if given position is higher than this */
  uint64_t last_one;

  int L;
  int block_size;
  int block_length;
  uint64_t block_size_mask;
  uint64_t lower_L_bits_mask;
};


void
EF64_free (T *old) {
  if (*old) {
    if ((*old)->selectz_upper != NULL) {
      Select64_zeroes_free(&(*old)->selectz_upper);
    }
    if ((*old)->select_upper != NULL) {
      Select64_ones_free(&(*old)->select_upper);
    }
    FREE((*old)->upper_bits);
    FREE((*old)->lower_bits);
    FREE(*old);
  }

  return;
}
  
  
/** Creates a new instance using an
 *  explicit list of positions for the ones in a bit vector.
 *
 *  Note that the list is read only at construction time.
 *
 *  In practice this constructor builds an Elias-Fano
 *  representation of the given list. select(const uint64_t rank) will retrieve
 *  an element of the list, and rank(const size_t pos) will return how many
 *  element of the list are smaller than the argument.
 */


static uint64_t
power (int base, int exponent) {
  uint64_t result = 1;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}

/* For repetitive oligos */
T
EF64_new_from_oligos (Oligospace_T *repetitive_oligos, uint64_t nrepetitive, int oligosize) {
  T new = (T) MALLOC(sizeof(*new));

  uint64_t oligo;

  uint64_t lower_bits_mask;
  int block_size;
  int L;
  uint64_t i;

  
  new->num_ones = nrepetitive;
  new->nbits = power(4,oligosize);

  if (new->num_ones == 0) {
    L = 0;
  } else if ((L = lambda_safe(new->nbits / new->num_ones)) < 0) {
    L = 0;
  }
  new->L = L;

#ifdef DEBUG
  printf("Number of ones: %lld.  L: %d\n", new->num_ones, L);
  printf("nbits: %u\n",new->nbits);
  printf("Upper bits: %lld\n", new->num_ones + (new->nbits >> L) + 1);
  printf("Lower bits: %lld\n", new->num_ones * L);
#endif

  lower_bits_mask = (1ULL << L) - 1;
  
  new->lower_bits = (uint64_t *) CALLOC(((new->num_ones * L + 63) / 64 + 2 * (L == 0)), sizeof(uint64_t));
  new->upper_bits = (uint64_t *) CALLOC((((new->num_ones + (new->nbits >> L) + 1) + 63) / 64), sizeof(uint64_t));

  for (i = 0; i < nrepetitive; i++) {
    oligo = (uint64_t) repetitive_oligos[i];
    debug2(printf("EF storing oligo %lu %016lx = %016lx %016lx\n",oligo,oligo,oligo>>L,oligo & lower_bits_mask));
    if (L != 0) {
      set_bits(new->lower_bits, /*start*/i * L, /*width*/L, /*value*/oligo & lower_bits_mask);
    }
    set(new->upper_bits, (oligo >> L) + i);
  }
  new->last_one = oligo;

#ifdef DEBUG
  printf("First lower: %016llx %016llx %016llx %016llx\n",
	 new->lower_bits[0], new->lower_bits[1], new->lower_bits[2], new->lower_bits[3]);
  for (i = 0; i < ((new->num_ones + (new->nbits >> L) + 1) + 63) / 64; i++) {
    printf("%d %016llx\n",i,new->upper_bits[i]);
  }
#endif

  if (new->num_ones == 0) {
    /* Otherwise, if num_ones == 0, these are time-consuming */
    new->select_upper = (Select64_ones_T) NULL;
    new->selectz_upper = (Select64_zeroes_T) NULL;
  } else {
    new->select_upper = Select64_ones_new(new->upper_bits, new->num_ones + (new->nbits >> L) + 1); /* Added + 1 from original code */
    new->selectz_upper = Select64_zeroes_new(new->upper_bits, new->num_ones + (new->nbits >> L) + 1); /* Added + 1 from original code */
  }

  block_size = 0;
  do
    ++block_size;
  while (block_size * L + block_size <= 64 && block_size <= L);
  block_size--;

  new->block_size = block_size;

#ifdef DEBUG
  printf("Block size: %d\n", block_size);
#endif
  new->block_size_mask = (1ULL << block_size) - 1;
  new->block_length = block_size * L;

  new->lower_L_bits_mask = (1ULL << L) - 1;

  return new;
}


/* positions: an array of positions of the ones in a bit vector, in ascending order */
/* nbits: the length (in bits) of the bit vector */
T
EF64_new_from_univcoordlist (Univcoordlist_T positions, const uint64_t nbits) {
  T new = (T) MALLOC(sizeof(*new));

  Univcoordlist_T p;
  uint64_t position;

  uint64_t lower_bits_mask;
  int block_size;
  int L;
  uint64_t i;

  
  new->num_ones = Univcoordlist_length(positions); /* num_ones */
  new->nbits = nbits;

  if (new->num_ones == 0) {
    L = 0;
  } else if ((L = lambda_safe(new->nbits / new->num_ones)) < 0) {
    L = 0;
  }
  new->L = L;

#ifdef DEBUG
  printf("Number of ones: %lld.  L: %d\n", new->num_ones, L);
  printf("nbits: %u\n",nbits);
  printf("Upper bits: %lld\n", new->num_ones + (nbits >> L) + 1);
  printf("Lower bits: %lld\n", new->num_ones * L);
#endif

  lower_bits_mask = (1ULL << L) - 1;
  
  new->lower_bits = (uint64_t *) CALLOC(((new->num_ones * L + 63) / 64 + 2 * (L == 0)), sizeof(uint64_t));
  new->upper_bits = (uint64_t *) CALLOC((((new->num_ones + (nbits >> L) + 1) + 63) / 64), sizeof(uint64_t));

  position = (uint64_t) 0;
  for (p = positions, i = 0; p != NULL; p = Univcoordlist_next(p), i++) {
    position = (uint64_t) Univcoordlist_head(p);
    debug2(printf("EF storing position %lu %016lx = %016lx %016lx\n",position,position,position>>L,position & lower_bits_mask));
    if (L != 0) {
      set_bits(new->lower_bits, /*start*/i * L, /*width*/L, /*value*/position & lower_bits_mask);
    }
    set(new->upper_bits, (position >> L) + i);
  }
  new->last_one = position;

#ifdef DEBUG
  printf("First lower: %016llx %016llx %016llx %016llx\n",
	 new->lower_bits[0], new->lower_bits[1], new->lower_bits[2], new->lower_bits[3]);
  for (i = 0; i < ((new->num_ones + (nbits >> L) + 1) + 63) / 64; i++) {
    printf("%d %016llx\n",i,new->upper_bits[i]);
  }
#endif

  if (new->num_ones == 0) {
    /* Otherwise, if num_ones == 0, these are time-consuming */
    new->select_upper = (Select64_ones_T) NULL;
    new->selectz_upper = (Select64_zeroes_T) NULL;
  } else {
    new->select_upper = Select64_ones_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code */
    new->selectz_upper = Select64_zeroes_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code */
  }

  block_size = 0;
  do
    ++block_size;
  while (block_size * L + block_size <= 64 && block_size <= L);
  block_size--;

  new->block_size = block_size;

#ifdef DEBUG
  printf("Block size: %d\n", block_size);
#endif
  new->block_size_mask = (1ULL << block_size) - 1;
  new->block_length = block_size * L;

  new->lower_L_bits_mask = (1ULL << L) - 1;

  return new;
}


T
EF64_new_from_chromosome_iit (Univ_IIT_T chromosome_iit) {
  T new = (T) MALLOC(sizeof(*new));

  uint64_t position;

  int nintervals, intervali;
  uint64_t lower_bits_mask;
  int block_size;
  int L;
  uint64_t nbits, i;


  new->num_ones = nintervals = Univ_IIT_total_nintervals(chromosome_iit);
  new->nbits = nbits = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);   /* genomelength */

  if (new->num_ones == 0) {
    L = 0;
  } else if ((L = lambda_safe(new->nbits / new->num_ones)) < 0) {
    L = 0;
  }
  new->L = L;

#ifdef DEBUG
  printf("Number of ones: %lld.  L: %d\n", new->num_ones, L);
  printf("nbits: %u\n",nbits);
  printf("Upper bits: %lld\n", new->num_ones + (nbits >> L) + 1);
  printf("Lower bits: %lld\n", new->num_ones * L);
#endif

  lower_bits_mask = (1ULL << L) - 1;
  
  new->lower_bits = (uint64_t *) CALLOC(((new->num_ones * L + 63) / 64 + 2 * (L == 0)), sizeof(uint64_t));
  new->upper_bits = (uint64_t *) CALLOC((((new->num_ones + (nbits >> L) + 1) + 63) / 64), sizeof(uint64_t));

  position = (uint64_t) 0;
  for (intervali = 0, i = 0; intervali < nintervals; intervali++, i++) {
    position = (uint64_t) Univ_IIT_interval_low(chromosome_iit,/*index*/intervali+1);
    debug2(printf("EF storing position %lu %016lx = %016lx %016lx\n",position,position,position>>L,position & lower_bits_mask));
    if (L != 0) {
      set_bits(new->lower_bits, /*start*/i * L, /*width*/L, /*value*/position & lower_bits_mask);
    }
    set(new->upper_bits, (position >> L) + i);
  }
  new->last_one = position;

#ifdef DEBUG
  printf("First lower: %016llx %016llx %016llx %016llx\n",
	 new->lower_bits[0], new->lower_bits[1], new->lower_bits[2], new->lower_bits[3]);
  for (i = 0; i < ((new->num_ones + (nbits >> L) + 1) + 63) / 64; i++) {
    printf("%d %016llx\n",i,new->upper_bits[i]);
  }
#endif

  if (new->num_ones == 0) {
    /* Otherwise, if num_ones == 0, these are time-consuming */
    new->select_upper = (Select64_ones_T) NULL;
    new->selectz_upper = (Select64_zeroes_T) NULL;
  } else {
    new->select_upper = Select64_ones_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code*/
    new->selectz_upper = Select64_zeroes_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code */
  }

  block_size = 0;
  do
    ++block_size;
  while (block_size * L + block_size <= 64 && block_size <= L);
  block_size--;

  new->block_size = block_size;

#ifdef DEBUG
  printf("Block size: %d\n", block_size);
#endif
  new->block_size_mask = (1ULL << block_size) - 1;
  new->block_length = block_size * L;

  new->lower_L_bits_mask = (1ULL << L) - 1;

  return new;
}


T
EF64_new_from_interleaved_univcoords (Univcoord_T *endpoints, uint64_t nintervals, const uint64_t nbits) {
  T new = (T) MALLOC(sizeof(*new));

  uint64_t position;

  uint64_t lower_bits_mask;
  int block_size;
  int L;
  uint64_t i, k;

  
  new->num_ones = nintervals;
  new->nbits = nbits;

  if (new->num_ones == 0) {
    L = 0;
  } else if ((L = lambda_safe(new->nbits / new->num_ones)) < 0) {
    L = 0;
  }
  new->L = L;

#ifdef DEBUG
  printf("Number of ones: %lld.  L: %d\n", new->num_ones, L);
  printf("nbits: %u\n",nbits);
  printf("Upper bits: %lld\n", new->num_ones + (nbits >> L) + 1);
  printf("Lower bits: %lld\n", new->num_ones * L);
#endif

  lower_bits_mask = (1ULL << L) - 1;
  
  new->lower_bits = (uint64_t *) CALLOC(((new->num_ones * L + 63) / 64 + 2 * (L == 0)), sizeof(uint64_t));
  new->upper_bits = (uint64_t *) CALLOC((((new->num_ones + (nbits >> L) + 1) + 63) / 64), sizeof(uint64_t));

  /* endpoints are interleaved start position and end position, so jump by 2, and convert to uint8 if necessary */
  position = (uint64_t) 0;
  for (i = 0, k = 0; i < nintervals; i++, k += 2) {
    position = (uint64_t) endpoints[k];
    debug2(printf("EF storing position %lu %016lx = %016lx %016lx\n",position,position,position>>L,position & lower_bits_mask));
    if (L != 0) {
      set_bits(new->lower_bits, /*start*/i * L, /*width*/L, /*value*/position & lower_bits_mask);
    }
    set(new->upper_bits, (position >> L) + i);
  }
  new->last_one = position;


#ifdef DEBUG
  printf("First lower: %016llx %016llx %016llx %016llx\n",
	 new->lower_bits[0], new->lower_bits[1], new->lower_bits[2], new->lower_bits[3]);
  printf("First upper: %016llx %016llx %016llx %016llx\n",
	 new->upper_bits[0], new->upper_bits[1], new->upper_bits[2], new->upper_bits[3]);
#endif

  if (new->num_ones == 0) {
    /* Otherwise, if num_ones == 0, these are time-consuming */
    new->select_upper = (Select64_ones_T) NULL;
    new->selectz_upper = (Select64_zeroes_T) NULL;
  } else {
    new->select_upper = Select64_ones_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code */
    new->selectz_upper = Select64_zeroes_new(new->upper_bits, new->num_ones + (nbits >> L) + 1); /* Added + 1 from original code */
  }

  block_size = 0;
  do
    ++block_size;
  while (block_size * L + block_size <= 64 && block_size <= L);
  block_size--;

  new->block_size = block_size;

#ifdef DEBUG
  printf("Block size: %d\n", block_size);
#endif
  new->block_size_mask = (1ULL << block_size) - 1;
  new->block_length = block_size * L;

  new->lower_L_bits_mask = (1ULL << L) - 1;

  return new;
}


uint64_t
EF64_rank (T this, const size_t k) {
  debug(printf("Entered EF64_rank with position %lu (%016lx %016lx)...\n",
	       k,k >> this->L, k & this->lower_L_bits_mask));

  if (this->num_ones == 0) {
    debug(printf("num_ones is 0, so returning 0\n"));
    return 0;
#if 1
  } else if (k >= this->nbits) { /* Was k >= this->nbits, but need to compare against this->last_one */
    debug(printf("value exceeds nbits %lu, so returning num_ones %lu\n",this->nbits,this->num_ones));
    return this->num_ones;
#else
  } else if (k >= this->last_one) {
    debug(printf("value exceeds last_one %lu, so returning num_ones %lu\n",this->last_one,this->num_ones));
    return this->num_ones;
#endif
  }
      
  const uint64_t k_shiftr_L = k >> this->L;

  int64_t pos = Select64_zeroes_select(this->selectz_upper,k_shiftr_L);
  uint64_t rank = pos - (k_shiftr_L);
  debug(printf("select zeroes yields position %ld and rank: %llu\n", pos, rank));

  /* Scan backward */
  uint64_t rank_times_L = rank * this->L;
  const uint64_t k_lower_bits = k & this->lower_L_bits_mask;

  do {
    rank--;
    rank_times_L -= this->L;
    pos--;
    debug(printf("At rank %d, upper bit %lu, scanning backward: comparing %016lx with %016lx\n",
		 rank,(this->upper_bits[pos / 64] & 1ULL << pos % 64),
		 get_bits(this->lower_bits,/*start*/rank_times_L,/*width*/this->L),k_lower_bits));
  } while (pos >= 0 && (this->upper_bits[pos / 64] & 1ULL << pos % 64) &&
	   get_bits(this->lower_bits, /*start*/rank_times_L, /*width*/this->L) > k_lower_bits);
  /* was >=, but that was giving chrhigh for a given position */

  debug(printf("EF64_rank returning %u\n",rank + 1));
  return ++rank;
}


/* We generally want ranks of two positions close to each other in the genome */
void
EF64_two_ranks (uint64_t *rank1, uint64_t *rank2, T this, const size_t k1, const size_t k2) {
  int64_t pos;
  uint64_t rank;
  uint64_t rank_times_L;

  debug(printf("Ranking %lu (%016lx) and %lu (%016lx)...\n", k1, k1, k2, k2));
  assert(k1 <= k2);

  if (this->num_ones == 0) {
    debug(printf("num_ones is zero\n"));
    *rank1 = *rank2 = 0;
  } else {
    const uint64_t k1_shiftr_L = k1 >> this->L;
    const uint64_t k2_shiftr_L = k2 >> this->L;
    const uint64_t k1_lower_bits = k1 & this->lower_L_bits_mask;
    const uint64_t k2_lower_bits = k2 & this->lower_L_bits_mask;

    if (k1_shiftr_L != k2_shiftr_L) {
      /* In different blocks, so need to have two separate calls for select */
      /* Find block for k2 */
      if (k2 >= this->last_one) { /* was k2 >= this->nbits */
	*rank2 = this->num_ones;
      } else {
	pos = Select64_zeroes_select(this->selectz_upper,k2_shiftr_L);
	rank = pos - (k2_shiftr_L);
	rank_times_L = rank * this->L;

	/* Scan backward for rank2 */
	do {
	  rank--;
	  rank_times_L -= this->L;
	  pos--;
	  debug(printf("At rank %d, upper bit %lu, comparing %016lx with %016lx\n",
		       rank,(this->upper_bits[pos / 64] & 1ULL << pos % 64),
		       get_bits(this->lower_bits,/*start*/rank_times_L,/*width*/this->L),k2_lower_bits));
	} while (pos >= 0 && (this->upper_bits[pos / 64] & 1ULL << pos % 64) &&
		 get_bits(this->lower_bits, /*start*/rank_times_L, /*width*/this->L) > k2_lower_bits);
	/* was >=, but changed to > to match with EF64_rank */
    
	*rank2 = rank + 1;
	debug(printf("Different blocks, got rank2 of %u\n",*rank2));
	debug15(assert(*rank2 == EF64_rank(this,k2)));
      }

      if (k1 >= this->last_one) { /* was k1 >= this->nbits */
	*rank1 = this->num_ones;
      } else {
	/* Find block for k1 */
	pos = Select64_zeroes_select(this->selectz_upper,k1_shiftr_L);
	rank = pos - (k1_shiftr_L);
	rank_times_L = rank * this->L;

	/* Scan backward for rank1 */
	do {
	  rank--;
	  rank_times_L -= this->L;
	  pos--;
	  debug(printf("At rank %d, upper bit %lu, comparing %016lx with %016lx\n",
		       rank,(this->upper_bits[pos / 64] & 1ULL << pos % 64),
		       get_bits(this->lower_bits,/*start*/rank_times_L,/*width*/this->L),k1_lower_bits));
	} while (pos >= 0 && (this->upper_bits[pos / 64] & 1ULL << pos % 64) &&
		 get_bits(this->lower_bits, /*start*/rank_times_L, /*width*/this->L) > k1_lower_bits);
	/* was >=, but changed to > to match with EF64_rank */
      
	*rank1 = rank + 1;
	debug(printf("Different blocks, got rank1 of %u\n",*rank1));
	debug15(assert(*rank1 == EF64_rank(this,k1)));
      }

    } else {
      /* In the same block, so can save and re-use the result (pos_orig, rank_orig) of the select command */
      /* Find block for k1 and k2 */
      pos = Select64_zeroes_select(this->selectz_upper,k2_shiftr_L);
      rank = pos - (k2_shiftr_L);
      rank_times_L = rank * this->L;

      if (k2 >= this->last_one) { /* was k2 >= this->nbits */
	*rank2 = this->num_ones;
      } else {
	/* Scan backward for rank2 */
	do {
	  rank--;
	  rank_times_L -= this->L;
	  pos--;
	  debug(printf("At rank %d, upper bit %lu, comparing %016lx with %016lx\n",
		       rank,(this->upper_bits[pos / 64] & 1ULL << pos % 64),
		       get_bits(this->lower_bits,/*start*/rank_times_L,/*width*/this->L),k2_lower_bits));
	} while (pos >= 0 && (this->upper_bits[pos / 64] & 1ULL << pos % 64) &&
		 get_bits(this->lower_bits, /*start*/rank_times_L, /*width*/this->L) > k2_lower_bits);
	/* was >=, but changed to > to match with EF64_rank */
    
	*rank2 = ++rank;		/* Need to increment pos, rank and rank_times_L so we test it again for k1 */
	rank_times_L += this->L;
	++pos;

	debug(printf("Same block, got rank2 of %u\n",*rank2));
	debug15(assert(*rank2 == EF64_rank(this,k2)));
      }

      if (k1 >= this->last_one) { /* was k1 >= this->nbits */
	*rank1 = this->num_ones;
      } else {
	/* Scan backward for rank1 */
	do {
	  rank--;
	  rank_times_L -= this->L;
	  pos--;
	  debug(printf("At rank %d, upper bit %lu, comparing %016lx with %016lx\n",
		       rank,(this->upper_bits[pos / 64] & 1ULL << pos % 64),
		       get_bits(this->lower_bits,/*start*/rank_times_L,/*width*/this->L),k1_lower_bits));
	} while (pos >= 0 && (this->upper_bits[pos / 64] & 1ULL << pos % 64) &&
		 get_bits(this->lower_bits, /*start*/rank_times_L, /*width*/this->L) > k1_lower_bits);
	/* was >=, but changed to > to match with EF64_rank */
      
	*rank1 = rank + 1;
	debug(printf("Same block, got rank1 of %u\n",*rank1));
	debug15(assert(*rank1 == EF64_rank(this,k1)));
      }
    }
  }

  return;
}


size_t
EF64_select (T this, const uint64_t rank) {
  debug(printf("Selecting %lld...\n", rank));
#ifdef DEBUG
  printf("Returning %lld = %llx << %d | %llx\n",
	 (Select64_ones_select(this->select_upper,rank) - rank) << this->L | get_bits(this->lower_bits, rank * this->L, this->L),
	 Select64_ones_select(this->select_upper,rank) - rank, this->L, get_bits(this->lower_bits, rank * this->L, this->L));
#endif
  return (Select64_ones_select(this->select_upper,rank) - rank) << this->L | get_bits(this->lower_bits, rank * this->L, this->L);
}


#if 1
/* Cannot call this with rank being the last chromosome */
uint64_t
EF64_next (T this, const uint64_t rank, size_t *next) {
  uint64_t s, t;

  debug(printf("EF64_next: comparing rank %u with num_ones %u\n",
	       rank,this->num_ones));

  s = Select64_ones_next(this->select_upper, rank, &t) - rank;
  t -= rank + 1;

  debug(printf("Select64_ones_next returns s %u and t %u\n",s,t));
  debug(printf("Select64_ones_select would have returned %u\n",
	       Select64_ones_select(this->select_upper, rank) - rank));

  const uint64_t position = rank * this->L;


  *next = t << this->L | get_bits(this->lower_bits, position + this->L, this->L);
  return s << this->L | get_bits(this->lower_bits, position, this->L);
}
#else
uint64_t
EF64_next (T this, const uint64_t rank) {
  uint64_t t;
  /* s = */ Select64_ones_next(this->select_upper, rank, &t) /*- rank (want to process t, not s) */;
  t -= rank + 1;

  const uint64_t position = rank * this->L;
  return t << this->L | get_bits(this->lower_bits, position + this->L, this->L);
}
#endif


#if 0
Chrnum_T
EF64_chrnum_old (Univcoord_T *chroffset, Univcoord_T *chrhigh, T this, Univcoord_T position) {
  uint64_t rank;
  size_t next_position;

  /* previously added 1 to position to handle the case where the
     position == 0, but that gives the wrong chromosome in some cases.
     Then we had a special case for position 0.  Now, we have changed
     the test from >= to > in EF64_rank and EF64_two_ranks */

#if 0
  if (position == 0) {
    *chroffset = 0;
    *chrhigh = (Univcoord_T) EF64_select(this,/*rank*/1);
    rank = 1;
  }
#endif

  if ((rank = EF64_rank(this,(size_t) position)) == this->num_ones) {
    *chroffset = (Univcoord_T) EF64_select(this,rank - 1);
    *chrhigh = this->nbits;

  } else {
    *chroffset = (Univcoord_T) EF64_next(this,rank - 1,&next_position);
    *chrhigh = (Univcoord_T) next_position;
  }

  assert(position >= *chroffset);
  assert(position < *chrhigh);
  assert(rank <= this->num_ones);

  debug(printf("EF64_chrnum called with %u => chrnum %d (out of %lu), chroffset %u, chrhigh %u\n",
	       position,rank,this->num_ones,*chroffset,*chrhigh));

  return (Chrnum_T) rank;		/* chrnum is 1-based */
}
#endif


Trnum_T
EF64_trnum (Trcoord_T *troffset, Trcoord_T *trhigh, T this,
	    Trcoord_T low_trcoord, Trcoord_T high_trcoord) {
  /* Trcoord_T trdiagonal, int querylength, int qstart, int qend */
  uint64_t rank;
  size_t next_position;
  Trcoord_T middle_trcoord;

  assert(low_trcoord <= high_trcoord);
  middle_trcoord = low_trcoord + (high_trcoord - low_trcoord)/2;

  /* previously added 1 to position to handle the case where the
     position == 0, but that gives the wrong chromosome in some cases.
     Then we had a special case for position 0.  Now, we have changed
     the test from >= to > in EF64_rank and EF64_two_ranks */

  debug0(printf("Entered EF64_trnum with %u..%u => %u\n",
		low_trcoord,high_trcoord,middle_trcoord));

#if 0
  if (position == 0) {
    *troffset = 0;
    *trhigh = (Univcoord_T) EF64_select(this,/*rank*/1);
    rank = 1;
  }
#endif

  if ((rank = EF64_rank(this,(size_t) middle_trcoord)) == this->num_ones) {
    /* Last transcript */
    *troffset = (Univcoord_T) EF64_select(this,rank - 1);
    *trhigh = this->nbits;

  } else {
    assert(rank < this->num_ones);

    *troffset = (Univcoord_T) EF64_next(this,rank - 1,&next_position);
    *trhigh = (Univcoord_T) next_position;

    assert(middle_trcoord >= *troffset);
    assert(middle_trcoord < *trhigh);

#if 0
    if (high_trcoord < *trhigh) {
      /* Both positions are within the same transcript */
    } else if (*trhigh - low_trcoord >= high_trcoord - *trhigh) {
      /* Most of sequence is in the lower transcript */
    } else {
      /* Advance to next transcript */
      *troffset = (Univcoord_T) EF64_next(this,rank,&next_position);
      *trhigh = (Univcoord_T) next_position;
    }
#endif
  }

  debug0(printf("EF64_trnum on %u returning trnum %d (out of %lu), troffset %u, trhigh %u\n",
		middle_trcoord,rank,this->num_ones,*troffset,*trhigh));

  return (Trnum_T) rank;		/* trnum is 1-based */
}


Chrnum_T
EF64_chrnum (Univcoord_T *chroffset, Univcoord_T *chrhigh, T this,
	     Univcoord_T low_univcoord, Univcoord_T high_univcoord) {
  uint64_t rank;
  size_t next_position;
  Univcoord_T middle_univcoord;

  assert(low_univcoord <= high_univcoord);
  middle_univcoord = low_univcoord + (high_univcoord - low_univcoord)/2;

  /* previously added 1 to position to handle the case where the
     position == 0, but that gives the wrong chromosome in some cases.
     Then we had a special case for position 0.  Now, we have changed
     the test from >= to > in EF64_rank and EF64_two_ranks */

  debug0(printf("Entered EF64_chrnum with %u..%u => %u\n",
		low_univcoord,high_univcoord,middle_univcoord));

#if 0
  if (position == 0) {
    *chroffset = 0;
    *chrhigh = (Univcoord_T) EF64_select(this,/*rank*/1);
    rank = 1;
  }
#endif

  if ((rank = EF64_rank(this,(size_t) middle_univcoord)) == this->num_ones) {
    /* Last chromosome */
    *chroffset = (Univcoord_T) EF64_select(this,rank - 1);
    *chrhigh = this->nbits;

  } else {
    assert(rank < this->num_ones);

    *chroffset = (Univcoord_T) EF64_next(this,rank - 1,&next_position);
    *chrhigh = (Univcoord_T) next_position;

    assert(middle_univcoord >= *chroffset);
    assert(middle_univcoord < *chrhigh);
  }

  debug0(printf("EF64_chrnum on %u returning chrnum %d (out of %lu), chroffset %u, chrhigh %u\n",
		middle_univcoord,rank,this->num_ones,*chroffset,*chrhigh));

  return (Chrnum_T) rank;		/* chrnum is 1-based */
}


void
EF64_chrbounds (Univcoord_T *chroffset, Univcoord_T *chrhigh, T this, Chrnum_T chrnum) {
  uint64_t rank;
  size_t next_position;

  rank = (uint64_t) chrnum;
  assert(rank > 0);

  if (rank == this->num_ones) {
    /* Last chromosome */
    *chroffset = (Univcoord_T) EF64_select(this,rank - 1);
    *chrhigh = this->nbits;
  } else {
    *chroffset = (Univcoord_T) EF64_next(this,rank - 1,&next_position);
    *chrhigh = (Univcoord_T) next_position;
  }

  return;
}


bool
EF64_presentp (size_t value, T this) {
  uint64_t rank;
  size_t closest;

  rank = EF64_rank(this,value);
  closest = EF64_select(this,rank - 1); /* Should be the one <= value */

  if (value == closest) {
    return true;
  } else {
    return false;
  }
}


#if 0
/* Faster to make two separate calls to EF64_rank */
void
EF64_ranks_in_range (uint64_t *low_rank, uint64_t *high_rank,
		     T this, const size_t low, const size_t high) {
  size_t position;

#if 0
  *low_rank = *high_rank = EF64_rank(this,low);
  debug(printf("Rank for %lu is %lu\n",low,*low_rank));
  if (*low_rank < this->num_ones && (position = EF64_select(this,*low_rank)) < high) {
    while (position < high && (*high_rank) + 1 < this->num_ones) {
      debug(printf("Rank %d in range because %u < %u\n",*high_rank,position,high));
      *high_rank += 1;
      position = EF64_select(this,*high_rank);
    }
  }
#else
  *low_rank = EF64_rank(this,low);
  *high_rank = EF64_rank(this,high);
#endif

  debug(printf("Range of ranks is %d..%d\n",*low_rank,*high_rank));

  return;
}
#endif




