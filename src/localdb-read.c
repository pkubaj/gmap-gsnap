#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "localdb-read.h"
#include "debug_groups.h"
#include "accel-profile-stats.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/mman.h>		/* For munmap */

#include "assert.h"
#include "mem.h"
#include "genomicpos.h"
#include "univdiagdef.h"

#include "genomebits_consec.h"
#include "genomebits_count.h"
#include "genomebits_trim.h"
/* #include "merge-uint4.h" -- for get_by_bisection */
#include "intersect-uint2.h"

#include "simd.h"

#define EPI16_MAX 0x8000


/* Print univdiagonals */
#if defined(DEBUG0) || defined(DEBUG_LOCALDB_DEEP)
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Method used */
#if defined(DEBUG) || defined(DEBUG_LOCALDB)
#define debug(x) x
#else
#define debug(x)
#endif

/* Details of suffix array search */
#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
#define debug1(x) x
#else
#define debug1(x)
#endif

/* get_exhaustive */
#if defined(DEBUG2) || defined(DEBUG_LOCALDB_DEEP)
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Localdb_get_one */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Dump each suffix array contents */
#if defined(DEBUG9) || defined(DEBUG_LOCALDB_DEEP)
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Dump saindex */
#ifdef DEBUG9A
#define debug9a(x) x
#else
#define debug9a(x)
#endif

/* Compare specialized Localdb_get against a general method */
#ifdef DEBUG15
#define debug15(x) x
#else
#define debug15(x)
#endif


/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif


/* It is critical to have the "U" at the end of each of these values.
   Otherwise, a compiler sometimes optimizes to give strange results */
#define SAINDEX16_LENGTH 256U	/* 2^8 */
#define SARRAY16_LENGTH 65536U	/* 2^16 */
#define SARRAY8_LENGTH 256U	/* 2^8 */


#define OLIGOSIZE 4 /* nucleotides for saindex to specify 256 pointers */


static Genome_T genomecomp;
static Univcoord_T genomelength;
static int index1part;
static int index1interval;

void
Localdb_setup (Genome_T genomecomp_in, Univcoord_T genomelength_in,
	       int index1part_in, int index1interval_in) {
  genomecomp = genomecomp_in;
  genomelength = genomelength_in;
  index1part = index1part_in;
  index1interval = index1interval_in;

  return;
}


#define T Localdb_T
struct T {
  uint16_t *saindex16;
  uint16_t *sarray16;
  uint8_t *sarray8;
  uint16_t *sasort16;

  Access_T localdb_access;

  int saindex16_shmid;
  key_t saindex16_key;
  int saindex16_fd;
  size_t saindex16_len;

  int sarray16_shmid;
  key_t sarray16_key;
  int sarray16_fd;
  size_t sarray16_len;

  int sarray8_shmid;
  key_t sarray8_key;
  int sarray8_fd;
  size_t sarray8_len;

  int sasort16_shmid;
  key_t sasort16_key;
  int sasort16_fd;
  size_t sasort16_len;
};

#if defined(__GNUC__)
#define LOCALDB_CACHE_TLS __thread
#else
#define LOCALDB_CACHE_TLS _Thread_local
#endif

#define LOCALDB_CACHE_INITIAL_NSLOTS 32768
/* Table growth ceiling: at the cap the table holds LOCALDB_CACHE_MAX_NSLOTS
   entries (~ nslots * sizeof(entry) of index overhead per worker thread). */
#define LOCALDB_CACHE_MAX_NSLOTS (1 << 17)
/* Result-payload ceiling per worker thread.  Exceeding it evicts the coldest
   half of the cache (approximate LRU) rather than disabling storage, so the
   pathological tail keeps its hot working set instead of running uncached. */
#define LOCALDB_CACHE_MAX_BYTES (64 * 1024 * 1024)
/* Single results too large to be worth their footprint are not cached; a miss
   just recomputes the identical value. */
#define LOCALDB_CACHE_MAX_RESULTS_PER_ENTRY 8192

/* Memoization is a pure per-pair value cache (Localdb_get is a pure function of
   its behavioral key; localdb_tally/localdb_overflowp are scratch, not outputs).
   It is unconditional in production and disabled only for the legacy
   request-graph statistics build, which must count every logical execution. */
#ifdef GSNAP_STATS_LEGACY_REQUEST_GRAPH
static const bool localdb_memoize_enabled_p = false;
#else
static const bool localdb_memoize_enabled_p = true;
#endif

struct Localdb_cache_entry_T {
  bool occupiedp;
  unsigned long long hash;
  unsigned long long last_used;
  T localdb;
  Genomebits_T genomebits;
  Compress_T query_compress;
  char *queryptr;
  int pos5;
  int pos3;
  int querylength;
  Univcoord_T low_univdiagonal;
  Univcoord_T high_univdiagonal;
  bool plusp;
  int genestrand;
  int nhits;
  int matchlength;
  bool sortedp;
  bool trimmedp;
  struct Univdiag_T *results;
  size_t result_bytes;
};

static LOCALDB_CACHE_TLS struct Localdb_cache_entry_T *localdb_cache_entries = NULL;
static LOCALDB_CACHE_TLS int *localdb_cache_used_slots = NULL;
static LOCALDB_CACHE_TLS int localdb_cache_nslots = 0;
static LOCALDB_CACHE_TLS int localdb_cache_nused = 0;
static LOCALDB_CACHE_TLS int localdb_cache_nentries = 0;
static LOCALDB_CACHE_TLS size_t localdb_cache_bytes = 0U;
static LOCALDB_CACHE_TLS unsigned long long localdb_cache_tick = 0ULL;

static unsigned long long
localdb_cache_mix (unsigned long long hash, unsigned long long value) {
  int i;

  for (i = 0; i < 8; i++) {
    hash ^= (value & 0xffU);
    hash *= 1099511628211ULL;
    value >>= 8;
  }

  return hash;
}

static unsigned long long
localdb_cache_hash (T this, Genomebits_T genomebits, Compress_T query_compress, char *queryptr,
		    int pos5, int pos3, int querylength,
		    Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		    bool plusp, int genestrand) {
  unsigned long long hash = 1469598103934665603ULL;
  int pos;

  hash = localdb_cache_mix(hash,(unsigned long long) (uintptr_t) this);
  hash = localdb_cache_mix(hash,(unsigned long long) (uintptr_t) genomebits);
  hash = localdb_cache_mix(hash,(unsigned long long) (uintptr_t) query_compress);
  hash = localdb_cache_mix(hash,(unsigned long long) pos5);
  hash = localdb_cache_mix(hash,(unsigned long long) pos3);
  hash = localdb_cache_mix(hash,(unsigned long long) querylength);
  hash = localdb_cache_mix(hash,(unsigned long long) low_univdiagonal);
  hash = localdb_cache_mix(hash,(unsigned long long) high_univdiagonal);
  hash = localdb_cache_mix(hash,(unsigned long long) plusp);
  hash = localdb_cache_mix(hash,(unsigned long long) genestrand);

  if (queryptr != NULL) {
    for (pos = pos5; pos < pos3; pos++) {
      hash ^= (unsigned char) queryptr[pos];
      hash *= 1099511628211ULL;
    }
  }

  return hash;
}

static bool
localdb_cache_key_matches_p (struct Localdb_cache_entry_T *entry, unsigned long long hash,
			     T this, Genomebits_T genomebits, Compress_T query_compress, char *queryptr,
			     int pos5, int pos3, int querylength,
			     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
			     bool plusp, int genestrand) {
  int span = pos3 - pos5;

  if (entry->hash != hash ||
      entry->localdb != this ||
      entry->genomebits != genomebits ||
      entry->query_compress != query_compress ||
      entry->pos5 != pos5 ||
      entry->pos3 != pos3 ||
      entry->querylength != querylength ||
      entry->low_univdiagonal != low_univdiagonal ||
      entry->high_univdiagonal != high_univdiagonal ||
      entry->plusp != plusp ||
      entry->genestrand != genestrand) {
    return false;

  } else if (span <= 0) {
    return true;

  } else if (entry->queryptr == NULL || queryptr == NULL) {
    return entry->queryptr == queryptr ? true : false;

  } else if (memcmp(&(entry->queryptr[pos5]),&(queryptr[pos5]),span) == 0) {
    return true;

  } else {
    return false;
  }
}

static bool
localdb_cache_allocated_p (void) {
  if (localdb_cache_entries != NULL) {
    return true;
  }

  localdb_cache_nslots = LOCALDB_CACHE_INITIAL_NSLOTS;
  localdb_cache_entries =
    (struct Localdb_cache_entry_T *) CALLOC(localdb_cache_nslots,
					   sizeof(struct Localdb_cache_entry_T));
  localdb_cache_used_slots = (int *) MALLOC(localdb_cache_nslots * sizeof(int));
  return true;
}

static void
localdb_cache_grow (void) {
  struct Localdb_cache_entry_T *old_entries = localdb_cache_entries;
  int *old_used_slots = localdb_cache_used_slots;
  int old_nused = localdb_cache_nused;
  int old_nslots = localdb_cache_nslots;
  unsigned int slot, mask;
  int i;

  localdb_cache_nslots *= 2;
  mask = (unsigned int) localdb_cache_nslots - 1U;
  localdb_cache_entries =
    (struct Localdb_cache_entry_T *) CALLOC(localdb_cache_nslots,
					   sizeof(struct Localdb_cache_entry_T));
  localdb_cache_used_slots = (int *) MALLOC(localdb_cache_nslots * sizeof(int));
  localdb_cache_nused = 0;

  for (i = 0; i < old_nused; i++) {
    struct Localdb_cache_entry_T entry = old_entries[old_used_slots[i]];
    slot = (unsigned int) entry.hash & mask;
    while (localdb_cache_entries[slot].occupiedp == true) {
      slot = (slot + 1U) & mask;
    }
    localdb_cache_entries[slot] = entry;
    localdb_cache_used_slots[localdb_cache_nused++] = (int) slot;
  }

  FREE(old_used_slots);
  FREE(old_entries);
  (void) old_nslots;
}

static int
localdb_cache_tick_cmp (const void *a, const void *b) {
  unsigned long long x = *(const unsigned long long *) a;
  unsigned long long y = *(const unsigned long long *) b;

  return x < y ? -1 : (x > y ? 1 : 0);
}

/* Approximate-LRU eviction: drop the coldest ~half of the cache by last-used
   tick and rehash the survivors into a fresh same-size table.  A cache miss
   only recomputes the identical value, so eviction never changes alignment
   output; it only trades a bounded memory footprint for some recomputation on
   the pathological tail. */
static void
localdb_cache_evict_half (void) {
  struct Localdb_cache_entry_T *old_entries = localdb_cache_entries;
  int *old_used_slots = localdb_cache_used_slots;
  int old_nused = localdb_cache_nused;
  unsigned long long *ticks, threshold;
  unsigned int slot, mask;
  int i, keep_count;

  if (old_nused <= 0) {
    return;
  }

  keep_count = old_nused / 2;
  if (keep_count == 0) {
    /* Degenerate (defensive): nothing worth keeping, drop everything. */
    for (i = 0; i < old_nused; i++) {
      struct Localdb_cache_entry_T *entry = &(old_entries[old_used_slots[i]]);
      if (entry->results != NULL) {
	FREE(entry->results);
      }
      memset(entry,0,sizeof(*entry));
    }
    localdb_cache_nused = 0;
    localdb_cache_nentries = 0;
    localdb_cache_bytes = 0U;
    return;
  }

  ticks = (unsigned long long *) MALLOC(old_nused * sizeof(unsigned long long));
  for (i = 0; i < old_nused; i++) {
    ticks[i] = old_entries[old_used_slots[i]].last_used;
  }
  qsort(ticks,old_nused,sizeof(unsigned long long),localdb_cache_tick_cmp);
  threshold = ticks[old_nused - keep_count]; /* keep the freshest keep_count entries */
  FREE(ticks);

  mask = (unsigned int) localdb_cache_nslots - 1U;
  localdb_cache_entries =
    (struct Localdb_cache_entry_T *) CALLOC(localdb_cache_nslots,
					   sizeof(struct Localdb_cache_entry_T));
  localdb_cache_used_slots = (int *) MALLOC(localdb_cache_nslots * sizeof(int));
  localdb_cache_nused = 0;
  localdb_cache_nentries = 0;
  localdb_cache_bytes = 0U;

  for (i = 0; i < old_nused; i++) {
    struct Localdb_cache_entry_T entry = old_entries[old_used_slots[i]];
    if (entry.last_used >= threshold) {
      slot = (unsigned int) entry.hash & mask;
      while (localdb_cache_entries[slot].occupiedp == true) {
	slot = (slot + 1U) & mask;
      }
      localdb_cache_entries[slot] = entry; /* moves results ownership to survivor */
      localdb_cache_used_slots[localdb_cache_nused++] = (int) slot;
      localdb_cache_nentries++;
      localdb_cache_bytes += entry.result_bytes;
    } else if (entry.results != NULL) {
      FREE(entry.results);
    }
  }

  FREE(old_used_slots);
  FREE(old_entries);
}

void
Localdb_cache_reset (void) {
  struct Localdb_cache_entry_T *entry;
  int i, slot;

  if (localdb_cache_entries == NULL) {
    return;
  }

  for (i = 0; i < localdb_cache_nused; i++) {
    slot = localdb_cache_used_slots[i];
    entry = &(localdb_cache_entries[slot]);
    if (entry->results != NULL) {
      FREE(entry->results);
    }
    memset(entry,0,sizeof(*entry));
  }

  localdb_cache_nused = 0;
  localdb_cache_nentries = 0;
  localdb_cache_bytes = 0U;
  localdb_cache_tick = 0ULL;
  return;
}

void
Localdb_cache_free (void) {
  Localdb_cache_reset();
  if (localdb_cache_used_slots != NULL) {
    FREE(localdb_cache_used_slots);
  }
  if (localdb_cache_entries != NULL) {
    FREE(localdb_cache_entries);
  }
  localdb_cache_used_slots = NULL;
  localdb_cache_entries = NULL;
  localdb_cache_nslots = 0;
  localdb_cache_nused = 0;
  localdb_cache_nentries = 0;
  localdb_cache_bytes = 0U;
  localdb_cache_tick = 0ULL;
  return;
}

static struct Localdb_cache_entry_T *
localdb_cache_lookup (unsigned long long hash, T this, Genomebits_T genomebits,
		      Compress_T query_compress, char *queryptr, int pos5, int pos3, int querylength,
		      Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		      bool plusp, int genestrand) {
  struct Localdb_cache_entry_T *entry;
  unsigned int slot;
  int probe;

  if (localdb_cache_entries == NULL) {
    return (struct Localdb_cache_entry_T *) NULL;
  }

  slot = (unsigned int) hash & ((unsigned int) localdb_cache_nslots - 1U);
  for (probe = 0; probe < localdb_cache_nslots; probe++) {
    entry = &(localdb_cache_entries[slot]);
    if (entry->occupiedp == false) {
      return (struct Localdb_cache_entry_T *) NULL;
    } else if (localdb_cache_key_matches_p(entry,hash,this,genomebits,query_compress,queryptr,
					   pos5,pos3,querylength,
					   low_univdiagonal,high_univdiagonal,
					   plusp,genestrand) == true) {
      entry->last_used = ++localdb_cache_tick;
      return entry;
    } else {
      slot = (slot + 1U) & ((unsigned int) localdb_cache_nslots - 1U);
    }
  }

  return (struct Localdb_cache_entry_T *) NULL;
}

static void
localdb_cache_replay (bool *sortedp, bool *trimmedp, int *matchlength,
		      Univdiag_T *univdiags_alloc,
		      struct Localdb_cache_entry_T *entry) {
  int i;

  *sortedp = entry->sortedp;
  *trimmedp = entry->trimmedp;
  *matchlength = entry->matchlength;
  for (i = 0; i < entry->nhits; i++) {
    Univdiag_transfer(univdiags_alloc[i],&(entry->results[i]));
  }

  return;
}

static void
localdb_cache_insert (unsigned long long hash, T this, Genomebits_T genomebits,
		      Compress_T query_compress, char *queryptr, int pos5, int pos3, int querylength,
		      Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		      bool plusp, int genestrand,
		      int nhits, int matchlength, bool sortedp, bool trimmedp,
		      Univdiag_T *univdiags_alloc) {
  struct Localdb_cache_entry_T *entry;
  size_t result_bytes;
  unsigned int slot;
  int i, probe;

  if (nhits > LOCALDB_CACHE_MAX_RESULTS_PER_ENTRY) {
    /* Oversized single result: cheaper to recompute than to hold. */
    return;
  }
  result_bytes = (size_t) nhits * sizeof(struct Univdiag_T);
  if (localdb_cache_allocated_p() == false) {
    return;
  }

  /* Bound result-payload memory: evict the coldest half until the new entry
     fits under the byte ceiling. */
  while (localdb_cache_nentries > 0 &&
	 localdb_cache_bytes + result_bytes > LOCALDB_CACHE_MAX_BYTES) {
    localdb_cache_evict_half();
  }

  /* Bound table size: grow while under the slot ceiling, otherwise evict to
     free slots. */
  if ((localdb_cache_nentries + 1) * 10 >= localdb_cache_nslots * 7) {
    if (localdb_cache_nslots * 2 <= LOCALDB_CACHE_MAX_NSLOTS) {
      localdb_cache_grow();
    } else {
      localdb_cache_evict_half();
    }
  }

  slot = (unsigned int) hash & ((unsigned int) localdb_cache_nslots - 1U);
  for (probe = 0; probe < localdb_cache_nslots; probe++) {
    entry = &(localdb_cache_entries[slot]);
    if (entry->occupiedp == false) {
      entry->occupiedp = true;
      entry->hash = hash;
      entry->localdb = this;
      entry->genomebits = genomebits;
      entry->query_compress = query_compress;
      entry->queryptr = queryptr;
      entry->pos5 = pos5;
      entry->pos3 = pos3;
      entry->querylength = querylength;
      entry->low_univdiagonal = low_univdiagonal;
      entry->high_univdiagonal = high_univdiagonal;
      entry->plusp = plusp;
      entry->genestrand = genestrand;
      entry->nhits = nhits;
      entry->matchlength = matchlength;
      entry->sortedp = sortedp;
      entry->trimmedp = trimmedp;
      entry->result_bytes = result_bytes;
      entry->last_used = ++localdb_cache_tick;
      if (nhits > 0) {
	entry->results = (struct Univdiag_T *) MALLOC(result_bytes);
	for (i = 0; i < nhits; i++) {
	  Univdiag_transfer(&(entry->results[i]),univdiags_alloc[i]);
	}
      }
      localdb_cache_bytes += result_bytes;
      localdb_cache_used_slots[localdb_cache_nused++] = (int) slot;
      localdb_cache_nentries++;
      return;

    } else if (localdb_cache_key_matches_p(entry,hash,this,genomebits,query_compress,queryptr,
					   pos5,pos3,querylength,
					   low_univdiagonal,high_univdiagonal,
					   plusp,genestrand) == true) {
      return;

    } else {
      slot = (slot + 1U) & ((unsigned int) localdb_cache_nslots - 1U);
    }
  }

  return;
}


void
Localdb_free (T *old) {

  if ((*old)->localdb_access == ALLOCATED_PRIVATE) {
    FREE_KEEP((*old)->saindex16);
    FREE_KEEP((*old)->sarray16);
    FREE_KEEP((*old)->sarray8);
    FREE_KEEP((*old)->sasort16);

  } else if ((*old)->localdb_access == ALLOCATED_SHARED) {
    Access_deallocate((*old)->saindex16,(*old)->saindex16_shmid,(*old)->saindex16_key);
    Access_deallocate((*old)->sarray16,(*old)->sarray16_shmid,(*old)->sarray16_key);
    Access_deallocate((*old)->sarray8,(*old)->sarray8_shmid,(*old)->sarray8_key);
    Access_deallocate((*old)->sasort16,(*old)->sasort16_shmid,(*old)->sasort16_key);

#ifdef HAVE_MMAP
  } else if ((*old)->localdb_access == MMAPPED) {
    munmap((void *) (*old)->saindex16,(*old)->saindex16_len);
    munmap((void *) (*old)->sarray16,(*old)->sarray16_len);
    munmap((void *) (*old)->sarray8,(*old)->sarray8_len);
    munmap((void *) (*old)->sasort16,(*old)->sasort16_len);
    close((*old)->saindex16_fd);
    close((*old)->sarray16_fd);
    close((*old)->sarray8_fd);
    close((*old)->sasort16_fd);
#endif
  }
  
  FREE(*old);
}


#if 0
/* Position 10470: CGGTACCC */
/* TODO: Write get char that works on genomebits */

/* Taken from Johannes Fischer, Advanced Text Indexing Techniques, Algorithm 1 */
/* Does not use LCP, so time is O(m * log(n)) */
static void
sarray_search_uint1 (int *matchlength, int *initptr, int *finalptr, uint8_t *SA, int n,
		     char *query, int querylength, int pos5, int pos3, int min_matchlength,
		     Univcoord_T region_term, Compress_T query_compress, bool plusp,
		     Genomebits_T genomebits) {
  int low, high, mid;
  Univcoord_T pos;
  int nmatches = 0;
  int seglength = pos3 - pos5;
  char mismatch_char;

#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
  char *Buffer;
  Buffer = (char *) MALLOC((seglength+2)*sizeof(char));
#endif

  debug1(printf("Entered sarray_search_uint1\n"));

  low = 0;
  high = n;

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    /* low 0, high 2, mid = 1.  low 0, high 3, mid = 1.  low 1, high 2, mid = 1.  low 1, high 3, mid = 2 */
    mid = low + (high - low + 0)/2;

    pos = region_term + SA[mid];

#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
    printf("low %d, high %d, mid %d => position %u\n",low,high,mid,pos);
    Genome_fill_buffer_simple(genomecomp,/*left*/pos,seglength,Buffer);
    if (SA[mid] + seglength >= n) {
      Buffer[n - SA[mid]] = '$';
      Buffer[n - SA[mid] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",mid,SA[mid],Buffer);
    printf("query\t\t%.*s\n",seglength,&(query[pos5]));
#endif

    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    if (SA[mid] + nmatches >= n) {
      /* Terminating char is always less than query */
      debug1(printf("nmatches %d, genomic char $\n",nmatches));
      low = mid + 1;
    } else if (nmatches >= seglength) {
      debug1(printf("nmatches %d >= seglength\n",nmatches));
      high = mid;
    } else if (mismatch_char > query[pos5+nmatches]) {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c > query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      high = mid;
    } else {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c <= query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      low = mid + 1;
    }
  }

  if ((*initptr = low) == n) {
    /* Query is higher than any suffix */
    *matchlength = 0;
    *finalptr = *initptr;
  } else {
    pos = region_term + SA[*initptr];
    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    debug1(printf("Setting initptr to be %d with %d matches\n",low,nmatches));

    if (nmatches < seglength) {
      *matchlength = nmatches;
    } else {
      *matchlength = seglength;
    }

    if (*matchlength < min_matchlength) {
      *finalptr = *initptr;
    } else {
      *finalptr = (*initptr) + 1;
      while (*finalptr < n &&
	     Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						/*univdiagonal*/region_term + SA[*finalptr] - pos5 + querylength,
						querylength,pos5,pos3,plusp,/*genestrand*/0) >= *matchlength) {
	*finalptr += 1;
      }
    }
  }


#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
  printf("sarray_search_uint1 returning %d..%d with matchlength %d\n",*initptr,*finalptr,*matchlength);
  for (mid = *initptr; mid < *finalptr; mid++) {
    printf("%u\n",region_term + SA[mid]);
  }
  printf("\n");
#endif

#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
  FREE(Buffer);
#endif

  return;
}
#endif



#if defined(DEBUG2) || defined(DEBUG_LOCALDB_DEEP)
static void
oligo_nt (char nt[], int oligo, int oligosize) {
  int i, j;

  j = oligosize-1;
  for (i = 0; i < oligosize; i++) {
    switch (oligo & 0x03) {
    case 0: nt[j] = 'A'; break;
    case 1: nt[j] = 'C'; break;
    case 2: nt[j] = 'G'; break;
    case 3: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  nt[oligosize] = '\0';
  return;
}
#endif


static int
nt_oligo (char *query, int oligosize) {
  int oligo = 0;
  int i;

  for (i = 0; i < oligosize; i++) {
    oligo *= 4;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default: oligo += 3; break;	/* Treat X as T */
    }
  }

  return oligo;
}

static unsigned char
nt_oligo_byte (char *query, int oligosize) {
  unsigned char oligo = 0;
  int i;

  for (i = 0; i < oligosize; i++) {
    oligo <<= 2;
    
    switch (query[i]) {
    case 'A': break;
    case 'C': oligo += 1; break;
    case 'G': oligo += 2; break;
    case 'T': oligo += 3; break;
    default: oligo += 3; break;	/* Treat X as T */
    }
  }

  return oligo;
}

static inline unsigned char
revise_oligo_byte (unsigned char oligo, char c) {
  switch (c) {
  case 'A': return (oligo << 2);
  case 'C': return ((oligo << 2) + 1);
  case 'G': return ((oligo << 2) + 2);
  case 'T': return ((oligo << 2) + 3);
  default: return ((oligo << 2) + 3);	/* Treat X as T */
  }
}


#if 0
static void
sarray_search_uint2 (int *matchlength, int *initptr, int *finalptr, uint16_t *saindex, uint16_t *SA, int n,
		     char *query, int querylength, int pos5, int pos3, int min_matchlength, Univcoord_T region_term,
		     Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int low, high, mid;
  int oligo;
  Univcoord_T pos;
  int nmatches = 0;
  int seglength = pos3 - pos5;
  char mismatch_char;

#if defined(DEBUG1) || defined(DEBUG9) || defined(DEBUG_LOCALDB_DEEP)
  char *suffix;
  suffix = (char *) MALLOC((seglength+2)*sizeof(char)); /* Add 2 so we can print $ at position n */
#endif

  debug1(printf("Entered sarray_search_uint2\n"));

  
#if defined(DEBUG9) || defined(DEBUG_LOCALDB_DEEP)
  int i;
  unsigned int k = 0;

#ifdef DEBUG9A
  /* Dump saindex pointers */
  for (k = 0; k < SAINDEX16_LENGTH; k++) {
    printf("%d\t%d\n",k,saindex[k]);
  }
  printf("\n");
#endif

  for (i = 0; i < n; i++) {
    Genome_fill_buffer_simple(genomecomp,/*left*/region_term+SA[i],seglength,suffix);
    if (SA[i] + seglength >= n) {
      suffix[n - SA[i]] = '$';
      suffix[n - SA[i] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",i,SA[i],suffix);
  }
#endif

#if 0
  for (int i = 0; i < n; i++) {
    if (SA[i] >= n) {
      printf("At i %d, SA[i] is %u, but n is %u\n",i,SA[i],n);
      abort();
    }
  }
#endif

  if (seglength < OLIGOSIZE) {
    low = 0;
    high = n;
  } else if ((oligo = nt_oligo(&(query[pos5]),OLIGOSIZE)) == SAINDEX16_LENGTH - 1) {
    low = saindex[oligo];
    high = n;
  } else {
    low = saindex[oligo];
    high = saindex[oligo+1];
  }
  debug1(printf("query: %.*s => oligo %d, initial low %d, initial high %d\n",
	       OLIGOSIZE,&(query[pos5]),oligo,low,high));

  if (low == high) {
    /* matchlength must be less than OLIGOSIZE */
    *matchlength = 0;
    *initptr = *finalptr = low;
    debug1(printf("sarray_search_uint2 returning %d..%d\n\n",*initptr,*finalptr));
#if defined(DEBUG1) || defined(DEBUG9) || defined(DEBUG_LOCALDB_DEEP)
    FREE(suffix);
#endif
    return;
  }

  while (low < high) {
    /* Compute mid for unsigned ints.  Want floor((low+high)/2). */
    /* low 0, high 2, mid = 1.  low 0, high 3, mid = 1.  low 1, high 2, mid = 1.  low 1, high 3, mid = 2 */
    mid = low + (high - low + 0)/2;

    pos = region_term + SA[mid];

#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
    printf("low %d, high %d, mid %d => position %u\n",low,high,mid,pos);
    Genome_fill_buffer_simple(genomecomp,/*left*/pos,seglength,suffix);
    if (SA[mid] + seglength >= n) {
      suffix[n - SA[mid]] = '$';
      suffix[n - SA[mid] + 1] = '\0';
    }
    printf("%d\t%d\t%s\n",mid,SA[mid],suffix);
    printf("query\t\t%.*s\n",seglength,&(query[pos5]));
#endif

    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);


    if (SA[mid] + nmatches >= n) {
      /* Terminating char is always less than query */
      debug1(printf("nmatches %d, genomic char $\n",nmatches));
      low = mid + 1;
    } else if (nmatches >= seglength) {
      debug1(printf("nmatches %d >= seglength\n",nmatches));
      high = mid;
    } else if (mismatch_char > query[pos5+nmatches]) {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c > query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      high = mid;
    } else {
      /* assert(mismatch_char == Genome_bits_get_char(pos+nmatches)); */
      debug1(printf("nmatches %d, genomic char %c <= query %c\n",nmatches,
		   mismatch_char,query[pos5+nmatches]));
      low = mid + 1;
    }
  }

  if ((*initptr = low) == n) {
    /* Query is higher than any suffix */
    *matchlength = 0;
    *finalptr = *initptr;
  } else {
    pos = region_term + SA[*initptr];
    nmatches = Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						  /*univdiagonal*/pos-pos5+querylength,querylength,
						  pos5,pos3,plusp,/*genestrand*/0);
    debug1(printf("Setting initptr to be %d with %d matches\n",low,nmatches));


    if (nmatches < seglength) {
      *matchlength = nmatches;
    } else {
      *matchlength = seglength;
    }

    if (*matchlength < min_matchlength) {
      *finalptr = *initptr;
    } else {
      *finalptr = (*initptr) + 1;
      while (*finalptr < n &&
	     Genomebits_consecutive_matches_wmm(&mismatch_char,genomebits,query_compress,
						/*univdiagonal*/(region_term + SA[*finalptr])-pos5+querylength,
						querylength,pos5,pos3,
						plusp,/*genestrand*/0) >= *matchlength) {
	*finalptr += 1;
      }
    }
  }


#if defined(DEBUG1) || defined(DEBUG_LOCALDB_DEEP)
  printf("sarray_search_uint2 returning %d..%d with matchlength %d\n",*initptr,*finalptr,*matchlength);
  for (mid = *initptr; mid < *finalptr; mid++) {
    printf("%u\n",region_term + SA[mid]);
  }
  printf("\n");
#endif
  
#if defined(DEBUG1) || defined(DEBUG9) || defined(DEBUG_LOCALDB_DEEP)
  FREE(suffix);
#endif

  return;
}
#endif


#if 0
static int
univcoord_compare (const void *a, const void *b) {
  Univcoord_T x = * (Univcoord_T *) a;
  Univcoord_T y = * (Univcoord_T *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}
#endif


#if 0
/* Hard to find maximum matchlength without scanning first */
static int
get_general (int *matchlength, Univdiag_T *univdiags, T this,
	     char *queryptr, int pos5, int pos3, int querylength, int min_matchlength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int region_matchlength;
  int nentries;
  Univdiag_T *out = univdiags, univdiag, last_univdiag, *p, *q;
  Univcoord_T low_position_start, high_position_start, diagonal;
  int initptr, finalptr, ptr, n;
  uint16_t *saindex16;
  uint16_t *sarray16;
  uint8_t *sarray8;
  Univcoord_T low_regioni, high_regioni, regioni;
  Univcoord_T large_region_term, small_region_term;

  *matchlength = 0;

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  debug1(printf("Entering get_general\n"));

  /* low region */
  regioni = low_regioni;
  large_region_term = regioni * SARRAY16_LENGTH;
  if (low_position_start < (small_region_term = large_region_term + SARRAY16_LENGTH - SARRAY8_LENGTH/2)) {
    /* Start with L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = SARRAY16_LENGTH;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = univdiags;		/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = large_region_term + sarray16[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + region_matchlength;
	  univdiag->nmismatches = 0;
	}
      }
    }
  }

  /* S */
  if (high_regioni > low_regioni) {
    sarray8 = &(this->sarray8[regioni * SARRAY8_LENGTH]);
    n = SARRAY8_LENGTH;
    debug1(printf("S* (region %u): %u..%u\n",regioni,small_region_term,small_region_term + n));
    sarray_search_uint1(&region_matchlength,&initptr,&finalptr,sarray8,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,small_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = univdiags;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = small_region_term + sarray8[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + region_matchlength;
	  univdiag->nmismatches = 0;
	}
      }
    }
  }
  
  /* Middle regions */
  for (regioni = low_regioni + 1; regioni < high_regioni; regioni++) {
    large_region_term += SARRAY16_LENGTH;
    small_region_term += SARRAY16_LENGTH;

    /* L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = SARRAY16_LENGTH;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = univdiags;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	diagonal = large_region_term + sarray16[ptr] - pos5 + querylength;
	assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	univdiag = *out++;
	univdiag->univdiagonal = diagonal;
	univdiag->qstart = pos5;
	univdiag->qend = pos5 + region_matchlength;
	univdiag->nmismatches = 0;
      }
    }
    
    /* S */
    sarray8 = &(this->sarray8[regioni * SARRAY8_LENGTH]);
    n = (small_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small_region_term;
    debug1(printf("S* (region %u): %u..%u\n",regioni,small_region_term,small_region_term + n));
    sarray_search_uint1(&region_matchlength,&initptr,&finalptr,sarray8,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,small_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = univdiags;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = small_region_term + sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + region_matchlength;
	  univdiag->nmismatches = 0;
	}
      }
    }
  }

  /* High region */
  large_region_term += SARRAY16_LENGTH;
  if (high_univdiagonal < large_region_term + SARRAY8_LENGTH/2) {
    /* No need for L */
  } else {      
    /* L */
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sarray16 = &(this->sarray16[large_region_term]);
    n = (large_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - large_region_term;
    debug1(printf("L* (region %u): %u..%u\n",regioni,large_region_term,large_region_term + n));
    sarray_search_uint2(&region_matchlength,&initptr,&finalptr,saindex16,sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,large_region_term,
			query_compress,plusp,genomebits);
    if (region_matchlength < *matchlength) {
      /* Skip */
    } else {
      if (region_matchlength > *matchlength) {
	*matchlength = region_matchlength;
	out = univdiags;	/* Reset output */
      }
      for (ptr = initptr; ptr < finalptr; ptr++) {
	if ((diagonal = large_region_term + sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + region_matchlength;
	  univdiag->nmismatches = 0;
	}
      }
    }
  }

  if ((nentries = out - univdiags) == 0) {
    return 0;

  } else {
    qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);

    /* Remove duplicates */
    q = univdiags;
    last_univdiag = *q;
    for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
      if (Univdiag_equal(*p,last_univdiag) == true) {
	/* Skip: repetitive */
      } else {
	Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
      }
      last_univdiag = *p;
    }
    nentries = (q - univdiags) + 1;
    return nentries;
  }
}
#endif


#if 0
static void
print_univdiagonals (Univcoord_T *diagonals, int nentries) {
  int i;

  printf("%d entries starting at address %p\n",nentries,diagonals);
  for (i = 0; i < nentries; i++) {
    printf("%u\n",diagonals[i]);
  }
  printf("\n");
    
  return;
}
#endif

#if defined(DEBUG0) || defined(DEBUG2) || defined(DEBUG_LOCALDB_DEEP)
static void
print_univdiags (Univdiag_T *univdiags, int nentries) {
  int i;

  printf("%d entries starting at address %p\n",nentries,univdiags);
  for (i = 0; i < nentries; i++) {
    printf("%u %d..%d\n",univdiags[i]->univdiagonal,univdiags[i]->qstart,univdiags[i]->qend);
  }
  printf("\n");
    
  return;
}
#endif


#ifdef DEBUG15
static void
check_univdiagonals_general (Univcoord_T *diagonals, int nentries, T this,
			     char *queryptr, int pos5, int pos3, int querylength,
			     Univcoord_T low_position, Univcoord_T high_position,
			     Compress_T query_compress, Genomebits_T genomebits) {
  Univcoord_T *diagonals_alt;
  int nentries_alt, i;
  bool errorp = false;


  diagonals_alt = get_general(&nentries_alt,this,queryptr,
					    pos5,pos3,querylength,
					    low_position,high_position,
					    query_compress,plusp,genomebits);

  if (nentries != nentries_alt) {
    errorp = true;
  } else {
    for (i = 0; i < nentries; i++) {
      if (diagonals[i] != diagonals_alt[i]) {
	errorp = true;
      }
    }
  }

  if (errorp == true) {
    printf("Standard method: %d entries\n",nentries);
    for (i = 0; i < nentries; i++) {
      printf("%u\n",diagonals[i]);
    }
    printf("\n");
    
    printf("General method: %d entries\n",nentries_alt);
    for (i = 0; i < nentries_alt; i++) {
      printf("%u\n",diagonals_alt[i]);
    }
    printf("\n");
    abort();
  }

  FREE(diagonals_alt);
  return;
}
#endif


#if 0
/* ? Do results need to be sorted */
/* This is time-consuming compared with get_by_ends and
   get_exhaustive, because of the calls to
   Genomebits_consecutive_matches_wmm */
static int
get_by_sarray (int *matchlength, Univdiag_T *univdiags, T this,
	       char *queryptr, int pos5, int pos3, int querylength, int min_matchlength,
	       Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	       Compress_T query_compress, bool plusp, Genomebits_T genomebits) {
  int nentries;
  Univdiag_T *out = univdiags, univdiag, last_univdiag, *p, *q;
  Univcoord_T low_position_start, high_position_start, diagonal;
  int low_initptr, low_finalptr, mid_initptr, mid_finalptr, high_initptr, high_finalptr,
    small0_initptr, small0_finalptr, small1_initptr, small1_finalptr, ptr, n;
  uint16_t *low_saindex16, *mid_saindex16, *high_saindex16;
  uint16_t *low_sarray16, *mid_sarray16, *high_sarray16;
  uint8_t *small0_sarray8, *small1_sarray8;
  Univcoord_T low_regioni, mid_regioni, high_regioni;
  Univcoord_T low_region_term, mid_region_term, high_region_term,
    small0_region_term, small1_region_term;
  int low_matchlength, mid_matchlength, high_matchlength, small0_matchlength, small1_matchlength;

  debug1(printf("Entered get_by_sarray with pos5 %d, pos3 %d, %u..%u\n",
	       pos5,pos3,low_univdiagonal,high_univdiagonal));

  assert(pos5 < pos3);

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  debug1(printf("low %u => region %u.  high %u => region %u\n",
	       low_position_start,low_regioni,high_position_start,high_regioni));

  if (high_regioni == low_regioni) {
    debug1(printf("Case L\n"));
    low_region_term = low_regioni * SARRAY16_LENGTH;
    low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
    low_sarray16 = &(this->sarray16[low_region_term]);
    n = (low_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH : genomelength - low_region_term;
    debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
    sarray_search_uint2(&(*matchlength),&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			n,queryptr,querylength,pos5,pos3,
			min_matchlength,low_region_term,
			query_compress,plusp,genomebits);

    for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
      if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal &&
	  diagonal <= high_univdiagonal) {
	univdiag = *out++;
	univdiag->univdiagonal = diagonal;
	univdiag->qstart = pos5;
	univdiag->qend = pos5 + (*matchlength);
	univdiag->nmismatches = 0;
      }
    }

    if ((nentries = out - univdiags) == 0) {
      debug0(print_univdiags(univdiags,nentries));
      return 0;

    } else {
      /* One region, so no need to remove duplicates */
#ifdef ALWAYS_SORT
      qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
#endif
      debug0(print_univdiags(univdiags,nentries));
      return nentries;
    }

  } else if (high_regioni == low_regioni + 1) {
    low_region_term = low_regioni * SARRAY16_LENGTH;
    high_region_term = low_region_term + SARRAY16_LENGTH;

    if (low_position_start >= high_region_term - SARRAY8_LENGTH/2 &&
	high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case S\n"));
      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = (small0_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small0_region_term;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&(*matchlength),&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal &&
	    diagonal <= high_univdiagonal) {
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;
      } else {
	/* One region, so no need to remove duplicates */
#ifdef ALWAYS_SORT
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
#endif
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else if (low_position_start >= high_region_term - SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L\n"));
      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH : genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else if (high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case L,S\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = (small0_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small0_region_term;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else {
      debug1(printf("Case L,S,L\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = high_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }
    }

  } else if (high_regioni == low_regioni + 2) {
    mid_regioni = low_regioni + 1;
    low_region_term = low_regioni * SARRAY16_LENGTH;
    mid_region_term = low_region_term + SARRAY16_LENGTH;
    high_region_term = mid_region_term + SARRAY16_LENGTH;

    if (low_position_start >= mid_region_term - SARRAY8_LENGTH/2 &&
	high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L,S\n"));
      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);

      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = (small1_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small1_region_term;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  if ((diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else if (low_position_start >= mid_region_term - SARRAY8_LENGTH/2) {
      debug1(printf("Case S,L,S,L\n"));
      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);

      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);


      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  if ((diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else if (high_position_start < high_region_term + SARRAY8_LENGTH/2) {
      debug1(printf("Case L,S,L,S\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);
      
      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = (small1_region_term + SARRAY8_LENGTH < genomelength) ? SARRAY8_LENGTH : genomelength - small1_region_term;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      
      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  if ((diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return nentries;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }

    } else {
      debug1(printf("Case L,S,L,S,L\n"));
      low_saindex16 = &(this->saindex16[low_regioni * SAINDEX16_LENGTH]);
      low_sarray16 = &(this->sarray16[low_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",low_regioni,low_region_term,low_region_term + n));
      sarray_search_uint2(&low_matchlength,&low_initptr,&low_finalptr,low_saindex16,low_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,low_region_term,
			  query_compress,plusp,genomebits);

      small0_region_term = mid_region_term - SARRAY8_LENGTH/2;
      small0_sarray8 = &(this->sarray8[low_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",low_regioni,small0_region_term,small0_region_term + n));
      sarray_search_uint1(&small0_matchlength,&small0_initptr,&small0_finalptr,small0_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small0_region_term,
			  query_compress,plusp,genomebits);

      mid_saindex16 = &(this->saindex16[mid_regioni * SAINDEX16_LENGTH]);
      mid_sarray16 = &(this->sarray16[mid_region_term]);
      n = SARRAY16_LENGTH;
      debug1(printf("L (region %u): %u..%u\n",mid_regioni,mid_region_term,mid_region_term + n));
      sarray_search_uint2(&mid_matchlength,&mid_initptr,&mid_finalptr,mid_saindex16,mid_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,mid_region_term,
			  query_compress,plusp,genomebits);
      
      small1_region_term = high_region_term - SARRAY8_LENGTH/2;
      small1_sarray8 = &(this->sarray8[mid_regioni * SARRAY8_LENGTH]);
      n = SARRAY8_LENGTH;
      debug1(printf("S (region %u): %u..%u\n",mid_regioni,small1_region_term,small1_region_term + n));
      sarray_search_uint1(&small1_matchlength,&small1_initptr,&small1_finalptr,small1_sarray8,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,small1_region_term,
			  query_compress,plusp,genomebits);

      high_saindex16 = &(this->saindex16[high_regioni * SAINDEX16_LENGTH]);
      high_sarray16 = &(this->sarray16[high_region_term]);
      n = (high_region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - high_region_term;
      debug1(printf("L (region %u): %u..%u\n",high_regioni,high_region_term,high_region_term + n));
      sarray_search_uint2(&high_matchlength,&high_initptr,&high_finalptr,high_saindex16,high_sarray16,
			  n,queryptr,querylength,pos5,pos3,
			  min_matchlength,high_region_term,
			  query_compress,plusp,genomebits);

      *matchlength = low_matchlength;
      if (small0_matchlength > *matchlength) *matchlength = small0_matchlength;
      if (mid_matchlength > *matchlength) *matchlength = mid_matchlength;
      if (small1_matchlength > *matchlength) *matchlength = small1_matchlength;
      if (high_matchlength > *matchlength) *matchlength = high_matchlength;

      if (low_matchlength == *matchlength) {
	for (ptr = low_initptr; ptr < low_finalptr; ptr++) {
	  if ((diagonal = low_region_term + low_sarray16[ptr] - pos5 + querylength) >= low_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }
      if (small0_matchlength == *matchlength) {
	for (ptr = small0_initptr; ptr < small0_finalptr; ptr++) {
	  diagonal = small0_region_term + small0_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (mid_matchlength == *matchlength) {
	for (ptr = mid_initptr; ptr < mid_finalptr; ptr++) {
	  diagonal = mid_region_term + mid_sarray16[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (small1_matchlength == *matchlength) {
	for (ptr = small1_initptr; ptr < small1_finalptr; ptr++) {
	  diagonal = small1_region_term + small1_sarray8[ptr] - pos5 + querylength;
	  assert(diagonal >= low_univdiagonal && diagonal <= high_univdiagonal);
	  univdiag = *out++;
	  univdiag->univdiagonal = diagonal;
	  univdiag->qstart = pos5;
	  univdiag->qend = pos5 + (*matchlength);
	  univdiag->nmismatches = 0;
	}
      }
      if (high_matchlength == *matchlength) {
	for (ptr = high_initptr; ptr < high_finalptr; ptr++) {
	  if ((diagonal = high_region_term + high_sarray16[ptr] - pos5 + querylength) <= high_univdiagonal) {
	    univdiag = *out++;
	    univdiag->univdiagonal = diagonal;
	    univdiag->qstart = pos5;
	    univdiag->qend = pos5 + (*matchlength);
	    univdiag->nmismatches = 0;
	  }
	}
      }

      if ((nentries = out - univdiags) == 0) {
	debug0(print_univdiagonals(/*diagonals*/NULL,nentries));
	return 0;

      } else {
	qsort(univdiags[0],nentries,sizeof(struct Univdiag_T),Univdiag_struct_cmp);
	/* Remove duplicates */
	q = univdiags;
	last_univdiag = *q;
	for (p = &(univdiags[1]); p < &(univdiags[nentries]); p++) {
	  if (Univdiag_equal(*p,last_univdiag) == true) {
	    /* Skip: repetitive */
	  } else {
	    Univdiag_transfer(*++q,*p); /* Not pointers: *++q = *p */
	  }
	  last_univdiag = *p;
	}
	nentries = (q - univdiags) + 1;
	debug0(print_univdiags(univdiags,nentries));
	return nentries;
      }
    }

  } else {
#if 0
    /* If the region is larger than 65536, then it is likely coming from end < start */
    fprintf(stderr,"get_by_sarray cannot handle a region this large (%llu) from %llu to %llu => regioni %d..%d\n",
	    (unsigned long long) (high_position_start - low_position_start),
	    (unsigned long long) low_position_start,(unsigned long long) high_position_start,
	    low_regioni,high_regioni);
    abort();
    return 0;
#else
    return get_general(&(*matchlength),&(*univdiags),this,
		       queryptr,pos5,pos3,querylength,min_matchlength,
		       low_univdiagonal,high_univdiagonal,query_compress,
		       plusp,genomebits);
#endif
  }
}
#endif



#if 0
/* Overwrites existing array */
static int
remove_duplicate_diagonals (Univcoord_T *diagonals, int ndiagonals) {
  Univcoord_T *out = diagonals;
  int i, j;
  
  i = 0;
  while (i < ndiagonals) {
    *out++ = diagonals[i];
    j = i + 1;
    while (j < ndiagonals && diagonals[j] == diagonals[i]) {
      j++;
    }
    i = j;
  }
  
  return (out - diagonals);
}
#endif


/* This method can lead to incorrect results if a mismatch occurs near
   an end, so allow it only with zero mismatches, or require that the
   mismatch is in the middle */
static int
get_by_ends (Univdiag_T *univdiags,
	     T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
	     char *queryptr, int pos5, int pos3, int querylength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, int genestrand,
	     Genomebits_T genomebits, int nmismatches_allowed) {
  int nentries;
  Univdiag_T *out, univdiag;
  Univcoord_T *ptr = diagonals_alloc;
  int lowi5, highi5, lowi3, highi3, initptr, finalptr, nptrs5, nptrs3, n;
  uint16_t *saindex16;
  uint16_t *sasort16;
  Univcoord_T low_position_start, high_position_start, low_regioni, high_regioni, regioni;
  Univcoord_T region_term;

  int i;
  unsigned char oligo5, oligo3;
  int ref_nmismatches, best_nmismatches, mm;
#ifndef HAVE_SSE2
  unsigned short *ignore_ptr;
#endif
  
#if defined(DEBUG2) || defined(DEBUG_LOCALDB_DEEP)
  char nt[5];
#endif


  debug2(printf("Entering get_by_ends with %.*s\n",pos3-pos5,&(queryptr[pos5])));  
  debug2(printf("low_univdiagonal %u, high_univdiagonal %u\n",low_univdiagonal,high_univdiagonal));
  assert(pos3 - pos5 >= OLIGOSIZE);

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;

  oligo5 = nt_oligo_byte(&(queryptr[pos5]),OLIGOSIZE);
  oligo3 = nt_oligo_byte(&(queryptr[pos3 - OLIGOSIZE]),OLIGOSIZE);

  /* Gather 4-mer positions for the two end oligos */
  region_term = low_regioni * SARRAY16_LENGTH;
  for (regioni = low_regioni; regioni <= high_regioni; regioni++) {
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sasort16 = &(this->sasort16[region_term]);
    n = (region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - region_term;

    initptr = saindex16[oligo5];
    finalptr = (oligo5 == SAINDEX16_LENGTH - 1) ? n : saindex16[oligo5+1];
    debug2(oligo_nt(nt,oligo5,4));
    debug2(printf("oligo5 %d (%s) => initptr %d, finalptr %d\n",oligo5,nt,initptr,finalptr));

    /* Filter diagonals */
    /* Bound by low_position and high_position */
    lowi5 = initptr;
    while (lowi5 < finalptr && (region_term + sasort16[lowi5]) - pos5 + querylength < low_univdiagonal) {
      lowi5++;
    }
    highi5 = finalptr - 1;
    while (highi5 >= lowi5 && (region_term + sasort16[highi5]) - pos5 + querylength > high_univdiagonal) {
      highi5--;
    }
    nptrs5 = highi5 - lowi5 + 1;


    initptr = saindex16[oligo3];
    finalptr = (oligo3 == SAINDEX16_LENGTH - 1) ? n : saindex16[oligo3+1];
    debug2(oligo_nt(nt,oligo3,4));
    debug2(printf("oligo3 %d (%s) => initptr %d, finalptr %d\n",oligo3,nt,initptr,finalptr));

    /* Filter diagonals */
    /* Bound by low_position and high_position */
    lowi3 = initptr;
    while (lowi3 < finalptr && (region_term + sasort16[lowi3]) - pos5 + querylength < low_univdiagonal) {
      lowi3++;
    }
    highi3 = finalptr - 1;
    while (highi3 >= lowi3 && (region_term + sasort16[highi3]) - pos5 + querylength > high_univdiagonal) {
      highi3--;
    }
    nptrs3 = highi3 - lowi3 + 1;

    if (nptrs5 > 0 && nptrs3 > 0) {
#ifdef HAVE_SSE2
      ptr += Intersect_uint2(ptr,localdb_alloc,
			     &(sasort16[lowi5]),nptrs5,/*diagterm1*/querylength - pos5,
			     &(sasort16[lowi3]),nptrs3,/*diagterm2*/querylength - (pos3 - OLIGOSIZE),
			     region_term);
#else
      ptr += Intersect_uint2_scalar(ptr,&ignore_ptr,
				    &(sasort16[lowi5]),nptrs5,/*diagterm1*/querylength - pos5,
				    &(sasort16[lowi3]),nptrs3,/*diagterm2*/querylength - (pos3 - OLIGOSIZE),
				    region_term);
#endif
    }

    region_term += SARRAY16_LENGTH;
  }

  if ((nentries = ptr - diagonals_alloc) == 0) {
    return 0;

  } else {
    best_nmismatches = querylength;
    out = univdiags;
    for (i = 0; i < nentries; i++) {
      if ((mm = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,/*alt*/NULL,query_compress,
						      /*univdiagonal*/diagonals_alloc[i],querylength,
						      pos5,pos3,plusp,genestrand)) < best_nmismatches) {
	/* Reset */
	out = univdiags;
	best_nmismatches = mm;
	
	univdiag = *out++;
	univdiag->univdiagonal = diagonals_alloc[i];
	univdiag->qstart = pos5;
	univdiag->qend = pos3;
	univdiag->nmismatches = best_nmismatches;

      } else if (mm == best_nmismatches) {
	univdiag = *out++;
	univdiag->univdiagonal = diagonals_alloc[i];
	univdiag->qstart = pos5;
	univdiag->qend = pos3;
	univdiag->nmismatches = best_nmismatches;
      }
    }

    if (best_nmismatches > nmismatches_allowed) {
      return 0;

    } else {
      nentries = out - univdiags;
      debug2(print_univdiags(univdiags,nentries));
      /* *nmismatches = best_nmismatches; */
      return nentries;
    }
  }
}


#define SUBOPT 2


/* Taken from segmentring.c */
/* Unlike other binary search algorithms in the code, these expect
   highi to be (length - 1).  They are a little faster because they
   remove one if condition */

/* Returns first value greater than or equal to the goal */
static int
binary_search_ge (int lowi, int highi, const uint16_t *positions, uint16_t goal) {
  int answer = highi + 1, middlei;

  debug20(printf("entered binary_search_ge with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi <= highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (positions[middlei] >= goal) {
      answer = middlei;
      highi = middlei - 1;
    } else {
      lowi = middlei + 1;
    }
  }

  debug20(printf("binary_search_ge returns %d (lowi %d, middlei %d, highi %d)\n",
		 answer,lowi,middlei,highi));
  return answer;
}

/* Returns first value less than or equal to the goal */
static int
binary_search_le (int lowi, int highi, const uint16_t *positions, uint16_t goal) {
  int answer = lowi - 1, middlei;

  debug20(printf("entered binary_search_le with lowi=%d, highi=%d, goal=%u\n",lowi,highi,goal));

  while (lowi <= highi) {
    middlei = lowi + ((highi - lowi) / 2);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,positions[lowi],middlei,positions[middlei],
		   highi,positions[highi],goal));
    if (positions[middlei] <= goal) {
      answer = middlei;
      lowi = middlei + 1;
    } else {
      highi = middlei - 1;
    }
  }

  debug20(printf("binary_search_le returns %d (lowi %d, middlei %d, highi %d)\n",
		 answer,lowi,middlei,highi));
  return answer;
}



/* Gets positions of 4-mers from low_position to high_position.  We do
   not worry about the small regions */

/* If trim5p or trim3p is true, this indicates where distal trimming
   can be performed.  Otherwise, medial trimming is not done, since
   we want to know what that solution yields. */
static int
get_exhaustive (bool *trimmedp, int *best_matchlength, Univdiag_T *univdiags,
		T this, int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		char *queryptr, int pos5, int pos3, int querylength,
		Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		struct Accel_profile_localdb_work_T *profile_work) {
  int nentries;
  Univdiag_T *q, univdiag;
  Univcoord_T low_position_start, high_position_start, univdiagonal;
  int rel_univdiagonal, ii, i;

  int goal;
  int initptr, finalptr, ptr, lowi, highi, n;
  uint16_t *saindex16;
  uint16_t *sasort16;
  Univcoord_T low_regioni, high_regioni, regioni;
  Univcoord_T region_term;

  bool overflow_exists_p = false;
  int16_t max_count = 0, subopt_count, count; /* signed so we can use SIMD _cmpgt_epi16 */
  int *qstarts, *qends;
  int nblocks, npositions;

  int pos;
  unsigned char oligo;
  
  int matchlength;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long profile_t0 = 0ULL;
#endif

#if !defined(HAVE_SSE2)
#elif !defined(HAVE_AVX2)
  __m128i _tally, _overflowp, _cmp, _threshold;
#elif !defined(HAVE_AVX512)
  __m256i _tally, _overflowp, _cmp, _threshold;
#else
  __m512i _tally, _overflowp, _cmp, _threshold;
  __mmask32 result;
#endif


#if defined(DEBUG2) || defined(DEBUG_LOCALDB_DEEP)
  char nt[OLIGOSIZE + 1];
#endif


  *trimmedp = false;

  /* This would allow medial trimming, but we want to know what the medial alignment yields */
  /* trim5p = trim3p = true; */

  debug2(printf("Entering get_exhaustive with %.*s\n",pos3-pos5,&(queryptr[pos5])));

  if (low_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    low_position_start = 0U;
  } else {
    low_position_start = low_univdiagonal - querylength + pos5;
  }
  if (high_univdiagonal < (Univcoord_T) (querylength - pos5)) {
    high_position_start = 0U;
  } else {
    high_position_start = high_univdiagonal - querylength + pos5;
  }

  npositions = high_univdiagonal - low_univdiagonal + 1;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->npositions = npositions;

  profile_t0 = Accel_profile_now_us();
#endif
#if !defined(HAVE_SSE2)
  /* localdb_tally = (int16_t *) CALLOC((high_univdiagonal - low_univdiagonal + 1),sizeof(int16_t)); */
  /* localdb_overflowp = (int16_t *) CALLOC((high_univdiagonal - low_univdiagonal + 1),sizeof(int16_t)); */
  nblocks = npositions;
  assert(nblocks <= max_nblocks);

  memset((void *) localdb_tally,0,nblocks * sizeof(int16_t));
  memset((void *) localdb_overflowp,0,nblocks * sizeof(int16_t));

#elif !defined(HAVE_AVX2)
  nblocks = (npositions + 7)/8;
  assert(nblocks <= max_nblocks);

  /* localdb_tally = (int16_t *) MALLOC_ALIGN(nblocks * 8 * sizeof(int16_t)); */
  /* localdb_overflowp = (int16_t *) MALLOC_ALIGN(nblocks * 8 * sizeof(int16_t)); */
  memset((void *) localdb_tally,0,nblocks * 8 * sizeof(int16_t));
  memset((void *) localdb_overflowp,0,nblocks * 8 * sizeof(int16_t));

#elif !defined(HAVE_AVX512)
  nblocks = (npositions + 15)/16;
  assert(nblocks <= max_nblocks);

  /* localdb_tally = (int16_t *) MALLOC_ALIGN(nblocks * 16 * sizeof(int16_t)); */
  /* localdb_overflowp = (int16_t *) MALLOC_ALIGN(nblocks * 16 * sizeof(int16_t)); */
  memset((void *) localdb_tally,0,nblocks * 16 * sizeof(int16_t));
  memset((void *) localdb_overflowp,0,nblocks * 16 * sizeof(int16_t));

#else
  nblocks = (npositions + 31)/32;
  if (nblocks > max_nblocks) {
    fprintf(stderr,"npositions %d, nblocks %d, max_nblocks %d\n",npositions,nblocks,max_nblocks);
    exit(9);
  }
  assert(nblocks <= max_nblocks);

  /* localdb_tally = (int16_t *) MALLOC_ALIGN(nblocks * 32 * sizeof(int16_t)); */
  /* localdb_overflowp = (int16_t *) MALLOC_ALIGN(nblocks * 32 * sizeof(int16_t)); */
  memset((void *) localdb_tally,0,nblocks * 32 * sizeof(int16_t));
  memset((void *) localdb_overflowp,0,nblocks * 32 * sizeof(int16_t));
#endif
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->nblocks = nblocks;
  profile_work->memset_us += Accel_profile_now_us() - profile_t0;

  profile_t0 = Accel_profile_now_us();
#endif
  qstarts = (int *) MALLOC(npositions * sizeof(int));
  qends = (int *) MALLOC(npositions * sizeof(int));
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->qalloc_us += Accel_profile_now_us() - profile_t0;
#endif


  low_regioni = low_position_start/SARRAY16_LENGTH;
  high_regioni = high_position_start/SARRAY16_LENGTH;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->region_count = (int) (high_regioni - low_regioni + 1);
#endif

  region_term = low_regioni * SARRAY16_LENGTH;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_t0 = Accel_profile_now_us();
#endif
  for (regioni = low_regioni; regioni <= high_regioni; regioni++) {
    saindex16 = &(this->saindex16[regioni * SAINDEX16_LENGTH]);
    sasort16 = &(this->sasort16[region_term]);
    n = (region_term + SARRAY16_LENGTH < genomelength) ? SARRAY16_LENGTH: genomelength - region_term;

    if (pos5 + OLIGOSIZE <= pos3) {
      oligo = nt_oligo_byte(&(queryptr[pos5]),OLIGOSIZE - 1);
    }
    for (pos = pos5; pos + OLIGOSIZE <= pos3; pos++) {
#ifdef GSNAP_ACCEL_PROFILE_STATS
      profile_work->oligo_probe_count++;
      profile_work->sa_range_count++;
#endif
      /* printf("Revising oligo with %c\n",queryptr[pos + OLIGOSIZE - 1]); */
      if ((oligo = revise_oligo_byte(oligo,queryptr[pos + OLIGOSIZE - 1])) == SAINDEX16_LENGTH - 1) {
	initptr = saindex16[oligo];
	finalptr = n;
      } else {
	initptr = saindex16[oligo];
	finalptr = saindex16[oligo+1];
      }

      /* printf("regioni %d compared with low_regioni %d and high_regioni %d\n",
	 regioni,low_regioni,high_regioni); */
      if (regioni > low_regioni + 1) { /* Check over two regions, in case we are near a boundary */
	lowi = initptr;
      } else if ((goal = (low_univdiagonal - querylength + pos) - region_term) < 0) {
	/* All of the array elements pass */
	lowi = initptr;
      } else if (goal > 65535) {
	/* None of the array elements pass */
	lowi = finalptr;
      } else {
	/* Filter diagonals */
	/* printf("Goal is low_univdiagonal %u - querylength %d + pos %d - region_term %u = %d\n",
	   low_univdiagonal,querylength,pos,region_term,(low_univdiagonal - querylength + pos) - region_term); */
	lowi = binary_search_ge(/*lowi*/initptr,/*highi*/finalptr - 1,sasort16,(uint16_t) goal);
      }

      if (high_regioni > 1 && regioni < high_regioni - 1) {
	highi = finalptr - 1;
      } else if ((goal = (high_univdiagonal - querylength + pos) - region_term) > 65535) {
	/* All of the array elements pass */
	highi = finalptr - 1;
      } else if (goal < 0) {
	/* None of the array elements pass */
	highi = initptr - 1;
      } else {
	/* Filter diagonals */
	/* printf("Goal is high_univdiagonal %u - querylength %d + pos %d - region_term %u = %d\n",
	   high_univdiagonal,querylength,pos,region_term,(high_univdiagonal - querylength + pos) - region_term); */
	highi = binary_search_le(lowi,/*highi*/finalptr - 1,sasort16,(uint16_t) goal);
      }

      debug2(oligo_nt(nt,oligo,4));
      debug2(printf("oligo %d (%s) => initptr %d, finalptr %d => lowi %d, highi %d\n",
		    oligo,nt,initptr,finalptr,lowi,highi));
      /* lowi and highi are inclusive */
#ifdef GSNAP_ACCEL_PROFILE_STATS
      if (highi >= lowi) {
	profile_work->tally_increment_count += (unsigned long long) (highi - lowi + 1);
      }
#endif
      
      for (ptr = lowi; ptr <= highi; ptr++) {
	univdiagonal = region_term + sasort16[ptr] - pos + querylength; /* diagonal */
	/* printf("univdiagonal %u = region_term %u + %d - pos %d + querylength %d, low_univdiagonal %u, high_univdiagonal %u\n", 
	   univdiagonal,region_term,sasort16[ptr],pos,querylength,low_univdiagonal,high_univdiagonal); */
	assert(univdiagonal >= low_univdiagonal && univdiagonal <= high_univdiagonal);

	rel_univdiagonal = (int) (univdiagonal - low_univdiagonal);
	if ((count = ++localdb_tally[rel_univdiagonal]) < 0) {
	  /* printf("Got overflow of %d at rel position %d\n",count,rel_univdiagonal); */
	  localdb_tally[rel_univdiagonal] = 0;
	  localdb_overflowp[rel_univdiagonal] = 1;
	  overflow_exists_p = true;
#ifdef GSNAP_ACCEL_PROFILE_STATS
	  profile_work->overflow_count++;
#endif

	} else if (localdb_overflowp[rel_univdiagonal] == true) {
	  /* Don't change qstart or max_count */

	} else if (count == 1) {
	  qstarts[rel_univdiagonal] = pos;

	} else if (count > max_count) {
	  max_count = count;
	}
	qends[rel_univdiagonal] = pos + OLIGOSIZE;
      }
    }

    region_term += SARRAY16_LENGTH;
  }
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->vote_us += Accel_profile_now_us() - profile_t0;
#endif
  debug2(printf("max_count is %d\n",max_count));

  if ((subopt_count = max_count - SUBOPT) < 1) {
    subopt_count = 1;
  }
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->max_count = max_count;
  profile_work->subopt_count = subopt_count;

  profile_t0 = Accel_profile_now_us();
#endif
  *best_matchlength = 0;
  q = univdiags;

#if !defined(HAVE_SSE2)
  if (overflow_exists_p == true) {
    for (i = 0; i < npositions; i++) {
      if (localdb_overflowp[i] > 0) {
	if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	  *best_matchlength = matchlength;
	}
	
	univdiag = *q++;
	univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	univdiag->qstart = qstarts[i];
	univdiag->qend = qends[i];
	univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
      }
    }

  } else {
    for (i = 0; i < npositions; i++) {
      if (localdb_tally[i] >= subopt_count) {
	if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	  *best_matchlength = matchlength;
	}
	
	univdiag = *q++;
	univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	univdiag->qstart = qstarts[i];
	univdiag->qend = qends[i];
	univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
      }
    }
  }

  /* FREE(localdb_tally); */
  /* FREE(localdb_overflowp); */

#elif !defined(HAVE_AVX2)
  /* Using SSE2 commands */
  /* _epi16_offset = _mm_set1_epi16(EPI16_MAX); -- Could subtract from both terms if we want epu16 arithmetic */

  if (overflow_exists_p == true) {
    for (ii = 0; ii < nblocks * 8; ii += 8) {
      _overflowp = _mm_load_si128((const __m128i *) &(localdb_overflowp[ii]));

#ifdef SSE4_1
      if (_mm_testz_si128(_overflowp, _overflowp)) {
#else
      if (!_mm_movemask_epi8(_overflowp)) {
#endif
      } else {
	for (i = ii; i < ii + 8; i++) {
	  if (localdb_overflowp[i] > 0) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }

  } else {
    _threshold = _mm_set1_epi16(subopt_count - 1); /* Subtract 1 because only cmpgt is available not cmpge */

    for (ii = 0; ii < nblocks * 8; ii += 8) {
      _tally = _mm_load_si128((const __m128i *) &(localdb_tally[ii]));
      _cmp = _mm_cmpgt_epi16(_tally, _threshold);

#ifdef SSE4_1
      if (_mm_testz_si128(_cmp, _cmp)) {
#else
      if (!_mm_movemask_epi8(_cmp)) {
#endif
      } else {
	for (i = ii; i < ii + 8; i++) {
	  if (localdb_tally[i] >= subopt_count) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }
  }

  /* FREE_ALIGN(localdb_tally); */
  /* FREE_ALIGN(localdb_overflowp); */

#elif !defined(HAVE_AVX512)
  /* Using AVX2 commands */
  /* _epi16_offset = _mm256_set1_epi16(EPI16_MAX); -- Could subtract from both terms if we want epu16 arithmetic */

  if (overflow_exists_p == true) {
    for (ii = 0; ii < nblocks * 16; ii += 16) {
      _overflowp = _mm256_load_si256((const __m256i *) &(localdb_overflowp[ii]));

      if (_mm256_testz_si256(_overflowp, _overflowp)) {
	/* Skip */
      } else {
	for (i = ii; i < ii + 16; i++) {
	  if (localdb_tally[i] >= subopt_count) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }

  } else {
    _threshold = _mm256_set1_epi16(subopt_count - 1); /* Subtract 1 because only cmpgt is available not cmpge */

    for (ii = 0; ii < nblocks * 16; ii += 16) {
      _tally = _mm256_load_si256((const __m256i *) &(localdb_tally[ii]));
      _cmp = _mm256_cmpgt_epi16(_tally, _threshold);

      if (_mm256_testz_si256(_cmp, _cmp)) {
	/* Skip */
      } else {
	for (i = ii; i < ii + 16; i++) {
	  if (localdb_tally[i] >= subopt_count) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }
  }

  /* FREE_ALIGN(localdb_tally); */
  /* FREE_ALIGN(localdb_overflowp); */

#else
  /* Using AVX512 commands */
  /* _epi16_offset = _mm512_set1_epi16(EPI16_MAX); -- Could subtract from both terms if we want epu16 arithmetic */

  if (overflow_exists_p == true) {
    for (ii = 0; ii < nblocks * 32; ii += 32) {
      _overflowp = _mm512_load_si512((const __m512i *) &(localdb_overflowp[ii]));

      if ((uint32_t) _mm512_test_epi16_mask(_overflowp, _overflowp) == 0) {
	/* Skip */
      } else {
	for (i = ii; i < ii + 32; i++) {
	  if (localdb_tally[i] >= subopt_count) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }

  } else {
    _threshold = _mm512_set1_epi32(subopt_count - 1); /* Subtract 1 because only cmpgt is available not cmpge */

    for (ii = 0; ii < nblocks * 32; ii += 32) {
      _tally = _mm512_load_si512((const __m512i *) &(localdb_tally[ii]));
      if ((uint32_t) _mm512_cmpgt_epi16_mask(_tally, _threshold) == 0) {
	/* Skip */
      } else {
	for (i = ii; i < ii + 32; i++) {
	  if (localdb_tally[i] >= subopt_count) {
	    if ((matchlength = qends[i] - qstarts[i]) > *best_matchlength) {
	      *best_matchlength = matchlength;
	    }

	    univdiag = *q++;
	    univdiag->univdiagonal = low_univdiagonal + (Univcoord_T) i;
	    univdiag->qstart = qstarts[i];
	    univdiag->qend = qends[i];
	    univdiag->nmismatches = matchlength - localdb_tally[i]; /* This is not quite correct */
	  }
	}
      }
    }
  }

  /* FREE_ALIGN(localdb_tally); */
  /* FREE_ALIGN(localdb_overflowp); */
#endif

  FREE(qstarts);
  FREE(qends);

  /* *nmismatches = (*best_matchlength) - best_nmatches; */
  nentries = q - univdiags;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  profile_work->nentries = nentries;
  profile_work->emit_us += Accel_profile_now_us() - profile_t0;
#endif

  debug2(print_univdiags(univdiags,nentries));
  return nentries;
}


#if 0
#ifdef CHECK_ASSERTIONS
static void
check_ascending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  prev_coord = coords[0];
  for (i = 1; i < n; i++) {
    if (coords[i] <= prev_coord) {
      printf("Expecting ascending, but at %d, got %u <= %u\n",
	     i,coords[i],prev_coord);
      abort();
    }
    prev_coord = coords[i];
  }
 
  return;
}
#endif
#endif



/* univdiags_alloc are needed for storing output */
/* diagonals_alloc and localdb_alloc are needed for the results of the first and last sasort entries in get_by_ends */
/* localdb_tally and localdb_overflowp are needed for get_exhaustive */
int
Localdb_get (bool *sortedp, bool *trimmedp, int *matchlength, Univdiag_T *univdiags_alloc, 
	     T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
	     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
	     char *queryptr, int pos5, int pos3, int querylength,
	     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
	     Compress_T query_compress, bool plusp, int genestrand,
	     Genomebits_T genomebits) {
  int nhits;
  unsigned long long cache_hash = 0;
  struct Localdb_cache_entry_T *cache_entry;
#ifdef GSNAP_ACCEL_PROFILE_STATS
  unsigned long long profile_start_us;
#endif
  struct Accel_profile_localdb_work_T profile_work;


  debug(printf("\nEntered Localdb_get with low_univdiagonal %u, high_univdiagonal %u, plusp %d, pos5 %d, pos3 %d, query %s, compress_fwdp %d\n",
	       low_univdiagonal,high_univdiagonal,plusp,pos5,pos3,queryptr,Compress_fwdp(query_compress)));
  assert(Compress_fwdp(query_compress) == plusp);

#ifdef GSNAP_ACCEL_PROFILE_STATS
  Accel_profile_localdb_work_init(&profile_work);
  profile_start_us = Accel_profile_now_us();
#endif
  *matchlength = pos3 - pos5;	/* Needed for get_by_ends */

  if (pos5 >= pos3) {
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_localdb_record("empty",queryptr,pos5,pos3,querylength,
				 low_univdiagonal,high_univdiagonal,plusp,genestrand,
				 /*nhits*/0,*matchlength,/*trimmedp*/false,/*sortedp*/false,
				 Accel_profile_now_us() - profile_start_us,&profile_work);
#endif
    return 0;

  } else if (high_univdiagonal < low_univdiagonal) {
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_localdb_record("empty",queryptr,pos5,pos3,querylength,
				 low_univdiagonal,high_univdiagonal,plusp,genestrand,
				 /*nhits*/0,*matchlength,/*trimmedp*/false,/*sortedp*/false,
				 Accel_profile_now_us() - profile_start_us,&profile_work);
#endif
    return 0;

#if 0
  } else if ((nhits = get_by_sarray(&(*matchlength),&(*univdiags_alloc),this,
				    queryptr,pos5,pos3,querylength,
 				    /*min_matchlength*/(pos3 - pos5),low_univdiagonal,high_univdiagonal,
				    query_compress,plusp,genomebits)) > 0) {
    /* Too time-consuming */
    /* Previously required pos3 - pos5 < index1part + index1interval because otherwise we would have found it from indexdb */
    /* assert(pos3 - pos5 < index1part + index1interval); */

    debug(printf("Localdb_get using get_by_sarray (and no mismatches allowed).  Returning %d hits\n",nhits));
    /* Not guaranteed to be ascending */
    *sortedp = false;
    *trimmedp = false;		/* Entire string matches */
    return nhits;
#endif

  } else if (localdb_memoize_enabled_p == true &&
	     (cache_hash = localdb_cache_hash(this,genomebits,query_compress,queryptr,pos5,pos3,querylength,
					     low_univdiagonal,high_univdiagonal,
					     plusp,genestrand),
	      cache_entry = localdb_cache_lookup(cache_hash,this,genomebits,query_compress,queryptr,
						 pos5,pos3,querylength,
						 low_univdiagonal,high_univdiagonal,
						 plusp,genestrand)) !=
	     (struct Localdb_cache_entry_T *) NULL) {
    localdb_cache_replay(&(*sortedp),&(*trimmedp),&(*matchlength),univdiags_alloc,cache_entry);
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_localdb_record("cache_hit",queryptr,pos5,pos3,querylength,
				 low_univdiagonal,high_univdiagonal,plusp,genestrand,
				 cache_entry->nhits,*matchlength,*trimmedp,*sortedp,
				 Accel_profile_now_us() - profile_start_us,&profile_work);
#endif
    return cache_entry->nhits;

  } else if (pos3 - pos5 > 2*OLIGOSIZE &&
	     (nhits = get_by_ends(&(*univdiags_alloc),
				  this,diagonals_alloc,localdb_alloc,queryptr,pos5,pos3,querylength,
				  low_univdiagonal,high_univdiagonal,query_compress,
				  plusp,genestrand,genomebits,/*nmismatches_allowed*/0)) > 0) {
    /* This method can lead to incorrect results if a mismatch is at
       the ends, so allow only if there are no mismatches allowed */
    /* If equal to OLIGOSIZE, we would have found it using get_by_sarray */
    /* If <= 2*OLIGOSIZE, and there is a mismatch we won't find it, since at least one end will have a mismatch */
    debug(printf("Localdb_get using get_by_ends and no mismatches allowed.  Returning %d hits\n",nhits));
    *sortedp = true;
    *trimmedp = false;		/* Both ends match */
    if (localdb_memoize_enabled_p == true) {
      localdb_cache_insert(cache_hash,this,genomebits,query_compress,queryptr,pos5,pos3,querylength,
			   low_univdiagonal,high_univdiagonal,plusp,genestrand,
			   nhits,*matchlength,*sortedp,*trimmedp,univdiags_alloc);
    }
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_localdb_record("by_ends",queryptr,pos5,pos3,querylength,
				 low_univdiagonal,high_univdiagonal,plusp,genestrand,
				 nhits,*matchlength,*trimmedp,*sortedp,
				 Accel_profile_now_us() - profile_start_us,&profile_work);
#endif
    return nhits;

  } else {
    nhits = get_exhaustive(&(*trimmedp),&(*matchlength),&(*univdiags_alloc),this,
			   localdb_tally,localdb_overflowp,max_nblocks,
			   queryptr,pos5,pos3,querylength,low_univdiagonal,high_univdiagonal,
			   &profile_work);
    debug(printf("Localdb_get using get_exhaustive and %d mismatches allowed.  Returning %d hits\n",
		 (pos3 - pos5)/2,nhits));
    *sortedp = true;
    if (localdb_memoize_enabled_p == true) {
      localdb_cache_insert(cache_hash,this,genomebits,query_compress,queryptr,pos5,pos3,querylength,
			   low_univdiagonal,high_univdiagonal,plusp,genestrand,
			   nhits,*matchlength,*sortedp,*trimmedp,univdiags_alloc);
    }
#ifdef GSNAP_ACCEL_PROFILE_STATS
    Accel_profile_localdb_record("exhaustive",queryptr,pos5,pos3,querylength,
				 low_univdiagonal,high_univdiagonal,plusp,genestrand,
				 nhits,*matchlength,*trimmedp,*sortedp,
				 Accel_profile_now_us() - profile_start_us,&profile_work);
#endif
    return nhits;
  }
}


Univcoord_T
Localdb_get_one_low (Univdiag_T *univdiags_alloc,
		     T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		     int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		     char *queryptr, int querylength,
		     Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		     Compress_T query_compress, bool plusp, int genestrand,
		     Genomebits_T genomebits) {
  Univcoord_T lowest;
  bool sortedp, trimmedp;
  int matchlength;
  int nunivdiags, i;

  Accel_profile_localdb_site("mate_rescue_low");
  if ((nunivdiags = Localdb_get(&sortedp,&trimmedp,&matchlength,univdiags_alloc,
				this,diagonals_alloc,localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				queryptr,/*pos5*/0,/*pos3*/querylength,querylength,
				/*anchor_univdiagonal:low_univdiagonal,*/
				low_univdiagonal,high_univdiagonal,
				query_compress,plusp,genestrand,genomebits)) == 0) {
    return (Univcoord_T) 0;

  } else if (nunivdiags == 1) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  univdiags_alloc[0]->univdiagonal,matchlength,nmismatches));
    return univdiags_alloc[0]->univdiagonal;

  } else if (sortedp == true) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  univdiags_alloc[0]->univdiagonal,matchlength,nmismatches));
    return univdiags_alloc[0]->univdiagonal;

  } else {
    lowest = univdiags_alloc[0]->univdiagonal;
    for (i = 1; i < nunivdiags; i++) {
      if (univdiags_alloc[i]->univdiagonal < lowest) {
	lowest = univdiags_alloc[i]->univdiagonal;
      }
    }

    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  lowest,matchlength,nmismatches));
    return lowest;
  }
}


Univcoord_T
Localdb_get_one_high (Univdiag_T *univdiags_alloc, 
		      T this, Univcoord_T *diagonals_alloc, unsigned short *localdb_alloc,
		      int16_t *localdb_tally, int16_t *localdb_overflowp, int max_nblocks,
		      char *queryptr, int querylength,
		      Univcoord_T low_univdiagonal, Univcoord_T high_univdiagonal,
		      Compress_T query_compress, bool plusp, int genestrand,
		      Genomebits_T genomebits) {
  Univcoord_T highest;
  bool sortedp, trimmedp;
  int matchlength;
  int nunivdiags, i;

  Accel_profile_localdb_site("mate_rescue_high");
  if ((nunivdiags = Localdb_get(&sortedp,&trimmedp,&matchlength,univdiags_alloc,
				this,diagonals_alloc,localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
				queryptr,/*pos5*/0,/*pos3*/querylength,querylength,
				/*anchor_univdiagonal:high_univdiagonal,*/
				low_univdiagonal,high_univdiagonal,
				query_compress,plusp,genestrand,genomebits)) == 0) {
    return (Univcoord_T) 0;

  } else if (nunivdiags == 1) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  univdiags_alloc[0]->univdiagonal,matchlength,nmismatches));
    return univdiags_alloc[0]->univdiagonal;

  } else if (sortedp == true) {
    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  univdiags_alloc[nunivdiags - 1]->univdiagonal,matchlength,nmismatches));
    return univdiags_alloc[nunivdiags - 1]->univdiagonal;

  } else {
    highest = univdiags_alloc[0]->univdiagonal;
    for (i = 1; i < nunivdiags; i++) {
      if (univdiags_alloc[i]->univdiagonal > highest) {
	highest = univdiags_alloc[i]->univdiagonal;
      }
    }

    debug5(printf("Returning %u with matchlength %d, nmismatches %d\n",
		  highest,matchlength,nmismatches));
    return highest;
  }
}



T
Localdb_new (char *genomesubdir, char *fileroot,
	     Access_mode_T localdb_access, bool sharedp,
	     bool multiple_sequences_p, bool preload_shared_memory_p, bool unload_shared_memory_p) {
  T new;

  char *saindex16_filename, *sarray16_filename, *sarray8_filename, *sasort16_filename;
  char *comma;
  double seconds0, seconds1, seconds2;

  saindex16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".saindex16")+1)*sizeof(char));
  sarray16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".sarray16")+1)*sizeof(char));
  sarray8_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".sarray8")+1)*sizeof(char));
  sasort16_filename = (char *) MALLOC((strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".sasort16")+1)*sizeof(char));
  sprintf(saindex16_filename,"%s/%s.saindex16",genomesubdir,fileroot);
  sprintf(sarray16_filename,"%s/%s.sarray16",genomesubdir,fileroot);
  sprintf(sarray8_filename,"%s/%s.sarray8",genomesubdir,fileroot);
  sprintf(sasort16_filename,"%s/%s.sasort16",genomesubdir,fileroot);

  if (Access_file_exists_p(saindex16_filename) == false ||
      Access_file_exists_p(sarray16_filename) == false ||
      Access_file_exists_p(sarray8_filename) == false ||
      Access_file_exists_p(sasort16_filename) == false) {
    fprintf(stderr,"Unable to find localdb files %s, %s, %s, or %s.  Will align without them.\n",
	    saindex16_filename,sarray16_filename,sarray8_filename,sasort16_filename);
    FREE(sasort16_filename);
    FREE(sarray8_filename);
    FREE(sarray16_filename);
    FREE(saindex16_filename);
    return (T) NULL;

  } else if (localdb_access == USE_ALLOCATE) {
    new = (T) MALLOC(sizeof(*new));
    fprintf(stderr,"Allocating memory for localdb...");

    if (
#ifndef HAVE_MMAP
	0 &&
#endif
	multiple_sequences_p == false && preload_shared_memory_p == false && unload_shared_memory_p == false) {
#ifdef HAVE_MMAP
      new->saindex16 = (uint16_t *) Access_mmap(&new->saindex16_fd,&new->saindex16_len,&seconds0,
					       saindex16_filename,/*randomp*/false);
      new->sarray16 = (uint16_t *) Access_mmap(&new->sarray16_fd,&new->sarray16_len,&seconds1,
					       sarray16_filename,/*randomp*/false);
      new->sarray8 = (uint8_t *) Access_mmap(&new->sarray8_fd,&new->sarray8_len,&seconds2,
					     sarray8_filename,/*randomp*/false);
      new->sasort16 = (uint16_t *) Access_mmap(&new->sasort16_fd,&new->sasort16_len,&seconds1,
					       sasort16_filename,/*randomp*/false);
      new->localdb_access = MMAPPED;
#endif
    } else if (sharedp == true) {
      new->saindex16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->saindex16_shmid,&new->saindex16_key,
							  &new->saindex16_fd,&new->saindex16_len,&seconds0,
							  saindex16_filename,sizeof(uint16_t));
      new->sarray16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->sarray16_shmid,&new->sarray16_key,
							  &new->sarray16_fd,&new->sarray16_len,&seconds1,
							  sarray16_filename,sizeof(uint16_t));
      new->sarray8 = (uint8_t *) Access_allocate_shared(&new->localdb_access,&new->sarray8_shmid,&new->sarray8_key,
							  &new->sarray8_fd,&new->sarray8_len,&seconds2,
							  sarray8_filename,sizeof(uint8_t));
      new->sasort16 = (uint16_t *) Access_allocate_shared(&new->localdb_access,&new->sasort16_shmid,&new->sasort16_key,
							  &new->sasort16_fd,&new->sasort16_len,&seconds1,
							  sasort16_filename,sizeof(uint16_t));
    } else {
      new->saindex16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->saindex16_len,&seconds0,
							   saindex16_filename,sizeof(uint16_t));
      new->sarray16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->sarray16_len,&seconds1,
							   sarray16_filename,sizeof(uint16_t));
      new->sarray8 = (uint8_t *) Access_allocate_private(&new->localdb_access,&new->sarray8_len,&seconds2,
							 sarray8_filename,sizeof(uint8_t));
      new->sasort16 = (uint16_t *) Access_allocate_private(&new->localdb_access,&new->sasort16_len,&seconds1,
							   sasort16_filename,sizeof(uint16_t));
    }

    /* Don't check for sarray8, since it can be empty when genomelength < 65536 */
    if (new->saindex16 == NULL || new->sarray16 == NULL ||
	/*new->sarray8 == NULL ||*/ new->sasort16 == NULL) {
      fprintf(stderr,"insufficient memory\n");
      exit(9);
    } else {
      comma = Genomicpos_commafmt(new->saindex16_len + new->sarray16_len + new->sarray8_len + new->sasort16_len);
      fprintf(stderr,"done (%s bytes",comma);
      FREE(comma);

      if (multiple_sequences_p == true) {
	fprintf(stderr,", %.2f sec",seconds0 + seconds1 + seconds2);
      }
      fprintf(stderr,")\n");
    }

#ifdef HAVE_MMAP
  } else if (localdb_access == USE_MMAP_PRELOAD || localdb_access == USE_MMAP_ONLY) {
    fprintf(stderr,"Memory mapping localdb files...");
    new = (T) MALLOC(sizeof(*new));
    new->saindex16 = (uint16_t *) Access_mmap(&new->saindex16_fd,&new->saindex16_len,&seconds0,
					     saindex16_filename,/*randomp*/true);
    new->sarray16 = (uint16_t *) Access_mmap(&new->sarray16_fd,&new->sarray16_len,&seconds1,
					     sarray16_filename,/*randomp*/true);
    new->sarray8 = (uint8_t *) Access_mmap(&new->sarray8_fd,&new->sarray8_len,&seconds2,
					     sarray8_filename,/*randomp*/true);
    new->sasort16 = (uint16_t *) Access_mmap(&new->sasort16_fd,&new->sasort16_len,&seconds1,
					     sasort16_filename,/*randomp*/true);

    if (new->saindex16 == NULL || new->sarray16 == NULL || new->sarray8 == NULL || new->sasort16 == NULL) {
      fprintf(stderr,"Insufficient memory for mmap (will use disk file instead, but program will be slow)\n");
      new->localdb_access = FILEIO;
    } else {
      comma = Genomicpos_commafmt(new->saindex16_len + new->sarray16_len + new->sarray8_len + new->sasort16_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds0 + seconds1 + seconds2);
      FREE(comma);
      new->localdb_access = MMAPPED;
    }
#endif

  } else {
    fprintf(stderr,"Don't recognize localdb_access %d\n",localdb_access);
    abort();
  }

  FREE(sasort16_filename);
  FREE(sarray8_filename);
  FREE(sarray16_filename);
  FREE(saindex16_filename);

  return new;
}
