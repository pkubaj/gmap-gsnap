static char rcsid[] = "$Id: 73e9e01f2bacd365ee285ebcc81d671ce96a81be $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "localdb-write.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_COMPARISON_SORT
#include <math.h>		/* For qsort */
#endif

#include "bool.h"
#include "assert.h"
#include "mem.h"
#include "fopen.h"
#include "littleendian.h"

#ifdef USE_SACA
/* One disadvantage of SACA algorithm is it needs a sentinel character
   at position n, leading to the first elt of the suffix array being
   n, which looks like 0 */
#include "saca-k.h"
#endif


#define MAXN_UINT1 256U
#define MAXN_UINT2 65536U

#define MONITOR_INTERVAL 100000000 /* 100 million nt */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif 

/* make_saindex_uint2 */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif 



#define OLIGOSIZE 4 /* nucleotides, occupying 8 bits, to specify 256 pointers */
#define OLIGOSPACE 256




static char *oligos[] =
  {"AAAA","AAAC","AAAG","AAAT","AACA","AACC","AACG","AACT","AAGA","AAGC","AAGG","AAGT","AATA","AATC","AATG","AATT",
   "ACAA","ACAC","ACAG","ACAT","ACCA","ACCC","ACCG","ACCT","ACGA","ACGC","ACGG","ACGT","ACTA","ACTC","ACTG","ACTT",
   "AGAA","AGAC","AGAG","AGAT","AGCA","AGCC","AGCG","AGCT","AGGA","AGGC","AGGG","AGGT","AGTA","AGTC","AGTG","AGTT",
   "ATAA","ATAC","ATAG","ATAT","ATCA","ATCC","ATCG","ATCT","ATGA","ATGC","ATGG","ATGT","ATTA","ATTC","ATTG","ATTT",

   "CAAA","CAAC","CAAG","CAAT","CACA","CACC","CACG","CACT","CAGA","CAGC","CAGG","CAGT","CATA","CATC","CATG","CATT",
   "CCAA","CCAC","CCAG","CCAT","CCCA","CCCC","CCCG","CCCT","CCGA","CCGC","CCGG","CCGT","CCTA","CCTC","CCTG","CCTT",
   "CGAA","CGAC","CGAG","CGAT","CGCA","CGCC","CGCG","CGCT","CGGA","CGGC","CGGG","CGGT","CGTA","CGTC","CGTG","CGTT",
   "CTAA","CTAC","CTAG","CTAT","CTCA","CTCC","CTCG","CTCT","CTGA","CTGC","CTGG","CTGT","CTTA","CTTC","CTTG","CTTT",

   "GAAA","GAAC","GAAG","GAAT","GACA","GACC","GACG","GACT","GAGA","GAGC","GAGG","GAGT","GATA","GATC","GATG","GATT",
   "GCAA","GCAC","GCAG","GCAT","GCCA","GCCC","GCCG","GCCT","GCGA","GCGC","GCGG","GCGT","GCTA","GCTC","GCTG","GCTT",
   "GGAA","GGAC","GGAG","GGAT","GGCA","GGCC","GGCG","GGCT","GGGA","GGGC","GGGG","GGGT","GGTA","GGTC","GGTG","GGTT",
   "GTAA","GTAC","GTAG","GTAT","GTCA","GTCC","GTCG","GTCT","GTGA","GTGC","GTGG","GTGT","GTTA","GTTC","GTTG","GTTT",

   "TAAA","TAAC","TAAG","TAAT","TACA","TACC","TACG","TACT","TAGA","TAGC","TAGG","TAGT","TATA","TATC","TATG","TATT",
   "TCAA","TCAC","TCAG","TCAT","TCCA","TCCC","TCCG","TCCT","TCGA","TCGC","TCGG","TCGT","TCTA","TCTC","TCTG","TCTT",
   "TGAA","TGAC","TGAG","TGAT","TGCA","TGCC","TGCG","TGCT","TGGA","TGGC","TGGG","TGGT","TGTA","TGTC","TGTG","TGTT",
   "TTAA","TTAC","TTAG","TTAT","TTCA","TTCC","TTCG","TTCT","TTGA","TTGC","TTGG","TTGT","TTTA","TTTC","TTTG","TTTT"};


#ifdef COMPUTE_OLIGOS
#define LOW_TWO_BITS  0x03

#define BYTE_RIGHT_A 0x00
#define BYTE_RIGHT_C 0x01
#define BYTE_RIGHT_G 0x02
#define BYTE_RIGHT_T 0x03

static void
oligo_nt (char *nt, uint8_t oligo, int oligolength) {
  int i, j;
  uint8_t lowbits;

  nt[oligolength] = '\0';

  j = oligolength - 1;
  for (i = 0; i < oligolength; i++) {
    lowbits = oligo & LOW_TWO_BITS;
    switch (lowbits) {
    case BYTE_RIGHT_A: nt[j] = 'A'; break;
    case BYTE_RIGHT_C: nt[j] = 'C'; break;
    case BYTE_RIGHT_G: nt[j] = 'G'; break;
    case BYTE_RIGHT_T: nt[j] = 'T'; break;
    }
    oligo >>= 2;
    j--;
  }

  return;
}
#endif


/* Creates a set of 256 pointers, each an unsigned short, having a
   value from 0 to 65535, indicating where the 4-nt string starts in
   the suffix array */
static void
make_saindex_uint2 (unsigned short *saindex16, unsigned short *sarray16, int n,
		  char *gbuffer16) {
#ifdef COMPUTE_OLIGOS
  char oligo[OLIGOSIZE+1];
#else
  char *oligo;
#endif
  int i, k;
  
#ifdef DEBUG9
  for (i = 0; i < n; i++) {
    printf("%d\t%d\t%.*s\n",i,sarray16[i],OLIGOSIZE,&(gbuffer16[sarray16[i]]));
  }
  printf("\n");
#endif

  i = 0;
  for (k = 0; k < OLIGOSPACE; k++) {
#ifdef COMPUTE_OLIGOS
    oligo_nt(oligo,k,OLIGOSIZE);
#else
    oligo = oligos[k];
#endif
    while (i < n && strcmp(/*suffix*/&(gbuffer16[sarray16[i]]),oligo) < 0) {
      i++;
    }
    debug9(printf("saindex at %d is %d\n",k,i));
    saindex16[k] = i;	   /* Should point to first suffix >= oligo */
  }

  return;
}


static int
ushort_cmp (const void *x, const void *y) {
  unsigned short a = * (unsigned short *) x;
  unsigned short b = * (unsigned short *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return +1;
  } else {
    return 0;
  }
}

/* For each of the 256 pointers in sarray16, sorts the values,
   essentially creating a hash table of 4-mers for the block of 65536
   positions */

static void
sarray_sort_uint2 (unsigned short *sarray16, unsigned short *saindex16, int n) {
  int initptr, finalptr, nitems;
  unsigned char oligoi;

  for (oligoi = 0; oligoi < OLIGOSPACE - 1; oligoi++) {
    initptr = saindex16[oligoi];
    finalptr = saindex16[oligoi+1];
    if ((nitems = finalptr - initptr) > 0) {
      qsort(&(sarray16[initptr]),nitems,sizeof(unsigned short),ushort_cmp);
    }
  }

  initptr = saindex16[OLIGOSPACE - 1];
  finalptr = n;
  if ((nitems = finalptr - initptr) > 0) {
    qsort(&(sarray16[initptr]),nitems,sizeof(unsigned short),ushort_cmp);
  }
  
  return;
}


#ifdef USE_SACA
static void
convert_to_bytes (uint8_t *dest, UINT4 *src, int n) {
  int i;

  for (i = 0; i < n; i++) {
    dest[i] = (uint8_t) src[i];
  }

  return;
}

static void
convert_to_shorts (uint16_t *dest, UINT4 *src, int n) {
  int i;

  for (i = 0; i < n; i++) {
    dest[i] = (uint16_t) src[i];
  }

  return;
}
#endif


typedef struct Rankpair_T *Rankpair_T;
struct Rankpair_T {
  int rank0;
  int rank1;
  int index;
  int order;
};

static void
Rankpair_initialize_ranks (Rankpair_T *rankpair, unsigned char *intstring, int n) {
  int *count = (int *) CALLOC(256+1,sizeof(int)); /* 256 possible values for a char */
  int i, k;

  for (i = 0; i < n; i++) {
    count[1 + (int) intstring[i]]++; /* Add 1, so count at the given character is the rank */
  }
  for (k = 1; k <= 256; k++) {
    count[k] += count[k-1];	/* cumsum */
  }
  for (i = 0; i < n; i++) {
    rankpair[i]->rank0 = count[(int) intstring[i]];
    rankpair[i]->rank1 = (i + 1 >= n) ? -1 : count[(int) intstring[i+1]];
    rankpair[i]->index = i;
    debug(rankpair[i]->order = 0);
  }

  FREE(count);
  return;
}


static int
Rankpair_rank_cmp (const void *a, const void *b) {
  Rankpair_T x = * (Rankpair_T *) a;
  Rankpair_T y = * (Rankpair_T *) b;

  if (x->rank0 < y->rank0) {
    return -1;
  } else if (y->rank0 < x->rank0) {
    return +1;
  } else if (x->rank1 < y->rank1) {
    return -1;
  } else if (y->rank1 < x->rank1) {
    return +1;
  } else {
    return 0;
  }
}

#ifndef USE_COMPARISON_SORT
/* Radix sort, O(n) */
static Rankpair_T *
Rankpair_sort_byranks (Rankpair_T *rankpair, int n) {
  Rankpair_T *sorted = (Rankpair_T *) MALLOC(n*sizeof(Rankpair_T));
  Rankpair_T *temp;
  int *count = (int *) CALLOC((n+1),sizeof(n));
  int i, k;

  /* Sort by rank1, which can go from -1 to (n-1), so we add 1 */
  for (i = 0; i < n; i++) {
    count[1 + rankpair[i]->rank1]++;
  }
  for (k = 1; k <= n; k++) {
    count[k] += count[k-1];	/* cumsum */
  }
  for (i = n-1; i >= 0; i--) {
    sorted[--count[1 + rankpair[i]->rank1]] = rankpair[i];
  }
  
  /* Swap rankpair and sorted */
  temp = rankpair;
  rankpair = sorted;
  sorted = temp;

  /* Sort by rank0, which can go from 0 to (n-1) */
  memset(count,0,n*sizeof(int));
  for (i = 0; i < n; i++) {
    count[0 + rankpair[i]->rank0]++;
  }
  for (k = 1; k <= n; k++) {
    count[k] += count[k-1];	/* cumsum */
  }
  for (i = n-1; i >= 0; i--) {
    sorted[--count[0 + rankpair[i]->rank0]] = rankpair[i];
  }

  FREE(count);
  FREE(rankpair);
  
  return sorted;
}
#endif


#ifdef USE_COMPARISON_SORT
/* Puts items back in their original order */
static int
Rankpair_index_cmp (const void *a, const void *b) {
  Rankpair_T x = * (Rankpair_T *) a;
  Rankpair_T y = * (Rankpair_T *) b;

  if (x->index < y->index) {
    return -1;
  } else if (y->index < x->index) {
    return +1;
  } else {
    return 0;
  }
}

#else
/* Radix sort, O(n) */
/* Assumes there is exactly one index from 0 to n-1 */
static Rankpair_T *
Rankpair_sort_byindex (Rankpair_T *rankpair, int n) {
  Rankpair_T *sorted = (Rankpair_T *) MALLOC(n*sizeof(Rankpair_T));
  int i;

  for (i = 0; i < n; i++) {
    sorted[rankpair[i]->index] = rankpair[i];
  }
  FREE(rankpair);
  
  return sorted;
}
#endif


#ifdef DEBUG
static void
print_rankpairs (Rankpair_T *rankpair, int n) {
  int i;

  for (i = 0; i < n; i++) {
    printf("%d,%d => index %d => order %d\n",
	   rankpair[i]->rank0,rankpair[i]->rank1,rankpair[i]->index,rankpair[i]->order);
  }
  printf("\n");

  return;
}
#endif


static Rankpair_T *
sarray_compute (Rankpair_T *rankpair, unsigned char *intstring, int n) {
  int order, kmer, i;
  bool donep = false;

  kmer = 1;
  Rankpair_initialize_ranks(rankpair,intstring,n);
  debug(print_rankpairs(rankpair,n));

#ifdef USE_COMPARISON_SORT
  qsort(rankpair,n,sizeof(Rankpair_T),Rankpair_rank_cmp);
#else
  rankpair = Rankpair_sort_byranks(rankpair,n);
#endif

  /* Assign order, accounting for duplicates */
  rankpair[0]->order = order = 0;
  for (i = 1; i < n; i++) {
    if (Rankpair_rank_cmp(&(rankpair[i]),&(rankpair[i-1])) == 0) {
      /* Tie */
      rankpair[i]->order = order;
    } else {
      rankpair[i]->order = ++order;
    }
  }
  debug(printf(">kmer %d\n",kmer));
  debug(print_rankpairs(rankpair,n));
#ifdef USE_COMPARISON_SORT
  qsort(rankpair,n,sizeof(Rankpair_T),Rankpair_index_cmp);
#else
  rankpair = Rankpair_sort_byindex(rankpair,n);
#endif
  debug(print_rankpairs(rankpair,n));

  kmer = 2;
  while (donep == false) {
    /* Assign ranks based on previous order */
    for (i = 0; i < n; i++) {
      rankpair[i]->rank0 = rankpair[i]->order;
      rankpair[i]->rank1 = (i + kmer >= n) ? -1 : rankpair[i + kmer]->order;
      rankpair[i]->index = i;
    }
#ifdef USE_COMPARISON_SORT
    qsort(rankpair,n,sizeof(Rankpair_T),Rankpair_rank_cmp);
#else
    rankpair = Rankpair_sort_byranks(rankpair,n);
#endif
  
    /* Assign order, accounting for duplicates */
    rankpair[0]->order = order = 0;
    for (i = 1; i < n; i++) {
      if (Rankpair_rank_cmp(&(rankpair[i]),&(rankpair[i-1])) == 0) {
	/* Tie */
	rankpair[i]->order = order;
      } else {
	rankpair[i]->order = ++order;
      }
    }
    if (order == n - 1) {
      donep = true;
    }

    debug(printf(">kmer %d\n",kmer));
    debug(print_rankpairs(rankpair,n));
#ifdef USE_COMPARISON_SORT
    qsort(rankpair,n,sizeof(Rankpair_T),Rankpair_index_cmp);
#else
    rankpair = Rankpair_sort_byindex(rankpair,n);
#endif
    debug(print_rankpairs(rankpair,n));

    kmer *= 2;
  }

  return rankpair;
}


static void
sarray_build_uint2 (unsigned short *sarray, unsigned char *intstring, int n,
		    struct Rankpair_T *space) {
  Rankpair_T *rankpair;
  int i;

  rankpair = (Rankpair_T *) MALLOC(MAXN_UINT2*sizeof(Rankpair_T));
  for (i = 0; i < n; i++) {
    rankpair[i] = &(space[i]);
  }

  rankpair = sarray_compute(rankpair,intstring,n);

  for (i = 0; i < n; i++) {
    sarray[rankpair[i]->order] = (unsigned short) i; /* rankpair[i].index; */
  }

  FREE(rankpair);
  return;
}


static void
sarray_build_uint1 (unsigned char *sarray, unsigned char *intstring, int n,
		    struct Rankpair_T *space) {
  Rankpair_T *rankpair;
  int i;

  rankpair = (Rankpair_T *) MALLOC(MAXN_UINT1*sizeof(Rankpair_T));
  for (i = 0; i < n; i++) {
    rankpair[i] = &(space[i]);
  }

  rankpair = sarray_compute(rankpair,intstring,n);

  for (i = 0; i < n; i++) {
    sarray[rankpair[i]->order] = (unsigned char) i; /* rankpair[i].index; */
  }

  FREE(rankpair);
  return;
}


void
Localdb_write (char *saindex16file, char *sarray16file,
	       char *sarray8file, char *sasort16file,
	       Genome_T genomecomp, Univcoord_T genomelength) {
  Univcoord_T left, n;
#ifdef USE_SACA
  UINT4 *SA8, *SA16;
#else
  struct Rankpair_T *space_uint1, *space_uint2;
#endif

  unsigned short *saindex16;
  unsigned short *sarray16;
  unsigned char *sarray8;
  unsigned char *intstring16, *intstring8;
  char *gbuffer16;
  FILE *fp0, *fp1, *fp2, *fp3;

  Univcoord_T monitor = MONITOR_INTERVAL;
  char *comma;


  if ((fp0 = FOPEN_WRITE_BINARY(saindex16file)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",saindex16file);
    exit(9);
  }
  if ((fp1 = FOPEN_WRITE_BINARY(sarray16file)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarray16file);
    exit(9);
  }
  if ((fp2 = FOPEN_WRITE_BINARY(sarray8file)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sarray8file);
    exit(9);
  }
  if ((fp3 = FOPEN_WRITE_BINARY(sasort16file)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",sasort16file);
    exit(9);
  }

  
#ifdef USE_SACA
  SA16 = (UINT4 *) MALLOC((MAXN_UINT2+1)*sizeof(UINT4));
  SA8 = (UINT4 *) MALLOC((MAXN_UINT1+1)*sizeof(UINT4));
#else
  space_uint1 = (struct Rankpair_T *) MALLOC(MAXN_UINT1*sizeof(struct Rankpair_T));
  space_uint2 = (struct Rankpair_T *) MALLOC(MAXN_UINT2*sizeof(struct Rankpair_T));
#endif

  saindex16 = (unsigned short *) MALLOC(OLIGOSPACE*sizeof(unsigned short));
  sarray16 = (unsigned short *) MALLOC(MAXN_UINT2*sizeof(unsigned short));
  sarray8 = (unsigned char *) MALLOC(MAXN_UINT1*sizeof(unsigned char));

  gbuffer16 = (char *) CALLOC(MAXN_UINT2+1,sizeof(char));
  intstring16 = (unsigned char *) CALLOC(MAXN_UINT2+1,sizeof(unsigned char));
  intstring8 = (unsigned char *) CALLOC(MAXN_UINT1+1,sizeof(unsigned char));

  left = 0;
  while (left + MAXN_UINT2 < genomelength) {
    if (left > monitor) {
      comma = Genomicpos_commafmt(monitor);
      fprintf(stderr,"Computed regiondb to position %s\n",comma);
      FREE(comma);
      monitor += MONITOR_INTERVAL;
    }

    /* Full large block */
    Genome_fill_buffer_int_string(genomecomp,left,/*length*/MAXN_UINT2,
				  intstring16,/*conversion*/NULL);
#ifdef USE_SACA
    intstring16[MAXN_UINT2] = 0;		       /* Terminator.  Tried N/X, but SACA_K fails */
    SACA_K(intstring16,SA16,MAXN_UINT2+/*virtual sentinel*/1,/*K, alphabet_size*/5,
	   /*m*/MAXN_UINT2+1,/*level*/0);
    convert_to_shorts(sarray16,SA16,MAXN_UINT2);
#else
    /* Suffix array */
    sarray_build_uint2(sarray16,intstring16,MAXN_UINT2,space_uint2);
#endif
    FWRITE_USHORTS(sarray16,MAXN_UINT2,fp1);

    /* Index of 4-mers */
    Genome_fill_buffer_simple(genomecomp,left,/*length*/MAXN_UINT2,gbuffer16);
    make_saindex_uint2(saindex16,sarray16,MAXN_UINT2,gbuffer16);
    FWRITE_USHORTS(saindex16,OLIGOSPACE,fp0);

    /* Sorted 4-mer array */
    sarray_sort_uint2(sarray16,saindex16,MAXN_UINT2);
    FWRITE_USHORTS(sarray16,MAXN_UINT2,fp3);


    left += MAXN_UINT2;
    if (left + 128 < genomelength) {
      /* Full small block */
      Genome_fill_buffer_int_string(genomecomp,left - 128,/*length*/MAXN_UINT1,
				    intstring8,/*conversion*/NULL);
#ifdef USE_SACA
      intstring8[MAXN_UINT1] = 0;		     /* Terminator */
      SACA_K(intstring8,SA8,MAXN_UINT1+/*virtual sentinel*/1,/*K, alphabet_size*/5,
	     /*m*/MAXN_UINT1+1,/*level*/0);
      convert_to_bytes(sarray8,SA8,MAXN_UINT1);
#else
      /* Suffix array */
      sarray_build_uint1(sarray8,intstring8,MAXN_UINT1,space_uint1);
#endif
      FWRITE_CHARS(sarray8,MAXN_UINT1,fp2);

    } else {
      /* Partial small block */
      n = genomelength - (left - 128);
      Genome_fill_buffer_int_string(genomecomp,left - 128,/*length*/n,
				    intstring8,/*conversion*/NULL);
#ifdef USE_SACA
      intstring8[n] = 0;		     /* Terminator */
      SACA_K(intstring8,SA8,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,
	     /*m*/n+1,/*level*/0);
      convert_to_bytes(sarray8,SA8,n);
#else
      /* Suffix array */
      sarray_build_uint1(sarray8,intstring8,n,space_uint1);
#endif
      FWRITE_CHARS(sarray8,n,fp2);
    }
  }

  if ((n = genomelength - left) > 0) {
    Genome_fill_buffer_int_string(genomecomp,left,/*length*/n,
				  intstring16,/*conversion*/NULL);
#ifdef USE_SACA
    intstring16[n] = 0;		       /* Terminator.  Tried N/X, but SACA_K fails */
    SACA_K(intstring16,SA16,n+/*virtual sentinel*/1,/*K, alphabet_size*/5,
	   /*m*/n+1,/*level*/0);
    convert_to_shorts(sarray16,SA16,n);
#else
    /* Suffix array */
    sarray_build_uint2(sarray16,intstring16,n,space_uint2);
#endif
    FWRITE_USHORTS(sarray16,n,fp1);

    /* Index of 4-mers */
    Genome_fill_buffer_simple(genomecomp,left,/*length*/n,gbuffer16);
    make_saindex_uint2(saindex16,sarray16,n,gbuffer16);
    FWRITE_USHORTS(saindex16,OLIGOSPACE,fp0);

    /* Sorted 4-mer array */
    sarray_sort_uint2(sarray16,saindex16,n);
    FWRITE_USHORTS(sarray16,n,fp3);
  }

  FREE(intstring8);
  FREE(intstring16);
  FREE(gbuffer16);

  FREE(sarray8);
  FREE(sarray16);
  FREE(saindex16);

#ifdef USE_SACA
  FREE(SA8);
  FREE(SA16);
#else
  FREE(space_uint1);
  FREE(space_uint2);
#endif

  fclose(fp3);
  fclose(fp2);
  fclose(fp1);
  fclose(fp0);

  return;
}


#if 0
int
main (int argc, char *argv[]) {
  char *string = "mississippi";
  int n = strlen(string), i;
  unsigned short *sarray = malloc(n*sizeof(unsigned short));

  sarray_build_uint2(sarray,string,n);
    
  for (i = 0; i < n; i++) {
    printf("%d ",sarray[i]);
  }
  printf("\n");

  return 0;
}
#endif



