static char rcsid[] = "$Id: compress-write.c 226315 2023-02-28 18:22:20Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "compress-write.h"

#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>		/* For isalpha, toupper */
#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For off_t */
#endif
#include "complement.h"
#include "mem.h"		/* For Compress_new */
#include "getline.h"


/* Compress_cat */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Another MONITOR_INTERVAL is in indexdb.c */
#define MONITOR_INTERVAL 10000000 /* 10 million nt */

#define MAX_BADCHAR_MESSAGES 10
#define BADCHAR_INTERVAL 1000000


static char uppercaseCode[128] = UPPERCASE_U2T;

/* We use int *, rather than char *, because we eventually return an int,
   and we see problems converting from char to int */
static void
fill_buffer (int *Buffer, Genomecomp_T high, Genomecomp_T low, Genomecomp_T flags, Univcoord_T position) {
  int i;

  /* printf("%08X %08X %08X => ",high,low,flags); */
  for (i = 0; i < 16; i++) {
    switch (low & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    low >>= 2;
  }
  for ( ; i < 32; i++) {
    switch (high & 3U) {
    case 0U: Buffer[i] = 'A'; break;
    case 1U: Buffer[i] = 'C'; break;
    case 2U: Buffer[i] = 'G'; break;
    case 3U: Buffer[i] = 'T'; break;
    default: abort();
    }
    high >>= 2;
  }
  for (i = 0; i < 32; i++) {
    if ((flags & 1U) == 1U) {
      if (Buffer[i] == 'A') {
	Buffer[i] = 'N';
      } else if (Buffer[i] == 'T') {
	Buffer[i] = 'X';
      } else {
#if defined(LARGE_GENOMES) || (defined(UTILITYP) && defined(HAVE_64_BIT))
	printf("Parsing error; saw non-ACGT flag plus %c at position %llu = %llu + %d\n",
	       Buffer[i],(unsigned long long) (position+i),(unsigned long long) position,i);
#else
	printf("Parsing error; saw non-ACGT flag plus %c at position %u = %u + %d\n",
	       Buffer[i],(unsigned int) (position+i),(unsigned int) position,i);
#endif
	printf("flags: %08X\n",flags);
	abort();
      }
    }
    flags >>= 1;
  }

  return;
}


/* Based on genomecomp */
int
Compress_get_char (FILE *sequence_fp, Univcoord_T position, bool uncompressedp) {
  Genomecomp_T high, low, flags;
  static int SAVEBUFFER[32];
  int ptr, c;

  if (uncompressedp == true) {
    while ((c = fgetc(sequence_fp)) != EOF && isspace(c)) {
    }
    if (c == EOF) {
      return EOF;
    } else {
      return c;
    }
  } else if ((ptr = position % 32) == 0) {
    if (FREAD_UINT(&high,sequence_fp) <= 0 ||
	FREAD_UINT(&low,sequence_fp) <= 0 ||
	FREAD_UINT(&flags,sequence_fp) <= 0) {
      return EOF;
    } else {
      fill_buffer(SAVEBUFFER,high,low,flags,position);
      return SAVEBUFFER[0];
    }
  } else {
    return SAVEBUFFER[ptr];
  }
}


/************************************************************************
 *   Compression and uncompression of the genome
 ************************************************************************/

/*                       87654321 */
#define UINT4_LEFT_A   0x00000000
#define UINT4_LEFT_C   0x40000000
#define UINT4_LEFT_G   0x80000000
#define UINT4_LEFT_T   0xC0000000
#define UINT4_LEFT_BIT 0x80000000

/*                     87654321 */
#define UINT4_LEFT_0 0x00000000
#define UINT4_LEFT_1 0x40000000
#define UINT4_LEFT_2 0x80000000
#define UINT4_LEFT_3 0xC0000000


/* Genomecomp format */
/* A = 000  Stored as 00 in first two bytes, 0 in flag byte
   C = 001            01                     0
   G = 010            10                     0
   T = 011            11                     0
   N = 100            00                     1
   X = 111            11                     1
*/


/*                         87654321 */
#define UINT4_LEFT_SET   0x80000000
#define UINT4_LEFT_CLEAR 0x00000000


/* Genome128 format */
/*          High bit   Low bit   Flag bit
   A = 000     0          0         0
   C = 001     0          1         0
   G = 010     1          0         0
   T = 011     1          1         1
   N = 100     0          0         1
   X = 111     1          1         1
*/




/************************************************************************/

static void
genomecomp_move_absolute (FILE *fp, Univcoord_T ptr) {
#ifdef HAVE_FSEEKO
  off_t offset = ptr*((off_t) sizeof(Genomecomp_T));

  if (fseeko(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#else
  long int offset = ptr*((long int) sizeof(Genomecomp_T));

  if (fseek(fp,offset,SEEK_SET) < 0) {
    perror("Error in gmapindex, genomecomp_move_absolute");
    exit(9);
  }
#endif

  return;
}

static void
genomecomp_read_current (Genomecomp_T *high, Genomecomp_T *low, Genomecomp_T *flags, FILE *fp,
			 int index1part) {
  char section[15];

  if (fread(section,sizeof(char),index1part,fp) < (unsigned int) index1part) {
    *high = 0xFFFFFFFF;
    *low = 0xFFFFFFFF;
    *flags = 0xFFFFFFFF;
    return;
  }

  *high = (section[3] & 0xff);
  *high <<= 8;
  *high |= (section[2] & 0xff);
  *high <<= 8;
  *high |= (section[1] & 0xff);
  *high <<= 8;
  *high |= (section[0] & 0xff);

  *low = (section[7] & 0xff);
  *low <<= 8;
  *low |= (section[6] & 0xff);
  *low <<= 8;
  *low |= (section[5] & 0xff);
  *low <<= 8;
  *low |= (section[4] & 0xff);

  *flags = (section[11] & 0xff);
  *flags <<= 8;
  *flags |= (section[10] & 0xff);
  *flags <<= 8;
  *flags |= (section[9] & 0xff);
  *flags <<= 8;
  *flags |= (section[8] & 0xff);

  return;
}


static void
write_compressed_one (FILE *fp, int *nbadchars, char Buffer[], Univcoord_T position) {
  Genomecomp_T high = 0U, low = 0U, flags = 0U, carry;
  int i;

  for (i = 0; i < 32; i++) {
    carry = high & 3U;
    high >>= 2;
    low >>= 2;
    flags >>= 1;
    switch (carry) {
    case 0U: break;
    case 1U: low |= UINT4_LEFT_C; break;
    case 2U: low |= UINT4_LEFT_G; break;
    case 3U: low |= UINT4_LEFT_T; break;
    default: abort();
    }

    switch (uppercaseCode[(int) Buffer[i]]) {
    case 'A': break;
    case 'C': high |= UINT4_LEFT_C; break;
    case 'G': high |= UINT4_LEFT_G; break;
    case 'T': high |= UINT4_LEFT_T; break;
    case 'N': flags |= UINT4_LEFT_BIT; break;
    case 'X': high |= UINT4_LEFT_T; flags |= UINT4_LEFT_BIT; break;
    default: 
      (*nbadchars) += 1;
      if (*nbadchars < MAX_BADCHAR_MESSAGES) {
#if defined(LARGE_GENOMES) || (defined(UTILITYP) && defined(HAVE_64_BIT))
	fprintf(stderr,"Don't recognize character %c at position %llu.  Using N instead\n",
		Buffer[i],(unsigned long long) (position+i));
#else
	fprintf(stderr,"Don't recognize character %c at position %u.  Using N instead\n",
		Buffer[i],(unsigned int) (position+i));
#endif
      } else if (*nbadchars == MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Too many non-recognizable characters.  Not reporting each individual occurrence anymore.\n");
      } else if ((*nbadchars) % BADCHAR_INTERVAL == 0) {
	fprintf(stderr,"A total of %d non-ACGTNX characters seen so far.\n",*nbadchars);
      }
      flags |= UINT4_LEFT_BIT;
      break;
    }
  }
  
  FWRITE_UINT(high,fp);
  FWRITE_UINT(low,fp);
  FWRITE_UINT(flags,fp);
  
  return;
}


static void
put_compressed_one (Genomecomp_T *sectioncomp, int *nbadchars, char Buffer[], Univcoord_T position) {
  Genomecomp_T high = 0U, low = 0U, flags = 0U, carry;
  int i;

  for (i = 0; i < 32; i++) {
    carry = high & 3U;
    high >>= 2;
    low >>= 2;
    flags >>= 1;
    switch (carry) {
    case 0U: break;
    case 1U: low |= UINT4_LEFT_C; break;
    case 2U: low |= UINT4_LEFT_G; break;
    case 3U: low |= UINT4_LEFT_T; break;
    default: abort();
    }

    switch (uppercaseCode[(int) Buffer[i]]) {
    case 'A': break;
    case 'C': high |= UINT4_LEFT_C; break;
    case 'G': high |= UINT4_LEFT_G; break;
    case 'T': high |= UINT4_LEFT_T; break;
    case 'N': flags |= UINT4_LEFT_BIT; break;
    case 'X': high |= UINT4_LEFT_T; flags |= UINT4_LEFT_BIT; break;
    default: 
      (*nbadchars) += 1;
      if (*nbadchars < MAX_BADCHAR_MESSAGES) {
#if defined(LARGE_GENOMES) || (defined(UTILITYP) && defined(HAVE_64_BIT))
	fprintf(stderr,"Don't recognize character %c at position %llu.  Using N instead\n",Buffer[i],position+i);
#else
	fprintf(stderr,"Don't recognize character %c at position %u.  Using N instead\n",Buffer[i],position+i);
#endif
      } else if (*nbadchars == MAX_BADCHAR_MESSAGES) {
	fprintf(stderr,"Too many non-recognizable characters.  Not reporting each individual occurrence anymore.\n");
      } else if ((*nbadchars) % BADCHAR_INTERVAL == 0) {
	fprintf(stderr,"A total of %d non-ACGTNX characters seen so far.\n",*nbadchars);
      }
      flags |= UINT4_LEFT_BIT;
      break;
    }
  }
  
  sectioncomp[0] = high;
  sectioncomp[1] = low;
  sectioncomp[2] = flags;
  
  return;
}


static char acgt[4] = {'A','C','G','T'};
static char non_acgt[4] = {'N','?','?','X'};

/* if gbuffer is NULL, then we fill with X's */
/* Based on genomecomp.  Version for genome128 not implemented yet */
int
Compress_update_file (int nbadchars, FILE *fp, char *gbuffer, Univcoord_T startpos,
		      Univcoord_T endpos, int index1part) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  unsigned int startdiscard, enddiscard, i;
  Genomecomp_T high, low, flags;
  char Buffer[32];


  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp,index1part);

    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,&nbadchars,Buffer,startpos);

  } else {

    genomecomp_move_absolute(fp,ptr);
    genomecomp_read_current(&high,&low,&flags,fp,index1part);

    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    genomecomp_move_absolute(fp,ptr);
    write_compressed_one(fp,&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      write_compressed_one(fp,&nbadchars,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      genomecomp_move_absolute(fp,ptr);
      genomecomp_read_current(&high,&low,&flags,fp,index1part);

      for (i = 0; i < 16; i++) {
	Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      genomecomp_move_absolute(fp,ptr);
      write_compressed_one(fp,&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}


int
Compress_update_memory (int nbadchars, Genomecomp_T *genomecomp, char *gbuffer, Univcoord_T startpos,
			Univcoord_T endpos) {
  /* Chrpos_T length = endpos - startpos; */
  Univcoord_T startblock, endblock, ptr;
  Genomecomp_T high, low, flags;
  char Buffer[32];
  unsigned int startdiscard, enddiscard, i;


  ptr = startblock = startpos/32U*3;
  endblock = endpos/32U*3;
  startdiscard = startpos % 32;
  enddiscard = endpos % 32;
  
  if (endblock == startblock) {
    /* Special case */
    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

    /* Fill Buffer with original contents */
    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < enddiscard; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);

  } else {

    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];

    /* Fill Buffer with original contents */
    for (i = 0; i < 16; i++) {
      Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
      low >>= 2;
      flags >>= 1;
    }
    for ( ; i < 32; i++) {
      Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
      high >>= 2;
      flags >>= 1;
    }
    for (i = startdiscard; i < 32; i++) {
      Buffer[i] = gbuffer ? *gbuffer++ : 'X';
    }
    put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,startpos);
    ptr += 3;
      
    while (ptr < endblock) {
      for (i = 0; i < 32; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
      ptr += 3;
    }

    if (enddiscard > 0) {
      high = genomecomp[ptr];
      low = genomecomp[ptr+1];
      flags = genomecomp[ptr+2];

      /* Fill Buffer with original contents */
      for (i = 0; i < 16; i++) {
	Buffer[i] = flags & 1U ? non_acgt[low & 3U] : acgt[low & 3U];
	low >>= 2;
	flags >>= 1;
      }
      for ( ; i < 32; i++) {
	Buffer[i] = flags & 1U ? non_acgt[high & 3U] : acgt[high & 3U];
	high >>= 2;
	flags >>= 1;
      }
      for (i = 0; i < enddiscard; i++) {
	Buffer[i] = gbuffer ? *gbuffer++ : 'X';
      }
      put_compressed_one(&(genomecomp[ptr]),&nbadchars,Buffer,ptr/3*32U);
    }
  }

  return nbadchars;
}


void
Compress_compress (char **files, int nfiles, bool stdin_p) {
  Genomecomp_T low = 0U, high = 0U, flags = 0U, carry;
  Univcoord_T position = 0UL;
  FILE *fp;
  char *line, *p;
  int c;
  int in_counter = 0, filei;


  fprintf(stderr,"Compressing genome");

  if (stdin_p == true) {
    nfiles = 1;
  }
  for (filei = 0; filei < nfiles; filei++) {
    if (stdin_p == true) {
      fp = stdin;
    } else if ((fp = fopen(files[filei],"r")) == NULL) {
      fprintf(stderr,"Could not open file %s\n",files[filei]);
      exit(9);
    }

    while ((line = Getline(fp)) != NULL) {
      if (line[0] == '>') {
	/* Skip header */
      } else {
	p = line;
	while ((c = (int) *p++) != '\0') {
	  if (isalpha(c)) {
	    in_counter++;

	    carry = high & 3U;
	    high >>= 2;
	    low >>= 2;
	    flags >>= 1;
	    switch (carry) {
	    case 0U: break;
	    case 1U: low |= UINT4_LEFT_C; break;
	    case 2U: low |= UINT4_LEFT_G; break;
	    case 3U: low |= UINT4_LEFT_T; break;
	    default: abort();
	    }

	    switch (uppercaseCode[c]) {
	    case 'A': break;
	    case 'C': high |= UINT4_LEFT_C; break;
	    case 'G': high |= UINT4_LEFT_G; break;
	    case 'T': high |= UINT4_LEFT_T; break;
	    case 'N': flags |= UINT4_LEFT_BIT; break;
	    case 'X': high |= UINT4_LEFT_T; flags |= UINT4_LEFT_BIT; break;
	    default: 
#if defined(LARGE_GENOMES) || (defined(UTILITYP) && defined(HAVE_64_BIT))
	      fprintf(stderr,"Non-standard nucleotide %c at position %llu.  Using N instead\n",c,position);
#else
	      fprintf(stderr,"Non-standard nucleotide %c at position %u.  Using N instead\n",c,position);
#endif
	      flags |= UINT4_LEFT_BIT;
	      break;
	    }
      
	    /* 8 is simply bits per byte */
	    if (in_counter == 8 * (int) sizeof(Genomecomp_T)) {
	      FWRITE_UINT(high,stdout);
	      FWRITE_UINT(low,stdout);
	      FWRITE_UINT(flags,stdout);
	    
	      low = high = flags = 0U;
	      in_counter = 0;
	    }

	    position++;
	    if (position % MONITOR_INTERVAL == 0) {
	      fprintf(stderr,".");
	    }
	  }
	}
      }

      FREE(line);
    }

    if (stdin_p == false) {
      fclose(fp);
    }
  }

  if (in_counter > 0) {
    while (in_counter < 8 * (int) sizeof(Genomecomp_T)) {
      carry = high & 3U;
      high >>= 2;
      low >>= 2;
      flags >>= 1;
      switch (carry) {
      case 0U: break;
      case 1U: low |= UINT4_LEFT_C; break;
      case 2U: low |= UINT4_LEFT_G; break;
      case 3U: low |= UINT4_LEFT_T; break;
      default: abort();
      }
      high |= UINT4_LEFT_T; flags |= UINT4_LEFT_BIT; /* Create X's at end */
      in_counter++;
    }

    FWRITE_UINT(high,stdout);
    FWRITE_UINT(low,stdout);
    FWRITE_UINT(flags,stdout);
  }

  /* Plus extra high and low */
  high = low = 0xFFFFFFFF;
  FWRITE_UINT(high,stdout);
  FWRITE_UINT(low,stdout);

  fprintf(stderr,"\n");

  return;
}



#ifdef HAVE_64_BIT
static inline void
nt_unshuffle (UINT4 *highbits, UINT4 *lowbits, UINT4 high, UINT4 low) {
  UINT8 x, t;

  x = (UINT8) high;
  x <<= 32;
  x |= low;

  t = (x ^ (x >> 1))  & 0x2222222222222222;  x = x ^ t ^ (t << 1);
  t = (x ^ (x >> 2))  & 0x0C0C0C0C0C0C0C0C;  x = x ^ t ^ (t << 2);
  t = (x ^ (x >> 4))  & 0x00F000F000F000F0;  x = x ^ t ^ (t << 4);
  t = (x ^ (x >> 8))  & 0x0000FF000000FF00;  x = x ^ t ^ (t << 8);
  t = (x ^ (x >> 16)) & 0x00000000FFFF0000;  x = x ^ t ^ (t << 16);

  *highbits = (UINT4) (x >> 32);
  *lowbits = (UINT4) x;

  return;
}

#else

static inline void
nt_unshuffle (UINT4 *highbits, UINT4 *lowbits, UINT4 high, UINT4 low) {
  UINT4 t;

  /* unshuffle high */
  t = (high ^ (high >> 1)) & 0x22222222;  high = high ^ t ^ (t << 1);
  t = (high ^ (high >> 2)) & 0x0C0C0C0C;  high = high ^ t ^ (t << 2);
  t = (high ^ (high >> 4)) & 0x00F000F0;  high = high ^ t ^ (t << 4);
  t = (high ^ (high >> 8)) & 0x0000FF00;  high = high ^ t ^ (t << 8);

  /* unshuffle low */
  t = (low ^ (low >> 1)) & 0x22222222;  low = low ^ t ^ (t << 1);
  t = (low ^ (low >> 2)) & 0x0C0C0C0C;  low = low ^ t ^ (t << 2);
  t = (low ^ (low >> 4)) & 0x00F000F0;  low = low ^ t ^ (t << 4);
  t = (low ^ (low >> 8)) & 0x0000FF00;  low = low ^ t ^ (t << 8);

  *highbits = (high & 0xFFFF0000) | (low >> 16);
  *lowbits = (high << 16) | (low & 0x0000FFFF);

  return;
}
#endif


void
Compress_unshuffle (FILE *out, FILE *in) {
  Genomecomp_T high, low, flags;
  Genomecomp_T highbits, lowbits;

  while (FREAD_UINT(&high,in) > 0 &&
	 FREAD_UINT(&low,in) > 0 &&
	 FREAD_UINT(&flags,in) > 0) {
    nt_unshuffle(&highbits,&lowbits,high,low);
    FWRITE_UINT(highbits,out);
    FWRITE_UINT(lowbits,out);
    FWRITE_UINT(flags,out);
  }

  return;
}


#define BLOCK_SIZE 4		/* 128 / 32 = 4, adequate for a 128-bit register */

void
Compress_unshuffle_bits (FILE *high_out, FILE *low_out, FILE *flags_out, FILE *in) {
  UINT4 high, low;
  UINT4 highbits[BLOCK_SIZE], lowbits[BLOCK_SIZE], flags[BLOCK_SIZE];
  int i;

  while (!feof(in)) {
    for (i = 0; i < BLOCK_SIZE; i++) {
      if (FREAD_UINT(&high,in) == 0 ||
	  FREAD_UINT(&low,in) == 0 ||
	  FREAD_UINT(&(flags[i]),in) == 0) {
	/* End of genome */
	highbits[i] = -1U;
	lowbits[i] = -1U;
	flags[i] = -1U;
      } else {
	nt_unshuffle(&(highbits[i]),&(lowbits[i]),high,low);
      }
    }
      
    FWRITE_UINTS(highbits,i,high_out);
    FWRITE_UINTS(lowbits,i,low_out);
    FWRITE_UINTS(flags,i,flags_out);
  }

  return;
}


/* Needed for user-provided segment in GMAP */
Genomecomp_T *
Compress_create_blocks_comp (char *genomicseg, Univcoord_T genomelength) {
  Genomecomp_T *genomecomp;
  size_t nuint4;

  nuint4 = ((genomelength + 31)/32U)*3;
  genomecomp = (Genomecomp_T *) CALLOC(nuint4+4,sizeof(Genomecomp_T));
  /* Add 4 because Oligoindex_hr procedures point to nextlow as ptr+4 */

  /* Creates X's at end */
  genomecomp[nuint4-3] = 0xFFFFFFFF;
  genomecomp[nuint4-2] = 0xFFFFFFFF;
  genomecomp[nuint4-1] = 0xFFFFFFFF;

  /* Plus extra 4 */
  genomecomp[nuint4]   = 0xFFFFFFFF;
  genomecomp[nuint4+1] = 0xFFFFFFFF;
  genomecomp[nuint4+2] = 0xFFFFFFFF;
  genomecomp[nuint4+3] = 0xFFFFFFFF;

  Compress_update_memory(/*nbadchars*/0,genomecomp,genomicseg,/*currposition*/0,genomelength);

  return genomecomp;
}


/* Needed for user-provided segment in GMAP */
Genomecomp_T *
Compress_create_blocks_bits (Genomecomp_T *genomecomp, Univcoord_T genomelength) {
  Genomecomp_T *genomebits, highbits, lowbits, high, low, flags;
  size_t nuint4, ptr;

  nuint4 = ((genomelength + 31)/32U)*3;
  genomebits = (Genomecomp_T *) CALLOC(nuint4+4,sizeof(Genomecomp_T));

  for (ptr = 0; ptr < nuint4; ptr += 3) {
    high = genomecomp[ptr];
    low = genomecomp[ptr+1];
    flags = genomecomp[ptr+2];
    
    nt_unshuffle(&highbits,&lowbits,high,low);
    genomebits[ptr] = highbits;
    genomebits[ptr+1] = lowbits;
    genomebits[ptr+2] = flags;
  }

  return genomebits;
}


#define GENOMECOMP_LENGTH 32	/* In bp */

void
Compress_cat (FILE *out, char **files, Univcoord_T *genomelengths, int nfiles) {
  FILE *fp;
  char *file;
  int current_pos = 0, bufferlen, shift;
  UINT4 Buffer[3];
#ifdef HAVE_64_BIT
  UINT8 current_highlow = 0, buffer_highlow = 0;
#else
  UINT4 buffer_high = 0, buffer_low = 0;
#endif
  UINT4 current_high = 0, current_low = 0, current_flags = 0, buffer_flags = 0, high, low;
#ifdef DEBUG1
  int k;
#endif

#ifdef DEBUG1
  printf("nfiles: %d\n",nfiles);
  printf("* Genomelengths");
  for (k = 0; k < nfiles; k++) {
    printf(" %llu",genomelengths[k]);
  }
  printf("\n");
#endif


  if ((fp = fopen((file = *files++),"rb")) == NULL) {
    fprintf(stderr,"Genomecomp file %s is not valid\n",file);
    exit(9);
  } else {
    debug1(printf("*** Opening genomecomp file %s\n",file));
  }
    
  while (fp != NULL) {
#ifdef DEBUG1
    printf("* Genomelengths");
    for (k = 0; k < nfiles; k++) {
      printf(" %llu",genomelengths[k]);
    }
    printf("\n");
#endif

    FREAD_UINTS(Buffer,3,fp);
    debug1(printf("Reading buffer: %08x %08x %08x\n",Buffer[0],Buffer[1],Buffer[2]));
#ifdef HAVE_64_BIT
    buffer_highlow = (((UINT8) Buffer[0]) << 32) | ((UINT8) Buffer[1]);
#else
    buffer_high = Buffer[0];
    buffer_low = Buffer[1];
#endif
    buffer_flags = Buffer[2];

    if (*genomelengths <= GENOMECOMP_LENGTH) {
      bufferlen = *genomelengths;
      *genomelengths = 0;

#ifdef HAVE_64_BIT
      buffer_highlow &= 0xFFFFFFFFFFFFFFFF >> (64 - 2*bufferlen); /* mask trailing FFFFFFFF */
#else
      if (bufferlen <= 16) {
	buffer_high = 0;
	buffer_low &= 0xFFFFFFFF >> (32 - 2*bufferlen);
      } else {
	buffer_high &= 0xFFFFFFFF >> (64 - 2*bufferlen);
      }
#endif
      buffer_flags &= 0xFFFFFFFF >> (32 - bufferlen); /* mask trailing FFFFFFFF */

      fclose(fp);
      if (--nfiles == 0) {
	fp = (FILE *) NULL;
      } else if ((fp = fopen((file = *files++),"rb")) == NULL) {
	fprintf(stderr,"Regiondb file %s is not valid\n",file);
	exit(9);
      } else {
	debug1(printf("*** Opening regiondb file %s\n",file));
	genomelengths++;
      }

    } else {
      bufferlen = GENOMECOMP_LENGTH;
      *genomelengths -= GENOMECOMP_LENGTH;
    }
#ifdef HAVE_64_BIT
    debug1(printf("Buffer with bufferlen %d: %016llx %08x\n",bufferlen,buffer_highlow,buffer_flags));
    debug1(printf("Current with current_pos %d: %016llx %08x\n",current_pos,current_highlow,current_flags));
#else
    debug1(printf("Buffer with bufferlen %d: %08x %08x %08x\n",bufferlen,buffer_high,buffer_low,buffer_flags));
    debug1(printf("Current with current_pos %d: %08x %08x %08x\n",current_pos,current_high,current_low,current_flags));
#endif

    if (current_pos + bufferlen < GENOMECOMP_LENGTH) {
      debug1(printf("Case 1\n"));
      debug1(printf("Appending buffer with bufferlen %d\n",bufferlen));
#ifdef HAVE_64_BIT
      current_highlow |= buffer_highlow << 2*current_pos; /* append */
#else
      /* append */
      if (current_pos == 0) {
        current_low = buffer_low;
        current_high = buffer_high;
      } else if (current_pos < 16) {
	current_low |= buffer_low << 2*current_pos;
	current_high = (buffer_low >> (32 - 2*current_pos)) | (buffer_high << 2*current_pos);
      } else if (current_pos == 16) {
	current_low = 0;
	current_high = buffer_low;
      } else {
	current_high |= buffer_low << (2*current_pos - 32);
      }
#endif
      current_flags |= buffer_flags << current_pos; /* append */
      current_pos += bufferlen;

    } else if (current_pos + bufferlen == GENOMECOMP_LENGTH) {
      debug1(printf("Case 2\n"));
      debug1(printf("Appending buffer with bufferlen %d\n",bufferlen));
#ifdef HAVE_64_BIT
      current_highlow |= buffer_highlow << 2*current_pos; /* append */
#else
      /* append */
      if (current_pos == 0) {
        current_low = buffer_low;
        current_high = buffer_high;
      } else if (current_pos < 16) {
	current_low |= buffer_low << 2*current_pos;
	current_high = (buffer_low >> (32 - 2*current_pos)) | (buffer_high << 2*current_pos);
      } else if (current_pos == 16) {
	current_low = 0;
	current_high = buffer_low;
      } else {
	current_high |= buffer_low << (2*current_pos - 32);
      }
#endif
      current_flags |= buffer_flags << current_pos; /* append */

#ifdef HAVE_64_BIT
      current_high = (UINT4) (current_highlow >> 32); current_low = (UINT4) (current_highlow & 0xFFFFFFFF); /* write */
#endif
      FWRITE_UINT(current_high,out); FWRITE_UINT(current_low,out); FWRITE_UINT(current_flags,out); /* write */
      debug1(printf("Writing current: %08x %08x %08x\n",current_high,current_low,current_flags));

#ifdef HAVE_64_BIT
      current_highlow = 0; current_flags = 0; /* clear */
#else
      current_high = current_low = current_flags = 0; /* clear */
#endif

      current_pos = 0;

    } else {
      debug1(printf("Case 3\n"));
      debug1(printf("Appending buffer with truncated bufferlen %d\n",GENOMECOMP_LENGTH-current_pos));
#ifdef HAVE_64_BIT
      current_highlow |= buffer_highlow << 2*current_pos; /* append */
#else
      /* append */
      if (current_pos == 0) {
        current_low = buffer_low;
        current_high = buffer_high;
      } else if (current_pos < 16) {
	current_low |= buffer_low << 2*current_pos;
	current_high = (buffer_low >> (32 - 2*current_pos)) | (buffer_high << 2*current_pos);
      } else if (current_pos == 16) {
	current_low = 0;
	current_high = buffer_low;
      } else {
	current_high |= buffer_low << (2*current_pos - 32);
      }
#endif
      current_flags |= buffer_flags << current_pos; /* append */

#ifdef HAVE_64_BIT
      current_high = (UINT4) (current_highlow >> 32); current_low = (UINT4) (current_highlow & 0xFFFFFFFF); /* write */
#endif
      FWRITE_UINT(current_high,out); FWRITE_UINT(current_low,out); FWRITE_UINT(current_flags,out); /* write */
      debug1(printf("Writing current: %08x %08x %08x\n",current_high,current_low,current_flags));

      shift = GENOMECOMP_LENGTH - current_pos;
      debug1(printf("Saving remainder with shift %d\n",shift));
#ifdef HAVE_64_BIT
      debug1(printf("Buffer: %016llx %08x\n",buffer_highlow,buffer_flags));
      current_highlow = buffer_highlow >> 2*shift; /* remainder */
#else
      debug1(printf("Buffer: %08x %08x %08x\n",buffer_high,buffer_low,buffer_flags));
      /* remainder */
      if (shift == 0) {
	current_high = buffer_high;
	current_low = buffer_low;
      } else if (shift < 16) {
	current_high = buffer_high >> 2*shift;
	current_low = (buffer_low >> 2*shift) | (buffer_high << (32 - 2*shift));
      } else if (shift == 16) {
	current_high = 0;
	current_low = buffer_high;
      } else {
	current_high = 0;
	current_low = buffer_high >> (2*shift - 32);
      }
#endif
      current_flags = buffer_flags >> shift; /* remainder */
#ifdef HAVE_64_BIT
      debug1(printf("Remainder: %016llx %08x\n",current_highlow,current_flags));
#else
      debug1(printf("Remainder: %08x %08x %08x\n",current_high,current_low,current_flags));
#endif

      current_pos = current_pos + bufferlen - GENOMECOMP_LENGTH;
    }
#ifdef HAVE_64_BIT
    debug1(printf("Now current with current_pos %d: %016llx %08x\n",current_pos,current_highlow,current_flags));
#else
    debug1(printf("Now current with current_pos %d: %08x %08x %08x\n",current_pos,current_high,current_low,current_flags));
#endif
  }

  if (current_pos > 0) {
    /* Introduce trailing FFFFFFFF.  Based on current_pos. */
    debug1(printf("Final block.  current_pos %d\n",current_pos));
#ifdef HAVE_64_BIT
    current_highlow |= 0xFFFFFFFFFFFFFFFF << 2*current_pos; /* add trailing FFFFFFFF */
    current_high = (UINT4) (current_highlow >> 32); current_low = (UINT4) (current_highlow & 0xFFFFFFFF); /* write */
#else
    /* add trailing FFFFFFFF */
    if (current_pos < 16) {
      current_high = 0xFFFFFFFF;
      current_low |= 0xFFFFFFFF << 2*current_pos;
    } else {
      current_high |= 0xFFFFFFFF << (2*current_pos - 32);
    }
#endif
    current_flags |= 0xFFFFFFFF << current_pos; /* add trailing FFFFFFFF */
    FWRITE_UINT(current_high,out); FWRITE_UINT(current_low,out); FWRITE_UINT(current_flags,out); /* write */
    debug1(printf("Writing current: %08x %08x %08x\n",current_high,current_low,current_flags));
  }

  /* Write an extra two UINT4 */
  high = low = 0xFFFFFFFF;
  FWRITE_UINT(high,out);
  FWRITE_UINT(low,out);

  return;
}

