static char rcsid[] = "$Id: 68aaf96c89bab30042fbecce02123d2f6a9bdb5b $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "genomebits.h"

#include <stdio.h>
#include <string.h>
#include <sys/mman.h>		/* For mmap */

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "genomicpos.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Genomebits_T

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


void
Genomebits_print (T this, Univcoord_T startpos, Univcoord_T endpos) {
  Genomecomp_T *high_ptr, *low_ptr, *flags_ptr, *end_ptr;
  Univcoord_T startblocki_32, endblocki_32;
  int startdiscard32, enddiscard32;
  Genomecomp_T high, low, flags;
  int i;
  bool firstp = true;

  startblocki_32 = startpos/32;
  endblocki_32 = endpos/32;
  startdiscard32 = startpos % 32;
  enddiscard32 = endpos % 32;

  high_ptr = &(this->high_blocks[startblocki_32]);
  low_ptr = &(this->low_blocks[startblocki_32]);
  flags_ptr = &(this->flags_blocks[startblocki_32]);
  end_ptr = &(this->high_blocks[endblocki_32]);

  while (high_ptr <= end_ptr) {
    if (firstp == true) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 -- spacing*/
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < startdiscard32; i++) {
	printf("v");
      }
      for ( ; i < 32; i++) {
	printf(" ");
      }
      printf("\n");
      firstp = false;
    }

    high = *high_ptr; low = *low_ptr; flags = *flags_ptr;
    printf("high: %08X  low: %08X  flags: %08X\t",high,low,flags);
    write_chars(high,low,flags);
    printf("\n");

    if (high_ptr == end_ptr) {
      /*      high: 9F61B62A  low: 6D68A157  flags: 00000000 -- spacing */
      printf("                                              \t");
      /* printf("%u\t",startblocki); */
      for (i = 0; i < enddiscard32; i++) {
	printf(" ");
      }
      for ( ; i < 32; i++) {
	printf("^");
      }
      printf("\n");
    }

    high_ptr += 1; low_ptr += 1; flags_ptr += 1;
  }

  printf("\n");
  return;
}



void
Genomebits_free (T *old) {
  if (*old) {
    if ((*old)->access == ALLOCATED_PRIVATE) {
      FREE_KEEP((*old)->flags_blocks);
      FREE_KEEP((*old)->low_blocks);
      FREE_KEEP((*old)->high_blocks);

    } else if ((*old)->access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->flags_blocks,(*old)->flags_blocks_shmid,(*old)->flags_blocks_key);
      Access_deallocate((*old)->low_blocks,(*old)->low_blocks_shmid,(*old)->low_blocks_key);
      Access_deallocate((*old)->high_blocks,(*old)->high_blocks_shmid,(*old)->high_blocks_key);

#ifdef HAVE_MMAP
    } else if ((*old)->access == MMAPPED) {
      munmap((void *) (*old)->flags_blocks,(*old)->flags_len);
      munmap((void *) (*old)->low_blocks,(*old)->low_len);
      munmap((void *) (*old)->high_blocks,(*old)->high_len);
      close((*old)->flags_fd);
      close((*old)->low_fd);
      close((*old)->high_fd);
#endif
    }

    FREE(*old);
  }
  return;
}


T
Genomebits_new (char *genomesubdir, char *fileroot, char *alt_suffix,
		Access_mode_T access, bool sharedp, bool revcompp) {
  T new = (T) MALLOC(sizeof(*new));
  char *high_filename, *low_filename, *flags_filename;
  double seconds0, seconds1, seconds2;
  int npages0, npages1, npages2;
#ifndef UTILITYP
  char *comma;
#endif

#if 0
  new->genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias*/true);
#endif

  if (alt_suffix != NULL) {
    if (revcompp == true) {
      high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".genomerevh.")+strlen(alt_suffix)+1,sizeof(char));
      low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				     strlen(".genomerevl.")+strlen(alt_suffix)+1,sizeof(char));
      flags_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".genomerevf.")+strlen(alt_suffix)+1,sizeof(char));
      sprintf(high_filename,"%s/%s.genomerevh.%s",genomesubdir,fileroot,alt_suffix);
      sprintf(low_filename,"%s/%s.genomerevl.%s",genomesubdir,fileroot,alt_suffix);
      sprintf(flags_filename,"%s/%s.genomerevf.%s",genomesubdir,fileroot,alt_suffix);
    } else {
      high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".genomefwdh.")+strlen(alt_suffix)+1,sizeof(char));
      low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				     strlen(".genomefwdl.")+strlen(alt_suffix)+1,sizeof(char));
      flags_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".genomefwdf.")+strlen(alt_suffix)+1,sizeof(char));
      sprintf(high_filename,"%s/%s.genomefwdh.%s",genomesubdir,fileroot,alt_suffix);
      sprintf(low_filename,"%s/%s.genomefwdl.%s",genomesubdir,fileroot,alt_suffix);
      sprintf(flags_filename,"%s/%s.genomefwdf.%s",genomesubdir,fileroot,alt_suffix);
    }
  } else {
    if (revcompp == true) {
      high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".genomerevh")+1,sizeof(char));
      low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				     strlen(".genomerevl")+1,sizeof(char));
      flags_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".genomerevf")+1,sizeof(char));
      sprintf(high_filename,"%s/%s.genomerevh",genomesubdir,fileroot);
      sprintf(low_filename,"%s/%s.genomerevl",genomesubdir,fileroot);
      sprintf(flags_filename,"%s/%s.genomerevf",genomesubdir,fileroot);
    } else {
      high_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				      strlen(".genomefwdh")+1,sizeof(char));
      low_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				     strlen(".genomefwdl")+1,sizeof(char));
      flags_filename = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(fileroot)+
				       strlen(".genomefwdf")+1,sizeof(char));
      sprintf(high_filename,"%s/%s.genomefwdh",genomesubdir,fileroot);
      sprintf(low_filename,"%s/%s.genomefwdl",genomesubdir,fileroot);
      sprintf(flags_filename,"%s/%s.genomefwdf",genomesubdir,fileroot);
    }
  }

  if (Access_file_exists_p(high_filename) == false) {
    fprintf(stderr,"Unable to detect new version of genome index: %s file not available.  This version of GSNAP is not backwards compatible.\n",high_filename);
    FREE(flags_filename);
    FREE(low_filename);
    FREE(high_filename);
    FREE(new);
    exit(9);
    return (T) NULL;
  }

  if (access == USE_ALLOCATE) {
#ifndef UTILITYP
    fprintf(stderr,"Allocating memory for compressed genome (bits)...");
#endif
    if (sharedp == true) {
      new->high_blocks = (Genomecomp_T *) Access_allocate_shared(&new->access,&new->high_blocks_shmid,&new->high_blocks_key,
								 &new->high_fd,&new->high_len,&seconds0,high_filename,sizeof(Genomecomp_T));
      new->low_blocks = (Genomecomp_T *) Access_allocate_shared(&new->access,&new->low_blocks_shmid,&new->low_blocks_key,
								&new->low_fd,&new->low_len,&seconds1,low_filename,sizeof(Genomecomp_T));
      new->flags_blocks = (Genomecomp_T *) Access_allocate_shared(&new->access,&new->flags_blocks_shmid,&new->flags_blocks_key,
								&new->flags_fd,&new->flags_len,&seconds2,flags_filename,sizeof(Genomecomp_T));
    } else {
      new->high_blocks = (Genomecomp_T *) Access_allocate_private(&new->access,&new->high_len,&seconds0,high_filename,sizeof(Genomecomp_T));
      new->low_blocks = (Genomecomp_T *) Access_allocate_private(&new->access,&new->low_len,&seconds1,low_filename,sizeof(Genomecomp_T));
      new->flags_blocks = (Genomecomp_T *) Access_allocate_private(&new->access,&new->flags_len,&seconds2,flags_filename,sizeof(Genomecomp_T));
    }
    if (new->high_blocks == NULL || new->low_blocks == NULL || new->flags_blocks == NULL) {
#ifndef UTILITYP
      fprintf(stderr,"insufficient memory (need to use a lower batch mode (-B))\n");
#endif
      exit(9);
    } else {
#ifndef UTILITYP
      comma = Genomicpos_commafmt(new->high_len + new->low_len + new->flags_len);
      fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds0 + seconds1 + seconds2);
      FREE(comma);
#endif
    }
    
#ifdef HAVE_MMAP
  } else if (access == USE_MMAP_PRELOAD) {
#ifndef UTILITYP
    fprintf(stderr,"Pre-loading compressed genome (bits)...");
#endif
    new->high_blocks = (Genomecomp_T *) Access_mmap_and_preload(&new->high_fd,&new->high_len,&npages0,&seconds0,
								high_filename,sizeof(Genomecomp_T));
    new->low_blocks = (Genomecomp_T *) Access_mmap_and_preload(&new->low_fd,&new->low_len,&npages1,&seconds1,
								low_filename,sizeof(Genomecomp_T));
    new->flags_blocks = (Genomecomp_T *) Access_mmap_and_preload(&new->flags_fd,&new->flags_len,&npages2,&seconds2,
								flags_filename,sizeof(Genomecomp_T));
    if (new->high_blocks == NULL || new->low_blocks == NULL || new->flags_blocks == NULL) {
#ifndef UTILITYP
      fprintf(stderr,"insufficient memory\n");
#endif
      exit(9);
    } else {
#ifndef UTILITYP
      comma = Genomicpos_commafmt(new->high_len + new->low_len + new->flags_len);
      fprintf(stderr,"done (%s bytes, %d pages, %.2f sec)\n",comma,npages0 + npages1 + npages2,seconds0 + seconds1 + seconds2);
      FREE(comma);
#endif
      new->access = MMAPPED;
    }
    
  } else if (access == USE_MMAP_ONLY) {
    new->high_blocks = (Genomecomp_T *) Access_mmap(&new->high_fd,&new->high_len,&seconds0,high_filename,/*randomp*/false);
    new->low_blocks = (Genomecomp_T *) Access_mmap(&new->low_fd,&new->low_len,&seconds1,low_filename,/*randomp*/false);
    new->flags_blocks = (Genomecomp_T *) Access_mmap(&new->flags_fd,&new->flags_len,&seconds2,flags_filename,/*randomp*/false);
    if (new->high_blocks == NULL || new->low_blocks == NULL || new->flags_blocks == NULL) {
      fprintf(stderr,"Insufficient memory for genome mmap\n");
      exit(9);
    } else {
      new->access = MMAPPED;
    }
#endif
    
  }
  
  FREE(flags_filename);
  FREE(low_filename);
  FREE(high_filename);

  return new;
}

static char CHARTABLE[] = "ACGTXXXX";

char
Genomebits_get_char (T this, Univcoord_T pos) {
  Univcoord_T blocki_32;
  int relpos;
  Genomecomp_T high, low, flags;
  int idx;

  blocki_32 = pos/32;
  relpos = pos % 32;

  debug(Genomebits_print(this,pos,pos+1));

  high = this->high_blocks[blocki_32]; low = this->low_blocks[blocki_32]; flags = this->flags_blocks[blocki_32];
#ifndef HAVE_SSE2
  idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#elif defined(SSE2_SLLI_CONST_IMM8)
  idx = ((flags >> relpos) & 0x1) << 2 |((high >> relpos) & 0x1) << 1 | ((low >> relpos) & 0x1);
#else
  idx = _mm_movemask_ps((__m128) (_mm_slli_epi32(_mm_set_epi32(0,flags,high,low),31 - relpos)));
#endif
  return CHARTABLE[idx];
}


/* Logic for genomebits processing:

1. Processing rightward from pos5 to pos3:

Case A. endblocki == startblocki, word [0]:
        <<<<<<<0
        FFFFFFFF
   pos3 pos5/ptr
        read 32@ptr
        startdiscard,enddiscard   

Case B. endblocki == startblocki, word [0..1]:
<<<<<<<1<<<<<<<0
FFFFFFFFFFFFFFFF
  pos3  pos5/ptr
read 64@ptr
startdiscard,enddiscard + 32

Case D. ptr == endptr, words [0..1],[2]:
C/D loop                   D terminate
<<<<<<<1<<<<<<<0           <<<<<<<2
FFFFFFFFFFFFFFFF           FFFFFFFF
        pos5/ptr           pos3/end
read 64@ptr                read 32@ptr
startdiscard               enddiscard

Case C. ptr + 1 == endptr [0..1],[2..3]:
C/D loop           C terminate
<<<<<<<1<<<<<<<0   <<<<<<<3<<<<<<<2
FFFFFFFFFFFFFFFF   FFFFFFFFFFFFFFFF
        pos5/ptr   pos3/end     ptr
read 64@ptr        read 64@ptr
startdiscard       enddiscard + 32


2. Processing leftward from pos3 to pos5:

Case A. startblocki == endblocki, word [3]:
                   <<<<<<<3
                   FFFFFFFF
                   pos3/ptr pos5
read 32@ptr
enddiscard,startdiscard

Case B. startblocki + 1 == endblocki, word [3..2]:
                   <<<<<<<3<<<<<<<2
                   FFFFFFFFFFFFFFFF
                   pos3/ptr  pos5
                   read 64@ptr - 1
                   enddiscard + 32,startdiscard

Case D. ptr == startptr, words [3..2],[1]:
D terminate        C/D loop
<<<<<<<1           <<<<<<<3<<<<<<<2
FFFFFFFF           FFFFFFFFFFFFFFFF
pos5/ptr           pos3/ptr
read 32@ptr        read 64@ptr - 1
                   enddiscard + 32

Case C. ptr == startptr + 1, words [3..2],[1..0]:
C terminate        C/D loop
<<<<<<<1<<<<<<<0   <<<<<<<3<<<<<<<2
FFFFFFFFFFFFFFFF   FFFFFFFFFFFFFFFF
        pos5/ptr   pos3/ptr
read 64@ptr-1      read 64@ptr - 1
startdiscard       enddiscard + 32

*/
