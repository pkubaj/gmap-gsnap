static char rcsid[] = "$Id: d2e7d015de963a8324654340fc1a4c7c388395a2 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "transcriptome.h"

#include <stdio.h>
#include <stddef.h>
#include <string.h>		/* Needed for strlen and strcpy */

#include "mem.h"
#include "access.h"
#include "list.h"
#include "genomicpos.h"
#include "bitpack64-readtwo.h"


#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif



#define T Transcriptome_T
struct T {

  int *dbindex;

  Chrnum_T *chrnums;

  Access_T offsetsmeta_access;
  int offsetsmeta_shmid;
  key_t offsetsmeta_key;
  int offsetsmeta_fd;
  size_t offsetsmeta_len;
  UINT4 *offsetsmeta;

  Access_T offsetsstrm_access;
  int offsetsstrm_shmid;
  key_t offsetsstrm_key;
  int offsetsstrm_fd;
  size_t offsetsstrm_len;
  UINT4 *offsetsstrm;

  Access_T exoninfo_access;
  int exoninfo_shmid;
  key_t exoninfo_key;
  int exoninfo_fd;
  size_t exoninfo_len;

  /* exoninfo has exonbounds (int) and exonstarts (unsigned int) interleaved */
  unsigned int *exoninfo;
};


#define BUFFERSIZE 1024000

void
Transcriptome_free (T *old) {
  if (*old) {

    if ((*old)->exoninfo_access == ALLOCATED_PRIVATE) {
      FREE_KEEP((*old)->exoninfo);
    } else if ((*old)->exoninfo_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->exoninfo,(*old)->exoninfo_shmid,(*old)->exoninfo_key);
    } else {
      abort();
    }

    if ((*old)->offsetsstrm_access == ALLOCATED_PRIVATE) {
      FREE_KEEP((*old)->offsetsstrm);
    } else if ((*old)->offsetsstrm_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsstrm,(*old)->offsetsstrm_shmid,(*old)->offsetsstrm_key);
    } else {
      abort();
    }
      
    if ((*old)->offsetsmeta_access == ALLOCATED_PRIVATE) {
      FREE_KEEP((*old)->offsetsmeta);
    } else if ((*old)->offsetsmeta_access == ALLOCATED_SHARED) {
      Access_deallocate((*old)->offsetsmeta,(*old)->offsetsmeta_shmid,(*old)->offsetsmeta_key);
    } else {
      abort();
    }

    FREE((*old)->dbindex);
    FREE((*old)->chrnums);
    FREE(*old);
  }
  return;
}


T
Transcriptome_new (char *genomesubdir, char *genome_fileroot, char *transcriptome_fileroot,
		   bool sharedp) {
  T new = (T) MALLOC(sizeof(*new));

  FILE *fp;
  char *info_root, *filename;

  char *comma;
  double seconds;
  int ntranscripts, nalignments, trnum, map_index;
  int divint;


  info_root = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(genome_fileroot)+strlen(".transcripts/")+strlen(transcriptome_fileroot)+1,
			      sizeof(char));
  sprintf(info_root,"%s/%s.transcripts/%s",genomesubdir,genome_fileroot,transcriptome_fileroot);

  /* Read dbindex file */
  filename = (char *) CALLOC(strlen(info_root)+strlen(".dbindex")+1,sizeof(char));
  sprintf(filename,"%s.dbindex",info_root);
  nalignments = Access_filesize(filename)/sizeof(int) - 1;

  if ((fp = fopen(filename,"rb")) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",filename);
    exit(9);
  }
  FREE(filename);

  new->dbindex = (int *) MALLOC((nalignments+1)*sizeof(int));
  for (map_index = 0; map_index <= nalignments; map_index++) {
    FREAD_INT(&trnum,fp);
    new->dbindex[map_index] = trnum;
  }
  fclose(fp);


  /* Read chrnums file */
  filename = (char *) CALLOC(strlen(info_root)+strlen(".chrnums")+1,sizeof(char));
  sprintf(filename,"%s.chrnums",info_root);
  ntranscripts = Access_filesize(filename)/sizeof(int);

  if ((fp = fopen(filename,"rb")) == NULL) {
    fprintf(stderr,"Cannot open file %s\n",filename);
    exit(9);
  }
  FREE(filename);

  new->chrnums = (Chrnum_T *) MALLOC((ntranscripts+1)*sizeof(Chrnum_T));
  new->chrnums[0] = 0;
  for (trnum = 1; trnum <= ntranscripts; trnum++) {
    FREAD_INT(&divint,fp);
    new->chrnums[trnum] = divint;
  }
  fclose(fp);


  /* Read bitpack info */
  filename = (char *) CALLOC(strlen(info_root)+strlen(".offsets64meta")+1,sizeof(char));
  sprintf(filename,"%s.offsets64meta",info_root);
  if (sharedp == true) {
    new->offsetsmeta = (UINT4 *) Access_allocate_shared(&new->offsetsmeta_access,&new->offsetsmeta_shmid,&new->offsetsmeta_key,
							&new->offsetsmeta_fd,&new->offsetsmeta_len,&seconds,
							filename,sizeof(UINT4));
  } else {
    new->offsetsmeta = (UINT4 *) Access_allocate_private(&new->offsetsmeta_access,&new->offsetsmeta_len,&seconds,
							filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->offsetsmeta_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);


  filename = (char *) CALLOC(strlen(info_root)+strlen(".offsets64strm")+1,sizeof(char));
  sprintf(filename,"%s.offsets64strm",info_root);
  if (sharedp == true) {
    new->offsetsstrm = (UINT4 *) Access_allocate_shared(&new->offsetsstrm_access,&new->offsetsstrm_shmid,&new->offsetsstrm_key,
							&new->offsetsstrm_fd,&new->offsetsstrm_len,&seconds,
							filename,sizeof(UINT4));
  } else {
    new->offsetsstrm = (UINT4 *) Access_allocate_private(&new->offsetsstrm_access,&new->offsetsstrm_len,&seconds,
							 filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->offsetsstrm_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);



  /* Read exoninfo file */
  filename = (char *) CALLOC(strlen(info_root)+strlen(".exoninfo")+1,sizeof(char));
  sprintf(filename,"%s.exoninfo",info_root);
  if (sharedp == true) {
    new->exoninfo = (UINT4 *) Access_allocate_shared(&new->exoninfo_access,&new->exoninfo_shmid,&new->exoninfo_key,
						     &new->exoninfo_fd,&new->exoninfo_len,&seconds,
						     filename,sizeof(UINT4));
  } else {
    new->exoninfo = (UINT4 *) Access_allocate_private(&new->exoninfo_access,&new->exoninfo_len,&seconds,
						      filename,sizeof(UINT4));
  }
  comma = Genomicpos_commafmt(new->exoninfo_len);
  fprintf(stderr,"done (%s bytes, %.2f sec)\n",comma,seconds);
  FREE(comma);
  FREE(filename);
    
  FREE(info_root);

  return new;
}


Chrnum_T
Transcriptome_chrnum (int *transcript_genestrand, T this, int trnum) {
  Chrnum_T chrnum;

  if ((chrnum = this->chrnums[trnum]) == 0) {
    return 0;
  } else if (chrnum > 0) {
    *transcript_genestrand = +1;
    return chrnum;
  } else {
    *transcript_genestrand = -1;
    return -chrnum;
  }
}


int
Transcriptome_nexons (T this, int trnum) {
  UINT4 offset0, offset1;

  offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
  return (int) (offset1 - offset0);
}


int
Transcriptome_exons (int **exonbounds, Chrpos_T **exonstarts, T this, Trnum_T trnum) {
  int nexons;
  UINT4 offset0, offset1;

  offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
  nexons = (int) (offset1 - offset0);
  *exonbounds = &(((int *) this->exoninfo)[2*offset0]);
  *exonstarts = &(this->exoninfo[2*offset0 + nexons]);

#if 0
  printf("Transcript %d has %d exons:\n",trnum,nexons);
  for (int i = 0; i < nexons; i++) {
    printf("transcript pos %d chrpos %u\n",(*exonbounds)[i],(*exonstarts)[i]);
  }
  printf("\n");
#endif

  return nexons;
}


/* Assumes that we have already checked that this->chrnums != 0 */
Trnum_T
Transcriptome_trnum (int *nexons, int **exonbounds, Chrpos_T **exonstarts, T this, int map_index) {
  Trnum_T trnum;
  UINT4 offset0, offset1;

  if ((trnum = this->dbindex[map_index]) < 1) {
    /* Not found in transcriptomedb */
    return -1;

  } else {
    offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
    *nexons = (int) (offset1 - offset0);
    *exonbounds = &(((int *) this->exoninfo)[2*offset0]);
    *exonstarts = &(this->exoninfo[2*offset0 + (*nexons)]);

#if 0
    printf("Transcript %d has %d exons:\n",trnum,*nexons);
    for (int i = 0; i < *nexons; i++) {
      printf("transcript pos %d chrpos %u\n",(*exonbounds)[i],(*exonstarts)[i]);
    }
    printf("\n");
#endif

    return trnum;
  }
}


#if 0
bool
Transcriptome_genomic_bounded_p (int trnum, Chrnum_T chrbound, Chrpos_T lowbound, Chrpos_T highbound, T this) {
  UINT4 offset0, offset1;
  int nexons;
  int *exonbounds;
  Chrnum_T chrnum;
  Chrpos_T *exonstarts;
  Univcoord_T start, end;

  if ((chrnum = this->chrnums[trnum]) == 0) {
    return false;
  } else if (chrnum > 0) {
    if (chrnum != chrbound) {
      return false;
    }
  } else {
    if (-chrnum != chrbound) {
      return false;
    }
  }

  offset0 = Bitpack64_read_two(&offset1,(Oligospace_T) (trnum - 1),this->offsetsmeta,this->offsetsstrm);
  nexons = (int) (offset1 - offset0);

  exonbounds = &(((int *) this->exoninfo)[2*offset0]);
  exonstarts = &(this->exoninfo[2*offset0 + nexons]);

  start = exonstarts[0];
  /* ? Shouldn't this be: ... + exonbounds[nexons - 1] - exonbounds[nexons - 1] */
  end = exonstarts[nexons - 1] + exonbounds[nexons - 1];
  if (start < end) {
    if (end < lowbound) {
      return false;
    } else if (start > highbound) {
      return false;
    } else {
      return true;
    }

  } else {
    if (start < lowbound) {
      return false;
    } else if (end > highbound) {
      return false;
    } else {
      return true;
    }
  }
}
#endif

