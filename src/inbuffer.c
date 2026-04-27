static char rcsid[] = "$Id: efecfcbf5f5f3bc753a0d69f85834f01f118b034 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "inbuffer.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif

#include "mem.h"


/* GMAP uses Sequence_T, while GSNAP, GEXACT, and GFILTER use
   Shortread_T.  GSNAP and GFILTER can handle single-end or paired-end
   reads, while GEXACT handles only single-end reads */

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static bool single_cell_p;
static bool filter_if_both_p;


void
Inbuffer_setup (bool single_cell_p_in, bool filter_if_both_p_in) {
  single_cell_p = single_cell_p_in;
  filter_if_both_p = filter_if_both_p_in;
  return;
}


#define T Inbuffer_T

struct T {
  Outbuffer_T outbuffer;
  FILE *input;
#if defined(GSNAP) || defined(GFILTER)
  FILE *input2;
#endif

#ifdef HAVE_ZLIB
  gzFile gzipped;
  gzFile gzipped2;
#else
  void *gzipped;
  void *gzipped2;
#endif

#ifdef HAVE_BZLIB
  Bzip2_T bzipped;
  Bzip2_T bzipped2;
#else
  void *bzipped;
  void *bzipped2;
#endif

  char *read_files_command;
  char **files;
  int nfiles;
  int nextchar;

#if defined(GSNAP) || defined(GFILTER)
  bool interleavedp;
#endif

#if defined(HAVE_PTHREAD)
  pthread_mutex_t lock;
#endif

  Request_T *buffer;
  unsigned int nspaces;
  int ptr;
  int nleft;
  unsigned int inputid;
  unsigned int requestid;

#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
  /* GMAP */
  bool user_pairalign_p;
  List_T user_genomes;
  List_T overflow_head;
  List_T overflow_tail;
#endif

  unsigned int part_modulus;
  unsigned int part_interval;
};


#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
T
Inbuffer_cmdline (char *queryseq, int querylength, Genome_T genome, Genome_T genomealt) {
  T new = (T) MALLOC(sizeof(*new));

  new->input = (FILE *) NULL;

  new->files = (char **) NULL;
  new->nfiles = 0;
  new->nextchar = '\0';

  new->buffer = (Request_T *) CALLOC(1,sizeof(Request_T));
  new->nspaces = 0;		/* Set to zero so we don't try to read from a file */

  new->ptr = 0;
  new->nleft = 1;
  new->inputid = 0;
  new->requestid = 0;

  new->user_pairalign_p = false;
  new->user_genomes = (List_T) NULL;

  new->part_modulus = 0;
  new->part_interval = 1;

  /* global_genome and global_genomealt were set by Inbuffer_setup */
  new->buffer[0] = Request_new(new->requestid++,genome,genomealt,
			       Sequence_genomic_new(queryseq,querylength,/*copyp*/true),
			       /*free_genome_p*/true);
  
#if defined(HAVE_PTHREAD)
  pthread_mutex_init(&new->lock,NULL);
#endif

  return new;
}
#endif


T
Inbuffer_new (int nextchar, FILE *input,

#if defined(GEXACT)
#ifdef HAVE_ZLIB
	      gzFile gzipped,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped,
#endif

#elif defined(GSNAP) || defined(GFILTER)
	      FILE *input2,
#ifdef HAVE_ZLIB
	      gzFile gzipped, gzFile gzipped2,
#endif
#ifdef HAVE_BZLIB
	      Bzip2_T bzipped, Bzip2_T bzipped2,
#endif
	      bool interleavedp,
#endif

	      char *read_files_command, char **files, int nfiles, unsigned int nspaces,
#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
	      bool user_pairalign_p, List_T user_genomes,
#endif
	      unsigned int part_modulus, unsigned int part_interval) {

  T new = (T) MALLOC(sizeof(*new));

  new->input = input;
#if defined(GEXACT)
#ifdef HAVE_ZLIB
  new->gzipped = gzipped;
#else
  new->gzipped = (void *) NULL;
#endif
#ifdef HAVE_BZLIB
  new->bzipped = bzipped;
#else
  new->bzipped = (void *) NULL;
#endif
  
#elif defined(GSNAP) || defined(GFILTER)
  new->input2 = input2;
#ifdef HAVE_ZLIB
  new->gzipped = gzipped;
  new->gzipped2 = gzipped2;
#else
  new->gzipped = (void *) NULL;
  new->gzipped2 = (void *) NULL;
#endif
#ifdef HAVE_BZLIB
  new->bzipped = bzipped;
  new->bzipped2 = bzipped2;
#else
  new->bzipped = (void *) NULL;
  new->bzipped2 = (void *) NULL;
#endif
  new->interleavedp = interleavedp;
#endif

  new->read_files_command = read_files_command;
  new->files = files;
  new->nfiles = nfiles;
  new->nextchar = nextchar;

#if defined(HAVE_PTHREAD)
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->buffer = (Request_T *) CALLOC(nspaces,sizeof(Request_T));
  new->nspaces = nspaces;
  new->ptr = 0;
  new->nleft = 0;
  new->inputid = 0;
  new->requestid = 0;

#if !defined(GFILTER) && !defined(GEXACT) && !defined(GSNAP)
  new->user_pairalign_p = user_pairalign_p;
  new->user_genomes = user_genomes;
  new->overflow_head = (List_T) NULL;
  new->overflow_tail = (List_T) NULL;
#endif

  new->part_modulus = part_modulus;
  new->part_interval = part_interval;

  return new;
}

void
Inbuffer_set_outbuffer (T this, Outbuffer_T outbuffer) {
  this->outbuffer = outbuffer;
  return;
}

void
Inbuffer_free (T *old) {
  if (*old) {
    /* No need to close input, since done by Shortread and Sequence read procedures */
    FREE((*old)->buffer);
    
#if defined(HAVE_PTHREAD)
    pthread_mutex_destroy(&(*old)->lock);
#endif

    FREE(*old);
  }
  return;
}


#if defined(GEXACT)

/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  Shortread_T queryseq1, queryseq2;
  bool skipp;
  int nchars1 = 0, nchars2 = 0;		/* Returned only because MPI master needs it.  Doesn't need to be saved as a field in Inbuffer_T. */

  FILE *input2 = NULL;
#ifdef HAVE_ZLIB
  gzFile gzipped2 = NULL;
#endif
#ifdef HAVE_BZLIB
  Bzip2_T bzipped2 = NULL;
#endif
  
  /* fprintf(stderr,"Entered fill_buffer\n"); */
  while (nread < this->nspaces &&
	 (queryseq1 = Shortread_read(&this->nextchar,&nchars1,&nchars2,&queryseq2,
				     &this->input,&input2,
#ifdef HAVE_ZLIB
				     &this->gzipped,&gzipped2,
#endif
#ifdef HAVE_BZLIB
				     &this->bzipped,&bzipped2,
#endif
				     /*interleavedp*/false,
				     this->read_files_command,&this->files,&this->nfiles,single_cell_p,
				     skipp = (this->inputid % this->part_interval != this->part_modulus))) != NULL) {
    if (skipp) {
#if 0
      /* Shortread procedures won't allocate in this situation */
      Shortread_free(&queryseq1);
#endif
      
    } else {
      this->buffer[nread++] = Request_new(this->requestid++,queryseq1);
    }
    this->inputid++;
  }
  /* fprintf(stderr,"Read %d reads\n",nread); */

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}


#elif defined(GSNAP) || defined(GFILTER)

/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  Shortread_T queryseq1, queryseq2;
  bool skipp;
  int nchars1 = 0, nchars2 = 0;		/* Returned only because MPI master needs it.  Doesn't need to be saved as a field in Inbuffer_T. */

  /* fprintf(stderr,"Entered fill_buffer\n"); */
  while (nread < this->nspaces &&
	 (queryseq1 = Shortread_read(&this->nextchar,&nchars1,&nchars2,&queryseq2,
				     &this->input,&this->input2,
#ifdef HAVE_ZLIB
				     &this->gzipped,&this->gzipped2,
#endif
#ifdef HAVE_BZLIB
				     &this->bzipped,&this->bzipped2,
#endif
				     this->interleavedp,
				     this->read_files_command,&this->files,&this->nfiles,single_cell_p,
				     skipp = (this->inputid % this->part_interval != this->part_modulus))) != NULL) {
    if (skipp) {
#if 0
      /* Shortread procedures won't allocate in this situation */
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
#endif
      
    } else if (filter_if_both_p == true &&
	       Shortread_filterp(queryseq1) == true && (queryseq2 == NULL || Shortread_filterp(queryseq2) == true)) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else if (filter_if_both_p == false &&
	       (Shortread_filterp(queryseq1) == true || (queryseq2 != NULL && Shortread_filterp(queryseq2) == true))) {
      Shortread_free(&queryseq1);
      if (queryseq2 != NULL) {
	Shortread_free(&queryseq2);
      }
      
    } else {
      this->buffer[nread++] = Request_new(this->requestid++,queryseq1,queryseq2);
    }
    this->inputid++;
  }
  /* fprintf(stderr,"Read %d reads\n",nread); */

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#else
	 
/* GMAP version */
/* Returns number of requests read */
static unsigned int
fill_buffer (T this) {
  unsigned int nread = 0;
  Request_T request;
  Sequence_T genomeseq, queryseq;
  Genome_T genome, genomealt;
  List_T p;

  /* Take care of overflow */
  while (nread < this->nspaces && this->overflow_head != NULL) {
    this->overflow_head = List_pop_in(this->overflow_head,(void **) &request);
    this->buffer[nread++] = request;
  }
  if (this->overflow_head == NULL) {
    this->overflow_tail = (List_T) NULL;
  }


  /* Handle new inputs */
  while (nread < this->nspaces &&
	 (this->user_pairalign_p == false ||
	  (genomeseq = Sequence_read_multifile(&this->nextchar,&this->input,
					       this->read_files_command,&this->files,&this->nfiles)) != NULL) &&
	 (queryseq = Sequence_read_multifile(&this->nextchar,&this->input,
					     this->read_files_command,&this->files,&this->nfiles)) != NULL) {
    if (this->inputid % this->part_interval != this->part_modulus) {
      Sequence_free(&genomeseq);
      Sequence_free(&queryseq);

    } else if (this->user_pairalign_p == true) {
      debug(printf("inbuffer creating request %d\n",this->requestid));
      genome = genomealt = Genome_from_sequence(genomeseq);
      this->buffer[nread++] = Request_new(this->requestid++,genome,genomealt,queryseq,
					  /*free_genome_p*/true);
      Sequence_free(&genomeseq);

    } else if (this->user_genomes == NULL) {
      debug(printf("inbuffer creating request %d\n",this->requestid));
      this->buffer[nread++] = Request_new(this->requestid++,/*genome*/NULL,/*genomealt*/NULL,queryseq,
					  /*free_genome_p*/false);
    } else {
      for (p = this->user_genomes; p != NULL; p = List_next(p)) {
	genome = genomealt = (Genome_T) List_head(p);
	debug(printf("inbuffer creating request %d with genome %p\n",this->requestid,genome));
	if (nread < this->nspaces) {
	  this->buffer[nread++] = Request_new(this->requestid++,genome,genomealt,queryseq,
					      /*free_genome_p*/false);
	} else {
	  /* Put into overflow */
	  this->overflow_tail = List_unshift_in(&this->overflow_head,this->overflow_tail,
						(void *) Request_new(this->requestid++,genome,genomealt,queryseq,
								     /*free_genome_p*/false));
	}
      }
    }

    this->inputid++;
  }

  this->nleft = nread;
  this->ptr = 0;

  return nread;
}

#endif


/* No need to lock, since only main thread calls */
/* Returns nread to give to Outbuffer_new */
unsigned int
Inbuffer_fill_init (T this) {
  unsigned int nread;

  debug(printf("inbuffer filling initially\n"));
  nread = fill_buffer(this);
  debug(printf("inbuffer read %d sequences\n",nread));

  return nread;
}
  

Request_T
Inbuffer_get_request (T this) {
  Request_T request;
  unsigned int nread;

  debug(printf("Calling Inbuffer_get_request\n"));

#if defined(HAVE_PTHREAD)
  pthread_mutex_lock(&this->lock);
#endif
  
  if (this->nleft > 0) {
    request = this->buffer[this->ptr++];
    this->nleft -= 1;

#if 0
  } else if (this->nextchar == EOF) {
    /* Causes --force-single-end to fail when reads in a file are a multiple of nspaces */
    /* Want to call fill_buffer to find out if the input is exhausted */
    Outbuffer_add_nread(this->outbuffer,/*nread*/0);
    request = NULL;
#endif

  } else {
    debug(printf("inbuffer filling with nextchar %c (%d)\n",this->nextchar,this->nextchar));

    nread = fill_buffer(this);

    Outbuffer_add_nread(this->outbuffer,nread);
    debug(printf("inbuffer read %d sequences\n",nread));
    
    if (nread == 0) {
      /* Still empty */
      request = NULL;
    } else {
      request = this->buffer[this->ptr++];
      this->nleft -= 1;
    }
  }

#if defined(HAVE_PTHREAD)
  pthread_mutex_unlock(&this->lock);
#endif

  return request;
}


#if 0
/* Previously used by GMAP for selfalign feature. */
/* Same as Inbuffer_get_request, but leaves sequence in buffer. */
#ifndef GSNAP
Request_T
Inbuffer_first_request (T this) {
  Request_T request;
  unsigned int nread;

  debug(printf("Calling Inbuffer_first_request with %d left\n",this->nleft));

#if defined(HAVE_PTHREAD)
  pthread_mutex_lock(&this->lock);
#endif
  
  if (this->nleft > 0) {
    request = this->buffer[this->ptr/*++*/];
    /* this->nleft -= 1; */

  } else {
    debug(printf("inbuffer filling\n"));
    nread = fill_buffer(this);
    Outbuffer_add_nread(this->outbuffer,nread);
    debug(printf("inbuffer read %d sequences\n",nread));
    
    if (nread == 0) {
      /* Still empty */
      request = NULL;
    } else {
      request = this->buffer[this->ptr/*++*/];
      /* this->nleft -= 1; */
    }
  }

#if defined(HAVE_PTHREAD)
  pthread_mutex_unlock(&this->lock);
#endif

  return request;
}
#endif
#endif


