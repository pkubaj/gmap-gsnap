static char rcsid[] = "$Id: 3529459b95fcdd128017e80db0e0ff6b525111b9 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "outbuffer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif

#include "assert.h"
#include "bool.h"
#include "mem.h"
#include "samheader.h"
#include "printbuffer.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* RRlist_dump */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* USE_WRITE use write instead of fwrite.  May not be portable, and is
   also not thread-safe, so it can cause a race condition where
   alignment output is printed before headers.  Therefore, do not
   use */
/* #define USE_WRITE 0 */

/* #define USE_SETVBUF 1 */


static bool paired_end_p;
static bool appendp;

#if !defined(GFILTER)
/* sam-to-bam conversions always need the headers */
#define SAM_HEADERS_ON_EMPTY_FILES 1

static bool any_circular_p;
static bool quiet_if_excessive_p;
#endif

#if defined(GFILTER)
static char *output_root;
#else
static bool split_simple_p;
static char *split_output_root;
static char *failedinput_root;
#endif

static char *write_mode;

#if defined(GFILTER)
static FILE *output1 = NULL;
static FILE *output2 = NULL;
static char *write_mode;

#elif defined(GSNAP)
static char *output_file;
static FILE **outputs = NULL;
static FILE *output_failedinput;
static FILE *output_failedinput_1;
static FILE *output_failedinput_2;
#else
/* GEXACT or GMAP */
static char *output_file;
static FILE **outputs = NULL;
static FILE *output_failedinput;
#endif


/* If we want to use setvbuf, we would need to coordinate with printbuffer.c */
#ifdef USE_SETVBUF
/* OUTPUTLEN defines the size for the fwrite buffer after a call to
   setvbuf.  The typical fwrite buffer is 8 K, so 64 K results in
   fewer system calls. */

#define OUTPUTLEN 65536

#ifdef (GFILTER)
static char *buffer_output1 = NULL;
static char *buffer_output2 = NULL;
#elif defined(GSNAP)
static char **buffer_outputs;
static char buffer_failedinput[OUTPUTLEN];
static char buffer_failedinput_1[OUTPUTLEN];
static char buffer_failedinput_2[OUTPUTLEN];
static char buffer_failedinput[OUTPUTLEN];
#else /* GEXACT or GMAP */
static char **buffer_outputs;
static char buffer_failedinput[OUTPUTLEN];
#endif

#endif	/* USE_SETVBUF */



#if 0
#ifndef GSNAP
/* Taken from Univ_IIT_dump_sam */
static void
dump_sam_genome (FILE *fp, Genome_T genome,
		 char *sam_read_group_id, char *sam_read_group_name,
		 char *sam_read_group_library, char *sam_read_group_platform) {

  fprintf(fp,"@SQ\tSN:%s",Genome_accession(genome));
  fprintf(fp,"\tLN:%u",Genome_genomelength(genome));
  fprintf(fp,"\n");

  if (sam_read_group_id != NULL) {
    fprintf(fp,"@RG\tID:%s",sam_read_group_id);
    if (sam_read_group_platform != NULL) {
      fprintf(fp,"\tPL:%s",sam_read_group_platform);
    }
    if (sam_read_group_library != NULL) {
      fprintf(fp,"\tLB:%s",sam_read_group_library);
    }
    fprintf(fp,"\tSM:%s",sam_read_group_name);
    fprintf(fp,"\n");
  }

  return;
}
#endif
#endif


#if defined(GFILTER)

#else

static void
failedinput_open (char *failedinput_root) {
  char *filename;

#if defined(GSNAP)
  /* GSNAP*/
  filename = (char *) MALLOC((strlen(failedinput_root)+strlen(".1")+1) * sizeof(char));
  sprintf(filename,"%s.1",failedinput_root);

  if ((output_failedinput_1 = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  } else {
#ifdef USE_SETVBUF
    setvbuf(output_failedinput_1,buffer_failedinput_1,_IOFBF,OUTPUTLEN);
#endif
  }

  /* Re-use filename, since it is the same length */
  sprintf(filename,"%s.2",failedinput_root);
  if ((output_failedinput_2 = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  } else {
#ifdef USE_SETVBUF
    setvbuf(output_failedinput_2,buffer_failedinput_2,_IOFBF,OUTPUTLEN);
#endif
  }

  /* Re-use filename, since it is shorter */
  sprintf(filename,"%s",failedinput_root);
  if ((output_failedinput = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  } else {
#ifdef USE_SETVBUF
    setvbuf(output_failedinput,buffer_failedinput,_IOFBF,OUTPUTLEN);
#endif
  }
  FREE(filename);

#else
  /* GEXACT or GMAP */
  filename = (char *) MALLOC((strlen(failedinput_root)+1) * sizeof(char));
  sprintf(filename,"%s",failedinput_root);
  if ((output_failedinput = fopen(filename,write_mode)) == NULL) {
    fprintf(stderr,"Cannot open file %s for writing\n",filename);
    exit(9);
  } else {
#ifdef USE_SETVBUF
    setvbuf(output_failedinput,buffer_failedinput,_IOFBF,OUTPUTLEN);
#endif
  }
  FREE(filename);
#endif  

  return;
}

#endif


#if defined(GFILTER)
void
Outbuffer_setup (bool paired_end_p_in, bool appendp_in, char *output_root_in) {

  paired_end_p = paired_end_p_in;
  appendp = appendp_in; 
  output_root = output_root_in;
  write_mode = (appendp == true) ? "a" : "w";

#ifdef USE_SETVBUF
  buffer_output1 = (char *) NULL;
  buffer_output2 = (char *) NULL;
#endif

  return;
}

#else  /* GSNAP, GEXACT, or GMAP */

void
Outbuffer_setup (bool any_circular_p_in, bool quiet_if_excessive_p_in,
		 bool paired_end_p_in, bool appendp_in, char *output_file_in,
		 bool split_simple_p_in, char *split_output_root_in, char *failedinput_root_in) {

  any_circular_p = any_circular_p_in;
  quiet_if_excessive_p = quiet_if_excessive_p_in;
  paired_end_p = paired_end_p_in;
  appendp = appendp_in; 
  split_simple_p = split_simple_p_in;
  split_output_root = split_output_root_in;
  output_file = output_file_in;


  /************************************************************************/
  /* Output files */
  /************************************************************************/

  /* Only the output thread needs to run fopen, and can open files when needed */
  if (appendp == true) {
    write_mode = "a";
  } else {
    write_mode = "w";
  }

  outputs = (FILE **) CALLOC_KEEP(1+N_SPLIT_OUTPUTS,sizeof(FILE *));
#ifdef USE_SETVBUF
  buffer_outputs = (char **) CALLOC_KEEP(1+N_SPLIT_OUTPUTS,sizeof(char *));
#endif


  /************************************************************************/
  /* Failed input files */
  /************************************************************************/

  failedinput_root = failedinput_root_in;
  if (failedinput_root == NULL) {
    output_failedinput = NULL;
#ifdef GSNAP
    output_failedinput_1 = output_failedinput_2 = NULL;
#endif
  } else {
    failedinput_open(failedinput_root);
  }

  return;
}

#endif


#if defined(GFILTER)

void
Outbuffer_cleanup () {

#ifdef USE_SETVBUF
  if (buffer_output2 != NULL) {
    FREE_KEEP(buffer_output2);
  }
  if (buffer_output1 != NULL) {
    FREE_KEEP(buffer_output1);
  }
#endif

  return;
}

#else

void
Outbuffer_cleanup () {
#ifdef USE_SETVBUF
  SAM_split_output_type split_output;
#endif
  
#ifdef USE_SETVBUF
  for (split_output = 0; split_output <= N_SPLIT_OUTPUTS; split_output++) {
    if (buffer_outputs[split_output] != NULL) {
      FREE_KEEP(buffer_outputs[split_output]);
    }
  }
  FREE_KEEP(buffer_outputs); /* Matches CALLOC_KEEP in Outbuffer_setup */
#endif
  
  FREE_KEEP(outputs);		/* Matches CALLOC_KEEP in Outbuffer_setup */
  return;
}

#endif


typedef struct RRlist_T *RRlist_T;
struct RRlist_T {
  unsigned int id;

#if defined(GFILTER)
  char *string1;
  char *string2;
#elif defined(GSNAP)
  SAM_split_output_type split_output;
  char *string;
  char *string_failedinput;
  char *string_failedinput_1;
  char *string_failedinput_2;
#else
  /* GEXACT or GMAP */
  SAM_split_output_type split_output;
  char *string;
  char *string_failedinput;
#endif

  RRlist_T next;
};


#ifdef DEBUG1
static void
RRlist_dump (RRlist_T head, RRlist_T tail) {
  RRlist_T this;

  printf("head %p\n",head);
  for (this = head; this != NULL; this = this->next) {
    printf("%p: next %p\n",this,this->next);
  }
  printf("tail %p\n",tail);
  printf("\n");
  return;
}
#endif


#if defined(GFILTER)

/* Returns new tail, which must be reassigned to tail */
/* Allows fp_reads1 or fp_reads2 to be NULL, which result in NULL strings stored */
static RRlist_T
RRlist_unshift (RRlist_T *head, RRlist_T tail, int request_id,
		Filestring_T fp_reads1, Filestring_T fp_reads2) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new)); /* Called by worker thread */
  new->id = request_id;
  if (fp_reads1 == NULL) {
    new->string1 = (char *) NULL;
    new->string2 = (char *) NULL;
  } else if (fp_reads2 == NULL) {
    new->string1 = Filestring_string(fp_reads1);
    new->string2 = (char *) NULL;
    Filestring_free(&fp_reads1,/*free_string_p*/false);
  } else {
    new->string1 = Filestring_string(fp_reads1);
    new->string2 = Filestring_string(fp_reads2);
    Filestring_free(&fp_reads1,/*free_string_p*/false);
    Filestring_free(&fp_reads2,/*free_string_p*/false);
  }
  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}


/* Returns new head */
static RRlist_T
RRlist_pop (RRlist_T head, unsigned int *id, char **string1, char **string2) {
  RRlist_T newhead;

  *id = head->id;
  *string1 = head->string1;
  *string2 = head->string2;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}

#elif defined(GSNAP)

static RRlist_T
RRlist_unshift_pass1 (RRlist_T *head, RRlist_T tail, int request_id) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new)); /* Called by worker thread */
  new->id = request_id;
  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}

/* Returns new tail, which must be reassigned to tail */
static RRlist_T
RRlist_unshift (RRlist_T *head, RRlist_T tail, int request_id,
		Filestring_T fp, Filestring_T fp_failedinput,
		Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new)); /* Called by worker thread */
  new->id = request_id;
  new->split_output = Filestring_split_output(fp);
  new->string = Filestring_string(fp);
  if (fp_failedinput == NULL) {
    new->string_failedinput = (char *) NULL;
  } else {
    new->string_failedinput = Filestring_string(fp_failedinput);
  }

  if (fp_failedinput_1 == NULL) {
    new->string_failedinput_1 = (char *) NULL;
  } else {
    new->string_failedinput_1 = Filestring_string(fp_failedinput_1);
  }
  if (fp_failedinput_2 == NULL) {
    new->string_failedinput_2 = (char *) NULL;
  } else {
    new->string_failedinput_2 = Filestring_string(fp_failedinput_2);
  }

  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}


static RRlist_T
RRlist_pop_pass1 (RRlist_T head, unsigned int *id) {
  RRlist_T newhead;

  *id = head->id;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}


/* Returns new head */
static RRlist_T
RRlist_pop (RRlist_T head, unsigned int *id, SAM_split_output_type *split_output, char **string,
	    char **string_failedinput, char **string_failedinput_1, char **string_failedinput_2) {
  RRlist_T newhead;

  *id = head->id;
  *split_output = head->split_output;
  *string = head->string;
  *string_failedinput = head->string_failedinput;
  *string_failedinput_1 = head->string_failedinput_1;
  *string_failedinput_2 = head->string_failedinput_2;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}

#else /* GEXACT or GMAP */

/* Returns new tail, which must be reassigned to tail */
static RRlist_T
RRlist_unshift (RRlist_T *head, RRlist_T tail, int request_id,
		Filestring_T fp, Filestring_T fp_failedinput) {
  RRlist_T new;

  new = (RRlist_T) MALLOC_OUT(sizeof(*new)); /* Called by worker thread */
  new->id = request_id;
  new->split_output = Filestring_split_output(fp);
  new->string = Filestring_string(fp);
  if (fp_failedinput == NULL) {
    new->string_failedinput = (char *) NULL;
  } else {
    new->string_failedinput = Filestring_string(fp_failedinput);
  }
  new->next = (RRlist_T) NULL;
  
  if (*head == NULL) {		/* Equivalent to tail == NULL, but using *head avoids having to set tail in RRlist_pop */
    *head = new;
  } else {
    tail->next = new;
  }

  return new;
}


/* Returns new head */
static RRlist_T
RRlist_pop (RRlist_T head, unsigned int *id, SAM_split_output_type *split_output, char **string,
	    char **string_failedinput) {
  RRlist_T newhead;

  *id = head->id;
  *split_output = head->split_output;
  *string = head->string;
  *string_failedinput = head->string_failedinput;

  newhead = head->next;

  FREE_OUT(head);		/* Called by outbuffer thread */
  return newhead;
}

#endif


/* The queue is normally stored in descending order to make additions
   faster since we expect IDs to arrive as ascending order */
static RRlist_T
RRlist_insert (RRlist_T queue, RRlist_T head, unsigned int id) {
  RRlist_T *p;
  
  p = &queue;
  while (*p != NULL && id < (*p)->id) { /* Use < to store in descending order */
    p = &(*p)->next;
  }

  head->next = *p;
  *p = head;

  return queue;
}


static RRlist_T
RRlist_reverse (RRlist_T queue) {
  RRlist_T head = NULL, next;

  for ( ; queue; queue = next) {
    next = queue->next;
    queue->next = head;
    head = queue;
  }
  return head;
}



#define T Outbuffer_T
struct T {

#ifdef HAVE_PTHREAD
  pthread_mutex_t lock;
#endif

  unsigned int output_buffer_size;
  unsigned int nread;  /* Total reads input so far */
  unsigned int ntotal; /* A flag set to -1 while reading input.  After
			  an increment of nread is 0, signaling input
			  is done, then this is set to be nread */

  unsigned int nprocessed;	/* Incremented with each RRlist_unshift */
  /* unsigned int nqueued; */
  /* unsigned int noutput; */

#if defined(GFILTER)
  unsigned int npassed;
#endif

  RRlist_T head;
  RRlist_T tail;
  
#ifdef HAVE_PTHREAD
  pthread_cond_t filestring_avail_p;
#endif
};


/************************************************************************
 *   File routines
 ************************************************************************/


T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif

  new->output_buffer_size = output_buffer_size;
  new->nread = nread;

  /* Set to infinity until all reads are input.  Used for Pthreads version */
  new->ntotal = (unsigned int) -1;

  new->nprocessed = 0;

#if defined(GFILTER)
  new->npassed = 0;
#endif

  new->head = (RRlist_T) NULL;
  new->tail = (RRlist_T) NULL;

#ifdef HAVE_PTHREAD
  pthread_cond_init(&new->filestring_avail_p,NULL);
#endif

  return new;
}


#if !defined(GFILTER)

/* Open empty files, and add SAM headers if SAM_HEADERS_ON_EMPTY_FILES is set */
static void
touch_all_single_outputs (FILE **outputs, char *split_output_root, bool appendp) {
  SAM_split_output_type split_output;

  split_output = 1;
  while (split_output <= N_SPLIT_OUTPUTS_SINGLE_STD) {
    if (outputs[split_output] == NULL) {
      outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
      SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
      buffer_outputs[split_output] = (char *) MALLOC_KEEP(OUTPUTLEN*sizeof(char));
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
    }
    split_output++;
  }

  if (any_circular_p == false) {
    split_output = N_SPLIT_OUTPUTS_SINGLE_TOCIRC + 1;
  } else {
    while (split_output <= N_SPLIT_OUTPUTS_SINGLE_TOCIRC) {
      if (outputs[split_output] == NULL) {
	outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
        SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
	buffer_outputs[split_output] = (char *) MALLOC_KEEP(OUTPUTLEN*sizeof(char));
	setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      }
      split_output++;
    }
  }

  if (quiet_if_excessive_p == true) {
    while (split_output <= N_SPLIT_OUTPUTS_SINGLE) {
      if (outputs[split_output] == NULL) {
	outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
        SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
	buffer_outputs[split_output] = (char *) MALLOC_KEEP(OUTPUTLEN*sizeof(char));
	setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      }
      split_output++;
    }
  }

  return;
}


/* Open empty files, and add SAM headers if SAM_HEADERS_ON_EMPTY_FILES is set */
static void
touch_all_paired_outputs_simple (FILE **outputs, char *split_output_root, bool appendp) {
  SAM_split_output_type split_output;

  split_output = N_SPLIT_OUTPUTS_SINGLE + 1;

  while (split_output <= N_SPLIT_OUTPUTS_SIMPLE) {
    if (outputs[split_output] == NULL) {
      outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
      SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
    }
    split_output++;
  }

  return;
}


/* Open empty files, and add SAM headers if SAM_HEADERS_ON_EMPTY_FILES is set */
static void
touch_all_paired_outputs (FILE **outputs, char *split_output_root, bool appendp) {
  SAM_split_output_type split_output;

  split_output = N_SPLIT_OUTPUTS_SINGLE + 1;

  while (split_output <= N_SPLIT_OUTPUTS_STD) {
    if (outputs[split_output] == NULL) {
      outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
      SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
    }
    split_output++;
  }

  if (any_circular_p == false) {
    split_output = N_SPLIT_OUTPUTS_TOCIRC + 1;
  } else {
    while (split_output <= N_SPLIT_OUTPUTS_TOCIRC) {
      if (outputs[split_output] == NULL) {
	outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
        SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
	setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      }
      split_output++;
    }
  }

  if (quiet_if_excessive_p == true) {
    while (split_output <= N_SPLIT_OUTPUTS) {
      if (outputs[split_output] == NULL) {
	outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef SAM_HEADERS_ON_EMPTY_FILES
        SAM_header_print_all(outputs[split_output]);
#endif
#ifdef USE_SETVBUF
	setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      }
      split_output++;
    }
  }

  return;
}


static bool
paired_outputs_p (FILE **outputs) {
  SAM_split_output_type split_output;

  split_output = N_SPLIT_OUTPUTS_SINGLE + 1;
  while (split_output <= N_SPLIT_OUTPUTS) {
    if (outputs[split_output] != NULL) {
      return true;
    }
    split_output++;
  }

  return false;
}


static void
touch_all_files (FILE **outputs, char *split_output_root, bool appendp) {
  if (split_simple_p == true) {
    touch_all_paired_outputs_simple(outputs,split_output_root,appendp);
  } else {
    touch_all_single_outputs(outputs,split_output_root,appendp);
    if (paired_outputs_p(outputs) == true) {
      touch_all_paired_outputs(outputs,split_output_root,appendp);
    }
  }
  return;
}

#endif	/* GFILTER */



void
Outbuffer_close_files () {
#if !defined(GFILTER)
  SAM_split_output_type split_output;
#endif

#if defined(GFILTER)

#elif defined(GSNAP)
  if (failedinput_root != NULL) {
    fclose(output_failedinput);
    fclose(output_failedinput_1);
    fclose(output_failedinput_2);
  }
#else  /* GEXACT or GMAP */
  if (failedinput_root != NULL) {
    fclose(output_failedinput);
  }
#endif

#if defined(GFILTER)
  if (output1 != NULL) {
    fclose(output1);
  } else {
    /* Used printbuffer in multi-threaded mode instead of outbuffer */
  }
  if (output2 != NULL) {
    fclose(output2);
  } else {
    /* Single-end reads, or used printbuffer in multi-threaded mode instead of outbuffer */
  }

#else
  if (split_output_root != NULL) {
    touch_all_files(outputs,split_output_root,appendp);

    for (split_output = 1; split_output <= N_SPLIT_OUTPUTS; split_output++) {
      if (outputs[split_output] != NULL) {
	fclose(outputs[split_output]);
      }
    }
  } else if (output_file != NULL && outputs[OUTPUT_FILE] != NULL) {
    fclose(outputs[OUTPUT_FILE]);
  } else {
    /* Wrote to stdout */
  }

  FREE_KEEP(outputs);
#endif

  return;
}


void
Outbuffer_free (T *old) {

  if (*old) {
#ifdef HAVE_PTHREAD
    pthread_cond_destroy(&(*old)->filestring_avail_p);
    pthread_mutex_destroy(&(*old)->lock);
#endif

    FREE_KEEP(*old);
  }

  return;
}


unsigned int
Outbuffer_nread (T this) {
  return this->nread;
}

#if defined(GFILTER)
unsigned int
Outbuffer_npassed (T this) {
  return this->npassed;
}
#endif

void
Outbuffer_add_nread (T this, unsigned int nread) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  if (nread == 0) {
    /* Finished reading, so able to determine total reads in input */
    this->ntotal = this->nread;
    debug(fprintf(stderr,"__Outbuffer_add_nread added 0 reads, so setting ntotal to be %u\n",this->ntotal));

#ifdef HAVE_PTHREAD
    pthread_cond_signal(&this->filestring_avail_p);
#endif

  } else {
    this->nread += nread;
    debug(fprintf(stderr,"__Outbuffer_add_nread added %d read, now %d\n",nread,this->nread));
  }

#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}


#if defined(GFILTER)
/* Handles the case where fp_reads1 is NULL, indicating the read did not pass */
/* RRlist_unshift can handle NULL values for filestrings.
   RRlist_unshift takes care of freeing filestrings */
void
Outbuffer_put_filestrings (T this, int request_id, Filestring_T fp_reads1, Filestring_T fp_reads2) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  this->nprocessed += 1;
  if (fp_reads1 != (Filestring_T) NULL) {
    this->npassed += 1;
  }

  this->tail = RRlist_unshift(&this->head,this->tail,request_id,fp_reads1,fp_reads2);
  debug1(RRlist_dump(this->head,this->tail));

#ifdef HAVE_PTHREAD
  debug(printf("Signaling that filestring is available\n"));
  pthread_cond_signal(&this->filestring_avail_p);
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}

#elif defined(GSNAP)
void
Outbuffer_put_pass1 (T this, int request_id) {
  debug(printf("Outbuffer_put_pass1 called\n"));

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  this->tail = RRlist_unshift_pass1(&this->head,this->tail,request_id);

  debug1(RRlist_dump(this->head,this->tail));
  this->nprocessed += 1;

#ifdef HAVE_PTHREAD
  debug(printf("Signaling that filestring is available\n"));
  pthread_cond_signal(&this->filestring_avail_p);
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}


void
Outbuffer_put_filestrings (T this, int request_id,
			   Filestring_T fp, Filestring_T fp_failedinput,
			   Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  this->tail = RRlist_unshift(&this->head,this->tail,request_id,
			      fp,fp_failedinput,fp_failedinput_1,fp_failedinput_2);
  Filestring_free(&fp,/*free_string_p*/false);
  Filestring_free(&fp_failedinput,/*free_string_p*/false);
  Filestring_free(&fp_failedinput_1,/*free_string_p*/false);
  Filestring_free(&fp_failedinput_2,/*free_string_p*/false);

  debug1(RRlist_dump(this->head,this->tail));
  this->nprocessed += 1;

#ifdef HAVE_PTHREAD
  debug(printf("Signaling that filestring is available\n"));
  pthread_cond_signal(&this->filestring_avail_p);
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}

#else  /* GEXACT or GMAP */
void
Outbuffer_put_filestrings (T this, int request_id, Filestring_T fp, Filestring_T fp_failedinput) {

#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif

  this->tail = RRlist_unshift(&this->head,this->tail,request_id,fp,fp_failedinput);
  Filestring_free(&fp,/*free_string_p*/false);
  Filestring_free(&fp_failedinput,/*free_string_p*/false);
  debug1(RRlist_dump(this->head,this->tail));
  this->nprocessed += 1;

#ifdef HAVE_PTHREAD
  debug(printf("Signaling that filestring is available\n"));
  pthread_cond_signal(&this->filestring_avail_p);
  pthread_mutex_unlock(&this->lock);
#endif

  return;
}

#endif



/* Outbuffer_print_filestrings is used in single-thread mode.  The
   outbuffer object is not used, except to increment npassed for
   gfilter. */
/* In contrast, Outbuffer_put_filestrings and Printbuffer_print are
   used on multi-threaded mode */
#if defined(GFILTER)
/* Handles the case where the filestrings are NULL */
void
Outbuffer_print_filestrings (T this, Filestring_T fp_reads1, Filestring_T fp_reads2) {
  char *filename;

  if (paired_end_p == false) {
    if (fp_reads1 == NULL) {
      /* No output */
    } else {
      if (output1 == NULL) {
	if ((output1 = fopen(output_root,write_mode)) == NULL) {
	  fprintf(stderr,"Cannot open file %s for writing\n",output_root);
	}
      }
      Filestring_print(output1,fp_reads1);
      Filestring_free(&fp_reads1,/*free_string_p*/true);
      this->npassed += 1;
    }

  } else {
    if (fp_reads1 == NULL) {
      /* No output */
    } else {
      if (output1 == NULL) {
	filename = (char *) CALLOC(strlen(output_root)+strlen("_1")+1,sizeof(char));
	sprintf(filename,"%s_1",output_root);
	if ((output1 = fopen(filename,write_mode)) == NULL) {
	  fprintf(stderr,"Cannot open file %s for writing\n",filename);
	}
	sprintf(filename,"%s_2",output_root);
	if ((output2 = fopen(filename,write_mode)) == NULL) {
	  fprintf(stderr,"Cannot open file %s for writing\n",filename);
	}
	FREE(filename);
      }

      Filestring_print(output1,fp_reads1);
      Filestring_print(output2,fp_reads2);
      Filestring_free(&fp_reads1,/*free_string_p*/true);
      Filestring_free(&fp_reads2,/*free_string_p*/true);
      this->npassed += 1;
    }
  }

  return;
}

#elif defined(GSNAP)
void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2) {
  SAM_split_output_type split_output;
  FILE *output;

  if (split_simple_p == true) {
    split_output = Filestring_split_output(fp);
    if ((output = outputs[split_output]) == NULL) {
      output = outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
    }

  } else if (split_output_root != NULL) {
    split_output = Filestring_split_output(fp);
    if ((output = outputs[split_output]) == NULL) {
      output = outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      if (split_output == OUTPUT_NONE) {
	/* Don't print file headers, since no output will go to
	   stdout.  Must be a nomapping when --nofails is specified */
      } else {
	SAM_header_print_all(output);
      }
    }

  } else if ((output = outputs[OUTPUT_NOT_SPLIT]) == NULL) {
    if (output_file == NULL) {
      output = outputs[OUTPUT_STDOUT] = stdout;
      SAM_header_print_all(stdout);
    } else {
      output = outputs[OUTPUT_FILE] = SAM_header_fopen(/*split_output*/OUTPUT_FILE,split_simple_p,output_file,paired_end_p,appendp);
      SAM_header_print_all(output);
#ifdef USE_SETVBUF
      setvbuf(outputs[OUTPUT_FILE],buffer_outputs[OUTPUT_FILE],_IOFBF,OUTPUTLEN);
#endif
    }
  }

  Filestring_print(output,fp);
  Filestring_free(&fp,/*free_string_p*/true);

  if (failedinput_root != NULL) {
    if (fp_failedinput != NULL) {
      Filestring_print(output_failedinput,fp_failedinput);
      Filestring_free(&fp_failedinput,/*free_string_p*/true);
    }
    if (fp_failedinput_1 != NULL) {
      Filestring_print(output_failedinput_1,fp_failedinput_1);
      Filestring_free(&fp_failedinput_1,/*free_string_p*/true);
    }
    if (fp_failedinput_2 != NULL) {
      Filestring_print(output_failedinput_2,fp_failedinput_2);
      Filestring_free(&fp_failedinput_2,/*free_string_p*/true);
    }
  }

  return;
}

#else  /* GEXACT or GMAP */
void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput) {
  SAM_split_output_type split_output;
  FILE *output;

  if (split_simple_p == true) {
    split_output = Filestring_split_output(fp);
    if ((output = outputs[split_output]) == NULL) {
      output = outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
    }

  } else if (split_output_root != NULL) {
    split_output = Filestring_split_output(fp);
    if ((output = outputs[split_output]) == NULL) {
      output = outputs[split_output] = SAM_header_fopen(split_output,split_simple_p,split_output_root,paired_end_p,appendp);
#ifdef USE_SETVBUF
      setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
      if (split_output == OUTPUT_NONE) {
	/* Don't print file headers, since no output will go to
	   stdout.  Must be a nomapping when --nofails is specified */
      } else {
	SAM_header_print_all(output);
      }
    }

  } else if ((output = outputs[OUTPUT_NOT_SPLIT]) == NULL) {
    if (output_file == NULL) {
      output = outputs[OUTPUT_STDOUT] = stdout;
      SAM_header_print_all(stdout);
#ifdef USE_SETVBUF
      setvbuf(outputs[OUTPUT_STDOUT],buffer_outputs[OUTPUT_STDOUT],_IOFBF,OUTPUTLEN);
#endif
    } else {
      output = outputs[OUTPUT_FILE] = SAM_header_fopen(/*split_output*/OUTPUT_FILE,split_simple_p,output_file,paired_end_p,appendp);
      SAM_header_print_all(output);
#ifdef USE_SETVBUF
      setvbuf(outputs[OUTPUT_FILE],buffer_outputs[OUTPUT_FILE],_IOFBF,OUTPUTLEN);
#endif
    }
  }

  Filestring_print(output,fp);
  Filestring_free(&fp,/*free_string_p*/true);

  if (failedinput_root != NULL) {
    if (fp_failedinput != NULL) {
      Filestring_print(output_failedinput,fp_failedinput);
      Filestring_free(&fp_failedinput,/*free_string_p*/true);
    }
  }

  return;
}
#endif



#ifdef GSNAP
void *
Outbuffer_thread_pass1 (void *data) {
  T this = (T) data;
  unsigned int noutput = 0;	/* Incremented with each RRlist_pop */
  unsigned int ntotal, id;

  
  /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif
  ntotal = this->ntotal;
#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  while (noutput < ntotal) {	/* Previously checked against this->ntotal */
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->lock);
    while (this->head == NULL && noutput < this->ntotal) {
      /* Wait for Outbuffer_put_pass1 to push a set of elts onto the end of the list  */
      debug(fprintf(stderr,"__outbuffer_thread_pass1 waiting for filestring_avail_p\n"));
      pthread_cond_wait(&this->filestring_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_anyorder woke up\n"));
#endif

    while (this->head != NULL) {
      this->head = RRlist_pop_pass1(this->head,&id);
      noutput += 1;
    }

    /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
    ntotal = this->ntotal;
#ifdef HAVE_PTHREAD
    pthread_mutex_unlock(&this->lock);
#endif
  }
    
  assert(this->head == NULL);
  return (void *) NULL;
}
#endif


void *
Outbuffer_thread_anyorder (void *data) {
  T this = (T) data;
  unsigned int output_buffer_size = this->output_buffer_size;
  unsigned int noutput = 0;	/* Incremented with each RRlist_pop */
  unsigned int ntotal, id;
  Printbuffer_T printbuffer;

#if defined(GFILTER)
  char *string1, *string2;
#elif defined(GSNAP)
  SAM_split_output_type split_output;
  char *string, *string_failedinput;
  char *string_failedinput_1, *string_failedinput_2;
#else  /* GEXACT or GMAP */
  SAM_split_output_type split_output;
  char *string, *string_failedinput;
#endif
  
#if defined(GFILTER)
  printbuffer = Printbuffer_new(output_root,paired_end_p,appendp);
#else
  printbuffer = Printbuffer_new(split_simple_p,split_output_root,output_file);
#endif

  /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif
  ntotal = this->ntotal;
#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  while (noutput < ntotal) {	/* Previously checked against this->ntotal */
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->lock);
    while (this->head == NULL && noutput < this->ntotal) {
      /* Wait for Outbuffer_put_filestrings to push a set of filestrings onto the end of the list  */
      debug(fprintf(stderr,"__outbuffer_thread_anyorder waiting for filestring_avail_p\n"));
      pthread_cond_wait(&this->filestring_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_anyorder woke up\n"));
#endif

    while (this->head != NULL) {
#if defined(GFILTER)
      this->head = RRlist_pop(this->head,&id,&string1,&string2);
      if (string1 != NULL) {
	Printbuffer_store(printbuffer,string1,string2);
      }
#elif defined(GSNAP)
      this->head = RRlist_pop(this->head,&id,&split_output,&string,&string_failedinput,
			      &string_failedinput_1,&string_failedinput_2);
      Printbuffer_store(printbuffer,split_output,string,string_failedinput,
			string_failedinput_1,string_failedinput_2);
#else  /* GEXACT or GMAP */
      this->head = RRlist_pop(this->head,&id,&split_output,&string,&string_failedinput);
      Printbuffer_store(printbuffer,split_output,string,string_failedinput);
#endif
      debug1(RRlist_dump(this->head,this->tail));
      noutput += 1;
    }

    /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
    ntotal = this->ntotal;
    

    if (Printbuffer_nlines(printbuffer) < output_buffer_size) {
      /* Do unlock here to allow worker threads to add results while we print */
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif

#if defined(GFILTER)
      Printbuffer_print(printbuffer);

#elif defined(GSNAP)
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			output_failedinput_1,output_failedinput_2,
			paired_end_p,appendp);
#else  /* GEXACT or GMAP */
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			paired_end_p,appendp);

#endif

    } else {
      /* Slow worker threads during printing so that results do not
	 accumulate */
#if defined(GFILTER)
      Printbuffer_print(printbuffer);

#elif defined(GSNAP)
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			output_failedinput_1,output_failedinput_2,
			paired_end_p,appendp);

#else  /* GEXACT or GMAP */
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			paired_end_p,appendp);
#endif

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif
    }

    debug(fprintf(stderr,"__outbuffer_thread_anyorder has noutput %d, ntotal %d\n",
		  noutput,ntotal));
  }

  Printbuffer_free(&printbuffer);

  assert(this->head == NULL);
  return (void *) NULL;
}



void *
Outbuffer_thread_ordered (void *data) {
  T this = (T) data;
  unsigned int output_buffer_size = this->output_buffer_size;
  unsigned int noutput = 0;	/* Incremented with each RRlist_pop */
  unsigned int nqueued = 0, ntotal, queue_min_id, id;
  Printbuffer_T printbuffer;
#if defined(GFILTER)
  char *string1, *string2;
#elif defined(GSNAP)
  SAM_split_output_type split_output;
  char *string, *string_failedinput;
  char *string_failedinput_1, *string_failedinput_2;
#else  /* GEXACT or GMAP */
  SAM_split_output_type split_output;
  char *string, *string_failedinput;
#endif
  RRlist_T queue = NULL, next;

#if defined(GFILTER)
  printbuffer = Printbuffer_new(output_root,paired_end_p,appendp);
#else
  printbuffer = Printbuffer_new(split_simple_p,split_output_root,output_file);
#endif

  /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
#ifdef HAVE_PTHREAD
  pthread_mutex_lock(&this->lock);
#endif
  ntotal = this->ntotal;
#ifdef HAVE_PTHREAD
  pthread_mutex_unlock(&this->lock);
#endif

  while (noutput < ntotal) {	/* Previously checked against this->ntotal */
#ifdef HAVE_PTHREAD
    pthread_mutex_lock(&this->lock);
    while (this->head == NULL && noutput < this->ntotal) {
      pthread_cond_wait(&this->filestring_avail_p,&this->lock);
    }
    debug(fprintf(stderr,"__outbuffer_thread_ordered woke up\n"));
#endif

    while (this->head != NULL) {
      /* Store results in ordered queue (in descending order) */
      next = this->head->next;
      id = this->head->id;
      if (queue == NULL || id < queue_min_id) {
	queue_min_id = id;
      }
      queue = RRlist_insert(queue,this->head,id);
      this->head = next;
      nqueued++;
    }
    /* printf("Queue has %d items, starting with %u\n",nqueued,queue->id); */

    /* Obtain this->ntotal while locked, to prevent race between output thread and input thread */
    ntotal = this->ntotal;

    /* Return control to worker threads */
    pthread_mutex_unlock(&this->lock);

    if (queue_min_id == noutput) {
      /* Print stored queue as applicable  */
      /* Reverse queue to be in ascending order */
      queue = RRlist_reverse(queue);
      while (queue != NULL && queue->id == noutput) {
#if defined(GFILTER)
	queue = RRlist_pop(queue,&id,&string1,&string2);
	if (string1 != NULL) {
	  Printbuffer_store(printbuffer,string1,string2);
	}
#elif defined(GSNAP)
	queue = RRlist_pop(queue,&id,&split_output,&string,&string_failedinput,
			   &string_failedinput_1,&string_failedinput_2);
	Printbuffer_store(printbuffer,split_output,string,string_failedinput,
			  string_failedinput_1,string_failedinput_2);
#else  /* GEXACT or GMAP */
	queue = RRlist_pop(queue,&id,&split_output,&string,&string_failedinput);
	Printbuffer_store(printbuffer,split_output,string,string_failedinput);
#endif
	debug1(RRlist_dump(this->head,this->tail));
	nqueued--;
	noutput += 1;
      }

      if (queue != NULL) {
	queue_min_id = queue->id;
	/* Restore queue to be in descending order */
	queue = RRlist_reverse(queue);
      }
    }

    if (Printbuffer_nlines(printbuffer) < output_buffer_size) {
      /* Allow worker threads to add results while we print */
#if defined(GFILTER)
      Printbuffer_print(printbuffer);

#elif defined(GSNAP)
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			output_failedinput_1,output_failedinput_2,
			paired_end_p,appendp);

#else  /* GEXACT or GMAP */
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			paired_end_p,appendp);
#endif


    } else {
      /* Slow worker threads during printing so that results do not
	 accumulate */
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&this->lock);
#endif

#if defined(GFILTER)
      Printbuffer_print(printbuffer);
#elif defined(GSNAP)
      Printbuffer_print(printbuffer,outputs,output_failedinput,
			output_failedinput_1,output_failedinput_2,
			paired_end_p,appendp);
#else  /* GEXACT or GMAP */
      Printbuffer_print(printbuffer,outputs,output_failedinput,paired_end_p,appendp);
#endif

#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&this->lock);
#endif
    }

    debug(fprintf(stderr,"__outbuffer_thread_ordered has noutput %d, ntotal %d\n",
		  noutput,ntotal));
  }

  Printbuffer_free(&printbuffer);

  assert(queue == NULL);
  return (void *) NULL;
}

