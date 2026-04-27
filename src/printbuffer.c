static char rcsid[] = "$Id: 2dd3ad12e9b7a144f7cd36167c6a4d1ed9b0416e $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "printbuffer.h"
#include <stdlib.h>
#include <string.h>

#include "assert.h"
#include "mem.h"
#include "list.h"
#include "samheader.h"


/* If we define USE_SETVBUF in outbuffer.c, this value may need to match */
#define USE_BUFFERED 1
#define OUTPUTLEN 65536


#define T Printbuffer_T
struct T {
  unsigned int nlines;

#if defined(GFILTER)
  FILE *output1;
  FILE *output2;

  List_T lines1, tails1;
  List_T lines2, tails2;

#elif defined(GSNAP)
  bool split_simple_p;
  char *split_output_root;
  char *output_file;

  List_T *lines, *tails;
  List_T lines_failedinput, tail_failedinput;
  List_T lines_failedinput_1, tail_failedinput_1;
  List_T lines_failedinput_2, tail_failedinput_2;

#else  /* GEXACT or GMAP */
  bool split_simple_p;
  char *split_output_root;
  char *output_file;

  List_T *lines, *tails;
  List_T lines_failedinput, tail_failedinput;
#endif
};


unsigned int
Printbuffer_nlines (T this) {
  return this->nlines;
}



#if defined(GFILTER)
T
Printbuffer_new (char *output_root, bool paired_end_p, bool appendp) {
  T new = (T) MALLOC_KEEP(sizeof(*new));
  char *filename;
  char *write_mode = (appendp == true) ? "a" : "w";

  new->nlines = 0;

  if (paired_end_p == false) {
    if ((new->output1 = fopen(output_root,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",output_root);
      exit(9);
    } else {
#ifdef USE_SETVBUF
      setvbuf(new->output1,buffer_output1,_IOFBF,OUTPUTLEN);
#endif
    }
    new->output2 = (FILE *) NULL;

  } else {
    filename = (char *) CALLOC(strlen(output_root)+strlen("_1")+1,sizeof(char));

    sprintf(filename,"%s_1",output_root);
    if ((new->output1 = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    } else {
#ifdef USE_SETVBUF
      setvbuf(new->output1,buffer_output1,_IOFBF,OUTPUTLEN);
#endif
    }
    
    sprintf(filename,"%s_2",output_root);
    if ((new->output2 = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
    } else {
#ifdef USE_SETVBUF
      setvbuf(new->output2,buffer_output2,_IOFBF,OUTPUTLEN);
#endif
    }

    FREE(filename);
  }

  new->lines1 = (List_T) NULL;
  new->tails1 = (List_T) NULL;

  new->lines2 = (List_T) NULL;
  new->tails2 = (List_T) NULL;

  return new;
}

#elif defined(GSNAP)
T
Printbuffer_new (bool split_simple_p, char *split_output_root, char *output_file) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlines = 0;
  new->split_simple_p = split_simple_p;
  new->split_output_root = split_output_root;
  new->output_file = output_file;

  if (split_simple_p == true) {
    new->lines = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS_SIMPLE+1),sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS_SIMPLE+1),sizeof(List_T));

  } else if (split_output_root != NULL) {
    new->lines = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS+1),sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS+1),sizeof(List_T));

  } else {
    new->lines = (List_T *) CALLOC_KEEP(1,sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP(1,sizeof(List_T));
  }
    
  new->lines_failedinput = new->tail_failedinput = (List_T) NULL;
  new->lines_failedinput_1 = new->tail_failedinput_1 = (List_T) NULL;
  new->lines_failedinput_2 = new->tail_failedinput_2 = (List_T) NULL;

  return new;
}

#else  /* GEXACT or GMAP */
T
Printbuffer_new (bool split_simple_p, char *split_output_root, char *output_file) {
  T new = (T) MALLOC_KEEP(sizeof(*new));

  new->nlines = 0;
  new->split_simple_p = split_simple_p;
  new->split_output_root = split_output_root;
  new->output_file = output_file;

  if (split_simple_p == true) {
    new->lines = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS_SIMPLE+1),sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS_SIMPLE+1),sizeof(List_T));

  } else if (split_output_root != NULL) {
    new->lines = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS+1),sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP((N_SPLIT_OUTPUTS+1),sizeof(List_T));

  } else {
    new->lines = (List_T *) CALLOC_KEEP(1,sizeof(List_T));
    new->tails = (List_T *) CALLOC_KEEP(1,sizeof(List_T));
  }
    
  new->lines_failedinput = new->tail_failedinput = (List_T) NULL;

  return new;
}
#endif

#if defined(GFILTER)
void
Printbuffer_free (T *old) {
  /* Individual List_T objects are freed by Outbuffer_lines_print */
  fclose((*old)->output1);
  if ((*old)->output2 != NULL) {
    fclose((*old)->output2);
  }
  FREE_KEEP(*old);
  return;
}

#else
void
Printbuffer_free (T *old) {
  /* Individual List_T objects are freed by Outbuffer_lines_print */
  FREE_KEEP((*old)->tails);
  FREE_KEEP((*old)->lines);
  FREE_KEEP(*old);
  return;
}
#endif


#if defined(GFILTER)
void
Printbuffer_store (T this, char *string1, char *string2) {

  this->tails1 = List_unshift_out(&(this->lines1),this->tails1,(void *) string1);
  this->tails2 = List_unshift_out(&(this->lines2),this->tails2,(void *) string2);
  this->nlines += 2;

  return;
}

#else

void
Printbuffer_store (T this, SAM_split_output_type split_output, char *string,
		   char *string_failedinput
#ifdef GSNAP
		   , char *string_failedinput_1, char *string_failedinput_2
#endif
		   ) {

  if (this->split_simple_p == true) {
    assert(split_output != 0);
    if (split_output == OUTPUT_NONE) {
      FREE(string);
    } else {
      this->tails[split_output] = List_unshift_out(&(this->lines[split_output]),this->tails[split_output],
						   (void *) string);
      this->nlines++;
    }

  } else if (this->split_output_root != NULL) {
    assert(split_output != 0);
    if (split_output == OUTPUT_NONE) {
      FREE(string);
    } else {
      this->tails[split_output] = List_unshift_out(&(this->lines[split_output]),this->tails[split_output],
						   (void *) string);
      this->nlines++;
    }

  } else if (split_output == OUTPUT_NONE) {
    FREE(string); /* ? No need to free string, since it is NULL */

  } else {
    this->tails[OUTPUT_NOT_SPLIT] = List_unshift_out(&(this->lines[OUTPUT_NOT_SPLIT]),
						     this->tails[OUTPUT_NOT_SPLIT],(void *) string);
    this->nlines++;
  }
  
  if (string_failedinput != NULL) {
    this->tail_failedinput = List_unshift_out(&(this->lines_failedinput),
					      this->tail_failedinput,(void *) string_failedinput);
  }

#ifdef GSNAP
  if (string_failedinput_1 != NULL) {
    this->tail_failedinput_1 = List_unshift_out(&(this->lines_failedinput_1),
						this->tail_failedinput_1,(void *) string_failedinput_1);
  }
  if (string_failedinput_2 != NULL) {
    this->tail_failedinput_2 = List_unshift_out(&(this->lines_failedinput_2),
						this->tail_failedinput_2,(void *) string_failedinput_2);
  }
#endif

  return;
}
#endif


#ifdef USE_BUFFERED
/* Taken from sam_sort */
static void
print_lines (FILE *fp_output, List_T *output_list, List_T *tail) {
  char buffer[OUTPUTLEN], *ptr;
  char *line, *p;
  int total_linelength, linelength, allowed;
  List_T l;

#ifdef USE_PUSH
  *output_list = List_reverse(*output_list);
#endif

  total_linelength = 0;
  ptr = &(buffer[0]);
  for (l = *output_list; l != NULL; l = List_next(l)) {
    p = line = (char *) List_head(l);
    if (line == NULL) {
      /* Possible with gff3 output */
    } else {
      linelength = strlen(line);

      while (total_linelength + linelength >= OUTPUTLEN) {
	allowed = OUTPUTLEN - total_linelength;
	strncpy(ptr,p,allowed);
#ifdef USE_WRITE
	write(fileno(fp_output),buffer,OUTPUTLEN*sizeof(char));
#else
	fwrite(buffer,sizeof(char),OUTPUTLEN,fp_output);
#endif

	ptr = &(buffer[0]);
	total_linelength = 0;

	p += allowed;
	linelength -= allowed;
      }
      
      strcpy(ptr,p);
      ptr += linelength;
      total_linelength += linelength;

      FREE_OUT(line);
    }
  }

  if (total_linelength > 0) {
#ifdef USE_WRITE
    write(fileno(fp_output),buffer,total_linelength*sizeof(char));
#else
    fwrite(buffer,sizeof(char),total_linelength,fp_output);
#endif
  }

  List_free_out(&(*output_list));
  *tail = (List_T) NULL;

  return;
}

#else
/* Prints all lines in a single system call */
static void
print_lines (FILE *fp_output, List_T *output_list, List_T *tail) {
  char *buffer, *ptr;
  char *line, *p;
  size_t total_linelength, linelength, allowed;
  List_T l;

#ifdef USE_PUSH
  *output_list = List_reverse(*output_list);
#endif

  total_linelength = 0;
  for (l = *output_list; l != NULL; l = List_next(l)) {
    line = (char *) List_head(l);
    if (line == NULL) {
      /* Possible with gff3 output */
    } else {
      total_linelength += strlen(line);
    }
  }

  ptr = buffer = MALLOC_OUT((total_linelength+1)*sizeof(char));
  for (l = *output_list; l != NULL; l = List_next(l)) {
    line = (char *) List_head(l);
    if (line == NULL) {
      /* Possible with gff3 output */
    } else {
      strcpy(ptr,line);
      ptr += strlen(line);
      FREE_OUT(line);
    }
  }
  
#ifdef USE_WRITE
  write(fileno(fp_output),buffer,total_linelength*sizeof(char));
#else
  fwrite(buffer,sizeof(char),total_linelength,fp_output);
#endif

  FREE_OUT(buffer);

  List_free_out(&(*output_list));
  *tail = (List_T) NULL;

  return;
}
#endif




#if defined(GFILTER)
void
Printbuffer_print (T this) {

  if (this->lines1 != NULL) {
    print_lines(this->output1,&(this->lines1),&(this->tails1));
  }

  if (this->lines2 != NULL) {
    print_lines(this->output2,&(this->lines2),&(this->tails2));
  }

  this->nlines = 0;

  return;
}

#else
void
Printbuffer_print (T this, FILE **outputs, FILE *output_failedinput,
#ifdef GSNAP
		   FILE *output_failedinput_1, FILE *output_failedinput_2,
#endif
		   bool paired_end_p, bool appendp) {
  SAM_split_output_type split_output;
  FILE *output;

  if (this->split_simple_p == true) {
    for (split_output = OUTPUT_NONE+1; split_output <= N_SPLIT_OUTPUTS_SIMPLE; split_output++) {
      if (this->lines[split_output] != NULL) {
	if ((output = outputs[split_output]) == NULL) {
	  output = outputs[split_output] = SAM_header_fopen(split_output,/*split_simple_p*/true,/*fileroot*/this->split_output_root,
							    paired_end_p,appendp);
#ifdef USE_SETVBUF
	  setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
	  SAM_header_print_all(output);
	}
	print_lines(output,&(this->lines[split_output]),&(this->tails[split_output]));
	this->nlines = 0;
      }
    }

  } else if (this->split_output_root != NULL) {
    /* Don't print file headers to OUTPUT_NONE (for --nofails, --omit-concordant-uniq, or --omit-concordant-mult */
    for (split_output = OUTPUT_NONE+1; split_output <= N_SPLIT_OUTPUTS; split_output++) {
      if (this->lines[split_output] != NULL) {
	if ((output = outputs[split_output]) == NULL) {
	  output = outputs[split_output] = SAM_header_fopen(split_output,/*split_simple_p*/false,/*fileroot*/this->split_output_root,
							    paired_end_p,appendp);
#ifdef USE_SETVBUF
	  setvbuf(outputs[split_output],buffer_outputs[split_output],_IOFBF,OUTPUTLEN);
#endif
	  SAM_header_print_all(output);
	}
	print_lines(output,&(this->lines[split_output]),&(this->tails[split_output]));
	this->nlines = 0;
      }
    }

  } else if (this->lines[OUTPUT_NOT_SPLIT] != NULL) {
    if ((output = outputs[OUTPUT_NOT_SPLIT]) == NULL) {
      /* Need to write file headers and assign output for future use */
      if (this->output_file == NULL) {
	output = outputs[OUTPUT_STDOUT] = stdout;
#ifdef USE_SETVBUF
	setvbuf(outputs[OUTPUT_STDOUT],buffer_outputs[OUTPUT_STDOUT],_IOFBF,OUTPUTLEN);
#endif
      } else {
	output = outputs[OUTPUT_FILE] = SAM_header_fopen(/*split_output*/OUTPUT_FILE,/*split_simple_p*/false,/*fileroot*/this->output_file,
							 paired_end_p,appendp);
#ifdef USE_SETVBUF
	setvbuf(outputs[OUTPUT_FILE],buffer_outputs[OUTPUT_FILE],_IOFBF,OUTPUTLEN);
#endif
      }
      SAM_header_print_all(output);
    }
    print_lines(output,&(this->lines[OUTPUT_NOT_SPLIT]),&(this->tails[OUTPUT_NOT_SPLIT]));
    this->nlines = 0;
  }

  if (this->lines_failedinput != NULL) {
    print_lines(output_failedinput,&(this->lines_failedinput),&(this->tail_failedinput));
  }

#ifdef GSNAP
  if (this->lines_failedinput_1 != NULL) {
    print_lines(output_failedinput_1,&(this->lines_failedinput_1),&(this->tail_failedinput_1));
  }
  if (this->lines_failedinput_2 != NULL) {
    print_lines(output_failedinput_2,&(this->lines_failedinput_2),&(this->tail_failedinput_2));
  }
#endif

  return;
}

#endif
