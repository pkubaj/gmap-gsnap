static char rcsid[] = "$Id: 6a208f7437a5956a06259c86a01a194548c02af6 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "single-cell.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* For index */
#include <math.h>		/* For pow */

#include "mem.h"
#include "tableuint.h"
#include "getline.h"
#include "fopen.h"
#include "shortread.h"


#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif


/* The first part of read1 is the cell barcode, and the second part is the molecular barcode (UMI) */
#define CELL_BARCODE_LEN 16

static Tableuint_T whitelistp = NULL;
static Tableuint_T whitelist_counts = NULL;
static double whitelist_total_double = 0.0;
static int wellpos = 4;


#define HEADERLEN 512
#define DISCARDLEN 8192

static char Header[HEADERLEN];
static char Discard[DISCARDLEN];


#define MAX_SC_READ1_LENGTH 150
static char Read1[MAX_SC_READ1_LENGTH+1];
static char Quality[MAX_SC_READ1_LENGTH+1];


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static bool
skip_header (int *nchars, FILE *fp, int nextchar) {
  char *p;

  if (nextchar == EOF) {
    return false;
  } else if ((p = fgets(&(Header[0]),HEADERLEN,fp)) == NULL) {
    /* File must terminate after > */
    return false;
  } else {
    *nchars += strlen(p);
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while ((p = fgets(&(Discard[0]),DISCARDLEN,fp)) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) {
      *nchars += strlen(p);
    }
  }

  return true;
} 

#ifdef HAVE_ZLIB
static bool
skip_header_gzip (gzFile fp, int nextchar) {

  if (nextchar == EOF) {	/* was gzeof(fp) */
    return false;
  } else if (gzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (gzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 
#endif

#ifdef HAVE_BZLIB
static bool
skip_header_bzip2 (Bzip2_T fp, int nextchar) {

  if (nextchar == EOF) {	/* was bzeof(fp) */
    return false;
  } else if (bzgets(fp,&(Header[0]),HEADERLEN) == NULL) {
    /* File must terminate after > */
    return false;
  }

  if (rindex(&(Header[0]),'\n') == NULL) {
    /* Eliminate rest of header from input */
    while (bzgets(fp,&(Discard[0]),DISCARDLEN) != NULL &&
	   rindex(&(Discard[0]),'\n') == NULL) ;
  }

  return true;
} 
#endif



#define SPACE 32

static char *
find_spaces (int *nspaces, char *line) {
  char *first, *p2;

  p2 = index(line,SPACE);

  if (p2 == NULL /* && p3 == NULL*/) {
    return NULL;
  } else {
    first = p2++;
    while (*p2 != '\0' && *p2 == SPACE) {
      p2++;
    }

    *nspaces = p2 - first;
    return first;
  }
}


/* Simplified from version in shortread.c */
static int
input_oneline (int *nextchar, int *nchars, char *Start, FILE *fp) {
  char *ptr, *p = NULL;
  int strlenp, nspaces;


  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  if (*nextchar == EOF) {
    debug(printf("nchars %d: EOF or > or +: Returning 0\n",*nchars));
    return 0;
  } else if (*nextchar == '\n') {
    debug(printf("nchars %d: Blank line: Returning 0\n",*nchars));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = fgets(ptr,MAX_SC_READ1_LENGTH+1,fp)) == NULL) {
      /* NULL if file ends with a blank line */
      debug(printf("Blank line. read %s.\n",ptr));
    } else {
      *nchars += strlen(p);
      debug(printf("nchars %d: Read %s.\n",*nchars,ptr));
      while ((p = find_spaces(&nspaces,ptr)) != NULL) {
	ptr = p;
	p += nspaces;
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
	ptr[strlenp] = '\0';
	debug(printf("Found %d spaces.  Did memmove of %d chars at %p to %p\n",nspaces,strlenp,p,ptr));
      }

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("nchars %d: Now string is %s.\n",*nchars,ptr));
      } else if (*ptr == EOF) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("nchars %d: End of file seen\n",*nchars));
      } else {
	/* No line feed, but not end of file.  Read too long, so using another method. */
	fprintf(stderr,"No line feed, but not end of file\n");
	exit(9);
      }
    }

    ptr += strlen(ptr);

    /* Peek at character after eoln */
    *nextchar = fgetc(fp);
    *nchars += 1;

    debug(printf("nchars %d: Returning %ld with nextchar %c\n",*nchars,(ptr - &(Start[0]))/sizeof(char),*nextchar));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}


/* Simplified from version in shortread.c */
#ifdef HAVE_ZLIB
static int
input_oneline_gzip (int *nextchar, char *Start, gzFile fp) {
  char *ptr, *p = NULL;
  int strlenp, nspaces;


  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  if (*nextchar == EOF) {
    debug(printf("EOF or > or +: Returning 0\n"));
    return 0;
  } else if (*nextchar == '\n') {
    debug(printf("Blank line: Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = gzgets(fp,ptr,MAX_SC_READ1_LENGTH+1)) == NULL) {
      /* NULL if file ends with a blank line */
      debug(printf("Blank line. read %s.\n",ptr));
    } else {
      debug(printf("Read %s.\n",ptr));
      while ((p = find_spaces(&nspaces,ptr)) != NULL) {
	ptr = p;
	p += nspaces;
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
	ptr[strlenp] = '\0';
	debug(printf("Found %d spaces.  Did memmove of %d chars at %p to %p\n",nspaces,strlenp,p,ptr));
      }

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (*ptr == EOF) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file.  Read too long, so using another method. */
	fprintf(stderr,"No line feed, but not end of file\n");
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    *nextchar = gzgetc(fp);

    debug(printf("Returning %ld with nextchar %c\n",(ptr - &(Start[0]))/sizeof(char),*nextchar));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}
#endif


/* Simplified from version in shortread.c */
#ifdef HAVE_BZLIB
static int
input_oneline_bzip2 (int *nextchar, char *Start, Bzip2_T fp) {
  char *ptr, *p = NULL;
  int strlenp, nspaces;

  debug(printf("Entering input_oneline with nextchar = %c\n",*nextchar));

  ptr = &(Start[0]);
  if (*nextchar == EOF) {
    debug(printf("EOF or > or +: Returning 0\n"));
    return 0;
  } else if (*nextchar == '\n') {
    debug(printf("Blank line: Returning 0\n"));
    return 0;
  } else {
    *ptr++ = (char) *nextchar;
    if ((p = bzgets(fp,ptr,MAX_SC_READ1_LENGTH+1)) == NULL) {
      /* NULL if file ends with a blank line */
      debug(printf("Blank line. read %s.\n",ptr));
    } else {
      debug(printf("Read %s.\n",ptr));
      while ((p = find_spaces(&nspaces,ptr)) != NULL) {
	ptr = p;
	p += nspaces;
	strlenp = strlen(p);
	memmove(ptr,p,strlenp);
	ptr[strlenp] = '\0';
	debug(printf("Found %d spaces.  Did memmove of %d chars at %p to %p\n",nspaces,strlenp,p,ptr));
      }

      if (*ptr == '\n') {
	*ptr = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if ((p = index(ptr,'\n')) != NULL) {
	if (p[-1] == '\r') {
	  p--;
	}
	*p = '\0';
	debug(printf("Now string is %s.\n",ptr));
      } else if (*ptr == EOF) {
	/* No line feed, but end of file.  Handle below. */
	debug(printf("End of file seen\n"));
      } else {
	/* No line feed, but not end of file.  Read too long, so using another method. */
	fprintf(stderr,"No line feed, but not end of file\n");
	exit(9);
      }
    }
    ptr += strlen(ptr);

    /* Peek at character after eoln */
    *nextchar = bzgetc(fp);

    debug(printf("Returning %ld with nextchar %c\n",(ptr - &(Start[0]))/sizeof(char),*nextchar));
    return (ptr - &(Start[0]))/sizeof(char);
  }
}
#endif



/* Modified from Shortread_read_fastq_text in shortread.c */
static bool
update_whitelist_counts_text (unsigned int *whitelist_total,
			      Tableuint_T whitelist_counts, Tableuint_T whitelistp,
			      int *nextchar, int *nchars1, FILE **input1,
			      char *read_files_command, char ***files, int *nfiles) {
  int fulllength;
  unsigned int count;
  
  while (1) {
    if (*input1 == NULL || *nextchar == EOF) { /* was feof(input1) */
      if (*input1 != NULL) {
	if (read_files_command != NULL) {
	  pclose(*input1);
	} else {
	  fclose(*input1);
	}
	*input1 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return false;

      } else if (*nfiles == 1) {
	fprintf(stderr,"Error: odd number of fastq files\n");
	(*files) += 1;
	(*nfiles) -= 1;
	*nextchar = EOF;
	return false;

      } else {
	while (*nfiles > 0 &&
	       (*input1 = Fopen_read_text(read_files_command,(*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it and its pair\n",
		  (*files)[0]);
	  (*files) += 2;
	  (*nfiles) -= 2;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return false;
	} else {
	  (*files) += 2;
	  (*nfiles) -= 2;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init(&(*nchars1),*input1)) == EOF) {
	return false;
      }
    }

    debug(printf("** Getting header for fastq text\n"));
    if (skip_header(&(*nchars1),*input1,*nextchar) == false) {
      /* Skip */
      *nextchar = EOF;
    } else {
      *nextchar = fgetc(*input1);
      *nchars1 += 1;
      if ((fulllength = input_oneline(&(*nextchar),&(*nchars1),&(Read1[0]),*input1)) == 0) {
	/* Skip */
      } else {
	if (*nextchar == '+') {
	  /* Ignore quality */
	  skip_header(&(*nchars1),*input1,*nextchar);
	  *nextchar = fgetc(*input1);
	  *nchars1 += 1;
	  input_oneline(&(*nextchar),&(*nchars1),&(Quality[0]),*input1);
	}

	/* Can use Read1 because of string_compare_cell_barcode and string_hash_cell_barcode functions */
	/* strncpy(barcode,Read1,CELL_BARCODE_LEN); */
	/* barcode[CELL_BARCODE_LEN] = '\0'; */
	if (Tableuint_get(whitelistp,(const void *) Read1) == (unsigned int) true) {
	  count = Tableuint_get(whitelist_counts,(const void *) Read1);
	  Tableuint_put(whitelist_counts,(const void *) Read1,count+1);
	  *whitelist_total += 1;
	}

	return true;
      }
    }
  }
}


/* Modified from Shortread_read_fastq_gzip in shortread.c */
#ifdef HAVE_ZLIB
static bool
update_whitelist_counts_gzip (unsigned int *whitelist_total,
			      Tableuint_T whitelist_counts, Tableuint_T whitelistp, 
			      int *nextchar, gzFile *input1, char ***files, int *nfiles) {
  int fulllength;
  unsigned int count;

  while (1) {
    if (*input1 == NULL || *nextchar == EOF) { /* was gzeof(*input1) */
      if (*input1 != NULL) {
	gzclose(*input1);
	*input1 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return false;

      } else if (*nfiles == 1) {
	fprintf(stderr,"Error: odd number of fastq files\n");
	(*files) += 1;
	(*nfiles) -= 1;
	*nextchar = EOF;
	return false;

      } else {
	while (*nfiles > 0 && (*input1 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it and its pair\n",
		  (*files)[0]);
	  (*files) += 2;
	  (*nfiles) -= 2;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return false;
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*input1,GZBUFFER_SIZE);
#endif
	  (*files) += 2;
	  (*nfiles) -= 2;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init_gzip(*input1)) == EOF) {
	return false;
      }
    }

    debug(printf("** Getting header for fastq gzip\n"));
    if (skip_header_gzip(*input1,*nextchar) == false) {
      /* Skip */
      *nextchar = EOF;
    } else {
      *nextchar = gzgetc(*input1);
      if ((fulllength = input_oneline_gzip(&(*nextchar),&(Read1[0]),*input1)) == 0) {
	/* Skip */
      } else {
	if (*nextchar == '+') {
	  /* Ignore quality */
	  skip_header_gzip(*input1,*nextchar);
	  *nextchar = gzgetc(*input1);
	  input_oneline_gzip(&(*nextchar),&(Quality[0]),*input1);
	}

	/* strncpy(barcode,Read1,CELL_BARCODE_LEN); */
	/* barcode[CELL_BARCODE_LEN] = '\0'; */
	if (Tableuint_get(whitelistp,(const void *) Read1) == (unsigned int) true) {
	  count = Tableuint_get(whitelist_counts,(const void *) Read1);
	  Tableuint_put(whitelist_counts,(const void *) Read1,count+1);
	  *whitelist_total += 1;
	}

	return true;
      }
    }
  }
}
#endif

/* Modified from Shortread_read_fastq_bzip2 in shortread.c */
#ifdef HAVE_BZLIB
static bool
update_whitelist_counts_bzip2 (unsigned int *whitelist_total,
			       Tableuint_T whitelist_counts, Tableuint_T whitelistp,
			       int *nextchar, Bzip2_T *input1, char ***files, int *nfiles) {
  int fulllength;
  unsigned int count;

  while (1) {
    if (*input1 == NULL || *nextchar == EOF) { /* Was bzeof(*input1) */
      if (*input1 != NULL) {
	Bzip2_free(&(*input1));
	*input1 = NULL;
      }

      if (*nfiles == 0) {
	*nextchar = EOF;
	return false;

      } else if (*nfiles == 1) {
	fprintf(stderr,"Error: odd number of fastq files\n");
	(*files) += 1;
	(*nfiles) -= 1;
	*nextchar = EOF;
	return false;

      } else {
	while (*nfiles > 0 && (*input1 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Can't open file %s => skipping it and its pair\n",
		  (*files)[0]);
	  (*files) += 2;
	  (*nfiles) -= 2;
	}
	if (*input1 == NULL) {
	  *nextchar = EOF;
	  return false;
	} else {
	  (*files) += 2;
	  (*nfiles) -= 2;
	  *nextchar = '\0';
	}
      }
    }

    if (*nextchar == '\0') {
      if ((*nextchar = Shortread_input_init_bzip2(*input1)) == EOF) {
	return false;
      }
    }

    debug(printf("** Getting header for fastq bzip2\n"));
    if (skip_header_bzip2(*input1,*nextchar) == false) {
      /* Skip */
      *nextchar = EOF;
    } else {
      *nextchar = bzgetc(*input1);
      if ((fulllength = input_oneline_bzip2(&(*nextchar),&(Read1[0]),*input1)) == 0) {
	/* Skip */
      } else {
	if (*nextchar == '+') {
	  /* Ignore quality */
	  skip_header_bzip2(*input1,*nextchar);
	  *nextchar = bzgetc(*input1);
	  input_oneline_bzip2(&(*nextchar),&(Quality[0]),*input1);
	}

	/* Can use Read1 because of string_compare_cell_barcode and string_hash_cell_barcode functions */
	/* strncpy(barcode,Read1,CELL_BARCODE_LEN); */
	/* barcode[CELL_BARCODE_LEN] = '\0'; */
	if (Tableuint_get(whitelistp,(const void *) Read1) == (unsigned int) true) {
	  count = Tableuint_get(whitelist_counts,(const void *) Read1);
	  Tableuint_put(whitelist_counts,(const void *) Read1,count+1);
	  *whitelist_total += 1;
	}

	return true;
      }
    }
  }
}
#endif


static bool
update_whitelist_counts (unsigned int *whitelist_total, Tableuint_T whitelist_counts,
			 Tableuint_T whitelistp, int *nextchar, int *nchars1, FILE **input1,
#ifdef HAVE_ZLIB
			 gzFile *gzipped1,
#endif
#ifdef HAVE_BZLIB
			 Bzip2_T *bzipped1,
#endif
			 char *read_files_command, char ***files, int *nfiles) {

#ifdef HAVE_ZLIB    
  if (*gzipped1 != NULL) {
    return update_whitelist_counts_gzip(&(*whitelist_total),whitelist_counts,whitelistp,
					&(*nextchar),&(*gzipped1),&(*files),&(*nfiles));
  }
#endif

#ifdef HAVE_BZLIB
  if (*bzipped1 != NULL) {
    return update_whitelist_counts_bzip2(&(*whitelist_total),whitelist_counts,whitelistp,
					 &(*nextchar),&(*bzipped1),&(*files),&(*nfiles));
  }
#endif

  return update_whitelist_counts_text(&(*whitelist_total),whitelist_counts,whitelistp,
				      &(*nextchar),&(*nchars1),&(*input1),
				      read_files_command,&(*files),&(*nfiles));
}


/* Compares only the first CELL_BARCODE_LEN characters of each */
static int
string_compare_cell_barcode (const void *x, const void *y) {
  char *a = (char *) x;
  char *b = (char *) y;

  return strncmp(a,b,CELL_BARCODE_LEN);
}

/* This is the X31 hash function */
/* Uses the first CELL_BARCODE_LEN characters */
static unsigned int
string_hash_cell_barcode (const void *x) {
  unsigned int h = 0U;
  const char *p = x;
  int i;
  
  for (i = 0; i < CELL_BARCODE_LEN; i++) {
    h = (h << 5) - h + *p++;
  }
  return h;
}


/* Modified from Inbuffer_new and fill_buffer in inbuffer.c */
void
Single_cell_compute_priors (char *whitelist_file,
			    char *read_files_command, bool gunzip_p, bool bunzip2_p,
			    char **files, int nfiles) {
  FILE *fp, *input1;
#ifdef HAVE_ZLIB
  gzFile gzipped1 = NULL;
#endif
#ifdef HAVE_BZLIB
  Bzip2_T bzipped1 = NULL;
#endif

  char *barcode;
  unsigned int whitelist_total;
  int nchars1 = 0; /* Returned only because MPI master needs it.  Doesn't need to be saved as a field in Inbuffer_T. */
  int nextchar = '\0';

  /* Read in whitelist */
  fp = fopen(whitelist_file,"r");
  whitelistp = Tableuint_new(/*hint*/800000,string_compare_cell_barcode,string_hash_cell_barcode);
  whitelist_counts = Tableuint_new(/*hint*/800000,string_compare_cell_barcode,string_hash_cell_barcode);

  while ((barcode = Getline(fp)) != NULL) {
    Tableuint_put(whitelistp,(const void *) barcode,(unsigned int) true);
    Tableuint_put(whitelist_counts,(const void *) barcode,0);
  }
  fclose(fp);

  /* Taken from open_streams_parser in gsnap.c */
  if (gunzip_p == true) {
#ifdef HAVE_ZLIB
    if ((gzipped1 = gzopen(files[0],"rb")) == NULL) {
      fprintf(stderr,"Cannot open gzipped file %s\n",files[0]);
      exit(9);
    } else {
#ifdef HAVE_ZLIB_GZBUFFER
      gzbuffer(gzipped1,GZBUFFER_SIZE);
#endif
      nextchar = Shortread_input_init_gzip(gzipped1);
    }
#endif

  } else if (bunzip2_p == true) {
    
#ifdef HAVE_BZLIB
    if ((bzipped1 = Bzip2_new(files[0])) == NULL) {
      fprintf(stderr,"Cannot open bzipped file %s\n",files[0]);
      exit(9);
    } else {
      nextchar = Shortread_input_init_bzip2(bzipped1);
    }
#endif

  } else {
    if ((input1 = Fopen_read_text(read_files_command,files[0])) == NULL) {
      fprintf(stderr,"Unable to read input\n");
      exit(9);
    } else {
      nextchar = Shortread_input_init(&nchars1,input1);
    }
  }

  if (nextchar == EOF) {
    fprintf(stderr,"Warning: input is empty\n");
    exit(9);
  } else {
    files += 2;
    nfiles -= 2;
  }

  /* Tally barcodes in whitelist from read1 fastq files */
  whitelist_total = 0;
  while (update_whitelist_counts(&whitelist_total,whitelist_counts,whitelistp,
				 &nextchar,&nchars1,&input1,
#ifdef HAVE_ZLIB
				 &gzipped1,
#endif
#ifdef HAVE_BZLIB
				 &bzipped1,
#endif
				 read_files_command,&files,&nfiles) == true) {
    /* Nothing to do */
  }
  whitelist_total_double = (double) whitelist_total;

  return;
}


/* From http://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam
 *
 *
 * CB: Chromium cellular barcode sequence that is error-corrected and
 * confirmed against a list of known-good barcode sequences.  Includes
 * a suffix with a dash separator followed by a number: GEM well
 *
 * CR: Chromium cellular barcode sequence as reported by the sequencer
 * CY: Chromium cellular barcode read quality.  Phred scores as reported by sequencer
 *
 * WE (GSNAP-specific): the GEM well, needed if we don't have a CB field and need to use CR
 *
 * UB: Chromium molecular barcode sequence that is error-corrected
 * among other molecular barcodes with the same cellular barcode and
 * gene alignment.  Cannot correct while we are aligning since we need
 * all alignments
 *
 * UR: Chromium molecular barcode sequence as reported by the sequencer
 * UY: Chromium molecular barcode read quality.  Phred scores as reported by sequencer
 *
 * TR: Trimmed sequence.  For the single cell 3' v1 chemistry, this is
 * trailing sequence following the UMI on Read 2.  For the single cell
 * 3' v2 chemistry, this is trailing sequence following the cell and
 * molecular barcodes on Read 1
 */

#define CONFIDENCE_THRESHOLD 0.975
#define ILLUMINA_QUAL_OFFSET 33

/* Based on correct_bc_error in cellranger-master/lib/python/cellranger/stats.py */
double
correct_barcode_error (char **new_barcode, char *contents, char *quality) {
  char *barcode, existing, new_char;
  double max_likelihood = 0.0, likelihood_sum = 0.0, likelihood, edit_qv, prior, p_edit;
  int qv;
  int new_pos, pos;
  unsigned int count;

  barcode = (char *) MALLOC((CELL_BARCODE_LEN+1)*sizeof(char));
  strncpy(barcode,contents,CELL_BARCODE_LEN);
  barcode[CELL_BARCODE_LEN] = '\0';

  for (pos = 0; pos < CELL_BARCODE_LEN; pos++) {
    qv = (int) quality[pos] - ILLUMINA_QUAL_OFFSET;
    edit_qv = (qv < 33 ? qv : 33);
    p_edit = pow(10.0,-edit_qv/10.0);

    existing = barcode[pos];
    if (existing != 'A') {
      barcode[pos] = 'A';
      if ((count = Tableuint_get(whitelist_counts,(const void *) barcode)) > 0) {
	prior = (double) count/whitelist_total_double;
	likelihood_sum += likelihood = prior * p_edit;
	if (likelihood > max_likelihood) {
	  new_pos = pos;
	  new_char = 'A';
	  max_likelihood = likelihood;
	}
      }
    }
    if (existing != 'C') {
      barcode[pos] = 'C';
      if ((count = Tableuint_get(whitelist_counts,(const void *) barcode)) > 0) {
	prior = (double) count/whitelist_total_double;
	likelihood_sum += likelihood = prior * p_edit;
	if (likelihood > max_likelihood) {
	  new_pos = pos;
	  new_char = 'C';
	  max_likelihood = likelihood;
	}
      }
    }
    if (existing != 'G') {
      barcode[pos] = 'G';
      if ((count = Tableuint_get(whitelist_counts,(const void *) barcode)) > 0) {
	prior = (double) count/whitelist_total_double;
	likelihood_sum += likelihood = prior * p_edit;
	if (likelihood > max_likelihood) {
	  new_pos = pos;
	  new_char = 'G';
	  max_likelihood = likelihood;
	}
      }
    }
    if (existing != 'T') {
      barcode[pos] = 'T';
      if ((count = Tableuint_get(whitelist_counts,(const void *) barcode)) > 0) {
	prior = (double) count/whitelist_total_double;
	likelihood_sum += likelihood = prior * p_edit;
	if (likelihood > max_likelihood) {
	  new_pos = pos;
	  new_char = 'T';
	  max_likelihood = likelihood;
	}
      }
    }
    barcode[pos] = existing;
  }

  if (max_likelihood > 0.0) {
    barcode[new_pos] = new_char;
  }
  *new_barcode = barcode;
  
  return max_likelihood/likelihood_sum;
}


static char *
get_well_string (Shortread_T infoseq, int wellpos) {
  char *well_string = NULL;
  char *header;
  char *well, *end, c;
  int well_length, i;

  if (wellpos == 0) {
    return (char *) NULL;

  } else {
    /* Try accession */
    well = Shortread_accession(infoseq);
    for (i = 1; i < wellpos; i++) {
      /* Typically, wellpos is 4, for Instrument ID, Run number (numeric), and Flowcell ID */
      while ((c = *well) != '\0' && c != ':') {
	well++;
      }
      if (c == ':') {
	well++;
      }
    }

    if (*well != '\0') {
      /* Found a well */
      end = well;
      while ((c = *end) != '\0' && c != ':') {
	end++;
      }
      if (c == ':') {
	end++;
      }
      if (c != '\0' && ((well_length = (end - 1 - well)/sizeof(char)) == 1)) {
	well_string = (char *) MALLOC((well_length+1)*sizeof(char));
	strncpy(well_string,well,well_length);
	well_string[well_length] = '\0';
      }
    }

    if (well_string == NULL && (header = Shortread_header(infoseq)) != NULL && header[0] != '\0') {
      /* Try rest of header */
      well = header;
      
      for (i = 1; i < wellpos; i++) {
	/* Typically, wellpos is 4, for Instrument ID, Run number (numeric), and Flowcell ID */
	while ((c = *well) != '\0' && c != ':') {
	  well++;
	}
	if (c == ':') {
	  well++;
	}
      }

      if (*well != '\0') {
	/* Found a well */
	end = well;
	while ((c = *end) != '\0' && c != ':') {
	  end++;
	}
	if (c == ':') {
	  end++;
	}
	if (c != '\0' && ((well_length = (end - 1 - well)/sizeof(char)) == 1)) {
	  well_string = (char *) MALLOC((well_length+1)*sizeof(char));
	  strncpy(well_string,well,well_length);
	  well_string[well_length] = '\0';
	}
      }
    }
  }

  return well_string;
}


/* Fields are defined at https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam */
void
Single_cell_print_fields (Filestring_T fp, Shortread_T infoseq) {
  char *new_barcode, *contents, *quality;
  char *well_string;
  /* int barcode_length; */

  /* Index */
  /* BC: Sample index.  Not necessary for downstream analysis.  Also,
     files from SRA do not have this index in the header */
#if 0
  if ((header = Shortread_header(infoseq)) == NULL || header[0] == '\0') {
    /* Don't print BC field */
  } else {
    FPRINTF(fp,"\tBC:Z:%s",header);
  }
  /* Not handling QT yet, which requires input of the I1 file */
#endif

  /* Chromium cellular (GEM) barcode */

  /* CB: Chromium cellular barcode sequence that is error-corrected
     and confirmed against a list of known-good barcode sequences.
     The cell barcode CB tag includes a suffix with a dash separator
     followed by a number, e.g., AAACCCAAGGAGAGTA-1 */

  /* WE (GSNAP-specific): Well */
  /* CR: Chromium cellular barcode sequence as reported by the sequencer. */
  /* CY: Chromium cellular barcode read quality */

  contents = Shortread_fullpointer(infoseq);
  quality = Shortread_quality_string(infoseq);

  if ((well_string = get_well_string(infoseq,wellpos)) != NULL) {
    FPRINTF(fp,"\tWE:i:%s",well_string);
  }

  FPRINTF(fp,"\tCR:Z:%.*s",CELL_BARCODE_LEN,contents);
  if (quality != NULL) {
    FPRINTF(fp,"\tCY:Z:%.*s",CELL_BARCODE_LEN,quality);
  }

  if (whitelistp == NULL) {
    /* No correction is possible */

  } else if (Tableuint_get(whitelistp,(const void *) contents) == (unsigned int) true) {
    /* Can use contents because of string_compare_cell_barcode and string_hash_cell_barcode functions */
    
    /* No need to correct */
    if (well_string != NULL) {
      /* Print with well */
      FPRINTF(fp,"\tCB:Z:%.*s-%s",CELL_BARCODE_LEN,contents,well_string);
    } else {
      /* Print without well */
      FPRINTF(fp,"\tCB:Z:%.*s",CELL_BARCODE_LEN,contents);
    }
      
  } else {
    if (correct_barcode_error(&new_barcode,contents,quality) > CONFIDENCE_THRESHOLD) {
      /* Print corrected barcode */
      if (well_string != NULL) {
	/* Print with well */
	FPRINTF(fp,"\tCB:Z:%s-%s",new_barcode,well_string);
      } else {
	/* Print without well */
	FPRINTF(fp,"\tCB:Z:%s",new_barcode);
      }
    } else {
      /* Print nothing for CB */
    }
    FREE(new_barcode);
  }

  FREE(well_string);


  /* Chromium molecular barcode */
  /* UB: Chromium molecular barcode sequence that is error-corrected
     among other molecular barcodes with the same cellular barcode and
     gene alignment.  (Not able to perform this correction at
     alignment time.) */
  /* UR: Chromium molecular barcode sequence as reported by the
     sequencer. */
  /* UY: Chromium molecular barcode read quality.  Phred scores as
     reported by sequencer. */
  FPRINTF(fp,"\tUR:Z:%s",&(contents[CELL_BARCODE_LEN]));
  if (quality != NULL) {
    FPRINTF(fp,"\tUY:Z:%s",&(quality[CELL_BARCODE_LEN]));
  }

  return;
}


void
Single_cell_setup (int wellpos_in) {
  wellpos = wellpos_in;
  return;
}


void
Single_cell_cleanup () {
  char **barcodes;
  int n, i;

  /* barcodes were allocated by Getline() */
  barcodes = (char **) Tableuint_keys(whitelist_counts,(void *) NULL);
  n = Tableuint_length(whitelist_counts);
  for (i = 0; i < n; i++) {
    FREE(barcodes[i]);
  }
  FREE(barcodes);

  /* Keys are in whitelist_counts */
  Tableuint_free(&whitelist_counts);
  Tableuint_free(&whitelistp);

  return;
}

