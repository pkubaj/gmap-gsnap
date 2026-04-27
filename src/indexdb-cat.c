static char rcsid[] = "$Id: indexdb-cat.c 224777 2021-12-17 22:27:38Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
# define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
# define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef WORDS_BIGENDIAN
#include "bigendian.h"
#else
#include "littleendian.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <libgen.h>		/* For basename */

#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"		/* For Positionsptr_T and Oligospace_T */
#include "filesuffix.h"

#include "assert.h"
#include "bool.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "indexdb.h"
#include "indexdbdef.h"
#include "indexdb-write.h"

#include "bitpack64-write.h"
#include "bitpack64-read.h"	/* for Bitpack64_block_offsets */

#include "datadir.h"
#include "getopt.h"


#define POSITIONS8_HIGH_SHIFT 32
#define POSITIONS8_LOW_MASK 0xFFFFFFFF

#define BLOCKSIZE 64
#define BUFFER_SIZE 1000000

#define MONITOR_INTERVAL 10000000
#define MONITOR_INTERVAL_OLIGO 1000000


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



static char *user_destdir = NULL;
static char *destdir = ".";
static char *fileroot = NULL;
static int compression_type;

static int index1part = 15;
static int required_index1part = 0;
static int index1interval;
static int required_index1interval = 0;

static char *snps_root = NULL;


static struct option long_options[] = {
  /* Input options */
  {"destdir", required_argument, 0, 'D'},	/* user_destdir */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part */
  {"sampling", required_argument, 0, 'q'}, /* required_index1interval */
  {"db", required_argument, 0, 'd'}, /* fileroot */
  {"usesnps", required_argument, 0, 'v'}, /* snps_root */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  fprintf(stdout,"\n");
  fprintf(stdout,"INDEXDB_CAT: Concatenates GMAP index files\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Default gmap directory: %s\n",GMAPDB);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

static void
print_program_usage ();


static Oligospace_T
power (int base, int exponent) {
  Oligospace_T result = 1UL;
  int i;

  for (i = 0; i < exponent; i++) {
    result *= base;
  }
  return result;
}


/*                      87654321 */
#define LOW_TWO_BITS  0x00000003



/* All indexdbs are not huge and positions are 4-byte */
static void
merge_indexdbs_uint4 (char *new_pointers_filename, char *new_offsets_filename,
		      char *new_positions_filename,
		      Indexdb_T *indexdbs, Univcoord_T *genome_univcoord, int ngenomes, 
		      Oligospace_T oligospace) {
  UINT4 *ascending4, npositions, k;
  UINT4 offset, increment;
  UINT4 **block_offsets;

  FILE *positions_fp;
  UINT4 **positions_ptrs, *ptr;
  UINT4 position4;

  Indexdb_T indexdb;
  Univcoord_T univcoord;
  Oligospace_T oligoi;
  int ii, genomei;
#ifdef DEBUG
  char *nt1, *nt2;
#endif

  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  } else {
    /* Allocate space */
    positions_ptrs = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      positions_ptrs[genomei] = &(indexdbs[genomei]->positions[0]);
    }
    ascending4 = MALLOC((oligospace+1)*sizeof(UINT4));

    block_offsets = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      block_offsets[genomei] = (UINT4 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT4));
    }
  }

  fprintf(stderr,"Merging offsets and positions");

  offset = 0;
  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    if (oligoi % MONITOR_INTERVAL_OLIGO == 0) {
      fprintf(stderr,".");
    }

    /* Obtain a block of BLOCKSIZE (64) offsets for each genome */
    for (genomei = 0; genomei < ngenomes; genomei++) {
      indexdb = indexdbs[genomei];
      Bitpack64_block_offsets(block_offsets[genomei],oligoi,
			      indexdb->offsetsmeta,indexdb->offsetsstrm);
    }

    /* Loop through each oligo in the block */
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      ascending4[oligoi+ii] = offset; /* ascending4[0] == 0 */

      increment = 0;
      for (genomei = 0; genomei < ngenomes; genomei++) {
	indexdb = indexdbs[genomei];
	univcoord = genome_univcoord[genomei];
	ptr = positions_ptrs[genomei];
	assert(ptr == &(indexdb->positions[block_offsets[genomei][ii]]));

	/* Determine number of positions for genomei for this oligo, and append them */
	increment += (npositions = block_offsets[genomei][ii+1] - block_offsets[genomei][ii]);
	for (k = 0; k < npositions; k++) {
	  position4 = (UINT4) (univcoord + (*ptr++));
	  FWRITE_UINT(position4,positions_fp);
	}

	positions_ptrs[genomei] = ptr; /* Update for next oligo */
      }
      offset += increment;
    }
  }

  /* Final offset */
  ascending4[oligospace] = offset;
  fclose(positions_fp);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Writing offsets...");
  Bitpack64_write_differential(new_pointers_filename,new_offsets_filename,
			       ascending4,oligospace);
  fprintf(stderr,"done\n");

  for (genomei = 0; genomei < ngenomes; genomei++) {
    FREE(block_offsets[genomei]);
  }
  FREE(block_offsets);
  FREE(ascending4);
  FREE(positions_ptrs);

  return;
}


#ifdef HAVE_64_BIT
/* All indexdbs are not huge and positions are 8-byte */
static void
merge_indexdbs_uint8 (char *new_pointers_filename, char *new_offsets_filename,
		      char *new_positions_high_filename, char *new_positions_filename,
		      Indexdb_T *indexdbs, Univcoord_T *genome_univcoord,
		      bool *genome_8p, int ngenomes, Oligospace_T oligospace) {
  UINT4 *ascending4, npositions, k;
  UINT4 offset, increment;
  UINT4 **block_offsets;

  FILE *positions_high_fp, *positions_fp;
  unsigned char **positions_high_ptrs, *high_ptr, position8_high;
  UINT4 **positions_ptrs, *ptr, position8_low;
  UINT8 position8;

  Indexdb_T indexdb;
  Univcoord_T univcoord;
  Oligospace_T oligoi;
  int ii, genomei;

  if ((positions_high_fp = FOPEN_WRITE_BINARY(new_positions_high_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_high_filename);
    exit(9);
  } else if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  } else {
    /* Allocate space */
    positions_high_ptrs = (unsigned char **) MALLOC(ngenomes*sizeof(unsigned char *));
    positions_ptrs = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      if (genome_8p[genomei] == false) {
	positions_high_ptrs[genomei] = (unsigned char *) NULL;
      } else {
	positions_high_ptrs[genomei] = &(indexdbs[genomei]->positions_high[0]);
      }
      positions_ptrs[genomei] = &(indexdbs[genomei]->positions[0]);
    }
    ascending4 = MALLOC((oligospace+1)*sizeof(UINT4));

    block_offsets = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      block_offsets[genomei] = (UINT4 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT4));
    }
  }

  fprintf(stderr,"Merging offsets and positions");

  offset = 0;
  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    if (oligoi % MONITOR_INTERVAL_OLIGO == 0) {
      fprintf(stderr,".");
    }

    /* Obtain a block of BLOCKSIZE (64) offsets for each genome */
    for (genomei = 0; genomei < ngenomes; genomei++) {
      indexdb = indexdbs[genomei];
      Bitpack64_block_offsets(block_offsets[genomei],oligoi,
			      indexdb->offsetsmeta,indexdb->offsetsstrm);
    }

    /* Loop through each oligo in the block */
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      ascending4[oligoi+ii] = offset; /* ascending4[0] == 0 */

      increment = 0;
      for (genomei = 0; genomei < ngenomes; genomei++) {
	indexdb = indexdbs[genomei];
	univcoord = genome_univcoord[genomei];
	ptr = positions_ptrs[genomei];
	assert(ptr == &(indexdb->positions[block_offsets[genomei][ii]]));

	/* Determine number of positions for genomei for this oligo, and append them */
	increment += (npositions = block_offsets[genomei][ii+1] - block_offsets[genomei][ii]);
	if ((high_ptr = positions_high_ptrs[genomei]) == NULL) {
	  for (k = 0; k < npositions; k++) {
	    position8 = univcoord + (Univcoord_T) (*ptr++);
	    position8_high = (unsigned char) (position8 >> POSITIONS8_HIGH_SHIFT);
	    position8_low = (UINT4) (position8 & POSITIONS8_LOW_MASK);

	    FWRITE_CHAR(position8_high,positions_high_fp);
	    FWRITE_UINT(position8_low,positions_fp);
	  }
	  positions_ptrs[genomei] = ptr; /* Update for next oligo */

	} else {
	  for (k = 0; k < npositions; k++) {
	    position8 = univcoord + ((Univcoord_T) (*high_ptr++) << 32) + (Univcoord_T) (*ptr++);
	    position8_high = (unsigned char) (position8 >> POSITIONS8_HIGH_SHIFT);
	    position8_low = (UINT4) (position8 & POSITIONS8_LOW_MASK);

	    FWRITE_CHAR(position8_high,positions_high_fp);
	    FWRITE_UINT(position8_low,positions_fp);
	  }
	  positions_high_ptrs[genomei] = high_ptr; /* Update for next oligo */
	  positions_ptrs[genomei] = ptr;
	}

      }
      offset += increment;
    }
  }

  /* Final offset */
  ascending4[oligospace] = offset;
  fclose(positions_fp);
  fclose(positions_high_fp);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Writing offsets...");
  Bitpack64_write_differential(new_pointers_filename,new_offsets_filename,
			       ascending4,oligospace);
  fprintf(stderr,"done\n");

  for (genomei = 0; genomei < ngenomes; genomei++) {
    FREE(block_offsets[genomei]);
  }
  FREE(block_offsets);
  FREE(ascending4);
  FREE(positions_ptrs);
  FREE(positions_high_ptrs);

  return;
}
#endif


#ifdef HAVE_64_BIT
/* All or some indexdbs are huge and positions are 4-byte */
static void
merge_indexdbs_huge_uint4 (char *new_pages_filename, char *new_pointers_filename,
			   char *new_offsets_filename, char *new_positions_filename,
			   Indexdb_T *indexdbs, bool *hugep, Univcoord_T *genome_univcoord,
			   int ngenomes, Oligospace_T oligospace) {
  UINT8 *ascending8, npositions, k;
  UINT8 offset, increment;
  UINT8 **block_offsets_uint8;
  UINT4 **block_offsets_uint4;

  FILE *positions_fp;
  UINT4 **positions_ptrs, *ptr;
  UINT4 position4;

  Indexdb_T indexdb;
  Univcoord_T univcoord;
  Oligospace_T oligoi;
  int ii, genomei;

  if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  } else {
    /* Allocate space */
    positions_ptrs = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      positions_ptrs[genomei] = &(indexdbs[genomei]->positions[0]);
    }
    ascending8 = MALLOC((oligospace+1)*sizeof(UINT8));

    block_offsets_uint8 = (UINT8 **) MALLOC(ngenomes*sizeof(UINT8 *));
    block_offsets_uint4 = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      if (hugep[genomei] == true) {
	block_offsets_uint8[genomei] = (UINT8 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT8));
      } else {
	block_offsets_uint4[genomei] = (UINT4 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT4));
      }
    }
  }

  fprintf(stderr,"Merging offsets and positions");

  offset = 0;
  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    if (oligoi % MONITOR_INTERVAL_OLIGO == 0) {
      fprintf(stderr,".");
    }

    /* Obtain a block of BLOCKSIZE (64) offsets for each genome */
    for (genomei = 0; genomei < ngenomes; genomei++) {
      indexdb = indexdbs[genomei];
      if (hugep[genomei] == true) {
	Bitpack64_block_offsets_huge(block_offsets_uint8[genomei],oligoi,indexdb->offsetspages,
				     indexdb->offsetsmeta,indexdb->offsetsstrm);
      } else {
	Bitpack64_block_offsets(block_offsets_uint4[genomei],oligoi,
				indexdb->offsetsmeta,indexdb->offsetsstrm);
      }
    }

    /* Loop through each oligo in the block */
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      ascending8[oligoi+ii] = offset; /* ascending4[0] == 0 */

      increment = 0;
      for (genomei = 0; genomei < ngenomes; genomei++) {
	indexdb = indexdbs[genomei];
	univcoord = genome_univcoord[genomei];
	ptr = positions_ptrs[genomei];

	/* Determine number of positions for genomei for this oligo, and append them */
	if (hugep[genomei] == true) {
	  assert(ptr == &(indexdb->positions[block_offsets_uint8[genomei][ii]]));
	  increment += (npositions = block_offsets_uint8[genomei][ii+1] - block_offsets_uint8[genomei][ii]);
	} else {
	  assert(ptr == &(indexdb->positions[block_offsets_uint4[genomei][ii]]));
	  increment += (npositions = (UINT8) (block_offsets_uint4[genomei][ii+1] - block_offsets_uint4[genomei][ii]));
	}

	for (k = 0; k < npositions; k++) {
	  position4 = (UINT4) (univcoord + (*ptr++));
	  FWRITE_UINT(position4,positions_fp);
	}
	positions_ptrs[genomei] = ptr; /* Update for next oligo */
      }

      offset += increment;
    }
  }

  /* Final offset */
  ascending8[oligospace] = offset;
  fclose(positions_fp);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Writing offsets...");
  Bitpack64_write_differential_huge(new_pages_filename,new_pointers_filename,new_offsets_filename,
				    ascending8,oligospace);
  fprintf(stderr,"done\n");

  for (genomei = 0; genomei < ngenomes; genomei++) {
    if (hugep[genomei] == true) {
      FREE(block_offsets_uint8[genomei]);
    } else {
      FREE(block_offsets_uint4[genomei]);
    }
  }
  FREE(block_offsets_uint4);
  FREE(block_offsets_uint8);
  FREE(ascending8);
  FREE(positions_ptrs);

  return;
}
#endif

#ifdef HAVE_64_BIT
/* All or some indexdbs are huge and positions are 8-byte */
static void
merge_indexdbs_huge_uint8 (char *new_pages_filename, char *new_pointers_filename,
			   char *new_offsets_filename,
			   char *new_positions_high_filename, char *new_positions_filename,
			   Indexdb_T *indexdbs, bool *hugep, Univcoord_T *genome_univcoord,
			   bool *genome_8p, int ngenomes, Oligospace_T oligospace) {
  UINT8 *ascending8, npositions, k;
  UINT8 offset, increment;
  UINT8 **block_offsets_uint8;
  UINT4 **block_offsets_uint4;

  FILE *positions_high_fp, *positions_fp;
  unsigned char **positions_high_ptrs, *high_ptr, position8_high;
  UINT4 **positions_ptrs, *ptr, position8_low;
  UINT8 position8;

  Indexdb_T indexdb;
  Univcoord_T univcoord;
  Oligospace_T oligoi;
  int ii, genomei;

  if ((positions_high_fp = FOPEN_WRITE_BINARY(new_positions_high_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_high_filename);
    exit(9);
  } else if ((positions_fp = FOPEN_WRITE_BINARY(new_positions_filename)) == NULL) {
    fprintf(stderr,"Can't write to file %s\n",new_positions_filename);
    exit(9);
  } else {
    /* Allocate space */
    positions_high_ptrs = (unsigned char **) MALLOC(ngenomes*sizeof(unsigned char *));
    positions_ptrs = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      if (genome_8p[genomei] == false) {
	positions_high_ptrs[genomei] = (unsigned char *) NULL;
      } else {
	positions_high_ptrs[genomei] = &(indexdbs[genomei]->positions_high[0]);
      }
      positions_ptrs[genomei] = &(indexdbs[genomei]->positions[0]);
    }
    ascending8 = MALLOC((oligospace+1)*sizeof(UINT8));

    block_offsets_uint8 = (UINT8 **) MALLOC(ngenomes*sizeof(UINT8 *));
    block_offsets_uint4 = (UINT4 **) MALLOC(ngenomes*sizeof(UINT4 *));
    for (genomei = 0; genomei < ngenomes; genomei++) {
      if (hugep[genomei] == true) {
	block_offsets_uint8[genomei] = (UINT8 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT8));
      } else {
	block_offsets_uint4[genomei] = (UINT4 *) MALLOC((BLOCKSIZE+1)*sizeof(UINT4));
      }
    }
  }

  fprintf(stderr,"Merging offsets and positions");

  offset = 0;
  for (oligoi = 0; oligoi < oligospace; oligoi += BLOCKSIZE) {
    if (oligoi % MONITOR_INTERVAL_OLIGO == 0) {
      fprintf(stderr,".");
    }

    /* Obtain a block of BLOCKSIZE (64) offsets for each genome */
    for (genomei = 0; genomei < ngenomes; genomei++) {
      indexdb = indexdbs[genomei];
      if (hugep[genomei] == true) {
	Bitpack64_block_offsets_huge(block_offsets_uint8[genomei],oligoi,indexdb->offsetspages,
				     indexdb->offsetsmeta,indexdb->offsetsstrm);
      } else {
	Bitpack64_block_offsets(block_offsets_uint4[genomei],oligoi,
				indexdb->offsetsmeta,indexdb->offsetsstrm);
      }
    }

    /* Loop through each oligo in the block */
    for (ii = 0; ii < BLOCKSIZE; ii++) {
      ascending8[oligoi+ii] = offset; /* ascending8[0] == 0 */

      increment = 0;
      for (genomei = 0; genomei < ngenomes; genomei++) {
	indexdb = indexdbs[genomei];
	univcoord = genome_univcoord[genomei];
	ptr = positions_ptrs[genomei];

	/* Determine number of positions for genomei for this oligo, and append them */
	if (hugep[genomei] == true) {
	  assert(ptr == &(indexdb->positions[block_offsets_uint8[genomei][ii]]));
	  increment += (npositions = block_offsets_uint8[genomei][ii+1] - block_offsets_uint8[genomei][ii]);
	} else {
	  assert(ptr == &(indexdb->positions[block_offsets_uint4[genomei][ii]]));
	  increment += (npositions = (UINT8) (block_offsets_uint4[genomei][ii+1] - block_offsets_uint4[genomei][ii]));
	}

	if ((high_ptr = positions_high_ptrs[genomei]) == NULL) {
	  for (k = 0; k < npositions; k++) {
	    position8 = univcoord + (Univcoord_T) (*ptr++);
	    position8_high = (unsigned char) (position8 >> POSITIONS8_HIGH_SHIFT);
	    position8_low = (UINT4) (position8 & POSITIONS8_LOW_MASK);

	    FWRITE_CHAR(position8_high,positions_high_fp);
	    FWRITE_UINT(position8_low,positions_fp);
	  }
	  positions_ptrs[genomei] = ptr; /* Update for next oligo */

	} else {
	  for (k = 0; k < npositions; k++) {
	    position8 = univcoord + ((Univcoord_T) (*high_ptr++) << 32) + (Univcoord_T) (*ptr++);
	    position8_high = (unsigned char) (position8 >> POSITIONS8_HIGH_SHIFT);
	    position8_low = (UINT4) (position8 & POSITIONS8_LOW_MASK);

	    FWRITE_CHAR(position8_high,positions_high_fp);
	    FWRITE_UINT(position8_low,positions_fp);
	  }
	  positions_high_ptrs[genomei] = high_ptr; /* Update for next oligo */
	  positions_ptrs[genomei] = ptr;
	}
      }

      offset += increment;
    }
  }

  /* Final offset */
  ascending8[oligospace] = offset;
  fclose(positions_fp);
  fclose(positions_high_fp);
  fprintf(stderr,"done\n");

  fprintf(stderr,"Writing offsets...");
  Bitpack64_write_differential_huge(new_pages_filename,new_pointers_filename,new_offsets_filename,
				    ascending8,oligospace);
  fprintf(stderr,"done\n");

  for (genomei = 0; genomei < ngenomes; genomei++) {
    if (hugep[genomei] == true) {
      FREE(block_offsets_uint8[genomei]);
    } else {
      FREE(block_offsets_uint4[genomei]);
    }
  }
  FREE(block_offsets_uint4);
  FREE(block_offsets_uint8);
  FREE(ascending8);
  FREE(positions_ptrs);
  FREE(positions_high_ptrs);

  return;
}
#endif


/* Usage: indexdb_cat -d <genome> <input_genomes...> */

/* Note: Depends on having gmapindex sampled on mod 3 bounds */
int
main (int argc, char *argv[]) {
  Indexdb_filenames_T ifilenames;
  char *source_genome, *source_dbroot, *filename;
  char *new_pages_filename = NULL, *new_pointers_filename = NULL, *new_offsets_filename = NULL,
    *new_positions_high_filename = NULL, *new_positions_filename = NULL;
  char interval_char;

  int ngenomes, genomei;
  Univ_IIT_T chromosome_iit;
  Univcoord_T *genome_univcoord, total_genomelength, total_noffsets;
  Oligospace_T oligospace;
  bool *hugep, huge_offsets_p;
  bool *genome_8p, coord_values_8p;

  Indexdb_T *indexdbs;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;

  while ((opt = getopt_long(argc,argv,"F:D:d:k:q:v:",
			    long_options,&long_option_index)) != -1) {
    switch (opt) {
    case 0: 
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	exit(0);
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'cmetindex --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_destdir = optarg; break;
    case 'd': fileroot = optarg; break;
    case 'k': required_index1part = atoi(optarg); break;
    case 'q': required_index1interval = atoi(optarg); break;
    case 'v': 
      snps_root = optarg;
      fprintf(stderr,"Combination of cmetindex and snps is not yet supported in 2018 versions\n");
      exit(9);
      break;
    default: fprintf(stderr,"Do not recognize flag %c\n",opt); exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  if (fileroot == NULL) {
    fprintf(stderr,"Missing name of destination genome database.  Must specify with -d flag.\n");
    fprintf(stderr,"Usage: indexdb-cat -d <genome> <input_genomes...>\n");
    exit(9);
  }

  destdir = Datadir_find_genomedir(user_destdir);
  ngenomes = argc;
  fprintf(stderr,"Concatenating %d genomes\n",ngenomes);

  /* Chromosome IIT files.  Need to determine coord_values_8p and genomelength */
  genome_8p = (bool *) MALLOC(ngenomes*sizeof(bool));
  genome_univcoord = (Univcoord_T *) MALLOC(ngenomes*sizeof(Univcoord_T));
  total_genomelength = 0;

  for (genomei = 0; genomei < ngenomes; genomei++) {
    genome_univcoord[genomei] = total_genomelength;
    source_genome = MALLOC((strlen(argv[genomei])+1)*sizeof(char));
    strcpy(source_genome,argv[genomei]);
    source_dbroot = basename(source_genome);

    filename = (char *) CALLOC(strlen(argv[genomei])+strlen("/")+
			       strlen(source_dbroot)+strlen(".chromosome.iit")+1,sizeof(char));
    sprintf(filename,"%s/%s.chromosome.iit",argv[genomei],source_dbroot);
    if ((chromosome_iit = Univ_IIT_read(filename,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      fprintf(stderr,"Could not find chromosome IIT file for the input genome %s\n",argv[genomei]);
      fprintf(stderr,"Expecting to find %s\n",filename);
      exit(9);
    } else {
      FREE(filename);
      FREE(source_genome);
      genome_8p[genomei] = Univ_IIT_coord_values_8p(chromosome_iit);
      total_genomelength += Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
      Univ_IIT_free(&chromosome_iit);
    }
  }

#ifdef HAVE_64_BIT
  if (total_genomelength > 4294967295) {
    coord_values_8p = true;
  } else {
    coord_values_8p = false;
  }
#else
  coord_values_8p = false;
#endif

  indexdbs = (Indexdb_T *) MALLOC(ngenomes*sizeof(Indexdb_T));

  for (genomei = 0; genomei < ngenomes; genomei++) {
    fprintf(stderr,"Reading genome %s\n",argv[genomei]);
    source_genome = MALLOC((strlen(argv[genomei])+1)*sizeof(char));
    strcpy(source_genome,argv[genomei]);
    source_dbroot = basename(source_genome);

    if ((ifilenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
					    /*genomesubdir*/source_genome,/*fileroot*/source_dbroot,
					    IDX_FILESUFFIX,snps_root,/*required_index1part*/0,/*required_index1interval*/0,
					    /*offsets_only_p*/false)) == NULL) {
      fprintf(stderr,"Could not find indexdb file for the input genome %s\n",argv[genomei]);
      exit(9);
    } else {
      if (required_index1part == 0) {
	required_index1part = index1part;
      }
      if (required_index1interval == 0) {
	required_index1interval = index1interval;
      } else if (index1interval > required_index1interval) {
	/* Allows q3 and q1 indices to be integrated */
	required_index1interval = index1interval;
      }
    }
  }

  fprintf(stderr,"Integrating at kmer %d and sampling interval %d\n",required_index1part,required_index1interval);

  for (genomei = 0; genomei < ngenomes; genomei++) {
    fprintf(stderr,"Reading genome %s\n",argv[genomei]);
    source_genome = MALLOC((strlen(argv[genomei])+1)*sizeof(char));
    strcpy(source_genome,argv[genomei]);
    source_dbroot = basename(source_genome);

    if ((ifilenames = Indexdb_get_filenames(&compression_type,&index1part,&index1interval,
					    /*genomesubdir*/source_genome,/*fileroot*/source_dbroot,
					    IDX_FILESUFFIX,snps_root,/*required_index1part*/0,/*required_index1interval*/0,
					    /*offsets_only_p*/false)) == NULL ||
	index1part != required_index1part || index1interval > required_index1interval) {
      fprintf(stderr,"Could not find indexdb file at kmer %d and sampling interval %d for %s\n",
	      required_index1part,required_index1interval,argv[genomei]);
      exit(9);
    } else if ((indexdbs[genomei] = Indexdb_new_genome(&index1part,&index1interval,
						       /*genomesubdir*/source_genome,/*snpsdir*/NULL,
						       /*fileroot*/source_dbroot,IDX_FILESUFFIX,/*snps_root*/NULL,
						       /*required_index1part*/0,/*required_index1interval*/0,
						       /*offsetsstrm_access*/USE_MMAP_ONLY,
						       /*positions_access*/USE_MMAP_ONLY,/*sharedp*/false,
						       /*multiple_sequences_p*/false,/*preload_shared_memory_p*/false,
						       /*unload_shared_memory_p*/false)) == NULL ||
	index1part != required_index1part || index1interval > required_index1interval) {
      fprintf(stderr,"Could not open genome file kmer %d and sampling interval %d for %s\n",
	      required_index1part,required_index1interval,argv[genomei]);
      exit(9);
    } else {
      Indexdb_filenames_free(&ifilenames);
      FREE(source_genome);
    }
  }


  index1part = required_index1part;
  index1interval = required_index1interval;

  Indexdb_setup(index1part);
  oligospace = power(4,index1part);
  if (index1interval == 3) {
    interval_char = '3';
  } else if (index1interval == 2) {
    interval_char = '2';
  } else if (index1interval == 1) {
    interval_char = '1';
  } else {
    fprintf(stderr,"Selected indexing interval %d is not allowed.  Only values allowed are 3, 2, or 1\n",
	    index1interval);
    exit(9);
  }

  total_noffsets = 0;
  hugep = (bool *) MALLOC(ngenomes*sizeof(bool));
  for (genomei = 0; genomei < ngenomes; genomei++) {
    if (indexdbs[genomei]->hugep == true) {
      hugep[genomei] = true;
      total_noffsets += Bitpack64_read_one_huge(oligospace,indexdbs[genomei]->offsetspages,
						indexdbs[genomei]->offsetsmeta,indexdbs[genomei]->offsetsstrm);
    } else {
      hugep[genomei] = false;
      total_noffsets += (Univcoord_T) Bitpack64_read_one(oligospace,indexdbs[genomei]->offsetsmeta,
							 indexdbs[genomei]->offsetsstrm);
    }
  }
    

  if (total_noffsets <= 4294967295) {
    fprintf(stderr,"Total number of offsets: %llu => pages file not required\n",total_noffsets);
    huge_offsets_p = false;
    new_pages_filename = (char *) NULL;
  } else {
    fprintf(stderr,"Total number of offsets: %llu => pages file required\n",total_noffsets);
    huge_offsets_p = true;
    new_pages_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					 strlen(".")+strlen(IDX_FILESUFFIX)+
					 /*for kmer*/2+/*for interval char*/1+
					 strlen("offsets64pages")+1,sizeof(char));
    sprintf(new_pages_filename,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,"offsets64pages");
  }

  new_pointers_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					  strlen(".")+strlen(IDX_FILESUFFIX)+
					  /*for kmer*/2+/*for interval char*/1+
					  strlen("offsets64meta")+1,sizeof(char));
  sprintf(new_pointers_filename,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,"offsets64meta");

  new_offsets_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					 strlen(".")+strlen(IDX_FILESUFFIX)+
					 /*for kmer*/2+/*for interval char*/1+
					 strlen("offsets64strm")+1,sizeof(char));
  sprintf(new_offsets_filename,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,"offsets64strm");

  if (coord_values_8p == true) {
    new_positions_high_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
						  strlen(".")+strlen(IDX_FILESUFFIX)+
						  /*for kmer*/2+/*for interval char*/1+
						  strlen(POSITIONS_HIGH_FILESUFFIX)+1,sizeof(char));
    sprintf(new_positions_high_filename,"%s/%s.%s%02d%c%s",
	    destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_HIGH_FILESUFFIX);
  }

  new_positions_filename = (char *) CALLOC(strlen(destdir)+strlen("/")+strlen(fileroot)+
					   strlen(".")+strlen(IDX_FILESUFFIX)+
					   /*for kmer*/2+/*for interval char*/1+
					   strlen(POSITIONS_FILESUFFIX)+1,sizeof(char));
  sprintf(new_positions_filename,"%s/%s.%s%02d%c%s",
	  destdir,fileroot,IDX_FILESUFFIX,index1part,interval_char,POSITIONS_FILESUFFIX);

#ifdef HAVE_64_BIT
  if (huge_offsets_p == true && coord_values_8p == true) {
    fprintf(stderr,"Genomic coords are 8-byte\n");
    merge_indexdbs_huge_uint8(new_pages_filename,new_pointers_filename,new_offsets_filename,
			      new_positions_high_filename,new_positions_filename,
			      indexdbs,hugep,genome_univcoord,genome_8p,ngenomes,
			      oligospace);

  } else if (huge_offsets_p == true) {
    fprintf(stderr,"Genomic coords are 4-byte\n");
    merge_indexdbs_huge_uint4(new_pages_filename,new_pointers_filename,new_offsets_filename,
			      new_positions_filename,indexdbs,hugep,genome_univcoord,ngenomes,
			      oligospace);

  } else if (coord_values_8p == true) {
    fprintf(stderr,"Genomic coords are 8-byte\n");
    merge_indexdbs_uint8(new_pointers_filename,new_offsets_filename,
			 new_positions_high_filename,new_positions_filename,
			 indexdbs,genome_univcoord,genome_8p,ngenomes,
			 oligospace);

  } else {
    fprintf(stderr,"Genomic coords are 4-byte\n");
    merge_indexdbs_uint4(new_pointers_filename,new_offsets_filename,
			 new_positions_filename,indexdbs,genome_univcoord,ngenomes,
			 oligospace);
  }
#else
  merge_indexdbs_uint4(new_pointers_filename,new_offsets_filename,
		       new_positions_filename,indexdbs,genome_univcoord,ngenomes,
		       oligospace);
#endif

  /* Print values similarly to "gmapindex -N", so gmap_cat knows whether the result will be large */
  printf("%llu %llu\n",total_noffsets,total_genomelength);


  FREE(genome_8p);
  FREE(genome_univcoord);
  FREE(hugep);
  if (huge_offsets_p == true) {
    FREE(new_pages_filename);
  }
  if (coord_values_8p == true) {
    FREE(new_positions_high_filename);
  }
  FREE(new_positions_filename);
  FREE(new_offsets_filename);
  FREE(new_pointers_filename);

  for (genomei = 0; genomei < ngenomes; genomei++) {
    Indexdb_free(&(indexdbs[genomei]));
  }
  FREE(indexdbs);

  return 0;
}



static void
print_program_usage () {
  fprintf(stdout,"\
Usage: indexdb_cat [OPTIONS...] -d </path/to/genome>\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Options (must include -d)\n");
  fprintf(stdout,"\
  -D, --destdir=directory        Directory where to write cmet index files (default is\n\
                                   value of -F, if provided; otherwise the value of the\n\
                                   GMAP genome directory specified at compile time)\n\
  -d, --db=STRING                Genome database\n\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less).\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  -q, --sampling=INT             Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected basesize and k-mer size\n\
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

