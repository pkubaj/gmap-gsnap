static char rcsid[] = "$Id: 61fe98feabdece46f142f4996992017da56f8e6c $";
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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>		/* For open */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>		/* For memset */
#include <ctype.h>		/* For toupper */
#include <sys/mman.h>		/* For munmap */
#include <math.h>		/* For qsort */
#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#include "assert.h"
#include "mem.h"
#include "fopen.h"
#include "access.h"
#include "types.h"
#include "iit-read.h"
#include "iit-read-univ.h"
#include "bitpack64-write.h"
#include "transcriptome.h"

#include "datadir.h"
#include "getopt.h"


/* Creates the following files:

   Note: <genomesubdir> is <gmapdb>/<genomename>
   and <transcriptomesubdir> is <gmapdb>/<transcriptomename>
   <info_root> is <genomesubdir>/<genomename.transcripts>/<transcriptomename>

   The following file (known as alignment_iit, genes_iit, or transcript_map_iit) is in genome order:
   <gmapdb>/<genomename>/<genomename>.transcripts/<transcriptomename>.genes.iit
   (Transcripts are genomic order, with genomic coordinates.  Needed
   by GSNAP for remapping, and by gcount programs for the gene symbol)

   The following are in transcriptome order:
   <gmapdb>/<genomename>/<genomename>.transcripts/<transcriptomename>.exoninfo
   <gmapdb>/<genomename>/<genomename>.transcripts/<transcriptomename>.chrnums
   <gmapdb>/<genomename>/<genomename>.transcripts/<transcriptomename>.offsets64meta
   <gmapdb>/<genomename>/<genomename>.transcripts/<transcriptomename>.offsets64strm
*/

static struct option long_options[] = {
  /* Input options */
  {"dir", required_argument, 0, 'D'}, /* user_gmapdb */
  {"genomedb", required_argument, 0, 'd'}, /* genome_dbroot */
  {"transcriptomedb", required_argument, 0, 'c'}, /* transcriptome_dbroot */

  /* Help options */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};

static void print_program_usage ();

static char *user_gmapdb = NULL;
static char *genome_dbroot = NULL;
static char *genome_fileroot = NULL;
static char *genomesubdir = NULL;

static char *transcriptome_dbroot = NULL;
static char *transcriptome_fileroot = NULL;
static char *transcriptomesubdir = NULL;


#define BUF_SIZE 1024

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


static void
check_data (IIT_T transcript_map_iit, Transcriptome_T transcriptome) {
  Trnum_T trnum;
  int map_index;
  int *exonbounds;
  Chrpos_T *exonstarts;
  int nexons;
  int transcript_genestrand;
#ifdef DEBUG
  bool alloc_label_p;
#endif

  for (map_index = 1; map_index <= IIT_total_nintervals(transcript_map_iit); map_index++) {
    if ((trnum = Transcriptome_trnum(&nexons,&exonbounds,&exonstarts,transcriptome,map_index)) >= 1) {
      Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      debug(printf("%d: %s\n",map_index,IIT_label(transcript_map_iit,index,&alloc_label_p)));

      if (transcript_genestrand > 0) {
	assert(IIT_interval_low(transcript_map_iit,map_index) == exonstarts[0]);
      } else {
	assert(IIT_interval_high(transcript_map_iit,map_index) == exonstarts[0]);
      }
    }
  }

  return;
}


int
main (int argc, char *argv[]) {
  char *dbversion;
  FILE *fp;
  char *pointersfile, *offsetsfile, *filename, *iitfile;
  char *divstring;
  int *divints, divint;
  unsigned int *offsets;
  
  char *info_root;
  char *alignment_file;
  IIT_T alignment_iit;		/* Same as transcript_map_iit */
  int in, out;
  char buf[BUF_SIZE];
  int nread;

  char *transcript_label;
  Univ_IIT_T chromosome_iit, transcript_iit;
  int ntranscripts, transcripti, nalignments, alignment_index;
  int *transcriptomedb_index;
  int transcript_genestrand;
  int nexons;
  int *exonbounds;
  Chrpos_T *exonstarts;
  bool allocp;

  Transcriptome_T transcriptome;

  int opt;
  extern int optind;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;


  while ((opt = getopt_long(argc,argv,"D:d:c:",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;

      if (!strcmp(long_name,"help")) {
	print_program_usage();
	exit(0);
	
      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'trindex --help'",long_name);
	exit(9);
      }
      break;

    case 'D': user_gmapdb = optarg; break;
    case 'd': genome_dbroot = optarg; break;
    case 'c': transcriptome_dbroot = optarg; break;

    case '?': fprintf(stderr,"For usage, run 'trindex --help'\n"); exit(9);
    default: exit(9);
    }
  }
  argc -= optind;
  argv += optind;

  
  if (genome_dbroot == NULL) {
    fprintf(stderr,"Need to specify the genome used for alignment with the -d flag\n");
    exit(9);
  } else {
    genomesubdir = Datadir_find_genomesubdir(&genome_fileroot,&dbversion,user_gmapdb,
					     genome_dbroot);
    FREE(dbversion);
  }

  if (transcriptome_dbroot == NULL) {
    fprintf(stderr,"Need to specify the transcriptome index name with the -c flag\n");
    exit(9);
  } else {
    transcriptomesubdir = Datadir_find_genomesubdir(&transcriptome_fileroot,&dbversion,user_gmapdb,
						    transcriptome_dbroot);
    FREE(dbversion);
  }

  info_root = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(genome_fileroot)+strlen(".transcripts/")+strlen(transcriptome_fileroot)+1,
			      sizeof(char));
  sprintf(info_root,"%s/%s.transcripts/%s",genomesubdir,genome_fileroot,transcriptome_fileroot);

  if (argc <= 0) {
    fprintf(stderr,"Need to specify an alignment IIT file\n");
    exit(9);
  } else {
    alignment_file = argv[0];
  }

  /* (1) Copy alignment (genes) IIT file */
  if ((alignment_iit = IIT_read(alignment_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true)) == NULL) {
    fprintf(stderr,"Unable to read %s as an IIT file\n",alignment_file);
    exit(9);
  } else {
    /* Copy alignment file back out */
    in = open(alignment_file, O_RDONLY);

    filename = (char *) CALLOC(strlen(info_root)+strlen(".genes.iit")+1,sizeof(char));
    sprintf(filename,"%s.genes.iit",info_root);
    out = open(filename,/*openFlags*/O_CREAT | O_WRONLY | O_TRUNC,
	       /*filePerms*/S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);

    while ((nread = read(in,buf,BUF_SIZE)) > 0) {
      if (write(out,buf,nread) != nread) {
	fprintf(stderr,"Could not write whole buffer to %s\n",filename);
	abort();
      }
    }

    close(out);
    close(in);
    FREE(filename);
  }


  /* Get chromosome_iit and transcript_iit from the genome and transcriptome, respectively */
  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(genome_fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,genome_fileroot);
  if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"Genome IIT file %s is not valid\n",iitfile);
    exit(9);
  } else {
    FREE(iitfile);
  }

  iitfile = (char *) CALLOC(strlen(transcriptomesubdir)+strlen("/")+
			    strlen(transcriptome_fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",transcriptomesubdir,transcriptome_fileroot);
  if ((transcript_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"Transcriptome IIT file %s is not valid\n",iitfile);
    exit(9);
  } else {
    FREE(iitfile);
    ntranscripts = Univ_IIT_total_nintervals(transcript_iit);
  }


  /* (2) Write exoninfo and compute offsets */
  /* Need to follow the transcript numbering in the existing
     chromosome_iit file from the transcriptome db, not in
     alignment_iit (or genes_iit) */
  nalignments = IIT_total_nintervals(alignment_iit);
  transcriptomedb_index = (int *) MALLOC((nalignments+1)*sizeof(int));
  for (alignment_index = 0; alignment_index <= nalignments; alignment_index++) {
    transcriptomedb_index[alignment_index] = -1;
  }

  divints = MALLOC(ntranscripts*sizeof(int));
  offsets = MALLOC((ntranscripts+1)*sizeof(unsigned int));
  offsets[0] = 0;

  filename = (char *) CALLOC(strlen(info_root)+strlen(".exoninfo")+1,sizeof(char));
  sprintf(filename,"%s.exoninfo",info_root);
  fp = fopen(filename,"wb");
  FREE(filename);

  for (transcripti = 1; transcripti <= ntranscripts; transcripti++) {
    transcript_label = Univ_IIT_label(transcript_iit,transcripti,&allocp);
    alignment_index = IIT_find_one(alignment_iit,transcript_label); /* See how transcript aligned to genome */
    /* printf("transcripti %d, label %s, alignment_index %d\n",
       transcripti,transcript_label,alignment_index); */

    if (alignment_index <= 0) {
      offsets[transcripti] = offsets[transcripti-1];
      divints[transcripti-1] = 0;
      fprintf(stderr,"No alignment found for transcript %s\n",transcript_label);

    } else {
      transcriptomedb_index[alignment_index] = transcripti;
      nexons = IIT_gene_exons_array(&transcript_genestrand,&exonbounds,&exonstarts,
				    alignment_iit,alignment_index);
      offsets[transcripti] = offsets[transcripti-1] + (unsigned int) nexons;
      /* printf("%d %s %d %d\n",transcripti,transcript_label,alignment_index,nexons); */

      FWRITE_INTS(exonbounds,nexons,fp);
      FWRITE_UINTS(exonstarts,nexons,fp);
      FREE(exonstarts);
      FREE(exonbounds);

      divstring = IIT_divstring_from_index(alignment_iit,alignment_index); /* Get chromosome name */
      divint = Univ_IIT_find_one(chromosome_iit,divstring); /* And look it up in the genome chromosome_iit file */
      if (transcript_genestrand > 0) {
	divints[transcripti-1] = divint;
      } else {
	divints[transcripti-1] = -divint;
      }
    }

    if (allocp == true) {
      FREE(transcript_label);
    }
  }
  fclose(fp);

  filename = (char *) CALLOC(strlen(info_root)+strlen(".dbindex")+1,sizeof(char));
  sprintf(filename,"%s.dbindex",info_root);
  fp = fopen(filename,"wb");
  FREE(filename);

  FWRITE_INTS(transcriptomedb_index,nalignments+1,fp);
  fclose(fp);
  FREE(transcriptomedb_index);


  /* (3) Write chrnums */
  filename = (char *) CALLOC(strlen(info_root)+strlen(".chrnums")+1,sizeof(char));
  sprintf(filename,"%s.chrnums",info_root);
  fp = fopen(filename,"wb");
  FREE(filename);

  FWRITE_INTS(divints,ntranscripts,fp);
  fclose(fp);
  FREE(divints);


  /* (4) Write offsets */
  pointersfile = (char *) CALLOC(strlen(info_root)+strlen(".offsets64meta")+1,sizeof(char));
  sprintf(pointersfile,"%s.offsets64meta",info_root);

  offsetsfile = (char *) CALLOC(strlen(info_root)+strlen(".offsets64strm")+1,sizeof(char));
  sprintf(offsetsfile,"%s.offsets64strm",info_root);
  
  fprintf(stderr,"Writing %d offsets compressed via bitpack64...",ntranscripts+1);
  Bitpack64_write_differential(pointersfile,offsetsfile,offsets,ntranscripts);
  fprintf(stderr,"done\n");
  FREE(offsetsfile);
  FREE(pointersfile);

  FREE(offsets);

  fprintf(stderr,"Wrote transcriptome files to %s/%s.transcripts/%s.*\n",
	  genomesubdir,genome_fileroot,transcriptome_fileroot);

  transcriptome = Transcriptome_new(genomesubdir,genome_fileroot,transcriptome_fileroot,
				    /*sharedp*/false);
  check_data(alignment_iit,transcriptome);
  Transcriptome_free(&transcriptome);

  IIT_free(&alignment_iit);
  FREE(info_root);
  FREE(transcriptome_fileroot);
  FREE(genomesubdir);

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: trindex [OPTIONS...] -D <gmapdb> -d <genome> -c <transcriptome> <alignment IIT file>\n\
\n\
Utility file used by gmap_build for linking transcriptome alignments to a genome\n\
\n\
");

  return;
}

/* 
trindex requires an alignment IIT file, which gives the transcript
coordinates for some genome index <genome>.  You can use the utility
programs provided in the GMAP package: ensembl_genes, gff3_genes,
gtf_genes, or psl_genes to generate this file, or you can align
transcripts using GMAP --format=map_exons.  You then cat this file to
iit_store -o <filename> to generate the genes IIT file

In the gmap_build program, we make this easy for the user by using
GMAP alignment.

You also need to create a GMAP index for the transcripts in the
alignment IIT file.  If you already have the transcripts in a FASTA
file, you can do

      gmap_build -d <transcriptome> -q 1 <fasta file>

Once you have built the GMAP index for the transcriptome, you can run
this program by doing

      trindex -d <genome> -c <transcriptome> <alignment IIT file>

Finally, you can then perform transcriptome-guided genomic alignment by doing

      gsnap -d <genome> -c <transcriptome> <reads>
*/
