static char rcsid[] = "$Id: samheader.c 224951 2022-10-08 20:03:19Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "samheader.h"

#include <stdlib.h>
#include <string.h>
#include "mem.h"


#define CHUNK 1024

static int argc;
static char **argv;
static int optind_save;
static int nworkers;
static bool orderedp;
static Univ_IIT_T chromosome_iit;

#if defined(GEXACT)
static Outputtype_T output_type = SAM_OUTPUT;
#elif defined(GSNAP)
static Outputtype_T output_type;
#else
static Printtype_T printtype;
static Genome_T global_genome;
#endif

static bool sam_headers_p;
static char *sam_read_group_id;
static char *sam_read_group_name;
static char *sam_read_group_library;
static char *sam_read_group_platform;



void
SAM_header_setup (int argc_in, char **argv_in, int optind_in,
		  int nworkers_in, bool orderedp_in,
		  Univ_IIT_T chromosome_iit_in,
#if defined(GEXACT)
#elif defined(GSNAP)
		  Outputtype_T output_type_in,
#else
		  Printtype_T printtype_in, Genome_T global_genome_in,
#endif
		  bool sam_headers_p_in, char *sam_read_group_id_in, char *sam_read_group_name_in,
		  char *sam_read_group_library_in, char *sam_read_group_platform_in) {

  argc = argc_in;
  argv = argv_in;
  optind_save = optind_in;
  nworkers = nworkers_in;
  orderedp = orderedp_in;

  chromosome_iit = chromosome_iit_in;

#if defined(GEXACT)
  output_type = SAM_OUTPUT;
#elif defined(GSNAP)
  output_type = output_type_in;
#else
  printtype = printtype_in;
  global_genome = global_genome_in;
#endif

  sam_headers_p = sam_headers_p_in;
  sam_read_group_id = sam_read_group_id_in;
  sam_read_group_name = sam_read_group_name_in;
  sam_read_group_library = sam_read_group_library_in;
  sam_read_group_platform = sam_read_group_platform_in;

  return;
}



FILE *
SAM_header_fopen (SAM_split_output_type split_output, bool split_simple_p,
		  char *fileroot, bool paired_end_p, bool appendp) {
  FILE *output;
  char *write_mode;
  char *filename, *suffix;

  if (split_output == OUTPUT_FILE) {
    if (appendp == true) {
      write_mode = "a";
    } else {
      write_mode = "w";
    }

    /* fileroot here is output_file */
    filename = (char *) CALLOC(strlen(fileroot)+1,sizeof(char));
    sprintf(filename,"%s",fileroot);

    if ((output = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
      return (FILE *) NULL;
    } else {
      FREE(filename);
      return output;
    }

  } else if (split_output == OUTPUT_NONE) {
    /* Don't open a file */
    return (FILE *) NULL;

  } else {
    /* split_output */
    if (split_simple_p == true) {
      switch (split_output) {
      case OUTPUT_FILE: /* Handled above */ abort(); break;
      case OUTPUT_NONE: /* Handled above */ abort(); break;
	
      case OUTPUT_OTHER: suffix = "other"; break;
      case OUTPUT_CONCORDANT_FULL: suffix = "concordant_full"; break;
	
      default:
	fprintf(stderr,"Cannot handle split output type %d\n",split_output);
	abort();
      }
      
    } else {
      switch (split_output) {
      case OUTPUT_FILE: /* Handled above */ abort(); break;
      case OUTPUT_NONE: /* Handled above */ abort(); break;

      case OUTPUT_NM: suffix = "nomapping"; break;

      case OUTPUT_UU: suffix = (paired_end_p ? "unpaired_uniq" : "uniq"); break;
      case OUTPUT_UT: suffix = (paired_end_p ? "unpaired_transloc" : "transloc"); break;
      case OUTPUT_UM: suffix = (paired_end_p ? "unpaired_mult" : "mult"); break;

      case OUTPUT_UC: suffix = (paired_end_p ? "unpaired_circular" : "circular"); break;
      case OUTPUT_UX: suffix = (paired_end_p ? "unpaired_mult_xs" : "mult_xs"); break;

      case OUTPUT_HU: suffix = "halfmapping_uniq"; break;
      case OUTPUT_HT: suffix = "halfmapping_transloc"; break;
      case OUTPUT_HM: suffix = "halfmapping_mult"; break;

      case OUTPUT_PI: suffix = "paired_uniq_inv"; break;
      case OUTPUT_PS: suffix = "paired_uniq_scr"; break;
      case OUTPUT_PL: suffix = "paired_uniq_long"; break;
      case OUTPUT_PM: suffix = "paired_mult"; break;

      case OUTPUT_CU: suffix = "concordant_uniq"; break;
      case OUTPUT_CT: suffix = "concordant_transloc"; break;
      case OUTPUT_CM: suffix = "concordant_mult"; break;

      case OUTPUT_HC: suffix = "halfmapping_circular"; break;
      case OUTPUT_PC: suffix = "paired_uniq_circular"; break;
      case OUTPUT_CC: suffix = "concordant_circular"; break;

      case OUTPUT_HX: suffix = "halfmapping_mult_xs"; break;
      case OUTPUT_PX: suffix = "paired_mult_xs"; break;
      case OUTPUT_CX: suffix = "concordant_mult_xs"; break;

      default:
	fprintf(stderr,"Cannot handle split output type %d\n",split_output);
	abort();
      }
    }

    /* fileroot here is split_output_root */
    filename = (char *) CALLOC(strlen(fileroot)+strlen(".")+strlen(suffix)+1,sizeof(char));
    sprintf(filename,"%s.%s",fileroot,suffix);

    if (appendp == true) {
      write_mode = "a";
    } else {
      write_mode = "w";
    }

    if ((output = fopen(filename,write_mode)) == NULL) {
      fprintf(stderr,"Cannot open file %s for writing\n",filename);
      exit(9);
      return (FILE *) NULL;
    } else {
      FREE(filename);
      return output;
    }
  }
}


/* Called only by sam_sort */
Filestring_T
SAM_header_change_HD_tosorted (FILE *input, int headerlen) {
  Filestring_T fp;
  char buffer[CHUNK], c, c0, c1, c2;

  fp = Filestring_new();

  /* @HD */
  while (headerlen > 0 && (c = fgetc(input)) != '\t') {
    PUTC(c,fp);
    headerlen--;
  }
  if (headerlen > 0) {
    PUTC('\t',fp);
    headerlen--;
  }

  /* VN */
  while (headerlen > 0 && (c = fgetc(input)) != '\t') {
    PUTC(c,fp);
    headerlen--;
  }
  if (headerlen > 0) {
    PUTC('\t',fp);
    headerlen--;
  }

  if (headerlen > 3) {
    /* SO: */
    c0 = fgetc(input);
    c1 = fgetc(input);
    c2 = fgetc(input);
    FPRINTF(fp,"%c%c%c",c0,c1,c2);
    headerlen -= 3;

    if (c0 == 'S' && c1 == 'O' && c2 == ':') {
      FPRINTF(fp,"coordinate\n");
      while (headerlen > 0 && fgetc(input) != '\n') {
	/* Skip given SO value */
	headerlen--;
      }
      headerlen--;
    }
  }

  while (headerlen > CHUNK) {
    fread(buffer,sizeof(char),CHUNK,input);
    Filestring_puts(fp,buffer,/*strlength*/CHUNK);
    headerlen -= CHUNK;
  }
  if (headerlen > 0) {
    fread(buffer,sizeof(char),headerlen,input);
    Filestring_puts(fp,buffer,/*strlength*/headerlen);
  }

  /* Need to call this once for all future writes */
  Filestring_stringify(fp);

  return fp;
}


static void
SAM_header_print_HD (FILE *fp, int nworkers, bool orderedp) {

  fprintf(fp,"@HD");
  fprintf(fp,"\tVN:1.0");	/* or 1.4 */
  if (nworkers > 1 && orderedp == false) {
    fprintf(fp,"\tSO:unsorted");
  } else {
    /* Picard does not recognize type unknown */
    /* fprintf(fp,"\tSO:unknown"); */
    fprintf(fp,"\tSO:unsorted");
  }
  fprintf(fp,"\n");

  return;
}


static void
SAM_header_print_PG (FILE *fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  fprintf(fp,"@PG");
#ifdef GEXACT
  fprintf(fp,"\tID:GEXACT");
  fprintf(fp,"\tPN:gexact");
#elif defined(GSNAP)
  fprintf(fp,"\tID:GSNAP");
  fprintf(fp,"\tPN:gsnap");
#elif defined(PMAP)
  fprintf(fp,"\tID:PMAP");
  fprintf(fp,"\tPN:pmap");
#else
  fprintf(fp,"\tID:GMAP");
  fprintf(fp,"\tPN:gmap");
#endif
  fprintf(fp,"\tVN:%s",PACKAGE_VERSION);

  fprintf(fp,"\tCL:");
  argstart = &(argv[-optind]);
  fprintf(fp,"%s",argstart[0]);
  for (c = 1; c < argc + optind; c++) {
    fprintf(fp," %s",argstart[c]);
  }
  fprintf(fp,"\n");

  return;
}


#if !defined(GEXACT) && !defined(GSNAP)
static void
print_gff_header (FILE *fp, int argc, char **argv, int optind) {
  char **argstart;
  int c;

  fprintf(fp,"##gff-version   3\n");
  fprintf(fp,"# Generated by GMAP version %s using call: ",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 0; c < argc + optind; c++) {
    fprintf(fp," %s",argstart[c]);
  }
  fprintf(fp,"\n");
  return;
}
#endif



void
SAM_header_print_all (FILE *output) {

  if (output == NULL) {
    /* Possible since we are no longer creating a file for OUTPUT_NONE */
    return;
  }

#if defined(GEXACT)
  if (sam_headers_p == true) {
    SAM_header_print_HD(output,nworkers,orderedp);
    SAM_header_print_PG(output,argc,argv,optind_save);
    Univ_IIT_dump_sam(output,chromosome_iit,sam_read_group_id,sam_read_group_name,
		      sam_read_group_library,sam_read_group_platform);
  }

#elif defined(GSNAP)
  if (output_type == SAM_OUTPUT && sam_headers_p == true) {
    SAM_header_print_HD(output,nworkers,orderedp);
    SAM_header_print_PG(output,argc,argv,optind_save);
    Univ_IIT_dump_sam(output,chromosome_iit,sam_read_group_id,sam_read_group_name,
		      sam_read_group_library,sam_read_group_platform);
  }

#else
  if (printtype == GFF3_GENE || printtype == GFF3_MATCH_CDNA || printtype == GFF3_MATCH_EST) {
    print_gff_header(output,argc,argv,optind_save);
      
#ifndef PMAP
  } else if (printtype == SAM && sam_headers_p == true) {
    if (chromosome_iit == NULL) {
      /* Does not make sense to print SAM headers */
    } else {
      SAM_header_print_HD(output,nworkers,orderedp);
      SAM_header_print_PG(output,argc,argv,optind_save);
      Univ_IIT_dump_sam(output,chromosome_iit,sam_read_group_id,sam_read_group_name,
			sam_read_group_library,sam_read_group_platform);
    }
#endif

  }
#endif

  return;
}



/* Called only by sam_sort */
int
SAM_header_length (int *lastchar, FILE *fp) {
  int headerlen = 0;
  int c;

  while (!feof(fp) && (c = getc(fp)) == '@') {
    headerlen++;
    while (!feof(fp) && (c = getc(fp)) != '\n') {
      headerlen++;
    }
    headerlen++;
  }
  /* headerlen++; -- Don't count in header, but as part of first SAM line */

  *lastchar = c;
  return headerlen;
}

