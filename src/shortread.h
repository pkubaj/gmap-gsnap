/* $Id: 7ea41ff790b6eb8231731168c8e1ed931eeaffc9 $ */
#ifndef SHORTREAD_INCLUDED
#define SHORTREAD_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_ZLIB, HAVE_BZLIB */
#endif

#include <stdio.h>
#include "bool.h"
#include "filestring.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif


#define T Shortread_T
typedef struct T *T;

extern void
Shortread_setup (int acc_fieldi_start_in, int acc_fieldi_end_in,
		 bool force_singled_end_p_in, bool filter_chastity_p_in, bool keep_chastity_p_in,
		 bool allow_paired_end_mismatch_p_in, bool fastq_format_p_in,
		 int barcode_length_in, int endtrim_length_in,
		 bool invert_first_p_in, bool invert_second_p_in,
		 bool chop_poly_at_first_p_in, bool chop_poly_at_second_p_in);

extern char *
Shortread_accession (T this);
extern char *
Shortread_header (T this);
extern bool
Shortread_filterp (T this);
extern bool
Shortread_invertedp (T this);

extern int
Shortread_input_init (int *nchars, FILE *fp);

#ifdef HAVE_ZLIB
extern int
Shortread_input_init_gzip (gzFile fp);
#endif

#ifdef HAVE_BZLIB
extern int
Shortread_input_init_bzip2 (Bzip2_T fp);
#endif

extern char *
Shortread_fullpointer (T this);
extern char *
Shortread_trimpointer (T this);

extern char *
Shortread_fullpointer_uc (T this);
extern char *
Shortread_contents_uc (T this);

extern char *
Shortread_queryuc_ptr (T this);
extern char *
Shortread_queryrc (T this);

extern int
Shortread_barcode_length (T this);
extern char *
Shortread_barcode (T this);

extern int
Shortread_left_choplength (T this);
extern int
Shortread_right_choplength (T this);

extern char *
Shortread_quality_string (T this);

extern int
Shortread_fulllength (T this);

extern void
Shortread_free (T *old);

extern bool
Shortread_chop_primers (T queryseq1, T queryseq2);
extern bool
Shortread_find_primers (T queryseq1, T queryseq2);
extern int
Shortread_max_overlap (T queryseq1, T queryseq2);
extern int
Shortread_find_overlap (T queryseq1, T queryseq2);

extern bool
Shortread_left_polya_p (T this, int overhang);
extern bool
Shortread_left_polyt_p (T this, int overhang);
extern bool
Shortread_right_polya_p (T this, int overhang);
extern bool
Shortread_right_polyt_p (T this, int overhang);

extern T
Shortread_read_fasta_text (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
			   FILE **input1, FILE **input2, char *read_files_command,
			   char ***files, int *nfiles, bool single_cell_p, bool skipp);
#ifdef HAVE_ZLIB
extern T
Shortread_read_fasta_gzip (int *nextchar, T *queryseq2,
			   gzFile *input1, gzFile *input2,
			   char ***files, int *nfiles, bool single_cell_p, bool skipp);
#endif

#ifdef HAVE_BZLIB
extern T
Shortread_read_fasta_bzip2 (int *nextchar, T *queryseq2,
			    Bzip2_T *input1, Bzip2_T *input2,
			    char ***files, int *nfiles, bool single_cell_p, bool skipp);
#endif


extern T
Shortread_read_fastq_text (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
			   FILE **input1, FILE **input2, char *read_files_command,
			   char ***files, int *nfiles, bool single_cell_p, bool skipp);

#ifdef HAVE_ZLIB
extern T
Shortread_read_fastq_gzip (int *nextchar, T *queryseq2,
			   gzFile *input1, gzFile *input2,
			   char ***files, int *nfiles, bool single_cell_p, bool skipp);
#endif

#ifdef HAVE_BZLIB
extern T
Shortread_read_fastq_bzip2 (int *nextchar, T *queryseq2,
			    Bzip2_T *input1, Bzip2_T *input2,
			    char ***files, int *nfiles, bool single_cell_p, bool skipp);
#endif

extern T
Shortread_read (int *nextchar, int *nchars1, int *nchars2, T *queryseq2,
		FILE **input1, FILE **input2,
#ifdef HAVE_ZLIB
		gzFile *gzipped1, gzFile *gzipped2,
#endif
#ifdef HAVE_BZLIB
		Bzip2_T *bzipped1, Bzip2_T *bzipped2,
#endif
		bool interleavedp,
		char *read_files_command, char ***files, int *nfiles,
		bool single_cell_p, bool skipp);

extern void
Shortread_print_header (Filestring_T fp, T queryseq1, T queryseq2);

extern void
Shortread_stderr_query_singleend_fasta (T queryseq, T headerseq);
extern void
Shortread_stderr_query_pairedend_fasta (T queryseq1, T queryseq2,
					bool invert_first_p, bool invert_second_p);

extern void
Shortread_print_query_singleend_fastq (Filestring_T fp, T queryseq, T headerseq);
extern void
Shortread_print_query_pairedend_fastq (Filestring_T fp1, Filestring_T fp2, T queryseq1, T queryseq2,
				      bool invert_first_p, bool invert_second_p);

extern void
Shortread_print_query_singleend (Filestring_T fp, T queryseq, T headerseq);
extern void
Shortread_print_query_pairedend (Filestring_T fp1, Filestring_T fp2, T queryseq1, T queryseq2);

extern void
Shortread_print_oneline (Filestring_T fp, T this);
extern void
Shortread_print_oneline_revcomp (Filestring_T fp, T this);

extern void
Shortread_print_chopped_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_revcomp_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);

extern void
Shortread_print_hardclipped_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_hardclipped_revcomp_sam (Filestring_T fp, T this, int hardclip_low, int hardclip_high);

extern void
Shortread_print_hardclipped_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
				     int shift);
extern void
Shortread_print_hardclipped_reverse_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
					     int shift);

extern void
Shortread_print_chopped_end (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end_revcomp (Filestring_T fp, T this, int hardclip_low, int hardclip_high);
extern void
Shortread_print_chopped_end_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
				     int shift);
extern void
Shortread_print_chopped_end_quality_reverse (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
					     int shift);

extern void
Shortread_print_barcode (Filestring_T fp, T this);
extern void
Shortread_print_chop (Filestring_T fp, T this, bool invertp);
extern void
Shortread_print_left_chop_symbols (Filestring_T fp, T this);
extern void
Shortread_print_right_chop_symbols (Filestring_T fp, T this);
extern void
Shortread_print_quality (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
			int shift, bool show_chopped_p);
extern void
Shortread_print_quality_revcomp (Filestring_T fp, T this, int hardclip_low, int hardclip_high,
				int shift, bool show_chopped_p);
extern void
Shortread_print_subseq_uc (Filestring_T fp, T this, int querystart, int queryend);
extern void
Shortread_print_subseq_revcomp_uc (Filestring_T fp, T this, int querystart, int queryend);

#undef T
#endif
