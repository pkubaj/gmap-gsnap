#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For strcpy */
#include <strings.h>		/* For rindex */
#include <ctype.h>
#include <math.h>		/* For rint */

#include "simd.h"

#if 0
#ifdef HAVE_SSE2
#include <emmintrin.h>
#endif
#ifdef HAVE_SSE4_1
#include <smmintrin.h>
#endif
#if !defined(HAVE_SSE4_2)
/* Skip popcnt */
#elif defined(HAVE_POPCNT)
#include <immintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip mm_popcnt */
#elif defined(HAVE_MM_POPCNT)
#include <nmmintrin.h>
#endif

#if !defined(HAVE_SSE4_2)
/* Skip lzcnt/tzcnt */
#elif defined(HAVE_LZCNT) || defined(HAVE_TZCNT)
#include <immintrin.h>
#endif
#endif

#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#ifdef HAVE_ZLIB
#include <zlib.h>
#define GZBUFFER_SIZE 131072
#endif

#ifdef HAVE_BZLIB
#include "bzip2.h"
#endif

#include <signal.h>

#include "assert.h"
#include "except.h"
#include "mem.h"
#include "bool.h"
#include "types.h"
#include "univcoord.h"
#include "fopen.h"
#include "getline.h"
#include "filesuffix.h"

#include "pass.h"
#include "mode.h"
#include "shortread.h"		/* For Shortread_setup */
#include "stopwatch.h"
#include "transcriptome.h"
#include "transcript.h"		/* For Transcript_setup */

#ifdef LARGE_GENOMES
#include "intersect-approx-uint8.h" /* For Intersect_approx_uint4_setup */
#include "intersect-approx-indices-uint8.h"
#else
#include "intersect-approx-indices-uint4.h"
#endif
#include "intersect-approx-uint4.h" /* For Intersect_approx_uint4_setup */
#include "intersect-wdups-indices.h"

#ifdef LARGE_GENOMES
#include "intersect-lower-large.h"
#include "intersect-higher-large.h"
#else
#include "intersect-lower-small.h"
#include "intersect-higher-small.h"
#endif

#ifdef LARGE_GENOMES
#include "intersect-indices2-large.h"
#include "intersect-indices-uint8.h"
#else
#include "intersect-indices2-small.h"
#endif
#include "intersect-indices-small.h"

#ifdef LARGE_GENOMES
#include "intersect-concordance-uint8.h"
#else
#include "intersect-concordance-uint4.h"
#endif

#include "genome.h"
#include "genomebits.h"
#include "genomebits_consec.h"
#include "genomebits_count.h"
#include "genomebits_mismatches.h"
#include "genomebits_trim.h"
#include "genomebits_kmer.h"
#include "genomebits_indel.h"
#include "genome_sites.h"	/* For Genome_sites_setup */
#include "maxent_hr.h"		/* For Maxent_hr_setup */
#include "knownsplicing.h"
#include "knownindels.h"
#include "mapq.h"
#include "concordance.h"	/* For Concordance_setup */
#include "splice.h"		/* For Splice_setup */
#include "spliceends.h"		/* For Spliceends_setup */
#include "oligo.h"		/* For Oligo_setup */

#include "intlistpool.h"
#include "uintlistpool.h"
#include "listpool.h"
#include "trdiagpool.h"
#include "univdiagpool.h"
#include "auxinfopool.h"
#include "pathstorepool.h"
#include "trpathpool.h"
#include "pathpool.h"
#include "hitlistpool.h"
#include "transcriptpool.h"
#include "vectorpool.h"
#include "spliceendsgen.h"

#include "repair.h"		/* For Repair_setup */
#include "path-solve.h"
#include "path-fusion.h"
#include "path-trim.h"
#include "path-eval.h"
#include "pathpair-eval.h"

#include "transcript-remap.h"
#include "path-print-alignment.h"
#include "path-print-m8.h"
#include "path-print-sam.h"
#include "path-learn.h"

#include "trpath-solve.h"
#include "trpath-eval.h"
#include "trpath-convert.h"
#include "transcriptome-search.h"
#include "tr-extension-search.h"
#include "kmer-search.h"
#include "extension-search.h"
/* #include "merge-search.h" */
/* #include "segment-search.h" */
#include "repetitive.h"

#include "indel.h"		/* For Indel_setup */
#include "pathstore.h"		/* For Pathstore_setup */
#include "stage1hr.h"
#include "stage1hr-single.h"
#include "stage1hr-paired.h"
#include "indexdb.h"
#include "localdb-read.h"
#include "resulthr.h"
#include "request.h"
#include "intlist.h"
#include "list.h"
#include "iit-read.h"
#include "iit-read-univ.h"
#include "ef64.h"
#include "datadir.h"

#include "tableuint.h"
#include "single-cell.h"

#include "filestring.h"
#include "outputtype.h"
#include "output.h"
#include "inbuffer.h"
#include "samheader.h"
#include "outbuffer.h"

#include "simplepair.h"
#include "getopt.h"


#define SEED 31415926
#define ONE_GB 1000000000

#define LOCALDB_REGION_SIZE 65536

#define MIN_INDEXDB_SIZE_THRESHOLD 1000

#define MAX_QUERYLENGTH_FOR_ALLOC    100000
#define MAX_GENOMICLENGTH_FOR_ALLOC 1000000


/* File open/close.  Want to turn on in shortread.c also. */
#ifdef DEBUGF
#define debugf(x) x
#else
#define debugf(x)
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/************************************************************************
 *   Global parameters
 ************************************************************************/

static bool genome_align_p = true;
static bool transcriptome_align_p = true;

/* Transcriptome */
static Transcriptome_T transcriptome = NULL;
static Univ_IIT_T transcript_iit = NULL;
static EF64_T transcript_ef64 = NULL;
static IIT_T transcript_map_iit = NULL;
static int *transcript_chrnum_crosstable = NULL;
static Trcoord_T transcriptomelength = 0;
static int ntranscripts = 0;
static int nalignments = 0;

/* Genome */
static Univ_IIT_T chromosome_iit = NULL;
static EF64_T chromosome_ef64 = NULL;
static int circular_typeint = -1;
static Univcoord_T genomelength;
static int nchromosomes = 0;
static bool *circularp = NULL;
static bool any_circular_p;

static char *chrsubset_file = NULL;
static bool *chrsubsetp = NULL;	/* Created if --chrsubset is provided */

static bool *altlocp = NULL;
static Univcoord_T *alias_starts = NULL;
static Univcoord_T *alias_ends = NULL;

static Indexdb_T indexdb = NULL;
static Indexdb_T indexdb_nonstd = NULL; /* For cmet or atoi */
static Indexdb_T tr_indexdb = NULL;


static bool two_pass_p = false;
static Pass_T pass = PASS2;	/* Use pass 2, unless user specifies two-pass mode */
#ifdef HAVE_PTHREAD
pthread_mutex_t pass1_lock;
#endif


/* default_localdb_p applies if user did not specify whether to use
   localdb.  For RNA-seq (for end splicing).  Could also use for
   DNA-seq (for end indels), but slows down program */
static bool default_localdb_p = true;
static bool user_localdb_p = false;
static bool use_localdb_p = false;

static Localdb_T localdb = NULL;

static Genome_T genome = NULL;
static Genome_T genomealt = NULL;
static Genomebits_T genomebits = NULL;
static Genomebits_T genomebits_alt = NULL;

static Genomebits_T transcriptomebits = NULL;
static Genomebits_T transcriptomebits_alt = NULL;

#if 0
static char STANDARD_CHARTABLE[4] = {'A','C','G','T'};
static char CMET_FWD_CHARTABLE[4] = {'A','T','G','T'}; /* CT */
static char CMET_REV_CHARTABLE[4] = {'A','C','A','T'}; /* GA */
static char ATOI_FWD_CHARTABLE[4] = {'G','C','G','T'};     /* AG */
static char ATOI_REV_CHARTABLE[4] = {'A','C','G','C'};     /* TC */
#endif


static bool fastq_format_p = false;
static bool resolve_inner_p = true;

#ifdef CHECK_ASSERTIONS
static bool want_random_p = false; /* order by genomic for testing and comparison */
#else
static bool want_random_p = true; /* randomize among equivalent scores */
#endif

static Stopwatch_T stopwatch = NULL;

/************************************************************************
 *   Program options
 ************************************************************************/

/* Input options */
static char *user_transcriptomedir = NULL;
static char *transcriptome_dbroot = NULL;

static char *user_genomedir = NULL;
static char *genome_dbroot = NULL;

static double skip_prob = 0.0;
static unsigned int part_modulus = 0;
static unsigned int part_interval = 1;
static int barcode_length = 0;
static int endtrim_length = 0;

static bool invert_first_p = false;
static bool invert_second_p = true;

static bool chop_poly_at_first_p = true;
static bool chop_poly_at_second_p = true;

static bool single_cell_p = false;
static char *whitelist_file = NULL;
static int wellpos = 4;

static int acc_fieldi_start = 0;
static int acc_fieldi_end = 0;
static bool force_single_end_p = false;
static bool filter_chastity_p = false;
static bool keep_chastity_p = false;
static bool allow_paired_end_mismatch_p = false;
static bool filter_if_both_p = false;

static char *read_files_command = NULL;
static bool gunzip_p = false;
static bool bunzip2_p = false;
static bool interleavedp = false;

static bool paired_end_p;	/* Determined by open_input_streams */


/* Compute options */
static bool chop_primers_p = false;

static bool query_unk_mismatch_p = false;
static bool genome_unk_mismatch_p = true; /* Needs to be false for path-eval assertions, but needs to be true for circular alignments */
static bool novelsplicingp = false;
static bool splicingp = false;

static bool user_find_dna_chimeras_p = false;
/* If user_specifies, then follow that.  Otherwise, on for
   novelsplicing RNA-seq, because improves those alignments */
static bool find_dna_chimeras_p; 

static bool multiple_sequences_p;
static bool sharedp = false;
static bool preload_shared_memory_p = false;
static bool unload_shared_memory_p = false;
static bool expand_offsets_p = false;

/* Level 5 is now default */
static Access_mode_T offsetsstrm_access = USE_ALLOCATE; /* About 0.5 GB */
static Access_mode_T positions_access = USE_ALLOCATE; /* About 1 GB for transcriptome, 3.5 GB for genome */
static Access_mode_T localdb_access = USE_ALLOCATE;	      /* About 10 GB */
static Access_mode_T locoffsetsstrm_access = USE_ALLOCATE;
static Access_mode_T locpositions_access = USE_ALLOCATE;
static Access_mode_T genome_access = USE_ALLOCATE; /* About 1 GB for oligos and 1 GB for bits */


static int min_intronlength = 9;

static int pairmax_transcriptome;
static int pairmax_linear;
static int pairmax_circular;
static int pairmax_dna = 2000;
static int pairmax_rna = 200000;

static int pass1_min_support = 20;
static double defect_rate = 0.01;
static int expected_pairlength = 200;
static int pairlength_deviation = 100;
static int min_insertlength;	/* expected_pairlength - 5*pairlength_deviation if >= 0 */
static int max_insertlength;	/* expected_pairlength + 5*pairlength_deviation */


#ifdef HAVE_PTHREAD
static pthread_t output_thread_id, *worker_thread_ids;
static pthread_key_t global_request_key;
static int nthreads = 1;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#else
static int nthreads = 0;	/* (int) sysconf(_SC_NPROCESSORS_ONLN) */
#endif

/* static Masktype_T masktype = MASK_REPETITIVE; */
static int subopt_levels = 0;

static int min_indel_end_matches = 4;

/* If negative, then hasn't been specified by user.  If between 0 and
   1, then treated as a fraction of the querylength.  Else, treated as
   an integer */
static double user_nmismatches_filter_float = -1.0;
static double user_mincoverage_filter_float = 0.5;
static int localdb_nmismatches_allowed = 2;

static int max_insertionlen = 6;
static int max_deletionlen = 9;
static Chrpos_T shortsplicedist = 200000;
/* static double min_spliceprob = 0.90; */

static Width_T index1part = 15;
static Width_T index1part_tr = 15;
static Width_T required_index1part = 0;
static Width_T index1interval;
static Width_T min_querylength;
static Width_T index1interval_tr;
static Width_T required_index1interval = 0;

static int indexdb_size_threshold;


#if 0
/* Genes IIT */
static char *genes_file = (char *) NULL;
static IIT_T genes_iit = NULL;	    /* For genome alignment */
static int *genes_divint_crosstable = NULL;
static int *genes_chrnum_crosstable = NULL;
static bool favor_multiexon_p = false;
#endif

static EF64_T repetitive_ef64;

/* Known splicing (from either splicing IIT or two-pass mode) */
static bool splice_evalp = true;

static char *dump_splices_file = NULL;
static bool dump_splices_include_known_p = true;
static char *read_splices_file = NULL;

static Knownsplicing_T knownsplicing = NULL;
static bool knownsplicingp = false;
static bool intron_level_p = false;

static char *user_splicingdir = (char *) NULL;
static char *splicing_file = (char *) NULL;
static IIT_T splicing_iit = NULL;

static int donor_typeint = -1;		/* for splicing_iit */
static int acceptor_typeint = -1;	/* for splicing_iit */

static int *splicing_divint_crosstable = NULL;


/* Known indels, from two-pass mode */
static FILE *dump_indels_fp = NULL;
static char *indels_file = NULL;
static Knownindels_T knownindels = NULL;


/* Need a mutex lock */
static unsigned long long total_mismatches = 0;
static unsigned long long total_querylength = 0;

static Support_T support_tables = NULL;


static Univcoordtableuint_T donor_table = NULL;
static Univcoordtableuint_T acceptor_table = NULL;
static Univcoordtableuint_T antidonor_table = NULL;
static Univcoordtableuint_T antiacceptor_table = NULL;

static Univcoordlist_T donor_startpoints = NULL;
static Univcoordlist_T donor_partners = NULL;
static Univcoordlist_T acceptor_startpoints = NULL;
static Univcoordlist_T acceptor_partners = NULL;
static Univcoordlist_T antidonor_startpoints = NULL;
static Univcoordlist_T antidonor_partners = NULL;
static Univcoordlist_T antiacceptor_startpoints = NULL;
static Univcoordlist_T antiacceptor_partners = NULL;

static Univcoordtable_T indel_table = NULL;

static Uintlist_T insertlengths = NULL;


/* Cmet and AtoI */
static char *user_modedir = NULL;  /* user_cmetdir, user_atoidir */
static Mode_T mode = STANDARD;

/* Masked (e.g., introns) genome */
static char *masked_suffix = (char *) NULL;
static bool maskedp = false;


/* SNPs IIT */
static char *user_snpsdir = NULL;
static char *snps_root = (char *) NULL;
static IIT_T snps_iit = NULL;
static int *snps_divint_crosstable = NULL;


/* Tally IIT */
static char *user_tallydir = NULL;
static char *tally_root = (char *) NULL;
static IIT_T tally_iit = NULL;
static int *tally_divint_crosstable = NULL;

static char *user_runlengthdir = NULL;
static char *runlength_root = (char *) NULL;
static IIT_T runlength_iit = NULL;
static int *runlength_divint_crosstable = NULL;


/* Output options */
static unsigned int output_buffer_size = 1000;
static Outputtype_T output_type = STD_OUTPUT;

/* For Illumina, subtract 64.  For Sanger, subtract 33. */
/* static int quality_score_adj = 64;  -- Stored in mapq.c */

static bool user_quality_score_adj = false;
static bool user_quality_shift = false;
static int quality_shift = 0;   /* For printing, may want -31 */

static bool exception_raise_p = true;
static bool add_paired_nomappers_p = false;
static bool paired_flag_means_concordant_p = false;
static bool quiet_if_excessive_p = false;
static int maxpaths_search = 1000;
static int maxpaths_report = 100;
static bool orderedp = false;
static bool failsonlyp = false;
static bool nofailsp = false;

static bool print_nsnpdiffs_p = false;
static bool print_snplabels_p = false;

static bool show_refdiff_p = false;
static bool clip_overlap_p = false;
static bool merge_overlap_p = false;
static bool merge_samechr_p = false;
static bool method_print_p = false;
static bool print_univdiagonal_p = false;

/* SAM */
static int sam_headers_batch = -1;
static bool sam_headers_p = true;
static bool sam_hardclip_use_S_p = false;
static bool sam_insert_0M_p = true;
static bool sam_cigar_extended_p = false;
static bool sam_multiple_primaries_p = false;
static bool sam_sparse_secondaries_p = false;
static char *sam_read_group_id = NULL;
static char *sam_read_group_name = NULL;
static char *sam_read_group_library = NULL;
static char *sam_read_group_platform = NULL;
static bool force_xs_direction_p = false;
static bool md_report_snps_p = false;
static bool allow_soft_clips_p = true;
static bool extend_soft_clips_p = false;
static Cigar_action_T cigar_action = CIGAR_ACTION_WARNING;

static bool only_concordant_p = false;
static bool omit_concordant_uniq_p = false;
static bool omit_concordant_mult_p = false;
static bool omit_softclipped_p = false;
static bool only_tr_consistent_p = false;


/* Input/output */
static bool split_simple_p = false;
static char *split_output_root = NULL;
static char *output_file = NULL;
static char *failedinput_root = NULL;
static bool appendp = false;
static Outbuffer_T outbuffer;
static Inbuffer_T inbuffer;
static unsigned int input_buffer_size = 10000; /* previously inbuffer_nspaces */
static bool timingp = false;
static bool unloadp = false;


/* getopt used alphabetically: AaBCcDdeGgiJjKklMmNnOoQqstVvwYyZz7 */

static struct option long_options[] = {
  /* Input options */
  {"transcriptdir", required_argument, 0, 'C'},	/* user_transcriptomedir */
  {"transcriptdb", required_argument, 0, 'c'}, /* transcriptome_dbroot */
  {"transcriptome-mode", required_argument, 0, 0}, /* genome_align_p, transcriptome_align_p */
  {"dir", required_argument, 0, 'D'},	/* user_genomedir */
  {"db", required_argument, 0, 'd'}, /* genome_dbroot */
  {"two-pass", no_argument, 0, 0},   /* two_pass_p */
  {"use-localdb", required_argument, 0, 0}, /* user_localdb_p, use_localdb_p */
  {"kmer", required_argument, 0, 'k'}, /* required_index1part, index1part */
  {"sampling", required_argument, 0, 0}, /* required_index1interval, index1interval */
  {"align-fraction", required_argument, 0, 0}, /* skip_prob */
  {"part", required_argument, 0, 'q'}, /* part_modulus, part_interval */
  {"orientation", required_argument, 0, 0}, /* single_cell_p, invert_first_p, invert_second_p */
  {"input-buffer-size", required_argument, 0, 0}, /* input_buffer_size */
  {"barcode-length", required_argument, 0, 0},	  /* barcode_length */
  {"endtrim-length", required_argument, 0, 0},	  /* endtrim_length */
  {"fastq-id-start", required_argument, 0, 0},	  /* acc_fieldi_start */
  {"fastq-id-end", required_argument, 0, 0},	  /* acc_fieldi_end */
  {"force-single-end", no_argument, 0, 0},	  /* force_single_end_p */
  {"filter-chastity", required_argument, 0, 0},	/* filter_chastity_p, filter_if_both_p */
  {"allow-pe-name-mismatch", no_argument, 0, 0}, /* allow_paired_end_mismatch_p */
  {"10x-whitelist", required_argument, 0, 0},	 /* whitelist_file */
  {"10x-well-position", required_argument, 0, 0}, /* wellpos */

  {"read-files-command", required_argument, 0, 0}, /* read_files_command */
#ifdef HAVE_ZLIB
  {"gunzip", no_argument, 0, 0}, /* gunzip_p */
#endif

#ifdef HAVE_BZLIB
  {"bunzip2", no_argument, 0, 0}, /* bunzip2_p */
#endif

  {"interleaved", no_argument, 0, 0}, /* interleavedp */

  /* Compute options */
  {"use-shared-memory", required_argument, 0, 0}, /* sharedp */
  {"preload-shared-memory", no_argument, 0, 0},	  /* preload_shared_memory_p */
  {"unload-shared-memory", no_argument, 0, 0},	  /* unload_shared_memory_p */
#ifdef HAVE_MMAP
  {"batch", required_argument, 0, 'B'}, /* offsetsstrm_access, positions_access, genome_access */
#endif
  {"expand-offsets", required_argument, 0, 0}, /* expand_offsets_p */
  {"pairmax-dna", required_argument, 0, 0}, /* pairmax_dna */
  {"pairmax-rna", required_argument, 0, 0}, /* pairmax_rna */

  {"resolve-inner", required_argument, 0, 0}, /* resolve_inner_p */
  {"pairexpect", required_argument, 0, 0},  /* expected_pairlength */
  {"pairdev", required_argument, 0, 0},  /* pairlength_deviation */

  {"pass1-min-support", required_argument, 0, 0}, /* pass1_min_support */

  {"nthreads", required_argument, 0, 't'}, /* nthreads */
  {"adapter-strip", required_argument, 0, 'a'},	/* chop_primers_p */

  {"query-unk-mismatch", required_argument, 0, 0}, /* query_unk_mismatch_p */
  {"genome-unk-mismatch", required_argument, 0, 0}, /* genome_unk_mismatch_p */

  {"trim-mismatch-score", required_argument, 0, 0}, /* trim_mismatch_score, user_trim_mismatch_score_p */
  {"novelsplicing", required_argument, 0, 'N'}, /* novelsplicingp */
  {"find-dna-chimeras", required_argument, 0, 0}, /* user_find_dna_chimeras_p */

  {"max-mismatches", required_argument, 0, 'm'}, /* user_nmismatches_filter_float */
  {"min-coverage", required_argument, 0, 0}, /* user_mincoverage_filter_float */

  {"indel-endlength", required_argument, 0, 0}, /* min_indel_end_matches */

  {"max-insertions", required_argument, 0, 0}, /* max_insertionlen */
  {"max-deletions", required_argument, 0, 0}, /* max_deletionlen */
  {"suboptimal-levels", required_argument, 0, 'M'}, /* subopt_levels */

  {"localsplicedist", required_argument, 0, 'w'}, /* shortsplicedist */
  {"splicingdir", required_argument, 0, 0},	  /* user_splicingdir */
  {"use-splicing", required_argument, 0, 's'}, /* splicing_iit, knownsplicingp, find_dna_chimeras_p */

  {"splices-noeval", no_argument, 0, 0},	     /* splice_evalp */
  {"splices-dump", required_argument, 0, 0},	       /* dump_splices_fp */
  {"splices-exclude-known", no_argument, 0, 0},	       /* dump_splices_include_known_p */
  {"splices-read", required_argument, 0, 0},    /* splices_file */

  {"indels-dump", required_argument, 0, 0},	       /* dump_indels_fp */
  {"indels-read", required_argument, 0, 0},    /* indels_file */

  {"merge-distant-samechr", no_argument, 0, 0},		    /* merge_samechr_p */

  {"cmetdir", required_argument, 0, 0}, /* user_modedir */
  {"atoidir", required_argument, 0, 0}, /* user_modedir */
  {"mode", required_argument, 0, 0}, /* mode */

  {"use-mask", required_argument, 0, 'e'}, /* masked_suffix */
  {"snpsdir", required_argument, 0, 'V'},   /* user_snpsdir */
  {"use-snps", required_argument, 0, 'v'}, /* snps_root */

  {"tallydir", required_argument, 0, 0},   /* user_tallydir */
  {"use-tally", required_argument, 0, 0}, /* tally_root */

  {"runlengthdir", required_argument, 0, 0},   /* user_runlengthdir */
  {"use-runlength", required_argument, 0, 0}, /* runlength_root */

  /* Output options */
  {"output-buffer-size", required_argument, 0, 0}, /* output_buffer_size */
  {"format", required_argument, 0, 'A'}, /* output_type */

  {"quality-protocol", required_argument, 0, 0}, /* quality_score_adj, quality_shift */
  {"quality-zero-score", required_argument, 0, 'J'}, /* quality_score_adj */
  {"quality-print-shift", required_argument, 0, 'j'}, /* quality_shift */
  {"sam-headers-batch", required_argument, 0, 0},	/* sam_headers_batch */
  {"no-sam-headers", no_argument, 0, 0},	/* sam_headers_p */
  {"sam-hardclip-use-S", no_argument, 0, 0},		/* sam_hardclip_use_S_p */
  {"sam-use-0M", required_argument, 0, 0},		/* sam_insert_0M_p */
  {"sam-extended-cigar", no_argument, 0, 0},	/* sam_cigar_extended_p */
  {"sam-multiple-primaries", no_argument, 0, 0}, /* sam_multiple_primaries_p */
  {"sam-sparse-secondaries", no_argument, 0, 0}, /* sam_sparse_secondaries_p */
  {"read-group-id", required_argument, 0, 0},	/* sam_read_group_id */
  {"read-group-name", required_argument, 0, 0},	/* sam_read_group_name */
  {"read-group-library", required_argument, 0, 0},	/* sam_read_group_library */
  {"read-group-platform", required_argument, 0, 0},	/* sam_read_group_platform */
  {"force-xs-dir", no_argument, 0, 0},			/* force_xs_direction_p */
  {"md-report-snps", no_argument, 0, 0},	    /* md_report_snps_p */
  {"no-soft-clips", no_argument, 0, 0},	    /* allow_soft_clips_p */
  {"extend-soft-clips", no_argument, 0, 0},		/* extend_soft_clips_p */
  {"action-if-cigar-error", required_argument, 0, 0},	/* cigar_action */

  {"noexceptions", no_argument, 0, '0'}, /* exception_raise_p */

  {"maxsearch", required_argument, 0, 0}, /* maxpaths_search */
  {"npaths", required_argument, 0, 'n'}, /* maxpaths_report */

  {"add-paired-nomappers", no_argument, 0, 0}, /* add_paired_nomappers_p */
  {"paired-flag-means-concordant", required_argument, 0, 0}, /* paired_flag_means_concordant_p */
  {"quiet-if-excessive", no_argument, 0, 'Q'}, /* quiet_if_excessive_p */
  {"ordered", no_argument, 0, 'O'}, /* orderedp */
  {"clip-overlap", no_argument, 0, 0},	     /* clip_overlap_p */
  {"merge-overlap", no_argument, 0, 0},	     /* merge_overlap_p */
  {"show-method", no_argument, 0, 0},	     /* method_print_p */
  {"show-univdiagonal", no_argument, 0, 0},	     /* print_univdiagonal_p */

  {"show-refdiff", no_argument, 0, 0},	       /* show_refdiff_p */
  {"print-snps", no_argument, 0, 0}, /* print_snplabels_p */
  {"failsonly", no_argument, 0, 0}, /* failsonlyp */
  {"nofails", no_argument, 0, 0}, /* nofailsp */
  {"chrsubset", required_argument, 0, 0}, /* chrsubsetp */
  {"output-file", required_argument, 0, 'o'}, /* output_file */

  {"split-output", required_argument, 0, 0}, /* split_output_root */
  {"split-simple", no_argument, 0, 0},	     /* split_simple_p */

  {"failed-input", required_argument, 0, 0}, /* failed_input_root */
  {"append-output", no_argument, 0, 0},	     /* appendp */

  {"order-among-best", required_argument, 0, 0}, /* want_random_p */

  {"only-concordant", no_argument, 0, 0}, /* only_concordant_p */
  {"omit-concordant-uniq", no_argument, 0, 0}, /* omit_concordant_uniq_p */
  {"omit-concordant-mult", no_argument, 0, 0}, /* omit_concordant_mult_p */
  {"omit-softclipped", no_argument, 0, 0}, /* omit_softclipped_p */
  {"only-tr-consistent", no_argument, 0, 0}, /* only_tr_consistent_p */

  /* Diagnostic options */
  {"time", no_argument, 0, 0},	/* timingp */
  {"unload", no_argument, 0, 0},	/* unloadp */

  /* Obsolete, but included for backward compatibility */
  {"use-sarray", required_argument, 0, 0},
  {"terminal-threshold", required_argument, 0, 0},

  /* Help options */
  {"check", no_argument, 0, 0}, /* check_compiler_assumptions */
  {"version", no_argument, 0, 0}, /* print_program_version */
  {"help", no_argument, 0, 0}, /* print_program_usage */
  {0, 0, 0, 0}
};


static void
print_program_version () {
  char *genomedir;

  fprintf(stdout,"\n");
  fprintf(stdout,"GSNAP: Genomic Short Nucleotide Alignment Program\n");
  fprintf(stdout,"Part of GMAP package, version %s\n",PACKAGE_VERSION);
  fprintf(stdout,"Build target: %s\n",TARGET);
  fprintf(stdout,"Features: ");
#ifdef HAVE_PTHREAD
  fprintf(stdout,"pthreads enabled, ");
#else
  fprintf(stdout,"no pthreads, ");
#endif
#ifdef HAVE_ALLOCA
  fprintf(stdout,"alloca available, ");
#else
  fprintf(stdout,"no alloca, ");
#endif
#ifdef HAVE_ZLIB
  fprintf(stdout,"zlib available, ");
#else
  fprintf(stdout,"no zlib, ");
#endif
#ifdef HAVE_MMAP
  fprintf(stdout,"mmap available, ");
#else
  fprintf(stdout,"no mmap, ");
#endif
#ifdef WORDS_BIGENDIAN
  fprintf(stdout,"bigendian, ");
#else
  fprintf(stdout,"littleendian, ");
#endif
#ifdef HAVE_SIGACTION
  fprintf(stdout,"sigaction available, ");
#else
  fprintf(stdout,"no sigaction, ");
#endif
#ifdef HAVE_64_BIT
  fprintf(stdout,"64 bits available");
#else
  fprintf(stdout,"64 bits not available");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Popcnt:");
#ifdef HAVE_POPCNT
  fprintf(stdout," popcnt/lzcnt/tzcnt");
#endif
#ifdef HAVE_MM_POPCNT
  fprintf(stdout," mm_popcnt");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"Builtin functions:");
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stdout," builtin_clz");
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stdout," builtin_ctz");
#endif
#ifdef HAVE_BUILTIN_POPCOUNT
  fprintf(stdout," builtin_popcount");
#endif
  fprintf(stdout,"\n");

  fprintf(stdout,"SIMD functions compiled:");
#ifdef HAVE_ALTIVEC
  fprintf(stdout," Altivec");
#endif
#ifdef HAVE_MMX
  fprintf(stdout," MMX");
#endif
#ifdef HAVE_SSE
  fprintf(stdout," SSE");
#endif
#ifdef HAVE_SSE2
  fprintf(stdout," SSE2");
#endif
#ifdef HAVE_SSE3
  fprintf(stdout," SSE3");
#endif
#ifdef HAVE_SSSE3
  fprintf(stdout," SSSE3");
#endif
#ifdef HAVE_SSE4_1
  fprintf(stdout," SSE4.1");
#endif
#ifdef HAVE_SSE4_2
  fprintf(stdout," SSE4.2");
#endif
#ifdef HAVE_AVX2
  fprintf(stdout," AVX2");
#endif
#ifdef HAVE_AVX512
  fprintf(stdout," AVX512");
#endif
#ifdef HAVE_AVX512BW
  fprintf(stdout," AVX512BW");
#endif
  fprintf(stdout,"\n");


  fprintf(stdout,"Sizes: off_t (%d), size_t (%d), unsigned int (%d), long int (%d), long long int (%d)\n",
	  (int) sizeof(off_t),(int) sizeof(size_t),(int) sizeof(unsigned int),(int) sizeof(long int),(int) sizeof(long long int));
  fprintf(stdout,"Default gmap directory (compiled): %s\n",GMAPDB);
  genomedir = Datadir_find_genomedir(/*user_genomedir*/NULL);
  fprintf(stdout,"Default gmap directory (environment): %s\n",genomedir);
  FREE(genomedir);
  fprintf(stdout,"Thomas D. Wu, Genentech, Inc.\n");
  fprintf(stdout,"Contact: twu@gene.com\n");
  fprintf(stdout,"\n");
  return;
}

/* This flag is not well-supported, and therefore hidden, but
   kept for backwards compatibility */
/*  -R, --rel=STRING               Release\n\ */


static void
print_program_usage ();


static void
check_compiler_assumptions () {
  unsigned int x = rand(), y = rand();
#ifdef HAVE_SSE2
  int z;
  __m128i a;
#ifdef HAVE_SSE4_1
  char negx, negy;
#endif
#endif

#ifdef HAVE_SSE2
  /* With -mavx, compiler may use assembly instructions for _mm_set1_epi32 that don't work on non-AVX machines */
  fprintf(stderr,"Checking compiler assumptions for SSE2: ");
  fprintf(stderr,"%08X %08X",x,y);
  a = _mm_xor_si128(_mm_set1_epi32(x),_mm_set1_epi32(y));
  z = _mm_cvtsi128_si32(a);
  fprintf(stderr," xor=%08X\n",z);
#endif


#ifdef HAVE_SSE4_1
  if ((negx = (char) x) > 0) {
    negx = -negx;
  }
  if ((negy = (char) y) > 0) {
    negy = -negy;
  }

  fprintf(stderr,"Checking compiler assumptions for SSE4.1: ");
  fprintf(stderr,"%d %d",negx,negy);
  a = _mm_max_epi8(_mm_set1_epi8(negx),_mm_set1_epi8(negy));
  z = _mm_extract_epi8(a,0);
  fprintf(stderr," max=%d => ",z);
  if (negx > negy) {
    if (z == (int) negx) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  } else {
    if (z == (int) negy) {
      fprintf(stderr,"compiler sign extends\n"); /* technically incorrect, but SIMD procedures behave properly */
    } else {
      fprintf(stderr,"compiler zero extends\n");
    }
  }
#endif


#ifdef HAVE_SSE4_2
  fprintf(stderr,"Checking compiler assumptions for SSE4.2 options: ");
  fprintf(stderr,"%08X ",x);
#ifdef HAVE_LZCNT
  fprintf(stderr,"_lzcnt_u32=%d ",_lzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CLZ
  fprintf(stderr,"__builtin_clz=%d ",__builtin_clz(x));
#endif
#ifdef HAVE_TZCNT
  fprintf(stderr,"_tzcnt_u32=%d ",_tzcnt_u32(x));
#endif
#ifdef HAVE_BUILTIN_CTZ
  fprintf(stderr,"__builtin_ctz=%d ",__builtin_ctz(x));
#endif

#ifdef HAVE_POPCNT
  fprintf(stderr,"_popcnt32=%d ",_popcnt32(x));
#endif
#if defined(HAVE_MM_POPCNT)
  fprintf(stderr,"_mm_popcnt_u32=%d ",_mm_popcnt_u32(x));
#endif
#if defined(HAVE_BUILTIN_POPCOUNT)
  fprintf(stderr,"__builtin_popcount=%d ",__builtin_popcount(x));
#endif

  fprintf(stderr,"\n");
#endif

  fprintf(stderr,"Finished checking compiler assumptions\n");

  return;
}


/************************************************************************/

/* Pass1 updates the global list of intron info */
static void
process_request_pass1 (Request_T request,

		       int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc,

		       Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool, 
		       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		       Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		       Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		       Spliceendsgen_T spliceendsgen,
		       Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3) {
  Stage1_T stage1_5, stage1_3;
  Shortread_T queryseq1, queryseq2;
  Path_T *patharray, *patharray5, *patharray3, path5, path3;
  Pathpair_T *pathpairarray, pathpair;

  int npaths_primary, npaths_altloc, npaths5_primary, npaths5_altloc, npaths3_primary, npaths3_altloc, i;
  int first_absmq, second_absmq, first_absmq5, second_absmq5, first_absmq3, second_absmq3;
  Pairtype_T final_pairtype;

  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  if (single_cell_p == true) {
    Spliceendsgen_reset(spliceendsgen);
    stage1_3 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_3,streamspace_max_alloc,streamspace_alloc,
				   queryseq2,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/true,/*first_read_p*/true,
				   /*pass*/PASS1);
    if (npaths_primary + npaths_altloc == 1) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&pass1_lock);
#endif
      Path_learn_defect_rate(patharray[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&pass1_lock);
#endif
    }

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }
    for (i = 0; i < npaths_primary + npaths_altloc; i++) {
      Path_free(&(patharray[i]),intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(patharray);

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

  } else if (queryseq2 == NULL) {
    Spliceendsgen_reset(spliceendsgen);
    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_5,streamspace_max_alloc,streamspace_alloc,
				   queryseq1,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/false,/*first_read_p*/true,
				   /*pass*/PASS1);
    if (npaths_primary + npaths_altloc == 1) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&pass1_lock);
#endif
      Path_learn_defect_rate(patharray[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&pass1_lock);
#endif
    }

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }

    for (i = 0; i < npaths_primary + npaths_altloc; i++) {
      Path_free(&(patharray[i]),intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(patharray);

    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

  } else if (Shortread_fulllength(queryseq1) < min_querylength &&
	     Shortread_fulllength(queryseq2) < min_querylength) {
    /* Nothing to learn */

  } else if (Shortread_fulllength(queryseq1) < min_querylength) {
    /* Solve only queryseq2 */
    Spliceendsgen_reset(spliceendsgen);
    stage1_3 = Stage1_new(/*first_read_p*/false);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_3,streamspace_max_alloc,streamspace_alloc,
				   queryseq2,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/false,/*first_read_p*/false,
				   /*pass*/PASS1);
    if (npaths_primary + npaths_altloc == 1) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&pass1_lock);
#endif
      Path_learn_defect_rate(patharray[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&pass1_lock);
#endif
    }

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }

    for (i = 0; i < npaths_primary + npaths_altloc; i++) {
      Path_free(&(patharray[i]),intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(patharray);

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

  } else if (Shortread_fulllength(queryseq2) < min_querylength) {
    /* Solve only queryseq1 */
    Spliceendsgen_reset(spliceendsgen);
    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_5,streamspace_max_alloc,streamspace_alloc,
				   queryseq1,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/false,/*first_read_p*/true,
				   /*pass*/PASS1);
    if (npaths_primary + npaths_altloc == 1) {
#ifdef HAVE_PTHREAD
      pthread_mutex_lock(&pass1_lock);
#endif
      Path_learn_defect_rate(patharray[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
      pthread_mutex_unlock(&pass1_lock);
#endif
    }

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }

    for (i = 0; i < npaths_primary + npaths_altloc; i++) {
      Path_free(&(patharray[i]),intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT(patharray);

    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

  } else {
    Spliceendsgen_reset(spliceendsgen5);
    Spliceendsgen_reset(spliceendsgen3);

    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    stage1_3 = Stage1_new(/*first_read_p*/false);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    if ((pathpairarray = Stage1_paired_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,&final_pairtype,
					    &patharray5,&npaths5_primary,&npaths5_altloc,&first_absmq5,&second_absmq5,
					    &patharray3,&npaths3_primary,&npaths3_altloc,&first_absmq3,&second_absmq3,
					    stage1_5,stage1_3,queryseq1,queryseq2,repetitive_ef64,
					    knownsplicing,knownindels,(Chrpos_T) pairmax_linear,
					    streamspace_max_alloc,streamspace_alloc,
					    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb_nmismatches_allowed,
					    trdiagpool,univdiagpool,auxinfopool,pathstorepool,
					    intlistpool,uintlistpool,univcoordlistpool,
					    listpool,trpathpool,pathpool,vectorpool,hitlistpool,
					    transcriptpool,spliceendsgen5,spliceendsgen3,
					    /*pass*/PASS1)) != NULL) {
      /* Paired or concordant hits found */
      if (npaths_primary + npaths_altloc == 1) {
	pathpair = pathpairarray[0];
	path5 = pathpair->path5;
	path3 = pathpair->path3;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pass1_lock);
#endif
	Pathpair_learn_insertlengths(pathpairarray[0],&insertlengths);
	Path_learn_defect_rate(path5,&total_mismatches,&total_querylength);
	Path_learn_defect_rate(path3,&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pass1_lock);
#endif
      }

      if (support_tables != NULL) {
	Support_revise_paired(support_tables,pathpairarray,npaths_primary + npaths_altloc);
      }

      for (i = 0; i < npaths_primary + npaths_altloc; i++) {
	Pathpair_free(&(pathpairarray[i]),intlistpool,univcoordlistpool,
		      listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(pathpairarray);

    } else {
      /* Process unpaired ends */
      if (npaths5_primary + npaths5_altloc == 1) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pass1_lock);
#endif
	Path_learn_defect_rate(patharray5[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pass1_lock);
#endif
      }

      if (npaths3_primary + npaths3_altloc == 1) {
#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pass1_lock);
#endif
	Path_learn_defect_rate(patharray3[0],&total_mismatches,&total_querylength);
#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pass1_lock);
#endif
      }

      if (support_tables != NULL) {
	Support_revise_single(support_tables,patharray5,npaths5_primary + npaths5_altloc);
	Support_revise_single(support_tables,patharray3,npaths3_primary + npaths3_altloc);
      }

      for (i = 0; i < npaths5_primary + npaths5_altloc; i++) {
	Path_free(&(patharray5[i]),intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      for (i = 0; i < npaths3_primary + npaths3_altloc; i++) {
	Path_free(&(patharray3[i]),intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray3);
      FREE_OUT(patharray5);
    }

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);
    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);
  }

  return;
}


static Filestring_T
process_request_pass2 (Filestring_T *fp_failedinput, Filestring_T *fp_failedinput_1, Filestring_T *fp_failedinput_2,
		       double *worker_runtime, Request_T request,
		       
		       int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		       Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		       unsigned short *localdb_alloc,

		       Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
		       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		       Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		       Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool,
		       Spliceendsgen_T spliceendsgen,
		       Spliceendsgen_T spliceendsgen5, Spliceendsgen_T spliceendsgen3,
		       Stopwatch_T worker_stopwatch) {
  Filestring_T fp;
  Result_T result;
  int jobid;

  Stage1_T stage1_5, stage1_3;
  Shortread_T queryseq1, queryseq2;
  Path_T *patharray, *patharray5, *patharray3;
  Pathpair_T *pathpairarray;

  int npaths_primary, npaths_altloc, npaths5_primary, npaths5_altloc, npaths3_primary, npaths3_altloc;
  int first_absmq, second_absmq, first_absmq5, second_absmq5, first_absmq3, second_absmq3;
  Pairtype_T final_pairtype;

  jobid = Request_id(request);
  queryseq1 = Request_queryseq1(request);
  queryseq2 = Request_queryseq2(request);

  /* printf("%s\n",Shortread_accession(queryseq1)); */

  if (worker_stopwatch != NULL) {
    Stopwatch_start(worker_stopwatch);
  }

  if (single_cell_p == true) {
    Spliceendsgen_reset(spliceendsgen);
    stage1_3 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_3,streamspace_max_alloc,streamspace_alloc,
				   queryseq2,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/true,/*first_read_p*/true,
				   /*pass*/PASS2);

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }
    result = Result_single_read_new(jobid,(void **) patharray,npaths_primary,npaths_altloc,first_absmq,second_absmq);
    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

    return fp;

  } else if (queryseq2 == NULL) {
    Spliceendsgen_reset(spliceendsgen);
    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    patharray = Stage1_single_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,
				   stage1_5,streamspace_max_alloc,streamspace_alloc,
				   queryseq1,repetitive_ef64,knownsplicing,knownindels,localdb,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   localdb_nmismatches_allowed,

				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				   spliceendsgen,/*single_cell_p*/false,/*first_read_p*/true,
				   /*pass*/PASS2);

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray,npaths_primary + npaths_altloc);
    }
    result = Result_single_read_new(jobid,(void **) patharray,npaths_primary,npaths_altloc,first_absmq,second_absmq);
    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

    return fp;

  } else if (Shortread_fulllength(queryseq1) < min_querylength &&
	     Shortread_fulllength(queryseq2) < min_querylength) {
    patharray3 = patharray5 = (Path_T *) NULL;
    /* No results for Support_revise */
    result = Result_paired_as_singles_new(jobid,(void **) patharray5,/*npaths5_primary*/0,/*npaths5_altloc*/0,
					  /*first_absmq5*/0,/*second_absmq5*/0,
					  (void **) patharray3,/*npaths3_primary*/0,/*npaths3_altloc*/0,
					  /*first_absmq3*/0,/*second_absmq3*/0);
    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    return fp;

  } else if (Shortread_fulllength(queryseq1) < min_querylength) {
    /* Solve only queryseq2 */
    Spliceendsgen_reset(spliceendsgen);
    stage1_3 = Stage1_new(/*first_read_p*/false);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    patharray5 = (Path_T *) NULL;
    patharray3 = Stage1_single_read(&npaths3_primary,&npaths3_altloc,&first_absmq3,&second_absmq3,
				    stage1_3,streamspace_max_alloc,streamspace_alloc,
				    queryseq2,repetitive_ef64,knownsplicing,knownindels,localdb,
				    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				    localdb_nmismatches_allowed,

				    trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				    intlistpool,uintlistpool,univcoordlistpool,
				    listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				    spliceendsgen,/*single_cell_p*/false,/*first_read_p*/false,
				    /*pass*/PASS2);

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray3,npaths3_primary + npaths3_altloc);
    }
    result = Result_paired_as_singles_new(jobid,(void **) patharray5,/*npaths5_primary*/0,/*npaths5_altloc*/0,
					  /*first_absmq5*/0,/*second_absmq5*/0,
					  (void **) patharray3,npaths3_primary,npaths3_altloc,
					  first_absmq3,second_absmq3);
    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

    return fp;

  } else if (Shortread_fulllength(queryseq2) < min_querylength) {
    /* Solve only queryseq1 */
    Spliceendsgen_reset(spliceendsgen);
    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    patharray3 = (Path_T *) NULL;
    patharray5 = Stage1_single_read(&npaths5_primary,&npaths5_altloc,&first_absmq5,&second_absmq5,
				    stage1_5,streamspace_max_alloc,streamspace_alloc,
				    queryseq1,repetitive_ef64,knownsplicing,knownindels,localdb,
				    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				    localdb_nmismatches_allowed,

				    trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				    intlistpool,uintlistpool,univcoordlistpool,
				    listpool,trpathpool,pathpool,transcriptpool,vectorpool,hitlistpool,
				    spliceendsgen,/*single_cell_p*/false,/*first_read_p*/true,
				    /*pass*/PASS2);

    if (support_tables != NULL) {
      Support_revise_single(support_tables,patharray5,npaths5_primary + npaths5_altloc);
    }
    result = Result_paired_as_singles_new(jobid,(void **) patharray5,npaths5_primary,npaths5_altloc,
					  first_absmq5,second_absmq5,
					  (void **) patharray3,/*npaths3_primary*/0,/*npaths3_altloc*/0,
					  /*first_absmq3*/0,/*second_absmq3*/0);
    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

    return fp;

  } else {
    Spliceendsgen_reset(spliceendsgen5);
    Spliceendsgen_reset(spliceendsgen3);

    stage1_5 = Stage1_new(/*first_read_p*/true);
    Stage1_set_sequence(stage1_5,queryseq1,localdb);

    stage1_3 = Stage1_new(/*first_read_p*/false);
    Stage1_set_sequence(stage1_3,queryseq2,localdb);

    if ((pathpairarray = Stage1_paired_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,&final_pairtype,
					    &patharray5,&npaths5_primary,&npaths5_altloc,&first_absmq5,&second_absmq5,
					    &patharray3,&npaths3_primary,&npaths3_altloc,&first_absmq3,&second_absmq3,
					    stage1_5,stage1_3,queryseq1,queryseq2,repetitive_ef64,
					    knownsplicing,knownindels,(Chrpos_T) pairmax_linear,
					    streamspace_max_alloc,streamspace_alloc,
					    novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb_nmismatches_allowed,
					    trdiagpool,univdiagpool,auxinfopool,pathstorepool,
					    intlistpool,uintlistpool,univcoordlistpool,
					    listpool,trpathpool,pathpool,vectorpool,hitlistpool,
					    transcriptpool,spliceendsgen5,spliceendsgen3,
					    /*pass*/PASS2)) != NULL) {

      /* Paired or concordant hits found */
      if (support_tables != NULL) {
	Support_revise_paired(support_tables,pathpairarray,npaths_primary + npaths_altloc);
      }
      result = Result_paired_read_new(jobid,(void **) pathpairarray,npaths_primary,npaths_altloc,first_absmq,second_absmq,
				      final_pairtype);

#if 0
    } else if (chop_primers_p == false || Shortread_chop_primers(queryseq1,queryseq2) == false) {
      /* SAM output is now set for for poly-A/T chop, rather than primer chop */
      /* No paired or concordant hits found, and no adapters found */
      /* Report ends as unpaired */
      if (support_tables != NULL) {
	Support_revise_single(support_tables,patharray5,npaths5_primary + npaths5_altloc);
	Support_revise_single(support_tables,patharray3,npaths3_primary + npaths3_altloc);
      }
      result = Result_paired_as_singles_new(jobid,(void **) patharray5,npaths5_primary,npaths5_altloc,first_absmq5,second_absmq5,
					    (void **) patharray3,npaths3_primary,npaths3_altloc,first_absmq3,second_absmq3);
#endif

#if 0
    } else if (0) {
      /* Try with potential primers chopped.  queryseq1 and queryseq2 altered by Shortread_chop_primers. */
      for (i = 0; i < npaths5_primary + npaths5_altloc; i++) {
	Path_free(&(patharray5[i]),intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray5);
      
      for (i = 0; i < npaths3_primary + npaths3_altloc; i++) {
	Path_free(&(patharray3[i]),intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }
      FREE_OUT(patharray3);
      
      if ((pathpairarray = Stage1_paired_read(&npaths_primary,&npaths_altloc,&first_absmq,&second_absmq,&final_pairtype,
					      &patharray5,&npaths5_primary,&npaths5_altloc,&first_absmq5,&second_absmq5,
					      &patharray3,&npaths3_primary,&npaths3_altloc,&first_absmq3,&second_absmq3,
					      stage1_5,stage1_3,queryseq1,queryseq2,repetitive_ef64,
					      knownsplicing,knownindels,(Chrpos_T) pairmax_linear,
					      streamspace_max_alloc,streamspace_alloc,
					      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,localdb_nmismatches_allowed,
					      trdiagpool,univdiagpool,auxinfopool,pathstorepool,
					      intlistpool,uintlistpool,univcoordlistpool,
					      listpool,trpathpool,pathpool,vectorpool,hitlistpool,
					      transcriptpool,spliceendsgen5,spliceendsgen3,
					      /*pass*/PASS2)) != NULL) {
	/* Paired or concordant hits found, after chopping adapters */
	if (support_tables != NULL) {
	  Support_revise_paired(support_tables,pathpairarray,npaths_primary + npaths_altloc);
	}
	result = Result_paired_read_new(jobid,(void **) pathpairarray,npaths_primary,npaths_altloc,first_absmq,second_absmq,
					final_pairtype);
	
      } else {
	/* No paired or concordant hits found, after chopping adapters */
	if (support_tables != NULL) {
	  Support_revise_single(support_tables,patharray5,npaths5_primary + npaths5_altloc);
	  Support_revise_single(support_tables,patharray3,npaths3_primary + npaths3_altloc);
	}
	result = Result_paired_as_singles_new(jobid,(void **) patharray5,npaths5_primary,npaths5_altloc,first_absmq5,second_absmq5,
					      (void **) patharray3,npaths3_primary,npaths3_altloc,first_absmq3,second_absmq3);
      }
#endif
      
    } else {
      if (support_tables != NULL) {
	Support_revise_single(support_tables,patharray5,npaths5_primary + npaths5_altloc);
	Support_revise_single(support_tables,patharray3,npaths3_primary + npaths3_altloc);
      }
      result = Result_paired_as_singles_new(jobid,(void **) patharray5,npaths5_primary,npaths5_altloc,first_absmq5,second_absmq5,
					    (void **) patharray3,npaths3_primary,npaths3_altloc,first_absmq3,second_absmq3);
    }

    fp = Output_filestring_fromresult(&(*fp_failedinput),&(*fp_failedinput_1),&(*fp_failedinput_2),
				      result,request,listpool);
    *worker_runtime = worker_stopwatch == NULL ? 0.00 : Stopwatch_stop(worker_stopwatch);
    Result_free(&result,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);

    Stage1_free(&stage1_3,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);
    Stage1_free(&stage1_5,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
		intlistpool,uintlistpool,univcoordlistpool,
		listpool,pathpool,trpathpool,transcriptpool,hitlistpool);

    return fp;
  }
}



#ifdef HAVE_SIGACTION
static const Except_T sigfpe_error = {"SIGFPE--arithmetic exception"};
static const Except_T sigsegv_error = {"SIGSEGV--segmentation violation"};
static const Except_T sigtrap_error = {"SIGTRAP--hardware fault"};
static const Except_T misc_signal_error = {"Miscellaneous signal"};

static void
signal_handler (int sig) {
  Request_T request;
  Shortread_T queryseq1, queryseq2;

  if (sig == SIGPIPE) {
    /* Allow pipe */
    return;
  }

  switch (sig) {
  case SIGABRT: fprintf(stderr,"Signal received: SIGABRT\n"); break;
  case SIGFPE: fprintf(stderr,"Signal received: SIGFPE\n"); break;
  case SIGHUP: fprintf(stderr,"Signal received: SIGHUP\n"); break;
  case SIGILL:
    fprintf(stderr,"Signal received: SIGILL\n");
    fprintf(stderr,"An illegal instruction means that this program is being run on a computer\n");
    fprintf(stderr,"  with different features than the computer used to compile the program\n");
    fprintf(stderr,"You may need to re-compile the program on the same computer type as the target machine\n");
    fprintf(stderr,"  or re-compile with fewer features by doing something like\n");
    fprintf(stderr,"  ./configure --disable-simd\n");
    break;
  case SIGINT: fprintf(stderr,"Signal received: SIGINT\n"); break;
  case SIGPIPE: fprintf(stderr,"Signal received: SIGPIPE\n"); break;
  case SIGQUIT: fprintf(stderr,"Signal received: SIGQUIT\n"); break;
  case SIGSEGV: fprintf(stderr,"Signal received: SIGSEGV\n"); break;
  case SIGSYS: fprintf(stderr,"Signal received: SIGSYS\n"); break;
  case SIGTERM: fprintf(stderr,"Signal received: SIGTERM\n"); break;
  case SIGTRAP: fprintf(stderr,"Signal received: SIGTRAP\n"); break;
  case SIGXCPU: fprintf(stderr,"Signal received: SIGXCPU\n"); break;
  case SIGXFSZ: fprintf(stderr,"Signal received: SIGXFSZ\n"); break;
  }

  Access_emergency_cleanup();


#ifdef HAVE_PTHREAD
  request = (Request_T) pthread_getspecific(global_request_key);
  if (request == NULL) {
    /* fprintf(stderr,"Unable to retrieve request for thread\n"); */
  } else {
    queryseq1 = Request_queryseq1(request);
    queryseq2 = Request_queryseq2(request);
    if (queryseq1 == NULL) {
      fprintf(stderr,"Unable to retrieve queryseq for request\n");
    } else {
      fprintf(stderr,"Problem sequence: ");
      fprintf(stderr,"%s (%d bp)\n",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      if (queryseq2 == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,queryseq2,invert_first_p,invert_second_p);
      }
    }
  }
#endif

  exit(9);

  return;
}
#endif


/* #define POOL_FREE_INTERVAL 200 */
#define POOL_FREE_INTERVAL 1


static void
single_thread () {
  Request_T request;

  Chrpos_T positive_gap_distance;
  int max_nintersections, max_localdb_nregions;

  int streamspace_max_alloc;
  Univcoord_T *streamspace_alloc;
  Univdiag_T *novel_univdiags_alloc;
  Univcoord_T *novel_diagonals_alloc;
  unsigned short *localdb_alloc;

  Filestring_T fp, fp_failedinput, fp_failedinput_1, fp_failedinput_2;
  Shortread_T queryseq1;
  Stopwatch_T worker_stopwatch;

  /* For GMAP */
  /* Oligoindex_array_T oligoindices_major, oligoindices_minor; */
  /* Dynprog_T dynprogL, dynprogM, dynprogR; */
  Trdiagpool_T trdiagpool;
  Univdiagpool_T univdiagpool;
  Auxinfopool_T auxinfopool;
  Pathstorepool_T pathstorepool;
  Intlistpool_T intlistpool;
  Uintlistpool_T uintlistpool;
  Univcoordlistpool_T univcoordlistpool;
  Listpool_T listpool;
  Trpathpool_T trpathpool;
  Pathpool_T pathpool;
  Hitlistpool_T hitlistpool;
  Transcriptpool_T transcriptpool;
  Vectorpool_T vectorpool;
  Spliceendsgen_T spliceendsgen, spliceendsgen5, spliceendsgen3;
  int jobid = 0;
  double worker_runtime;

#ifdef MEMUSAGE
  long int overall_max = 0, memusage_constant = 0, memusage;
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20], comma9[20];
#endif

  /* Intersection indices are paired */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  if (localdb == NULL) {
    max_localdb_nregions = 0;
    max_nintersections = 0;

    streamspace_max_alloc = 0;
    streamspace_alloc = (Univcoord_T *) NULL;

    novel_univdiags_alloc = (Univdiag_T *) NULL;
    novel_diagonals_alloc = (Univcoord_T *) NULL;
    localdb_alloc = (unsigned short *) NULL;

  } else {
    /* Used by Localdb_get */
    /* Because we re-use these, use the keep memory pool */
    max_localdb_nregions = (positive_gap_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
    max_nintersections = max_localdb_nregions * LOCALDB_REGION_SIZE;

    streamspace_max_alloc = max_localdb_nregions * LOCALDB_REGION_SIZE; /* same as max_nintersections */
    streamspace_alloc = (Univcoord_T *) MALLOC_KEEP_ALIGN(streamspace_max_alloc * sizeof(Univcoord_T));

    novel_univdiags_alloc = Univdiag_new_array_keep(2 * max_nintersections);
    novel_diagonals_alloc = (Univcoord_T *) MALLOC_KEEP(2 * max_nintersections * sizeof(Univcoord_T));
    localdb_alloc = (unsigned short *) MALLOC_KEEP_ALIGN(max_nintersections * sizeof(unsigned short));
  }

  trdiagpool = Trdiagpool_new(); Trdiagpool_init(trdiagpool);
  univdiagpool = Univdiagpool_new(); Univdiagpool_init(univdiagpool);
  auxinfopool = Auxinfopool_new(); Auxinfopool_init(auxinfopool);
  pathstorepool = Pathstorepool_new(); Pathstorepool_init(pathstorepool);
  intlistpool = Intlistpool_new(); Intlistpool_init(intlistpool);
  uintlistpool = Uintlistpool_new(); Uintlistpool_init(uintlistpool);
  univcoordlistpool = Univcoordlistpool_new(); Univcoordlistpool_init(univcoordlistpool);
  listpool = Listpool_new(); Listpool_init(listpool);
  trpathpool = Trpathpool_new(); Trpathpool_init(trpathpool);
  pathpool = Pathpool_new(); Pathpool_init(pathpool);
  hitlistpool = Hitlistpool_new(); Hitlistpool_init(hitlistpool);
  transcriptpool = Transcriptpool_new(); Transcriptpool_init(transcriptpool);
  vectorpool = Vectorpool_new(); Vectorpool_init(vectorpool);

  spliceendsgen = Spliceendsgen_new();
  spliceendsgen5 = Spliceendsgen_new();
  spliceendsgen3 = Spliceendsgen_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  /* Except_stack_create(); -- requires pthreads */

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  /* Mem_usage_reset_heap_baseline(0); */
#endif


  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(printf("single_thread got request %d\n",Request_id(request)));

#ifdef MEMUSAGE
    queryseq1 = Request_queryseq1(request);
    /* fprintf(stderr,"Single thread starting %s\n",Shortread_accession(queryseq1)); */
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      if (pass == PASS1) {
	process_request_pass1(request,streamspace_max_alloc,streamspace_alloc,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      trdiagpool,univdiagpool,auxinfopool,pathstorepool,
			      intlistpool,uintlistpool,univcoordlistpool,
			      listpool,trpathpool,pathpool,vectorpool,hitlistpool,
			      transcriptpool,spliceendsgen,spliceendsgen5,spliceendsgen3);
      } else if (pass == PASS2) {
	fp = process_request_pass2(&fp_failedinput,&fp_failedinput_1,&fp_failedinput_2,&worker_runtime,
				   request,streamspace_max_alloc,streamspace_alloc,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,vectorpool,hitlistpool,
				   transcriptpool,spliceendsgen,spliceendsgen5,spliceendsgen3,
				   worker_stopwatch);
      } else {
	fprintf(stderr,"Unknown pass %d\n",pass);
	abort();
      }

      if (timingp == true) {
        queryseq1 = Request_queryseq1(request);
	/* Report time in microseconds */
        fprintf(stderr,"%s\t%.6f\n",Shortread_accession(queryseq1),worker_runtime * 1e6);
      }

    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,Request_queryseq2(request),
					       invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    /* Only a single thread, but this way Outbuffer_print_filestrings can assume Filestring_stringify has been called */
    if (pass == PASS1) {
      /* Nothing to do */
    } else if (pass == PASS2) {
      Filestring_stringify(fp);
      if (fp_failedinput != NULL) {
	Filestring_stringify(fp_failedinput);
      }
      if (fp_failedinput_1 != NULL) {
	Filestring_stringify(fp_failedinput_1);
      }
      if (fp_failedinput_2 != NULL) {
	Filestring_stringify(fp_failedinput_2);
      }
      Outbuffer_print_filestrings(fp,fp_failedinput,fp_failedinput_1,fp_failedinput_2);
    }

    if (jobid % POOL_FREE_INTERVAL == 0) {
      Trdiagpool_reset_memory(trdiagpool);
      Univdiagpool_reset_memory(univdiagpool);
      Auxinfopool_reset_memory(auxinfopool);
      Pathstorepool_reset_memory(pathstorepool);
      Intlistpool_reset_memory(intlistpool);
      Uintlistpool_reset_memory(uintlistpool);
      Univcoordlistpool_reset_memory(univcoordlistpool);
      Listpool_reset_memory(listpool);
      Trpathpool_reset_memory(trpathpool);
      Pathpool_reset_memory(pathpool);
      Hitlistpool_reset_memory(hitlistpool);
      Transcriptpool_reset_memory(transcriptpool);
      Vectorpool_reset_memory(vectorpool);

      Spliceendsgen_free_memory(spliceendsgen);
      Spliceendsgen_free_memory(spliceendsgen5);
      Spliceendsgen_free_memory(spliceendsgen3);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq1 = Request_queryseq1(request);
    strncpy(acc,Shortread_accession(queryseq1),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    if (Mem_usage_report_std_heap_max() > overall_max) {
      overall_max = Mem_usage_report_std_heap_max();
    }
    Genomicpos_commafmt_fill(comma9,overall_max);
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s: constant %s  overall_max %s  query_max %s  std %s  keep %s  in %s  out %s\n",
	    acc,comma0,comma9,comma1,comma2,comma3,comma4,comma5);

    /* We check only the std memory pool, and not the keep memory pool, which can be used across queries */
    if ((memusage = Mem_usage_report_std_heap()) > memusage_constant) {
      fprintf(stderr,"Acc %s: Memory leak in single thread of %ld bytes\n",acc,memusage);
      fflush(stdout);
      /* exit(9); */
    } else if (Mem_usage_report_std_heap_max() > ONE_GB) {
      fprintf(stderr,"Acc %s: Excessive memory usage in single thread of %s bytes\n",acc,comma1);
      fflush(stdout);
      /* exit(9); */
    }
#endif
  }


#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant,__FILE__,__LINE__);
#endif

  /* Except_stack_destroy(); -- requires pthreads */

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }

  Spliceendsgen_free(&spliceendsgen3);
  Spliceendsgen_free(&spliceendsgen5);
  Spliceendsgen_free(&spliceendsgen);

  Trpathpool_free(&trpathpool);
  Pathpool_free(&pathpool);
  Hitlistpool_free(&hitlistpool);
  Transcriptpool_free(&transcriptpool);
  Vectorpool_free(&vectorpool);
  Listpool_free(&listpool);
  Univcoordlistpool_free(&univcoordlistpool);
  Uintlistpool_free(&uintlistpool);
  Intlistpool_free(&intlistpool);
  Auxinfopool_free(&auxinfopool);
  Pathstorepool_free(&pathstorepool);
  Univdiagpool_free(&univdiagpool);
  Trdiagpool_free(&trdiagpool);

  FREE_KEEP_ALIGN(localdb_alloc);
  FREE_KEEP(novel_diagonals_alloc);
  Univdiag_free_array_keep(&novel_univdiags_alloc);
  FREE_KEEP_ALIGN(streamspace_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  return;
}


#ifdef HAVE_PTHREAD
static void *
worker_thread (void *data) {
  Request_T request;

  Chrpos_T positive_gap_distance;
  int max_nintersections, max_localdb_nregions;

  int streamspace_max_alloc;
  Univcoord_T *streamspace_alloc;
  Univdiag_T *novel_univdiags_alloc;
  Univcoord_T *novel_diagonals_alloc;
  unsigned short *localdb_alloc;

  Filestring_T fp, fp_failedinput, fp_failedinput_1, fp_failedinput_2;
  Shortread_T queryseq1;
  Stopwatch_T worker_stopwatch;

  Trdiagpool_T trdiagpool;
  Univdiagpool_T univdiagpool;
  Auxinfopool_T auxinfopool;
  Pathstorepool_T pathstorepool;
  Intlistpool_T intlistpool;
  Uintlistpool_T uintlistpool;
  Univcoordlistpool_T univcoordlistpool;
  Listpool_T listpool;
  Trpathpool_T trpathpool;
  Pathpool_T pathpool;
  Hitlistpool_T hitlistpool;
  Transcriptpool_T transcriptpool;
  Vectorpool_T vectorpool;
  Spliceendsgen_T spliceendsgen, spliceendsgen5, spliceendsgen3;

  int worker_jobid = 0;
  double worker_runtime;
#if defined(DEBUG) || defined(MEMUSAGE)
  long int worker_id = (long int) data;
#endif

#ifdef MEMUSAGE
  long int overall_max = 0, memusage_constant = 0, memusage;
  char threadname[12];
  char acc[100+1], comma0[20], comma1[20], comma2[20], comma3[20], comma4[20], comma5[20], comma9[20];
  sprintf(threadname,"thread-%ld",worker_id);
  Mem_usage_set_threadname(threadname);
#endif

  debug(fprintf(stderr,"worker_thread %ld starting\n",worker_id));

  /* Thread-specific data and storage */

  /* Intersection indices are paired */
  positive_gap_distance = (shortsplicedist > (Chrpos_T) max_deletionlen) ? shortsplicedist : (Chrpos_T) max_deletionlen;

  if (localdb == NULL) {
    max_localdb_nregions = 0;
    max_nintersections = 0;

    streamspace_max_alloc = 0;
    streamspace_alloc = (Univcoord_T *) NULL;

    novel_univdiags_alloc = (Univdiag_T *) NULL;
    novel_diagonals_alloc = (Univcoord_T *) NULL;
    localdb_alloc = (unsigned short *) NULL;

  } else {
    /* Used by Localdb_get */
    /* Because we re-use these, use the keep memory pool */
    max_localdb_nregions = (positive_gap_distance + LOCALDB_REGION_SIZE) / LOCALDB_REGION_SIZE + 1;
    max_nintersections = max_localdb_nregions * LOCALDB_REGION_SIZE;

    streamspace_max_alloc = max_localdb_nregions * LOCALDB_REGION_SIZE; /* same as max_nintersections */
    streamspace_alloc = (Univcoord_T *) MALLOC_KEEP_ALIGN(streamspace_max_alloc * sizeof(Univcoord_T));

    novel_univdiags_alloc = Univdiag_new_array_keep(2 * max_nintersections);
    novel_diagonals_alloc = (Univcoord_T *) MALLOC_KEEP(2 * max_nintersections * sizeof(Univcoord_T));
    localdb_alloc = (unsigned short *) MALLOC_KEEP_ALIGN(max_nintersections * sizeof(unsigned short));
  }

  trdiagpool = Trdiagpool_new(); Trdiagpool_init(trdiagpool);
  univdiagpool = Univdiagpool_new(); Univdiagpool_init(univdiagpool);
  auxinfopool = Auxinfopool_new(); Auxinfopool_init(auxinfopool);
  pathstorepool = Pathstorepool_new(); Pathstorepool_init(pathstorepool);
  intlistpool = Intlistpool_new(); Intlistpool_init(intlistpool);
  uintlistpool = Uintlistpool_new(); Uintlistpool_init(uintlistpool);
  univcoordlistpool = Univcoordlistpool_new(); Univcoordlistpool_init(univcoordlistpool);
  listpool = Listpool_new(); Listpool_init(listpool);
  trpathpool = Trpathpool_new(); Trpathpool_init(trpathpool);
  pathpool = Pathpool_new(); Pathpool_init(pathpool);
  hitlistpool = Hitlistpool_new(); Hitlistpool_init(hitlistpool);
  transcriptpool = Transcriptpool_new(); Transcriptpool_init(transcriptpool);
  vectorpool = Vectorpool_new(); Vectorpool_init(vectorpool);

  spliceendsgen = Spliceendsgen_new();
  spliceendsgen5 = Spliceendsgen_new();
  spliceendsgen3 = Spliceendsgen_new();
  worker_stopwatch = (timingp == true) ? Stopwatch_new() : (Stopwatch_T) NULL;

  Except_stack_create();

#ifdef MEMUSAGE
  memusage_constant += Mem_usage_report_std_heap();
  Genomicpos_commafmt_fill(comma0,memusage_constant);
  /* Mem_usage_reset_heap_baseline(0); */
#endif

  while ((request = Inbuffer_get_request(inbuffer)) != NULL) {
    debug(fprintf(stderr,"worker_thread %ld starting request %d (%s)\n",
		  worker_id,Request_id(request),Shortread_accession(Request_queryseq1(request))));
    pthread_setspecific(global_request_key,(void *) request);

#ifdef MEMUSAGE
    queryseq1 = Request_queryseq1(request);
    /* fprintf(stderr,"Thread %d starting %s\n",worker_id,Shortread_accession(queryseq1)); */
    Mem_usage_reset_stack_max();
    Mem_usage_reset_heap_max();
#endif

    TRY
      if (pass == PASS1) {
	process_request_pass1(request,streamspace_max_alloc,streamspace_alloc,
			      novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
			      trdiagpool,univdiagpool,auxinfopool,pathstorepool,
			      intlistpool,uintlistpool,univcoordlistpool,
			      listpool,trpathpool,pathpool,vectorpool,hitlistpool,
			      transcriptpool,spliceendsgen,spliceendsgen5,spliceendsgen3);
      } else if (pass == PASS2) {
	fp = process_request_pass2(&fp_failedinput,&fp_failedinput_1,&fp_failedinput_2,&worker_runtime,
				   request,streamspace_max_alloc,streamspace_alloc,
				   novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				   trdiagpool,univdiagpool,auxinfopool,pathstorepool,
				   intlistpool,uintlistpool,univcoordlistpool,
				   listpool,trpathpool,pathpool,vectorpool,hitlistpool,
				   transcriptpool,spliceendsgen,spliceendsgen5,spliceendsgen3,
				   worker_stopwatch);
      } else {
	fprintf(stderr,"Unknown pass %d\n",pass);
	abort();
      }

      if (timingp == true) {
        queryseq1 = Request_queryseq1(request);
	/* Report time in microseconds */
        fprintf(stderr,"%s\t%.6f\n",Shortread_accession(queryseq1),worker_runtime * 1e6);
      }

    ELSE
      queryseq1 = Request_queryseq1(request);
      if (queryseq1 == NULL) {
	fprintf(stderr,"NULL");
      } else if (Shortread_accession(queryseq1) == NULL) {
	fprintf(stderr,"unnamed (%d bp)",Shortread_fulllength(queryseq1));
      } else {
	fprintf(stderr,"Problem sequence: ");
	fprintf(stderr,"%s (%d bp)",Shortread_accession(queryseq1),Shortread_fulllength(queryseq1));
      }
      fprintf(stderr,"\n");
      if (Request_queryseq2(request) == NULL) {
	Shortread_stderr_query_singleend_fasta(queryseq1,/*headerseq*/queryseq1);
      } else {
	Shortread_stderr_query_pairedend_fasta(queryseq1,Request_queryseq2(request),
					       invert_first_p,invert_second_p);
      }
      fprintf(stderr,"\n");
      fprintf(stderr,"To obtain a core dump, re-run program on problem sequence with the -0 [zero] flag\n");

      fprintf(stderr,"Exiting...\n");
      exit(9);
    RERAISE;
    END_TRY;

    debug(fprintf(stderr,"worker_thread %ld finished with request %d (%s)\n",
		  worker_id,Request_id(request),Shortread_accession(Request_queryseq1(request))));

    if (pass == PASS1) {
      Outbuffer_put_pass1(outbuffer,Request_id(request));

    } else if (pass == PASS2) {
      /* Parallelize the stringify operation by performing by worker thread and not the output thread */
      Filestring_stringify(fp);
      if (fp_failedinput != NULL) {
	Filestring_stringify(fp_failedinput);
      }
      if (fp_failedinput_1 != NULL) {
	Filestring_stringify(fp_failedinput_1);
      }
      if (fp_failedinput_2 != NULL) {
	Filestring_stringify(fp_failedinput_2);
      }

      Outbuffer_put_filestrings(outbuffer,Request_id(request),fp,fp_failedinput,
				fp_failedinput_1,fp_failedinput_2);
    }

    if (worker_jobid % POOL_FREE_INTERVAL == 0) {
      Trdiagpool_reset_memory(trdiagpool);
      Univdiagpool_reset_memory(univdiagpool);
      Auxinfopool_reset_memory(auxinfopool);
      Pathstorepool_reset_memory(pathstorepool);
      Intlistpool_reset_memory(intlistpool);
      Uintlistpool_reset_memory(uintlistpool);
      Univcoordlistpool_reset_memory(univcoordlistpool);
      Listpool_reset_memory(listpool);
      Trpathpool_reset_memory(trpathpool);
      Pathpool_reset_memory(pathpool);
      Hitlistpool_reset_memory(hitlistpool);
      Transcriptpool_reset_memory(transcriptpool);
      Vectorpool_reset_memory(vectorpool);

      Spliceendsgen_free_memory(spliceendsgen);
      Spliceendsgen_free_memory(spliceendsgen5);
      Spliceendsgen_free_memory(spliceendsgen3);
    }

#ifdef MEMUSAGE
    /* Copy acc before we free the request */
    queryseq1 = Request_queryseq1(request);
    strncpy(acc,Shortread_accession(queryseq1),100);
    acc[100] = '\0';
#endif

    Request_free(&request);

#ifdef MEMUSAGE
    if (Mem_usage_report_std_heap_max() > overall_max) {
      overall_max = Mem_usage_report_std_heap_max();
    }
    Genomicpos_commafmt_fill(comma9,overall_max);
    Genomicpos_commafmt_fill(comma1,Mem_usage_report_std_heap_max());
    Genomicpos_commafmt_fill(comma2,Mem_usage_report_std_heap());
    Genomicpos_commafmt_fill(comma3,Mem_usage_report_keep());
    Genomicpos_commafmt_fill(comma4,Mem_usage_report_in());
    Genomicpos_commafmt_fill(comma5,Mem_usage_report_out());

    fprintf(stderr,"Acc %s, thread %d: constant %s  overall_max %s  query_max %s  std %s  keep %s  in %s  out %s\n",
	    acc,worker_id,comma0,comma9,comma1,comma2,comma3,comma4,comma5);

    /* We check only the std memory pool, and not the keep memory pool, which can be used across queries */
    if ((memusage = Mem_usage_report_std_heap()) > memusage_constant) {
      fprintf(stderr,"Acc %s: Memory leak in worker thread %ld of %ld bytes\n",acc,worker_id,memusage);
      fflush(stdout);
      /* exit(9); */
    } else if (Mem_usage_report_std_heap_max() > ONE_GB) {
      fprintf(stderr,"Acc %s: Excessive memory usage in worker thread %ld of %s bytes\n",acc,worker_id,comma1);
      fflush(stdout);
      /* exit(9); */
    }
#endif
  }

#ifdef MEMUSAGE
  Mem_usage_std_heap_add(memusage_constant,__FILE__,__LINE__);
#endif

  Except_stack_destroy();

  if (worker_stopwatch != NULL) {
    Stopwatch_free(&worker_stopwatch);
  }

  Spliceendsgen_free(&spliceendsgen3);
  Spliceendsgen_free(&spliceendsgen5);
  Spliceendsgen_free(&spliceendsgen);

  Trpathpool_free(&trpathpool);
  Pathpool_free(&pathpool);
  Hitlistpool_free(&hitlistpool);
  Transcriptpool_free(&transcriptpool);
  Vectorpool_free(&vectorpool);
  Listpool_free(&listpool);
  Univcoordlistpool_free(&univcoordlistpool);
  Uintlistpool_free(&uintlistpool);
  Intlistpool_free(&intlistpool);
  Auxinfopool_free(&auxinfopool);
  Pathstorepool_free(&pathstorepool);
  Univdiagpool_free(&univdiagpool);
  Trdiagpool_free(&trdiagpool);

  FREE_KEEP_ALIGN(localdb_alloc);
  FREE_KEEP(novel_diagonals_alloc);
  Univdiag_free_array_keep(&novel_univdiags_alloc);
  FREE_KEEP_ALIGN(streamspace_alloc);

#ifdef MEMUSAGE
  Mem_usage_set_threadname("main");
#endif

  debug(fprintf(stderr,"worker_thread %ld finished\n",worker_id));

  return (void *) NULL;
}
#endif


static void
parse_part (unsigned int *part_modulus, unsigned int *part_interval, char *string) {
  char *p = string;

  if (sscanf(p,"%u",&(*part_modulus)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  while (*p != '\0' && !isdigit(*p)) {
    p++;
  }
  if (sscanf(p,"%d",&(*part_interval)) < 1) {
    fprintf(stderr,"Cannot parse first integer from %s\n",string);
    exit(9);
  }
  if ((*part_modulus) >= (*part_interval)) {
    fprintf(stderr,"In %s, batch number %u must be less than the number of batches %u\n",
	    string,*part_modulus,*part_interval);
    exit(9);
  }
  if (*part_interval == 0) {
    fprintf(stderr,"Bad batch specification %s.  Batch interval cannot be 0.\n",string);
    exit(9);
  }

  return;
}


#if 0
static int
add_gmap_mode (char *string) {
  if (!strcmp(string,"none")) {
    gmap_mode = 0;
    return 0;
  } else if (!strcmp(string,"all")) {
    gmap_mode = (GMAP_IMPROVEMENT | GMAP_ENDS | GMAP_PAIRSEARCH);
    return 1;
  } else {
    if (!strcmp(string,"improve")) {
      gmap_mode |= GMAP_IMPROVEMENT;
    } else if (!strcmp(string,"ends")) {
      gmap_mode |= GMAP_ENDS;
    } else if (!strcmp(string,"pairsearch")) {
      gmap_mode |= GMAP_PAIRSEARCH;
    } else {
      fprintf(stderr,"Don't recognize gmap-mode type %s\n",string);
      fprintf(stderr,"Allowed values are: none, all, improve, ends, pairsearch\n");
      exit(9);
    }
    return 1;
  }
}
#endif


static char *
check_valid_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+') {
      p++;
    }
    if (!isdigit(*p)) {
      return false;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid int\n",string);
    exit(9);
    return NULL;
  }
}


static double
check_valid_float (char *string, const char *option) {
  double value;
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
      exit(9);
      return 0.0;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    if ((value = atof(string)) > 1.0 || value < 0.0) {
      fprintf(stderr,"Value for option %s should be between 0.0 and 1.0\n",option);
      exit(9);
    } else {
      return value;
    }
  } else {
    fprintf(stderr,"Value %s for option %s is not a valid float\n",string,option);
    exit(9);
    return 0.0;
  }
}


static char *
check_valid_float_or_int (char *string) {
  char *p = string;

  if (*p == '+' || *p == '-') {
    p++;
  }

  while (*p != '\0' && isdigit(*p)) {
    p++;
  }
  if (*p == '\0') {
    return string;
  }

  if (*p == '.') {
    p++;
  }

  if (!isdigit(*p)) {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
  while (*p != '\0' && isdigit(*p)) {
    p++;
  }

  if (*p == 'e') {
    p++;
    if (*p == '+' || *p == '-') {
      p++;
    }
    if (!isdigit(*p)) {
      fprintf(stderr,"value %s is not a valid float\n",string);
      exit(9);
      return NULL;
    }
    while (*p != '\0' && isdigit(*p)) {
      p++;
    }
  }

  if (*p == '\0') {
    return string;
  } else {
    fprintf(stderr,"value %s is not a valid float\n",string);
    exit(9);
    return NULL;
  }
}


static int
parse_command_line (int argc, char *argv[], int optind) {
  int opt, c;
  extern char *optarg;
  int long_option_index = 0;
  const char *long_name;
  char **argstart;
  /* char *string; */

  fprintf(stderr,"GSNAP version %s called with args:",PACKAGE_VERSION);
  argstart = &(argv[-optind]);
  for (c = 1; c < argc + optind; c++) {
      fprintf(stderr," %s",argstart[c]);
  }
  fprintf(stderr,"\n");

  while ((opt = getopt_long(argc,argv,
			    "C:c:D:d:k:q:o:a:N:M:m:w:e:J:s:V:v:B:t:A:j:0n:QO",
			    long_options, &long_option_index)) != -1) {
    switch (opt) {
    case 0:
      long_name = long_options[long_option_index].name;
      if (!strcmp(long_name,"version")) {
	print_program_version();
	return 1;
      } else if (!strcmp(long_name,"check")) {
	check_compiler_assumptions();
	return 1;
      } else if (!strcmp(long_name,"help")) {
	print_program_usage();
	return 1;

      } else if (!strcmp(long_name,"use-sarray")) {
	fprintf(stderr,"Ignoring the --use-sarray flag, since this version of GSNAP does not use suffix arrays\n");

      } else if (!strcmp(long_name,"terminal-threshold")) {
	fprintf(stderr,"Ignoring the --terminal-threshold flag, which is obsolete\n");

      } else if (!strcmp(long_name,"two-pass")) {
	two_pass_p = true;

      } else if (!strcmp(long_name,"use-localdb")) {
	if (!strcmp(optarg,"1")) {
	  user_localdb_p = true;
	  use_localdb_p = true;
	} else if (!strcmp(optarg,"0")) {
	  user_localdb_p = true;
	  use_localdb_p = false;
	} else {
	  fprintf(stderr,"--use-localdb flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"transcriptome-mode")) {
	if (!strcmp(optarg,"assist")) {
	  genome_align_p = true;
	  transcriptome_align_p = true;
	} else if (!strcmp(optarg,"only")) {
	  genome_align_p = false;
	  transcriptome_align_p = true;
	} else if (!strcmp(optarg,"annotate")) {
	  genome_align_p = true;
	  transcriptome_align_p = false;
	} else {
	  fprintf(stderr,"--transcriptome-mode flag must be assist, only, or annotate\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"use-shared-memory")) {
	if (!strcmp(optarg,"1")) {
	  sharedp = true;
	} else if (!strcmp(optarg,"0")) {
	  sharedp = false;
	} else {
	  fprintf(stderr,"--use-shared-memory flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"preload-shared-memory")) {
	preload_shared_memory_p = true;

      } else if (!strcmp(long_name,"unload-shared-memory")) {
	unload_shared_memory_p = true;

      } else if (!strcmp(long_name,"expand-offsets")) {
	fprintf(stderr,"Note: --expand-offsets flag is no longer supported.  With the latest algorithms, it doesn't improve speed much.  Ignoring this flag");

      } else if (!strcmp(long_name,"sampling")) {
	required_index1interval = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"align-fraction")) {
	skip_prob = 1.0 - check_valid_float(optarg,long_name);

      } else if (!strcmp(long_name,"time")) {
	timingp = true;

      } else if (!strcmp(long_name,"unload")) {
	unloadp = true;

      } else if (!strcmp(long_name,"maxsearch")) {
	maxpaths_search = atoi(check_valid_int(optarg)); break;

      } else if (!strcmp(long_name,"mode")) {
	if (!strcmp(optarg,"standard")) {
	  mode = STANDARD;
	} else if (!strcmp(optarg,"cmet-stranded")) {
	  mode = CMET_STRANDED;
	} else if (!strcmp(optarg,"cmet-nonstranded")) {
	  mode = CMET_NONSTRANDED;
	} else if (!strcmp(optarg,"atoi-stranded")) {
	  mode = ATOI_STRANDED;
	} else if (!strcmp(optarg,"atoi-nonstranded")) {
	  mode = ATOI_NONSTRANDED;
	} else if (!strcmp(optarg,"ttoc-stranded")) {
	  mode = TTOC_STRANDED;
	} else if (!strcmp(optarg,"ttoc-nonstranded")) {
	  mode = TTOC_NONSTRANDED;
	} else {
	  fprintf(stderr,"--mode must be standard, cmet-stranded, cmet-nonstranded, atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"cmetdir")) {
	user_modedir = optarg;

      } else if (!strcmp(long_name,"atoidir")) {
	user_modedir = optarg;

      } else if (!strcmp(long_name,"splicingdir")) {
	user_splicingdir = optarg;

      } else if (!strcmp(long_name,"splices-noeval")) {
	splice_evalp = false;

      } else if (!strcmp(long_name,"splices-dump")) {
	dump_splices_file = optarg;

      } else if (!strcmp(long_name,"splices-exclude-known")) {
	dump_splices_include_known_p = false;

      } else if (!strcmp(long_name,"splices-read")) {
	read_splices_file = optarg;

      } else if (!strcmp(long_name,"indels-dump")) {
	if ((dump_indels_fp = fopen(optarg,"w")) == NULL) {
	  fprintf(stderr,"Cannot write to --indels-dump file %s\n",optarg);
	  exit(9);
	}

      } else if (!strcmp(long_name,"indels-read")) {
	indels_file = optarg;

      } else if (!strcmp(long_name,"find-dna-chimeras")) {
	if (!strcmp(optarg,"1")) {
	  user_find_dna_chimeras_p = true;
	  find_dna_chimeras_p = true;
	} else if (!strcmp(optarg,"0")) {
	  user_find_dna_chimeras_p = true;
	  find_dna_chimeras_p = false;
	} else {
	  fprintf(stderr,"--find-dna-chimeras flag must be 0 or 1\n");
	  exit(9);
	}

      } else if (!strcmp(long_name,"tallydir")) {
	user_tallydir = optarg;

      } else if (!strcmp(long_name,"use-tally")) {
	tally_root = optarg;

      } else if (!strcmp(long_name,"runlengthdir")) {
	user_runlengthdir = optarg;

      } else if (!strcmp(long_name,"use-runlength")) {
	runlength_root = optarg;

      } else if (!strcmp(long_name,"input-buffer-size")) {
	input_buffer_size = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"output-buffer-size")) {
	output_buffer_size = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"barcode-length")) {
	barcode_length = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"endtrim-length")) {
	endtrim_length = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"fastq-id-start")) {
	acc_fieldi_start = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_start < 0) {
	  fprintf(stderr,"Value for fastq-id-start must be 1 or greater\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"fastq-id-end")) {
	acc_fieldi_end = atoi(check_valid_int(optarg)) - 1;
	if (acc_fieldi_end < 0) {
	  fprintf(stderr,"Value for fastq-id-end must be 1 or greater\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"force-single-end")) {
	force_single_end_p = true;

      } else if (!strcmp(long_name,"filter-chastity")) {
	if (!strcmp(optarg,"off")) {
	  filter_chastity_p = false;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"either")) {
	  filter_chastity_p = true;
	  filter_if_both_p = false;
	} else if (!strcmp(optarg,"both")) {
	  filter_chastity_p = true;
	  filter_if_both_p = true;
	} else {
	  fprintf(stderr,"--filter-chastity values allowed: off, either, both\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"allow-pe-name-mismatch")) {
	allow_paired_end_mismatch_p = true;

      } else if (!strcmp(long_name,"10x-whitelist")) {
	whitelist_file = optarg;

      } else if (!strcmp(long_name,"10x-well-position")) {
	wellpos = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"read-files-command")) {
	read_files_command = optarg;

#ifdef HAVE_ZLIB
      } else if (!strcmp(long_name,"gunzip")) {
	gunzip_p = true;
#endif

#ifdef HAVE_BZLIB
      } else if (!strcmp(long_name,"bunzip2")) {
	bunzip2_p = true;
#endif

      } else if (!strcmp(long_name,"interleaved")) {
	interleavedp = true;

      } else if (!strcmp(long_name,"orientation")) {
	if (!strcmp(optarg,"FR")) {
	  invert_first_p = false;
	  invert_second_p = true;

	} else if (!strcmp(optarg,"RF")) {
	  invert_first_p = true;
	  invert_second_p = false;

	} else if (!strcmp(optarg,"FF")) {
	  invert_first_p = invert_second_p = false;

	} else if (!strcmp(optarg,"10X") || !strcmp(optarg,"10x")) {
	  single_cell_p = true;
	  invert_first_p = false;
	  invert_second_p = true;
	  keep_chastity_p = false;

	} else {
	  fprintf(stderr,"Currently allowed values for orientation: FR (fwd-rev), RF (rev-fwd), FF (fwd-fwd), or 10x (barcodes-rev)\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"split-output")) {
	split_output_root = optarg;

      } else if (!strcmp(long_name,"split-simple")) {
	split_simple_p = true;

      } else if (!strcmp(long_name,"failed-input")) {
	failedinput_root = optarg;

      } else if (!strcmp(long_name,"append-output")) {
	appendp = true;

      } else if (!strcmp(long_name,"order-among-best")) {
	if (!strcmp(optarg,"genomic")) {
	  want_random_p = false;
	} else if (!strcmp(optarg,"random")) {
	  want_random_p = true;
	} else {
	  fprintf(stderr,"--order-among-best values allowed: genomic, random (default)\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"pairmax-dna")) {
	pairmax_dna = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairmax-rna")) {
	pairmax_rna = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"resolve-inner")) {
	if (!strcmp(optarg,"1")) {
	  resolve_inner_p = true;
	} else if (!strcmp(optarg,"0")) {
	  resolve_inner_p = false;
	} else {
	  fprintf(stderr,"--resolve-inner flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"pairexpect")) {
	expected_pairlength = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"pairdev")) {
	pairlength_deviation = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"min-coverage")) {
	user_mincoverage_filter_float = atof(check_valid_float_or_int(optarg));
	if (user_mincoverage_filter_float > 1.0 && user_mincoverage_filter_float != rint(user_mincoverage_filter_float)) {
	  fprintf(stderr,"Cannot specify fractional value %f for --min-coverage except between 0.0 and 1.0\n",
		  user_mincoverage_filter_float);
	  return 9;
	}

      } else if (!strcmp(long_name,"pass1-min-support")) {
	pass1_min_support = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"indel-endlength")) {
	min_indel_end_matches = atoi(check_valid_int(optarg));
#if 0
	if (min_indel_end_matches > 14) {
	  allow_end_indels_p = false;
	}
#endif

      } else if (!strcmp(long_name,"max-insertions")) {
	max_insertionlen = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"max-deletions")) {
	max_deletionlen = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"merge-distant-samechr")) {
	merge_samechr_p = true;

      } else if (!strcmp(long_name,"query-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  query_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  query_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--query-unk-mismatch flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"genome-unk-mismatch")) {
	if (!strcmp(optarg,"1")) {
	  genome_unk_mismatch_p = true;
	} else if (!strcmp(optarg,"0")) {
	  genome_unk_mismatch_p = false;
	} else {
	  fprintf(stderr,"--genome-unk-mismatch flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"trim-mismatch-score")) {
	fprintf(stderr,"--trim-mismatch-score is now deprecated.  GSNAP is now tuned to trim well at the ends\n");
	/* trim_mismatch_score = atoi(check_valid_int(optarg)); */
	/* user_trim_mismatch_score_p = true; */

      } else if (!strcmp(long_name,"force-xs-dir")) {
	force_xs_direction_p = true;

      } else if (!strcmp(long_name,"show-refdiff")) {
	show_refdiff_p = true;

      } else if (!strcmp(long_name,"clip-overlap")) {
	clip_overlap_p = true;

      } else if (!strcmp(long_name,"merge-overlap")) {
	merge_overlap_p = true;

      } else if (!strcmp(long_name,"show-method")) {
	method_print_p = true;

      } else if (!strcmp(long_name,"show-univdiagonal")) {
	print_univdiagonal_p = true;

      } else if (!strcmp(long_name,"no-sam-headers")) {
	sam_headers_p = false;

      } else if (!strcmp(long_name,"add-paired-nomappers")) {
	add_paired_nomappers_p = true;

      } else if (!strcmp(long_name,"paired-flag-means-concordant")) {
	if (!strcmp(optarg,"1")) {
	  paired_flag_means_concordant_p = true;
	} else if (!strcmp(optarg,"0")) {
	  paired_flag_means_concordant_p = false; /* Default */
	} else {
	  fprintf(stderr,"--paired-flag-means-concordant flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"sam-headers-batch")) {
	sam_headers_batch = atoi(check_valid_int(optarg));

      } else if (!strcmp(long_name,"sam-hardclip-use-S")) {
	sam_hardclip_use_S_p = true;

      } else if (!strcmp(long_name,"sam-use-0M")) {
	if (!strcmp(optarg,"1")) {
	  sam_insert_0M_p = true;
	} else if (!strcmp(optarg,"0")) {
	  sam_insert_0M_p = false;
	} else {
	  fprintf(stderr,"--sam-use-0M flag must be 0 or 1\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"sam-extended-cigar")) {
	sam_cigar_extended_p = true;

      } else if (!strcmp(long_name,"sam-multiple-primaries")) {
	sam_multiple_primaries_p = true;

      } else if (!strcmp(long_name,"sam-sparse-secondaries")) {
	sam_sparse_secondaries_p = true;

      } else if (!strcmp(long_name,"quality-protocol")) {
	if (user_quality_score_adj == true) {
	  fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	  return 9;
	} else if (user_quality_shift == true) {
	  fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	  return 9;
	} else if (!strcmp(optarg,"illumina")) {
	  MAPQ_init(/*quality_score_adj*/64);
	  user_quality_score_adj = true;
	  quality_shift = -31;
	  user_quality_shift = true;
	} else if (!strcmp(optarg,"sanger")) {
	  MAPQ_init(/*quality_score_adj*/33);
	  user_quality_score_adj = true;
	  quality_shift = 0;
	  user_quality_shift = true;
	} else {
	  fprintf(stderr,"The only values allowed for --quality-protocol are illumina or sanger\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"md-report-snps")) {
	md_report_snps_p = true;

      } else if (!strcmp(long_name,"no-soft-clips")) {
	allow_soft_clips_p = false;

      } else if (!strcmp(long_name,"extend-soft-clips")) {
	extend_soft_clips_p = true;

      } else if (!strcmp(long_name,"action-if-cigar-error")) {
	if (!strcmp(optarg,"ignore")) {
	  cigar_action = CIGAR_ACTION_IGNORE;
	} else if (!strcmp(optarg,"warning")) {
	  cigar_action = CIGAR_ACTION_WARNING;
	} else if (!strcmp(optarg,"noprint")) {
	  cigar_action = CIGAR_ACTION_NOPRINT;
	} else if (!strcmp(optarg,"abort")) {
	  cigar_action = CIGAR_ACTION_ABORT;
	} else {
	  fprintf(stderr,"The only values allowed for --action-if-cigar-error are ignore, warning, noprint, abort\n");
	  return 9;
	}

      } else if (!strcmp(long_name,"read-group-id")) {
	sam_read_group_id = optarg;

      } else if (!strcmp(long_name,"read-group-name")) {
	sam_read_group_name = optarg;

      } else if (!strcmp(long_name,"read-group-library")) {
	sam_read_group_library = optarg;

      } else if (!strcmp(long_name,"read-group-platform")) {
	sam_read_group_platform = optarg;

      } else if (!strcmp(long_name,"print-snps")) {
	print_snplabels_p = true;

      } else if (!strcmp(long_name,"failsonly")) {
	if (nofailsp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  failsonlyp = true;
	}

      } else if (!strcmp(long_name,"nofails")) {
	if (failsonlyp == true) {
	  fprintf(stderr,"Cannot specify both --nofails and --failsonly\n");
	  return 9;
	} else {
	  nofailsp = true;
	}

      } else if (!strcmp(long_name,"chrsubset")) {
	chrsubset_file = optarg;

      } else if (!strcmp(long_name,"only-concordant")) {
	only_concordant_p = true;
      } else if (!strcmp(long_name,"omit-concordant-uniq")) {
	omit_concordant_uniq_p = true;
      } else if (!strcmp(long_name,"omit-concordant-mult")) {
	omit_concordant_mult_p = true;

      } else if (!strcmp(long_name,"omit-softclipped")) {
	omit_softclipped_p = false;

      } else if (!strcmp(long_name,"only-tr-consistent")) {
	only_tr_consistent_p = true;

      } else {
	/* Shouldn't reach here */
	fprintf(stderr,"Don't recognize option %s.  For usage, run 'gsnap --help'",long_name);
	return 9;
      }
      break;

    case 'C': user_transcriptomedir = optarg; break;
    case 'c': transcriptome_dbroot = optarg; break;

    case 'D': user_genomedir = optarg; break;
    case 'd': genome_dbroot = optarg; break;

    case 'k':
      required_index1part = atoi(check_valid_int(optarg));
      if (required_index1part > MAXIMUM_KMER) {
	fprintf(stderr,"The value for k-mer size must be %d or less\n",MAXIMUM_KMER);
	return 9;
      }
      break;

    case 'q': parse_part(&part_modulus,&part_interval,optarg); break;
    case 'o': output_file = optarg; break;

    case 'a': 
      if (!strcmp(optarg,"paired")) {
	chop_primers_p = true;
      } else if (!strcmp(optarg,"off")) {
	chop_primers_p = false;
      } else {
	fprintf(stderr,"Currently allowed values for adapter stripping (-a): off, paired\n");
	return 9;
      }
      break;

    case 'N':
      if (!strcmp(optarg,"1")) {
	novelsplicingp = true;
      } else if (!strcmp(optarg,"0")) {
	novelsplicingp = false;
      } else {
	fprintf(stderr,"Novel splicing (-N flag) must be 0 or 1\n");
	return 9;
      }
      break;

#if 0
    case 'R': 
      if (!strcmp(optarg,"0")) {
	masktype = MASK_NONE;
      } else if (!strcmp(optarg,"1")) {
	masktype = MASK_FREQUENT;
      } else if (!strcmp(optarg,"2")) {
	masktype = MASK_REPETITIVE;
      } else if (!strcmp(optarg,"3")) {
	masktype = MASK_GREEDY_FREQUENT;
      } else if (!strcmp(optarg,"4")) {
	masktype = MASK_GREEDY_REPETITIVE;
      } else {
	fprintf(stderr,"Masking mode %s not recognized.\n",optarg);
	fprintf(stderr,"Mode 0 means no masking, mode 1 masks frequent oligomers;\n");
	fprintf(stderr,"  mode 2 masks frequent and repetitive oligomers;\n");
	fprintf(stderr,"  mode 3 does greedy masking of frequent oligomers,\n");
	fprintf(stderr,"    then no masking if necessary;\n");
	fprintf(stderr,"  mode 4 does greedy masking of frequent and repetitive oligomers,\n");
	fprintf(stderr,"    then no masking if necessary.\n");
	return 9;
      }
      break;
#endif

    case 'M': subopt_levels = atoi(check_valid_int(optarg)); break;
    case 'm':
      user_nmismatches_filter_float = atof(check_valid_float_or_int(optarg));
      if (user_nmismatches_filter_float > 1.0 && user_nmismatches_filter_float != rint(user_nmismatches_filter_float)) {
	fprintf(stderr,"Cannot specify fractional value %f for --max-mismatches except between 0.0 and 1.0\n",user_nmismatches_filter_float);
	return 9;
      }
      break;

    case 'w': shortsplicedist = strtoul(optarg,NULL,10); break;

    case 's':
      splicing_file = optarg;
      knownsplicingp = true;
      break;

    case 'e':
      masked_suffix = optarg;
      genome_unk_mismatch_p = true; /* Want introns to count as mismatches so we can count exonic matches */
      maskedp = true;
      print_nsnpdiffs_p = true;
      break;

    case 'V': user_snpsdir = optarg; break;
    case 'v': snps_root = optarg; break;

    case 'B':
      if (!strcmp(optarg,"5")) {
#if 0
	/* Not true.  -B 5 allocates suffix array and suffix aux files */
	fprintf(stderr,"Note: Batch mode 5 is now the same as batch mode 4.\n");
	fprintf(stderr,"Expansion of offsets is now controlled separately by --expand-offsets (default=0).\n");
#endif
	offsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE; /* Doesn't matter */
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_ALLOCATE;
	localdb_access = USE_ALLOCATE;

#ifdef HAVE_MMAP
      } else if (!strcmp(optarg,"4")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE;
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_ALLOCATE;
	localdb_access = USE_MMAP_ONLY;

      } else if (!strcmp(optarg,"3")) {
	offsetsstrm_access = USE_ALLOCATE;
	positions_access = USE_ALLOCATE;
	locoffsetsstrm_access = USE_ALLOCATE;
	locpositions_access = USE_ALLOCATE;

	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
	localdb_access = USE_MMAP_ONLY;

      } else if (!strcmp(optarg,"2")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */

	genome_access = USE_MMAP_PRELOAD; /* was batch_genome_p = true */
	localdb_access = USE_MMAP_ONLY;

      } else if (!strcmp(optarg,"1")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_PRELOAD; /* was batch_positions_p = true */

	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
	localdb_access = USE_MMAP_ONLY;

      } else if (!strcmp(optarg,"0")) {
	offsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	positions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */
	locoffsetsstrm_access = USE_ALLOCATE; /* was batch_offsets_p = true */
	locpositions_access = USE_MMAP_ONLY; /* was batch_positions_p = false */

	genome_access = USE_MMAP_ONLY; /* was batch_genome_p = false */
	localdb_access = USE_MMAP_ONLY;

#endif
      } else {
#ifdef HAVE_MMAP
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 0-5.  Run 'gsnap --help' for more information.\n",optarg);
#else
	fprintf(stderr,"Batch mode %s not recognized.  Only allow 4-5, since mmap is disabled.  Run 'gsnap --help' for more information.\n",optarg);
#endif
	return 9;
      }
      break;

#if defined(HAVE_PTHREAD)
    case 't': nthreads = atoi(check_valid_int(optarg)); break;
#else
    case 't': fprintf(stderr,"This version of GSNAP has pthreads disabled, so ignoring the value of %s for -t\n",optarg); break;
#endif

    case 'A':
      if (!strcmp(optarg,"standard")) {
	output_type = STD_OUTPUT;
      } else if (!strcmp(optarg,"sam")) {
	output_type = SAM_OUTPUT;
      } else if (!strcmp(optarg,"m8")) {
	output_type = M8_OUTPUT;
      } else {
	fprintf(stderr,"Output format %s not recognized.  Allowed values: standard (default), sam, m8\n",optarg);
	return 9;
      }
      break;

    case 'j':
      if (user_quality_shift == true) {
	fprintf(stderr,"Cannot specify both -j (--quality-print-shift) and --quality-protocol\n");
	return 9;
      } else {
	quality_shift = atoi(check_valid_int(optarg));
	user_quality_shift = true;
      }
      break;

    case 'J':
      if (user_quality_score_adj == true) {
	fprintf(stderr,"Cannot specify both -J (--quality-zero-score) and --quality-protocol\n");
	return 9;
      } else {
	MAPQ_init(/*quality_score_adj*/atoi(check_valid_int(optarg)));
	user_quality_score_adj = true;
      }
      break;

    case '0': exception_raise_p = false; break; /* Allows signals to pass through */
    case 'n': maxpaths_report = atoi(check_valid_int(optarg)); break;
    case 'Q': quiet_if_excessive_p = true; break;

    case 'O': orderedp = true; break;

    case '?': fprintf(stderr,"For usage, run 'gsnap --help'\n"); return 9;
    default: return 9;
    }
  }

  /* Make inferences */
  if (genome_dbroot == NULL) {
    fprintf(stderr,"Need to specify the -d flag.  For usage, run 'gsnap --help'\n");
    /* print_program_usage(); */
    return 9;
  }

  if (single_cell_p == true) {
    chop_poly_at_first_p = false;
    chop_poly_at_second_p = true;
  } else {
    chop_poly_at_first_p = true;
    chop_poly_at_second_p = true;
  }

  if (whitelist_file != NULL && single_cell_p == false) {
    fprintf(stderr,"--10x-whitelist was given, so assuming that orientation is 10x\n");
    single_cell_p = true;
    invert_first_p = false;
    invert_second_p = true;
    keep_chastity_p = false;
  }

  if (acc_fieldi_end < acc_fieldi_start) {
    fprintf(stderr,"--fastq-id-end must be equal to or greater than --fastq-id-start\n");
    return 9;
  }

  if (clip_overlap_p == true && merge_overlap_p == true) {
    fprintf(stderr,"Cannot specify both --clip-overlap and --merge-overlap.  Please choose one.\n");
    return 9;
  }

#if 0
  if (masked_suffix != NULL && genome_unk_mismatch_p == true) {
    fprintf(stderr,"Unexpected combination of --use-mask and --genome-unk-mismatch=true.  Are you sure?\n");
  }
#endif

#if 0
  /* Need to work on precise coordinates in Pathpair_resolve */
  if (expected_pairlength < 5*pairlength_deviation) {
    min_insertlength = 0;
  } else {
    min_insertlength = expected_pairlength - 5*pairlength_deviation;
  }
#else
  min_insertlength = 0;
#endif
  max_insertlength = expected_pairlength + 5*pairlength_deviation;

  if (transcriptome_dbroot == NULL) {
    transcriptome_align_p = false;
  }

  if (transcriptome_dbroot != NULL) {
    knownsplicingp = true;
    if (novelsplicingp == false) {
      fprintf(stderr,"Transcriptome specified => assume reads are RNA-Seq.  Turning on novel splicing.\n");
      novelsplicingp = true;
    } else {
      fprintf(stderr,"Transcriptome specified => assume reads are RNA-Seq\n");
    }
    if (user_find_dna_chimeras_p == false) {
      find_dna_chimeras_p = true;
    }

    pairmax_transcriptome = pairmax_dna;
    pairmax_linear = pairmax_rna;
    pairmax_circular = pairmax_dna;

    if (user_localdb_p == false) {
      use_localdb_p = true;	/* Use for backup when a transcript is not found */
    }

#if 0
  } else if (find_dna_chimeras_p == true) {
    /* DNA-Seq, but need to trim to find chimeras */
    fprintf(stderr,"Neither novel splicing (-N) nor known splicing (-s) turned on => assume reads are DNA-Seq (genomic)\n");
    pairmax_transcriptome = 0;
    pairmax_linear = pairmax_dna;
    pairmax_circular = pairmax_dna;
    shortsplicedist = 0U;
#endif

  } else if (novelsplicingp == true && knownsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) and known splicing (-s) both turned on => assume reads are RNA-Seq\n");
    if (user_find_dna_chimeras_p == false) {
      find_dna_chimeras_p = true;
    }

    pairmax_transcriptome = 0;
    pairmax_linear = pairmax_rna;
    pairmax_circular = pairmax_dna;

    if (user_localdb_p == false) {
      use_localdb_p = default_localdb_p;
    }

  } else if (knownsplicingp == true) {
    fprintf(stderr,"Known splicing (-s) turned on => assume reads are RNA-Seq\n");
    if (user_find_dna_chimeras_p == false) {
      find_dna_chimeras_p = false;
    }

    pairmax_transcriptome = 0;
    pairmax_linear = pairmax_rna;
    pairmax_circular = pairmax_dna;

    if (user_localdb_p == false) {
      use_localdb_p = default_localdb_p;
    }

  } else if (novelsplicingp == true) {
    fprintf(stderr,"Novel splicing (-N) turned on => assume reads are RNA-Seq\n");
    if (user_find_dna_chimeras_p == false) {
      find_dna_chimeras_p = true;
    }

    pairmax_transcriptome = 0;
    pairmax_linear = pairmax_rna;
    pairmax_circular = pairmax_dna;

    if (user_localdb_p == false) {
      use_localdb_p = default_localdb_p;
    }

  } else {
    /* Straight DNA-Seq */
    fprintf(stderr,"Neither novel splicing (-N 1) nor known splicing (-s) turned on => assume reads are DNA-Seq (genomic)\n");
    if (dump_splices_file != NULL) {
      fprintf(stderr,"--splices-dump given, but neither novel splicing (-N 1) nor transcriptome mode (-c) is given\n");
      exit(9);
    }

    if (user_find_dna_chimeras_p == false) {
      find_dna_chimeras_p = true;
    }

    pairmax_transcriptome = 0;
    pairmax_linear = pairmax_dna;
    pairmax_circular = pairmax_dna;
    shortsplicedist = 0;

    if (user_localdb_p == false) {
      use_localdb_p = true;  	/* Previously did not need this for DNA-seq, but need it for end indels */
    }
    use_localdb_p = false;	/* For faster speed */

    /* Genomic reads do not have poly-adenylation */
    chop_poly_at_first_p = false;
    chop_poly_at_second_p = false;

#if 0
    /* Large soft-clips not expected for DNA-Seq, so require good alignments */
    if (user_mismatches_refalt_float < 0.0) {
      user_mismatches_refalt_float = 0.10;
    }
#endif
  }


  if (split_simple_p == true && split_output_root == NULL) {
    fprintf(stderr,"Cannot specify --split-simple without specifying --split-output\n");
    return 9;
  }

  if (sam_headers_batch >= 0) {
    if ((int) part_modulus == sam_headers_batch) {
      sam_headers_p = true;
    } else {
      sam_headers_p = false;
    }
  }

  if (sam_read_group_id == NULL && sam_read_group_name != NULL) {
    sam_read_group_id = sam_read_group_name;
  } else if (sam_read_group_id != NULL && sam_read_group_name == NULL) {
    sam_read_group_name = sam_read_group_id;
  }

  if (chop_primers_p == true) {
    if (invert_first_p == false && invert_second_p == true) {
      /* orientation FR */
    } else {
      fprintf(stderr,"Adapter stripping not currently implemented for given orientation\n");
      return 9;
    }
  }

  if (maxpaths_search < maxpaths_report) {
    fprintf(stderr,"Value for --maxsearch (%d) is less than --npaths (%d), so raising --maxsearch to be %d also\n",
	    maxpaths_search,maxpaths_report,maxpaths_report);
    maxpaths_search = maxpaths_report;
  }

  return 0;
}


static bool
open_input_streams_parser (int *nextchar, int *nchars1, int *nchars2, bool *paired_end_p,
			   char ***files, int *nfiles, FILE **input, FILE **input2,
#ifdef HAVE_ZLIB
			   gzFile *gzipped, gzFile *gzipped2,
#endif
#ifdef HAVE_BZLIB
			   Bzip2_T *bzipped, Bzip2_T *bzipped2,
#endif
			   char *read_files_command, bool gunzip_p, bool bunzip2_p, bool interleavedp,
			   int argc, char **argv) {
  bool fastq_format_p = false;
  
  *input = *input2 = NULL;
  *paired_end_p = true;
#ifdef HAVE_ZLIB
  *gzipped = *gzipped2 = NULL;
#endif
#ifdef HAVE_BZLIB
  *bzipped = *bzipped2 = NULL;
#endif

  /* Open input stream and peek at first char */
  if (preload_shared_memory_p == true || unload_shared_memory_p == true) {
    /* Ignore any input */
    *input = stdin;
    *files = (char **) NULL;
    *nfiles = 0;
    *nextchar = EOF;
    *paired_end_p = false;
    return false;

  } else if (argc == 0) {
    fprintf(stderr,"Reading from stdin\n");
    *input = stdin;
    *files = (char **) NULL;
    *nfiles = 0;
    if (interleavedp == true) {
      *nextchar = '\0';
      /* I believe interleaved is paired-end by definition */
    } else {
      *nextchar = Shortread_input_init(&(*nchars1),*input);
      *paired_end_p = false;
    }

  } else {
    *files = argv;
    *nfiles = argc;

    if (gunzip_p == true) {
#ifdef HAVE_ZLIB
      if ((*gzipped = gzopen((*files)[0],"rb")) == NULL) {
	fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	exit(9);
      } else {
#ifdef HAVE_ZLIB_GZBUFFER
	gzbuffer(*gzipped,GZBUFFER_SIZE);
#endif
	if (interleavedp == true) {
	  *nextchar = '\0';
	} else {
	  *nextchar = Shortread_input_init_gzip(*gzipped);
	}
      }
#endif

    } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
      if ((*bzipped = Bzip2_new((*files)[0])) == NULL) {
	fprintf(stderr,"Cannot open bzipped file %s\n",(*files)[0]);
	exit(9);
      } else {
	if (interleavedp == true) {
	  *nextchar = '\0';
	} else {
	  *nextchar = Shortread_input_init_bzip2(*bzipped);
	}
      }
#endif

    } else {
      if ((*input = Fopen_read_text(read_files_command,(*files)[0])) == NULL) {
	fprintf(stderr,"Unable to read input\n");
	exit(9);
      } else {
	debugf(fprintf(stderr,"Master opening file %s using fopen for input\n",(*files)[0]));
	if (interleavedp == true) {
	  *nextchar = '\0';
	} else {
	  *nextchar = Shortread_input_init(&(*nchars1),*input);
	}
      }
    }

    (*files)++;
    (*nfiles)--;
  }

  if (interleavedp == true) {
    /* Read one line at a time */
    /* I believe interleaved is paired-end by definition */
    return /*fastq_format_p*/false;
    
  } else if (*nextchar == EOF) {
    fprintf(stderr,"Warning: input is empty\n");
    exit(0);

  } else if (*nextchar == '@') {
    /* Looks like a FASTQ file */
    if (*nfiles == 0 || force_single_end_p == true) {
#ifdef HAVE_ZLIB
      *gzipped2 = (gzFile) NULL;
#endif
#ifdef HAVE_BZLIB
      *bzipped2 = (Bzip2_T) NULL;
#endif
      *input2 = (FILE *) NULL;
      *paired_end_p = false;

    } else {
      if (gunzip_p == true) {
#ifdef HAVE_ZLIB
	if ((*gzipped2 = gzopen((*files)[0],"rb")) == NULL) {
	  fprintf(stderr,"Cannot open gzipped file %s\n",(*files)[0]);
	  exit(9);
	} else {
#ifdef HAVE_ZLIB_GZBUFFER
	  gzbuffer(*gzipped2,GZBUFFER_SIZE);
#endif
	  /* nextchar2 = */ Shortread_input_init_gzip(*gzipped2);
	}
#endif

      } else if (bunzip2_p == true) {
#ifdef HAVE_BZLIB
	if ((*bzipped2 = Bzip2_new((*files)[0])) == NULL) {
	  fprintf(stderr,"Cannot open bzip2 file %s\n",(*files)[0]);
	  exit(9);
	} else {
	  /* nextchar2 = */ Shortread_input_init_bzip2(*bzipped2);
	}
#endif

      } else {
	if ((*input2 = Fopen_read_text(read_files_command,(*files)[0])) == NULL) {
	  fprintf(stderr,"Unable to read input\n");
	  exit(9);
	} else {
	  debugf(fprintf(stderr,"Master opening file %s using fopen for input2\n",(*files)[0]));
	  /* nextchar2 = */ Shortread_input_init(&(*nchars2),*input2);
	}
      }
      /* Keep paired_end_p as true */
      (*files)++;
      (*nfiles)--;
    }
    fastq_format_p = true;

  } else if (*nextchar == '>') {
    /* Looks like a FASTA file */
    *paired_end_p = false;

  } else {
    fprintf(stderr,"First char is %c.  Expecting either '>' for FASTA or '@' for FASTQ format.  If file is gzipped, specify --gunzip\n",*nextchar);
    exit(9);
  }

  return fastq_format_p;
}


static Univ_IIT_T
transcript_iit_setup (Trcoord_T *transcriptomelength, int *ntranscripts,
		      char *transcriptomesubdir, char *fileroot) {
  Univ_IIT_T transcript_iit = NULL;
  char *iitfile = NULL;

  /* Prepare transcript data */

  iitfile = (char *) CALLOC(strlen(transcriptomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",transcriptomesubdir,fileroot);
  if ((transcript_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",iitfile);
    exit(9);
  } else {
    FREE(iitfile);
    *ntranscripts = Univ_IIT_total_nintervals(transcript_iit);
  }

  *transcriptomelength = (Trcoord_T) Univ_IIT_genomelength(transcript_iit,/*with_circular_alias_p*/false);
  return transcript_iit;
}



static Univ_IIT_T
chromosome_iit_setup (Univcoord_T *genomelength, int *nchromosomes, int *circular_typeint,
		      bool *any_circular_p, bool **circularp, char *genomesubdir, char *fileroot) {
  Univ_IIT_T chromosome_iit = NULL, altscaffold_iit = NULL;
  char *iitfile = NULL;


  /* Prepare genomic data */

  iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			    strlen(fileroot)+strlen(".chromosome.iit")+1,sizeof(char));
  sprintf(iitfile,"%s/%s.chromosome.iit",genomesubdir,fileroot);
  if ((chromosome_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
    fprintf(stderr,"IIT file %s is not valid\n",iitfile);
    exit(9);
#ifdef LARGE_GENOMES
  } else if (Univ_IIT_coord_values_8p(chromosome_iit) == false) {
    fprintf(stderr,"This program gsnapl is designed for large genomes.\n");
    fprintf(stderr,"For small genomes of less than 2^32 (4 billion) bp, please run gsnap instead.\n");
    exit(9);
#endif
  } else {
    FREE(iitfile);
    *nchromosomes = Univ_IIT_total_nintervals(chromosome_iit);
    *circular_typeint = Univ_IIT_typeint(chromosome_iit,"circular");
    *circularp = Univ_IIT_circularp(&(*any_circular_p),chromosome_iit);

    iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
			      strlen(fileroot)+strlen(".altscaffold.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.altscaffold.iit",genomesubdir,fileroot);
    if ((altscaffold_iit = Univ_IIT_read(iitfile,/*readonlyp*/true,/*add_iit_p*/false)) == NULL) {
      /* fprintf(stderr,"No altscaffold file found\n"); */
      altlocp = (bool *) CALLOC((*nchromosomes)+1,sizeof(bool));
      alias_starts = (Univcoord_T *) CALLOC((*nchromosomes)+1,sizeof(Univcoord_T));
      alias_ends = (Univcoord_T *) CALLOC((*nchromosomes)+1,sizeof(Univcoord_T));

    } else {
      fprintf(stderr,"Found altscaffold file found\n");
      altlocp = Univ_IIT_altlocp(&alias_starts,&alias_ends,chromosome_iit,altscaffold_iit);
      Univ_IIT_free(&altscaffold_iit);
    }

    FREE(iitfile);
  }

  *genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
  return chromosome_iit;
}


static void
worker_setup (char *transcriptomesubdir, char *transcriptome_fileroot,
	      char *genomesubdir, char *genome_fileroot, Univ_IIT_T chromosome_iit) {
  char *snpsdir = NULL, *modedir = NULL, *mapdir = NULL, *iitfile = NULL;
  /* Univcoord_T genome_totallength; */
  char *idx_filesuffix1, *idx_filesuffix2;
  FILE *dump_fp;


  if (snps_root == NULL && maskedp == false) {
    if (transcriptome_fileroot != NULL) {
      transcriptomebits = Genomebits_new(transcriptomesubdir,transcriptome_fileroot,/*snps_root*/NULL,
					 genome_access,sharedp,/*revcompp*/false);
    }

    genome = genomealt = Genome_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
				    chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    genomebits = Genomebits_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
				genome_access,sharedp,/*revcompp*/false);

  } else if (maskedp == true) {
    genome = Genome_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
			chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    genomebits = Genomebits_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
				genome_access,sharedp,/*revcompp*/false);
    genomealt = Genome_new(genomesubdir,genome_fileroot,/*alt_root*/masked_suffix,
			   chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    genomebits_alt = Genomebits_new(genomesubdir,genome_fileroot,/*alt_root*/masked_suffix,
				    genome_access,sharedp,/*revcompp*/false);
  } else {
    if (user_snpsdir == NULL) {
      snpsdir = genomesubdir;
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,genome_fileroot);
    } else {
      snpsdir = user_snpsdir;
      mapdir = user_snpsdir;
    }

    /* SNPs */
    if (transcriptome_fileroot != NULL) {
      transcriptomebits = Genomebits_new(transcriptomesubdir,transcriptome_fileroot,/*snps_root*/NULL,
					 genome_access,sharedp,/*revcompp*/false);
      transcriptomebits_alt = Genomebits_new(snpsdir,transcriptome_fileroot,snps_root,
					     genome_access,sharedp,/*revcompp*/false);
    }

    genome = Genome_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
			chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    genomebits = Genomebits_new(genomesubdir,genome_fileroot,/*alt_root*/NULL,
				genome_access,sharedp,/*revcompp*/false);
    genomealt = Genome_new(snpsdir,genome_fileroot,snps_root,
			   chromosome_iit,genome_access,sharedp,/*revcompp*/false);
    genomebits_alt = Genomebits_new(snpsdir,genome_fileroot,snps_root,
				    genome_access,sharedp,/*revcompp*/false);
  }

  /* Must be done before Knownsplicing_retrieve_via_splicesites */
  Genome_setup(genome,genomealt,circular_typeint);

  if (user_modedir == NULL) {
    modedir = genomesubdir;
  } else {
    modedir = user_modedir;
  }

  if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
    idx_filesuffix1 = "metct";
    idx_filesuffix2 = "metga";
  } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED) {
    idx_filesuffix1 = "a2iag";
    idx_filesuffix2 = "a2itc";
  } else if (mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
    idx_filesuffix1 = "a2itc";
    idx_filesuffix2 = "a2iag";
  } else {
    idx_filesuffix1 = IDX_FILESUFFIX; /* "ref" */
    idx_filesuffix2 = (char *) NULL;
  }

  if ((indexdb = Indexdb_new_genome(&index1part,&index1interval,
				    /*genomesubdir*/modedir,snpsdir,
				    genome_fileroot,idx_filesuffix1,snps_root,
				    required_index1part,required_index1interval,
				    offsetsstrm_access,positions_access,sharedp,
				    multiple_sequences_p,preload_shared_memory_p,unload_shared_memory_p)) == NULL) {

    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      fprintf(stderr,"Cannot find %s index file.  Need to run cmetindex first\n",idx_filesuffix1);
    } else if (mode == ATOI_STRANDED || mode == ATOI_NONSTRANDED ||
	       mode == TTOC_STRANDED || mode == TTOC_NONSTRANDED) {
      fprintf(stderr,"Cannot find %s index file.  Need to run atoiindex first\n",idx_filesuffix1);
    } else {
      fprintf(stderr,"Cannot find %s index file\n",idx_filesuffix1);
    }
    exit(9);
  }

  if (idx_filesuffix2 == NULL) {
    indexdb_nonstd = indexdb;
  } else if ((indexdb_nonstd = Indexdb_new_genome(&index1part,&index1interval,
						  /*genomesubdir*/modedir,snpsdir,
						  genome_fileroot,idx_filesuffix2,snps_root,
						  required_index1part,required_index1interval,
						  offsetsstrm_access,positions_access,sharedp,
						  multiple_sequences_p,preload_shared_memory_p,unload_shared_memory_p)) == NULL) {
    if (mode == CMET_STRANDED || mode == CMET_NONSTRANDED) {
      fprintf(stderr,"Cannot find %s index file.  Need to run cmetindex first\n",idx_filesuffix2);
    } else {
      fprintf(stderr,"Cannot find %s index file.  Need to run atoiindex first\n",idx_filesuffix2);
    }
    exit(9);
  }

  if (snps_root != NULL) {
    iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
			      strlen(snps_root)+1,sizeof(char));
    sprintf(iitfile,"%s/%s",mapdir,snps_root);
    fprintf(stderr,"Reading SNPs file %s/%s...",mapdir,snps_root);
    if ((snps_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			     /*divstring*/NULL,/*add_iit_p*/true)) == NULL) {
      fprintf(stderr,"SNPs file %s.iit not found in %s.\n",snps_root,mapdir);
      if (user_snpsdir == NULL) {
	fprintf(stderr,"Available files:\n");
	Datadir_list_directory(stderr,genomesubdir);
	fprintf(stderr,"Either install file %s.iit or specify a directory for the IIT file\n",snps_root);
	fprintf(stderr,"using the -M flag.\n");
	exit(9);
      }
    }

    print_nsnpdiffs_p = true;
    snps_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,snps_iit);

    fprintf(stderr,"done\n");
    FREE(iitfile);
    if (user_snpsdir == NULL) {
      FREE(mapdir);
    }
  }


  indexdb_size_threshold = (int) (10*Indexdb_mean_size(indexdb,mode,index1part));
  debug(printf("Size threshold is %d\n",indexdb_size_threshold));
  if (indexdb_size_threshold < MIN_INDEXDB_SIZE_THRESHOLD) {
    indexdb_size_threshold = MIN_INDEXDB_SIZE_THRESHOLD;
  }

#if 0
  if (genes_file != NULL) {
    if ((genes_iit = IIT_read(genes_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
			      /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
      fprintf(stderr,"Reading genes file %s locally...",genes_file);
    } else {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,genome_fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(genes_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,genes_file);
      if ((genes_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading genes file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Genes file %s.iit not found locally or in %s.  Available files:\n",genes_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",genes_file);
	exit(9);
      }
    }
    genes_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,genes_iit);
    genes_chrnum_crosstable = Univ_IIT_chrnum_crosstable(chromosome_iit,genes_iit);
  }
#endif

  if (splicing_file != NULL) {
    if (user_splicingdir == NULL) {
      if ((splicing_iit = IIT_read(splicing_file,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s locally...",splicing_file);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_splicingdir)+strlen("/")+strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_splicingdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (splicing_iit == NULL) {
      mapdir = Datadir_find_mapdir(/*user_mapdir*/NULL,genomesubdir,genome_fileroot);
      iitfile = (char *) CALLOC(strlen(mapdir)+strlen("/")+
				strlen(splicing_file)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",mapdir,splicing_file);
      if ((splicing_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				   /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading splicing file %s...",iitfile);
	FREE(iitfile);
	FREE(mapdir);
      } else {
	fprintf(stderr,"Splicing file %s.iit not found locally or in %s.  Available files:\n",splicing_file,mapdir);
	Datadir_list_directory(stderr,mapdir);
	fprintf(stderr,"Either install file %s or specify a full directory path\n",splicing_file);
	exit(9);
      }
    }

    splicing_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,splicing_iit);
    if ((donor_typeint = IIT_typeint(splicing_iit,"donor")) >= 0 && 
	(acceptor_typeint = IIT_typeint(splicing_iit,"acceptor")) >= 0) {
      fprintf(stderr,"found donor and acceptor tags, so treating as splicesites file.  Note: intron-level splicing information gives better results\n");
      intron_level_p = false;
      knownsplicing = Knownsplicing_from_splicing_iit(splicing_iit,splicing_divint_crosstable,
						      donor_typeint,acceptor_typeint,chromosome_iit,
						      /*intron_level_p*/false);
      if (Knownsplicing_nintervals(knownsplicing) == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		genome_dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);
      }

    } else {
      fprintf(stderr,"did not find donor and acceptor tags, so treating as an intron file\n");
      intron_level_p = true;
      knownsplicing = Knownsplicing_from_splicing_iit(splicing_iit,splicing_divint_crosstable,
						      donor_typeint,acceptor_typeint,chromosome_iit,
						      /*intron_level_p*/true);
      if (Knownsplicing_nintervals(knownsplicing) == 0) {
	fprintf(stderr,"\nWarning: No splicesites observed for genome %s.  Are you sure this splicesite file was built for this genome?  Please compare chromosomes below:\n",
		genome_dbroot);
	fprintf(stderr,"Chromosomes in the genome: ");
	Univ_IIT_dump_labels(stderr,chromosome_iit);
	fprintf(stderr,"Chromosomes in the splicesites IIT file: ");
	IIT_dump_divstrings(stderr,splicing_iit);
	exit(9);
      }
    }

    /* For benchmarking purposes.  Can spend time/memory to load
       splicesites, but then not use them. */
    if (unloadp == true) {
      fprintf(stderr,"unloading...");

      FREE(splicing_divint_crosstable);
      IIT_free(&splicing_iit);
      splicing_iit = NULL;
      knownsplicingp = false;
      splicing_file = (char *) NULL;
    }

    fprintf(stderr,"done\n");
  }

#if 0
  if (splices_file == NULL) {
    /* Skip */
  } else if ((dump_fp = fopen(splices_file,"r")) == NULL) {
    fprintf(stderr,"Cannot open splices file %s\n",splices_file);
    exit(9);
  } else {
    knownsplicing = Knownsplicing_new_from_dump(dump_fp,genomelength);
    fclose(dump_fp);
  }
#endif

  if (indels_file == NULL) {
    /* Skip */
  } else if ((dump_fp = fopen(indels_file,"r")) == NULL) {
    fprintf(stderr,"Cannot open indels file %s\n",indels_file);
    exit(9);
  } else {
    knownindels = Knownindels_new_from_dump(dump_fp,genomelength);
    fclose(dump_fp);
  }


  if (novelsplicingp == true || knownsplicingp == true) {
    /* RNA-seq */
    splicingp = true;
  } else {
    /* DNA-seq */
    splicingp = false;
  }


  if (tally_root != NULL) {
    if (user_tallydir == NULL) {
      if ((tally_iit = IIT_read(tally_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s.iit locally...",tally_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_tallydir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_tallydir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (tally_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(tally_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,tally_root);
      if ((tally_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				/*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading tally file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Tally file %s.iit not found locally",tally_root);
	if (user_tallydir != NULL) {
	  fprintf(stderr," or in %s",user_tallydir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    tally_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,tally_iit);
    fprintf(stderr,"done\n");
  }


  if (runlength_root != NULL) {
    if (user_runlengthdir == NULL) {
      if ((runlength_iit = IIT_read(runlength_root,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s.iit locally...",runlength_root);
      }
    } else {
      iitfile = (char *) CALLOC(strlen(user_runlengthdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",user_runlengthdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      }
    }

    if (runlength_iit == NULL) {
      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+strlen(runlength_root)+1,sizeof(char));
      sprintf(iitfile,"%s/%s",genomesubdir,runlength_root);
      if ((runlength_iit = IIT_read(iitfile,/*name*/NULL,/*readonlyp*/true,/*divread*/READ_ALL,
				    /*divstring*/NULL,/*add_iit_p*/true)) != NULL) {
	fprintf(stderr,"Reading runlength file %s...",iitfile);
	FREE(iitfile);
      } else {
	fprintf(stderr,"Runlength file %s.iit not found locally",runlength_root);
	if (user_runlengthdir != NULL) {
	  fprintf(stderr," or in %s",user_runlengthdir);
	}
	fprintf(stderr," or in %s\n",genomesubdir);
	exit(9);
      }
    }

    runlength_divint_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,runlength_iit);
    fprintf(stderr,"done\n");
  }


  /* Approx procedures used by transcriptome-search and kmer-search */
#ifdef LARGE_GENOMES
  Intersect_approx_uint8_setup();
  Intersect_approx_indices_uint8_setup();
#else
  Intersect_approx_indices_uint4_setup();
#endif
  Intersect_approx_uint4_setup();
  Intersect_wdups_indices_setup();

  /* Lower/higher and indices procedures used by segment-search */
  Intersect_lower_setup();
  Intersect_higher_setup();
#ifdef LARGE_GENOMES
  Intersect_indices2_large_setup();
  Intersect_indices_uint8_setup();
#else
  Intersect_indices2_small_setup();
#endif
  Intersect_indices_small_setup();

  Intersect_concordance_setup();

  Altsplice_setup(max_insertlength,genomebits,genomebits_alt);
  Path_setup(circularp,altlocp);
  Pathpair_setup(output_type,genomelength,shortsplicedist);
  Repair_setup(genomebits,genomebits_alt);
  Path_solve_setup(circularp,transcriptome,chromosome_ef64,
		   genomebits,genomebits_alt,genomelength,
		   index1part,min_intronlength,
		   max_insertionlen,max_deletionlen,
		   novelsplicingp,knownsplicingp);
  Path_fusion_setup(circularp,chromosome_ef64,genomelength,shortsplicedist,transcriptome,
		    genomebits);
  Path_trim_setup(genomebits,genomebits_alt);
  Path_eval_setup(genomebits,genomebits_alt,transcriptome,
		  circularp,chrsubsetp,altlocp,index1part,index1interval,
		  output_type,md_report_snps_p,want_random_p,allow_soft_clips_p);
  Pathpair_eval_setup(expected_pairlength,pairlength_deviation,transcriptome,
		      circularp,chrsubsetp,altlocp,
		      output_type,splicingp,resolve_inner_p,want_random_p,
		      allow_soft_clips_p);
  Path_print_alignment_setup(chromosome_iit,method_print_p,print_univdiagonal_p,print_nsnpdiffs_p);
  Path_print_m8_setup(chromosome_iit);
  Path_print_sam_setup(add_paired_nomappers_p,paired_flag_means_concordant_p,sam_insert_0M_p,
		       quiet_if_excessive_p,maxpaths_report,
		       failedinput_root,fastq_format_p,extend_soft_clips_p,method_print_p,
		       only_concordant_p,omit_concordant_uniq_p,omit_concordant_mult_p,
		       only_tr_consistent_p,circularp,
		       clip_overlap_p,merge_overlap_p,merge_samechr_p,
		       sam_multiple_primaries_p,sam_sparse_secondaries_p,
		       chromosome_iit,transcript_iit,snps_iit,maskedp);

  Trpath_solve_setup(transcriptome,transcriptomebits,transcript_ef64,
		     max_insertionlen,max_deletionlen,index1part_tr);
  Trpath_eval_setup(transcriptomebits);
  Trpath_convert_setup(transcriptome,chromosome_ef64);
  Transcriptome_search_setup(index1part_tr,index1interval_tr,tr_indexdb,
			     transcriptome,transcript_ef64,
			     genomebits,genomebits_alt,transcriptomebits);
  Tr_extension_search_setup(transcriptome,transcriptomelength,transcript_ef64,
			    transcriptomebits,tr_indexdb,index1part_tr,
			    max_insertionlen,max_deletionlen,maxpaths_search);

  Kmer_search_setup(index1part,index1interval,indexdb,chromosome_ef64,
		    genomebits,genomebits_alt,genomelength,
		    max_insertionlen,max_deletionlen,shortsplicedist,splicingp);
  Extension_search_setup(mode,genomelength,
			 circular_typeint,circularp,chromosome_ef64,
			 genomebits,genomebits_alt,indexdb,indexdb_nonstd,
			 index1part,index1interval,maxpaths_search,
			 max_insertionlen,max_deletionlen,shortsplicedist);
  /* Merge_search_setup(index1part,chromosome_ef64); */
  /* Segment_search_setup(index1part,index1interval,chromosome_iit,nchromosomes,
                          circular_typeint,chromosome_ef64,mode); */


  if (genomebits == NULL) {
    fprintf(stderr,"This version of GSNAP requires the genomefwdh file for the genome\n");
    exit(9);
  } else {
    Genomebits_consec_setup(query_unk_mismatch_p,genome_unk_mismatch_p,mode);
    Genomebits_count_setup(genomebits,/*snp_blocks*/genomebits_alt ? genomebits_alt : NULL,
			   query_unk_mismatch_p,genome_unk_mismatch_p,
			   mode,md_report_snps_p,maskedp);
    Genomebits_mismatches_setup(query_unk_mismatch_p,genome_unk_mismatch_p,mode,maskedp);
    Genomebits_trim_setup(mode,maskedp);
    Genomebits_kmer_setup(genome_unk_mismatch_p,mode);
    Genomebits_indel_setup(max_insertionlen,max_deletionlen,
			   query_unk_mismatch_p,genome_unk_mismatch_p,mode,maskedp);
  }

  Genome_sites_setup(genome,genomealt);
  Maxent_hr_setup(genome,genomealt);

  Simplepair_setup(splicingp,transcript_iit,sam_insert_0M_p,
		   /*md_lowercase_variant_p*/false,/*snps_p*/snps_iit ? true : false,
		   sam_cigar_extended_p,cigar_action);

  /* Oligoindex_hr_setup(Genome_blocks(genome),mode); */
  /* Oligoindex_localdb_setup(chromosome_iit,circular_typeint,localdb,local1part); */

  min_querylength = index1part + index1interval - 1;
  Indexdb_setup(index1part);
  /* Localdb_setup(genome,genomelength,index1part,index1interval); */

  Splice_setup(genomebits,genomebits_alt,max_insertionlen,max_deletionlen,
	       novelsplicingp,splice_evalp);
  Spliceends_setup(circularp,genomebits,genomebits_alt,genomelength,
		   indexdb,index1part,index1interval,
		   max_insertionlen,max_deletionlen,shortsplicedist,
		   allow_soft_clips_p,novelsplicingp,knownsplicingp);
  Indel_setup(genomebits,genomebits_alt,
	      min_indel_end_matches,maskedp);
  Transcript_setup(pairmax_transcriptome,output_type,transcriptome,transcript_iit);
  Transcript_remap_setup(circularp,transcriptome,transcript_map_iit,transcript_chrnum_crosstable);
  Single_cell_setup(wellpos);

  repetitive_ef64 = Repetitive_setup(index1part);
  Pathstore_setup(chromosome_ef64);
  Stage1hr_setup(indexdb,indexdb_nonstd,tr_indexdb,repetitive_ef64,
		 index1part,index1interval,index1part_tr,index1interval_tr,
		 max_deletionlen,shortsplicedist,transcriptome);
  Stage1hr_single_setup(mode,index1part,index1interval,index1part_tr,
			transcriptome,genome_align_p,transcriptome_align_p,
			user_nmismatches_filter_float,user_mincoverage_filter_float,splicingp);
  Stage1hr_paired_setup(mode,index1part,index1interval,index1part_tr,
			transcriptome,genome_align_p,transcriptome_align_p,
			genomebits,chromosome_ef64,
			user_nmismatches_filter_float,user_mincoverage_filter_float,
			max_deletionlen,max_insertlength,shortsplicedist,
			splicingp,maxpaths_search,maxpaths_report,
			circularp,pairmax_linear,pairmax_circular);

  Concordance_setup(subopt_levels,pairmax_transcriptome,pairmax_linear,pairmax_circular,
		    circularp,merge_samechr_p,two_pass_p);

#if 0
  SAM_setup(add_paired_nomappers_p,paired_flag_means_concordant_p,
	    only_concordant_p,omit_concordant_uniq_p,omit_concordant_mult_p,quiet_if_excessive_p,
	    maxpaths_report,failedinput_root,fastq_format_p,extend_soft_clips_p,method_print_p,
	    clip_overlap_p,merge_overlap_p,merge_samechr_p,sam_multiple_primaries_p,sam_sparse_secondaries_p,
	    force_xs_direction_p,snps_iit,maskedp,find_dna_chimeras_p,
	    splicingp,donor_typeint,acceptor_typeint,chromosome_iit,transcript_iit);
  Cigar_setup(genome,sam_cigar_extended_p,extend_soft_clips_p,merge_samechr_p,md_lowercase_variant_p,
	      sam_hardclip_use_S_p,sam_insert_0M_p);
#endif

  Output_setup(chromosome_iit,nofailsp,failsonlyp,quiet_if_excessive_p,maxpaths_report,
	       failedinput_root,quality_shift,
	       output_type,invert_first_p,invert_second_p,single_cell_p,
	       only_concordant_p,only_tr_consistent_p,sam_read_group_id);

  return;
}


static void
worker_localdb_setup (char *genomesubdir, char *genome_fileroot) {

  /* Note: localdb based on suffix arrays, so it is not SNP-tolerant */
  if (use_localdb_p == false) { 
    fprintf(stderr,"Not loading localdb files because --use-local=0\n");
    localdb = (Localdb_T) NULL; 

  } else if (splicingp == false) {
    fprintf(stderr,"Not loading localdb files because reads are DNA-seq (without -N, -s, or -c flags)\n");
    localdb = (Localdb_T) NULL; 

  } else {
    localdb = Localdb_new(genomesubdir,genome_fileroot,localdb_access,sharedp,
			  multiple_sequences_p,preload_shared_memory_p,unload_shared_memory_p);
#if 0
    if (localdb == NULL) {
      fprintf(stderr,"Warning: Cannot find localdb files.  GSNAP will run on this older index, but there is a newer index format available\n");
    }
#endif
  }
  
  Localdb_setup(genome,genomelength,index1part,index1interval);
  Path_solve_localdb_setup(localdb);
  Spliceends_localdb_setup(localdb);
  Stage1hr_paired_localdb_setup(localdb);

  return;
}


static void
worker_cleanup () {

  if (whitelist_file != NULL) {
    Single_cell_cleanup();
  }

  /* Segment_search_cleanup(); */

  EF64_free(&repetitive_ef64);

  if (localdb != NULL) {
    Localdb_free(&localdb);
  }

  if (knownindels != NULL) {
    Knownindels_free(&knownindels);
  }

  if (knownsplicing != NULL) {
    Knownsplicing_free(&knownsplicing);
  }

  if (indexdb_nonstd != indexdb) {
    Indexdb_free(&indexdb_nonstd);
  }
  if (indexdb != NULL) {
    Indexdb_free(&indexdb);
  }

  if (transcriptomebits_alt != NULL) {
    Genomebits_free(&transcriptomebits_alt);
  }
  if (transcriptomebits != NULL) {
    Genomebits_free(&transcriptomebits);
  }


  if (genomealt != NULL && genomealt != genome) {
    Genome_free(&genomealt);
    Genomebits_free(&genomebits_alt);
  }
  if (genomebits != NULL) {
    Genomebits_free(&genomebits);
  }
  if (genome != NULL) {
    Genome_free(&genome);
  }

  if (runlength_iit != NULL) {
    FREE(runlength_divint_crosstable);
    IIT_free(&runlength_iit);
  }

  if (tally_iit != NULL) {
    FREE(tally_divint_crosstable);
    IIT_free(&tally_iit);
  }

  if (splicing_iit != NULL) {
    FREE(splicing_divint_crosstable);
    IIT_free(&splicing_iit);
  }

#if 0
  if (genes_iit != NULL) {
    FREE(genes_chrnum_crosstable);
    FREE(genes_divint_crosstable);
    IIT_free(&genes_iit);
  }
#endif

  if (snps_iit != NULL) {
    FREE(snps_divint_crosstable);
    IIT_free(&snps_iit);
  }

  if (altlocp != NULL) {
    FREE(alias_ends);
    FREE(alias_starts);
    FREE(altlocp);
  }

  if (chrsubsetp != NULL) {
    FREE(chrsubsetp);
  }

  if (circularp != NULL) {
    FREE(circularp);
  }

  if (tr_indexdb != NULL) {
    Indexdb_free(&tr_indexdb);
  }

  if (transcriptome != NULL) {
    Transcriptome_free(&transcriptome);
  }

  if (transcript_map_iit != NULL) {
    IIT_free(&transcript_map_iit);
    FREE(transcript_chrnum_crosstable);
  }

  if (transcript_ef64 != NULL) {
    EF64_free(&transcript_ef64);
  }

  if (transcript_iit != NULL) {
    Univ_IIT_free(&transcript_iit);
  }

  if (chromosome_ef64 != NULL) {
    EF64_free(&chromosome_ef64);
  }

  if (chromosome_iit != NULL) {
    Univ_IIT_free(&chromosome_iit);
  }

  Access_controlled_cleanup();

  return;
}


int
main (int argc, char *argv[]) {
  int nchars1 = 0, nchars2 = 0;
  int cmdline_status;

  char *transcriptomesubdir = NULL, *genomesubdir, *transcriptome_fileroot = NULL, *genome_fileroot, *dbversion;
  char *iitfile;
  char **files;
  int nfiles;

  char *line, *chr, *p, c;
  int length;
  Chrnum_T chrnum;
  FILE *fp, *input, *input2;

#ifdef HAVE_ZLIB
  gzFile gzipped, gzipped2;
#endif

#ifdef HAVE_BZLIB
  Bzip2_T bzipped, bzipped2;
#endif

  long int worker_id;

  int nread;
  int nextchar = '\0';
  double runtime;

#ifdef HAVE_PTHREAD
  int ret;
  pthread_attr_t thread_attr_join;
#endif

#ifdef HAVE_SIGACTION
  struct sigaction signal_action;
#endif

  extern int optind;

#ifdef MEMUSAGE
  Mem_usage_init();
  Mem_usage_set_threadname("main");
#endif


  cmdline_status = parse_command_line(argc,argv,optind);
  argc -= optind;
  argv += optind;

  if (cmdline_status == 0) {
    /* okay to continue */
  } else if (cmdline_status == 1) {
    exit(0);
  } else {
    exit(cmdline_status);
  }

  check_compiler_assumptions();

  if (exception_raise_p == false) {
    fprintf(stderr,"Allowing signals and exceptions to pass through.  If using shared memory, need to remove segments manually.\n");
    Except_inactivate();
  } else {
#ifdef HAVE_SIGACTION
    signal_action.sa_handler = signal_handler;
    signal_action.sa_flags = 0;
    sigfillset(&signal_action.sa_mask); /* After first signal, block all other signals */

    /* Note: SIGKILL and SIGSTOP cannot be caught */

    sigaction(SIGABRT,&signal_action,NULL); /* abnormal termination (abort) */
    sigaction(SIGBUS,&signal_action,NULL);  /* bus error */
    sigaction(SIGFPE,&signal_action,NULL);  /* arithmetic exception */
    sigaction(SIGHUP,&signal_action,NULL);  /* hangup */
    sigaction(SIGILL,&signal_action,NULL);  /* illegal hardware instruction */
    sigaction(SIGINT,&signal_action,NULL);  /* terminal interruption (control-C) */
    sigaction(SIGPIPE,&signal_action,NULL);  /* write to pipe with no readers */
    sigaction(SIGQUIT,&signal_action,NULL);  /* terminal quit (control-backslash) */
    sigaction(SIGSEGV,&signal_action,NULL);  /* invalid memory reference */
    sigaction(SIGSYS,&signal_action,NULL);  /* invalid system call */
    sigaction(SIGTERM,&signal_action,NULL);  /* Unix kill command */
    sigaction(SIGTRAP,&signal_action,NULL);  /* hardware fault */
    sigaction(SIGXCPU,&signal_action,NULL);  /* CPU limit exceeded */
    sigaction(SIGXFSZ,&signal_action,NULL);  /* file size limit exceeded */
#endif
  }

  if (whitelist_file != NULL) {
    Single_cell_compute_priors(whitelist_file,
			       read_files_command,gunzip_p,bunzip2_p,
			       /*files*/argv,/*nfiles*/argc);
  }

  fastq_format_p = open_input_streams_parser(&nextchar,&nchars1,&nchars2,&paired_end_p,
					     &files,&nfiles,&input,&input2,
#ifdef HAVE_ZLIB
					     &gzipped,&gzipped2,
#endif
#ifdef HAVE_BZLIB
					     &bzipped,&bzipped2,
#endif
					     read_files_command,gunzip_p,bunzip2_p,interleavedp,argc,argv);


  /* If we are processing oligos (e.g., for repetitiveness), this
     needs to occur before the initial call to Inbuffer_fill_init */
  Oligo_setup(mode);

  /* Needs to be after open_input_streams_parser to get correct fastq_format_p */
  Shortread_setup(acc_fieldi_start,acc_fieldi_end,force_single_end_p,
		  filter_chastity_p,keep_chastity_p,
		  allow_paired_end_mismatch_p,fastq_format_p,barcode_length,endtrim_length,
		  invert_first_p,invert_second_p,chop_poly_at_first_p,chop_poly_at_second_p);

  Inbuffer_setup(single_cell_p,filter_if_both_p);

  inbuffer = Inbuffer_new(nextchar,input,input2,
#ifdef HAVE_ZLIB
			  gzipped,gzipped2,
#endif
#ifdef HAVE_BZLIB
			  bzipped,bzipped2,
#endif
			  interleavedp,read_files_command,files,nfiles,input_buffer_size,
			  skip_prob,/*seed*/SEED,part_modulus,part_interval);

  if ((nread = Inbuffer_fill_init(inbuffer)) > 1) {
    multiple_sequences_p = true;
  } else {
    multiple_sequences_p = false;
  }

  if (preload_shared_memory_p == true || unload_shared_memory_p == true) {
    /* No need for pre-load */

  } else if (multiple_sequences_p == true) {
    /* Expected mode */

  } else {
    /* Faster processing of a single sequence by using mmap instead of allocate */
    /* fprintf(stderr,"Note: only 1 sequence detected.  Ignoring batch (-B) command\n"); */
    expand_offsets_p = false;
#ifdef HAVE_MMAP
    offsetsstrm_access = USE_MMAP_ONLY;
    positions_access = USE_MMAP_ONLY;
    genome_access = USE_MMAP_ONLY;
    locoffsetsstrm_access = USE_MMAP_ONLY;
    locpositions_access = USE_MMAP_ONLY;
    localdb_access = USE_MMAP_ONLY;
#else
    /* No choice, since mmap is not available */
    offsetsstrm_access = USE_ALLOCATE;
    positions_access = USE_ALLOCATE;
    genome_access = USE_ALLOCATE;
    locoffsetsstrm_access = USE_ALLOCATE;
    locpositions_access = USE_ALLOCATE;
    localdb_access = USE_ALLOCATE;
#endif
  }

  genomesubdir = Datadir_find_genomesubdir(&genome_fileroot,&dbversion,user_genomedir,genome_dbroot);
  FREE(dbversion);

  chromosome_iit = chromosome_iit_setup(&genomelength,&nchromosomes,&circular_typeint,&any_circular_p,&circularp,
					genomesubdir,genome_fileroot);
  chromosome_ef64 = EF64_new_from_chromosome_iit(chromosome_iit);

  if (chrsubset_file == NULL) {
    /* Skip */
  } else if ((fp = fopen(chrsubset_file,"r")) == NULL) {
    fprintf(stderr,"Unable to open --chrsubset file %s\n",chrsubset_file);
    exit(9);
  } else {
    chrsubsetp = (bool *) CALLOC(nchromosomes+1,sizeof(bool)); /* Chrnum_T is 1-based */
    while ((line = Getline(fp)) != NULL) {
      p = line;
      while ((c = *p++) != '\0' && !isspace(c)) ;

      length = (p - 1) - line;
      chr = (char *) MALLOC((length+1)*sizeof(char));
      strncpy(chr,line,length);
      chr[length] = '\0';
	
      if ((chrnum = Univ_IIT_find_one(chromosome_iit,chr)) <= 0) {
	fprintf(stderr,"Unable to find chromosome %s in --chrsubset file\n",line);
	exit(9);
      } else {
	chrsubsetp[chrnum] = true;
      }
      FREE(chr);
      FREE(line);
    }
    fclose(fp);
  }

  if (transcriptome_dbroot == NULL) {
    transcriptome_align_p = false;

  } else {
    if (user_transcriptomedir == NULL && user_genomedir != NULL) {
      user_transcriptomedir = user_genomedir;
    }
    transcriptomesubdir = Datadir_find_genomesubdir(&transcriptome_fileroot,&dbversion,
						    user_transcriptomedir,transcriptome_dbroot);
    FREE(dbversion);

    /* transcript_iit relates transcripts to transcriptome coordinates and is a Univ_IIT_T */
    transcript_iit = transcript_iit_setup(&transcriptomelength,&ntranscripts,transcriptomesubdir,
					  transcriptome_fileroot);
    transcript_ef64 = EF64_new_from_chromosome_iit(transcript_iit);


    /* transcript_map_iit relates transcripts to genomic coordinates and is an IIT_T */
    /* New format */
    iitfile = (char *) CALLOC(strlen(transcriptomesubdir)+strlen("/")+
			      strlen(transcriptome_fileroot)+strlen(".genomes/")+
			      strlen(genome_fileroot)+strlen(".genes.iit")+1,sizeof(char));
    sprintf(iitfile,"%s/%s.genomes/%s.genes.iit",transcriptomesubdir,transcriptome_fileroot,genome_fileroot);

    if ((transcript_map_iit = IIT_read(iitfile,/*name*/transcriptome_fileroot,/*readonlyp*/true,/*divread*/READ_ALL,
				       /*divstring*/NULL,/*add_iit_p*/false)) != NULL) {
      FREE(iitfile);
    } else {
      fprintf(stderr,"Transcript map file %s not found.  Trying old format.\n",iitfile);
      FREE(iitfile);

      iitfile = (char *) CALLOC(strlen(genomesubdir)+strlen("/")+
				strlen(genome_fileroot)+strlen(".transcripts/")+
				strlen(transcriptome_fileroot)+strlen(".genes.iit")+1,sizeof(char));
      sprintf(iitfile,"%s/%s.transcripts/%s.genes.iit",genomesubdir,genome_fileroot,transcriptome_fileroot);
      
      if ((transcript_map_iit = IIT_read(iitfile,/*name*/transcriptome_fileroot,/*readonlyp*/true,/*divread*/READ_ALL,
				       /*divstring*/NULL,/*add_iit_p*/false)) != NULL) {
	FREE(iitfile);
      } else {
	fprintf(stderr,"Transcript map file %s not found in new format either.\n",iitfile);
	exit(9);
      }
    }

    /* transcript_map_iit (generated by iit_store) probably has genes
       in alphanumeric order But chromosome_iit may not have them in
       that order, and transcript_map_iit may be missing some
       chromosomes */
    transcript_chrnum_crosstable = Univ_IIT_divint_crosstable(chromosome_iit,transcript_map_iit);

    transcriptome = Transcriptome_new(genomesubdir,genome_fileroot,transcriptomesubdir,transcriptome_fileroot,sharedp);
    tr_indexdb = Indexdb_new_transcriptome(&index1part_tr,&index1interval_tr,transcriptomesubdir,
					   transcriptome_fileroot,/*idx_filesuffix*/"ref",/*snps_root*/NULL,
					   /*required_index1part*/0,/*required_index1interval*/0,
					   offsetsstrm_access,positions_access,sharedp,
					   multiple_sequences_p,preload_shared_memory_p,unload_shared_memory_p);
    nalignments = IIT_total_nintervals(transcript_map_iit);
    knownsplicing = Knownsplicing_from_transcriptome(transcriptome,nalignments,
						     chromosome_ef64,genomelength,/*intron_level_p*/true);
    knownsplicingp = true;
  }


  Access_setup(preload_shared_memory_p,unload_shared_memory_p);

  worker_setup(transcriptomesubdir,transcriptome_fileroot,genomesubdir,genome_fileroot,chromosome_iit);

  FREE(transcriptomesubdir);
  FREE(transcriptome_fileroot);
  /* FREE(genomesubdir); -- Needed for worker_localdb_setup */
  /* FREE(genome_fileroot); -- Needed for worker_localdb_setup */
  if (preload_shared_memory_p == true || unload_shared_memory_p == true) {
    worker_cleanup();
    return 0;
  }

  stopwatch = Stopwatch_new();

  if (two_pass_p == true) {
    /* Pass 1 */

    Stopwatch_start(stopwatch);
    outbuffer = Outbuffer_new(output_buffer_size,nread);
    Inbuffer_set_outbuffer(inbuffer,outbuffer);

#ifdef HAVE_PTHREAD
    pthread_mutex_init(&pass1_lock,NULL);
#endif

    pass = PASS1;		/* Causes worker or single thread to run process_request_pass1 */
    fprintf(stderr,"Starting pass 1.  Learning defect rate, splice sites, introns, and insert lengths (alignments are being analyzed internally, without any output)\n");
    donor_table = Univcoordtableuint_new(/*hint*/500000);
    acceptor_table = Univcoordtableuint_new(/*hint*/500000);
    antidonor_table = Univcoordtableuint_new(/*hint*/500000);
    antiacceptor_table = Univcoordtableuint_new(/*hint*/500000);

#if 0
    indel_table = Univcoordtable_new(/*hint*/500000);
#endif

#if !defined(HAVE_PTHREAD)
    /* Serial version */
    single_thread();

#else
    /* Pthreads version */
    if (nthreads == 0) {
      single_thread();

#if 0
    } else if (multiple_sequences_p == false) {
      single_thread();
#endif

    } else {
      pthread_attr_init(&thread_attr_join);
      if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
	fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
	exit(1);
      }

      worker_thread_ids = (pthread_t *) CALLOC(nthreads,sizeof(pthread_t));
      Except_init_pthread();
      pthread_key_create(&global_request_key,NULL);

      /* Create output thread */
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_pass1,
		     (void *) outbuffer);

      for (worker_id = 0; worker_id < nthreads; worker_id++) {
	/* Need to have worker threads finish before we call Inbuffer_free() */
	pthread_create(&(worker_thread_ids[worker_id]),&thread_attr_join,worker_thread,(void *) worker_id);
      }
    
      pthread_join(output_thread_id,NULL);
      for (worker_id = 0; worker_id < nthreads; worker_id++) {
	pthread_join(worker_thread_ids[worker_id],NULL);
      }

      pthread_key_delete(global_request_key);
      /* Do not delete global_except_key, because worker threads might still need it */
      /* Except_term_pthread(); */
      
      FREE(worker_thread_ids);
    }
#endif

    runtime = Stopwatch_stop(stopwatch);
    nread = Outbuffer_nread(outbuffer);
    fprintf(stderr,"Pass 1: Processed %u queries in %.2f seconds (%.2f queries/sec)\n",
	  nread,runtime,(double) nread/runtime);

#ifdef HAVE_PTHREAD
    pthread_mutex_destroy(&pass1_lock);
#endif

    defect_rate = (double) total_mismatches/(double) total_querylength;
    fprintf(stderr,"%llu mismatches/%llu aligned bp => %f defect rate\n",
	    total_mismatches,total_querylength,defect_rate);
    if (knownsplicing != NULL) {
      Knownsplicing_free(&knownsplicing);
    }
    knownsplicing = Knownsplicing_new(donor_startpoints,donor_partners,
				      acceptor_startpoints,acceptor_partners,
				      antidonor_startpoints,antidonor_partners,
				      antiacceptor_startpoints,antiacceptor_partners,
				      donor_table,acceptor_table,antidonor_table,antiacceptor_table,
				      genomelength,chromosome_iit,chromosome_ef64,
				      /*intron_level_p*/true);

#if 0
    /* Not working */
    knownindels = Knownindels_new(indel_table,genomelength,dump_indels_fp,
				  chromosome_iit,chromosome_ef64);
#endif
    max_insertlength = expected_pairlength + 5*pairlength_deviation;

    if (dump_indels_fp != NULL) {
      fclose(dump_indels_fp);
    }

    Pathpair_analyze_insertlengths(&expected_pairlength,&pairlength_deviation,insertlengths);
    max_insertlength = (int) (expected_pairlength + 5*pairlength_deviation);

    /* Stage3hr_pass2_setup(expected_pairlength,pairlength_deviation); */
    /* Pathpair_pass2_setup(expected_pairlength,pairlength_deviation); */
    Altsplice_setup(max_insertlength,genomebits,genomebits_alt);
    Stage1hr_paired_pass2_setup(max_insertlength,shortsplicedist);
    Uintlist_free(&insertlengths);

    /* TODO: Free values of indel_table */
    /* Univcoordtable_free(&indel_table); */

    Univcoordtableuint_free(&antiacceptor_table);
    Univcoordtableuint_free(&antidonor_table);
    Univcoordtableuint_free(&acceptor_table);
    Univcoordtableuint_free(&donor_table);

    Outbuffer_free(&outbuffer);
    Inbuffer_free(&inbuffer);
    pass = PASS2;

    /* Reset inbuffer and outbuffer for pass 2 */
    fastq_format_p = open_input_streams_parser(&nextchar,&nchars1,&nchars2,&paired_end_p,
					       &files,&nfiles,&input,&input2,
#ifdef HAVE_ZLIB
					       &gzipped,&gzipped2,
#endif
#ifdef HAVE_BZLIB
					       &bzipped,&bzipped2,
#endif
					       read_files_command,gunzip_p,bunzip2_p,interleavedp,argc,argv);
    Inbuffer_setup(single_cell_p,filter_if_both_p);

    inbuffer = Inbuffer_new(nextchar,input,input2,
#ifdef HAVE_ZLIB
			    gzipped,gzipped2,
#endif
#ifdef HAVE_BZLIB
			    bzipped,bzipped2,
#endif
			    interleavedp,read_files_command,files,nfiles,input_buffer_size,
			    skip_prob,/*seed*/SEED,part_modulus,part_interval);

    if ((nread = Inbuffer_fill_init(inbuffer)) > 1) {
      multiple_sequences_p = true;
    } else {
      multiple_sequences_p = false;
    }
  }
    

  /* Pass 2 (for either 1-pass or 2-pass mode) */
  worker_localdb_setup(genomesubdir,genome_fileroot);	/* localdb not used in pass 1 */
  FREE(genomesubdir);
  FREE(genome_fileroot);

  if (read_splices_file == NULL) {
    /* Without support_tables, this branch and subsequent ones must be 1-pass mode */
    /* No externally-provided splices */

    if (dump_splices_file != NULL) {
      /* Also want to record and write splice counts */
      support_tables = Support_new(knownsplicing);
    }

  } else if ((fp = fopen(read_splices_file,"r")) == NULL) {
    /* Externally-provided splices not readable.  Treat as no externally-provided splices */
    fprintf(stderr,"Unable to open --splices-read file %s\n",read_splices_file);

    if (dump_splices_file != NULL) {
      /* Also want to record and write splice counts */
      support_tables = Support_new(knownsplicing);
    }

  } else {
    support_tables = Support_read(fp,chromosome_ef64,chromosome_iit);
    fclose(fp);
	
    /* External program has already filtered splices, so convert only  */
    Support_convert(support_tables,genomelength); /* Prepare for pass 2 */
#if 0
    Splice_pass2_setup(/*for read*/support_tables);
    Trpath_convert_pass2_setup(/*for read*/support_tables);
#endif
  }


  Filestring_setup(split_simple_p);
  SAM_header_setup(argc,argv,optind,nthreads,orderedp,chromosome_iit,output_type,
		   sam_headers_p,sam_read_group_id,sam_read_group_name,
		   sam_read_group_library,sam_read_group_platform);
  Outbuffer_setup(any_circular_p,quiet_if_excessive_p,
		  paired_end_p,appendp,output_file,
		  split_simple_p,split_output_root,failedinput_root);

  Stopwatch_start(stopwatch);
  outbuffer = Outbuffer_new(output_buffer_size,nread);
  Inbuffer_set_outbuffer(inbuffer,outbuffer);

  if (output_file != NULL) {
    fprintf(stderr,"Starting alignment.  Writing results to %s\n",output_file);
  } else if (split_output_root != NULL) {
    fprintf(stderr,"Starting alignment.  Writing results to %s.*\n",split_output_root);
  } else {
    fprintf(stderr,"Starting alignment\n");
  }

#if !defined(HAVE_PTHREAD)
  /* Serial version */
  single_thread();

#else
  /* Pthreads version */
  if (nthreads == 0) {
    single_thread();

  } else if (two_pass_p == false && multiple_sequences_p == false) {
    single_thread();

  } else {
    pthread_attr_init(&thread_attr_join);
    if ((ret = pthread_attr_setdetachstate(&thread_attr_join,PTHREAD_CREATE_JOINABLE)) != 0) {
      fprintf(stderr,"ERROR: pthread_attr_setdetachstate returned %d\n",ret);
      exit(1);
    }

    worker_thread_ids = (pthread_t *) CALLOC(nthreads,sizeof(pthread_t));
    Except_init_pthread();
    pthread_key_create(&global_request_key,NULL);

    if (orderedp == true) {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_ordered,
		     (void *) outbuffer);
    } else {
      pthread_create(&output_thread_id,&thread_attr_join,Outbuffer_thread_anyorder,
		     (void *) outbuffer);
    }

    for (worker_id = 0; worker_id < nthreads; worker_id++) {
      /* Need to have worker threads finish before we call Inbuffer_free() */
      pthread_create(&(worker_thread_ids[worker_id]),&thread_attr_join,worker_thread,(void *) worker_id);
    }
    
    pthread_join(output_thread_id,NULL);
    for (worker_id = 0; worker_id < nthreads; worker_id++) {
      pthread_join(worker_thread_ids[worker_id],NULL);
    }

    pthread_key_delete(global_request_key);
    /* Do not delete global_except_key, because worker threads might still need it */
    /* Except_term_pthread(); */

    FREE(worker_thread_ids);

  }
#endif /* HAVE_PTHREAD */


  /* Note: Shortread and Sequence procedures should close their own input files */
  /* Single CPU or Pthreads version */
  runtime = Stopwatch_stop(stopwatch);
  nread = Outbuffer_nread(outbuffer);
  /* nbeyond = Outbuffer_nbeyond(outbuffer); */

  if (two_pass_p == true) {
    fprintf(stderr,"Pass 2: Aligned %u queries in %.2f seconds (%.2f queries/sec)\n",
	    nread,runtime,(double) nread/runtime);
  } else {
    fprintf(stderr,"Aligned %u queries in %.2f seconds (%.2f queries/sec)\n",
	    nread,runtime,(double) nread/runtime);
  }
  
  Stopwatch_free(&stopwatch);

  if (dump_splices_file == NULL) {
    /* Skip */
  } else if ((fp = fopen(dump_splices_file,"w")) == NULL) {
    fprintf(stderr,"Cannot write to --splices-dump file %s\n",dump_splices_file);
  } else {
    fprintf(stderr,"Dumping splices to %s\n",dump_splices_file);
    Support_dump(fp,support_tables,chromosome_ef64,chromosome_iit,
		 genomebits,dump_splices_include_known_p);
  fclose(fp);
  }

  if (support_tables != NULL) {
    Support_free(&support_tables);
  }

  Outbuffer_free(&outbuffer);
  Inbuffer_free(&inbuffer);

  Outbuffer_close_files();
  Outbuffer_cleanup();

  worker_cleanup();

  return 0;
}


static void
print_program_usage () {
  fprintf(stdout,"\
Usage: gsnap [OPTIONS...] <FASTA file>, or\n\
       cat <FASTA file> | gmap [OPTIONS...]\n\
\n\
");

  /* Input options */
  fprintf(stdout,"Input options (must include -d)\n");
  fprintf(stdout,"\
  -D, --dir=directory            Genome directory.  Default (as specified by --with-gmapdb to the configure program) is\n\
                                   %s\n\
",GMAPDB);
  fprintf(stdout,"\
  -d, --db=STRING                Genome database\n\
  --two-pass                     Two-pass mode, in which the sequences are processed first to identify splice sites\n\
                                   and introns, and then aligned using this splicing information\n\
  --use-localdb=INT              Whether to use the local suffix arrays, which help with finding extensions to the ends\n\
                                   of alignments in the presence of splicing or indels (0=no, 1=yes if available (default))\n\
\n\
");

  /* Transcriptome-guided options */
  fprintf(stdout,"Transcriptome-guided options (optional)\n");
  fprintf(stdout,"\
  -C, --transcriptdir=directory  Transcriptome directory.  Default is the value for --dir above\n\
  -c, --transcriptdb=STRING      Transcriptome database\n\
  --transcriptome-mode=STRING    Options: assist, only, annotate (default).  The option assist means\n\
                                   to try transcriptome alignment first, but then use genomic alignment\n\
                                   if nothing is found.  The option only means to try transcriptome\n\
                                   alignment only.  The option annotate means to try only genomic\n\
                                   alignment, to use the transcriptome only for annotation; this is\n\
                                   the fastest option.  In the other two options, annotation is also\n\
                                   performed\n\
\n\
");

  /* Computation options */
  fprintf(stdout,"Computation options\n");
  fprintf(stdout,"\
  -k, --kmer=INT                 kmer size to use in genome database (allowed values: 16 or less)\n\
                                   If not specified, the program will find the highest available\n\
                                   kmer size in the genome database\n\
  --sampling=INT                 Sampling to use in genome database.  If not specified, the program\n\
                                   will find the smallest available sampling value in the genome database\n\
                                   within selected k-mer size\n\
  --align-fraction=FLOAT         Process only the given fraction of reads, selected at random\n\
                                   If --align-fraction and --part are given, --align-fraction takes precedence\n\
  -q, --part=INT/INT             Process only the i-th out of every n sequences\n\
                                   e.g., 0/100 or 99/100 (useful for distributing jobs\n\
                                   to a computer farm).\n\
");
  fprintf(stdout,"\
  --input-buffer-size=INT        Size of input buffer (program reads this many sequences\n\
                                   at a time for efficiency) (default %d)\n\
",input_buffer_size);

  fprintf(stdout,"\
  --barcode-length=INT           Amount of barcode to remove from start of every read before alignment\n\
                                   (default %d)\n\
",barcode_length);

  fprintf(stdout,"\
  --endtrim-length=INT           Amount of trim to remove from the end of every read before alignment\n\
                                   (default %d)\n\
",endtrim_length);

  fprintf(stdout,"\
  --orientation=STRING           Orientation of paired-end reads\n\
                                   Allowed values: FR (fwd-rev, or typical Illumina; default),\n\
                                   RF (rev-fwd, for circularized inserts), or FF (fwd-fwd, same strand),\n\
                                   or 10X (single-cell where read 1 has barcode information; read 2 is rev)\n\
  --10x-whitelist=FILE           Whitelist of 10X Genomics GEM bead barcodes, needed to perform correction of\n\
                                   cellular barcodes.  This file can be obtained at\n\
                                   cellranger-x.y.z/lib/python/cellranger/barcodes (for Cell Ranger version >= 4)\n\
                                   cellranger-x.y.z/lib/cellranger-cs/x.y.z/lib/python/cellranger/barcodes (<= 3)\n\
  --10x-well-position=INT        Position of well information in the accession, when separated by colons\n\
                                   If set to 0, then no well information will be printed in the CB field (default: 4)\n\
  --fastq-id-start=INT           Starting position of identifier in FASTQ header, space-delimited (>= 1)\n\
  --fastq-id-end=INT             Ending position of identifier in FASTQ header, space-delimited (>= 1)\n\
                                 Examples:\n\
                                   @HWUSI-EAS100R:6:73:941:1973#0/1\n\
                                      start=1, end=1 (default) => identifier is HWUSI-EAS100R:6:73:941:1973#0\n\
                                   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36\n\
                                      start=1, end=1  => identifier is SRR001666.1\n\
                                      start=2, end=2  => identifier is 071112_SLXA-EAS1_s_7:5:1:817:345\n\
                                      start=1, end=2  => identifier is SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345\n\
  --force-single-end             When multiple FASTQ files are provided on the command line, GSNAP assumes\n\
                                    they are matching paired-end files.  This flag treats each file as single-end.\n\
  --filter-chastity=STRING       Skips reads marked by the Illumina chastity program.  Expecting a string\n\
                                   after the accession having a 'Y' after the first colon, like this:\n\
                                         @accession 1:Y:0:CTTGTA\n\
                                   where the 'Y' signifies filtering by chastity.\n\
                                   Values: off (default), either, both.  For 'either', a 'Y' on either end\n\
                                   of a paired-end read will be filtered.  For 'both', a 'Y' is required\n\
                                   on both ends of a paired-end read (or on the only end of a single-end read).\n\
  --allow-pe-name-mismatch       Allows accession names of reads to mismatch in paired-end files\n\
  --interleaved                  Input is in interleaved format (one read per line, tab-delimited\n\
");
#ifdef HAVE_ZLIB
  fprintf(stdout,"\
  --gunzip                       Uncompress gzipped input files\n\
");
#endif
#ifdef HAVE_BZLIB
  fprintf(stdout,"\
  --bunzip2                      Uncompress bzip2-compressed input files\n\
");
#endif
  fprintf(stdout,"\n");

  /* Computation options */
  fprintf(stdout,"Computation options\n");
#if 0
  fprintf(stdout,"\
\n\
  Note: GSNAP has an ultrafast algorithm for calculating mismatches up to and including\n\
((readlength+2)/kmer - 2) (\"ultrafast mismatches\").  The program will run fastest if\n\
max-mismatches (plus suboptimal-levels) is within that value.\n\
Also, indels, especially end indels, take longer to compute, although the algorithm\n\
is still designed to be fast.\n\
\n\
");
#endif

#ifdef HAVE_MMAP
  fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 5)\n\
                                 Mode  Hash offsets  Hash positions  Genome          Local hash offsets  Local hash positions  Localdb\n\
                                   0   allocate      mmap            mmap            allocate            mmap                  mmap\n\
                                   1   allocate      mmap & preload  mmap            allocate            mmap & preload        mmap\n\
                                   2   allocate      mmap & preload  mmap & preload  allocate            mmap & preload        mmap\n\
                                   3   allocate      allocate        mmap & preload  allocate            allocate              mmap\n\
                                   4   allocate      allocate        allocate        allocate            allocate              mmap\n\
                      (default)    5   allocate      allocate        allocate        allocate            allocate              allocate\n\
                           Note: For a single sequence, all data structures use mmap\n\
                           A batch level of 5 means the same as 4, and is kept only for backward compatibility\n\
");
#else
  fprintf(stdout,"\
  -B, --batch=INT                Batch mode (default = 4, modes 0-3 disallowed because program configured without mmap)\n\
                                 Mode  Hash offsets  Hash positions  Genome          Local hash offsets  Local hash positions  Localdb\n\
                      (default)    5   allocate      allocate        allocate        allocate            allocate              allocate\n\
");
#endif
  fprintf(stdout,"\
  --use-shared-memory=INT        If 1, then allocated memory is shared among all processes on this node\n\
                                   If 0 (default), then each process has private allocated memory\n\
  --preload-shared-memory        Load files indicated by --batch mode into shared memory for use by other\n\
                                   GMAP/GSNAP processes on this node, and then exit.  Ignore any input files.\n\
  --unload-shared-memory         Unload files indicated by --batch mode into shared memory, or allow them\n\
                                   to be unloaded when existing GMAP/GSNAP processes on this node are finished\n\
                                   with them.  Ignore any input files.\n\
");

  fprintf(stdout,"\
  -m, --max-mismatches=FLOAT     Maximum number of mismatches allowed (if not specified, then\n\
                                   GSNAP tries to find the best possible match in the genome)\n\
                                   If specified between 0.0 and 1.0, then treated as a fraction\n\
                                   of each read length.  Otherwise, treated as an integral number\n\
                                   of mismatches (including indel and splicing penalties).\n\
                                   Default is 0.3\n\
  --query-unk-mismatch=INT       Whether to count unknown (N) characters in the query as a mismatch\n\
                                   (0=no (default), 1=yes)\n\
  --genome-unk-mismatch=INT      Whether to count unknown (N) characters in the genome as a mismatch\n\
                                   (0=no, 1=yes).  If --use-mask is specified, default is no, otherwise yes.\n\
");

  fprintf(stdout,"\
  --maxsearch=INT                Maximum number of alignments to find (default %d).\n\
                                   Should be larger than --npaths, which is the number to report.\n\
                                   Keeping this number large will allow for random selection among multiple alignments.\n\
                                   Reducing this number can speed up the program.\n\
",maxpaths_search);

  fprintf(stdout,"\
  --indel-endlength=INT          Minimum length at end required for indel alignments (default %d)\n\
",min_indel_end_matches);

  fprintf(stdout,"\
  --max-insertions=INT           Maximum number of insertions allowed (default %d)\n\
",max_insertionlen);
  fprintf(stdout,"\
  --max-deletions=INT            Maximum number of deletions allowed (default %d)\n\
",max_deletionlen);
  fprintf(stdout,"\
  -M, --suboptimal-levels=INT    Report suboptimal hits beyond best hit (default %d)\n\
                                   All hits with best score plus suboptimal-levels are reported\n\
                                   (Note: Not currently implemented)\n\
",subopt_levels);
  fprintf(stdout,"\
  -a, --adapter-strip=STRING     Method for removing adapters from reads.  Currently allowed values: off, paired.\n\
                                   Default is \"off\".  To turn on, specify \"paired\", which removes adapters\n\
                                   from paired-end reads if they appear to be present.\n\
");

  fprintf(stdout,"\
  -e, --use-mask=STRING          Use genome containing masks (e.g. for non-exons) for scoring preference\n\
  -V, --snpsdir=STRING           Directory for SNPs index files (created using snpindex) (default is\n\
                                   location of genome index files specified using -D and -d)\n \
  -v, --use-snps=STRING          Use database containing known SNPs (in <STRING>.iit, built\n\
                                   previously using snpindex) for tolerance to SNPs\n\
  --cmetdir=STRING               Directory for methylcytosine index files (created using cmetindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --atoidir=STRING               Directory for A-to-I RNA editing index files (created using atoiindex)\n\
                                   (default is location of genome index files specified using -D, -V, and -d)\n\
  --mode=STRING                  Alignment mode: standard (default), cmet-stranded, cmet-nonstranded,\n\
                                    atoi-stranded, atoi-nonstranded, ttoc-stranded, or ttoc-nonstranded.\n\
                                    Non-standard modes requires you to have previously run the cmetindex\n\
                                    or atoiindex programs (which also cover the ttoc modes) on the genome\n\
");


#if 0
  fprintf(stdout,"\
  --tallydir=STRING              Directory for tally IIT file to resolve concordant multiple alignments (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-tally instead.\n\
  --use-tally=STRING             Use this tally IIT file to resolve concordant multiple alignments\n\
  --runlengthdir=STRING          Directory for runlength IIT file to resolve concordant multiple alignments (default is\n\
                                   location of genome index files specified using -D and -d).  Note: can\n\
                                   just give full path name to --use-runlength instead.\n\
  --use-runlength=STRING         Use this runlength IIT file to resolve concordant multiple alignments\n\
");
#endif


#ifdef HAVE_PTHREAD
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads\n\
");
#else
  fprintf(stdout,"\
  -t, --nthreads=INT             Number of worker threads.  Flag is ignored in this version of GSNAP, which has pthreads disabled\n\
");
#endif


#if 0
  fprintf(stdout,"Genes options for RNA-Seq\n");
  fprintf(stdout,"\
  -g, --genes=STRING             Look for known genes in <STRING>.iit, to be used for resolving\n\
                                   multiple mapping reads.  See README instructions for the correct\n\
                                   formatting of a genes IIT file.\n\
  --favor-multiexon              In resolving multiple mapping reads, overlaps with known\n\
                                   multi-exon genes are favored over those with known single-exon\n\
                                   genes.  This favors spliced genes over psuedogenes.\n\
");
  fprintf(stdout,"\n");
#endif


  /* Splicing options */
  fprintf(stdout,"Splicing options for DNA-Seq\n");
  fprintf(stdout,"\
  --find-dna-chimeras=INT              Look for distant splicing involving poor splice sites (0=no, 1=yes)\n\
                                         If not specified, then default is to be on unless only known splicing\n\
                                         is desired (--use-splicing is specified and --novelsplicing is off)\n\
");
  fprintf(stdout,"\n");

  /* Splicing options */
  fprintf(stdout,"Splicing options for RNA-Seq\n");
  fprintf(stdout,"\
  -N, --novelsplicing=INT              Look for novel splicing (0=no (default), 1=yes)\n\
  --splicingdir=STRING                 Directory for splicing involving known sites or known introns,\n\
                                         as specified by the -s or --use-splicing flag (default is\n\
                                         directory computed from -D and -d flags).  Note: can\n\
                                         just give full pathname to the -s flag instead.\n\
  -s, --use-splicing=STRING            Look for splicing involving known sites or known introns\n\
                                         (in <STRING>.iit), at short or long distances\n\
                                         See README instructions for the distinction between known sites\n\
                                         and known introns\n\
  --splices-noeval                     Do not evaluate splices for probability or intron length, but depend\n\
                                         only on sequence alignment\n\
  --splices-dump=FILE                  Write splice junction information to FILE, in the same format as for STAR\n\
                                         plus MaxEnt probabilities for the two intron positions.  Note that in this\n\
                                         dump file, the annotation column is reserved strictly for known introns,\n\
                                         and not novel introns that passed some criterion from a first pass.\n\
  --splices-exclude-known              In the file for --splices-dump, exclude all known introns\n\
");

  fprintf(stdout,"\
  -w, --localsplicedist=INT            Definition of local novel splicing event (default %d)\n\
",shortsplicedist);

  fprintf(stdout,"\
  --merge-distant-samechr              Report distant splices on the same chromosome as a single splice, if possible.\n\
                                         Will produce a single SAM line instead of two SAM lines, which is also done\n\
                                         for translocations, inversions, and scramble events\n\
");
  fprintf(stdout,"\n");


  /* Paired-end options */
  fprintf(stdout,"Options for paired-end reads\n");
  fprintf(stdout,"\
  --pairmax-dna=INT              Max total genomic length for DNA-Seq paired reads, or other reads\n\
                                   without splicing (default %d).  Used if -N or -s is not specified.\n\
                                   This value is also used for circular chromosomes when splicing in\n\
                                   linear chromosomes is allowed\n\
",pairmax_dna);
  fprintf(stdout,"\
  --pairmax-rna=INT              Max total genomic length for RNA-Seq paired reads, or other reads\n\
                                   that could have a splice (default %d).  Used if -N or -s is specified.\n\
                                   Should probably match the value for -w, --localsplicedist.\n\
",pairmax_rna);

  fprintf(stdout,"\
  --resolve-inner=INT            Whether to resolve soft-clipping on the insides of paired-end reads (default %d)\n\
",resolve_inner_p);

  fprintf(stdout,"\
  --pairexpect=INT               Expected paired-end length, used for resolving soft-clipping on the insides\n\
                                   of paired-end reads, and for pairing DNA-seq reads (default %d)\n\
",expected_pairlength);

  fprintf(stdout,"\
  --pairdev=INT                  Allowable deviation from expected paired-end length, used for\n\
                                   resolving soft-clipping on the insides of paired-end reads (default %d).\n\
",pairlength_deviation);

#if 0
  fprintf(stdout,"\n");
#endif
  fprintf(stdout,"\
  --pass1-min-support=INT        Threshold read support for learning an intron during pass 1 of --two-pass mode\n\
                                   (default %d)\n\
",pass1_min_support);


  /* Quality score options */
  fprintf(stdout,"Options for quality scores\n");
  fprintf(stdout,"\
  --quality-protocol=STRING      Protocol for input quality scores.  Allowed values:\n\
                                   illumina (ASCII 64-126) (equivalent to -J 64 -j -31)\n\
                                   sanger   (ASCII 33-126) (equivalent to -J 33 -j 0)\n\
                                 Default is sanger (no quality print shift)\n\
                                 SAM output files should have quality scores in sanger protocol\n\
\n\
                                 Or you can customize this behavior with these flags:\n\
  -J, --quality-zero-score=INT   FASTQ quality scores are zero at this ASCII value\n\
                                   (default is 33 for sanger protocol; for Illumina, select 64)\n\
  -j, --quality-print-shift=INT  Shift FASTQ quality scores by this amount in output\n\
                                   (default is 0 for sanger protocol; to change Illumina input\n\
                                   to Sanger output, select -31)\n\
");

  /* Output options */
  fprintf(stdout,"Output options\n");
  fprintf(stdout,"\
  -n, --npaths=INT               Maximum number of paths to print (default %d).\n\
",maxpaths_report);
  fprintf(stdout,"\
  -Q, --quiet-if-excessive       If more than maximum number of paths are found,\n\
                                   then nothing is printed.\n\
  -O, --ordered                  Print output in same order as input (relevant\n\
                                   only if there is more than one worker thread)\n\
  --show-refdiff                 For GSNAP output in SNP-tolerant alignment, shows all differences\n\
                                   relative to the reference genome as lower case (otherwise, it shows\n\
                                   all differences relative to both the reference and alternate genome)\n\
  --clip-overlap                 For paired-end reads whose alignments overlap, clip the overlapping region.\n\
  --merge-overlap                For paired-end reads whose alignments overlap, merge the two ends into a single end (beta implementation)\n\
  --print-snps                   Print detailed information about SNPs in reads (works only if -v also selected)\n\
                                   (not fully implemented yet)\n\
  --failsonly                    Print only failed alignments, those with no results\n\
  --nofails                      Exclude printing of failed alignments\n\
  --only-concordant              Print only concordant alignments (concordant_uniq, concordant_mult, concordant_circular)\n\
  --omit-concordant-uniq         Do not print any concordant_uniq alignments\n\
  --omit-concordant-mult         Do not print any concordant_mult alignments\n\
  --omit-softclipped             Do not allow any alignments with soft clips\n\
  --only-tr-consistent           Print only alignments with consistent transcripts (XX field present, identical if paired-end)\n\
");

  fprintf(stdout,"\
  -A, --format=STRING            Another format type, other than default.\n\
                                   Currently implemented: sam, m8 (BLAST tabular format)\n\
");

  fprintf(stdout,"\
  --split-output=STRING          Basename for multiple-file output, separately for nomapping,\n\
                                   halfmapping_uniq, halfmapping_mult, unpaired_uniq, unpaired_mult,\n\
                                   paired_uniq, paired_mult, concordant_uniq, and concordant_mult results\n\
  -o, --output-file=STRING       File name for a single stream of output results.\n\
  --failed-input=STRING          Print completely failed alignments as input FASTA or FASTQ format,\n\
                                    to the given file, appending .1 or .2, for paired-end data.\n\
                                    If the --split-output flag is also given, this file is generated\n\
                                    in addition to the output in the .nomapping file.\n\
  --append-output                When --split-output or --failed-input is given, this flag will append output\n\
                                    to the existing files.  Otherwise, the default is to create new files.\n\
  --order-among-best=STRING      Among alignments tied with the best score, order those alignments in this order.\n\
                                    Allowed values: genomic, random (default)\n\
");
  fprintf(stdout,"\
  --output-buffer-size=INT       Buffer size, in queries, for output thread (default %d).  When the number\n\
                                   of results to be printed exceeds this size, worker threads wait\n\
                                   until the backlog is cleared\n\
",output_buffer_size);
  fprintf(stdout,"\n");

  /* SAM options */
  fprintf(stdout,"Options for SAM output\n");
  fprintf(stdout,"\
  --no-sam-headers               Do not print headers beginning with '@'\n\
  --add-paired-nomappers         Add nomapper lines as needed to make all paired-end results alternate\n\
                                   between first end and second end\n\
  --paired-flag-means-concordant=INT  Whether the paired bit in the SAM flags means concordant only (1)\n\
                                 or paired plus concordant (0, default)\n\
  --sam-headers-batch=INT        Print headers only for this batch, as specified by -q\n\
  --sam-hardclip-use-S           Use S instead of H for hardclips\n\
  --sam-use-0M=INT               If 1 (default), then insert 0M in CIGAR between adjacent indels and introns\n\
                                   If 0, do not allow 0M.  Picard disallows 0M, but other tools may require it\n\
  --sam-extended-cigar           Use extended CIGAR format (using X and = symbols instead of M,\n\
                                   to indicate matches and mismatches, respectively\n\
  --sam-multiple-primaries       Allows multiple alignments to be marked as primary if they\n\
                                   have equally good mapping scores\n\
  --sam-sparse-secondaries       For secondary alignments (in multiple mappings), uses '*' for SEQ\n\
                                   and QUAL fields, to give smaller file sizes.  However, the output\n\
                                   will give warnings in Picard to give warnings and may not work\n\
                                   with downstream tools\n\
  --force-xs-dir                 For RNA-Seq alignments, disallows XS:A:? when the sense direction\n\
                                   is unclear, and replaces this value arbitrarily with XS:A:+.\n\
                                   May be useful for some programs, such as Cufflinks, that cannot\n\
                                   handle XS:A:?.  However, if you use this flag, the reported value\n\
                                   of XS:A:+ in these cases will not be meaningful.\n\
  --md-report-snps               In MD string, when known SNPs are given by the -v flag,\n\
                                   prints difference nucleotides when they differ\n\
                                   from reference but match a known alternate allele\n\
  --no-soft-clips                Does not allow soft clips at ends.  Mismatches will be counted over the entire query\n\
  --extend-soft-clips            Extends alignments through soft clipped regions.  CIGAR string and coordinates\n\
                                   will be revised, but mismatches and the MD string will reflect the clipped CIGAR\n\
  --action-if-cigar-error        Action to take if there is a disagreement between CIGAR length and sequence length\n\
                                   Allowed values: ignore, warning (default), noprint, abort\n\
                                   Note that the noprint option does not print the CIGAR string at all if there\n\
                                   is an error, so it may break a SAM parser\n\
  --read-group-id=STRING         Value to put into read-group id (RG-ID) field\n\
  --read-group-name=STRING       Value to put into read-group name (RG-SM) field\n\
  --read-group-library=STRING    Value to put into read-group library (RG-LB) field\n\
  --read-group-platform=STRING   Value to put into read-group library (RG-PL) field\n\
");
  fprintf(stdout,"\n");

  /* Help options */
  fprintf(stdout,"Help options\n");
  fprintf(stdout,"\
  --check                        Check compiler assumptions\n\
  --version                      Show version\n\
  --help                         Show this help message\n\
");
  return;
}

