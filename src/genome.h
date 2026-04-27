/* $Id: genome.h 225310 2022-11-09 17:01:06Z twu $ */
#ifndef GENOME_INCLUDED
#define GENOME_INCLUDED

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For size_t, off_t */
#endif

#include "bool.h"
#include "access.h"
#include "types.h"
#include "genomicpos.h"
#include "iit-read-univ.h"
#include "chrnum.h"
#include "mode.h"
#include "univcoord.h"

#ifndef GSNAP
#include "sequence.h"
#endif

#define OUTOFBOUNDS '*'

/* typedef enum {GENOME_OLIGOS, GENOME_BITS} Genometype_T; */

#define T Genome_T
struct T {
  Access_T access;		/* For blocks: ALLOCATED_PRIVATE, ALLOCATED_SHARED, or MMAPPED */
  int blocks_shmid;
  key_t blocks_key;

  int fd;
  size_t len;			/* filesize */
  Univcoord_T genomelength;	/* from chromosome_iit */

  Genomecomp_T *blocks;
  char *accession;		/* for a usersegment */
};

typedef struct T *T;

extern void
Genome_free (T *old);
extern Genomecomp_T *
Genome_blocks (T this);
extern Univcoord_T
Genome_genomelength (T this);
extern char *
Genome_accession (T this);
#if 0
extern void
Genome_shmem_remove (char *genomesubdir, char *fileroot, char *alt_suffix, Genometype_T genometype,
		     bool genome_lc_p);
#endif
extern T
Genome_new (char *genomesubdir, char *fileroot, char *alt_suffix,
	    Univ_IIT_T chromosome_iit, Access_mode_T access, bool sharedp, bool revcompp);

#ifndef GSNAP
extern T
Genome_from_sequence (Sequence_T sequence);
#endif

extern void
Genome_setup (
#ifdef GSNAP
	      Genome_T genome_in, Genome_T genomealt_in,
#endif
	      int circular_typeint_in);

#ifdef UTILITYP
extern void
Genome_dump (T genome, Univcoord_T startpos, Univcoord_T endpos, int wraplength);
#endif

extern void
Genome_uncompress_mmap (char *gbuffer1, Genomecomp_T *blocks, Univcoord_T startpos, 
			Univcoord_T endpos);
extern void
Genome_uncompress_memory (char *gbuffer1, Genomecomp_T *blocks, Univcoord_T startpos, 
			  Univcoord_T endpos);
#if 0
extern bool
Genome_fill_buffer_old (Chrnum_T *chrnum, int *nunknowns, T this, Univcoord_T left, Chrpos_T length, char *gbuffer1,
			Univ_IIT_T chromosome_iit);
#endif

extern void
Genome_fill_buffer (
#if !defined(GSNAP)
		    T genome,
#endif
		    Univcoord_T left, Chrpos_T length, char *gbuffer1);

/* Used by sarray-write.c */
extern void
Genome_fill_buffer_simple (T this, Univcoord_T left, Chrpos_T length, char *gbuffer1);


#if !defined(GEXACT) && !defined(GSNAP) && !defined(UTILITYP)
extern void
Genome_fill_buffer_blocks_noterm (T genome, T genomealt, Univcoord_T left, Chrpos_T length, char *gbuffer1, char *gbuffer2);
#endif

extern void
Genome_fill_buffer_alt (
#ifndef GSNAP
			T genome, T genomealt,
#endif
			Univcoord_T left, Chrpos_T length, char *gbuffer1);

/* Used by sarray-write.c */
extern void
Genome_fill_buffer_int_string (T this, Univcoord_T left, Chrpos_T length, unsigned char *intstring,
			       unsigned char *conversion);
/* Used by sarray-write.c */
extern char
Genome_get_char_lex (T this, Univcoord_T left, Univcoord_T genomelength, char chartable[]);

extern char
Genome_get_char (char *charalt,
#if !defined(GSNAP)
		 T genome, T genomealt,
#endif
		 Univcoord_T left);

#if !defined(GEXACT) && !defined(GSNAP)
extern void
Genome_get_segment_right (char *segment, char *segmentalt, T genome, T genomealt,
			  Univcoord_T left, Chrpos_T length, Univcoord_T chrhigh, bool revcomp);
extern void
Genome_get_segment_left (char *segment, char *segmentalt, T genome, T genomealt,
			 Univcoord_T left, Chrpos_T length, Univcoord_T chroffset, bool revcomp);
#endif

#if !defined(GSNAP)
extern Sequence_T
Genome_get_segment (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
		    bool revcomp);
extern Sequence_T
Genome_get_segment_alt (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
			bool revcomp);
extern Sequence_T
Genome_get_segment_snp (T this, Univcoord_T left, Chrpos_T length, Univ_IIT_T chromosome_iit,
			bool revcomp);
#endif

#undef T
#endif
