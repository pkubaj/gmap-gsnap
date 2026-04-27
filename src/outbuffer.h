/* $Id: aabce3f4e51654ac045975d907a462c8d5876ad6 $ */
#ifndef OUTBUFFER_INCLUDED
#define OUTBUFFER_INCLUDED

#include "types.h"
#include "bool.h"
#include "genomicpos.h"
#include "filestring.h"

#define T Outbuffer_T
typedef struct T *T;

#if defined(GFILTER)
extern void
Outbuffer_setup (bool paired_end_p_in, bool appendp_in, char *output_root_in);
#else  /* GSNAP, GEXACT, or GMAP */
extern void
Outbuffer_setup (bool any_circular_p_in, bool quiet_if_excessive_p_in,
		 bool paired_end_p_in, bool appendp_in,
		 char *output_file_in, bool split_simple_p_in, char *split_output_root_in,
		 char *failedinput_root_in);
#endif


extern void
Outbuffer_cleanup ();

extern T
Outbuffer_new (unsigned int output_buffer_size, unsigned int nread);

extern void
Outbuffer_close_files ();

extern void
Outbuffer_free (T *old);

extern unsigned int
Outbuffer_nread (T this);

#if defined(GFILTER)
extern unsigned int
Outbuffer_npassed (T this);
#endif

extern void
Outbuffer_add_nread (T this, unsigned int nread);

#if defined(GFILTER)
extern void
Outbuffer_put_filestrings (T this, int request_id, Filestring_T fp_reads_1, Filestring_T fp_reads_2);

extern void
Outbuffer_print_filestrings (T this, Filestring_T fp_reads_1, Filestring_T fp_reads_2);

#elif defined(GSNAP)
extern void
Outbuffer_put_pass1 (T this, int request_id);
extern void
Outbuffer_put_filestrings (T this, int request_id, Filestring_T fp, Filestring_T fp_failedinput,
			   Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2);

extern void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput, Filestring_T fp_failedinput_1, Filestring_T fp_failedinput_2);

#else  /* GEXACT or GMAP */
extern void
Outbuffer_put_filestrings (T this, int request_id, Filestring_T fp, Filestring_T fp_failedinput);

extern void
Outbuffer_print_filestrings (Filestring_T fp, Filestring_T fp_failedinput_1);

#endif

#ifdef GSNAP
extern void *
Outbuffer_thread_pass1 (void *data);
#endif

extern void *
Outbuffer_thread_anyorder (void *data);

extern void *
Outbuffer_thread_ordered (void *data);


#undef T
#endif

