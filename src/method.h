/* $Id: e988fc99b690afe19695fc2fc8bbcafd4509c126 $ */
#ifndef METHOD_INCLUDED
#define METHOD_INCLUDED

#include "filestring.h"


/* Keep in order of calls to single_read in stage1hr-single.c and paired_read in stage1hr-paired.c */
/* Separate categories for tr, RNA, and DNA alignment allow us to order each task separately */
typedef enum {METHOD_INIT,
	      TR_EXACT1, /* TR   AUX              Transcriptome_exact1 */
	      TR_EXACT2, /* TR   AUX              Transcriptome_exact2 */
	      TR_ANYPAIR,	/* TR   AUX              Transcriptome_anypair */
	      TR_PREVALENT,	/* TR   AUX              Transcriptome_prevalent */

	      TR_END,		/* TR   AUX              Transcriptome_search_end */
	      TR_EXT,         	/* TR   AUX              Tr_extension_search */
	      TR_CHAIN,		/* TR   AUX              Trchain_search */
	      TR_CHAIN_K13,	/* TR   AUX              Trchain_search, k13 seeds */
	      TR_CHAIN_K11,	/* TR   AUX              Trchain_search, k11 rescue seeds */
	      TR_CHAIN_FM_TRDIAG,	/* TR   AUX              Trchain_search, FM seeds + TR diag chain */
	      TR_CHAIN_FM_LCHAIN,	/* TR   AUX              Trchain_search, FM seeds + FM chain */

	      TR_METHOD,

	      KMER_EXACT1,	/* TR         RNA   DNA  Kmer_exact1 */
	      KMER_ANYPAIR,	/* TR         RNA   DNA  Kmer_anypair */	
	      KMER_EXACT2,	/*            RNA   DNA  Kmer_exact2 */

	      EXT,  	        /* TR         RNA   DNA  Extension_search (for concordance) */
	      EXT_NOSPLICE,	/* TR         RNA   DNA  Extension_search */
	      EXT_SPLICE,	/* TR         RNA   DNA  Extension_search */

	      SEGMENT1,		/* TR         RNA   DNA  Kmer_segment */

	      KMER_COMPLETE,	/*                  DNA  */
	      KMER_APPROX,	/*            RNA   DNA  Kmer_search_approx */

	      KMER_PREVALENT,	/* TR         RNA   DNA  Kmer_prevalent */
	      KMER_ANCHORED,	/* TR         RNA   DNA  Kmer_anchored_5, Kmer_anchored_3 */
	      KMER_PAIRED,	/* TR         RNA   DNA  paired search */
	      KMER_WIDEST,	/* TR         RNA   DNA  Kmer_widest */
	      END,	        /* TR         RNA   DNA  Kmer_search_end */
	      TR_COMPLETE,	/* TR   AUX              Transcriptome_search_complete */
	      TR_APPROX,	/* TR   AUX              Transcriptome_search_ends */

	      LOCAL_MATE,
	      EXHAUSTIVE,

	      FUSION,	      /* fusion, from path-fusion */

	      NMETHODS} Method_T;

extern void
Method_dump (bool *methodp);
extern char *
Method_string (Method_T method);
extern void
Method_samprint (Filestring_T fp, Method_T method);
extern void
Method_print (Filestring_T fp, Method_T method);

#endif
