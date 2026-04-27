/* $Id: 5e2b52b93c49e3c499d2ab51795af2b6218a4d37 $ */
#ifndef PATH_TRIM_INCLUDED
#define PATH_TRIM_INCLUDED

#include "path.h"
#include "compress.h"
#include "genomebits.h"

#include "intlistpool.h"
#include "univcoord.h"
#include "listpool.h"
#include "pathpool.h"
#include "transcriptpool.h"

#define T Path_T

extern bool
Path_trim_qstart_trimdiag (T this,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool,
			   Transcriptpool_T transcriptpool, Univcoord_T trimdiag);

extern bool
Path_trim_qend_trimdiag (T this,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool,
			 Transcriptpool_T transcriptpool, int querylength,
			 Univcoord_T trimdiag);

extern void
Path_trim_chrbounds (T this,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,
		     Transcriptpool_T transcriptpool);

extern void
Path_trim_circular_unalias (T this);
extern void
Path_trim_circular_unalias_pair (T pathL, T pathH);
extern void
Path_trim_circular (T this, Compress_T query_compress,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool);
extern void
Path_trim_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in);


#undef T
#endif


