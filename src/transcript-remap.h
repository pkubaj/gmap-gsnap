/* $Id: cc06ef8fefcd502a3d2b92394506fd2be8b568b8 $ */
#ifndef TRANSCRIPT_REMAP_INCLUDED
#define TRANSCRIPT_REMAP_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "bool.h"
#include "types.h"
#include "chrnum.h"
#include "intlist.h"
#include "uintlist.h"
#include "univcoord.h"
#include "list.h"
#include "path.h"
#include "shortread.h"

#include "uintlistpool.h"
#include "transcriptpool.h"
#include "listpool.h"
#include "vectorpool.h"


extern Uintlist_T
Transcript_remap_endpoint_coords (Intlist_T endpoints, Univcoordlist_T univdiagonals, List_T junctions,
				  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, int querylength,
				  Uintlistpool_T uintlistpool, bool extend_qstart_p, bool extend_qend_p);

extern List_T
Transcript_remap_geneplus (int *overall_trstart, int *overall_trend,
			   int *trstart_overhang, Chrpos_T *trstart_splice_distance,
			   int *trend_overhang, Chrpos_T *trend_splice_distance,
			   Shortread_T queryseq, bool plusp,
			   Uintlist_T coords, int *exonbounds, Chrpos_T *exonstarts,
			   int transcript_nexons, Transcriptpool_T transcriptpool, Listpool_T listpool);
extern List_T
Transcript_remap_geneminus (int *overall_trstart, int *overall_trend,
			    int *trstart_overhang, Chrpos_T *trstart_splice_distance,
			    int *trend_overhang, Chrpos_T *trend_splice_distance,
			    Shortread_T queryseq, bool plusp,
			    Uintlist_T coords, int *exonbounds, Chrpos_T *exonstarts,
			    int transcript_nexons, Transcriptpool_T transcriptpool, Listpool_T listpool);

extern bool
Transcript_remap_invalid (Transcript_T transcript, Path_T path, Transcriptome_T transcriptome,
			  Shortread_T queryseq, Uintlistpool_T uintlistpool, Listpool_T listpool,
			  Transcriptpool_T transcriptpool);

extern List_T
Transcript_remap_all (List_T *transcripts, List_T *invalid_transcripts,
		      Intlist_T endpoints, Univcoordlist_T univdiagonals, List_T junctions,
		      Shortread_T queryseq, int querylength, bool plusp,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      Uintlistpool_T uintlistpool, Listpool_T listpool,
		      Transcriptpool_T transcriptpool,
		      int desired_genestrand, bool extend_qstart_p, bool extend_qend_p,
		      bool repairp);

extern bool
Transcript_remap_matchp (Chrpos_T low_chrpos, Chrpos_T high_chrpos, Chrnum_T chrnum);

extern void
Transcript_remap_setup (bool *circularp_in,
			Transcriptome_T transcriptome_in,
			IIT_T transcript_map_iit_in,
			int *transcript_chrnum_crosstable_in);

#endif


