/* $Id: e7b440f9b937f73f5d6b5e1db1a3c38d30ddd816 $ */
#ifndef SPLICE_INCLUDED
#define SPLICE_INCLUDED

typedef struct Spliceinfo_T *Spliceinfo_T;

#include "bool.h"
#include "types.h"
#include "genomicpos.h"
#include "compress.h"
#include "genomebits.h"
#include "knownsplicing.h"
#include "univcoord.h"
#include "univdiag.h"
#include "indel.h"
#include "stage1hr.h"
#include "localdb-read.h"
#include "support.h"


#define KNOWN_SITE 0
#define TYPE1_SITE -1		/* GT|GC - AG */
#define TYPE2_SITE -2		/* AT - AC */


#define SMALL_SPLICE 1000	/* For these distances and smaller, accept splice with only one good prob */
#define MEDIUM_SPLICE 10000
#define HUGE_SPLICE 30000

#define SPLICEINDEL_PROB 0.8
/* #define SPLICE_PROB_HIGH 0.9 */
#define SPLICE_PROB_LOW 0.2 /* Used by Splice_fusion_sense and Splice_fusion_antisense */

#define T Spliceinfo_T
struct T {
  int *mismatch_positions_left1;
  int *mismatch_positions_right1;
  int *mismatch_positions_left2;
  int *mismatch_positions_right2;

  int *segmenti_sites_alloc1;
  int *segmenti_types_alloc1;
  int *segmenti_sites_alloc2;
  int *segmenti_types_alloc2;
  double *segmenti_probs_alloc;

  int *segmentk1_sites_alloc1;
  int *segmentk1_types_alloc1;
  int *segmentk1_sites_alloc2;
  int *segmentk1_types_alloc2;
  double *segmentk1_probs_alloc;

  int *segmentk2_sites_alloc1;
  int *segmentk2_types_alloc1;
  int *segmentk2_sites_alloc2;
  int *segmentk2_types_alloc2;
  double *segmentk2_probs_alloc;

  int *segmentj_sites_alloc1;
  int *segmentj_types_alloc1;
  int *segmentj_sites_alloc2;
  int *segmentj_types_alloc2;
  double *segmentj_probs_alloc;
};


extern void
Splice_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in,
	      int max_insertionlen_in, int max_deletionlen_in,
	      bool novelsplicingp_in, bool splice_evalp_in);

extern void
Spliceinfo_free (T *old);

extern T
Spliceinfo_new (int querylength);

extern bool
Splice_accept_p (int splice_qpos, int querylength,
		 double splice_prob_1, double splice_prob_2,
		 Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		 int supporti, int nmismatchesi, int supportj, int nmismatchesj,
		 bool plusp, bool sense_forward_p, bool innerp);

extern int
Splice_resolve (int *best_type, int *typei, int *typej, int *trimpos5, int *trimpos3,
		Univcoord_T *middle_univdiagonal, int *splice_qpos_i, int *splice_qpos_j,
		int *best_nindels, int *best_indel_pos,
		int *best_nmismatches_i, int *nmismatches_middle, int *best_nmismatches_j,
		int *best_nmismatches_indel,

		int *best_ref_nmismatches_i, int *ref_nmismatches_middle, int *best_ref_nmismatches_j,
		int *best_ref_nmismatches_indel,
		double *best_donor1_prob, double *best_acceptor1_prob,
		double *best_donor2_prob, double *best_acceptor2_prob,

		Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		Stage1_T stage1, Compress_T query_compress, char *queryptr,
		bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		Localdb_T localdb, int localdb_nmismatches_allowed,

		int pos5, int pos3, int querylength,
		Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		Intlistpool_T intlistpool,
		bool sense_forward_p, int genestrand, bool trim5p, bool trim3p, bool innerp,
		bool allow_indels_p);

extern int
Splice_resolve_qstart (Univcoordlist_T *univdiagonals, Intlist_T *nmismatches, Intlist_T *ref_nmismatches,
		       Intlist_T *endpoints, List_T *junctions,

		       Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		       Stage1_T stage1, Compress_T query_compress, char *queryptr,
		       bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		       int streamspace_max_alloc, Univcoord_T *streamspace_alloc,

		       Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		       Localdb_T localdb, int localdb_nmismatches_allowed,

		       int pos5, int pos3, int querylength, int cont_trimpos, int cont_nmismatches,
		       Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		       Univcoordlistpool_T univcoordlistpool, Intlistpool_T intlistpool,
		       Listpool_T listpool, Pathpool_T pathpool,

		       int sensedir, int genestrand, bool trim5p, bool trim3p, bool innerp,
		       bool allow_indels_p, bool find_middle_p);

extern int
Splice_resolve_qend (Univcoordlist_T *univdiagonals, Intlist_T *nmismatches, Intlist_T *ref_nmismatches,
		     Intlist_T *endpoints, List_T *junctions,

		     Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		     Stage1_T stage1, Compress_T query_compress, char *queryptr,
		     bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,
		     
		     int streamspace_max_alloc, Univcoord_T *streamspace_alloc,
		     Univdiag_T *univdiags_alloc, Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		     Localdb_T localdb, int localdb_nmismatches_allowed,

		     int pos5, int pos3, int querylength, int cont_trimpos, int cont_nmismatches,
		     Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     Univcoordlistpool_T univcoordlistpool, Intlistpool_T intlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,

		     int sensedir, int genestrand, bool trim5p, bool trim3p, bool innerp,
		     bool allow_indels_p, bool find_middle_p);

extern int
Splice_nomiddle (int *best_type, int *trimpos5, int *trimpos3, int *supporti, int *supportj,
		 int *best_nindels, int *best_indel_pos,
		 int *best_nmismatches_i, int *best_nmismatches_j,
		 int *best_nmismatches_indel,

		 int *best_ref_nmismatches_i, int *best_ref_nmismatches_j,
		 int *best_ref_nmismatches_indel,
		 double *best_donor_prob, double *best_acceptor_prob,

		 Univcoord_T univdiagonal_i, Univcoord_T univdiagonal_j,
		 Compress_T query_compress,
		 bool plusp, Univcoord_T chroffset, Univcoord_T chrhigh,

		 int pos5, int pos3, int querylength,
		 Indelinfo_T indelinfo, Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		 Intlistpool_T intlistpool,
		 bool sense_forward_p, int genestrand, bool trim5p, bool trim3p, bool innerp,
		 bool allow_indels_p);

extern int
Splice_fusion_sense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
		     double *best_donor_prob, double *best_acceptor_prob,

		     int *best_nmismatches_D, int *best_nmismatches_A,
		     int *best_ref_nmismatches_D, int *best_ref_nmismatches_A,
		       
		     Univcoord_T univdiagonalD, Univcoord_T univdiagonalA,
		     Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
		     Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
		     
		     int queryposD, int queryposA, int querylength,
		     Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
		     Intlistpool_T intlistpool, int genestrand);


extern int
Splice_fusion_antisense (int *best_type, char *donor1, char *donor2, char *acceptor1, char *acceptor2,
			 double *best_donor_prob, double *best_acceptor_prob,

			 int *best_nmismatches_A, int *best_nmismatches_D,
			 int *best_ref_nmismatches_A, int *best_ref_nmismatches_D,
		       
			 Univcoord_T univdiagonalA, Univcoord_T univdiagonalD,
			 Compress_T query_compress_A, bool plusAp, Univcoord_T chroffset_A,
			 Compress_T query_compress_D, bool plusDp, Univcoord_T chroffset_D,
			 
			 int queryposA, int queryposD, int querylength,
			 Spliceinfo_T spliceinfo, Knownsplicing_T knownsplicing,
			 Intlistpool_T intlistpool, int genestrand);

#undef T
#endif

