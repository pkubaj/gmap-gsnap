/* $Id: cf2289f50ac16cc01c3309b6b84f4b65b5f0f640 $ */
#ifndef PATH_INCLUDED
#define PATH_INCLUDED

typedef struct Path_T *Path_T;

#include "bool.h"
#include "types.h"		/* For Splicetype_T */
#include "sense.h"
#include "method.h"

#include "auxinfo.h"
#include "pathstore.h"

#include "list.h"
#include "intlist.h"
#include "univcoord.h"
#include "chrnum.h"
#include "compress.h"
#include "altsplice.h"
#include "junction.h"

#include "pathstorepool.h"
#include "pathpool.h"
#include "transcriptpool.h"
#include "vectorpool.h"
#include "intlistpool.h"
#include "listpool.h"
#include "hitlistpool.h"


/* A qstart path goes has endpoints with increasing qpos, whereas a
   qend path has them in decreasing qpos.  Therefore, we can add
   segments to the head of qstart path and to the head of a qend
   path */

/* Possible fields to add: first_read_p */

#define T Path_T
struct T {
  /* MAPQ scores */
  float mapq_loglik;
  int mapq_score;
  int absmq_score;		/* Absolute MAPQ, for XQ and X2 flags */

  /* Evaluation metrics */
  int nmatches;			/* Includes all completed alignments, including altsplices, but excludes ambig splice ends */
  int ref_nmatches;
  double junction_splice_prob;	/* Excludes ambig splice ends.  Use for comparing at each locus */
  double total_splice_prob;	/* Includes ambig splice ends.  Comes into play when there are no splice junctions */

  int found_score;		/* Trims count toward score */
  int score_within_trims;	/* Ignore trims when scoring */
  char *genomic_diff;

  bool plusp;
  int genestrand;
  int sensedir;
  int querylength;

  Chrnum_T chrnum;
  Univcoord_T chroffset;
  Univcoord_T chrhigh;

 /* Used for concordance.  Avoids issues where spliceends could have
    yielded a wrong distal univdiagonal, which could prevent a
    concordance from being found */
  Univcoord_T main_univdiagonal;

  Intlist_T endpoints;
  Univcoordlist_T univdiagonals;
  Intlist_T nmismatches;
  Intlist_T ref_nmismatches;
  List_T junctions;
  
  bool splice5p;
  bool splice3p;
  Splicetype_T splicetype5;
  Splicetype_T splicetype3;
  double ambig_prob_5;
  double ambig_prob_3;

  Altsplice_T qstart_alts;
  Altsplice_T qend_alts;


  /* For circular alignments */
  bool circular_high_p;
  Intlist_T circular_endpoints;
  Univcoordlist_T circular_univdiagonals;
  Intlist_T circular_nmismatches;
  Intlist_T circular_ref_nmismatches;
  List_T circular_junctions;


  /* For a gene fusion */
  Junction_T fusion_querystart_junction;
  Junction_T fusion_queryend_junction;
  Chrnum_T fusion_chrnum;
  Univcoord_T fusion_chroffset;
  Univcoord_T fusion_chrhigh;
  bool fusion_plusp;

  Intlist_T fusion_endpoints;
  Univcoordlist_T fusion_univdiagonals;
  Intlist_T fusion_nmismatches;
  Intlist_T fusion_ref_nmismatches;
  List_T fusion_junctions;

  Altsplice_T fusion_alts;

  bool fusion_splicep;
  Splicetype_T fusion_splicetype;
  double fusion_ambig_prob;

  List_T transcripts;
  List_T invalid_transcripts;
  List_T fusion_transcripts;
  List_T fusion_invalid_transcripts;

  bool completep;
  bool childp;		     /* Result of Path_extend run on parent */
  bool extendedp;    /* Path_extend was called on this path */

  Method_T method;
  bool transcriptome_method_p;

  /* Needed for Path_eval_nmatches_chopped */
  Compress_T query_compress_fwd;
  Compress_T query_compress_rev;
  int chop_qstart;
  int chop_qend;
};


static inline Univcoord_T
Path_chroffset (T this) {
  return this->chroffset;
}

static inline Univcoord_T
Path_chrhigh (T this) {
  return this->chrhigh;
}

extern int
Path_nbadsplices (T this);

extern double
Path_avg_junction_splice_prob (T this);

extern double
Path_end_splice_prob (T this);

extern int
Path_effective_sensedir (T this);

extern int
Path_qstart (T this);

extern int
Path_qend (T this);

extern int
Path_coverage (T this);

extern Chrpos_T
Path_chrlength (T this);


extern Univcoord_T
Path_genomiclow (T this);

extern Univcoord_T
Path_genomichigh (T this);

extern int
Path_querystart (T this);

extern int
Path_queryend (T this);


#ifdef USE_HIGHLOW_UNIVDIAGONALS
static inline Univcoord_T
Path_low_univdiagonal (T this) {
  return Univcoordlist_head(this->univdiagonals);
}

static inline Univcoord_T
Path_high_univdiagonal (T this) {
  return Univcoordlist_last_value(this->univdiagonals);
}
#else
static inline Univcoord_T
Path_main_univdiagonal (T this) {
  return this->main_univdiagonal;
}
#endif


/* Does include soft-clipped region */
static inline Univcoord_T
Path_genomicstart (T this) {
  if (this->plusp == true) {
    return Univcoordlist_head(this->univdiagonals) - this->querylength;
  } else {
    return Univcoordlist_last_value(this->univdiagonals);
  }
}

/* Does include soft-clipped region */
static inline Univcoord_T
Path_genomicend (T this) {
  if (this->plusp == true) {
    return Univcoordlist_last_value(this->univdiagonals);
  } else {
    return Univcoordlist_head(this->univdiagonals) - this->querylength;
  }
}

/* Does not include soft-clipped region */
static inline Univcoord_T
Path_genomiclow_softclipped (T this) {
  return Univcoordlist_head(this->univdiagonals) - this->querylength + Intlist_head(this->endpoints);
}

static inline Univcoord_T
Path_genomiclow_fusion_softclipped (T this) {
  return Univcoordlist_head(this->fusion_univdiagonals) - this->querylength + Intlist_head(this->fusion_endpoints);
}

static inline Univcoord_T
Path_genomiclow_circular_softclipped (T this) {
  return Univcoordlist_head(this->circular_univdiagonals) - this->querylength + Intlist_head(this->circular_endpoints);
}


/* Does not include soft-clipped region */
static inline Univcoord_T
Path_genomichigh_softclipped (T this) {
  return Univcoordlist_last_value(this->univdiagonals) - this->querylength + Intlist_last_value(this->endpoints);
}


extern int
Path_interval_cmp (const void *x, const void *y);

extern int
Path_structure_cmp (const void *x, const void *y);
extern int
Path_structure_ignore_sense_cmp (const void *x, const void *y);


extern bool
Path_overlap_p (T x, T y);

static inline int
Path_qstart_trim (T this) {
  return Intlist_head(this->endpoints);
}

static inline int
Path_qend_trim (T this) {
  return this->querylength - Intlist_last_value(this->endpoints);
}

extern List_T
Path_filter (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool);

extern int
Path_setdiff_univdiagonals_auxinfo (Univcoord_T *old_univdiagonals, Auxinfo_T *old_auxinfo,
				    Pathstore_T *old_pathstores, int nold,
				    Univcoord_T *new_univdiagonals, Auxinfo_T *new_auxinfo,
				    Pathstore_T *new_pathstores, int nnew,
				    Pathstorepool_T pathstorepool, Intlistpool_T intlistpool);

extern void
Path_merge_univdiagonals_auxinfo (Univcoord_T **_all_univdiagonals, Auxinfo_T **all_auxinfo,
				  Pathstore_T **all_pathstores, int *nall,
				  Univcoord_T *_new_univdiagonals, Auxinfo_T *new_auxinfo,
				  Pathstore_T *new_pathstores, int nnew,
				  Pathstorepool_T pathstorepool, Intlistpool_T intlistpool);

extern List_T
Path_unique (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool);

extern List_T
Path_optimal_nmatches (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool);

extern List_T
Path_array_to_list (Path_T *paths, int n, Hitlistpool_T hitlistpool);

extern bool
Path_softclippedp (T this);
extern bool
Path_has_distant_splice_p (T this);
extern bool
Path_resolved_qstart_p (T this);
extern bool
Path_resolved_qend_p (T this);
extern bool
Path_unextended_qstart_p (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_unextended_qend_p (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_unextended_querystart_p (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_unextended_queryend_p (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_unextendedp (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_completep (T this, int endtrim_allowed, bool allow_ambig_p);
extern bool
Path_unsolvedp (T this);

extern int
Path_ndiffs (T this);

extern int
Path_nindels (T this);

extern unsigned int
Path_trnum_low (T this);
unsigned int
Path_trnum_high (T this);

extern int
Path_trnum_low_cmp (const void *x, const void *y);
extern int
Path_trnum_high_cmp (const void *x, const void *y);

#ifdef USE_HIGHLOW_UNIVDIAGONALS
extern int
Path_low_univdiagonal_cmp (const void *x, const void *y);
extern int
Path_high_univdiagonal_cmp (const void *x, const void *y);
#else
extern int
Path_main_univdiagonal_cmp (const void *x, const void *y);
#endif

#ifdef USE_HIGHLOW_UNIVDIAGONALS
extern int
Path_fill_low_univdiagonal (Univcoord_T *coords, int *cumsum, T *paths, int n);
extern int
Path_fill_high_univdiagonal (Univcoord_T *coords, int *cumsum, T *paths, int n);
#else
extern int
Path_fill_main_univdiagonal (Univcoord_T *coords, int *cumsum, T *paths, int n);
#endif

extern void
Path_print (T this);

#ifdef CHECK_ASSERTIONS
static inline T
Path_expect_fwd (T this) {
  int prev_endpoint;
  Intlist_T q;

  prev_endpoint = Intlist_head(this->endpoints);
  for (q = Intlist_next(this->endpoints); q != NULL; q = Intlist_next(q)) {
    /* Need to use < instead of <= because we allow repeated endpoints in an unsolved path */
    if (Intlist_head(q) < prev_endpoint) {
      printf("Expecting forward, but got\n");
      Path_print(this);
      abort();
    }
    prev_endpoint = Intlist_head(q);
  }
 
  return this;
}
#else
static inline T
Path_expect_fwd (T this) {
  return this;
}
#endif

#ifdef CHECK_ASSERTIONS
static inline T
Path_expect_rev (T this) {
  int prev_endpoint;
  Intlist_T q;

  prev_endpoint = Intlist_head(this->endpoints);
  for (q = Intlist_next(this->endpoints); q != NULL; q = Intlist_next(q)) {
    /* Need to use > instead of >= because we allow repeated endpoints in an unsolved path */
    if (Intlist_head(q) > prev_endpoint) {
      printf("Expecting reverse, but got\n");
      Path_print(this);
      abort();
    }
    prev_endpoint = Intlist_head(q);
  }

  return this;
}
#else
static inline T
Path_expect_rev (T this) {
  return this;
}
#endif

extern T
Path_reverse (T this, bool expect_fwd_p);
extern T
Path_reverse_nocheck (T this);
extern bool
Path_endpoints_acceptable_p (Intlist_T endpoints, List_T junctions);

/* Called by Trpath_convert procedures */
extern T
Path_convert_simple (Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Intlist_T endpoints, Univcoordlist_T univdiagonals, Intlist_T nmismatches,
		     Intlist_T ref_nmismatches, List_T junctions, bool evalp,
		     bool plusp, bool first_read_p, int genestrand, int sensedir, int querylength,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     Listpool_T listpool, Pathpool_T pathpool, Method_T method);
extern T
Path_create (Univcoord_T main_univdiagonal,
	     Intlist_T endpoints, Univcoordlist_T univdiagonals, Intlist_T nmismatches,
	     Intlist_T ref_nmismatches, List_T junctions,
	     bool plusp, bool first_read_p, int genestrand,
	     int sensedir, int querylength, Method_T method,
	     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
	     bool splice5p, Splicetype_T splicetype5, double ambig_prob_5,
	     bool splice3p, Splicetype_T splicetype3, double ambig_prob_3,
	     Altsplice_T qstart_alts, Altsplice_T qend_alts,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     Pathpool_T pathpool, Vectorpool_T vectorpool);
extern T
Path_copy (T old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
	   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool);
extern T
Path_copy_5 (T old, bool splice5p, Splicetype_T splicetype5, double ambig_prob_5,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool);
extern T
Path_copy_3 (T old,  bool splice3p, Splicetype_T splicetype3, double ambig_prob_3,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool);
extern int
Path_exon_origin (T this);
extern void
Path_count (int *npaths_primary, int *npaths_altloc, List_T paths);

/* Almost components of path, including junctions, altsplice, and
   genomic_diff are allocated by Pathpool.  However, components of altsplice are allocated. */
extern void
Path_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	   Hitlistpool_T hitlistpool);

/* All components of path, including junctions, altsplice, and genomic_diff are allocated by Pathpool */
extern void
Path_gc (List_T *list, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	 Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	 Hitlistpool_T hitlistpool);

extern void
Path_free_parent (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool);

extern void
Path_array_gc (T *paths, int n, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool);

extern T
Path_new_from_ends (Univcoord_T univdiagonal5, int qstart5, int qend5,
		    Univcoord_T univdiagonal3, int qstart3, int qend3,
		    bool plusp, bool first_read_p, int genestrand,
		    int sensedir, int querylength, Method_T method,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool);

extern T
Path_new_exact (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches, int ref_nmismatches,
		bool plusp, bool first_read_p, int genestrand, int querylength, int found_score,
		Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		Compress_T query_compress_fwd, Compress_T query_compress_rev,
		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Pathpool_T pathpool,
		Method_T method);

extern T
Path_new_for_qstart_extension (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches,
			       bool plusp, bool first_read_p, int genestrand,
			       int sensedir, int querylength, Method_T method,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			       bool splice5p, Splicetype_T splicetype5, double ambig_prob_5, 
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			       Pathpool_T pathpool);
extern T
Path_new_for_qend_extension (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches,
			     bool plusp, bool first_read_p, int genestrand,
			     int sensedir, int querylength, Method_T method,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool splice3p, Splicetype_T splicetype3, double ambig_prob_3, 
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Pathpool_T pathpool);

extern bool
Path_segments_good_p (T this);

extern void
Path_check_valid (T this);

extern void
Path_setup (bool *circularp_in, bool *altlocp_in);

#undef T
#endif

