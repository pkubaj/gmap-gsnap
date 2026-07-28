/* $Id: 54f0ad3e1cd222a84798316000428b44bd9842a6 $ */
#ifndef ALTSPLICE_INCLUDED
#define ALTSPLICE_INCLUDED

typedef struct Altsplice_T *Altsplice_T;

#include "bool.h"
#include "univcoord.h"
#include "univdiag.h"
#include "pathpool.h"
#include "vectorpool.h"
#include "genomebits.h"
#include "compress.h"

#define GOOD_DISTAL_PROB 0.8

#define T Altsplice_T
struct T {

  bool boundedp;
  bool innerp;
  bool nmismatches_known_p;

  /* Because of combine_leftright_paths, the medial_nmismatches might not be valid for the combined path */
  int splice_qpos;
  int anchor_qpos;		/* The qpos more medial than splice_qpos; beginning of the segment */

  /* Needed by Path_eval_nmatches */
  int best_distali;

  int nunivdiagonals;
  Univcoord_T *univdiagonals;

  /* Assume a single medial splice_qpos */
  Univcoord_T medial_univdiagonal;
  int medial_nmismatches;
  double medial_prob;
  int medial_support;

  /* Multiple distal solutions */
  int *distal_lengths;
  int *distal_qpos;
  int *distal_trimpos;
  int *distal_nmismatches;
  double *distal_probs;
  int *distal_support;

  /* int *ref_nmismatches; */
};


static inline bool
Altsplice_boundedp (T this) {
  return this->boundedp;
}

/* Needed by Path_genomiclow and Path_genomichigh */
static inline Univcoord_T
Altsplice_best_univdiagonal (T this) {
  return this->univdiagonals[this->best_distali];
}

static inline double
Altsplice_best_distal_prob (T this) {
  return this->distal_probs[this->best_distali];
}

static inline int
Altsplice_best_nmatches (T this) {
  return this->distal_support[this->best_distali] - this->distal_nmismatches[this->best_distali];
}

static inline int
Altsplice_found_score (T this) {
  return this->distal_nmismatches[this->best_distali];
}

extern void
Altsplice_free (T *old, Pathpool_T pathpool);

extern void
Altsplice_print (T this);

extern T
Altsplice_copy (T old, Pathpool_T pathpool, Vectorpool_T vectorpool);

extern Univcoord_T
Altsplice_qstart_best (int *distal_trimpos, int splice_qpos, int querylength,
		       Univcoord_T medial_univdiagonal,
		       double medial_prob, int medial_support, int medial_nmismatches,
		       int *distal_support, Univdiag_T *distal_univdiags,
		       double *distal_probs, int npartners,
		       bool plusp, bool sense_forward_p, bool innerp);

extern Univcoord_T
Altsplice_qend_best (int *distal_trimpos, int splice_qpos, int querylength,
		     Univcoord_T medial_univdiagonal,
		     double medial_prob, int medial_support, int medial_nmismatches,
		     int *distal_support, Univdiag_T *distal_univdiags,
		     double *distal_probs, int npartners,
		     bool plusp, bool sense_forward_p, bool innerp);

extern T
Altsplice_qstart_new (int *endpoint, bool boundedp, int splice_qpos, int anchor_qpos, int querylength,
		      Univcoord_T medial_univdiagonal, int medial_nmismatches,
		      double medial_prob, int medial_support, Univcoord_T *distal_positions,
		      int *distal_support, Univdiag_T *distal_univdiags,
		      double *distal_probs, int *order_medial_qpos, int npartners,
		      bool plusp, bool sense_forward_p, bool innerp,
		      Pathpool_T pathpool, Vectorpool_T vectorpool, bool sort_bydistal_p);

extern T
Altsplice_qend_new (int *endpoint, bool boundedp, int splice_qpos, int anchor_qpos, int querylength,
		    Univcoord_T medial_univdiagonal, int medial_nmismatches,
		    double medial_prob, int medial_support, Univcoord_T *distal_positions,
		    int *distal_support, Univdiag_T *distal_univdiags,
		    double *distal_probs, int *order_medial_qpos, int npartners,
		    bool plusp, bool sense_forward_p, bool innerp,
		    Pathpool_T pathpool, Vectorpool_T vectorpool, bool sort_bydistal_p);

extern Univcoord_T
Altsplice_first_univdiagonal (T this);

extern Univcoord_T
Altsplice_last_univdiagonal (T this);

extern bool
Altsplice_trim_qstart_chrbounds (T this, Univcoord_T chroffset, int querylength);

extern bool
Altsplice_trim_qend_chrbounds (T this, Univcoord_T chrhigh);

extern bool
Altsplice_resolve_qend (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L, 
			double *medial_prob_L, double *distal_prob_L,
			T this, int anchor_qpos_L, int querylengthL, int querylengthH,
			Univcoord_T genomicstartH, Compress_T query_compress,
			bool plusp, int genestrand, bool sense_forward_p);
extern bool
Altsplice_resolve_qstart (Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			  int *medial_nmismatches_L, int *distal_nmismatches_H,
			  double *medial_prob_H, double *distal_prob_H,
			  T this, int anchor_qpos_H, int querylengthL, int querylengthH,
			  Univcoord_T genomicendL, Compress_T query_compress,
			  bool plusp, int genestrand, bool sense_forward_p);
extern bool
Altsplice_resolve_both (Univcoord_T *univdiagonal_L, int *splice_qpos_L, int *distal_trimpos_L,
			int *medial_nmismatches_L, int *distal_nmismatches_L,
			double *medial_prob_L, double *distal_prob_L,

			Univcoord_T *univdiagonal_H, int *splice_qpos_H, int *distal_trimpos_H,
			int *medial_nmismatches_H, int *distal_nmismatches_H,
			double *medial_prob_H, double *distal_prob_H,

			T thisL, int anchor_qpos_L, T thisH, int anchor_qpos_H,
			int querylengthL, int querylengthH,
			Compress_T query_compress_L, Compress_T query_compress_H,
			bool plusp, int genestrand, bool sense_forward_p);

extern int
Altsplice_select_qend (T this, int querylength, Compress_T query_compress, bool plusp, int genestrand);

extern int
Altsplice_select_qstart (T this, int querylength, Compress_T query_compress, bool plusp, int genestrand);

extern void
Altsplice_setup (int max_insertlength_in, Genomebits_T genomebits_in,
		 Genomebits_T genomebits_alt_in);

#undef T
#endif


