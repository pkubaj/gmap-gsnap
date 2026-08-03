#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-trim.h"
#include "debug_groups.h"
#include "path-eval.h"

#include <stdio.h>
#include <math.h>		/* For rint */

#include "assert.h"
#include "intlist.h"
#include "univcoord.h"
#include "list.h"
#include "genomebits_count.h"

#include "altsplice.h"
#include "junction.h"


static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;


/* #define MIN_END 4 */

#if defined(DEBUG) || defined(DEBUG_PATH_CORE)
#define debug(x) x
#else
#define debug(x)
#endif

/* Path_trim_circular */
#if defined(DEBUG2) || defined(DEBUG_PATH_CORE)
#define debug2(x) x
#else
#define debug2(x)
#endif

/* unalias_path */
#if defined(DEBUG3) || defined(DEBUG_PATH_CORE)
#define debug3(x) x
#else
#define debug3(x)
#endif


/* TODO: Handle alt splices that go outside of chromosome bounds */
/* TODO: Trim fusion paths */

#define T Path_T

#if 0
static int
count_outside_qstart_byposition (int *nsegments, int *remainder,
				 Intlist_T endpoints, Univcoordlist_T univdiagonals,
				 Univcoord_T trimbounds, int querylength) {
  int noutside = 0;
  Univcoord_T left, univdiagonal;
  int qstart, qend;

  *nsegments = 0;
  *remainder = 0;

  debug(printf("Entered count_outside_qstart_byposition\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    left = univdiagonal - querylength;
    qstart = Intlist_head(endpoints);
    qend = Intlist_second_value(endpoints);
    debug(printf("univdiagonal %u, left %u (%u..%u) vs trimbounds %u\n",
		 univdiagonal,left,left+qstart,left+qend,trimbounds));

    if (left + qstart >= trimbounds) {
      *remainder = 0;
      return noutside;

    } else if (left + qend >= trimbounds) {
      *remainder = (int) (trimbounds - (left + qstart));
      return noutside + (*remainder);

    } else {
      *nsegments += 1;
      noutside += qend - qstart;
      endpoints = Intlist_next(endpoints);
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return noutside;
}
#endif


static int
count_outside_qstart_bydiagonal (Univcoordlist_T univdiagonals, Univcoord_T trimdiag) {
  int nsegments = 0;
  Univcoord_T univdiagonal;

  debug(printf("Entered count_outside_qstart_bydiagonal\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    /* Allow the same univdiagonal */
    if (univdiagonal >= trimdiag) {
      return nsegments;

    } else {
      nsegments += 1;
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return nsegments;
}


/* nmatches_trimmed is valid only if return value (evalp) is false */
static bool
trim_qstart (int *nmatches_trimmed, int nsegments, int remainder, T this,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  bool evalp = false;
  int i;
  /* int qstart, qend; */
  int ignore_int, nmismatches, endpoint, ninserts;
  Univcoord_T ignore_univcoord;
  Junction_T junction;
  Junctiontype_T type;

    
  debug(printf("Path needs trimming at qstart (%d segments, %d remainder)\n",
	       nsegments,remainder));
  debug(Path_print(this));

  *nmatches_trimmed = 0;

  if (nsegments > 0 && this->qstart_alts != NULL) {
    Altsplice_free(&this->qstart_alts,pathpool);
    this->qstart_alts = (Altsplice_T) NULL;
    debug(printf("(1) Setting splice5p to be false\n"));
    this->splice5p = false;
    this->splicetype5 = NO_SPLICE;
    this->ambig_prob_5 = 0.0;
    evalp = true;
  }

  ninserts = 0;
  for (i = 0; i < nsegments; i++) {
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&endpoint
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&nmismatches
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    *nmatches_trimmed += (Intlist_head(this->endpoints) - endpoint) - nmismatches;

    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));

    ninserts = Junction_ninserts(junction);
    if ((type = Junction_type(junction)) == SPLICE_JUNCTION) {
      this->splice5p = true;
      if (this->sensedir == SENSE_FORWARD) {
	if (this->plusp == true) {
	  this->splicetype5 = DONOR;
	  this->ambig_prob_5 = junction->donor_prob;
	} else {
	  this->splicetype5 = ANTIACCEPTOR;
	  this->ambig_prob_5 = junction->acceptor_prob;
	}
      } else {
	if (this->plusp == true) {
	  this->splicetype5 = ACCEPTOR;
	  this->ambig_prob_5 = junction->acceptor_prob;
	} else {
	  this->splicetype5 = ANTIDONOR;
	  this->ambig_prob_5 = junction->donor_prob;
	}
      }

    } else if (type == INS_JUNCTION) {
      Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) + Junction_ninserts(junction));
      if (Intlist_head(this->nmismatches) == 0) {
	/* Keep */
      } else {
	Intlist_head_set(this->nmismatches,-1);
      }

      debug(printf("(2) Setting splice5p to be false\n"));
      this->splice5p = false;
      this->splicetype5 = NO_SPLICE;
      this->ambig_prob_5 = 0.0;

    } else {
      debug(printf("(3) Setting splice5p to be false\n"));
      this->splice5p = false;
      this->splicetype5 = NO_SPLICE;
      this->ambig_prob_5 = 0.0;
    }

    Junction_free(&junction,pathpool);
    evalp = true;		/* Possibly a way to avoid calling Path_eval_nmatches, but need to handle ambig_prob_5 */
  }

  if (remainder == 0) {
    /* Keep endpoint and nmismatches */
  } else {
    remainder -= ninserts;	/* Undo advance of head endpoint */
    debug(printf("trim_qstart remainder is %d\n",remainder));
    debug(printf("(4) Setting splice5p to be false, as a result of trimming\n"));
    this->splice5p = false;
    this->splicetype5 = NO_SPLICE;
    this->ambig_prob_5 = 0.0;

    debug(printf("trim_qstart changing first endpoint from %d to %d\n",
		 Intlist_head(this->endpoints),Intlist_head(this->endpoints) + remainder));
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) + remainder);
    if (Intlist_head(this->nmismatches) == 0) {
      *nmatches_trimmed += remainder;

    } else {
      Intlist_head_set(this->nmismatches,-1);
      Intlist_head_set(this->ref_nmismatches,-1);
      evalp = true;
    }
    debug(printf("trim_qstart endpoints is now %s\n",Intlist_to_string(this->endpoints)));
  }
    
#if 0
  /* See if end segment is too short near the chromosome bound */
  qstart = Intlist_head(this->endpoints);
  qend = Intlist_second_value(this->endpoints);
  if (qend - qstart < MIN_END && this->junctions != NULL) {
    /* If so, then trim one more segment */
    debug(printf("Trimming one short segment at qstart\n"));
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
      
    /* Re-align to the chromosome start */
    Intlist_head_set(this->endpoints,qstart);
    Intlist_head_set(this->nmismatches,-1);
    Intlist_head_set(this->ref_nmismatches,-1);
  }
#endif

  if (nsegments > 0) {
    /* Alignment no longer follows the transcripts */
    /* TODO: Handle the case where nsegments is zero, but remainder is non-zero */
    Transcript_list_gc(&this->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&this->invalid_transcripts,listpool,transcriptpool);
  }
      
  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);

  return evalp;
}


/* Previously used when using Path_common_structure_trim_p in Path_consolidate */
/* Now called by Path_qstart_resolve */
void
Path_trim_qstart_n (int noutside, T this,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  int nsegments = 0, remainder, nmatches_trimmed;
  Intlist_T p;
  List_T j;
  int overall_qstart, last_ninserts, ninserts;
  int ignore_found_score = this->querylength;
  int ignore_end_score = this->querylength;

  
  debug(printf("Entering Path_trim_qstart_n with noutside %d against %s\n",
	       noutside,Intlist_to_string(this->endpoints)));
  debug(Path_print(this));

  p = this->endpoints;
  overall_qstart = Intlist_head(p);

  j = this->junctions;
  last_ninserts = 0;
  if (j == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts((Junction_T) List_head(j));
  }
  debug(printf("Comparing endpoint %d - ninserts %d against overall_qstart %d + noutside %d\n",
	       Intlist_second_value(p),ninserts,overall_qstart,noutside));
  while (Intlist_next(p) != NULL && Intlist_second_value(p) - ninserts <= overall_qstart + noutside) {
    nsegments++;
    p = Intlist_next(p);
    last_ninserts = 0;
    if (j == NULL) {
      ninserts = 0;
    } else if ((j = List_next(j)) == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }
  
  debug(printf("Now endpoints is %s\n",Intlist_to_string(p)));
  if (Intlist_next(p) == NULL) {
    debug(printf("Trimming would remove all segments, so unable to comply\n"));

  } else {
    /* Previously checked qstart = Intlist_head(p)) - overall_qstart >= noutside */
    remainder = noutside - (Intlist_head(p) - overall_qstart) + last_ninserts;
    debug(printf("Computed %d segments, remainder %d (using last_ninserts %d)\n",
		 nsegments,remainder,last_ninserts));
    
    if (trim_qstart(&nmatches_trimmed,nsegments,remainder,this,
		    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
      /* Call even if nsegments == 0 and remainder == 0 */
      this->nmatches -= nmatches_trimmed;
    } else {
      this->nmatches = -1;
      Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
    }
    
    /* Used to always call Path_eval_nmatches */
    debug(printf("After trimming: "));
    debug(Path_print(this));
    debug(printf("\n"));
  }

  return;
}


#if 0
/* Will need to call Path_eval_nmatches again after this, so don't
   worry about trimmedp */
/* Replaced by Path_trim_chrbounds */
void
Path_trim_qstart_trimbounds (T this,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			     Univcoord_T trimbounds) {
  int noutside, nsegments;
  int remainder;
  bool trimmedp;

  debug(printf("Path_trim_qstart_trimbounds with trimbounds %u\n",trimbounds));
  debug(Path_print(this));

  if ((noutside = count_outside_qstart_byposition(&nsegments,&remainder,
						  this->endpoints,this->univdiagonals,
						  trimbounds,this->querylength)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    fprintf(stderr,"Attempt to trim entire path\n");
    abort();

  } else if (trim_qstart(&nmatches_trimmed,nsegments,remainder,this,
			 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
    this->nmatches -= nmatches_trimmed;
  } else {
    this->nmatches = -1;
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
  }

  return;
}
#endif


/* Returns false is result is not valid, usually because of the mate
   being inside of the path */
bool
Path_trim_qstart_trimdiag (T this,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
			   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Univcoord_T trimdiag) {
  int nsegments;
  int ignore_found_score = this->querylength;
  int ignore_end_score = this->querylength;
  int nmatches_trimmed;

  debug(printf("Entering Path_trim_qstart_trimdiag with trimdiag %u\n",trimdiag));
  debug(Path_print(this));

  if ((nsegments = count_outside_qstart_bydiagonal(this->univdiagonals,trimdiag)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));
    return true;

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    /* fprintf(stderr,"Attempt to trim entire path\n"); */
    return false;

  } else if (trim_qstart(&nmatches_trimmed,nsegments,/*remainder*/0,this,
			 intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
    /* Call even if nsegments == 0 and remainder == 0 */
    this->nmatches -= nmatches_trimmed;
#ifdef CHECK_ASSERTIONS
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
#endif
    return true;

  } else {
    this->nmatches = -1;
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
    return true;
  }
}


/* Needed by divide_path */
/* Expecting path to be reversed */
static int
count_outside_qend_byposition (int *nsegments, int *remainder,
			       Intlist_T endpoints, Univcoordlist_T univdiagonals,
			       Univcoord_T trimbounds, int querylength) {
  int noutside = 0;
  Univcoord_T left, univdiagonal;
  int qstart, qend;

  *nsegments = 0;
  *remainder = 0;

  debug(printf("Entered count_outside_qend_byposition\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);
    left = univdiagonal - querylength;
    qend = Intlist_head(endpoints);
    qstart = Intlist_second_value(endpoints);
    debug(printf("univdiagonal %u, left %u (%u..%u) vs trimbounds %u\n",
		 univdiagonal,left,left+qstart,left+qend,trimbounds));

    if (left + qend < trimbounds) {
      *remainder = 0;
      return noutside;

    } else if (left + qstart < trimbounds) {
      *remainder = (int) ((left + qend) - trimbounds);
      return noutside + (*remainder);

    } else {
      *nsegments += 1;
      noutside += qend - qstart;
      endpoints = Intlist_next(endpoints);
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return noutside;
}


static int
count_outside_qend_bydiagonal (Univcoordlist_T univdiagonals, int querylength,
			       Univcoord_T trimdiag) {
  int nsegments = 0;
  Univcoord_T univdiagonal;

  debug(printf("Entered count_outside_qend_bydiagonal\n"));
  while (univdiagonals != NULL) {
    univdiagonal = Univcoordlist_head(univdiagonals);

    /* Allow the same univdiagonal */
    if (univdiagonal - querylength <= trimdiag) {
      return nsegments;

    } else {
      nsegments += 1;
      univdiagonals = Univcoordlist_next(univdiagonals);
    }
  }

  return nsegments;
}


/* nmatches_trimmed is valid only if return value (evalp) is false */
/* Assumes path has been reversed */
static bool
trim_qend (int *nmatches_trimmed, int nsegments, int remainder, T this,
	   Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  bool evalp = false;
  int i;
  /* int qstart, qend; */
  int ignore_int = 0, nmismatches, endpoint;
  Univcoord_T ignore_univcoord;
  Junction_T junction;

  debug(printf("Path needs trimming at qend (%d segments, %d remainder)\n",
	       nsegments,remainder));
  debug(Path_print(this));

  *nmatches_trimmed = 0;

  if (nsegments > 0 && this->qend_alts != NULL) {
    Altsplice_free(&this->qend_alts,pathpool);
    this->qend_alts = (Altsplice_T) NULL;
    debug(printf("(1) Setting splice3p to be false\n"));
    this->splice3p = false;
    this->splicetype3 = NO_SPLICE;
    this->ambig_prob_3 = 0.0;
    evalp = true;
  }

  for (i = 0; i < nsegments; i++) {
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&endpoint
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&nmismatches
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    *nmatches_trimmed += (endpoint - Intlist_head(this->endpoints)) - nmismatches;

    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));

    if (Junction_type(junction) == SPLICE_JUNCTION) {
      this->splice3p = true;
      if (this->sensedir == SENSE_FORWARD) {
	if (this->plusp == true) {
	  this->splicetype3 = ACCEPTOR;
	  this->ambig_prob_3 = junction->acceptor_prob;
	} else {
	  this->splicetype3 = ANTIDONOR;
	  this->ambig_prob_3 = junction->donor_prob;
	}
      } else {
	if (this->plusp == true) {
	  this->splicetype3 = DONOR;
	  this->ambig_prob_3 = junction->donor_prob;
	} else {
	  this->splicetype3 = ANTIACCEPTOR;
	  this->ambig_prob_3 = junction->acceptor_prob;
	}
      }

    } else {
      debug(printf("(2) Setting splice3p to be false\n"));
      this->splice3p = false;
      this->splicetype3 = NO_SPLICE;
      this->ambig_prob_3 = 0.0;
    }

    Junction_free(&junction,pathpool);

    evalp = true;		/* Possibly a way to avoid calling Path_eval_nmatches, but need to handle ambig_prob_3 */
  }

  if (remainder == 0) {
    /* Keep endpoint and nmismatches */
  } else {
    debug(printf("trim_qend remainder is %d\n",remainder));
    debug(printf("(3) Setting splice3p to be false, as a result of trimming\n"));
    this->splice3p = false;
    this->splicetype3 = NO_SPLICE;
    this->ambig_prob_3 = 0.0;

    debug(printf("trim_qend changing first endpoint from %d to %d\n",
		 Intlist_head(this->endpoints),Intlist_head(this->endpoints) - remainder));
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) - remainder);
    if (Intlist_head(this->nmismatches) == 0) {
      *nmatches_trimmed += remainder;
    } else {
      Intlist_head_set(this->nmismatches,-1);
      Intlist_head_set(this->ref_nmismatches,-1);
      evalp = true;
    }
    debug(printf("trim_qend endpoints is now %s\n",Intlist_to_string(this->endpoints)));
  }
    
#if 0
  /* See if end segment is too short near the chromosome bound */
  qend = Intlist_head(this->endpoints);
  qstart = Intlist_second_value(this->endpoints);
  if (qend - qstart < MIN_END && this->junctions != NULL) {
    /* If so, then trim one more segment */
    debug(printf("Trimming one short segment at qend\n"));
    this->endpoints = Intlistpool_pop(this->endpoints,intlistpool,&ignore_int
				      intlistpool_trace(__FILE__,__LINE__));
    this->nmismatches = Intlistpool_pop(this->nmismatches,intlistpool,&ignore_int
					intlistpool_trace(__FILE__,__LINE__));
    this->ref_nmismatches = Intlistpool_pop(this->ref_nmismatches,intlistpool,&ignore_int
					    intlistpool_trace(__FILE__,__LINE__));
    this->univdiagonals = Univcoordlistpool_pop(this->univdiagonals,univcoordlistpool,&ignore_univcoord
						univcoordlistpool_trace(__FILE__,__LINE__));
    this->junctions = Listpool_pop(this->junctions,listpool,(void **) &junction
				   listpool_trace(__FILE__,__LINE__));
    Junction_free(&junction,pathpool);
      
    /* Re-align to the chromosome end */
    Intlist_head_set(this->endpoints,qend);
    Intlist_head_set(this->nmismatches,-1);
    Intlist_head_set(this->ref_nmismatches,-1);
  }
#endif

  if (nsegments > 0) {
    /* Alignment no longer follows the transcripts */
    /* TODO: Handle the case where nsegments is zero, but remainder is non-zero */
    Transcript_list_gc(&this->transcripts,listpool,transcriptpool);
    Transcript_list_gc(&this->invalid_transcripts,listpool,transcriptpool);
  }
    
  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);

  return evalp;
}


/* Previously used when using Path_common_structure_trim_p in Path_consolidate */
/* Now called by Path_qend_resolve */
void
Path_trim_qend_n (int noutside, T this,
		  Compress_T query_compress_fwd, Compress_T query_compress_rev,
		  Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
  int nsegments = 0, remainder, nmatches_trimmed;
  Intlist_T p;
  List_T j;
  int overall_qend, last_ninserts, ninserts;
  int ignore_found_score = this->querylength;
  int ignore_end_score = this->querylength;

  debug(printf("Entering Path_trim_qend_n with noutside %d against %s\n",
	       noutside,Intlist_to_string(this->endpoints)));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  p = this->endpoints;
  overall_qend = Intlist_head(p);

  j = this->junctions;
  last_ninserts = 0;
  if (j == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts((Junction_T) List_head(j));
  }
  debug(printf("Comparing endpoint %d + ninserts %d against overall_qend %d - noutside %d\n",
	       Intlist_second_value(p),ninserts,overall_qend,noutside));
  while (Intlist_next(p) != NULL && Intlist_second_value(p) + ninserts >= overall_qend - noutside) {
    nsegments++;
    p = Intlist_next(p);
    last_ninserts = ninserts;
    if (j == NULL) {
      ninserts = 0;
    } else if ((j = List_next(j)) == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }

  debug(printf("Now endpoints is %s\n",Intlist_to_string(p)));
  if (Intlist_next(p) == NULL) {
    debug(printf("Trimming would remove all segments, so unable to comply\n"));
    Path_reverse(this,/*expect_fwd_p*/true);

  } else {
    /* Previously checked overall_qend - (qend = Intlist_head(p)) >= noutside */
    remainder = noutside - (overall_qend - Intlist_head(p)) + last_ninserts;
    debug(printf("Computed %d segments, remainder %d (using last_ninserts %d)\n",
		 nsegments,remainder,last_ninserts));
    
    if (trim_qend(&nmatches_trimmed,nsegments,remainder,this,
		  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
      /* Call even if nsegments == 0 and remainder == 0 */
      this->nmatches -= nmatches_trimmed;
      Path_reverse(this,/*expect_fwd_p*/true);
#ifdef CHECK_ASSERTIONS
      Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
#endif
      
    } else {
      Path_reverse(this,/*expect_fwd_p*/true);
      this->nmatches = -1;
      Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
    }

    debug(printf("After trimming: "));
    debug(Path_print(this));
    debug(printf("\n"));
  }

  return;
}



#if 0
/* Replaced by Path_trim_chrbounds */
/* Will need to call Path_eval_nmatches again after this, so don't
   worry about trimmedp */
void
Path_trim_qend_trimbounds (T this, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Univcoord_T trimbounds) {
  int noutside, nsegments;
  int remainder;
  bool trimmedp;

  debug(printf("Path_trim_qend_trimbounds with trimbounds %u\n",trimbounds));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  if ((noutside = count_outside_qend_byposition(&nsegments,&remainder,
						this->endpoints,this->univdiagonals,
						trimbounds,this->querylength)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    fprintf(stderr,"Attempt to trim entire path\n");
    abort();

  } else if (trim_qend(&nmatches_trimmed,nsegments,remainder,this,
		       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
    /* Call even if nsegments == 0 and remainder == 0 */
    this->nmatches -= nmatches_trimmed;
#ifdef CHECK_ASSERTIONS
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
#endif

  } else {
    this->nmatches = -1;
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
  }

  Path_reverse(this,/*expect_fwd_p*/true);

  return;
}
#endif


/* Returns false is result is not valid, usually because of the mate
   being inside of the path */
bool
Path_trim_qend_trimdiag (T this,
			 Compress_T query_compress_fwd, Compress_T query_compress_rev,
			 Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			 Listpool_T listpool, Pathpool_T pathpool,
			 Transcriptpool_T transcriptpool, int querylength, Univcoord_T trimdiag) {
  int nsegments;
  int ignore_found_score = this->querylength;
  int ignore_end_score = this->querylength;
  int nmatches_trimmed;

  debug(printf("Entering Path_trim_qend_trimdiag with trimdiag %u\n",trimdiag));
  debug(Path_print(this));

  Path_reverse(this,/*expect_fwd_p*/false);

  if ((nsegments = count_outside_qend_bydiagonal(this->univdiagonals,querylength,trimdiag)) == 0) {
    /* Skip */
    debug(printf("No trimming needed\n"));
    Path_reverse(this,/*expect_fwd_p*/true);
    return true;

  } else if (nsegments >= Univcoordlist_length(this->univdiagonals)) {
    /* fprintf(stderr,"Attempt to trim entire path\n"); */
    Path_reverse(this,/*expect_fwd_p*/true);
    return false;

  } else if (trim_qend(&nmatches_trimmed,nsegments,/*remainder*/0,this,
		       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == false) {
    /* Call even if nsegments == 0 and remainder == 0 */
    this->nmatches -= nmatches_trimmed;
    Path_reverse(this,/*expect_fwd_p*/true);
#ifdef CHECK_ASSERTIONS
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
#endif
    return true;

  } else {
    this->nmatches = -1;
    Path_reverse(this,/*expect_fwd_p*/true);
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
    return true;
  }
}


bool
Path_trim_chrbounds (int *found_score, int *end_score, T this,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		     Listpool_T listpool, Pathpool_T pathpool,
		     Transcriptpool_T transcriptpool) {
  int nsegments, remainder, qstart, qend, ninserts;
  Univcoord_T univdiagonal;
  Univcoordlist_T p;
  Intlist_T q;
  List_T j;
  int ignore_found_score = this->querylength;
  int ignore_end_score = this->querylength;
  bool qstart_evalp = false, qend_evalp = false;
  int qstart_nmatches_trimmed = 0, qend_nmatches_trimmed = 0;
  
  
  debug(printf("Entering Path_trim_chrbounds with chroffset %u and chrhigh %u\n",
	       this->chroffset,this->chrhigh));
  debug(Path_print(this));

  /* qstart */
  if (this->qstart_alts != NULL &&
      Altsplice_trim_qstart_chrbounds(this->qstart_alts,this->chroffset,this->querylength) == false) {
    Altsplice_free(&this->qstart_alts,pathpool);
    this->qstart_alts = (Altsplice_T) NULL;
    this->splice5p = false;
    this->splicetype5 = NO_SPLICE;
    this->ambig_prob_5 = 0.0;
    qstart_evalp = true;
  }

  nsegments = 0;
  p = this->univdiagonals;
  q = this->endpoints;
  while (p != NULL &&
	 Univcoordlist_head(p) + /*qend*/Intlist_second_value(q) < this->chroffset + this->querylength) {
    nsegments++;
    p = Univcoordlist_next(p);
    q = Intlist_next(q);
  }
  
  if (p == NULL) {
    /* All segments are before chromosomal start */
    debug(printf("All segments are before chromosomal start, so returning false\n"));
    return false;

  } else {
    univdiagonal = Univcoordlist_head(p);
    qstart = Intlist_head(q);
    debug(printf("Now univdiagonals is %s, qstart %d\n",Univcoordlist_to_string(p),qstart));
    debug(printf("Comparing position %u = left %u + qstart %d, with chroffset %u\n",
		 univdiagonal - this->querylength + qstart,univdiagonal - this->querylength,qstart,
		 this->chroffset));
    
    /* remainder = chroffset - (univdiagonal - querylength + endpoint) */
    if (univdiagonal + qstart >= this->chroffset + this->querylength) {
      remainder = 0;
      debug(printf("Computed %d segments, remainder 0\n",nsegments));
    } else {
      remainder = this->chroffset + this->querylength - univdiagonal - qstart;
      debug(printf("Computed %d segments, remainder %d = chroffset %u - (left %u + qstart %d)\n",
		   nsegments,remainder,this->chroffset,univdiagonal - this->querylength,qstart));
    }
    /* Call even if nsegments == 0 and remainder == 0 */
    if (trim_qstart(&qstart_nmatches_trimmed,nsegments,remainder,this,
		    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == true) {
      qstart_evalp = true;
    }
  }

  /* qend */
  if (this->qend_alts != NULL &&
      Altsplice_trim_qend_chrbounds(this->qend_alts,this->chrhigh) == false) {
    debug(printf("Removing qend_alts\n"));
    Altsplice_free(&this->qend_alts,pathpool);
    this->qend_alts = (Altsplice_T) NULL;
    this->splice3p = false;
    this->splicetype3 = NO_SPLICE;
    this->ambig_prob_3 = 0.0;
    qend_evalp = true;
  }

  if (Path_endpoints_acceptable_p(this->endpoints,this->junctions) == false) {
    debug(printf("Path endpoints not acceptable\n"));
    return false;
  } else {
    Path_reverse(this,/*expect_fwd_p*/false);
  }

  nsegments = 0;
  ninserts = 0;
  p = this->univdiagonals;
  q = this->endpoints;
  j = this->junctions;
  while (p != NULL &&
	 Univcoordlist_head(p) + /*qstart*/Intlist_second_value(q) + ninserts >= this->chrhigh + this->querylength) {
    nsegments++;
    p = Univcoordlist_next(p);
    q = Intlist_next(q);
    if (j == NULL) {
      ninserts = 0;
    } else {
      ninserts = Junction_ninserts((Junction_T) List_head(j));
      j = List_next(j);
    }
  }

  if (p == NULL) {
    /* All segments are after chromosomal end */
    debug(printf("All segments are after chromosomal end, so returning false\n"));
    return false;

  } else {
    univdiagonal = Univcoordlist_head(p);
    qend = Intlist_head(q);
    debug(printf("Now univdiagonals is %s, qend %d\n",Univcoordlist_to_string(p),qend));
    if (univdiagonal + qend < this->chrhigh + this->querylength) {
      remainder = 0;
      debug(printf("Computed %d segments, remainder 0\n",nsegments));
    } else {
      remainder = univdiagonal + qend - this->querylength - this->chrhigh;
      debug(printf("Computed %d segments, remainder %d = (left %u + qend %d) - chrhigh %u\n",
		   nsegments,remainder,univdiagonal - this->querylength,qend,this->chrhigh));
    }

    /* Call even if nsegments == 0 and remainder == 0 */
    if (trim_qend(&qend_nmatches_trimmed,nsegments,remainder,this,
		  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool) == true) {
      qend_evalp = true;
    }
  }

  /* Want nocheck in case endpoints are not acceptable */
  Path_reverse_nocheck(this);

  if (Path_endpoints_acceptable_p(this->endpoints,this->junctions) == false) {
    debug(printf("Path endpoints not acceptable\n"));
    return false;
  }

  if (this->nmatches < 0) {
    /* Path was changed by caller, e.g., Pathpair_resolve */
    Path_eval_nmatches(&(*found_score),&(*end_score),this,query_compress_fwd,query_compress_rev);

  } else if (qstart_evalp == false && qend_evalp == false) {
    debug(printf("Starting with nmatches %d, but subtracting %d from start and %d from end\n",
		 this->nmatches,qstart_nmatches_trimmed,qend_nmatches_trimmed));
    this->nmatches -= qstart_nmatches_trimmed;
    this->nmatches -= qend_nmatches_trimmed;
#ifdef CHECK_ASSERTIONS
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);
#endif
  } else {
    this->nmatches = -1;
    Path_eval_nmatches(&(*found_score),&(*end_score),this,query_compress_fwd,query_compress_rev);
  }

  debug(printf("After trimming: "));
  debug(Path_print(this));
  debug(printf("\n"));

  return true;
}



static int
count_total (Intlist_T endpoints, Univcoordlist_T univdiagonals) {
  int ntotal = 0;
  int qstart, qend;

  while (univdiagonals != NULL) {
    qstart = Intlist_head(endpoints);
    qend = Intlist_second_value(endpoints);

    ntotal += qend - qstart;
    endpoints = Intlist_next(endpoints);
    univdiagonals = Univcoordlist_next(univdiagonals);
  }

  return ntotal;
}


static int
divide_path (int *nhigh, int *qend_nsegments, int *qend_remainder,
	     Univcoord_T *linear_chrhigh, Univcoord_T *linear_chrlength, T this) {
  int ntotal;

  /* Determine the point of division */
  *linear_chrlength = (this->chrhigh - this->chroffset)/2;
  *linear_chrhigh = this->chroffset + (*linear_chrlength);
  debug2(printf("divide_path\n"));
  debug2(Path_print(this));
  debug2(printf("chrhigh %u, chroffset %u, chrlength %u\n",
		this->chrhigh,this->chroffset,*linear_chrlength));

  ntotal = count_total(this->endpoints,this->univdiagonals);

  Path_reverse(this,/*expect_fwd_p*/false);
  *nhigh = count_outside_qend_byposition(&(*qend_nsegments),&(*qend_remainder),
					 this->endpoints,this->univdiagonals,
					 /*trimbounds*/*linear_chrhigh,this->querylength);
  Path_reverse(this,/*expect_fwd_p*/true);
  debug2(printf("nlow %d, nhigh %d\n",ntotal - (*nhigh),*nhigh));

  return (ntotal - (*nhigh));
}


/* Moves path from higher coordinates to lower ones */
/* Need to check that this alignment has no splices (although it could have indels) */
static void
unalias_path (T this, Univcoord_T linear_chrlength, Compress_T query_compress,
	      Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  Univcoord_T circular_origin, new_univdiagonal, univdiagonal;
  Univcoordlist_T u;
  Altsplice_T altsplice;
  Intlist_T q, r, s;
  List_T j;
  int qstart, qend, last_ninserts;
  int nmismatches, ref_nmismatches, ignore_found_score, ignore_end_score;
  bool cross_origin_p = false;
  int i;

  debug3(printf("Entered unalias path with "));
  debug3(Path_print(this));

  if (this->main_univdiagonal < linear_chrlength) {
    this->main_univdiagonal = 0;
  } else {
    this->main_univdiagonal -= linear_chrlength;
  }
  circular_origin = this->chroffset + linear_chrlength;

  last_ninserts = 0;
  for (u = this->univdiagonals, q = this->endpoints, s = this->nmismatches, r = this->ref_nmismatches, j = this->junctions;
       u != NULL; u = Univcoordlist_next(u), q = Intlist_next(q), s = Intlist_next(s), r = Intlist_next(r), j = List_next(j)) {
    univdiagonal = Univcoordlist_head(u);
    new_univdiagonal = univdiagonal - linear_chrlength;
    Univcoordlist_head_set(u,new_univdiagonal);

    qstart = Intlist_head(q) + last_ninserts;
    qend = Intlist_second_value(q);

    if (new_univdiagonal + qstart < circular_origin + this->querylength &&
	new_univdiagonal + qend >= circular_origin + this->querylength) {
      /* If alignment crosses the circular origin, need to recompute
	 nmismatches and genomic_diff because genomic sequence around
	 the origin is different from the previous and subsequent
	 chromosome */
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  /*univdiagonal*/new_univdiagonal,
							  this->querylength,/*pos5*/qstart,/*pos3*/qend,
							  this->plusp,/*genestrand*/0);
      debug3(printf("Recomputing nmismatches as %d for univdiagonal %u, %d..%d\n",
		    nmismatches,new_univdiagonal,qstart,qend));

      Intlist_head_set(s,nmismatches);
      Intlist_head_set(r,ref_nmismatches);
      cross_origin_p = true;
    }

    if (j != NULL) {
      last_ninserts = Junction_ninserts((Junction_T) List_head(j));
    }
  }
  
  if (cross_origin_p == true) {
    this->genomic_diff = (char *) NULL;

    /* Need to re-evaluate path, since matches and mismatches might have gone past chrhigh */
    this->nmatches = -1;
    ignore_found_score = ignore_end_score = this->querylength;
    Path_eval_nmatches(&ignore_found_score,&ignore_end_score,this,query_compress_fwd,query_compress_rev);

    /* Pathpool_free_string(&this->genomic_diff); */
  }

  /* Should not have splicing in a circular chromosome */
  /* It is possible that some qstart alts_coords are above linear_chrhigh and some below */
  /* However, if we don't unalias all of them, we could get a negative splice distance */
  if ((altsplice = this->qstart_alts) != NULL) {
    altsplice->medial_univdiagonal -= linear_chrlength;
    for (i = 0; i < altsplice->nunivdiagonals; i++) {
      /* if (altsplice->univdiagonals[i] >= linear_chrhigh) { */
      altsplice->univdiagonals[i] -= linear_chrlength;
      /* } */
    }
  }

  if ((altsplice = this->qend_alts) != NULL) {
    altsplice->medial_univdiagonal -= linear_chrlength;
    for (i = 0; i < altsplice->nunivdiagonals; i++) {
      altsplice->univdiagonals[i] -= linear_chrlength;
    }
  }

  debug3(printf("Exited unalias path with "));
  debug3(Path_print(this));

  return;
}


/* Returns true if the entire path is in the upper half, meaning that
   it must have an alias in the lower half */
void
Path_trim_circular_unalias (T this, Compress_T query_compress,
			    Compress_T query_compress_fwd, Compress_T query_compress_rev) {
  int nlow, nhigh, qend_nsegments, qend_remainder;
  Univcoord_T linear_chrhigh, linear_chrlength;

  if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			  &linear_chrhigh,&linear_chrlength,this)) > 0) {
    /* Do not alter the path */
  } else {
    unalias_path(this,linear_chrlength,query_compress,query_compress_fwd,query_compress_rev);
  }

  return;
}


void
Path_trim_circular_unalias_pair (T path5, T path3, Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
				 Compress_T query3_compress_fwd, Compress_T query3_compress_rev) {
  int nlow, nhigh, qend_nsegments, qend_remainder;
  Univcoord_T linear_chrhigh, linear_chrlength;

  if (path5->plusp == true) {
    if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			    &linear_chrhigh,&linear_chrlength,/*pathL*/path5)) > 0) {
      /* Do not alter either path */
    } else {
      /* Alter both paths */
      unalias_path(path5,linear_chrlength,query5_compress_fwd,
		   query5_compress_fwd,query5_compress_rev);
      unalias_path(path3,linear_chrlength,query3_compress_fwd,
		   query3_compress_rev,query3_compress_rev);
    }
  } else {
    if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			    &linear_chrhigh,&linear_chrlength,/*pathL*/path3)) > 0) {
      /* Do not alter either path */
    } else {
      /* Alter both paths */
      unalias_path(path5,linear_chrlength,query5_compress_rev,
		   query5_compress_fwd,query5_compress_rev);
      unalias_path(path3,linear_chrlength,query3_compress_rev,
		   query3_compress_fwd,query3_compress_rev);
    }
  }

  return;
}


/* For circular alignments: (1) if the entire path is in the lower
   path, we keep it; (2) if the entire path is in the upper half,
   then we unalias it.  Otherwise, if it straddles the circular
   origin, we do not touch the main endpoints and univdiagonals, but
   create circlow and circhigh endpoints */

void
Path_trim_circular (T this, Compress_T query_compress,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool) {
  int nlow, nhigh, qend_nsegments, qend_remainder, i;
  Univcoord_T linear_chrhigh, linear_chrlength, univdiagonal;
  int nmismatches, ref_nmismatches;

  Intlist_T high_endpoints = NULL, high_nmismatches = NULL, high_ref_nmismatches = NULL;
  Univcoordlist_T high_univdiagonals = NULL;
  List_T high_junctions = NULL;
  Junction_T junction;
  Altsplice_T altsplice;

  debug2(printf("Entering Path_trim_circular\n"));

  if ((nlow = divide_path(&nhigh,&qend_nsegments,&qend_remainder,
			  &linear_chrhigh,&linear_chrlength,this)) == 0) {
    /* Entire alignment is in the top half of the double chromosome, so unalias all univdiagonals */
    debug2(printf("Entire alignment is in the top half\n"));
    Path_trim_circular_unalias(this,query_compress,
			       query_compress_fwd,query_compress_rev);

  } else if (nhigh == 0) {
    /* Entire alignment is in the bottom half of the double chromosome, so nothing to do */
    debug2(printf("Entire alignment is in the lower half\n"));

  } else {
    /* Divide up the components */
    debug2(printf("Dividing alignment: qend_nsegments %d, qend_remainder %d\n",
		  qend_nsegments,qend_remainder));
    Path_reverse(this,/*expect_fwd_p*/false);

    for (i = 0; i < qend_nsegments; i++) {
      debug2(printf("Pushing %d onto high_endpoints\n",Intlist_head(this->endpoints)));
      high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints)
					intlistpool_trace(__FILE__,__LINE__));
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,Intlist_head(this->nmismatches)
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,Intlist_head(this->ref_nmismatches)
					      intlistpool_trace(__FILE__,__LINE__));
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,
						  Univcoordlist_head(this->univdiagonals) - linear_chrlength
						  univcoordlistpool_trace(__FILE__,__LINE__));

      junction = (Junction_T) List_head(this->junctions);
      high_junctions = Listpool_push(high_junctions,listpool,(void *) junction
				     listpool_trace(__FILE__,__LINE__));

      this->endpoints = Intlist_next(this->endpoints);
      this->nmismatches = Intlist_next(this->nmismatches);
      this->ref_nmismatches = Intlist_next(this->ref_nmismatches);
      this->univdiagonals = Univcoordlist_next(this->univdiagonals);
      this->junctions = List_next(this->junctions);
    }
			    
    debug2(printf("Pushing %d onto high_endpoints\n",Intlist_head(this->endpoints)));
    high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints)
				      intlistpool_trace(__FILE__,__LINE__));

    high_endpoints = Intlistpool_push(high_endpoints,intlistpool,Intlist_head(this->endpoints) - qend_remainder
				      intlistpool_trace(__FILE__,__LINE__));
    Intlist_head_set(this->endpoints,Intlist_head(this->endpoints) - qend_remainder);

    univdiagonal = Univcoordlist_head(this->univdiagonals);
    if (univdiagonal < this->chroffset + linear_chrlength) {
      /* Can happen if there is splicing across the circular origin */
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,univdiagonal
						  univcoordlistpool_trace(__FILE__,__LINE__));
    } else {
      high_univdiagonals = Univcoordlistpool_push(high_univdiagonals,univcoordlistpool,
						  univdiagonal - linear_chrlength
						  univcoordlistpool_trace(__FILE__,__LINE__));
    }

    if (Intlist_head(this->nmismatches) == 0 && Intlist_head(this->ref_nmismatches) == 0) {
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,0
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,0
					      intlistpool_trace(__FILE__,__LINE__));
      
    } else {
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,this->querylength,
							  /*pos5*/Intlist_second_value(this->endpoints),
							  /*pos3*/Intlist_head(this->endpoints),
							  this->plusp,this->genestrand);
      Intlist_head_set(this->nmismatches,nmismatches);
      Intlist_head_set(this->ref_nmismatches,ref_nmismatches);
      
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  univdiagonal,this->querylength,
							  /*pos5*/Intlist_head(high_endpoints),
							  /*pos3*/Intlist_second_value(high_endpoints),
							  this->plusp,this->genestrand);
      high_nmismatches = Intlistpool_push(high_nmismatches,intlistpool,nmismatches
					  intlistpool_trace(__FILE__,__LINE__));
      high_ref_nmismatches = Intlistpool_push(high_ref_nmismatches,intlistpool,nmismatches
					      intlistpool_trace(__FILE__,__LINE__));
    }
    

    /* Unalias qend_alts, which are above chrlength */
    if ((altsplice = this->qend_alts) != NULL) {
      altsplice->medial_univdiagonal -= linear_chrlength;
      for (i = 0; i < altsplice->nunivdiagonals; i++) {
	altsplice->univdiagonals[i] -= linear_chrlength;
      }
    }

    /* Restore low path */
    Path_reverse(this,/*expect_fwd_p*/true);

    if (nlow >= nhigh) {
      /* Make high path the circular part */
      this->circular_high_p = true;
      this->circular_endpoints = high_endpoints;
      this->circular_nmismatches = high_nmismatches;
      this->circular_ref_nmismatches = high_ref_nmismatches;
      this->circular_univdiagonals = high_univdiagonals;
      this->circular_junctions = high_junctions;

    } else {
      /* Make low path the circular part */
      this->circular_high_p = false;
      this->circular_endpoints = this->endpoints;
      this->circular_nmismatches = this->nmismatches;
      this->circular_ref_nmismatches = this->ref_nmismatches;
      this->circular_univdiagonals = this->univdiagonals;
      this->circular_junctions = this->junctions;
      
      this->endpoints = high_endpoints;
      this->nmismatches = high_nmismatches;
      this->ref_nmismatches = high_ref_nmismatches;
      this->univdiagonals = high_univdiagonals;
      this->junctions = high_junctions;
    }

    /* TODO: Fix transcripts (or print them only with the primary entry) */
  }

  debug2(printf("Result of dividing path\n"));
  debug2(Path_print(this));

  return;
}


void
Path_trim_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in) {
  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;
  return;
}
