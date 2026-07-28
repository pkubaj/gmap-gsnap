static char rcsid[] = "$Id: 1dd83c68fef05cae9d673fae0bca77867d24bd17 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "transcript.h"

#include "sedgesort.h"
#include "genomebits_count.h"


static bool *circularp;
static bool *altlocp;


#define ENDTRIM_ALLOWED 4


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Path_print */
/* Also, need to define DEBUG1 in junction.c */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Path_endpoints_acceptable_p */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Path_filter */
#ifdef DEBUG4
#define debug4(x) x
#else
#define debug4(x)
#endif

/* Getting unique genomic coords */
#ifdef DEBUG6
#define debug6(x) x
#else
#define debug6(x)
#endif

/* Setdiff univdiagonals auxinfo */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* Merge univdiagonals auxinfo */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif


#define T Path_T


int
Path_nbadsplices (T this) {
  int nbadsplices = 0;
  List_T p;
  Junction_T junction;

  for (p = this->junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      if (Junction_donor_prob(junction) * Junction_acceptor_prob(junction) < 0.6) {
	nbadsplices++;
      }
    }
  }

  return nbadsplices;
}


double
Path_avg_junction_splice_prob (T this) {
  int nsplices = 0;
  List_T p;
  Junction_T junction;

  for (p = this->junctions; p != NULL; p = List_next(p)) {
    junction = (Junction_T) List_head(p);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      nsplices++;
    }
  }

  if (nsplices == 0) {
    return 0.0;
  } else {
    return this->junction_splice_prob/(double) nsplices;
  }
}

double
Path_end_splice_prob (T this) {
  return this->ambig_prob_5 + this->ambig_prob_3;
}



/* TODO: Consider whether we need to check ambig_prob_5 and
   ambig_prob_3, since they are supposed to be incorporated into
   splice_prob in Path_nmatches */
int
Path_effective_sensedir (T this) {
  List_T p;

  if (this->qstart_alts != NULL) {
    return this->sensedir;
  } else if (this->qend_alts != NULL) {
    return this->sensedir;
  } else if (this->splicetype5 != NO_SPLICE) {
    return this->sensedir;
  } else if (this->splicetype3 != NO_SPLICE) {
    return this->sensedir;
  } else {
    for (p = this->junctions; p != NULL; p = List_next(p)) {
      if (Junction_type((Junction_T) List_head(p)) == SPLICE_JUNCTION) {
	return this->sensedir;
      }
    }
    return SENSE_NULL;
  }
}


int
Path_coverage (T this) {
  int pos5, pos3, fusion_coverage;

  if (this->qstart_alts != NULL) {
    pos5 = 0;
  } else {
    pos5 = Intlist_head(this->endpoints);
  }
  if (this->qend_alts != NULL) {
    pos3 = this->querylength;
  } else {
    pos3 = Intlist_last_value(this->endpoints);
  }
  
  if (this->fusion_querystart_junction != NULL) {
    fusion_coverage = Intlist_last_value(this->fusion_endpoints) - Intlist_head(this->fusion_endpoints);
  } else if (this->fusion_queryend_junction != NULL) {
    fusion_coverage = Intlist_last_value(this->fusion_endpoints) - Intlist_head(this->fusion_endpoints);
  } else {
    fusion_coverage = 0;
  }

  return pos3 - pos5 + fusion_coverage;
}


bool
Path_softclippedp (T this) {
  if (Intlist_head(this->endpoints) > 0) {
    return true;
  } else if (this->querylength - Intlist_last_value(this->endpoints) > 0) {
    return true;
  } else {
    return false;
  }
}


bool
Path_resolved_qstart_p (T this) {
  /* Presence of qstart_alts means not resolved, and implies endpoint is not 0 */
  if (Intlist_head(this->endpoints) == 0) {
    return true;
  } else if (this->plusp == true && this->fusion_querystart_junction != NULL) {
    return true;
  } else if (this->plusp == false && this->fusion_queryend_junction != NULL) {
    return true;
  } else {
    return false;
  }
}


bool
Path_resolved_qend_p (T this) {
  /* Presence of qend_alts means not resolved, and implies endpoint is not querylength */
  if (Intlist_last_value(this->endpoints) == this->querylength) {
    return true;
  } else if (this->plusp == true && this->fusion_queryend_junction != NULL) {
    return false;
  } else if (this->plusp == false && this->fusion_querystart_junction != NULL) {
    return false;
  } else {
    return false;
  }
}


bool
Path_unextended_qstart_p (T this, int endtrim_allowed, bool allow_ambig_p) {
  /* printf("Path_unextended_qstart_p => "); */
  if (this->qstart_alts != NULL) {
    /* printf("qstart_alts => false\n"); */
    return false;
  } else if (this->plusp == true && this->fusion_querystart_junction != NULL) {
    /* printf("fusion_querystart => false\n"); */
    return false;
  } else if (this->plusp == false && this->fusion_queryend_junction != NULL) {
    /* printf("fusion_queryend => false\n"); */
    return false;
  } else if (allow_ambig_p == true && this->splice5p == true) {
    /* printf("ambig5 => false\n"); */
    return false;
  } else if (Intlist_head(this->endpoints) <= endtrim_allowed) {
    /* printf("endpoints %d => false\n",Intlist_head(this->endpoints)); */
    return false;
  } else {
    /* printf("true\n"); */
    return true;
  }
}


bool
Path_unextended_qend_p (T this, int endtrim_allowed, bool allow_ambig_p) {
  /* printf("Path_unextended_qend_p => "); */
  if (this->qend_alts != NULL) {
    /* printf("qend_alts => false\n"); */
    return false;
  } else if (this->plusp == true && this->fusion_queryend_junction != NULL) {
    /*printf("fusion_queryend => false\n"); */
    return false;
  } else if (this->plusp == false && this->fusion_querystart_junction != NULL) {
    /* printf("fusion_querystart => false\n"); */
    return false;
  } else if (allow_ambig_p == true && this->splice3p == true) {
    /* printf("ambig3 => false\n"); */
    return false;
  } else if (this->querylength - Intlist_last_value(this->endpoints) <= endtrim_allowed) {
    /* printf("endpoints %d => false\n",Intlist_last_value(this->endpoints)); */
    return false;
  } else {
    /* printf("true\n"); */
    return true;
  }
}


bool
Path_unextended_querystart_p (T this, int endtrim_allowed, bool allow_ambig_p) {
  if (this->plusp == true) {
    return Path_unextended_qstart_p(this,endtrim_allowed,allow_ambig_p);
  } else {
    return Path_unextended_qend_p(this,endtrim_allowed,allow_ambig_p);
  }
}


bool
Path_unextended_queryend_p (T this, int endtrim_allowed, bool allow_ambig_p) {
  if (this->plusp == true) {
    return Path_unextended_qend_p(this,endtrim_allowed,allow_ambig_p);
  } else {
    return Path_unextended_qstart_p(this,endtrim_allowed,allow_ambig_p);
  }
}


bool
Path_unextendedp (T this, int endtrim_allowed, bool allow_ambig_p) {
  return false;

  if (this->transcriptome_method_p == true) {
    return false;
  } else if (Path_unextended_qstart_p(this,endtrim_allowed,allow_ambig_p) == true) {
    return true;
  } else if (Path_unextended_qend_p(this,endtrim_allowed,allow_ambig_p) == true) {
    return true;
  } else {
    return false;
  }
}

/* Checks for good splice sites */
bool
Path_completep (T this, int endtrim_allowed, bool allow_ambig_p) {
  return true;

  if (this->transcriptome_method_p == true) {
    return true;
  } else if (Path_unextended_qstart_p(this,endtrim_allowed,allow_ambig_p) == true) {
    return false;
  } else if (Path_unextended_qend_p(this,endtrim_allowed,allow_ambig_p) == true) {
    return false;
  } else if (Path_nbadsplices(this) > 0) {
    return false;
  } else {
    return true;
  }
}

bool
Path_unsolvedp (T this) {
  List_T p;

  for (p = this->junctions; p != NULL; p = List_next(p)) {
    if ((Junction_T) List_head(p) == JUNCTION_UNSOLVED) {
      return true;
    }
  }

  return false;
}


/* Used for NM:i in SAM output */
int
Path_ndiffs (T this) {
  int nmismatches, nindels = 0;
  List_T j;

  nmismatches = this->score_within_trims;
  for (j = this->junctions; j != NULL; j = List_next(j)) {
    /* Insertions and deletions count as one mismatch per base */
    nindels += Junction_nindels((Junction_T) List_head(j));
  }

  return nmismatches + nindels;
}


/* Used for Path_local_cmp and Pathpair_local_cmp */
int
Path_nindels (T this) {
  int nindels = 0;
  List_T j;

  for (j = this->junctions; j != NULL; j = List_next(j)) {
    /* Insertions and deletions count as one mismatch per base */
    nindels += Junction_nindels((Junction_T) List_head(j));
  }

  return nindels;
}


Chrpos_T
Path_chrlength (T this) {
  if (this == NULL) {
    /* Can happen if we call upon a mate in a halfmapping */
    return 0;
  } else if (circularp[this->chrnum] == true) {
    return (this->chrhigh - this->chroffset)/2;
  } else {
    return (this->chrhigh - this->chroffset);
  }
}


Univcoord_T
Path_genomiclow (T this) {
  Altsplice_T altsplice;

  if ((altsplice = this->qstart_alts) != NULL) {
    return Altsplice_best_univdiagonal(altsplice) - (Univcoord_T) this->querylength;
  } else if (Univcoordlist_head(this->univdiagonals) < (Univcoord_T) this->querylength) {
    return 0;
  } else {
    return Univcoordlist_head(this->univdiagonals) - this->querylength;
  }
}

Univcoord_T
Path_genomichigh (T this) {
  Altsplice_T altsplice;

  if ((altsplice = this->qend_alts) != NULL) {
    return Altsplice_best_univdiagonal(altsplice);
  } else {
    return Univcoordlist_last_value(this->univdiagonals);
  }
}


#if 0
/* To consider splice junctions, sites, or alts, use Path_effective_sensedir */
int
Path_sensedir (T this) {
  List_T j;
  Junction_T junction;

  if (this->qstart_alts != NULL) {
    return this->sensedir;
  } else if (this->qend_alts != NULL) {
    return this->sensedir;
  } else {
    for (j = this->junctions; j != NULL; j = List_next(j)) {
      junction = (Junction_T) List_head(j);
      if (Junction_type(junction) == SPLICE_JUNCTION) {
	return this->sensedir;
      }
    }
    return SENSE_NULL;
  }
}
#endif
    

unsigned int
Path_trnum_low (T this) {
  if (this->transcripts != NULL) {
    return ((Transcript_T) List_head(this->transcripts))->num;
  } else {
    /* Can occur after Trpath_convert calls Path_solve_junctions */
    return ((Transcript_T) List_head(this->invalid_transcripts))->num;
  }
}

unsigned int
Path_trnum_high (T this) {
  if (this->transcripts != NULL) {
    return ((Transcript_T) List_last_value(this->transcripts,NULL))->num;
  } else {
    /* Can occur after Trpath_convert calls Path_solve_junctions */
    return ((Transcript_T) List_last_value(this->invalid_transcripts,NULL))->num;
  }
}


int
Path_trnum_low_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  unsigned int trnum_low_x;
  unsigned int trnum_low_y;

  if (a->transcripts != NULL) {
    trnum_low_x = ((Transcript_T) List_head(a->transcripts))->num;
  } else {
    trnum_low_x = ((Transcript_T) List_head(a->invalid_transcripts))->num;
  }
    
  if (b->transcripts != NULL) {
    trnum_low_y = ((Transcript_T) List_head(b->transcripts))->num;
  } else {
    trnum_low_y = ((Transcript_T) List_head(b->invalid_transcripts))->num;
  }
    

  if (trnum_low_x < trnum_low_y) {
    return -1;
  } else if (trnum_low_y < trnum_low_x) {
    return +1;
  } else {
    return 0;
  }
}


int
Path_trnum_high_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  unsigned int trnum_high_x;
  unsigned int trnum_high_y;

  if (a->transcripts != NULL) {
    trnum_high_x = ((Transcript_T) List_last_value(a->transcripts,NULL))->num;
  } else {
    trnum_high_x = ((Transcript_T) List_last_value(a->invalid_transcripts,NULL))->num;
  }
    
  if (b->transcripts != NULL) {
    trnum_high_y = ((Transcript_T) List_last_value(b->transcripts,NULL))->num;
  } else {
    trnum_high_y = ((Transcript_T) List_last_value(b->invalid_transcripts,NULL))->num;
  }
  
  if (trnum_high_x < trnum_high_y) {
    return -1;
  } else if (trnum_high_y < trnum_high_x) {
    return +1;
  } else {
    return 0;
  }
}



#ifdef USE_HIGHLOW_UNIVDIAGONALS
int
Path_low_univdiagonal_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  Univcoord_T genomiclow_x = Path_low_univdiagonal(a);
  Univcoord_T genomiclow_y = Path_low_univdiagonal(b);

  if (genomiclow_x < genomiclow_y) {
    return -1;
  } else if (genomiclow_y < genomiclow_x) {
    return +1;
  } else {
    return 0;
  }
}


int
Path_high_univdiagonal_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  Univcoord_T genomichigh_x = Path_high_univdiagonal(a);
  Univcoord_T genomichigh_y = Path_high_univdiagonal(b);

  if (genomichigh_x < genomichigh_y) {
    return -1;
  } else if (genomichigh_y < genomichigh_x) {
    return +1;
  } else {
    return 0;
  }
}

#else
int
Path_main_univdiagonal_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->main_univdiagonal < b->main_univdiagonal) {
    return -1;
  } else if (b->main_univdiagonal < a->main_univdiagonal) {
    return +1;
  } else {
    return 0;
  }
}

#endif



#ifdef CHECK_ASSERTIONS
static void
check_ascending (Univcoord_T *coords, int n) {
  Univcoord_T prev_coord;
  int i;

  prev_coord = coords[0];
  for (i = 1; i < n; i++) {
    if (coords[i] <= prev_coord) {
      printf("Expecting forward, but at %d, got %u <= %u\n",
	     i,coords[i],prev_coord);
      abort();
    }
    prev_coord = coords[i];
  }
 
  return;
}
#endif



#ifdef USE_HIGHLOW_UNIVDIAGONALS
#ifdef CHECK_ASSERTIONS
static void
check_unique_coords_low_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int nunique) {
  int last_cumsum = 0;
  int k, i;

  for (k = 0; k < nunique; k++) {
    assert(indices[k] > last_cumsum);
    for (i = last_cumsum; i < indices[k]; i++) {
      assert(coords[k] == Path_low_univdiagonal(paths[i]));
    }
    last_cumsum = indices[k];
  }

  return;
}
#endif


#ifdef CHECK_ASSERTIONS
static void
check_unique_coords_high_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int nunique) {
  int last_cumsum = 0;
  int k, i;

  for (k = 0; k < nunique; k++) {
    assert(indices[k] > last_cumsum);
    for (i = last_cumsum; i < indices[k]; i++) {
      assert(coords[k] == Path_high_univdiagonal(paths[i]));
    }
    last_cumsum = indices[k];
  }

  return;
}
#endif
#else
#ifdef CHECK_ASSERTIONS
static void
check_unique_coords_main_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int nunique) {
  int last_cumsum = 0;
  int k, i;

  for (k = 0; k < nunique; k++) {
    assert(indices[k] > last_cumsum);
    for (i = last_cumsum; i < indices[k]; i++) {
      assert(coords[k] == Path_main_univdiagonal(paths[i]));
    }
    last_cumsum = indices[k];
  }

  return;
}
#endif
#endif


#ifdef USE_HIGHLOW_UNIVDIAGONALS
/* Handles duplicate coords by storing a cumsum of the counts */
int
Path_fill_low_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int n) {
  Univcoord_T last_coord, coord;
  int cumsum = 0, i, nunique = 0;

  last_coord = Path_low_univdiagonal(paths[0]);
  /* cumsum++; */

  for (i = 0; i < n; i++) {
    if ((coord = Path_low_univdiagonal(paths[i])) != last_coord) {
      coords[nunique] = last_coord;
      indices[nunique++] = cumsum;
      last_coord = coord;
    }
    cumsum++;
  }

  /* Last entry */
  coords[nunique] = last_coord;
  indices[nunique++] = cumsum;

#ifdef CHECK_ASSERTIONS
  check_unique_coords_low_univdiagonal(coords,indices,paths,nunique);
  check_ascending(coords,nunique);
#endif

  return nunique;
}


/* Handles duplicate coords by storing a cumsum of the counts */
int
Path_fill_high_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int n) {
  Univcoord_T last_coord, coord;
  int cumsum = 0, i, nunique = 0;

  last_coord = Path_high_univdiagonal(paths[0]);
  /* cumsum++; */

  for (i = 0; i < n; i++) {
    if ((coord = Path_high_univdiagonal(paths[i])) != last_coord) {
      coords[nunique] = last_coord;
      indices[nunique++] = cumsum;
      last_coord = coord;
    }
    cumsum++;
  }

  /* Last entry */
  coords[nunique] = last_coord;
  indices[nunique++] = cumsum;

#ifdef CHECK_ASSERTIONS
  check_unique_coords_high_univdiagonal(coords,indices,paths,nunique);
  check_ascending(coords,nunique);
#endif

  return nunique;
}

#else
/* Handles duplicate coords by storing a cumsum of the counts */
int
Path_fill_main_univdiagonal (Univcoord_T *coords, int *indices, T *paths, int n) {
  Univcoord_T last_coord, coord;
  int cumsum = 0, i, nunique = 0;

  last_coord = Path_main_univdiagonal(paths[0]);
  /* cumsum++; */

  for (i = 0; i < n; i++) {
    if ((coord = Path_main_univdiagonal(paths[i])) != last_coord) {
      coords[nunique] = last_coord;
      indices[nunique++] = cumsum;
      last_coord = coord;
    }
    cumsum++;
  }

  /* Last entry */
  coords[nunique] = last_coord;
  indices[nunique++] = cumsum;

#ifdef CHECK_ASSERTIONS
  check_unique_coords_main_univdiagonal(coords,indices,paths,nunique);
  check_ascending(coords,nunique);
#endif

  return nunique;
}
#endif



int
Path_interval_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  Univcoord_T genomiclow_x, genomiclow_y,
    genomichigh_x, genomichigh_y;

  genomiclow_x = Path_genomiclow(a);
  genomiclow_y = Path_genomiclow(b);

  if (genomiclow_x < genomiclow_y) {
    return -1;
  } else if (genomiclow_y < genomiclow_x) {
    return +1;
  } else {
    genomichigh_x = Path_genomichigh(a);
    genomichigh_y = Path_genomichigh(b);
    if (genomichigh_x > genomichigh_y) {
      return -1;
    } else if (genomichigh_y > genomichigh_x) {
      return +1;
    } else {
      return 0;
    }
  }
}


int
Path_structure_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  Univcoordlist_T p, q;
  Intlist_T i, j;
  List_T r, s;
  int nindels_x, nindels_y;
  int sensedir_a, sensedir_b;

  /* Need to use effective sensedirs, not sensedir field */
  if ((sensedir_a = Path_effective_sensedir(a)) > (sensedir_b = Path_effective_sensedir(b))) {
    return -1;
  } else if (sensedir_b > sensedir_a) {
    return +1;
  } else {
    p = a->univdiagonals;
    q = b->univdiagonals;

    while (p != NULL && q != NULL) {
      if (Univcoordlist_head(p) < Univcoordlist_head(q)) {
	return -1;
      } else if (Univcoordlist_head(q) < Univcoordlist_head(p)) {
	return +1;
      } else {
	p = Univcoordlist_next(p);
	q = Univcoordlist_next(q);
      }
    }

    if (p == NULL && q != NULL) {
      return -1;
    } else if (p != NULL && q == NULL) {
      return +1;
    } else {

      /* Check for internal endpoint structure, which could mean
	 different introns, but ignore the first and last endpoints,
	 which could mean different trimming */

      i = Intlist_next(a->endpoints);
      j = Intlist_next(b->endpoints);
      r = a->junctions;
      s = b->junctions;
      
      while (Intlist_next(i) != NULL && Intlist_next(j) != NULL) {
	if ((nindels_x = Junction_nindels((Junction_T) List_head(r))) < (nindels_y = Junction_nindels((Junction_T) List_head(s)))) {
	  return -1;
	} else if (nindels_y < nindels_x) {
	  return +1;
	} else if (nindels_x > 0) {
	  /* Allow differences in endpoints for an indel */
	  i = Intlist_next(i);
	  j = Intlist_next(j);
	  r = List_next(r);
	  s = List_next(s);
	    
	} else if (Intlist_head(i) < Intlist_head(j)) {
	  /* Difference in splice endpoints */
	  return -1;
	} else if (Intlist_head(j) < Intlist_head(i)) {
	  /* Difference in splice endpoints */
	  return +1;
	  /* Same splice endpoints */
	} else {
	  i = Intlist_next(i);
	  j = Intlist_next(j);
	  r = List_next(r);
	  s = List_next(s);
	}
      }

      if (Intlist_next(i) == NULL && Intlist_next(j) != NULL) {
	return -1;
      } else if (Intlist_next(i) != NULL && Intlist_next(j) == NULL) {
	return +1;
      } else {
	return 0;
      }
    }
  }
}


int
Path_structure_ignore_sense_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;
  Univcoordlist_T p, q;
  Intlist_T i, j;

  p = a->univdiagonals;
  q = b->univdiagonals;

  while (p != NULL && q != NULL) {
    if (Univcoordlist_head(p) < Univcoordlist_head(q)) {
      return -1;
    } else if (Univcoordlist_head(q) < Univcoordlist_head(p)) {
      return +1;
    } else {
      p = Univcoordlist_next(p);
      q = Univcoordlist_next(q);
    }
  }

  if (p == NULL && q != NULL) {
    return -1;
  } else if (p != NULL && q == NULL) {
    return +1;
  } else {

    /* Check for internal endpoint structure, which could mean
       different introns, but ignore the first and last endpoints,
       which could mean different trimming */

    i = Intlist_next(a->endpoints);
    j = Intlist_next(b->endpoints);
    
    while (Intlist_next(i) != NULL && Intlist_next(j) != NULL) {
      if (Intlist_head(i) < Intlist_head(j)) {
	return -1;
      } else if (Intlist_head(j) < Intlist_head(i)) {
	return +1;
      } else {
	i = Intlist_next(i);
	j = Intlist_next(j);
      }
    }
    
    if (Intlist_next(i) == NULL && Intlist_next(j) != NULL) {
      return -1;
    } else if (Intlist_next(i) != NULL && Intlist_next(j) == NULL) {
      return +1;
    } else {
      return 0;
    }
  }
}


/* Called only by stage1hr-paired */
static int
optimal_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->nmatches > b->nmatches) {
    return -1;
  } else if (b->nmatches > a->nmatches) {
    return +1;
  } else if (a->junction_splice_prob > b->junction_splice_prob) {
    return -1;
  } else if (b->junction_splice_prob > a->junction_splice_prob) {
    return +1;
  } else if (a->method > b->method) {
    return -1;
  } else if (b->method > a->method) {
    return +1;
  } else {
    return 0;
  }
}


static bool
Path_identical_mind_sensedir_p (T a, T b) {
  if (a->sensedir != b->sensedir) {
    return false;
  } else if (a->nmatches != b->nmatches) {
    return false;
  } else if (a->junction_splice_prob != b->junction_splice_prob) {
    return false;
  } else if (a->method != b->method) {
    return false;
  } else {
    return true;
  }
}


bool
Path_overlap_p (T x, T y) {
  Univcoord_T genomiclow_x, genomiclow_y,
    genomichigh_x, genomichigh_y;

  genomiclow_x = Path_genomiclow(x);
  genomiclow_y = Path_genomiclow(y);
  genomichigh_x = Path_genomichigh(x);
  genomichigh_y = Path_genomichigh(y);
  
  if (genomichigh_x < genomiclow_y) {
    return false;
  } else if (genomichigh_y < genomiclow_x) {
    return false;
  } else {
    return true;
  }
}


int
Path_max_trim (T this) {
  int overall_qstart, overall_qend;

  overall_qstart = Intlist_head(this->endpoints);
  overall_qend = Intlist_last_value(this->endpoints);
  if (overall_qstart > this->querylength - overall_qend) {
    return overall_qstart;
  } else {
    return this->querylength - overall_qend;
  }
}


/* Called only by stage1hr-paired */
/* Removes duplicates and subsumed paths.  Keeps both sensedirs, but
   called only on lists of a single sensedir anyway */
/* Frees paths, so needs to be called on a stage1 list and returned to a stage1 list */
/* Not sure if this works with extended paths */
List_T
Path_filter (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool) {
  List_T list;
  T *array, path;
  Univcoord_T genomichigh;
  int n, i, j, k;

  if ((n = List_length(paths)) == 0) {
    return (List_T) NULL;

  } else {
    list = (List_T) NULL;

    array = (T *) List_to_array(paths,NULL);
    qsort(array,n,sizeof(T),Path_interval_cmp);

    k = 0;
    i = 0;
    while (i < n) {
      genomichigh = Path_genomichigh(array[i]);
      j = i + 1;
      while (j < n && Path_genomiclow(array[j]) <= genomichigh) {
	j++;
      }

      if (j - i > 1) {
	qsort(&(array[i]),j - i,sizeof(T),optimal_cmp);
      }

      list = Hitlist_push(list,hitlistpool,(void *) array[i]
			  hitlistpool_trace(__FILE__,__LINE__));
      for (k = i + 1; k < j; k++) {
	if (Path_identical_mind_sensedir_p(array[k],array[i]) == true) {
	  path = array[k];
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	} else if (optimal_cmp(&(array[k]),&(array[i])) == 0) {
	  list = Hitlist_push(list,hitlistpool,(void *) array[k]
			      hitlistpool_trace(__FILE__,__LINE__));
	} else {
	  path = array[k];
	  Path_free(&path,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
	}
      }

      i = j;
    }
    FREE(array);
    Hitlistpool_free_list(&paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));  /* Assigned by Hitlist_push */

    debug0(printf("\n"));
    return List_reverse(list);
  }
}


static int
setdiff_scalar (Univcoord_T *old, Univcoord_T *end_old,
		Auxinfo_T *old_auxinfo, Pathstore_T *old_pathstores,
		Univcoord_T *new, Univcoord_T *end_new,
		Auxinfo_T *new_auxinfo, Pathstore_T *new_pathstores,
		Pathstorepool_T pathstorepool, Intlistpool_T intlistpool) {

  Univcoord_T *new_start, *ptr = new;
  Auxinfo_T *ptr_auxinfo = new_auxinfo;
  Pathstore_T *ptr_pathstores = new_pathstores;

  new_start = new;
  while (old < end_old && new < end_new) {
    debug9(printf("Comparing %u and %u: ",(*old),(*new)));
    if ((*old) < (*new)) {
      debug9(printf("Old is first.  Not keeping\n"));
      old++; old_auxinfo++; old_pathstores++;

    } else if ((*new) < (*old)) {
      /* Overwrite at top of the new list */
      debug9(printf("New is first.  Keeping\n"));
      *ptr++ = *new++;
      *ptr_auxinfo++ = *new_auxinfo++;
      *ptr_pathstores++ = *new_pathstores++;

    } else {
      /* Equal: Do not put into the list, but merge new auxinfo into old */
      old++; new++;
      *old_auxinfo = Auxinfo_append(*old_auxinfo,*new_auxinfo);
      old_auxinfo++; new_auxinfo++;
      *old_pathstores = Pathstore_merge(*old_pathstores,&(*new_pathstores),pathstorepool,intlistpool);
      old_pathstores++; new_pathstores++;
    }
  }
      
  while (new < end_new) {
    debug9(printf("Observed new %u: Keeping\n"));
    *ptr++ = *new++;
    *ptr_auxinfo++ = *new_auxinfo++;
    *ptr_pathstores++ = *new_pathstores++;
  }

  return (ptr - new_start);
}



/* Modifies old_auxinfo by appending the new_auxinfo to its end.
   Overwrites new_univdiagonals and new_auxinfo.  Returns nnovel */
int
Path_setdiff_univdiagonals_auxinfo (Univcoord_T *old_univdiagonals, Auxinfo_T *old_auxinfo,
				    Pathstore_T *old_pathstores, int nold,
				    Univcoord_T *new_univdiagonals, Auxinfo_T *new_auxinfo,
				    Pathstore_T *new_pathstores, int nnew,
				    Pathstorepool_T pathstorepool, Intlistpool_T intlistpool) {
  int nnovel;

#ifdef DEBUG9
  printf("Entering Path_setdiff_univdiagonals_auxinfo with\n");
  for (int i = 0; i < nold; i++) {
    printf("old univdiagonal %u, auxinfo ",old_univdiagonals[i]);
    Auxinfo_print(old_auxinfo[i]);
    Pathstore_print(old_pathstores[i]);
    printf("\n");
  }
  for (int i = 0; i < nnew; i++) {
    printf("new univdiagonal %u, auxinfo ",new_univdiagonals[i]);
    Auxinfo_print(new_auxinfo[i]);
    Pathstore_print(new_pathstores[i]);
    printf("\n");
  }
#endif

  if (nold == 0) {
    /* No change */
#ifdef DEBUG9
    printf("(1) Exiting Path_setdiff_univdiagonals_auxinfo with\n");
    for (int i = 0; i < nnew; i++) {
      printf("new univdiagonal %u, auxinfo ",new_univdiagonals[i]);
      Auxinfo_print(new_auxinfo[i]);
      Pathstore_print(new_pathstores[i]);
      printf("\n");
    }
#endif

    return nnew;

  } else if (nnew == 0) {
    /* No change */
#ifdef DEBUG9
    printf("(2) Exiting Path_setdiff_univdiagonals_auxinfo with\n");
    for (int i = 0; i < nnew; i++) {
      printf("new univdiagonal %u, auxinfo ",new_univdiagonals[i]);
      Auxinfo_print(new_auxinfo[i]);
      Pathstore_print(new_pathstores[i]);
      printf("\n");
    }
#endif

    return 0;
    
  } else {
#ifdef DEBUG9
    for (int i = 0; i < nold; i++) {
      printf("Old %u",old_univdiagonals[i]);
      for (Auxinfo_T p = old_auxinfo[i]; p != NULL; p = p->rest) {
	printf(" %p %s",p,Method_string(p->method));
      }
      printf("\n");
    }
    for (int i = 0; i < nnew; i++) {
      printf("New %u",new_univdiagonals[i]);
      for (Auxinfo_T p = new_auxinfo[i]; p != NULL; p = p->rest) {
	printf(" %p %s",p,Method_string(p->method));
      }
      printf("\n");
    }
#endif

    nnovel = setdiff_scalar(old_univdiagonals,&(old_univdiagonals[nold]),
			    old_auxinfo,old_pathstores,
			    new_univdiagonals,&(new_univdiagonals[nnew]),
			    new_auxinfo,new_pathstores,
			    pathstorepool,intlistpool);

#ifdef DEBUG9
    printf("(3) Exiting Path_setdiff_univdiagonals_auxinfo with\n");
    for (int i = 0; i < nnovel; i++) {
      printf("new univdiagonal %u, auxinfo ",new_univdiagonals[i]);
      Auxinfo_print(new_auxinfo[i]);
      Pathstore_print(new_pathstores[i]);
      printf("\n");
    }
#endif

    return nnovel;
  }
}



static int
merge_scalar (Univcoord_T *out, Auxinfo_T *out_auxinfo, Pathstore_T *out_pathstores,
	      Univcoord_T *A, Univcoord_T *end_A,
	      Auxinfo_T *A_auxinfo, Pathstore_T *A_pathstores,
	      Univcoord_T *B, Univcoord_T *end_B,
	      Auxinfo_T *B_auxinfo, Pathstore_T *B_pathstores,
	      Pathstorepool_T pathstorepool, Intlistpool_T intlistpool) {

  Univcoord_T *out_start;
  
  out_start = out;
  while (A < end_A && B < end_B) {
    debug11(printf("Comparing %u and %u\n",(*A),(*B)));
    if ((*A) < (*B)) {
      *out++ = *A++;
      *out_auxinfo++ = *A_auxinfo++;
      *out_pathstores++ = *A_pathstores++;
    } else if ((*B) < (*A)) {
      *out++ = *B++;
      *out_auxinfo++ = *B_auxinfo++;
      *out_pathstores++ = *B_pathstores++;
    } else {
      *out++ = *A++;
      *out_auxinfo++ = Auxinfo_append(*A_auxinfo,*B_auxinfo);
      *out_pathstores++ = Pathstore_merge(*A_pathstores,&(*B_pathstores),
					  pathstorepool,intlistpool);
      A_auxinfo++; B_auxinfo++;
      A_pathstores++; B_pathstores++;
      B++;
    }
  }
      
  while (A < end_A) {
    debug11(printf("Moving %u from A\n",(*A)));
    *out++ = *A++;
    *out_auxinfo++ = *A_auxinfo++;
    *out_pathstores++ = *A_pathstores++;
  }

  while (B < end_B) {
    debug11(printf("Moving %u from B\n",(*B)));
    *out++ = *B++;
    *out_auxinfo++ = *B_auxinfo++;
    *out_pathstores++ = *B_pathstores++;
  }

  return (out - out_start);
}


static int
merge_rare_freq (Univcoord_T *out, Auxinfo_T *out_auxinfo,
		 Univcoord_T *rare, Univcoord_T *end_rare, Auxinfo_T *rare_auxinfo,
		 Univcoord_T *freq, Univcoord_T *end_freq, Auxinfo_T *freq_auxinfo) {

  Univcoord_T val_rare, *freq_start, *out_start;
  int ncopy;
  
  out_start = out;
  while (rare < end_rare && freq < end_freq) {
    val_rare = *rare;
    freq_start = freq;
    while (freq < end_freq && (*freq) < val_rare) {
      freq++;
    }

    ncopy = freq - freq_start;
    debug11(printf("Copying %d from freq\n",ncopy));

    memcpy(out,freq_start,ncopy*sizeof(Univcoord_T));
    memcpy(out_auxinfo,freq_auxinfo,ncopy*sizeof(Auxinfo_T));
    freq_auxinfo += ncopy;
    out += ncopy;
    out_auxinfo += ncopy;
    
    if (freq >= end_freq) {
      debug11(printf("Reached end of freq, so copying rare\n"));
      *out++ = val_rare;
      *out_auxinfo++ = *rare_auxinfo++;

    } else if (val_rare < (*freq)) {
      debug11(printf("Passed rare, so copying rare.  Not advancing freq\n"));
      *out++ = val_rare;
      *out_auxinfo++ = *rare_auxinfo++;

    } else {
      /* Equal */
      debug11(printf("Equal\n"));
      *out++ = val_rare;
      *out_auxinfo++ = Auxinfo_append(*rare_auxinfo++,*freq_auxinfo++);
      freq++;
    }

    rare++;
  }

  if (rare < end_rare) {
    ncopy = end_rare - rare;
    debug11(printf("Finishing %d elts on rare\n",ncopy));
    memcpy(out,rare,ncopy*sizeof(Univcoord_T));
    memcpy(out_auxinfo,rare_auxinfo,ncopy*sizeof(Auxinfo_T));
    out += ncopy;

  } else if (freq < end_freq) {
    ncopy = end_freq - freq;
    debug11(printf("Finishing %d elts on freq\n",ncopy));
    memcpy(out,freq,ncopy*sizeof(Univcoord_T));
    memcpy(out_auxinfo,freq_auxinfo,ncopy*sizeof(Auxinfo_T));
    out += ncopy;
  }

  return (out - out_start);
}


/* new_univdiagonals are aligned, coming from Kmer_exact1, Extension_search, or Kmer_segment */
/* Because we transfer pointer to all_univdiagonals, it is also aligned.  And therefore _old_univdiagonals are aligned */
void
Path_merge_univdiagonals_auxinfo (Univcoord_T **_all_univdiagonals, Auxinfo_T **all_auxinfo,
				  Pathstore_T **all_pathstores, int *nall,
				  Univcoord_T *_new_univdiagonals, Auxinfo_T *new_auxinfo,
				  Pathstore_T *new_pathstores, int nnew,
				  Pathstorepool_T pathstorepool, Intlistpool_T intlistpool) {
  Univcoord_T *_old_univdiagonals;
  Auxinfo_T *old_auxinfo;
  Pathstore_T *old_pathstores;
  int nold;

  debug11(printf("Entering Path_merge_univdiagonals_auxinfo\n"));

  if ((*nall) + nnew == 0) {
    /* printf("(1) KEEPING *nall to be ZERO\n"); */
    /* It is possible that nnew (really nnovel) is 0, but new_univdiagonals is not NULL */
    /* printf("FREEING %p\n",_new_univdiagonals); */
    FREE_ALIGN(_new_univdiagonals);
    FREE(new_auxinfo);
    FREE(new_pathstores);
    return;

  } else if (nnew == 0) {
    /* printf("(2) KEEPING all to be same with %d entries\n",*nall); */
    /* It is possible that nnew (really nnovel) is 0, but new_univdiagonals is not NULL */
    /* printf("FREEING %p\n",_new_univdiagonals); */
    FREE_ALIGN(_new_univdiagonals);
    FREE(new_auxinfo);
    FREE(new_pathstores);
    return;
    
  } else if ((*nall) == 0) {
    /* printf("(3) TRANSFERRING new to all with %d entries\n",nnew); */
    /* printf("Re-assigning %p to all\n",new_univdiagonals); */
    /* printf("TRANSFERRING %p to ALL\n",_new_univdiagonals); */
    FREE_ALIGN(*_all_univdiagonals);
    FREE(*all_auxinfo);
    FREE(*all_pathstores);
    *_all_univdiagonals = _new_univdiagonals;
    *all_auxinfo = new_auxinfo;
    *all_pathstores = new_pathstores;
    *nall = nnew;
    return;

  } else {
    /* Keep pointers to _all_univdiagonals and all_auxinfo */
    _old_univdiagonals = *_all_univdiagonals;
    old_auxinfo = *all_auxinfo;
    old_pathstores = *all_pathstores;
    nold = *nall;

#ifdef DEBUG11
    for (int i = 0; i < nold; i++) {
      printf("Old %u",_old_univdiagonals[i]);
      for (Auxinfo_T p = old_auxinfo[i]; p != NULL; p = p->rest) {
	printf(" %p %s",p,Method_string(p->method));
      }
      printf("\n");
    }
    for (int i = 0; i < nnew; i++) {
      printf("New %u",_new_univdiagonals[i]);
      for (Auxinfo_T p = new_auxinfo[i]; p != NULL; p = p->rest) {
	printf(" %p %s",p,Method_string(p->method));
      }
      printf("\n");
    }
#endif

    MALLOC_ALIGN(*_all_univdiagonals,(nold+nnew)*sizeof(Univcoord_T));
    *all_auxinfo = (Auxinfo_T *) MALLOC((nold+nnew)*sizeof(Auxinfo_T));
    *all_pathstores = (Pathstore_T *) MALLOC((nold+nnew)*sizeof(Pathstore_T));

#if 1
    *nall = merge_scalar(/*out*/*_all_univdiagonals,*all_auxinfo,*all_pathstores,
  		         /*A*/_old_univdiagonals,&(_old_univdiagonals[nold]),
			 old_auxinfo,old_pathstores,
			 /*B*/_new_univdiagonals,&(_new_univdiagonals[nnew]),
			 new_auxinfo,new_pathstores,pathstorepool,intlistpool);
#else
    if (nold < 8 && nnew < 8) {
      *nall = merge_scalar(/*out*/*_all_univdiagonals,*all_auxinfo,
			   /*A*/_old_univdiagonals,&(_old_univdiagonals[nold]),old_auxinfo,
			   /*B*/_new_univdiagonals,&(_new_univdiagonals[nnew]),new_auxinfo);
    } else if (nold <= nnew) {
      /* Results in a memory leak */
      *nall = merge_rare_freq(/*out*/*_all_univdiagonals,*all_auxinfo,
			      /*rare*/_old_univdiagonals,&(_old_univdiagonals[nold]),old_auxinfo,
			      /*freq*/_new_univdiagonals,&(_new_univdiagonals[nnew]),new_auxinfo);
    } else {
      /* Results in a memory leak */
      *nall = merge_rare_freq(/*out*/*_all_univdiagonals,*all_auxinfo,
			      /*rare*/_new_univdiagonals,&(_new_univdiagonals[nnew]),new_auxinfo,
			      /*freq*/_old_univdiagonals,&(_old_univdiagonals[nold]),old_auxinfo);
    }
#endif

#ifdef DEBUG11
    printf("(4) Path_merge returning %d univdiagonals/auxinfo\n",*nall);
    for (int i = 0; i < *nall; i++) {
      printf("%u",(*_all_univdiagonals)[i]);
      for (Auxinfo_T p = (*all_auxinfo)[i]; p != NULL; p = p->rest) {
	printf(" %p",p);
      }
      printf("\n");
      /* printf("%u %s\n",(*_all_univdiagonals)[i],Method_string((*all_auxinfo)[i]->method)); */
    }
    printf("\n");
#endif

    FREE_ALIGN(_new_univdiagonals);
    FREE(new_auxinfo);
    FREE(new_pathstores);
    FREE_ALIGN(_old_univdiagonals);
    FREE(old_auxinfo);
    FREE(old_pathstores);
    
    return;
  }
}



static int
Path_method_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->method > b->method) {
    return -1;
  } else if (b->method > a->method) {
    return +1;
  } else {
    return 0;
  }
}


/* Using instead of Path_filter.  Removes only reads that are structurally identical */
List_T
Path_unique (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool) {
  List_T list;
  T *patharray, path;
  int n, i, j, l;

  if ((n = List_length(paths)) == 0) {
    return (List_T) NULL;

  } else if (n == 1) {
    return paths;

  } else {
    list = (List_T) NULL;
    patharray = (T *) List_to_array(paths,NULL);
    qsort(patharray,n,sizeof(T),Path_structure_cmp);
    
    i = 0;
    while (i < n) {
      j = i + 1;
      while (j < n && Path_structure_cmp(&(patharray[j]),&(patharray[i])) == 0) {
	j++;
      }
      debug4(printf("Found an identical group by structure (except sensedir) of %d paths => Re-sorting by method_cmp\n",j - i));
      
      qsort(&(patharray[i]),j - i,sizeof(T),Path_method_cmp);
      debug4(printf("(0) Keeping by method_cmp\n")); debug4(Path_print(patharray[i]));
      list = Hitlist_push(list,hitlistpool,(void *) patharray[i]
			  hitlistpool_trace(__FILE__,__LINE__));

      for (l = i + 1; l < j; l++) {
	debug4(printf("(0) Eliminating by method_cmp\n")); debug4(Path_print(patharray[l]));
	path = patharray[l];
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      }

      i = j;
    }
    FREE(patharray);

    debug4(printf("After removing duplicates (except for sensedir), have %d paths\n",List_length(list)));

    return List_reverse(list);
  }
}


#if 0
/* Removes duplicates and subsumed paths */
int
Path_filter_array (List_T *duplicates, T *paths, int npaths,
		   Hitlistpool_T hitlistpool) {
  T *out, path;
  Univcoord_T genomichigh;
  int i, j, k;

  if (npaths == 0) {
    return 0;

  } else {
    /* Caller, which is Ladder_paths_for_score, already allocates paths, so just re-use */
    /* out = unique = (T *) MALLOC((npaths+1)*sizeof(T)); */

    qsort(paths,npaths,sizeof(T),Path_interval_cmp);

    out = paths;
    k = 0;
    i = 0;
    while (i < npaths) {
      genomichigh = Path_genomichigh(paths[i]);

      j = i + 1;
      while (j < npaths && Path_genomiclow(paths[j]) <= genomichigh) {
	j++;
      }

      if (j - i > 1) {
	qsort(&(paths[i]),j - i,sizeof(T),optimal_cmp);
      }

      *out++ = paths[i];
      for (k = i + 1; k < j; k++) {
	if (Path_identical_mind_sensedir_p(paths[k],paths[i]) == true) {
	  path = paths[k];
	  /* Path_free(&path); -- newladder and ladder share paths during Concordance_byscore, so store for freeing later */
	  /* printf("Pushing path %p into duplicates\n",path); */
	  *duplicates = Hitlist_push(*duplicates,hitlistpool,(void *) path
				     hitlistpool_trace(__FILE__,__LINE__));
	} else if (optimal_cmp(&(paths[k]),&(paths[i])) == 0) {
	  *out++ = paths[k];
	} else {
	  path = paths[k];
	  /* Path_free(&path); -- newladder and ladder share paths during Concordance_byscore, so store for freeing later */
	  /* printf("Pushing path %p into duplicates\n",path); */
	  *duplicates = Hitlist_push(*duplicates,hitlistpool,(void *) path
				     hitlistpool_trace(__FILE__,__LINE__));
	}
      }

      i = j;
    }

    /* FREE(paths); */

    return out - paths;
  }
}
#endif


/* Compares across all loci.  Relies upon effective_sensedir rather than try_sensedir */
List_T
Path_optimal_nmatches (List_T paths, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Hitlistpool_T hitlistpool) {
  List_T list, p;
  T path;
  int max_nmatches;

  if (List_length(paths) == 0) {
    return (List_T) NULL;

  } else {
    max_nmatches = 0;

    for (p = paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (path->nmatches > max_nmatches) {
	max_nmatches = path->nmatches;
      }
    }

    list = NULL;
    for (p = paths; p != NULL; p = List_next(p)) {
      path = (T) List_head(p);
      if (path->nmatches < max_nmatches) {
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
      } else {
	list = Hitlist_push(list,hitlistpool,(void *) path
			    hitlistpool_trace(__FILE__,__LINE__));
      }
    }

    return List_reverse(list);
  }
}


List_T
Path_array_to_list (Path_T *paths, int n, Hitlistpool_T hitlistpool) {
  List_T list = NULL;
  int i;

  for (i = 0; i < n; i++) {
    list = Hitlist_push(list,hitlistpool,(void *) paths[i]
			hitlistpool_trace(__FILE__,__LINE__));
  }

  return list;
}


/* expect_fwd_p holds for the output, not the input */
T
Path_reverse (T this, bool expect_fwd_p) {
  Intlist_T q;
#ifdef CHECK_ASSERTIONS
  int prev_endpoint;
#endif

  this->endpoints = Intlist_reverse(this->endpoints);
  this->univdiagonals = Univcoordlist_reverse(this->univdiagonals);
  this->nmismatches = Intlist_reverse(this->nmismatches);
  this->ref_nmismatches = Intlist_reverse(this->ref_nmismatches);
  this->junctions = List_reverse(this->junctions);
  
#ifdef CHECK_ASSERTIONS
  if (expect_fwd_p == true) {
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
 
 } else {
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
  }
#endif

  return this;
}


/* Sometimes merging of left and right paths can result in anomalies */
/* Same as Trpath_endpoints_acceptable_p */
bool
Path_endpoints_acceptable_p (Intlist_T endpoints, List_T junctions) {
  Intlist_T p;
  List_T q;
  Junction_T junction;
  int last_endpoint;

  debug2(printf("Evaluating endpoints for acceptability: %s\n",Intlist_to_string(endpoints)));

  /* last_endpoint = 0; */
  /* Skip first endpoint */
  for (p = Intlist_next(endpoints), q = junctions; Intlist_next(p) != NULL; p = Intlist_next(p), q = List_next(q)) {
    last_endpoint = Intlist_head(p);
    junction = (Junction_T) List_head(q);
    /* Previously used >=, but a deletion can yield duplicated endpoints */
    if (last_endpoint + Junction_ninserts(junction) > Intlist_head(Intlist_next(p))) {
      debug2(printf("Endpoint %d + %d >= %d, so unacceptable\n",
		    last_endpoint,Junction_ninserts(junction),Intlist_head(Intlist_next(p))));
      return false;
    } else {
      debug2(printf("Endpoint %d + %d < %d, so acceptable\n",
		    last_endpoint,Junction_ninserts(junction),Intlist_head(Intlist_next(p))));
    }
  }

  return true;
}


/* Called by conversion procedures in trpath-convert.c */
T
Path_convert_simple (Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		     Intlist_T endpoints, Univcoordlist_T univdiagonals, Intlist_T nmismatches,
		     Intlist_T ref_nmismatches, List_T junctions, bool evalp,
		     bool plusp, bool first_read_p, int genestrand, int sensedir, int querylength,
		     Compress_T query_compress_fwd, Compress_T query_compress_rev,
		     Listpool_T listpool, Pathpool_T pathpool, Method_T method) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
  Intlist_T q;
#ifdef CHECK_ASSERTIONS
  int prev_endpoint;
#endif

  assert(sensedir == SENSE_FORWARD || sensedir == SENSE_ANTI);

  assert(Univcoordlist_length(univdiagonals) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(nmismatches) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(ref_nmismatches) == Intlist_length(endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(endpoints) - 2);

  /* Avoids having to call Path_eval_nmatches */

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = querylength;
  new->score_within_trims = querylength;
  new->genomic_diff = (char *) NULL;
  if (evalp == false) {
    new->nmatches = new->ref_nmatches = Intlist_last_value(endpoints) - Intlist_head(endpoints) - Intlist_sum(nmismatches) - Junction_total_ninserts(junctions);
#if 0
    printf("endpoints: %s\n",Intlist_to_string(endpoints));
    printf("Setting nmatches to be endpoints (%d - %d) - nmismatches %d - ninserts %d\n",
	   Intlist_last_value(endpoints),Intlist_head(endpoints),Intlist_sum(nmismatches),Junction_total_ninserts(junctions));
#endif
  } else {
    new->nmatches = -1;
  }

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = sensedir;
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  /* Use furthest ends, to allow cases where only the tails of each read overlap */
  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = Univcoordlist_last_value(univdiagonals);
    } else if (Univcoordlist_head(univdiagonals) < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = Univcoordlist_head(univdiagonals) - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = Univcoordlist_last_value(univdiagonals);
    } else if (Univcoordlist_head(univdiagonals) < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = Univcoordlist_head(univdiagonals) - querylength;
    }
  }

  debug0(printf("%d: Creating path %p, %u, by Path_convert_simple, plusp %d\n",
		++call_i,new,new->main_univdiagonal,plusp));
  
  new->endpoints = endpoints;
  new->univdiagonals = univdiagonals;
  new->nmismatches = nmismatches;
  new->ref_nmismatches = ref_nmismatches;
  new->junctions = junctions;

  new->splice5p = false;
  new->splicetype5 = NO_SPLICE;
  new->ambig_prob_5 = 0.0;

  new->splice3p = false;
  new->splicetype3 = NO_SPLICE;
  new->ambig_prob_3 = 0.0;

  new->qstart_alts = (Altsplice_T) NULL;
  new->qend_alts = (Altsplice_T) NULL;

  new->circular_endpoints = (Intlist_T) NULL;
#if 0
  new->circular_high_p = false;
  new->circular_nmismatches = (Intlist_T) NULL;
  new->circular_ref_nmismatches = (Intlist_T) NULL;
  new->circular_univdiagonals = (Univcoordlist_T) NULL;
  new->circular_junctions = (List_T) NULL;
#endif

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->fusion_chrnum = -1;
  new->fusion_chroffset = 0;
  new->fusion_chrhigh = 0;

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* Determined by Path_solve procedures */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = true;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

#ifdef CHECK_ASSERTIONS
  prev_endpoint = Intlist_head(new->endpoints);
  for (q = Intlist_next(new->endpoints); q != NULL; q = Intlist_next(q)) {
    /* Previously used <=, but a deletion can yield duplicated endpoints */
    if (Intlist_head(q) < prev_endpoint) {
      printf("Path_create_from_transcript expected forward, but got\n");
      Path_print(new);
      abort();
    }
    prev_endpoint = Intlist_head(q);
  }
#endif

  debug0(Path_print(new));
  return new;
}



/* Called by combine_leftright_paths */
T
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
	     Pathpool_T pathpool, Vectorpool_T vectorpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
  Intlist_T q;
#ifdef CHECK_ASSERTIONS
  int prev_endpoint;
#endif

  assert(sensedir == SENSE_FORWARD || sensedir == SENSE_ANTI);

  assert(Univcoordlist_length(univdiagonals) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(nmismatches) == Intlist_length(endpoints) - 1);
  assert(Intlist_length(ref_nmismatches) == Intlist_length(endpoints) - 1);
  assert(List_length(junctions) == Intlist_length(endpoints) - 2);

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = querylength;
  new->score_within_trims = querylength;
  new->genomic_diff = (char *) NULL;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = sensedir;
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = main_univdiagonal;
    } else if (main_univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = main_univdiagonal - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = main_univdiagonal;
    } else if (main_univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = main_univdiagonal - querylength;
    }
  }

  debug0(printf("%d: Creating path %p, %u, by Path_create\n",
		++call_i,new,new->main_univdiagonal));
  
  new->endpoints = endpoints;
  new->univdiagonals = univdiagonals;
  new->nmismatches = nmismatches;
  new->ref_nmismatches = ref_nmismatches;
  new->junctions = junctions;

  new->splice5p = splice5p;
  new->splicetype5 = splicetype5;
  new->ambig_prob_5 = ambig_prob_5;

  new->splice3p = splice3p;
  new->splicetype3 = splicetype3;
  new->ambig_prob_3 = ambig_prob_3;

  /* Check that trimming has not invalidated altsplice */
  if (qstart_alts != NULL && qstart_alts->splice_qpos < Intlist_second_value(endpoints)) {
    new->qstart_alts = Altsplice_copy(qstart_alts,pathpool,vectorpool);
  } else {
    new->qstart_alts = (Altsplice_T) NULL;
  }
  if (qend_alts != NULL && qend_alts->splice_qpos > Intlist_penultimate_value(endpoints)) {
    new->qend_alts = Altsplice_copy(qend_alts,pathpool,vectorpool);
  } else {
    new->qend_alts = (Altsplice_T) NULL;
  }

  new->circular_endpoints = (Intlist_T) NULL;

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* Determined by caller */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = false;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

#ifdef CHECK_ASSERTIONS
  prev_endpoint = Intlist_head(new->endpoints);
  for (q = Intlist_next(new->endpoints); q != NULL; q = Intlist_next(q)) {
    if (Intlist_head(q) <= prev_endpoint) {
      printf("Path_create expected forward, but got\n");
      Path_print(new);
      abort();
    }
    prev_endpoint = Intlist_head(q);
  }
#endif

  return new;
}


T
Path_copy (T old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
	   Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));

  assert(old->sensedir == SENSE_FORWARD || old->sensedir == SENSE_ANTI);

  assert(Univcoordlist_length(old->univdiagonals) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->ref_nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(List_length(old->junctions) == Intlist_length(old->endpoints) - 2);

  new->query_compress_fwd = old->query_compress_fwd;
  new->query_compress_rev = old->query_compress_rev;

  new->nmatches = old->nmatches;
  new->ref_nmatches = old->ref_nmatches;
  new->junction_splice_prob = old->junction_splice_prob;
  new->total_splice_prob = old->total_splice_prob;
  new->found_score = old->found_score;
  new->score_within_trims = old->score_within_trims;

  if (old->genomic_diff == NULL) {
    new->genomic_diff = (char *) NULL;
  } else {
    new->genomic_diff = Pathpool_new_string(pathpool,old->querylength+1);
    strcpy(new->genomic_diff,old->genomic_diff);
  }

  new->plusp = old->plusp;
  new->genestrand = old->genestrand;
  new->sensedir = old->sensedir;
  new->querylength = old->querylength;

  new->chrnum = old->chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;

  new->main_univdiagonal = old->main_univdiagonal;

  debug0(printf("%d: Creating path %p, %u, by Path_copy from %p\n",
		++call_i,new,new->main_univdiagonal,old));

  new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
  new->univdiagonals = Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool);
  new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
  new->ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);
  new->junctions = Junction_copy_list(old->junctions,listpool,pathpool);
  
  new->splice5p = old->splice5p;
  new->splicetype5 = old->splicetype5;
  new->ambig_prob_5 = old->ambig_prob_5;

  new->splice3p = old->splice3p;
  new->splicetype3 = old->splicetype3;
  new->ambig_prob_3 = old->ambig_prob_3;

  new->qstart_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
  new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
  
  if (old->circular_endpoints == NULL) {
    new->circular_endpoints = (Intlist_T) NULL;
  } else {
    new->circular_high_p = old->circular_high_p;
    new->circular_endpoints = Intlistpool_copy(old->circular_endpoints,intlistpool);
    new->circular_univdiagonals = Univcoordlistpool_copy(old->circular_univdiagonals,univcoordlistpool);
    new->circular_nmismatches = Intlistpool_copy(old->circular_nmismatches,intlistpool);
    new->circular_ref_nmismatches = Intlistpool_copy(old->circular_ref_nmismatches,intlistpool);
    new->circular_junctions = Junction_copy_list(old->circular_junctions,listpool,pathpool);
  }

  new->fusion_querystart_junction = Junction_copy(old->fusion_querystart_junction,pathpool);
  new->fusion_queryend_junction = Junction_copy(old->fusion_queryend_junction,pathpool);
  if (new->fusion_querystart_junction == NULL && new->fusion_queryend_junction == NULL) {
    new->fusion_endpoints = (Intlist_T) NULL;
    /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */
  } else {
    new->fusion_chrnum = old->fusion_chrnum;
    new->fusion_chroffset = old->fusion_chroffset;
    new->fusion_chrhigh = old->fusion_chrhigh;
    new->fusion_plusp = old->fusion_plusp;
    
    new->fusion_endpoints = Intlistpool_copy(old->fusion_endpoints,intlistpool);
    new->fusion_univdiagonals = Univcoordlistpool_copy(old->fusion_univdiagonals,univcoordlistpool);
    new->fusion_nmismatches = Intlistpool_copy(old->fusion_nmismatches,intlistpool);
    new->fusion_ref_nmismatches = Intlistpool_copy(old->fusion_ref_nmismatches,intlistpool);
    new->fusion_junctions = Junction_copy_list(old->fusion_junctions,listpool,pathpool);

#if 0
    /* This computation was already done when constructing the fusion path */
    if (new->fusion_querystart_junction != NULL) {
      if (new->fusion_plusp == true) {
	new->fusion_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
      } else {
	new->fusion_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
      }
    } else {
      if (new->fusion_plusp == true) {
	new->fusion_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
      } else {
	new->fusion_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
      }
    }
#else
    new->fusion_alts = Altsplice_copy(old->fusion_alts,pathpool,vectorpool);
#endif
  }

  new->transcripts = Transcript_copy_list(old->transcripts,transcriptpool,listpool);
  new->invalid_transcripts = Transcript_copy_list(old->invalid_transcripts,transcriptpool,listpool);
  new->fusion_transcripts = Transcript_copy_list(old->fusion_transcripts,transcriptpool,listpool);
  new->fusion_invalid_transcripts = Transcript_copy_list(old->fusion_invalid_transcripts,transcriptpool,listpool);

  new->completep = old->completep;
  new->childp = old->childp;
  new->extendedp = old->extendedp;

  new->method = old->method;
  new->transcriptome_method_p = old->transcriptome_method_p;

  /* new->spliceends_qstart = old->spliceends_qstart; */
  /* new->spliceends_qend = old->spliceends_qend; */

  return new;
}


#if 0
List_T
Path_copy_list (List_T old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		Hitlistpool_T hitlistpool, Transcriptpool_T transcriptpool) {
  List_T new = NULL, p;

  for (p = old; p != NULL; p = List_next(p)) {
    new = Hitlist_push(new,hitlistpool,
		       (void *) Path_copy((T) List_head(p),intlistpool,univcoordlistpool,
					  listpool,pathpool,vectorpool,transcriptpool)
		       hitlistpool_trace(__FILE__,__LINE__));
  }
  return List_reverse(new);
}
#endif


T
Path_copy_5 (T old, bool splice5p, Splicetype_T splicetype5, double ambig_prob_5,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
  
  assert(old->sensedir == SENSE_FORWARD || old->sensedir == SENSE_ANTI);

  assert(Univcoordlist_length(old->univdiagonals) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->ref_nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(List_length(old->junctions) == Intlist_length(old->endpoints) - 2);

  new->query_compress_fwd = old->query_compress_fwd;
  new->query_compress_rev = old->query_compress_rev;

 /* Need to call Path_eval_nmatches on this copy */
  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = old->querylength;
  new->score_within_trims = old->querylength;

  if (old->genomic_diff == NULL) {
    new->genomic_diff = (char *) NULL;
  } else {
    new->genomic_diff = Pathpool_new_string(pathpool,old->querylength+1);
    strcpy(new->genomic_diff,old->genomic_diff);
  }

  new->plusp = old->plusp;
  new->genestrand = old->genestrand;
  new->sensedir = old->sensedir;
  new->querylength = old->querylength;

  new->chrnum = old->chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;

  new->main_univdiagonal = old->main_univdiagonal;

  debug0(printf("%d: Creating path %p, %u, by Path_copy_5 from %p\n",
		++call_i,new,new->main_univdiagonal,old));

  new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
  new->univdiagonals = Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool);
  new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
  new->ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);
  new->junctions = Junction_copy_list(old->junctions,listpool,pathpool);
  
  new->splice5p = splice5p;
  new->splicetype5 = splicetype5;
  new->ambig_prob_5 = ambig_prob_5;

  new->splice3p = old->splice3p;
  new->splicetype3 = old->splicetype3;
  new->ambig_prob_3 = old->ambig_prob_3;

  new->qstart_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
  new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
  
  assert(old->circular_endpoints == NULL);
  new->circular_endpoints = (Intlist_T) NULL;

  assert(old->fusion_querystart_junction == NULL);
  assert(old->fusion_queryend_junction == NULL);
  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = old->completep;
  new->childp = false;
  new->extendedp = false;

  new->method = old->method;
  new->transcriptome_method_p = old->transcriptome_method_p;

  /* new->spliceends_qstart = old->spliceends_qstart; */
  /* new->spliceends_qend = old->spliceends_qend; */

  return new;
}

T
Path_copy_3 (T old,  bool splice3p, Splicetype_T splicetype3, double ambig_prob_3,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
  
  assert(old->sensedir == SENSE_FORWARD || old->sensedir == SENSE_ANTI);

  assert(Univcoordlist_length(old->univdiagonals) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(Intlist_length(old->ref_nmismatches) == Intlist_length(old->endpoints) - 1);
  assert(List_length(old->junctions) == Intlist_length(old->endpoints) - 2);

  new->query_compress_fwd = old->query_compress_fwd;
  new->query_compress_rev = old->query_compress_rev;

 /* Need to call Path_eval_nmatches on this copy */
  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = old->querylength;
  new->score_within_trims = old->querylength;

  if (old->genomic_diff == NULL) {
    new->genomic_diff = (char *) NULL;
  } else {
    new->genomic_diff = Pathpool_new_string(pathpool,old->querylength+1);
    strcpy(new->genomic_diff,old->genomic_diff);
  }

  new->plusp = old->plusp;
  new->genestrand = old->genestrand;
  new->sensedir = old->sensedir;
  new->querylength = old->querylength;

  new->chrnum = old->chrnum;
  new->chroffset = old->chroffset;
  new->chrhigh = old->chrhigh;

  new->main_univdiagonal = old->main_univdiagonal;

  debug0(printf("%d: Creating path %p, %u, by Path_copy_3 from %p\n",
		++call_i,new,new->main_univdiagonal,old));

  new->endpoints = Intlistpool_copy(old->endpoints,intlistpool);
  new->univdiagonals = Univcoordlistpool_copy(old->univdiagonals,univcoordlistpool);
  new->nmismatches = Intlistpool_copy(old->nmismatches,intlistpool);
  new->ref_nmismatches = Intlistpool_copy(old->ref_nmismatches,intlistpool);
  new->junctions = Junction_copy_list(old->junctions,listpool,pathpool);
  
  new->splice5p = old->splice5p;
  new->splicetype5 = old->splicetype5;
  new->ambig_prob_5 = old->ambig_prob_5;

  new->splice3p = splice3p;
  new->splicetype3 = splicetype3;
  new->ambig_prob_3 = ambig_prob_3;

  new->qstart_alts = Altsplice_copy(old->qstart_alts,pathpool,vectorpool);
  new->qend_alts = Altsplice_copy(old->qend_alts,pathpool,vectorpool);
  
  assert(old->circular_endpoints == NULL);
  new->circular_endpoints = (Intlist_T) NULL;

  assert(old->fusion_querystart_junction == NULL);
  assert(old->fusion_queryend_junction == NULL);
  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = old->completep;
  new->childp = false;
  new->extendedp = false;

  new->method = old->method;
  new->transcriptome_method_p = old->transcriptome_method_p;

  /* new->spliceends_qstart = old->spliceends_qstart; */
  /* new->spliceends_qend = old->spliceends_qend; */

  return new;
}


int
Path_exon_origin (T this) {
  int exon_origin;
  Intlist_T p = this->endpoints;
  List_T j = this->junctions;

  p = Intlist_next(p);
  exon_origin = Intlist_head(p);

  while (j != NULL && Junction_type((Junction_T) List_head(j)) != SPLICE_JUNCTION) {
    p = Intlist_next(p);
    exon_origin = Intlist_head(p);

    j = List_next(j);
  }

  return exon_origin;
}


void
Path_count (int *npaths_primary, int *npaths_altloc, List_T paths) {
  T path;

  *npaths_primary = *npaths_altloc = 0;

  while (paths != NULL) {
    path = (T) List_head(paths);
    if (altlocp[path->chrnum] == true) {
      *npaths_altloc += 1;
    } else {
      *npaths_primary += 1;
    }
    paths = List_next(paths);
  }

  return;
}


void
Path_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	   Hitlistpool_T hitlistpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  debug0(printf("%d: Freeing path %p\n",++call_i,*old));

  assert((*old)->sensedir == SENSE_FORWARD || (*old)->sensedir == SENSE_ANTI);
  /* Note: assertions about lengths do not hold because procedures,
     such as Repair_path, can call Path_free before completing the
     path */

  /* Recursive calls */
  Intlistpool_free_list(&(*old)->endpoints,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Univcoordlistpool_free_list(&(*old)->univdiagonals,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
  Intlistpool_free_list(&(*old)->nmismatches,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Intlistpool_free_list(&(*old)->ref_nmismatches,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Junction_list_gc(&(*old)->junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */

  if ((*old)->circular_endpoints != NULL) {
    Intlistpool_free_list(&(*old)->circular_endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Univcoordlistpool_free_list(&(*old)->circular_univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
    Intlistpool_free_list(&(*old)->circular_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Intlistpool_free_list(&(*old)->circular_ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Junction_list_gc(&(*old)->circular_junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */
  }

  if ((*old)->fusion_endpoints != NULL) {
    Intlistpool_free_list(&(*old)->fusion_endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Univcoordlistpool_free_list(&(*old)->fusion_univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
    Intlistpool_free_list(&(*old)->fusion_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Intlistpool_free_list(&(*old)->fusion_ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Junction_list_gc(&(*old)->fusion_junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */

    if ((*old)->fusion_alts != NULL) {
      Altsplice_free(&(*old)->fusion_alts,pathpool);
    }
  }

  if ((*old)->fusion_querystart_junction != NULL) {
    Pathpool_free_junction(&(*old)->fusion_querystart_junction,pathpool
			   pathpool_trace(__FILE__,__LINE__));
  }
  if ((*old)->fusion_queryend_junction != NULL) {
    Pathpool_free_junction(&(*old)->fusion_queryend_junction,pathpool
			   pathpool_trace(__FILE__,__LINE__));
  }

  /* altsplices are allocated by Pathpool_new_altsplice */
  if ((*old)->qstart_alts != NULL) {
    Altsplice_free(&(*old)->qstart_alts,pathpool);
  }
  if ((*old)->qend_alts != NULL) {
    Altsplice_free(&(*old)->qend_alts,pathpool);
  }
  
#if 0
  /* genomic_diff is allocated by Pathpool_new_string */
  if ((*old)->genomic_diff != NULL) {
    FREE((*old)->genomic_diff);
  }
#endif

  Transcript_list_gc(&(*old)->transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->invalid_transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->fusion_transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->fusion_invalid_transcripts,listpool,transcriptpool);

  Pathpool_free_path(&(*old),pathpool
		     pathpool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */

  return;
}



void
Path_gc (List_T *list, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	 Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	 Hitlistpool_T hitlistpool) {
  List_T p;
  T old;
  
  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Path_free(&old,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }
  Hitlistpool_free_list(&(*list),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* allocated by Hitlistpool_push */
  return;
}


void
Path_array_gc (T *paths, int n, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	       Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	       Hitlistpool_T hitlistpool) {
  int i;
  T old;
  
  for (i = 0; i < n; i++) {
    old = paths[i];
    Path_free(&old,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
  }
  FREE(paths);

  return;
}



void
Path_free_parent (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  debug0(printf("%d: Freeing parent path %p\n",++call_i,*old));

  assert((*old)->sensedir == SENSE_FORWARD || (*old)->sensedir == SENSE_ANTI);
  /* Note: assertions about lengths do not hold because procedures,
     such as Repair_path, can call Path_free before completing the
     path */


  Intlistpool_free_list(&(*old)->endpoints,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Univcoordlistpool_free_list(&(*old)->univdiagonals,univcoordlistpool
			      univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
  Intlistpool_free_list(&(*old)->nmismatches,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Intlistpool_free_list(&(*old)->ref_nmismatches,intlistpool
			intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
  Junction_list_gc(&(*old)->junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */

  if ((*old)->circular_endpoints != NULL) {
    Intlistpool_free_list(&(*old)->circular_endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Univcoordlistpool_free_list(&(*old)->circular_univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
    Intlistpool_free_list(&(*old)->circular_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Intlistpool_free_list(&(*old)->circular_ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Junction_list_gc(&(*old)->circular_junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */
  }

  if ((*old)->fusion_endpoints != NULL) {
    Intlistpool_free_list(&(*old)->fusion_endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Univcoordlistpool_free_list(&(*old)->fusion_univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* allocated by Univcoordlistpool_push */
    Intlistpool_free_list(&(*old)->fusion_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Intlistpool_free_list(&(*old)->fusion_ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* allocated by Intlistpool_push */
    Junction_list_gc(&(*old)->fusion_junctions,listpool,pathpool); /* junctions allocated by Pathpool_new_junction, and list allocated by Listpool_push */

    if ((*old)->fusion_alts != NULL) {
      Altsplice_free(&(*old)->fusion_alts,pathpool);
    }
  }

  if ((*old)->fusion_querystart_junction != NULL) {
    Pathpool_free_junction(&(*old)->fusion_querystart_junction,pathpool
			   pathpool_trace(__FILE__,__LINE__));
  }
  if ((*old)->fusion_queryend_junction != NULL) {
    Pathpool_free_junction(&(*old)->fusion_queryend_junction,pathpool
			   pathpool_trace(__FILE__,__LINE__));
  }

  /* altsplices are allocated by Pathpool_new_altsplice */
  if ((*old)->qstart_alts != NULL) {
    Altsplice_free(&(*old)->qstart_alts,pathpool);
  }
  if ((*old)->qend_alts != NULL) {
    Altsplice_free(&(*old)->qend_alts,pathpool);
  }
  
#if 0
  /* genomic_diff is allocated by Pathpool_new_string */
  if ((*old)->genomic_diff != NULL) {
    FREE((*old)->genomic_diff);
  }
#endif

  Transcript_list_gc(&(*old)->transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->invalid_transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->fusion_transcripts,listpool,transcriptpool);
  Transcript_list_gc(&(*old)->fusion_invalid_transcripts,listpool,transcriptpool);

  Pathpool_free_path(&(*old),pathpool
		     pathpool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */

  return;
}


#if 0
static int
Path_nsegments (T this) {
  int nsegments;

  nsegments = Univcoordlist_length(this->univdiagonals);
  if (this->qstart_alts != NULL) {
    nsegments += 1;
  }
  if (this->qend_alts != NULL) {
    nsegments += 1;
  }
  return nsegments;
}
#endif


#if 0
static int
Path_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;
  Univcoordlist_T p, q;

  p = x->univdiagonals;
  q = y->univdiagonals;
  while (p != NULL && q != NULL) {

    if (Univcoordlist_head(p) < Univcoordlist_head(q)) {
      return -1;
    } else if (Univcoordlist_head(q) < Univcoordlist_head(p)) {
      return +1;
    } else {
      p = Univcoordlist_next(p);
      q = Univcoordlist_next(q);
    }
  }

  if (p == NULL && q == NULL) {
    return 0;
  } else if (p == NULL) {
    return -1;
  } else if (q == NULL) {
    return +1;
  } else {
    /* Not possible */
    return 0;
  }
}
#endif


T
Path_new_from_ends (Univcoord_T univdiagonal5, int qstart5, int qend5,
		    Univcoord_T univdiagonal3, int qstart3, int qend3,
		    bool plusp, bool first_read_p, int genestrand,
		    int sensedir, int querylength, Method_T method,
		    Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool) {
  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));

  assert(sensedir == SENSE_FORWARD || sensedir == SENSE_ANTI);

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = querylength;
  new->score_within_trims = querylength;
  new->genomic_diff = (char *) NULL;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = sensedir;
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  /* Use furthest ends, to allow cases where only the tails of each read overlap */
  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = univdiagonal3;
    } else if (univdiagonal5 < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal5 - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = univdiagonal3;
    } else if (univdiagonal5 < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal5 - querylength;
    }
  }

  debug0(printf("Creating path %p, %u, by Path_new_from_ends\n",
		new,new->main_univdiagonal));
  
  new->endpoints = Intlistpool_push(NULL,intlistpool,qend3
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,qstart3
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,qstart5
				    intlistpool_trace(__FILE__,__LINE__));

  new->univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal3
					      univcoordlistpool_trace(__FILE__,__LINE__));
  new->univdiagonals = Univcoordlistpool_push(new->univdiagonals,univcoordlistpool,univdiagonal5
					      univcoordlistpool_trace(__FILE__,__LINE__));

  new->nmismatches = Intlistpool_push(NULL,intlistpool,-1
				      intlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(new->nmismatches,intlistpool,-1
				      intlistpool_trace(__FILE__,__LINE__));

  new->ref_nmismatches = Intlistpool_push(NULL,intlistpool,-1
					  intlistpool_trace(__FILE__,__LINE__));
  new->ref_nmismatches = Intlistpool_push(new->ref_nmismatches,intlistpool,-1
					  intlistpool_trace(__FILE__,__LINE__));

#ifdef ALLOCATE_UNSOLVED_JUNCTION
  new->junctions = Listpool_push(NULL,listpool,Junction_new_unsolved(pathpool)
				 listpool_trace(__FILE__,__LINE__));
#else
  new->junctions = Listpool_push(NULL,listpool,(void *) JUNCTION_UNSOLVED
				 listpool_trace(__FILE__,__LINE__));
#endif
  
  new->splice5p = false;
  new->splicetype5 = NO_SPLICE;
  new->ambig_prob_5 = 0.0;

  new->splice3p = false;
  new->splicetype3 = NO_SPLICE;
  new->ambig_prob_3 = 0.0;

  new->qstart_alts = (Altsplice_T) NULL;
  new->qend_alts = (Altsplice_T) NULL;
  
  new->circular_endpoints = (Intlist_T) NULL;

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* Determined by Path_solve procedures */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = false;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

  return new;
}


T
Path_new_exact (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches, int ref_nmismatches,
		bool plusp, bool first_read_p, int genestrand, int querylength, int found_score,
		Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		Compress_T query_compress_fwd, Compress_T query_compress_rev,
		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool, Pathpool_T pathpool,
		Method_T method) {
#ifdef DEBUG0
  static int call_i = 0;
#endif

  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));

  assert(found_score >= 0);

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->nmatches = querylength - found_score;
  new->ref_nmatches = querylength - found_score;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = found_score;
  new->score_within_trims = nmismatches;
  new->genomic_diff = (char *) NULL;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = SENSE_FORWARD; /* antisense generated after Path_copy */
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  }

  debug0(printf("%d: Creating path %p, %u %d..%d nmismatches=%d by Path_new_exact\n",
		++call_i,new,new->main_univdiagonal,qstart,qend,nmismatches));
  
#if 0
  /* Should be handled by Path_trim_chrbound */
  if (univdiagonal + qend >= chrhigh + querylength) {
    new->endpoints = Intlistpool_push(NULL,intlistpool,univdiagonal - chrhigh
				      intlistpool_trace(__FILE__,__LINE__));
  } else {
    new->endpoints = Intlistpool_push(NULL,intlistpool,qend
				      intlistpool_trace(__FILE__,__LINE__));
  }
  if (univdiagonal + qstart < chroffset + querylength) {
    new->endpoints = Intlistpool_push(new->endpoints,intlistpool,chroffset + querylength - univdiagonal
				      intlistpool_trace(__FILE__,__LINE__));
  } else {
    new->endpoints = Intlistpool_push(new->endpoints,intlistpool,qstart
				      intlistpool_trace(__FILE__,__LINE__));
  }
#else
  new->endpoints = Intlistpool_push(NULL,intlistpool,qend
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,qstart
				    intlistpool_trace(__FILE__,__LINE__));
#endif

  new->univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));

  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__));
  new->ref_nmismatches = Intlistpool_push(NULL,intlistpool,ref_nmismatches
					  intlistpool_trace(__FILE__,__LINE__));

  new->junctions = (List_T) NULL;
  
  new->splice5p = false;
  new->splicetype5 = NO_SPLICE;
  new->ambig_prob_5 = 0.0;

  new->splice3p = false;
  new->splicetype3 = NO_SPLICE;
  new->ambig_prob_3 = 0.0;

  new->qstart_alts = (Altsplice_T) NULL;
  new->qend_alts = (Altsplice_T) NULL;
  
  new->circular_endpoints = (Intlist_T) NULL;

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* Determined by Path_solve procedures */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = false;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

  return new;
}



/* qstart and qend are the genome-normalized coordinates, so qstart
   marks the left coordinate and qend marks the right coordinate.  For
   a plus-strand alignment, qstart = querystart and qend = queryend.
   For a minus-strand alignment qstart = querylength - querystart and
   qend = querylength - queryend. */

/* Need to convert from qstart and qend to querystart and queryend
   when creating Altsplice_T objects */

T
Path_new_for_qstart_extension (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches,
			       bool plusp, bool first_read_p, int genestrand,
			       int sensedir, int querylength, Method_T method,
			       Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			       bool splice5p, Splicetype_T splicetype5, double ambig_prob_5, 
			       Compress_T query_compress_fwd, Compress_T query_compress_rev,
			       Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			       Pathpool_T pathpool) {
  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));
  
  assert(sensedir == SENSE_FORWARD || sensedir == SENSE_ANTI);

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = querylength;
  new->score_within_trims = querylength;
  new->genomic_diff = (char *) NULL;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = sensedir;
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  }

  debug0(printf("Creating path %p, %u, by Path_new_for_qstart_extension\n",new,new->main_univdiagonal));

  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,qend
						     intlistpool_trace(__FILE__,__LINE__)),intlistpool,qstart
				    intlistpool_trace(__FILE__,__LINE__));
  new->univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__)); /* qstart..qend found by consecutive matches */
  new->ref_nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
					  intlistpool_trace(__FILE__,__LINE__));
  new->junctions = (List_T) NULL;
  
  new->splice5p = splice5p;
  new->splicetype5 = splicetype5;
  new->ambig_prob_5 = ambig_prob_5;

  new->splice3p = false;
  new->splicetype3 = NO_SPLICE;
  new->ambig_prob_3 = 0.0;

  new->qstart_alts = (Altsplice_T) NULL;
  new->qend_alts = (Altsplice_T) NULL;
  
  new->circular_endpoints = (Intlist_T) NULL;

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* To be determined by Path_solve procedures */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = false;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

  return new;
}


T
Path_new_for_qend_extension (Univcoord_T univdiagonal, int qstart, int qend, int nmismatches,
			     bool plusp, bool first_read_p, int genestrand,
			     int sensedir, int querylength, Method_T method,
			     Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			     bool splice3p, Splicetype_T splicetype3, double ambig_prob_3, 
			     Compress_T query_compress_fwd, Compress_T query_compress_rev,
			     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
			     Pathpool_T pathpool) {
  T new = Pathpool_new_path(pathpool
			    pathpool_trace(__FILE__,__LINE__));

  assert(sensedir == SENSE_FORWARD || sensedir == SENSE_ANTI);

  new->query_compress_fwd = query_compress_fwd;
  new->query_compress_rev = query_compress_rev;

  new->nmatches = -1;
  new->ref_nmatches = -1;
  new->junction_splice_prob = 0.0;
  new->total_splice_prob = 0.0;
  new->found_score = querylength;
  new->score_within_trims = querylength;
  new->genomic_diff = (char *) NULL;

  new->plusp = plusp;
  new->genestrand = genestrand;
  new->sensedir = sensedir;
  new->querylength = querylength;

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  if (plusp == true) {
    if (first_read_p == false) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  } else {
    if (first_read_p == true) {
      new->main_univdiagonal = univdiagonal;
    } else if (univdiagonal < (Univcoord_T) querylength) {
      new->main_univdiagonal = 0;
    } else {
      new->main_univdiagonal = univdiagonal - querylength;
    }
  }

  debug0(printf("Creating path %p, %u, by Path_new_for_qend_extension\n",new,new->main_univdiagonal));

  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,qstart
						     intlistpool_trace(__FILE__,__LINE__)),intlistpool,qend
				    intlistpool_trace(__FILE__,__LINE__));
  new->univdiagonals = Univcoordlistpool_push(NULL,univcoordlistpool,univdiagonal
					      univcoordlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__));
  new->ref_nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
					  intlistpool_trace(__FILE__,__LINE__));
  new->junctions = (List_T) NULL;

  new->splice5p = false;
  new->splicetype5 = NO_SPLICE;
  new->ambig_prob_5 = 0.0;

  new->splice3p = splice3p;
  new->splicetype3 = splicetype3;
  new->ambig_prob_3 = ambig_prob_3;

  new->qstart_alts = (Altsplice_T) NULL;
  new->qend_alts = (Altsplice_T) NULL;

  new->circular_endpoints = (Intlist_T) NULL;

  new->fusion_querystart_junction = (Junction_T) NULL;
  new->fusion_queryend_junction = (Junction_T) NULL;
  new->fusion_endpoints = (Intlist_T) NULL;
  /* Obviates need to set fusion_univdiagonals, fusion_nmismatches, fusion_ref_nmismatches, fusion_junctions and fusion_alts */

  new->transcripts = (List_T) NULL;
  new->invalid_transcripts = (List_T) NULL;
  new->fusion_transcripts = (List_T) NULL;
  new->fusion_invalid_transcripts = (List_T) NULL;

  new->completep = false;	/* To be determined by Path_solve procedures */
  new->childp = false;
  new->extendedp = false;

  new->method = method;
  new->transcriptome_method_p = false;

  /* new->spliceends_qstart = -1; */
  /* new->spliceends_qend = -1; */

  return new;
}


#if 0
void
Path_fill_adj_genomiclows (Univcoord_T *positions, T *paths, int n, int querylength) {
  int i;

  for (i = 0; i < n; i++) {
    positions[i] = Path_genomiclow(paths[i]) + querylength;
  }

  return;
}
#endif

#if 0
void
Path_fill_adj_genomichighs (Univcoord_T *positions, T *paths, int n, int querylength) {
  int i;

  for (i = 0; i < n; i++) {
    positions[i] = Path_genomichigh(paths[i]) - querylength;
  }

  return;
}
#endif



static char *
sensedir_string (int sensedir) {
  if (sensedir == SENSE_NULL) {
    return "sense:null";
  } else if (sensedir == SENSE_ANTI) {
    return "sense:anti";
  } else if (sensedir == SENSE_FORWARD) {
    return "sense:fwd";
  } else {
    abort();
  }
}


static char *
splicetype_string (Splicetype_T splicetype) {
  switch (splicetype) {
  case DONOR: return "donor";
  case ACCEPTOR: return "acceptor";
  case ANTIDONOR: return "antidonor";
  case ANTIACCEPTOR: return "antiacceptor";
  default: abort();
  }
}


#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
void
Path_print (T this) {
  Junction_T fusion_junction;

  if (this != NULL) {
    printf(">> Path %p (%u): ",this,this->main_univdiagonal);
    printf("%s  ",Intlist_to_string(this->endpoints));
    printf("%s  ",Univcoordlist_to_string(this->univdiagonals));
    printf("%d:%s  ",this->chrnum,Univcoordlist_to_string_offset(this->univdiagonals,this->chroffset));

    if (this->nmatches >= 0) {
      /* Info added for new version */
      printf("found_score:%d, within_trims:%d; nmatches:%d splice_prob:%f->%f  ",
	     this->found_score,this->score_within_trims,this->nmatches,
	     this->junction_splice_prob,this->total_splice_prob);
    }

    /* Info added for new version */
    printf("plusp:%d try_%s eff_%s method:%s  ",
	   this->plusp,sensedir_string(this->sensedir),sensedir_string(Path_effective_sensedir(this)),
	   Method_string(this->method));

    printf("nmismatches:%s  ",Intlist_to_string(this->nmismatches));

    if (this->splice5p == false) {
      printf("splice5p:0 ");
    } else {
      printf("splice5p:1(%s,%f) ",splicetype_string(this->splicetype5),this->ambig_prob_5);
    }

    if (this->splice3p == false) {
      printf("splice3p:0  ");
    } else {
      printf("splice3p:1(%s,%f)  ",splicetype_string(this->splicetype3),this->ambig_prob_3);
    }

    printf("jcns:");
    Junction_print_list(this->junctions);
    printf("  qstart_alts:");
    Altsplice_print(this->qstart_alts);
    printf("  qend_alts:");
    Altsplice_print(this->qend_alts);

    printf("  transcripts:");
    Transcript_print_nums(this->transcripts);
    printf("  invalid:");
    Transcript_print_nums(this->invalid_transcripts);

    if (this->circular_endpoints != NULL) {
      printf(" CIRCULAR (high:%d): %s %s ",
	     this->circular_high_p,
	     Intlist_to_string(this->circular_endpoints),
	     Univcoordlist_to_string(this->circular_univdiagonals));
      Junction_print_list(this->circular_junctions);
    }

    if (this->fusion_querystart_junction != NULL || this->fusion_queryend_junction != NULL) {
      if (this->fusion_querystart_junction != NULL) {
	fusion_junction = this->fusion_querystart_junction;
	printf(" QUERYSTART_FUSION: plusp:%d,%c%c-%c%c,%f,%f ",
	       this->fusion_plusp,fusion_junction->donor1,fusion_junction->donor2,
	       fusion_junction->acceptor2,fusion_junction->acceptor1,
	       fusion_junction->donor_prob,fusion_junction->acceptor_prob);
      } else {
	fusion_junction = this->fusion_queryend_junction;
	printf(" QUERYEND_FUSION: plusp:%d,%c%c-%c%c,%f,%f ",
	       this->fusion_plusp,fusion_junction->donor1,fusion_junction->donor2,
	       fusion_junction->acceptor2,fusion_junction->acceptor1,
	       fusion_junction->donor_prob,fusion_junction->acceptor_prob);
      }
      printf("%s  %d:%s  %s  nmismatches:%s  ref_nmismatches:%s  ",
	     Univcoordlist_to_string(this->fusion_univdiagonals),this->fusion_chrnum,
	     Univcoordlist_to_string_offset(this->fusion_univdiagonals,this->fusion_chroffset),
	     Intlist_to_string(this->fusion_endpoints),
	     Intlist_to_string(this->fusion_nmismatches),Intlist_to_string(this->fusion_ref_nmismatches));
      printf("jcns:");
      Junction_print_list(this->fusion_junctions);
    }

    printf(" completep:%d",this->completep);

    if (this->genomic_diff != NULL) {
      printf(" %s",this->genomic_diff);
    }

    printf("\n");
  }

  return;
}
#endif


void
Path_check_valid (T this) {

  assert(Univcoordlist_length(this->univdiagonals) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(Intlist_length(this->ref_nmismatches) == Intlist_length(this->endpoints) - 1);
  assert(List_length(this->junctions) == Intlist_length(this->endpoints) - 2);

  return;
}


void
Path_setup (bool *circularp_in, bool *altlocp_in) {
  circularp = circularp_in;
  altlocp = altlocp_in;
  return;
}


