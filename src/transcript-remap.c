static char rcsid[] = "$Id: d20f49dc738fb3527d24196da4db83d8701ae09b $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "transcript-remap.h"

#include <stdio.h>
#include "assert.h"
#include "uintlist.h"
#include "junction.h"
#include "chrnum.h"
#include "transcript.h"
#include "repair.h"


/* Needed for speed when unspliced transcripts are possible */
#define INVALID_TRANSCRIPTS_SUFFICIENT 1

/* A 3' UTR extension is needed because poly-A tail is added to
   transcript, not represented in genome */
/* TODO: Look at A's or T's past UTR3 extension */
#define DISALLOW_UTR5_EXTENSION 1
/* #define DISALLOW_UTR3_EXTENSION 1 */


/* Transcript_remap_matchp */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Transcript_remap_all */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif

/* Transcript_remap_invalid */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif


/* binary search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif



static bool *circularp;
static Transcriptome_T transcriptome;
static IIT_T transcript_map_iit;
static int *transcript_chrnum_crosstable;


#define T Transcript_T


static inline Chrpos_T
exonend_geneplus (int exoni, int *exonbounds, Chrpos_T *exonstarts) {
  if (exoni == 0) {
    return exonstarts[exoni] + /*exonlength*/exonbounds[0] - 1;
  } else {
    return exonstarts[exoni] + /*exonlength*/(exonbounds[exoni] - exonbounds[exoni - 1]) - 1;
  }
}

/* exonstarts are in ascending order */
static int
binary_search_geneplus (int lowi, int highi, Chrpos_T *exonstarts, Univcoord_T chr_alignstart) {
  int middlei;
#ifdef CHECK_ASSERTIONS
  int n = highi;
#endif

  debug20(printf("entered binary search_geneplus with lowi=%d, highi=%d, alignstart=%u\n",
		 lowi,highi,chr_alignstart));

  while (lowi < highi) {
    /* middlei = lowi + ((highi - lowi) / 2); */
    middlei = (lowi + highi)/2;
    assert(middlei < n);

    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,exonstarts[lowi],middlei,exonstarts[middlei],highi,exonstarts[highi],
		   chr_alignstart));
    if (chr_alignstart < exonstarts[middlei]) {
      highi = middlei;
    } else if (chr_alignstart > exonstarts[middlei]) {
      lowi = middlei + 1;
    } else {
      debug20(printf("(1) Equality with alignstart: binary_search_geneplus returns middlei %d\n",middlei));
      return middlei;
    }
  }

  debug20(printf("(2) binary_search_geneplus returns highi %d - 1\n",highi));
  return highi - 1;
}



static inline Chrpos_T
exonend_geneminus (int exoni, int *exonbounds, Chrpos_T *exonstarts) {
  if (exoni == 0) {
    return exonstarts[exoni] - /*exonlength*/exonbounds[0] + 1;
  } else {
    return exonstarts[exoni] - /*exonlength*/(exonbounds[exoni] - exonbounds[exoni - 1]) + 1;
  }
}

/* exonstarts are in descending order */
static int
binary_search_geneminus (int lowi, int highi, Chrpos_T *exonstarts, Univcoord_T chr_alignstart) {
  int middlei;
#ifdef CHECK_ASSERTIONS
  int n = highi;
#endif

  debug20(printf("entered binary search_geneminus with lowi=%d, highi=%d, alignstart=%u\n",
		 lowi,highi,chr_alignstart));

  while (lowi < highi) {
    /* middlei = lowi + ((highi - lowi) / 2); */
    middlei = (lowi + highi)/2;
    assert(middlei < n);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,exonstarts[lowi],middlei,exonstarts[middlei],highi,exonstarts[highi],
		   chr_alignstart));
    if (chr_alignstart > exonstarts[middlei]) {
      highi = middlei;
    } else if (chr_alignstart < exonstarts[middlei]) {
      lowi = middlei + 1;
    } else {
      debug20(printf("(1) Equality with alignstart: binary_search_geneminus returns middlei %d\n",middlei));
      return middlei;
    }
  }

  debug20(printf("(2) binary_search_geneminus returns highi %d - 1\n",highi));
  return highi - 1;
}


/* For geneminus, endpoints/univdiagonals/junctions need to be in reverse genomic order */
Uintlist_T
Transcript_remap_endpoint_coords (Intlist_T endpoints, Univcoordlist_T univdiagonals, List_T junctions,
				  Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh, int querylength,
				  Uintlistpool_T uintlistpool, bool extend_qstart_p, bool extend_qend_p) {
  Uintlist_T coords = NULL;
  Chrpos_T chr_alignstart, chr_alignend, chr_left;
  Univcoord_T univdiagonal, chrmidpoint;
  Junction_T junction;
  int type;
  int endpoint, qstart, qend;

  debug10(printf("Entering Transcript_remap_endpoint_coords with univdiagonals %s (%s)\n",
		 Univcoordlist_to_string(univdiagonals),Univcoordlist_to_string_offset(univdiagonals,chroffset)));

  /* Start first exon */
  univdiagonal = Univcoordlist_head(univdiagonals);
  if (circularp[chrnum] == true &&
      univdiagonal > (chrmidpoint = chroffset + (chrhigh - chroffset)/2) + querylength) {
    chr_left = (univdiagonal - chrmidpoint) - querylength;
  } else {    
    chr_left = (univdiagonal - chroffset) - querylength;
  }
  debug10(printf("First exon: chr_left is %u\n",chr_left));
  
  qstart = extend_qstart_p == true ? 0: Intlist_head(endpoints);
  chr_alignstart = chr_left + qstart + 1; /* 1-based */

  endpoints = Intlist_next(endpoints);
  univdiagonals = Univcoordlist_next(univdiagonals);

  while (univdiagonals != NULL) {
    junction = (Junction_T) List_head(junctions);
    if ((type = Junction_type(junction)) == INS_JUNCTION || type == DEL_JUNCTION) {
      /* Update chr_left for next iteration */
      univdiagonal = Univcoordlist_head(univdiagonals);
      if (circularp[chrnum] == true && univdiagonal > chrmidpoint + querylength) {
	chr_left = (univdiagonal - chrmidpoint) - querylength;
      } else {
	chr_left = (univdiagonal - chroffset) - querylength;
      }

    } else {
      /* Emit exon */
      endpoint = Intlist_head(endpoints);
      chr_alignend = chr_left + endpoint; /* Use previous value of chr_left; 1-based */
      coords = Uintlistpool_push(coords,uintlistpool,chr_alignstart
				 uintlistpool_trace(__FILE__,__LINE__));
      coords = Uintlistpool_push(coords,uintlistpool,chr_alignend
				  uintlistpool_trace(__FILE__,__LINE__));

      /* Start a new exon.  Revise chr_alignstart */
      univdiagonal = Univcoordlist_head(univdiagonals);
      if (circularp[chrnum] == true && univdiagonal > chrmidpoint + querylength) {
	chr_left = (univdiagonal - chrmidpoint) - querylength;
      } else {
	chr_left = (univdiagonal - chroffset) - querylength;
      }
      chr_alignstart = chr_left + endpoint + 1; /* 1-based */
    }

    endpoints = Intlist_next(endpoints);
    univdiagonals = Univcoordlist_next(univdiagonals);
    junctions = List_next(junctions);
  }

  /* Emit last exon */
  /* Intlist_head(endpoints) should be the same as Intlist_last_value(endpoints) */
  qend = extend_qend_p == true ? querylength : Intlist_head(endpoints);
  chr_alignend = chr_left + qend; /* 1-based */

  coords = Uintlistpool_push(coords,uintlistpool,chr_alignstart
			      uintlistpool_trace(__FILE__,__LINE__));
  coords = Uintlistpool_push(coords,uintlistpool,chr_alignend
			      uintlistpool_trace(__FILE__,__LINE__));

#ifdef DEBUG10
  coords = Uintlist_reverse(coords);
  debug10(printf("Returning coords %s\n",Uintlist_to_string(coords)));
  return coords;
#else
  return Uintlist_reverse(coords);
#endif
}
	

#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
static Chrpos_T *
compute_exonends_geneplus (int *exonbounds, Chrpos_T *exonstarts, int transcript_nexons,
			   Vectorpool_T vectorpool) {
  Chrpos_T *exonends;
  int last_bound, exonlength;
  int exoni;

  exonends = Vectorpool_new_uintvector(vectorpool,transcript_nexons);

  last_bound = 0;
  for (exoni = 0; exoni < transcript_nexons; exoni++) {
    exonlength = exonbounds[exoni] - last_bound;
    exonends[exoni] = exonstarts[exoni] + exonlength - 1;
    last_bound = exonbounds[exoni];
  }

  return exonends;
}
#endif

#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
static Chrpos_T *
compute_exonends_geneminus (int *exonbounds, Chrpos_T *exonstarts, int transcript_nexons,
			    Vectorpool_T vectorpool) {
  Chrpos_T *exonends;
  int last_bound, exonlength;
  int exoni;

  exonends = Vectorpool_new_uintvector(vectorpool,transcript_nexons);

  last_bound = 0;
  for (exoni = 0; exoni < transcript_nexons; exoni++) {
    exonlength = exonbounds[exoni] - last_bound;
    exonends[exoni] = exonstarts[exoni] - exonlength + 1;
    last_bound = exonbounds[exoni];
  }

  return exonends;
}
#endif


static void
bound_segment_geneplus (int *exoni, int *exonj, Chrpos_T chr_alignstart, Chrpos_T chr_alignend,
			int *exonbounds, Chrpos_T *exonstarts,
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
			Chrpos_T *exonends,
#endif
			int transcript_nexons) {
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  int i;
#endif

  debug10(printf("bound_segment_geneplus with chr_alignstart %u, chr_alignend %u\n",
		 chr_alignstart,chr_alignend));

  if (chr_alignstart <= exonstarts[0]) {
    *exoni = -1;
    *exonj = 0;
  } else {
    *exoni = binary_search_geneplus(0,transcript_nexons,exonstarts,chr_alignstart);
    *exonj = *exoni;
    assert(chr_alignstart >= exonstarts[*exoni]);
#ifdef CHECK_ASSERTIONS
    if (*exoni < transcript_nexons - 1) {
      assert(chr_alignstart < exonstarts[(*exoni)+1]);
    }
#endif
  }
  /* chr_alignstart could be anywhere from the start of exoni to the start of (exoni + 1) */
  debug20(printf("Got exoni %d\n",*exoni));

  while (*exonj < transcript_nexons &&
	 chr_alignend >= exonend_geneplus(*exonj,exonbounds,exonstarts)) {
    (*exonj)++;
    debug20(printf("Advanced exonj to %d\n",*exonj));
  }
  debug20(printf("Comparing chr_alignend %u against %u\n",chr_alignend,exonend_geneplus(*exonj,exonbounds,exonstarts)));
#ifdef CHECK_ASSERTIONS
  if (*exonj < transcript_nexons - 1) {
    assert(chr_alignend <= exonend_geneplus(*exonj,exonbounds,exonstarts));
  }
#endif

  debug10(printf("bound_segment_geneplus returning exoni %d, exonj %d\n\n",*exoni,*exonj));

  return;
}


static void
bound_segment_geneminus (int *exoni, int *exonj, Chrpos_T chr_alignstart, Chrpos_T chr_alignend,
			 int *exonbounds, Chrpos_T *exonstarts,
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
			 Chrpos_T *exonends,
#endif
			 int transcript_nexons) {
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  int i;
#endif

  debug10(printf("bound_segment_geneminus with chr_alignstart %u, chr_alignend %u\n",
		 chr_alignstart,chr_alignend));

  if (chr_alignstart >= exonstarts[0]) {
    *exoni = -1;
    *exonj = 0;
  } else {
    *exoni = binary_search_geneminus(0,transcript_nexons,exonstarts,chr_alignstart);
    *exonj = *exoni;
    assert(chr_alignstart <= exonstarts[*exoni]);
#ifdef CHECK_ASSERTIONS
    if (*exoni < transcript_nexons - 1) {
      assert(chr_alignstart > exonstarts[(*exoni)+1]);
    }
#endif
  }
  /* chr_alignstart could be anywhere from the start of exoni to the start of (exoni + 1) */
  debug20(printf("Got exoni %d\n",*exoni));

  while (*exonj < transcript_nexons &&
	 chr_alignend <= exonend_geneminus(*exonj,exonbounds,exonstarts)) {
    (*exonj)++;
    debug20(printf("Advanced exonj to %d\n",*exonj));
  }
  debug20(printf("Comparing chr_alignend %u against %u\n",chr_alignend,exonend_geneminus(*exonj,exonbounds,exonstarts)));
#ifdef CHECK_ASSERTIONS
  if (*exonj < transcript_nexons - 1) {
    assert(chr_alignend >= exonend_geneminus(*exonj,exonbounds,exonstarts));
  }
#endif

  debug10(printf("bound_segment_geneminus returning exoni %d, exonj %d\n\n",*exoni,*exonj));

  return;
}


List_T
Transcript_remap_geneplus (int *overall_trstart, int *overall_trend,
			   int *trstart_overhang, Chrpos_T *trstart_splice_distance,
			   int *trend_overhang, Chrpos_T *trend_splice_distance,
			   Shortread_T queryseq, bool plusp,
			   Uintlist_T coords, int *exonbounds, Chrpos_T *exonstarts,
			   int transcript_nexons, Transcriptpool_T transcriptpool, Listpool_T listpool) {
  List_T exons = NULL;
  Uintlist_T p = coords;
  char firstchar, lastchar;
  Chrpos_T chr_alignstart, chr_alignend;
  Chrpos_T exonend_exoni, exonend_exonj, exonend_prev;
  int exoni, exonj, i;
  int trstart, trlength;
  bool intronp = false, utrp = false, inconsistent_polya = false;
  int choplength;
  
  trlength = exonbounds[transcript_nexons-1];

  *overall_trstart = -1;
  *trstart_overhang = *trend_overhang = 0;
  *trstart_splice_distance = *trend_splice_distance = 0;

  if (plusp == true) {
    choplength = Shortread_right_choplength(queryseq);
  } else {
    choplength = Shortread_left_choplength(queryseq);
  }
  debug10(printf("Transcript_remap_geneplus: plusp %d, choplength %d\n",plusp,choplength));

#ifdef DEBUG10
  printf("exon 0: exonbound %d, exonstart %u, exonend %u\n",
	 exonbounds[0],exonstarts[0],exonend_geneplus(0,exonbounds,exonstarts));
  for (exoni = 1; exoni < transcript_nexons; exoni++) {
    printf("intron: %u\n",exonstarts[exoni] - exonend_geneplus(exoni-1,exonbounds,exonstarts) - 1);
    printf("exon %d: exonbound %d, exonstart %u, exonend %u\n",
	   exoni,exonbounds[exoni],exonstarts[exoni],exonend_geneplus(exoni,exonbounds,exonstarts));
  }
#endif

  while (p != NULL) {
    /* Map one segment */
    chr_alignstart = Uintlist_head(p);
    chr_alignend = Uintlist_head(Uintlist_next(p));

    bound_segment_geneplus(&exoni,&exonj,chr_alignstart,chr_alignend,
			   exonbounds,exonstarts,transcript_nexons);
    debug10(printf("bounds geneplus: exoni %d, exonj %d\n",exoni,exonj));

    if (exoni < 0) {
      firstchar = 'u';
      trstart = 1;
      debug10(printf("Start case 0: Segment starts before transcript\n"));

    } else {
      exonend_exoni = exonend_geneplus(exoni,exonbounds,exonstarts);
      debug10(printf("First char: chr_alignstart %u versus exonstarts %u, exoni %d vs transcript_nexons %d\n",
		     chr_alignstart,exonstarts[exoni],exoni,transcript_nexons));

      if (chr_alignstart == exonstarts[exoni]) {
	debug10(printf("Start case 1: Segment starts at beginning of exon\n"));
	trstart = exonbounds[exoni] - (/*exonends[exoni]*/exonend_exoni - chr_alignstart);
	firstchar = (exoni == 0) ? '.' : 's';

      } else if (chr_alignstart <= exonend_exoni) {
	debug10(printf("Start case 2: Segment starts within exon, p p%, coords %p\n",p,coords));
	trstart = exonbounds[exoni] - (/*exonends[exoni]*/exonend_exoni - chr_alignstart);
	firstchar = (p == coords) ? '.' : 'y';
      
      } else if (exoni >= transcript_nexons) {
	debug10(printf("Start case N: Segment starts after transcript\n"));
	trstart = trlength + 1;
	firstchar = 'u';

      } else if (chr_alignend < exonstarts[exoni+1]) {
	debug10(printf("Start case 3: Segment is entirely within intron\n"));
	trstart = exonbounds[exoni] + 1; /* First bp of next exon */
	firstchar = 'i';
      
      } else {
	trstart = exonbounds[exoni] + 1 /* First bp of next exon */;
	*trstart_overhang = exonstarts[exoni+1] - chr_alignstart;
	*trstart_splice_distance = exonstarts[exoni+1] - exonend_exoni - 1;
	debug10(printf("Start case 4: Segment straddles intron and next exon: trstart_overhang %d = %u - %u, trstart_splice_distance %u\n",
		       *trstart_overhang,exonstarts[exoni+1],chr_alignstart,*trstart_splice_distance));
	firstchar = 'x';
      }
    }

    if (*overall_trstart >= 0) {
      /* Skip.  Already computed */
    } else {
      *overall_trstart = trstart;
      debug10(printf("Computing trstart as %d\n",*overall_trstart));
    }

    i = exoni;
    while (i < exonj) {
      debug10(printf("Emitting exon: %c%d%c\n",firstchar,i+1,'i'));
      exons = Listpool_push(exons,listpool,
			    Exon_new(firstchar,exoni,'i',transcriptpool)
			    listpool_trace(__FILE__,__LINE__));
      if (firstchar == 'i') {
	intronp = true;
      }
      
      firstchar = 'i';
      i++;
    }

    p = Uintlist_next(Uintlist_next(p));

    if (exonj == transcript_nexons) {
      lastchar = 'u';
      *overall_trend = trlength;
      debug10(printf("End case 0: Segment ends after transcript\n"));

    } else {
      exonend_exonj = exonend_geneplus(exonj,exonbounds,exonstarts);
      debug10(printf("Last char: chr_alignend %u versus exonends %u, exonj %d vs transcript_nexons %d\n",
		     chr_alignend,exonend_exonj,exonj,transcript_nexons));

      if (chr_alignend == /*exonends[exonj]*/exonend_exonj) {
	*overall_trend = exonbounds[exonj];
	debug10(printf("End case 1: Segment ends at end of exon. trend is %d\n",*overall_trend));
	if (p != NULL) {
	  lastchar = (exonj == transcript_nexons - 1) ? '.' : 's';
	} else if (choplength == 0) {
	  lastchar = (exonj == transcript_nexons - 1) ? '.' : 's';
	} else if (exonj == transcript_nexons - 1) {
	  lastchar = '.';
	} else {
	  lastchar = 's';
	  inconsistent_polya = true;
	}

      } else if (chr_alignend >= exonstarts[exonj]) {
	*overall_trend = exonbounds[exonj] - (/*exonends[exonj]*/exonend_exonj - chr_alignend);
	debug10(printf("End case 2: Segment ends within exon: trend is %d\n",*overall_trend));
	if (p != NULL) {
	  lastchar = 'y';
	} else if (choplength == 0) {
	  lastchar = '.';
	} else {
	  lastchar = '.';
	  inconsistent_polya = true;
	}

      } else if (exonj == 0) {
        *overall_trend = 0;
	debug10(printf("End case N: Segment ends before transcript\n"));
	lastchar = 'u';

      } else if (chr_alignstart > (exonend_prev = exonend_geneplus(exonj-1,exonbounds,exonstarts))) {
	*overall_trend = exonbounds[exonj-1]; /* Last bp of previous exon */
	debug10(printf("End case 3: Segment is entirely within intron\n"));
	lastchar = 'i';

      } else {
	*overall_trend = exonbounds[exonj-1]; /* Last bp of previous exon */
	*trend_overhang = chr_alignend - exonend_prev;
	*trend_splice_distance = exonstarts[exonj] - exonend_prev - 1;
	debug10(printf("End case 4: Segment straddles intron and previous exon: trend_overhang %d = %u - %u, trend_splice_distance %u\n",
		       *trend_overhang,chr_alignend,exonend_prev,*trend_splice_distance));
	lastchar = 'x';
      }
    }

    if (lastchar == 'i') {
      /* Skip.  Already generated intron previously */
    } else {
      debug10(printf("Emitting final exon: %c%d%c\n",firstchar,i+1,lastchar));
      exons = Listpool_push(exons,listpool,
			    Exon_new(firstchar,exoni,lastchar,transcriptpool)
			    listpool_trace(__FILE__,__LINE__));
    }
  }

#if 0
  if (extend_qend_p == true && *trend_overhang == 0) {
    *overall_trend -= querylength - Intlist_last_value(endpoints);
  }
#endif

  if (*overall_trstart >= *overall_trend) {
    debug10(printf("Endpoints %d..%d don't make sense\n",*overall_trstart,*overall_trend));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (*overall_trstart <= 0) {
    debug10(printf("Begins before transcript\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (*overall_trend > trlength) {
    /* No overlap with any exon or intron */
    debug10(printf("Ends after transcript\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (Exon_list_consecutivep(exons) == false) {
    /* Not consistent with transcript */
    debug10(printf("Non-consecutive exons\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (inconsistent_polya == true) {
    /* Not consistent with transcript */
    debug10(printf("inconsistent polya\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (utrp == true) {
    /* Not consistent with transcript */
    debug10(printf("UTR\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (intronp == true && List_length(exons) > 1) {
    /* Intron not allowed with other exons */
    debug10(printf("Intron with other exons\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else {
    return List_reverse(exons);
  }
}


List_T
Transcript_remap_geneminus (int *overall_trstart, int *overall_trend,
			    int *trstart_overhang, Chrpos_T *trstart_splice_distance,
			    int *trend_overhang, Chrpos_T *trend_splice_distance,
			    Shortread_T queryseq, bool plusp,
			    Uintlist_T coords, int *exonbounds, Chrpos_T *exonstarts,
			    int transcript_nexons, Transcriptpool_T transcriptpool, Listpool_T listpool) {
  List_T exons = NULL;
  Uintlist_T p = coords;
  char firstchar, lastchar;
  Chrpos_T chr_alignstart, chr_alignend;
  Chrpos_T exonend_exoni, exonend_exonj, exonend_prev;
  int exoni, exonj, i;
  int trstart, trlength;
  bool intronp = false, utrp = false, inconsistent_polya = false;
  int choplength;

  trlength = exonbounds[transcript_nexons-1];

  *overall_trstart = -1;
  *trstart_overhang = *trend_overhang = 0;
  *trstart_splice_distance = *trend_splice_distance = 0;

  if (plusp == true) {
    choplength = Shortread_left_choplength(queryseq);
  } else {
    choplength = Shortread_right_choplength(queryseq);
  }
  debug10(printf("Transcript_remap_geneminus: plusp %d, choplength %d\n",plusp,choplength));

#ifdef DEBUG10
    printf("exon 0: exonbound %d, exonstart %u, exonend %u\n",
	   exonbounds[0],exonstarts[0],exonend_geneminus(0,exonbounds,exonstarts));
    for (exoni = 1; exoni < transcript_nexons; exoni++) {
      printf("intron: %u\n",exonend_geneminus(exoni - 1,exonbounds,exonstarts) - exonstarts[exoni] - 1);
      printf("exon %d: exonbound %d, exonstart %u, exonend %u\n",
	     exoni,exonbounds[exoni],exonstarts[exoni],exonend_geneminus(exoni,exonbounds,exonstarts));
    }
#endif

  while (p != NULL) {
    /* Map one segment */
    chr_alignstart = Uintlist_head(p);
    chr_alignend = Uintlist_head(Uintlist_next(p));

    bound_segment_geneminus(&exoni,&exonj,chr_alignstart,chr_alignend,
			    exonbounds,exonstarts,transcript_nexons);
    debug10(printf("bounds geneminus: exoni %d, exonj %d\n",exoni,exonj));

    debug10(printf("First char: chr_alignstart %u versus exonstarts %u, exoni %d vs transcript_nexons %d\n",
		   chr_alignstart,exonstarts[exoni],exoni,transcript_nexons));
    if (exoni < 0) {
      firstchar = 'u';
      trstart = 1;
      debug10(printf("Start case 0: Segment starts before transcript\n"));

    } else {
      exonend_exoni = exonend_geneminus(exoni,exonbounds,exonstarts);
      if (chr_alignstart == exonstarts[exoni]) {
	debug10(printf("Start case 1: Segment starts at beginning of exon\n"));
	trstart = exonbounds[exoni] - (chr_alignstart - /*exonends[exoni]*/exonend_exoni);
	firstchar = (exoni == 0) ? '.' : 's';
      
      } else if (chr_alignstart >= exonend_exoni) {
	debug10(printf("Start case 2: Segment starts within exon\n"));
	trstart = exonbounds[exoni] - (chr_alignstart - /*exonends[exoni]*/exonend_exoni);
	firstchar = (p == coords) ? '.' : 'y';
      
      } else if (exoni >= transcript_nexons) {
	debug10(printf("Start case N: Segment starts after transcript\n"));
	trstart = trlength + 1;
	firstchar = 'u';

      } else if (chr_alignend > exonstarts[exoni+1]) {
	debug10(printf("Start case 3: Segment is entirely within intron\n"));
	trstart = exonbounds[exoni] + 1; /* First bp of next exon */
	firstchar = 'i';

      } else {
	trstart = exonbounds[exoni] + 1 /* First bp of next exon */;
	*trstart_overhang = chr_alignstart - exonstarts[exoni+1];
	*trstart_splice_distance = exonend_exoni - exonstarts[exoni+1] - 1;
	debug10(printf("Start case 4: Segment straddles intron and next exon: trstart_overhang %d = %u - %u, trstart_splice_distance %u\n",
		       *trstart_overhang,chr_alignstart,exonstarts[exoni+1],*trstart_splice_distance));
	firstchar = 'x';
      }
    }
    
    if (*overall_trstart >= 0) {
      /* Skip */
    } else {
      *overall_trstart = trstart;
      debug10(printf("Computing trstart as %d\n",*overall_trstart));
    }
    
    i = exoni;
    while (i < exonj) {
      debug10(printf("Emitting exon: %c%d%c\n",firstchar,i+1,'i'));
      exons = Listpool_push(exons,listpool,
			    Exon_new(firstchar,exoni,'i',transcriptpool)
			    listpool_trace(__FILE__,__LINE__));
      if (firstchar == 'i') {
	intronp = true;
      }
      
      firstchar = 'i';
      i++;
    }
    
    p = Uintlist_next(Uintlist_next(p));
      
    if (exonj == transcript_nexons) {
      lastchar = 'u';
      *overall_trend = trlength;
      debug10(printf("End case 0: Segment ends after transcript\n"));

    } else {
      exonend_exonj = exonend_geneminus(exonj,exonbounds,exonstarts);
      debug10(printf("Last char: chr_alignend %u versus exonend %u, exonj %d vs transcript_nexons %d\n",
		     chr_alignend,exonend_exonj,exonj,transcript_nexons));

      if (chr_alignend == /*exonends[exonj]*/exonend_exonj) {
	*overall_trend = exonbounds[exonj];
	debug10(printf("End case 1: Segment ends at end of exon. trend is %d\n",*overall_trend));
	if (p != NULL) {
	  lastchar = (exonj == transcript_nexons - 1) ? '.' : 's';
	} else if (choplength == 0) {
	  lastchar = (exonj == transcript_nexons - 1) ? '.' : 's';
	} else if (exonj == transcript_nexons - 1) {
	  lastchar = '.';
	} else {
	  lastchar = 's';
	  inconsistent_polya = true;
	}

      } else if (chr_alignend <= exonstarts[exonj]) {
	*overall_trend = exonbounds[exonj] - (chr_alignend - /*exonends[exonj]*/exonend_exonj);
	debug10(printf("End case 2: Segment ends within exon. trend is %d\n",*overall_trend));
	if (p != NULL) {
	  lastchar = 'y';
	} else if (choplength == 0) {
	  lastchar = '.';
	} else {
	  lastchar = '.';
	  inconsistent_polya = true;
	}
      
      } else if (exonj == 0) {
        *overall_trend = 0;
	debug10(printf("End case N: Segment ends before transcript\n"));
	lastchar = 'u';

      } else if (chr_alignstart < (exonend_prev = exonend_geneminus(exonj-1,exonbounds,exonstarts))) {
	*overall_trend = exonbounds[exonj];
	debug10(printf("End case 3: Segment is entirely within intron\n"));
	lastchar = 'i';
      
      } else {
	*overall_trend = exonbounds[exonj-1]; /* Last bp of previous exon */
	*trend_overhang = /*exonends[exonj]*/exonend_prev - chr_alignend;
	*trend_splice_distance = exonend_prev - exonstarts[exonj] - 1;
	debug10(printf("End case 4: Segment straddles intron and previous exon: trend_overhang %d = %u - %u, trend_splice_distance %u\n",
		       *trend_overhang,chr_alignend,exonend_prev,*trend_splice_distance));
	lastchar = 'x';
      }
    }
    
    if (lastchar == 'i') {
      /* Skip.  Already generated intron previously */
    } else {
      debug10(printf("Emitting final_exon: %c%d%c\n",firstchar,i+1,lastchar));
      exons = Listpool_push(exons,listpool,
			    Exon_new(firstchar,exoni,lastchar,transcriptpool)
			    listpool_trace(__FILE__,__LINE__));
    }
  }

#if 0
  if (extrend_qstart_p == true && *trend_overhang == 0) {
    *overall_trend -= Intlist_head(endpoints);
  }
#endif

  if (*overall_trstart >= *overall_trend) {
    debug10(printf("Endpoints %d..%d don't make sense\n",*overall_trstart,*overall_trend));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (*overall_trstart <= 0) {
    debug10(printf("Begins before transcript\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (*overall_trend > trlength) {
    debug10(printf("Ends after transcript\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (Exon_list_consecutivep(exons) == false) {
    /* Not consistent with transcript */
    debug10(printf("Non-consecutive exons\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (inconsistent_polya == true) {
    /* Not consistent with transcript */
    debug10(printf("inconsistent polya\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (utrp == true) {
    /* Not consistent with transcript */
    debug10(printf("UTR\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else if (intronp == true && List_length(exons) > 1) {
    /* Intron not allowed with other exons */
    debug10(printf("Intron with other exons\n\n"));
    Exon_list_gc(&exons,listpool,transcriptpool);
    return (List_T) NULL;

  } else {
    return List_reverse(exons);
  }
}


bool
Transcript_remap_invalid (Transcript_T transcript, Path_T path, Transcriptome_T transcriptome,
			  Shortread_T queryseq, Uintlistpool_T uintlistpool, Listpool_T listpool,
			  Transcriptpool_T transcriptpool) {
  bool assignedp;

  int transcript_nexons;
  int *exonbounds;
  Chrpos_T *exonstarts;
  Uintlist_T coords;
  List_T exons;

  int trstart, trend, trstart_overhang, trend_overhang;
  Chrpos_T trstart_splice_distance, trend_splice_distance;
  bool repairablep;

  debug11(printf("Entered Transcript_remap_invalid with path: "));
  debug11(Path_print(path));

  transcript_nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,transcript->num);
  if (transcript->genestrand > 0) {
    coords = Transcript_remap_endpoint_coords(path->endpoints,path->univdiagonals,path->junctions,
					      path->chrnum,path->chroffset,path->chrhigh,path->querylength,
					      uintlistpool,/*extend_qstart_p*/false,/*extend_qend_p*/false);
    if ((exons = Transcript_remap_geneplus(&trstart,&trend,&trstart_overhang,&trstart_splice_distance,
					   &trend_overhang,&trend_splice_distance,
					   queryseq,path->plusp,coords,exonbounds,exonstarts,
					   transcript_nexons,transcriptpool,listpool)) == NULL) {
      debug11(printf("exons could not be remapped, so transcript remains invalid\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      /* path->transcriptome_method_p = false; */
      assignedp = false;

    } else if (Exon_list_validp(&repairablep,exons) == false) {
      debug11(printf("Exon list is not valid\n\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      transcript = Transcript_new(/*num*/transcript->num,/*transcript_genestrand*/+1,
				  trstart,trend,trstart_overhang,trend_overhang,
				  exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				  transcriptpool);
      path->invalid_transcripts = Listpool_push(path->invalid_transcripts,listpool,(void *) transcript
						listpool_trace(__FILE__,__LINE__));
      /* path->transcriptome_method_p = false; */
      assignedp = true;
      
    } else {
      debug11(printf("Exon list is valid\n\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      transcript = Transcript_new(/*num*/transcript->num,/*transcript_genestrand*/+1,
				  trstart,trend,trstart_overhang,trend_overhang,
				  exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				  transcriptpool);
      path->transcripts = Listpool_push(path->transcripts,listpool,(void *) transcript
					listpool_trace(__FILE__,__LINE__));
      /* path->transcriptome_method_p = true; */
      /* validp = true; */
      assignedp = true;
    }

  } else {
    path = Path_reverse(path,/*expect_fwd_p*/false);  /* Needed for geneminus */

    coords = Transcript_remap_endpoint_coords(path->endpoints,path->univdiagonals,path->junctions,
					      path->chrnum,path->chroffset,path->chrhigh,path->querylength,
					      uintlistpool,/*extend_qstart_p*/false,/*extend_qend_p*/false);
    if ((exons = Transcript_remap_geneminus(&trstart,&trend,&trstart_overhang,&trstart_splice_distance,
					    &trend_overhang,&trend_splice_distance,
					    queryseq,path->plusp,coords,exonbounds,exonstarts,
					    transcript_nexons,transcriptpool,listpool)) == NULL) {
      debug11(printf("exons could not be remapped, so transcript remains invalid\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      /* path->transcriptome_method_p = false; */
      assignedp = false;
      
    } else if (Exon_list_validp(&repairablep,exons) == false) {
      debug11(printf("Exon list is not valid\n\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      transcript = Transcript_new(/*num*/transcript->num,/*transcript_genestrand*/-1,
				  trstart,trend,trstart_overhang,trend_overhang,
				  exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				  transcriptpool);
      path->invalid_transcripts = Listpool_push(path->invalid_transcripts,listpool,(void *) transcript
						listpool_trace(__FILE__,__LINE__));
      /* path->transcriptome_method_p = false; */
      assignedp = true;
      
    } else {
      debug11(printf("Exon list is valid\n\n"));
      /* Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool); */
      transcript = Transcript_new(/*num*/transcript->num,/*transcript_genestrand*/-1,
				  trstart,trend,trstart_overhang,trend_overhang,
				  exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				  transcriptpool);
      path->transcripts = Listpool_push(path->transcripts,listpool,(void *) transcript
					listpool_trace(__FILE__,__LINE__));
      /* path->transcriptome_method_p = true; */
      /* validp = true; */
      assignedp = true;
    }
    
    path = Path_reverse(path,/*expect_fwd_p*/true);  /* Needed for geneminus */
  }

  return assignedp;
}



/* Returns repairs */
/* Need queryseq so we can access the chopped part of the query to evaluate for poly-A */
List_T
Transcript_remap_all (List_T *transcripts, List_T *invalid_transcripts,
		      Intlist_T endpoints, Univcoordlist_T univdiagonals, List_T junctions,
		      Shortread_T queryseq, int querylength, bool plusp,
		      Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
		      Uintlistpool_T uintlistpool, Listpool_T listpool,
		      Transcriptpool_T transcriptpool,
		      int desired_genestrand, bool extend_qstart_p, bool extend_qend_p,
		      bool repairp) {

  List_T repairs = NULL;
  bool repairablep;
  Repair_T repair;
  int trstart_overhang, trend_overhang;
  Chrpos_T trstart_splice_distance, trend_splice_distance, chrlength;

  Chrnum_T map_chrnum;
  Univcoord_T genomicstart, genomicend, univdiagonal;
  Trnum_T trnum;
  int *matches, map_index;
  int nmatches, matchi;

  Uintlist_T coords;
  int transcript_nexons;
  
  int transcript_genestrand;
  Chrpos_T *exonstarts;
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  Chrpos_T *exonends;
#endif
  int *exonbounds;

  T transcript;
  List_T exons;
  int trstart, trend;
#ifdef DEBUG10
  bool alloc_label_p;
#endif


  assert(transcript_map_iit != NULL);

  *transcripts = (List_T) NULL;
  *invalid_transcripts = (List_T) NULL;
  if (endpoints == NULL) {
    /* Can happen for fusion endpoints */
    return (List_T) NULL;
  }

  debug10(printf("Entered Transcript_remap_all with chrnum %d, extend_qstart_p %d and extend_qend_p %d\n",
		 chrnum,extend_qstart_p,extend_qend_p));

  if ((map_chrnum = transcript_chrnum_crosstable[chrnum]) <= 0) {
    /* No corresponding transcript in map_iit file for this chrnum */
    debug10(printf("chrnum %d => map_chrnum %d, so returning NULL\n",chrnum,map_chrnum));
    return (List_T) NULL;
  } else {
    debug10(printf("chrnum %d => map_chrnum %d\n",chrnum,map_chrnum));
  }
  debug10(printf("endpoints %s\n",Intlist_to_string(endpoints)));
  debug10(printf("univdiagonals %s (%s)\n",
		 Univcoordlist_to_string(univdiagonals),Univcoordlist_to_string_offset(univdiagonals,chroffset)));

  coords = Transcript_remap_endpoint_coords(endpoints,univdiagonals,junctions,
					    chrnum,chroffset,chrhigh,querylength,
					    uintlistpool,extend_qstart_p,extend_qend_p);
  debug10(printf("coords %s\n",Uintlist_to_string(coords)));

  univdiagonal = Univcoordlist_head(univdiagonals);
  if (circularp[chrnum] == true &&
      univdiagonal > chroffset + (chrlength = (chrhigh - chroffset)/2) + querylength) {
    genomicstart = univdiagonal - chrlength - querylength;
  } else {
    genomicstart = univdiagonal - querylength;
  }

  univdiagonal = Univcoordlist_last_value(univdiagonals);
  if (circularp[chrnum] == true && univdiagonal > chroffset + chrlength + querylength) {
    genomicend = univdiagonal - chrlength;
  } else {
    genomicend = univdiagonal;
  }

  /* Note: Specifying sortp == true would ensure that map_indices are
     in ascending order, but trnums may not be. */
#if 0
  /* Should ignore desired_genestrand, which can be incorrect due to a faulty splice end */
  /* But leads to repairs by transcripts in the wrong sense */
  matches = IIT_get_with_divno(&nmatches,transcript_map_iit,map_chrnum,
			       genomicstart - chroffset + 1 /* 1-based */,
			       genomicend - chroffset,/*sortp*/false);
#else
  matches = IIT_get_signed_with_divno(&nmatches,transcript_map_iit,map_chrnum,
				      genomicstart - chroffset + 1 /* 1-based */,
				      genomicend - chroffset,/*sortp*/false,
				      /*sign*/desired_genestrand);
  debug10(printf("With desired genestrand %d, got %d matches\n",desired_genestrand,nmatches));
#endif

  for (matchi = 0; matchi < nmatches; matchi++) {
    map_index = matches[matchi];
    if ((trnum = Transcriptome_trnum(&transcript_nexons,&exonbounds,&exonstarts,transcriptome,map_index)) < 1) {
      /* Skip.  Not in transcriptome */
    } else {
      Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      debug10(printf("%d: %s with %d transcript exons\n",
		   trnum,IIT_label(transcript_map_iit,map_index,&alloc_label_p),transcript_nexons));

      if (transcript_genestrand > 0) {
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
	exonends = compute_exonends_geneplus(exonbounds,exonstarts,transcript_nexons,vectorpool);
	assert(IIT_interval_low(transcript_map_iit,map_index) == exonstarts[0]);
	assert(IIT_interval_high(transcript_map_iit,map_index) == exonends[transcript_nexons-1]);
#endif

	if ((exons = Transcript_remap_geneplus(&trstart,&trend,&trstart_overhang,&trstart_splice_distance,
					       &trend_overhang,&trend_splice_distance,
					       queryseq,plusp,coords,exonbounds,exonstarts,
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
					       exonends,
#endif
					       transcript_nexons,transcriptpool,listpool)) == NULL) {
	  /* Skip */
	} else if (trstart >= trend) {
	  /* Skip.  Entirely within introns */

	} else if (Exon_list_validp(&repairablep,exons) == true) {
	  debug10(printf("Exon list is valid\n\n"));
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *transcripts = Listpool_push(*transcripts,listpool,(void *) transcript
				       listpool_trace(__FILE__,__LINE__));

	} else if (repairp == false || repairablep == false) {
	  debug10(printf("Exon list is not valid\n\n"));
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *invalid_transcripts = 
	    Listpool_push(*invalid_transcripts,listpool,(void *) transcript
			  listpool_trace(__FILE__,__LINE__));

	} else {
	  debug10(printf("(1) Exon list is repairable."));
	  debug10(printf("  trstart overhang %d, trstart splice distance %u, trend overhang %d, trend splice_distance %u\n\n",
		       trstart_overhang,trstart_splice_distance,trend_overhang,trend_splice_distance));
	  /* Store transcript in case repair fails */
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *invalid_transcripts = 
	    Listpool_push(*invalid_transcripts,listpool,(void *) transcript
			  listpool_trace(__FILE__,__LINE__));

	  transcript = Transcript_copy(transcript,transcriptpool,listpool);
	  repair = Repair_new(transcript_genestrand,trstart_overhang,trstart_splice_distance,
			      trend_overhang,trend_splice_distance,transcript,listpool);
	  repairs = Listpool_push(repairs,listpool,(void *) repair
				  listpool_trace(__FILE__,__LINE__));
	}

      } else {
	coords = Uintlist_reverse(coords);
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
	exonends = compute_exonends_geneminus(exonbounds,exonstarts,transcript_nexons,vectorpool);
	assert(IIT_interval_high(transcript_map_iit,map_index) == exonstarts[0]);
	assert(IIT_interval_low(transcript_map_iit,map_index) == exonends[transcript_nexons-1]);
#endif

	if ((exons = Transcript_remap_geneminus(&trstart,&trend,&trstart_overhang,&trstart_splice_distance,
						&trend_overhang,&trend_splice_distance,
						queryseq,plusp,coords,exonbounds,exonstarts,
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
						exonends,
#endif
						transcript_nexons,transcriptpool,listpool)) == NULL) {
	  /* Skip */
	} else if (trstart >= trend) {
	  /* Skip.  Entirely within introns */

	} else if (Exon_list_validp(&repairablep,exons) == true) {
	  debug10(printf("Exon list is valid\n\n"));
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *transcripts = Listpool_push(*transcripts,listpool,(void *) transcript
				       listpool_trace(__FILE__,__LINE__));

	} else if (repairp == false || repairablep == false) {
	  debug10(printf("Exon list is not valid\n\n"));
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *invalid_transcripts = 
	    Listpool_push(*invalid_transcripts,listpool,(void *) transcript
			  listpool_trace(__FILE__,__LINE__));

	} else {
	  debug10(printf("(2) Exon list is repairable."));
	  debug10(printf("  trstart overhang %d, trstart splice distance %u, trend overhang %d, trend splice_distance %u\n\n",
		       trstart_overhang,trstart_splice_distance,trend_overhang,trend_splice_distance));
	  /* Store transcript in case repair fails */
	  transcript = Transcript_new(/*num*/trnum,transcript_genestrand,
				      trstart,trend,trstart_overhang,trend_overhang,
				      exons,transcript_nexons,/*trlength*/exonbounds[transcript_nexons-1],
				      transcriptpool);
	  *invalid_transcripts = 
	    Listpool_push(*invalid_transcripts,listpool,(void *) transcript
			  listpool_trace(__FILE__,__LINE__));

	  transcript = Transcript_copy(transcript,transcriptpool,listpool);
	  repair = Repair_new(transcript_genestrand,trstart_overhang,trstart_splice_distance,
			      trend_overhang,trend_splice_distance,transcript,listpool);
	  repairs = Listpool_push(repairs,listpool,(void *) repair
				  listpool_trace(__FILE__,__LINE__));
	}

	coords = Uintlist_reverse(coords);
      }
    }
  }

  FREE(matches);

  Uintlistpool_free_list(&coords,uintlistpool
			 uintlistpool_trace(__FILE__,__LINE__));


  *transcripts = Transcript_list_sort(*transcripts,listpool,transcriptpool);
  *invalid_transcripts = Transcript_list_sort(*invalid_transcripts,listpool,transcriptpool);
  repairs = Repair_make_unique(repairs,listpool,transcriptpool);

#ifdef DEBUG10
  printf("Returning transcripts ");
  Transcript_print_nums(*transcripts);
  printf("\n");
  printf("Returning invalid transcripts ");
  Transcript_print_nums(*invalid_transcripts);
  printf("\n");
#endif
  
#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(*transcripts);
  Transcript_list_ascendingp(*invalid_transcripts);
#endif

  return repairs;
}




static bool
matchp_exonends_geneplus (Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos,
			  int *exonbounds, Chrpos_T *exonstarts, int transcript_nexons) {
  Chrpos_T last_exonend;
  int last_bound, exonlength;
  int exoni;

  last_bound = 0;
  exonlength = exonbounds[0] - last_bound;
  last_exonend = exonstarts[0] + exonlength - 1;

  for (exoni = 1; exoni < transcript_nexons; exoni++) {
    debug2(printf("exon %d: comparing donor %u, acceptor %u with exonend %u, exonstart %u\n",
		  exoni,donor_chrpos,acceptor_chrpos,last_exonend,exonstarts[exoni]));
    if (last_exonend == donor_chrpos && exonstarts[exoni] == acceptor_chrpos) {
      return true;
    } else {
      exonlength = exonbounds[exoni] - last_bound;
      last_exonend = exonstarts[exoni] + exonlength - 1;
      last_bound = exonbounds[exoni];
    }
  }

  return false;
}

static bool
matchp_exonends_geneminus (Chrpos_T donor_chrpos, Chrpos_T acceptor_chrpos,
			   int *exonbounds, Chrpos_T *exonstarts, int transcript_nexons) {
  Chrpos_T last_exonend;
  int last_bound, exonlength;
  int exoni;

  last_bound = 0;
  exonlength = exonbounds[0] - last_bound;
  last_exonend = exonstarts[0] - exonlength + 1;

  for (exoni = 1; exoni < transcript_nexons; exoni++) {
    debug2(printf("exon %d: comparing donor %u, acceptor %u with exonend %u, exonstart %u\n",
		  exoni,donor_chrpos,acceptor_chrpos,last_exonend,exonstarts[exoni]));
    if (last_exonend == donor_chrpos && exonstarts[exoni] == acceptor_chrpos) {
      return true;
    } else {
      exonlength = exonbounds[exoni] - last_bound;
      last_exonend = exonstarts[exoni] - exonlength + 1;
      last_bound = exonbounds[exoni];
    }
  }

  return false;
}


/* Called by querystart_fusion_p and queryend_fusion_p*/
bool
Transcript_remap_matchp (Chrpos_T low_chrpos, Chrpos_T high_chrpos, Chrnum_T chrnum) {

  Chrnum_T map_chrnum;
  Trnum_T trnum;
  int *matches, map_index;
  int nmatches, matchi;

  int transcript_genestrand;
  Chrpos_T *exonstarts;
  int *exonbounds;

  int transcript_nexons;
#ifdef DEBUG
  bool alloc_label_p;
#endif

  if (transcript_map_iit == NULL) {
    return false;
  } else {
    debug2(printf("Looking for transcripts that match chrnum %d, low_chrpos %u, high_chrpos %u\n",
		  chrnum,low_chrpos,high_chrpos));
    if ((map_chrnum = transcript_chrnum_crosstable[chrnum]) <= 0) {
      return false;
    } else {
      matches = IIT_get_with_divno(&nmatches,transcript_map_iit,map_chrnum,
				   low_chrpos,high_chrpos,/*sortp*/false);
      
      for (matchi = 0; matchi < nmatches; matchi++) {
	map_index = matches[matchi];
	debug2(printf("%s\n",IIT_label(transcript_map_iit,map_index,&alloc_label_p)));
	
	if ((trnum = Transcriptome_trnum(&transcript_nexons,&exonbounds,&exonstarts,transcriptome,map_index)) < 1) {
	  /* Skip.  Not in transcriptome */
	} else {
	  Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	  if (transcript_genestrand > 0) {
	    if (matchp_exonends_geneplus(/*donor_chrpos*/low_chrpos,
					 /*acceptor_chrpos*/high_chrpos + 1/*1-based*/,
					 exonbounds,exonstarts,transcript_nexons) == true) {
	      FREE(matches);
	      return true;
	    }
	    
	  } else {
	    if (matchp_exonends_geneminus(/*donor_chrpos*/high_chrpos,
					  /*acceptor_chrpos*/low_chrpos + 1/*1-based*/,
					  exonbounds,exonstarts,transcript_nexons) == true) {
	      FREE(matches);
	      return true;
	    }
	  }
	}
      }

      FREE(matches);
    }

    return false;
  }
}
   
    
void
Transcript_remap_setup (bool *circularp_in,
			Transcriptome_T transcriptome_in,
			IIT_T transcript_map_iit_in,
			int *transcript_chrnum_crosstable_in) {
  circularp = circularp_in;
  transcriptome = transcriptome_in;
  transcript_map_iit = transcript_map_iit_in;
  transcript_chrnum_crosstable = transcript_chrnum_crosstable_in;
  return;
}

  
