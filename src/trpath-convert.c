#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath.h"
#include "trpath-convert.h"
#include "debug_groups.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "parameters.h"
#include "path-eval.h"
#include "path-solve.h"


#if defined(DEBUG0) || defined(DEBUG_TRPATH)
#define debug0(x) x
#else
#define debug0(x)
#endif

/* endpoints_acceptable_p */
#if defined(DEBUG2) || defined(DEBUG_TRPATH)
#define debug2(x) x
#else
#define debug2(x)
#endif


/* Conversion to Path_T */
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
#define debug10(x) x
#else
#define debug10(x)
#endif


/* binary_search */
#if defined(DEBUG20) || defined(DEBUG_TRPATH)
#define debug20(x) x
#else
#define debug20(x)
#endif


static Transcriptome_T transcriptome;
static EF64_T chromosome_ef64;

#define MAX_TRANSCRIPT_SPLICE_BOUNDARY_DELETION 10
#define MIN_TRANSCRIPT_SPLICE_BOUNDARY_DELETION_INTRON 1000


#define T Trpath_T

static bool
mismatches_all_zero_p (Intlist_T nmismatches, Intlist_T ref_nmismatches) {
  Intlist_T p;

  for (p = nmismatches; p != NULL; p = Intlist_next(p)) {
    if (Intlist_head(p) != 0) {
      return false;
    }
  }

  for (p = ref_nmismatches; p != NULL; p = Intlist_next(p)) {
    if (Intlist_head(p) != 0) {
      return false;
    }
  }

  return true;
}

static Univcoord_T
advance_univdiagonal (Univcoord_T univdiagonal, int direction, int adj) {
  if (direction > 0) {
    return univdiagonal + (Univcoord_T) adj;
  } else {
    return univdiagonal - (Univcoord_T) adj;
  }
}

/* Emits a known intron adjacent to a short transcriptome deletion as two explicit
   junctions instead of a middle JUNCTION_UNSOLVED.  The caller has already
   emitted the aligned segment ending at endpoint; this routine adds the
   zero-length segment between the splice and deletion. */
static Univcoord_T
push_splice_boundary_deletion (List_T *junctions, Intlist_T *endpoints,
			       Univcoordlist_T *univdiagonals,
			       Intlist_T *nmismatches, Intlist_T *ref_nmismatches,
			       Univcoord_T univdiagonal, int endpoint,
			       Chrpos_T splice_distance, int deletionlen,
			       int sensedir, bool deletion_first_p, int direction,
			       Listpool_T listpool, Intlistpool_T intlistpool,
			       Univcoordlistpool_T univcoordlistpool,
			       Pathpool_T pathpool) {
  Junction_T first_junction, second_junction;
  Univcoord_T middle_univdiagonal;

  if (deletion_first_p == true) {
    first_junction = Junction_new_deletion(deletionlen,pathpool);
    second_junction = Junction_new_splice_annotated(splice_distance,sensedir,
						    /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
  } else {
    first_junction = Junction_new_splice_annotated(splice_distance,sensedir,
						   /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
    second_junction = Junction_new_deletion(deletionlen,pathpool);
  }

  middle_univdiagonal = advance_univdiagonal(univdiagonal,direction,Junction_adj(first_junction));
  *junctions = Listpool_push(*junctions,listpool,(void *) first_junction
			     listpool_trace(__FILE__,__LINE__));
  *endpoints = Intlistpool_push(*endpoints,intlistpool,endpoint
				intlistpool_trace(__FILE__,__LINE__));
  *univdiagonals = Univcoordlistpool_push(*univdiagonals,univcoordlistpool,middle_univdiagonal
					  univcoordlistpool_trace(__FILE__,__LINE__));
  *nmismatches = Intlistpool_push(*nmismatches,intlistpool,0
				  intlistpool_trace(__FILE__,__LINE__));
  *ref_nmismatches = Intlistpool_push(*ref_nmismatches,intlistpool,0
				      intlistpool_trace(__FILE__,__LINE__));
  *junctions = Listpool_push(*junctions,listpool,(void *) second_junction
			     listpool_trace(__FILE__,__LINE__));

  return advance_univdiagonal(middle_univdiagonal,direction,Junction_adj(second_junction));
}

/* exonbounds does not include the initial value of 0, so caller needs a special case when the initial exoni is 0 */
static int
binary_search (int lowi, int highi, int *exonbounds, int trstart) {
  int middlei;

  debug20(printf("entered binary search with lowi=%d, highi=%d, trstart=%d (but testing at %d)\n",
		 lowi,highi,trstart,trstart+1));
  trstart += 1;			/* Because when trstart == exonbound, we want the next exon */

  while (lowi < highi) {
    /* middlei = lowi + ((highi - lowi) / 2); */
    middlei = (lowi + highi)/2;
    assert(middlei != 0);
    debug20(printf("  binary: %d:%u %d:%u %d:%u   vs. %u\n",
		   lowi,exonbounds[lowi],middlei,exonbounds[middlei],
		   highi,exonbounds[highi],trstart));
    if (trstart < exonbounds[middlei]) {
      highi = middlei;
    } else if (trstart > exonbounds[middlei]) {
      lowi = middlei + 1;
    } else {
      debug20(printf("binary search returns middlei %d - 1\n",middlei));
      return middlei - 1;
    }
  }

  debug20(printf("(2) binary search returns highi %d - 1\n",highi));
  return highi - 1;
}


static List_T
compute_exons (int trstart, int trend, int *exonbounds, int nexons,
	       Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T exons;
  char firstchar, lastchar;
  int exoni, exonj, i;

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("compute_exons with trstart %d, trend %d\n",trstart,trend);
  for (i = 0; i < nexons; i++) {
    printf("exon %d: exonbound %d\n",i,exonbounds[i]);
  }
#endif

#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  exoni = nexons;
  /* Test: >= trstart + 1 or > trstart */
  while (exoni >= 1 && exonbounds[exoni - 1] > trstart) {
    exoni--;
  }
#ifdef TEST_BINARY_SEARCH
  if (trstart < exonbounds[0]) {
    assert(exoni == 0);
  } else {
    assert(exoni == binary_search(0,nexons,&(exonbounds[-1]),trstart));
  }
#endif
#else
  if (trstart < exonbounds[0]) {
    exoni = 0;
  } else {
    exoni = binary_search(0,nexons,&(exonbounds[-1]),trstart);
  }
#endif


  debug10(printf("exoni is %d\n",exoni));

  exonj = exoni;
  while (exonj < nexons - 1 && exonbounds[exonj] < trend) { /* Not <= */
    exonj++;
  }
  debug10(printf("exonj is %d\n",exonj));

  if (exoni == 0) {
    firstchar = '.';
  } else if (trstart == exonbounds[exoni - 1]) {
    firstchar = 's';
  } else {
    firstchar = '.';
  }

  if (exonj == nexons - 1) {
    lastchar = '.';
  } else if (trend == exonbounds[exonj]) {
    lastchar = 's';
  } else {
    lastchar = '.';
  }

  if (exonj == exoni) {
    exons = Listpool_push(NULL,listpool,
			  Exon_new(firstchar,exoni,lastchar,transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
  } else {
    i = exoni;
    exons = Listpool_push(NULL,listpool,
			  Exon_new(firstchar,i,'s',transcriptpool)
			  listpool_trace(__FILE__,__LINE__));

    i++;
    while (i < exonj) {
      exons = Listpool_push(exons,listpool,
			    Exon_new('s',i,'s',transcriptpool)
			    listpool_trace(__FILE__,__LINE__));
      i++;
    }

    exons = Listpool_push(exons,listpool,
			  Exon_new('s',i,lastchar,transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
    
  }

  debug10(Exon_print_list_stdout(exons));

  assert(exons != NULL);
  return List_reverse(exons);
}


static bool
trim_initial_transcript_start (int *tstart, int *trstart, int tend,
			       Intlist_T tr_endpoints) {
  int trim;

  if (*trstart >= 0) {
    return true;
  } else if ((trim = -(*trstart)) + (*tstart) >= tend) {
    return false;
  } else {
    *tstart += trim;
    *trstart = 0;
    Intlist_head_set(tr_endpoints,*tstart);
    return true;
  }
}



static Path_T
convert_to_path_geneplus (int *nintrons, int *found_score, int *end_score,
			  Trpath_T trpath, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  bool plusp, bool first_read_p, int sensedir,
			  Shortread_T queryseq, char *queryptr, int querylength,
			  Stage1_T this, Knownsplicing_T knownsplicing, Compress_T query_compress,
			  Compress_T query_compress_fwd, Compress_T query_compress_rev,

			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			  int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool,
			  Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			  Hitlistpool_T hitlistpool, Method_T method) {

  Path_T path, newpath;
  bool unsolvedp = false, evalp = false;

  int *exonbounds;
  Chrpos_T *exonstarts;
  int nexons;
  Trnum_T trnum;

  List_T exons = NULL;
#ifdef COMPUTE_EXONS_INCREMENTALLY
  Exon_T exon;
  char firstchar, lastchar;
#endif
  int exoni;
#ifdef TEST_EXONS
  List_T exons_simple;
#endif

  Intlist_T endpoints;
  Univcoordlist_T univdiagonals = NULL;
  Intlist_T nmismatches = NULL, ref_nmismatches = NULL;
  List_T junctions = NULL;

  Univcoord_T univdiagonal;
  Chrpos_T splice_distance;
  Junction_T junction;

  int exstart, exend;
  int trstart, trend;		/* From beginning of transcript */
  int tstart, tend;		/* From beginning of read */
  int exonpos, exonlength, uselength, adj;

  Transcript_T transcript;
  int overall_trstart, overall_trend;

  /* Procedure modifies tr_endpoints and tr_nmismatches, which prevenst a second call to convert trpath */
  Intlist_T tr_endpoints = Intlistpool_copy(trpath->endpoints,intlistpool);
  Uintlist_T trdiagonals = trpath->trdiagonals;
  Intlist_T tr_nmismatches = Intlistpool_copy(trpath->nmismatches,intlistpool);
  List_T tr_junctions = trpath->junctions;
  Trcoord_T troffset = trpath->troffset;
  Trcoord_T trdiagonal;

  int local_trimpos, local_nmismatches;
  int ntrim_qstart, ntrim_qend, pos5, pos3;


  *nintrons = 0;

  trnum = trpath->trnum;
  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);

  debug0(printf("Entered convert_to_path_geneplus with trpath %p, trnum %d, chrnum %d from method %s\n",
		trpath,trnum,chrnum,Method_string(method)));

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("\nEntered convert_to_path_geneplus with trpath %p, trnum %d, chrnum %d from method %s\n",
	 trpath,trnum,chrnum,Method_string(method));
  printf("tr_endpoints: %s (querylength %d)\n",Intlist_to_string(tr_endpoints),querylength);
  printf("trdiagonals: %s\n",Uintlist_to_string(trdiagonals));
  printf("tr_nmismatches: %s\n",Intlist_to_string(tr_nmismatches));
  printf("tr_junctions: ");
  Junction_print_list(tr_junctions);
  printf("\n");
  for (exoni = 0; exoni < nexons; exoni++) {
    printf("exon %d: exonbound %d, exonstart %u\n",exoni,exonbounds[exoni],exonstarts[exoni]);
  }
#endif

  tstart = Intlist_head(tr_endpoints);
  tend = Intlist_second_value(tr_endpoints);

  trdiagonal = Uintlist_head(trdiagonals);
  trstart = (trdiagonal - troffset - querylength) + tstart;
  trend = (trdiagonal - troffset - querylength) + tend;
  debug10(printf("Initial tstart %d => trstart %d.  tend %d => trend %d.\n",
		 tstart,trstart,tend,trend));
  assert(trstart < trend);
  if (trim_initial_transcript_start(&tstart,&trstart,tend,tr_endpoints) == false) {
    debug10(printf("starting before coordinate 0 consumes the segment => aborting\n"));
    return (Path_T) NULL;
  } else if (trend >= exonbounds[nexons - 1]) {
    debug10(printf("ending past the last exonbound => aborting\n"));
    return (Path_T) NULL;
  }


  /* Find initial exon */
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  exoni = 0;
  exstart = 0;
  exend = exonbounds[0];
  while (exoni < nexons && trstart >= exend) {
    exstart = exend;
    exoni += 1;
    debug10(printf(" trstart %d >= exend %d => Advancing exoni to be %d\n",trstart,exend,exoni));
    exend = exonbounds[exoni];
  }
#ifdef TEST_BINARY_SEARCH
  if (trstart < exonbounds[0]) {
    assert(exoni == 0);
  } else {
    assert(exoni == binary_search(0,nexons,&(exonbounds[-1]),trstart));
  }
#endif
#else
  if (trstart < exonbounds[0]) {
    exoni = 0;
    exstart = 0;
    exend = exonbounds[0];
  } else {
    exoni = binary_search(0,nexons,&(exonbounds[-1]),trstart);
    exstart = exonbounds[exoni-1];
    exend = exonbounds[exoni];
  }
#endif

  exonpos = trstart - exstart;
  exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
  if (exoni == 0) {
    firstchar = '.';
  } else if (exonpos == 0) {
    firstchar = 's';
  } else {
    firstchar = '.';
  }
#endif
  debug10(printf("Initial exon %d\n",exoni));

  /* Initial endpoint */
  endpoints = Intlistpool_push(NULL,intlistpool,tstart
			       intlistpool_trace(__FILE__,__LINE__)); /* gplus */

  /* Need to subtract 1 from from exonstarts[exoni] only for geneplus */
  univdiagonal = chroffset + (Univcoord_T) (exonstarts[exoni] - 1/*1-based*/ + exonpos - tstart + querylength); /* geneplus */
  trdiagonal = Uintlist_head(trdiagonals);
  /* trend = (trdiagonal - troffset - querylength) + tend; -- computed above */
  
  while (exoni < nexons && trdiagonals != NULL) {
    exonpos = trend - exstart;
    debug10(printf(">>Segment univdiagonal %u, tcoords %d..%d, trcoords %d..%d\n",
		   univdiagonal,tstart,tend,trstart,trend));
    debug10(printf("Exoni %d, exstart %d, exend %d, exonpos %d/%d, univdiagonal %u\n",
		   exoni,exstart,exend,exonpos,exonlength,univdiagonal));
    debug10(printf(" tend %d => Exonpos is %d = trend %d - exstart %d, out of an exonlength of %d\n",
		   tend,exonpos,trend,exstart,exonlength));

    if (exonpos <= exonlength) {
      /* Transcript segment is in the exon */
      debug10(printf(">>Output within-exon: univdiagonal %u, tcoords %d..%d, trcoords %d..%d\n",
		     univdiagonal,tstart,tend,trstart,trend));
      endpoints = Intlistpool_push(endpoints,intlistpool,tend
				   intlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(tr_nmismatches)
				     intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(tr_nmismatches)
					 intlistpool_trace(__FILE__,__LINE__));

      if ((trdiagonals = Uintlist_next(trdiagonals)) == NULL) {
	/* End of trpath */

      } else if (exoni < nexons - 1 &&
		 exonpos + Junction_adj((Junction_T) List_head(tr_junctions)) >= exonlength - 3) {
	/* gplus pre-intron: Indel is at or near donor splice site */
	debug10(printf("gplus pre-intron\n"));
	++exoni;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	lastchar = (exonpos < exonlength) ? 'y' : 'x';
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp before intron of length %u\n",
		       adj,exonlength - (exonpos + adj),splice_distance));

	if (Parameters_fix_splice_boundary_indels_p == true &&
	    adj > 0 && adj <= MAX_TRANSCRIPT_SPLICE_BOUNDARY_DELETION &&
	    splice_distance >= MIN_TRANSCRIPT_SPLICE_BOUNDARY_DELETION_INTRON) {
	  univdiagonal = push_splice_boundary_deletion(&junctions,&endpoints,&univdiagonals,
						       &nmismatches,&ref_nmismatches,
						       univdiagonal,tend,splice_distance,adj,
						       sensedir,/*deletion_first_p*/exonlength - (exonpos + adj) > 0,
						       /*direction*/+1,listpool,intlistpool,
						       univcoordlistpool,pathpool);
	  tr_junctions = List_next(tr_junctions);
	  (*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf(">>Output splice-boundary deletion.  Advancing tr_junctions\n");
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	} else if ((splice_distance += adj) == 0) {
	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Junction_list_gc(&junctions,listpool,pathpool);
	  Exon_list_gc(&exons,listpool,transcriptpool);

	  Intlistpool_free_list(&tr_endpoints,intlistpool);
	  Intlistpool_free_list(&tr_nmismatches,intlistpool);

	  return (Path_T) NULL;

	} else {
	  junctions = Listpool_push(junctions,listpool,(void *) JUNCTION_UNSOLVED
				    listpool_trace(__FILE__,__LINE__));
	  debug10(printf(">>Output UNSOLVED splice junction.  Advancing tr_junctions\n"));
	  tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	  Intlist_head_set(nmismatches,-1);
	  Intlist_head_set(tr_nmismatches->rest,-1);
	  univdiagonal += splice_distance; /* with Junction_adj.  geneplus */
	  evalp = true;
	  unsolvedp = true;
	}

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);

	tstart = tend; trstart = trend;
	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

      } else if (exoni < nexons - 1 && exonpos >= exonlength - 3) {
	/* gplus post-intron: Indel is at or near acceptor splice site */
	debug10(printf("gplus post-intron\n"));
	++exoni;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	if (exonpos < exonlength) {
	  lastchar = 'y';
	} else if (exonpos > exonlength) {
	  lastchar = 'x';
	} else {
	  lastchar = 's';
	}
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp after intron of length %u\n",
		       adj,exonlength - exonpos,splice_distance));

	if (Parameters_fix_splice_boundary_indels_p == true &&
	    adj > 0 && adj <= MAX_TRANSCRIPT_SPLICE_BOUNDARY_DELETION &&
	    splice_distance >= MIN_TRANSCRIPT_SPLICE_BOUNDARY_DELETION_INTRON) {
	  univdiagonal = push_splice_boundary_deletion(&junctions,&endpoints,&univdiagonals,
						       &nmismatches,&ref_nmismatches,
						       univdiagonal,tend,splice_distance,adj,
						       sensedir,/*deletion_first_p*/false,
						       /*direction*/+1,listpool,intlistpool,
						       univcoordlistpool,pathpool);
	  tr_junctions = List_next(tr_junctions);
	  (*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf(">>Output splice-boundary deletion.  Advancing tr_junctions\n");
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	} else if ((splice_distance += adj) == 0) {
	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Junction_list_gc(&junctions,listpool,pathpool);
	  Exon_list_gc(&exons,listpool,transcriptpool);

	  Intlistpool_free_list(&tr_endpoints,intlistpool);
	  Intlistpool_free_list(&tr_nmismatches,intlistpool);

	  return (Path_T) NULL;

	} else {
	  junctions = Listpool_push(junctions,listpool,(void *) JUNCTION_UNSOLVED
				    listpool_trace(__FILE__,__LINE__));
	  debug10(printf(">>Output UNSOLVED splice junction.  Advancing tr_junctions\n"));
	  tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	  Intlist_head_set(nmismatches,-1);
	  Intlist_head_set(tr_nmismatches->rest,-1);
	  univdiagonal += splice_distance; /* with Junction_adj.  geneplus */
	  evalp = true;
	  unsolvedp = true;
	}

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);

	tstart = tend; trstart = trend;
	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

      } else {
	/* Indel is within exon */

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);
	tstart = tend; trstart = trend;

	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

	junction = Junction_copy((Junction_T) List_head(tr_junctions),pathpool);
	/* Needed only for geneminus, to move genomic endpoint before insertion */
	/* Intlist_head_incr(endpoints,-Junction_ninserts(junction)); */
#if 0
	if (Junction_type(junction) == DEL_JUNCTION) {
	  /* Need to translate transcriptome-based deletionpos to genomic deletionpos */
	  deletionpos = univdiagonal - querylength + tstart;
	  debug10(printf("At geneplus deletion, computing deletionpos to be %u\n",deletionpos));
	  Junction_set_deletionpos(junction,deletionpos);
	}
#endif
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf(">>(1) Output indel junction: ");
	Junction_print(junction);
	printf("\n");
#endif
	univdiagonal += Junction_adj(junction); /* geneplus */
	tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf("Head of junctions: ");
	if (tr_junctions != NULL) {
	  Junction_print((Junction_T) List_head(tr_junctions));
	}
	printf("\n");
#endif	

	/* Do not advance exon information */
      }

    } else {
      /* Transcript segment extends to the next exon */

      /* Revise transcript segment */
      trstart = (trdiagonal - troffset - querylength) + tstart;
      uselength = exend - trstart;
      debug10(printf(">>Output next-exon: univdiagonal %u, %d..%d\n",
		     univdiagonal,tstart,tstart + uselength));

      tstart += uselength; trstart += uselength;
      Intlist_head_set(tr_endpoints,tstart);
      endpoints = Intlistpool_push(endpoints,intlistpool,tstart
				   intlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));

      if (Intlist_head(tr_nmismatches) == 0) {
	nmismatches = Intlistpool_push(nmismatches,intlistpool,0
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,0
					   intlistpool_trace(__FILE__,__LINE__));
      } else {
	/* Don't know how to divide a non-zero value of nmismatches between the exons */
	nmismatches = Intlistpool_push(nmismatches,intlistpool,-1
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,-1
					   intlistpool_trace(__FILE__,__LINE__));
	Intlist_head_set(tr_nmismatches,-1);
	evalp = true;
      }
      
      if (++exoni < nexons) {
#ifdef COMPUTE_EXONS_INCREMENTALLY
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,/*lastchar*/'s',transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	junction = Junction_new_splice_annotated(splice_distance,sensedir,
						 /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
	(*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf(">>Output splice junction: ");
	Junction_print(junction);
	printf("\n");
#endif

	univdiagonal += splice_distance; /* geneplus */

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif
      }
    }
  }

  debug10(printf(">>Final exoni %d, tcoords %d..%d, trcoords %d..%d\n",
		 exoni,tstart,tend,trstart,trend));
#ifdef COMPUTE_EXONS_INCREMENTALLY
  if (exoni >= nexons - 1) {
    exons = Listpool_push(exons,listpool,
			  Exon_new(firstchar,nexons-1,/*lastchar*/'.',transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
  } else {
    exonpos = trend - exstart;
    if (exonpos == exonlength) {
      lastchar = 's';
    } else {
      lastchar = '.';
    }
    exons = Listpool_push(exons,listpool,
			  Exon_new(firstchar,exoni,lastchar,transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
  }
  exons = List_reverse(exons);

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("New exons: ");
  Exon_print_list_stdout(exons);
  printf("\n");
#endif
#endif


  /* Reverse lists to put them in gplus order */
  endpoints = Intlist_reverse(endpoints);
  univdiagonals = Univcoordlist_reverse(univdiagonals);
  nmismatches = Intlist_reverse(nmismatches);
  ref_nmismatches = Intlist_reverse(ref_nmismatches);
  junctions = List_reverse(junctions);

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("endpoints: %s\n",Intlist_to_string(endpoints));
  printf("univdiagonals: %s (%s)\n",
	 Univcoordlist_to_string(univdiagonals),Univcoordlist_to_string_offset(univdiagonals,chroffset));
  printf("nmismatches: %s\n",Intlist_to_string(nmismatches));
  printf("junctions: ");
  Junction_print_list(junctions);
  printf("\n");
#endif

  if (Path_endpoints_acceptable_p(endpoints,junctions) == false) {
    debug10(printf("Endpoints not acceptable\n"));
    Intlistpool_free_list(&endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
    Intlistpool_free_list(&nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Intlistpool_free_list(&ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Junction_list_gc(&junctions,listpool,pathpool);
#ifdef COMPUTE_EXONS_INCREMENTALLY
    Exon_list_gc(&exons,listpool,transcriptpool);
#endif

    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);

    return (Path_T) NULL;
    
  } else if (unsolvedp == true) {
    path = Path_convert_simple(chrnum,chroffset,chrhigh,
			       endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,evalp,
			       plusp,first_read_p,/*genestrand*/0,sensedir,querylength,
			       query_compress_fwd,query_compress_rev,listpool,pathpool,method);
    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);

    if ((newpath = Path_solve_junctions(&(*found_score),&(*end_score),path,sensedir,/*genestrand*/0,
					query_compress,query_compress_fwd,query_compress_rev,
					queryseq,queryptr,querylength,

					novel_univdiags_alloc,novel_diagonals_alloc,
					localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					localdb,localdb_nmismatches_allowed,

					this,knownsplicing,
					uintlistpool,intlistpool,univcoordlistpool,listpool,
					pathpool,transcriptpool)) == NULL) {
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      return (Path_T) NULL;

    } else {
      /* Path_solve_junctions called Path_eval_nmatches */
      debug10(printf("Have solved path: "));
      debug10(Path_print(newpath));

      newpath->completep = true;

      overall_trstart = (trdiagonal - troffset - querylength) + Intlist_head(path->endpoints);
      overall_trend = (trdiagonal - troffset - querylength) + Intlist_last_value(path->endpoints);
      exons = compute_exons(overall_trstart,overall_trend,exonbounds,nexons,
			    listpool,transcriptpool);

      transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/+1,
				  overall_trstart+/*1-based*/1,overall_trend,
				  /*trstart_overhang*/0,/*trend_overhang*/0,
				  exons,nexons,/*trlength*/exonbounds[nexons-1],
				  transcriptpool);
      /* Path_consolidate and Concordance_tr need to handle paths without valid transcripts */
      newpath->invalid_transcripts = Listpool_push(NULL,listpool,(void *) transcript
						   listpool_trace(__FILE__,__LINE__));
      return newpath;
    }

  } else {
#ifdef TEST_EXONS
    exons_simple = compute_exons_simple(overall_trstart,overall_trend,exonbounds,nexons,
					listpool,transcriptpool);
    printf("Correct exons: ");
    Exon_print_list_stdout(exons_simple);
    printf("\n");
    assert(Exon_list_equal(exons,exons_simple));
    Exon_list_gc(&exons_simple,listpool,transcriptpool);
#endif

    path = Path_convert_simple(chrnum,chroffset,chrhigh,
			       endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,evalp,
			       plusp,first_read_p,/*genestrand*/0,sensedir,querylength,
			       query_compress_fwd,query_compress_rev,listpool,pathpool,method);

    if (trpath->found_score == 0 && path->found_score == 0 &&
	mismatches_all_zero_p(path->nmismatches,path->ref_nmismatches) == true) {
      /* Can skip Path_eval_nmatches only when genomic conversion is also exact. */
      path->found_score = path->score_within_trims = 0;
      path->nmatches = path->ref_nmatches = querylength /*- trpath->found_score*/ - trpath->total_ninserts;
      path->total_splice_prob = path->junction_splice_prob = (*nintrons) * 4.0;
      *found_score = 0;
      *end_score = 0;
    } else {
      debug10(printf("(1) Calling Path_eval_nmatches\n"));
      Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
    }

#if 0
    /* Not effective.  Pruning transcriptome paths early yields to partial
       trpaths, which can obscure complete genomic paths */
    if (Path_segments_good_p(path) == false) {
      debug10(printf("(1) Path segments are not good, so returning NULL\n"));
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      return (Path_T) NULL;
    }
#endif

    overall_trstart = (trdiagonal - troffset - querylength) + Intlist_head(path->endpoints);
    overall_trend = (trdiagonal - troffset - querylength) + Intlist_last_value(path->endpoints);
    exons = compute_exons(overall_trstart,overall_trend,exonbounds,nexons,
			  listpool,transcriptpool);
    
    transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/+1,
				overall_trstart+/*1-based*/1,overall_trend,
				/*trstart_overhang*/0,/*trend_overhang*/0,
				exons,nexons,/*trlength*/exonbounds[nexons-1],
				transcriptpool);
    path->transcripts = Listpool_push(NULL,listpool,(void *) transcript
				      listpool_trace(__FILE__,__LINE__));
    debug10(printf("Have path: "));
    debug10(Path_print(path));

    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);
    
    return path;
  }
}


static Path_T
convert_to_path_geneminus (int *nintrons, int *found_score, int *end_score,
			   Trpath_T trpath, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, bool first_read_p, int sensedir,
			   Shortread_T queryseq, char *queryptr, int querylength,
			   Stage1_T this, Knownsplicing_T knownsplicing, Compress_T query_compress,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,

			   Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			   unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			   int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
			   Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
			   Hitlistpool_T hitlistpool, Method_T method) {
  Path_T path, newpath;
  bool unsolvedp = false, evalp = false;

  int *exonbounds;
  Chrpos_T *exonstarts;
  int nexons;
  Trnum_T trnum;

  List_T exons = NULL;
#ifdef COMPUTE_EXONS_INCREMENTALLY
  Exon_T exon;
  char firstchar, lastchar;
#endif
  int exoni;
#ifdef TEST_EXONS
  List_T exons_simple;
#endif

  Intlist_T endpoints;
  Univcoordlist_T univdiagonals = NULL;
  Intlist_T nmismatches = NULL, ref_nmismatches = NULL;
  List_T junctions = NULL;

  Univcoord_T univdiagonal;
  Chrpos_T splice_distance;
  Junction_T junction;

  int exstart, exend;
  int trstart, trend;		/* From beginning of transcript */
  int tstart, tend;		/* From beginning of read */
  int exonpos, exonlength, uselength, adj;

  Transcript_T transcript;
  int overall_trstart, overall_trend;

  /* Procedure modifies tr_endpoints and tr_nmismatches, which prevents a second call to convert trpath */
  Intlist_T tr_endpoints = Intlistpool_copy(trpath->endpoints,intlistpool);
  Uintlist_T trdiagonals = trpath->trdiagonals;
  Intlist_T tr_nmismatches = Intlistpool_copy(trpath->nmismatches,intlistpool);
  List_T tr_junctions = trpath->junctions;
  Trcoord_T troffset = trpath->troffset;
  Trcoord_T trdiagonal;

  int local_trimpos, local_nmismatches;
  int ntrim_qstart, ntrim_qend, pos5, pos3;


  *nintrons = 0;

  trnum = trpath->trnum;
  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);

  debug0(printf("Entered convert_to_path_geneminus with trpath %p, trnum %d, chrnum %d from method %s\n",
		trpath,trnum,chrnum,Method_string(method)));

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("\nEntered convert_to_path_geneminus with trpath %p, trnum %d, chrnum %d from method %s\n",
	 trpath,trnum,chrnum,Method_string(method));
  printf("tr_endpoints: %s (querylength %d)\n",Intlist_to_string(tr_endpoints),querylength);
  printf("trdiagonals: %s\n",Uintlist_to_string(trdiagonals));
  printf("tr_nmismatches: %s\n",Intlist_to_string(tr_nmismatches));
  printf("tr_junctions: ");
  Junction_print_list(tr_junctions);
  printf("\n");
  for (exoni = 0; exoni < nexons; exoni++) {
    printf("exon %d: exonbound %d, exonstart %u\n",exoni,exonbounds[exoni],exonstarts[exoni]);
  }
#endif

  tstart = Intlist_head(tr_endpoints);
  tend = Intlist_second_value(tr_endpoints);

  trdiagonal = Uintlist_head(trdiagonals);
  trstart = (trdiagonal - troffset - querylength) + tstart;
  trend = (trdiagonal - troffset - querylength) + tend;
  debug10(printf("Initial tstart %d => trstart %d.  tend %d => trend %d.\n",
		 tstart,trstart,tend,trend));
  assert(trstart < trend);
  if (trim_initial_transcript_start(&tstart,&trstart,tend,tr_endpoints) == false) {
    debug10(printf("starting before coordinate 0 consumes the segment => aborting\n"));
    return (Path_T) NULL;
  } else if (trend>= exonbounds[nexons - 1]) {
    debug10(printf("ending past the last exonbound => aborting\n"));
    return (Path_T) NULL;
  }

  /* Find initial exon */
#if defined(LINEAR_SEARCH) || defined(TEST_BINARY_SEARCH)
  exoni = 0;
  exstart = 0;
  exend = exonbounds[0];
  while (exoni < nexons && trstart >= exend) {
    exstart = exend;
    exoni += 1;
    debug10(printf(" trstart %d >= exend %d => Advancing exoni to be %d\n",trstart,exend,exoni));
    exend = exonbounds[exoni];
  }
#ifdef TEST_BINARY_SEARCH
  if (trstart < exonbounds[0]) {
    assert(exoni == 0);
  } else {
    assert(exoni == binary_search(0,nexons,&(exonbounds[-1]),trstart));
  }
#endif
#else
  if (trstart < exonbounds[0]) {
    exoni = 0;
    exstart = 0;
    exend = exonbounds[0];
  } else {
    exoni = binary_search(0,nexons,&(exonbounds[-1]),trstart);
    exstart = exonbounds[exoni-1];
    exend = exonbounds[exoni];
  }
#endif

  exonpos = trstart - exstart;
  exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
  if (exoni == 0) {
    firstchar = '.';
  } else if (exonpos == 0) {
    firstchar = 's';
  } else {
    firstchar = '.';
  }
#endif
  debug10(printf("Initial exon %d\n",exoni));

  /* Initial endpoint */
  endpoints = Intlistpool_push(NULL,intlistpool,querylength - tstart
			       intlistpool_trace(__FILE__,__LINE__)); /* gminus */

  /* Need to subtract 1 from from exonstarts[exoni] only for geneplus */
  univdiagonal = chroffset + (Univcoord_T) (exonstarts[exoni] - exonpos + tstart); /* geneminus */
  trdiagonal = Uintlist_head(trdiagonals);
  /* trend = (trdiagonal - troffset - querylength) + tend; -- computed above */

  while (exoni < nexons && trdiagonals != NULL) {
    exonpos = trend - exstart;
    debug10(printf(">>Segment univdiagonal %u, tcoords %d..%d, trcoords %d..%d\n",
		   univdiagonal,tstart,tend,trstart,trend));
    debug10(printf("Exoni %d, exstart %d, exend %d, exonpos %d/%d, univdiagonal %u\n",
		   exoni,exstart,exend,exonpos,exonlength,univdiagonal));
    debug10(printf(" => Exonpos is %d = trend %d - exstart %d, out of an exonlength of %d\n",
		   exonpos,trend,exstart,exonlength));

    if (exonpos <= exonlength) {
      /* Transcript segment is in the exon */
      debug10(printf(">>Output within-exon: univdiagonal %u, tcoords %d..%d, trcoords %d..%d\n",
		     univdiagonal,tstart,tend,trstart,trend));
      endpoints = Intlistpool_push(endpoints,intlistpool,querylength - tend
				   intlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(tr_nmismatches)
				     intlistpool_trace(__FILE__,__LINE__));
      ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,Intlist_head(tr_nmismatches)
					 intlistpool_trace(__FILE__,__LINE__));

      if ((trdiagonals = Uintlist_next(trdiagonals)) == NULL) {
	/* End of trpath */
	
      } else if (exoni < nexons - 1 &&
		 exonpos + Junction_adj((Junction_T) List_head(tr_junctions)) >= exonlength - 3) {
	/* gminus pre-intron: Indel is at or near donor splice site */
	debug10(printf("gplus pre-intron\n"));
	++exoni;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	lastchar = (exonpos < exonlength) ? 'y' : 'x';
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp before intron of length %u\n",
		       adj,exonlength - (exonpos + adj),splice_distance));

	if (Parameters_fix_splice_boundary_indels_p == true &&
	    adj > 0 && adj <= MAX_TRANSCRIPT_SPLICE_BOUNDARY_DELETION &&
	    splice_distance >= MIN_TRANSCRIPT_SPLICE_BOUNDARY_DELETION_INTRON) {
	  univdiagonal = push_splice_boundary_deletion(&junctions,&endpoints,&univdiagonals,
						       &nmismatches,&ref_nmismatches,
						       univdiagonal,querylength - tend,splice_distance,adj,
						       sensedir,/*deletion_first_p*/exonlength - (exonpos + adj) > 0,
						       /*direction*/-1,listpool,intlistpool,
						       univcoordlistpool,pathpool);
	  tr_junctions = List_next(tr_junctions);
	  (*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf(">>Output splice-boundary deletion.  Advancing tr_junctions\n");
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	} else if ((splice_distance += adj) == 0) {
	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Junction_list_gc(&junctions,listpool,pathpool);
	  Exon_list_gc(&exons,listpool,transcriptpool);

	  Intlistpool_free_list(&tr_endpoints,intlistpool);
	  Intlistpool_free_list(&tr_nmismatches,intlistpool);

	  return (Path_T) NULL;

	} else {
	  junctions = Listpool_push(junctions,listpool,(void *) JUNCTION_UNSOLVED
				    listpool_trace(__FILE__,__LINE__));
	  debug10(printf(">>Output UNSOLVED splice junction.  Advancing tr_junctions\n"));
	  tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	  Intlist_head_set(nmismatches,-1);
	  Intlist_head_set(tr_nmismatches->rest,-1);
	  univdiagonal -= splice_distance; /* with Junction_adj.  geneminus */
	  evalp = true;
	  unsolvedp = true;
	}

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);

	tstart = tend; trstart = trend;
	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

      } else if (exoni < nexons - 1 && exonpos >= exonlength - 3) {
	/* gminus post-intron: Indel is at or near acceptor splice site */
	debug10(printf("gplus post-intron\n"));
	++exoni;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	if (exonpos < exonlength) {
	  lastchar = 'y';
	} else if (exonpos > exonlength) {
	  lastchar = 'x';
	} else {
	  lastchar = 's';
	}
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp after intron of length %u\n",
		       adj,exonlength - exonpos,splice_distance));

	if (Parameters_fix_splice_boundary_indels_p == true &&
	    adj > 0 && adj <= MAX_TRANSCRIPT_SPLICE_BOUNDARY_DELETION &&
	    splice_distance >= MIN_TRANSCRIPT_SPLICE_BOUNDARY_DELETION_INTRON) {
	  univdiagonal = push_splice_boundary_deletion(&junctions,&endpoints,&univdiagonals,
						       &nmismatches,&ref_nmismatches,
						       univdiagonal,querylength - tend,splice_distance,adj,
						       sensedir,/*deletion_first_p*/false,
						       /*direction*/-1,listpool,intlistpool,
						       univcoordlistpool,pathpool);
	  tr_junctions = List_next(tr_junctions);
	  (*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf(">>Output splice-boundary deletion.  Advancing tr_junctions\n");
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	} else if ((splice_distance += adj) == 0) {
	  Intlistpool_free_list(&endpoints,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				      univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
	  Intlistpool_free_list(&nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Intlistpool_free_list(&ref_nmismatches,intlistpool
				intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
	  Junction_list_gc(&junctions,listpool,pathpool);
	  Exon_list_gc(&exons,listpool,transcriptpool);

	  Intlistpool_free_list(&tr_endpoints,intlistpool);
	  Intlistpool_free_list(&tr_nmismatches,intlistpool);

	  return (Path_T) NULL;

	} else {
	  junctions = Listpool_push(junctions,listpool,(void *) JUNCTION_UNSOLVED
				    listpool_trace(__FILE__,__LINE__));
	  debug10(printf(">>Output UNSOLVED splice junction.  Advancing tr_junctions\n"));
	  tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	  printf("Head of junctions: ");
	  if (tr_junctions != NULL) {
	    Junction_print((Junction_T) List_head(tr_junctions));
	  }
	  printf("\n");
#endif

	  Intlist_head_set(nmismatches,-1);
	  Intlist_head_set(tr_nmismatches->rest,-1);
	  univdiagonal -= splice_distance; /* with Junction_adj.  geneminus */
	  evalp = true;
	  unsolvedp = true;
	}

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);

	tstart = tend; trstart = trend;
	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

      } else {
	/* Indel is within exon */

	/* Advance to next transcript segment */
	trdiagonal = Uintlist_head(trdiagonals);
	tr_endpoints = Intlist_next(tr_endpoints);
	tr_nmismatches = Intlist_next(tr_nmismatches);
	tstart = tend; trstart = trend;

	tend = Intlist_second_value(tr_endpoints);
	trend = (trdiagonal - troffset - querylength) + tend;

	junction = Junction_copy((Junction_T) List_head(tr_junctions),pathpool);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	if (Junction_ninserts(junction) != 0) {
	  printf("For ninserts, revising endpoint %d to %d\n",
		 Intlist_head(endpoints),Intlist_head(endpoints) - Junction_ninserts(junction));
	}
#endif
	/* Needed only for geneminus, to move genomic endpoint before insertion */
	Intlist_head_incr(endpoints,-Junction_ninserts(junction));

#if 0
	if (Junction_type(junction) == DEL_JUNCTION) {
	  /* Need to translate transcriptome-based deletionpos to genomic deletionpos */
	  /* univdiagonal - querylength + (querylength - tstart + Junction_nindels(junction)) */
	  deletionpos = univdiagonal - tstart - Junction_nindels(junction);
	  debug10(printf("At geneminus deletion, computing deletionpos to be %u\n",deletionpos));
	  Junction_set_deletionpos(junction,deletionpos);
	}
#endif
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf(">>(2) Output indel junction: ");
	Junction_print(junction);
	printf("\n");
#endif
	univdiagonal -= Junction_adj(junction); /* geneminus */
	tr_junctions = List_next(tr_junctions);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf("Head of junctions: ");
	if (tr_junctions != NULL) {
	  Junction_print((Junction_T) List_head(tr_junctions));
	}
	printf("\n");
#endif	
	
	/* Do not advance exon information */
      }

    } else {
      /* Transcript segment extends to the next exon */

      /* Revise transcript segment */
      trstart = (trdiagonal - troffset - querylength) + tstart;
      uselength = exend - trstart;
      debug10(printf(">>Output next-exon: univdiagonal %u, %d..%d\n",
		     univdiagonal,tstart,tstart + uselength));

      tstart += uselength; trstart += uselength;
      Intlist_head_set(tr_endpoints,tstart);
      endpoints = Intlistpool_push(endpoints,intlistpool,querylength - tstart
				   intlistpool_trace(__FILE__,__LINE__));
      univdiagonals = Univcoordlistpool_push(univdiagonals,univcoordlistpool,univdiagonal
					     univcoordlistpool_trace(__FILE__,__LINE__));

      if (Intlist_head(tr_nmismatches) == 0) {
	nmismatches = Intlistpool_push(nmismatches,intlistpool,0
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,0
					   intlistpool_trace(__FILE__,__LINE__));
      } else {
	/* Don't know how to divide a non-zero value of nmismatches between the exons */
	nmismatches = Intlistpool_push(nmismatches,intlistpool,-1
				       intlistpool_trace(__FILE__,__LINE__));
	ref_nmismatches = Intlistpool_push(ref_nmismatches,intlistpool,-1
					   intlistpool_trace(__FILE__,__LINE__));
	Intlist_head_set(tr_nmismatches,-1);
	evalp = true;
      }

      if (++exoni < nexons) {
#ifdef COMPUTE_EXONS_INCREMENTALLY
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,/*lastchar*/'s',transcriptpool)
			      listpool_trace(__FILE__,__LINE__));
#endif

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	junction = Junction_new_splice_annotated(splice_distance,sensedir,
						 /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
	(*nintrons) += 1;
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
	printf(">>Output splice junction: ");
	Junction_print(junction);
	printf("\n");
#endif

	univdiagonal -= splice_distance; /* geneminus */

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
#ifdef COMPUTE_EXONS_INCREMENTALLY
	firstchar = 's';
#endif
      }
    }
  }

  debug10(printf(">>Final exoni %d, tcoords %d..%d, trcoords %d..%d\n",
		 exoni,tstart,tend,trstart,trend));
#ifdef COMPUTE_EXONS_INCREMENTALLY
  if (exoni >= nexons - 1) {
    exons = Listpool_push(exons,listpool,
			  Exon_new(firstchar,nexons-1,/*lastchar*/'.',transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
  } else {
    exonpos = trend - exstart;
    if (exonpos == exonlength) {
      lastchar = 's';
    } else {
      lastchar = '.';
    }
    exons = Listpool_push(exons,listpool,
			  Exon_new(firstchar,exoni,lastchar,transcriptpool)
			  listpool_trace(__FILE__,__LINE__));
  }

  exons = List_reverse(exons);
#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("New exons: ");
  Exon_print_list_stdout(exons);
  printf("\n");
#endif
#endif


  /*Do not reverse lists, since they are already in gplus order */

#if defined(DEBUG10) || defined(DEBUG_TRANSCRIPTOME) || defined(DEBUG_TRPATH)
  printf("endpoints: %s\n",Intlist_to_string(endpoints));
  printf("univdiagonals: %s (%s)\n",
	 Univcoordlist_to_string(univdiagonals),Univcoordlist_to_string_offset(univdiagonals,chroffset));
  printf("nmismatches: %s\n",Intlist_to_string(nmismatches));
  printf("junctions: ");
  Junction_print_list(junctions);
  printf("\n");
#endif

  if (Path_endpoints_acceptable_p(endpoints,junctions) == false) {
    debug10(printf("Endpoints not acceptable\n"));
    Intlistpool_free_list(&endpoints,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Univcoordlistpool_free_list(&univdiagonals,univcoordlistpool
				univcoordlistpool_trace(__FILE__,__LINE__)); /* Allocated by Univcoordlistpool_T */
    Intlistpool_free_list(&nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Intlistpool_free_list(&ref_nmismatches,intlistpool
			  intlistpool_trace(__FILE__,__LINE__)); /* Allocated by Intlistpool_T */
    Junction_list_gc(&junctions,listpool,pathpool);
    Exon_list_gc(&exons,listpool,transcriptpool);

    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);

    return (Path_T) NULL;

  } else if (unsolvedp == true) {
    path = Path_convert_simple(chrnum,chroffset,chrhigh,
			       endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,evalp,
			       plusp,first_read_p,/*genestrand*/0,sensedir,querylength,
			       query_compress_fwd,query_compress_rev,listpool,pathpool,method);
    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);

    if ((newpath = Path_solve_junctions(&(*found_score),&(*end_score),path,sensedir,/*genestrand*/0,
					query_compress,query_compress_fwd,query_compress_rev,
					queryseq,queryptr,querylength,

					novel_univdiags_alloc,novel_diagonals_alloc,
					localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					localdb,localdb_nmismatches_allowed,

					this,knownsplicing,
					uintlistpool,intlistpool,univcoordlistpool,listpool,
					pathpool,transcriptpool)) == NULL) {
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      return (Path_T) NULL;

    } else {
      /* Path_solve_junctions called Path_eval_nmatches */
      debug10(printf("Have solved path: "));
      debug10(Path_print(newpath));

      newpath->completep = true;

      overall_trstart = (trdiagonal - troffset - querylength) +
	(querylength - Intlist_last_value(path->endpoints));
      overall_trend = (trdiagonal - troffset - querylength) +
	(querylength - Intlist_head(path->endpoints));
      exons = compute_exons(overall_trstart,overall_trend,exonbounds,nexons,
			    listpool,transcriptpool);

      transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/-1,
				  overall_trstart+/*1-based*/1,overall_trend,
				  /*trstart_overhang*/0,/*trend_overhang*/0,
				  exons,nexons,/*trlength*/exonbounds[nexons-1],
				  transcriptpool);
      /* Path_consolidate and Concordance_tr need to handle paths without valid transcripts */
      newpath->invalid_transcripts = Listpool_push(NULL,listpool,(void *) transcript
						   listpool_trace(__FILE__,__LINE__));
      return newpath;
    }
    
  } else {
#ifdef TEST_EXONS
    exons_simple = compute_exons_simple(overall_trstart,overall_trend,exonbounds,nexons,
					listpool,transcriptpool);
    printf("Correct exons: ");
    Exon_print_list_stdout(exons_simple);
    printf("\n");
    assert(Exon_list_equal(exons,exons_simple));
    Exon_list_gc(&exons_simple,listpool,transcriptpool);
#endif

    path = Path_convert_simple(chrnum,chroffset,chrhigh,
			       endpoints,univdiagonals,nmismatches,ref_nmismatches,junctions,evalp,
			       plusp,first_read_p,/*genestrand*/0,sensedir,querylength,
			       query_compress_fwd,query_compress_rev,listpool,pathpool,method);

    if (trpath->found_score == 0 && path->found_score == 0 &&
	mismatches_all_zero_p(path->nmismatches,path->ref_nmismatches) == true) {
      /* Can skip Path_eval_nmatches only when genomic conversion is also exact. */
      path->found_score = path->score_within_trims = 0;
      path->nmatches = path->ref_nmatches = querylength /*- trpath->found_score*/ - trpath->total_ninserts;
      path->total_splice_prob = path->junction_splice_prob = (*nintrons) * 4.0;
      *found_score = 0;
      *end_score = 0;
    } else {
      debug10(printf("(2) Calling Path_eval_nmatches\n"));
      Path_eval_nmatches(&(*found_score),&(*end_score),path,query_compress_fwd,query_compress_rev);
    }

#if 0    
    /* Not effective.  Pruning transcriptome paths early yields to partial
       trpaths, which can obscure complete genomic paths */
    if (Path_segments_good_p(path) == false) {
      debug10(printf("(2) Path segments are not good, so returning NULL\n"));
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
      return (Path_T) NULL;
    }
#endif  

    overall_trstart = (trdiagonal - troffset - querylength) +
      (querylength - Intlist_last_value(path->endpoints));
    overall_trend = (trdiagonal - troffset - querylength) +
      (querylength - Intlist_head(path->endpoints));
    exons = compute_exons(overall_trstart,overall_trend,exonbounds,nexons,
			  listpool,transcriptpool);

    transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/-1,
				overall_trstart+/*1-based*/1,overall_trend,
				/*trstart_overhang*/0,/*trend_overhang*/0,
				exons,nexons,/*trlength*/exonbounds[nexons-1],
				transcriptpool);
    path->transcripts = Listpool_push(NULL,listpool,(void *) transcript
				      listpool_trace(__FILE__,__LINE__));

    debug10(printf("Have path: "));
    debug10(Path_print(path));

    Intlistpool_free_list(&tr_endpoints,intlistpool);
    Intlistpool_free_list(&tr_nmismatches,intlistpool);
    
    return path;
  }
}


/* tplusp (direction of read along transcript) and
   transcript_genestrand (geneplus/geneminus) (direction of
   transcript along genome) combine to give plusp (direction of read
   along genome) */

void
Trpath_convert_sense (int *found_score, int *end_score,
		      List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		      List_T sense_trpaths, bool first_read_p,
		      Shortread_T queryseq, int querylength,
		      Stage1_T this, Knownsplicing_T knownsplicing,

		      Compress_T query_compress_fwd, Compress_T query_compress_rev,

		      Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
		      unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
		      int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		      Hitlistpool_T hitlistpool, bool crossover_sense_p) {
  List_T p;
  Path_T path, copy;
  T trpath;
  int nintrons;
  char *queryptr;
  /* int ignore_found_score = 0; */

  Univcoord_T chroffset, chrhigh;

  for (p = sense_trpaths; p != NULL; p = List_next(p)) {
    trpath = (T) List_head(p);

    EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,trpath->chrnum);
    /* nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trpath->trnum); */

    if (trpath->transcript_genestrand > 0) {
      queryptr = Shortread_queryuc_ptr(queryseq);
      if ((path = convert_to_path_geneplus(&nintrons,&(*found_score),&(*end_score),
					   trpath,trpath->chrnum,chroffset,chrhigh,
					   /*plusp*/true,first_read_p,/*sensedir*/SENSE_FORWARD,
					   queryseq,queryptr,querylength,this,knownsplicing,
					   /*query_compress*/query_compress_fwd,
					   query_compress_fwd,query_compress_rev,

					   novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					   localdb,localdb_nmismatches_allowed,

					   intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
					   transcriptpool,hitlistpool,trpath->method)) != NULL) {
	/* Path_eval_nmatches(&ignore_found_score,path,query_compress_fwd,query_compress_rev); */
	debug10(Path_print(path));

	path->completep = true;
	*sense_paths_gplus = Hitlist_push(*sense_paths_gplus,hitlistpool,(void *) path
					  hitlistpool_trace(__FILE__,__LINE__));

	if (crossover_sense_p == true && nintrons == 0) {
	  /* Add a copy, without transcript info, to the other sense direction */
	  copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,transcriptpool,hitlistpool);
	  copy->sensedir = SENSE_ANTI;
	  Transcript_list_gc(&copy->transcripts,listpool,transcriptpool);
	  Transcript_list_gc(&copy->invalid_transcripts,listpool,transcriptpool);
	  *antisense_paths_gplus = Hitlist_push(*antisense_paths_gplus,hitlistpool,(void *) copy
						hitlistpool_trace(__FILE__,__LINE__));
	}
      }

    } else {
      queryptr = Shortread_queryrc(queryseq);
      if ((path = convert_to_path_geneminus(&nintrons,&(*found_score),&(*end_score),
					    trpath,trpath->chrnum,chroffset,chrhigh,
					    /*plusp*/false,first_read_p,/*sensedir*/SENSE_FORWARD,
					    queryseq,queryptr,querylength,this,knownsplicing,
					    /*query_compress*/query_compress_rev,
					    query_compress_fwd,query_compress_rev,
					    
					    novel_univdiags_alloc,novel_diagonals_alloc,
					    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					    localdb,localdb_nmismatches_allowed,

					    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
					    transcriptpool,hitlistpool,trpath->method)) != NULL) {
	/* Path_eval_nmatches(&ignore_found_score,path,query_compress_fwd,query_compress_rev); */
	debug10(Path_print(path));

	path->completep = true;
	*sense_paths_gminus = Hitlist_push(*sense_paths_gminus,hitlistpool,(void *) path
					   hitlistpool_trace(__FILE__,__LINE__));

	if (crossover_sense_p == true && nintrons == 0) {
	  /* Add a copy, without transcript info, to the other sense direction */
	  copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,transcriptpool,hitlistpool);
	  copy->sensedir = SENSE_ANTI;
	  Transcript_list_gc(&copy->transcripts,listpool,transcriptpool);
	  Transcript_list_gc(&copy->invalid_transcripts,listpool,transcriptpool);
	  *antisense_paths_gminus = Hitlist_push(*antisense_paths_gminus,hitlistpool,(void *) copy
						hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  return;
}


void
Trpath_convert_antisense (int *found_score, int *end_score,
			  List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,
			  List_T *sense_paths_gplus, List_T *sense_paths_gminus,

			  List_T antisense_trpaths, bool first_read_p,
			  Shortread_T queryseq, int querylength,
			  Stage1_T this, Knownsplicing_T knownsplicing,

			  Compress_T query_compress_fwd, Compress_T query_compress_rev,
			  
			  Univdiag_T *novel_univdiags_alloc, Univcoord_T *novel_diagonals_alloc,
			  unsigned short *localdb_alloc, int16_t *localdb_tally, int16_t *localdb_overflowp,
			  int max_nblocks, Localdb_T localdb, int localdb_nmismatches_allowed,

			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			  Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			  Hitlistpool_T hitlistpool, bool crossover_sense_p) {
  List_T p;
  Path_T path, copy;
  T trpath;
  int nintrons;
  char *queryptr;
  /* int ignore_found_score = 0; */

  Univcoord_T chroffset, chrhigh;

  for (p = antisense_trpaths; p != NULL; p = List_next(p)) {
    trpath = (T) List_head(p);

    EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,trpath->chrnum);
    /* nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trpath->trnum); */

    if (trpath->transcript_genestrand > 0) {
      queryptr = Shortread_queryrc(queryseq);
      if ((path = convert_to_path_geneplus(&nintrons,&(*found_score),&(*end_score),
					   trpath,trpath->chrnum,chroffset,chrhigh,
					   /*plusp*/false,first_read_p,/*sensedir*/SENSE_ANTI,
					   queryseq,queryptr,querylength,this,knownsplicing,
					   /*query_compress*/query_compress_rev,
					   query_compress_fwd,query_compress_rev,

					   novel_univdiags_alloc,novel_diagonals_alloc,
					   localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					   localdb,localdb_nmismatches_allowed,

					   intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
					   transcriptpool,hitlistpool,trpath->method)) != NULL) {
	/* Path_eval_nmatches(&ignore_found_score,path,query_compress_fwd,query_compress_rev); */
	debug10(Path_print(path));

	path->completep = true;
	*antisense_paths_gminus = Hitlist_push(*antisense_paths_gminus,hitlistpool,(void *) path
					      hitlistpool_trace(__FILE__,__LINE__));
	if (crossover_sense_p == true && nintrons == 0) {
	  /* Add a copy, without transcript info, to the other sense direction */
	  copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,transcriptpool,hitlistpool);
	  copy->sensedir = SENSE_FORWARD;
	  Transcript_list_gc(&copy->transcripts,listpool,transcriptpool);
	  Transcript_list_gc(&copy->invalid_transcripts,listpool,transcriptpool);
	  *sense_paths_gminus = Hitlist_push(*sense_paths_gminus,hitlistpool,(void *) copy
					     hitlistpool_trace(__FILE__,__LINE__));
	}
      }

    } else {
      queryptr = Shortread_queryuc_ptr(queryseq);
      if ((path = convert_to_path_geneminus(&nintrons,&(*found_score),&(*end_score),
					    trpath,trpath->chrnum,chroffset,chrhigh,
					    /*plusp*/true,first_read_p,/*sensedir*/SENSE_ANTI,
					    queryseq,queryptr,querylength,this,knownsplicing,
					    /*query_compress*/query_compress_fwd,
					    query_compress_fwd,query_compress_rev,

					    novel_univdiags_alloc,novel_diagonals_alloc,
					    localdb_alloc,localdb_tally,localdb_overflowp,max_nblocks,
					    localdb,localdb_nmismatches_allowed,

					    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
					    transcriptpool,hitlistpool,trpath->method)) != NULL) {
	/* Path_eval_nmatches(&ignore_found_score,path,query_compress_fwd,query_compress_rev); */
	debug10(Path_print(path));

	path->completep = true;
	*antisense_paths_gplus = Hitlist_push(*antisense_paths_gplus,hitlistpool,(void *) path
					       hitlistpool_trace(__FILE__,__LINE__));
	if (crossover_sense_p == true && nintrons == 0) {
	  /* Add a copy, without transcript info, to the other sense direction */
	  copy = Path_copy(path,intlistpool,univcoordlistpool,listpool,pathpool,
			   vectorpool,transcriptpool,hitlistpool);
	  copy->sensedir = SENSE_FORWARD;
	  Transcript_list_gc(&copy->transcripts,listpool,transcriptpool);
	  Transcript_list_gc(&copy->invalid_transcripts,listpool,transcriptpool);
	  *sense_paths_gplus = Hitlist_push(*sense_paths_gplus,hitlistpool,(void *) copy
					    hitlistpool_trace(__FILE__,__LINE__));
	}
      }
    }
  }

  return;
}


void
Trpath_convert_setup (Transcriptome_T transcriptome_in,
		      EF64_T chromosome_ef64_in) {

  transcriptome = transcriptome_in;
  chromosome_ef64 = chromosome_ef64_in;

  return;
}
