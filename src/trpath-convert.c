static char rcsid[] = "$Id: 2d1589d3621db3a00e22fdb481d7ebb666d5c143 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath.h"
#include "trpath-convert.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"
#include "path-eval.h"
#include "path-solve.h"


#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* endpoints_acceptable_p */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif


/* Conversion to Path_T */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif


/* binary_search */
#ifdef DEBUG20
#define debug20(x) x
#else
#define debug20(x)
#endif


static Transcriptome_T transcriptome;
static EF64_T chromosome_ef64;


#define T Trpath_T

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


#ifdef TEST_EXONS
static List_T
compute_exons_simple (int trstart, int trend, int *exonbounds, int nexons,
		      Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T exons;
  char firstchar, lastchar;
  int exoni, exonj, i;

#ifdef DEBUG10
  printf("compute_exons_simple with trstart %d, trend %d\n",trstart,trend);
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
#endif



static Path_T
convert_to_path_geneplus (int *nintrons, int *found_score,
			  Trpath_T trpath, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			  bool plusp, bool first_read_p, int sensedir,
			  Shortread_T queryseq, int querylength,
			  Stage1_T this, Knownsplicing_T knownsplicing, Compress_T query_compress,
			  Compress_T query_compress_fwd, Compress_T query_compress_rev,
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
  char firstchar, lastchar;
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
  int overall_trstart;

  /* Procedure modifies tr_endpoints and tr_nmismatches, which prevenst a second call to convert trpath */
  Intlist_T tr_endpoints = Intlistpool_copy(trpath->endpoints,intlistpool);
  Uintlist_T trdiagonals = trpath->trdiagonals;
  Intlist_T tr_nmismatches = Intlistpool_copy(trpath->nmismatches,intlistpool);
  List_T tr_junctions = trpath->junctions;
  Trcoord_T troffset = trpath->troffset;
  Trcoord_T trdiagonal;


  *nintrons = 0;

  trnum = trpath->trnum;
  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);

  debug0(printf("Entered convert_to_path_geneplus with trpath %p, trnum %d, chrnum %d from method %s\n",
		trpath,trnum,chrnum,Method_string(method)));

#ifdef DEBUG10
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
  overall_trstart = trstart = (trdiagonal - troffset - querylength) + tstart;
  trend = (trdiagonal - troffset - querylength) + tend;
  debug10(printf("Initial tstart %d => trstart %d.  tend %d => trend %d.\n",
		 tstart,trstart,tend,trend));
  assert(trstart < trend);
  if (trstart < 0) {
    debug10(printf("starting before coordinate 0 => aborting\n"));
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
  if (exoni == 0) {
    firstchar = '.';
  } else if (exonpos == 0) {
    firstchar = 's';
  } else {
    firstchar = '.';
  }
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
	lastchar = (exonpos < exonlength) ? 'y' : 'x';
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp before intron of length %u\n",
		       adj,exonlength - (exonpos + adj),splice_distance));
	splice_distance += adj;

	if (splice_distance == 0) {
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
#ifdef DEBUG10
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
	firstchar = 's';

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

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp after intron of length %u\n",
		       adj,exonlength - exonpos,splice_distance));
	splice_distance += adj;

	if (splice_distance == 0) {
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
#ifdef DEBUG10
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
	firstchar = 's';

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
#ifdef DEBUG10
	printf(">>(1) Output indel junction: ");
	Junction_print(junction);
	printf("\n");
#endif
	univdiagonal += Junction_adj(junction); /* geneplus */
	tr_junctions = List_next(tr_junctions);
#ifdef DEBUG10
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
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,/*lastchar*/'s',transcriptpool)
			      listpool_trace(__FILE__,__LINE__));

	splice_distance = exonstarts[exoni] - (exonstarts[exoni-1] + exonlength); /* gplus */
	junction = Junction_new_splice(splice_distance,sensedir,
				       /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
	(*nintrons) += 1;
#ifdef DEBUG10
	printf(">>Output splice junction: ");
	Junction_print(junction);
	printf("\n");
#endif

	univdiagonal += splice_distance; /* geneplus */

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
	firstchar = 's';
      }
    }
  }

  debug10(printf(">>Final exoni %d, tcoords %d..%d, trcoords %d..%d\n",
		 exoni,tstart,tend,trstart,trend));
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
#ifdef DEBUG10
  printf("New exons: ");
  Exon_print_list_stdout(exons);
  printf("\n");
#endif


  /* Reverse lists to put them in gplus order */
  endpoints = Intlist_reverse(endpoints);
  univdiagonals = Univcoordlist_reverse(univdiagonals);
  nmismatches = Intlist_reverse(nmismatches);
  ref_nmismatches = Intlist_reverse(ref_nmismatches);
  junctions = List_reverse(junctions);

#ifdef DEBUG10
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

    if ((newpath = Path_solve_junctions(&(*found_score),path,sensedir,/*genestrand*/0,
					query_compress,query_compress_fwd,query_compress_rev,
					queryseq,querylength,this,knownsplicing,
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
      transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/+1,
				  overall_trstart+/*1-based*/1,trend,/*trstart_overhang*/0,/*trend_overhang*/0,
				  exons,nexons,/*trlength*/exonbounds[nexons-1],
				  transcriptpool);
      /* Path_consolidate and Concordance_tr need to handle paths without valid transcripts */
      newpath->invalid_transcripts = Listpool_push(NULL,listpool,(void *) transcript
						   listpool_trace(__FILE__,__LINE__));
      return newpath;
    }

  } else {
#ifdef TEST_EXONS
    exons_simple = compute_exons_simple(overall_trstart,/*overall_trend*/trend,exonbounds,nexons,
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

    if (trpath->found_score == 0) {
      /* Can skip Path_eval_nmatches */
      path->found_score = path->score_within_trims = 0;
      path->nmatches = path->ref_nmatches = querylength /*- trpath->found_score*/ - trpath->total_ninserts;
      path->total_splice_prob = path->junction_splice_prob = (*nintrons) * 4.0;
      if (*found_score > 0) {
	*found_score = 0;
      }
    } else {
      Path_eval_nmatches(&(*found_score),path,query_compress_fwd,query_compress_rev);
    }

    transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/+1,
				overall_trstart+/*1-based*/1,trend,/*trstart_overhang*/0,/*trend_overhang*/0,
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
convert_to_path_geneminus (int *nintrons, int *found_score,
			   Trpath_T trpath, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh,
			   bool plusp, bool first_read_p, int sensedir,
			   Shortread_T queryseq, int querylength,
			   Stage1_T this, Knownsplicing_T knownsplicing, Compress_T query_compress,
			   Compress_T query_compress_fwd, Compress_T query_compress_rev,
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
  char firstchar, lastchar;
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
  int overall_trstart;

  /* Procedure modifies tr_endpoints and tr_nmismatches, which prevents a second call to convert trpath */
  Intlist_T tr_endpoints = Intlistpool_copy(trpath->endpoints,intlistpool);
  Uintlist_T trdiagonals = trpath->trdiagonals;
  Intlist_T tr_nmismatches = Intlistpool_copy(trpath->nmismatches,intlistpool);
  List_T tr_junctions = trpath->junctions;
  Trcoord_T troffset = trpath->troffset;
  Trcoord_T trdiagonal;


  *nintrons = 0;

  trnum = trpath->trnum;
  nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trnum);

  debug0(printf("Entered convert_to_path_geneminus with trpath %p, trnum %d, chrnum %d from method %s\n",
		trpath,trnum,chrnum,Method_string(method)));

#ifdef DEBUG10
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
  overall_trstart = trstart = (trdiagonal - troffset - querylength) + tstart;
  trend = (trdiagonal - troffset - querylength) + tend;
  debug10(printf("Initial tstart %d => trstart %d.  tend %d => trend %d.\n",
		 tstart,trstart,tend,trend));
  assert(trstart < trend);
  if (trstart < 0) {
    debug10(printf("starting before coordinate 0 => aborting\n"));
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
  if (exoni == 0) {
    firstchar = '.';
  } else if (exonpos == 0) {
    firstchar = 's';
  } else {
    firstchar = '.';
  }
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
	lastchar = (exonpos < exonlength) ? 'y' : 'x';
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,lastchar,transcriptpool)
			      listpool_trace(__FILE__,__LINE__));

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp before intron of length %u\n",
		       adj,exonlength - (exonpos + adj),splice_distance));
	splice_distance += adj;

	if (splice_distance == 0) {
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
#ifdef DEBUG10
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
	firstchar = 's';

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

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	adj = Junction_adj((Junction_T) List_head(tr_junctions));
	debug10(printf("Indel of adj %d is %d bp after intron of length %u\n",
		       adj,exonlength - exonpos,splice_distance));
	splice_distance += adj;

	if (splice_distance == 0) {
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
#ifdef DEBUG10
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
	firstchar = 's';

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
#ifdef DEBUG10
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
#ifdef DEBUG10
	printf(">>(2) Output indel junction: ");
	Junction_print(junction);
	printf("\n");
#endif
	univdiagonal -= Junction_adj(junction); /* geneminus */
	tr_junctions = List_next(tr_junctions);
#ifdef DEBUG10
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
	exons = Listpool_push(exons,listpool,
			      Exon_new(firstchar,exoni-1,/*lastchar*/'s',transcriptpool)
			      listpool_trace(__FILE__,__LINE__));

	splice_distance = (exonstarts[exoni-1] - exonlength) - exonstarts[exoni]; /* gminus */
	junction = Junction_new_splice(splice_distance,sensedir,
				       /*donor_prob*/2.0,/*acceptor_prob*/2.0,pathpool);
	junctions = Listpool_push(junctions,listpool,(void *) junction
				  listpool_trace(__FILE__,__LINE__));
	(*nintrons) += 1;
#ifdef DEBUG10
	printf(">>Output splice junction: ");
	Junction_print(junction);
	printf("\n");
#endif

	univdiagonal -= splice_distance; /* geneminus */

	/* Advance exon information */
	exstart = exend;
	exend = exonbounds[exoni];
	exonlength = exend - exstart;
	firstchar = 's';
      }
    }
  }

  debug10(printf(">>Final exoni %d, tcoords %d..%d, trcoords %d..%d\n",
		 exoni,tstart,tend,trstart,trend));
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
#ifdef DEBUG10
  printf("New exons: ");
  Exon_print_list_stdout(exons);
  printf("\n");
#endif


  /*Do not reverse lists, since they are already in gplus order */

#ifdef DEBUG10
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

    if ((newpath = Path_solve_junctions(&(*found_score),path,sensedir,/*genestrand*/0,
					query_compress,query_compress_fwd,query_compress_rev,
					queryseq,querylength,this,knownsplicing,
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
      transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/-1,
				  overall_trstart+/*1-based*/1,trend,/*trstart_overhang*/0,/*trend_overhang*/0,
				  exons,nexons,/*trlength*/exonbounds[nexons-1],
				  transcriptpool);
      /* Path_consolidate and Concordance_tr need to handle paths without valid transcripts */
      newpath->invalid_transcripts = Listpool_push(NULL,listpool,(void *) transcript
						   listpool_trace(__FILE__,__LINE__));
      return newpath;
    }
    
  } else {
#ifdef TEST_EXONS
    exons_simple = compute_exons_simple(overall_trstart,/*overall_trend*/trend,exonbounds,nexons,
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

    if (trpath->found_score == 0) {
      /* Can skip Path_eval_nmatches */
      path->found_score = path->score_within_trims = 0;
      path->nmatches = path->ref_nmatches = querylength /*- trpath->found_score*/ - trpath->total_ninserts;
      path->total_splice_prob = path->junction_splice_prob = (*nintrons) * 4.0;
    } else {
      Path_eval_nmatches(&(*found_score),path,query_compress_fwd,query_compress_rev);
    }

    transcript = Transcript_new(/*num*/trnum,/*transcript_genestrand*/-1,
				overall_trstart+/*1-based*/1,trend,/*trstart_overhang*/0,/*trend_overhang*/0,
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
Trpath_convert_sense (int *found_score,
		      List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		      List_T sense_trpaths, bool first_read_p,
		      Shortread_T queryseq, int querylength,
		      Stage1_T this, Knownsplicing_T knownsplicing,

		      Compress_T query_compress_fwd, Compress_T query_compress_rev,

		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		      Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
		      Hitlistpool_T hitlistpool, bool crossover_sense_p) {
  List_T p;
  Path_T path, copy;
  T trpath;
  int nintrons;
  /* int ignore_found_score = 0; */

  Univcoord_T chroffset, chrhigh;

  for (p = sense_trpaths; p != NULL; p = List_next(p)) {
    trpath = (T) List_head(p);

    EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,trpath->chrnum);
    /* nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trpath->trnum); */

    if (trpath->transcript_genestrand > 0) {
      if ((path = convert_to_path_geneplus(&nintrons,&(*found_score),
					   trpath,trpath->chrnum,chroffset,chrhigh,
					   /*plusp*/true,first_read_p,/*sensedir*/SENSE_FORWARD,
					   queryseq,querylength,this,knownsplicing,
					   /*query_compress*/query_compress_fwd,
					   query_compress_fwd,query_compress_rev,
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
      if ((path = convert_to_path_geneminus(&nintrons,&(*found_score),
					    trpath,trpath->chrnum,chroffset,chrhigh,
					    /*plusp*/false,first_read_p,/*sensedir*/SENSE_FORWARD,
					    queryseq,querylength,this,knownsplicing,
					    /*query_compress*/query_compress_rev,
					    query_compress_fwd,query_compress_rev,
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
Trpath_convert_antisense (int *found_score,
			  List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,
			  List_T *sense_paths_gplus, List_T *sense_paths_gminus,

			  List_T antisense_trpaths, bool first_read_p,
			  Shortread_T queryseq, int querylength,
			  Stage1_T this, Knownsplicing_T knownsplicing,

			  Compress_T query_compress_fwd, Compress_T query_compress_rev,
			  
			  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			  Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
			  Pathpool_T pathpool, Vectorpool_T vectorpool, Transcriptpool_T transcriptpool,
			  Hitlistpool_T hitlistpool, bool crossover_sense_p) {
  List_T p;
  Path_T path, copy;
  T trpath;
  int nintrons;
  /* int ignore_found_score = 0; */

  Univcoord_T chroffset, chrhigh;

  for (p = antisense_trpaths; p != NULL; p = List_next(p)) {
    trpath = (T) List_head(p);

    EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,trpath->chrnum);
    /* nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,trpath->trnum); */

    if (trpath->transcript_genestrand > 0) {
      if ((path = convert_to_path_geneplus(&nintrons,&(*found_score),
					   trpath,trpath->chrnum,chroffset,chrhigh,
					   /*plusp*/false,first_read_p,/*sensedir*/SENSE_ANTI,
					   queryseq,querylength,this,knownsplicing,
					   /*query_compress*/query_compress_rev,
					   query_compress_fwd,query_compress_rev,
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
      if ((path = convert_to_path_geneminus(&nintrons,&(*found_score),
					    trpath,trpath->chrnum,chroffset,chrhigh,
					    /*plusp*/true,first_read_p,/*sensedir*/SENSE_ANTI,
					    queryseq,querylength,this,knownsplicing,
					    /*query_compress*/query_compress_fwd,
					    query_compress_fwd,query_compress_rev,
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
