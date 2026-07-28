#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath-solve.h"
#include "trpath-eval.h"
#include "path-eval.h"

#include "assert.h"
#include "mem.h"
#include "spliceends.h"
#include "genomebits_indel.h"
#include "genomebits_count.h"
#include "genomebits_trim.h"

static Transcriptome_T transcriptome;
static Genomebits_T transcriptomebits;
static EF64_T transcript_ef64;

static int max_insertionlen;
static int max_deletionlen;
static int index1part_tr;


#define MIN_SUPPORT_INDEL 6	/* Also defined in kmer-search.c */

/* Turns out to worsen results */
/* #define END_SUPPORT 18		-- Require transcript alignments to go this far to ends */


#ifdef CHECK_ASSERTIONS
#define CHECK_NMISMATCHES 1
#endif


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* endpoints_acceptable_p */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

#ifdef DEBUG13
#define debug13(x) x
#else
#define debug13(x)
#endif


#define T Trpath_T


T
Trpath_solve_from_trdiagonal (int *found_score, bool *partialp, int sufficient_score,
			      Trcoord_T trdiagonal, int tstart, int tend,
			      
			      Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			      Compress_T query_compress_tr, bool tplusp, int querylength,
			      int *mismatch_positions_alloc, Chrnum_T chrnum, int transcript_genestrand,

			      bool paired_end_p, bool first_read_p, Indelinfo_T indelinfo,
			      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			      Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			      Method_T method) {

  int score, nmismatches, ref_nmismatches;
  int trimmed_tstart, trimmed_tend;
  /* int nmismatches_to_trimpos; */

  debug(printf("Entered Trpath_solve_from_trdiagonal with trnum %d, plusp %d: %d..%d\n",
	       trnum,tplusp,tstart,tend));

#if 0
  if (tstart == 0) {
    trimmed_tstart = 0;
  } else {
    trimmed_tstart = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_tr,
					    /*bits*/transcriptomebits,(Univcoord_T) trdiagonal,querylength,
					    /*pos5*/0,/*pos3*/tend,tplusp,/*genestrand*/0);
  }

  if (tend == querylength) {
    trimmed_tend = querylength;
  } else {
    trimmed_tend = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_tr,
					/*bits*/transcriptomebits,(Univcoord_T) trdiagonal,querylength,
					/*pos5*/tstart,/*pos3*/querylength,tplusp,/*genestrand*/0);
  }
#else
  /* Genomebits_trim already called by caller */
  trimmed_tstart = tstart;
  trimmed_tend = tend;
#endif

  nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,transcriptomebits,query_compress_tr,
						      trdiagonal,querylength,
						      /*pos5*/trimmed_tstart,/*pos3*/trimmed_tend,tplusp,/*genestrand*/0);

  if (trimmed_tstart <= 8 && querylength - trimmed_tend <= 8) {
    score = trimmed_tstart + (querylength - trimmed_tend) + nmismatches;
    if (0 && score > sufficient_score) {
      /* Previously discarded these, but we need them for cases with low sequence quality */
      return (Trpath_T) NULL;
    } else {
      /* Reasonably good alignment */
      if (score < *found_score) {
	*found_score = score;
      }
      
      /* Makes no difference if we use index1part_tr or 1 */
      if (trimmed_tstart > index1part_tr) {
	*partialp = true;
      } else if (trimmed_tend < querylength - index1part_tr) {
	*partialp = true;
      } else {
	*partialp = false;
      }
      *partialp = false;
	 
      return Trpath_new_exact(trdiagonal,trimmed_tstart,trimmed_tend,
			      nmismatches,tplusp,trnum,troffset,trhigh,
			      chrnum,transcript_genestrand,
			      intlistpool,uintlistpool,trpathpool,score,method);
    }

  } else {
    return Trpath_solve_from_diagonals(&(*found_score),&(*partialp),sufficient_score,/*middle_trdiagonal*/trdiagonal,
				       trimmed_tstart,trimmed_tend,nmismatches,
				       /*qstart_trdiag*/NULL,/*qend_trdiag*/NULL,
				       tplusp,querylength,query_compress_tr,
				       mismatch_positions_alloc,trnum,troffset,trhigh,
				       chrnum,transcript_genestrand,paired_end_p,first_read_p,
				       indelinfo,intlistpool,uintlistpool,listpool,trpathpool,
				       pathpool,method);
  }
}



/* Modifies path */
static void
attach_unknown_tstart (T trpath, Trcoord_T low_trdiagonal, int low_tstart,
		       Trcoord_T trhigh, int querylength, Indelinfo_T indelinfo,
		       Compress_T query_compress_tr, bool want_lowest_coordinate_tr_p,
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Listpool_T listpool,
		       Pathpool_T pathpool) {
  Trcoord_T trdiagonal;
  int tend;
  int nindels, indel_pos;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j;
  int nmismatches, ref_nmismatches;
#ifdef DEBUG13
  int tstart;
#endif
  
  /* Do not need to call Spliceends_qstart_trim, because
     Genomebits_indel_solve_low has already computed trimpos, provided
     as low_tstart */

  trdiagonal = Uintlist_head(trpath->trdiagonals);

  /* Assume that left+qend gives a coordinate within genome */
  tend = Intlist_second_value(trpath->endpoints) /*+ ninserts*/;

#ifdef DEBUG13
  if (trpath->junctions == NULL) {
    /* ninserts = 0; */
  } else {
    /* ninserts = Junction_ninserts(List_head(trpath->junctions)); */
  }

  tstart = Intlist_head(trpath->endpoints) /*+ ninserts*/;
  printf("Entering attach_unknown_tstart with low_trdiagonal %u, low_tstart %d, and trdiagonal %u %d..%d (diff %d)\n",
	 low_trdiagonal - trpath->troffset,low_tstart,trdiagonal - trpath->troffset,
	 tstart,tend,trdiagonal - low_trdiagonal);
#endif

  if (low_tstart >= tend) {
    debug13(printf("Does not add to start of path: low_qstart %d >= qend %d\n",low_tstart,tend));

  } else if (low_trdiagonal == trdiagonal) {
    if (low_tstart >= Intlist_head(trpath->endpoints)) {
      debug13(printf("Mismatch fails, since new endpoint %d >= old endpoint %d\n",low_tstart,Intlist_head(trpath->endpoints)));
    } else {
      /* Mismatch: Revise the endpoint */
      debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(trpath->endpoints),low_tstart));

      /* Determine nmismatches */
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,transcriptomebits,query_compress_tr,
							  trdiagonal,querylength,/*pos5*/low_tstart,/*pos3*/tend,
							  trpath->tplusp,/*genestrand*/0);
      /* debug13(printf("Counting mismatches from %d to %d => %d (%d ref)\n",
	 low_tstart,tend,nmismatches,ref_nmismatches)); */
      
      Intlist_head_set(trpath->nmismatches,nmismatches);
      Intlist_head_set(trpath->endpoints,low_tstart);
    }

  } else if (low_trdiagonal > trdiagonal + max_insertionlen) {
    /* Impossible */
    debug13(printf("Impossible\n"));

  } else if (low_trdiagonal > trdiagonal) {
    /* (A) Insertion */
    nindels = low_trdiagonal - trdiagonal;
    debug13(printf("Insertion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal_i*/(Univcoord_T) low_trdiagonal,/*indels*/+nindels,
						    trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						    /*pos5*/low_tstart,/*pos3*/tend,querylength,
						    indelinfo,trpath->tplusp,/*genestrand*/0,
						    want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Insertion fails\n"));
      
    } else {
#if 0
      supporti = indel_pos - low_tstart;
      supportj = qend - (indel_pos + nindels);
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      }
#endif
      debug13(printf("nmismatches: %d and %d\n",nmismatches_i,nmismatches_j));
      debug13(printf("(3) attach_indel_tstart is modifying trpath %p\n",trpath));
      Intlist_head_set(trpath->endpoints,indel_pos);
      trpath->endpoints = Intlistpool_push(trpath->endpoints,intlistpool,low_tstart
					   intlistpool_trace(__FILE__,__LINE__));
      trpath->junctions = Listpool_push(trpath->junctions,listpool,
					(void *) Junction_new_insertion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
      /* For tstart, push j first, then push i */
      Intlist_head_set(trpath->nmismatches,nmismatches_j);
      trpath->nmismatches = Intlistpool_push(trpath->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
	
      trpath->trdiagonals = Uintlistpool_push(trpath->trdiagonals,uintlistpool,low_trdiagonal
					      uintlistpool_trace(__FILE__,__LINE__));
      debug13(printf("Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		     low_tstart,tend,indel_pos,nindels,nmismatches_i,nmismatches_j));
    }
    
  } else if (low_trdiagonal + max_deletionlen >= trdiagonal) {
    /* (B) Deletion (or short intron) */
    nindels = trdiagonal - low_trdiagonal;
    debug13(printf("Deletion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/(Univcoord_T) low_trdiagonal,/*indels*/-nindels,
						   trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						   /*pos5*/low_tstart,/*pos3*/tend,querylength,
						   indelinfo,trpath->tplusp,/*genestrand*/0,
						   want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Deletion fails\n"));
	  
    } else {
#if 0
      supporti = indel_pos - low_tstart;
      supportj = qend - indel_pos;
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      }
#endif

      assert(nindels >= 0);
      /* deletionpos = (low_trdiagonal - querylength) + indel_pos; */
      trpath->junctions = Listpool_push(trpath->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));

      debug13(printf("nmismatches: %d and %d\n",nmismatches_i,nmismatches_j));
      debug13(printf("(4) attach_unknown_tstart is modifying trpath %p\n",trpath));
      Intlist_head_set(trpath->endpoints,indel_pos);
      trpath->endpoints = Intlistpool_push(trpath->endpoints,intlistpool,low_tstart
					   intlistpool_trace(__FILE__,__LINE__));
	
      /* For tstart, push j first, then push i */
      Intlist_head_set(trpath->nmismatches,nmismatches_j);
      trpath->nmismatches = Intlistpool_push(trpath->nmismatches,intlistpool,nmismatches_i
					     intlistpool_trace(__FILE__,__LINE__));
      
      trpath->trdiagonals = Uintlistpool_push(trpath->trdiagonals,uintlistpool,low_trdiagonal
					       uintlistpool_trace(__FILE__,__LINE__));
      debug13(printf("Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		     low_tstart,tend,indel_pos,nindels,nmismatches_i,nmismatches_j));
    }
  }

  return;
}


/* Modifies trpath */
static void
attach_unknown_tend (T trpath, Trcoord_T high_trdiagonal, int high_tend,
		     Trcoord_T trhigh, int querylength, Indelinfo_T indelinfo,
		     Compress_T query_compress_tr, bool want_lowest_coordinate_tr_p, 
		     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Listpool_T listpool,
		     Pathpool_T pathpool) {
  Trcoord_T trdiagonal;
  int tstart, ninserts;
  int nindels, indel_pos;
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j;
  int nmismatches, ref_nmismatches;
#ifdef DEBUG13
  int tend;
#endif

  /* Do not need to call Spliceends_qend_trim, because
     Genomebits_indel_solve_high has already computed trimpos,
     provided as high_qend*/

  trdiagonal = Uintlist_head(trpath->trdiagonals);

#if 0
  ninserts = Junction_total_ninserts(trpath->junctions);
#else
  if (trpath->junctions == NULL) {
    ninserts = 0;
  } else {
    ninserts = Junction_ninserts(List_head(trpath->junctions));
  }
#endif

  /* Assume that left+qstart gives a coordinate within genome */
  tstart = Intlist_second_value(trpath->endpoints) + ninserts;

#ifdef DEBUG13
  tend = Intlist_head(trpath->endpoints) /*+ ninserts*/;
  printf("Entering attach_unknown_tend with trdiagonal %u %d..%d and high_trdiagonal %u, high_tend %d (diff %d)\n",
	 trdiagonal - trpath->troffset,tstart,tend,high_trdiagonal - trpath->troffset,
	 high_tend,high_trdiagonal - trdiagonal);
#endif


  if (tstart >= high_tend) {
    debug13(printf("Does not add to end of path: qstart %d >= high_qend %d\n",tstart,high_tend));
    
  } else if (high_trdiagonal == trdiagonal) {
    if (high_tend <= Intlist_head(trpath->endpoints)) {
      debug13(printf("Mismatch fails, since new endpoint %d <= old endpoint %d\n",high_tend,Intlist_head(trpath->endpoints)));

    } else {
      /* Mismatch: Revise the endpoint */
      debug13(printf("Mismatch extends from %d to %d\n",Intlist_head(trpath->endpoints),high_tend));

      /* Determine nmismatches */
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,transcriptomebits,transcriptomebits,query_compress_tr,
							  trdiagonal,querylength,
							  /*pos5*/Intlist_second_value(trpath->endpoints) + ninserts,
							  /*pos3*/high_tend,trpath->tplusp,/*genestrand*/0);
      /* debug13(printf("Counting mismatches from %d to %d => %d (%d ref)\n",
	 Intlist_head(Intlist_next(trpath->endpoints)),high_tend,nmismatches,ref_nmismatches)); */
      Intlist_head_set(trpath->nmismatches,nmismatches);
      Intlist_head_set(trpath->endpoints,high_tend);
    }

  } else if (high_trdiagonal + max_insertionlen < trdiagonal) {
    /* Impossible */
    debug13(printf("Impossible\n"));

  } else if (high_trdiagonal < trdiagonal) {
    /* (A) Insertion */
    nindels = trdiagonal - high_trdiagonal;
    debug13(printf("Insertion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal_i*/(Univcoord_T) trdiagonal,/*indels*/+nindels,
						    trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						    /*pos5*/tstart,/*pos3*/high_tend,querylength,
						    indelinfo,trpath->tplusp,/*genestrand*/0,
						    want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Insertion fails\n"));

    } else {
#if 0
      supporti = indel_pos - tstart;
      supportj = high_qend - (indel_pos + nindels);
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      }
#endif
      debug13(printf("nmismatches: %d and %d\n",nmismatches_i,nmismatches_j));
      debug13(printf("(3) attach_unknown_tend is modifying trpath %p\n",trpath));
      Intlist_head_set(trpath->endpoints,indel_pos);
      trpath->endpoints = Intlistpool_push(trpath->endpoints,intlistpool,high_tend
					   intlistpool_trace(__FILE__,__LINE__));
      trpath->junctions = Listpool_push(trpath->junctions,listpool,
					(void *) Junction_new_insertion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	
      /* For tend, push i first, then push j */
      Intlist_head_set(trpath->nmismatches,nmismatches_i);
      trpath->nmismatches = Intlistpool_push(trpath->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	
      trpath->trdiagonals = Uintlistpool_push(trpath->trdiagonals,uintlistpool,high_trdiagonal
					      uintlistpool_trace(__FILE__,__LINE__));
      debug13(printf("Insertion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		     tstart,high_tend,indel_pos,nindels,nmismatches_i,nmismatches_j));
    }

  } else if (high_trdiagonal <= trdiagonal + max_deletionlen) {
    /* (B) Deletion (or short intron) */
    nindels = high_trdiagonal - trdiagonal;
    debug13(printf("Deletion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/(Univcoord_T) trdiagonal,/*indels*/-nindels,
						   trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						   /*pos5*/tstart,/*pos3*/high_tend,querylength,
						   indelinfo,trpath->tplusp,/*genestrand*/0,
						   want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Deletion fails\n"));
      
    } else {
#if 0
      supporti = indel_pos - tstart;
      supportj = high_qend - indel_pos;
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      }
#endif

      assert(nindels >= 0);
      /* deletionpos = (trdiagonal - querylength) + indel_pos; */
      trpath->junctions = Listpool_push(trpath->junctions,listpool,
					(void *) Junction_new_deletion(nindels,pathpool)
					listpool_trace(__FILE__,__LINE__));
	  
      debug13(printf("nmismatches: %d and %d\n",nmismatches_i,nmismatches_j));
      debug13(printf("(4) attach_unknown_tend is modifying trpath %p\n",trpath));
      Intlist_head_set(trpath->endpoints,indel_pos);
      trpath->endpoints = Intlistpool_push(trpath->endpoints,intlistpool,high_tend
					   intlistpool_trace(__FILE__,__LINE__));
	  
      /* For qend, push i first, then push j */
      Intlist_head_set(trpath->nmismatches,nmismatches_i);
      trpath->nmismatches = Intlistpool_push(trpath->nmismatches,intlistpool,nmismatches_j
					     intlistpool_trace(__FILE__,__LINE__));
	
      trpath->trdiagonals = Uintlistpool_push(trpath->trdiagonals,uintlistpool,high_trdiagonal
					      uintlistpool_trace(__FILE__,__LINE__));
      debug13(printf("Deletion in range %d..%d is at %d with %d indels and nmismatches %d+%d\n",
		     tstart,high_tend,indel_pos,nindels,nmismatches_i,nmismatches_j));
    }
  }

  return;
}


/* The same as Path_endpoints_acceptable_p */
static bool
Trpath_endpoints_acceptable_p (Intlist_T endpoints, List_T junctions) {
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
    if (last_endpoint + Junction_ninserts(junction) >= Intlist_second_value(p)) {
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

static T
combine_leftright_trpaths (T tstart_trpath, T tend_trpath, Compress_T query_compress_tr, bool tplusp,
			   int sufficient_score, int querylength,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Listpool_T listpool,
			   Trpathpool_T trpathpool, Pathpool_T pathpool, Method_T method) {

  T trpath = NULL;
  Intlist_T tr_endpoints, q;
  Uintlist_T trdiagonals, u;
  Intlist_T nmismatches, s;
  List_T junctions, j;

  int tstart1, tend1, tstart2, tend2, ninserts1, ninserts2;
  int middle_nmismatches;
  int score = 0, total_ninserts = 0;
  bool score_knownp = true;


#ifdef DEBUG13
  printf("\n");
  printf("*** Entered combine_leftright_trpaths\n");
#endif

  debug13(printf("++ Tstart/left trpath %p: ",tstart_trpath));
  debug13(Trpath_print(tstart_trpath));
  debug13(printf("++ Tend/right path %p: ",tend_trpath));
  debug13(Trpath_print(tend_trpath));
  debug13(printf("\n"));

  /* Combine tstart_trpath with tend_trpath */
  /* If either list is NULL, must have obtained an unacceptable result */
	
  /* tr_endpoints = (Intlist_T) NULL; -- Initialized with first push */
  trdiagonals = (Uintlist_T) NULL;
  nmismatches = (Intlist_T) NULL;
  junctions = (List_T) NULL;
	
  tend1 = Intlist_last_value(tstart_trpath->endpoints);
  tstart1 = Intlist_penultimate_value(tstart_trpath->endpoints);
	
  tstart2 = Intlist_last_value(tend_trpath->endpoints);
  tend2 = Intlist_penultimate_value(tend_trpath->endpoints);
	
  q = tstart_trpath->endpoints;
  u = tstart_trpath->trdiagonals;
  s = tstart_trpath->nmismatches;
  j = tstart_trpath->junctions;
	
  tr_endpoints = Intlistpool_push(NULL,intlistpool,Intlist_head(q)
				  intlistpool_trace(__FILE__,__LINE__));
  q = Intlist_next(q);
  ninserts1 = 0;
  while (j != NULL) {
    tr_endpoints = Intlistpool_push(tr_endpoints,intlistpool,Intlist_head(q)
				    intlistpool_trace(__FILE__,__LINE__));
    trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,Uintlist_head(u)
				     uintlistpool_trace(__FILE__,__LINE__));
    nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(s)
				   intlistpool_trace(__FILE__,__LINE__));
    score += Intlist_head(nmismatches);

    junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(j),pathpool)
			      listpool_trace(__FILE__,__LINE__));
    ninserts1 = Junction_ninserts((Junction_T) List_head(junctions));
    total_ninserts += ninserts1;

    q = Intlist_next(q);
    u = Uintlist_next(u);
    s = Intlist_next(s);
    j = List_next(j);
  }
	
  /* Reached middle_trdiagonal */
  tstart1 = Intlist_head(tr_endpoints);
  tend1 = Intlist_head(q);
	
  Trpath_reverse(tend_trpath,/*expect_fwd_p*/true);
	
  tstart2 = Intlist_head(tend_trpath->endpoints);
  tend2 = Intlist_second_value(tend_trpath->endpoints);
  ninserts2 = 0;
  if (tend2 <= tstart1) {
    /* No overlap */
    debug13(printf("++ Combined trpath not possible, due to lack of overlap between %d..%d and %d..%d\n",
		   tstart1,tend1,tstart2,tend2));
    Junction_list_gc(&junctions,listpool,pathpool);

  } else {
    tr_endpoints = Intlistpool_push(tr_endpoints,intlistpool,tend2
				    intlistpool_trace(__FILE__,__LINE__));
    trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,Uintlist_head(u)
				    uintlistpool_trace(__FILE__,__LINE__));

    if (tstart1 + ninserts1 == tstart2 + ninserts2 && tend1 == tend2) {
      /* Take nmismatches from either tstart_trpath or tend_trpath (if available) */
      if ((middle_nmismatches = Intlist_head(s)) == -1) {
	middle_nmismatches = Intlist_head(tend_trpath->nmismatches);
      }
	  
    } else if (tstart1 + ninserts1 == tstart2 + ninserts2) {
      /* Take nmismatches from either tend_trpath */
      middle_nmismatches = Intlist_head(tend_trpath->nmismatches);
    } else if (tend1 == tend2) {
      /* Take nmismatches from either tstart_trpath */
      middle_nmismatches = Intlist_head(s);
    } else {
      middle_nmismatches = -1;
      score_knownp = false;
    }
    nmismatches = Intlistpool_push(nmismatches,intlistpool,middle_nmismatches
				   intlistpool_trace(__FILE__,__LINE__));
    score += Intlist_head(nmismatches);
	
    q = Intlist_next(tend_trpath->endpoints);
    u = tend_trpath->trdiagonals;
    s = tend_trpath->nmismatches;
    j = tend_trpath->junctions;
	
    while (j != NULL) {
      junctions = Listpool_push(junctions,listpool,(void *) Junction_copy((Junction_T) List_head(j),pathpool)
				listpool_trace(__FILE__,__LINE__));
      ninserts2 = Junction_ninserts((Junction_T) List_head(junctions));
      total_ninserts += ninserts2;

      q = Intlist_next(q);
      u = Uintlist_next(u);
      s = Intlist_next(s);
      tr_endpoints = Intlistpool_push(tr_endpoints,intlistpool,Intlist_head(q)
				      intlistpool_trace(__FILE__,__LINE__));
      trdiagonals = Uintlistpool_push(trdiagonals,uintlistpool,Uintlist_head(u)
				      uintlistpool_trace(__FILE__,__LINE__));
      nmismatches = Intlistpool_push(nmismatches,intlistpool,Intlist_head(s)
				     intlistpool_trace(__FILE__,__LINE__));
      score += Intlist_head(nmismatches);

      j = List_next(j);
    }

    trdiagonals = Uintlist_reverse(trdiagonals);
    nmismatches = Intlist_reverse(nmismatches);
    tr_endpoints = Intlist_reverse(tr_endpoints);
    junctions = List_reverse(junctions);

    if (Trpath_endpoints_acceptable_p(tr_endpoints,junctions) == false) {
      debug13(printf("++ Combined trpath not possible, due to unacceptable endpoints\n"));
      Intlistpool_free_list(&nmismatches,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      Uintlistpool_free_list(&trdiagonals,uintlistpool
			     uintlistpool_trace(__FILE__,__LINE__));
      Intlistpool_free_list(&tr_endpoints,intlistpool
			    intlistpool_trace(__FILE__,__LINE__));
      Junction_list_gc(&junctions,listpool,pathpool);

    } else {
      if (score_knownp == false) {
	/* Also modifies nmismatches */
	score = Trpath_eval_score(trdiagonals,tr_endpoints,nmismatches,junctions,
				  query_compress_tr,querylength,tplusp);
      }
      
      if (0 && score > sufficient_score) {
	/* Previously discarded these, but we need them for cases with low sequence quality */
	debug13(printf("++ Combined trpath has poor score %d > sufficient score %d\n",score,sufficient_score));
	Intlistpool_free_list(&nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Uintlistpool_free_list(&trdiagonals,uintlistpool
			       uintlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&tr_endpoints,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Junction_list_gc(&junctions,listpool,pathpool);
	
#if 0
      } else if (Intlist_head(tr_endpoints) > END_SUPPORT || querylength - Intlist_last_value(tr_endpoints) > END_SUPPORT) {
	/* Want pretty good alignments, especially since we pair up
	   ends by transcript, which can preclude finding better
	   genomic alignments */
	debug13(printf("++ Combined trpath has poor endpoints %s\n",Intlist_to_string(tr_endpoints)));
	Intlistpool_free_list(&nmismatches,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Uintlistpool_free_list(&trdiagonals,uintlistpool
			       uintlistpool_trace(__FILE__,__LINE__));
	Intlistpool_free_list(&tr_endpoints,intlistpool
			      intlistpool_trace(__FILE__,__LINE__));
	Junction_list_gc(&junctions,listpool,pathpool);
#endif

      } else {
	trpath = Trpath_create(tr_endpoints,trdiagonals,nmismatches,junctions,
			       tplusp,tstart_trpath->trnum,tstart_trpath->troffset,tstart_trpath->trhigh,
			       tstart_trpath->chrnum,tstart_trpath->transcript_genestrand,
			       trpathpool,score,total_ninserts,method);

	debug13(printf("++ Combined trpath %p: ",trpath));
	debug13(Trpath_print(trpath));
	debug13(printf("\n"));
      }
    }
  }
 
  /* Undo reversal */
  Trpath_reverse(tend_trpath,/*expect_fwd_p*/false);

  return trpath;
}


T
Trpath_solve_from_diagonals (int *found_score, bool *partialp, int sufficient_score, Trcoord_T middle_trdiagonal,
			     int middle_trdiagonal_qstart, int middle_trdiagonal_qend,
			     int middle_nmismatches,
			     Trdiag_T qstart_trdiag, Trdiag_T qend_trdiag,
			     bool tplusp, int querylength, Compress_T query_compress_tr,
			     int *mismatch_positions_alloc, 
			     Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			     Chrnum_T chrnum, int transcript_genestrand,
			     bool paired_end_p, bool first_read_p, Indelinfo_T indelinfo,
			     Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			     Method_T method) {

  T trpath = NULL, tstart_path, tend_path;
  Trcoord_T trdiagonal;

  int nmismatches_to_trimpos;
  int tstart, tend, trimpos, pos5, pos3;
  int adj;
  int found_score_ignore = querylength;

  bool want_lowest_coordinate_tr_p;

  if (middle_trdiagonal < (Trcoord_T) querylength) {
    return (T) NULL;
  } else if (transcript_genestrand > 0) {
    /* gplus */
    want_lowest_coordinate_tr_p = true;
  } else {
    /* gminus */
    want_lowest_coordinate_tr_p = false;
  }

  debug13(printf("tplusp %d.  transcript_genestrand %d => want_lowest_coordinate_tr_p %d\n",
		 tplusp,transcript_genestrand,want_lowest_coordinate_tr_p));
  debug13(printf("Entered Trpath_solve_from_diagonals with middle_diagonal %u, %d..%d, qstart_trdiag %p, qend_trdiag %p, want_lowest_coordinate_tr_p %d\n",
		 middle_trdiagonal,middle_trdiagonal_qstart,middle_trdiagonal_qend,
		 qstart_trdiag,qend_trdiag,want_lowest_coordinate_tr_p));

  if (qstart_trdiag != NULL) {
    tstart_path = Trpath_new_for_tstart_extension(middle_trdiagonal,middle_trdiagonal_qstart,middle_trdiagonal_qend,
						  middle_nmismatches,tplusp,trnum,troffset,trhigh,
						  chrnum,transcript_genestrand,
						  intlistpool,uintlistpool,trpathpool,method);

    attach_unknown_tstart(tstart_path,/*low_diagonal*/qstart_trdiag->trdiagonal,
			  /*low_tstart*/qstart_trdiag->qstart,
			  trhigh,querylength,indelinfo,query_compress_tr,
			  want_lowest_coordinate_tr_p,
			  intlistpool,uintlistpool,listpool,pathpool);

    debug13(printf("TSTART PATH FROM ATTACHED DIAGONAL:\n"));
    debug13(Trpath_print(tstart_path));

  } else if ((tstart =
	      Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_tr,
				     /*bits*/transcriptomebits,
				     /*univdiagonal*/(Univcoord_T) middle_trdiagonal,querylength,
				     /*pos5*/0,/*pos3*/middle_trdiagonal_qend,tplusp,/*genestrand*/0))
	     < middle_trdiagonal_qend) {
    tstart_path = Trpath_new_for_tstart_extension(middle_trdiagonal,tstart,/*tend*/middle_trdiagonal_qend,
						  nmismatches_to_trimpos,tplusp,trnum,troffset,trhigh,
						  chrnum,transcript_genestrand,
						  intlistpool,uintlistpool,trpathpool,method);
    debug13(printf("tstart %d, tend %d\n",tstart,middle_trdiagonal_qend));

    if (tstart <= 1) {
      /* Already extended to the start */
      debug13(printf("Not extending because tstart is %d\n",tstart));
      
    } else if ((adj = Genomebits_indel_solve_low(&trimpos,&nmismatches_to_trimpos,
						 /*univdiagonal*/(Univcoord_T) middle_trdiagonal,querylength,
						 /*pos5*/0,/*pos3*/tstart,
						 query_compress_tr,mismatch_positions_alloc,
						 transcriptomebits,/*bits_alt*/NULL,tplusp,/*genestrand*/0)) != 0) {
      debug13(printf("Genomebits_indel_solve_low succeeds with adj %d and %d mismatches_to_trimpos => trimpos %d\n",
		     adj,nmismatches_to_trimpos,trimpos));
      
      /* Subtract adj to get low diagonal, but add adj to get high diagonal */
      attach_unknown_tstart(tstart_path,/*low_diagonal*/middle_trdiagonal - adj,
			    /*low_tstart*/trimpos,
			    trhigh,querylength,indelinfo,query_compress_tr,
			    want_lowest_coordinate_tr_p,
			    intlistpool,uintlistpool,listpool,pathpool);
    }

    debug13(printf("TSTART PATH FROM COMPUTED INDEL:\n"));
    debug13(Trpath_print(tstart_path));

  } else {
    tstart_path = (T) NULL;
  }
  

  /* Final extension to qstart */
  if (tstart_path == (T) NULL) {
    /* Skip */
  } else if ((tstart = Intlist_head(tstart_path->endpoints)) <= 1) {
    /* Already extended to the start */
    debug13(printf("Not extending because tstart is %d\n",tstart));
  } else {
    trdiagonal = Uintlist_head(tstart_path->trdiagonals);
    pos3 = Intlist_second_value(tstart_path->endpoints);
    if ((tstart = Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_tr,
					 /*bits*/transcriptomebits,
					 /*univdiagonal*/(Univcoord_T) trdiagonal,querylength,
					 /*pos5*/0,pos3,tplusp,/*genestrand*/0)) < pos3) {
      Intlist_head_set(tstart_path->endpoints,tstart);
      Intlist_head_set(tstart_path->nmismatches,nmismatches_to_trimpos);
    }
    debug13(printf("TSTART PATH AFTER FINAL EXTENSION:\n"));
    debug13(Trpath_print(tstart_path));
  }


  if (qend_trdiag != NULL) {
    tend_path = Trpath_new_for_tend_extension(middle_trdiagonal,middle_trdiagonal_qstart,middle_trdiagonal_qend,
					      middle_nmismatches,tplusp,trnum,troffset,trhigh,
					      chrnum,transcript_genestrand,
					      intlistpool,uintlistpool,trpathpool,method);

    attach_unknown_tend(tend_path,/*high_diagonal*/qend_trdiag->trdiagonal,
			/*high_qend*/qend_trdiag->qend,
			trhigh,querylength,indelinfo,query_compress_tr,
			want_lowest_coordinate_tr_p,
			intlistpool,uintlistpool,listpool,pathpool);

    debug13(printf("TEND PATH FROM ATTACHED DIAGONAL:\n"));
    debug13(Trpath_print(tend_path));

  } else if ((tend =
	      Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_tr,
				   /*bits*/transcriptomebits,
				   /*univdiagonal*/(Univcoord_T) middle_trdiagonal,querylength,
				   /*pos5*/middle_trdiagonal_qstart,/*pos3*/querylength,tplusp,/*genestrand*/0))
	     > middle_trdiagonal_qstart) {
    tend_path = Trpath_new_for_tend_extension(middle_trdiagonal,/*tstart*/middle_trdiagonal_qstart,tend,
					      nmismatches_to_trimpos,tplusp,trnum,troffset,trhigh,
					      chrnum,transcript_genestrand,
					      intlistpool,uintlistpool,trpathpool,method);
    debug13(printf("tstart %d, tend %d\n",middle_trdiagonal_qstart,tend));

    if (tend >= querylength - 1) {
      /* Already extended to the end */
      debug13(printf("Not extending because tend is %d\n",tend));
      
    } else if ((adj = Genomebits_indel_solve_high(&trimpos,&nmismatches_to_trimpos,
						  /*univdiagonal*/(Univcoord_T) middle_trdiagonal,querylength,
						  /*pos5*/tend,/*pos3*/querylength,
						  query_compress_tr,mismatch_positions_alloc,
						  transcriptomebits,/*bits_alt*/NULL,tplusp,/*genestrand*/0)) != 0) {
      debug13(printf("Genomebits_indel_solve_high succeeds with adj %d and %d mismatches_to_trimpos => trimpos %d\n",
		     adj,nmismatches_to_trimpos,trimpos));
      
      /* Subtract adj to get low diagonal, but add adj to get high diagonal */
      attach_unknown_tend(tend_path,/*high_diagonal*/middle_trdiagonal + adj,
			  /*high_qend*/trimpos,
			  trhigh,querylength,indelinfo,query_compress_tr,
			  want_lowest_coordinate_tr_p,
			  intlistpool,uintlistpool,listpool,pathpool);
    }

    debug13(printf("TEND PATH FROM COMPUTED INDEL:\n"));
    debug13(Trpath_print(tend_path));

  } else {
    tend_path = (T) NULL;
  }


  /* Final extension to qend */
  if (tend_path == (T) NULL) {
    /* Skip */
  } else if ((tend = Intlist_head(tend_path->endpoints)) >= querylength - 1) {
    /* Already extended to the end */
    debug13(printf("Not extending because tend is %d\n",tend));
  } else {
    trdiagonal = Uintlist_head(tend_path->trdiagonals);
    pos5 = Intlist_second_value(tend_path->endpoints);
    if (tend_path->junctions != NULL) {
      pos5 += Junction_ninserts((Junction_T) List_head(tend_path->junctions));
    }
    if ((tend = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_tr,
				     /*bits*/transcriptomebits,
				     /*univdiagonal*/(Univcoord_T) trdiagonal,querylength,
				     pos5,/*pos3*/querylength,tplusp,/*genestrand*/0)) > pos5) {
      Intlist_head_set(tend_path->endpoints,tend);
      Intlist_head_set(tend_path->nmismatches,nmismatches_to_trimpos);
    }
    debug13(printf("TEND PATH AFTER FINAL EXTENSION:\n"));
    debug13(Trpath_print(tend_path));
  }


  if (tstart_path == NULL && tend_path == NULL) {
    return (T) NULL;
  } else if (tstart_path != NULL && tend_path == NULL) {
    trpath = tstart_path;
  } else if (tstart_path == NULL && tend_path != NULL) {
    trpath = tend_path;
  } else if ((trpath = combine_leftright_trpaths(tstart_path,tend_path,query_compress_tr,tplusp,
						 sufficient_score,querylength,
						 intlistpool,uintlistpool,listpool,
						 trpathpool,pathpool,method)) == NULL) {
    /* Must have exceeded nmismatches_allowed */
    debug13(printf("Result of combine_leftright_trpaths is NULL\n"));
    Trpath_free(&tstart_path,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
    Trpath_free(&tend_path,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
    return (T) NULL;
  } else {
    Trpath_free(&tstart_path,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
    Trpath_free(&tend_path,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
    /* Continue below */
  }

  /* Final filtering.  Requires attempts at final extension to qstart and qend */
  /* Sometimes a transcript is partially correct and the inner part can be finished by Pathpair_resolve */

  /* Testing shows slight improvement in using index1part_tr compared with 1 */
  if (Intlist_head(trpath->endpoints) > index1part_tr) {
    /* Doesn't extend to the end, and not due to a single mismatch */
    *partialp = true;
    Trpath_eval_nmatches(&found_score_ignore,trpath,querylength,query_compress_tr);

  } else if (Intlist_last_value(trpath->endpoints) < querylength - index1part_tr) {
    /* Doesn't extend to the end, and not due to a single mismatch */
    *partialp = true;
    Trpath_eval_nmatches(&found_score_ignore,trpath,querylength,query_compress_tr);

  } else {
    *partialp = false;
    Trpath_eval_nmatches(&(*found_score),trpath,querylength,query_compress_tr);
  }

  *partialp = false;

  debug13(printf("Returning result of combine_leftright_trpaths\n"));
  debug13(Trpath_print(trpath));
  /* Need to set found_score */
  return trpath;
}


T
Trpath_solve_from_trstart (Trcoord_T trdiagonal,
			   bool tplusp, int querylength, Compress_T query_compress_tr,
			   Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			   Trpathpool_T trpathpool, Method_T method) {
  /* T trpath; */
  /* int trimpos, pos5; */
  /* int adj; */

  int nmismatches_to_trimpos, found_score;
  int tend;

  Trnum_T trnum;
  Trcoord_T troffset, trhigh;

  Chrnum_T chrnum;
  int transcript_genestrand;

  debug13(printf("Entered Trpath_solve_from_trstart with trdiagonal %u, querylength %d\n",
		 trdiagonal,querylength));

  if ((tend = Genomebits_trim_qend(&nmismatches_to_trimpos,query_compress_tr,
				   /*bits*/transcriptomebits,
				   /*univdiagonal*/(Univcoord_T) trdiagonal,querylength,
				   /*pos5*/0,/*pos3*/querylength,tplusp,/*genestrand*/0)) < querylength - 1) {
    return (T) NULL;

  } else {
    found_score = (querylength - tend) + nmismatches_to_trimpos;
    trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonal - querylength,trdiagonal);
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    return Trpath_new_exact(trdiagonal,/*tstart*/0,tend,
			    nmismatches_to_trimpos,tplusp,trnum,troffset,trhigh,
			    chrnum,transcript_genestrand,
			    intlistpool,uintlistpool,trpathpool,found_score,method);
  }
}


T
Trpath_solve_from_trend (Trcoord_T trdiagonal,
			 bool tplusp, int querylength, Compress_T query_compress_tr,
			 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			 Trpathpool_T trpathpool, Method_T method) {
  /* T trpath; */
  /* int trimpos, pos3; */
  /* int adj; */
  int nmismatches_to_trimpos, found_score;
  int tstart;

  Trnum_T trnum;
  Trcoord_T troffset, trhigh;

  Chrnum_T chrnum;
  int transcript_genestrand;

  debug13(printf("Entered Trpath_solve_from_trend with trdiagonal %u, querylength %d\n",
		 trdiagonal,querylength));

  if ((tstart =
       Genomebits_trim_qstart(&nmismatches_to_trimpos,query_compress_tr,
			      /*bits*/transcriptomebits,
			      /*univdiagonal*/(Univcoord_T) trdiagonal,querylength,
			      /*pos5*/0,/*pos3*/querylength,tplusp,/*genestrand*/0)) > 1) {
    return (T) NULL;
  } else {
    found_score = tstart + nmismatches_to_trimpos;
    trnum = EF64_trnum(&troffset,&trhigh,transcript_ef64,trdiagonal - querylength,trdiagonal);
    chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
    return Trpath_new_exact(trdiagonal,tstart,/*tend*/querylength,
			    nmismatches_to_trimpos,tplusp,trnum,troffset,trhigh,
			    chrnum,transcript_genestrand,
			    intlistpool,uintlistpool,trpathpool,found_score,method);
  }
}


T
Trpath_solve_from_ends (int *found_score,
			Trcoord_T trdiagonal_i, int pos5_0, int pos3_0,
			Trcoord_T trdiagonal_j, int pos5_1, int pos3_1,
			bool tplusp, int querylength, Compress_T query_compress_tr,
			Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			Chrnum_T chrnum, int transcript_genestrand, Indelinfo_T indelinfo,
			Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			Method_T method) {
			
  T trpath = NULL;
  Junction_T junction;
  
  int nmismatches_i, nmismatches_j, ref_nmismatches_i, ref_nmismatches_j;
  int nindels;
  int supporti, supportj;

  int indel_pos, tstart, tend;
  /* Trcoord_T segmenti_left; */

  List_T j;
  Intlist_T q;
  Intlist_T s;

  bool want_lowest_coordinate_tr_p;

  if (transcript_genestrand > 0) {
    /* gplus */
    want_lowest_coordinate_tr_p = true;
  } else {
    /* gminus */
    want_lowest_coordinate_tr_p = false;
  }

  debug13(printf("Entered Trpath_solve_from_ends, with low_trdiagonal %u, %d..%d, and trdiagonal_j %u, %d..%d, want_lowest_coordinate_tr_p %d\n",
		 trdiagonal_i - troffset,pos5_0,pos3_0,trdiagonal_j - troffset,pos5_1,pos3_1,want_lowest_coordinate_tr_p));

  assert(trdiagonal_i != trdiagonal_j); /* Caller should handle this case */

  /* segmenti_left = trdiagonal_i - querylength; */
  /* segmentj_left = trdiagonal_j - querylength; */

  tstart = pos5_0;
  tend = pos3_1;

  /* Follows attach_unknown_tstart and attach_unknown_tend */
  if (trdiagonal_i > trdiagonal_j + max_insertionlen) {
    /* Impossible */
    debug13(printf("Impossible\n"));

  } else if (trdiagonal_i > trdiagonal_j) {
    /* (A) Insertion */
    nindels = trdiagonal_i - trdiagonal_j;
    debug13(printf("Insertion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_insertion(&nmismatches_i,&nmismatches_j,
						    &ref_nmismatches_i,&ref_nmismatches_j,
						    /*univdiagonal_i*/(Univcoord_T) trdiagonal_i,/*indels*/+nindels,
						    trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						    /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						    /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						    /*pos5*/tstart,/*pos3*/tend,querylength,
						    indelinfo,tplusp,/*genestrand*/0,
						    want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Insertion fails\n"));

    } else {
      supporti = indel_pos - tstart;
      supportj = tend - (indel_pos + nindels);
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
	
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	trpath = Trpath_new_from_ends(trdiagonal_i,pos5_0,pos3_0,trdiagonal_j,pos5_1,pos3_1,
				      tplusp,trnum,troffset,trhigh,chrnum,transcript_genestrand,
				      intlistpool,uintlistpool,listpool,trpathpool,
				      /*found_score*/nmismatches_i + nmismatches_j,/*total_ninserts*/nindels,
				      method);
	j = trpath->junctions; q = trpath->endpoints; s = trpath->nmismatches;

	junction = (Junction_T) List_head(j);
	Junction_free(&junction,pathpool);
	List_head_set(j,(void *) Junction_new_insertion(nindels,pathpool));

	/* No need to change trdiagonals */
	Intlist_head_set(q->rest,indel_pos);
	Intlist_head_set(s->rest,nmismatches_j);
	/* Intlist_head_set(r->rest,ref_nmismatches_j); */
	Intlist_head_set(s,nmismatches_i);
	/* Intlist_head_set(r,ref_nmismatches_i); */
      }
    }

  } else if (trdiagonal_i + max_deletionlen >= trdiagonal_j) {
    /* (B) Deletion (or short intron) */
    nindels = trdiagonal_j - trdiagonal_i;
    debug13(printf("Deletion with want_lowest_coordinate_tr_p %d\n",want_lowest_coordinate_tr_p));
    if ((indel_pos = Indel_resolve_middle_deletion(&nmismatches_i,&nmismatches_j,
						   &ref_nmismatches_i,&ref_nmismatches_j,
						   /*univdiagonal_i*/(Univcoord_T) trdiagonal_i,/*indels*/-nindels,
						   trhigh,/*mismatch_positions_left*/NULL,/*nmismatches_left*/0,
						   /*mismatch_positions_right*/NULL,/*nmismatches_right*/0,
						   /*ome*/transcriptomebits,/*ome_alt*/NULL,query_compress_tr,
						   /*pos5*/tstart,/*pos3*/tend,querylength,
						   indelinfo,tplusp,/*genestrand*/0,
						   want_lowest_coordinate_tr_p)) <= 0) {
      debug13(printf("Deletion fails\n"));
	  
    } else {
      supporti = indel_pos - tstart;
      supportj = tend - indel_pos;
      if (supporti - 3*nmismatches_i < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supporti,nmismatches_i));
      } else if (supportj - 3*nmismatches_j < MIN_SUPPORT_INDEL) {
	debug13(printf("Not enough support for indel: supporti %d and mismatches %d\n",supportj,nmismatches_j));
      } else {
	trpath = Trpath_new_from_ends(trdiagonal_i,/*tstart5*/pos5_0,/*tend5*/pos3_0,
				      trdiagonal_j,/*tstart3*/pos5_1,/*tend3*/pos3_1,
				      tplusp,trnum,troffset,trhigh,chrnum,transcript_genestrand,
				      intlistpool,uintlistpool,listpool,trpathpool,
				      /*found_score*/nmismatches_i + nmismatches_j,/*total_ninserts*/0,
				      method);
	j = trpath->junctions; q = trpath->endpoints; s = trpath->nmismatches;

	assert(nindels >= 0);
	/* deletionpos = segmenti_left + indel_pos; */

	junction = (Junction_T) List_head(j);
	Junction_free(&junction,pathpool);
	List_head_set(j,(void *) Junction_new_deletion(nindels,pathpool));

	/* No need to change trdiagonals */
	Intlist_head_set(q->rest,indel_pos);
	Intlist_head_set(s->rest,nmismatches_j);
	/* Intlist_head_set(r->rest,ref_nmismatches_j); */
	Intlist_head_set(s,nmismatches_i);
	/* Intlist_head_set(r,ref_nmismatches_i); */
      }
    }
  }

  if (trpath == NULL) {
    debug13(printf("Could not be solved\n"));
    return (T) NULL;

  } else if (Trpath_endpoints_acceptable_p(trpath->endpoints,trpath->junctions) == false) {
    debug13(printf("Endpoints were not acceptable\n"));
    Trpath_free(&trpath,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
    return (T) NULL;

  } else {
    /* Need to set found_score */
    Trpath_eval_nmatches(&(*found_score),trpath,querylength,query_compress_tr);
    return trpath;
  }
}


void
Trpath_solve_setup (Transcriptome_T transcriptome_in,
		    Genomebits_T transcriptomebits_in, EF64_T transcript_ef64_in,
		    int max_insertionlen_in, int max_deletionlen_in, int index1part_tr_in) {

  transcriptome = transcriptome_in;
  transcriptomebits = transcriptomebits_in;
  transcript_ef64 = transcript_ef64_in;
  max_insertionlen = max_insertionlen_in;
  max_deletionlen = max_deletionlen_in;
  index1part_tr = index1part_tr_in;

  return;
}
