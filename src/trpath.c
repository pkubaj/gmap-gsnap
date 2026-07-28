static char rcsid[] = "$Id: bca0ed5ea758bf30f802e36d3bb04dc048136821 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "trpath.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mem.h"
#include "assert.h"


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Allocation and freeing of objects */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


/* Printing */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Trpath_T

T
Trpath_create (Intlist_T endpoints, Uintlist_T trdiagonals, Intlist_T nmismatches,
	       List_T junctions, bool tplusp,
	       Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
	       Chrnum_T chrnum, int transcript_genestrand,
	       Trpathpool_T trpathpool, int found_score, int total_ninserts, Method_T method) {

  T new = Trpathpool_new_trpath(trpathpool
				trpathpool_trace(__FILE__,__LINE__));
  Intlist_T q;
#ifdef CHECK_ASSERTIONS
  int prev_endpoint;
#endif

  debug0(printf("Creating trpath %p by Trpath_create with trnum %d, found score %d, mismatches %s\n",
		new,trnum,found_score,Intlist_to_string(nmismatches)));
  
  new->tplusp = tplusp;

  new->trnum = trnum;
  new->troffset = troffset;
  new->trhigh = trhigh;

  new->chrnum = chrnum;
  new->transcript_genestrand = transcript_genestrand;

  new->endpoints = endpoints;
  new->trdiagonals = trdiagonals;
  new->nmismatches = nmismatches;
  new->junctions = junctions;

  new->found_score = found_score;
  new->total_ninserts = total_ninserts;

  new->method = method;

#ifdef CHECK_ASSERTIONS
  prev_endpoint = Intlist_head(new->endpoints);
  for (q = Intlist_next(new->endpoints); q != NULL; q = Intlist_next(q)) {
    if (Intlist_head(q) <= prev_endpoint) {
      printf("Trpath_create expected forward, but got\n");
      Trpath_print(new);
      abort();
    }
    prev_endpoint = Intlist_head(q);
  }
#endif

  return new;
}


void
Trpath_free (T *old, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	     Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool) {
  Intlistpool_free_list(&(*old)->endpoints,intlistpool
			intlistpool_trace(__FILE__,__LINE__));
  Uintlistpool_free_list(&(*old)->trdiagonals,uintlistpool
			  uintlistpool_trace(__FILE__,__LINE__));
  Intlistpool_free_list(&(*old)->nmismatches,intlistpool
			intlistpool_trace(__FILE__,__LINE__));
  Junction_list_gc(&(*old)->junctions,listpool,pathpool);
  Trpathpool_free_trpath(&(*old),trpathpool
			 trpathpool_trace(__FILE__,__LINE__));
  return;
}


void
Trpath_gc (List_T *list, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
	   Listpool_T listpool, Pathpool_T pathpool, Trpathpool_T trpathpool,
	   Hitlistpool_T hitlistpool) {
  List_T p;
  T old;
  
  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Trpath_free(&old,intlistpool,uintlistpool,listpool,trpathpool,pathpool);
  }
  Hitlistpool_free_list(&(*list),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__)); /* allocated by Hitlistpool_push */
  return;
}


int
Trpath_trnum_cmp (const void *x, const void *y) {
  T a = * (T *) x;
  T b = * (T *) y;

  if (a->trnum < b->trnum) {
    return -1;
  } else if (b->trnum < a->trnum) {
    return +1;
  } else {
    return 0;
  }
}


/* fwd direction */
T
Trpath_new_for_tstart_extension (Trcoord_T trdiagonal, int tstart, int tend, int nmismatches,
				 bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
				 Chrnum_T chrnum, int transcript_genestrand,
				 Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
				 Trpathpool_T trpathpool, Method_T method) {
  T new = Trpathpool_new_trpath(trpathpool
				trpathpool_trace(__FILE__,__LINE__));

  new->tplusp = tplusp;

  new->trnum = trnum;
  new->troffset = troffset;
  new->trhigh = trhigh;

  new->chrnum = chrnum;
  new->transcript_genestrand = transcript_genestrand;

  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,tend
						     intlistpool_trace(__FILE__,__LINE__)),intlistpool,tstart
				    intlistpool_trace(__FILE__,__LINE__));
  new->trdiagonals = Uintlistpool_push(NULL,uintlistpool,trdiagonal
				        uintlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__)); /* nmismatches_to_trimpos */
  new->junctions = (List_T) NULL;
  
  new->found_score = -1;
  new->total_ninserts = -1;

  new->method = method;

  return new;
}


/* rev direction */
T
Trpath_new_for_tend_extension (Trcoord_T trdiagonal, int tstart, int tend, int nmismatches,
			       bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
			       Chrnum_T chrnum, int transcript_genestrand,
			       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			       Trpathpool_T trpathpool, Method_T method) {
  T new = Trpathpool_new_trpath(trpathpool
				trpathpool_trace(__FILE__,__LINE__));

  new->tplusp = tplusp;

  new->trnum = trnum;
  new->troffset = troffset;
  new->trhigh = trhigh;

  new->chrnum = chrnum;
  new->transcript_genestrand = transcript_genestrand;

  new->endpoints = Intlistpool_push(Intlistpool_push(NULL,intlistpool,tstart
						     intlistpool_trace(__FILE__,__LINE__)),intlistpool,tend
				    intlistpool_trace(__FILE__,__LINE__));
  new->trdiagonals = Uintlistpool_push(NULL,uintlistpool,trdiagonal
				        uintlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__));
  new->junctions = (List_T) NULL;

  new->found_score = -1;
  new->total_ninserts = -1;

  new->method = method;

  return new;
}


T
Trpath_new_from_ends (Trcoord_T trdiagonal5, int tstart5, int tend5,
		      Trcoord_T trdiagonal3, int tstart3, int tend3,
		      bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
		      Chrnum_T chrnum, int transcript_genestrand,
		      Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Listpool_T listpool, Trpathpool_T trpathpool,
		      int found_score, int total_ninserts, Method_T method) {
  (void)(tend5);

  T new = Trpathpool_new_trpath(trpathpool
				trpathpool_trace(__FILE__,__LINE__));
  debug0(printf("Creating trpath %p by Trpath_new_from_ends\n",new));
  
  new->tplusp = tplusp;
  new->trnum = trnum;
  new->troffset = troffset;
  new->trhigh = trhigh;

  new->chrnum = chrnum;
  new->transcript_genestrand = transcript_genestrand;

  new->endpoints = Intlistpool_push(NULL,intlistpool,tend3
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,tstart3
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,tstart5
				    intlistpool_trace(__FILE__,__LINE__));

  new->trdiagonals = Uintlistpool_push(NULL,uintlistpool,trdiagonal3
				        uintlistpool_trace(__FILE__,__LINE__));
  new->trdiagonals = Uintlistpool_push(new->trdiagonals,uintlistpool,trdiagonal5
				        uintlistpool_trace(__FILE__,__LINE__));

  new->nmismatches = Intlistpool_push(NULL,intlistpool,-1
				      intlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(new->nmismatches,intlistpool,-1
				      intlistpool_trace(__FILE__,__LINE__));

#ifdef ALLOCATE_UNSOLVED_JUNCTION
  new->junctions = Listpool_push(NULL,listpool,Junction_new_unsolved(pathpool)
				 listpool_trace(__FILE__,__LINE__));
#else
  new->junctions = Listpool_push(NULL,listpool,(void *) JUNCTION_UNSOLVED
				 listpool_trace(__FILE__,__LINE__));
#endif
  
  new->found_score = found_score;
  new->total_ninserts = total_ninserts;

  new->method = method;

  return new;
}


T
Trpath_new_exact (Trcoord_T trdiagonal, int tstart, int tend, int nmismatches,
		  bool tplusp, Trnum_T trnum, Trcoord_T troffset, Trcoord_T trhigh,
		  Chrnum_T chrnum, int transcript_genestrand,
		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		  Trpathpool_T trpathpool, int found_score, Method_T method) {

  T new = Trpathpool_new_trpath(trpathpool
				trpathpool_trace(__FILE__,__LINE__));
  debug0(printf("Creating trpath %p by Trpath_new_exact with trnum %d, found score %d, nmismatches %d\n",
		new,trnum,found_score,nmismatches));
  
  new->tplusp = tplusp;
  new->trnum = trnum;
  new->troffset = troffset;
  new->trhigh = trhigh;

  new->chrnum = chrnum;
  new->transcript_genestrand = transcript_genestrand;

  new->endpoints = Intlistpool_push(NULL,intlistpool,tend
				    intlistpool_trace(__FILE__,__LINE__));
  new->endpoints = Intlistpool_push(new->endpoints,intlistpool,tstart
				    intlistpool_trace(__FILE__,__LINE__));

  new->trdiagonals = Uintlistpool_push(NULL,uintlistpool,trdiagonal
				        uintlistpool_trace(__FILE__,__LINE__));
  new->nmismatches = Intlistpool_push(NULL,intlistpool,nmismatches
				      intlistpool_trace(__FILE__,__LINE__));
  new->junctions = (List_T) NULL;
  
  new->found_score = found_score;
  new->total_ninserts = 0;

  new->method = method;

  return new;
}


T
Trpath_reverse (T this, bool expect_fwd_p) {
  Intlist_T q;
#ifdef CHECK_ASSERTIONS
  int prev_endpoint;
#endif

  this->endpoints = Intlist_reverse(this->endpoints);
  this->trdiagonals = Uintlist_reverse(this->trdiagonals);
  this->nmismatches = Intlist_reverse(this->nmismatches);
  this->junctions = List_reverse(this->junctions);

#ifdef CHECK_ASSERTIONS
  if (expect_fwd_p == true) {
    prev_endpoint = Intlist_head(this->endpoints);
    for (q = Intlist_next(this->endpoints); q != NULL; q = Intlist_next(q)) {
      if (Intlist_head(q) <= prev_endpoint) {
	printf("Expecting forward, but got\n");
	Trpath_print(this);
	abort();
      }
      prev_endpoint = Intlist_head(q);
    }
 
 } else {
    prev_endpoint = Intlist_head(this->endpoints);
    for (q = Intlist_next(this->endpoints); q != NULL; q = Intlist_next(q)) {
      if (Intlist_head(q) >= prev_endpoint) {
	printf("Expecting reverse, but got\n");
	Trpath_print(this);
	abort();
      }
      prev_endpoint = Intlist_head(q);
    }
  }
#endif

  return this;
}


List_T
Trpath_filter_best (List_T *subopt, List_T trpaths, int found_score) {
  List_T best = NULL, p, pnext;
  T trpath;

  p = trpaths;
  while (p != NULL) {
    trpath = (T) List_head(trpaths);
    if (trpath->found_score == found_score) {
      /* Transfer onto best */
      pnext = p->rest;
      p->rest = best;
	best = p;
    } else {
      /* Transfer onto subopt */
      pnext = p->rest;
      p->rest = *subopt;
      *subopt = p;
    }
    p = pnext;
  }

  return best;
}



#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
void
Trpath_print (T this) {

  if (this != NULL) {
    printf(">> Trpath %p: ",this);

    /* Info added for new version */
    printf("tplusp:%d ",this->tplusp);

    printf("%s  %d:%s  %s  nmismatches:%s  ",
	   Uintlist_to_string(this->trdiagonals),this->trnum,
	   Uintlist_to_string_offset(this->trdiagonals,this->troffset),
	   Intlist_to_string(this->endpoints),Intlist_to_string(this->nmismatches));

    printf("jcns:");
    Junction_print_list(this->junctions);

    printf("\n");
  }

  return;
}
#endif


