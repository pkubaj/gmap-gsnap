#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "auxinfo.h"

#include <stdio.h>
#include "assert.h"
#include "mem.h"
#include "path-eval.h"		/* For Path_local_cmp */

/* Auxinfo_new */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


#define T Auxinfo_T


void
Auxinfo_print (T this) {
  while (this != NULL) {
    printf("%p.%s ",this,Method_string(this->method));
    this = this->rest;
  }

  return;
}


T
Auxinfo_pop (T old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Listpool_T listpool) {
  T prev;
  List_T chain, p;
  int nsplices;

  prev = old;
  old = prev->rest;

  for (nsplices = 0; nsplices <= prev->max_nsplices; nsplices++) {
    for (p = prev->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
      chain = (List_T) List_head(p);
      Univdiagpool_gc(&chain,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
    Listpool_free_list(&(prev->chains_by_nsplices[nsplices]),listpool
		       listpool_trace(__FILE__,__LINE__));
  }
  FREE(prev->chains_by_nsplices);

  Auxinfopool_free_auxinfo(&prev,auxinfopool
			   auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */

  return old;
}


/* Frees the entire list */
void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Listpool_T listpool,
	      bool free_chains_p) {
  T prev;
  List_T chain, p;
  int nsplices;

  debug0(printf("Freeing auxinfo %p\n",*old));

  while ((prev = *old) != NULL) {
    *old = prev->rest;

    if (free_chains_p == true) {
      for (nsplices = 0; nsplices <= prev->max_nsplices; nsplices++) {
	for (p = prev->chains_by_nsplices[nsplices]; p != NULL; p = List_next(p)) {
	  chain = (List_T) List_head(p);
	  Univdiagpool_gc(&chain,univdiagpool
			  univdiagpool_trace(__FILE__,__LINE__));
	}
	Listpool_free_list(&(prev->chains_by_nsplices[nsplices]),listpool
			   listpool_trace(__FILE__,__LINE__));
      }
    }
    FREE(prev->chains_by_nsplices);
    
    Auxinfopool_free_auxinfo(&prev,auxinfopool
			     auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */
  }

  return;
}


void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool, Listpool_T listpool) {
  int i;
  T this;


  for (i = 0; i < n; i++) {
    this = array[i];
    Auxinfo_free(&this,univdiagpool,auxinfopool,listpool,/*free_chains_p*/true);
  }

  FREE(array);
  return;
}  



/* Standard version */
T
Auxinfo_new (Method_T method, int qstart, int qend, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  debug0(printf("Created auxinfo %p using Auxinfo_new\n",new));

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;

  new->chains_by_nsplices = (List_T *) NULL;
  new->nsplices = 0;
  new->max_nsplices = -1;
  new->chaindir = CHAIN_NONE;

  new->elti = -1;

  new->rest = (T) NULL;

  return new;
}


T
Auxinfo_new_tr (Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  debug0(printf("Created auxinfo %p using Auxinfo_new_tr\n",new));

  new->method = TR_METHOD;	/* Some tr method, e.g., TR_EXACT1 */
  new->qstart = -1;
  new->qend = -1;

  new->chains_by_nsplices = (List_T *) NULL;
  new->nsplices = 0;
  new->max_nsplices = -1;
  new->chaindir = CHAIN_NONE;

  new->elti = -1;

  new->rest = (T) NULL;

  return new;
}




T
Auxinfo_new_chains (Method_T method, int qstart, int qend,
		    List_T chains, int chaindir, int nsplices, int max_nsplices,
		    Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  debug0(printf("Created auxinfo %p using Auxinfo_new_chains\n",new));

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;

  new->chains_by_nsplices = (List_T *) CALLOC(max_nsplices + 1,sizeof(List_T));
  new->chains_by_nsplices[nsplices] = chains;
  new->chaindir = chaindir;
  new->nsplices = nsplices;
  new->max_nsplices = max_nsplices;

  new->elti = -1;

  new->rest = (T) NULL;

  return new;
}


T
Auxinfo_new_elt (Method_T method, int qstart, int qend, int elti,
		 Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  debug0(printf("Created auxinfo %p using Auxinfo_new\n",new));

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;

  new->chains_by_nsplices = (List_T *) NULL;
  new->nsplices = 0;
  new->max_nsplices = -1;
  new->chaindir = CHAIN_NONE;

  new->elti = elti;

  new->rest = (T) NULL;

  return new;
}


int
Auxinfo_length (T list) {
  int n;
  
  for (n = 0; list; list = list->rest) {
    n++;
  }

  return n;
}


/* Modified from List_append */
T
Auxinfo_append (T first, T second) {
  T *p = &first;

  /* Doesn't hold when first and second are NULL */
  /* assert(first != second); */

  while (*p) {
    p = &(*p)->rest;
  }
  *p = second;

  return first;
}


