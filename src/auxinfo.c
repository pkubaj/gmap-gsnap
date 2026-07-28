static char rcsid[] = "$Id: b387ae682f2a4f454698d47d4aa790e3598e558e $";
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
Auxinfo_pop (T old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool) {
  T prev;

  prev = old;
  old = prev->rest;

  if (prev->right_univdiags != NULL) {
    Univdiagpool_gc(&prev->right_univdiags,univdiagpool
		    univdiagpool_trace(__FILE__,__LINE__));
  }
  if (prev->left_univdiags != NULL) {
    Univdiagpool_gc(&prev->left_univdiags,univdiagpool
		    univdiagpool_trace(__FILE__,__LINE__));
  }
  
  Auxinfopool_free_auxinfo(&prev,auxinfopool
			   auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */

  return old;
}


/* Frees the entire list */
void
Auxinfo_free (T *old, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool) {
  T prev;

  debug0(printf("Freeing auxinfo %p\n",*old));

  while ((prev = *old) != NULL) {
    *old = prev->rest;

    if (prev->right_univdiags != NULL) {
      Univdiagpool_gc(&prev->right_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
    if (prev->left_univdiags != NULL) {
      Univdiagpool_gc(&prev->left_univdiags,univdiagpool
		      univdiagpool_trace(__FILE__,__LINE__));
    }
      
    Auxinfopool_free_auxinfo(&prev,auxinfopool
			     auxinfopool_trace(__FILE__,__LINE__)); /* Allocated by Pathpool_new_path */
  }

  return;
}


void
Auxinfo_gc (T *array, int n, Univdiagpool_T univdiagpool, Auxinfopool_T auxinfopool) {
  int i;
  T this;


  for (i = 0; i < n; i++) {
    this = array[i];
    Auxinfo_free(&this,univdiagpool,auxinfopool);
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
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

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
  new->nmismatches = -1;

  new->right_univdiags = (List_T) NULL;
  new->left_univdiags = (List_T) NULL;

  new->rest = (T) NULL;

  return new;
}

T
Auxinfo_new_univdiags (Method_T method, int qstart, int qend, int nmismatches,
		       List_T right_univdiags, List_T left_univdiags, Auxinfopool_T auxinfopool) {
  T new = Auxinfopool_new_auxinfo(auxinfopool
				  auxinfopool_trace(__FILE__,__LINE__));

  debug0(printf("Created auxinfo %p using Auxinfo_new_univdiags\n",new));

  new->method = method;
  new->qstart = qstart;
  new->qend = qend;
  new->nmismatches = nmismatches;

  new->right_univdiags = right_univdiags;
  new->left_univdiags = left_univdiags;

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

  assert(first != second);

  while (*p) {
    p = &(*p)->rest;
  }
  *p = second;

  return first;
}


