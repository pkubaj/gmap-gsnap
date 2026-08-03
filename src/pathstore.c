static char rcsid[] = "$Id: d24c359b18df7b27503e69d0599bd24967eceb0f $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "auxinfo.h"

#include <stdio.h>
#include "assert.h"
#include "mem.h"
#include "path-eval.h"		/* For Path_local_cmp */

/* Pathstore_new */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Pathstore_set_best */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif

/* Auxinfo_assign_chrinfo */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

static EF64_T chromosome_ef64;


#define T Pathstore_T


void
Pathstore_print (T this) {
  List_T p;

  printf("%p best:(%d,%d) complete:(%d,%d) ",
	 this,List_length(this->best_sense_paths),List_length(this->best_antisense_paths),
	 List_length(this->complete_sense_paths),List_length(this->complete_antisense_paths));

  for (p = this->complete_sense_paths; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  for (p = this->complete_antisense_paths; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }

  return;
}


void
Pathstore_free (T *old, Pathstorepool_T pathstorepool,
		Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, Pathpool_T pathpool,
		Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Freeing pathstore %p\n",++call_i,*old));

  Intlistpool_free_list(&(*old)->best_sense_partners,intlistpool
			intlistpool_trace(__FILE__,__LINE__));
  Intlistpool_free_list(&(*old)->best_antisense_partners,intlistpool
			intlistpool_trace(__FILE__,__LINE__));
  
  Hitlistpool_free_list(&(*old)->best_sense_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  Hitlistpool_free_list(&(*old)->best_antisense_paths,hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  
  Path_gc(&(*old)->unextended_sense_paths,intlistpool,univcoordlistpool,
	  listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&(*old)->unextended_antisense_paths,intlistpool,univcoordlistpool,
	  listpool,pathpool,transcriptpool,hitlistpool);
  
  Path_gc(&(*old)->complete_sense_paths,intlistpool,univcoordlistpool,
	  listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&(*old)->complete_antisense_paths,intlistpool,univcoordlistpool,
	  listpool,pathpool,transcriptpool,hitlistpool);
  
  Pathstorepool_free_pathstore(&(*old),pathstorepool
			       pathstorepool_trace(__FILE__,__LINE__));

  return;
}


void
Pathstore_gc (T *array, int n, Pathstorepool_T pathstorepool,
	      Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	      Listpool_T listpool, Pathpool_T pathpool,
	      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  int i;
  T this;

  for (i = 0; i < n; i++) {
    this = array[i];
    Pathstore_free(&this,pathstorepool,intlistpool,univcoordlistpool,
		   listpool,pathpool,transcriptpool,hitlistpool);
  }

  FREE(array);
  return;
}  



static List_T
best_paths (List_T paths, Hitlistpool_T hitlistpool) {
  List_T best;

  int npaths, i, j, k;
  Path_T *patharray, path0, path;
  int chop_qstart, chop_qend, qstart, qend;


  debug8(printf("Finding best paths\n"));
  if ((npaths = List_length(paths)) == 0) {
    return (List_T) NULL;

  } else if (npaths == 1) {
    path = (Path_T) List_head(paths);
    debug8(printf("Single: ")); debug8(Path_print(path)); debug8(printf("\n"));
    return Hitlist_push(NULL,hitlistpool,(void *) path
			hitlistpool_trace(__FILE__,__LINE__));

  } else {
    patharray = (Path_T *) List_to_array(paths,NULL);

    /* Eliminate duplicates */
    qsort(patharray,npaths,sizeof(Path_T),Path_structure_ignore_sense_cmp);

    k = 0;
    i = 0;
    while (i < npaths) {
      j = i + 1;
      while (j < npaths && Path_structure_ignore_sense_cmp(&(patharray[j]),&(patharray[i])) == 0) {
	debug8(printf("Identical: ")); debug8(Path_print(patharray[j]));
#if 0
	/* Don't want to free path, since it is included in a complete paths list */
	Path_free(&(patharray[j]),intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
#endif
	j++;
      }

      /* Found an identical group by structure.  Keep just the first one */
      patharray[k++] = patharray[i];
	
      i = j;
    }
    npaths = k;

    /* Find best ones */
    if (npaths > 1) {
      chop_qstart = Intlist_head(patharray[0]->endpoints);
      chop_qend = Intlist_last_value(patharray[0]->endpoints);

      for (i = 1; i < npaths; i++) {
	if ((qstart = Intlist_head(patharray[i]->endpoints)) > chop_qstart) {
	  chop_qstart = qstart;
	}
	if ((qend = Intlist_last_value(patharray[i]->endpoints)) < chop_qend) {
	  chop_qend = qend;
	}
      }

      for (i = 0; i < npaths; i++) {
	patharray[i]->chop_qstart = chop_qstart;
	patharray[i]->chop_qend = chop_qend;
      }

      qsort(patharray,npaths,sizeof(Path_T),Path_best_cmp);
    }

    path0 = patharray[0];
    debug8(printf("Best: ")); debug8(Path_print(path0));
    best = Hitlist_push(NULL,hitlistpool,(void *) path0
			hitlistpool_trace(__FILE__,__LINE__));
    for (i = 1; i < npaths; i++) {
      path = patharray[i];
      if (Path_best_cmp(&path,&path0) <= 0) {
	debug8(printf("Tie: ")); debug8(Path_print(path));
	best = Hitlist_push(best,hitlistpool,(void *) path
			    hitlistpool_trace(__FILE__,__LINE__));
      } else {
	debug8(printf("Not best: ")); debug8(Path_print(path));
#if 0
	/* Don't want to free path, since it is included in a complete paths list */
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
#endif
      }
    }
    debug8(printf("\n"));

    FREE(patharray);
    return List_reverse(best);
  }
}


void
Pathstore_set_best_paths (T this, Hitlistpool_T hitlistpool) {

  debug8(printf("Setting best paths for pathstore %p\n",this));
  if (this->complete_sense_paths != NULL) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else {
    this->best_sense_paths = best_paths(this->unextended_sense_paths,hitlistpool);
  }

  if (this->complete_antisense_paths != NULL) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else {
    this->best_antisense_paths = best_paths(this->unextended_antisense_paths,hitlistpool);
  }

  return;
}


void
Pathstore_set_best_sense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p) {

  debug8(printf("Setting best sense paths for pathstore %p\n",this));
  if (0 && only_complete_p == true) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else if (this->complete_sense_paths != NULL) {
    this->complete_sense_p = true;
    this->best_sense_paths = best_paths(this->complete_sense_paths,hitlistpool);
  } else {
    this->best_sense_paths = best_paths(this->unextended_sense_paths,hitlistpool);
  }

  return;
}

void
Pathstore_set_best_antisense_paths (T this, Hitlistpool_T hitlistpool, bool only_complete_p) {

  debug8(printf("Setting best antisense paths for pathstore %p\n",this));
  if (0 && only_complete_p == true) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else if (this->complete_antisense_paths != NULL) {
    this->complete_antisense_p = true;
    this->best_antisense_paths = best_paths(this->complete_antisense_paths,hitlistpool);
  } else {
    this->best_antisense_paths = best_paths(this->unextended_antisense_paths,hitlistpool);
  }

  return;
}


T
Pathstore_merge (T old, T *new, Pathstorepool_T pathstorepool, Intlistpool_T intlistpool) {

  debug0(printf("Merging pathstore %p into %p\n",*new,old));

  if ((*new)->solvedp == true) {
    old->solvedp = true;
  }
  if ((*new)->complete_sense_p == true || (*new)->complete_sense_paths != NULL) {
    old->complete_sense_p = true;
  }
  if ((*new)->complete_antisense_p == true || (*new)->complete_antisense_paths != NULL) {
    old->complete_antisense_p = true;
  }

  old->best_sense_paths = List_append(old->best_sense_paths,(*new)->best_sense_paths);
  old->best_antisense_paths = List_append(old->best_antisense_paths,(*new)->best_antisense_paths);

  old->unextended_sense_paths = List_append(old->unextended_sense_paths,(*new)->unextended_sense_paths);
  old->unextended_antisense_paths = List_append(old->unextended_antisense_paths,(*new)->unextended_antisense_paths);

  old->complete_sense_paths = List_append(old->complete_sense_paths,(*new)->complete_sense_paths);
  old->complete_antisense_paths = List_append(old->complete_antisense_paths,(*new)->complete_antisense_paths);
  
  /* Do not free paths */
  Intlistpool_free_list(&(*new)->best_sense_partners,intlistpool
			intlistpool_trace(__FILE__,__LINE__));
  Intlistpool_free_list(&(*new)->best_antisense_partners,intlistpool
			intlistpool_trace(__FILE__,__LINE__));

  Pathstorepool_free_pathstore(&(*new),pathstorepool
			       pathstorepool_trace(__FILE__,__LINE__));

  return old;
}


/* Standard version */
T
Pathstore_new (Pathstorepool_T pathstorepool) {
  T new = Pathstorepool_new_pathstore(pathstorepool
				      pathstorepool_trace(__FILE__,__LINE__));

#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Created pathstore %p using Pathstore_new\n",++call_i,new));

  new->chrnum = 0;

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = (List_T) NULL;

  return new;
}


T
Pathstore_new_tr (Pathstorepool_T pathstorepool, Chrnum_T chrnum, Univcoord_T chroffset, Univcoord_T chrhigh) {
  T new = Pathstorepool_new_pathstore(pathstorepool
				      pathstorepool_trace(__FILE__,__LINE__));

#ifdef DEBUG0
  static int call_i = 0;
#endif
  
  debug0(printf("%d: Created pathstore %p using Pathstore_new_tr\n",++call_i,new));

  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->chrhigh = chrhigh;

  new->solvedp = true;
  new->complete_sense_p = false;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = (List_T) NULL;

  return new;
};



T
Pathstore_new_sense_path (Path_T path, Hitlistpool_T hitlistpool, Pathstorepool_T pathstorepool) {
  T new = Pathstorepool_new_pathstore(pathstorepool
				      pathstorepool_trace(__FILE__,__LINE__));

#ifdef DEBUG0
  static int call_i = 0;
#endif

  debug0(printf("%d: Created pathstore %p using Pathstore_new_sense_path\n",++call_i,new));

  new->solvedp = false;
  new->complete_sense_p = true;
  new->complete_antisense_p = false;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
				       hitlistpool_trace(__FILE__,__LINE__));
  new->best_antisense_paths = (List_T) NULL;

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
					   hitlistpool_trace(__FILE__,__LINE__));
  new->complete_antisense_paths = (List_T) NULL;

  return new;
}


T
Pathstore_new_antisense_path (Path_T path, Hitlistpool_T hitlistpool, Pathstorepool_T pathstorepool) {
  T new = Pathstorepool_new_pathstore(pathstorepool
				      pathstorepool_trace(__FILE__,__LINE__));

#ifdef DEBUG0
  static int call_i = 0;
#endif

  debug0(printf("%d: Created pathstore %p using Pathstore_new_antisense_path\n",++call_i,new));

  new->solvedp = false;
  new->complete_sense_p = false;
  new->complete_antisense_p = true;

  new->best_sense_partners = (Intlist_T) NULL;
  new->best_antisense_partners = (Intlist_T) NULL;

  new->best_sense_paths = (List_T) NULL;
  new->best_antisense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
					   hitlistpool_trace(__FILE__,__LINE__));

  new->unextended_sense_paths = (List_T) NULL;
  new->unextended_antisense_paths = (List_T) NULL;
  new->complete_sense_paths = (List_T) NULL;
  new->complete_antisense_paths = Hitlist_push(NULL,hitlistpool,(void *) path
					       hitlistpool_trace(__FILE__,__LINE__));

  return new;
}


#if 0
void
Pathstore_collect_unextended_paths (List_T *unextended_sense_paths,
				  List_T *unextended_antisense_paths,
				  T *pathstore_array, int nunivdiagonals,
				  Hitlistpool_T hitlistpool) {
  T pathstore;
  Path_T path;
  List_T p;
  int i;

  *unextended_sense_paths = *unextended_antisense_paths = (List_T) NULL;

  for (i = 0; i < nunivdiagonals; i++) {
    pathstore = pathstore_array[i];
    for (p = pathstore->unextended_sense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      *unextended_sense_paths =
	Hitlist_push(*unextended_sense_paths,hitlistpool,(void *) path
		     hitlistpool_trace(__FILE__,__LINE__));
    }

    for (p = pathstore->unextended_antisense_paths; p != NULL; p = List_next(p)) {
      path = (Path_T) List_head(p);
      *unextended_antisense_paths =
	Hitlist_push(*unextended_antisense_paths,hitlistpool,(void *) path
		     hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  *unextended_sense_paths = List_reverse(*unextended_sense_paths);
  *unextended_antisense_paths = List_reverse(*unextended_antisense_paths);

  return;
}
#endif


void
Pathstore_collect_paths (bool *foundp, List_T *sense_paths, List_T *antisense_paths,
		       T *pathstore_array, int nunivdiagonals,
		       Hitlistpool_T hitlistpool) {

  T this;
  Path_T path;
  List_T p;
  int i;

  *sense_paths = *antisense_paths = (List_T) NULL;

  for (i = 0; i < nunivdiagonals; i++) {
    this = pathstore_array[i];

    /* Collect the best ones */
    for (p = this->best_sense_paths; p != NULL; p = List_next(p)) {
      *foundp = true;
      path = (Path_T) List_head(p);
      *sense_paths = Hitlist_push(*sense_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    }
    for (p = this->best_antisense_paths; p != NULL; p = List_next(p)) {
      *foundp = true;
      path = (Path_T) List_head(p);
      *antisense_paths = Hitlist_push(*antisense_paths,hitlistpool,(void *) path
				  hitlistpool_trace(__FILE__,__LINE__));
    }
  }

  *sense_paths = List_reverse(*sense_paths);
  *antisense_paths = List_reverse(*antisense_paths);

  return;
}


void
Pathstore_assign_chrinfo (Univcoord_T * univdiagonals, T *pathstore, int n, int querylength) {
  int i, j;
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  T pathstore_i, pathstore_j;

  debug9(printf("Entering Pathstore_assign_chrinfo with %d univdiagonals/pathstore\n",n));
  i = 0;
  while (i < n && univdiagonals[i] < (Univcoord_T) querylength) {
    pathstore_i = pathstore[i];
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 /*univdiagonals[i] - querylength*/0,univdiagonals[i]);
    pathstore_i->chrnum = chrnum;
    pathstore_i->chroffset = chroffset;
    pathstore_i->chrhigh = chrhigh;

    debug9(printf("i* %d %u %d %u..%u\n",i,univdiagonals[i],chrnum,chroffset,chrhigh));
    
    j = i + 1;
    while (j < n && univdiagonals[j] < chrhigh) {
      pathstore_j = pathstore[j];
      pathstore_j->chrnum = chrnum;
      pathstore_j->chroffset = chroffset;
      pathstore_j->chrhigh = chrhigh;
      debug9(printf("j* %d %u %d %u..%u\n",j,univdiagonals[j],chrnum,chroffset,chrhigh));

      j++;
    }
    
    i = j;
  }

  while (i < n) {
    pathstore_i = pathstore[i];
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonals[i] - querylength,univdiagonals[i]);
    pathstore_i->chrnum = chrnum;
    pathstore_i->chroffset = chroffset;
    pathstore_i->chrhigh = chrhigh;

    debug9(printf("i %d %u %d %u..%u\n",i,univdiagonals[i],chrnum,chroffset,chrhigh));
    
    j = i + 1;
    while (j < n && univdiagonals[j] < chrhigh) {
      pathstore_j = pathstore[j];
      pathstore_j->chrnum = chrnum;
      pathstore_j->chroffset = chroffset;
      pathstore_j->chrhigh = chrhigh;
      debug9(printf("j %d %u %d %u..%u\n",j,univdiagonals[j],chrnum,chroffset,chrhigh));

      j++;
    }
    
    i = j;
  }

  return;
}


void
Pathstore_setup (EF64_T chromosome_ef64_in) {

  chromosome_ef64 = chromosome_ef64_in;

  return;
}
