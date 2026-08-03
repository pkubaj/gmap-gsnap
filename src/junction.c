#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "junction.h"
#include "mem.h"
#include "assert.h"
#include "maxent_hr.h"
#include "sense.h"


#define MIN_INTRONLEN 30	/* Same as in splice.c */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Debugging procedures */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Junction_T


#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
void
Junction_print (T this) {
  if (this == NULL) {
#ifdef ALLOCATE_UNSOLVED_JUNCTION
    printf("No junction\n");
#else
    printf("Unsolved junction\n");
#endif
  } else if (this->type == INS_JUNCTION) {
    printf("Insertion of %d\n",this->nindels);
  } else if (this->type == DEL_JUNCTION) {
    printf("Deletion of %d\n",this->nindels);
  } else if (this->type == SPLICE_JUNCTION) {
    if (this->splice_distance == 0) {
      printf("Splice ambiguous with sense %d, prob %f and %f\n",
	     this->sensedir,this->donor_prob,this->acceptor_prob);
    } else {
      printf("Splice with sense %d of %u, prob %f and %f\n",
	     this->sensedir,this->splice_distance,this->donor_prob,this->acceptor_prob);
    }
  } else if (this->type == UNSOLVED_JUNCTION) {
    printf("Unsolved junction\n");
  }
  return;
}
#endif

#if defined(CHECK_ASSERTIONS) || defined(DEBUG1)
void
Junction_print_list (List_T list) {
  T this;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this == NULL) {
#ifdef ALLOCATE_UNSOLVED_JUNCTION
      printf("None,");
#else
      printf("Unsolved,");
#endif
    } else if (this->type == INS_JUNCTION) {
      printf("Ins:%d,",this->nindels);
    } else if (this->type == DEL_JUNCTION) {
      printf("Del:%d,",this->nindels);
    } else if (this->type == SPLICE_JUNCTION) {
      if (0 && this->splice_distance == 0) {
	/* Should not happen with current representation */
	printf("Amb:%f-%f,",this->donor_prob,this->acceptor_prob);
      } else if (this->sensedir == SENSE_FORWARD) {
	printf("Splice(sense,%f-%f):%u,",
	       this->donor_prob,this->acceptor_prob,this->splice_distance);
      } else if (this->sensedir == SENSE_ANTI) {
	printf("Splice(antisense,%f-%f):%u,",
	       this->donor_prob,this->acceptor_prob,this->splice_distance);
      } else if (this->sensedir == SENSE_NULL) {
	printf("Splice(null,%f-%f):%u,",
	       this->donor_prob,this->acceptor_prob,this->splice_distance);
      }
    } else if (this->type == UNSOLVED_JUNCTION) {
      printf("Unsolved,");
    }
  }

  return;
}
#endif

void
Junction_free (T *old, Pathpool_T pathpool) {
  if (*old) {
    Pathpool_free_junction(&(*old),pathpool
			   pathpool_trace(__FILE__,__LINE__));
  }
  return;
}

void
Junction_list_gc (List_T *list, Listpool_T listpool, Pathpool_T pathpool) {
  List_T p;
  T old;

  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Junction_free(&old,pathpool);
  }
  Listpool_free_list(&(*list),listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */
  return;
}

T
Junction_new_insertion (int nindels, Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  assert(nindels > 0);

  new->type = INS_JUNCTION;
  new->nindels = nindels;
  new->annotatedp = false;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_deletion (int nindels, Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  assert(nindels > 0);

  new->type = DEL_JUNCTION;
  new->nindels = nindels;
  new->annotatedp = false;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}

T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob,
		     Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  assert((int) splice_distance > 0);

  if (splice_distance < MIN_INTRONLEN && donor_prob <= 1.0 && acceptor_prob <= 1.00) {
    /* Short splices are reported as deletions, unless they are known introns */
    return Junction_new_deletion(/*nindels*/(int) splice_distance,pathpool);

  } else {
    new->type = SPLICE_JUNCTION;
    new->nindels = 0;
    new->annotatedp = false;

    new->splice_distance = splice_distance;
    new->sensedir = sensedir;
    new->donor_prob = donor_prob;
    new->acceptor_prob = acceptor_prob;

    /* printf("DONOR_PROB %f, ACCEPTOR_PROB %f\n",donor_prob,acceptor_prob); */

    return new;
  }
}


T
Junction_new_splice_annotated (Chrpos_T splice_distance, int sensedir,
			       double donor_prob, double acceptor_prob,
			       Pathpool_T pathpool) {
  T new = Junction_new_splice(splice_distance,sensedir,donor_prob,acceptor_prob,pathpool);

  if (Junction_type(new) == SPLICE_JUNCTION) {
    new->annotatedp = true;
  }

  return new;
}


T
Junction_new_ambig_splice (int sensedir, double donor_prob, double acceptor_prob,
			   Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  new->type = SPLICE_JUNCTION;
  new->nindels = 0;
  new->annotatedp = false;

  new->splice_distance = 0;    /* A zero splice distance is created for ambiguous splices */
  new->sensedir = sensedir;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}


T
Junction_new_chimera (char donor1, char donor2, char acceptor1, char acceptor2,
		      double donor_prob, double acceptor_prob, Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  new->type = CHIMERA_JUNCTION;
  new->nindels = 0;
  new->annotatedp = false;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor1 = donor1;
  new->donor2 = donor2;
  new->acceptor1 = acceptor1;
  new->acceptor2 = acceptor2;
  new->donor_prob = donor_prob;
  new->acceptor_prob = acceptor_prob;

  return new;
}


#ifdef ALLOCATE_UNSOLVED_JUNCTION
T
Junction_new_unsolved (Pathpool_T pathpool) {
  T new = Pathpool_new_junction(pathpool
				pathpool_trace(__FILE__,__LINE__));

  new->type = UNSOLVED_JUNCTION;
  new->nindels = 0;
  new->annotatedp = false;

  new->splice_distance = 0;
  new->sensedir = 0;
  new->donor_prob = 0.0;
  new->acceptor_prob = 0.0;

  return new;
}
#endif


T
Junction_copy (T old, Pathpool_T pathpool) {
  if (old == JUNCTION_UNSOLVED) {
    return JUNCTION_UNSOLVED;
  } else {
    T new = Pathpool_new_junction(pathpool
				  pathpool_trace(__FILE__,__LINE__));

    new->type = old->type;
    new->nindels = old->nindels;
    new->annotatedp = old->annotatedp;

    new->splice_distance = old->splice_distance;
    new->sensedir = old->sensedir;
    new->donor1 = old->donor1;
    new->donor2 = old->donor2;
    new->acceptor1 = old->acceptor1;
    new->acceptor2 = old->acceptor2;
    new->donor_prob = old->donor_prob;
    new->acceptor_prob = old->acceptor_prob;

    return new;
  }
}


List_T
Junction_copy_list (List_T old, Listpool_T listpool, Pathpool_T pathpool) {
  List_T new = NULL, p;

  for (p = old; p != NULL; p = List_next(p)) {
    new = Listpool_push(new,listpool,(void *) Junction_copy((T) List_head(p),pathpool)
			listpool_trace(__FILE__,__LINE__));
  }
  return List_reverse(new);
}


Junctiontype_T
Junction_type (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return UNSOLVED_JUNCTION;
  } else {
    return this->type;
  }
}

bool
Junction_annotatedp (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return false;
  } else {
    return this->annotatedp;
  }
}

char *
Junction_typestring (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return "Unsolved";
  } else {
    switch (this->type) {
    case NO_JUNCTION: return "None";
    case INS_JUNCTION: return "Insertion";
    case DEL_JUNCTION: return "Deletion";
    case SPLICE_JUNCTION: return "Splice";
    case CHIMERA_JUNCTION: return "Chimera";
    case AMB_JUNCTION: return "Amb";
    case END_JUNCTION: return "End";
    case UNSOLVED_JUNCTION: return "Unsolved";
    }
  }
  return (char *) NULL;
}

double
Junction_prob (T this) {
  if (this == NULL) {
    return 0.0;
  } else {
    return this->donor_prob + this->acceptor_prob;
  }
}

int
Junction_sensedir (T this) {
  return this->sensedir;
}

double
Junction_donor_prob (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0.0;
  } else {
    return this->donor_prob;
  }
}

double
Junction_acceptor_prob (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0.0;
  } else {
    return this->acceptor_prob;
  }
}

double
Junction_splice_score (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0.0;
  } else {
    return this->donor_prob + this->acceptor_prob;
  }
}

int
Junction_nindels (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0;
  } else {
    return this->nindels;
  }
}

int
Junction_adj (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0;
  } else if (this->type == DEL_JUNCTION) {
    return +this->nindels;
  } else if (this->type == INS_JUNCTION) {
    return -this->nindels;
  } else if (this->type == SPLICE_JUNCTION) {
    return this->splice_distance;
  } else {
    fprintf(stderr,"Unexpected Junction type %d\n",this->type);
    abort();
    return 0;
  }
}

int
Junction_ninserts (T this) {
  if (this == JUNCTION_UNSOLVED) {
    return 0;
  } else if (this->type == INS_JUNCTION) {
    return this->nindels;
  } else {
    return 0;
  }
}

int
Junction_total_ninserts (List_T list) {
  int ninserts = 0;
  T this;
  List_T p;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    if (this == JUNCTION_UNSOLVED) {
      /* Skip */
    } else if (this->type == INS_JUNCTION) {
      ninserts += this->nindels;
    }
  }

  return ninserts;
}



#if 0
static char complCode[128] = COMPLEMENT_LC;

static char *
make_complement_inplace (char *sequence, unsigned int length) {
  char temp;
  unsigned int i, j;

  for (i = 0, j = length-1; i < length/2; i++, j--) {
    temp = complCode[(int) sequence[i]];
    sequence[i] = complCode[(int) sequence[j]];
    sequence[j] = temp;
  }
  if (i == j) {
    sequence[i] = complCode[(int) sequence[i]];
  }

  return sequence;
}
#endif


/* Called only by Path_print_sam and Path_print_alignment, so we need to treat only the plusp case */
char *
Junction_deletion_string (T this, Univcoord_T deletionpos) {
  char *deletion_string;
  
  /* printf("Entered Junction_deletion_string\n"); */
  /* printf("deletionpos = %u\n",this->deletionpos); */

  deletion_string = (char *) MALLOC((this->nindels+1)*sizeof(char));
  Genome_fill_buffer(deletionpos,this->nindels,deletion_string);
#if 0
  if (plusp == false) {
    make_complement_inplace(deletion_string,this->nindels);
  }
#endif

  /* printf("string = %s\n",deletion_string); */
  return deletion_string;
}


Chrpos_T
Junction_splice_distance (T this) {
  return this->splice_distance;
}

void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob) {
  assert(distance != 0);
  this->splice_distance = distance;
  this->donor_prob = donor_prob;
  this->acceptor_prob = acceptor_prob;
  this->annotatedp = false;

  return;
}

void
Junction_set_ambiguous (T this) {
  this->splice_distance = 0;
  this->annotatedp = false;

  return;
}
