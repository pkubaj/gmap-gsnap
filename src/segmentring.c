#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "segmentring.h"

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"
#include "univcoord.h"


/* Segmentring_smooth */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif



#define T Segmentring_T
struct T {
  Segmentlist_T head;
  Univcoord_T middle_univdiagonal;
  int length;
};


T
Segmentring_new (Univcoord_T middle_univdiagonal) {
  T new = (T) MALLOC(sizeof(*new));

  new->head = (Segmentlist_T) NULL;
  new->middle_univdiagonal = middle_univdiagonal;
  new->length = 0;

  return new;
}


void
Segmentring_free (T *old) {
  Segmentlist_T p, q;
  int n;

  if ((p = (*old)->head) != NULL) {
    n = (*old)->length;
    while (n-- > 0) {
      q = p->rlink;
      FREE(p);
      p = q;
    }
  }

  FREE(*old);
  return;
}


Segmentlist_T
Segmentring_head (T this) {
  return this->head;
}

Segmentlist_T
Segmentring_tail (T this) {
  return this->head->llink;
}

Univcoord_T
Segmentring_middle_univdiagonal (T this) {
  return this->middle_univdiagonal;
}

int
Segmentring_length (T this) {
  return this->length;
}


void
Segmentring_print (T this) {
  int i;
  Segmentlist_T p = this->head;

  for (i = 0; i < this->length; i++) {
    Segmentlist_print(p);
    p = p->rlink;
  }

  printf("\n");

  return;
}


void
Segmentring_get (int *index, int *qstart, int *qend, Univcoord_T *univdiagonal,
		 int *count, T this, int i) {
  Segmentlist_T q;
  int n;

  /* Get the ith node */
  q = this->head;
  if (i <= this->length/2) {
    n = i;
    while (n-- > 0) {
      q = q->rlink;
    }
  } else {
    n = this->length - i;
    while (n-- > 0) {
      q = q->llink;
    }
  }

  *index = q->index;
  *qstart = q->qstart;
  *qend = q->qend;
  *univdiagonal = q->univdiagonal;
  *count = q->count;

  return;
}


void
Segmentring_put (T this, int i, int index, int qstart, int qend, Univcoord_T univdiagonal,
		 int count) {
  Segmentlist_T q;
  int n;

  /* Get the ith node */
  q = this->head;
  if (i <= this->length/2) {
    n = i;
    while (n-- > 0) {
      q = q->rlink;
    }
  } else {
    n = this->length - i;
    while (n-- > 0) {
      q = q->llink;
    }
  }
  
  q->index = index;
  q->qstart = qstart;
  q->qend = qend;
  q->univdiagonal = univdiagonal;
  q->count = count;

  return;
}


void
Segmentring_addhi (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
		   int count) {
  Segmentlist_T p, q;

  p = (Segmentlist_T) MALLOC(sizeof(*p));
  if ((q = this->head) == NULL) {
    /* Make p the only value of the ring */
    this->head = p->llink = p->rlink = p;
  } else {
    /* Insert p to the left of q */
    p->llink = q->llink;
    q->llink->rlink = p;
    p->rlink = q;
    q->llink = p;
  }

  this->length++;
  p->index = index;
  p->qstart = qstart;
  p->qend = qend;
  p->univdiagonal = univdiagonal;
  p->count = count;

  return;
}


void
Segmentring_addlo (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
		   int count) {
  Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);
  this->head = this->head->llink;
  return;
}


void
Segmentring_add (T this, int pos, int index, int qstart, int qend, Univcoord_T univdiagonal,
		 int count) {
  Segmentlist_T p, q;
  int n, i;

  if (pos == 1 || pos == -this->length) {
    Segmentring_addlo(this,index,qstart,qend,univdiagonal,count);

  } else if (pos == 0 || pos == this->length + 1) {
    Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);

  } else {
    if (pos < 0) {
      i = pos + this->length;
    } else {
      i = pos - 1;
    }

    /* Get the ith node */
    q = this->head;
    if (i <= this->length/2) {
      n = i;
      while (n-- > 0) {
	q = q->rlink;
      }
    } else {
      n = this->length - i;
      while (n-- > 0) {
	q = q->llink;
      }
    }

    p = (Segmentlist_T) MALLOC(sizeof(*p));
    
    /* Insert p to the left of q */
    p->llink = q->llink;
    q->llink->rlink = p;
    p->rlink = q;
    q->llink = p;

    this->length++;
    p->index = index;
    p->qstart = qstart;
    p->qend = qend;
    p->univdiagonal = univdiagonal;
    p->count = count;
  }

  return;
}


bool
Segmentring_insert_qstart (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
			   int count, int max_overlap) {
  Segmentlist_T prev, next, head, new;
  int overlap_start, overlap_end;
  int i;

  if (this->length == 0) {
    Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);
    return true;

  } else if (this->length == 1) {
    head = this->head;
    
    if (qstart < head->qstart && qend < head->qend) {
      Segmentring_addlo(this,index,qstart,qend,univdiagonal,count);
      return true;
    } else {
      return false;
    }

  } else {
    prev = this->head->llink;	/* tail */
    if (qstart >= prev->qstart || qend >= prev->qend) {
      /* On qend side of middle_univdiagonal */
      return false;

    } else {
      for (i = this->length - 2; i >= 0; i--) {
	next = prev->llink;
      
	if (qstart > next->qstart && qend < prev->qend) {
	  /* Fits between next and prev.  Add if the overlaps are reasonable */
	  overlap_start = next->qend - qstart;
	  overlap_end = qend - prev->qstart;
	
	  if (overlap_start <= max_overlap && overlap_end <= max_overlap) {
	    new = (Segmentlist_T) MALLOC(sizeof(*new));
	
	    /* Insert new to the left of prev */
	    new->llink = prev->llink;
	    prev->llink->rlink = new;
	    new->rlink = prev;
	    prev->llink = new;
	    
	    this->length++;
	    new->index = index;
	    new->qstart = qstart;
	    new->qend = qend;
	    new->univdiagonal = univdiagonal;
	    new->count = count;
	    
	    return true;
	  }
	}

	prev = next;
      }

      if (qstart >= prev->qstart) {
	return false;
      } else if (qend >= prev->qend) {
	return false;
      } else {
	overlap_end = qend - prev->qstart;
	if (overlap_end <= max_overlap) {
	  Segmentring_addlo(this,index,qstart,qend,univdiagonal,count);
	  return true;
	} else {
	  return false;
	}
      }
    }
  }
}


bool
Segmentring_insert_qend (T this, int index, int qstart, int qend, Univcoord_T univdiagonal,
			 int count, int max_overlap) {
  Segmentlist_T prev, next, head, new;
  int overlap_start, overlap_end;
  int i;

  if (this->length == 0) {
    Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);
    return true;

  } else if (this->length == 1) {
    head = this->head;
    
    if (qstart > head->qstart && qend > head->qend) {
      Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);
      return true;
    } else {
      return false;
    }

  } else {
    prev = this->head;
    if (qstart <= prev->qstart || qend <= prev->qend) {
      /* On qstart side of middle_univdiagonal */
      return false;
      
    } else {
      for (i = 0; i < this->length - 1; i++) {
	next = prev->rlink;
      
	if (qstart > prev->qstart && qend < next->qend) {
	  /* Fits between prev and next.  Add if the overlaps are reasonable */
	  overlap_start = prev->qend - qstart;
	  overlap_end = qend - next->qstart;
	
	  if (overlap_start <= max_overlap && overlap_end <= max_overlap) {
	    new = (Segmentlist_T) MALLOC(sizeof(*new));
	
	    /* Insert new to the left of next */
	    new->llink = next->llink;
	    next->llink->rlink = new;
	    new->rlink = next;
	    next->llink = new;
	
	    this->length++;
	    new->index = index;
	    new->qstart = qstart;
	    new->qend = qend;
	    new->univdiagonal = univdiagonal;
	    new->count = count;
	
	    return true;
	  }
	}

	prev = next;
      }

      if (qstart <= prev->qstart) {
	return false;
      } else if (qend <= prev->qend) {
	return false;
      } else {
	overlap_start = prev->qend - qstart;
	if (overlap_start <= max_overlap) {
	  Segmentring_addhi(this,index,qstart,qend,univdiagonal,count);
	  return true;
	} else {
	  return false;
	}
      }
    }
  }
}


/* Dynamic programming */
List_T
Segmentring_smooth (bool *usedp, T this, Univdiagpool_T univdiagpool) {
  List_T chain = NULL;
  Segmentlist_T global_ptr, best_ptr, p, q;
  int i, j;
  int global_score, best_score, score;

  /* Initialize */
  j = 0;
  q = this->head;
  global_score = q->score = q->count;
  q->ptr = (Segmentlist_T) NULL;
  global_ptr = q;
  

  for (j++, q = q->rlink; j < this->length; j++, q = q->rlink) {
    /* Handle the first one */
    i = j - 1;
    p = q->llink;
    best_score = p->score;
    if (p->univdiagonal > q->univdiagonal) {
      /* Insertion: Take penalty */
      best_score -= (int) (p->univdiagonal - q->univdiagonal);
    } else {
      /* Deletion or splice: Ignore penalty as long as we are monotonically increasing */
      /* best_score -= 0; */
    }
    best_ptr = p;

    /* Handle the rest */
    for (i--, p = p->llink; i >= 0; i--, p = p->llink) {
      score = p->score;
      if (p->univdiagonal > q->univdiagonal) {
	/* Insertion: Take penalty */
	score -= (int) (p->univdiagonal - q->univdiagonal);
      } else {
	/* Deletion or splice: Ignore penalty as long as we are monotonically increasing */
	/* best_score -= 0; */
      }
      
      if (score > best_score) {
	best_ptr = p;
	best_score = score;
      }
    }

    /* Store best results */
    q->score = q->count + best_score;
    q->ptr = best_ptr;
    debug9(printf("At %d..%d %u, score is %d, best ptr is %d..%d %u\n",
		  q->qstart,q->qend,q->univdiagonal,q->score,
		  best_ptr->qstart,best_ptr->qend,best_ptr->univdiagonal));

    if (q->score > global_score) {
      global_score = q->score;
      global_ptr = q;
    }
  }

  /* Backtrack */
  debug9(printf("Best path\n"));
  for (q = global_ptr; q != NULL; q = q->ptr) {
    debug9(printf("%d..%d %u score:%d\n",q->qstart,q->qend,q->univdiagonal,q->score));
    chain = Univdiagpool_push(chain,univdiagpool,
			      q->qstart,q->qend,/*nmismatches*/-1,q->univdiagonal
			      univdiagpool_trace(__FILE__,__LINE__));
    usedp[q->index] = true;
  }

  return chain;
}

