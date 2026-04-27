static char rcsid[] = "$Id: d98901be34c5550c1fd076939e82ac032889b2f2 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "exon.h"
#include <stdio.h>
#include <string.h>		/* For memcpy */

#include "mem.h"


/* Exon_list_validp */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


/* Printing routines */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


#define T Exon_T

#if 0
/* For TEST_EXONS in trpath.c */
bool
Exon_list_equal (List_T a, List_T b) {
  T exon1, exon2;

  while (a != NULL && b != NULL) {
    exon1 = (T) List_head(a);
    exon2 = (T) List_head(b);
    if (exon1->exoni != exon2->exoni) {
      return false;
    } else if (exon1->firstchar != exon2->firstchar) {
      return false;
    } else if (exon1->lastchar != exon2->lastchar) {
      return false;
    } else {
      a = List_next(a);
      b = List_next(b);
    }
  }

  if (a != NULL) {
    return false;
  } else if (b != NULL) {
    return false;
  } else {
    return true;
  }
}
#endif


void
Exon_free (T *old, Transcriptpool_T transcriptpool) {
  Transcriptpool_free_exon(&(*old),transcriptpool
			   transcriptpool_trace(__FILE__,__LINE__)); /* Allocated by Transcriptpool_new_exon */
  return;
}

void
Exon_list_gc (List_T *exons, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T p;
  T exon;

  for (p = *exons; p != NULL; p = List_next(p)) {
    exon = (T) List_head(p);
    Exon_free(&exon,transcriptpool);
  }
  Listpool_free_list(&(*exons),listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */

  return;
}


T
Exon_new (char firstchar, int exoni, char lastchar,
	  Transcriptpool_T transcriptpool) {
  T new = Transcriptpool_new_exon(transcriptpool
				  transcriptpool_trace(__FILE__,__LINE__));

  new->exoni = exoni;
  new->firstchar = firstchar;
  new->lastchar = lastchar;

  return new;
}

static T
Exon_copy (T old, Transcriptpool_T transcriptpool) {
  T new = Transcriptpool_new_exon(transcriptpool
				  transcriptpool_trace(__FILE__,__LINE__));
  
#if 0
  new->exoni = old->exoni;
  new->firstchar = old->firstchar;
  new->lastchar = old->lastchar;
#else
  memcpy(new,old,sizeof(*new));
#endif

  return new;
}

List_T
Exon_list_copy (List_T old, Transcriptpool_T transcriptpool,
		Listpool_T listpool) {
  List_T new = NULL;

  while (old != NULL) {
    new = Listpool_push(new,listpool,
			(void *) Exon_copy((T) List_head(old),transcriptpool)
			listpool_trace(__FILE__,__LINE__));
    old = List_next(old);
  }

  return List_reverse(new);
}


/* Assumes exons are in reverse order */
bool
Exon_list_consecutivep (List_T exons) {
  T exon;
  int next_exoni;

  debug(printf("Entering Exon_list_consecutivep\n"));

  if (exons == NULL) {
    return true;

  } else {
    exon = (T) List_head(exons);
    next_exoni = exon->exoni;
    exons = List_next(exons);

    while (exons != NULL) {
      exon = (T) List_head(exons);
      if (exon->exoni != next_exoni - 1) {
	/* Non-consecutive exons */
	return false;
      } else {
	next_exoni = exon->exoni;
	exons = List_next(exons);
      }
    }

    return true;
  }
}



bool
Exon_list_validp (bool *repairablep, List_T exons) {
  bool validp;
  T exon;
  int last_exoni;

  debug(printf("Entering Exon_list_validp\n"));

  if (exons == NULL) {
    debug(printf("(1) exons is NULL\n"));
    *repairablep = false;
    debug(printf("Returning false\n"));
    return false;

  } else if (List_length(exons) == 1) {
    /* Single exon */
    *repairablep = true;

    exon = (T) List_head(exons);
    if (exon->firstchar == 'i') {
      debug(printf("(2) firstchar is i\n"));
      *repairablep = false;
      return false;

    } else if (exon->firstchar == 'x') {
      debug(printf("(3) firstchar is x\n"));
      /* Still repairable */
      return false;

    } else if (exon->lastchar == 'i') {
      debug(printf("(4) lastchar is i\n"));
      *repairablep = false;
      return false;

    } else if (exon->lastchar == 'x') {
      debug(printf("(5) lastchar is x\n"));
      /* Still repairable */
      return false;

    } else if (exon->firstchar == 'u' && exon->lastchar == 'u') {
      debug(printf("(6) firstchar is u and lastchar is u\n"));
      *repairablep = false;
      return false;

    } else {
      return true;
    }

  } else {
    *repairablep = true;
    validp = true;

    /* First exon */
    exon = (T) List_head(exons);
    if (exon->firstchar == 'i') {
      debug(printf("(6) firstchar is i\n"));
      *repairablep = false;
      validp = false;

    } else if (exon->firstchar == 'x') {
      debug(printf("(7) firstchar is x\n"));
      /* Still repairable */
      validp = false;

#if 0
    } else if (exon->firstchar == 'u') {
      /* Allow extension from 5' UTR */
#endif      

    } else if (exon->lastchar == 'i' || exon->lastchar == 'x' || exon->lastchar == 'u') {
      debug(printf("(8) lastchar is %c\n",exon->lastchar));
      *repairablep = false;
      validp = false;
    }

    last_exoni = exon->exoni;
    exons = List_next(exons);

    /* Middle exons */
    while (List_next(exons) != NULL) {
      exon = (T) List_head(exons);
      if (exon->firstchar == 'i' || exon->firstchar == 'x' || exon->firstchar == 'u') {
	debug(printf("(9) firstchar is %c\n",exon->firstchar));
	*repairablep = false;
	validp = false;

      } else if (exon->lastchar == 'i' || exon->lastchar == 'x' || exon->lastchar == 'u') {
	debug(printf("(10) lastchar is %c\n",exon->lastchar));
	*repairablep = false;
	validp = false;

      } else if (exon->exoni != last_exoni + 1) {
	/* Non-consecutive exons */
	debug(printf("(11) Non-consecutive exons\n"));
	*repairablep = false;
	validp = false;
      }

      last_exoni = exon->exoni;
      exons = List_next(exons);
    }

    /* Last exon */
    exon = (T) List_head(exons);
    if (exon->firstchar == 'i' || exon->firstchar == 'x' || exon->firstchar == 'u') {
      debug(printf("(12) firstchar is %c\n",exon->firstchar));
      *repairablep = false;
      validp = false;

    } else if (exon->lastchar == 'i') {
      debug(printf("(13) lastchar is i\n"));
      *repairablep = false;
      validp = false;

    } else if (exon->lastchar == 'x') {
      debug(printf("(14) lastchar is x\n"));
      /* Still repairable */
      validp = false;

#if 0
    } else if (exon->lastchar == 'u') {
      /* Allow extension into 3' UTR */
#endif      

    } else if (exon->exoni != last_exoni + 1) {
      /* Non-consecutive exons */
      debug(printf("(15) Non-consecutive exons\n"));
      *repairablep = false;
      validp = false;
    }

    debug(printf("Returning validp %d\n",validp));
    return validp;
  }
}


void
Exon_print_list (Filestring_T fp, List_T exons) {
  List_T p;
  T exon;

  if (exons != NULL) {
    exon = (T) List_head(exons);
    FPRINTF(fp,"%c%d%c",exon->firstchar,exon->exoni + /*1-based*/1,exon->lastchar);
  
    for (p = List_next(exons); p != NULL; p = List_next(p)) {
      exon = (T) List_head(p);
      FPRINTF(fp,"|%c%d%c",exon->firstchar,exon->exoni + /*1-based*/1,exon->lastchar);
    }
  }

  return;
}


void
Exon_print_list_stdout (List_T exons) {
  List_T p;
  T exon;

  if (exons != NULL) {
    exon = (T) List_head(exons);
    printf("%c%d%c",exon->firstchar,exon->exoni + /*1-based*/1,exon->lastchar);
  
    for (p = List_next(exons); p != NULL; p = List_next(p)) {
      exon = (T) List_head(p);
      printf("|%c%d%c",exon->firstchar,exon->exoni + /*1-based*/1,exon->lastchar);
    }
  }

  return;
}


