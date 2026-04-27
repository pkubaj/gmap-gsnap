/* $Id: d69ca8ddb89d064622d58ffc91e71f69f70b67de $ */
#ifndef JUNCTION_INCLUDED
#define JUNCTION_INCLUDED

typedef enum {NO_JUNCTION, INS_JUNCTION, DEL_JUNCTION, SPLICE_JUNCTION,
	      CHIMERA_JUNCTION, AMB_JUNCTION, END_JUNCTION, UNSOLVED_JUNCTION} Junctiontype_T;

typedef struct Junction_T *Junction_T;

#include "types.h"
#include "genomicpos.h"
#include "univcoord.h"
#include "bool.h"
#include "genome.h"
#include "list.h"
#include "listpool.h"
#include "pathpool.h"


#define T Junction_T
struct T {
  Junctiontype_T type;
  int nindels;			/* Should be positive */
  /* Univcoord_T deletionpos; -- Computed now at print time */

  Chrpos_T splice_distance;
  int sensedir;

  char donor1;
  char donor2;
  char acceptor1;
  char acceptor2;
  
  double donor_prob;
  double acceptor_prob;
};


extern void
Junction_print (T this);
extern void
Junction_print_list (List_T list);

extern void
Junction_free (T *old, Pathpool_T pathpool);
extern void
Junction_list_gc (List_T *list, Listpool_T listpool, Pathpool_T pathpool);

extern T
Junction_new_insertion (int nindels, Pathpool_T pathpool);
extern T
Junction_new_deletion (int nindels, Pathpool_T pathpool);
extern T
Junction_new_splice (Chrpos_T splice_distance, int sensedir, double donor_prob, double acceptor_prob,
		     Pathpool_T pathpool);
extern T
Junction_new_ambig_splice (int sensedir, double donor_prob, double acceptor_prob, Pathpool_T pathpool);

extern T
Junction_new_chimera (char donor1, char donor2, char acceptor1, char acceptor2,
		      double donor_prob, double acceptor_prob, Pathpool_T pathpool);

/* More efficient to use NULL to represent an unsolved junction.
   Saves on numerous calls to malloc */
#ifdef ALLOCATE_UNSOLVED_JUNCTION
extern T
Junction_new_unsolved (Pathpool_T pathpool);
#else
#define JUNCTION_UNSOLVED (Junction_T) NULL
#endif

extern T
Junction_copy (T old, Pathpool_T pathpool);
extern List_T
Junction_copy_list (List_T old, Listpool_T listpool, Pathpool_T pathpool);


extern Junctiontype_T
Junction_type (T this);
extern char *
Junction_typestring (T this);
extern int
Junction_sensedir (T this);
extern double
Junction_prob (T this);
extern double
Junction_donor_prob (T this);
extern double
Junction_acceptor_prob (T this);
extern double
Junction_splice_score (T this);

extern int
Junction_nindels (T this);
extern int
Junction_adj (T this);
extern int
Junction_ninserts (T this);
extern int
Junction_total_ninserts (List_T list);

extern char *
Junction_deletion_string (T this, Univcoord_T deletionpos);
extern Chrpos_T
Junction_splice_distance (T this);
extern void
Junction_set_unambiguous (T this, Chrpos_T distance, double donor_prob, double acceptor_prob);
extern void
Junction_set_ambiguous (T this);

#undef T
#endif

