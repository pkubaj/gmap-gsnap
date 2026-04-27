static char rcsid[] = "$Id: bede0610b8a97ae142ad31be7ef2cc9235ebcf73 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "resulthr.h"
#include <stdlib.h>
#include "assert.h"
#include "mem.h"
#include "path.h"
#include "pathpair.h"


/* Assignment of resulttype */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif



#define T Result_T
struct T {
  Resulttype_T resulttype;
  int id;
  bool tr_consistent_p;

  void **array;
  int npaths_primary;
  int npaths_altloc;
  int first_absmq;
  int second_absmq;
  void **array2;
  int npaths2_primary;
  int npaths2_altloc;
  int first_absmq2;
  int second_absmq2;
  double worker_runtime;

  SAM_split_output_type split_output;
  int strlength;
  /* char *string; */
};


Resulttype_T
Result_resulttype (T this) {
  return this->resulttype;
}


bool
Result_tr_consistent_p (T this) {
  return this->tr_consistent_p;
}


char *
Pairtype_string (Pairtype_T pairtype) {
  switch (pairtype) {

  case CONCORDANT: return "concordant";
  case PAIRED_UNSPECIFIED: return "paired_unspecified";
  case PAIRED_INVERSION: return "paired_scramble";
  case PAIRED_SCRAMBLE: return "paired_scramble";
  case PAIRED_TOOLONG: return "paired_toolong";
  case CONCORDANT_TRANSLOCATIONS: return "concordant_translocations";
  case UNPAIRED: return "unpaired";
  case UNSPECIFIED: return "unspecified";
  default: 
    fprintf(stderr,"Unknown pairtype %d\n",pairtype);
    abort();
  }
  return "";
}


char *
Resulttype_string (Resulttype_T resulttype) {
  switch (resulttype) {
  case SINGLEEND_NOMAPPING: return "singleend_nomapping";
  case PAIREDEND_NOMAPPING: return "pairedend_nomapping";
  case SINGLEEND_UNIQ: return "singleend_uniq";
  case SINGLEEND_TRANSLOC: return "singleend_transloc";
  case SINGLEEND_MULT: return "singleend_mult";
  case PAIRED_UNIQ_INV: return "paired_uniq_inv";
  case PAIRED_UNIQ_SCR: return "paired_uniq_scr";
  case PAIRED_UNIQ_TOOLONG: return "paired_uniq_long";
  case PAIRED_MULT: return "paired_mult";
  case CONCORDANT_UNIQ: return "concordant_uniq";
  case CONCORDANT_TRANSLOC: return "concordant_transloc";
  case CONCORDANT_MULT: return "concordant_mult";
  case HALFMAPPING_UNIQ: return "halfmapping_uniq";
  case HALFMAPPING_TRANSLOC: return "halfmapping_transloc";
  case HALFMAPPING_MULT: return "halfmapping_mult";
  case UNPAIRED_UNIQ: return "unpaired_uniq";
  case UNPAIRED_TRANSLOC: return "unpaired_transloc";
  case UNPAIRED_MULT: return "unpaired_mult";
  default: 
    fprintf(stderr,"Unknown resulttype %d\n",resulttype);
    abort();
  }
  return "";
}


int
Result_id (T this) {
  return this->id;
}


void **
Result_array (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, T this) {
  *npaths_primary = this->npaths_primary;
  *npaths_altloc = this->npaths_altloc;
  *first_absmq = this->first_absmq;
  *second_absmq = this->second_absmq;
  return this->array;
}


/* For second end, when not paired with first end */
void **
Result_array2 (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq, T this) {
  *npaths_primary = this->npaths2_primary;
  *npaths_altloc = this->npaths2_altloc;
  *first_absmq = this->first_absmq2;
  *second_absmq = this->second_absmq2;
  return this->array2;
}

double
Result_worker_runtime (T this) {
  return this->worker_runtime;
}


/* resulttype can be SINGLEEND_NOMAPPING, SINGLEEND_TRANSLOC, SINGLEEND_UNIQ, SINGLEEND_MULT */
T
Result_single_read_new (int id, void **resultarray, int npaths_primary, int npaths_altloc,
			int first_absmq, int second_absmq) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Path_T path;
  int pathi;

  new->tr_consistent_p = false;

  if (npaths_primary + npaths_altloc == 0) {
    new->resulttype = SINGLEEND_NOMAPPING;
    
  } else {
    if (npaths_primary + npaths_altloc > 1) {
      new->resulttype = SINGLEEND_MULT;
    } else {
      new->resulttype = SINGLEEND_UNIQ;
    }
    
    for (pathi = 0; pathi < npaths_primary + npaths_altloc; pathi++) {
      path = ((Path_T *) resultarray)[pathi];
      if (path->fusion_querystart_junction != NULL ||
	  path->fusion_queryend_junction != NULL) {
	new->resulttype = SINGLEEND_TRANSLOC;

      } else if (path->transcripts != NULL) {
	new->tr_consistent_p = true;
      }
    }
  }

  new->id = id;
  new->array = resultarray;
  new->npaths_primary = npaths_primary;
  new->npaths_altloc = npaths_altloc;
  new->first_absmq = first_absmq;
  new->second_absmq = second_absmq;

  return new;
}


/* resulttype can be CONCORDANT_TRANSLOC, CONCORDANT_UNIQ, CONCORDANT_MULT, PAIRED_UNIQ_INV, PAIRED_UNIQ_SCR, PAIRED_UNIQ_TOOLONG, PAIRED_MULT */
/* final_pairtype can be CONCORDANT_TRANSLOCATIONS, CONCORDANT, PAIRED_INVERSION, PAIRED_SCRAMBLE, PAIRED_TOOLONG */
/* final_pairtype cannot be UNPAIRED */
T
Result_paired_read_new (int id, void **resultarray, int npaths_primary, int npaths_altloc,
			int first_absmq, int second_absmq, Pairtype_T final_pairtype) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Pathpair_T pathpair;
  Path_T path5, path3;
  int pathi;

  /* printf("Entered Result_paired_read_new with %d primary, %d altloc\n",
     npaths_primary,npaths_altloc); */
  /* printf("final_pairtype %d\n",final_pairtype); */

  new->tr_consistent_p = false;

  if (final_pairtype == CONCORDANT_TRANSLOCATIONS) {
    new->resulttype = CONCORDANT_TRANSLOC;

  } else if (final_pairtype == CONCORDANT) {
    if (npaths_primary + npaths_altloc > 1) {
      new->resulttype = CONCORDANT_MULT;
    } else {
      new->resulttype = CONCORDANT_UNIQ;
    }

    for (pathi = 0; pathi < npaths_primary + npaths_altloc; pathi++) {
      pathpair = ((Pathpair_T *) resultarray)[0];
      path5 = pathpair->path5;
      path3 = pathpair->path3;
      
      if (path5->fusion_querystart_junction != NULL ||
	  path5->fusion_queryend_junction != NULL ||
	  path3->fusion_querystart_junction != NULL ||
	  path3->fusion_queryend_junction != NULL) {
	new->resulttype = CONCORDANT_TRANSLOC;

      } else if (pathpair->transcript_concordant_p == true) {
	new->tr_consistent_p = true;
      }
    }

  } else if (npaths_primary + npaths_altloc > 1) {
    new->resulttype = PAIRED_MULT;

  } else if (final_pairtype == PAIRED_INVERSION) {
    new->resulttype = PAIRED_UNIQ_INV;

  } else if (final_pairtype == PAIRED_SCRAMBLE) {
    new->resulttype = PAIRED_UNIQ_SCR;

  } else if (final_pairtype == PAIRED_TOOLONG) {
    new->resulttype = PAIRED_UNIQ_TOOLONG;

  } else {
    fprintf(stderr,"final_pairtype %d not recognized\n",final_pairtype);
    abort();
  }

  new->id = id;
  new->array = resultarray;
  new->npaths_primary = npaths_primary;
  new->npaths_altloc = npaths_altloc;
  new->first_absmq = first_absmq;
  new->second_absmq = second_absmq;

  /* printf("resulttype %d\n",new->resulttype); */

  return new;
}


/* resulttype can be PAIREDEND_NOMAPPING, HALFMAPPING_TRANSLOC, HALFMAPPING_UNIQ, HALFMAPPING_MULT,
   UNPAIRED_TRANSLOC, UNPAIRED_UNIQ, UNPAIRED_MULT */
T
Result_paired_as_singles_new (int id, void **hits5, int npaths5_primary, int npaths5_altloc,
			      int first_absmq5, int second_absmq5,
			      void **hits3, int npaths3_primary, int npaths3_altloc,
			      int first_absmq3, int second_absmq3) {
  T new = (T) MALLOC_OUT(sizeof(*new));
  Path_T path5, path3;
  int npaths5, npaths3, pathi;

  npaths5 = npaths5_primary + npaths5_altloc;
  npaths3 = npaths3_primary + npaths3_altloc;
  debug(printf("npaths5 = %d, npaths3 = %d\n",npaths5,npaths3));

  new->tr_consistent_p = false;

  if (npaths5 == 0 && npaths3 == 0) {
    new->resulttype = PAIREDEND_NOMAPPING;

  } else if (npaths5 == 0) {
    if (npaths3 > 1) {
      new->resulttype = HALFMAPPING_MULT;
    } else {
      new->resulttype = HALFMAPPING_UNIQ;
    }

    for (pathi = 0; pathi < npaths3; pathi++) {
      path3 = (Path_T) hits3[pathi];
      if (path3->fusion_querystart_junction != NULL ||
	  path3->fusion_queryend_junction != NULL) {
	new->resulttype = HALFMAPPING_TRANSLOC;
      }
    }

  } else if (npaths3 == 0) {
    if (npaths5 > 1) {
      new->resulttype = HALFMAPPING_MULT;
    } else {
      new->resulttype = HALFMAPPING_UNIQ;
    }
    
    for (pathi = 0; pathi < npaths5; pathi++) {
      path5 = (Path_T) hits5[pathi];
      if (path5->fusion_querystart_junction != NULL ||
	  path5->fusion_queryend_junction != NULL) {
	new->resulttype = HALFMAPPING_TRANSLOC;
      }
    }

  } else {
    if (npaths5 > 1 || npaths3 > 1) {
      new->resulttype = UNPAIRED_MULT;
    } else {
      new->resulttype = UNPAIRED_UNIQ;
    }
    
    for (pathi = 0; pathi < npaths5; pathi++) {
      path5 = (Path_T) hits5[pathi];
      if (path5->fusion_querystart_junction != NULL ||
	  path5->fusion_queryend_junction != NULL) {
	new->resulttype = UNPAIRED_TRANSLOC;
      }
    }

    for (pathi = 0; pathi < npaths3; pathi++) {
      path3 = (Path_T) hits3[pathi];
      if (path3->fusion_querystart_junction != NULL ||
	  path3->fusion_queryend_junction != NULL) {
	new->resulttype = UNPAIRED_TRANSLOC;
      }
    }
  }

  new->id = id;
  new->array = hits5;
  new->npaths_primary = npaths5_primary;
  new->npaths_altloc = npaths5_altloc;
  new->first_absmq = first_absmq5;
  new->second_absmq = second_absmq5;
  new->array2 = hits3;
  new->npaths2_primary = npaths3_primary;
  new->npaths2_altloc = npaths3_altloc;
  new->first_absmq2 = first_absmq3;
  new->second_absmq2 = second_absmq3;

  return new;
}

void
Result_free (T *old, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
	     Hitlistpool_T hitlistpool) {
  int i;
  Path_T path;
  Pathpair_T pathpair;

  /* printf("Entered Result_free with resulttype %s, %d primary, %d altloc\n",
     Resulttype_string((*old)->resulttype),(*old)->npaths_primary,(*old)->npaths_altloc); */

  switch ((*old)->resulttype) {
  case SINGLEEND_NOMAPPING: case PAIREDEND_NOMAPPING:
    /* It is possible to have results that are not printed because of various flags */
    for (i = 0; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
      path = (*old)->array[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT((*old)->array);
    break;

  case SINGLEEND_UNIQ: case SINGLEEND_TRANSLOC: case SINGLEEND_MULT:
    for (i = 0; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
      path = (*old)->array[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT((*old)->array);
    break;

  case PAIRED_UNIQ_INV: case PAIRED_UNIQ_SCR: case PAIRED_UNIQ_TOOLONG: case PAIRED_MULT:
  case CONCORDANT_UNIQ: case CONCORDANT_TRANSLOC: case CONCORDANT_MULT:
    for (i = 0; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
      pathpair = (*old)->array[i];
      Pathpair_free(&pathpair,intlistpool,univcoordlistpool,
		    listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT((*old)->array);
    break;

  default:
    for (i = 0; i < (*old)->npaths2_primary + (*old)->npaths2_altloc; i++) {
      path = (*old)->array2[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT((*old)->array2);

    for (i = 0; i < (*old)->npaths_primary + (*old)->npaths_altloc; i++) {
      path = (*old)->array[i];
      Path_free(&path,intlistpool,univcoordlistpool,
		listpool,pathpool,transcriptpool,hitlistpool);
    }
    FREE_OUT((*old)->array);
  }

  FREE_OUT(*old);

  return;
}



