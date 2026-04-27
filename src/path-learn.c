static char rcsid[] = "$Id: stage3hr.c 225533 2022-12-21 00:40:38Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "path-learn.h"

#include <math.h>
#include "junction.h"
#include "orderstat.h"


#define MIN_TWOPASS_SPLICING_SUPPORT 10
#define MIN_TWOPASS_INDEL_SUPPORT 10


/* In multi-threaded mode, assumes that a lock has been obtained */
void
Path_learn_defect_rate (Path_T this, unsigned long long *total_mismatches,
			unsigned long long *total_querylength) {
  *total_mismatches += this->querylength - this->nmatches;
  *total_querylength += this->querylength;
  return;
}


#ifdef TO_FIX
void
Path_learn_splicesites (T this, Univcoordtableuint_T donor_table, Univcoordtableuint_T acceptor_table,
			Univcoordtableuint_T antidonor_table, Univcoordtableuint_T antiacceptor_table) {
  Univcoord_T splicesite_position;
  Endtype_T endtype;
  unsigned int count;

  substring1 = (Substring_T) List_head(this->substrings_1toN);
  if ((endtype = Substring_start_endtype(substring1)) == DON) {
    splicesite_position = Substring_alignstart_trim(substring1);
    if (Substring_plusp(substring1) == true) {
      /* printf("antidonor at %u\n",splicesite_position) */
      count = Univcoordtableuint_get(antidonor_table,splicesite_position);
      Univcoordtableuint_put(antidonor_table,splicesite_position,count + 1);
    } else {
      /* printf("donor at %u\n",splicesite_position) */
      count = Univcoordtableuint_get(donor_table,splicesite_position);
      Univcoordtableuint_put(donor_table,splicesite_position,count + 1);
    }
  } else if (endtype == ACC) {
    splicesite_position = Substring_alignstart_trim(substring1);
    if (Substring_plusp(substring1) == true) {
      /* printf("acceptor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(acceptor_table,splicesite_position);
      Univcoordtableuint_put(acceptor_table,splicesite_position,count + 1);
    } else {
      /* printf("antiacceptor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(antiacceptor_table,splicesite_position);
      Univcoordtableuint_put(antiacceptor_table,splicesite_position,count + 1);
    }
  }

  substringN = (Substring_T) List_head(this->substrings_Nto1);
  if ((endtype = Substring_end_endtype(substringN)) == DON) {
    splicesite_position = Substring_alignend_trim(substringN);
    if (Substring_plusp(substringN) == true) {
      /* printf("donor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(donor_table,splicesite_position);
      Univcoordtableuint_put(donor_table,splicesite_position,count + 1);
    } else {
      /* printf("antidonor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(antidonor_table,splicesite_position);
      Univcoordtableuint_put(antidonor_table,splicesite_position,count + 1);
    }
  } else if (endtype == ACC) {
    splicesite_position = Substring_alignend_trim(substringN);
    if (Substring_plusp(substringN) == true) {
      /* printf("antiacceptor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(antiacceptor_table,splicesite_position);
      Univcoordtableuint_put(antiacceptor_table,splicesite_position,count + 1);
    } else {
      /* printf("acceptor at %u\n",splicesite_position); */
      count = Univcoordtableuint_get(acceptor_table,splicesite_position);
      Univcoordtableuint_put(acceptor_table,splicesite_position,count + 1);
    }
  }
  
  return;
}
#endif


void
Path_learn_introns (Path_T this, Univcoordlist_T *donor_startpoints, Univcoordlist_T *donor_partners,
		    Univcoordlist_T *acceptor_startpoints, Univcoordlist_T *acceptor_partners,
		    Univcoordlist_T *antidonor_startpoints, Univcoordlist_T *antidonor_partners,
		    Univcoordlist_T *antiacceptor_startpoints, Univcoordlist_T *antiacceptor_partners) {
  Intlist_T q;
  Univcoordlist_T p;
  List_T j;
  Junction_T junction;
  Univcoord_T donor_splicesite, acceptor_splicesite, antidonor_splicesite, antiacceptor_splicesite,
    pre_univdiagonal, post_univdiagonal;
  int splicesite_qpos, support;

  q = this->endpoints;
  p = this->univdiagonals;
  for (j = this->junctions; j != NULL; j = List_next(j)) {
    junction = (Junction_T) List_head(j);

    if (Junction_type(junction) == SPLICE_JUNCTION) {
      splicesite_qpos = Intlist_head(Intlist_next(q));

      if (splicesite_qpos < this->querylength - splicesite_qpos) {
	support = splicesite_qpos;
      } else {
	support = this->querylength - splicesite_qpos;
      }

      if (support < MIN_TWOPASS_SPLICING_SUPPORT) {
	/* Don't use */

      } else if (Junction_sensedir(junction) == SENSE_FORWARD) {
	pre_univdiagonal = Univcoordlist_head(p);
	post_univdiagonal = Univcoordlist_head(Univcoordlist_next(p));

	if (this->plusp == true) {
	  donor_splicesite = pre_univdiagonal - this->querylength + splicesite_qpos;
	  acceptor_splicesite = post_univdiagonal - this->querylength + splicesite_qpos;

	  *donor_startpoints = Univcoordlist_push(*donor_startpoints,donor_splicesite);
	  *donor_partners = Univcoordlist_push(*donor_partners,acceptor_splicesite);
	  *acceptor_startpoints = Univcoordlist_push(*acceptor_startpoints,acceptor_splicesite);
	  *acceptor_partners = Univcoordlist_push(*acceptor_partners,donor_splicesite);

	} else {
	  antidonor_splicesite = pre_univdiagonal - this->querylength + splicesite_qpos;
	  antiacceptor_splicesite = post_univdiagonal - this->querylength + splicesite_qpos;

	  *antidonor_startpoints = Univcoordlist_push(*antidonor_startpoints,antidonor_splicesite);
	  *antidonor_partners = Univcoordlist_push(*antidonor_partners,antiacceptor_splicesite);
	  *antiacceptor_startpoints = Univcoordlist_push(*antiacceptor_startpoints,antiacceptor_splicesite);
	  *antiacceptor_partners = Univcoordlist_push(*antiacceptor_partners,antidonor_splicesite);
	}

      } else if (Junction_sensedir(junction) == SENSE_ANTI) {
	pre_univdiagonal = Univcoordlist_head(p);
	post_univdiagonal = Univcoordlist_head(Univcoordlist_next(p));

	if (this->plusp == true) {
	  antidonor_splicesite = post_univdiagonal - this->querylength + splicesite_qpos;
	  antiacceptor_splicesite = pre_univdiagonal - this->querylength + splicesite_qpos;

	  *antidonor_startpoints = Univcoordlist_push(*antidonor_startpoints,antidonor_splicesite);
	  *antidonor_partners = Univcoordlist_push(*antidonor_partners,antiacceptor_splicesite);
	  *antiacceptor_startpoints = Univcoordlist_push(*antiacceptor_startpoints,antiacceptor_splicesite);
	  *antiacceptor_partners = Univcoordlist_push(*antiacceptor_partners,antidonor_splicesite);

	} else {
	  donor_splicesite = post_univdiagonal - this->querylength + splicesite_qpos;
	  acceptor_splicesite = pre_univdiagonal - this->querylength + splicesite_qpos;

	  *donor_startpoints = Univcoordlist_push(*donor_startpoints,donor_splicesite);
	  *donor_partners = Univcoordlist_push(*donor_partners,acceptor_splicesite);
	  *acceptor_startpoints = Univcoordlist_push(*acceptor_startpoints,acceptor_splicesite);
	  *acceptor_partners = Univcoordlist_push(*acceptor_partners,donor_splicesite);
	}
      }
    }
  }

  return;
}


void
Path_learn_indels (Path_T this, Univcoordtable_T indel_table) {
  Intlist_T q;
  Univcoordlist_T p;
  List_T j;
  Junction_T junction;
  Univcoord_T indel_position, univdiagonal;
  int indel_qpos, support;
  Intlist_T indel_adjs;
  int adj;

  q = this->endpoints;
  p = this->univdiagonals;
  for (j = this->junctions; j != NULL; j = List_next(j)) {
    junction = (Junction_T) List_head(j);

    if ((adj = Junction_adj(junction)) != 0) {
      indel_qpos = Intlist_head(Intlist_next(q));

      if (indel_qpos < this->querylength - indel_qpos) {
	support = indel_qpos;
      } else {
	support = this->querylength - indel_qpos;
      }
      
      if (support < MIN_TWOPASS_INDEL_SUPPORT) {
	/* Don't use */
      } else {
	univdiagonal = Univcoordlist_head(p);
	indel_position = univdiagonal - this->querylength + indel_qpos;
	indel_adjs = (Intlist_T) Univcoordtable_get(indel_table,indel_position);
	indel_adjs = Intlist_push(indel_adjs,adj);
	Univcoordtable_put(indel_table,indel_position,indel_adjs);
      }
    }

    q = Intlist_next(q);
    p = Univcoordlist_next(p);
  }

  return;
}


void
Pathpair_learn_insertlengths (Pathpair_T pathpair, Uintlist_T *insertlengths) {
  *insertlengths = Uintlist_push(*insertlengths,pathpair->insertlength);
  return;
}

void
Pathpair_analyze_insertlengths (int *expected_pairlength, int *pairlength_deviation,
				Uintlist_T insertlengths) {
  Chrpos_T *array, median, mad;
  int n, i;

  if ((n = Uintlist_length(insertlengths)) > 0) {
    array = Uintlist_to_array(insertlengths,0);
    median = Orderstat_uint_pct(array,n,/*pct*/0.50);
  
    for (i = 0; i < n; i++) {
      if (array[i] > median) {
	array[i] = array[i] - median;
      } else {
	array[i] = median - array[i];
      }
    }
    
    mad = Orderstat_uint_pct(array,n,/*pct*/0.50);
    fprintf(stderr,"Insert lengths have a median of %u and mad of %u\n",
	    median,mad);
    
    /* Store information for pass 2 */
    *expected_pairlength = (int) median;
    *pairlength_deviation = rint(1.4826 * (double) mad);
    if (*pairlength_deviation < 1) {
      *pairlength_deviation = 1;
    }

    FREE(array);
  }

  return;
}
