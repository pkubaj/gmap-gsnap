static char rcsid[] = "$Id: 3c6cc7085a888103177cf6ee40c84cc31b656e45 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "knownsplicing.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"

#include "iitdef.h"
#include "interval.h"
#include "sedgesort.h"

#include "nr-x.h"
#include "getline.h"

#include "genomebits_count.h"
#include "sense.h"


#define PVALUE_THRESHOLD 1e-6


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* process_intervals.  Creation of EF64 object */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Knownsplicing_start_univdiagonals and Knownsplicing_end_univdiagonals */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif


/* Previous data structure:
   Genomecomp_T splicecomp (bit vector over the genome)
   Univcoord_T *splicesites (over splicesites)
   Splicetypes *splicetypes (over nsplicesites)
   Chrpos_T *splicedists (over nsplicesites)
*/

/* New data structure:
   EF_T donor_ef (a bit vector over the genome using Elias Fano compression)
   Univcoord_T *donor_endpoints (2*nsplicesites, alternating donor position and then acceptor position)

   Same for acceptor, antidonor, and antiacceptor */


#define T Knownsplicing_T
struct T {
  EF64_T donor_ef;
  Univcoord_T *donor_endpoints;
  int donor_nintervals;

  EF64_T acceptor_ef;
  Univcoord_T *acceptor_endpoints;
  int acceptor_nintervals;

  EF64_T antidonor_ef;
  Univcoord_T *antidonor_endpoints;
  int antidonor_nintervals;

  EF64_T antiacceptor_ef;
  Univcoord_T *antiacceptor_endpoints;
  int antiacceptor_nintervals;
};

void
Knownsplicing_free (T *old) {

  FREE((*old)->antiacceptor_endpoints);
  EF64_free(&(*old)->antiacceptor_ef);
  FREE((*old)->antidonor_endpoints);
  EF64_free(&(*old)->antidonor_ef);
  FREE((*old)->acceptor_endpoints);
  EF64_free(&(*old)->acceptor_ef);
  FREE((*old)->donor_endpoints);
  EF64_free(&(*old)->donor_ef);
  FREE(*old);

  return;
}


static int
univcoord_cmp (const void *x, const void *y) {
  Univcoord_T a = * (Univcoord_T *) x;
  Univcoord_T b = * (Univcoord_T *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return +1;
  } else {
    return 0;
  }
}


/* TODO: Currently splices are rejected if they are suboptimal on either
   donor and acceptor sides, but we should get better performance if
   the rejected on both */

/* endpoints are position pairs of splice sites */
static int
process_intervals (Univcoord_T **endpoints, Univcoord_T **rejections, int *nrejections,
		   Univcoordlist_T startpoints, Univcoordlist_T partners,
		   Univcoordtableuint_T start_table, Univcoordtableuint_T partner_table) {
  int nintervals = 0;
  Univcoord_T *array1, *array2, *subarray2;
  Univcoordlist_T reject_list = NULL, reject_partners = NULL;
  double pvalue;
  int length, n, m, i, ii, j, jj, k;
  int max_count, count;
  unsigned int start_count, partner_count;
  int *order;


  debug0(printf("Entered process_intervals with %d startpoints and partners\n",
		Univcoordlist_length(startpoints)));

  if ((length = Univcoordlist_length(startpoints)) == 0) {
    *endpoints = (Univcoord_T *) NULL;
    *rejections = (Univcoord_T *) NULL;
    *nrejections = 0;
    return 0;
      
  } else {
    array1 = Univcoordlist_to_array(startpoints,/*end*/(Univcoord_T) 0); /* Creates n+1 elts */
    array2 = Univcoordlist_to_array(partners,/*end*/(Univcoord_T) 0);
    Univcoordlist_free(&startpoints);
    Univcoordlist_free(&partners);

#ifdef LARGE_GENOMES
    order = Sedgesort_order_uint8(array1,length);	     /* Requires n+1 elts */
#else
    order = Sedgesort_order_uint4(array1,length);	     /* Requires n+1 elts */
#endif

    /* Remove duplicates */
    *endpoints = (Univcoord_T *) MALLOC(2*length*sizeof(Univcoord_T));

    n = 0;
    i = 0;
    while (i < length) {
      ii = i + 1;
      while (ii < length && array1[order[ii]] == array1[order[i]]) {
	ii++;
      }

      /* Have a unique startpoint, and sort the partners */
      m = ii - i;
      subarray2 = (Univcoord_T *) MALLOC((m + 1)*sizeof(Univcoord_T));
      for (k = i, j = 0; k < ii; k++, j++) {
	subarray2[j] = array2[order[k]];
      }
      qsort(subarray2,m,sizeof(Univcoord_T),univcoord_cmp);


      max_count = 0;		/* Max for this startpoint */
      j = 0;
      while (j < m) {
	jj = j + 1;
	while (jj < m && subarray2[jj] == subarray2[j]) {
	  jj++;
	}
	if (jj - j > max_count) {
	  max_count = jj - j;
	}
	j = jj;
      }
	
      /* Check all partners */
      j = 0;
      while (j < m) {
	jj = j + 1;
	while (jj < m && subarray2[jj] == subarray2[j]) {
	  jj++;
	}

	if ((count = jj - j) == 1) {
	  reject_list = Univcoordlist_push(reject_list,subarray2[j]);
	  reject_partners = Univcoordlist_push(reject_partners,array1[order[i]]);

	} else if (count == max_count) {
	  /* Accept */
	  (*endpoints)[n++] = (Univcoord_T) array1[order[i]];
	  (*endpoints)[n++] = (Univcoord_T) subarray2[j];
	    
	} else if ((pvalue = NR_ppois((double) count,/*lambda*/(double) max_count)) >= PVALUE_THRESHOLD) {
	  /* Accept as a relatively high count */
	  (*endpoints)[n++] = (Univcoord_T) array1[order[i]];
	  (*endpoints)[n++] = (Univcoord_T) subarray2[j];

	} else if (start_table != NULL && partner_table != NULL &&
		   (start_count = Univcoordtableuint_get(start_table,array1[order[i]])) > 0 &&
		   (partner_count = Univcoordtableuint_get(partner_table,subarray2[j])) > 0) {

	  /* Acct as a suboptimal intron supported by other reads */
	  fprintf(stderr,"Accepting %llu..%llu with count %d relative to max %d because of splicesite counts %u and %u\n",
		  (unsigned long long) array1[order[i]],(unsigned long long) subarray2[j],
		  count,max_count,start_count,partner_count);

	  (*endpoints)[n++] = (Univcoord_T) array1[order[i]];
	  (*endpoints)[n++] = (Univcoord_T) subarray2[j];

	} else {
	  debug(printf("Rejecting %llu..%llu with count of 1\n",
		       (unsigned long long) array1[order[i]],(unsigned long long) subarray2[j]));
	  reject_list = Univcoordlist_push(reject_list,subarray2[j]);
	  reject_partners = Univcoordlist_push(reject_partners,array1[order[i]]);
	}

	j = jj;
      }
      FREE(subarray2);

      i = ii;
    }

    FREE(array2);
    FREE(array1);
    FREE(order);
    nintervals = n/2;


    /* Process rejections to be used in store_intervals */
    if ((length = Univcoordlist_length(reject_list)) == 0) {
      *rejections = (Univcoord_T *) NULL;
      *nrejections = 0;

    } else {
      array1 = Univcoordlist_to_array(reject_list,/*end*/(Univcoord_T) 0); /* Creates n+1 elts */
      array2 = Univcoordlist_to_array(reject_partners,/*end*/(Univcoord_T) 0);
      Univcoordlist_free(&reject_list);
      Univcoordlist_free(&reject_partners);

#ifdef LARGE_GENOMES
      order = Sedgesort_order_uint8(array1,length);	     /* Requires n+1 elts */
#else
      order = Sedgesort_order_uint4(array1,length);	     /* Requires n+1 elts */
#endif

      *rejections = (Univcoord_T *) MALLOC(2*length*sizeof(Univcoord_T));
      n = 0;
      i = 0;
      while (i < length) {
	ii = i + 1;
	while (ii < length && array1[order[ii]] == array1[order[i]]) {
	  ii++;
	}

	/* Sort the partners */
	m = ii - i;
	subarray2 = (Univcoord_T *) MALLOC((m + 1)*sizeof(Univcoord_T));
	for (k = i, j = 0; k < ii; k++, j++) {
	  subarray2[j] = array2[order[k]];
	}
	qsort(subarray2,m,sizeof(Univcoord_T),univcoord_cmp);

	for (j = 0; j < m; j++) {
	  /* Write in interleaved format */
	  (*rejections)[n++] = (Univcoord_T) array1[order[i]];
	  (*rejections)[n++] = (Univcoord_T) subarray2[j];
	}
	FREE(subarray2);

	i = ii;
      }
      
      FREE(array2);
      FREE(array1);
      FREE(order);
      *nrejections = n/2;
    }

    return nintervals;
  }
}


static EF64_T
store_intervals (Univcoord_T **endpoints, int *nintervals,
		 Univcoord_T *candidates, int ncandidates,
		 Univcoord_T *rejections, int nrejections,
		 Univcoord_T genomelength, FILE *dump_splices_fp,
		 Univ_IIT_T chromosome_iit, EF64_T chromosome_ef64) {
  Univcoord_T startpoint, partner;
  Chrnum_T chrnum;
  char *chr = NULL;
  Univcoord_T chroffset, chrhigh = 0;
  int i, j, k;
  bool allocp;

  if (candidates == NULL) {
    *endpoints = (Univcoord_T *) NULL;
    *nintervals = 0;

  } else {
    k = 0;
    i = j = 0;
    while (i < 2*ncandidates) {
      startpoint = candidates[i];
      partner = candidates[i+1];
      
      /* Advance to the corresponding rejection */
      while (j < 2*nrejections && (rejections[j] < startpoint ||
				   (rejections[j] == startpoint && rejections[j+1] < partner))) {
	j += 2;
      }
    
      if (j < 2*nrejections && rejections[j] == startpoint && rejections[j+1] == partner) {
	debug(printf("Rejecting %u..%u\n",startpoint,partner));
      } else {
	candidates[k++] = startpoint;
	candidates[k++] = partner;
	debug(printf("Writing %u..%u\n",startpoint,partner));
	if (dump_splices_fp != NULL) {
	  if (startpoint > partner) {
	    fprintf(dump_splices_fp,"%llu %llu %u",
		    (unsigned long long) startpoint,(unsigned long long) partner,
		    /*splice distance*/(Chrpos_T) (startpoint - partner));
	  } else {
	    fprintf(dump_splices_fp,"%llu %llu %u",
		    (unsigned long long) startpoint,(unsigned long long) partner,
		    /*splice distance*/(Chrpos_T) (partner - startpoint));
	  }

	  if (startpoint > chrhigh) {
	    /* Update chromosome */
	    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
				 /*alignstart*/startpoint,/*alignend*/startpoint);
	    if (chr != NULL && allocp == true) {
	      FREE(chr);
	    }
	    chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	  }
	  fprintf(dump_splices_fp," %s:%u..%u\n",
		  chr,(Chrpos_T) (startpoint - chroffset),(Chrpos_T) (partner - chroffset));
	}
      }
      
      i += 2;
    }

    if ((*nintervals = k/2) == 0) {
      *endpoints = (Univcoord_T *) NULL;
    } else {
      /* Reduce memory allocated, which was an overestimate */
      *endpoints = (Univcoord_T *) MALLOC(2*(*nintervals)*sizeof(Univcoord_T));
      memcpy(*endpoints,candidates,2*(*nintervals)*sizeof(Univcoord_T));

#ifdef DEBUG0
      printf("%d: %u %u\n",0,(*endpoints)[2*0],(*endpoints)[2*0+1]);
      for (i = 1; i < *nintervals; i++) {
	printf("%d: %u %u\n",i,(*endpoints)[2*i],(*endpoints)[2*i+1]);
	assert((*endpoints)[2*i] >= (*endpoints)[2*(i-1)]);
	if ((*endpoints)[2*i] == (*endpoints)[2*(i-1)]) {
	  assert((*endpoints)[2*i+1] > (*endpoints)[2*(i-1)+1]);
	}
      }
      printf("\n");
#endif
    }

    FREE(candidates);
    if (rejections != NULL) {
      FREE(rejections);
    }
  }

  if (chr != NULL && allocp == true) {
    FREE(chr);
  }

  return EF64_new_from_interleaved_univcoords(*endpoints,*nintervals,genomelength);
}



/* Frees all uintlists */
T
Knownsplicing_new (Univcoordlist_T donor_startpoints, Univcoordlist_T donor_partners,
		   Univcoordlist_T acceptor_startpoints, Univcoordlist_T acceptor_partners,
		   Univcoordlist_T antidonor_startpoints, Univcoordlist_T antidonor_partners,
		   Univcoordlist_T antiacceptor_startpoints, Univcoordlist_T antiacceptor_partners,
		   Univcoordtableuint_T donor_table, Univcoordtableuint_T acceptor_table, 
		   Univcoordtableuint_T antidonor_table, Univcoordtableuint_T antiacceptor_table, 
		   Univcoord_T genomelength, FILE *dump_splices_fp, Univ_IIT_T chromosome_iit,
		   EF64_T chromosome_ef64, bool intron_level_p) {
  T new;
  int from_donor_nintervals, from_acceptor_nintervals,
    from_antidonor_nintervals, from_antiacceptor_nintervals;
  int from_donor_nrejections, from_acceptor_nrejections,
    from_antidonor_nrejections, from_antiacceptor_nrejections;
  Univcoord_T *from_donor_endpoints, *from_acceptor_endpoints,
    *from_antidonor_endpoints, *from_antiacceptor_endpoints;
  Univcoord_T *from_donor_rejections, *from_acceptor_rejections,
    *from_antidonor_rejections, *from_antiacceptor_rejections;

  new = (T) MALLOC(sizeof(*new));
  debug(printf("Plus genestrand\n"));
  /* We have a bipartite graph */
  from_donor_nintervals = process_intervals(&from_donor_endpoints,
					    &from_acceptor_rejections,&from_acceptor_nrejections,
					    donor_startpoints,donor_partners,
					    donor_table,/*partner_table*/acceptor_table);
  from_acceptor_nintervals = process_intervals(&from_acceptor_endpoints,
					       &from_donor_rejections,&from_donor_nrejections,
					       acceptor_startpoints,acceptor_partners,
					       acceptor_table,/*partner_table*/donor_table);
  
  if (dump_splices_fp != NULL) {
    fprintf(dump_splices_fp,">Donors\n");
  }
  
  new->donor_ef = store_intervals(&new->donor_endpoints,&new->donor_nintervals,
				  from_donor_endpoints,from_donor_nintervals,
				  from_donor_rejections,from_donor_nrejections,
				  genomelength,dump_splices_fp,chromosome_iit,chromosome_ef64);
  if (dump_splices_fp != NULL) {
    fprintf(dump_splices_fp,">Acceptors\n");
  }
  new->acceptor_ef = store_intervals(&new->acceptor_endpoints,&new->acceptor_nintervals,
				     from_acceptor_endpoints,from_acceptor_nintervals,
				     from_acceptor_rejections,from_acceptor_nrejections,
				     genomelength,dump_splices_fp,chromosome_iit,chromosome_ef64);
  debug(printf("\n"));

  debug(printf("Minus genestrand\n"));
  from_antidonor_nintervals = process_intervals(&from_antidonor_endpoints,
						&from_antiacceptor_rejections,&from_antiacceptor_nrejections,
						antidonor_startpoints,antidonor_partners,
						antidonor_table,/*partner_table*/antiacceptor_table);
  from_antiacceptor_nintervals = process_intervals(&from_antiacceptor_endpoints,
						   &from_antidonor_rejections,&from_antidonor_nrejections,
						   antiacceptor_startpoints,antiacceptor_partners,
						   antiacceptor_table,/*partner_table*/antidonor_table);
  if (dump_splices_fp != NULL) {
    fprintf(dump_splices_fp,">Antidonors\n");
  }
  new->antidonor_ef = store_intervals(&new->antidonor_endpoints,&new->antidonor_nintervals,
				      from_antidonor_endpoints,from_antidonor_nintervals,
				      from_antidonor_rejections,from_antidonor_nrejections,
				      genomelength,dump_splices_fp,chromosome_iit,chromosome_ef64);
  if (dump_splices_fp != NULL) {
    fprintf(dump_splices_fp,">Antiacceptors\n");
  }
  new->antiacceptor_ef = store_intervals(&new->antiacceptor_endpoints,&new->antiacceptor_nintervals,
					 from_antiacceptor_endpoints,from_antiacceptor_nintervals,
					 from_antiacceptor_rejections,from_antiacceptor_nrejections,
					 genomelength,dump_splices_fp,chromosome_iit,chromosome_ef64);
  debug(printf("\n"));

  if (intron_level_p == true && new->donor_nintervals != new->acceptor_nintervals) {
    fprintf(stderr,"Knownsplicing_new: donor nintervals %d != acceptor nintervals %d.  Please report to twu@gene.com\n",
	    new->donor_nintervals,new->acceptor_nintervals);
    exit(9);
  } else if (intron_level_p == true && new->antidonor_nintervals != new->antiacceptor_nintervals) {
    fprintf(stderr,"Knownsplicing_new: antidonor nintervals %d != antiacceptor nintervals %d.  Please report to twu@gene.com\n",
	    new->antidonor_nintervals,new->antiacceptor_nintervals);
    exit(9);
  } else {
    fprintf(stderr,"Observed %d distinct introns on genome plus strand and %d distinct introns on genome minus strand\n",
	    new->donor_nintervals,new->antidonor_nintervals);
  }

  return new;
}


T
Knownsplicing_new_from_dump (FILE *fp, Univcoord_T genomelength) {
  T new = (T) MALLOC(sizeof(*new));
  Univcoordlist_T endpoint_list;
  Univcoord_T startpoint, partner;
  char *line;

  /* Donors */
  endpoint_list = (Univcoordlist_T) NULL;
  if ((line = Getline(fp)) == NULL || line[0] != '>') {
    fprintf(stderr,"Knownsplicing expecting file to start with >Donors.  Got %s\n",line);
    exit(9);
  } else {
    FREE(line);
  }
  while ((line = Getline(fp)) != NULL && line[0] != '>') {
    if (
#ifdef LARGE_GENOMES
	sscanf(line,"%llu %llu",&startpoint,&partner)
#else
	sscanf(line,"%u %u",&startpoint,&partner)
#endif
	< 2) {
      fprintf(stderr,"Knownsplicing cannot parse line %s\n",line);
      exit(9);
    } else {
      endpoint_list = Univcoordlist_push(endpoint_list,startpoint);
      endpoint_list = Univcoordlist_push(endpoint_list,partner);
    }
    FREE(line);
  }

  endpoint_list = Univcoordlist_reverse(endpoint_list);
  new->donor_nintervals = Univcoordlist_length(endpoint_list)/2;
  new->donor_endpoints = Univcoordlist_to_array(endpoint_list,/*end*/0);
  Univcoordlist_free(&endpoint_list);
  new->donor_ef = EF64_new_from_interleaved_univcoords(new->donor_endpoints,new->donor_nintervals,genomelength);
  

  /* Acceptors */
  endpoint_list = (Univcoordlist_T) NULL;
  if (line == NULL || line[0] != '>') {
    fprintf(stderr,"Knownsplicing expecting file to have >Acceptors.  Got %s\n",line);
    exit(9);
  } else {
    FREE(line);
  }
  while ((line = Getline(fp)) != NULL && line[0] != '>') {
    if (
#ifdef LARGE_GENOMES
	sscanf(line,"%llu %llu",&startpoint,&partner)
#else
	sscanf(line,"%u %u",&startpoint,&partner)
#endif
	< 2) {
      fprintf(stderr,"Knownsplicing cannot parse line %s\n",line);
      exit(9);
    } else {
      endpoint_list = Univcoordlist_push(endpoint_list,startpoint);
      endpoint_list = Univcoordlist_push(endpoint_list,partner);
    }
    FREE(line);
  }

  endpoint_list = Univcoordlist_reverse(endpoint_list);
  new->acceptor_nintervals = Univcoordlist_length(endpoint_list)/2;
  new->acceptor_endpoints = Univcoordlist_to_array(endpoint_list,/*end*/0);
  Univcoordlist_free(&endpoint_list);
  new->acceptor_ef = EF64_new_from_interleaved_univcoords(new->acceptor_endpoints,new->acceptor_nintervals,genomelength);


  /* Antidonors */
  endpoint_list = (Univcoordlist_T) NULL;
  if (line == NULL || line[0] != '>') {
    fprintf(stderr,"Knownsplicing expecting file to have >Antidonors.  Got %s\n",line);
    exit(9);
  } else {
    FREE(line);
  }
  while ((line = Getline(fp)) != NULL && line[0] != '>') {
    if (
#ifdef LARGE_GENOMES
	sscanf(line,"%llu %llu",&startpoint,&partner)
#else
	sscanf(line,"%u %u",&startpoint,&partner)
#endif
	< 2) {
      fprintf(stderr,"Knownsplicing cannot parse line %s\n",line);
      exit(9);
    } else {
      endpoint_list = Univcoordlist_push(endpoint_list,startpoint);
      endpoint_list = Univcoordlist_push(endpoint_list,partner);
    }
    FREE(line);
  }

  endpoint_list = Univcoordlist_reverse(endpoint_list);
  new->antidonor_nintervals = Univcoordlist_length(endpoint_list)/2;
  new->antidonor_endpoints = Univcoordlist_to_array(endpoint_list,/*end*/0);
  Univcoordlist_free(&endpoint_list);
  new->antidonor_ef = EF64_new_from_interleaved_univcoords(new->antidonor_endpoints,new->antidonor_nintervals,genomelength);


  /* Antiacceptors */
  endpoint_list = (Univcoordlist_T) NULL;
  if (line == NULL || line[0] != '>') {
    fprintf(stderr,"Knownsplicing expecting file to have >Antiacceptors.  Got %s\n",line);
    exit(9);
  } else {
    FREE(line);
  }
  while ((line = Getline(fp)) != NULL) {
    if (
#ifdef LARGE_GENOMES
	sscanf(line,"%llu %llu",&startpoint,&partner)
#else
	sscanf(line,"%u %u",&startpoint,&partner)
#endif
	< 2) {
      fprintf(stderr,"Knownsplicing cannot parse line %s\n",line);
      exit(9);
    } else {
      endpoint_list = Univcoordlist_push(endpoint_list,startpoint);
      endpoint_list = Univcoordlist_push(endpoint_list,partner);
    }
    FREE(line);
  }

  endpoint_list = Univcoordlist_reverse(endpoint_list);
  new->antiacceptor_nintervals = Univcoordlist_length(endpoint_list)/2;
  new->antiacceptor_endpoints = Univcoordlist_to_array(endpoint_list,/*end*/0);
  Univcoordlist_free(&endpoint_list);
  new->antiacceptor_ef = EF64_new_from_interleaved_univcoords(new->antiacceptor_endpoints,new->antiacceptor_nintervals,genomelength);

  return new;
}


int
Knownsplicing_nintervals (T this) {
  return this->donor_nintervals + this->acceptor_nintervals +
    this->antidonor_nintervals + this->antiacceptor_nintervals;
}
    

Univcoord_T *
Knownsplicing_donors (uint64_t *low_rank, uint64_t *high_rank, T this,
		      Univcoord_T univdiagonal, int querylength, int pos5, int pos3) {
  Univcoord_T left;

  if (this->donor_ef == NULL) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;

  } else if (pos5 >= pos3) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;
    
  } else {
    /* *low_rank = EF64_rank(this->donor_ef,left + pos5); */
    /* *high_rank = EF64_rank(this->donor_ef,left + pos3 - 1); */
    /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
    assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
    left = univdiagonal - querylength;
    EF64_two_ranks(&(*low_rank),&(*high_rank),this->donor_ef,
		   /*low*/left + pos5,/*high*/left + pos3 - 1);
    return this->donor_endpoints;
  }
}

Univcoord_T *
Knownsplicing_acceptors (uint64_t *low_rank, uint64_t *high_rank, T this,
			 Univcoord_T univdiagonal, int querylength, int pos5, int pos3) {
  Univcoord_T left;

  if (this->acceptor_ef == NULL) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;

  } else if (pos5 >= pos3) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;
    
  } else {
    /* *low_rank = EF64_rank(this->acceptor_ef,left + pos5); */
    /* *high_rank = EF64_rank(this->acceptor_ef,left + pos3 - 1); */
    /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
    assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
    left = univdiagonal - querylength;
    EF64_two_ranks(&(*low_rank),&(*high_rank),this->acceptor_ef,
		   /*low*/left + pos5,/*high*/left + pos3 - 1);
    return this->acceptor_endpoints;
  }
}

Univcoord_T *
Knownsplicing_antidonors (uint64_t *low_rank, uint64_t *high_rank, T this,
			  Univcoord_T univdiagonal, int querylength, int pos5, int pos3) {
  Univcoord_T left;

  if (this->antidonor_ef == NULL) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;

  } else if (pos5 >= pos3) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;
    
  } else {
    /* *low_rank = EF64_rank(this->antidonor_ef,left + pos5); */
    /* *high_rank = EF64_rank(this->antidonor_ef,left + pos3 - 1); */
    /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
    assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
    left = univdiagonal - querylength;
    EF64_two_ranks(&(*low_rank),&(*high_rank),this->antidonor_ef,
		   /*low*/left + pos5,/*high*/left + pos3 - 1);
    return this->antidonor_endpoints;
  }
}

Univcoord_T *
Knownsplicing_antiacceptors (uint64_t *low_rank, uint64_t *high_rank, T this,
			     Univcoord_T univdiagonal, int querylength, int pos5, int pos3) {
  Univcoord_T left;

  if (this->antiacceptor_ef == NULL) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;

  } else if (pos5 >= pos3) {
    *low_rank = *high_rank = 0;
    return (Univcoord_T *) NULL;
    
  } else {
    /* *low_rank = EF64_rank(this->antiacceptor_ef,left + pos5); */
    /* *high_rank = EF64_rank(this->antiacceptor_ef,left + pos3 - 1); */
    /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
    assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
    left = univdiagonal - querylength;
    EF64_two_ranks(&(*low_rank),&(*high_rank),this->antiacceptor_ef,
		   /*low*/left + pos5,/*high*/left + pos3 - 1);
    return this->antiacceptor_endpoints;
  }
}


#if 0
/* Interface designed to match that of Regiondb_get_univdiagonals */
Univcoord_T *
Knownsplicing_qstart_univdiagonals (int *nentries, T this,
				    Univcoord_T univdiagonal, Compress_T query_compress,
				    int pos5, int pos3, int querylength,
				    Genomebits_T genomebits, Genomebits_T genomebits_alt,
				    bool plusp, int try_sensedir) {
  Univcoord_T *result;
  Univcoord_T *endpoints, splice_dist, distal_univdiagonal;
  uint64_t low_rank, high_rank, rank;

  int splice_pos;
  int best_nmismatches, distal_nmismatches, distal_ref_nmismatches;
  Univcoordlist_T best_univdiagonals = NULL;

  debug1(printf("Entered Knownsplicing_qstart_univdiagonals with univdiagonal %u, pos5 %d, pos3 %d, plusp %d, try_sensedir %d\n",
		univdiagonal,pos5,pos3,plusp,try_sensedir));

  assert(pos5 < pos3);
  if (plusp == true) {
    if (try_sensedir == SENSE_FORWARD) {
      endpoints = Knownsplicing_acceptors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing acceptors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else if (try_sensedir == SENSE_ANTI) {
      endpoints = Knownsplicing_antidonors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing antidonors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else {
      fprintf(stderr,"Unexpected value for try_sensedir\n");
      abort();
    }

  } else {
    if (try_sensedir == SENSE_FORWARD) {
      endpoints = Knownsplicing_antidonors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing antidonors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else if (try_sensedir == SENSE_ANTI) {
      endpoints = Knownsplicing_acceptors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing acceptors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else {
      fprintf(stderr,"Unexpected value for try_sensedir\n");
      abort();
    }
  }

  best_nmismatches = querylength;
  for (rank = low_rank; rank < high_rank; rank++) {
    debug1(printf("Splice %u..%u relative to left %u\n",endpoints[2*rank],endpoints[2*rank+1],left));
    assert(endpoints[2*rank] >= left + pos5);
    assert(endpoints[2*rank] < left + pos3);
    
    splice_pos = endpoints[2*rank] - left;
    splice_dist = endpoints[2*rank+1] - endpoints[2*rank];
    distal_univdiagonal = univdiagonal + splice_dist;
    distal_nmismatches = Genomebits_count_mismatches_substring(&distal_ref_nmismatches,
							       genomebits,genomebits_alt,query_compress,
							       /*univdiagonal*/distal_univdiagonal,querylength,
							       pos5,/*pos3*/splice_pos,
							       plusp,/*genestrand*/0);
    debug1(printf("Rank #%lu at qpos %d, %u..%u => %d nmismatches\n",
		  rank,splice_pos,endpoints[2*rank],endpoints[2*rank+1],distal_nmismatches));

    if (distal_nmismatches < best_nmismatches) {
      Univcoordlist_free(&best_univdiagonals);
      best_univdiagonals = Univcoordlist_push(NULL,/*univdiagonal*/distal_univdiagonal);
      best_nmismatches = distal_nmismatches;
    } else if (distal_nmismatches == best_nmismatches) {
      best_univdiagonals = Univcoordlist_push(best_univdiagonals,/*univdiagonal*/distal_univdiagonal);
    }
  }

  if (best_univdiagonals == (Univcoordlist_T) NULL) {
    *nentries = 0;
    return (Univcoord_T *) NULL;
  } else {
    *nentries = Univcoordlist_length(best_univdiagonals);
    result = Univcoordlist_to_array(best_univdiagonals,/*end*/0U);
    Univcoordlist_free(&best_univdiagonals);
    return result;
  }
}
#endif


#if 0
/* Interface designed to match that of Regiondb_get_univdiagonals */
Univcoord_T *
Knownsplicing_qend_univdiagonals (int *nentries, T this,
				  Univcoord_T univdiagonal, Compress_T query_compress,
				  int pos5, int pos3, int querylength,
				  Genomebits_T genomebits, Genomebits_T genomebits_alt,
				  bool plusp, int try_sensedir) {
  Univcoord_T *result;
  Univcoord_T *endpoints, splice_dist, distal_univdiagonal;
  uint64_t low_rank, high_rank, rank;

  int splice_pos;
  int best_nmismatches, distal_nmismatches, distal_ref_nmismatches;
  Univcoordlist_T best_univdiagonals = NULL;

  debug1(printf("Entered Knownsplicing_qend_univdiagonals with left %u, pos5 %d, pos3 %d, plusp %d, try_sensedir %d\n",
		left,pos5,pos3,plusp,try_sensedir));

  assert(pos5 < pos3);
  if (plusp == true) {
    if (try_sensedir == SENSE_FORWARD) {
      endpoints = Knownsplicing_donors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing donors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else if (try_sensedir == SENSE_ANTI) {
      endpoints = Knownsplicing_antiacceptors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing antiacceptors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else {
      fprintf(stderr,"Unexpected value for try_sensedir\n");
      abort();
    }

  } else {
    /* Not sure how to handle pos5 and pos3 */
    if (try_sensedir == SENSE_FORWARD) {
      endpoints = Knownsplicing_antiacceptors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing antiacceptors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else if (try_sensedir == SENSE_ANTI) {
      endpoints = Knownsplicing_donors(&low_rank,&high_rank,this,univdiagonal,querylength,pos5,pos3);
      debug1(printf("Knownsplicing donors at %u + %d..%d yields low_rank %lu to high_rank %lu\n",
		    left,pos5,pos3,low_rank,high_rank));
    } else {
      fprintf(stderr,"Unexpected value for try_sensedir\n");
      abort();
    }
  }

  best_nmismatches = querylength;
  for (rank = low_rank; rank < high_rank; rank++) {
    splice_pos = endpoints[2*rank] - left;
    splice_dist = endpoints[2*rank+1] - endpoints[2*rank];
    distal_univdiagonal = univdiagonal + splice_dist;
    distal_nmismatches = Genomebits_count_mismatches_substring(&distal_ref_nmismatches,
							       genomebits,genomebits_alt,query_compress,
							       /*univdiagonal*/distal_univdiagonal,querylength,
							       /*pos5*/splice_pos,pos3,
							       plusp,/*genestrand*/0);
    debug1(printf("Rank #%lu at qpos %d, %u..%u => %d nmismatches\n",
		  rank,splice_pos,endpoints[2*rank],endpoints[2*rank+1],distal_nmismatches));

    if (distal_nmismatches < best_nmismatches) {
      Univcoordlist_free(&best_univdiagonals);
      best_univdiagonals = Univcoordlist_push(NULL,/*univdiagonal*/distal_univdiagonal);
      best_nmismatches = distal_nmismatches;
    } else if (distal_nmismatches == best_nmismatches) {
      best_univdiagonals = Univcoordlist_push(best_univdiagonals,/*univdiagonal*/distal_univdiagonal);
    }
  }

  if (best_univdiagonals == (Univcoordlist_T) NULL) {
    *nentries = 0;
    return (Univcoord_T *) NULL;
  } else {
    *nentries = Univcoordlist_length(best_univdiagonals);
    result = Univcoordlist_to_array(best_univdiagonals,/*end*/0U);
    Univcoordlist_free(&best_univdiagonals);
    return result;
  }
}
#endif



T
Knownsplicing_from_splicing_iit (IIT_T splicing_iit, int *splicing_divint_crosstable,
				 int donor_typeint, int acceptor_typeint, Univ_IIT_T chromosome_iit,
				 bool intron_level_p) {
  Univcoord_T genomelength, chroffset, chrhigh, low_position, high_position;
  char *chr;
  bool allocp;
  Chrpos_T chrlength;
  struct Interval_T *intervals;
  Interval_T interval;
  int divno, nintervals, i;
  Chrnum_T chrnum;
  Univcoordlist_T donor_startpoints = NULL, donor_partners = NULL, acceptor_startpoints = NULL, acceptor_partners = NULL,
    antidonor_startpoints = NULL, antidonor_partners = NULL, antiacceptor_startpoints = NULL, antiacceptor_partners = NULL;
  /* char donor1, donor2, acceptor1, acceptor2; */

  for (chrnum = 1; chrnum <= Univ_IIT_total_nintervals(chromosome_iit); chrnum++) {
    if ((divno = splicing_divint_crosstable[chrnum]) > 0) {
      Univ_IIT_interval_bounds(&chroffset,&chrhigh,&chrlength,chromosome_iit,chrnum,/*circular_typeint*/-1);
      intervals = splicing_iit->intervals[divno];
      nintervals = splicing_iit->nintervals[divno];

      for (i = 0; i < nintervals; i++) {
	interval = &(intervals[i]);
	low_position = chroffset + Interval_low(interval);
	high_position = chroffset + Interval_high(interval) - 1; /* To get to 0-based */

	if (low_position >= chrhigh || high_position >= chrhigh) {
	  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	  fprintf(stderr,"Splice site %s:%u..%u extends beyond chromosome length %u.  Discarding\n",
		  chr,(Chrpos_T) (low_position - chroffset),(Chrpos_T) (high_position - chroffset),chrlength);
	  if (allocp) FREE(chr);

	} else if (intron_level_p == true) {
	  /* Treat all intervals as donor..acceptor */
	  if (Interval_sign(interval) > 0) {
#if 0
	    donor1 = Genomebits_get_char(ref,low_position+1);
	    donor2 = Genomebits_get_char(ref,low_position+2);
	    acceptor2 = Genomebits_get_char(ref,high_position-2);
	    acceptor1 = Genomebits_get_char(ref,high_position-1);
#endif
	    donor_startpoints = Univcoordlist_push(donor_startpoints,low_position);
	    donor_partners = Univcoordlist_push(donor_partners,high_position);
	    acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,high_position);
	    acceptor_partners = Univcoordlist_push(acceptor_partners,low_position);
	  } else {
#if 0
	    acceptor1 = Genomebits_get_char(ref,low_position+1);
	    acceptor2 = Genomebits_get_char(ref,low_position+2);
	    donor2 = Genomebits_get_char(ref,high_position-2);
	    donor1 = Genomebits_get_char(ref,high_position-1);
#endif
	    antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,high_position);
	    antidonor_partners = Univcoordlist_push(antidonor_partners,low_position);
	    antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,low_position);
	    antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,high_position);
	  }

	} else if (Interval_type(interval) == donor_typeint) {
	  if (Interval_sign(interval) > 0) {
	    donor_startpoints = Univcoordlist_push(donor_startpoints,low_position);
	    donor_partners = Univcoordlist_push(donor_partners,0);
	  } else {
	    antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,low_position);
	    antidonor_partners = Univcoordlist_push(antidonor_partners,0);
	  }
	    
	} else if (Interval_type(interval) == acceptor_typeint) {
	  if (Interval_sign(interval) > 0) {
	    acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,low_position);
	    acceptor_partners = Univcoordlist_push(acceptor_partners,0);
	  } else {
	    antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,low_position);
	    antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,0);
	  }
	}
      }
    }
  }

  genomelength = Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/true);
  return Knownsplicing_new(donor_startpoints,donor_partners,
			   acceptor_startpoints,acceptor_partners,
			   antidonor_startpoints,antidonor_partners,
			   antiacceptor_startpoints,antiacceptor_partners,
			   /*donor_table*/NULL,/*acceptor_table*/NULL,
			   /*antidonor_table*/NULL,/*antiacceptor_table*/NULL,
			   genomelength,/*dump_splices_fp*/NULL,/*chromosome_iit*/NULL,
			   /*chromosome_ef64*/NULL,intron_level_p);
}


/* Modified from transcript-remap.c */
static Chrpos_T *
compute_exonends_geneplus (int *exonbounds, Chrpos_T *exonstarts, int nexons) {
  Chrpos_T *exonends;
  int last_bound, exonlength;
  int exoni;

  exonends = (Chrpos_T *) MALLOC(nexons*sizeof(Chrpos_T));

  last_bound = 0;
  for (exoni = 0; exoni < nexons; exoni++) {
    exonlength = exonbounds[exoni] - last_bound;
    exonends[exoni] = exonstarts[exoni] + exonlength - 1;
    last_bound = exonbounds[exoni];
  }

  return exonends;
}


/* Modified from transcript-remap.c */
static Chrpos_T *
compute_exonends_geneminus (int *exonbounds, Chrpos_T *exonstarts, int nexons) {
  Chrpos_T *exonends;
  int last_bound, exonlength;
  int exoni;

  exonends = (Chrpos_T *) MALLOC(nexons*sizeof(Chrpos_T));

  last_bound = 0;
  for (exoni = 0; exoni < nexons; exoni++) {
    exonlength = exonbounds[exoni] - last_bound;
    exonends[exoni] = exonstarts[exoni] - exonlength + 1;
    last_bound = exonbounds[exoni];
  }

  return exonends;
}


T
Knownsplicing_from_transcriptome (Transcriptome_T transcriptome, int nalignments,
				  EF64_T chromosome_ef64, Univcoord_T genomelength, bool intron_level_p) {
  Trnum_T trnum;
  int map_index;

  int transcript_genestrand;
  int nexons, exoni;
  Chrpos_T *exonstarts, *exonends;
  int *exonbounds;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh, donor_position, acceptor_position;
  Univcoordlist_T donor_startpoints = NULL, donor_partners = NULL, acceptor_startpoints = NULL, acceptor_partners = NULL,
    antidonor_startpoints = NULL, antidonor_partners = NULL, antiacceptor_startpoints = NULL, antiacceptor_partners = NULL;

  for (map_index = 1; map_index <= nalignments; map_index++) {
    if ((trnum = Transcriptome_trnum(&nexons,&exonbounds,&exonstarts,transcriptome,map_index)) < 1) {
      /* Skip.  Not in transcriptome */
    } else {
      chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
      EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,chrnum);

      if (transcript_genestrand > 0) {
	exonends = compute_exonends_geneplus(exonbounds,exonstarts,nexons);

	for (exoni = 1; exoni < nexons; exoni++) {
	  donor_position = chroffset + (Univcoord_T) exonends[exoni - 1] /*+ 1 (after exon) - 1*/;
	  acceptor_position = chroffset + (Univcoord_T) exonstarts[exoni] - 1; /* change 1-based to 0-based */

	  if (intron_level_p == true) {
	    donor_startpoints = Univcoordlist_push(donor_startpoints,donor_position);
	    donor_partners = Univcoordlist_push(donor_partners,acceptor_position);
	    acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,acceptor_position);
	    acceptor_partners = Univcoordlist_push(acceptor_partners,donor_position);
	  } else {
	    donor_startpoints = Univcoordlist_push(donor_startpoints,donor_position);
	    donor_partners = Univcoordlist_push(donor_partners,0);
	    acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,acceptor_position);
	    acceptor_partners = Univcoordlist_push(acceptor_partners,0);
	  }
	}

	FREE(exonends);

      } else {
	exonends = compute_exonends_geneminus(exonbounds,exonstarts,nexons);

	for (exoni = 1; exoni < nexons; exoni++) {
	  donor_position = chroffset + (Univcoord_T) exonends[exoni - 1] - 1; /* change 1-based to 0-based */
	  acceptor_position = chroffset + (Univcoord_T) exonstarts[exoni] /*+ 1 (after exon) - 1*/;

	  if (intron_level_p == true) {
	    antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,donor_position);
	    antidonor_partners = Univcoordlist_push(antidonor_partners,acceptor_position);
	    antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,acceptor_position);
	    antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,donor_position);
	  } else {
	    antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,donor_position);
	    antidonor_partners = Univcoordlist_push(antidonor_partners,0);
	    antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,acceptor_position);
	    antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,0);
	  }
	}

	FREE(exonends);
      }
    }
  }

  return Knownsplicing_new(donor_startpoints,donor_partners,
			   acceptor_startpoints,acceptor_partners,
			   antidonor_startpoints,antidonor_partners,
			   antiacceptor_startpoints,antiacceptor_partners,
			   /*donor_table*/NULL,/*acceptor_table*/NULL,
			   /*antidonor_table*/NULL,/*antiacceptor_table*/NULL,
			   genomelength,/*dump_splices_fp*/NULL,/*chromosome_iit*/NULL,
			   /*chromosome_ef64*/NULL,intron_level_p);
}


