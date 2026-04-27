static char rcsid[] = "$Id: knownsplicing.c 225459 2022-12-15 01:18:40Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "knownindels.h"

#include <stdlib.h>
#include <string.h>		/* For memcpy */
#include <math.h>		/* For qsort */

#include "assert.h"
#include "mem.h"

#include "getline.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Knownindels_T
struct T {
  EF64_T ef;
  Univcoord_T *positions;
  char *adjs;
};

void
Knownindels_free (T *old) {

  FREE((*old)->adjs);
  FREE((*old)->positions);
  EF64_free(&(*old)->ef);
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

static int
int_cmp (const void *x, const void *y) {
  int a = * (int *) x;
  int b = * (int *) y;

  if (a < b) {
    return -1;
  } else if (a > b) {
    return +1;
  } else {
    return 0;
  }
}


/* endpoints are position pairs of splice sites */
static Univcoord_T *
process_indels (char **adjs, int **counts, int *npositions, Univcoordtable_T indel_table) {
  Univcoord_T *positions;
  Univcoord_T *array, indel_position;
  Intlist_T indel_adjs;
  int *subarray;
  int length, n, nitems, i, j, jj;
  int max_count;
  int best_adj, adj;
  bool tiep;

  if ((length = Univcoordtable_length(indel_table)) == 0) {
    *adjs = (char *) NULL;
    *counts = (int *) NULL;
    *npositions = 0;
    return (Univcoord_T *) NULL;
      
  } else {
    array = Univcoordtable_keys(indel_table,/*end*/(Univcoord_T) 0); /* Creates n+1 elts */
    qsort(array,length,sizeof(Univcoord_T),univcoord_cmp);

    positions = (Univcoord_T *) MALLOC(length*sizeof(Univcoord_T));
    *adjs = (char *) MALLOC(length*sizeof(char));
    *counts = (int *) MALLOC(length*sizeof(int));

    n = 0;
    for (i = 0; i < length; i++) {
      indel_position = array[i];
      indel_adjs = (Intlist_T) Univcoordtable_get(indel_table,indel_position);
      subarray = Intlist_to_array(&nitems,indel_adjs);
      Intlist_free(&indel_adjs);

      qsort(subarray,nitems,sizeof(int),int_cmp);

      debug(printf("%u:",indel_position));
      max_count = 0;
      best_adj = 0;

      j = 0;
      while (j < nitems) {
	adj = subarray[j];

	jj = j + 1;
	while (jj < nitems && subarray[jj] == adj) {
	  jj++;
	}

	if (adj < -127 || adj > 127) {
	  fprintf(stderr,"At position %llu, saw an indel of length %d.  Ignoring.\n",
		  (unsigned long long) indel_position,adj);
	} else {
	  debug(printf(" %d:%d",adj,jj-j));

	  if (jj - j > max_count) {
	    max_count = jj - j;
	    best_adj = subarray[j];
	    tiep = false;
	  } else if (jj - j == max_count) {
	    tiep = true;
	  }

	  j = jj;
	}
      }
      debug(printf("\n"));
	
      if (max_count == 1) {
	/* Ignore single counts, which are probably due to sequencing error and not the genome */

      } else if (tiep == true) {
	/* Do not store info */
	fprintf(stderr,"At position %llu, ambiguous indels seen\n",
		(unsigned long long) indel_position);

      } else {
	positions[n] = indel_position;
	(*adjs)[n] = (char) best_adj;
	(*counts)[n] = max_count;
	n++;
      }

      FREE(subarray);
    }

    FREE(array);

    *npositions = n;
    return positions;
  }
}


#define abs(x) ((x) < 0) ? -(x) : +(x)

static Univcoordlist_T
filter_indels (char **filtered_adjs, int *n_indel_positions,
	       Univcoord_T *raw_positions, char *raw_adjs, int *raw_counts,
	       int raw_npositions, FILE *dump_indels_fp,
	       Univ_IIT_T chromosome_iit, EF64_T chromosome_ef64) {
  Univcoordlist_T filtered_positions = NULL;
  bool *rejectp;
  int starti, endi, i, k;
  Univcoord_T indel_position;
  int nvalid;
  char nindels;
  int max_count;
  
  Chrnum_T chrnum;
  char *chr = NULL;
  Univcoord_T chroffset, chrhigh = 0;
  bool allocp;

  rejectp = (bool *) CALLOC(raw_npositions,sizeof(bool));

  for (i = 0; i < raw_npositions; i++) {
    indel_position = raw_positions[i];
    nindels = abs(raw_adjs[i]);

    starti = i - 1;
    while (starti >= 0 &&
	   raw_positions[starti] + (Univcoord_T) (abs(raw_adjs[starti])) >= indel_position) {
      starti--;
    }

    endi = i + 1;
    while (endi < raw_npositions &&
	   raw_positions[endi] <= indel_position + (Univcoord_T) nindels) {
      endi++;
    }

    max_count = raw_counts[i];
    for (k = starti + 1; k < i; k++) {
      if (raw_counts[k] < max_count) {
	rejectp[k] = true;
      }
    }
    for (k = i + 1; k < endi; k++) {
      if (raw_counts[k] < max_count) {
	rejectp[k] = true;
      }
    }
  }

  nvalid = 0;
  for (i = 0; i < raw_npositions; i++) {
    if (rejectp[i] == true) {
      fprintf(stderr,"Filtering overlapping indel position %llu with adj %d and count %d\n",
	      (unsigned long long) raw_positions[i],(int) raw_adjs[i],raw_counts[i]);
    } else {
      nvalid++;
    }
  }
  *n_indel_positions = nvalid;


  *filtered_adjs = (char *) MALLOC(nvalid*sizeof(char));
  nvalid = 0;
  for (i = 0; i < raw_npositions; i++) {
    if (rejectp[i] == false) {
      indel_position = raw_positions[i];
      filtered_positions = Univcoordlist_push(filtered_positions,indel_position);
      (*filtered_adjs)[nvalid++] = raw_adjs[i];
      
      if (dump_indels_fp != NULL) {
	fprintf(dump_indels_fp,"%llu %d %d",
		(unsigned long long) indel_position,raw_adjs[i],raw_counts[i]);
	if (indel_position > chrhigh) {
	  /* Update chromosome */
	  chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			       /*alignstart*/indel_position,/*alignend*/indel_position);
	  if (chr != NULL && allocp == true) {
	    FREE(chr);
	  }
	  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
	}
	fprintf(dump_indels_fp," %s:%u\n",chr,(Chrpos_T) (indel_position - chroffset));
      }
    }
  }

  if (chr != NULL && allocp == true) {
    FREE(chr);
  }

  FREE(rejectp);

  return Univcoordlist_reverse(filtered_positions);
}


T
Knownindels_new (Univcoordtable_T indel_table, Univcoord_T genomelength,
		 FILE *dump_indels_fp, Univ_IIT_T chromosome_iit, EF64_T chromosome_ef64) {
  T new = (T) MALLOC(sizeof(*new));
  Univcoordlist_T positions_list;
  int n_indel_positions;

  Univcoord_T *raw_positions;
  char *raw_adjs;
  int *raw_counts;
  int raw_npositions;
  

  raw_positions = process_indels(&raw_adjs,&raw_counts,&raw_npositions,indel_table);
  if (raw_npositions == 0) {
    fprintf(stderr,"Observed 0 distinct indel positions\n");
    return (T) NULL;

  } else {
    positions_list = filter_indels(&new->adjs,&n_indel_positions,
				   raw_positions,raw_adjs,raw_counts,raw_npositions,
				   dump_indels_fp,chromosome_iit,chromosome_ef64);
    FREE(raw_counts);
    FREE(raw_adjs);
    FREE(raw_positions);

    new->positions = Univcoordlist_to_array(positions_list,/*end*/0);
    new->ef = EF64_new_from_univcoordlist(positions_list,genomelength);
    Univcoordlist_free(&positions_list);

    fprintf(stderr,"Observed %d distinct indel positions\n",n_indel_positions);

    return new;
  }
}


T
Knownindels_new_from_dump (FILE *fp, Univcoord_T genomelength) {
  T new = (T) MALLOC(sizeof(*new));
  Univcoordlist_T positions_list = NULL;
  Intlist_T adjs = NULL, p;
  Univcoord_T indel_position;
  int adj;
  int n_indel_positions, k = 0;
  char *line;

  while ((line = Getline(fp)) != NULL) {
    if (
#ifdef LARGE_GENOMES
	sscanf(line,"%llu %d",&indel_position,&adj)
#else
	sscanf(line,"%u %d",&indel_position,&adj)
#endif
	< 2) {
      fprintf(stderr,"Knownindels cannot parse line %s\n",line);
      exit(9);
    } else {
      positions_list = Univcoordlist_push(positions_list,indel_position);
      adjs = Intlist_push(adjs,adj);
    }
    FREE(line);
  }

  positions_list = Univcoordlist_reverse(positions_list);
  adjs = Intlist_reverse(adjs);
  
  n_indel_positions = Univcoordlist_length(positions_list);
  
  new->positions = Univcoordlist_to_array(positions_list,/*end*/0);
  new->ef = EF64_new_from_univcoordlist(positions_list,genomelength);
  Univcoordlist_free(&positions_list);
  
  new->adjs = (char *) MALLOC(n_indel_positions*sizeof(char));
  for (p = adjs; p != NULL; p = Intlist_next(p)) {
    new->adjs[k++] = (char) Intlist_head(p);
  }
  Intlist_free(&adjs);

  return new;
}


int
Knownindels_find_lowest (int *indel_pos, T this,
			 Univcoord_T univdiagonal, int querylength,
			 int pos5, int pos3) {
  Univcoord_T left;
  uint64_t low_rank, high_rank;

  /* *low_rank = EF64_rank(this->ef,left + pos5); */
  /* *high_rank = EF64_rank(this->ef,left + pos3 - 1); */
  /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
  assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
  left = univdiagonal - querylength;
  EF64_two_ranks(&low_rank,&high_rank,this->ef,
		 /*low*/left + pos5,/*high*/left + pos3 - 1);
  if (high_rank == low_rank) {
    return 0;
  } else {
    *indel_pos = this->positions[low_rank] - left;

    /* Do not need to modify indel_pos when adding to qend, because insertion changes qend, but not qstart */

    return (int) this->adjs[low_rank];
  }
}


int
Knownindels_find_highest (int *indel_pos, T this,
			  Univcoord_T univdiagonal, int querylength,
			  int pos5, int pos3) {
  Univcoord_T left;
  int adj;
  uint64_t low_rank, high_rank;

  /* *low_rank = EF64_rank(this->ef,left + pos5); */
  /* *high_rank = EF64_rank(this->ef,left + pos3 - 1); */
  /* Subtract 1 from pos3 because EF64_rank counts a match at the given position */
  assert(univdiagonal + pos5 >= (Univcoord_T) querylength);
  left = univdiagonal - querylength;
  EF64_two_ranks(&low_rank,&high_rank,this->ef,
		 /*low*/left + pos5,/*high*/left + pos3 - 1);
  if (high_rank == low_rank) {
    return 0;
  } else {
    *indel_pos = this->positions[high_rank - 1] - left;

    /* Need to modify indel_pos when adding to qstart, because insertion changes qend, but not qstart */
    if ((adj = (int) this->adjs[high_rank - 1]) < 0) {
      *indel_pos += adj;
      if (*indel_pos <= 0) {
	return 0;
      }
    }

    return adj;
  }
}
