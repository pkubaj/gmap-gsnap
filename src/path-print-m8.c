static char rcsid[] = "$Id: d3702aaaafcb6421b65b07e02907dd9d86072aae $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#include "path-print-m8.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>		/* For tolower */
#include <math.h>

#include "bool.h"
#include "assert.h"
#include "univcoord.h"
#include "junction.h"


/* making substrings */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif


static Univ_IIT_T chromosome_iit;
static double genomesize;	/* For BLAST E-value */


/* Taken from NCBI Blast 2.2.29, algo/blast/core/blast_stat.c */
/* Karlin-Altschul formula: m n exp(-lambda * S + log k) = k m n exp(-lambda * S) */
/* Also in pair.c */

static double
blast_evalue (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);

  return k * (double) alignlength * genomesize * exp(-lambda * score);
}

static double
blast_bitscore (int alignlength, int nmismatches) {
  double k = 0.1;
  double lambda = 1.58;		/* For a +1, -1 scoring scheme */
  double score;
  
  score = (double) ((alignlength - nmismatches) /* scored as +1 */ - nmismatches /* scored as -1 */);
  return (score * lambda - log(k)) / log(2.0);
}


typedef struct Substring_T *Substring_T;
struct Substring_T {
  Univcoord_T univdiagonal;
  int querystart;
  int queryend;
  int nmismatches;
  Chrnum_T chrnum;
  Univcoord_T chroffset;
  bool plusp;
};


static void
Substring_free (Substring_T *old) {
  FREE(*old);
  return;
}


static Substring_T
Substring_new (Univcoord_T univdiagonal, int querystart, int queryend,
	       int nmismatches, Chrnum_T chrnum, Univcoord_T chroffset, bool plusp) {
  Substring_T new = (Substring_T) MALLOC(sizeof(*new));

  new->univdiagonal = univdiagonal;
  new->querystart = querystart;
  new->queryend = queryend;
  new->nmismatches = nmismatches;
  new->chrnum = chrnum;
  new->chroffset = chroffset;
  new->plusp = plusp;

  return new;
}


/* Modified from Substring_print_m8 */
static void
print_substring (Filestring_T fp, Substring_T substring, char *accession,
		 char *acc_suffix, int querylength, bool invertp) {
  char *chr;
  double identity;
  int qstart, qend;
  int alignlength;
  bool allocp;


  chr = Univ_IIT_label(chromosome_iit,substring->chrnum,&allocp);

  FPRINTF(fp,"%s%s",accession,acc_suffix); /* field 0: accession */

  FPRINTF(fp,"\t%s",chr);	/* field 1: chr */

  /* field 2: identity */
  alignlength = substring->queryend - substring->querystart;
  identity = (double) (alignlength - substring->nmismatches)/(double) alignlength;
  FPRINTF(fp,"\t%.1f",100.0*identity);


  FPRINTF(fp,"\t%d",alignlength); /* field 3: query length */

  FPRINTF(fp,"\t%d",substring->nmismatches); /* field 4: nmismatches */

  FPRINTF(fp,"\t0");		/* field 5: gap openings */

  FPRINTF(fp,"\t%d",substring->querystart + 1); /* field 6: query start */

  FPRINTF(fp,"\t%d",substring->queryend); /* field 7: query end */

  /* fields 8 and 9: chr start and end */
  /* Modified from print_coordinates in path-print-alignment.c */
  if (substring->plusp != invertp) {
    qstart = substring->querystart;
    qend = substring->queryend;
    FPRINTF(fp,"\t%u\t%u",substring->univdiagonal - querylength + qstart - substring->chroffset + 1,
	    substring->univdiagonal - querylength + qend - substring->chroffset);
  } else {
    qstart = querylength - substring->queryend;
    qend = querylength - substring->querystart;
    FPRINTF(fp,"\t%u\t%u",substring->univdiagonal - querylength + qend - substring->chroffset,
	    substring->univdiagonal - querylength + qstart - substring->chroffset + 1);
  }

  /* field 10: E value */
  FPRINTF(fp,"\t%.2g",blast_evalue(alignlength,substring->nmismatches));

 /* field 11: bit score */
  FPRINTF(fp,"\t%.1f",blast_bitscore(alignlength,substring->nmismatches));
  
  FPRINTF(fp,"\n");

  return;
}


/* Modified from Path_print_alignment */
void
Path_print_m8 (Filestring_T fp, Path_T path, char *accession, char *acc_suffix,
	       bool invertp, Listpool_T listpool) {
  List_T substrings = NULL, fusion_substrings = NULL;
  List_T j;
  Intlist_T q;			/* for endpoints */
  Intlist_T m;			/* for nmismatches */
  Univcoordlist_T p;
  Univcoord_T univdiagonal;
  Junction_T fusion_junction;
  int querystart, queryend;
  Substring_T substring;
  List_T a;


  debug0(Path_print(path));

  /* Process along watson strand */

  if (path->fusion_querystart_junction != NULL) {
    fusion_junction = path->fusion_querystart_junction;
  } else if (path->fusion_queryend_junction != NULL) {
    fusion_junction = path->fusion_queryend_junction;
  } else {
    fusion_junction = (Junction_T) NULL;
  }


  /* 1a.  Make fusion_substrings.  Process along watson strand, but convert
     qstart and qend to querystart and queryend along query
     direction */
  if (fusion_junction != NULL) {
    /* pre_junction = (Junction_T) NULL; */
    for (p = path->fusion_univdiagonals, q = path->fusion_endpoints,
	   m = path->fusion_nmismatches, j = path->fusion_junctions;
	 p != NULL; p = Univcoordlist_next(p), q = Intlist_next(q), m = Intlist_next(m), j = List_next(j)) {
      univdiagonal = Univcoordlist_head(p);
      /* post_junction = (j == NULL) ? (Junction_T) NULL : (Junction_T) List_head(j); */
      if (path->fusion_plusp != invertp) {
	querystart = Intlist_head(q);
	queryend = Intlist_second_value(q);
      } else {
	querystart = path->querylength - Intlist_second_value(q);
	queryend = path->querylength - Intlist_head(q);
      }
      debug0(printf("Creating pre-fusion substring at %u, %d:%u, %d..%d, plusp:%d\n",
		    univdiagonal,path->fusion_chrnum,univdiagonal - path->fusion_chroffset,
		    querystart,queryend,path->fusion_plusp));
      fusion_substrings = Listpool_push(fusion_substrings,listpool,
					(void *) Substring_new(univdiagonal,querystart,queryend,
							       /*nmismatches*/Intlist_head(m),path->fusion_chrnum,
							       path->fusion_chroffset,path->fusion_plusp)
					listpool_trace(__FILE__,__LINE__));
      /* pre_junction = post_junction; */
    }
			   
    /* 1b.  Order fusion substrings */
    if (path->fusion_plusp != invertp) {
      debug0(printf("Putting fusion substrings back into original order\n"));
      fusion_substrings = List_reverse(fusion_substrings);

    } else {
      debug0(printf("Keeping fusion substrings reversed\n"));
#if 0
      for (a = fusion_substrings; a != NULL; a = List_next(a)) {
	substring = (Substring_T) List_head(a);
	/* temp_junction = substring->pre_junction; */
	/* substring->pre_junction = substring->post_junction; */
	/* substring->post_junction = temp_junction; */
      }
#endif
    }
  }


  /* 2a.  Make main substrings.  Process along watson strand, but convert
     qstart and qend to querystart and queryend along query
     direction */
  /* pre_junction = (Junction_T) NULL; */
  for (p = path->univdiagonals, q = path->endpoints, m = path->nmismatches, j = path->junctions;
       p != NULL; p = Univcoordlist_next(p), q = Intlist_next(q), m = Intlist_next(m), j = List_next(j)) {
    univdiagonal = Univcoordlist_head(p);
    /* post_junction = (j == NULL) ? (Junction_T) NULL : (Junction_T) List_head(j); */
    if (path->plusp != invertp) {
      querystart = Intlist_head(q);
      queryend = Intlist_second_value(q);
    } else {
      querystart = path->querylength - Intlist_second_value(q);
      queryend = path->querylength - Intlist_head(q);
    }
    debug0(printf("Creating main substring at %u, %d:%u, %d..%d, plusp:%d\n",
		  univdiagonal,path->chrnum,univdiagonal - path->chroffset,
		  querystart,queryend,path->plusp));
    substrings = Listpool_push(substrings,listpool,
			       (void *) Substring_new(univdiagonal,querystart,queryend,
						      /*nmismatches*/Intlist_head(m),path->chrnum,
						      path->chroffset,path->plusp)
			       listpool_trace(__FILE__,__LINE__));
    /* pre_junction = post_junction; */
  }

  /* 2b.  Order main substrings */
  if (path->plusp != invertp) {
    debug0(printf("Putting substrings back into original order\n"));
    substrings = List_reverse(substrings);

  } else {
    debug0(printf("Keeping substrings reversed\n"));
#if 0
    for (a = substrings; a != NULL; a = List_next(a)) {
      substring = (Substring_T) List_head(a);
      /* temp_junction = substring->pre_junction; */
      /* substring->pre_junction = substring->post_junction; */
      /* substring->post_junction = temp_junction; */
    }
#endif
  }


  /* 3. Merge fusion substrings and main substrings */
  if (fusion_substrings == NULL) {
    /* Skip */
  } else if (invertp == false) {
    if (path->fusion_querystart_junction != NULL) {
      /* Put fusion at querystart */
      debug0(printf("Putting fusion at querystart\n"));
      substring = List_last_value(fusion_substrings,NULL);
      /* substring->post_junction = fusion_junction; */
      substring = List_head(substrings);
      /* substring->pre_junction = fusion_junction; */
      substrings = List_append(fusion_substrings,substrings);
      
    } else {
      /* Put fusion at queryend */
      debug0(printf("Putting fusion at queryend\n"));
      substring = List_last_value(substrings,NULL);
      /* substring->post_junction = fusion_junction; */
      substring = List_head(fusion_substrings);
      /* substring->pre_junction = fusion_junction; */
      substrings = List_append(substrings,fusion_substrings);
    }
  } else {
    if (path->fusion_querystart_junction != NULL) {
      /* Put fusion at queryend */
      debug0(printf("Putting fusion at queryend\n"));
      substring = List_last_value(substrings,NULL);
      /* substring->post_junction = fusion_junction; */
      substring = List_head(fusion_substrings);
      /* substring->pre_junction = fusion_junction; */
      substrings = List_append(substrings,fusion_substrings);
    } else {
      /* Put fusion at querystart */
      debug0(printf("Putting fusion at querystart\n"));
      substring = List_last_value(fusion_substrings,NULL);
      /* substring->post_junction = fusion_junction; */
      substring = List_head(substrings);
      /* substring->pre_junction = fusion_junction; */
      substrings = List_append(fusion_substrings,substrings);
    }
  }

  /* 4.  Print the substrings */
  for (a = substrings; a != NULL; a = List_next(a)) {
    substring = (Substring_T) List_head(a);
    print_substring(fp,substring,accession,acc_suffix,path->querylength,invertp);
    Substring_free(&substring);
  }

  /* fusion_substrings are appended to substrings */
  Listpool_free_list(&substrings,listpool
		     listpool_trace(__FILE__,__LINE__)); /* Allocated by Listpool_push */

  return;
}



void
Path_print_m8_setup (Univ_IIT_T chromosome_iit_in) {

  chromosome_iit = chromosome_iit_in;
  genomesize = (double) Univ_IIT_genomelength(chromosome_iit,/*with_circular_alias_p*/false);

  return;
}
