static char rcsid[] = "$Id: 264fa2f6d8d3060b6e1d13dd1f840ec1599c8685 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "transcript.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>		/* For qsort */
#include <ctype.h>		/* For isalpha */

#include "assert.h"
#include "mem.h"
#include "sense.h"
#include "exon.h"


/* Transcript_new */
/* Also need to turn on DEBUG1 in exon.c */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

/* Printing and debugging procedures */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

/* Concordance */
#ifdef DEBUG2
#define debug2(x) x
#else
#define debug2(x)
#endif

/* Transcript_list_trim */
#ifdef DEBUG5
#define debug5(x) x
#else
#define debug5(x)
#endif

/* Transcript intersectp */
#ifdef DEBUG11
#define debug11(x) x
#else
#define debug11(x)
#endif


static int pairmax_transcriptome;
static Outputtype_T output_type;
static Transcriptome_T transcriptome;

/* For debugging */
static Univ_IIT_T transcript_iit;


#define T Transcript_T


void
Transcript_free (T *old, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  Exon_list_gc(&(*old)->exons,listpool,transcriptpool);
  Transcriptpool_free_transcript(&(*old),transcriptpool
				 transcriptpool_trace(__FILE__,__LINE__));
  return;
}


void
Transcript_list_gc (List_T *list, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T p;
  T old;

  for (p = *list; p != NULL; p = List_next(p)) {
    old = (T) List_head(p);
    Exon_list_gc(&old->exons,listpool,transcriptpool);
    Transcriptpool_free_transcript(&old,transcriptpool
				   transcriptpool_trace(__FILE__,__LINE__));
  }
  Listpool_free_list(&(*list),listpool
		     listpool_trace(__FILE__,__LINE__));
  return;
}


T
Transcript_new (unsigned int num, int transcript_genestrand, int trstart, int trend,
		int trstart_overhang, int trend_overhang,
		List_T exons, int nexons, int trlength,
		Transcriptpool_T transcriptpool) {

  T new = Transcriptpool_new_transcript(transcriptpool
					transcriptpool_trace(__FILE__,__LINE__));

#ifdef DEBUG0
  static int call_i = 0;
  char *label;
  bool allocp;
#endif

  assert(trstart < trend);

#ifdef DEBUG0
  label = Univ_IIT_label(transcript_iit,num,&allocp);
  printf("%d: %p: Making transcript %s with trnum %d, transcript_genestrand %d, start %d, end %d\n",
	 ++call_i,new,label,num,transcript_genestrand,trstart,trend);
  if (allocp) {
    FREE(label);
  }
#endif
  
  new->num = num;
  new->genestrand = transcript_genestrand;

  new->trstart = trstart;
  new->trend = trend;
  
  new->trstart_overhang = trstart_overhang;
  new->trend_overhang = trend_overhang;

  new->exons = exons;
  new->nexons = nexons;
  new->trlength = trlength;

  new->velocity = -1;

  return new;
}


T
Transcript_copy (T old, Transcriptpool_T transcriptpool, Listpool_T listpool) {
  T new = Transcriptpool_new_transcript(transcriptpool
					transcriptpool_trace(__FILE__,__LINE__));

  debug0(printf("Copying %p from %p\n",new,old));

  new->num = old->num;
  new->genestrand = old->genestrand;
  
  new->trstart = old->trstart;
  new->trend = old->trend;

  new->trstart_overhang = old->trstart_overhang;
  new->trend_overhang = old->trend_overhang;

  new->exons = Exon_list_copy(old->exons,transcriptpool,listpool);
  new->nexons = old->nexons;
  new->trlength = old->trlength;

  new->velocity = old->velocity;

  return new;
}

List_T
Transcript_copy_list (List_T old, Transcriptpool_T transcriptpool, Listpool_T listpool) {
  List_T new = NULL, p;
  T this;

  for (p = old; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    new = Listpool_push(new,listpool,(void *) Transcript_copy(this,transcriptpool,listpool)
			listpool_trace(__FILE__,__LINE__));
  }
  return List_reverse(new);
}


void
Transcript_list_ascendingp (List_T list) {
  T this;
  List_T q;
  unsigned int prev_trnum;

  if (list != NULL) {
    this = (T) List_head(list);
    prev_trnum = this->num;
    for (q = List_next(list); q != NULL; q = List_next(q)) {
      this = (T) List_head(q);
      if (this->num < prev_trnum) {
	printf("Expecting ascending, but got\n");
	Transcript_print_nums(list);
	abort();
      }
      prev_trnum = this->num;
    }
  }

  return;
}


static bool
Transcript_equal (T x, T y) {
  if (x->num != y->num) {
    return false;
  } else if (x->trstart != y->trstart) {
    return false;
  } else if (x->trend != y->trend) {
    return false;
  } else {
    return true;
  }
}


#if 0
bool
Transcript_in_list_p (T x, List_T list) {
  List_T p;
  T y;

  for (p = list; p != NULL; p = List_next(p)) {
    y = (T) List_head(p);
    if (x->num == y->num &&
	x->trstart == y->trstart &&
	x->trend == y->trend) {
      return true;
    }
  }
  return false;
}
#endif


static int
transcript_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->num < y->num) {
    return -1;
  } else if (y->num < x->num) {
    return +1;
  } else {
    return 0;
  }
}


List_T
Transcript_list_sort (List_T transcripts, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T p, q, prev;
  T *array, last_transcript;
  int n, i;

  if ((n = List_length(transcripts)) == 0) {
    return (List_T) NULL;
  } else {
    array = (T *) List_to_array(transcripts,NULL);
    qsort(array,n,sizeof(T),transcript_cmp);

    List_head_set(transcripts,(void *) array[0]);
    last_transcript = array[0];
    prev = transcripts;

    for (q = p = List_next(transcripts), i = 1; p != NULL; p = List_next(p), i++) {
      if (Transcript_equal(array[i],last_transcript) == true) {
	Transcript_free(&(array[i]),listpool,transcriptpool);
      } else {
	List_head_set(q,(void *) array[i]);
	prev = q;
	q = List_next(q);
	last_transcript = array[i];	
      }
    }
    prev->rest = (List_T) NULL;

    /* Transcript_print_nums(transcripts); */
    /* printf("\n"); */

    FREE(array);
    return transcripts;
  }
}


#if 0
/* Previously took trim5 and trim3 as parameters, when we used
   Path_common_structure_trim_p in Path_consolidate */
void
Transcript_list_trim (List_T transcripts, int qstart, int qend, Transcriptome_T transcriptome,
		      Listpool_T listpool, Transcriptpool_T transcriptpool) {
  T this;
  int *exonbounds;
  Chrpos_T *exonstarts;
  List_T exons, p;
  Exon_T exon;
  char firstchar, lastchar;
#ifdef DEBUG5
  int nexons, exoni;
#endif

  debug5(printf("Entered Transcript_list_trim with qstart %d and qend %d\n",qstart,qend));

  for (p = transcripts; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    exons = this->exons;

#ifdef DEBUG5
    nexons = Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,/*trnum*/this->num);
    printf("Trimming transcript %p %d with trstart %d, trend %d and qstart %d, qend %d: ",
	   this,this->num,this->trstart,this->trend,qstart,qend);
    Exon_print_list_stdout(exons);
    printf("\n");
    for (exoni = 0; exoni < nexons; exoni++) {
      printf("Exon %d, exonbound %d, exonstart %u\n",
	     exoni,exonbounds[exoni],exonstarts[exoni]);
    }
#else
    Transcriptome_exons(&exonbounds,&exonstarts,transcriptome,/*trnum*/this->num);
#endif

    if (trim5 > 0) {
      exon = (Exon_T) List_head(exons);
      lastchar = exon->lastchar;
      while (exons != NULL && this->trstart >= exonbounds[((Exon_T) List_head(exons))->exoni]) {
	exon = (Exon_T) List_head(exons);
	lastchar = exon->lastchar;
#ifdef DEBUG5
	printf("Trimming start exon %d because trstart %d >= exonbound %d\n",
	       exon->exoni,this->trstart,exonbounds[exon->exoni]);
#endif
	exons = Listpool_pop(exons,listpool,(void **) &exon
			     listpool_trace(__FILE__,__LINE__));
	Exon_free(&exon,transcriptpool);
      }
      if (lastchar != 'x') {
	this->trstart += trim5;
      }

      /* Fix firstchar */
      exon = (Exon_T) List_head(exons);
      debug5(printf("Comparing trstart %d with exonbound %d\n",
		    this->trstart,exonbounds[exon->exoni - 1]));
      if (exon->exoni > 0 && this->trstart > exonbounds[exon->exoni - 1]) {
	exon->firstchar = '.';
      }
    }

    if (trim3 > 0) {
      exons = List_reverse(exons);
      exon = (Exon_T) List_head(exons);
      firstchar = exon->firstchar;
      while (exons != NULL && ((Exon_T) List_head(exons))->exoni > 0 &&
	     this->trend <= exonbounds[((Exon_T) List_head(exons))->exoni - 1]) {
	exon = (Exon_T) List_head(exons);
	firstchar = exon->firstchar;
#ifdef DEBUG5
	printf("Trimming end exon %d because trend %d <= exonbound %d\n",
	       exon->exoni,this->trend,exonbounds[exon->exoni - 1]);
#endif
	exons = Listpool_pop(exons,listpool,(void **) &exon
			     listpool_trace(__FILE__,__LINE__));
	Exon_free(&exon,transcriptpool);
      }
      if (firstchar != 'x') {
	this->trend -= trim3;
      }

      /* Fix lastchar */
      exon = (Exon_T) List_head(exons);
      debug5(printf("Comparing trend %d with exonbound %d\n",
		    this->trend,exonbounds[exon->exoni]));
      if (this->trend < exonbounds[exon->exoni]) {
	exon->lastchar = '.';
      }

      exons = List_reverse(exons);
    }

#ifdef DEBUG5
    printf("Result of trimming for transcript %p %d [%d..%d]: ",
	   this,this->num,this->trstart,this->trend);
    Exon_print_list_stdout(exons);
    printf("\n");
#endif
    
    assert(exons != NULL);
    this->exons = exons;
  }
    
  return;
}
#endif


void
Transcript_print_nums (List_T list) {
  List_T p;
  T this;

  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %u[%d|%d..%d|%d]:",
	   this->num,this->trstart_overhang,this->trstart,this->trend,this->trend_overhang);
    assert(this->exons != NULL);
    Exon_print_list_stdout(this->exons);
  }
  return;
}


#ifdef DEBUG1
void
Transcript_print_list_debug (List_T list) {
  List_T p;
  T this;

  printf(" Trnums:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %u",this->num);
  }
  printf("\n");

  printf(" Trstarts:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d",this->trstart);
  }
  printf("\n");

  printf(" Trends:");
  for (p = list; p != NULL; p = List_next(p)) {
    this = (T) List_head(p);
    printf(" %d",this->trend);
  }
  printf("\n");

  return;
}
#endif


static void
print_one (Filestring_T fp, T this, Univ_IIT_T transcript_iit) {
  char *label;
  bool allocp;
  int trstart, trend;

  label = Univ_IIT_label(transcript_iit,this->num,&allocp);

  if ((trstart = this->trstart - this->trstart_overhang) < 1) {
    trstart = 1;
  }
  if ((trend = this->trend + this->trend_overhang) > this->trlength) {
    trend = this->trlength;
  }

#ifdef CHECK_ASSERTIONS
  FPRINTF(fp,"%d:",this->num);
#endif
  FPRINTF(fp,"%s:%d..%d:",label,trstart,trend); /* Do not need to add 1 for 1-based */
  switch (this->velocity) {
  case SPLICED: FPRINTF(fp,"S:"); break;
  case UNSPLICED: FPRINTF(fp,"U:"); break;
  case BOTH: FPRINTF(fp,"B:"); break;
  default: FPRINTF(fp,"*:"); break;
  }
  Exon_print_list(fp,this->exons);

  if (allocp) {
    FREE(label);
  }

  return;
}


void
Transcript_print_list (Filestring_T fp, List_T transcripts,
		       Univ_IIT_T transcript_iit, char *header) {
  T transcript;
  List_T p;

  if (transcripts != NULL) {
    transcript = List_head(transcripts);
    FPRINTF(fp,"%s",header);
    print_one(fp,transcript,transcript_iit);

    for (p = List_next(transcripts); p != NULL; p = List_next(p)) {
      transcript = List_head(p);
      PUTC(',',fp);
      print_one(fp,transcript,transcript_iit);
    }
  }

  return;
}


bool
Transcript_intersectp (List_T transcripts5, List_T transcripts3) {
  List_T p, q;
  T transcript5, transcript3;

#ifdef DEBUG11
  printf("Looking for intersection between ");
  Transcript_print_nums(transcripts5);
  printf(" and ");
  Transcript_print_nums(transcripts3);
  printf("\n");
#endif

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts5);
  Transcript_list_ascendingp(transcripts3);
#endif

  p = transcripts5;
  q = transcripts3;

  while (p != NULL && q != NULL) {
    transcript5 = (T) List_head(p);
    transcript3 = (T) List_head(q);
    if (transcript5->num < transcript3->num) {
      p = List_next(p);
    } else if (transcript3->num < transcript5->num) {
      q = List_next(q);
    } else {
      return true;
    }
  }
  
  return false;
}


#if 0
int
Transcript_intersection (List_T *transcripts5, List_T *transcripts3,
			 Path_T path5, Path_T path3,
			 Listpool_T listpool, Transcriptpool_T transcriptpool) {
  int min_fragment_length = 0, fragment_length;
  List_T intersection5 = NULL, intersection3 = NULL, p, q;
  T transcript5, transcript3;

  for (p = *transcripts5; p != NULL; p = List_next(p)) {
    transcript5 = (T) List_head(p);
    for (q = *transcripts3; q != NULL; q = List_next(q)) {
      transcript3 = (T) List_head(q);
      if (transcript5->num == transcript3->num) {
	intersection5 = Listpool_push(intersection5,listpool,(void *) transcript5);
	intersection3 = Listpool_push(intersection3,listpool,(void *) transcript3);
	if (transcript5->trstart < transcript3->trend) {
	  fragment_length = transcript3->trend - transcript5->trstart + 1;
	  if (transcript5->genestrand > 0) {
	    fragment_length += Path_qstart_trim(path5) + Path_qend_trim(path3);
	  } else {
	    fragment_length += Path_qend_trim(path5) + Path_qstart_trim(path3);
	  }
	  /* printf("Have an intersection of %d..%d and %d..%d => fragment length %d\n",
	     transcript5->trstart,transcript5->trend,
	     transcript3->trstart,transcript3->trend,fragment_length); */
	  if (min_fragment_length == 0 || fragment_length < min_fragment_length) {
	    min_fragment_length = fragment_length;
	  }
	} else if (transcript3->trstart < transcript5->trend) {
	  fragment_length = transcript5->trend - transcript3->trstart + 1;
	  if (transcript5->genestrand > 0) {
	    fragment_length += Path_qend_trim(path5) + Path_qstart_trim(path3);
	  } else {
	    fragment_length += Path_qstart_trim(path5) + Path_qend_trim(path3);
	  }
	  /* printf("Have an intersection of %d..%d and %d..%d => fragment length %d\n",
		 transcript5->trstart,transcript5->trend,
		 transcript3->trstart,transcript3->trend,fragment_length); */
	  if (min_fragment_length == 0 || fragment_length < min_fragment_length) {
	    min_fragment_length = fragment_length;
	  }
	}
      }
    }
  }

  if (intersection5 != NULL && intersection3 != NULL) {
    Transcript_list_gc(&(*transcripts5),listpool,transcriptpool);
    Transcript_list_gc(&(*transcripts3),listpool,transcriptpool);
    *transcripts5 = List_reverse(intersection5);
    *transcripts3 = List_reverse(intersection3);
  } else if (intersection5 != NULL) {
    Transcript_list_gc(&intersection5,listpool,transcriptpool);
  } else if (intersection3 != NULL) {
    Transcript_list_gc(&intersection3,listpool,transcriptpool);
  }

  return min_fragment_length;
}
#endif
  

static List_T
compute_intersection5_valid (List_T transcripts5, Path_T path3, Listpool_T listpool) {
  List_T intersection5 = NULL, p, q;
  T transcript5, transcript3;
  unsigned int transcript3_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts5);
  Transcript_list_ascendingp(path3->transcripts);
#endif

  p = transcripts5;
  q = path3->transcripts;

  while (p != NULL && q != NULL) {
    if (q == NULL) {
      transcript3_num = (unsigned int) -1U;
    } else {
      transcript3 = (T) List_head(q);
      transcript3_num = transcript3->num;
    }

    transcript5 = (T) List_head(p);
    if (transcript5->num < transcript3_num) {
      p = List_next(p);
    } else if (transcript3_num < transcript5->num) {
      q = List_next(q);
    } else {
      intersection5 = Listpool_push(intersection5,listpool,(void *) transcript5
				    listpool_trace(__FILE__,__LINE__));
      p = List_next(p);
      q = List_next(q);
    }
  }

  return List_reverse(intersection5);
}

static List_T
compute_intersection5_invalid (List_T transcripts5, Path_T path3, Listpool_T listpool) {
  List_T intersection5 = NULL, p, r;
  T transcript5, invalid_transcript3;
  unsigned int invalid_transcript3_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts5);
  Transcript_list_ascendingp(path3->invalid_transcripts);
#endif

  p = transcripts5;
  r = path3->invalid_transcripts;

  while (p != NULL && r != NULL) {
    if (r == NULL) {
      invalid_transcript3_num = (unsigned int) -1U;
    } else {
      invalid_transcript3 = (T) List_head(r);
      invalid_transcript3_num = invalid_transcript3->num;
    }
    
    transcript5 = (T) List_head(p);
    if (transcript5->num < invalid_transcript3_num) {
      p = List_next(p);
    } else if (invalid_transcript3_num < transcript5->num) {
      r = List_next(r);
    } else {
      intersection5 = Listpool_push(intersection5,listpool,(void *) transcript5
				    listpool_trace(__FILE__,__LINE__));
      p = List_next(p);
      r = List_next(r);
    }
  }

  return List_reverse(intersection5);
}


static List_T
compute_intersection3_valid (List_T transcripts3, Path_T path5, Listpool_T listpool) {
  List_T intersection3 = NULL, p, q;
  T transcript5, transcript3;
  unsigned int transcript5_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts3);
  Transcript_list_ascendingp(path5->transcripts);
#endif

  p = transcripts3;
  q = path5->transcripts;

  while (p != NULL && q != NULL) {
    if (q == NULL) {
      transcript5_num = (unsigned int) -1U;
    } else {
      transcript5 = (T) List_head(q);
      transcript5_num = transcript5->num;
    }

    transcript3 = (T) List_head(p);
    if (transcript3->num < transcript5_num) {
      p = List_next(p);
    } else if (transcript5_num < transcript3->num) {
      q = List_next(q);
    } else {
      intersection3 = Listpool_push(intersection3,listpool,(void *) transcript3
				    listpool_trace(__FILE__,__LINE__));
      p = List_next(p);
      q = List_next(q);
    }
  }

  return List_reverse(intersection3);
}


static List_T
compute_intersection3_invalid (List_T transcripts3, Path_T path5, Listpool_T listpool) {
  List_T intersection3 = NULL, p, r;
  T transcript3, invalid_transcript5;
  unsigned int invalid_transcript5_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts3);
  Transcript_list_ascendingp(path5->invalid_transcripts);
#endif

  p = transcripts3;
  r = path5->invalid_transcripts;

  while (p != NULL && r != NULL) {
    if (r == NULL) {
      invalid_transcript5_num = (unsigned int) -1U;
    } else {
      invalid_transcript5 = (T) List_head(r);
      invalid_transcript5_num = invalid_transcript5->num;
    }
    
    transcript3 = (T) List_head(p);
    if (transcript3->num < invalid_transcript5_num) {
      p = List_next(p);
    } else if (invalid_transcript5_num < transcript3->num) {
      r = List_next(r);
    } else {
      intersection3 = Listpool_push(intersection3,listpool,(void *) transcript3
				    listpool_trace(__FILE__,__LINE__));
      p = List_next(p);
      r = List_next(r);
    }
  }

  return List_reverse(intersection3);
}


static int
compute_fragment_length (int min_fragment_length, List_T transcripts5, List_T transcripts3,
			 Path_T path5, Path_T path3, Shortread_T queryseq5, Shortread_T queryseq3) {
  int fragment_length;
  List_T p, q;
  T transcript5, transcript3;
  int insert_end, insert_start;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts5);
  Transcript_list_ascendingp(transcripts3);
#endif

  p = transcripts5;
  q = transcripts3;

  while (p != NULL && q != NULL) {
    transcript5 = (T) List_head(p);
    transcript3 = (T) List_head(q);
    if (transcript5->num < transcript3->num) {
      p = List_next(p);
    } else if (transcript3->num < transcript5->num) {
      q = List_next(q);
    } else {
      if (Transcript_genestrand(transcript5) > 0) {
	if (path5->plusp == true) {
	  insert_end = transcript3->trstart + path3->querylength;
	  insert_start = transcript5->trend - path5->querylength;

	  insert_end -= Shortread_left_choplength(queryseq3);
	  insert_end += Shortread_right_choplength(queryseq3);
	  insert_start -= Shortread_left_choplength(queryseq5);
	  insert_start += Shortread_right_choplength(queryseq5);

	} else {
	  insert_end = transcript5->trstart + path5->querylength;
	  insert_start = transcript3->trend - path3->querylength;

	  insert_end -= Shortread_right_choplength(queryseq5);
	  insert_end += Shortread_left_choplength(queryseq5);
	  insert_start -= Shortread_right_choplength(queryseq3);
	  insert_start += Shortread_left_choplength(queryseq3);
	}

      } else {
	if (path5->plusp == true) {
	  insert_end = transcript3->trend - path3->querylength;
	  insert_start = transcript5->trstart + path5->querylength;

	  insert_end -= Shortread_left_choplength(queryseq3);
	  insert_end += Shortread_right_choplength(queryseq3);
	  insert_start -= Shortread_left_choplength(queryseq5);
	  insert_start += Shortread_right_choplength(queryseq5);

	} else {
	  insert_end = transcript5->trend - path5->querylength;
	  insert_start = transcript3->trstart + path3->querylength;

	  insert_end -= Shortread_right_choplength(queryseq5);
	  insert_end += Shortread_left_choplength(queryseq5);
	  insert_start -= Shortread_right_choplength(queryseq3);
	  insert_start += Shortread_left_choplength(queryseq3);
	}
      }

      fragment_length = insert_end - insert_start;
      if (min_fragment_length == 0 || fragment_length < min_fragment_length) {
	min_fragment_length = fragment_length;
      }

      p = List_next(p);
      q = List_next(q);
    }
  }
    
  return min_fragment_length;
}


static List_T
free_non_intersection (List_T all_transcripts, List_T intersection,
		       Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T keep = NULL;
  List_T p, q;
  T transcript, elt;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(all_transcripts);
  Transcript_list_ascendingp(intersection);
#endif

  p = all_transcripts;
  q = intersection;
  
  while (p != NULL && q != NULL) {
    transcript = (T) List_head(p);
    elt = (T) List_head(q);

    if (transcript->num < elt->num) {
      Transcript_free(&transcript,listpool,transcriptpool);
      p = List_next(p);

    } else if (elt->num < transcript->num) {
      fprintf(stderr,"Transcript %u is in intersection but not superset\n",
	      elt->num);
      abort();
      
    } else {
      keep = Listpool_push(keep,listpool,(void *) transcript
			   listpool_trace(__FILE__,__LINE__));
      p = List_next(p);
      q = List_next(q);
    }
  }
  
  while (p != NULL) {
    transcript = (T) List_head(p);
    Transcript_free(&transcript,listpool,transcriptpool);
    p = List_next(p);
  }

  Listpool_free_list(&all_transcripts,listpool
		     listpool_trace(__FILE__,__LINE__));

  return List_reverse(keep);
}


/* If an intersection is found, then uses it, and returns true.
   Otherwise, keeps the original transcripts, and returns false */
bool
Transcript_intersection (Path_T path5, Path_T path3,
			 Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T new_transcripts5, new_transcripts3,
    new_invalid_transcripts5, new_invalid_transcripts3;

#ifdef DEBUG11
  printf("Before intersection:\n");
  printf("Valid 5: "); Transcript_print_nums(path5->transcripts); printf("\n");
  printf("Invalid 5: "); Transcript_print_nums(path5->invalid_transcripts); printf("\n");
  printf("Valid 3: "); Transcript_print_nums(path3->transcripts); printf("\n");
  printf("Invalid 3: "); Transcript_print_nums(path3->invalid_transcripts); printf("\n");
  printf("\n");
#endif

  new_transcripts5 = compute_intersection5_valid(path5->transcripts,path3,listpool);
  new_transcripts3 = compute_intersection3_valid(path3->transcripts,path5,listpool);

  if (new_transcripts5 != NULL && new_transcripts3 != NULL) {
    /* Found intersection among valid transcripts, so use them.  No need to report invalid transcripts */
    Transcript_list_gc(&path5->invalid_transcripts,listpool,transcriptpool);
    Transcript_list_gc(&path3->invalid_transcripts,listpool,transcriptpool);
    path5->invalid_transcripts = (List_T) NULL;
    path3->invalid_transcripts = (List_T) NULL;
    /* Same as new_transcripts5 */
    path5->transcripts = free_non_intersection(path5->transcripts,new_transcripts5,
					       listpool,transcriptpool);
    /* Same as new_transcripts3 */
    path3->transcripts = free_non_intersection(path3->transcripts,new_transcripts3,
					       listpool,transcriptpool);
    debug11(printf("(1) Returning true for transcript_concordant_p\n"));
    return true;

  } else {
    /* Try intersection with valid to invalid transcripts */
    /* assert(new_transcripts5 == NULL); */
    /* assert(new_transcripts3 == NULL); */
    Transcript_list_gc(&new_transcripts5,listpool,transcriptpool);
    Transcript_list_gc(&new_transcripts3,listpool,transcriptpool);

    new_transcripts5 = compute_intersection5_invalid(path5->transcripts,path3,listpool);
    new_transcripts3 = compute_intersection3_invalid(path3->transcripts,path5,listpool);

    if (new_transcripts5 != NULL || new_transcripts3 != NULL) {
      /* Found an intersection, so use them */
      path5->invalid_transcripts = free_non_intersection(path5->invalid_transcripts,new_transcripts3,
							 listpool,transcriptpool);
      path3->invalid_transcripts = free_non_intersection(path3->invalid_transcripts,new_transcripts5,
							 listpool,transcriptpool);
      path5->transcripts = free_non_intersection(path5->transcripts,new_transcripts5,
						 listpool,transcriptpool);
      path3->transcripts = free_non_intersection(path3->transcripts,new_transcripts3,
						 listpool,transcriptpool);

      /* transcript_concordant_p means both ends have concordant transcripts that are valid */
      debug11(printf("(2) Returning false for transcript_concordant_p\n"));
      return false;

    } else {
      /* Try intersection with invalid to invalid transcripts */
      new_invalid_transcripts5 = compute_intersection5_invalid(path5->invalid_transcripts,path3,
							       listpool);
      new_invalid_transcripts3 = compute_intersection3_invalid(path3->invalid_transcripts,path5,
							       listpool);

      if (new_invalid_transcripts5 != NULL && new_invalid_transcripts3 != NULL) {
	/* Found intersection of invalid transcripts, so use them.  No need to report valid transcripts */
	Transcript_list_gc(&path5->transcripts,listpool,transcriptpool);
	Transcript_list_gc(&path3->transcripts,listpool,transcriptpool);
	path5->transcripts = (List_T) NULL;
	path3->transcripts = (List_T) NULL;
    
	path5->invalid_transcripts = free_non_intersection(path5->invalid_transcripts,new_invalid_transcripts5,
							   listpool,transcriptpool);
	path3->invalid_transcripts = free_non_intersection(path3->invalid_transcripts,new_invalid_transcripts3,
							   listpool,transcriptpool);
	debug11(printf("(3) Returning false for transcript_concordant_p\n"));
	return false;

      } else {
	/* No intersection found, so keep original transcripts */
	debug11(printf("(4) Returning false for transcript_concordant_p\n"));
	return false;
      }
    }
  }

#if 0
  printf("After intersection: ");
  Transcript_print_nums(path5->transcripts);
  Transcript_print_nums(path5->invalid_transcripts);
  Transcript_print_nums(path3->transcripts);
  Transcript_print_nums(path3->invalid_transcripts);
  printf("\n");
#endif
  
  return false;
}


int
Transcript_fragment_length (Path_T path5, Path_T path3, Shortread_T queryseq5, Shortread_T queryseq3) {
  int min_fragment_length = 0;

  min_fragment_length = compute_fragment_length(min_fragment_length,path5->transcripts,path3->transcripts,
						path5,path3,queryseq5,queryseq3);
  if (min_fragment_length > 0) {
    return min_fragment_length;
  }
  
  min_fragment_length = compute_fragment_length(min_fragment_length,path5->transcripts,path3->invalid_transcripts,
						path5,path3,queryseq5,queryseq3);
  min_fragment_length = compute_fragment_length(min_fragment_length,path5->invalid_transcripts,path3->transcripts,
						path5,path3,queryseq5,queryseq3);
  if (min_fragment_length > 0) {
    return min_fragment_length;
  }

  min_fragment_length = compute_fragment_length(min_fragment_length,path5->invalid_transcripts,path3->invalid_transcripts,
						path5,path3,queryseq5,queryseq3);
  return min_fragment_length;
}


void
Transcript_repair_trstart (T this, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  Exon_T exon;

  this->trstart -= this->trstart_overhang;
  this->trstart_overhang = 0;

  exon = (Exon_T) List_head(this->exons);
  exon->firstchar = 's';

  exon = Exon_new('.',exon->exoni - 1,'s',transcriptpool);
  this->exons = Listpool_push(this->exons,listpool,(void *) exon
			      listpool_trace(__FILE__,__LINE__));

  /* this->trstart -= trstart_overhang; */

  return;
}

void
Transcript_repair_trend (T this, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  Exon_T exon;

  this->trend += this->trend_overhang;
  this->trend_overhang = 0;

  this->exons = List_reverse(this->exons);

  exon = (Exon_T) List_head(this->exons);
  exon->lastchar = 's';

  exon = Exon_new('s',exon->exoni + 1,'.',transcriptpool);
  this->exons = Listpool_push(this->exons,listpool,(void *) exon
			      listpool_trace(__FILE__,__LINE__));

  this->exons = List_reverse(this->exons);

  /* this->trend += trend_overhang; */

  return;
}


List_T
Transcript_remove_subset (List_T transcripts, List_T subset,
			  Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T keep = NULL;
  List_T p, q;
  T transcript, remove;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts);
  Transcript_list_ascendingp(subset);
#endif

  p = transcripts;
  q = subset;
  
  while (q != NULL) {
    transcript = (T) List_head(p);
    remove = (T) List_head(q);

    if (transcript->num < remove->num) {
      keep = Listpool_push(keep,listpool,(void *) transcript
				listpool_trace(__FILE__,__LINE__));
      p = List_next(p);

    } else if (remove->num < transcript->num) {
      fprintf(stderr,"Transcript %u is not in the subset of transcripts\n",
	      remove->num);
      abort();
      
    } else {
      Transcript_free(&transcript,listpool,transcriptpool);
      Transcript_free(&remove,listpool,transcriptpool);
      p = List_next(p);
      q = List_next(q);
    }
  }
  
  while (p != NULL) {
    transcript = (T) List_head(p);
    keep = Listpool_push(keep,listpool,(void *) transcript
			      listpool_trace(__FILE__,__LINE__));
    p = List_next(p);
  }

  return List_reverse(keep);
}


void
Transcript_setup (int pairmax_transcriptome_in, 
		  Outputtype_T output_type_in, Transcriptome_T transcriptome_in,
		  Univ_IIT_T transcript_iit_in) {

  pairmax_transcriptome = pairmax_transcriptome_in;
  output_type = output_type_in;
  transcriptome = transcriptome_in;
  
  /* For debugging */
  transcript_iit = transcript_iit_in;

  return;
}


