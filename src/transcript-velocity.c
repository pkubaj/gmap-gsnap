static char rcsid[] = "$Id: 8ccc0e678ff3165c79ecdcb98e213a66a145f196 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "transcript-velocity.h"

#include "assert.h"
#include "list.h"
#include "exon.h"
#include "transcript.h"

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Transcript_T


/* For equivalence class, we need to consider both valid and invalid
   transcripts.  A read might be consistent with both the exon of one
   transcript and the intron of another, and both need to be
   considered in the equivalence class. */

/* Indicators of splicing: (1) Gene has only a single exon, (2)
   Multiple exons seen (including between paired-ends) */

/* Indicators of retained intron: (1) A singleton iNi, (2) xNs at the
   start or sNx at the end */

/* Retained false, Spliced false: Velocity B
   Retained false, Spliced true: Velocity S
   Retained true, Spliced false: Velocity U
   Retained true, Spliced true: Velocity U */


static Velocity_T
velocity_single (Transcript_T transcript) {
  Exon_T first_exon, last_exon, exon;
  bool splicedp, retainedp = false;

  if (transcript->nexons == 1) {
    splicedp = true;
    retainedp = false;

  } else if (List_length(transcript->exons) == 1) {
    splicedp = false;
    exon = (Exon_T) List_head(transcript->exons);
    if (exon->firstchar == 'i' && exon->lastchar == 'i') {
      retainedp = true;
    } else if (exon->firstchar == 'x' || exon->lastchar == 'x') {
      retainedp = true;
    }

  } else {
    splicedp = true;		/* From multiple exons */
    first_exon = (Exon_T) List_head(transcript->exons);
    last_exon = (Exon_T) List_last_value(transcript->exons,NULL);
  
    if (first_exon->firstchar == 'x') {
      retainedp = true;
    } else if (last_exon->lastchar == 'x') {
      retainedp = true;
    }
  }

  if (retainedp == true) {
    /* Must be pre-mRNA */
    return UNSPLICED;
  } else if (splicedp == true) {
    /* Must be processed mRNA */
    return SPLICED;
  } else {
    /* Could be either */
    return BOTH;
  }
}

void
Transcript_velocity_single (Path_T path) {
  List_T p;
  Transcript_T transcript;

  for (p = path->transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  for (p = path->invalid_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }


  /* Handle fusion transcripts as single */
  for (p = path->fusion_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  for (p = path->fusion_invalid_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  return;
}


static Velocity_T
velocity_paired (Transcript_T transcript5, Transcript_T transcript3) {
  Exon_T first_exon, last_exon, exon5, exon3;
  bool splicedp, retainedp = false;

  assert(transcript5->num == transcript3->num);
  debug(printf("Evaluating paired velocity of trnum %u with %d exons\n",transcript5->num,transcript5->nexons));

  if (transcript5->nexons == 1) {
    splicedp = true;
    retainedp = false;

  } else if (List_length(transcript5->exons) == 1 && List_length(transcript3->exons) == 1) {
    debug(printf("Have one exon on each paired end\n"));
    exon5 = (Exon_T) List_head(transcript5->exons);
    exon3 = (Exon_T) List_head(transcript3->exons);
    if (exon5->exoni == exon3->exoni) {
      splicedp = false;
    } else {
      splicedp = true;
    }

    if (exon5->firstchar == 'i' && exon5->lastchar == 'i') {
      retainedp = true;
    } else if (exon3->firstchar == 'i' && exon3->lastchar == 'i') {
      retainedp = true;
    } else if (exon5->firstchar == 'x' || exon5->lastchar == 'x') {
      retainedp = true;
    } else if (exon3->firstchar == 'x' || exon3->lastchar == 'x') {
      retainedp = true;
    }

  } else if (List_length(transcript5->exons) == 1) {
    splicedp = true;		/* From transcript3 */

    exon5 = (Exon_T) List_head(transcript5->exons);
    if (exon5->firstchar == 'i' && exon5->lastchar == 'i') {
      retainedp = true;
    } else if (exon5->firstchar == 'x' || exon5->lastchar == 'x') {
      retainedp = true;
    } else {
      first_exon = (Exon_T) List_head(transcript3->exons);
      last_exon = (Exon_T) List_last_value(transcript3->exons,NULL);
      
      if (first_exon->firstchar == 'x') {
	retainedp = true;
      } else if (last_exon->lastchar == 'x') {
	retainedp = true;
      }
    }

  } else if (List_length(transcript3->exons) == 1) {
    splicedp = true;		/* From transcript5 */

    exon3 = (Exon_T) List_head(transcript3->exons);
    if (exon3->firstchar == 'i' && exon3->lastchar == 'i') {
      retainedp = true;
    } else if (exon3->firstchar == 'x' || exon3->lastchar == 'x') {
      retainedp = true;
    } else {
      first_exon = (Exon_T) List_head(transcript5->exons);
      last_exon = (Exon_T) List_last_value(transcript5->exons,NULL);

      if (first_exon->firstchar == 'x') {
	retainedp = true;
      } else if (last_exon->lastchar == 'x') {
	retainedp = true;
      }
    }

  } else {
    splicedp = true;		/* From transcript5 and transcript3 */

    first_exon = (Exon_T) List_head(transcript5->exons);
    last_exon = (Exon_T) List_last_value(transcript5->exons,NULL);
    if (first_exon->firstchar == 'x') {
      retainedp = true;
    } else if (last_exon->lastchar == 'x') {
      retainedp = true;
    }

    first_exon = (Exon_T) List_head(transcript3->exons);
    last_exon = (Exon_T) List_last_value(transcript3->exons,NULL);
    if (first_exon->firstchar == 'x') {
      retainedp = true;
    } else if (last_exon->lastchar == 'x') {
      retainedp = true;
    }
  }
    
  if (retainedp == true) {
    /* Must be pre-mRNA */
    return UNSPLICED;
  } else if (splicedp == true) {
    /* Must be processed mRNA */
    return SPLICED;
  } else {
    /* Could be either */
    return BOTH;
  }
}


static void
compute_velocity5 (List_T transcripts5, Path_T path3) {
  List_T p, q, r;
  T transcript5, transcript3, invalid_transcript3;
  unsigned int transcript3_num, invalid_transcript3_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts5);
  Transcript_list_ascendingp(path3->transcripts);
  Transcript_list_ascendingp(path3->invalid_transcripts);
#endif

  p = transcripts5;
  q = path3->transcripts;
  r = path3->invalid_transcripts;

  while (p != NULL && (q != NULL || r != NULL)) {
    if (q == NULL) {
      transcript3_num = (unsigned int) -1U;
    } else {
      transcript3 = (T) List_head(q);
      transcript3_num = transcript3->num;
    }

    if (r == NULL) {
      invalid_transcript3_num = (unsigned int) -1U;
    } else {
      invalid_transcript3 = (T) List_head(r);
      invalid_transcript3_num = invalid_transcript3->num;
    }
    
    transcript5 = (T) List_head(p);
    if (transcript3_num < invalid_transcript3_num) {
      if (transcript5->num < transcript3_num) {
	transcript5->velocity = velocity_single(transcript5);
	p = List_next(p);
      } else if (transcript3_num < transcript5->num) {
	q = List_next(q);
      } else {
	transcript5->velocity = velocity_paired(transcript5,transcript3);
	p = List_next(p);
	q = List_next(q);
      }

    } else if (invalid_transcript3_num < transcript3_num) {
      if (transcript5->num < invalid_transcript3_num) {
	transcript5->velocity = velocity_single(transcript5);
	p = List_next(p);
      } else if (invalid_transcript3_num < transcript5->num) {
	r = List_next(r);
      } else {
	transcript5->velocity = velocity_paired(transcript5,invalid_transcript3);
	p = List_next(p);
	r = List_next(r);
      }

    } else {
      fprintf(stderr,"Same transcript %d is in both transcripts3 and invalid_transcripts3\n",
	      transcript3_num);
      /* Path_print(path3); */
      abort();
    }
  }

  while (p != NULL) {
    transcript5 = (T) List_head(p);
    transcript5->velocity = velocity_single(transcript5);
    p = List_next(p);
  }

  return;
}


static void
compute_velocity3 (List_T transcripts3, Path_T path5) {
  List_T p, q, r;
  T transcript5, invalid_transcript5, transcript3;
  unsigned int transcript5_num, invalid_transcript5_num;

#ifdef CHECK_ASSERTIONS
  Transcript_list_ascendingp(transcripts3);
  Transcript_list_ascendingp(path5->transcripts);
  Transcript_list_ascendingp(path5->invalid_transcripts);
#endif

  p = transcripts3;
  q = path5->transcripts;
  r = path5->invalid_transcripts;

  while (p != NULL && (q != NULL || r != NULL)) {
    if (q == NULL) {
      transcript5_num = (unsigned int) -1U;
    } else {
      transcript5 = (T) List_head(q);
      transcript5_num = transcript5->num;
    }

    if (r == NULL) {
      invalid_transcript5_num = (unsigned int) -1U;
    } else {
      invalid_transcript5 = (T) List_head(r);
      invalid_transcript5_num = invalid_transcript5->num;
    }
    
    transcript3 = (T) List_head(p);
    if (transcript5_num < invalid_transcript5_num) {
      if (transcript3->num < transcript5_num) {
	transcript3->velocity = velocity_single(transcript3);
	p = List_next(p);
      } else if (transcript5_num < transcript3->num) {
	q = List_next(q);
      } else {
	transcript3->velocity = velocity_paired(transcript5,transcript3);
	p = List_next(p);
	q = List_next(q);
      }

    } else if (invalid_transcript5_num < transcript5_num) {
      if (transcript3->num < invalid_transcript5_num) {
	transcript3->velocity = velocity_single(transcript3);
	p = List_next(p);
      } else if (invalid_transcript5_num < transcript3->num) {
	r = List_next(r);
      } else {
	transcript3->velocity = velocity_paired(invalid_transcript5,transcript3);
	p = List_next(p);
	r = List_next(r);
      }

    } else {
      fprintf(stderr,"Same transcript %d is in both transcripts5 and invalid_transcripts5\n",
	      transcript5_num);
      /* Path_print(path5); */
      abort();
    }
  }

  while (p != NULL) {
    transcript3 = (T) List_head(p);
    transcript3->velocity = velocity_single(transcript3);
    p = List_next(p);
  }

  return;
}

void
Transcript_velocity_paired (Path_T path5, Path_T path3) {
  List_T p;
  T transcript;

  compute_velocity5(path5->transcripts,path3);
  compute_velocity5(path5->invalid_transcripts,path3);
  compute_velocity3(path3->transcripts,path5);
  compute_velocity3(path3->invalid_transcripts,path5);

  /* Handle fusion transcripts as single */
  for (p = path5->fusion_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  for (p = path5->fusion_invalid_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  for (p = path3->fusion_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  for (p = path3->fusion_invalid_transcripts; p != NULL; p = List_next(p)) {
    transcript = (T) List_head(p);
    transcript->velocity = velocity_single(transcript);
  }

  return;
}


