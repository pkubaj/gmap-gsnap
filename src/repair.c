static char rcsid[] = "$Id: d26adf8a21da8418a0fbb8e0e3642f2293ef0a13 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "repair.h"
#include <math.h>		/* For qsort */
#include "mem.h"
#include "path-eval.h"
#include "genomebits_count.h"


/* Repair_path */
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Repair_unique */
#ifdef DEBUG1
#define debug1(x) x
#else
#define debug1(x)
#endif

static Genomebits_T genomebits;
static Genomebits_T genomebits_alt;


#define T Repair_T


void
Repair_free (T *old, Listpool_T listpool, Transcriptpool_T transcriptpool,
	     bool free_transcripts_p) {

  if (free_transcripts_p == true) {
    Transcript_list_gc(&(*old)->transcripts,listpool,transcriptpool);
  }
  FREE(*old);
  return;
}


T
Repair_new (int transcript_genestrand,
	    int tstart_overhang, Chrpos_T tstart_splice_distance,
	    int trend_overhang, Chrpos_T trend_splice_distance,
	    Transcript_T transcript, Listpool_T listpool) {
  T new = (T) MALLOC(sizeof(*new));

  new->transcript_genestrand = transcript_genestrand;
  new->trstart_overhang = tstart_overhang;
  new->trstart_splice_distance = tstart_splice_distance;
  new->trend_overhang = trend_overhang;
  new->trend_splice_distance = trend_splice_distance;
  new->transcripts = Listpool_push(NULL,listpool,(void *) transcript
				   listpool_trace(__FILE__,__LINE__));
 
  return new;
}
  

static int
repair_cmp (const void *a, const void *b) {
  T x = * (T *) a;
  T y = * (T *) b;

  if (x->transcript_genestrand > y->transcript_genestrand) {
    return -1;
  } else if (y->transcript_genestrand > x->transcript_genestrand) {
    return +1;
  } else if (x->trstart_overhang < y->trstart_overhang) {
    return -1;
  } else if (y->trstart_overhang < x->trstart_overhang) {
    return +1;
  } else if (x->trstart_splice_distance < y->trstart_splice_distance) {
    return -1;
  } else if (y->trstart_splice_distance < x->trstart_splice_distance) {
    return +1;
  } else if (x->trend_overhang < y->trend_overhang) {
    return -1;
  } else if (y->trend_overhang < x->trend_overhang) {
    return +1;
  } else if (x->trend_splice_distance < y->trend_splice_distance) {
    return -1;
  } else if (y->trend_splice_distance < x->trend_splice_distance) {
    return +1;
  } else {
    return 0;
  }
}


List_T
Repair_make_unique (List_T repairs, Listpool_T listpool, Transcriptpool_T transcriptpool) {
  List_T unique = NULL;
  T *array;
  int n, i, j;

  debug1(printf("Entered Repair_unique with %d repairs\n",
		List_length(repairs)));

  if ((n = List_length(repairs)) == 0) {
    return (List_T) NULL;

  } else {
    array = (T *) List_to_array(repairs,NULL);
    qsort(array,n,sizeof(T),repair_cmp);

    i = 0;
    while (i < n) {
      debug1(printf("Storing repair %d/%u, %d/%u\n",
		    array[i]->trstart_overhang,array[i]->trstart_splice_distance,
		    array[i]->trend_overhang,array[i]->trend_splice_distance));

      unique = Listpool_push(unique,listpool,(void *) array[i]
			     listpool_trace(__FILE__,__LINE__));
      j = i + 1;
      while (j < n && repair_cmp(&(array[j]),&(array[i])) == 0) {
	array[i]->transcripts = List_append(array[j]->transcripts,array[i]->transcripts);
	Repair_free(&(array[j]),listpool,transcriptpool,/*free_transcripts_p*/false);
	j++;
      }
      /* TODO: See if transcripts are in order */

      i = j;
    }
    FREE(array);

    Listpool_free_list(&repairs,listpool
		       listpool_trace(__FILE__,__LINE__));

    return List_reverse(unique);
  }
}


Path_T
Repair_path (int *found_score, T this, Path_T oldpath, int sensedir,
	     Compress_T query_compress,
	     Compress_T query_compress_fwd, Compress_T query_compress_rev,
	     Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
	     Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
	     Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {
  Path_T path;
  List_T p;
  Transcript_T transcript;
  int querypos;
  Univcoord_T distal_univdiagonal;
  int nmismatches, ref_nmismatches;


  path = Path_copy(oldpath,intlistpool,univcoordlistpool,listpool,
		   pathpool,vectorpool,transcriptpool,hitlistpool);

#ifdef DEBUG
  printf("Entered Repair_path, transcript genestrand %d, with %d/%u and %d/%u\n",
	 this->transcript_genestrand,this->trstart_overhang,this->trstart_splice_distance,
	 this->trend_overhang,this->trend_splice_distance);
  Transcript_print_nums(this->transcripts);
  Path_print(oldpath);
#endif

  if (this->transcript_genestrand > 0) {
    if (this->trstart_overhang > 0) {
      /* Don't know if the overhang actually aligns to the previous exon */
      distal_univdiagonal = Univcoordlist_head(path->univdiagonals) - this->trstart_splice_distance;
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  distal_univdiagonal,path->querylength,
							  /*pos5*/0,/*pos3*/this->trstart_overhang,
							  path->plusp,path->genestrand);
      debug(printf("(1) On a trstart_overhang of %d, got %d nmismatches\n",this->trstart_overhang,nmismatches));
      if (/*support*/this->trstart_overhang - 4*nmismatches < 0) {
	debug(printf("Insufficient support for this repair, so returning NULL\n"));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (Path_T) NULL;
      } else {
	if (Intlist_head(path->endpoints) == 0 && Intlist_head(path->nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->nmismatches,-1);
	}
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches
					     intlistpool_trace(__FILE__,__LINE__));

	if (Intlist_head(path->endpoints) == 0 && Intlist_head(path->ref_nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->ref_nmismatches,-1);
	}
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
      }
      
      querypos = Intlist_second_value(path->endpoints);
      /* No need to correct for Junction_ninserts on qstart side */

      if (this->trstart_overhang >= querypos) {
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	debug(printf("tstart_overhang %d >= querypos %d, so returning NULL\n",this->trstart_overhang,querypos));
	return (Path_T) NULL;
      } else {
	Intlist_head_set(path->endpoints,this->trstart_overhang);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,0
					   intlistpool_trace(__FILE__,__LINE__));
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
      }

      path->univdiagonals =
	Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,distal_univdiagonal
			       univcoordlistpool_trace(__FILE__,__LINE__));

      path->junctions = Listpool_push(path->junctions,listpool,
				      (void *) Junction_new_splice(this->trstart_splice_distance,
								   sensedir,/*donor_prob*/2.0,
								   /*acceptor_prob*/2.0,pathpool)
				      listpool_trace(__FILE__,__LINE__));

      Altsplice_free(&path->qstart_alts,pathpool);
      path->qstart_alts = (Altsplice_T) NULL;
    }

    if (this->trend_overhang > 0) {
      path->endpoints = Intlist_reverse(path->endpoints);
      path->univdiagonals = Univcoordlist_reverse(path->univdiagonals);
      path->nmismatches = Intlist_reverse(path->nmismatches);
      path->ref_nmismatches = Intlist_reverse(path->ref_nmismatches);
      path->junctions = List_reverse(path->junctions);

      /* Don't know if the overhang actually aligns to the next exon */
      distal_univdiagonal = Univcoordlist_head(path->univdiagonals) + this->trend_splice_distance;
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  distal_univdiagonal,path->querylength,
							  /*pos5*/path->querylength - this->trend_overhang,
							  /*pos3*/path->querylength,path->plusp,path->genestrand);
      debug(printf("(2) On a trend_overhang of %d, got %d nmismatches\n",this->trend_overhang,nmismatches));
      if (/*support*/this->trstart_overhang - 4*nmismatches < 0) {
	debug(printf("Insufficient support for this repair, so returning NULL\n"));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (Path_T) NULL;
      } else {
	if (Intlist_head(path->endpoints) == path->querylength && Intlist_head(path->nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->nmismatches,-1);
	}
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches
					     intlistpool_trace(__FILE__,__LINE__));

	if (Intlist_head(path->endpoints) == path->querylength && Intlist_head(path->ref_nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->ref_nmismatches,-1);
	}
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
      }

      querypos = Intlist_second_value(path->endpoints);
      if (path->junctions != NULL) {
	querypos += Junction_ninserts((Junction_T) List_head(path->junctions));
      }

      if (path->querylength - this->trend_overhang <= querypos) {
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	debug(printf("trend_overhang %d <= querypos %d, so returning NULL\n",this->trend_overhang,querypos));
	return (Path_T) NULL;
      } else {
	Intlist_head_set(path->endpoints,path->querylength - this->trend_overhang);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,path->querylength
					   intlistpool_trace(__FILE__,__LINE__));
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
      }

      path->univdiagonals =
	Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,distal_univdiagonal
			       univcoordlistpool_trace(__FILE__,__LINE__));

      path->junctions = Listpool_push(path->junctions,listpool,
				      (void *) Junction_new_splice(this->trend_splice_distance,
								   sensedir,/*donor_prob*/2.0,
								   /*acceptor_prob*/2.0,pathpool)
				      listpool_trace(__FILE__,__LINE__));

      path->endpoints = Intlist_reverse(path->endpoints);
      path->univdiagonals = Univcoordlist_reverse(path->univdiagonals);
      path->nmismatches = Intlist_reverse(path->nmismatches);
      path->ref_nmismatches = Intlist_reverse(path->ref_nmismatches);
      path->junctions = List_reverse(path->junctions);

      Altsplice_free(&path->qend_alts,pathpool);
      path->qend_alts = (Altsplice_T) NULL;
    }

  } else {
    if (this->trstart_overhang > 0) {
      path->endpoints = Intlist_reverse(path->endpoints);
      path->univdiagonals = Univcoordlist_reverse(path->univdiagonals);
      path->nmismatches = Intlist_reverse(path->nmismatches);
      path->ref_nmismatches = Intlist_reverse(path->ref_nmismatches);
      path->junctions = List_reverse(path->junctions);
      
      /* Don't know if the overhang actually aligns to the previous exon */
      distal_univdiagonal = Univcoordlist_head(path->univdiagonals) + this->trstart_splice_distance;
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  distal_univdiagonal,path->querylength,
							  /*pos5*/path->querylength - this->trstart_overhang,
							  /*pos3*/path->querylength,path->plusp,path->genestrand);
      debug(printf("(3) On a trstart_overhang of %d, got %d nmismatches\n",this->trstart_overhang,nmismatches));
      if (/*support*/this->trstart_overhang - 4*nmismatches < 0) {
	debug(printf("Insufficient support for this repair, so returning NULL\n"));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (Path_T) NULL;
      } else {
	if (Intlist_head(path->endpoints) == path->querylength && Intlist_head(path->nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->nmismatches,-1);
	}
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches
					     intlistpool_trace(__FILE__,__LINE__));

	if (Intlist_head(path->endpoints) == path->querylength && Intlist_head(path->ref_nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->ref_nmismatches,-1);
	}
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
      }

      querypos = Intlist_second_value(path->endpoints);
      if (path->junctions != NULL) {
	querypos += Junction_ninserts((Junction_T) List_head(path->junctions));
      }

      if (path->querylength - this->trstart_overhang <= querypos) {
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	debug(printf("tstart_overhang %d <= querypos %d, so returning NULL\n",this->trstart_overhang,querypos));
	return (Path_T) NULL;
      } else {
	Intlist_head_set(path->endpoints,path->querylength - this->trstart_overhang);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,path->querylength
					   intlistpool_trace(__FILE__,__LINE__));
	path->splice3p = false;
	path->splicetype3 = NO_SPLICE;
	path->ambig_prob_3 = 0.0;
      }

      path->univdiagonals =
	Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,distal_univdiagonal
			       univcoordlistpool_trace(__FILE__,__LINE__));

      path->junctions = Listpool_push(path->junctions,listpool,
				      (void *) Junction_new_splice(this->trstart_splice_distance,
								   sensedir,/*donor_prob*/2.0,
								   /*acceptor_prob*/2.0,pathpool)
				      listpool_trace(__FILE__,__LINE__));

      path->endpoints = Intlist_reverse(path->endpoints);
      path->univdiagonals = Univcoordlist_reverse(path->univdiagonals);
      path->nmismatches = Intlist_reverse(path->nmismatches);
      path->ref_nmismatches = Intlist_reverse(path->ref_nmismatches);
      path->junctions = List_reverse(path->junctions);

      Altsplice_free(&path->qend_alts,pathpool);
      path->qend_alts = (Altsplice_T) NULL;
    }

    if (this->trend_overhang > 0) {
      /* Don't know if the overhang actually aligns to the previous exon */
      distal_univdiagonal = Univcoordlist_head(path->univdiagonals) - this->trend_splice_distance;
      nmismatches = Genomebits_count_mismatches_substring(&ref_nmismatches,genomebits,genomebits_alt,query_compress,
							  distal_univdiagonal,path->querylength,
							  /*pos5*/0,/*pos3*/this->trend_overhang,
							  path->plusp,path->genestrand);
      debug(printf("(4) On a trend_overhang of %d, got %d nmismatches\n",this->trend_overhang,nmismatches));
      if (/*support*/this->trend_overhang - 4*nmismatches < 0) {
	debug(printf("Insufficient support for this repair, so returning NULL\n"));
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	return (Path_T) NULL;
      } else {
	if (Intlist_head(path->endpoints) == 0 && Intlist_head(path->nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->nmismatches,-1);
	}
	path->nmismatches = Intlistpool_push(path->nmismatches,intlistpool,nmismatches
					     intlistpool_trace(__FILE__,__LINE__));
      
	if (Intlist_head(path->endpoints) && Intlist_head(path->ref_nmismatches) == 0) {
	  /* Medial portion is still an exact match */
	} else {
	  Intlist_head_set(path->ref_nmismatches,-1);
	}
	path->ref_nmismatches = Intlistpool_push(path->ref_nmismatches,intlistpool,ref_nmismatches
						 intlistpool_trace(__FILE__,__LINE__));
      }
	
      querypos = Intlist_second_value(path->endpoints);
      /* No need to correct for Junction_ninserts on qstart side */

      if (this->trend_overhang >= querypos) {
	Path_free(&path,intlistpool,univcoordlistpool,
		  listpool,pathpool,transcriptpool,hitlistpool);
	debug(printf("trend_overhang %d >= querypos %d, so returning NULL\n",this->trend_overhang,querypos));
	return (Path_T) NULL;
      } else {
	Intlist_head_set(path->endpoints,this->trend_overhang);
	path->endpoints = Intlistpool_push(path->endpoints,intlistpool,0
					   intlistpool_trace(__FILE__,__LINE__));
	path->splice5p = false;
	path->splicetype5 = NO_SPLICE;
	path->ambig_prob_5 = 0.0;
      }

      path->univdiagonals =
	Univcoordlistpool_push(path->univdiagonals,univcoordlistpool,distal_univdiagonal
			       univcoordlistpool_trace(__FILE__,__LINE__));

      path->junctions = Listpool_push(path->junctions,listpool,
				      (void *) Junction_new_splice(this->trend_splice_distance,
								   sensedir,/*donor_prob*/2.0,
								   /*acceptor_prob*/2.0,pathpool)
				      listpool_trace(__FILE__,__LINE__));

      Altsplice_free(&path->qstart_alts,pathpool);
      path->qstart_alts = (Altsplice_T) NULL;
    }
  }
    
  Path_eval_nmatches(&(*found_score),path,query_compress_fwd,query_compress_rev);
#if 0
  if (path->score_within_trims > nmismatches_allowed) {
    Path_free(&path,intlistpool,univcoordlistpool,
	      listpool,pathpool,transcriptpool,hitlistpool);
    debug(printf("score_within_trims %d > nmismatches_allowed %d, so returning NULL\n",
		 path->score_within_trims,nmismatches_allowed));

    /* These transcripts are not applicable at all */
    oldpath->invalid_transcripts = Transcript_remove_subset(oldpath->invalid_transcripts,this->transcripts,
							    listpool,transcriptpool);

    path = (Path_T) NULL;
  } else {
#endif

    /* Extend transcripts */
    debug(printf("score_within_trims %d\n",path->score_within_trims));
    for (p = this->transcripts; p != NULL; p = List_next(p)) {
      transcript = (Transcript_T) List_head(p);
      if (this->trstart_overhang > 0) {
	Transcript_repair_trstart(transcript,listpool,transcriptpool);
      }
      if (this->trend_overhang > 0) {
	Transcript_repair_trend(transcript,listpool,transcriptpool);
      }
    }

    /* Make repair->transcripts valid ones for this path */
    path->transcripts = this->transcripts;

    /* Overkill if more than one repair is applicable to a path, since
       this removes all invalid transcripts */
    Transcript_list_gc(&path->invalid_transcripts,listpool,transcriptpool);

#ifdef DEBUG
    Path_print(path);
    printf("Exiting Repair_path\n");
#endif

#if 0
  }
#endif

  return path;
}


void
Repair_setup (Genomebits_T genomebits_in, Genomebits_T genomebits_alt_in) {

  genomebits = genomebits_in;
  genomebits_alt = genomebits_alt_in;

  return;
}
