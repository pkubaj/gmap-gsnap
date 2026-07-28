#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif

#ifdef HAVE_PTHREAD
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* Needed to define pthread_t on Solaris */
#endif
#include <pthread.h>
#endif


#include "support.h"

#include <ctype.h>
#include <string.h>

#include "assert.h"
#include "mem.h"
#include "intlist.h"
#include "list.h"
#include "junction.h"
#include "splice-count.h"
#include "maxent_hr.h"
#include "getline.h"


#if 0
#include "complement.h"
static char complCode[128] = COMPLEMENT_UC;
#endif


#define INTRON_ADJ 1		/* For reporting intron ends in Support_dump, as STAR does */


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Support_dump univcoords */
#ifdef DEBUG8
#define debug8(x) x
#else
#define debug8(x)
#endif


/* Support_check */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


#define T Support_T
struct T {

#ifdef HAVE_PTHREAD
  pthread_mutex_t lock;
#endif
  
  /* For updating counts and support, and for final dump */
  Univcoordpairtable_T plus_table;
  Univcoordpairtable_T minus_table;

  /* For faster lookup in pass 2 */
  Knownsplicing_T pass1_splicing;
};


T
Support_new (Knownsplicing_T knownsplicing) {
  T new = (T) MALLOC_KEEP(sizeof(*new));
  int i, k;

#ifdef HAVE_PTHREAD
  pthread_mutex_init(&new->lock,NULL);
#endif
  
  new->plus_table = Univcoordpairtable_new(/*hint*/100000);
  new->minus_table = Univcoordpairtable_new(/*hint*/100000);

  if (knownsplicing != NULL) {
    k = 0;
    for (i = 0; i < knownsplicing->donor_nintervals; i++) {
      /* printf("Putting known plus intron at %u..%u\n",
	 knownsplicing->donor_endpoints[k],knownsplicing->donor_endpoints[k+1]); */
      Univcoordpairtable_put(new->plus_table,
			     /*donor_position*/knownsplicing->donor_endpoints[k],
			     /*acceptor_position*/knownsplicing->donor_endpoints[k+1],
			     (void *) Splice_count_new(/*knownp*/true));
      k += 2;
    }

    k = 0;
    for (i = 0; i < knownsplicing->antiacceptor_nintervals; i++) {
      /*printf("Putting known minus intron at %u..%u\n",
	knownsplicing->antiacceptor_endpoints[k],knownsplicing->antiacceptor_endpoints[k+1]); */
      Univcoordpairtable_put(new->minus_table,
			     /*acceptor_position*/knownsplicing->antiacceptor_endpoints[k],
			     /*donor_position*/knownsplicing->antiacceptor_endpoints[k+1],
			     (void *) Splice_count_new(/*knownp*/true));
      k += 2;
    }
  }

  new->pass1_splicing = (Knownsplicing_T) NULL;

  return new;
}


void
Support_free (T *old) {
  int nintrons, i;
  Splice_count_T *splice_counts;

  if (*old) {
#ifdef HAVE_PTHREAD
    pthread_mutex_destroy(&(*old)->lock);
#endif

    if ((nintrons = Univcoordpairtable_length((*old)->plus_table)) > 0) {
      splice_counts = (Splice_count_T *) Univcoordpairtable_values((*old)->plus_table);
      for (i = 0; i < nintrons; i++) {
	Splice_count_free(&(splice_counts[i]));
      }
      FREE(splice_counts);
    }
    Univcoordpairtable_free(&(*old)->plus_table);

    if ((nintrons = Univcoordpairtable_length((*old)->minus_table)) > 0) {
      splice_counts = (Splice_count_T *) Univcoordpairtable_values((*old)->minus_table);
      for (i = 0; i < nintrons; i++) {
	Splice_count_free(&(splice_counts[i]));
      }
      FREE(splice_counts);
    }
    Univcoordpairtable_free(&(*old)->minus_table);

    if ((*old)->pass1_splicing != NULL) {
      Knownsplicing_free(&(*old)->pass1_splicing);
    }

    FREE_KEEP(*old);
  }

  return;
}

    
static void
Support_revise (T this, Path_T path, int npaths) {
  Intlist_T p;
  List_T q;
  Univcoordlist_T u;
  Univcoord_T donor_position, acceptor_position, last_univdiagonal, next_univdiagonal;
  Junction_T junction;
  int support_length, splice_qpos;
  Splice_count_T splice_count;
  bool sense_forward_p;
  int overall_qstart, overall_qend;

  debug(printf("Support_revise with npaths %d and path ",npaths));
  debug(Path_print(path));

  if (path->sensedir == SENSE_FORWARD) {
    sense_forward_p = true;
  } else {
    sense_forward_p = false;
  }

  overall_qstart = Intlist_head(path->endpoints);
  overall_qend = Intlist_last_value(path->endpoints);

  next_univdiagonal = Univcoordlist_head(path->univdiagonals);
  for (p = Intlist_next(path->endpoints), q = path->junctions, u = path->univdiagonals;
       Intlist_next(p) != NULL; p = Intlist_next(p), q = List_next(q), u = Univcoordlist_next(u)) {
    last_univdiagonal = next_univdiagonal;
    next_univdiagonal = Univcoordlist_head(Univcoordlist_next(u));

    junction = (Junction_T) List_head(q);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      splice_qpos = Intlist_head(p);
      if (splice_qpos - overall_qstart < overall_qend - splice_qpos) {
	support_length = splice_qpos - overall_qstart;
      } else {
	support_length = overall_qend - splice_qpos;
      }

      if (path->plusp == sense_forward_p) {
	/* plus */
	donor_position = last_univdiagonal - path->querylength + splice_qpos;
	acceptor_position = next_univdiagonal - path->querylength + splice_qpos;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&this->lock);
#endif
	if ((splice_count = Univcoordpairtable_get(this->plus_table,
						   /*key1*/donor_position,/*key2*/acceptor_position)) == NULL) {
	  splice_count = Splice_count_new(/*knownp*/false);
	  Univcoordpairtable_put(this->plus_table,donor_position,acceptor_position,splice_count);
	}

	if (npaths > 1) {
	  splice_count->nmulti += 1;
	} else {
	  splice_count->nunique += 1;
	}

	/* Need to count support for multimappers, in case genome has patches */
	if (support_length > splice_count->support) {
	  debug(printf("Updating support at %u..%u to be %d\n",donor_position,acceptor_position,support_length));
	  splice_count->support = support_length;
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&this->lock);
#endif

      } else {
	/* minus */
	acceptor_position = last_univdiagonal - path->querylength + splice_qpos;
	donor_position = next_univdiagonal - path->querylength + splice_qpos;

#ifdef HAVE_PTHREAD
	pthread_mutex_lock(&this->lock);
#endif
	if ((splice_count = Univcoordpairtable_get(this->minus_table,
						   /*key1*/acceptor_position,/*key2*/donor_position)) == NULL) {
	  splice_count = Splice_count_new(/*knownp*/false);
	  Univcoordpairtable_put(this->minus_table,acceptor_position,donor_position,splice_count);
	}

	if (npaths > 1) {
	  splice_count->nmulti += 1;
	} else {
	  splice_count->nunique += 1;
	}

	/* Need to count support for multimappers, in case genome has patches */
	if (support_length > splice_count->support) {
	  debug(printf("Updating support at %u..%u to be %d\n",donor_position,acceptor_position,support_length));
	  splice_count->support = support_length;
	}

#ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&this->lock);
#endif
      }
    }
  }

  return;
}


void
Support_revise_single (T this, Path_T *patharray, int npaths) {
  int i;

  for (i = 0; i < npaths; i++) {
    Support_revise(this,patharray[i],npaths);
  }

  return;
}


void
Support_revise_paired (T this, Pathpair_T *pathpairarray, int npaths) {
  int i;

  for (i = 0; i < npaths; i++) {
    Support_revise(this,pathpairarray[i]->pathL,npaths);
    Support_revise(this,pathpairarray[i]->pathH,npaths);
  }

  return;
}


/* Needed by Trpath_convert procedures, if USE_SUPPORT_CHECK is defined */
bool
Support_check (T this, Path_T path) {
  Intlist_T p;
  List_T q;
  Univcoordlist_T u;
  Univcoord_T donor_position, acceptor_position, last_univdiagonal, next_univdiagonal;
  Junction_T junction;
  int splice_qpos;
  Splice_count_T splice_count;
  bool sense_forward_p;


  debug9(printf("Support_check with path "));
  debug9(Path_print(path));

  if (path->sensedir == SENSE_FORWARD) {
    sense_forward_p = true;
  } else {
    sense_forward_p = false;
  }

  next_univdiagonal = Univcoordlist_head(path->univdiagonals);
  for (p = Intlist_next(path->endpoints), q = path->junctions, u = path->univdiagonals;
       Intlist_next(p) != NULL; p = Intlist_next(p), q = List_next(q), u = Univcoordlist_next(u)) {
    last_univdiagonal = next_univdiagonal;
    next_univdiagonal = Univcoordlist_head(Univcoordlist_next(u));

    junction = (Junction_T) List_head(q);
    if (Junction_type(junction) == SPLICE_JUNCTION) {
      splice_qpos = Intlist_head(p);

      if (path->plusp == sense_forward_p) {
	/* plus */
	donor_position = last_univdiagonal - path->querylength + splice_qpos;
	acceptor_position = next_univdiagonal - path->querylength + splice_qpos;

#ifdef HAVE_PTHREAD
	/* pthread_mutex_lock(&this->lock); -- read-only */
#endif
	if ((splice_count = Univcoordpairtable_get(this->plus_table,
						   /*key1*/donor_position,/*key2*/acceptor_position)) == NULL) {
	  debug9(printf("For plus intron %u..%u, no splice_count found\n",
			donor_position,acceptor_position));
	  return false;

	} else if (splice_count->support < 20) {
	  /* Same as in Splice_accept_p */
	  debug9(printf("For plus intron %u..%u, support is only %d\n",
			donor_position,acceptor_position,splice_count->support));
	  return false;
	}

#ifdef HAVE_PTHREAD
	/* pthread_mutex_unlock(&this->lock); -- read-only */
#endif

      } else {
	/* minus */
	acceptor_position = last_univdiagonal - path->querylength + splice_qpos;
	donor_position = next_univdiagonal - path->querylength + splice_qpos;

#ifdef HAVE_PTHREAD
	/* pthread_mutex_lock(&this->lock); -- read-only */
#endif
	if ((splice_count = Univcoordpairtable_get(this->minus_table,
						   /*key1*/acceptor_position,/*key2*/donor_position)) == NULL) {
	  debug9(printf("For minus intron %u..%u, no splice_count found\n",
			donor_position,acceptor_position));
	  return false;

	} else if (splice_count->support < 20) {
	  /* Same as in Splice_accept_p */
	  debug9(printf("For minus intron %u..%u, support is only %d\n",
			donor_position,acceptor_position,splice_count->support));
	  return false;
	}

#ifdef HAVE_PTHREAD
	/* pthread_mutex_unlock(&this->lock); -- read-only */
#endif
      }
    }
  }

  return true;
}


void
Support_reset_counts (T this) {
  Splice_count_T *splice_counts;
  int nintrons, i;

  /* plus */
  if ((nintrons = Univcoordpairtable_length(this->plus_table)) > 0) {
    splice_counts = (Splice_count_T *) Univcoordpairtable_values(this->plus_table);
    for (i = 0; i < nintrons; i++) {
      splice_counts[i]->nunique = 0;
      splice_counts[i]->nmulti = 0;
    }
    FREE(splice_counts);
  }

  /* minus */
  if ((nintrons = Univcoordpairtable_length(this->minus_table)) > 0) {
    splice_counts = (Splice_count_T *) Univcoordpairtable_values(this->minus_table);
    for (i = 0; i < nintrons; i++) {
      splice_counts[i]->nunique = 0;
      splice_counts[i]->nmulti = 0;
    }
    FREE(splice_counts);
  }
  
  return;
}



/* Filters introns from pass 1, based on unique mappers, support, and
   splice site probabilities.  Also converts Univcoordpairtable
   information into an EF64_T format */
void
Support_filter (T this, EF64_T chromosome_ef64, Univcoord_T genomelength, int threshold) {
  Univcoordlist_T donor_startpoints = NULL, donor_partners = NULL,
    acceptor_startpoints = NULL, acceptor_partners = NULL,
    antidonor_startpoints = NULL, antidonor_partners = NULL,
    antiacceptor_startpoints = NULL, antiacceptor_partners = NULL;

  Univcoord_T *donor_positions, *acceptor_positions, donor_position, acceptor_position;
  int nintrons, i;
  Splice_count_T splice_count;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  double donor_prob, acceptor_prob;


  /* plus */
  Univcoordpairtable_keys(&donor_positions,&acceptor_positions,this->plus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->plus_table);

  chrnum = 1;
  EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,/*chrnum*/1);

  for (i = 0; i < nintrons; i++) {
    donor_position = donor_positions[i];
    acceptor_position = acceptor_positions[i];

    splice_count = Univcoordpairtable_get(this->plus_table,donor_position,acceptor_position);

#if 0
    if (splice_count->nunique == 0) {
      /* Skip: requiring unique mappers may be counterproductive in small
	 datasets, but --two-pass is likely to be run on a large
	 dataset.  However, genome could have duplicate regions. */
    }
#endif

    if (splice_count->support < threshold) {
      /* Skip */
    } else {
      if (donor_position >= chrhigh) {
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,donor_position,donor_position);
      }

      donor_prob = Maxent_hr_donor_prob(donor_position,chroffset);
      acceptor_prob = Maxent_hr_acceptor_prob(acceptor_position,chroffset);

      if (donor_prob > 0.5 && acceptor_prob > 0.5) {
	donor_startpoints = Univcoordlist_push(donor_startpoints,donor_position);
	donor_partners = Univcoordlist_push(donor_partners,acceptor_position);
    
	acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,acceptor_position);
	acceptor_partners = Univcoordlist_push(acceptor_partners,donor_position);
      }
    }
  }

  FREE(acceptor_positions);
  FREE(donor_positions);


  /* minus */
  Univcoordpairtable_keys(&acceptor_positions,&donor_positions,this->minus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->minus_table);

  chrnum = 1;
  EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,/*chrnum*/1);

  for (i = 0; i < nintrons; i++) {
    acceptor_position = acceptor_positions[i];
    donor_position = donor_positions[i];

    splice_count = Univcoordpairtable_get(this->minus_table,acceptor_position,donor_position);

#if 0
    if (splice_count->nunique == 0) {
      /* Skip: requiring unique mappers may be counterproductive in
	 small datasets, but --two-pass is likely to be run on a large
	 dataset.  However, genome could have duplicate regions. */
    }
#endif

    if (splice_count->support < threshold) {
      /* Skip */
    } else {
      if (acceptor_position >= chrhigh) {
	chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,acceptor_position,acceptor_position);
      }
      donor_prob = Maxent_hr_antidonor_prob(donor_position,chroffset);
      acceptor_prob = Maxent_hr_antiacceptor_prob(acceptor_position,chroffset);
      
      if (donor_prob > 0.5 && acceptor_prob > 0.5) {
	antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,donor_position);
	antidonor_partners = Univcoordlist_push(antidonor_partners,acceptor_position);
    
	antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,acceptor_position);
	antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,donor_position);
      }
    }
  }

  FREE(acceptor_positions);
  FREE(donor_positions);


  /* Knownsplicing_new sorts startpoints and partners */
  this->pass1_splicing = Knownsplicing_new(donor_startpoints,donor_partners,
					   acceptor_startpoints,acceptor_partners,
					   antidonor_startpoints,antidonor_partners,
					   antiacceptor_startpoints,antiacceptor_partners,
					   /*donor_table*/NULL,/*acceptor_table*/NULL,
					   /*antidonor_table*/NULL,/*antiacceptor_table*/NULL,
					   genomelength,/*chromosome_iit*/NULL,
					   /*chromosome_ef64*/NULL,/*intron_level_p*/true);

  return;
}


void
Support_convert (T this, Univcoord_T genomelength) {
  Univcoordlist_T donor_startpoints = NULL, donor_partners = NULL,
    acceptor_startpoints = NULL, acceptor_partners = NULL,
    antidonor_startpoints = NULL, antidonor_partners = NULL,
    antiacceptor_startpoints = NULL, antiacceptor_partners = NULL;

  Univcoord_T *donor_positions, *acceptor_positions, donor_position, acceptor_position;
  int nintrons, i;


  /* plus */
  Univcoordpairtable_keys(&donor_positions,&acceptor_positions,this->plus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->plus_table);

  for (i = 0; i < nintrons; i++) {
    donor_position = donor_positions[i];
    acceptor_position = acceptor_positions[i];

    donor_startpoints = Univcoordlist_push(donor_startpoints,donor_position);
    donor_partners = Univcoordlist_push(donor_partners,acceptor_position);
    
    acceptor_startpoints = Univcoordlist_push(acceptor_startpoints,acceptor_position);
    acceptor_partners = Univcoordlist_push(acceptor_partners,donor_position);
  }

  FREE(acceptor_positions);
  FREE(donor_positions);


  /* minus */
  Univcoordpairtable_keys(&acceptor_positions,&donor_positions,this->minus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->minus_table);

  for (i = 0; i < nintrons; i++) {
    acceptor_position = acceptor_positions[i];
    donor_position = donor_positions[i];

    antidonor_startpoints = Univcoordlist_push(antidonor_startpoints,donor_position);
    antidonor_partners = Univcoordlist_push(antidonor_partners,acceptor_position);
    
    antiacceptor_startpoints = Univcoordlist_push(antiacceptor_startpoints,acceptor_position);
    antiacceptor_partners = Univcoordlist_push(antiacceptor_partners,donor_position);
  }

  FREE(acceptor_positions);
  FREE(donor_positions);


  /* Knownsplicing_new sorts startpoints and partners */
  this->pass1_splicing = Knownsplicing_new(donor_startpoints,donor_partners,
					   acceptor_startpoints,acceptor_partners,
					   antidonor_startpoints,antidonor_partners,
					   antiacceptor_startpoints,antiacceptor_partners,
					   /*donor_table*/NULL,/*acceptor_table*/NULL,
					   /*antidonor_table*/NULL,/*antiacceptor_table*/NULL,
					   genomelength,/*chromosome_iit*/NULL,
					   /*chromosome_ef64*/NULL,/*intron_level_p*/true);

  return;
}


#if 0
/* Old method based on Univcoordpairtable_T, slower */
int
Support_get (T this,
	     Univcoord_T splice_position_i, Univcoord_T splice_position_j,
	     bool plusp, bool sense_forward_p) {
  Splice_count_T splice_count;

  if (plusp == sense_forward_p) {
    splice_count = (Splice_count_T) Univcoordpairtable_get(this->plus_table,
							   /*donor_position*/splice_position_i,
							   /*acceptor_position*/splice_position_j);
  } else {
    splice_count = (Splice_count_T) Univcoordpairtable_get(this->minus_table,
							   /*acceptor_position*/splice_position_i,
							   /*donor_position*/splice_position_j);
  }

  if (splice_count == NULL) {
    return 0;
  } else {
    return splice_count->support;
  }
}
#endif


static bool
same_chr_p (char *last_chr, char *line) {
  char *p = line, *q, c;
  int length;
  
  q = p;
  while ((c = *p++) != '\0' && !isspace(c)) ;
  length = (p - 1) - q;

  if (strncmp(q,last_chr,length) == 0) {
    return true;
      
  } else {
    return false;
  }
}

static char *
get_chr (char *line) {
  char *chr, *p = line, c;
  int length;

  while ((c = *p++) != '\0' && !isspace(c)) ; /* chr name */
  length = (p - 1) - line;

  chr = (char *) MALLOC((length+1)*sizeof(char));
  strncpy(chr,line,length);
  chr[length] = '\0';

  return chr;
}


static Splice_count_T
parse_line (bool *watsonp, Chrpos_T *low_chrpos, Chrpos_T *high_chrpos, char *line) {
  Splice_count_T splice_count;
  char *p = line, c;
  int nchars;
  int nunique, nmulti, support;
  bool knownp;

  while ((c = *p++) != '\0' && !isspace(c)) ; /* chr name */
  
  if (sscanf(p,"%u%n",&(*low_chrpos),&nchars) < 1) {
    fprintf(stderr,"Cannot parse coords from line %s => %s\n",line,p);
    exit(9);
  } else {
    p += nchars + 1;	/* Skip over low_chrpos */
  }

  if (sscanf(p,"%u%n",&(*high_chrpos),&nchars) < 1) {
    fprintf(stderr,"Cannot parse coords from line %s => %s\n",line,p);
    exit(9);
  } else {
    p += nchars + 1;	/* Skip over high_chrpos */
  }

  if ((c = *p++) == '+') {
    *watsonp = true;
  } else if (c == '-') {
    *watsonp = false;
  } else {
    fprintf(stderr,"Observed %c for strand.  Expecting + or -.  Treating as +\n",c);
    *watsonp = true;
  }
  if (*p++ != '\t') {
    fprintf(stderr,"Expecting tab after strand\n");
    exit(9);
  }

  while ((c = *p++) != '\0' && !isspace(c)) ; /* intron type */

  if ((c = *p++) == '1') {
    knownp = true;
  } else if (c == '0') {
    knownp = false;
  } else {
    fprintf(stderr,"Observed %c for annotated.  Expecting 0 or 1.  Treating as 0\n",c);
    knownp = false;
  }
  if (*p++ != '\t') {
    fprintf(stderr,"Expecting tab after annotated\n");
    exit(9);
  }

  if (sscanf(p,"%u%n",&nunique,&nchars) < 1) {
    fprintf(stderr,"Cannot parse coords from line %s => %s\n",line,p);
    exit(9);
  } else {
    p += nchars + 1;	/* Skip over nunique */
  }

  if (sscanf(p,"%u%n",&nmulti,&nchars) < 1) {
    fprintf(stderr,"Cannot parse coords from line %s => %s\n",line,p);
    exit(9);
  } else {
    p += nchars + 1;	/* Skip over nmulti */
  }

  if (sscanf(p,"%u%n",&support,&nchars) < 1) {
    fprintf(stderr,"Cannot parse coords from line %s => %s\n",line,p);
    exit(9);
  } else {
    p += nchars + 1;	/* Skip over support */
  }

  splice_count = Splice_count_new(knownp);
  splice_count->nunique = nunique;
  splice_count->nmulti = nmulti;
  splice_count->support = support;

  return splice_count;
}


/* Reads a Support_dump file */
T
Support_read (FILE *fp, EF64_T chromosome_ef64, Univ_IIT_T chromosome_iit) {
  T new = Support_new(/*knownsplicing*/false);

  char *line, *last_chr;
  Splice_count_T splice_count;
  bool watsonp;

  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
  Chrpos_T low_chrpos, high_chrpos;


  line = Getline(fp);
  while (line != NULL) {
    last_chr = get_chr(line);
    if ((chrnum = Univ_IIT_find_one(chromosome_iit,last_chr)) < 0) {
      fprintf(stderr,"Cannot find chromosome %s in chromosome IIT file, so skipping %s",
	      last_chr,line);

      /* First iteration */
      FREE(line);
      line = Getline(fp);

      while (line != NULL && same_chr_p(last_chr,line) == true) {
	/* Subsequent iterations */
	FREE(line);
	line = Getline(fp);
      }

    } else {
      /* Note: file has 1-based coordinates of intron ends.  Need to invert the statements in Support_dump */
      EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,chrnum);
    
      /* First iteration */
      splice_count = parse_line(&watsonp,&low_chrpos,&high_chrpos,line);
      if (watsonp == true) {
	Univcoordpairtable_put(new->plus_table,
			       /*donor_position*/chroffset + (Univcoord_T) (low_chrpos - INTRON_ADJ),
			       /*acceptor_position*/chroffset + (Univcoord_T) (high_chrpos - 1 + INTRON_ADJ),
			       splice_count);
      } else {
	Univcoordpairtable_put(new->minus_table,
			       /*acceptor_position*/chroffset + (Univcoord_T) (low_chrpos - INTRON_ADJ),
			       /*donor_position*/chroffset + (Univcoord_T) (high_chrpos - 1 + INTRON_ADJ),
			       splice_count);
      }

      FREE(line);
      line = Getline(fp);

      while (line != NULL && same_chr_p(last_chr,line) == true) {
	/* Subsequent iterations */
	splice_count = parse_line(&watsonp,&low_chrpos,&high_chrpos,line);
	if (watsonp == true) {
	  Univcoordpairtable_put(new->plus_table,
				 /*donor_position*/chroffset + (Univcoord_T) (low_chrpos - INTRON_ADJ),
				 /*acceptor_position*/chroffset + (Univcoord_T) (high_chrpos - 1 + INTRON_ADJ),
				 splice_count);
	} else {
	  Univcoordpairtable_put(new->minus_table,
				 /*acceptor_position*/chroffset + (Univcoord_T) (low_chrpos - INTRON_ADJ),
				 /*donor_position*/chroffset + (Univcoord_T) (high_chrpos - 1 + INTRON_ADJ),
				 splice_count);
	}

	FREE(line);
	line = Getline(fp);
      }
    }

    FREE(last_chr);
    last_chr = (char *) NULL;
  }

  FREE(last_chr);

  return new;
}


bool
Support_presentp (T this,
		  Univcoord_T splice_position_i, Univcoord_T splice_position_j,
		  bool plusp, bool sense_forward_p) {

  assert(splice_position_i < splice_position_j);
  if (plusp == sense_forward_p) {
    return Knownsplicing_intron_p(this->pass1_splicing,
				  /*donor_position*/splice_position_i,
				  /*acceptor_position*/splice_position_j);
  } else {
    return Knownsplicing_antiintron_p(this->pass1_splicing,
				      /*acceptor_position*/splice_position_i,
				      /*donor_position*/splice_position_j);
  }
}


static inline int
dinucleotide_coded (Univcoord_T low_position, Univcoord_T high_position,
		    Genomebits_T ref, bool watsonp) {
  char low1, low2, high2, high1;

  low1 = Genomebits_get_char(ref,low_position);
  low2 = Genomebits_get_char(ref,low_position + 1);
  high2 = Genomebits_get_char(ref,high_position - 2);
  high1 = Genomebits_get_char(ref,high_position - 1);

  if (watsonp == true) {
    /* Check expected nucleotide pairs first */
    if (low1 == 'G' && high2 == 'A' && high1 == 'G') {
      if (low2 == 'T') {
	return 1;		/* GT-AG watson */
      } else if (low2 == 'C') {
	return 3;		/* GC-AG watson  */
      } else {
	return 0;		/* GX-AG watson */
      }

    } else if (low1 == 'A' && low2 == 'T' &&
	       high2 == 'A' && high1 == 'C') {
      return 5;			/* AT-AC watson */

    } else if (low1 == 'C' && low2 == 'T' && high1 == 'C') {
      if (high2 == 'A') {
	return 2;		/* GT-AG crick */
      } else if (high2 == 'G') {
	return 4;		/* GC-AG crick */
      } else {
	return 0;		/* GX-AG crick */
      }
    } else if (low1 == 'G' && low2 == 'T' &&
	       high2 == 'A' && high1 == 'T') {
      return 6;			/* AT-AC crick */

    } else {
      return 0;			/* Non-canonical */
    }

  } else {
    /* Check expected nucleotide pairs first */
    if (low1 == 'C' && low2 == 'T' && high1 == 'C') {
      if (high2 == 'A') {
	return 2;		/* GT-AG crick */
      } else if (high2 == 'G') {
	return 4;		/* GC-AG crick */
      } else {
	return 0;		/* GX-AG crick */
      }
    } else if (low1 == 'G' && low2 == 'T' &&
	       high2 == 'A' && high1 == 'T') {
      return 6;			/* AT-AC crick */

    } else if (low1 == 'G' && high2 == 'A' && high1 == 'G') {
      if (low2 == 'T') {
	return 1;		/* GT-AG watson */
      } else if (low2 == 'C') {
	return 3;		/* GC-AG watson  */
      } else {
	return 0;		/* GX-AG watson */
      }

    } else if (low1 == 'A' && low2 == 'T' &&
	       high2 == 'A' && high1 == 'C') {
      return 5;			/* AT-AC watson */

    } else {
      return 0;			/* Non-canonical */
    }
  }
}


void
Support_dump (FILE *fp, T this, EF64_T chromosome_ef64, Univ_IIT_T chromosome_iit,
	      Genomebits_T ref, bool include_known_p) {
  Chrnum_T chrnum;
  char *chr;
#if 0
  char donor1, donor2, acceptor2, acceptor1;
#endif

  Univcoord_T chroffset, chrhigh;
  Univcoord_T *donor_positions, *acceptor_positions, donor_position, acceptor_position;
  int nintrons, i;
  Splice_count_T splice_count;

  bool allocp = false;
  double donor_prob, acceptor_prob;


  /* plus */
  Univcoordpairtable_keys(&donor_positions,&acceptor_positions,this->plus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->plus_table);
  fprintf(stderr,"Have %d plus introns\n",nintrons);

  chrnum = 1;
  EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,/*chrnum*/1);
  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);

  for (i = 0; i < nintrons; i++) {
    donor_position = donor_positions[i];
    acceptor_position = acceptor_positions[i];

    if (donor_position >= chrhigh) {
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,donor_position,donor_position);
      if (allocp == true) {
	FREE(chr);
      }
      chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
    }

    splice_count = Univcoordpairtable_get(this->plus_table,donor_position,acceptor_position);
    if ((include_known_p == true && splice_count->knownp == true) ||
	(splice_count->nunique > 0 || splice_count->nmulti > 0)) {

      fprintf(fp,"%s\t",chr);
      fprintf(fp,"%u\t",donor_position - chroffset + INTRON_ADJ); /* low position */
      fprintf(fp,"%u\t",acceptor_position - chroffset + 1 - INTRON_ADJ); /* high position */

#ifdef DEBUG8
      fprintf(fp,"%u\t",donor_position);
      fprintf(fp,"%u\t",acceptor_position);
#endif      

      fprintf(fp,"+\t");	/* strand */
#if 0
      donor1 = Genomebits_get_char(ref,donor_position); /* (donor_position - 1) + 1 */
      donor2 = Genomebits_get_char(ref,donor_position+1);   /* (donor_position - 1) + 2 */
      acceptor2 = Genomebits_get_char(ref,acceptor_position-2);
      acceptor1 = Genomebits_get_char(ref,acceptor_position-1);
      fprintf(fp,"%c%c-%c%c\t",donor1,donor2,acceptor2,acceptor1); /* dinucleotides */
#else
      fprintf(fp,"%d\t",dinucleotide_coded(/*low*/donor_position,/*high*/acceptor_position,ref,/*watsonp*/true));
#endif

      if (splice_count->knownp == false) {
	fprintf(fp,"0\t");		/* annotatedp */
      } else {
	fprintf(fp,"1\t");		/* annotatedp */
      }

      fprintf(fp,"%d\t",splice_count->nunique);
      fprintf(fp,"%d\t",splice_count->nmulti);
      fprintf(fp,"%d\t",splice_count->support);

      donor_prob = Maxent_hr_donor_prob(donor_position,chroffset);
      acceptor_prob = Maxent_hr_acceptor_prob(acceptor_position,chroffset);
	      
      fprintf(fp,"%f\t",donor_prob);
      fprintf(fp,"%f\n",acceptor_prob);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  FREE(acceptor_positions);
  FREE(donor_positions);


  /* minus */
  Univcoordpairtable_keys(&acceptor_positions,&donor_positions,this->minus_table,/*sortp*/true,/*end*/0);
  nintrons = Univcoordpairtable_length(this->minus_table);
  fprintf(stderr,"Have %d minus introns\n",nintrons);

  chrnum = 1;
  EF64_chrbounds(&chroffset,&chrhigh,chromosome_ef64,/*chrnum*/1);
  chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);

  for (i = 0; i < nintrons; i++) {
    acceptor_position = acceptor_positions[i];
    donor_position = donor_positions[i];

    if (acceptor_position >= chrhigh) {
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,acceptor_position,acceptor_position);
      if (allocp == true) {
	FREE(chr);
      }
      chr = Univ_IIT_label(chromosome_iit,chrnum,&allocp);
    }

    splice_count = Univcoordpairtable_get(this->minus_table,acceptor_position,donor_position);
    if ((include_known_p == true && splice_count->knownp == true) ||
	(splice_count->nunique > 0 || splice_count->nmulti > 0)) {

      fprintf(fp,"%s\t",chr);
      fprintf(fp,"%u\t",acceptor_position - chroffset + INTRON_ADJ); /* low position */
      fprintf(fp,"%u\t",donor_position - chroffset + 1 - INTRON_ADJ); /* high position */
      
#ifdef DEBUG8
      fprintf(fp,"%u\t",acceptor_position);
      fprintf(fp,"%u\t",donor_position);
#endif      

      fprintf(fp,"-\t");	/* strand */

#if 0
      /* (acceptor_position - 1) + 1 */
      acceptor1 = complCode[(int) Genomebits_get_char(ref,acceptor_position)]; 
      /* (acceptor_position - 1) + 2 */
      acceptor2 = complCode[(int) Genomebits_get_char(ref,acceptor_position+1)];
      donor2 = complCode[(int) Genomebits_get_char(ref,donor_position-2)];
      donor1 = complCode[(int) Genomebits_get_char(ref,donor_position-1)];
      fprintf(fp,"%c%c-%c%c\t",donor1,donor2,acceptor2,acceptor1); /* dinucleotides */
#else
      fprintf(fp,"%d\t",dinucleotide_coded(/*low*/acceptor_position,/*high*/donor_position,ref,/*watsonp*/false));
#endif

      if (splice_count->knownp == false) {
	fprintf(fp,"0\t");		/* annotatedp */
      } else {
	fprintf(fp,"1\t");		/* annotatedp */
      }

      fprintf(fp,"%d\t",splice_count->nunique);
      fprintf(fp,"%d\t",splice_count->nmulti);
      fprintf(fp,"%d\t",splice_count->support);

      donor_prob = Maxent_hr_antidonor_prob(donor_position,chroffset);
      acceptor_prob = Maxent_hr_antiacceptor_prob(acceptor_position,chroffset);
	      
      fprintf(fp,"%f\t",acceptor_prob);
      fprintf(fp,"%f\n",donor_prob);
    }
  }

  if (allocp == true) {
    FREE(chr);
  }

  FREE(acceptor_positions);
  FREE(donor_positions);

  return;
}



