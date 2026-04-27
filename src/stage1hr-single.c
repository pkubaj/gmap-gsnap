static char rcsid[] = "$Id: 254bf8e98d592b24b209ae7013bf82a5182ef4bb $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifndef HAVE_MEMCPY
#define memcpy(d,s,n) bcopy((s),(d),(n))
#endif
#ifndef HAVE_MEMMOVE
#define memmove(d,s,n) bcopy((s),(d),(n))
#endif

#include "stage1hr.h"
#include "stage1hr-single.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>		/* For rint */
#include <string.h>		/* For memset */

#include "assert.h"
#include "mem.h"
#include "types.h"		/* Needed for HAVE_64_BIT */
#include "univcoord.h"

#include "list.h"
#include "compress.h"
#include "record.h"

#include "genomebits_mismatches.h" /* For MISMATCH_EXTRA */
#include "orderstat.h"
#include "transcriptome-search.h"
#include "tr-extension-search.h"
#include "kmer-search.h"
#include "extension-search.h"
/* #include "merge-search.h" */
/* #include "segment-search.h" */

#include "trpath-solve.h"
#include "trpath-convert.h"

#include "auxinfo.h"
#include "path.h"
#include "path-solve.h"
#include "path-fusion.h"

#include "transcript-remap.h"
#include "transcript-velocity.h"
#include "path-eval.h"


#define MIN_SIZELIMIT 100
#define MAX_HITS_EXACT 100	/* Excessive exact paths indicate a repetitive sequence */
#define MAX_END_POSITIONS 100	/* For END method.  Could be different for DNA-seq compared with RNA-seq */


static Mode_T mode;
static int index1part;
static int index1interval;
static int index1part_tr;

static Transcriptome_T transcriptome;
static bool transcriptome_align_p;
static bool genome_align_p;

static double user_nmismatches_filter_float;
static double user_mincoverage_filter_float;

static bool splicingp;


#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


#define T Stage1_T

#if 0
static int
determine_sizelimit (T this, int querylength) {
  int cutoff, *set, count;
  int n;
  int query_lastpos, querypos;

  assert(querylength >= index1part);

  query_lastpos = querylength - index1part;
  set = (int *) MALLOC(2*(query_lastpos+1)*sizeof(int));
  n = 0;
  for (querypos = 0; querypos <= query_lastpos; querypos++) {
    if (this->validp[querypos] == true) {
      set[n++] = count = this->plus_npositions[querypos];
      set[n++] = count = this->minus_npositions[querypos];
    }
  }

  if (n < 5) {
    cutoff = MIN_SIZELIMIT;
  } else if ((cutoff = Orderstat_int_pct_inplace(set,n,/*pct*/0.60)) < MIN_SIZELIMIT) {
    cutoff = MIN_SIZELIMIT;
  }
  FREE(set);

  return cutoff;
}
#endif


#ifdef DEBUG
static void
list_paths (List_T list, char *destination, bool expected_sensedir) {
  List_T p;
  Path_T path;

  for (p = list; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    printf("Destination %s: ",destination);
    Path_print(path);
    assert(path->sensedir == expected_sensedir);
  }

  return;
}
#endif


/* Benchmarks show slight speed improvement with TR_EXACT2 before TR_ANYPAIR */
Method_T
single_read_next_method_trdiagonal (Method_T last_method, T this, int querylength,
				    Compress_T query_compress_fwd, Compress_T query_compress_rev,
				    bool first_read_p) {

  if (last_method < TR_EXACT1) {
    /* 1.1.  TR_EXACT1 */
    debug(printf("%s Read: 1.1.  Running Transcriptome_exact1 for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_init_end_tr(this,querylength);
    Transcriptome_exact1(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			 &this->_sense_trdiagonals,&this->n_sense_trdiagonals,
			 &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			 &this->_antisense_trdiagonals,&this->n_antisense_trdiagonals,
			 this,querylength);
    return TR_EXACT1;

#if 0
  } else if (last_method < TR_EXACT2) {
    /* Can lead to ends that don't match the transcript */
    Stage1_trdiagonals_gc(this);
    
    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_init_end2_positions_tr(this,querylength);
    Transcriptome_exact2(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			 &this->_sense_trdiagonals,&this->n_sense_trdiagonals,
			 &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			 &this->_antisense_trdiagonals,&this->n_antisense_trdiagonals,
			 this,querylength);
    return TR_EXACT2;
#endif

#if 0
  } else if (last_method < TR_ANYPAIR) {
    /* Can lead to incorrect read loci */
    debug(printf("%s Read: 1.2.  Running Transcriptome_anypair for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_trdiagonals_gc(this);

    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_fill_all_positions_tr(this,querylength);

    Transcriptome_anypair(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			  &this->_sense_trdiagonals,&this->sense_tstarts,&this->sense_tends,
			  &this->n_sense_trdiagonals,
			  &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			  &this->_antisense_trdiagonals,&this->antisense_tstarts,&this->antisense_tends,
			  &this->n_antisense_trdiagonals,this,querylength);

    return TR_ANYPAIR;
#endif

  } else if (last_method < TR_PREVALENT) {
    debug(printf("%s Read: 1.3.  Running Transcriptome_prevalent for trdiagonals\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_trdiagonals_gc(this);

    Stage1_fill_all_oligos_tr(this,querylength);
    Stage1_fill_all_positions_tr(this,querylength);

    Transcriptome_prevalent(&this->sense_trnums,&this->sense_troffsets,&this->sense_trhighs,
			    &this->_sense_trdiagonals,&this->sense_tstarts,&this->sense_tends,
			    &this->n_sense_trdiagonals,
			    &this->antisense_trnums,&this->antisense_troffsets,&this->antisense_trhighs,
			    &this->_antisense_trdiagonals,&this->antisense_tstarts,&this->antisense_tends,
			    &this->n_antisense_trdiagonals,
			    this,querylength,query_compress_fwd,query_compress_rev);
    
    return TR_PREVALENT;

  } else {
    fprintf(stderr,"No trdiagonal method after TR_PREVALENT\n");
    abort();
  }
}


Method_T
single_read_next_method_tr (int *found_score, int sufficient_score, Method_T last_method,

			    List_T *partial_sense_trpaths, List_T *partial_antisense_trpaths,
			    List_T *sense_trpaths, List_T *antisense_trpaths,

			    T this, int genestrand, int querylength,
			    int *mismatch_positions_alloc,
			    Compress_T query_compress_fwd, Compress_T query_compress_rev,

			    int nmismatches_allowed,
			 
			    Trdiagpool_T trdiagpool,
			    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
			    Listpool_T listpool, Trpathpool_T trpathpool, Pathpool_T pathpool,
			    Hitlistpool_T hitlistpool, bool paired_end_p, bool first_read_p,
			    bool appendp) {

  if (appendp == false) {
    *partial_sense_trpaths = *partial_antisense_trpaths = (List_T) NULL;
    *sense_trpaths = *antisense_trpaths = (List_T) NULL;
  }

  if (last_method < TR_EXT) {

#if 0
    /* Supplanted by new Tr_extension_search,which also uses
       Genomebits_trim procedures, but still offers some speed
       advantage */
    /* 1.2.  TR_END */
    debug(printf("%s Read: 1.2.  Running Transcriptome_search_end\n",
		 first_read_p ? "5'" : "3'"));
    Transcriptome_search_end(&(*found_score),&(*sense_trpaths),&(*antisense_trpaths),

			     this,querylength,query_compress_fwd,query_compress_rev, 
			     intlistpool,uintlistpool,trpathpool,hitlistpool,
			     /*method*/TR_END);
    return TR_END;

  } else if (last_method == TR_END) {
#endif

    /* 1.3.  TR_EXT */
    debug(printf("%s Read: 1.4.  Running Tr_extension_search for trpaths\n",
		 first_read_p ? "5'" : "3'"));
    Stage1_fill_all_oligos_tr(this,querylength);
    Tr_extension_search(&(*found_score),sufficient_score,
			&(*partial_sense_trpaths),&(*partial_antisense_trpaths),
			&(*sense_trpaths),&(*antisense_trpaths),

			this,querylength,mismatch_positions_alloc,
			query_compress_fwd,query_compress_rev, 
			paired_end_p,first_read_p,

			trdiagpool,intlistpool,uintlistpool,listpool,
			trpathpool,pathpool,hitlistpool,
			nmismatches_allowed,genestrand,/*method*/TR_EXT);
    debug(printf(" => %d sense and %d antisense trpaths\n",
		 List_length(*sense_trpaths),List_length(*antisense_trpaths)));
    return TR_EXT;
    
  } else {
    fprintf(stderr,"No method after TR_EXT\n");
    abort();
  }
}


static void
single_read_extend (int *found_score, T this,

		    List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		    List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		    List_T *unextended_sense_paths_gplus, List_T *unextended_sense_paths_gminus,
		    List_T *unextended_antisense_paths_gplus, List_T *unextended_antisense_paths_gminus,

		    Shortread_T queryseq, char *queryuc_ptr, char *queryrc, int querylength,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		    int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		    Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		    Compress_T query_compress_fwd, Compress_T query_compress_rev, int nmismatches_allowed,

		    Intlistpool_T intlistpool,
		    Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Univdiagpool_T univdiagpool, Vectorpool_T vectorpool,
		    Hitlistpool_T hitlistpool, Spliceendsgen_T spliceendsgen) {

  List_T filtered_paths, unextended_paths, complete_paths, p;
  Path_T path;

  /* Attempt to extend paths */
  debug(printf("Entered single_read_extend\n"));

#if 0
  /* Path_filter appears to select incorrect splice lengths */
  debug(printf("Before Path_filter: %d unextended_sense_paths_gplus\n",
	       List_length(this->unextended_sense_paths_gplus)));
  debug(printf("Before Path_filter: %d unextended_antisense_paths_gplus\n",
	       List_length(this->unextended_antisense_paths_gplus)));
  debug(printf("Before Path_filter: %d unextended_sense_paths_gminus\n",
	       List_length(this->unextended_sense_paths_gminus)));
  debug(printf("Before Path_filter: %d unextended_antisense_paths_gminus\n",
	       List_length(this->unextended_antisense_paths_gminus)));

  this->unextended_sense_paths_gplus = Path_filter(this->unextended_sense_paths_gplus,
						   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_sense_paths_gminus = Path_filter(this->unextended_sense_paths_gminus,
						    intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_antisense_paths_gplus = Path_filter(this->unextended_antisense_paths_gplus,
						       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  this->unextended_antisense_paths_gminus = Path_filter(this->unextended_antisense_paths_gminus,
							intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#elif 0
  /* Path_unique gets rid of sense/antisense, needed for finding fusions */
  *unextended_sense_paths_gplus = Path_unique(*unextended_sense_paths_gplus,
					      intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_sense_paths_gminus = Path_unique(*unextended_sense_paths_gminus,
					       intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_antisense_paths_gplus = Path_unique(*unextended_antisense_paths_gplus,
						  intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  *unextended_antisense_paths_gminus = Path_unique(*unextended_antisense_paths_gminus,
						   intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
#endif
      
  debug(printf("Have %d unextended_sense_paths_gplus\n",List_length(*unextended_sense_paths_gplus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryuc_ptr,querylength,
			    mismatch_positions_alloc,novel_univdiags_alloc,
			    novel_diagonals_alloc,localdb_alloc,
			    this,knownsplicing,knownindels,
			    /*query_compress*/query_compress_fwd,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    nmismatches_allowed,/*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_sense_paths_gplus\n"));
      path->completep = true;
      *sense_paths_gplus = List_append(complete_paths,*sense_paths_gplus);

    } else {
      /* TODO: Path_gc(&unextended_paths); */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_sense_paths_gplus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_sense_paths_gplus = filtered_paths;


  debug(printf("Have %d unextended_sense_paths_gminus\n",List_length(*unextended_sense_paths_gminus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths =
		Path_extend(&(*found_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryrc,querylength,
			    mismatch_positions_alloc,novel_univdiags_alloc,
			    novel_diagonals_alloc,localdb_alloc,
			    this,knownsplicing,knownindels,
			    /*query_compress*/query_compress_rev,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    nmismatches_allowed,/*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_sense_paths_gminus\n"));
      path->completep = true;
      *sense_paths_gminus = List_append(complete_paths,*sense_paths_gminus);

    } else {
      /* TODO: Path_gc(&unextended_paths); */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_sense_paths_gminus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_sense_paths_gminus = filtered_paths;


  debug(printf("Have %d unextended_antisense_paths_gplus\n",List_length(*unextended_antisense_paths_gplus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryuc_ptr,querylength,
			    mismatch_positions_alloc,novel_univdiags_alloc,
			    novel_diagonals_alloc,localdb_alloc,
			    this,knownsplicing,knownindels,
			    /*query_compress*/query_compress_fwd,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    nmismatches_allowed,/*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_antisense_paths_gplus\n"));
      path->completep = true;
      *antisense_paths_gplus = List_append(complete_paths,*antisense_paths_gplus);
      
    } else {
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_antisense_paths_gplus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_antisense_paths_gplus = filtered_paths;


  debug(printf("Have %d unextended_antisense_paths_gminus\n",List_length(*unextended_antisense_paths_gminus)));
  filtered_paths = (List_T) NULL;
  for (p = *unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);
    if (path->extendedp == true) {
      /* Skip */
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    } else if ((complete_paths = 
		Path_extend(&(*found_score),&unextended_paths,/*original_path*/path,
			    queryseq,/*queryptr*/queryrc,querylength,
			    mismatch_positions_alloc,novel_univdiags_alloc,
			    novel_diagonals_alloc,localdb_alloc,
			    this,knownsplicing,knownindels,
			    /*query_compress*/query_compress_rev,
			    query_compress_fwd,query_compress_rev,/*genestrand*/0,
			    nmismatches_allowed,/*paired_end_p*/false,/*lowp*/true,
			    intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
			    univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			    /*extend_qstart_p*/true,/*extend_qend_p*/true)) != NULL) {
      debug(printf("Found extended_antisense_paths_gminus\n"));
      path->completep = true;
      *antisense_paths_gminus = List_append(complete_paths,*antisense_paths_gminus);
    } else {
      filtered_paths = Hitlist_push(filtered_paths,hitlistpool,(void *) path
				    hitlistpool_trace(__FILE__,__LINE__));
    }
  }
  Hitlistpool_free_list(&(*unextended_antisense_paths_gminus),hitlistpool
			hitlistpool_trace(__FILE__,__LINE__));
  *unextended_antisense_paths_gminus = filtered_paths;

  return;
}



static List_T
single_read_fusion (int *found_score, T this, int querylength,

		    List_T unextended_sense_paths_gplus, List_T unextended_sense_paths_gminus,
		    List_T unextended_antisense_paths_gplus, List_T unextended_antisense_paths_gminus,

		    Compress_T query_compress_fwd, Compress_T query_compress_rev,
		    Shortread_T queryseq, Knownsplicing_T knownsplicing, int nmismatches_allowed,

		    Univdiagpool_T univdiagpool, Intlistpool_T intlistpool,
		    Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		    Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Vectorpool_T vectorpool, Hitlistpool_T hitlistpool) {

  List_T paths = NULL, p;
  Path_T path;

  /* Look for possible fusions (or combinations of existing paths) */
  /* Find fusions.  Use code similar to finding outer fusions in Pathpair_eval_and_sort */

  debug(printf("Have %d unextended_sense_paths_gplus\n",List_length(unextended_sense_paths_gplus)));
  for (p = unextended_sense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_plus(&(*found_score),/*main*/path,this,

						      query_compress_fwd,query_compress_rev,
						      queryseq,querylength,knownsplicing,/*genestrand*/0,
						      nmismatches_allowed,
						      intlistpool,uintlistpool,univcoordlistpool,
						      listpool,univdiagpool,pathpool,vectorpool,
						      transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_plus(&(*found_score),/*main*/path,this,

						    query_compress_fwd,query_compress_rev,
						    queryseq,querylength,knownsplicing,/*genestrand*/0,
						    nmismatches_allowed,
						    intlistpool,uintlistpool,univcoordlistpool,
						    listpool,univdiagpool,pathpool,vectorpool,
						    transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_antisense_paths_gplus\n",List_length(unextended_antisense_paths_gplus)));
  for (p = unextended_antisense_paths_gplus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_plus(&(*found_score),/*main*/path,this,

						      query_compress_fwd,query_compress_rev,
						      queryseq,querylength,knownsplicing,/*genestrand*/0,
						      nmismatches_allowed,
						      intlistpool,uintlistpool,univcoordlistpool,
						      listpool,univdiagpool,pathpool,vectorpool,
						      transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_plus(&(*found_score),/*main*/path,this,

						    query_compress_fwd,query_compress_rev,
						    queryseq,querylength,knownsplicing,/*genestrand*/0,
						    nmismatches_allowed,
						    intlistpool,uintlistpool,univcoordlistpool,
						    listpool,univdiagpool,pathpool,vectorpool,
						    transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_sense_paths_gminus\n",List_length(unextended_sense_paths_gminus)));
  for (p = unextended_sense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_minus(&(*found_score),/*main*/path,this,

						       query_compress_fwd,query_compress_rev,
						       queryseq,querylength,knownsplicing,/*genestrand*/0,
						       nmismatches_allowed,
						       intlistpool,uintlistpool,univcoordlistpool,
						       listpool,univdiagpool,pathpool,vectorpool,
						       transcriptpool,hitlistpool),paths);
    }

    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_minus(&(*found_score),/*main*/path,this,

						     query_compress_fwd,query_compress_rev,
						     queryseq,querylength,knownsplicing,/*genestrand*/0,
						     nmismatches_allowed,
						     intlistpool,uintlistpool,univcoordlistpool,
						     listpool,univdiagpool,pathpool,vectorpool,
						     transcriptpool,hitlistpool),paths);
    }
  }
    
  debug(printf("Have %d unextended_antisense_paths_gminus\n",List_length(unextended_antisense_paths_gminus)));
  for (p = unextended_antisense_paths_gminus; p != NULL; p = List_next(p)) {
    path = (Path_T) List_head(p);

    if (Path_unextended_qend_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qstart_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_querystart_minus(&(*found_score),/*main*/path,this,

						       query_compress_fwd,query_compress_rev,
						       queryseq,querylength,knownsplicing,/*genestrand*/0,
						       nmismatches_allowed,
						       intlistpool,uintlistpool,univcoordlistpool,
						       listpool,univdiagpool,pathpool,vectorpool,
						       transcriptpool,hitlistpool),paths);
    }
      
    if (Path_unextended_qstart_p(path,/*endtrim_allowed*/8,/*allow_ambig_p*/false) == false &&
	Path_unextended_qend_p(path,/*endtrim_allowed*/25,/*allow_ambig_p*/false) == true) {
      paths = List_append(Path_fusion_queryend_minus(&(*found_score),/*main*/path,this,

						     query_compress_fwd,query_compress_rev,
						     queryseq,querylength,knownsplicing,/*genestrand*/0,
						     nmismatches_allowed,
						     intlistpool,uintlistpool,univcoordlistpool,
						     listpool,univdiagpool,pathpool,vectorpool,
						     transcriptpool,hitlistpool),paths);
    }
  }

  return paths;
}


static void
convert_trpaths (int *found_score,

		 List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		 List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		 List_T sense_trpaths, List_T antisense_trpaths, bool first_read_p,

		 T this, Knownsplicing_T knownsplicing,
		 Shortread_T queryseq, int querylength,
		 Compress_T query_compress_fwd, Compress_T query_compress_rev,

		 Uintlistpool_T uintlistpool, Intlistpool_T intlistpool, Univcoordlistpool_T univcoordlistpool,
		 Listpool_T listpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		 Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  *sense_paths_gplus = *sense_paths_gminus = (List_T) NULL;
  Trpath_convert_sense(&(*found_score),
		       &(*sense_paths_gplus),&(*sense_paths_gminus),
		       /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,
		       sense_trpaths,first_read_p,
		       queryseq,querylength,this,knownsplicing,
		       query_compress_fwd,query_compress_rev,
		       intlistpool,uintlistpool,univcoordlistpool,listpool,
		       pathpool,vectorpool,transcriptpool,hitlistpool,
		       /*crossover_sense_p*/false);

  *sense_paths_gplus = *sense_paths_gminus = (List_T) NULL;
  Trpath_convert_antisense(&(*found_score),
			   &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			   /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			   antisense_trpaths,first_read_p,
			   queryseq,querylength,this,knownsplicing,
			   query_compress_fwd,query_compress_rev,
			   intlistpool,uintlistpool,univcoordlistpool,listpool,
			   pathpool,vectorpool,transcriptpool,hitlistpool,
			   /*crossover_sense_p*/false);

  *sense_paths_gplus = Path_consolidate(*sense_paths_gplus,queryseq,
					query_compress_fwd,query_compress_rev,
					uintlistpool,intlistpool,univcoordlistpool,
					listpool,pathpool,transcriptpool,hitlistpool);
  *sense_paths_gminus = Path_consolidate(*sense_paths_gminus,queryseq,
					query_compress_fwd,query_compress_rev,
					uintlistpool,intlistpool,univcoordlistpool,
					listpool,pathpool,transcriptpool,hitlistpool);
  *antisense_paths_gplus = Path_consolidate(*antisense_paths_gplus,queryseq,
					   query_compress_fwd,query_compress_rev,
					   uintlistpool,intlistpool,univcoordlistpool,
					   listpool,pathpool,transcriptpool,hitlistpool);
  *antisense_paths_gminus = Path_consolidate(*antisense_paths_gminus,queryseq,
					    query_compress_fwd,query_compress_rev,
					    uintlistpool,intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,hitlistpool);

  return;
}
    

static bool
single_read_tr_paths (int *found_score, int sufficient_score, Method_T *last_method,

		      List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		      List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		      T this, int genestrand,
		      
		      Shortread_T queryseq, int querylength,
		      Knownsplicing_T knownsplicing, int *mismatch_positions_alloc, 
		      Compress_T query_compress_fwd, Compress_T query_compress_rev,

		      int nmismatches_allowed,
			 
		      Trdiagpool_T trdiagpool, Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		      Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		      Trpathpool_T trpathpool, Pathpool_T pathpool, Vectorpool_T vectorpool,
		      Transcriptpool_T transcriptpool, Hitlistpool_T hitlistpool) {

  Trpath_T trpath;
  int tstart, tend;
  int index;

  Trnum_T trnum;
  Chrnum_T chrnum;
  int transcript_genestrand;
  bool partialp;


  /* Univcoord_T *univdiagonals_gplus, *univdiagonals_gminus; */
  /* int nunivdiagonals_gplus, nunivdiagonals_gminus, i; */


  *sense_paths_gplus = NULL;
  *sense_paths_gminus = NULL;
  *antisense_paths_gplus = NULL;
  *antisense_paths_gminus = NULL;
    
  if (transcriptome_align_p == true) {
    /* A.  Transcriptome search using univdiagonals */
    while (*found_score > sufficient_score && *last_method < TR_PREVALENT) {
      /* Append results to lists in Stage1_T object */
      *last_method = single_read_next_method_trdiagonal(*last_method,this,querylength,
							query_compress_fwd,query_compress_rev,
							/*first_read_p*/true);
      tstart = 0;
      tend = querylength;
      for (index = 0; index < this->n_sense_trdiagonals; index++) {
	if (this->sense_tstarts != NULL) {
	  tstart = this->sense_tstarts[index];
	  tend = this->sense_tends[index];
	}

	trnum = this->sense_trnums[index];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&partialp,sufficient_score,
						   /*trdiagonal*/this->_sense_trdiagonals[index],
						   tstart,tend,trnum,/*troffset*/this->sense_troffsets[index],
						   /*trhigh*/this->sense_trhighs[index],
						   /*query_compress_tr*/query_compress_fwd,/*tplusp*/true,querylength,
						   mismatch_positions_alloc,chrnum,transcript_genestrand,
						   /*paired_end_p*/false,/*first_read_p*/true,this->indelinfo,
						   intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						   *last_method)) != NULL) {
	  if (partialp == true) {
	    this->partial_sense_trpaths = Hitlist_push(this->partial_sense_trpaths,hitlistpool,(void *) trpath
					       hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this->sense_trpaths = Hitlist_push(this->sense_trpaths,hitlistpool,(void *) trpath
					       hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
      }

      tstart = 0;
      tend = querylength;
      for (index = 0; index < this->n_antisense_trdiagonals; index++) {
	if (this->antisense_tstarts != NULL) {
	  tstart = this->antisense_tstarts[index];
	  tend = this->antisense_tends[index];
	}

	trnum = this->antisense_trnums[index];
	chrnum = Transcriptome_chrnum(&transcript_genestrand,transcriptome,trnum);
	if ((trpath = Trpath_solve_from_trdiagonal(&(*found_score),&partialp,sufficient_score,
						   /*trdiagonal*/this->_antisense_trdiagonals[index],
						   tstart,tend,trnum,/*troffset*/this->antisense_troffsets[index],
						   /*trhigh*/this->antisense_trhighs[index],
						   /*query_compress_tr*/query_compress_rev,/*tplusp*/false,querylength,
						   mismatch_positions_alloc,chrnum,transcript_genestrand,
						   /*paired_end_p*/false,/*first_read_p*/true,this->indelinfo,
						   intlistpool,uintlistpool,listpool,trpathpool,pathpool,
						   *last_method)) != NULL) {
	  if (partialp == true) {
	    this->partial_antisense_trpaths = Hitlist_push(this->partial_antisense_trpaths,hitlistpool,(void *) trpath
							   hitlistpool_trace(__FILE__,__LINE__));
	  } else {
	    this->antisense_trpaths = Hitlist_push(this->antisense_trpaths,hitlistpool,(void *) trpath
						   hitlistpool_trace(__FILE__,__LINE__));
	  }
	}
      }
    }

    /* Should be handled by Stage1hr_free */
    /* FREE_ALIGN(this->_sense_trdiagonals); */
    /* FREE_ALIGN(this->_antisense_trdiagonals); */

    if (*found_score <= sufficient_score) {
      convert_trpaths(&(*found_score),
		      
		      &(*sense_paths_gplus),&(*sense_paths_gminus),
		      &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		      this->sense_trpaths,this->antisense_trpaths,/*first_read_p*/true,
		      this,knownsplicing,queryseq,querylength,
		      query_compress_fwd,query_compress_rev,
		      uintlistpool,intlistpool,univcoordlistpool,listpool,
		      pathpool,vectorpool,transcriptpool,hitlistpool);
      return true;
    }



    /* B.  Transcriptome search using trpaths */
    while (*found_score > sufficient_score && *last_method < TR_EXT) {
      /* Append results to lists in Stage1_T object */
      *last_method = single_read_next_method_tr(&(*found_score),sufficient_score,*last_method,
					     
						&this->partial_sense_trpaths,&this->partial_antisense_trpaths,
						&this->sense_trpaths,&this->antisense_trpaths,
					     
						this,genestrand,querylength,
						mismatch_positions_alloc,
						query_compress_fwd,query_compress_rev,
						nmismatches_allowed,
						trdiagpool,intlistpool,uintlistpool,
						listpool,trpathpool,pathpool,hitlistpool,
						/*paired_end_p*/false,/*first_read_p*/true,/*appendp*/true);
    }

    if (*found_score <= sufficient_score) {
      convert_trpaths(&(*found_score),

		      &(*sense_paths_gplus),&(*sense_paths_gminus),
		      &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		      this->sense_trpaths,this->antisense_trpaths,/*first_read_p*/true,
		      this,knownsplicing,queryseq,querylength,
		      query_compress_fwd,query_compress_rev,
		      uintlistpool,intlistpool,univcoordlistpool,listpool,
		      pathpool,vectorpool,transcriptpool,hitlistpool);
      return true;
    }
    

    /* C.  Prep for genome search.  Convert all partial and complete trpaths to paths */
    *sense_paths_gplus = *sense_paths_gminus = (List_T) NULL;
    Trpath_convert_sense(&(*found_score),
			 &(*sense_paths_gplus),&(*sense_paths_gminus),
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 this->partial_sense_trpaths,/*first_read_p*/true,
			 queryseq,querylength,this,knownsplicing,
			 query_compress_fwd,query_compress_rev,
			 intlistpool,uintlistpool,univcoordlistpool,listpool,
			 pathpool,vectorpool,transcriptpool,hitlistpool,
			 /*crossover_sense_p*/false);

    Trpath_convert_sense(&(*found_score),
			 &(*sense_paths_gplus),&(*sense_paths_gminus),
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 this->sense_trpaths,/*first_read_p*/true,
			 queryseq,querylength,this,knownsplicing,
			 query_compress_fwd,query_compress_rev,
			 intlistpool,uintlistpool,univcoordlistpool,listpool,
			 pathpool,vectorpool,transcriptpool,hitlistpool,
			 /*crossover_sense_p*/false);

    *antisense_paths_gplus = *antisense_paths_gminus = (List_T) NULL;
    Trpath_convert_antisense(&(*found_score),
			     &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     this->partial_antisense_trpaths,/*first_read_p*/true,
			     queryseq,querylength,this,knownsplicing,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool,
			     /*crossover_sense_p*/false);

    Trpath_convert_antisense(&(*found_score),
			     &(*antisense_paths_gplus),&(*antisense_paths_gminus),
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     this->antisense_trpaths,/*first_read_p*/true,
			     queryseq,querylength,this,knownsplicing,
			     query_compress_fwd,query_compress_rev,
			     intlistpool,uintlistpool,univcoordlistpool,listpool,
			     pathpool,vectorpool,transcriptpool,hitlistpool,
			     /*crossover_sense_p*/false);


    if (*found_score <= sufficient_score) {
      *sense_paths_gplus = Path_consolidate(*sense_paths_gplus,queryseq,
					    query_compress_fwd,query_compress_rev,
					    uintlistpool,intlistpool,univcoordlistpool,
					    listpool,pathpool,transcriptpool,hitlistpool);
      *sense_paths_gminus = Path_consolidate(*sense_paths_gminus,queryseq,
					     query_compress_fwd,query_compress_rev,
					     uintlistpool,intlistpool,univcoordlistpool,
					     listpool,pathpool,transcriptpool,hitlistpool);
      *antisense_paths_gplus = Path_consolidate(*antisense_paths_gplus,queryseq,
						query_compress_fwd,query_compress_rev,
						uintlistpool,intlistpool,univcoordlistpool,
						listpool,pathpool,transcriptpool,hitlistpool);
      *antisense_paths_gminus = Path_consolidate(*antisense_paths_gminus,queryseq,
						 query_compress_fwd,query_compress_rev,
						 uintlistpool,intlistpool,univcoordlistpool,
						 listpool,pathpool,transcriptpool,hitlistpool);
      return true;
    }

    /* Previously solved unsolved paths here, but now Trpath_convert procedures solve them */
  }

  return false;
}


/* For single-end reads, we will just call
   Path_solve_from_univdiagonals on most prevalent univdiagonals.
   For paired-end reads, we will return all univdiagonals (and most
   prevalent) and try to anchor */

/* TODO: Need to return Auxinfo_T arrays, which have the paths */
static bool
single_read_gen_paths (int *found_score, Method_T *last_method,

		       List_T *unextended_sense_paths_gplus, List_T *unextended_sense_paths_gminus,
		       List_T *unextended_antisense_paths_gplus, List_T *unextended_antisense_paths_gminus,

		       List_T *sense_paths_gplus, List_T *sense_paths_gminus,
		       List_T *antisense_paths_gplus, List_T *antisense_paths_gminus,

		       T this, int genestrand,
		 
		       Shortread_T queryseq, char *queryuc_ptr, char *queryrc, int querylength,
		       Knownsplicing_T knownsplicing, Knownindels_T knownindels,
		       int *mismatch_positions_alloc, Univdiag_T *novel_univdiags_alloc,
		       Univcoord_T *novel_diagonals_alloc, unsigned short *localdb_alloc,
		       Compress_T query_compress_fwd, Compress_T query_compress_rev,
		 
		       int localdb_nmismatches_allowed, EF64_T repetitive_ef64,
		       
		       Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		       Univdiagpool_T univdiagpool, 
		       Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		       Univcoordlistpool_T univcoordlistpool, Listpool_T listpool, 
		       Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		       Vectorpool_T vectorpool, Hitlistpool_T hitlistpool,
		       Spliceendsgen_T spliceendsgen) {
  
  int sufficient_score = querylength/20;
  bool any_imperfect_ends_p = false;
  int kmer_querystart, kmer_queryend;
  int i;

  int nunivdiagonals_gplus, nunivdiagonals_gminus;
  int total_npositions_plus, total_npositions_minus;
  Univcoord_T *_univdiagonals_gplus, *_univdiagonals_gminus, univdiagonal;
  Auxinfo_T *auxinfo_gplus = NULL, *auxinfo_gminus = NULL, auxinfo;
  Pathstore_T *pathstores_gplus = NULL, *pathstores_gminus = NULL, pathstore;

#ifdef INDIVIDUAL_CHRINFO
  Chrnum_T chrnum;
  Univcoord_T chroffset, chrhigh;
#endif

  
  *sense_paths_gplus = NULL;
  *sense_paths_gminus = NULL;
  *antisense_paths_gplus = NULL;
  *antisense_paths_gminus = NULL;


  /* 1. Exact */
  debug(printf("Single Read: 1.  Running Kmer_exact1\n"));
  Stage1_init_end_gen(&kmer_querystart,&kmer_queryend,this,querylength,genestrand);
  Kmer_exact1(&_univdiagonals_gplus,&auxinfo_gplus,
	      &pathstores_gplus,&nunivdiagonals_gplus,
	      &_univdiagonals_gminus,&auxinfo_gminus,
	      &pathstores_gminus,&nunivdiagonals_gminus,
	      this,kmer_querystart,kmer_queryend,querylength,
	      auxinfopool,pathstorepool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf(" => %d plus and %d minus univdiagonals\n",nunivdiagonals_gplus,nunivdiagonals_gminus));
  debug(printf("Kmer_exact1 returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = KMER_EXACT1;
  any_imperfect_ends_p = false;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

    if (univdiagonal < (Univcoord_T) querylength) {
      /* Skip */
    } else {
#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength,univdiagonal);
#endif
      if (Path_solve_exact(&(*found_score),
			   
			   &(*sense_paths_gplus),&(*antisense_paths_gplus),

			   univdiagonal,auxinfo,querylength,
			   /*plusp*/true,/*first_read_p*/true,genestrand,
			   /*query_compress*/query_compress_fwd,
			   query_compress_fwd,query_compress_rev,
			   queryseq,queryuc_ptr,queryrc,
			   /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			   /*chrhigh*/pathstore->chrhigh,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			   /*method*/KMER_EXACT1) == false) {
	debug(printf("Imperfect end\n"));
	any_imperfect_ends_p = true;
      }
    }
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

    if (univdiagonal < (Univcoord_T) querylength) {
      /* Skip */
    } else {
#ifdef INDIVIDUAL_CHRINFO
      chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			   univdiagonal - querylength,univdiagonal);
#endif
      if (Path_solve_exact(&(*found_score),

			   &(*sense_paths_gminus),&(*antisense_paths_gminus),

			   univdiagonal,auxinfo,querylength,
			   /*plusp*/false,/*first_read_p*/true,genestrand,
			   /*query_compress*/query_compress_rev,
			   query_compress_fwd,query_compress_rev,
			   queryseq,queryuc_ptr,queryrc,
			   /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			   /*chrhigh*/pathstore->chrhigh,
			   intlistpool,uintlistpool,univcoordlistpool,
			   listpool,pathpool,vectorpool,hitlistpool,transcriptpool,
			   /*method*/KMER_EXACT1) == false) {
	debug(printf("Imperfect end\n"));
	any_imperfect_ends_p = true;
      }
    }
  }
    
  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);


  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    return true;
  }


  /* 2. Extension search */
  debug(printf("Single Read: 2.  Running Extension search\n"));
  Stage1_fill_all_oligos_gen(this,querylength,genestrand);
  Extension_search(&_univdiagonals_gplus,&auxinfo_gplus,
		   &pathstores_gplus,&nunivdiagonals_gplus,
		   &_univdiagonals_gminus,&auxinfo_gminus,
		   &pathstores_gminus,&nunivdiagonals_gminus,

		   this,query_compress_fwd,query_compress_rev,querylength,
		   auxinfopool,pathstorepool,univdiagpool,univcoordlistpool,
		   listpool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf("Extension_search returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = EXT;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

#ifdef INDIVIDUAL_CHRINFO
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonal - querylength + auxinfo->qstart,
			 univdiagonal - querylength + auxinfo->qend);
#endif
    Path_solve_from_diagonals(&(*found_score),

			      &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
			      &(*sense_paths_gplus),&(*antisense_paths_gplus),

			      univdiagonal,auxinfo,queryseq,/*queryptr*/queryuc_ptr,querylength,
			      mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
			      localdb_alloc,/*stage1*/this,knownsplicing,knownindels,
			      /*query_compress*/query_compress_fwd,query_compress_fwd,query_compress_rev,
			      /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			      /*chrhigh*/pathstore->chrhigh,
			      /*plusp*/true,genestrand,
			      localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,
			      intlistpool,uintlistpool,univcoordlistpool,listpool,
			      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			      /*method*/EXT,/*find_splices_p*/true);
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

#ifdef INDIVIDUAL_CHRINFO
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonal - querylength + auxinfo->qstart,
			 univdiagonal - querylength + auxinfo->qend);
#endif
    Path_solve_from_diagonals(&(*found_score),

			      &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
			      &(*sense_paths_gminus),&(*antisense_paths_gminus),

			      univdiagonal,auxinfo,queryseq,/*queryptr*/queryrc,querylength,
			      mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
			      localdb_alloc,/*stage1*/this,knownsplicing,knownindels,
			      /*query_compress*/query_compress_rev,query_compress_fwd,query_compress_rev,
			      /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			      /*chrhigh*/pathstore->chrhigh,
			      /*plusp*/false,genestrand,
			      localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,
			      intlistpool,uintlistpool,univcoordlistpool,listpool,
			      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			      /*method*/EXT,/*find_splices_p*/true);
  }

  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);

  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    return true;
  }


  /* 3. Segment search */
  debug(printf("Single Read: 3.  Running Kmer segment\n"));
  assert(this->all_oligos_gen_filledp == true);
  Stage1_fill_all_positions_gen(&total_npositions_plus,&total_npositions_minus,
				this,querylength,genestrand);
  Kmer_segment(&_univdiagonals_gplus,&auxinfo_gplus,&pathstores_gplus,&nunivdiagonals_gplus,
	       &_univdiagonals_gminus,&auxinfo_gminus,&pathstores_gminus,&nunivdiagonals_gminus,
	       this,querylength,repetitive_ef64,univdiagpool,auxinfopool,pathstorepool);

  Pathstore_assign_chrinfo(_univdiagonals_gplus,pathstores_gplus,nunivdiagonals_gplus,querylength);
  Pathstore_assign_chrinfo(_univdiagonals_gminus,pathstores_gminus,nunivdiagonals_gminus,querylength);

  debug(printf("Kmer_segment returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = SEGMENT1;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

#ifdef INDIVIDUAL_CHRINFO
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonal - querylength + auxinfo->qstart,
			 univdiagonal - querylength + auxinfo->qend);
#endif
    Path_solve_from_diagonals(&(*found_score),

			      &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
			      &(*sense_paths_gplus),&(*antisense_paths_gplus),

			      univdiagonal,auxinfo,queryseq,/*queryptr*/queryuc_ptr,querylength,
			      mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
			      localdb_alloc,/*stage1*/this,knownsplicing,knownindels,
			      /*query_compress*/query_compress_fwd,query_compress_fwd,query_compress_rev,
			      /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			      /*chrhigh*/pathstore->chrhigh,
			      /*plusp*/true,genestrand,
			      localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,
			      intlistpool,uintlistpool,univcoordlistpool,listpool,
			      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			      /*method*/SEGMENT1,/*find_splices_p*/true);
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

#ifdef INDIVIDUAL_CHRINFO
    chrnum = EF64_chrnum(&chroffset,&chrhigh,chromosome_ef64,
			 univdiagonal - querylength + auxinfo->qstart,
			 univdiagonal - querylength + auxinfo->qend);
#endif

    Path_solve_from_diagonals(&(*found_score),

			      &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
			      &(*sense_paths_gminus),&(*antisense_paths_gminus),

			      univdiagonal,auxinfo,queryseq,/*queryptr*/queryrc,querylength,
			      mismatch_positions_alloc,novel_univdiags_alloc,novel_diagonals_alloc,
			      localdb_alloc,/*stage1*/this,knownsplicing,knownindels,
			      /*query_compress*/query_compress_rev,query_compress_fwd,query_compress_rev,
			      /*chrnum*/pathstore->chrnum,/*chroffset*/pathstore->chroffset,
			      /*chrhigh*/pathstore->chrhigh,
			      /*plusp*/false,genestrand,
			      localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,
			      intlistpool,uintlistpool,univcoordlistpool,listpool,
			      pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
			      /*method*/SEGMENT1,/*find_splices_p*/true);
  }

  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);

#if 0
  /* 3. Prevalent (merging).  Equivalent of segment search */
  debug(printf("Single Read: 3.  Running Kmer prevalent\n"));
  assert(this->all_oligos_gen_filledp == true); /* From Extension_search */
  Stage1_fill_all_positions_gen(this,querylength,genestrand);
  Kmer_prevalent(&_univdiagonals_gplus,&auxinfo_gplus,&nunivdiagonals_gplus,
		 &_univdiagonals_gminus,&auxinfo_gminus,&nunivdiagonals_gminus,
		 this,querylength);
  debug(printf(" => %d plus and %d minus univdiagonals\n",nunivdiagonals_gplus,nunivdiagonals_gminus));
  debug(printf("Kmer_prevalent returning %d plus and %d minus univdiagonals\n",
	       nunivdiagonals_gplus,nunivdiagonals_gminus));
  *last_method = KMER_PREVALENT;

  for (i = 0; i < nunivdiagonals_gplus; i++) {
    univdiagonal = _univdiagonals_gplus[i];
    auxinfo = auxinfo_gplus[i];
    pathstore = pathstores_gplus[i];

    Path_solve_from_univdiagonal(&(*found_score),

				 &(*unextended_sense_paths_gplus),&(*unextended_antisense_paths_gplus),
				 &(*sense_paths_gplus),&(*antisense_paths_gplus),

				 univdiagonal,auxinfo,queryseq,/*queryptr*/queryuc_ptr,
				 /*query_compress*/query_compress_fwd,
				 query_compress_fwd,query_compress_rev,
				 /*plusp*/true,querylength,mismatch_positions_alloc,
				   
				 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				 this,knownsplicing,knownindels,
				 localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,

				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				 /*method*/KMER_PREVALENT,/*find_splices_p*/true);
  }
  
  for (i = 0; i < nunivdiagonals_gminus; i++) {
    univdiagonal = _univdiagonals_gminus[i];
    auxinfo = auxinfo_gminus[i];
    pathstore = pathstores_gminus[i];

    Path_solve_from_univdiagonal(&(*found_score),

				 &(*unextended_sense_paths_gminus),&(*unextended_antisense_paths_gminus),
				 &(*sense_paths_gminus),&(*antisense_paths_gminus),

				 univdiagonal,auxinfo,queryseq,/*queryptr*/queryrc,
				 /*query_compress*/query_compress_rev,
				 query_compress_fwd,query_compress_rev,
				 /*plusp*/false,querylength,mismatch_positions_alloc,
				 
				 novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				 this,knownsplicing,knownindels,
				 localdb_nmismatches_allowed,/*paired_end_p*/false,/*first_read_p*/true,
				 
				 intlistpool,uintlistpool,univcoordlistpool,listpool,
				 pathpool,transcriptpool,univdiagpool,vectorpool,hitlistpool,spliceendsgen,
				 /*method*/KMER_PREVALENT,/*find_splices_p*/true);
  }

  FREE_ALIGN(_univdiagonals_gplus);
  FREE_ALIGN(_univdiagonals_gminus);
  Auxinfo_gc(auxinfo_gplus,nunivdiagonals_gplus,univdiagpool,auxinfopool);
  Auxinfo_gc(auxinfo_gminus,nunivdiagonals_gminus,univdiagpool,auxinfopool);
  Pathstore_gc(pathstores_gplus,nunivdiagonals_gplus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
  Pathstore_gc(pathstores_gminus,nunivdiagonals_gminus,pathstorepool,
	       intlistpool,univcoordlistpool,listpool,pathpool,
	       transcriptpool,hitlistpool);
#endif



  debug(printf("found score %d vs sufficient score %d\n",*found_score,sufficient_score));
  if (*found_score <= sufficient_score) {
    return true;

  } else {
    /* For single-end reads, perform extension, since no paired-end read
       to help with extension */
    debug(printf("Single Read: 4.  Performing extension\n"));

    single_read_extend(&(*found_score),this,

		       &(*sense_paths_gplus),&(*sense_paths_gminus),
		       &(*antisense_paths_gplus),&(*antisense_paths_gminus),

		       &(*unextended_sense_paths_gplus),&(*unextended_sense_paths_gminus),
		       &(*unextended_antisense_paths_gplus),&(*unextended_antisense_paths_gminus),

		       queryseq,queryuc_ptr,queryrc,querylength,
		       knownsplicing,knownindels,mismatch_positions_alloc,
		       novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
		       query_compress_fwd,query_compress_rev,
		       localdb_nmismatches_allowed,
		       intlistpool,uintlistpool,univcoordlistpool,
		       listpool,pathpool,transcriptpool,univdiagpool,vectorpool,
		       hitlistpool,spliceendsgen);
  }

  if (*found_score <= sufficient_score) {
    return true;
  } else {
    return false;
  }
}


Path_T *
Stage1_single_read (int *npaths_primary, int *npaths_altloc, int *first_absmq, int *second_absmq,
		    Shortread_T queryseq, EF64_T repetitive_ef64,
		    Knownsplicing_T knownsplicing, Knownindels_T knownindels, Localdb_T localdb,
		    Trdiagpool_T trdiagpool, Univdiagpool_T univdiagpool,
		    Auxinfopool_T auxinfopool, Pathstorepool_T pathstorepool,
		    Intlistpool_T intlistpool, Uintlistpool_T uintlistpool,
		    Univcoordlistpool_T univcoordlistpool, Listpool_T listpool,
		    Trpathpool_T trpathpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		    Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, 
		    Spliceendsgen_T spliceendsgen, bool single_cell_p, bool first_read_p,
		    Pass_T pass) {
  Path_T *patharray;
  T this;
  List_T paths;

  bool any_imperfect_ends_p;

  List_T unextended_sense_paths_gplus = NULL, unextended_sense_paths_gminus = NULL,
    unextended_antisense_paths_gplus = NULL, unextended_antisense_paths_gminus = NULL;

  List_T sense_paths_gplus, sense_paths_gminus,
    antisense_paths_gplus, antisense_paths_gminus;

  int *mismatch_positions_alloc;
  Univdiag_T *novel_univdiags_alloc;
  Univcoord_T *novel_diagonals_alloc;
  unsigned short *localdb_alloc;

  int nmismatches_filter, mincoverage_filter;
  int nmismatches_allowed;

  int sufficient_score, querylength;
  char *queryuc_ptr, *queryrc;
  Compress_T query_compress_fwd, query_compress_rev;

  int found_score;
  Method_T last_method;
  
#if 0
  bool first_read_p;
  if (single_cell_p == true) {
    first_read_p = false;
  } else {
    first_read_p = true;
  }
#endif

  if ((querylength = Shortread_fulllength(queryseq)) < index1part + index1interval - 1) {
    *npaths_primary = *npaths_altloc = 0;
    return (Path_T *) NULL;
  } else {
    queryuc_ptr = Shortread_queryuc_ptr(queryseq);
    queryrc = Shortread_queryrc(queryseq);
    this = Stage1_new(queryuc_ptr,querylength,first_read_p);
  }

  /* nmismatches_allowed means nmismatches_search and is not specified
     by the user.  The user-specified value for -m represents
     nmismatches_filter */
  /* TODO: make this dependent upon the defect_rate */
  nmismatches_allowed = querylength/20; /* was querylength/index1part */

  if (user_nmismatches_filter_float < 0.0) {
    /* Not specified, so don't filter */
    nmismatches_filter = querylength;
  } else if (user_nmismatches_filter_float < 1.0) {
    nmismatches_filter = (int) rint(user_nmismatches_filter_float * (double) querylength);
  } else {
    nmismatches_filter = (int) user_nmismatches_filter_float;
  }

  if (user_mincoverage_filter_float <= 0.0) {
    mincoverage_filter = 0;
  } else if (user_mincoverage_filter_float <= 1.0) {
    /* Assuming that --min-coverage=1 must mean 1.0 and not a coverage of 1 bp */
    mincoverage_filter = (int) rint(user_mincoverage_filter_float * (double) querylength);
  } else {
    mincoverage_filter = (int) user_mincoverage_filter_float;
  }

#if 0
  if (max_insertionlen > querylength) {
    max_insertionlen = querylength;
  }
#endif

  mismatch_positions_alloc = (int *) MALLOC((querylength+MISMATCH_EXTRA)*sizeof(int));

  /* 2 localdb regions possible if shortsplicedist_novelend < 65536 */
  /* 65536 represents the worst possible case where every position in the localdb region matches the query */
  novel_univdiags_alloc = Univdiag_new_array(2 * 65536);
  novel_diagonals_alloc = (Univcoord_T *) MALLOC(2 * 65536 * sizeof(Univcoord_T));
  MALLOC_ALIGN(localdb_alloc,65536 * sizeof(unsigned short)); /* Maximum number of intersections in a localdb region */

  query_compress_fwd = Compress_new_fwd(queryuc_ptr,querylength);
  query_compress_rev = Compress_new_rev(queryuc_ptr,querylength);

  if (mode == STANDARD || mode == CMET_STRANDED || mode == ATOI_STRANDED || mode == TTOC_STRANDED) {
    found_score = querylength;
    sufficient_score = querylength/20;
    last_method = METHOD_INIT;

    if (transcriptome_align_p == true &&
	single_read_tr_paths(&found_score,sufficient_score,&last_method,

			     &sense_paths_gplus,&sense_paths_gminus,
			     &antisense_paths_gplus,&antisense_paths_gminus,
		       
			     this,/*genestrand*/0,queryseq,querylength,
			     knownsplicing,mismatch_positions_alloc,
			     query_compress_fwd,query_compress_rev,
		       
			     nmismatches_allowed,
			     trdiagpool,intlistpool,uintlistpool,univcoordlistpool,
			     listpool,trpathpool,pathpool,vectorpool,transcriptpool,
			     hitlistpool) == true) {

      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else if (genome_align_p == false) {
      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else if (single_read_gen_paths(&found_score,&last_method,
			       
				     &unextended_sense_paths_gplus,&unextended_sense_paths_gminus,
				     &unextended_antisense_paths_gplus,&unextended_antisense_paths_gminus,

				     &sense_paths_gplus,&sense_paths_gminus,
				     &antisense_paths_gplus,&antisense_paths_gminus,

				     this,/*genestrand*/0,
				     queryseq,queryuc_ptr,queryrc,querylength,
				     knownsplicing,knownindels,mismatch_positions_alloc,
				     novel_univdiags_alloc,novel_diagonals_alloc,localdb_alloc,
				     query_compress_fwd,query_compress_rev,nmismatches_allowed,

				     repetitive_ef64,auxinfopool,pathstorepool,univdiagpool,
				     intlistpool,uintlistpool,univcoordlistpool,
				     listpool,pathpool,transcriptpool,vectorpool,
				     hitlistpool,spliceendsgen) == true) {

      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else if (splicingp == false) {
      paths = List_append(sense_paths_gplus,
			  List_append(sense_paths_gminus,
				      List_append(antisense_paths_gplus,
						  antisense_paths_gminus)));

    } else {
      /* Perform fusions, if single-end and splicingp is true */
      debug(printf("Single Read: 6.  Performing fusions\n"));

      debug(Stage1_list_extension(this));

      paths = single_read_fusion(&found_score,this,querylength,

				 unextended_sense_paths_gplus,unextended_sense_paths_gminus,
				 unextended_antisense_paths_gplus,unextended_antisense_paths_gminus,

				 query_compress_fwd,query_compress_rev,
				 queryseq,knownsplicing,nmismatches_allowed,
				 univdiagpool,intlistpool,uintlistpool,univcoordlistpool,
				 listpool,pathpool,transcriptpool,vectorpool,hitlistpool);
    }
  }

  if (paths != NULL) {
    Path_gc(&unextended_sense_paths_gplus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_sense_paths_gminus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_antisense_paths_gplus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);
    Path_gc(&unextended_antisense_paths_gminus,intlistpool,univcoordlistpool,listpool,
	    pathpool,transcriptpool,hitlistpool);

  } else {
    /* As last resort, use unextended paths */
     /* Should have called single_read_extend, which means no further extensions are possible */
    paths = List_append(unextended_sense_paths_gplus,
			List_append(unextended_sense_paths_gminus,
				    List_append(unextended_antisense_paths_gplus,
						unextended_antisense_paths_gminus)));
    /* unextended_sense_paths_gplus = (List_T) NULL; */
    /* unextended_sense_paths_gminus = (List_T) NULL; */
    /* unextended_antisense_paths_gplus = (List_T) NULL;*/
    /* unextended_antisense_paths_gminus = (List_T) NULL; */
  }
    
  if (paths == NULL) {
    *npaths_primary = *npaths_altloc = 0;
    patharray = (Path_T *) NULL;

  } else {
    patharray = (Path_T *) List_to_array_out(paths,NULL);
    patharray = Path_eval_and_sort(&(*npaths_primary),&(*npaths_altloc),
				   &(*first_absmq),&(*second_absmq),patharray,
				   /*npaths*/List_length(paths),
				   query_compress_fwd,query_compress_rev,queryuc_ptr,queryrc,
				   Shortread_quality_string(queryseq),nmismatches_filter,mincoverage_filter,
				   intlistpool,univcoordlistpool,listpool,
				   pathpool,transcriptpool,hitlistpool,/*filterp*/true);


    Hitlistpool_free_list(&paths,hitlistpool
			  hitlistpool_trace(__FILE__,__LINE__));
  }
  
  Compress_free(&query_compress_fwd);
  Compress_free(&query_compress_rev);
  FREE_ALIGN(localdb_alloc);
  FREE(novel_diagonals_alloc);
  FREE(novel_univdiags_alloc[0]);
  FREE(novel_univdiags_alloc);
  FREE(mismatch_positions_alloc);

  /* Do not free paths, since they are now appended to paths */
  Stage1_free(&this,trdiagpool,univdiagpool,auxinfopool,pathstorepool,
	      intlistpool,uintlistpool,univcoordlistpool,
	      listpool,pathpool,trpathpool,transcriptpool,hitlistpool,
	      /*free_paths_p*/false);

  /* FREE(queryrc); -- Now taken from Shortread */

  return patharray;
}


void
Stage1hr_single_setup (Mode_T mode_in, int index1part_in, int index1interval_in, int index1part_tr_in,
		       Transcriptome_T transcriptome_in, bool genome_align_p_in, bool transcriptome_align_p_in,
		       double user_nmismatches_filter_float_in, double user_mincoverage_filter_float_in,
		       bool splicingp_in) {

  mode = mode_in;
  index1part = index1part_in;
  index1interval = index1interval_in;
  index1part_tr = index1part_tr_in;

  transcriptome = transcriptome_in;
  genome_align_p = genome_align_p_in;
  transcriptome_align_p = transcriptome_align_p_in;

  user_nmismatches_filter_float = user_nmismatches_filter_float_in;
  user_mincoverage_filter_float = user_mincoverage_filter_float_in;

  splicingp = splicingp_in;

  return;
}
