static char rcsid[] = "$Id: 4bd99054d3295ac2a5d388d66e1c7b6bf57fafa9 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "concordance.h"

#include <stdlib.h>		/* For qsort */
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "assert.h"
#include "mem.h"
#include "univdiag.h"
#include "univdiagdef.h"
#include "transcript.h"

#include "trpath.h"
#include "trpath-convert.h"

#include "path.h"
#include "pathpair.h"
#include "path-solve.h"
#include "path-eval.h"		/* For Path_consolidate */

#ifdef LARGE_GENOMES
#include "intersect-concordance-uint8.h"
#else
#include "intersect-concordance-uint4.h"
#endif


static int subopt_levels;

static Chrpos_T pairmax_transcriptome;
static Chrpos_T pairmax_linear;	/* For two ends that both lack a splice */
static Chrpos_T pairmax_circular;

static bool two_pass_p;
static Chrpos_T adjacent_pairlength;  /* For two ends, one of which has a splice */

static bool *circularp;
static bool merge_samechr_p;


/* #define MAX_HITS 10000 */


#define T Path_T

/* Concordant paths */
#ifdef DEBUG0
#define debug0(x) x
#else
#define debug0(x)
#endif

#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif

/* Exhaustive listing of jj loop */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif

/* binary search */
#ifdef DEBUG10
#define debug10(x) x
#else
#define debug10(x)
#endif




#ifdef DEBUG
static char *
print_sense (int sense) {
  if (sense == SENSE_NULL) {
    return "sense:null";
  } else if (sense == SENSE_ANTI) {
    return "sense:anti";
  } else if (sense == SENSE_FORWARD) {
    return "sense:fwd";
  } else {
    abort();
  }
}
#endif



static List_T
do_transcriptome (int *found_score_paired, int *found_score_5, int *found_score_3, List_T pathpairs,
		  Trpath_T *trpaths5, int ntrpaths5, Trpath_T *trpaths3, int ntrpaths3,
		  Stage1_T stage1_5, Stage1_T stage1_3, Knownsplicing_T knownsplicing,
		  Shortread_T queryseq5, Shortread_T queryseq3,
		  
		  Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		  Compress_T query3_compress_fwd, Compress_T query3_compress_rev,
		  int querylength5, int querylength3,

		  int nmismatches_filter_5, int nmismatches_filter_3,
		  int mincoverage_filter_5, int mincoverage_filter_3,

		  Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		  Listpool_T listpool,  Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		  Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, int sensedir) {
  
  List_T concordant_trpaths5 = NULL, concordant_trpaths3 = NULL;
  List_T unresolved_pathpairs = NULL;

  List_T paths5_gplus, paths5_gminus, paths3_gplus, paths3_gminus;
  int npaths5, npaths3;
  Trpath_T trpath5, trpath3;
  unsigned int trnum_low_5, trnum_high_5;

  T *allpaths5, *allpaths3;

  int i, j;
  T path5, path3;
  Pathpair_T pathpair;


  debug(printf("Entered do_transcriptome with %d trpaths5 and %d trpaths3\n",ntrpaths5,ntrpaths3));

  i = j = 0;
  while (i < ntrpaths5 && j < ntrpaths3) {
    trpath5 = trpaths5[i];
    trpath3 = trpaths3[j];
    debug(Trpath_print(trpath5));
    debug(Trpath_print(trpath3));

    if (trpath5->trnum < trpath3->trnum) {
      i++;
    } else if (trpath3->trnum < trpath5->trnum) {
      j++;
    } else {
      concordant_trpaths5 = Listpool_push(concordant_trpaths5,listpool,(void *) trpath5
					  listpool_trace(__FILE__,__LINE__));
      concordant_trpaths3 = Listpool_push(concordant_trpaths3,listpool,(void *) trpath3
					  listpool_trace(__FILE__,__LINE__));
      i++;
      j++;
    }
  }

  if (sensedir == SENSE_FORWARD) {
    paths5_gplus = paths5_gminus = (List_T) NULL;
    Trpath_convert_sense(&(*found_score_5),&paths5_gplus,&paths5_gminus,
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 concordant_trpaths5,/*first_read_p*/true,
			 queryseq5,querylength5,stage1_5,knownsplicing,
			 query5_compress_fwd,query5_compress_rev,
			 intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
			 vectorpool,transcriptpool,hitlistpool,/*crossover_sense_p*/false);

    paths3_gplus = paths3_gminus = (List_T) NULL;
    Trpath_convert_sense(&(*found_score_3),&paths3_gplus,&paths3_gminus,
			 /*antisense_paths_gplus*/NULL,/*antisense_paths_gminus*/NULL,

			 concordant_trpaths3,/*first_read_p*/false,
			 queryseq3,querylength3,stage1_3,knownsplicing,
			 query3_compress_fwd,query3_compress_rev,
			 intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
			 vectorpool,transcriptpool,hitlistpool,/*crossover_sense_p*/false);
  } else {
    paths5_gplus = paths5_gminus = (List_T) NULL;
    Trpath_convert_antisense(&(*found_score_5),&paths5_gplus,&paths5_gminus,
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     concordant_trpaths5,/*first_read_p*/true,
			     queryseq5,querylength5,stage1_5,knownsplicing,
			     query5_compress_fwd,query5_compress_rev,
			     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
			     vectorpool,transcriptpool,hitlistpool,/*crossover_sense_p*/false);

    paths3_gplus = paths3_gminus = (List_T) NULL;
    Trpath_convert_antisense(&(*found_score_3),&paths3_gplus,&paths3_gminus,
			     /*sense_paths_gplus*/NULL,/*sense_paths_gminus*/NULL,

			     concordant_trpaths3,/*first_read_p*/false,
			     queryseq3,querylength3,stage1_3,knownsplicing,
			     query3_compress_fwd,query3_compress_rev,
			     intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,
			     vectorpool,transcriptpool,hitlistpool,/*crossover_sense_p*/false);
  }


  /* Analogous to do_genome, for gplus */
  paths5_gplus = Path_consolidate(paths5_gplus,queryseq5,
				  query5_compress_fwd,query5_compress_rev,
				  uintlistpool,intlistpool,univcoordlistpool,
				  listpool,pathpool,transcriptpool,hitlistpool);
  paths3_gplus = Path_consolidate(paths3_gplus,queryseq3,
				  query3_compress_fwd,query3_compress_rev,
				  uintlistpool,intlistpool,univcoordlistpool,
				  listpool,pathpool,transcriptpool,hitlistpool);

  if ((npaths5 = List_length(paths5_gplus)) > 0 &&
      (npaths3 = List_length(paths3_gplus)) > 0) {
    debug(printf("gplus: %d consolidated paths5 and %d consolidated paths3\n",
		 npaths5,npaths3));

    allpaths5 = (T *) List_to_array(paths5_gplus,NULL);
    allpaths3 = (T *) List_to_array(paths3_gplus,NULL);
    qsort(allpaths5,npaths5,sizeof(T),Path_trnum_high_cmp);
    qsort(allpaths3,npaths3,sizeof(T),Path_trnum_low_cmp);

    i = j = 0;
    while (i < npaths5) {
      path5 = allpaths5[i];
#ifdef DEBUG
      printf("i=%d/%d ",i,npaths5);
      Path_print(path5);
#endif
      trnum_low_5 = Path_trnum_low(path5);
      trnum_high_5 = Path_trnum_high(path5);

      while (j >= 0 && Path_trnum_high(allpaths3[j]) >= trnum_low_5) {
#ifdef DEBUG
	printf("  backup: j=%d/%d ",j,npaths3);
	Path_print(allpaths3[j]);
#endif
	j--;
      }
      j++;		/* Finish backup */

      while (j < npaths3 && Path_trnum_low(allpaths3[j]) <= trnum_high_5) {
	path3 = allpaths3[j];
#ifdef DEBUG
	printf("  forward: j=%d/%d ",j,npaths3);
	Path_print(path3);
#endif
	if (Pathpair_transcript_intersectp(path5,path3) == true &&
	    (pathpair =
	     Pathpair_new_concordant(&unresolved_pathpairs,
				     path5,path3,queryseq5,queryseq3,/*plusp*/true,

				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,

				     intlistpool,univcoordlistpool,listpool,
				     pathpool,vectorpool,transcriptpool,hitlistpool,
				     /*check_inner_p*/false,
				     /*copyLp*/true,/*copyHp*/true)) != NULL) {
	  debug(printf("(1) Pathpair found, low/high: %u\n",
		       Pathpair_insertlength(pathpair)));
	  debug(Path_print(path5));
	  debug(Path_print(path3));
	  if (Pathpair_found_score(pathpair) < *found_score_paired) {
	    *found_score_paired = Pathpair_found_score(pathpair);
	  }
	  pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				   hitlistpool_trace(__FILE__,__LINE__));
	}

	j++;
      }
      
      j--;		/* Finish advance */
      
      i++;
    }

    FREE(allpaths3);
    FREE(allpaths5);
  }

  Path_gc(&paths3_gplus,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&paths5_gplus,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);


  /* Analogous to do_genome, for gminus */
  paths5_gminus = Path_consolidate(paths5_gminus,queryseq5,
				   query5_compress_fwd,query5_compress_rev,
				   uintlistpool,intlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,hitlistpool);
  paths3_gminus = Path_consolidate(paths3_gminus,queryseq3,
				   query3_compress_fwd,query3_compress_rev,
				   uintlistpool,intlistpool,univcoordlistpool,
				   listpool,pathpool,transcriptpool,hitlistpool);

  if ((npaths5 = List_length(paths5_gminus)) > 0 &&
      (npaths3 = List_length(paths3_gminus)) > 0) {
    debug(printf("gminus: %d consolidated paths5 and %d consolidated paths3\n",
		 npaths5,npaths3));

    allpaths5 = (T *) List_to_array(paths5_gminus,NULL);
    allpaths3 = (T *) List_to_array(paths3_gminus,NULL);
    qsort(allpaths5,npaths5,sizeof(T),Path_trnum_high_cmp);
    qsort(allpaths3,npaths3,sizeof(T),Path_trnum_low_cmp);

    i = j = 0;
    while (i < npaths5) {
      path5 = allpaths5[i];
#ifdef DEBUG
      printf("i=%d/%d ",i,npaths5);
      Path_print(path5);
#endif
      trnum_low_5 = Path_trnum_low(path5);
      trnum_high_5 = Path_trnum_high(path5);

      while (j >= 0 && Path_trnum_high(allpaths3[j]) >= trnum_low_5) {
#ifdef DEBUG
	printf("  backup: j=%d/%d ",j,npaths3);
	Path_print(allpaths3[j]);
#endif
	j--;
      }
      j++;		/* Finish backup */

      while (j < npaths3 && Path_trnum_low(allpaths3[j]) <= trnum_high_5) {
	path3 = allpaths3[j];
#ifdef DEBUG
	printf("  forward: j=%d/%d ",j,npaths3);
	Path_print(path3);
#endif
	if (Pathpair_transcript_intersectp(path5,path3) == true &&
	    (pathpair =
	     Pathpair_new_concordant(&unresolved_pathpairs,
				     path3,path5,queryseq3,queryseq5,/*plusp*/false,
				     nmismatches_filter_5,nmismatches_filter_3,
				     mincoverage_filter_5,mincoverage_filter_3,
				     intlistpool,univcoordlistpool,listpool,
				     pathpool,vectorpool,transcriptpool,hitlistpool,
				     /*check_inner_p*/false,
				     /*copyLp*/true,/*copyHp*/true)) != NULL) {
	  debug(printf("(2) Pathpair found, low/high, %u:\n",
		       Pathpair_insertlength(pathpair)));
	  debug(Path_print(path3));
	  debug(Path_print(path5));
	  if (Pathpair_found_score(pathpair) < *found_score_paired) {
	    *found_score_paired = Pathpair_found_score(pathpair);
	  }
	  pathpairs = Hitlist_push(pathpairs,hitlistpool,(void *) pathpair
				   hitlistpool_trace(__FILE__,__LINE__));
	}

	j++;
      }
      
      j--;		/* Finish advance */
      
      i++;
    }

    FREE(allpaths3);
    FREE(allpaths5);
  }


  Path_gc(&paths3_gminus,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);
  Path_gc(&paths5_gminus,intlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,hitlistpool);

  return pathpairs;
}


List_T
Concordance_tr (int *found_score_paired, int *found_score_5, int *found_score_3, List_T pathpairs,

		List_T newtrpaths5, List_T newtrpaths3,
		List_T trpaths5, List_T trpaths3,
		    
		Compress_T query5_compress_fwd, Compress_T query5_compress_rev,
		Compress_T query3_compress_fwd, Compress_T query3_compress_rev,

		Shortread_T queryseq5, Shortread_T queryseq3,
		int querylength5, int querylength3,
		Stage1_T stage1_5, Stage1_T stage1_3, Knownsplicing_T knownsplicing,

		int nmismatches_filter_5, int nmismatches_filter_3,
		int mincoverage_filter_5, int mincoverage_filter_3,

		Intlistpool_T intlistpool, Uintlistpool_T uintlistpool, Univcoordlistpool_T univcoordlistpool,
		Listpool_T listpool, Pathpool_T pathpool, Transcriptpool_T transcriptpool,
		Vectorpool_T vectorpool, Hitlistpool_T hitlistpool, int sensedir) {

  Trpath_T *alltrpaths5, *alltrpaths3;
  int nnewtrpaths5, nnewtrpaths3, ntrpaths5, ntrpaths3;

#ifdef DEBUG
  List_T p;
  printf("Entered Concordance_tr with %d pathpairs, newtrpaths5 %d, newtrpaths3 %d, trpaths5 %d, trpaths3 %d\n",
	 List_length(pathpairs),List_length(newtrpaths5),List_length(newtrpaths3),
	 List_length(trpaths5),List_length(trpaths3));
  printf("newtrpaths5:");
  for (p = newtrpaths5; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");

  printf("trpaths5:");
  for (p = trpaths5; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");

  printf("newtrpaths3:");
  for (p = newtrpaths3; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");

  printf("trpaths3:");
  for (p = trpaths3; p != NULL; p = List_next(p)) {
    printf(" %p",List_head(p));
  }
  printf("\n");
#endif

  nnewtrpaths5 = List_length(newtrpaths5);
  nnewtrpaths3 = List_length(newtrpaths3);
  ntrpaths5 = List_length(trpaths5);
  ntrpaths3 = List_length(trpaths3);

  debug(printf("%d new + %d old paths5 and %d new + %d old paths3\n",
	       nnewtrpaths5,ntrpaths5,nnewtrpaths3,ntrpaths3));

  if (nnewtrpaths5 + ntrpaths5 > 0 && nnewtrpaths3 + ntrpaths3 > 0) {
    debug(printf("Allocating %d+%d paths for allpaths5 and %d+%d paths for allpaths3\n",
		 nnewtrpaths5,ntrpaths5,nnewtrpaths3,ntrpaths3));
    alltrpaths5 = (Trpath_T *) MALLOC((nnewtrpaths5 + ntrpaths5)*sizeof(Path_T));
    alltrpaths3 = (Trpath_T *) MALLOC((nnewtrpaths3 + ntrpaths3)*sizeof(Path_T));

    /* (1) New trpaths 5 vs New trpaths 3 */
    List_fill_array((void **) alltrpaths5,newtrpaths5);
    List_fill_array((void **) alltrpaths3,newtrpaths3);

    qsort(alltrpaths5,nnewtrpaths5,sizeof(Trpath_T),Trpath_trnum_cmp);
    qsort(alltrpaths3,nnewtrpaths3,sizeof(Trpath_T),Trpath_trnum_cmp);

    pathpairs = do_transcriptome(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 alltrpaths5,nnewtrpaths5,
				 alltrpaths3,nnewtrpaths3,
				 stage1_5,stage1_3,knownsplicing,
				 /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,
				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 querylength5,querylength3,
				 
				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,sensedir);

    /* (2) New trpaths 5 vs Old trpaths 3 */
    List_fill_array((void **) &(alltrpaths3[nnewtrpaths3]),trpaths3);
    qsort(&(alltrpaths3[nnewtrpaths3]),ntrpaths3,sizeof(Trpath_T),Trpath_trnum_cmp);
    
    pathpairs = do_transcriptome(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 alltrpaths5,nnewtrpaths5,
				 &(alltrpaths3[nnewtrpaths3]),ntrpaths3,
				 stage1_5,stage1_3,knownsplicing,
				 /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,
				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 querylength5,querylength3,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,sensedir);

    /* (3) Old trpaths 5 vs New trpaths 3 */
    List_fill_array((void **) &(alltrpaths5[nnewtrpaths5]),trpaths5);
    qsort(&(alltrpaths5[nnewtrpaths5]),ntrpaths5,sizeof(T),Trpath_trnum_cmp);

    pathpairs = do_transcriptome(&(*found_score_paired),&(*found_score_5),&(*found_score_3),pathpairs,
				 &(alltrpaths5[nnewtrpaths5]),ntrpaths5,
				 alltrpaths3,nnewtrpaths3,
				 stage1_5,stage1_3,knownsplicing,
				 /*queryseqL*/queryseq5,/*queryseqH*/queryseq3,
				 query5_compress_fwd,query5_compress_rev,
				 query3_compress_fwd,query3_compress_rev,
				 querylength5,querylength3,

				 nmismatches_filter_5,nmismatches_filter_3,
				 mincoverage_filter_5,mincoverage_filter_3,
				 intlistpool,uintlistpool,univcoordlistpool,listpool,pathpool,transcriptpool,
				 vectorpool,hitlistpool,sensedir);

    FREE(alltrpaths3);
    FREE(alltrpaths5);
  }

  return pathpairs;
}


void
Concordance_setup (int subopt_levels_in, Chrpos_T pairmax_transcriptome_in,
		   Chrpos_T pairmax_linear_in, Chrpos_T pairmax_circular_in,
		   bool *circularp_in, bool merge_samechr_p_in, bool two_pass_p_in) {
  subopt_levels = subopt_levels_in;

  pairmax_transcriptome = pairmax_transcriptome_in;
  pairmax_linear = pairmax_linear_in;
  pairmax_circular = pairmax_circular_in;

  adjacent_pairlength = 0;

  circularp = circularp_in;
  merge_samechr_p = merge_samechr_p_in;
  two_pass_p = two_pass_p_in;

  return;
}


