static char rcsid[] = "$Id: 73edd6d2fcccbb00f12d4e6312bbfc39b8f61603 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "method.h"

#include <stdlib.h>


void
Method_dump (bool *methodp) {
  int i;

  for (i = 0; i < NMETHODS; i++) {
    printf("Method %d: %d\n",i,methodp[i]);
  }
  return;
}


char *
Method_string (Method_T method) {
  /* May want to turn on NO_COMPARE in substring.c also */
  switch (method) {
  case METHOD_INIT: return "method.init"; /* For debugging */
  case TR_EXACT1: return "tr.exact1";
  case TR_EXACT2: return "tr.exact2";
  case TR_ANYPAIR: return "tr.anypair";
  case TR_PREVALENT: return "tr.prev";
  case TR_END: return "tr.end";
  case TR_EXT: return "tr.ext";
  case TR_METHOD: return "tr.method"; /* For debugging */
  case KMER_EXACT1: return "exact1";
  case KMER_ANYPAIR: return "anypair";
  case KMER_PREVALENT: return "kmer.prev";
  case KMER_ANCHORED: return "kmer.anchor";
  case KMER_PAIRED: return "paired";
  case KMER_WIDEST: return "kmer.wide";
  case END: return "end";
  case EXT: return "ext";	/* Level for concordance */
  case EXT_NOSPLICE: return "ext.nosplice";
  case EXT_SPLICE: return "ext.splice";
  case KMER_COMPLETE: return "complete";
  case TR_COMPLETE: return "tr.complete";
  case TR_APPROX: return "tr.approx";
  case SEGMENT1: return "seg1";
  case KMER_EXACT2: return "exact2";
  case KMER_APPROX: return "approx";
  case LOCAL_MATE: return "mate";
  case EXHAUSTIVE: return "exhaustive";
  case FUSION: return "fusion";
  default: abort();
  }
}


void
Method_samprint (Filestring_T fp, Method_T method) {
  switch (method) {
  case TR_EXACT1: FPRINTF(fp,"\tXG:Z:tr.exact1"); break;
  case TR_EXACT2: FPRINTF(fp,"\tXG:Z:tr.exact2"); break;
  case TR_ANYPAIR: FPRINTF(fp,"\tXG:Z:tr.anypair"); break;
  case TR_PREVALENT: FPRINTF(fp,"\tXG:Z:tr.prev"); break;
  case TR_END: FPRINTF(fp,"\tXG:Z:tr.end"); break;
  case TR_EXT: FPRINTF(fp,"\tXG:Z:tr.ext"); break;
  case KMER_EXACT1: FPRINTF(fp,"\tXG:Z:exact1"); break;
  case KMER_ANYPAIR: FPRINTF(fp,"\tXG:Z:anypair"); break;
  case KMER_PREVALENT: FPRINTF(fp,"\tXG:Z:kmer.prev"); break;
  case KMER_ANCHORED: FPRINTF(fp,"\tXG:Z:kmer.anchor"); break;
  case KMER_PAIRED: FPRINTF(fp,"\tXG:Z:paired"); break;
  case KMER_WIDEST: FPRINTF(fp,"\tXG:Z:kmer.wide"); break;
  case END: FPRINTF(fp,"\tXG:Z:end"); break;
  case EXT: FPRINTF(fp,"\tXG:Z:ext"); break; 
  case EXT_NOSPLICE: FPRINTF(fp,"\tXG:Z:ext.nosplice"); break; 
  case EXT_SPLICE: FPRINTF(fp,"\tXG:Z:ext.splice"); break; 
  case KMER_COMPLETE: FPRINTF(fp,"\tXG:Z:complete"); break;
  case TR_COMPLETE: FPRINTF(fp,"\tXG:Z:tr.complete"); break;
  case TR_APPROX: FPRINTF(fp,"\tXG:Z:tr.approx"); break;
  case SEGMENT1: FPRINTF(fp,"\tXG:Z:seg1"); break;
  case KMER_EXACT2: FPRINTF(fp,"\tXG:Z:exact2"); break;
  case KMER_APPROX: FPRINTF(fp,"\tXG:Z:approx"); break;
  case LOCAL_MATE: FPRINTF(fp,"\tXG:Z:mate"); break;
  case EXHAUSTIVE: FPRINTF(fp,"\tXG:Z:exhaustive"); break;
  case FUSION: FPRINTF(fp,"\tXG:Z:fusion"); break;
  default: abort();
  }

  return;
}

void
Method_print (Filestring_T fp, Method_T method) {
  switch (method) {
  case TR_EXACT1: FPRINTF(fp,"\tmethod:tr.exact1"); break;
  case TR_EXACT2: FPRINTF(fp,"\tmethod:tr.exact2"); break;
  case TR_ANYPAIR: FPRINTF(fp,"\tmethod:tr.anypair"); break;
  case TR_PREVALENT: FPRINTF(fp,"\tmethod:tr.prev"); break;
  case TR_END: FPRINTF(fp,"\tmethod:tr.end"); break;
  case TR_EXT: FPRINTF(fp,"\tmethod:tr.ext"); break;
  case KMER_EXACT1: FPRINTF(fp,"\tmethod:exact1"); break;
  case KMER_ANYPAIR: FPRINTF(fp,"\tmethod:anypair"); break;
  case KMER_PREVALENT: FPRINTF(fp,"\tmethod:kmer.prev"); break;
  case KMER_ANCHORED: FPRINTF(fp,"\tmethod:kmer.anchor"); break;
  case KMER_PAIRED: FPRINTF(fp,"\tmethod:paired"); break;
  case KMER_WIDEST: FPRINTF(fp,"\tmethod:kmer.wide"); break;
  case END: FPRINTF(fp,"\tmethod:end"); break;
  case EXT: FPRINTF(fp,"\tmethod:ext"); break;
  case EXT_NOSPLICE: FPRINTF(fp,"\tmethod:ext.nosplice"); break;
  case EXT_SPLICE: FPRINTF(fp,"\tmethod:ext.splice"); break;
  case KMER_COMPLETE: FPRINTF(fp,"\tmethod:complete"); break;
  case TR_COMPLETE: FPRINTF(fp,"\tmethod:tr.complete"); break;
  case TR_APPROX: FPRINTF(fp,"\tmethod:tr.approx"); break;
  case SEGMENT1: FPRINTF(fp,"\tmethod:seg1"); break;
  case KMER_EXACT2: FPRINTF(fp,"\tmethod:exact2"); break;
  case KMER_APPROX: FPRINTF(fp,"\tmethod:approx"); break;
  case LOCAL_MATE: FPRINTF(fp,"\tmethod:mate"); break;
  case EXHAUSTIVE: FPRINTF(fp,"\tmethod:exhaustive"); break;
  case FUSION: FPRINTF(fp,"\tmethod:fusion"); break;
  default: abort();
  }

  return;
}

