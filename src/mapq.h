/* $Id: mapq.h 226345 2023-03-04 14:44:29Z twu $ */
#ifndef MAPQ_INCLUDED
#define MAPQ_INCLUDED

#include "types.h"
#include "compress.h"
#include "genomicpos.h"
#include "genomebits.h"


#define MAX_QUALITY_SCORE_INPUT 96	/* Was 40 */
#define MAX_QUALITY_SCORE 40
#define MAPQ_MAXIMUM_SCORE 40

extern void
MAPQ_init (int quality_score_adj_in);
extern int
MAPQ_max_quality_score (char *quality_string, int querylength);
extern float
MAPQ_loglik_string (char *genomic_diff, char *quality_string, int querylength, bool plusp);


#endif

