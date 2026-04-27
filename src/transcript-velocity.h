/* $Id: 00cf557b1e7af81a11ff3c6c2421a0c00fba807b $ */
#ifndef TRANSCRIPT_VELOCITY_INCLUDED
#define TRANSCRIPT_VELOCITY_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

typedef enum {SPLICED, UNSPLICED, BOTH} Velocity_T;

#include "path.h"

extern void
Transcript_velocity_single (Path_T path);

extern void
Transcript_velocity_paired (Path_T path5, Path_T path3);

#endif


