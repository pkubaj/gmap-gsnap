static char rcsid[] = "$Id: stopwatch.c 223349 2020-10-28 02:49:25Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


/* #define STANDALONE 1 */

#include "stopwatch.h"

/* #define USE_POSIX_C_TIME 1 -- Now using clock_gettime() */
#define BILLION 1000000000.0

#ifdef STANDALONE
#include <unistd.h>		/* For sysconf */
#elif defined(HAVE_UNISTD_H)
#include <unistd.h>		/* For sysconf */
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>		/* For sys/times.h under AT&T System V Interface */
#endif

#ifdef USE_POSIX_C_TIME
#include <sys/times.h>
#else
#include <time.h>		/* For clock(), clock_gettime() */
#endif

#ifdef STANDALONE
#include <stdio.h>
#include <stdlib.h>
#define MALLOC malloc
#define FREE free
#else
#include "mem.h"
#endif


#define T Stopwatch_T
struct T {
#ifdef USE_POSIX_C_TIME
  struct tms start;
  struct tms stop;
#endif
#ifdef USE_CLOCK
  clock_t start_elapsed;
  clock_t stop_elapsed;
#else
  struct timespec start_elapsed;
  struct timespec stop_elapsed;
#endif
};

T
Stopwatch_new () {
  T new = (T) MALLOC(sizeof(*new));
  return new;
}

void
Stopwatch_free (T *old) {
  if (*old) {
    FREE(*old);
  }
  return;
}

void
Stopwatch_start (T this) {
  if (this != NULL) {
#ifdef USE_POSIX_C_TIME
    this->start_elapsed = times(&this->start);
#elif defined(USE_CLOCK)
    this->start_elapsed = clock();
#else
    clock_gettime(CLOCK_MONOTONIC,&this->start_elapsed);
#endif
  }
  return;
}

/* Returns user seconds elapsed since start of process. */
/* struct tms is defined in <sys/times.h> (see man times); tms values are
   in "clock ticks per second", CLK_TCK */
double 
Stopwatch_stop (T this) {
#ifdef USE_POSIX_C_TIME
#ifndef _SC_CLK_TCK
  long clk_tck = 100;
#else
  /* For some reason, this isn't defined in STANDALONE MODE */
  long clk_tck = sysconf(_SC_CLK_TCK);
#endif
#endif

  if (this == NULL) {
    return 0.0;
  } else {
#ifdef USE_POSIX_C_TIME
    this->stop_elapsed = times(&this->stop);
    /* user time is in stop.tms_utime */
    return (double) (this->stop_elapsed - this->start_elapsed)/(double) clk_tck;
#elif defined(USE_CLOCK)
    this->stop_elapsed = clock();
    return (double) (this->stop_elapsed - this->start_elapsed)/(double) CLOCKS_PER_SEC;
#else
    clock_gettime(CLOCK_MONOTONIC,&this->stop_elapsed);
#if 0
    printf("%u (%u) -> %u sec (%u nsec) => \n",
	   this->start_elapsed.tv_sec,this->start_elapsed.tv_nsec,
	   this->stop_elapsed.tv_sec,this->stop_elapsed.tv_nsec,
	   (this->stop_elapsed.tv_sec - this->start_elapsed.tv_sec) + (double) (this->stop_elapsed.tv_nsec - this->start_elapsed.tv_nsec)/BILLION);
#endif
    
    return (this->stop_elapsed.tv_sec - this->start_elapsed.tv_sec) + (double) (this->stop_elapsed.tv_nsec - this->start_elapsed.tv_nsec)/BILLION;
#endif

  }
}

#ifdef STANDALONE
int
main (int argc, char *argv[]) {
  Stopwatch_T stopwatch = Stopwatch_new();

  Stopwatch_start(stopwatch);
  sleep(5);
  printf("%f\n",Stopwatch_stop(stopwatch));

  return 0;
}
#endif


