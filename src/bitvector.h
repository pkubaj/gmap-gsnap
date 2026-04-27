/* $Id: cd447da18b5b9c7da830403a9c0517f2d56f2b3d $ */
#ifndef BITVECTOR_INCLUDED
#define BITVECTOR_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stddef.h>
#include "bool.h"
#include "mem.h"
#include "list.h"
#include "filestring.h"
#include "uintlist.h"
#include "doublelist.h"


#define T Bitvector_T
typedef struct T *T;

typedef unsigned int Bitvector_elt_T;
#define BITS_PER_ELT 32		/* 8 bits in an unsigned int */
struct T {
  int nelts;			/* No need to store nbits */
  Bitvector_elt_T *elts;
};

typedef struct Assoc_T *Assoc_T;
struct Assoc_T {
  T bitvector;
  void *value;
};

extern void
Assoc_free (Assoc_T *old);

extern Assoc_T
Assoc_new (T bitvector, void *value);

extern Assoc_T
Assoc_find (T this, List_T assocs);

extern bool
Bitvector_find (T this, List_T bitvectors);

extern double
Bitvector_double_find (T this, List_T bitvectors, Doublelist_T values);

extern Doublelist_T
Bitvector_double_find_ptr (T this, List_T bitvectors, Doublelist_T values);

extern void
Bitvector_free (T *old);

extern T
Bitvector_new (int nbits);

extern T
Bitvector_copy (T old);

extern int
Bitvector_cmp (T x, T y);

extern bool
Bitvector_null_p (T this);

extern void
Bitvector_clear_all (T this);

extern bool
Bitvector_intersect_p (T x, T y);

extern bool
Bitvector_diff_p (T x, T y);

extern T
Bitvector_union (T dest, T source, bool copyp);

extern unsigned int
Bitvector_high_bit (T this);

extern unsigned int
Bitvector_low_bit (T this);

extern Uintlist_T
Bitvector_to_uintlist (T this);

extern void
Bitvector_print (Filestring_T fp, T this, char sep);

extern void
Bitvector_print_file (FILE *out, T this, char sep);


#if !defined(HAVE_INLINE)

extern void Bitvector_set (T this, unsigned int i);
extern void Bitvector_clear (T this, unsigned int i);
extern bool Bitvector_get (T this, unsigned int i);

#else

static inline void
Bitvector_set (T this, unsigned int i) {
  this->elts[i / BITS_PER_ELT] |= (1 << (i % BITS_PER_ELT));
}

static inline void
Bitvector_clear (T this, unsigned int i) {
  this->elts[i / BITS_PER_ELT] &= ~(1 << (i % BITS_PER_ELT));
}

static inline bool
Bitvector_get (T this, unsigned int i) {
  return this->elts[i / BITS_PER_ELT] & (1 << (i % BITS_PER_ELT)) ? true : false;
}

#endif

#undef T
#endif

