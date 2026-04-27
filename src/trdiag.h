/* $Id: 44f040c84c9667d7e9493a8f51c1904581c009a3 $ */
#ifndef TRDIAG_INCLUDED
#define TRDIAG_INCLUDED

typedef struct Trdiag_T *Trdiag_T;

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "types.h"

#define T Trdiag_T

#if 0
/* Now in trdiagpool.c */
extern T
Trdiag_new (int qstart, int qend, Trcoord_T trdiagonal);
#endif

extern void
Trdiag_free (T *old);
extern void
Trdiag_gc (List_T *list);
extern void
Trdiag_list_gc (List_T *paths);
extern int
Trdiag_list_length (List_T path);

extern int
Trdiag_ascending_cmp (const void *a, const void *b);
extern int
Trdiag_descending_cmp (const void *a, const void *b);
extern int
Trdiag_diagonal_cmp (const void *a, const void *b);

#undef T
#endif


