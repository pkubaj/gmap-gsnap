#ifndef UNIVDIAG_INCLUDED
#define UNIVDIAG_INCLUDED

typedef struct Univdiag_T *Univdiag_T;

#include "bool.h"
#include "list.h"
#include "genomicpos.h"
#include "types.h"

#define T Univdiag_T

#if 0
/* Now in univdiagpool.c */
extern T
Univdiag_new (int qstart, int qend, Univcoord_T univdiagonal);
#endif

extern void
Univdiag_free_array (Univdiag_T **array);
extern T *
Univdiag_new_array (int n);
extern void
Univdiag_free_array_keep (Univdiag_T **array);
extern T *
Univdiag_new_array_keep (int n);

extern void
Univdiag_transfer (T new, T old);
extern void
Univdiag_transfer_list_to_array (T *array, List_T p);
extern void
Univdiag_reverse_inplace (T *array, int starti, int endi);

extern void
Univdiag_free (T *old);
extern void
Univdiag_gc (List_T *list);
extern void
Univdiag_list_gc (List_T *paths);
extern int
Univdiag_list_length (List_T path);

extern bool
Univdiag_equal (T x, T y);
extern int
Univdiag_ascending_cmp (const void *a, const void *b);
extern int
Univdiag_descending_cmp (const void *a, const void *b);

extern int
Univdiag_struct_diagonal_fwd_cmp (const void *a, const void *b);
extern int
Univdiag_struct_diagonal_rev_cmp (const void *a, const void *b);

extern int
Univdiag_struct_cmp (const void *a, const void *b);

#undef T
#endif


