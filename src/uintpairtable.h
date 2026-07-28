/* $Id$ */
#ifndef UINTPAIRTABLE_INCLUDED
#define UINTPAIRTABLE_INCLUDED
#include "bool.h"

#define T Uintpairtable_T
typedef struct T *T;

extern T
Uintpairtable_new (int hint);
extern void 
Uintpairtable_free (T *table);
extern int   
Uintpairtable_length (T table);
extern void *
Uintpairtable_put (T table, const unsigned int key1, const unsigned int key2, void *value);
extern void *
Uintpairtable_get (T table, const unsigned int key1, const unsigned int key2);
extern void
Uintpairtable_keys (unsigned int **keys1, unsigned int **keys2,
		    T table, bool sortp, unsigned int end);
extern void **
Uintpairtable_values (T table);

#undef T
#endif
