/* $Id$ */
#ifndef UINT8PAIRTABLE_INCLUDED
#define UINT8PAIRTABLE_INCLUDED
#include "bool.h"
#include "types.h"

#ifdef HAVE_64_BIT

#define T Uint8pairtable_T
typedef struct T *T;

extern T
Uint8pairtable_new (int hint);
extern void 
Uint8pairtable_free (T *table);
extern int   
Uint8pairtable_length (T table);
extern void *
Uint8pairtable_put (T table, const UINT8 key1, const UINT8 key2, void *value);
extern void *
Uint8pairtable_get (T table, const UINT8 key1, const UINT8 key2);
extern void
Uint8pairtable_keys (UINT8 **keys1, UINT8 **keys2,
		     T table, bool sortp, UINT8 end);
extern void **
Uint8pairtable_values (T table);

#undef T

#endif /*HAVE_64_BIT */

#endif

