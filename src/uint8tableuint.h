/* $Id: ed0a5a2a4943a2e277f0e2fdcf6f17ee32b619f0 $ */
#ifndef UINT8TABLEUINT_INCLUDED
#define UINT8TABLEUINT_INCLUDED
#include "bool.h"
#include "types.h"

#define T Uint8tableuint_T
typedef struct T *T;

extern T
Uint8tableuint_new (int hint);
extern void 
Uint8tableuint_free (T *table);
extern int   
Uint8tableuint_length (T table);
extern unsigned int
Uint8tableuint_put (T table, const UINT8 key, unsigned int value);
extern unsigned int
Uint8tableuint_get (T table, const UINT8 key);
extern unsigned int
Uint8tableuint_remove (T table, const UINT8 key);
extern void   
Uint8tableuint_map (T table,
	       void (*apply)(const UINT8 key, unsigned int *value, void *cl),
	       void *cl);
extern UINT8 *
Uint8tableuint_keys (T table, bool sortp, UINT8 end);
extern UINT8 *
Uint8tableuint_keys_by_timeindex (T table);
extern unsigned int *
Uint8tableuint_values (T table, unsigned int end);

#undef T
#endif
