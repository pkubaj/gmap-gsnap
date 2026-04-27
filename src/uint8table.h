/* $Id: uint8table.h 223349 2020-10-28 02:49:25Z twu $ */
#ifndef UINT8TABLE_INCLUDED
#define UINT8TABLE_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"		/* For HAVE_64_BIT */
#endif

#include "bool.h"
#include "types.h"		/* Needed for HAVE_64_BIT */

#ifdef HAVE_64_BIT

#define T Uint8table_T
typedef struct T *T;

extern T
Uint8table_new (int hint);
extern void 
Uint8table_free (T *table);
extern int   
Uint8table_length (T table);
extern void *
Uint8table_put (T table, const UINT8 key, void *value);
extern void *
Uint8table_get (T table, const UINT8 key);
extern void *
Uint8table_remove (T table, const UINT8 key);
extern void   
Uint8table_map (T table,
	       void (*apply)(const UINT8 key, void **value, void *cl),
	       void *cl);
extern UINT8 *
Uint8table_keys (T table, bool sortp);
extern void
Uint8table_fill_keys (UINT8 *keyarray, T table, bool sortp);
extern UINT8 *
Uint8table_keys_by_timeindex (T table);
extern void **
Uint8table_values (T table);

#undef T

#endif /*HAVE_64_BIT */

#endif
