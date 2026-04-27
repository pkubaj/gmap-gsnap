static char rcsid[] = "$Id: uint8table.c 223349 2020-10-28 02:49:25Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "uint8table.h"
#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>		/* For qsort */
#include <string.h>		/* For strcmp */
#include "mem.h"
#include "assert.h"

#ifdef HAVE_64_BIT		/* Defined in types.h, included in uint8table.h */

#define T Uint8table_T
struct T {
  int size;
  int length;
  unsigned int timestamp;
  struct binding {
    struct binding *link;
    UINT8 key;
    void *value;
    unsigned int timeindex;
  } **buckets;
};



T 
Uint8table_new (int hint) {
  T table;
  int i;
  static int primes[] = { 509, 509, 1021, 2053, 4093,
			  8191, 16381, 32771, 65521, INT_MAX };

  assert(hint >= 0);
  for (i = 1; primes[i] < hint; i++) {
  }
  table = (T) MALLOC(sizeof(*table) +
		     primes[i-1]*sizeof(table->buckets[0]));
  table->size = primes[i-1];
  table->buckets = (struct binding **)(table + 1);
  for (i = 0; i < table->size; i++) {
    table->buckets[i] = NULL;
  }
  table->length = 0;
  table->timestamp = 0;
  return table;
}

void *
Uint8table_get (T table, const UINT8 key) {
  int i;
  struct binding *p;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = key % table->size;
  /* printf("Doing Uint8table_get on %s at bucket %d\n",(char *) key, i); */
  for (p = table->buckets[i]; p; p = p->link) {
    /* printf("  Comparing %s with %s at %p, key = %p\n",(char *) key, (char *) p->key, p, p->key); */
    if (key == p->key) {
      break;
    }
  }
  return p ? p->value : NULL;
}

void *
Uint8table_put (T table, const UINT8 key, void *value) {
  int i;
  struct binding *p;
  void *prev;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  i = key % table->size;
  for (p = table->buckets[i]; p; p = p->link) {
    if (key == p->key) {
      break;
    }
  }
  if (p == NULL) {
    NEW(p);
    p->key = key;
    /* printf("Doing Uint8table_put at %p, key = %p\n",p,p->key); */
    p->link = table->buckets[i];
    table->buckets[i] = p;
    table->length++;
    prev = NULL;
  } else {
    prev = p->value;
  }
  p->value = value;
  p->timeindex = table->timestamp;
  table->timestamp++;
  return prev;
}

int 
Uint8table_length (T table) {
  assert(table);
  return table->length;
}

void 
Uint8table_map (T table,
	       void (*apply)(const UINT8 key, void **value, void *cl),
	       void *cl) {
  int i;
  struct binding *p;

  assert(table);
  assert(apply);
  for (i = 0; i < table->size; i++)
    for (p = table->buckets[i]; p; p = p->link) {
      apply(p->key, &p->value, cl);
    }
}

void *
Uint8table_remove (T table, const UINT8 key) {
  int i;
  struct binding **pp;

  assert(table);
  /* assert(key); -- Doesn't hold for atomic 0 */
  table->timestamp++;
  i = key % table->size;
  for (pp = &table->buckets[i]; *pp; pp = &(*pp)->link) {
    if (key == (*pp)->key) {
      struct binding *p = *pp;
      void *value = p->value;
      *pp = p->link;
      FREE(p);
      table->length--;
      return value;
    }
  }
  return NULL;
}

static int
uint8_compare (const void *a, const void *b) {
  UINT8 x = * (UINT8 *) a;
  UINT8 y = * (UINT8 *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


UINT8 *
Uint8table_keys (T table, bool sortp) {
  UINT8 *keyarray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  keyarray = (UINT8 *) CALLOC(table->length+1,sizeof(UINT8));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      keyarray[j++] = p->key;
    }
  }

  if (sortp == true) {
    qsort(keyarray,table->length,sizeof(UINT8),uint8_compare);
  }

  return keyarray;
}


void
Uint8table_fill_keys (UINT8 *keyarray, T table, bool sortp) {
  int i, j = 0;
  struct binding *p;

  assert(table);
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      keyarray[j++] = p->key;
    }
  }

  if (sortp == true) {
    qsort(keyarray,table->length,sizeof(UINT8),uint8_compare);
  }

  return;
}

static int
timeindex_cmp (const void *x, const void *y) {
  struct binding *a = * (struct binding **) x;
  struct binding *b = * (struct binding **) y;

  if (a->timeindex < b->timeindex) {
    return -1;
  } else if (a->timeindex > b->timeindex) {
    return +1;
  } else {
    return 0;
  }
}


UINT8 *
Uint8table_keys_by_timeindex (T table) {
  UINT8 *keyarray;
  int i, j = 0;
  struct binding **buckets, *p;

  assert(table);
  buckets = (struct binding **) CALLOC(table->length+1,sizeof(struct binding *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      buckets[j++] = p;
    }
  }
  qsort(buckets,table->length,sizeof(struct binding *),timeindex_cmp);

  keyarray = (UINT8 *) CALLOC(table->length,sizeof(UINT8));
  for (j = 0; j < table->length; j++) {
    p = buckets[j];
    keyarray[j] = p->key;
  }
  FREE(buckets);

  return keyarray;
}


void **
Uint8table_values (T table) {
  void **valuearray;
  int i, j = 0;
  struct binding *p;

  assert(table);
  valuearray = (void **) CALLOC(table->length,sizeof(void *));
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      valuearray[j++] = (void *) p->value;
    }
  }
  return valuearray;
}

void 
Uint8table_free (T *table) {
  assert(table && *table);
  if ((*table)->length > 0) {
    int i;
    struct binding *p, *q;
    for (i = 0; i < (*table)->size; i++) {
      for (p = (*table)->buckets[i]; p; p = q) {
	q = p->link;
	FREE(p);
      }
    }
  }
  FREE(*table);
}


#endif	/* HAVE_64_BIT */

