static char rcsid[] = "$Id$";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "uintpairtable.h"
#include <stdio.h>
#include <limits.h>
#include <stddef.h>
#include <stdlib.h>		/* For qsort */
#include <string.h>		/* For strcmp */
#include "mem.h"
#include "assert.h"
#include "sedgesort.h"

#define T Uintpairtable_T
struct T {
  int size;
  int length;
  unsigned int timestamp;
  struct binding {
    struct binding *link;
    unsigned int key1;
    unsigned int key2;
    void *value;
    unsigned int timeindex;
  } **buckets;
};



T 
Uintpairtable_new (int hint) {
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
Uintpairtable_get (T table, const unsigned int key1, const unsigned int key2) {
  int i;
  struct binding *p;

  assert(table);
  /* assert(key1); -- Doesn't hold for atomic 0 */
  /* assert(key2); -- Doesn't hold for atomic 0 */
  i = (key1 + key2) % table->size;
  /* printf("Doing Uintpairtable_get on %s at bucket %d\n",(char *) key, i); */
  for (p = table->buckets[i]; p; p = p->link) {
    /* printf("  Comparing %s with %s at %p, key = %p\n",(char *) key, (char *) p->key, p, p->key); */
    if (key1 == p->key1 && key2 == p->key2) {
      break;
    }
  }
  return p ? p->value : NULL;
}

void *
Uintpairtable_put (T table, const unsigned int key1, const unsigned int key2, void *value) {
  int i;
  struct binding *p;
  void *prev;

  assert(table);
  /* assert(key1); -- Doesn't hold for atomic 0 */
  /* assert(key2); -- Doesn't hold for atomic 0 */
  i = (key1 + key2) % table->size;
  for (p = table->buckets[i]; p; p = p->link) {
    if (key1 == p->key1 && key2 == p->key2) {
      break;
    }
  }
  if (p == NULL) {
    NEW(p);
    p->key1 = key1;
    p->key2 = key2;
    /* printf("Doing Uintpairtable_put at %p, key = %p\n",p,p->key); */
    p->link = table->buckets[i];
    table->buckets[i] = p;
    table->length++;
    prev = 0;
  } else {
    prev = p->value;
  }
  p->value = value;
  p->timeindex = table->timestamp;
  table->timestamp++;
  return prev;
}

int 
Uintpairtable_length (T table) {
  assert(table);
  return table->length;
}


static int
uint_compare (const void *a, const void *b) {
  unsigned int x = * (unsigned int *) a;
  unsigned int y = * (unsigned int *) b;

  if (x < y) {
    return -1;
  } else if (y < x) {
    return 1;
  } else {
    return 0;
  }
}


void
Uintpairtable_keys (unsigned int **keys1, unsigned int **keys2,
		    T table, bool sortp, unsigned int end) {
  unsigned int *key1array, *key2array;
  int *order;
  int i, j, k;
  struct binding *p;

  assert(table);
  key1array = (unsigned int *) CALLOC(table->length+1,sizeof(unsigned int));
  key2array = (unsigned int *) CALLOC(table->length+1,sizeof(unsigned int));

  k = 0;
  for (i = 0; i < table->size; i++) {
    for (p = table->buckets[i]; p; p = p->link) {
      key1array[k] = p->key1;
      key2array[k] = p->key2;
      k++;
    }
  }

  if (sortp == false) {
    *keys1 = key1array;
    *keys2 = key2array;
  } else {
    *keys1 = (unsigned int *) CALLOC(table->length+1,sizeof(unsigned int));
    *keys2 = (unsigned int *) CALLOC(table->length+1,sizeof(unsigned int));
    order = Sedgesort_order_uint4(key1array,table->length);

    i = 0;
    while (i < table->length) {
      j = i + 1;
      while (j < table->length && key1array[order[j]] == key1array[order[i]]) {
	j++;
      }

      for (k = i; k < j; k++) {
	(*keys1)[k] = key1array[order[k]];
	(*keys2)[k] = key2array[order[k]];
      }
      if (j - i > 1) {
	qsort(&((*keys2)[i]),j - i,sizeof(unsigned int),uint_compare);
      }
	
      i = j;
    }

    FREE(order);
    FREE(key2array);
    FREE(key1array);
  }

  (*keys1)[table->length] = end;
  (*keys2)[table->length] = end;

  return;
}


void **
Uintpairtable_values (T table) {
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
Uintpairtable_free (T *table) {
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
  return;
}
