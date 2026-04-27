static char rcsid[] = "$Id: 74b107acd9c0bf275440ab84a57be16e48f4b386 $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>		/* For memset */
#include "assert.h"
#include "mem.h"
#include "bitvector.h"
#include "popcount.h"

/* Need SSE4_2 or lzcnt and tzcnt */
#include "simd.h"

#if !defined(HAVE_SSE4_2)
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#elif defined(HAVE_LZCNT)
#define count_leading_zeroes_32(diff) _lzcnt_u32(diff)
#elif defined(HAVE_BUILTIN_CLZ)
#define count_leading_zeroes_32(diff) __builtin_clz(diff)
#else
#define count_leading_zeroes_32(diff) ((diff >> 16) ? clz_table[diff >> 16] : 16 + clz_table[diff])
#endif

#if !defined(HAVE_SSE4_2)
#define count_trailing_zeroes_32(elt) mod_37_bit_position[(-elt & elt) % 37]
#elif defined(HAVE_TZCNT)
#define count_trailing_zeroes_32(elt) _tzcnt_u32(elt)
#elif defined(HAVE_BUILTIN_CTZ)
#define count_trailing_zeroes_32(elt) __builtin_ctz(elt)
#else
/* lowbit = -elt & elt */
#define count_trailing_zeroes_32(elt) mod_37_bit_position[(-elt & elt) % 37]
#endif

#define clear_lowbit_32(elt,relpos) (elt & (elt - 1))


#define T Bitvector_T


void
Assoc_free (Assoc_T *old) {
  FREE(*old);
  return;
}

Assoc_T
Assoc_new (T bitvector, void *value) {
  Assoc_T new = (Assoc_T) MALLOC(sizeof(*new));

  new->bitvector = bitvector;
  new->value = value;

  return new;
}

Assoc_T
Assoc_find (T this, List_T assocs) {
  List_T p;
  Assoc_T assoc;

  for (p = assocs; p != NULL; p = List_next(p)) {
    assoc = (Assoc_T) List_head(p);
    if (Bitvector_diff_p(assoc->bitvector,this) == false) {
      return assoc;
    }
  }
  return (Assoc_T) NULL;
}


bool
Bitvector_find (T this, List_T bitvectors) {
  List_T p;

  for (p = bitvectors; p != NULL; p = List_next(p)) {
    if (Bitvector_diff_p(this,(T) List_head(p)) == false) {
      return true;
    }
  }

  return false;
}
  

double
Bitvector_double_find (T this, List_T bitvectors, Doublelist_T values) {
  List_T p;
  Doublelist_T q;

  for (p = bitvectors, q = values; p != NULL; p = List_next(p), q = Doublelist_next(q)) {
    if (Bitvector_diff_p(this,(T) List_head(p)) == false) {
      return Doublelist_head(q);
    }
  }

  fprintf(stderr,"In Bitvector_double_find, could not find bitvector\n");
  abort();
  return 0.0;
}


Doublelist_T
Bitvector_double_find_ptr (T this, List_T bitvectors, Doublelist_T values) {
  List_T p;
  Doublelist_T q;

  for (p = bitvectors, q = values; p != NULL; p = List_next(p), q = Doublelist_next(q)) {
    if (Bitvector_diff_p(this,(T) List_head(p)) == false) {
      return q;
    }
  }

  fprintf(stderr,"In Bitvector_double_find_ptr, could not find bitvector\n");
  return NULL;
}


void
Bitvector_free (T *old) {
  if (*old) {
    FREE((*old)->elts);
    FREE(*old);
  }

  return;
}


T
Bitvector_new (int nbits) {
  T new;

  if (nbits == 0) {
    return (T) NULL;
  } else {
    new = MALLOC(sizeof(*new));      

    new->nelts = (nbits + BITS_PER_ELT - 1) / BITS_PER_ELT;
    new->elts = (Bitvector_elt_T *) CALLOC(new->nelts,sizeof(Bitvector_elt_T));
    return new;
  }
}


T
Bitvector_copy (T old) {
  T new;

  if (old == NULL) {
    return (T) NULL;

  } else {
    new = (T) MALLOC(sizeof(*new));

    new->nelts = old->nelts;
    new->elts = (Bitvector_elt_T *) CALLOC(new->nelts,sizeof(Bitvector_elt_T));
    memcpy(new->elts,old->elts,old->nelts * sizeof(Bitvector_elt_T));
    return new;
  }
}


int
Bitvector_cmp (T x, T y) {
  Bitvector_elt_T *p, *q, elt_x, elt_y;
  int i;

  if (x == NULL && y == NULL) {
    return 0;
  } else if (x == NULL) {
    return -1;
  } else if (y == NULL) {
    return +1;
  } else {
    assert(x->nelts == y->nelts);
    p = &(x->elts[x->nelts]);
    q = &(y->elts[y->nelts]);
    for (i = x->nelts - 1; i >= 0; --i) {
      elt_x = *--p;
      elt_y = *--q;
      if (elt_x < elt_y) {
	return -1;
      } else if (elt_y < elt_x) {
	return +1;
      }
    }

    return 0;
  }
}  


bool
Bitvector_null_p (T this) {
  Bitvector_elt_T *p;
  int i;

  if (this != NULL) {
    p = this->elts;
    for (i = 0; i < this->nelts; i++) {
      if (*p++ /* != 0*/) {
	return false;
      }
    }
  }

  return true;
}


void
Bitvector_clear_all (T this) {
  if (this != NULL) {
    memset(this->elts,0,this->nelts*sizeof(Bitvector_elt_T));
  }

  return;
}


bool
Bitvector_intersect_p (T x, T y) {
  int i;

  if (x == NULL || y == NULL) {
    return false;
  } else {
    for (i = 0; i < x->nelts; i++) {
      if ((x->elts[i] & y->elts[i]) != 0) {
	return true;
      }
    }
    return false;
  }
}


bool
Bitvector_diff_p (T x, T y) {
  int i;

  if (x == NULL && y == NULL) {
    return false;
  } else if (x == NULL) {
    return true;
  } else if (y == NULL) {
    return true;
  } else {
    for (i = 0; i < x->nelts; i++) {
      if ((x->elts[i] ^ y->elts[i]) != 0) {
	return true;
      }
    }
    return false;
  }
}


T
Bitvector_union (T dest, T source, bool copyp) {
  int i;

  if (source == NULL) {
    if (copyp == true) {
      return Bitvector_copy(dest);
    } else {
      return dest;
    }

  } else if (dest == NULL) {
    return Bitvector_copy(source);

  } else {
    if (copyp == true) {
      dest = Bitvector_copy(dest);
    }
    for (i = 0; i < dest->nelts; i++) {
      dest->elts[i] |= source->elts[i];
    }
    return dest;
  }
}


unsigned int
Bitvector_high_bit (T this) {
  Bitvector_elt_T *p, elt;
  unsigned int offset;
  int i;
  
  offset = this->nelts * BITS_PER_ELT;
  p = &(this->elts[this->nelts]);
  for (i = this->nelts - 1; i >= 0; --i) {
    elt = *--p;
    if (elt) {
      return offset - count_leading_zeroes_32(elt);
    }
    offset -= BITS_PER_ELT;
  }

  return (unsigned int) -1;		 
}


unsigned int
Bitvector_low_bit (T this) {
  Bitvector_elt_T *p, elt;
  unsigned int offset;
  int i;
  
  offset = 0;
  p = this->elts;
  for (i = 0; i < this->nelts; i++) {
    elt = *p++;
    if (elt) {
      return offset + count_trailing_zeroes_32(elt);
    }
    offset += BITS_PER_ELT;
  }

  return (unsigned int) -1;		 
}


Uintlist_T
Bitvector_to_uintlist (T this) {
  Uintlist_T list = NULL;
  Bitvector_elt_T *p, elt;
  unsigned int offset, relpos;
  int i;
  
  if (this != NULL) {
    p = this->elts;
    offset = 0;

    for (i = 0; i < this->nelts; i++) {
      elt = *p++;
      while (elt) {
	list = Uintlist_push(list,offset + (relpos = count_trailing_zeroes_32(elt)));
	elt = clear_lowbit_32(elt,relpos);
      }
      offset += BITS_PER_ELT;
    }
  }

  return list;
}


/* Prints indices as 1-based */
void
Bitvector_print (Filestring_T fp, T this, char sep) {
  Bitvector_elt_T *p, elt;
  unsigned int offset, relpos;
  int i;
  bool firstp = true;
  
  if (this != NULL) {
    p = this->elts;
    offset = 1;			/* 1-based */
    for (i = 0; i < this->nelts; i++) {
      elt = *p++;
      while (elt) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  PUTC(sep,fp);
	}
	FPRINTF(fp,"%u",offset + (relpos = count_trailing_zeroes_32(elt)));
	elt = clear_lowbit_32(elt,relpos);
      }
      offset += BITS_PER_ELT;
    }
  }

  return;
}


void
Bitvector_print_file (FILE *out, T this, char sep) {
  Bitvector_elt_T *p, elt;
  unsigned int offset, relpos;
  int i;
  bool firstp = true;
  
  if (this != NULL) {
    p = this->elts;
    offset = 1;			/* 1-based */
    for (i = 0; i < this->nelts; i++) {
      elt = *p++;
      while (elt) {
	if (firstp == true) {
	  firstp = false;
	} else {
	  fputc(sep,out);
	}
	fprintf(out,"%u",offset + (relpos = count_trailing_zeroes_32(elt)));
	elt = clear_lowbit_32(elt,relpos);
      }
      offset += BITS_PER_ELT;
    }
  }

  return;
}


#if !defined(HAVE_INLINE)
void
Bitvector_set (T this, unsigned int i) {
  this->elts[i / BITS_PER_ELT] |= (1 << (i % BITS_PER_ELT));
}

void
Bitvector_clear (Bitvector_T this, unsigned int i) {
  this->elts[i / BITS_PER_ELT] &= ~(1 << (i % BITS_PER_ELT));
}

bool
Bitvector_get (Bitvector_T this, unsigned int i) {
  return this->elts[i / BITS_PER_ELT] & (1 << (i % BITS_PER_ELT)) ? true : false;
}
#endif


