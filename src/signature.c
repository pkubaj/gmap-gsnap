#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "signature.h"

#include <stdio.h>
#include <stdint.h>
#include "assert.h"
#include "bool.h"
#include "oligo.h"


/* Signatures */
#ifdef DEBUG9
#define debug9(x) x
#else
#define debug9(x)
#endif


static Width_T ksize;
static Width_T msize;
static Oligospace_T mmer_mask;


void
Signature_setup (Width_T ksize_in, Width_T msize_in) {

  ksize = ksize_in;
  msize = msize_in;
  mmer_mask = ~(~ (Oligospace_T) 0U << (2*msize_in));

  return;
}


/* Signature_T must be 64 bits, to handle a 15-mer * more than 4 positions */

Signature_T
Signature_get (Oligospace_T kmer) {
  Signature_T signature;
  Oligospace_T original_kmer = kmer;
  Oligospace_T index, best_mmer = 0, mmer;
  Oligospace_T top, bottom, bottom_mask;
  int mmer_i = 0, i, npos = ksize - msize + 1;

  debug9(printf("Finding signature for %s\n",Oligo_one_nt(kmer,ksize)));

  /* First mmer is the rightmost one, in the least significant bits */
  mmer_i = 0;
  best_mmer = (uint32_t) (kmer & mmer_mask);
  debug9(printf("Testing at %d: best_mmer %s %lu\n",mmer_i,Oligo_one_nt(best_mmer,msize),best_mmer));
  kmer >>= 2;

  for (i = 1; i < npos; i++) {
    mmer = (uint32_t) (kmer & mmer_mask);
    debug9(printf("Testing at %d: mmer %s %lu",i,Oligo_one_nt(mmer,msize),mmer));
    if (mmer < best_mmer) {
      debug9(printf(" => best mmer so far\n"));
      best_mmer = mmer;
      mmer_i = i;
    } else {
      debug9(printf("\n"));
    }
    kmer >>= 2;
  }

  debug9(printf("minimizer position %d\n",mmer_i));

  top = (original_kmer >> (2*(mmer_i + msize))) << (2*mmer_i);
    
  bottom_mask = ~(~ (Oligospace_T) 0U << (2*mmer_i));
  bottom = original_kmer & bottom_mask;

  index = best_mmer*npos + mmer_i;
  signature = ((Signature_T) index << 2*(ksize - msize)) | (Signature_T) top | (Signature_T) bottom;

  debug9(printf("index %u * %u + %u = %u => %16lX, top %08X, bottom %08X, signature %16lX (%lu)\n",
		best_mmer,npos,mmer_i,index,(Signature_T) index << 2*(ksize - msize),
		top,bottom,signature,signature));

  return signature;
}


Signature_T
Signature_revise_right (Oligospace_T kmer, Signature_T prev_signature) {
  Signature_T signature;
  Oligospace_T index, prev_best_mmer, new_mmer;
  Oligospace_T top, bottom, bottom_mask;
  int mmer_i, npos = ksize - msize + 1;

  debug9(printf("Revising signature right for %s from prev_signature %016lX (%lu)\n",
		Oligo_one_nt(kmer,ksize),prev_signature,prev_signature));

  prev_signature >>= 2*(ksize - msize); /* now prev_index */
  debug9(printf("prev_index is %u, which should yield a remainder of %d (mod %d)\n",
		prev_signature,prev_signature % npos,npos));

  if ((mmer_i = prev_signature % npos) == npos - 1) {
    return Signature_get(kmer);

  } else {
    new_mmer = (uint32_t) (kmer & mmer_mask);
    prev_best_mmer = prev_signature/npos;

    if (new_mmer <= prev_best_mmer) { /* Use <= when revising right to handle ties properly */
      /* mmer_i = 0; */
      index = new_mmer*npos /*+ mmer_i 0*/;

      top = kmer >> (2*(/*mmer_i 0 +*/ msize));
      bottom = 0U;

    } else {
      mmer_i += 1;
      index = prev_signature + 1;

      top = (kmer >> (2*(mmer_i + msize))) << (2*mmer_i);
    
      bottom_mask = ~(~ (Oligospace_T) 0U << (2*mmer_i));
      bottom = kmer & bottom_mask;

    }

    debug9(printf("minimizer position %d\n",mmer_i));

    signature = ((Signature_T) index << 2*(ksize - msize)) | (Signature_T) top | (Signature_T) bottom;

    debug9(printf("index %u => %lu, top %08X, bottom %08X, signature %16lX (%lu)\n",
		  index,(Signature_T) index << 2*(ksize - msize),
		  top,bottom,signature,signature));

    assert(signature == Signature_get(kmer));
    return signature;
  }
}


Signature_T
Signature_revise_left (Oligospace_T kmer, Signature_T prev_signature) {
  Signature_T signature;
  Oligospace_T index, prev_best_mmer, new_mmer;
  Oligospace_T top, bottom, bottom_mask;
  int mmer_i, npos = ksize - msize + 1;

  debug9(printf("Revising signature left for %s from prev_signature %016lX (%lu)\n",
		Oligo_one_nt(kmer,ksize),prev_signature,prev_signature));

  prev_signature >>= 2*(ksize - msize); /* now prev_index */
  debug9(printf("prev_index is %u, which should yield a remainder of %d (mod %d)\n",
		prev_signature,prev_signature % npos,npos));

  if ((mmer_i = prev_signature % npos) == 0) {
    return Signature_get(kmer);

  } else {
    new_mmer = (uint32_t) ((kmer >> 2*(ksize - msize)) & mmer_mask);
    prev_best_mmer = prev_signature/npos;

    if (new_mmer < prev_best_mmer) { /* Use < when revising left to handle ties properly */
      mmer_i = npos - 1;
      index = new_mmer*npos + mmer_i;

      top = 0U;
    
    } else {
      mmer_i -= 1;
      index = prev_signature - 1;

      top = (kmer >> (2*(mmer_i + msize))) << (2*mmer_i);
    }

    bottom_mask = ~(~ (Oligospace_T) 0U << (2*mmer_i));
    bottom = kmer & bottom_mask;

    debug9(printf("minimizer position %d\n",mmer_i));

    signature = ((Signature_T) index << 2*(ksize - msize)) | (Signature_T) top | (Signature_T) bottom;

    debug9(printf("index %u => %lu, top %08X, bottom %08X, signature %16lX (%lu)\n",
		  index,(Signature_T) index << 2*(ksize - msize),
		  top,bottom,signature,signature));

    assert(signature == Signature_get(kmer));
    return signature;
  }
}


