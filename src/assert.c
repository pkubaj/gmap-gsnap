static char rcsid[] = "$Id: assert.c 223349 2020-10-28 02:49:25Z twu $";
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "assert.h"

const Except_T Assert_Failed = { "Assertion Failed" };

/*
void
(assert) (int e) {
  assert(e);
}
*/

