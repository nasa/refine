
#ifndef REF_MALLOC_H
#define REF_MALLOC_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

#define ref_malloc( ptr, n, ptr_type )			\
  (ptr) = (ptr_type *)malloc( (n) * sizeof(ptr_type) );	\
  RNS((ptr),"malloc " #ptr " NULL");

/* realloc of size zero with return NULL */

#define ref_realloc( ptr, n, ptr_type )					\
  if ( 0 < (n))								\
    (ptr) = (ptr_type *)realloc( (ptr), (n) * sizeof(ptr_type) );	\
  RNS((ptr),"realloc " #ptr " NULL");

#define ref_free(ptr) if ( NULL != (ptr) ) free((ptr));

END_C_DECLORATION

#endif /* REF_MALLOC_H */
