
#include <stdlib.h>
#include <stdio.h>

#include "ref_front.h"
#include "ref_malloc.h"

REF_STATUS ref_front_create( REF_FRONT *ref_front_ptr, REF_INT node_per )
{
  REF_FRONT ref_front;

  ref_malloc( *ref_front_ptr, 1, REF_FRONT_STRUCT );
  ref_front = (*ref_front_ptr);

  ref_front_n(ref_front) = 0;
  ref_front_node_per(ref_front) = node_per;

  return REF_SUCCESS;
}

REF_STATUS ref_front_free( REF_FRONT ref_front )
{
  if ( NULL == (void *)ref_front ) return REF_NULL;
  ref_free( ref_front );
  return REF_SUCCESS;
}

