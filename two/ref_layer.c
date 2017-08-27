
#include <stdlib.h>
#include <stdio.h>

#include "ref_layer.h"

#include "ref_malloc.h"
#include "ref_mpi.h"

REF_STATUS ref_layer_create( REF_LAYER *ref_layer_ptr )
{
  REF_LAYER ref_layer;

  ref_malloc( *ref_layer_ptr, 1, REF_LAYER_STRUCT );

  ref_layer = *ref_layer_ptr;

  ref_layer_n(ref_layer) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_layer_free( REF_LAYER ref_layer )
{
  if ( NULL == (void *)ref_layer ) return REF_NULL;

  ref_free( ref_layer );

  return REF_SUCCESS;
}

