
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

  RSS(ref_list_create(&(ref_layer_list(ref_layer))),"create list");

  return REF_SUCCESS;
}

REF_STATUS ref_layer_free( REF_LAYER ref_layer )
{
  if ( NULL == (void *)ref_layer ) return REF_NULL;

  ref_list_free( ref_layer_list(ref_layer) );
  ref_free( ref_layer );

  return REF_SUCCESS;
}

REF_STATUS ref_layer_attach( REF_LAYER ref_layer,
			     REF_GRID ref_grid, REF_DICT faceids )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell;
  
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      if (ref_dict_has_key(faceids,ref_cell_c2n( ref_cell, 3, cell )))
	{
	  RSS( ref_list_add(ref_layer_list(ref_layer), cell ), "mark tri" );
	}
    }
  
  return REF_SUCCESS;
}


