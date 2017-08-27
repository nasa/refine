
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
  RSS(ref_node_create(&(ref_layer_node(ref_layer))),"create node");
  RSS(ref_cell_create(&(ref_layer_cell(ref_layer)),6,REF_FALSE),"create cell");

  return REF_SUCCESS;
}

REF_STATUS ref_layer_free( REF_LAYER ref_layer )
{
  if ( NULL == (void *)ref_layer ) return REF_NULL;

  ref_cell_free( ref_layer_cell(ref_layer) );
  ref_node_free( ref_layer_node(ref_layer) );
  ref_list_free( ref_layer_list(ref_layer) );
  ref_free( ref_layer );

  return REF_SUCCESS;
}

REF_STATUS ref_layer_attach( REF_LAYER ref_layer,
			     REF_GRID ref_grid, REF_INT faceid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  /* copy nodes into local copy that provides compact index */
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( faceid == nodes[ref_cell_node_per(ref_cell)] )
      RSS( ref_list_add( ref_layer_list(ref_layer), cell ), "parent" );
  
  return REF_SUCCESS;
}

REF_STATUS ref_layer_puff( REF_LAYER ref_layer, REF_GRID ref_grid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_DICT node_dict;
  REF_INT item, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, local, global, i;
  
  RSS(ref_dict_create(&node_dict),"create nodes");

  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	if ( ! ref_dict_has_key( node_dict, nodes[cell_node] ) )
	  {
	    RSS( ref_dict_store( node_dict, nodes[cell_node], REF_EMPTY ),
		 "store");
	  }
    }
  printf(" layer ntri %d nnode %d\n",
	 ref_list_n(ref_layer_list(ref_layer)),
	 ref_dict_n(node_dict));

  each_ref_dict_key_value( node_dict, item, global, local )
    {
      RSS( ref_node_add(ref_layer_node(ref_layer), global, &node), "add");
      for (i=0;i<3;i++)
	ref_node_xyz(ref_layer_node(ref_layer), i, node) =
	  ref_node_xyz(ref_grid_node(ref_grid), i, global);
    }
  
  for (local = 0;local<ref_dict_n(node_dict);local++)
    {
      global = local+ref_dict_n(node_dict);
      RSS( ref_node_add(ref_layer_node(ref_layer), global, &node), "add");
      for (i=0;i<3;i++)
	ref_node_xyz(ref_layer_node(ref_layer), i, node) =
	  ref_node_xyz(ref_grid_node(ref_grid), i, global);
    }
  
  ref_dict_free(node_dict);
  
  return REF_SUCCESS;
}


