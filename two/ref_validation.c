
#include <stdlib.h>
#include <stdio.h>

#include "ref_validation.h"

REF_STATUS ref_validation_validate( REF_GRID ref_grid )
{
  REF_INT node;
  REF_BOOL problem;
  REF_ADJ ref_adj;
  REF_CELL ref_cell;
  REF_INT group, cell, node_per, *nodes;

  problem = REF_FALSE;
  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      if ( ref_adj_empty( ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_hex(ref_grid)), node) )
	{
	  problem = REF_TRUE;
	  printf(" hanging node %d\n",node);
	}
    }

  ref_adj_create( &ref_adj );

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      node_per = ref_cell_node_per(ref_cell);
      nodes = (REF_INT *) malloc( node_per * sizeof(REF_INT) );
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  for ( node = 0; node < node_per; node++ )
	    RSS( ref_adj_add( ref_adj, nodes[node], group+4*cell ), "add");
	}
      free(nodes);
    }

  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      if ( ref_adj_empty( ref_adj, node ) )
	{
	  problem = REF_TRUE;
	  printf(" hanging node %d\n",node);	  
	}
    }

  ref_adj_free(ref_adj);

  return (problem?REF_FAILURE:REF_SUCCESS);
}

