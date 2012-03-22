
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_collapse.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_mpi.h"
#include "ref_sort.h"
#include "ref_malloc.h"

#define MAX_CELL_COLLAPSE (100)

REF_STATUS ref_collapse_edge( REF_GRID ref_grid, 
			      REF_INT node0, REF_INT node1 )
{
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];

  ref_cell = ref_grid_tet(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_COLLAPSE, &ncell, cell_to_collapse ),"list");

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_collapse[cell_in_list];
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );
    }
  RSS( ref_cell_replace_node( ref_cell, node1, node0 ), "replace node" );

  ref_cell = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_COLLAPSE, &ncell, cell_to_collapse ),"list");

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_collapse[cell_in_list];
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );
    }
  RSS( ref_cell_replace_node( ref_cell, node1, node0 ), "replace node" );

  RSS( ref_node_remove(ref_grid_node(ref_grid),node1), "rm");

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_geometry( REF_GRID ref_grid, 
				       REF_INT node0, REF_INT node1,
				       REF_BOOL *allowed )
{

  SUPRESS_UNUSED_COMPILER_WARNING(node0);

  *allowed = ( ref_adj_empty( ref_cell_adj(ref_grid_tri(ref_grid)),
			      node1 ) &&
	       ref_adj_empty( ref_cell_adj(ref_grid_qua(ref_grid)),
			      node1 ) );

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_mixed( REF_GRID ref_grid, 
				    REF_INT node0, REF_INT node1,
				    REF_BOOL *allowed )
{
  SUPRESS_UNUSED_COMPILER_WARNING(node0);

  *allowed = ( ref_adj_empty( ref_cell_adj(ref_grid_pyr(ref_grid)),
			      node1 ) &&
	       ref_adj_empty( ref_cell_adj(ref_grid_pri(ref_grid)),
			      node1 ) &&
	       ref_adj_empty( ref_cell_adj(ref_grid_hex(ref_grid)),
			      node1 ) );

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_edge_local_tets( REF_GRID ref_grid, 
					 REF_INT node0, REF_INT node1,
					 REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, node;

  *allowed =  REF_FALSE;

  ref_cell = ref_grid_tet(ref_grid);

  each_ref_cell_having_node( ref_cell, node1, item, cell )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
      if ( ref_mpi_id != ref_node_part(ref_node,
				       ref_cell_c2n(ref_cell,node,cell)) )
	return REF_SUCCESS;

  /* may be able to relax node0 local if geom constraint is o.k. */
  each_ref_cell_having_node( ref_cell, node0, item, cell )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
      if ( ref_mpi_id != ref_node_part(ref_node,
				       ref_cell_c2n(ref_cell,node,cell)) )
	return REF_SUCCESS;



  *allowed =  REF_TRUE;

  return REF_SUCCESS;
}
