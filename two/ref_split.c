
#include <stdlib.h>
#include <stdio.h>

#include "ref_split.h"
#include "ref_cell.h"

#define MAX_CELL_SPLIT (100)

REF_STATUS ref_split_edge( REF_GRID ref_grid, 
			   REF_INT node0, REF_INT node1,
			   REF_INT new_node )
{
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_split[MAX_CELL_SPLIT];
  REF_INT node, new_cell;

  ref_cell = ref_grid_tet(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_SPLIT, &ncell, cell_to_split ), "get list" );

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node0 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node0 version");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node0;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node1 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node1 version");
    }

  ref_cell = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_SPLIT, &ncell, cell_to_split ), "get list" );

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node0 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node0 version");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node0;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node1 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node1 version");
    }

  return REF_SUCCESS;
}
