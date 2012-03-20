
#include <stdlib.h>
#include <stdio.h>

#include "ref_split.h"
#include "ref_cell.h"
#include "ref_mpi.h"

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

REF_STATUS ref_split_edge_mixed( REF_GRID ref_grid, 
				 REF_INT node0, REF_INT node1,
				 REF_BOOL *allowed )
{
  REF_BOOL pyr_side, pri_side, hex_side;

  RSS(ref_cell_has_side(ref_grid_pyr(ref_grid), node0, node1, &pyr_side),"pyr");
  RSS(ref_cell_has_side(ref_grid_pri(ref_grid), node0, node1, &pri_side),"pri");
  RSS(ref_cell_has_side(ref_grid_hex(ref_grid), node0, node1, &hex_side),"hex");

  *allowed = ( !pyr_side && !pri_side && !hex_side );

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_local_tets( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT node1,
				      REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, search_node, test_node;

  *allowed = REF_TRUE;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node( ref_cell, node0, item, cell )
    for ( search_node = 0 ; 
	  search_node < ref_cell_node_per(ref_cell); 
	  search_node++ )
      if ( node1 == ref_cell_c2n(ref_cell,search_node,cell) )
	for ( test_node = 0 ; 
	      test_node < ref_cell_node_per(ref_cell); 
	      test_node++ )
	  if ( ref_mpi_id != ref_node_part(ref_node,test_node) )
	    {
	      *allowed = REF_FALSE;
	      return REF_SUCCESS;
	    }

  return REF_SUCCESS;
}
