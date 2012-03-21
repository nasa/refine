
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_split.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_mpi.h"
#include "ref_sort.h"
#include "ref_malloc.h"

#define MAX_CELL_SPLIT (100)

REF_STATUS ref_split_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *edges, *order;
  REF_INT i, n, edge;
  REF_DBL ratio_limit;
  REF_BOOL allowed;
  REF_INT global, new_node;

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc( ratio, ref_edge_n(ref_edge), REF_DBL );
  ref_malloc( order, ref_edge_n(ref_edge), REF_INT );
  ref_malloc( edges, ref_edge_n(ref_edge), REF_INT );
  
  ratio_limit = sqrt(2.0);

  n=0;
  for(edge=0;edge<ref_edge_n(ref_edge);edge++)
    {
      RSS( ref_node_ratio( ref_node, 
			   ref_edge_e2n( ref_edge, edge, 0 ),
			   ref_edge_e2n( ref_edge, edge, 1 ),
			   &(ratio[n]) ), "ratio");
      if ( ratio[n] > ratio_limit)
	{
	  edges[n] = edge;
	  n++;
	}
    }

  RSS( ref_sort_heap_dbl( n, ratio, order), "sort lengths" );

  for ( i = n-1; i>= 0; i-- )
    {
      edge = edges[i];
      RSS( ref_split_edge_mixed( ref_grid,
				 ref_edge_e2n( ref_edge, edge, 0 ),
				 ref_edge_e2n( ref_edge, edge, 1 ),
				 &allowed ), "mixed" );
      if ( !allowed) continue;
      RSS( ref_split_edge_local_tets( ref_grid,
				      ref_edge_e2n( ref_edge, edge, 0 ),
				      ref_edge_e2n( ref_edge, edge, 1 ),
				      &allowed ), "local tet" );
      if ( !allowed) continue;

      RSS( ref_node_next_global( ref_node, &global ), "next global");
      RSS( ref_node_add( ref_node, global, &new_node ), "new node");
      RSS( ref_node_interpolate_edge( ref_node, 
				      ref_edge_e2n( ref_edge, edge, 0 ),
				      ref_edge_e2n( ref_edge, edge, 1 ),
				      new_node ), "new node");

      /* test potential quality */

      RSS( ref_split_edge( ref_grid,
			   ref_edge_e2n( ref_edge, edge, 0 ),
			   ref_edge_e2n( ref_edge, edge, 1 ),
			   new_node ), "split" );

    }
  
  ref_free( edges );
  ref_free( order );
  ref_free( ratio );

  ref_edge_free( ref_edge );

  return REF_SUCCESS;
}

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
