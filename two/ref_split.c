
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_split.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_mpi.h"
#include "ref_sort.h"
#include "ref_malloc.h"

#include "ref_adapt.h"

#include "ref_gather.h"
#include "ref_twod.h"

#define MAX_CELL_SPLIT (100)

REF_STATUS ref_split_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *edges, *order;
  REF_INT i, n, edge;
  REF_BOOL allowed;
  REF_INT global, new_node;

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc( ratio, ref_edge_n(ref_edge), REF_DBL );
  ref_malloc( order, ref_edge_n(ref_edge), REF_INT );
  ref_malloc( edges, ref_edge_n(ref_edge), REF_INT );
  
  n=0;
  for(edge=0;edge<ref_edge_n(ref_edge);edge++)
    {
      RSS( ref_node_ratio( ref_node, 
			   ref_edge_e2n( ref_edge, 0, edge ),
			   ref_edge_e2n( ref_edge, 1, edge ),
			   &(ratio[n]) ), "ratio");
      if ( ratio[n] > ref_adapt_split_ratio )
	{
	  edges[n] = edge;
	  n++;
	}
    }

  RSS( ref_sort_heap_dbl( n, ratio, order), "sort lengths" );

  for ( i = n-1; i>= 0; i-- )
    {
      edge = edges[order[i]];
      RSS( ref_split_edge_mixed( ref_grid,
				 ref_edge_e2n( ref_edge, 0, edge ),
				 ref_edge_e2n( ref_edge, 1, edge ),
				 &allowed ), "mixed" );
      if ( !allowed) continue;

      RSS( ref_node_next_global( ref_node, &global ), "next global");
      RSS( ref_node_add( ref_node, global, &new_node ), "new node");
      RSS( ref_node_interpolate_edge( ref_node, 
				      ref_edge_e2n( ref_edge, 0, edge ),
				      ref_edge_e2n( ref_edge, 1, edge ),
				      new_node ), "interp new node");

      RSS( ref_split_edge_quality( ref_grid,
				   ref_edge_e2n( ref_edge, 0, edge ),
				   ref_edge_e2n( ref_edge, 1, edge ),
				   new_node,
				   &allowed ), "edge qual" );
      if ( !allowed) 
	{
	  RSS( ref_node_remove( ref_node, new_node ), "remove new node");
	  continue;
	}

      RSS( ref_split_edge_local_tets( ref_grid,
				      ref_edge_e2n( ref_edge, 0, edge ),
				      ref_edge_e2n( ref_edge, 1, edge ),
				      &allowed ), "local tet" );
      if ( !allowed) 
	{
	  ref_node_age(ref_node,ref_edge_e2n( ref_edge, 0, edge ))++;
	  ref_node_age(ref_node,ref_edge_e2n( ref_edge, 1, edge ))++;
	  RSS( ref_node_remove( ref_node, new_node ), "remove new node");
	  continue;
	}

      RSS( ref_split_edge( ref_grid,
			   ref_edge_e2n( ref_edge, 0, edge ),
			   ref_edge_e2n( ref_edge, 1, edge ),
			   new_node ), "split" );

      ref_node_age(ref_node,ref_edge_e2n( ref_edge, 0, edge )) = 0;
      ref_node_age(ref_node,ref_edge_e2n( ref_edge, 1, edge )) = 0;

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
  RSS( ref_cell_list_with2(ref_cell,node0,node1,
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
  RSS( ref_cell_list_with2(ref_cell,node0,node1,
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

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node( ref_cell, node0, item, cell )
    for ( search_node = 0 ; 
	  search_node < ref_cell_node_per(ref_cell); 
	  search_node++ )
      if ( node1 == ref_cell_c2n(ref_cell,search_node,cell) )
	for ( test_node = 0 ; 
	      test_node < ref_cell_node_per(ref_cell); 
	      test_node++ )
	  if ( ref_mpi_id != ref_node_part( ref_node,
					    ref_cell_c2n(ref_cell,
							 test_node,cell) ) )
	    {
	      return REF_SUCCESS;
	    }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_quality( REF_GRID ref_grid, 
				   REF_INT node0, REF_INT node1,
				   REF_INT new_node,
				   REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_split[MAX_CELL_SPLIT];
  REF_INT node;
  REF_DBL quality, quality0, quality1;
  REF_DBL min_existing_quality;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tet(ref_grid);
  RSS( ref_cell_list_with2(ref_cell,node0,node1,
			  MAX_CELL_SPLIT, &ncell, cell_to_split ), "get list" );

  min_existing_quality = 1.0;
  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");

      RSS( ref_node_tet_quality( ref_node,nodes,&quality ), "q");
      min_existing_quality = MIN( min_existing_quality, quality );
    }

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node0 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_node_tet_quality( ref_node,nodes,&quality0 ), "q0");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node0;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node1 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_node_tet_quality( ref_node,nodes,&quality1 ), "q1");

      if ( quality0 < ref_adapt_split_quality_absolute ||
	   quality1 < ref_adapt_split_quality_absolute ||
	   quality0 < ref_adapt_split_quality_relative*min_existing_quality ||
	   quality1 < ref_adapt_split_quality_relative*min_existing_quality ) 
	return REF_SUCCESS;
    }

  /* FIXME check tris too */

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_twod_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *edges, *order;
  REF_INT edge, n, i;
  REF_BOOL active, allowed;
  REF_INT node0, node1, node2, node3, new_node0, new_node1;
  REF_INT global;

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc( ratio, ref_edge_n(ref_edge), REF_DBL );
  ref_malloc( order, ref_edge_n(ref_edge), REF_INT );
  ref_malloc( edges, ref_edge_n(ref_edge), REF_INT );
  
  n=0;
  for(edge=0;edge<ref_edge_n(ref_edge);edge++)
    {
      RSS( ref_node_edge_twod( ref_node, 
			       ref_edge_e2n( ref_edge, 0, edge ),
			       ref_edge_e2n( ref_edge, 1, edge ),
			       &active ), "act" );
      if ( !active ) continue;

      RSS( ref_node_ratio( ref_node, 
			   ref_edge_e2n( ref_edge, 0, edge ),
			   ref_edge_e2n( ref_edge, 1, edge ),
			   &(ratio[n]) ), "ratio");
      if ( ratio[n] > ref_adapt_split_ratio )
	{
	  edges[n] = edge;
	  n++;
	}
    }

  RSS( ref_sort_heap_dbl( n, ratio, order), "sort lengths" );

  for ( i = n-1; i>= 0; i-- )
    {
      edge = edges[order[i]];
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );

      RSS( ref_node_next_global( ref_node, &global ), "next global");
      RSS( ref_node_add( ref_node, global, &new_node0 ), "new node");
      RSS( ref_node_interpolate_edge( ref_node, node0, node1,
				      new_node0 ), "interp new node");

      RSS( ref_split_prism_tri_quality( ref_grid, node0, node1, new_node0,
					&allowed ), "quality of new tri" );
      if ( !allowed ) 
	{
	  RSS( ref_node_remove( ref_node, new_node0 ), "remove new node");
	  continue;
	}

      RSS( ref_split_edge_local_prisms( ref_grid, node0, node1,
					&allowed ), "local pri" );
      if ( !allowed ) 
	{
	  ref_node_age(ref_node,node0)++;
	  ref_node_age(ref_node,node1)++;
	  RSS( ref_node_remove( ref_node, new_node0 ), "remove new node");
	  continue;
	}

      RSS(ref_twod_opposite_edge( ref_grid_pri(ref_grid),
				  node0,node1,&node2,&node3),"opp");

      RSS( ref_node_next_global( ref_node, &global ), "next global");
      RSS( ref_node_add( ref_node, global, &new_node1 ), "new node");
      RSS( ref_node_interpolate_edge( ref_node, node2, node3,
				      new_node1 ), "interp new node");

      RSS( ref_split_face( ref_grid, node0, node1, new_node0,
			   node2, node3, new_node1 ), "split face");

      ref_node_age(ref_node,node0) = 0;
      ref_node_age(ref_node,node1) = 0;

    }

  ref_free( edges );
  ref_free( order );
  ref_free( ratio );

  ref_edge_free( ref_edge );

  return REF_SUCCESS;
}

REF_STATUS ref_split_face( REF_GRID ref_grid, 
			   REF_INT node0, REF_INT node1, REF_INT new_node0, 
			   REF_INT node2, REF_INT node3, REF_INT new_node1 )
{

  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL qua = ref_grid_qua(ref_grid);

  REF_INT ncell, cell_in_list;
  REF_INT cell_to_split[2];
  REF_INT orig_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell, new_cell, node;

  RSS( ref_cell_list_with2( pri, node0, node1,
			   2, &ncell, cell_to_split), "more than two" );

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes( pri, cell, orig_nodes),"cell nodes");
      RSS( ref_cell_remove( pri, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(pri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node0 == orig_nodes[node] ) new_nodes[node] = new_node0;
	  if ( node2 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      RSS( ref_cell_add(pri,new_nodes,&new_cell),
	   "add node0-node2 version");

      for ( node = 0 ; node < ref_cell_node_per(pri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node1 == orig_nodes[node] ) new_nodes[node] = new_node0;
	  if ( node3 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      RSS( ref_cell_add(pri,new_nodes,&new_cell),
	   "add node1-node3 version");
    }

  tri = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with2(tri,node0,node1,
			  2, &ncell, cell_to_split ), "more then two" );

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(tri, cell, orig_nodes),"cell nodes");
      RSS( ref_cell_remove( tri, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(tri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node0 == orig_nodes[node] ) new_nodes[node] = new_node0;
	}
      new_nodes[ref_cell_node_per(tri)] =  orig_nodes[ref_cell_node_per(tri)];
      RSS( ref_cell_add(tri,new_nodes,&new_cell),"add node0 version");

      for ( node = 0 ; node < ref_cell_node_per(tri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node1 == orig_nodes[node] ) new_nodes[node] = new_node0;
	}
      new_nodes[ref_cell_node_per(tri)] =  orig_nodes[ref_cell_node_per(tri)];
      RSS( ref_cell_add(tri,new_nodes,&new_cell),"add node1 version");

    }

  tri = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with2(tri,node2,node3,
			  2, &ncell, cell_to_split ), "more then two" );

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(tri, cell, orig_nodes),"cell nodes");
      RSS( ref_cell_remove( tri, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(tri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node2 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      new_nodes[ref_cell_node_per(tri)] =  orig_nodes[ref_cell_node_per(tri)];
      RSS( ref_cell_add(tri,new_nodes,&new_cell),"add node0 version");

      for ( node = 0 ; node < ref_cell_node_per(tri); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node3 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      new_nodes[ref_cell_node_per(tri)] =  orig_nodes[ref_cell_node_per(tri)];
      RSS( ref_cell_add(tri,new_nodes,&new_cell),"add node1 version");

    }

  qua = ref_grid_qua(ref_grid);
  orig_nodes[0] = node0;
  orig_nodes[1] = node1;
  orig_nodes[2] = node3;
  orig_nodes[3] = node2;
  RXS( ref_cell_with(qua, orig_nodes, &cell ), REF_NOT_FOUND, "qua with" );

  if ( REF_EMPTY != cell )
    {
      RSS( ref_cell_nodes(qua, cell, orig_nodes),"cell nodes");
      RSS( ref_cell_remove( qua, cell ), "remove" );

      for ( node = 0 ; node < ref_cell_node_per(qua); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node0 == orig_nodes[node] ) new_nodes[node] = new_node0;
	  if ( node2 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      new_nodes[ref_cell_node_per(qua)] =  orig_nodes[ref_cell_node_per(qua)];
      RSS( ref_cell_add(qua,new_nodes,&new_cell),"add node0-node3 version");

      for ( node = 0 ; node < ref_cell_node_per(qua); node++ )
	{
	  new_nodes[node] = orig_nodes[node];
	  if ( node1 == orig_nodes[node] ) new_nodes[node] = new_node0;
	  if ( node3 == orig_nodes[node] ) new_nodes[node] = new_node1;
	}
      new_nodes[ref_cell_node_per(qua)] =  orig_nodes[ref_cell_node_per(qua)];
      RSS( ref_cell_add(qua,new_nodes,&new_cell),"add node1-node2 version");

    }

  return REF_SUCCESS;
}

REF_STATUS ref_split_edge_local_prisms( REF_GRID ref_grid, 
					REF_INT node0, REF_INT node1,
					REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, search_node, test_node;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_having_node( ref_cell, node0, item, cell )
    for ( search_node = 0 ; 
	  search_node < ref_cell_node_per(ref_cell); 
	  search_node++ )
      if ( node1 == ref_cell_c2n(ref_cell,search_node,cell) )
	for ( test_node = 0 ; 
	      test_node < ref_cell_node_per(ref_cell); 
	      test_node++ )
	  if ( ref_mpi_id != ref_node_part( ref_node,
					    ref_cell_c2n(ref_cell,
							 test_node,cell) ) )
	    {
	      return REF_SUCCESS;
	    }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_split_prism_tri_quality( REF_GRID ref_grid, 
					REF_INT node0, REF_INT node1,
					REF_INT new_node,
					REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_split[MAX_CELL_SPLIT];
  REF_INT node;
  REF_DBL quality, quality0, quality1;
  REF_DBL min_existing_quality;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with2(ref_cell,node0,node1,
			  MAX_CELL_SPLIT, &ncell, cell_to_split ), "get list" );

  min_existing_quality = 1.0;
  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");
      
      RSS( ref_node_tri_quality( ref_node,nodes,&quality ), "q");
      min_existing_quality = MIN( min_existing_quality, quality );
    }

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_split[cell_in_list];
      RSS( ref_cell_nodes(ref_cell, cell, nodes),"cell nodes");

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node0 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_node_tri_quality( ref_node,nodes,&quality0 ), "q0");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node0;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node1 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_node_tri_quality( ref_node,nodes,&quality1 ), "q1");

      if ( quality0 < ref_adapt_split_quality_absolute ||
	   quality1 < ref_adapt_split_quality_absolute ||
	   quality0 < ref_adapt_split_quality_relative*min_existing_quality ||
	   quality1 < ref_adapt_split_quality_relative*min_existing_quality ) 
	return REF_SUCCESS;
    }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

