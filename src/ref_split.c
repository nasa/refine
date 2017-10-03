
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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

#include "ref_geom.h"
#include "ref_cavity.h"

#define MAX_CELL_SPLIT (100)

REF_STATUS ref_split_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *edges, *order;
  REF_INT i, n, edge;
  REF_BOOL allowed, allowed_quality, allowed_local, geom_support, valid_cavity;
  REF_INT global, new_node;
  REF_CAVITY ref_cavity = (REF_CAVITY)NULL;
  REF_STATUS status;

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
      if ( ratio[n] > ref_grid_adapt(ref_grid,split_ratio) )
	{
	  edges[n] = edge;
	  n++;
	}
    }

  RSS( ref_sort_heap_dbl( n, ratio, order), "sort lengths" );

  for ( i = n-1; i>= 0; i-- )
    {
      edge = edges[order[i]];
      
      RSS( ref_cell_has_side( ref_grid_tet(ref_grid),
			      ref_edge_e2n( ref_edge, 0, edge ),
			      ref_edge_e2n( ref_edge, 1, edge ),
			      &allowed ), "has side" );
      if ( !allowed) continue;

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
      RSS( ref_geom_add_between( ref_grid, 
				 ref_edge_e2n( ref_edge, 0, edge ),
				 ref_edge_e2n( ref_edge, 1, edge ),
				 new_node ), "geom new node");
      RSS( ref_geom_constrain( ref_grid, new_node ), "geom constraint");
      RSS( ref_geom_supported( ref_grid_geom(ref_grid), new_node,
			     &geom_support ), "geom support");

      RSS( ref_split_edge_quality( ref_grid,
				   ref_edge_e2n( ref_edge, 0, edge ),
				   ref_edge_e2n( ref_edge, 1, edge ),
				   new_node,
				   &allowed_quality ), "edge qual" );
      valid_cavity = REF_FALSE;
      if ( !allowed_quality && geom_support && ref_mpi_n == 1 )
	{
	  RSS( ref_cavity_create( &ref_cavity, 3 ), "cav create" );
	  RSS( ref_cavity_add_edge( ref_cavity, ref_grid,
				    ref_edge_e2n( ref_edge, 0, edge ),
				    ref_edge_e2n( ref_edge, 1, edge ) ),
	       "cav add" );
	  RSS( ref_cavity_split_edge( ref_cavity,
				      ref_edge_e2n( ref_edge, 0, edge ),
				      ref_edge_e2n( ref_edge, 1, edge ),
				      new_node ), "cav split" );
	  status = ref_cavity_enlarge_visible( ref_cavity, ref_grid,
					       new_node );
	  if ( REF_SUCCESS == status )
	    {
	      RSS( ref_cavity_local( ref_cavity, ref_grid,
				     new_node, &allowed_local ), "cav local");
	      valid_cavity = allowed_local;	      
	    }
	  RSS( ref_cavity_free( ref_cavity ), "cav free" );
	  ref_cavity = (REF_CAVITY)NULL;
	}
      
      if ( !valid_cavity && !allowed_quality )
	{
	  RSS( ref_node_remove( ref_node, new_node ), "remove new node");
	  RSS( ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
	  continue;
	}

      RSS( ref_split_edge_local_tets( ref_grid,
				      ref_edge_e2n( ref_edge, 0, edge ),
				      ref_edge_e2n( ref_edge, 1, edge ),
				      &allowed_local ), "local tet" );
      if ( !allowed_local )
	{
	  ref_node_age(ref_node,ref_edge_e2n( ref_edge, 0, edge ))++;
	  ref_node_age(ref_node,ref_edge_e2n( ref_edge, 1, edge ))++;
	  RSS( ref_node_remove( ref_node, new_node ), "remove new node");
	  RSS( ref_geom_remove_all(ref_grid_geom(ref_grid), new_node), "rm");
	  continue;
	}

      RSS( ref_split_edge( ref_grid,
			   ref_edge_e2n( ref_edge, 0, edge ),
			   ref_edge_e2n( ref_edge, 1, edge ),
			   new_node ), "split" );
      if ( valid_cavity )
	{
	  RSS( ref_cavity_create( &ref_cavity, 3 ), "cav create" );
	  RSS( ref_cavity_add_ball( ref_cavity, ref_grid,
				    new_node ), "cav split" );
	  RSS( ref_cavity_enlarge_visible( ref_cavity, ref_grid,
					   new_node ), "cav enlarge");
	  RSS( ref_cavity_replace_tet( ref_cavity, ref_grid,
				       new_node ), "cav replace" );
	  RSS( ref_cavity_free( ref_cavity ), "cav free" );
	  ref_cavity = (REF_CAVITY)NULL;
	}

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

  ref_cell = ref_grid_edg(ref_grid);
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

REF_STATUS ref_split_face( REF_GRID ref_grid,
			   REF_INT node0, REF_INT node1, REF_INT node2,
			   REF_INT new_node )
{
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER], face_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell0, cell1;
  REF_INT node, new_cell;

  face_nodes[0] = node0;
  face_nodes[1] = node1;
  face_nodes[2] = node2;
  face_nodes[3] = node0;

  ref_cell = ref_grid_tet(ref_grid);
  RSS( ref_cell_with_face(ref_cell,face_nodes,
			  &cell0, &cell1 ), "get tet(2)" );
  if ( REF_EMPTY != cell0 )
    {
      cell = cell0;
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
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node1;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node2 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node2 version");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node2;
    }
  if ( REF_EMPTY != cell1 )
    {
      cell = cell1;
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
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node1;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node2 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node2 version");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node2;
    }

  ref_cell = ref_grid_tri(ref_grid);
  RXS( ref_cell_with(ref_cell,face_nodes,&cell),
       REF_NOT_FOUND, "find tri" );
  if ( REF_EMPTY != cell )
    {
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
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node1;

      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( node2 == nodes[node] ) nodes[node] = new_node;
      RSS( ref_cell_add(ref_cell,nodes,&new_cell),"add node2 version");
      for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
	if ( new_node == nodes[node] ) nodes[node] = node2;
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

      if ( quality0 < ref_grid_adapt(ref_grid,split_quality_absolute) ||
	   quality1 < ref_grid_adapt(ref_grid,split_quality_absolute) ||
	   quality0 < ref_grid_adapt(ref_grid,split_quality_relative)*min_existing_quality ||
	   quality1 < ref_grid_adapt(ref_grid,split_quality_relative)*min_existing_quality ) 
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
      if ( ratio[n] > ref_grid_adapt(ref_grid,split_ratio) )
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
      RSS( ref_geom_add_between( ref_grid, node0, node1,
				 new_node0 ), "geom new node");
      RSS( ref_geom_constrain( ref_grid, new_node0 ), "geom constraint");

      RSS( ref_split_prism_tri_quality( ref_grid, node0, node1, new_node0,
					&allowed ), "quality of new tri" );
      if ( !allowed ) 
	{
	  RSS( ref_node_remove( ref_node, new_node0 ), "remove new node");
	  RSS( ref_geom_remove_all(ref_grid_geom(ref_grid), new_node0), "rm");
	  continue;
	}

      RSS( ref_split_edge_local_prisms( ref_grid, node0, node1,
					&allowed ), "local pri" );
      if ( !allowed ) 
	{
	  ref_node_age(ref_node,node0)++;
	  ref_node_age(ref_node,node1)++;
	  RSS( ref_node_remove( ref_node, new_node0 ), "remove new node");
	  RSS( ref_geom_remove_all(ref_grid_geom(ref_grid), new_node0), "rm");
	  continue;
	}

      RSS(ref_twod_opposite_edge( ref_grid_pri(ref_grid),
				  node0,node1,&node2,&node3),"opp");

      RSS( ref_node_next_global( ref_node, &global ), "next global");
      RSS( ref_node_add( ref_node, global, &new_node1 ), "new node");
      RSS( ref_node_interpolate_edge( ref_node, node2, node3,
				      new_node1 ), "interp new node");
      RSS( ref_geom_add_between( ref_grid, node2, node3,
				 new_node1 ), "geom new node");
      RSS( ref_geom_constrain( ref_grid, new_node1 ), "geom constraint");

      RSS( ref_split_twod_edge( ref_grid, node0, node1, new_node0,
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

REF_STATUS ref_split_twod_edge( REF_GRID ref_grid, 
				REF_INT node0, REF_INT node1,
				REF_INT new_node0, 
				REF_INT node2, REF_INT node3,
				REF_INT new_node1 )
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

      if ( quality0 < ref_grid_adapt(ref_grid,split_quality_absolute) ||
	   quality1 < ref_grid_adapt(ref_grid,split_quality_absolute) ||
	   quality0 < ref_grid_adapt(ref_grid,split_quality_relative)*min_existing_quality ||
	   quality1 < ref_grid_adapt(ref_grid,split_quality_relative)*min_existing_quality ) 
	return REF_SUCCESS;
    }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

