
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
#define MAX_NODE_LIST (1000)

REF_STATUS ref_collapse_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *order;
  REF_INT ntarget, *target, *node2target;
  REF_INT node, node0, node1;
  REF_INT i, edge;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL edge_ratio;
  REF_DBL ratio_limit;

  ratio_limit = 0.5 * sqrt(2.0);

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc_init( ratio, ref_node_max(ref_node), REF_DBL, 2.0*ratio_limit );

  for(edge=0;edge<ref_edge_n(ref_edge);edge++)
    {
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );
      RSS( ref_node_ratio( ref_node, node0, node1, &edge_ratio ), "ratio");
      ratio[node0] = MIN( ratio[node0], edge_ratio );
      ratio[node1] = MIN( ratio[node1], edge_ratio );
    }

  ref_malloc( target, ref_node_n(ref_node), REF_INT );
  ref_malloc_init( node2target, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  ntarget=0;
  for ( node=0 ; node < ref_node_max(ref_node) ; node++ )
    if ( ratio[node] < ratio_limit )
      {
	node2target[node] = ntarget;
	target[ntarget] = node;
	ratio[ntarget] = ratio[node];
	ntarget++;
      }

  ref_malloc( order, ntarget, REF_INT );

  RSS( ref_sort_heap_dbl( ntarget, ratio, order), "sort lengths" );

  for ( i = 0; i < ntarget; i++ )
    {
      if ( ratio[order[i]] > ratio_limit ) continue; 
      node1 = target[order[i]];
      RSS( ref_collapse_to_remove_node1( ref_grid, &node0, node1 ), 
	   "collapse rm" );
      if ( !ref_node_valid(ref_node,node1) )
	{
	  each_ref_cell_having_node( ref_cell, node1, item, cell )
	    {
	      RSS( ref_cell_nodes( ref_cell, cell, nodes), "cell nodes");
	      for (node=0;node<4;node++)
		if ( REF_EMPTY != node2target[nodes[node]] )
		  ratio[node2target[nodes[node]]] = 1.0;
	    }
	}
    }
  
  ref_free( order );
  ref_free( node2target );
  ref_free( target );
  ref_free( ratio );

  ref_edge_free( ref_edge );

  return REF_SUCCESS;
}

REF_STATUS ref_collapse_to_remove_node1( REF_GRID ref_grid, 
					 REF_INT *actual_node0, REF_INT node1 )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT nnode, node;
  REF_INT node_to_collapse[MAX_NODE_LIST];
  REF_INT order[MAX_NODE_LIST];
  REF_DBL ratio_to_collapse[MAX_NODE_LIST];
  REF_INT node0;
  REF_BOOL allowed;

  *actual_node0 = REF_EMPTY;

  RSS( ref_cell_node_list_around( ref_cell, node1, MAX_NODE_LIST,
				  &nnode, node_to_collapse ), "da hood");
  for ( node=0 ; node < nnode ; node++ )
    RSS( ref_node_ratio( ref_node, node_to_collapse[node], node1, 
			 &(ratio_to_collapse[node]) ), "ratio");

  RSS( ref_sort_heap_dbl( nnode, ratio_to_collapse, order), "sort lengths" );

  for ( node=0 ; node < nnode ; node++ )
    {
      node0 = node_to_collapse[order[node]];
  
      RSS(ref_collapse_edge_mixed(ref_grid,node0,node1,&allowed),"col mixed");
      if ( !allowed ) continue;

      RSS(ref_collapse_edge_geometry(ref_grid,node0,node1,&allowed),"col geom");
      if ( !allowed ) continue;

      RSS(ref_collapse_edge_local_tets(ref_grid,node0,node1,&allowed),"colloc");
      if ( !allowed ) continue;

      RSS(ref_collapse_edge_quality(ref_grid,node0,node1,&allowed),"qual");
      if ( !allowed ) continue;

      *actual_node0 = node0;
      RSS( ref_collapse_edge( ref_grid, node0, node1 ), "col!");

      break;

    }

  return REF_SUCCESS;
}

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
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT deg, degree1;
  REF_INT id, ids1[3]; 
  REF_BOOL already_have_it;

  REF_INT ncell;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];
  REF_INT id0, id1;

  degree1 = 0;
  each_ref_cell_having_node( ref_cell, node1, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      id = nodes[3];
      already_have_it = REF_FALSE;
      for (deg=0;deg<degree1;deg++)
	if ( id == ids1[deg] ) already_have_it = REF_TRUE;
      if ( !already_have_it )
	{
	  ids1[degree1] = id;
	  degree1++;
	  if ( 3 == degree1 ) break;
	}
    }

  *allowed = REF_FALSE;

  switch ( degree1 )
    {
    case 3: /* geometry node never allowed to move */
      *allowed = REF_FALSE;
      break;
    case 2: /* geometery edge allowed if collapse is on edge */
      RSS( ref_cell_list_with(ref_cell,node0,node1,
			      MAX_CELL_COLLAPSE, &ncell, 
			      cell_to_collapse ),"list");
      if ( 2 != ncell ) return REF_SUCCESS;
      RSS( ref_cell_nodes( ref_cell, cell_to_collapse[0], nodes ), "nodes" );
      id0 = nodes[3];
      RSS( ref_cell_nodes( ref_cell, cell_to_collapse[1], nodes ), "nodes" );
      id1 = nodes[3];
      if ( ( id0 == ids1[0] && id1 == ids1[1] ) ||
	   ( id1 == ids1[0] && id0 == ids1[1] ) ) *allowed = REF_TRUE;
      break;
    case 1: /* geometry face allowed if on that face */
      RSS( ref_cell_has_side( ref_cell, node0, node1, allowed ),
	   "allowed if a side of a triangle" );
      break;
    case 0: /* volume node always allowed */
      *allowed = REF_TRUE;
      break;
    }

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

REF_STATUS ref_collapse_edge_quality( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT node1,
				      REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL quality;
  REF_DBL quality_tolerence = 1.0e-3;
  REF_BOOL will_be_collapsed;

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node( ref_cell, node1, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      will_be_collapsed = REF_FALSE;
      for ( node = 0; node < ref_cell_node_per(ref_cell) ; node++ )
	if ( node0 == nodes[node] ) will_be_collapsed = REF_TRUE;
      if ( will_be_collapsed ) continue;
      for ( node = 0; node < ref_cell_node_per(ref_cell) ; node++ )
	if ( node1 == nodes[node] ) nodes[node] = node0;
      RSS( ref_node_tet_quality( ref_node,nodes,&quality ), "qual");
      if ( quality < quality_tolerence ) return REF_SUCCESS;
    }

  /* FIXME check tris too */

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

