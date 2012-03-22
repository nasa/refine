
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
