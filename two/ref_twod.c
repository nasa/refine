
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_twod.h"

REF_STATUS ref_twod_opposite_node( REF_CELL pri,
				   REF_INT node, REF_INT *opposite)
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell;

  *opposite = REF_EMPTY;
  
  cell = ref_cell_first_with( pri, node );

  RSS( ref_cell_nodes( pri, cell, nodes ), "get first prism about node" );

  if ( node == nodes[0] ) *opposite = nodes[3];
  if ( node == nodes[3] ) *opposite = nodes[0];

  if ( node == nodes[1] ) *opposite = nodes[4];
  if ( node == nodes[4] ) *opposite = nodes[1];

  if ( node == nodes[2] ) *opposite = nodes[5];
  if ( node == nodes[5] ) *opposite = nodes[2];

  return ((REF_EMPTY==(*opposite))?REF_NOT_FOUND:REF_SUCCESS);
}

REF_STATUS ref_twod_opposite_edge( REF_CELL pri, 
				   REF_INT node0, REF_INT node1,  
				   REF_INT *node2, REF_INT *node3 )
{
  REF_INT ncell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell_to_split[2];

  *node2 = REF_EMPTY;
  *node3 = REF_EMPTY;

  RSS( ref_cell_list_with( pri, node0, node1,
			   2, &ncell, cell_to_split), "more than two" );

  RAS( ncell > 0, "no cells found");

  RSS( ref_cell_nodes(pri, cell_to_split[0], nodes), "nodes" );

  if ( node0 == nodes[0] && node1 == nodes[1] )
    { *node2=nodes[3]; *node3=nodes[4]; return REF_SUCCESS; }
  if ( node0 == nodes[1] && node1 == nodes[0] )
    { *node2=nodes[4]; *node3=nodes[3]; return REF_SUCCESS; }

  if ( node0 == nodes[1] && node1 == nodes[2] )
    { *node2=nodes[4]; *node3=nodes[5]; return REF_SUCCESS; }
  if ( node0 == nodes[2] && node1 == nodes[1] )
    { *node2=nodes[5]; *node3=nodes[4]; return REF_SUCCESS; }

  if ( node0 == nodes[2] && node1 == nodes[0] )
    { *node2=nodes[5]; *node3=nodes[3]; return REF_SUCCESS; }
  if ( node0 == nodes[0] && node1 == nodes[2] )
    { *node2=nodes[3]; *node3=nodes[5]; return REF_SUCCESS; }

  if ( node0 == nodes[3] && node1 == nodes[4] )
    { *node2=nodes[0]; *node3=nodes[1]; return REF_SUCCESS; }
  if ( node0 == nodes[4] && node1 == nodes[3] )
    { *node2=nodes[1]; *node3=nodes[0]; return REF_SUCCESS; }

  if ( node0 == nodes[4] && node1 == nodes[5] )
    { *node2=nodes[1]; *node3=nodes[2]; return REF_SUCCESS; }
  if ( node0 == nodes[5] && node1 == nodes[4] )
    { *node2=nodes[2]; *node3=nodes[1]; return REF_SUCCESS; }

  if ( node0 == nodes[5] && node1 == nodes[3] )
    { *node2=nodes[2]; *node3=nodes[0]; return REF_SUCCESS; }
  if ( node0 == nodes[3] && node1 == nodes[5] )
    { *node2=nodes[0]; *node3=nodes[2]; return REF_SUCCESS; }

   return REF_FAILURE;
}

