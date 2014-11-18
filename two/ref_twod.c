
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
