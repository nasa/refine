
#include <stdlib.h>
#include <stdio.h>

#include "ref_axi.h"

#include "ref_node.h"

REF_STATUS ref_axi_wedge( REF_GRID ref_grid )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;

  ref_node = ref_grid_node(ref_grid);

  ref_node_n(ref_node)=0;
  ref_node_n_global(ref_node)=0;

  ref_cell = ref_grid_qua(ref_grid);
  ref_cell_n(ref_cell) = 0;

  return REF_SUCCESS;
}

