
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_smooth.h"

REF_STATUS ref_smooth_twod( REF_GRID ref_grid, REF_INT node )
{
  REF_DBL f, d[3];
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      RSS( ref_node_tri_dquality_dnode0(ref_node, nodes, 
					&f, d), "qual deriv" );
      printf("cost %10.8f : %12.8f %12.8f\n",f,d[0],d[2]);
    }

  return REF_SUCCESS;
}
