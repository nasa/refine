
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_smooth.h"

REF_STATUS ref_smooth_twod( REF_GRID ref_grid, REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL f, d[3];
  REF_DBL dcost, dcost_dl, dl;
  REF_BOOL verbose = REF_FALSE;

  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      RSS( ref_node_tri_dquality_dnode0(ref_node, nodes, 
					&f, d), "qual deriv" );
      if (verbose)
	printf("cost %10.8f : %12.8f %12.8f %12.8f\n",f,d[0],d[1],d[2]);
    }

  dcost = 1.0-f;
  dcost_dl = sqrt(d[0]*d[0]+d[1]*d[1]*d[2]*d[2]);
  dl = dcost/dcost_dl;

  ref_node_xyz(ref_node,0,node) += dl*d[0];
  ref_node_xyz(ref_node,1,node) += dl*d[1];
  ref_node_xyz(ref_node,2,node) += dl*d[2];
 
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      RSS( ref_node_tri_dquality_dnode0(ref_node, nodes, 
					&f, d), "qual deriv" );
    }

  if (verbose)
    printf("rate %12.8f dcost %12.8f dl %12.8f\n",
	   (1.0-f)/dcost,dcost,dl);
  

  return REF_SUCCESS;
}
