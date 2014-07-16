
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

REF_STATUS ref_smooth_tri_quality_around( REF_GRID ref_grid, 
					  REF_INT node,
					  REF_DBL *min_quality )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL quality;

  *min_quality = 1.0;
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
      none_found = REF_FALSE;
      RSS( ref_node_tri_quality( ref_node, 
				 nodes,  
				 &quality ), "qual" );
      *min_quality = MIN( *min_quality, quality );
    }
      
  if ( none_found )
    {
      *min_quality = -2.0;
      return REF_NOT_FOUND;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_ideal_tri( REF_GRID ref_grid, 
				 REF_INT node,
				 REF_INT tri,
				 REF_DBL *ideal_location )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT n0, n1;
  REF_INT ixyz;
  RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri, nodes ), "get tri");
  n0=REF_EMPTY;n1=REF_EMPTY;
  if ( node == nodes[0])
    n0=nodes[1];n1=nodes[2];
  if ( node == nodes[1])
    n0=nodes[2];n1=nodes[0];
  if ( node == nodes[2])
    n0=nodes[0];n1=nodes[1];
  if ( n0==REF_EMPTY || n1==REF_EMPTY)
    THROW("empty triangle side");
  for (ixyz=0;ixyz<3;ixyz++)
    ideal_location[ixyz] = 0.5*( ref_node_xyz(ref_node,ixyz,n0) +
				 ref_node_xyz(ref_node,ixyz,n1) );
  SKIP_BLOCK("finish ref_smooth_ideal_tri")
    return REF_IMPLEMENT;
  return REF_SUCCESS;
}
