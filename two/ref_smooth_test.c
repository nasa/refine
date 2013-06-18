
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_smooth.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_list.h"
#include "ref_node.h"
#include "ref_matrix.h"

#include "ref_mpi.h"

#include "ref_malloc.h"

static REF_STATUS ref_smooth_tri_twod( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS( ref_grid_create( ref_grid_ptr ), "grid" );
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /*
   0    z
   |\   |
   1-2  +-x
   */

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[0] = node;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[1] = node;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[2] = node;

  nodes[3] = 10;

  for ( node=0;node<4;node++)
    {
      ref_node_metric(ref_node,0,node) = 1.0;
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,2,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0;
      ref_node_metric(ref_node,4,node) = 0.0;
      ref_node_metric(ref_node,5,node) = 1.0;
    }

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}

int main( void )
{

  { /*  */
    REF_GRID ref_grid;
    RSS( ref_smooth_tri_twod( &ref_grid ), "2d fixture" );

    RSS( ref_smooth_twod( ref_grid, 0 ), "smooth" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
