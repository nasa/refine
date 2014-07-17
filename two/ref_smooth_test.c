
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
#include "ref_math.h"

#include "ref_import.h"
#include "ref_export.h"
#include "ref_part.h"
#include "ref_histogram.h"

#include "ref_edge.h"
#include "ref_dict.h"
#include "ref_migrate.h"
#include "ref_gather.h"
#include "ref_adapt.h"
#include "ref_collapse.h"
#include "ref_split.h"

#include "ref_metric.h"

#include "ref_mpi.h"

#include "ref_malloc.h"

static REF_STATUS ref_smooth_tri_twod( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL ideal[3];

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

  RSS( ref_metric_unit_node( ref_node ), "unit node" )

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  RSS(ref_smooth_ideal_tri(ref_grid,nodes[0], cell, ideal),"ideal");

  RWDS(0.5,ideal[0],-1,"ideal x");
  RWDS(0.0,ideal[1],-1,"ideal y");
  RWDS(0.5*sqrt(3.0),ideal[2],-1,"ideal z");
  return REF_SUCCESS;
}

/*
ev ~/bibtex-refs/Ref/alauzet-ec-2013-topology-moving-mesh.pdf
*/


int main( int argc, char *argv[] )
{

  if ( argc > 2 )
    {
      REF_GRID ref_grid;

      RSS( ref_mpi_start( argc, argv ), "start" );

      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );

      RSS( ref_part_metric( ref_grid_node(ref_grid), argv[2] ), "get metric");

      RSS( ref_histogram_quality( ref_grid ), "qual");

      RSS( ref_export_tec_surf( ref_grid, "ref_smooth_test_0.tec" ), 
	   "surf");

      RSS( ref_grid_free( ref_grid ), "free");

      RSS( ref_mpi_stop( ), "stop" );
    }

  { /*  */
    REF_GRID ref_grid;
    RSS( ref_smooth_tri_twod( &ref_grid ), "2d fixture" );

    RSS( ref_smooth_twod( ref_grid, 0 ), "smooth" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
