#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_edge.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_sort.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  {  /* make edges shared by two elements */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell;

    RSS(ref_grid_create(&ref_grid),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
    nodes[3] = 3; nodes[4] = 4; nodes[5] = 5;
    RSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

    nodes[0] = 3; nodes[1] = 4; nodes[2] = 5; nodes[3] = 6;
    RSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add tet");

    RSS(ref_edge_create(&ref_edge,ref_grid),"create");

    REIS( 12, ref_edge_n(ref_edge), "check total edges");

    RSS(ref_edge_free(ref_edge),"edge");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* find edge with nodes */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell, edge;

    RSS(ref_grid_create(&ref_grid),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
    nodes[3] = 3; nodes[4] = 4; nodes[5] = 5;
    RSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

    RSS(ref_edge_create(&ref_edge,ref_grid),"create");

    RSS( ref_edge_with( ref_edge, 0, 1, &edge ), "find" );
    REIS( 0, edge, "right one" );
    REIS( REF_NOT_FOUND,ref_edge_with( ref_edge, 0, 4, &edge ), "not found" );
    REIS( REF_EMPTY, edge, "not found" );

    RSS(ref_edge_free(ref_edge),"edge");
    RSS(ref_grid_free(ref_grid),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
