#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_edge.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_sort.h"
#include "ref_node.h"
#include "ref_adj.h"
#include "ref_metric.h"

int main( void )
{

  {  /* make edges shared by two elements */
    REF_EDGE ref_edge;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell;

    TSS(ref_grid_create(&ref_grid),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
    nodes[3] = 3; nodes[4] = 4; nodes[5] = 5;
    TSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

    nodes[0] = 3; nodes[1] = 4; nodes[2] = 5; nodes[3] = 6;
    TSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add tet");

    TSS(ref_edge_create(&ref_edge,ref_grid),"create");

    TEIS( 12, ref_edge_n(ref_edge), "check total edges");

    TSS(ref_edge_free(ref_edge),"edge");
    TSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
