#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_test.h"

#include "ref_adj.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_CELL ref_cell;
  REF_INT node_per;
  REF_INT group;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT cell;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  /* init */

  TFS(ref_grid_free(NULL),"dont free NULL");

  TSS(ref_grid_create(&ref_grid),"create");
  TSS(ref_grid_free(ref_grid),"free");

  /* each element */

  TSS(ref_grid_create(&ref_grid),"create");

  node_per = 3;
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      node_per += 1;
      if ( 7 == node_per ) node_per = 8;
      TES( node_per, ref_cell_node_per( ref_cell), "cells in order" );
    }

  TSS(ref_grid_free(ref_grid),"free");

  /* make edges */

  TSS(ref_grid_create(&ref_grid),"create");

  TES( 0, ref_grid_nedge(ref_grid), "check initial edges");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
  nodes[3] = 3; nodes[4] = 4; nodes[5] = 5;
  TSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

  nodes[0] = 3; nodes[1] = 4; nodes[2] = 5; nodes[3] = 6;
  TSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add tet");

  TSS( ref_grid_make_edges(ref_grid), "make edges");

  TEIS( 12, ref_grid_nedge(ref_grid), "check total edges");

  TSS(ref_grid_free(ref_grid),"free");
  

  return 0;
}
