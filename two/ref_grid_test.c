#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_CELL ref_cell;
  REF_INT node_per;
  REF_INT group;

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

  return 0;
}
