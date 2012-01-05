#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"

#include "ref_sort.h"

int main( void )
{
  REF_GRID ref_grid;

  {  /* init */
    TFS(ref_grid_free(NULL),"dont free NULL");

    TSS(ref_grid_create(&ref_grid),"create");
    TSS(ref_grid_free(ref_grid),"free");
  }

  {  /* each element */

    REF_CELL ref_cell;
    REF_INT node_per;
    REF_INT group;

    TSS(ref_grid_create(&ref_grid),"create");

    node_per = 3;
    each_ref_grid_ref_cell( ref_grid, group, ref_cell )
      {
	node_per += 1;
	if ( 7 == node_per ) node_per = 8;
	TES( node_per, ref_cell_node_per( ref_cell), "cells in order" );
      }

    TSS(ref_grid_free(ref_grid),"free"); 
  }

  return 0;
}
