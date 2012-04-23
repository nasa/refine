#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"

#include "ref_list.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_mpi.h"

int main( void )
{

  {  /* init */
    REF_GRID ref_grid;
    REIS(REF_NULL,ref_grid_free(NULL),"dont free NULL");

    RSS(ref_grid_create(&ref_grid),"create");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* each element */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;
    REF_INT group;

    RSS(ref_grid_create(&ref_grid),"create");

    node_per = 3;
    each_ref_grid_ref_cell( ref_grid, group, ref_cell )
      {
	node_per += 1;
	if ( 7 == node_per ) node_per = 8;
	RES( node_per, ref_cell_node_per( ref_cell), "cells in order" );
      }

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* cell with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid),"create");

    node_per = 4;
    RSS( ref_grid_cell_with( ref_grid, node_per, &ref_cell ), "with" );
    REIS( node_per, ref_cell_node_per( ref_cell ), "match" );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* face with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid),"create");

    node_per = 4;
    RSS( ref_grid_face_with( ref_grid, node_per, &ref_cell ), "with" );
    REIS( node_per, ref_cell_node_per( ref_cell ), "match" );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  return 0;
}
