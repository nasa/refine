#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"

#include "ref_list.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_fixture.h"
#include "ref_malloc.h"

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

  {  /* imply metric tet */
    REF_DBL tol = -1.0;
    REF_GRID ref_grid;
    REF_DBL *metric;
    REF_INT node;

    RSS( ref_fixture_tet_grid( &ref_grid ), "tet" );

    ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

    RSS( ref_grid_imply_metric( ref_grid, metric ), "imply" );

    each_ref_node_valid_node( ref_grid_node(ref_grid), node )
      {
	RWDS( 1.0, metric[0+6*node], tol, "m[0]");
	RWDS( 0.5, metric[1+6*node], tol, "m[1]");
	RWDS( 0.5, metric[2+6*node], tol, "m[2]");
	RWDS( 1.0, metric[3+6*node], tol, "m[3]");
	RWDS( 0.5, metric[4+6*node], tol, "m[4]");
	RWDS( 1.0, metric[5+6*node], tol, "m[5]");
      }

    ref_free( metric );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  return 0;
}
