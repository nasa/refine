#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_matrix.h"

#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_gather.h"
#include "ref_adapt.h"

#include "ref_collapse.h"
#include "ref_split.h"
#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 < argc )
    {
      REF_GRID ref_grid;

      ref_mpi_stopwatch_start();
      RSS(ref_part_b8_ugrid( &ref_grid, argv[1] ), "import" );
      ref_mpi_stopwatch_stop("read");
      RSS(ref_migrate_new_part(ref_grid),"new part");
      ref_mpi_stopwatch_stop("new part");
      RSS( ref_migrate_shufflin( ref_grid ), "shufflin");
      ref_mpi_stopwatch_stop("shufflin");

      RSS( ref_adapt_pass( ref_grid ), "pass");

      ref_mpi_stopwatch_start();
      RSS( ref_gather_b8_ugrid( ref_grid, "ref_adapt_test.b8.ugrid" ), 
	   "gather");
      ref_mpi_stopwatch_stop("gather");

      RSS( ref_grid_free( ref_grid ), "free");
    }

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
