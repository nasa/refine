#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_test.h"


int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 == argc )
    {
      REF_GRID import_grid;
      char grid_file[] = "ref_migrate_test.b8.ugrid";

      if ( ref_mpi_master ) 
	{
	  REF_GRID export_grid;
	  RSS(ref_fixture_pri_stack_grid( &export_grid ), "set up tet" );
	  RSS(ref_export_b8_ugrid( export_grid, grid_file ), "export" );
	  RSS(ref_grid_free(export_grid),"free" );
	}

      RSS(ref_part_b8_ugrid( &import_grid, grid_file ), "import" );
      RSS(ref_migrate_new_part(import_grid),"create");

      RSS( ref_migrate_shufflin( import_grid ), "shufflin");

      RSS( ref_grid_free( import_grid ), "free");
      if ( ref_mpi_master ) REIS(0, remove( grid_file ), "test clean up");
    }

  if ( 1 < argc )
    {
      REF_GRID import_grid;

      if ( ref_mpi_master ) 
	printf("%d procs, read %s\n",ref_mpi_n,argv[1]);

      ref_mpi_stopwatch_start();
      RSS(ref_part_b8_ugrid( &import_grid, argv[1] ), "import" );
      ref_mpi_stopwatch_stop("read");

      RSS(ref_migrate_new_part(import_grid),"new part");
      ref_mpi_stopwatch_stop("new part");

      RSS( ref_migrate_shufflin( import_grid ), "shufflin");
      ref_mpi_stopwatch_stop("shufflin");

      RSS( ref_migrate_part_viz( import_grid ), "part_viz");

      RSS( ref_grid_free( import_grid ), "free");
    }

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
