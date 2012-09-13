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
#include "ref_gather.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_mpi.h"
#include "ref_part.h"
#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 == ref_mpi_n )
    {
      REF_GRID ref_grid;
      REF_MIGRATE ref_migrate;
      REF_INT keep, lose;

      RSS(ref_fixture_pri_grid(&ref_grid),"set up grid");
      RSS(ref_migrate_create(&ref_migrate,ref_grid),"set up mig");

      /* fake 2 proc */
      ref_node_part(ref_grid_node(ref_grid),3) = 1;
      ref_node_part(ref_grid_node(ref_grid),4) = 1;
      ref_node_part(ref_grid_node(ref_grid),5) = 1;

      keep = 0; lose = 3; 
      RSS( ref_migrate_2d_agglomeration_keep( ref_migrate, keep, lose), "0-3");

      REIS( 0,         ref_migrate_global( ref_migrate, 0 ), "mark" );
      REIS( REF_EMPTY, ref_migrate_global( ref_migrate, 3 ), "mark" );

      RSS( ref_migrate_free( ref_migrate ), "free migrate");
      RSS( ref_grid_free( ref_grid ), "free gride");
    }

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

      RSS( ref_gather_tec_part( import_grid, "ref_migrate_test.tec" ), 
	   "part_viz");

      RSS( ref_grid_free( import_grid ), "free");
    }

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
