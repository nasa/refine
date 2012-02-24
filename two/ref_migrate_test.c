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

static REF_STATUS tear_down( REF_MIGRATE ref_migrate )
{
  REF_GRID ref_grid;

  ref_grid = ref_migrate_grid(ref_migrate);

  TSS(ref_migrate_free(ref_migrate),"free");

  TSS( ref_grid_free(ref_grid),"free" );

  return REF_SUCCESS;
}

int main( int argc, char *argv[] )
{

  TSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 == argc )
    {
      REF_GRID import_grid;
      REF_MIGRATE ref_migrate;
      char grid_file[] = "ref_migrate_test.b8.ugrid";

      if ( ref_mpi_master ) 
	{
	  REF_GRID export_grid;
	  TSS(ref_fixture_pri_stack_grid( &export_grid ), "set up tet" );
	  TSS(ref_export_b8_ugrid( export_grid, grid_file ), "export" );
	  TSS(ref_grid_free(export_grid),"free" );
	}

      TSS(ref_part_b8_ugrid( &import_grid, grid_file ), "import" );
      TSS(ref_migrate_create(&ref_migrate,import_grid),"create");

      TSS( tear_down( ref_migrate ), "tear down");
      if ( ref_mpi_master ) TEIS(0, remove( grid_file ), "test clean up");
    }

  if ( 1 < argc )
    {
      REF_GRID import_grid;
      REF_MIGRATE ref_migrate;

      TSS(ref_part_b8_ugrid( &import_grid, argv[1] ), "import" );
      TSS(ref_migrate_create(&ref_migrate,import_grid),"create");

      TSS( tear_down( ref_migrate ), "tear down");
    }

  TSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
