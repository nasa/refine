#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_test.h"

#include "ref_part.h"
#include "ref_mpi.h"
#include "ref_fixture.h"
#include "ref_export.h"

#include "ref_dict.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"


int main( int argc, char *argv[] )
{

  TSS( ref_mpi_start( argc, argv ), "start" );

  { /* part */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";

    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    if ( ref_mpi_master ) 
      TSS(ref_export_b8_ugrid( export_grid, file ), "export" );
    
    TSS(ref_part_b8_ugrid( &import_grid, file ), "import" );

    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    if ( ref_mpi_master ) TEIS(0, remove( file ), "test clean up");
  }

  TSS( ref_mpi_stop( ), "stop" );

  return 0;
}
