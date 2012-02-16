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
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"


int main( int argc, char *argv[] )
{

  TSS( ref_mpi_start( argc, argv ), "start" );

  { /* one even */
    TEIS( 0, ref_part_first( 4, 4, 0 ), "first");
    TEIS( 1, ref_part_first( 4, 4, 1 ), "first");
    TEIS( 2, ref_part_first( 4, 4, 2 ), "first");
    TEIS( 3, ref_part_first( 4, 4, 3 ), "first");
    TEIS( 4, ref_part_first( 4, 4, 4 ), "first");
  }

  { /* one run out */
    TEIS( 0, ref_part_first( 2, 4, 0 ), "first");
    TEIS( 1, ref_part_first( 2, 4, 1 ), "first");
    TEIS( 2, ref_part_first( 2, 4, 2 ), "first");
    TEIS( 2, ref_part_first( 2, 4, 3 ), "first");
    TEIS( 2, ref_part_first( 2, 4, 4 ), "first");
  }

  { /* two run out */
    TEIS( 0, ref_part_first( 7, 4, 0 ), "first");
    TEIS( 2, ref_part_first( 7, 4, 1 ), "first");
    TEIS( 4, ref_part_first( 7, 4, 2 ), "first");
    TEIS( 6, ref_part_first( 7, 4, 3 ), "first");
    TEIS( 7, ref_part_first( 7, 4, 4 ), "first");
  }

  { /* two even */
    TEIS( 0, ref_part_first( 4, 2, 0 ), "first");
    TEIS( 2, ref_part_first( 4, 2, 1 ), "first");
    TEIS( 4, ref_part_first( 4, 2, 2 ), "first");
  }

  { /* single */
    TEIS( 0, ref_part_implicit( 4, 1, 0 ), "part");
    TEIS( 0, ref_part_implicit( 4, 1, 1 ), "part");
    TEIS( 0, ref_part_implicit( 4, 1, 2 ), "part");
    TEIS( 0, ref_part_implicit( 4, 1, 3 ), "part");
  }

  { /* one even */
    TEIS( 0, ref_part_implicit( 4, 4, 0 ), "part");
    TEIS( 1, ref_part_implicit( 4, 4, 1 ), "part");
    TEIS( 2, ref_part_implicit( 4, 4, 2 ), "part");
    TEIS( 3, ref_part_implicit( 4, 4, 3 ), "part");
  }

  { /* one run out */
    TEIS( 0, ref_part_implicit( 2, 4, 0 ), "part");
    TEIS( 1, ref_part_implicit( 2, 4, 1 ), "part");
  }

  { /* two run out */
    TEIS( 0, ref_part_implicit( 7, 4, 0 ), "part");
    TEIS( 0, ref_part_implicit( 7, 4, 1 ), "part");
    TEIS( 1, ref_part_implicit( 7, 4, 2 ), "part");
    TEIS( 1, ref_part_implicit( 7, 4, 3 ), "part");
    TEIS( 2, ref_part_implicit( 7, 4, 4 ), "part");
    TEIS( 2, ref_part_implicit( 7, 4, 5 ), "part");
    TEIS( 3, ref_part_implicit( 7, 4, 6 ), "part");
  }

  { /* two even */
    TEIS( 0, ref_part_implicit( 4, 2, 0 ), "part");
    TEIS( 0, ref_part_implicit( 4, 2, 1 ), "part");
    TEIS( 1, ref_part_implicit( 4, 2, 2 ), "part");
    TEIS( 1, ref_part_implicit( 4, 2, 3 ), "part");
  }

  { /* part */
    REF_GRID export_grid, import_grid;
    char grid_file[] = "ref_part_test.b8.ugrid";
    char viz_file[256];
    
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    if ( ref_mpi_master ) 
      {
	TSS(ref_export_b8_ugrid( export_grid, grid_file ), "export" );
      }
    TSS(ref_part_b8_ugrid( &import_grid, grid_file ), "import" );

    if ( ref_mpi_n > 1 ) 
      {
	sprintf(viz_file, "ref_part_p%d.tec", ref_mpi_id);
	RSS( ref_export_by_extension( import_grid, viz_file ), "export");
	if ( ref_mpi_master ) ref_grid_inspect( export_grid );
	ref_grid_inspect( import_grid );  
      }

    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    if ( ref_mpi_master ) TEIS(0, remove( grid_file ), "test clean up");
  }

  TSS( ref_mpi_stop( ), "stop" );

  return 0;
}
