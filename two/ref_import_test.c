#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"
#include "ref_test.h"
#include "ref_fixture.h"

#include "ref_adj.h"

int main( void )
{

  { /* export import .fgrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.fgrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_fgrid( export_grid, file ), "export" );
    TSS(ref_import_fgrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export import .ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_ugrid( export_grid, file ), "export" );
    TSS(ref_import_ugrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  return 0;
}
