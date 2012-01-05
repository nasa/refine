#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_fixture.h"

#include "ref_adj.h"
#include "ref_sort.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      RSS( ref_import_examine_header( argv[1] ), "examine header" );
      return 0;
    }


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

  { /* export import .b8.ugrid tet */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_b8_ugrid( export_grid, file ), "export" );
    TSS(ref_import_b8_ugrid( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TEIS( ref_cell_n(ref_grid_tri(export_grid)),
	  ref_cell_n(ref_grid_tri(import_grid)), "tri count" );
    TEIS( ref_cell_n(ref_grid_qua(export_grid)),
	  ref_cell_n(ref_grid_qua(import_grid)), "qua count" );
    TEIS( ref_cell_n(ref_grid_tet(export_grid)),
	  ref_cell_n(ref_grid_tet(import_grid)), "tet count" );
    TES( ref_node_xyz( ref_grid_node(export_grid),0,1),
	 ref_node_xyz( ref_grid_node(import_grid),0,1), "x 1" );
    TEIS( ref_cell_c2n(ref_grid_tet(export_grid),0,0),
	  ref_cell_c2n(ref_grid_tet(import_grid),0,0), "tet node0" );
    TEIS( ref_cell_c2n(ref_grid_tet(export_grid),1,0),
	  ref_cell_c2n(ref_grid_tet(import_grid),1,0), "tet node 1" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export import .b8.ugrid tet by extension */
    REF_GRID export_grid, import_grid;
    char file[] = "ref_import_test.b8.ugrid";
    TSS(ref_fixture_tet_grid( &export_grid ), "set up tet" );
    TSS(ref_export_b8_ugrid( export_grid, file ), "export" );
    TSS(ref_import_by_extension( &import_grid, file ), "import" );
    TEIS( ref_node_n(ref_grid_node(export_grid)),
	  ref_node_n(ref_grid_node(import_grid)), "node count" );
    TSS(ref_grid_free(import_grid),"free");
    TSS(ref_grid_free(export_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  return 0;
}
