#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_test.h"
#include "ref_fixture.h"

int main( void )
{

  { /* export .vtk tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.vtk";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_vtk( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export .tec tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_tec( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export .fgrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.fgrid";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_fgrid( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export .ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.ugrid";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_ugrid( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  { /* export .b8.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.b8.ugrid";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_b8_ugrid( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    /* TEIS(0, remove( file ), "test clean up"); */
  }

  return 0;
}
