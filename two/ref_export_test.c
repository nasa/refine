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

  { /* export .ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.ugrid";
    TSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    TSS(ref_export_ugrid( ref_grid, file ),"export" );
    TSS(ref_grid_free(ref_grid),"free");
    TEIS(0, remove( file ), "test clean up");
  }

  return 0;
}
