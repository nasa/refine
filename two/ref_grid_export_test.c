#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid_export.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TSS(ref_grid_create( &ref_grid ), "create" );
  TSS(ref_grid_export_vtk( ref_grid, "ref_grid_export_test.vtk" ),"export" );
  TSS(ref_grid_free(ref_grid),"free");

  return 0;
}
