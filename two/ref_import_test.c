#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"
#include "ref_test.h"

#include "ref_adj.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  /* AFLR3 ugrid*/

  /*
  TSS(ref_import_ugrid("../test/gbumpf_MX.ugrid",&ref_grid),"from ugrid");
  TSS(ref_export_vtk(ref_grid,"gbumpf_MX.vtk"),"vtk");
  TSS(ref_grid_free(ref_grid),"free");
  */

  /* FAST fgrid */

  TSS(ref_import_fgrid("../test/gbumpn.fgrid",&ref_grid),"from fgrid");
  TSS(ref_export_vtk(ref_grid,"gbumpn.vtk"),"vtk");
  TSS(ref_export_fgrid(ref_grid,"gbumpn.fgrid"),"fgrid");
  TSS(ref_grid_free(ref_grid),"free");

  TSS(ref_import_fgrid("../test/om6_inv08.fgrid",&ref_grid),"from fgrid");
  TSS(ref_export_vtk(ref_grid,"om6_inv08.vtk"),"vtk");
  TSS(ref_export_fgrid(ref_grid,"om6_inv08.fgrid"),"fgrid");
  TSS(ref_grid_free(ref_grid),"free");

  return 0;
}
