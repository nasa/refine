#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid_import.h"
#include "ref_grid_export.h"
#include "ref_test.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc<2) 
    {
      printf("usage: %s filename.ugrid\n",argv[0]);
      return 0;
    }

  TSS(ref_grid_import_ugrid( argv[1], &ref_grid ),"from ugrid");

  TSS(ref_grid_export_vtk(ref_grid, "ref2.vtk"),"to vtk");

  TSS(ref_grid_free(ref_grid),"free");

  return 0;
}
