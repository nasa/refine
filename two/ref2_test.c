#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"
#include "ref_test.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

#include "ref_edge.h"

#include "ref_face.h"
#include "ref_sort.h"

#include "ref_quality.h"
#include "ref_hexdiv.h"

#include "ref_validation.h"

#include "ref_math.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc<2) 
    {
      printf("usage: %s filename.ugrid\n",argv[0]);
      return 0;
    }

  printf("reading %s\n",argv[1]);
  TSS(ref_import_ugrid( argv[1], &ref_grid ),"from ugrid");
  printf("complete.\n");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("hex quality.\n");
  TSS(ref_quality_hex( ref_grid ),"quality");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("vtk.\n");
  TSS(ref_export_vtk(ref_grid, "ref2.vtk"),"to vtk");

  printf("tec.\n");
  TSS(ref_export_tec(ref_grid, "ref2.tec"),"to tec");

  printf("validate.\n");
  TSS( ref_validation_all( ref_grid ), "invalid grid" );

  printf("ugrid.\n");
  TSS(ref_export_ugrid(ref_grid, "ref2.ugrid"),"to ugrid");

  printf("free.\n");
  TSS(ref_grid_free(ref_grid),"free");

  printf("done.\n");
  return 0;
}
