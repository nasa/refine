#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_validation.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

#include "ref_face.h"
#include "ref_sort.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {
      printf("validating\n");

      printf("reading %s\n",argv[1]);
      TSS(ref_import_ugrid( &ref_grid, argv[1] ),"from ugrid");
      printf("complete.\n");
      
      printf("vtk.\n");
      TSS( ref_export_vtk( ref_grid, "validate.vtk" ), "vtk" );
      printf("tec.\n");
      TSS( ref_export_tec( ref_grid, "validate.tec" ), "tec" );

      printf("validate.\n");
      TSS( ref_validation_all( ref_grid ), "invalid grid" );

      printf("done.\n");
      return 0;
    }

  return 0;
}
