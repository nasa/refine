#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_grid.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_validation.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"

#include "ref_cell.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {
      printf("validating\n");

      printf("reading %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from ugrid");
      printf("complete.\n");
      
      RSS(ref_grid_inspect( ref_grid ), "inspection");

      printf("validate.\n");
      RSS( ref_validation_all( ref_grid ), "invalid grid" );

      printf("vtk.\n");
      RSS( ref_export_vtk( ref_grid, "validate.vtk" ), "vtk" );

      printf("tec.\n");
      RSS( ref_export_tec( ref_grid, "validate.tec" ), "tec" );

      printf("done.\n");
      return 0;
    }

  return 0;
}
