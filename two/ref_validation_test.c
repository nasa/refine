#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_grid_import.h"
#include "ref_validation.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {
      printf("validating\n");

      printf("reading %s\n",argv[1]);
      TSS(ref_grid_import_ugrid( argv[1], &ref_grid ),"from ugrid");
      printf("complete.\n");
      
      printf("validate.\n");
      TSS( ref_validation_check( ref_grid ), "invalid grid" );

      printf("done.\n");
      return 0;
    }

  return 0;
}
