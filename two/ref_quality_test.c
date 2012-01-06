#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_quality.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_hexdiv.h"
#include "ref_subdiv.h"
#include "ref_edge.h"
#include "ref_math.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {

      printf("reading %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from ugrid");
      printf("complete.\n");
      
      RSS(ref_grid_inspect( ref_grid ), "inspection");

      printf("validate.\n");
      RSS( ref_quality_multiple_face_cell( ref_grid ), "invalid grid" );

      printf("done.\n");
      return 0;
    }

  return 0;
}
