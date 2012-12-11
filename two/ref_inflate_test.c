#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_inflate.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  if ( 3 == argc )
    {
      REF_GRID ref_grid;
      REF_INT faceid;

      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "read grid" );

      faceid = atoi( argv[2] );

      printf("inflating face %d\n",faceid);

      RSS( ref_inflate_face( ref_grid, faceid ), "inflate" );

      RSS( ref_export_by_extension( ref_grid, "ref_inflate_test.tec" ), "tec" );
      RSS( ref_export_by_extension( ref_grid, "ref_inflate_test.b8.ugrid" ), "b8" );

      RSS(ref_grid_free(ref_grid),"free");
    }

  return 0;
}
