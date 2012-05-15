#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_plot3d.h"
#include "ref_grid.h"
#include "ref_import.h"
#include "ref_dict.h"
#include "ref_node.h"
#include "ref_cell.h"
#include "ref_adj.h"
#include "ref_mpi.h"
#include "ref_matrix.h"
#include "ref_sort.h"
#include "ref_list.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_PLOT3D ref_plot3d;
      RSS( ref_plot3d_from_file( &ref_plot3d, argv[1] ), "from file" );
      RSS( ref_plot3d_tec( ref_plot3d, "ref_plot3d_test.tec" ), "tec" );
      RSS( ref_plot3d_free( ref_plot3d ), "free" );
      return 0;
    }

  if (3 == argc) 
    {
      REF_PLOT3D ref_plot3d;
      REF_GRID ref_grid;
      RSS( ref_plot3d_from_file( &ref_plot3d, argv[1] ), "from file" );
      RSS( ref_import_by_extension( &ref_grid, argv[2] ), "by ext" );

      RSS( ref_plot3d_free( ref_plot3d ), "free" );
      RSS( ref_grid_free( ref_grid ), "free" );
      return 0;
    }

  return 0;
}
