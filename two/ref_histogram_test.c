#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_histogram.h"
#include "ref_mpi.h"
#include "ref_grid.h"
#include "ref_migrate.h"
#include "ref_part.h"

#include "ref_node.h"
#include "ref_cell.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_list.h"
#include "ref_matrix.h"
#include "ref_math.h"

#include "ref_edge.h"
#include "ref_metric.h"

#include "ref_malloc.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 2 == argc )
    {
      REF_GRID ref_grid;

      RSS(ref_part_b8_ugrid( &ref_grid, argv[1] ), "part grid" );
      RSS( ref_metric_unit_node( ref_grid_node(ref_grid)), "id metric" );

      {
	REF_DBL *metric;
	ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
	RSS( ref_metric_imply_from( metric, ref_grid ), 
	     "from");
	RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid)), "to");
	RSS( ref_node_ghost_real( ref_grid_node(ref_grid) ), "ghost real");
	ref_free( metric );
      }

      RSS( ref_histogram_ratio( ref_grid ), "gram");

      RSS( ref_grid_free( ref_grid ), "free");
      RSS( ref_mpi_stop(  ), "stop" );
      return 0;
    }

  {
    REF_HISTOGRAM ref_histogram;
    REIS(REF_NULL, ref_histogram_free(NULL),"dont free NULL");
    RSS(ref_histogram_create(&ref_histogram),"create");
    RSS(ref_histogram_free(ref_histogram),"free");
  }

  RSS( ref_mpi_stop(  ), "stop" );
  return 0;
}
