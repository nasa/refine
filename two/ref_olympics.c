#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_inflate.h"

#include "ref_math.h"

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

#include "ref_fixture.h"
#include "ref_metric.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if ( 1 == argc )
    {
      RSS( ref_fixture_tet_brick_grid( &ref_grid ), "brick");
      RSS( ref_export_by_extension( ref_grid, "ref_olympics.meshb" ), "out" );

      return 0;
    }

  if ( 3 == argc )
    {
      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "in");
      RSS( ref_metric_olympic_node( ref_grid_node(ref_grid), atof(argv[2]) ),
	   "oly" );

      return 0;
    }

  return 0;
}
