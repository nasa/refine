#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_clump.h"

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_sort.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include "ref_fixture.h"
#include "ref_export.h"
#include  "ref_dict.h"
#include  "ref_edge.h"
#include "ref_split.h"
#include  "ref_adapt.h"
#include   "ref_collapse.h"
#include    "ref_math.h"
#include   "ref_smooth.h"
#include    "ref_twod.h"
#include   "ref_gather.h"
#include   "ref_metric.h"

int main( int argc, char *argv[] )
{

  if ( 2 == argc )
    { /* ball */
      REF_GRID ref_grid;
      REF_INT node;

      RSS( ref_fixture_twod_brick_grid( &ref_grid ), "2d brick" );
      node = 10;
      RSS(ref_clump_tri_around(ref_grid, node, argv[1]),"clump");
      RSS( ref_grid_free(ref_grid),"free");

      RSS( ref_fixture_tet_brick_grid( &ref_grid ), "tet brick" );
      RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "unit metric" );
      node = 10;
      RSS(ref_clump_around(ref_grid, node, argv[1]),"clump");
      RSS( ref_grid_free(ref_grid),"free");

    }

  return 0;
}
