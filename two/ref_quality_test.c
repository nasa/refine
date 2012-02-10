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
#include "ref_shard.h"
#include "ref_subdiv.h"
#include "ref_edge.h"
#include "ref_math.h"
#include "ref_dict.h"

#include "ref_swap.h"

#include "ref_test.h"
#include "ref_fixture.h"

int main( void )
{

  { /* find mark */
    REF_GRID ref_grid;
    REF_DBL vol;

    TSS( ref_fixture_tet_grid( &ref_grid ), "tet fixture" );
    TSS( ref_quality_tet_vol( ref_grid, 0, &vol), "get vol");

    TWDS(1.0/6.0,vol,-1.0,"expected vol");

    TSS( ref_grid_free( ref_grid ), "free" );
  }

  return 0;
}
