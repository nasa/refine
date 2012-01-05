#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"
#include "ref_fixture.h"

#include "ref_export.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_metric.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_GRID ref_grid;
      RSS( ref_fixture_pyr_grid( &ref_grid ), "pyr" );
      RSS( ref_export_b8_ugrid( ref_grid, argv[1] ), "export b8 ugrid" );
      return 0;
    }

  return 0;
}
