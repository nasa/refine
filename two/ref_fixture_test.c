#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "ref_fixture.h"

#include "ref_export.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_metric.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_GRID ref_grid;
      RSS( ref_fixture_pri_grid( &ref_grid ), "fix" );
      RSS( ref_export_tec( ref_grid, argv[1] ), "export" );
      return 0;
    }

  return 0;
}
