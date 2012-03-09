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
#include "ref_face.h"
#include "ref_mpi.h"

#include "ref_migrate.h"
#include "ref_validation.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 < ref_mpi_n ) 
    {
      REF_GRID ref_grid;

      RSS( ref_fixture_pri_stack_grid( &ref_grid ), "fix" );
      RSS( ref_export_tec_part( ref_grid, "orig_part" ), "see" );

      RSS( ref_validation_cell_node( ref_grid ), "valid nodes" );

      RSS( ref_migrate_to_balance( ref_grid ), "bal" );
      RSS( ref_export_tec_part( ref_grid, "bal_part" ), "see" );

    }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
