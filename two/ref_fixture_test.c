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

  {
    REF_GRID ref_grid;

    RSS( ref_fixture_pri_grid( &ref_grid ), "fix" );

    RSS( ref_validation_cell_node( ref_grid ), "invalid pri" );

    RSS( ref_grid_free(ref_grid),"free");
  }

  {
    REF_GRID ref_grid;

    RSS( ref_fixture_pri_stack_grid( &ref_grid ), "fix" );

    if ( 1 < ref_mpi_n ) 
      RSS( ref_export_tec_part( ref_grid, "ref_fixture_orig_stack" ), "see" );

    RSS( ref_validation_cell_node( ref_grid ), "invalid stack" );

    RSS( ref_migrate_to_balance( ref_grid ), "bal" );

    if ( 1 < ref_mpi_n ) 
      RSS( ref_export_tec_part( ref_grid, "ref_fixture_bal_stack" ), "see" );

    RSS( ref_grid_free(ref_grid),"free");
  }

  {
    REF_GRID ref_grid;

    RSS( ref_fixture_pri_grid( &ref_grid ), "fix" );

    if ( 1 < ref_mpi_n ) 
      RSS( ref_export_tec_part( ref_grid, "ref_fixture_orig_pri" ), "see" );

    RSS( ref_validation_cell_node( ref_grid ), "invalid pri" );

    RSS( ref_migrate_to_balance( ref_grid ), "bal" );

    if ( 1 < ref_mpi_n ) 
      RSS( ref_export_tec_part( ref_grid, "ref_fixture_bal_pri" ), "see" );

    RSS( ref_grid_free(ref_grid),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
