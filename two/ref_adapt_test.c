#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_matrix.h"

#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_fixture.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_gather.h"
#include "ref_adapt.h"

#include "ref_collapse.h"
#include "ref_split.h"
#include "ref_edge.h"

#include "ref_subdiv.h"
#include "ref_validation.h"
#include "ref_face.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 2 < argc )
    {
      REF_GRID ref_grid;
      REF_NODE ref_node;
      REF_INT i, passes;

      ref_mpi_stopwatch_start();
      RSS(ref_part_b8_ugrid( &ref_grid, argv[1] ), "part grid" );
      ref_node = ref_grid_node(ref_grid);
      ref_mpi_stopwatch_stop("read grid");

      RSS(ref_migrate_to_balance(ref_grid),"balance");
      ref_mpi_stopwatch_stop("balance");

      RSS(ref_part_metric( ref_node, argv[2] ), "part metric" );
      ref_mpi_stopwatch_stop("read metric");

      RSS( ref_export_tec_ratio(ref_grid,"ref_adapt_orig"),"export m");

      RSS(ref_validation_cell_volume(ref_grid),"vol");

      passes = 1;
      if ( 3 < argc ) passes = atoi(argv[3]);

      {
	REF_SUBDIV ref_subdiv;

        RSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

	RSS(ref_subdiv_mark_prism_by_metric(ref_subdiv),"mark metric");
	ref_mpi_stopwatch_stop("subdiv mark");

	RSS(ref_subdiv_split(ref_subdiv),"split");
	ref_mpi_stopwatch_stop("subdiv split");

	RSS(ref_subdiv_free(ref_subdiv),"free");

	RSS(ref_validation_cell_volume(ref_grid),"vol");
	RSS(ref_migrate_to_balance(ref_grid),"balance");
	ref_mpi_stopwatch_stop("balance");
      }

      for (i = 0; i<passes ; i++ )
	{
	  RSS( ref_adapt_pass( ref_grid ), "pass");
	  ref_mpi_stopwatch_stop("pass");
	  RSS(ref_validation_cell_volume(ref_grid),"vol");
	  RSS(ref_migrate_to_balance(ref_grid),"balance");
	  ref_mpi_stopwatch_stop("balance");
	}

      ref_mpi_stopwatch_start();
      RSS( ref_gather_b8_ugrid( ref_grid, "ref_adapt_test.b8.ugrid" ), 
	   "gather");
      ref_mpi_stopwatch_stop("gather");

      RSS( ref_export_tec_ratio(ref_grid,"ref_adapt_post"),"export m");

      RSS( ref_grid_free( ref_grid ), "free");

      if ( ref_mpi_master )
	{
	  RSS(ref_import_by_extension( &ref_grid, 
				       "ref_adapt_test.b8.ugrid" ), "imp" );
	  RSS(ref_export_by_extension( ref_grid, "ref_adapt_test.tec" ),"ex" );
	  RSS( ref_grid_free( ref_grid ), "free");
	}
      ref_mpi_stopwatch_stop("post");
    }

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
