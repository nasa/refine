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

#include "ref_smooth.h"
#include "ref_collapse.h"
#include "ref_split.h"
#include "ref_edge.h"

#include "ref_subdiv.h"
#include "ref_validation.h"
#include "ref_face.h"

#include "ref_malloc.h"
#include "ref_metric.h"
#include "ref_math.h"

#include "ref_histogram.h"

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

      RSS(ref_validation_cell_volume(ref_grid),"vol");
      RSS( ref_histogram_ratio( ref_grid ), "gram");

      if ( 3 < argc )
      {
	REF_SUBDIV ref_subdiv;
	REF_DBL *node_ratio;

        RSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

	ref_malloc( node_ratio, ref_node_max(ref_node), REF_DBL );

	RSS(ref_part_ratio( ref_node, node_ratio, argv[3] ), "part metric" );
	ref_mpi_stopwatch_stop("read ratio");
	
	RSS(ref_subdiv_mark_prism_by_ratio(ref_subdiv, node_ratio),"mark rat");
	ref_mpi_stopwatch_stop("subdiv mark");
	
	ref_free( node_ratio );

	RSS(ref_subdiv_split(ref_subdiv),"split");
	ref_mpi_stopwatch_stop("subdiv split");

	RSS(ref_subdiv_free(ref_subdiv),"free");

	RSS(ref_validation_cell_volume(ref_grid),"vol");
	RSS(ref_migrate_to_balance(ref_grid),"balance");
	ref_mpi_stopwatch_stop("balance");
      }

      if ( REF_FALSE )
	{
	  REF_INT gradation;
	  RSS(ref_export_tec( ref_grid, "ref_adapt_ellipseg.tec" ),"ex" );
	  RSS(ref_export_tec_metric_ellipse( ref_grid, "ref_adapt_0" ),"ex" );
	  for ( gradation =0 ; gradation<10 ; gradation++ )
	    {
	      printf("gradation %d\n",gradation);
	      RSS( ref_metric_gradation( ref_grid ), "grad");
	      RSS( ref_node_ghost_real( ref_node ), "ghost real");
	    }
	  RSS(ref_export_tec_metric_ellipse( ref_grid, "ref_adapt_1" ),"ex" );
	}

      if (REF_FALSE)
	{
	  printf("limit with existing grid\n");
	  RSS(ref_export_tec( ref_grid, "ref_adapt_ellipseg.tec" ),"ex" );
	  RSS(ref_export_tec_metric_ellipse( ref_grid, "ref_adapt_0" ),"ex" );
	  RSS( ref_metric_sanitize(ref_grid),"sant");
	  RSS( ref_node_ghost_real( ref_node ), "ghost real");
	  RSS(ref_export_tec_metric_ellipse( ref_grid, "ref_adapt_1" ),"ex" );
	}

      RSS( ref_gather_tec_movie_record_button( REF_TRUE ), "rec" );

      RSS(ref_validation_cell_volume(ref_grid),"vol");
      RSS( ref_histogram_ratio( ref_grid ), "gram");

      passes = 10;
      for (i = 0; i<passes ; i++ )
	{
	  RSS( ref_adapt_pass( ref_grid ), "pass");
	  ref_mpi_stopwatch_stop("pass");
	  RSS(ref_validation_cell_volume(ref_grid),"vol");
	  RSS( ref_histogram_quality( ref_grid ), "gram");
	  RSS( ref_histogram_ratio( ref_grid ), "gram");
	  RSS(ref_migrate_to_balance(ref_grid),"balance");
	  ref_mpi_stopwatch_stop("balance");
	}

      ref_mpi_stopwatch_start();
      RSS( ref_gather_b8_ugrid( ref_grid, "ref_adapt_test.b8.ugrid" ), 
	   "gather");
      ref_mpi_stopwatch_stop("gather");

      RSS(ref_export_tec_metric_ellipse( ref_grid, "ref_adapt_f" ),"ex" );

      RSS( ref_grid_free( ref_grid ), "free");

      if ( ref_mpi_master )
	{
	  RSS(ref_import_by_extension( &ref_grid, 
				       "ref_adapt_test.b8.ugrid" ), "imp" );
	  RSS(ref_export_tec_surf( ref_grid, "ref_adapt_test.tec" ),"ex" );
	  RSS( ref_grid_free( ref_grid ), "free");
	}
      ref_mpi_stopwatch_stop("post");
    }

  if ( 2 == argc )
    {
      REF_GRID ref_grid;
      REF_INT i, passes;

      ref_migrate_method = REF_MIGRATE_RCB;

      ref_mpi_stopwatch_start();
      RSS(ref_part_b8_ugrid( &ref_grid, argv[1] ), "part grid" );
      ref_mpi_stopwatch_stop("read grid");

      RSS(ref_migrate_to_balance(ref_grid),"balance");
      ref_mpi_stopwatch_stop("balance");

      RSS(ref_validation_cell_volume(ref_grid),"vol");

      {
	REF_DBL *metric;
	REF_INT node;
	ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
	RSS( ref_metric_imply_from( metric, ref_grid ), 
	     "from");
	RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid)), "to");
	RSS( ref_node_ghost_real( ref_grid_node(ref_grid) ), "ghost real");
	ref_free( metric );

	each_ref_node_valid_node( ref_grid_node(ref_grid), node )
  	  {
	    REF_DBL scale = 0.25;
	    ref_node_metric(ref_grid_node(ref_grid),0,node) /= (scale*scale);
	    ref_node_metric(ref_grid_node(ref_grid),1,node)  = 0.0;
	    ref_node_metric(ref_grid_node(ref_grid),2,node) /= (scale*scale);
	    ref_node_metric(ref_grid_node(ref_grid),3,node)  = 1.0;
	    ref_node_metric(ref_grid_node(ref_grid),4,node)  = 0.0;
	    ref_node_metric(ref_grid_node(ref_grid),5,node) /= (scale*scale);
	  }
      }

      RSS( ref_gather_tec_movie_record_button( REF_TRUE ), "rec" );

      RSS(ref_validation_cell_volume(ref_grid),"vol");
      RSS( ref_histogram_ratio( ref_grid ), "gram");
      RSS( ref_histogram_quality( ref_grid ), "qual");

      passes = 10;
      for (i = 0; i<passes ; i++ )
	{
	  RSS( ref_adapt_pass( ref_grid ), "pass");
	  ref_mpi_stopwatch_stop("pass");
	  RSS(ref_validation_cell_volume(ref_grid),"vol");
	  RSS( ref_histogram_ratio( ref_grid ), "gram");
	  RSS( ref_histogram_quality( ref_grid ), "qual");
	  RSS(ref_migrate_to_balance(ref_grid),"balance");
	  ref_mpi_stopwatch_stop("balance");
	}

      ref_mpi_stopwatch_start();
      RSS( ref_gather_b8_ugrid( ref_grid, "ref_adapt_test.b8.ugrid" ), 
	   "gather");
      ref_mpi_stopwatch_stop("gather");

      RSS( ref_gather_tec_part( ref_grid, "ref_adapt_test_part.tec" ), 
	   "part_viz");

      RSS( ref_grid_free( ref_grid ), "free");

      if ( ref_mpi_master )
	{
	  RSS(ref_import_by_extension( &ref_grid, 
				       "ref_adapt_test.b8.ugrid" ), "imp" );
	  RSS(ref_export_tec_surf( ref_grid, "ref_adapt_test.tec" ),"ex" );
	  RSS( ref_grid_free( ref_grid ), "free");
	}
      ref_mpi_stopwatch_stop("post");
    }

  if ( 1 == argc )
    {
      REF_GRID ref_grid;
      REF_INT i, passes;

      RSS(ref_fixture_pri_grid(&ref_grid),"set up grid");
      ref_grid_twod(ref_grid) = REF_TRUE;

      RSS(ref_migrate_to_balance(ref_grid),"balance");

      {
	REF_DBL *metric;
	ref_malloc( metric, 6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
	RSS( ref_metric_imply_from( metric, ref_grid ), 
	     "from");
	RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid)), "to");
	RSS( ref_node_ghost_real( ref_grid_node(ref_grid) ), "ghost real");
	ref_free( metric );
      }

      passes = 5;
      for (i = 0; i<passes ; i++ )
	{
	  RSS( ref_adapt_pass( ref_grid ), "pass");
	  RSS(ref_migrate_to_balance(ref_grid),"balance");
	}

      RSS( ref_grid_free( ref_grid ), "free");

    }

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
