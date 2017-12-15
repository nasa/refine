
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_interp.h"

#include "ref_mpi.h"
#include "ref_fixture.h"
#include "ref_grid.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include  "ref_node.h"
#include   "ref_matrix.h"
#include   "ref_math.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_geom.h"
#include   "ref_dict.h"
#include   "ref_export.h"
#include  "ref_gather.h"
#include  "ref_adapt.h"
#include   "ref_collapse.h"
#include    "ref_edge.h"
#include   "ref_split.h"
#include    "ref_cavity.h"
#include    "ref_twod.h"
#include   "ref_smooth.h"
#include    "ref_clump.h"
#include "ref_part.h"
#include  "ref_import.h"
#include  "ref_migrate.h"
#include "ref_metric.h"

REF_STATUS ref_interp_setup( REF_INTERP *ref_interp_ptr, REF_MPI ref_mpi )
{
  REF_GRID from, to;
  RSS(ref_grid_create(&from,ref_mpi),"create");
  RSS(ref_grid_create(&to,ref_mpi),"create");
  RSS( ref_interp_create( ref_interp_ptr, from, to ), "make interp" );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_teardown( REF_INTERP ref_interp )
{
  RSS(ref_grid_free(ref_interp_from_grid(ref_interp)),"free");
  RSS(ref_grid_free(ref_interp_to_grid(ref_interp)),"free");
  RSS( ref_interp_free( ref_interp ), "interp free" );
  return REF_SUCCESS;
}

int main( int argc, char *argv[] )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_start( argc, argv ), "start" );
  RSS( ref_mpi_create( &ref_mpi ), "make mpi" );

  if ( argc == 3 )
    {
      REF_GRID from, to;
      REF_INTERP ref_interp;

      RSS( ref_mpi_stopwatch_start( ref_mpi ), "sw start");
      RSS(ref_part_by_extension( &from, ref_mpi, argv[1] ), "import" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "from grid" ), "sw start");
      RSS(ref_part_by_extension( &to, ref_mpi, argv[2] ), "import" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "to grid" ), "sw start");

      RSS(ref_export_tec_surf( to, "ref_interp_test_to.tec" ),"export" );
      RSS(ref_export_tec_surf( from, "ref_interp_test_from.tec" ),"export" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "export viz" ), "sw start");
      
      RSS( ref_interp_create( &ref_interp, from, to ), "make interp" );
      ref_interp->instrument = REF_TRUE;
      RSS( ref_interp_locate(ref_interp), "map" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "locate" ), "sw start");
      RSS(ref_interp_tec( ref_interp, "ref_interp_test_exhaust.tec" ),"export" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "tec" ), "sw start");
      RSS( ref_interp_stats(ref_interp ), "err" );
      RSS( ref_mpi_stopwatch_stop( ref_mpi, "stats" ), "sw start");

      
      RSS( ref_interp_free( ref_interp ), "interp free" );
      RSS( ref_grid_free(to),"free");
      RSS( ref_grid_free(from),"free");

      RSS( ref_mpi_free( ref_mpi ), "mpi free" );
      RSS( ref_mpi_stop( ), "stop" );
      return 0;
    }

  {
    REF_INTERP ref_interp;
    RSS( ref_interp_setup( &ref_interp, ref_mpi ), "setup");
    RSS( ref_interp_teardown( ref_interp ), "teardown");
  }

  { /* bricks */
    REF_GRID from, to;
    REF_INT node;
    char file[] = "ref_interp_test.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if ( ref_mpi_once(ref_mpi) )
      {
	RSS( ref_fixture_tet_brick_grid( &from, ref_mpi ), "brick" );
	RSS(ref_export_by_extension( from, file ),"export" );
	RSS( ref_grid_free(from),"free");
      }
    RSS(ref_part_by_extension( &from, ref_mpi, file ), "import" );
    RSS(ref_part_by_extension( &to, ref_mpi, file ), "import" );
    if ( ref_mpi_once(ref_mpi) )
      REIS(0, remove( file ), "test clean up");

    each_ref_node_valid_node( ref_grid_node(to), node )
      if ( ( 0.01 < ref_node_xyz(ref_grid_node(to),0,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),0,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),1,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),1,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),2,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),2,node) ) )
	{
	  ref_node_xyz(ref_grid_node(to),0,node) += 1.0e-2;
	  ref_node_xyz(ref_grid_node(to),1,node) += 2.0e-2;
	  ref_node_xyz(ref_grid_node(to),2,node) += 4.0e-2;
	}

    RSS( ref_interp_create( &ref_interp, from, to ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
    if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 26, ref_interp->n_walk, "walk count" );
	REIS( 30, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.121 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 7.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_interp_create( &ref_interp, to, from ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
      if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 26, ref_interp->n_walk, "walk count" );
	REIS( 30, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.121 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 7.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_grid_free(to),"free");
    RSS( ref_grid_free(from),"free");
  }
  
  { /* odd split one brick */
    REF_GRID from, to;
    REF_INT node;
    char even[] = "ref_interp_test_even.meshb";
    char odd[] = "ref_interp_test_odd.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if ( ref_mpi_once(ref_mpi) )
      {
	REF_GRID ref_grid;

	RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "brick" );
	RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "m" );
	RSS( ref_export_by_extension( ref_grid, even ),"export" );
	RSS( ref_grid_free(ref_grid),"free");

	RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "brick" );
	RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "m" );
	RSS( ref_split_edge_pattern( ref_grid, 1, 2 ), "split" );
	RSS( ref_export_by_extension( ref_grid, odd ),"export" );
	RSS( ref_grid_free(ref_grid),"free");
      }
    RSS(ref_part_by_extension( &from, ref_mpi, even ), "import" );
    RSS(ref_part_by_extension( &to, ref_mpi, odd ), "import" );
    if ( ref_mpi_once(ref_mpi) )
      {
	REIS(0, remove( even ), "test clean up");
	REIS(0, remove( odd ), "test clean up");
      }

    each_ref_node_valid_node( ref_grid_node(to), node )
      if ( ( 0.01 < ref_node_xyz(ref_grid_node(to),0,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),0,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),1,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),1,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),2,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),2,node) ) )
	{
	  ref_node_xyz(ref_grid_node(to),0,node) += 1.0e-2;
	  ref_node_xyz(ref_grid_node(to),1,node) += 2.0e-2;
	  ref_node_xyz(ref_grid_node(to),2,node) += 4.0e-2;
	}

    RSS( ref_interp_create( &ref_interp, from, to ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
    if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 129, ref_interp->n_walk, "walk count" );
	REIS( 66, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.241 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 9.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_interp_create( &ref_interp, to, from ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
    if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 26, ref_interp->n_walk, "walk count" );
	REIS( 30, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.241 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 7.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_grid_free(to),"free");
    RSS( ref_grid_free(from),"free");
  }

  { /* odd/even split bricks */
    REF_GRID from, to;
    REF_INT node;
    char even[] = "ref_interp_test_even.meshb";
    char odd[] = "ref_interp_test_odd.meshb";
    REF_INTERP ref_interp;
    REF_DBL max_error, min_bary;

    if ( ref_mpi_once(ref_mpi) )
      {
	REF_GRID ref_grid;

	RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "brick" );
	RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "m" );
	RSS( ref_split_edge_pattern( ref_grid, 0, 2 ), "split" );
	RSS( ref_export_by_extension( ref_grid, even ),"export" );
	RSS( ref_grid_free(ref_grid),"free");

	RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "brick" );
	RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "m" );
	RSS( ref_split_edge_pattern( ref_grid, 1, 2 ), "split" );
	RSS( ref_export_by_extension( ref_grid, odd ),"export" );
	RSS( ref_grid_free(ref_grid),"free");
      }
    RSS(ref_part_by_extension( &from, ref_mpi, even ), "import" );
    RSS(ref_part_by_extension( &to, ref_mpi, odd ), "import" );
    if ( ref_mpi_once(ref_mpi) )
      {
	REIS(0, remove( even ), "test clean up");
	REIS(0, remove( odd ), "test clean up");
      }

    each_ref_node_valid_node( ref_grid_node(to), node )
      if ( ( 0.01 < ref_node_xyz(ref_grid_node(to),0,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),0,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),1,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),1,node) ) ||
	   ( 0.01 < ref_node_xyz(ref_grid_node(to),2,node) &&
	     0.99 > ref_node_xyz(ref_grid_node(to),2,node) ) )
	{
	  ref_node_xyz(ref_grid_node(to),0,node) += 1.0e-2;
	  ref_node_xyz(ref_grid_node(to),1,node) += 2.0e-2;
	  ref_node_xyz(ref_grid_node(to),2,node) += 4.0e-2;
	}

    RSS( ref_interp_create( &ref_interp, from, to ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
    if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 129, ref_interp->n_walk, "walk count" );
	REIS( 66, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.241 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 9.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_interp_create( &ref_interp, to, from ), "make interp" );
    RSS( ref_interp_locate(ref_interp), "map" );
    REIS( 8, ref_interp->n_geom, "geom missing" );
    REIS( 0, ref_interp->n_geom_fail, "geom fail" );
    if ( !ref_mpi_para(ref_mpi) )
      {
	REIS( 121, ref_interp->n_walk, "walk count" );
	REIS( 75, ref_interp->n_tree, "tree count" );	
      }
    RSS( ref_interp_min_bary(ref_interp, &min_bary), "min bary" );
    RAS( -0.241 < min_bary, "large extrapolation" );
    RSS( ref_interp_max_error(ref_interp, &max_error), "err" );
    RAS( 7.0e-16 > max_error, "large interp error" );
    RSS( ref_interp_free( ref_interp ), "interp free" );

    RSS( ref_grid_free(to),"free");
    RSS( ref_grid_free(from),"free");
  }

  RSS( ref_mpi_free( ref_mpi ), "mpi free" );
  RSS( ref_mpi_stop( ), "stop" );
  return 0;
}
