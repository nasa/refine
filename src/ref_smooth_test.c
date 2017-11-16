
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
#include <math.h>

#include "ref_smooth.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_list.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_math.h"
#include "ref_geom.h"

#include "ref_import.h"
#include "ref_export.h"
#include "ref_part.h"
#include "ref_histogram.h"

#include "ref_edge.h"
#include "ref_dict.h"
#include "ref_migrate.h"
#include "ref_gather.h"
#include "ref_adapt.h"
#include "ref_collapse.h"
#include "ref_split.h"

#include "ref_metric.h"

#include "ref_fixture.h"

#include "ref_mpi.h"

#include "ref_malloc.h"

#include "ref_twod.h"

static REF_STATUS ref_smooth_tri_single_fixture( REF_GRID *ref_grid_ptr, 
						 REF_INT *target_node, 
						 REF_INT *target_cell )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS( ref_grid_create( ref_grid_ptr ), "grid" );
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /*
   0    z
   |\   |
   1-2  +-x
   */

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[0] = node;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[1] = node;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[2] = node;

  nodes[3] = 10;

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  *target_node = nodes[0];
  *target_cell = cell;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[3] = node;

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[4] = node;

  RSS(ref_node_add(ref_node,5,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[5] = node;

  RSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add pri");

  RSS( ref_metric_unit_node( ref_node ), "unit node" );

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_tri_two_fixture( REF_GRID *ref_grid_ptr, 
					      REF_INT *target_node )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS( ref_grid_create( ref_grid_ptr ), "grid" );
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /*
   1-0  z
   |/|  |
   2-3  +-x
   */

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[0] = node;

  *target_node = node;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[1] = node;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[2] = node;

  nodes[3] = 10;

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[3] = node;

  RSS(ref_node_add(ref_node,5,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[4] = node;

  RSS(ref_node_add(ref_node,6,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[5] = node;

  RSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add pri");

  nodes[1] = nodes[2];
  nodes[4] = nodes[5];

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[2] = node;

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  RSS(ref_node_add(ref_node,7,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[5] = node;

  RSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add pri");

  RSS( ref_metric_unit_node( ref_node ), "unit node" );

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_tet_two_fixture( REF_GRID *ref_grid_ptr, 
					      REF_INT *target_node,
					      REF_INT *top_node )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT node;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS( ref_grid_create( ref_grid_ptr ), "grid" );
  ref_grid = *ref_grid_ptr;
  ref_node = ref_grid_node(ref_grid);

  /* 0-1-2 base, 0-2-4 wall, 3 free */

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[0] = node;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[1] = node;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;
  nodes[2] = node;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[3] = node;

  *target_node = node;

  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tri");

  nodes[1] = nodes[2];

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;
  nodes[2] = node;

  *top_node = node;

  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tri");

  RSS( ref_metric_unit_node( ref_node ), "unit node" );

  return REF_SUCCESS;
}

static REF_STATUS ref_smooth_tet_tri_fixture( REF_GRID *ref_grid_ptr, 
					      REF_INT *target_node )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, cell;
  
  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  *target_node = node;
  
  RSS(ref_node_add(ref_node,1,&node),"node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,2,&node),"node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"node");
  ref_node_xyz(ref_node,0,node) =-1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,4,&node),"node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) =-1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,5,&node),"node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  /*  2

  3   01  5

      4  */

  nodes[0] = 0; nodes[3] = 1; /* common tet top node and face id */

  nodes[1] = 2; nodes[2] = 3;
  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tet");

  nodes[1] = 3; nodes[2] = 4;
  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tet");

  nodes[1] = 4; nodes[2] = 5;
  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tet");

  nodes[1] = 5; nodes[2] = 2;
  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tet");

  RSS( ref_metric_unit_node( ref_node ), "unit node" );

  return REF_SUCCESS;
}

int main( int argc, char *argv[] )
{

  if ( argc > 2 )
    {
      REF_GRID ref_grid;

      RSS( ref_mpi_start( argc, argv ), "start" );

      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );

      RSS( ref_part_metric( ref_grid_node(ref_grid), argv[2] ), "get metric");

      RSS( ref_histogram_quality( ref_grid ), "qual");

      RSS( ref_export_tec_surf( ref_grid, "ref_smooth_test_0.tec" ),
           "surf");

      RSS( ref_grid_free( ref_grid ), "free");

      RSS( ref_mpi_stop( ), "stop" );
    }

  if ( argc == 2 )
    {
      REF_GRID ref_grid;
      REF_INT target_node;

      RSS( ref_mpi_start( argc, argv ), "start" );

      RSS( ref_smooth_tet_tri_fixture( &ref_grid, &target_node ), "fix" );
      
      RSS( ref_export_by_extension( ref_grid, argv[1] ),
           "fixture");

      RSS( ref_grid_free( ref_grid ), "free");

      RSS( ref_mpi_stop( ), "stop" );
    }

  {   
    REF_GRID ref_grid;
    REF_INT node, cell;
    RSS( ref_smooth_tri_single_fixture( &ref_grid, &node, &cell ), "2d fix" );

    RSS( ref_smooth_tri_steepest_descent( ref_grid, node ), "smooth" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* ideal in unit metric  */
    REF_GRID ref_grid;
    REF_INT node, cell;
    REF_DBL ideal[3];
    RSS( ref_smooth_tri_single_fixture( &ref_grid, &node, &cell ), "2d fix" );

    RSS(ref_smooth_tri_ideal(ref_grid, node, cell, ideal),"ideal");
    RWDS(0.5,ideal[0],-1,"ideal x");
    RWDS(1.0,ideal[1],-1,"ideal y");
    RWDS(0.5*sqrt(3.0),ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* ideal in hx=2.0 metric  */
    REF_GRID ref_grid;
    REF_INT node, cell;
    REF_DBL ideal[3];

    RSS( ref_smooth_tri_single_fixture( &ref_grid, &node, &cell ), "2d fix" );

    ref_node_metric( ref_grid_node(ref_grid), 0, node ) = 0.25;

    RSS(ref_smooth_tri_ideal(ref_grid, node, cell, ideal),"ideal");
    RWDS(0.5,ideal[0],-1,"ideal x");
    RWDS(1.0,ideal[1],-1,"ideal y");
    RWDS(0.5*sqrt(3.0),ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* ideal in hz=2.0 metric  */
    REF_GRID ref_grid;
    REF_INT node, cell;
    REF_DBL ideal[3];

    RSS( ref_smooth_tri_single_fixture( &ref_grid, &node, &cell ), "2d fix" );

    ref_node_metric( ref_grid_node(ref_grid), 5, node ) = 0.25;

    RSS(ref_smooth_tri_ideal(ref_grid, node, cell, ideal),"ideal");
    RWDS(0.5,ideal[0],-1,"ideal x");
    RWDS(1.0,ideal[1],-1,"ideal y");
    RWDS(2.0*0.5*sqrt(3.0),ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

 { /* weighted ideal */
    REF_GRID ref_grid;
    REF_INT node;
    REF_DBL ideal[3];

    RSS( ref_smooth_tri_two_fixture( &ref_grid, &node ), "2d fix" );

    RSS(ref_smooth_tri_weighted_ideal(ref_grid, node, ideal),"ideal");
    RWDS(0.5*(0.5+0.5*sqrt(3.0)),ideal[0],-1,"ideal x");
    RWDS(1.0,ideal[1],-1,"ideal y");
    RWDS(0.5*(0.5+0.5*sqrt(3.0)),ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

 { /* set to new ideal */
    REF_GRID ref_grid;
    REF_INT node, opposite;
    REF_DBL quality0, quality1;
    REF_BOOL allowed;

    RSS( ref_smooth_tri_two_fixture( &ref_grid, &node ), "2d fix" );

    ref_node_xyz( ref_grid_node(ref_grid), 2, node ) = 0.0000001;
    RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid), 
				 node, &opposite), "opp");
    ref_node_xyz( ref_grid_node(ref_grid), 2, opposite ) =
      ref_node_xyz( ref_grid_node(ref_grid), 2, node );

    ref_node_xyz( ref_grid_node(ref_grid), 0, 1 ) = 1.0;
    ref_node_xyz( ref_grid_node(ref_grid), 2, 1 ) = 0.5;

    RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid), 
				 1, &opposite), "opp");
    ref_node_xyz( ref_grid_node(ref_grid), 0, opposite ) =
      ref_node_xyz( ref_grid_node(ref_grid), 0, 1 );
    ref_node_xyz( ref_grid_node(ref_grid), 2, opposite ) =
      ref_node_xyz( ref_grid_node(ref_grid), 2, 1 );

    RSS( ref_smooth_tri_quality_around( ref_grid, node, &quality0),"q");

    RSS( ref_smooth_tri_improve( ref_grid, node ),"imp");

    RSS( ref_smooth_tri_quality_around( ref_grid, node, &quality1),"q");

    RAS( quality1 > quality0, "expected improvment");

    RSS( ref_smooth_outward_norm( ref_grid, node, &allowed ),"outward allowed");
    RAS( allowed, "expected validity");

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* ideal tet in unit metric  */
    REF_GRID ref_grid;
    REF_INT node, cell;
    REF_DBL ideal[3];
    RSS( ref_fixture_tet_grid( &ref_grid ), "2d fix" );
    node = 3;
    cell = 0;
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "unit node" );

    ref_node_xyz( ref_grid_node(ref_grid), 0, 2 ) = 0.5;
    ref_node_xyz( ref_grid_node(ref_grid), 1, 2 ) = 0.5*sqrt(3.0);

    RSS(ref_smooth_tet_ideal(ref_grid, node, cell, ideal),"ideal");
    RWDS(0.5,               ideal[0],-1,"ideal x");
    RWDS(1.0/6.0*sqrt(3.0), ideal[1],-1,"ideal y");
    RWDS(1.0/3.0*sqrt(6.0), ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* ideal tet in hz=2.0 metric  */
    REF_GRID ref_grid;
    REF_INT node, cell;
    REF_DBL ideal[3];
    REF_DBL hz;
    RSS( ref_fixture_tet_grid( &ref_grid ), "2d fix" );
    node = 3;
    cell = 0;
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "unit node" );

    hz = 2.0;
    ref_node_metric( ref_grid_node(ref_grid), 5, node ) = 1.0/hz/hz;

    ref_node_xyz( ref_grid_node(ref_grid), 0, 2 ) = 0.5;
    ref_node_xyz( ref_grid_node(ref_grid), 1, 2 ) = 0.5*sqrt(3.0);

    RSS(ref_smooth_tet_ideal(ref_grid, node, cell, ideal),"ideal");
    RWDS(0.5,               ideal[0],-1,"ideal x");
    RWDS(1.0/6.0*sqrt(3.0), ideal[1],-1,"ideal y");
    RWDS(hz*1.0/3.0*sqrt(6.0), ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

 { /* weighted tet ideal */
    REF_GRID ref_grid;
    REF_INT node, top_node;
    REF_DBL ideal[3];

    RSS( ref_smooth_tet_two_fixture( &ref_grid, &node, &top_node ), "3d 2tet" );

    RSS(ref_smooth_tet_weighted_ideal(ref_grid, node, ideal),"ideal");
    RWDS(0.574914957130530,ideal[0],-1,"ideal x");
    RWDS(1.0/3.0,ideal[1],-1,"ideal y");
    RWDS(0.574914957130530,ideal[2],-1,"ideal z");

    RSS(ref_grid_free(ref_grid),"free");
  }

 { /* weighted tet ideal */
    REF_GRID ref_grid;
    REF_INT node, top_node;
    REF_DBL quality0, quality1;

    RSS( ref_smooth_tet_two_fixture( &ref_grid, &node, &top_node ), "3d 2tet" );

    ref_node_xyz( ref_grid_node(ref_grid), 0, node ) = 1.0;
    ref_node_xyz( ref_grid_node(ref_grid), 1, node ) = 0.0;
    ref_node_xyz( ref_grid_node(ref_grid), 2, node ) = 0.0001;

    ref_node_xyz( ref_grid_node(ref_grid), 0, top_node ) = 0.8;
    ref_node_xyz( ref_grid_node(ref_grid), 1, top_node ) = 0.0;
    ref_node_xyz( ref_grid_node(ref_grid), 2, top_node ) = 0.5;

    RSS( ref_smooth_tet_quality_around( ref_grid, node, &quality0),"q");

    RSS( ref_smooth_tet_improve( ref_grid, node ),"imp");

    RSS( ref_smooth_tet_quality_around( ref_grid, node, &quality1),"q");

    RAS( quality1 > quality0, "expected improvment");
    RAS( quality1 > 0, "expected validity");

    RSS(ref_grid_free(ref_grid),"free");
  }

 { /* planar tri face on tet grid */
   REF_GRID ref_grid;
   REF_INT target_node;
   REF_INT cell;
   REF_DBL ideal[3];
   RSS( ref_smooth_tet_tri_fixture( &ref_grid, &target_node ), "fix" );

   ref_node_xyz(ref_grid_node(ref_grid),0,target_node) = 0.5;

   cell = 0;
   RSS(ref_smooth_tri_ideal(ref_grid, target_node, cell, ideal),"ideal");

   RWDS( 0.1123724356957946, ideal[0],-1,"ideal x");
   RWDS(-0.1123724356957946, ideal[1],-1,"ideal y");
   RWDS( 0.0,                ideal[2],-1,"ideal z");

   
   RSS( ref_grid_free( ref_grid ), "free");
 }

  return 0;
}
