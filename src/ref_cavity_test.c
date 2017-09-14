#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_cavity.h"

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
#include "ref_metric.h"
#include "ref_clump.h"
#include "ref_geom.h"

int main( int argc, char *argv[] )
{

  { /* init 2 */
    REF_CAVITY ref_cavity;
    REIS(REF_NULL, ref_cavity_free(NULL),"dont free NULL");
    RSS(ref_cavity_create(&ref_cavity,2),"create");
    REIS( 0, ref_cavity_n(ref_cavity), "init no cavity");
    REIS( 2, ref_cavity_node_per(ref_cavity), "init per");
    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add face increments count */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert");
    REIS( 1, ref_cavity_n(ref_cavity), "init no cavity");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add faces, force realloc */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];
    REF_INT f, n;

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    n = ref_cavity_max(ref_cavity) + 3;
    for (f = 0; f<n; f++)
      {
        nodes[0] = f; nodes[1] = f+1;
        RSS(ref_cavity_insert(ref_cavity,nodes),"insert");
        REIS( f+1, ref_cavity_n(ref_cavity), "init no cavity");
      }

    REIS( n, ref_cavity_n(ref_cavity), "count");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add same face 2, raise error */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    REIS(REF_INVALID,ref_cavity_insert(ref_cavity,nodes),"insert second");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add same face 3, raise error */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];

    RSS(ref_cavity_create(&ref_cavity,3),"create");

    nodes[0] = 1; nodes[1] = 2; nodes[2] = 3;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    REIS(REF_INVALID,ref_cavity_insert(ref_cavity,nodes),"insert second");

    nodes[0] = 2; nodes[1] = 3; nodes[2] = 1;
    REIS(REF_INVALID,ref_cavity_insert(ref_cavity,nodes),"insert second");

    nodes[0] = 3; nodes[1] = 1; nodes[2] = 2;
    REIS(REF_INVALID,ref_cavity_insert(ref_cavity,nodes),"insert second");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add opposite face 2, mutual destruction */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    nodes[0] = 2; nodes[1] = 1;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert opposite");

    REIS( 0, ref_cavity_n(ref_cavity), "cancel");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add opposite face 3, mutual destruction */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];

    RSS(ref_cavity_create(&ref_cavity,3),"create");
    nodes[0] = 1; nodes[1] = 2; nodes[2] = 3;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");
    nodes[0] = 1; nodes[1] = 3; nodes[2] = 2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert opposite");

    REIS( 0, ref_cavity_n(ref_cavity), "cancel");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* find face 2 */
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");

    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0] = 2; nodes[1] = 1;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");

    nodes[0] = 3; nodes[1] = 4;
    REIS(REF_NOT_FOUND,ref_cavity_find(ref_cavity,nodes,
                                       &face,&reversed),"missing");
    REIS(REF_EMPTY,face,"found");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* find face 3 */
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT face;
    REF_BOOL reversed;

    RSS(ref_cavity_create(&ref_cavity,3),"create");
    nodes[0] = 1; nodes[1] = 2; nodes[2] = 3;
    RSS(ref_cavity_insert(ref_cavity,nodes),"insert first");

    nodes[0] = 1; nodes[1] = 2; nodes[2] = 3;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");
    nodes[0] = 3; nodes[1] = 1; nodes[2] = 2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");
    nodes[0] = 2; nodes[1] = 3; nodes[2] = 1;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find same");
    REIS(0,face,"found");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0] = 2; nodes[1] = 1; nodes[2] = 3;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");
    nodes[0] = 3; nodes[1] = 2; nodes[2] = 1;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");
    nodes[0] = 1; nodes[1] = 3; nodes[2] = 2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find reversed");
    REIS(0,face,"found");
    REIS(REF_TRUE,reversed,"not rev");

    nodes[0] = 3; nodes[1] = 4; nodes[2] = 5;
    REIS(REF_NOT_FOUND,ref_cavity_find(ref_cavity,nodes,
                                       &face,&reversed),"missing");
    REIS(REF_EMPTY,face,"found");

    RSS(ref_cavity_free(ref_cavity),"free");
  }

  { /* add triangle adds faces*/
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT nodes[2];
    REF_INT face;
    REF_BOOL reversed;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,1),"insert first");

    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find 0");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0] = 2; nodes[1] = 0;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find 1");
    REIS(REF_FALSE,reversed,"not rev");

    nodes[0] = 1; nodes[1] = 2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find 2");
    REIS(REF_FALSE,reversed,"not rev");

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* add rm triangle counts */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,1),"insert first");
    REIS( 3, ref_cavity_n(ref_cavity), "n" );
    REIS( 1, ref_list_n(ref_cavity_list(ref_cavity)), "l" );
    RSS(ref_cavity_rm_tri(ref_cavity,ref_grid,1),"insert first");
    REIS( 0, ref_cavity_n(ref_cavity), "n" );
    REIS( 0, ref_list_n(ref_cavity_list(ref_cavity)), "l" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* add tet */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;
    REF_INT nodes[3];
    REF_INT face;
    REF_BOOL reversed;

    RSS( ref_fixture_tet_grid( &ref_grid ), "pri" );

    RSS(ref_cavity_create(&ref_cavity,3),"create");

    RSS(ref_cavity_add_tet(ref_cavity,ref_grid,0),"insert first");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
    RSS(ref_cavity_find(ref_cavity,nodes,&face,&reversed),"find 0");
    REIS(REF_FALSE,reversed,"not rev");

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* insert tri node */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT global, node, clone;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,0),"insert first");

    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "new node");
    ref_node_xyz(ref_node,0,node) = 0.2;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 0.3;
    RSS( ref_node_twod_clone( ref_node, node, &clone ), "new node");

    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, clone ),"free");

    REIS( 8, ref_node_n(ref_grid_node(ref_grid)), "nodes" );
    REIS( 6, ref_cell_n(ref_grid_tri(ref_grid)), "nodes" );
    REIS( 3, ref_cell_n(ref_grid_pri(ref_grid)), "nodes" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* insert tet node */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT global, node;

    RSS( ref_fixture_tet_grid( &ref_grid ), "pri" );
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_cavity_create(&ref_cavity,3),"create");
    RSS(ref_cavity_add_tet(ref_cavity,ref_grid,0),"insert first");

    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "new node");
    ref_node_xyz(ref_node,0,node) = 0.1;
    ref_node_xyz(ref_node,1,node) = 0.2;
    ref_node_xyz(ref_node,2,node) = 0.3;

    RSS(ref_cavity_replace_tet(ref_cavity, ref_grid, node ),"free");

    REIS( 5, ref_node_n(ref_grid_node(ref_grid)), "nodes" );
    REIS( 4, ref_cell_n(ref_grid_tet(ref_grid)), "nodes" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* visible tet */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT global, node, face;
    REF_BOOL visible;

    RSS( ref_fixture_tet_grid( &ref_grid ), "pri" );
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_cavity_create(&ref_cavity,3),"create");
    RSS(ref_cavity_add_tet(ref_cavity,ref_grid,0),"insert first");

    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "new node");

    ref_node_xyz(ref_node,0,node) = 0.1;
    ref_node_xyz(ref_node,1,node) = 0.2;
    ref_node_xyz(ref_node,2,node) = 0.3;
    face = 0;
    RSS(ref_cavity_visible(ref_cavity, ref_node, node, face, &visible ),"free");
    REIS( REF_TRUE, visible, "vis" );

    ref_node_xyz(ref_node,0,node) = 1.0;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 1.0;
    face = 0;
    RSS(ref_cavity_visible(ref_cavity, ref_node, node, face, &visible ),"free");
    REIS( REF_FALSE, visible, "vis" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* visible twod */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT global, node, face;
    REF_BOOL visible;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,0),"insert first");

    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "new node");
    RSS(ref_metric_unit_node( ref_node ), "unit metric");

    ref_node_xyz(ref_node,0,node) = 0.2;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 0.3;
    face = 0;
    RSS(ref_cavity_visible(ref_cavity, ref_node, node, face, &visible ),"free");
    REIS( REF_TRUE, visible, "vis" );

    ref_node_xyz(ref_node,0,node) = 0.6;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 0.6;
    face = 0;
    RSS(ref_cavity_visible(ref_cavity, ref_node, node, face, &visible ),"free");
    REIS( REF_FALSE, visible, "vis" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* gobble */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT node, opp;

    RSS( ref_fixture_twod_brick_grid( &ref_grid ), "brick" );
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_metric_unit_node( ref_node ), "unit metric");
    REIS( 32, ref_node_n(ref_grid_node(ref_grid)), "nodes" );

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    node = 1;
    RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");

    RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");

    ref_node_xyz(ref_node,2,node) = 0.5;
    ref_node_xyz(ref_node,2,opp ) = 0.5;

    RSS(ref_cavity_enlarge_visible(ref_cavity,ref_grid,node),"insert first");
    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, opp ),"free");

    REIS( 30, ref_node_n(ref_grid_node(ref_grid)), "nodes" );
    REIS( 32, ref_cell_n(ref_grid_tri(ref_grid)), "nodes" );
    REIS( 16, ref_cell_n(ref_grid_pri(ref_grid)), "nodes" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS( ref_grid_free(ref_grid),"free");
  }

  { /* gobble gooble*/
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT node, opp;

    RSS( ref_fixture_twod_brick_grid( &ref_grid ), "brick" );
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_metric_unit_node( ref_node ), "unit metric");

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    node = 1;
    RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");
    RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");
    ref_node_xyz(ref_node,2,node) = 0.5;
    ref_node_xyz(ref_node,2,opp ) = 0.5;
    RSS(ref_cavity_enlarge_visible(ref_cavity,ref_grid,node),"insert first");
    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, opp ),"free");
    RSS(ref_cavity_free(ref_cavity),"free");

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    node = 2;
    RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");
    RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");
    ref_node_xyz(ref_node,2,node) = 0.5;
    ref_node_xyz(ref_node,2,opp ) = 0.5;
    RSS(ref_cavity_enlarge_visible(ref_cavity,ref_grid,node),"insert first");
    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, opp ),"free");
    RSS(ref_cavity_free(ref_cavity),"free");

    REIS( 28, ref_node_n(ref_grid_node(ref_grid)), "nodes" );
    REIS( 28, ref_cell_n(ref_grid_tri(ref_grid)), "nodes" );
    REIS( 14, ref_cell_n(ref_grid_pri(ref_grid)), "nodes" );

    RSS( ref_grid_free(ref_grid),"free");
  }

  { /* insert and remove tri node */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT global, node, clone, opp;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri" );
    ref_node = ref_grid_node(ref_grid);

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,0),"insert first");

    RSS( ref_node_next_global( ref_node, &global ), "next global");
    RSS( ref_node_add( ref_node, global, &node ), "new node");
    ref_node_xyz(ref_node,0,node) = 0.1;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 0.1;
    RSS( ref_node_twod_clone( ref_node, node, &clone ), "new node");
    RSS(ref_metric_unit_node( ref_node ), "unit metric");

    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, clone ),"free");
    RSS(ref_cavity_free(ref_cavity),"free");

    RSS(ref_cavity_create(&ref_cavity,2),"create");

    node = 0;
    RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");

    RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");

    ref_node_xyz(ref_node,0,node) = 0.3;
    ref_node_xyz(ref_node,0,opp ) = 0.3;
    ref_node_xyz(ref_node,2,node) = 0.3;
    ref_node_xyz(ref_node,2,opp ) = 0.3;

    RSS(ref_cavity_enlarge_visible(ref_cavity,ref_grid,node),"insert first");
    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, opp ),"free");

    RSS(ref_cavity_free(ref_cavity),"free");

    REIS( 6, ref_node_n(ref_grid_node(ref_grid)), "nodes" );
    REIS( 2, ref_cell_n(ref_grid_tri(ref_grid)), "nodes" );
    REIS( 1, ref_cell_n(ref_grid_pri(ref_grid)), "nodes" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* enlarge shrink twod face */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS( ref_fixture_twod_brick_grid( &ref_grid ), "brick" );

    RSS(ref_cavity_create(&ref_cavity,2),"create");
    RSS(ref_cavity_add_tri(ref_cavity,ref_grid,8),"insert first tri");
    REIS( 3, ref_cavity_n(ref_cavity), "n" );
    RSS(ref_cavity_enlarge_face(ref_cavity,ref_grid,1),"enl face 1");
    REIS( 4, ref_cavity_n(ref_cavity), "n" );
    RSS(ref_cavity_shrink_face(ref_cavity,ref_grid,3),"insert first tri");
    REIS( 3, ref_cavity_n(ref_cavity), "n" );
    REIS( 1, ref_list_n(ref_cavity_list(ref_cavity)), "l" );

    RSS(ref_cavity_free(ref_cavity),"free");
    RSS( ref_grid_free(ref_grid),"free");
  }

  { /* enlarge shrink threed face */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS( ref_fixture_tet2_grid( &ref_grid ), "brick" );

    RSS(ref_cavity_create(&ref_cavity,3),"create");
    RSS(ref_cavity_add_tet(ref_cavity,ref_grid,0),"insert first tri");
    REIS( 4, ref_cavity_n(ref_cavity), "n" );
    REIS( 1, ref_list_n(ref_cavity_list(ref_cavity)), "l" );
    RSS(ref_cavity_enlarge_face(ref_cavity,ref_grid,0),"enl face 1");
    REIS( 6, ref_cavity_n(ref_cavity), "n" );
    REIS( 2, ref_list_n(ref_cavity_list(ref_cavity)), "l" );
    RSS(ref_cavity_shrink_face(ref_cavity,ref_grid,5),"insert first tri");
    REIS( 4, ref_cavity_n(ref_cavity), "n" );
    REIS( 1, ref_list_n(ref_cavity_list(ref_cavity)), "l" );
    RSS(ref_cavity_free(ref_cavity),"free");
    RSS( ref_grid_free(ref_grid),"free");
  }

  { /* small edge */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_CAVITY ref_cavity;
    REF_INT node, opp;

    RSS( ref_fixture_twod_brick_grid( &ref_grid ), "brick" );
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_metric_unit_node( ref_node ), "unit metric");

    node = 10;
    RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");
    ref_node_xyz(ref_node,2,node) = 0.6;
    ref_node_xyz(ref_node,2,opp ) = 0.6;
    RSS(ref_cavity_create(&ref_cavity,2),"create");
    RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");
    RSS(ref_cavity_enlarge_metric(ref_cavity,ref_grid,node),"enlarge short");
    RSS(ref_cavity_replace_tri(ref_cavity, ref_grid, node, opp ),"free");
    RSS(ref_cavity_free(ref_cavity),"free");

    if ( 2 == argc )
      RSS( ref_export_by_extension( ref_grid, argv[1] ), "export" );

    RSS( ref_grid_free(ref_grid),"free");
  }

  { /* split edge of tet */
    REF_GRID ref_grid;
    REF_CAVITY ref_cavity;

    RSS( ref_fixture_tet_grid( &ref_grid ), "pri" );
    RSS(ref_cavity_create(&ref_cavity,3),"create");

    RSS(ref_cavity_add_edge(ref_cavity,ref_grid,1,2),"insert edge");
    REIS( 4, ref_cavity_n(ref_cavity), "n" );
    REIS( 1, ref_list_n(ref_cavity_list(ref_cavity)), "l" );
    
    RSS(ref_cavity_free(ref_cavity),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  
  return 0;
}
