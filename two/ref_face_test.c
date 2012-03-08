#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"

#include "ref_face.h"
#include "ref_mpi.h"

#include "ref_test.h"

int main( void )
{

  {  /* make faces shared by two elements */
    REF_FACE ref_face;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell;
    REF_INT node;

    RSS(ref_grid_create(&ref_grid),"create");

    for( node=0; node<6; node++)
      nodes[node] = node;

    RSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

    nodes[0] = 3; nodes[1] = 4; nodes[2] = 5; nodes[3] = 6;
    RSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add tet");

    RSS(ref_face_create(&ref_face,ref_grid),"create");

    REIS( 8, ref_face_n(ref_face), "check total faces");

    RSS(ref_face_free(ref_face),"face");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* find face with rotation */
    REF_FACE ref_face;
    REF_GRID ref_grid;
    REF_INT nodes[4];
    REF_INT cell;
    REF_INT node;
    REF_INT face;

    RSS(ref_grid_create(&ref_grid),"create");

    for( node=0; node<4; node++)
      nodes[node] = node;

    RSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add pri");

    RSS(ref_face_create(&ref_face,ref_grid),"create");

    REIS( 4, ref_face_n(ref_face), "check total faces");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 0;
    RSS( ref_face_with( ref_face, nodes, &face ), "with");
    REIS( 3, face, "missing face");	 

    nodes[0] = 1; nodes[1] = 2; nodes[2] = 0; nodes[3] = 1;
    RSS( ref_face_with( ref_face, nodes, &face ), "with rotation");
    REIS( 3, face, "missing face");	 

    RSS(ref_face_free(ref_face),"face");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* find face spanning two nodes */
    REF_FACE ref_face;
    REF_GRID ref_grid;
    REF_INT nodes[8];
    REF_INT cell;
    REF_INT node, node0, node1, face;

    RSS(ref_grid_create(&ref_grid),"create");

    for( node=0; node<8; node++)
      nodes[node] = node;
    
    RSS(ref_cell_add( ref_grid_hex(ref_grid), nodes, &cell ), "add pri");

    RSS(ref_face_create(&ref_face,ref_grid),"create");

    REIS( 6, ref_face_n(ref_face), "check total faces");

    node0 = 1; node1 = 6;
    RSS( ref_face_spanning( ref_face, node0, node1, &face ), "span" );
    REIS( 1, face, "wrong face");

    node0 = 4; node1 = 3;
    RSS( ref_face_spanning( ref_face, node0, node1, &face ), "span" );
    REIS( 3, face, "wrong face");

    node0 = 0; node1 = 2;
    RSS( ref_face_spanning( ref_face, node0, node1, &face ), "span" );
    REIS( 4, face, "wrong face");

    RSS(ref_face_free(ref_face),"face");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* quad normal */
    REF_DBL xyz0[3] = {0.0,0.0,0.0};
    REF_DBL xyz1[3] = {1.0,0.0,0.0};
    REF_DBL xyz2[3] = {1.0,1.0,0.0};
    REF_DBL xyz3[3] = {0.0,1.0,0.0};
    REF_DBL normal[3];

    RSS( ref_face_normal( xyz0, xyz1, xyz2, xyz3, normal ), "norm");
    
    TWDS( 0.0, normal[0], -1.0, "x-norm");
    TWDS( 0.0, normal[1], -1.0, "y-norm");
    TWDS( 1.0, normal[2], -1.0, "z-norm");
  }

  {  /* tri normal */
    REF_DBL xyz0[3] = {0.0,0.0,0.0};
    REF_DBL xyz1[3] = {1.0,0.0,0.0};
    REF_DBL xyz2[3] = {1.0,1.0,0.0};
    REF_DBL xyz3[3] = {0.0,0.0,0.0};
    REF_DBL normal[3];

    RSS( ref_face_normal( xyz0, xyz1, xyz2, xyz3, normal ), "norm");
    
    TWDS( 0.0, normal[0], -1.0, "x-norm");
    TWDS( 0.0, normal[1], -1.0, "y-norm");
    TWDS( 0.5, normal[2], -1.0, "z-norm");
  }

  {  /* quad open node 1 */
    REF_DBL xyz0[3] = {0.0,0.0,0.0};
    REF_DBL xyz1[3] = {0.5,0.0,0.0};
    REF_DBL xyz2[3] = {1.0,0.0,0.0};
    REF_DBL xyz3[3] = {0.0,1.0,0.0};
    REF_INT open_node;

    RSS( ref_face_open_node( xyz0, xyz1, xyz2, xyz3, &open_node ), "open node");
    
    REIS( 1, open_node, "expect 180 deg node");
  }

  {  /* quad open node 2 */
    REF_DBL xyz0[3] = {0.0,0.0,0.0};
    REF_DBL xyz1[3] = {0.1,0.0,0.0};
    REF_DBL xyz2[3] = {0.5,0.5,0.0};
    REF_DBL xyz3[3] = {0.0,1.0,0.0};
    REF_INT open_node;

    RSS( ref_face_open_node( xyz0, xyz1, xyz2, xyz3, &open_node ), "open node");
    
    REIS( 2, open_node, "expect 180 deg node");
  }

  return 0;
}
