
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

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include  "ref_math.h"

#include   "ref_fixture.h"

#include  "ref_export.h"
#include   "ref_dict.h"
#include   "ref_edge.h"

#include "ref_malloc.h"

int main( void )
{
  REF_MPI ref_mpi;
  RSS( ref_mpi_create( &ref_mpi ), "create" );
 
  {  /* init */
    REF_GRID ref_grid;
    REIS(REF_NULL,ref_grid_free(NULL),"dont free NULL");

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* each element */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;
    REF_INT group;

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");

    node_per = 3;
    each_ref_grid_ref_cell( ref_grid, group, ref_cell )
      {
	node_per += 1;
	if ( 7 == node_per ) node_per = 8;
	RES( node_per, ref_cell_node_per( ref_cell), "cells in order" );
      }

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* cell with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");

    node_per = 4;
    RSS( ref_grid_cell_with( ref_grid, node_per, &ref_cell ), "with" );
    REIS( node_per, ref_cell_node_per( ref_cell ), "match" );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* face with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");

    node_per = 4;
    RSS( ref_grid_face_with( ref_grid, node_per, &ref_cell ), "with" );
    REIS( node_per, ref_cell_node_per( ref_cell ), "match" );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  { /* unique nodes of one tri */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT nnode, nface, *g2l, *l2g;

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");
    ref_cell = ref_grid_tri(ref_grid);

    nodes[0] = 5; nodes[1] = 8; nodes[2] = 6; nodes[3] = 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_grid_boundary_nodes(ref_grid,10,&nnode,&nface,&g2l,&l2g),"no list");
    REIS(1,nface, "mis count");
    REIS(3,nnode, "mis count");
    REIS(5,l2g[0], "not in list");
    REIS(6,l2g[1], "not in list");
    REIS(8,l2g[2], "not in list");

    ref_free( l2g );
    ref_free( g2l );

    RSS(ref_grid_free(ref_grid),"cleanup");
  }

  { /* unique nodes of two tri */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT nnode, nface, *g2l, *l2g;

    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");
    ref_cell = ref_grid_tri(ref_grid);

    nodes[0] = 5; nodes[1] = 8; nodes[2] = 6; nodes[3] = 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    nodes[0] = 6; nodes[1] = 8; nodes[2] = 9; nodes[3] = 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_grid_boundary_nodes(ref_grid,10,&nnode,&nface,&g2l,&l2g),"no list");
    REIS(2,nface, "mis count");
    REIS(4,nnode, "mis count");
    REIS(5,l2g[0], "not in list");
    REIS(6,l2g[1], "not in list");
    REIS(8,l2g[2], "not in list");
    REIS(9,l2g[3], "not in list");

    ref_free( l2g );
    ref_free( g2l );

    RSS(ref_grid_free(ref_grid),"cleanup");
  }

  { /* orient outward */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add inward tri");

    RSS( ref_grid_outward_boundary_orientation( ref_grid ), "get out" );

    RSS(ref_cell_nodes(ref_grid_tri(ref_grid),cell,nodes),"add tet");
    REIS( 1, nodes[0], "n0" );
    REIS( 0, nodes[1], "n1" );
    REIS( 2, nodes[2], "n2" );
    REIS(10, nodes[3], "id" );
    
    RSS(ref_grid_free(ref_grid),"cleanup");
  }
  
  { /* orient inward */
    REF_GRID ref_grid;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    RSS(ref_grid_create(&ref_grid,ref_mpi),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");

    nodes[0] = 1; nodes[1] = 0; nodes[2] = 2; nodes[3] = 10;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add inward tri");

    RSS( ref_grid_inward_boundary_orientation( ref_grid ), "get out" );

    RSS(ref_cell_nodes(ref_grid_tri(ref_grid),cell,nodes),"add tet");
    REIS( 0, nodes[0], "n0" );
    REIS( 1, nodes[1], "n1" );
    REIS( 2, nodes[2], "n2" );
    REIS(10, nodes[3], "id" );
    
    RSS(ref_grid_free(ref_grid),"cleanup");
  }
  
  { /* single tri enclosing */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[3];
    REF_INT tri;
    RSS( ref_fixture_pri_grid( &ref_grid, ref_mpi ), "fix" );

    xyz[0]= 0.2;
    xyz[1]= 0.0;
    xyz[2]= 0.3;
    tri = 0;
    RSS( ref_grid_enclosing_tri( ref_grid, xyz,
				 &tri, bary ), "enclose");

    REIS( 0, tri, "tri" );
    RWDS( 0.5, bary[0], -1, "b0" );
    RWDS( 0.2, bary[1], -1, "b1" );
    RWDS( 0.3, bary[2], -1, "b2" );

    xyz[0]= 0.2;
    xyz[1]= 0.0;
    xyz[2]= 0.3;
    tri = REF_EMPTY;
    RSS( ref_grid_enclosing_tri( ref_grid, xyz,
				 &tri, bary ), "enclose");

    REIS( 0, tri, "tri" );
    RWDS( 0.5, bary[0], -1, "b0" );
    RWDS( 0.2, bary[1], -1, "b1" );
    RWDS( 0.3, bary[2], -1, "b2" );

    RSS(ref_grid_free(ref_grid),"cleanup");
  }
  
  { /* walk to find enclosing tri */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[3];
    REF_INT tri;
    RSS( ref_fixture_twod_brick_grid( &ref_grid, ref_mpi ), "fix" );

    xyz[0]= 0.5;
    xyz[1]= 0.0;
    xyz[2]= 0.5;
    tri = REF_EMPTY;
    RSS( ref_grid_enclosing_tri( ref_grid, xyz,
				 &tri, bary ), "enclose");

    REIS( 8, tri, "tri" );
    RWDS( 0.0, bary[0], -1, "b0" );
    RWDS( 0.5, bary[1], -1, "b1" );
    RWDS( 0.5, bary[2], -1, "b2" );
    
    RSS(ref_grid_free(ref_grid),"cleanup");
  }

  { /* single tet enclosing */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS( ref_fixture_tet_grid( &ref_grid, ref_mpi ), "fix" );

    xyz[0]= 0.2;
    xyz[1]= 0.1;
    xyz[2]= 0.3;
    tet = 0;
    RSS( ref_grid_enclosing_tet( ref_grid, xyz,
				 &tet, bary ), "enclose");

    REIS( 0, tet, "tet" );
    RWDS( 0.4, bary[0], -1, "b0" );
    RWDS( 0.2, bary[1], -1, "b1" );
    RWDS( 0.1, bary[2], -1, "b2" );
    RWDS( 0.3, bary[3], -1, "b3" );

    xyz[0]= 0.2;
    xyz[1]= 0.1;
    xyz[2]= 0.3;
    tet = REF_EMPTY;
    RSS( ref_grid_enclosing_tet( ref_grid, xyz,
				 &tet, bary ), "enclose");

    REIS( 0, tet, "tet" );
    RWDS( 0.4, bary[0], -1, "b0" );
    RWDS( 0.2, bary[1], -1, "b1" );
    RWDS( 0.1, bary[2], -1, "b2" );
    RWDS( 0.3, bary[3], -1, "b3" );

    RSS(ref_grid_free(ref_grid),"cleanup");
  }
  
  { /* walk to find enclosing tet */
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "fix" );

    xyz[0]= 0.5;
    xyz[1]= 0.5;
    xyz[2]= 0.5;
    tet = REF_EMPTY;
    RSS( ref_grid_enclosing_tet( ref_grid, xyz,
				 &tet, bary ), "enclose");

    REIS( 79, tet, "tet" );
    RWDS( 0.0, bary[0], -1, "b0" );
    RWDS( 0.5, bary[1], -1, "b1" );
    RWDS( 0.0, bary[2], -1, "b2" );
    RWDS( 0.5, bary[3], -1, "b3" );
    
    RSS(ref_grid_free(ref_grid),"cleanup");
  }

  { /* walk to find enclosing tet falls outside*/
    REF_GRID ref_grid;
    REF_DBL xyz[3], bary[4];
    REF_INT tet;
    RSS( ref_fixture_tet_brick_grid( &ref_grid, ref_mpi ), "fix" );

    xyz[0]= 0.5;
    xyz[1]= 0.5;
    xyz[2]= -0.01;
    tet = REF_EMPTY;
    RSS( ref_grid_enclosing_tet( ref_grid, xyz,
				 &tet, bary ), "enclose");

    REIS( 25, tet, "tet" );
    RWDS( 0.53, bary[0], -1, "b0" );
    RWDS( 0.50, bary[1], -1, "b1" );
    RWDS( 0.00, bary[2], -1, "b2" );
    RWDS(-0.03, bary[3], -1, "b3" );
    
    RSS(ref_grid_free(ref_grid),"cleanup");
  }

  RSS( ref_mpi_free( ref_mpi ), "free" );
  return 0;
}
