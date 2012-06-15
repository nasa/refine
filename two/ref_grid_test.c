#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"

#include "ref_list.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_mpi.h"
#include "ref_malloc.h"

int main( void )
{

  {  /* init */
    REF_GRID ref_grid;
    REIS(REF_NULL,ref_grid_free(NULL),"dont free NULL");

    RSS(ref_grid_create(&ref_grid),"create");
    RSS(ref_grid_free(ref_grid),"free");
  }

  {  /* each element */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;
    REF_INT group;

    RSS(ref_grid_create(&ref_grid),"create");

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

    RSS(ref_grid_create(&ref_grid),"create");

    node_per = 4;
    RSS( ref_grid_cell_with( ref_grid, node_per, &ref_cell ), "with" );
    REIS( node_per, ref_cell_node_per( ref_cell ), "match" );

    RSS(ref_grid_free(ref_grid),"free"); 
  }

  {  /* face with these many nodes */
    REF_GRID ref_grid;
    REF_CELL ref_cell;
    REF_INT node_per;

    RSS(ref_grid_create(&ref_grid),"create");

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

    RSS(ref_grid_create(&ref_grid),"create");
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

    RSS(ref_grid_create(&ref_grid),"create");
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

  return 0;
}
