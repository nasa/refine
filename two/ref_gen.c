
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_gen.h"

#include "ref_node.h"
#include "ref_cell.h"

#include "ref_math.h"

/*
                                                     5
                                                  . /|
                                               .   / |
                                            .     /  |
                                         .       /   |
                                      .         /    |
                                   .           /     |
                                .             /      |
                             .               /       |
                          .                 /        |
                         3--q3---------q2--4         2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                   y     |               . |     /
                   ^     |            .    |    /
                   |     |         .       |   /
                   |     |      .          |  /
                         |   .             | /
                         |.                |/
                         0--q0---------q1--1

 */

/*
                                                     5
                                                  . /|
                                               .   / |
                                            .     /  |
                                         .       /   |
                                      .         /    |
                                   .           /     |
                                .             /      |
                             .               /       |
                          .                 /        |
                         3-----------------4         2
                         |                 |      . /
                         |                 |   .   /
                         |                 |.     /
                         |               . |     /
                         |            .    |    /
                         |         .       |   /
                         |      .          |  /
                         |   .             | /
                         |.                |/
                         0-----------------1

 */

REF_STATUS ref_gen_fill( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL front;
  REF_CELL ref_cell;
  REF_INT cell, new_cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pri[REF_CELL_MAX_SIZE_PER];
  REF_INT tri[REF_CELL_MAX_SIZE_PER];
  REF_INT qua[REF_CELL_MAX_SIZE_PER];
  REF_INT global;

  RSS( ref_cell_create( &front, 4, REF_TRUE ), "qua front" );

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    RSS( ref_cell_add( front, nodes, &new_cell ), "front add" );

  each_ref_cell_valid_cell_with_nodes( front, cell, nodes )
    {
      REF_DBL h = 0.4;
      REF_DBL norm[3];

      norm[0] = -(ref_node_xyz(ref_node,2,nodes[1]) -
		  ref_node_xyz(ref_node,2,nodes[0]));
      norm[1] = 0.0;
      norm[2] =  (ref_node_xyz(ref_node,0,nodes[1]) -
		  ref_node_xyz(ref_node,0,nodes[0]));
      RSS(ref_math_normalize(norm),"norm");
      norm[0] *= h; norm[1] *= h; norm[2] *= h;

      pri[0] = nodes[0];
      pri[1] = nodes[1];
      RSS( ref_node_next_global( ref_node, &global ), "glob");
      RSS( ref_node_add( ref_node, global, &(pri[2]) ), "nn");
      ref_node_xyz(ref_node,0,pri[2]) = norm[0] + 
	0.5 * ( ref_node_xyz(ref_node,0,pri[0]) +
		ref_node_xyz(ref_node,0,pri[1]) ); 
      ref_node_xyz(ref_node,1,pri[2]) = norm[1] + 
	0.5 * ( ref_node_xyz(ref_node,1,pri[0]) +
		ref_node_xyz(ref_node,1,pri[1]) ); 
      ref_node_xyz(ref_node,2,pri[2]) = norm[2] + 
	0.5 * ( ref_node_xyz(ref_node,2,pri[0]) +
		ref_node_xyz(ref_node,2,pri[1]) ); 

      tri[0] = pri[0];
      tri[1] = pri[1];
      tri[2] = pri[2];
      tri[3] = 2;
      RSS( ref_cell_add( ref_grid_tri(ref_grid), tri, &new_cell ), "tri add" );

      pri[4] = nodes[2];
      pri[3] = nodes[3];
      RSS( ref_node_next_global( ref_node, &global ), "glob");
      RSS( ref_node_add( ref_node, global, &(pri[5]) ), "nn");
      ref_node_xyz(ref_node,0,pri[5]) = norm[0] + 
	0.5 * ( ref_node_xyz(ref_node,0,pri[0]) +
		ref_node_xyz(ref_node,0,pri[1]) ); 
      ref_node_xyz(ref_node,1,pri[5]) = norm[1] + 
	0.5 * ( ref_node_xyz(ref_node,1,pri[0]) +
		ref_node_xyz(ref_node,1,pri[1]) ); 
      ref_node_xyz(ref_node,2,pri[5]) = norm[2] + 
	0.5 * ( ref_node_xyz(ref_node,2,pri[0]) +
		ref_node_xyz(ref_node,2,pri[1]) ); 

      tri[0] = pri[4];
      tri[1] = pri[3];
      tri[2] = pri[5];
      tri[3] = 1;
      RSS( ref_cell_add( ref_grid_tri(ref_grid), tri, &new_cell ), "tri add" );

      RSS( ref_cell_add( ref_grid_pri(ref_grid), pri, &new_cell ), "pri add" );

      qua[0] = pri[2];
      qua[1] = pri[1];

      RSS(ref_cell_remove(front,cell),"rm");
    }

  RSS( ref_cell_free( front ), "free front" );
  
  return REF_SUCCESS;
}

