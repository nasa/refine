
#include <stdlib.h>
#include <stdio.h>

#include "ref_inflate.h"

#include "ref_cell.h"

REF_STATUS ref_inflate_face( REF_GRID ref_grid, REF_INT faceid )
{
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT cell, tri_side, node0, node1;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ntri, tris[2];

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( faceid == nodes[3] )
      {
	for(tri_side=0;tri_side<3;tri_side++)
	  {
	    node0 = ref_cell_e2n(tri,0,tri_side,cell);
	    node1 = ref_cell_e2n(tri,1,tri_side,cell);
	    RSS( ref_cell_list_with( tri, 
				     node0, node1,
				     2, &ntri,
				     tris ),"bad tri count");
	    if ( 2 != ntri ) THROW("not manifold");
	    if ( ( ref_cell_c2n(tri,3,tris[0]) == faceid &&
		   ref_cell_c2n(tri,3,tris[1]) != faceid  ) || 
		 ( ref_cell_c2n(tri,3,tris[0]) != faceid &&
		   ref_cell_c2n(tri,3,tris[1]) == faceid )  )
	      {
		printf("side\n");
	      }
	  }
      }

  return REF_SUCCESS;
}

