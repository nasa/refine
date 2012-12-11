
#include <stdlib.h>
#include <stdio.h>

#include "ref_inflate.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_malloc.h"
#include "ref_math.h"

REF_STATUS ref_inflate_face( REF_GRID ref_grid, REF_INT faceid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT cell, tri_side, node0, node1;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ntri, tris[2];
  REF_INT tri_node;
  REF_INT *o2n;
  REF_INT global, new_node;
  REF_DBL thickness, radius, scale;
  ref_malloc_init( o2n, ref_node_max(ref_node), 
		   REF_INT, REF_EMPTY );

  thickness = 0.025;

  printf("orig nnode %d\n",ref_node_n(ref_node));
  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( faceid == nodes[3] )
      for(tri_node=0;tri_node<3;tri_node++)
	{
	  node0 = nodes[tri_node];
	  if ( REF_EMPTY == o2n[node0] )
	    {
	      RSS( ref_node_next_global( ref_node, &global ), "global" );
	      RSS( ref_node_add( ref_node, global, &new_node ), "add node" );
	      o2n[node0] = new_node;
	      radius = sqrt( ref_node_xyz(ref_node,1,node0) *
			     ref_node_xyz(ref_node,1,node0) +
			     ref_node_xyz(ref_node,2,node0) *
			     ref_node_xyz(ref_node,2,node0) );
	      if ( !ref_math_divisible(radius+thickness,radius) )
		THROW("div zero");
	      scale = (radius+thickness)/radius;
	      ref_node_xyz(ref_node,0,new_node) = 
		ref_node_xyz(ref_node,0,node0);
	      ref_node_xyz(ref_node,1,new_node) = 
		scale*ref_node_xyz(ref_node,1,node0);
	      ref_node_xyz(ref_node,2,new_node) = 
		scale*ref_node_xyz(ref_node,2,node0);
	    }
	}
  printf("new  nnode %d\n",ref_node_n(ref_node));	

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
		
	      }
	  }
      }

  ref_free( o2n );

  return REF_SUCCESS;
}

