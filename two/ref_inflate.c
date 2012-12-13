
#include <stdlib.h>
#include <stdio.h>

#include "ref_inflate.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_malloc.h"
#include "ref_math.h"

REF_STATUS ref_inflate_normal( REF_GRID ref_grid, 
			       REF_DICT faceids, 
			       REF_DBL thickness, REF_DBL xshift )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL projection[2];
  REF_DBL normal[3];
  REF_DBL rx[2];
  REF_DBL len, dot;

  projection[0] = xshift;
  projection[1] = thickness;
  len = sqrt( projection[0]*projection[0] + projection[1]*projection[1] );
  if ( !ref_math_divisible(projection[0],len) ||
       !ref_math_divisible(projection[1],len))
    THROW("div zero");
  projection[0] /= len;
  projection[1] /= len;

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	RSS( ref_node_tri_normal( ref_node, nodes, normal ), "norm" );
	rx[0] = normal[0];
	rx[1] = sqrt( normal[1]*normal[1] + normal[2]*normal[2] );

	len = sqrt( rx[0]*rx[0] + rx[1]*rx[1] );
	if ( !ref_math_divisible(rx[0],len) ||
	     !ref_math_divisible(rx[1],len))
	  THROW("div zero");
	rx[0] /= len;
	rx[1] /= len;

	dot = rx[0]*projection[0] + rx[1]*projection[1];
	if ( dot <= 0.0 )
	  {
	    printf("neg pri\n");
	  }
	if ( dot <= 0.1 )
	  {
	    printf("proj [%f,%f] rx [%f,%f] dot %f\n",
		   projection[0],projection[1],rx[0],rx[1],dot);
	  }
      }

  return REF_SUCCESS;
}

REF_STATUS ref_inflate_face( REF_GRID ref_grid, 
			     REF_DICT faceids, 
			     REF_DBL thickness, REF_DBL xshift )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL qua = ref_grid_qua(ref_grid);
  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_INT cell, tri_side, node0, node1;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ntri, tris[2], nquad, quads[2];
  REF_INT tri_node;
  REF_INT *o2n;
  REF_INT global, new_node;
  REF_DBL radius, scale;
  REF_INT new_cell;
  REF_DBL min_dot;

  ref_malloc_init( o2n, ref_node_max(ref_node), 
		   REF_INT, REF_EMPTY );

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
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
		xshift + ref_node_xyz(ref_node,0,node0);
	      ref_node_xyz(ref_node,1,new_node) = 
		scale * ref_node_xyz(ref_node,1,node0);
	      ref_node_xyz(ref_node,2,new_node) = 
		scale * ref_node_xyz(ref_node,2,node0);
	    }
	}

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	for(tri_side=0;tri_side<3;tri_side++)
	  {
	    node0 = ref_cell_e2n(tri,0,tri_side,cell);
	    node1 = ref_cell_e2n(tri,1,tri_side,cell);
	    RSS( ref_cell_list_with( tri, 
				     node0, node1,
				     2, &ntri,
				     tris ),"bad tri count");
	    if ( 1 == ntri ) 
	      {
		RSS( ref_cell_list_with( qua, 
					 node0, node1,
					 2, &nquad,
					 quads ),"bad quad count");
		if ( 1 != nquad ) THROW("tri without quad");
		new_nodes[4] = ref_cell_c2n(qua,4,quads[0]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	    if (  ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[0]) ) &&
		 !ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[1]) ) )
	      {
		new_nodes[4] = ref_cell_c2n(tri,3,tris[1]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	    if ( !ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[0]) ) &&
		  ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[1]) ) )
	      {
		new_nodes[4] = ref_cell_c2n(tri,3,tris[0]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	  }
      }

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	new_nodes[0] = nodes[0];
	new_nodes[1] = nodes[2];
	new_nodes[2] = nodes[1];
	new_nodes[3] = o2n[nodes[0]];
	new_nodes[4] = o2n[nodes[2]];
	new_nodes[5] = o2n[nodes[1]];
	
	RSS( ref_node_pri_min_dot( ref_node, new_nodes, &min_dot ), "md");
	if ( min_dot <= 0.0 ) printf("min_dot %f\n",min_dot);
	
	RSS( ref_cell_add( pri, new_nodes, &new_cell ), "pri");
      }

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	nodes[0] = o2n[nodes[0]];
	nodes[1] = o2n[nodes[1]];
	nodes[2] = o2n[nodes[2]];
	RSS( ref_cell_replace_whole( tri, cell, nodes ), "repl");
      }

  ref_free( o2n );

  return REF_SUCCESS;
}

