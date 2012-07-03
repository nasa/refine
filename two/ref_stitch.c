
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_stitch.h"

#include "ref_math.h"

#include "ref_node.h"
#include "ref_cell.h"

#include "ref_sort.h"
#include "ref_malloc.h"

REF_STATUS ref_stitch_together( REF_GRID ref_grid, 
				REF_INT tri_boundary, REF_INT qua_boundary )
{
  REF_NODE ref_node;
  REF_INT tri_nnode, tri_nface, *tri_g2l, *tri_l2g;
  REF_INT qua_nnode, qua_nface, *qua_g2l, *qua_l2g;
  REF_INT tri_node, qua_node;
  REF_DBL d, dist2, tol, tol2;
  REF_INT *t2q;
  REF_CELL ref_cell, hex, tri;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT hex_cell, hex_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_node, global;
  REF_INT node;
  REF_INT cell_face, base[4], tri_nodes[4], tri_cell;
  REF_INT top_face;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_grid_boundary_nodes(ref_grid,tri_boundary,
			      &tri_nnode,&tri_nface,&tri_g2l,&tri_l2g),"l2g");

  printf("tri bound %d : %d nnode %d nface\n",
	 tri_boundary,tri_nnode,tri_nface);

  RSS(ref_grid_boundary_nodes(ref_grid,qua_boundary,
			      &qua_nnode,&qua_nface,&qua_g2l,&qua_l2g),"l2g");

  printf("qua bound %d : %d nnode %d nface\n",
	 qua_boundary,qua_nnode,qua_nface);

  if ( tri_nnode != qua_nnode ||
       tri_nface != 2*qua_nface )
    THROW("face sizes do not match. stop.");

  ref_malloc_init( t2q, tri_nnode, REF_INT, REF_EMPTY );

  tol = 1.0e-8;
  printf("same point tol %e\n",tol);
  tol2 = tol*tol;

  printf("start %d sized n-squared search\n",tri_nnode);

  for( tri_node = 0 ; tri_node < tri_nnode ; tri_node++ )
    {
      for( qua_node = 0 ; qua_node < qua_nnode ; qua_node++ )
	{
	  dist2 = 0.0;
	  d = ref_node_xyz(ref_node,0,tri_l2g[tri_node]) - 
	    ref_node_xyz(ref_node,0,qua_l2g[qua_node]);
	  dist2 += d*d;
	  d = ref_node_xyz(ref_node,1,tri_l2g[tri_node]) - 
	    ref_node_xyz(ref_node,1,qua_l2g[qua_node]);
	  dist2 += d*d;
	  d = ref_node_xyz(ref_node,2,tri_l2g[tri_node]) - 
	    ref_node_xyz(ref_node,2,qua_l2g[qua_node]);
	  dist2 += d*d;
	  if ( dist2 < tol2 ) {
	    t2q[tri_node] = qua_l2g[qua_node];
	    break;
	  }
	}
      if ( REF_EMPTY == t2q[tri_node] )
	THROW("point not matched. stop.");
    }

  printf("collapsing duplicate nodes.\n");

  for( tri_node = 0 ; tri_node < tri_nnode ; tri_node++ )
    {
      RSS( ref_grid_replace_node( ref_grid, tri_l2g[tri_node], t2q[tri_node] ), 
	   "repl");
    }

  printf("removing duplicate nodes.\n");
  
  for( tri_node = 0 ; tri_node < tri_nnode ; tri_node++ )
    RSS( ref_node_remove_requiring_rebuild( ref_node, tri_l2g[tri_node] ), 
	 "rm node");

  printf("rebuild sorted globals.\n");
  RSS( ref_node_rebuild_sorted_global( ref_node ), "rebuild");

  printf("removing iterior quads and splitting near hexes.\n");

  hex = ref_grid_hex(ref_grid);
  tri = ref_grid_tri(ref_grid);

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( qua_boundary == nodes[4] )
      {
	RSS( ref_cell_with_face( hex, nodes, &hex_cell ), "missing hex" );
	RSS( ref_cell_nodes( hex, hex_cell, hex_nodes),"hex nodes");

	RSS( ref_node_next_global(ref_node,&global), "get glob");
	RSS( ref_node_add(ref_node,global,&new_node), "add node");

	ref_node_xyz(ref_node,0,new_node) = 0.0;
	ref_node_xyz(ref_node,1,new_node) = 0.0;
	ref_node_xyz(ref_node,2,new_node) = 0.0;
	for( node = 0 ; node < 8 ; node++ )
	  {
	    ref_node_xyz(ref_node,0,new_node) += 
	      0.125 * ref_node_xyz(ref_node,0,hex_nodes[node]);
	    ref_node_xyz(ref_node,1,new_node) += 
	      0.125 * ref_node_xyz(ref_node,1,hex_nodes[node]);
	    ref_node_xyz(ref_node,2,new_node) += 
	      0.125 * ref_node_xyz(ref_node,2,hex_nodes[node]);
	  }

	top_face = REF_EMPTY;
	for (cell_face=0;cell_face<ref_cell_face_per(hex);cell_face++)
	  {
	    for (node=0;node<4; node++)
	      base[node] = ref_cell_f2n(hex,node,cell_face,hex_cell);

	    for (node=0;node<3; node++)
	      tri_nodes[node] = base[node];
	    tri_nodes[3]=tri_nodes[0];

	    RXS( ref_cell_with_face( tri, tri_nodes, &tri_cell),
		 REF_NOT_FOUND,
		 "cell has face" );

	    if ( REF_EMPTY != tri_cell )
	      printf("hex_cell %d face %d\n",hex_cell,cell_face);

	    for (node=0;node<3; node++)
	      tri_nodes[node] = base[node+1];
	    tri_nodes[3]=tri_nodes[0];

	    RXS( ref_cell_with_face( tri, tri_nodes, &tri_cell),
		 REF_NOT_FOUND,
		 "cell has face" );

	    if ( REF_EMPTY != tri_cell )
	      printf("hex_cell %d face %d\n",hex_cell,cell_face);

	    
	    /* see if there is a triagle -> make tets */
	    /* no tri, make pyrimid */

	  }

	if ( REF_EMPTY == top_face ) THROW( "top tris not found" );

	RSS( ref_cell_remove( ref_cell, cell ), "rm qua" );
      }

  printf("removing iterior triangles.\n");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( tri_boundary == nodes[3] )
      RSS( ref_cell_remove( ref_cell, cell ), "rm tri" );

  ref_free( t2q );
  ref_free( qua_l2g );
  ref_free( qua_g2l );
  ref_free( tri_l2g );
  ref_free( tri_g2l );

  return REF_SUCCESS;
}

