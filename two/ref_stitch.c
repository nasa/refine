
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
  REF_CELL ref_cell;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

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
	    t2q[tri_node] = qua_node;
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

  printf("removing iterior boundary faces.\n");

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( tri_boundary == nodes[3] )
      RSS( ref_cell_remove( ref_cell, cell ), "rm tri" );

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    if ( qua_boundary == nodes[4] )
      RSS( ref_cell_remove( ref_cell, cell ), "rm qua" );

  ref_free( t2q );
  ref_free( qua_l2g );
  ref_free( qua_g2l );
  ref_free( tri_l2g );
  ref_free( tri_g2l );

  return REF_SUCCESS;
}

