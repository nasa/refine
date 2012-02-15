
#include <stdlib.h>
#include <stdio.h>

#include "ref_axi.h"

#include "ref_node.h"

REF_STATUS ref_axi_wedge( REF_GRID ref_grid )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_DBL pole_tol;
  REF_INT *o2n, node;

  ref_node = ref_grid_node(ref_grid);

  pole_tol = 1.0e-7;
  o2n = (REF_INT *)malloc( ref_node_n(ref_node) * sizeof(REF_INT));
  RNS(o2n,"malloc o2n NULL");

  each_ref_node_valid_node( ref_node, node )
    {
      if ( ABS(ref_node_xyz(ref_node,2,node)) < pole_tol &&
	   ABS(ref_node_xyz(ref_node,1,node)) > 0.5 )
	{
	  o2n[node] = node-ref_node_n(ref_node)/2;
	  RSS( ref_node_remove( ref_node, node ), "remove" );
	}
      else
	{
	  o2n[node] = node;
	}
    }

  /* FIXME, hack until globals are handled right */
  ref_node_n_global(ref_node) = ref_node_n(ref_node);

  ref_cell = ref_grid_qua(ref_grid);
  ref_cell_n(ref_cell) = 0;

  free(o2n);

  return REF_SUCCESS;
}

