
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_metric.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"

#include "ref_malloc.h"
#include "ref_matrix.h"

REF_STATUS ref_metric_imply_from( REF_DBL *metric, REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT *hits;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  REF_CELL ref_cell = ref_grid_tet(ref_grid);

  ref_malloc_init( hits, ref_node_max(ref_node),
		   REF_INT, 0 );

  for( node=0; node<ref_node_max(ref_node); node++ )
    for( im=0; im<6; im++ )
      metric[im+6*node] = 0.0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,nodes[0]), 
			       ref_node_xyz_ptr(ref_node,nodes[1]), 
			       ref_node_xyz_ptr(ref_node,nodes[2]), 
			       ref_node_xyz_ptr(ref_node,nodes[3])), "imply" );
      RSS( ref_matrix_log_m( m, log_m ), "log" );

      for( node=0; node<4; node++ )
	{
	  hits[nodes[node]]++;
	  for( im=0; im<6; im++ )
	    metric[im+6*nodes[node]] += log_m[im];
	}
    }

  each_ref_node_valid_node( ref_node, node )
    {
      RAS( 0 != hits[node], "zero metric contributions" );
      for( im=0; im<6; im++ )
	log_m[im] = metric[im+6*node] / ((REF_DBL)hits[node]);
      RSS( ref_matrix_exp_m( log_m, m ), "exp" );
      for( im=0; im<6; im++ )
	metric[im+6*node] = m[im];
    }

  ref_free( hits );

  return REF_SUCCESS;
}
