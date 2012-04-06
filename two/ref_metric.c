
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
  REF_DBL *total_node_volume, tet_volume;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], pri_nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init( total_node_volume, ref_node_max(ref_node),
		   REF_DBL, 0.0 );

  for( node=0; node<ref_node_max(ref_node); node++ )
    for( im=0; im<6; im++ )
      metric[im+6*node] = 0.0;

  each_ref_cell_valid_cell_with_nodes( ref_grid_tet(ref_grid), cell, tet_nodes)
    {
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<4; node++ )
	{
	  total_node_volume[tet_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*tet_nodes[node]] += tet_volume*log_m[im];
	}
    }

  each_ref_cell_valid_cell_with_nodes( ref_grid_pri(ref_grid), cell, pri_nodes)
    {
      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[4];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[3];
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[4];
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[2];
      tet_nodes[3] = pri_nodes[5];
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]]++;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += log_m[im];
	}
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

    }

  each_ref_node_valid_node( ref_node, node )
    {
      RAS( 0.0 < total_node_volume[node], "zero metric contributions" );
      for( im=0; im<6; im++ )
	{
	  log_m[im] = metric[im+6*node] / total_node_volume[node];
	}
      RSS( ref_matrix_exp_m( log_m, m ), "exp" );
      for( im=0; im<6; im++ )
	metric[im+6*node] = m[im];
    }

  ref_free( total_node_volume );

  return REF_SUCCESS;
}
REF_STATUS ref_metric_imply_non_tet( REF_DBL *metric, REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *total_node_volume, tet_volume;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], pri_nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init( total_node_volume, ref_node_max(ref_node),
		   REF_DBL, 0.0 );

  each_ref_cell_valid_cell_with_nodes( ref_grid_pri(ref_grid), cell, pri_nodes)
    for( node=0; node<6; node++ )
      for( im=0; im<6; im++ )
	metric[im+6*pri_nodes[node]] = 0.0;

  each_ref_cell_valid_cell_with_nodes( ref_grid_pri(ref_grid), cell, pri_nodes)
    {
      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[4];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[3];
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[4];
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[1];
      tet_nodes[2] = pri_nodes[2];
      tet_nodes[3] = pri_nodes[5];
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]]++;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += log_m[im];
	}
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,tet_nodes[0]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[1]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[2]), 
			       ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl");
      RSS( ref_matrix_log_m( m, log_m ), "log" );
      RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );
      for( node=0; node<6; node++ )
	{
	  total_node_volume[pri_nodes[node]] += tet_volume;
	  for( im=0; im<6; im++ )
	    metric[im+6*pri_nodes[node]] += tet_volume*log_m[im];
	}

    }

  each_ref_cell_valid_cell_with_nodes( ref_grid_pri(ref_grid), cell, pri_nodes)
    for( node=0; node<6; node++ )
      if ( 0.0 < total_node_volume[pri_nodes[node]] )
	{
	  for( im=0; im<6; im++ )
	    {
	      log_m[im] = 
		metric[im+6*node] / total_node_volume[pri_nodes[node]];
	    }
	  RSS( ref_matrix_exp_m( log_m, m ), "exp" );
	  for( im=0; im<6; im++ )
	    metric[im+6*node] = m[im];
	  total_node_volume[pri_nodes[node]] = 0.0;
	}

  ref_free( total_node_volume );

  return REF_SUCCESS;
}
