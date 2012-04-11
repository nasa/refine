
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_metric.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"

#include "ref_malloc.h"
#include "ref_matrix.h"
#include "ref_math.h"

REF_STATUS ref_metric_show( REF_DBL *m )
{
  printf(" %18.10e %18.10e %18.10e\n",m[0],m[1],m[2]);
  printf(" %18.10e %18.10e %18.10e\n",m[1],m[3],m[4]);
  printf(" %18.10e %18.10e %18.10e\n",m[2],m[4],m[5]);
  return REF_SUCCESS;
}

REF_STATUS ref_metric_from_node( REF_DBL *metric, REF_NODE ref_node )
{
  REF_INT node, im;

  each_ref_node_valid_node( ref_node, node )
    for(im=0;im<6;im++)
      metric[im+6*node] = 
	ref_node_metric(ref_node,im,node);

  return REF_SUCCESS;
}

REF_STATUS ref_metric_to_node( REF_DBL *metric, REF_NODE ref_node )
{
  REF_INT node, im;

  each_ref_node_valid_node( ref_node, node )
    for(im=0;im<6;im++)
      ref_node_metric(ref_node,im,node) = 
	metric[im+6*node];

  return REF_SUCCESS;
}

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

REF_STATUS ref_metric_smr( REF_DBL *metric0, REF_DBL *metric1, REF_DBL *metric,
			   REF_GRID ref_grid )
{
  REF_INT node;
  REF_DBL metric_inv[6];
  REF_DBL inv_m1_m2[9];
  REF_DBL n_values[3], n_vectors[9], inv_n_vectors[9];
  REF_DBL diagonal_system[12];
  REF_DBL h0, h1, h, hmax, hmin, h2;
  REF_DBL eig;
  REF_INT i;

  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      RSS( ref_matrix_inv_m( &(metric0[6*node]), metric_inv), "inv" );
      RSS( ref_matrix_mult_m( metric_inv, &(metric1[6*node]), inv_m1_m2 ), 
	   "mult" );
      RSS( ref_matrix_diag_gen( 3, inv_m1_m2, n_values, n_vectors ), "gen eig");
      for (i=0;i<3;i++)
	{
	  h0 = ref_matrix_sqrt_vt_m_v( &(metric0[6*node]), &(n_vectors[i*3]) );
	  if ( !ref_math_divisible( 1.0, h0 ) ) RSS( REF_DIV_ZERO, "inf h0");
	  h0 = 1.0/h0;
	  h1 = ref_matrix_sqrt_vt_m_v( &(metric1[6*node]), &(n_vectors[i*3]) );
	  if ( !ref_math_divisible( 1.0, h1 ) ) RSS( REF_DIV_ZERO, "inf h1");
	  h1 = 1.0/h1;
	  hmax = 4.00*h0; 
	  hmin = 0.25*h0; 
	  h = MIN( hmax, MAX( hmin, h1 ));
	  h2 = h*h;
	  if ( !ref_math_divisible( 1.0, h2 ) ) RSS( REF_DIV_ZERO, "zero h^2");
	  eig = 1.0/h2;
	  ref_matrix_eig( diagonal_system, i ) = eig;
	}
      RSS( ref_matrix_inv_gen( 3, n_vectors, inv_n_vectors ), "gen eig");
      RSS( ref_matrix_transpose_gen( 3, inv_n_vectors, 
				     &(diagonal_system[3]) ), "gen eig");
      RSS( ref_matrix_form_m( diagonal_system, &(metric[6*node]) ), "reform m");
    }

  return REF_SUCCESS;
}

