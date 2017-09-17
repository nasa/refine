
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_metric.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_cell.h"
#include "ref_edge.h"

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

REF_STATUS ref_metric_inspect( REF_NODE ref_node )
{
  REF_INT node;
  each_ref_node_valid_node( ref_node, node )
    RSS( ref_metric_show( ref_node_metric_ptr(ref_node,node) ), "show it" );
  
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

REF_STATUS ref_metric_unit_node( REF_NODE ref_node )
{
  REF_INT node;

  each_ref_node_valid_node( ref_node, node )
    {
      ref_node_metric(ref_node,0,node) = 1.0;
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,2,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0;
      ref_node_metric(ref_node,4,node) = 0.0;
      ref_node_metric(ref_node,5,node) = 1.0;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_olympic_node( REF_NODE ref_node, REF_DBL h )
{
  REF_INT node;
  REF_DBL hh;

  each_ref_node_valid_node( ref_node, node )
    {
      ref_node_metric(ref_node,0,node) = 1.0/(0.1*0.1);
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,2,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0/(0.1*0.1);
      ref_node_metric(ref_node,4,node) = 0.0;
      hh = h + (0.1-h)*ABS(ref_node_xyz(ref_node,2,node)-0.5)/0.5;
      ref_node_metric(ref_node,5,node) = 1.0/(hh*hh);
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_ring_node( REF_NODE ref_node )
{
  REF_INT node;
  REF_DBL hh;
  REF_DBL h =0.01;
  REF_DBL x;
  each_ref_node_valid_node( ref_node, node )
    {
      x= ref_node_xyz(ref_node,0,node);
      hh = h + (0.1-h)*MIN(2*ABS(x-1.0),1);
      ref_node_metric(ref_node,0,node) = 1.0/(hh*hh);
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,2,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0/(0.1*0.1);
      ref_node_metric(ref_node,4,node) = 0.0;
      ref_node_metric(ref_node,5,node) = 1.0/(0.1*0.1);
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_polar2d_node( REF_NODE ref_node )
{
  REF_INT node;
  REF_DBL x, z, r, t;
  REF_DBL h_y, h_t, h_r, h0;
  REF_DBL d[12], m[6];
  
  each_ref_node_valid_node( ref_node, node )
    {
      x = ref_node_xyz(ref_node,0,node);
      z = ref_node_xyz(ref_node,2,node);
      r = sqrt(x*x+z*z);
      t = atan2(z,x);
      h_y = 1.0;
      h_t = 0.1;
      h0 = 0.001;
      h_r = h0 + 2*(0.1-h0)*ABS(r-0.5);
      ref_matrix_eig( d, 0 ) = 1.0/(h_r*h_r);
      ref_matrix_eig( d, 1 ) = 1.0/(h_t*h_t);
      ref_matrix_eig( d, 2 ) = 1.0/(h_y*h_y);
      ref_matrix_vec( d, 0, 0 ) = cos( t );
      ref_matrix_vec( d, 1, 0 ) = 0.0;
      ref_matrix_vec( d, 2, 0 ) = sin( t );
      ref_matrix_vec( d, 0, 1 ) =-sin( t );
      ref_matrix_vec( d, 1, 1 ) = 0.0;
      ref_matrix_vec( d, 2, 1 ) = cos( t );
      ref_matrix_vec( d, 0, 2 ) = 0.0;
      ref_matrix_vec( d, 1, 2 ) = 1.0;
      ref_matrix_vec( d, 2, 2 ) = 0.0;
      ref_matrix_form_m( d, m );
      ref_node_metric(ref_node,0,node) = m[0];
      ref_node_metric(ref_node,1,node) = m[1];
      ref_node_metric(ref_node,2,node) = m[2];
      ref_node_metric(ref_node,3,node) = m[3];
      ref_node_metric(ref_node,4,node) = m[4];
      ref_node_metric(ref_node,5,node) = m[5];
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_ugawg_node( REF_NODE ref_node, REF_INT version )
{
  REF_INT node;
  REF_DBL x, y, r, t;
  REF_DBL h_z, h_t, h_r, h0, d0;
  REF_DBL d[12], m[6];
  
  each_ref_node_valid_node( ref_node, node )
    {
      x = ref_node_xyz(ref_node,0,node);
      y = ref_node_xyz(ref_node,1,node);
      r = sqrt(x*x+y*y);
      t = atan2(y,x);
      h_z = 0.1;
      h_t = 0.1;
      h0 = 0.001;
      h_r = h0 + 2*(0.1-h0)*ABS(r-0.5);
      if (2 == version)
	{ 
	  d0 = MIN( 10.0*ABS(r-0.5), 1.0 );
	  h_t = 0.1 * d0 + 0.025 * (1.0-d0);
	}
      ref_matrix_eig( d, 0 ) = 1.0/(h_r*h_r);
      ref_matrix_eig( d, 1 ) = 1.0/(h_t*h_t);
      ref_matrix_eig( d, 2 ) = 1.0/(h_z*h_z);
      ref_matrix_vec( d, 0, 0 ) = cos( t );
      ref_matrix_vec( d, 1, 0 ) = sin( t );
      ref_matrix_vec( d, 2, 0 ) = 0.0;
      ref_matrix_vec( d, 0, 1 ) =-sin( t );
      ref_matrix_vec( d, 1, 1 ) = cos( t );
      ref_matrix_vec( d, 2, 1 ) = 0.0;
      ref_matrix_vec( d, 0, 2 ) = 0.0;
      ref_matrix_vec( d, 1, 2 ) = 0.0;
      ref_matrix_vec( d, 2, 2 ) = 1.0;
      ref_matrix_form_m( d, m );
      ref_node_metric(ref_node,0,node) = m[0];
      ref_node_metric(ref_node,1,node) = m[1];
      ref_node_metric(ref_node,2,node) = m[2];
      ref_node_metric(ref_node,3,node) = m[3];
      ref_node_metric(ref_node,4,node) = m[4];
      ref_node_metric(ref_node,5,node) = m[5];
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_masabl_node( REF_NODE ref_node )
{
  REF_INT node;
  REF_DBL hx, hz, c,k1;

  each_ref_node_valid_node( ref_node, node )
    {
      hx = 0.01+0.2*cos(ref_math_pi*(ref_node_xyz(ref_node,0,node)-0.5));
      ref_node_metric(ref_node,0,node) = 1.0/(hx*hx);
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,2,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0/(0.1*0.1);
      ref_node_metric(ref_node,4,node) = 0.0;
      c = 0.001;
      k1 = 6.0;
      hz = c*exp(k1*ref_node_xyz(ref_node,2,node));
      ref_node_metric(ref_node,5,node) = 1.0/(hz*hz);
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_twod_node( REF_NODE ref_node )
{
  REF_INT node;

  each_ref_node_valid_node( ref_node, node )
    {
      ref_node_metric(ref_node,1,node) = 0.0;
      ref_node_metric(ref_node,3,node) = 1.0;
      ref_node_metric(ref_node,4,node) = 0.0;
    }

  return REF_SUCCESS;
}

static REF_STATUS ref_metric_interpolate_twod( REF_GRID ref_grid, REF_GRID parent_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_NODE parent_node = ref_grid_node(parent_grid);
  REF_INT node, tri, ixyz, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[3], interpolated_xyz[3], bary[3];
  REF_DBL tol = 1.0e-11;
  REF_DBL log_parent_m[3][6];
  REF_DBL log_interpolated_m[6];

  if (!ref_grid_twod(ref_grid))
    RSS(REF_IMPLEMENT,"ref_metric_interpolate only implemented for twod");

  each_ref_node_valid_node( ref_node, node )
    {
      for (ixyz=0; ixyz<3; ixyz++)
	xyz[ixyz] = ref_node_xyz(ref_node,ixyz,node); 
      tri = ref_node_guess(ref_node,node);
      RSS( ref_grid_enclosing_tri( parent_grid, xyz,
				   &tri, bary ), "enclosing tri" );
      if ( ref_node_guess_allocated(ref_node) )
	ref_node_raw_guess(ref_node,node) = tri;
      RSS( ref_cell_nodes( ref_grid_tri(parent_grid), tri, nodes ), "c2n");
      for (ixyz=0; ixyz<3; ixyz++)
	{
	  interpolated_xyz[ixyz] = 0.0;
	  for (ibary=0; ibary<3; ibary++)
	    interpolated_xyz[ixyz] += 
	      bary[ibary]*ref_node_real(parent_node,ixyz,nodes[ibary]);
	}
      /* override y for fake twod */
      interpolated_xyz[1] = ref_node_xyz(ref_node,1,node);  
      for (ixyz=0; ixyz<3; ixyz++)
	RWDS( xyz[ixyz], interpolated_xyz[ixyz], tol, "xyz check");
      for (ibary=0; ibary<3; ibary++)
	RSS( ref_matrix_log_m( ref_node_metric_ptr(parent_node,nodes[ibary]),
			       log_parent_m[ibary] ), "log(parentM)");
      for (im=0; im<6; im++)
	{
	  log_interpolated_m[im] = 0.0;
	  for (ibary=0; ibary<3; ibary++)
	    log_interpolated_m[im] += 
	      bary[ibary]*log_parent_m[ibary][im];
	}
      RSS(ref_matrix_exp_m( log_interpolated_m, 
			    ref_node_metric_ptr(ref_node,node) ),"exp(intrpM)");
    }

  return REF_SUCCESS;
}


REF_STATUS ref_metric_interpolate( REF_GRID ref_grid, REF_GRID parent_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_NODE parent_node = ref_grid_node(parent_grid);
  REF_INT node, tet, ixyz, ibary, im;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL xyz[3], interpolated_xyz[3], bary[4];
  REF_DBL tol = 1.0e-11;
  REF_DBL log_parent_m[4][6];
  REF_DBL log_interpolated_m[6];
  
  if (ref_grid_twod(ref_grid))
    {
      RSS( ref_metric_interpolate_twod( ref_grid, parent_grid ), "2d version");
      return REF_SUCCESS;
    }
  
  each_ref_node_valid_node( ref_node, node )
    {
      for (ixyz=0; ixyz<3; ixyz++)
	xyz[ixyz] = ref_node_xyz(ref_node,ixyz,node); 
      tet = ref_node_guess(ref_node,node);
      RSS( ref_grid_enclosing_tet( parent_grid, xyz,
				   &tet, bary ), "enclosing tet" );
      if ( ref_node_guess_allocated(ref_node) )
	ref_node_raw_guess(ref_node,node) = tet;
      RSS( ref_cell_nodes( ref_grid_tet(parent_grid), tet, nodes ), "c2n");
      for (ixyz=0; ixyz<3; ixyz++)
	{
	  interpolated_xyz[ixyz] = 0.0;
	  for (ibary=0; ibary<4; ibary++)
	    interpolated_xyz[ixyz] += 
	      bary[ibary]*ref_node_real(parent_node,ixyz,nodes[ibary]);
	}
      for (ixyz=0; ixyz<3; ixyz++)
	RWDS( xyz[ixyz], interpolated_xyz[ixyz], tol, "xyz check");
      for (ibary=0; ibary<4; ibary++)
	RSS( ref_matrix_log_m( ref_node_metric_ptr(parent_node,nodes[ibary]),
			       log_parent_m[ibary] ), "log(parentM)");
      for (im=0; im<6; im++)
	{
	  log_interpolated_m[im] = 0.0;
	  for (ibary=0; ibary<4; ibary++)
	    log_interpolated_m[im] += 
	      bary[ibary]*log_parent_m[ibary][im];
	}
      RSS(ref_matrix_exp_m( log_interpolated_m, 
			    ref_node_metric_ptr(ref_node,node) ),"exp(intrpM)");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_metric_gradation( REF_DBL *metric, REF_GRID ref_grid, REF_DBL r )
{
  REF_EDGE ref_edge;
  REF_DBL *metric_orig;
  REF_DBL *metric_limit;
  REF_DBL ratio, lr;
  REF_DBL l0[6], l1[6];
  REF_DBL m0[6], m1[6];
  REF_INT node, i;
  REF_INT edge, node0, node1;

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc( metric_orig, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  ref_malloc( metric_limit, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  
  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      for (i=0;i<6;i++) metric_orig[i+6*node] = metric[i+6*node];
    }
  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      for (i=0;i<6;i++) metric_limit[i+6*node] = metric[i+6*node]*(1.0/r/r);
    }

  each_ref_edge( ref_edge, edge )
    {
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );
      RSS( ref_node_ratio( ref_grid_node(ref_grid), 
			   node0, node1, &ratio),"ratio");
      lr = pow(r,ratio);
      for (i=0;i<6;i++) l0[i] = metric_limit[i+6*node0] * (1.0/lr/lr);
      for (i=0;i<6;i++) l1[i] = metric_limit[i+6*node1] * (1.0/lr/lr);
      RSS( ref_matrix_intersect( &(metric_orig[6*node0]), 
				 l1,
				 m0 ), "m0" );  
      RSS( ref_matrix_intersect( &(metric_orig[6*node1]), 
				 l0,
				 m1 ), "m1" );
      RSS( ref_matrix_intersect( m0,
				 &(metric[6*node0]), 
				 &(metric[6*node0]) ), "m0" );  
      RSS( ref_matrix_intersect( m1,
				 &(metric[6*node1]), 
				 &(metric[6*node1]) ), "m0" );  
    }

  ref_free( metric_limit );
  ref_free( metric_orig );

  ref_edge_free( ref_edge );

  return REF_SUCCESS;
}

REF_STATUS ref_metric_sanitize( REF_GRID ref_grid )
{
  if (ref_grid_twod(ref_grid)) 
    {
      RSS( ref_metric_sanitize_twod( ref_grid ), "threed" );
    }
  else
    {
      RSS( ref_metric_sanitize_threed( ref_grid ), "threed" );
    }
 
  return REF_SUCCESS;
}

REF_STATUS ref_metric_sanitize_threed( REF_GRID ref_grid )
{
  REF_DBL *metric_orig;
  REF_DBL *metric_imply;
  REF_DBL *metric;

  ref_malloc( metric_orig, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  ref_malloc( metric_imply, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  ref_malloc( metric, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

  RSS( ref_metric_from_node( metric_orig, ref_grid_node(ref_grid)), "from");
  
  RSS( ref_metric_imply_from( metric_imply, ref_grid ), "imply" );

  RSS( ref_metric_smr( metric_imply, metric_orig, metric, ref_grid ), "smr" );

  RSS( ref_metric_imply_non_tet( metric, ref_grid ), "imply non tet");

  RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid)), "to");

  ref_free( metric );
  ref_free( metric_imply );
  ref_free( metric_orig );

  return REF_SUCCESS;
}
REF_STATUS ref_metric_sanitize_twod( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *metric_orig;
  REF_DBL *metric_imply;
  REF_DBL *metric;
  REF_INT node;

  ref_malloc( metric_orig, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  ref_malloc( metric_imply, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );
  ref_malloc( metric, 
	      6*ref_node_max(ref_grid_node(ref_grid)), REF_DBL );

  RSS( ref_metric_from_node( metric_orig, ref_grid_node(ref_grid)), "from");
  for( node=0; node<ref_node_max(ref_node); node++ )
    {
      metric_orig[1+6*node] = 0.0;
      metric_orig[3+6*node] = 1.0;
      metric_orig[4+6*node] = 0.0;
    }
  
  RSS( ref_metric_imply_from( metric_imply, ref_grid ), "imply" );
  for( node=0; node<ref_node_max(ref_node); node++ )
    {
      metric_imply[1+6*node] = 0.0;
      metric_imply[3+6*node] = 1.0;
      metric_imply[4+6*node] = 0.0;
    }

  RSS( ref_metric_smr( metric_imply, metric_orig, metric, ref_grid ), "smr" );

  RSS( ref_metric_to_node( metric, ref_grid_node(ref_grid)), "to");

  ref_free( metric );
  ref_free( metric_imply );
  ref_free( metric_orig );

  return REF_SUCCESS;
}

REF_STATUS ref_metric_interpolated_curvature( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *metric;
  REF_INT gradation;

  ref_malloc( metric, 6*ref_node_max(ref_node), REF_DBL );
  RSS( ref_metric_from_curvature( metric, ref_grid ), "curve" );  
  for ( gradation =0 ; gradation<10 ; gradation++ )
    {
      RSS( ref_metric_gradation( metric, ref_grid, 1.2 ), "grad");
    }
  RSS(ref_metric_to_node( metric, ref_node ), "to node");
  ref_free( metric );

  return REF_SUCCESS;
}

REF_STATUS ref_metric_constrain_curvature( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *curvature_metric;
  REF_DBL m[6];
  REF_INT node, im;
  
  if ( !ref_geom_model_loaded(ref_grid_geom(ref_grid)) )
    {
      printf("No geometry model loaded, skipping curvature constraint.\n");
      return REF_SUCCESS;
    }

  ref_malloc( curvature_metric, 6*ref_node_max(ref_node), REF_DBL );

  RSS( ref_metric_from_curvature( curvature_metric, ref_grid ), "curve" );

  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_matrix_intersect( &(curvature_metric[6*node]), 
				 ref_node_metric_ptr(ref_node,node),
				 m ), "intersect" );  
      for(im=0;im<6;im++) ref_node_metric(ref_node,im,node) = m[im];
    }
  
  ref_free( curvature_metric );

  return REF_SUCCESS;
}

REF_STATUS ref_metric_from_curvature( REF_DBL *metric, REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_GEOM ref_geom = ref_grid_geom(ref_grid);
  REF_INT geom, node;
  REF_DBL kr, r[3], ks, s[3], n[3];
  REF_DBL diagonal_system[12];
  REF_INT i;
  REF_DBL drad;
  REF_DBL hmax;
  REF_DBL rlimit;
  REF_DBL h, hr, hs, hn;
  REF_DBL curvature_ratio, norm_ratio;

  if ( !ref_geom_model_loaded(ref_geom) )
    {
      printf("\nNo geometry model, did you forget to load it?\n\n");
      RSS(REF_IMPLEMENT,"...or implement non-CAD curvature estimate")
    }

  drad = 1.0/10.0; /* 1/segments per radian */
  RSS( ref_geom_egads_diagonal( ref_geom, &hmax ), "bbox diag");
  hmax *= 0.1; /* normal spacing and max tangential spacing */
  rlimit = hmax/drad; /* h = r*drad, r = h/drad */
  curvature_ratio = 1.0/20.0;
  norm_ratio = 5.0;

  each_ref_node_valid_node( ref_node, node )
    {
      h = hmax;
      metric[0+6*node] = 1.0/(h*h);
      metric[1+6*node] = 0.0;
      metric[2+6*node] = 0.0;
      metric[3+6*node] = 1.0/(h*h);
      metric[4+6*node] = 0.0;
      metric[5+6*node] = 1.0/(h*h);
    }

  each_ref_geom_face( ref_geom, geom )
    {
      RSS( ref_geom_curvature( ref_geom, geom, &kr, r, &ks, s ), "curve" );
      kr = ABS(kr);
      ks = ABS(ks);
      kr = MAX(kr,curvature_ratio*ks);
      ks = MAX(ks,curvature_ratio*kr);
      ref_math_cross_product( r, s, n );
      node = ref_geom_node(ref_geom,geom);
      for ( i=0 ; i<3 ; i++ )
	ref_matrix_vec(diagonal_system, i, 0 ) = r[i];
      hr = hmax;
      if ( 1.0/rlimit < kr )
	hr = drad/kr;
      ref_matrix_eig(diagonal_system, 0 ) = 1.0/hr/hr;
      for ( i=0 ; i<3 ; i++ )
	ref_matrix_vec(diagonal_system, i, 1 ) = s[i];
      hs = hmax;
      if ( 1.0/rlimit < ks )
	hs = drad/ks;
      ref_matrix_eig(diagonal_system, 1 ) = 1.0/hs/hs;
      for ( i=0 ; i<3 ; i++ )
	ref_matrix_vec(diagonal_system, i, 2 ) = n[i];
      hn = hmax;
      hn = MIN(hn, norm_ratio*hr);
      hn = MIN(hn, norm_ratio*hs);
      ref_matrix_eig(diagonal_system, 2 ) = 1.0/hn/hn;
      RSS( ref_matrix_form_m( diagonal_system, &(metric[6*node]) ), "reform m");
    }
  
  return REF_SUCCESS;
}

#define sub_tet_contribution(n0,n1,n2,n3)	\
  {						\
    tet_nodes[0] = nodes[(n0)];			\
    tet_nodes[1] = nodes[(n1)];			\
    tet_nodes[2] = nodes[(n2)];			\
    tet_nodes[3] = nodes[(n3)];			\
    RSS( ref_matrix_imply_m( m,						\
                             ref_node_xyz_ptr(ref_node,tet_nodes[0]),   \
                             ref_node_xyz_ptr(ref_node,tet_nodes[1]),   \
                             ref_node_xyz_ptr(ref_node,tet_nodes[2]),	\
                             ref_node_xyz_ptr(ref_node,tet_nodes[3])),"impl"); \
    RSS( ref_matrix_log_m( m, log_m ), "log" );				\
    RSS( ref_node_tet_vol( ref_node, tet_nodes, &tet_volume ), "vol" );	\
    for( node=0; node<ref_cell_node_per(ref_cell); node++ )	       	\
      {									\
        total_node_volume[nodes[node]] += tet_volume;			\
        for( im=0; im<6; im++ )						\
          metric[im+6*nodes[node]] += tet_volume*log_m[im];		\
      }									\
  }

REF_STATUS ref_metric_imply_from( REF_DBL *metric, REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_DBL *total_node_volume, tet_volume;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell;
  REF_CELL ref_cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init( total_node_volume, ref_node_max(ref_node),
		   REF_DBL, 0.0 );

  for( node=0; node<ref_node_max(ref_node); node++ )
    for( im=0; im<6; im++ )
      metric[im+6*node] = 0.0;

  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      sub_tet_contribution(0,1,2,3);
    }

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,4,5,3);
      sub_tet_contribution(0,1,5,4);
      sub_tet_contribution(0,1,2,5);
    }

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,4,1,2);
      sub_tet_contribution(0,3,4,2);
    }

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,5,7,4);
      sub_tet_contribution(0,1,7,5);
      sub_tet_contribution(1,6,7,5);

      sub_tet_contribution(0,7,2,3);
      sub_tet_contribution(0,7,1,2);
      sub_tet_contribution(1,7,6,2);
    }

  each_ref_node_valid_node( ref_node, node )
    {
      RAS( 0.0 < total_node_volume[node], "zero metric contributions" );
      for( im=0; im<6; im++ )
	{
	  if ( !ref_math_divisible( metric[im+6*node], 
				    total_node_volume[node]) ) 
	    RSS( REF_DIV_ZERO, "zero volume");
	  log_m[im] = metric[im+6*node] / total_node_volume[node];
	}
      RSS( ref_matrix_exp_m( log_m, m ), "exp" );
      for( im=0; im<6; im++ )
	metric[im+6*node] = m[im];
      total_node_volume[node] = 0.0;
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
  REF_CELL ref_cell;
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER], nodes[REF_CELL_MAX_SIZE_PER];

  ref_malloc_init( total_node_volume, ref_node_max(ref_node),
		   REF_DBL, 0.0 );

  
  each_ref_node_valid_node( ref_node, node )
    if ( ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_pyr(ref_grid)), 
				       node ) ) ||
	 ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_pri(ref_grid)), 
				       node ) ) ||
	 ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_hex(ref_grid)), 
				       node ) ) )
      for( im=0; im<6; im++ )
	metric[im+6*node] = 0.0;

  ref_cell = ref_grid_pri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,4,5,3);
      sub_tet_contribution(0,1,5,4);
      sub_tet_contribution(0,1,2,5);
    }

  ref_cell = ref_grid_pyr(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,4,1,2);
      sub_tet_contribution(0,3,4,2);
    }

  /*
VI1 VI6 VI8 VI5  VI1 VI2 VI8 VI6  VI2 VI7 VI8 VI6  
VI1 VI8 VI3 VI4  VI1 VI8 VI2 VI3  VI2 VI8 VI7 VI3
  */

  ref_cell = ref_grid_hex(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      sub_tet_contribution(0,5,7,4);
      sub_tet_contribution(0,1,7,5);
      sub_tet_contribution(1,6,7,5);

      sub_tet_contribution(0,7,2,3);
      sub_tet_contribution(0,7,1,2);
      sub_tet_contribution(1,7,6,2);
    }

  each_ref_node_valid_node( ref_node, node )
    if ( ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_pyr(ref_grid)), 
				       node ) ) ||
	 ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_pri(ref_grid)), 
				       node ) ) ||
	 ref_adj_valid( ref_adj_first( ref_cell_adj(ref_grid_hex(ref_grid)), 
				       node ) ) )
      {
	RAS( 0.0 < total_node_volume[node], "zero metric contributions" );
	for( im=0; im<6; im++ )
	  {
	    if ( !ref_math_divisible( metric[im+6*node], 
				      total_node_volume[node]) ) 
	      RSS( REF_DIV_ZERO, "zero volume");
	    log_m[im] = metric[im+6*node] / total_node_volume[node];
	  }
	RSS( ref_matrix_exp_m( log_m, m ), "exp" );
	for( im=0; im<6; im++ )
	  metric[im+6*node] = m[im];
	total_node_volume[node] = 0.0;
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
      if ( REF_SUCCESS != ref_matrix_inv_gen( 3, n_vectors, inv_n_vectors ) )
	{
	  printf(" unable to invert eigenvectors:\n");
	  printf(" %f %f %f\n",n_vectors[0],n_vectors[1],n_vectors[2]);
	  printf(" %f %f %f\n",n_vectors[3],n_vectors[4],n_vectors[5]);
	  printf(" %f %f %f\n",n_vectors[6],n_vectors[7],n_vectors[8]);
	  RSS( REF_FAILURE, "gen eig" );
	}
      RSS( ref_matrix_transpose_gen( 3, inv_n_vectors, 
				     &(diagonal_system[3]) ), "gen eig");
      RSS( ref_matrix_form_m( diagonal_system, &(metric[6*node]) ), "reform m");
    }

  return REF_SUCCESS;
}

