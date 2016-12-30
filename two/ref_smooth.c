
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_smooth.h"
#include "ref_adapt.h"
#include "ref_math.h"
#include "ref_matrix.h"
#include "ref_cell.h"
#include "ref_mpi.h"
#include "ref_twod.h"

REF_STATUS ref_smooth_tri_steepest_descent( REF_GRID ref_grid, REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL f, d[3];
  REF_DBL dcost, dcost_dl, dl;
  REF_BOOL verbose = REF_FALSE;

  each_ref_cell_having_node( ref_cell, node, item, cell )
  {
    RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
    RSS( ref_node_tri_dquality_dnode0(ref_node, nodes,
                                      &f, d), "qual deriv" );
    if (verbose)
      printf("cost %10.8f : %12.8f %12.8f %12.8f\n",f,d[0],d[1],d[2]);
  }

  dcost = 1.0-f;
  dcost_dl = sqrt(d[0]*d[0]+d[1]*d[1]*d[2]*d[2]);
  dl = dcost/dcost_dl;

  ref_node_xyz(ref_node,0,node) += dl*d[0];
  ref_node_xyz(ref_node,1,node) += dl*d[1];
  ref_node_xyz(ref_node,2,node) += dl*d[2];

  each_ref_cell_having_node( ref_cell, node, item, cell )
  {
    RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
    RSS( ref_node_tri_dquality_dnode0(ref_node, nodes,
                                      &f, d), "qual deriv" );
  }

  if (verbose)
    printf("rate %12.8f dcost %12.8f dl %12.8f\n",
           (1.0-f)/dcost,dcost,dl);


  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_quality_around( REF_GRID ref_grid,
                                          REF_INT node,
                                          REF_DBL *min_quality )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL quality;

  *min_quality = 1.0;
  each_ref_cell_having_node( ref_cell, node, item, cell )
  {
    RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
    none_found = REF_FALSE;
    RSS( ref_node_tri_quality( ref_node,
                               nodes,
                               &quality ), "qual" );
    *min_quality = MIN( *min_quality, quality );
  }

  if ( none_found )
    {
      *min_quality = -2.0;
      return REF_NOT_FOUND;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_outward_norm( REF_GRID ref_grid, 
				    REF_INT node,
				    REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normal[3];

  *allowed = REF_FALSE;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );

      RSS( ref_node_tri_normal( ref_node,nodes,normal ), "norm");

      if ( ( ref_node_xyz(ref_node,1,nodes[0]) > 0.5 &&
	     normal[1] >= 0.0 ) ||
	   ( ref_node_xyz(ref_node,1,nodes[0]) < 0.5 &&
	     normal[1] <= 0.0 ) ) return REF_SUCCESS;
    }

  *allowed = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_ideal( REF_GRID ref_grid,
                                 REF_INT node,
                                 REF_INT tri,
                                 REF_DBL *ideal_location )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT n0, n1;
  REF_INT ixyz;
  REF_DBL dn[3];
  REF_DBL dt[3];
  REF_DBL tangent_length, projection, scale, length_in_metric;

  RSS(ref_cell_nodes(ref_grid_tri(ref_grid), tri, nodes ), "get tri");
  n0 = REF_EMPTY; n1 = REF_EMPTY;
  if ( node == nodes[0])
    {
      n0 = nodes[1];
      n1 = nodes[2];
    }
  if ( node == nodes[1])
    {
      n0 = nodes[2];
      n1 = nodes[0];
    }
  if ( node == nodes[2])
    {
      n0 = nodes[0];
      n1 = nodes[1];
    }
  if ( n0==REF_EMPTY || n1==REF_EMPTY)
    THROW("empty triangle side");

  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] = 0.5*( ref_node_xyz(ref_node,ixyz,n0) +
                                 ref_node_xyz(ref_node,ixyz,n1) );
  for (ixyz = 0; ixyz<3; ixyz++)
    dn[ixyz] = ref_node_xyz(ref_node,ixyz,node) - ideal_location[ixyz];
  for (ixyz = 0; ixyz<3; ixyz++)
    dt[ixyz] = ref_node_xyz(ref_node,ixyz,n1) - ref_node_xyz(ref_node,ixyz,n0);

  tangent_length = ref_math_dot(dt,dt);
  projection = ref_math_dot(dn,dt);

  if ( ref_math_divisible(projection,tangent_length) )
    {
      for (ixyz = 0; ixyz<3; ixyz++)
        dn[ixyz] -=  (projection/tangent_length) * dt[ixyz];
    }
  else
    {
      printf("projection = %e tangent_length = %e\n",
             projection,tangent_length);
      return REF_DIV_ZERO;
    }

  RSS( ref_math_normalize( dn ), "normalize direction" );
  /* would an averaged metric be more appropriate? */
  length_in_metric =
    ref_matrix_sqrt_vt_m_v( ref_node_metric_ptr(ref_node,node), dn );

  scale = 0.5*sqrt(3.0); /* altitude of equilateral triangle */
  if ( ref_math_divisible(scale,length_in_metric) )
    {
      scale = scale/length_in_metric;
    }
  else
    {
      printf(" length_in_metric = %e, not invertable\n", length_in_metric);
      return REF_DIV_ZERO;
    }

  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] += scale*dn[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_weighted_ideal( REF_GRID ref_grid,
					  REF_INT node,
					  REF_DBL *ideal_location )
{
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ixyz;
  REF_DBL tri_ideal[3];
  REF_DBL quality, weight, normalization;

  normalization = 0.0;
  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] = 0.0;

  each_ref_cell_having_node( ref_grid_tri(ref_grid), node, item, cell )
    {
      RSS( ref_smooth_tri_ideal( ref_grid, node, cell, 
				 tri_ideal ), "tri ideal");
      RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes ), "nodes" );
      RSS( ref_node_tri_quality( ref_grid_node(ref_grid), 
				 nodes,  
				 &quality ), "tri qual");
      quality = MAX(quality,ref_adapt_smooth_min_quality);
      weight = 1.0/quality;
      normalization += weight;
      for (ixyz = 0; ixyz<3; ixyz++)
	ideal_location[ixyz] += weight*tri_ideal[ixyz];
    }

  if ( ref_math_divisible(1.0,normalization) )
    {
      for (ixyz = 0; ixyz<3; ixyz++)
        ideal_location[ixyz] =  (1.0/normalization) * ideal_location[ixyz];
    }
  else
    {
      printf("normalization = %e\n",normalization);
      return REF_DIV_ZERO;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tri_improve( REF_GRID ref_grid,
				   REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality0, quality;
  REF_INT ixyz, opposite;
  REF_BOOL allowed;

  /* can't handle boundaries yet */
  if ( !ref_cell_node_empty( ref_grid_qua( ref_grid ), node ) )
    return REF_SUCCESS;

  for (ixyz = 0; ixyz<3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node,ixyz,node);

  RSS( ref_smooth_tri_weighted_ideal( ref_grid, node, ideal ), "ideal" );

  RSS( ref_smooth_tri_quality_around( ref_grid, node, &quality0),"q");

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++)
    {
      for (ixyz = 0; ixyz<3; ixyz++)
	ref_node_xyz(ref_node,ixyz,node) = backoff*ideal[ixyz] +
	  (1.0 - backoff) * original[ixyz];
      RSS( ref_smooth_outward_norm( ref_grid, node, &allowed ), "normals" );
      if ( allowed )
	{
	  RSS( ref_smooth_tri_quality_around( ref_grid, node, &quality),"q");
	  if ( quality > quality0 )
	    {
	      /* update opposite side: X and Z only */
	      RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid), 
					   node, &opposite ),"opp");
	      ref_node_xyz(ref_node,0,opposite) = ref_node_xyz(ref_node,0,node);
	      ref_node_xyz(ref_node,2,opposite) = ref_node_xyz(ref_node,2,node);
	      return REF_SUCCESS;
	    }
	}
      backoff *= 0.5;
    }

  for (ixyz = 0; ixyz<3; ixyz++)
    ref_node_xyz(ref_node,ixyz,node) = original[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_local_pris_about( REF_GRID ref_grid, 
					REF_INT about_node, 
					REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, node;

  *allowed =  REF_FALSE;

  ref_cell = ref_grid_pri(ref_grid);

  each_ref_cell_having_node( ref_cell, about_node, item, cell )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
      if ( ref_mpi_id != ref_node_part(ref_node,
				       ref_cell_c2n(ref_cell,node,cell)) )
	return REF_SUCCESS;

  *allowed =  REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_twod_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node;
  REF_BOOL allowed;

  each_ref_node_valid_node( ref_node, node )
    {
      RSS( ref_node_node_twod( ref_node, node, &allowed ), "twod" );
      if ( !allowed ) continue;

      /* can't handle boundaries yet */
      allowed = ref_cell_node_empty( ref_grid_qua( ref_grid ), node );
      if ( !allowed ) continue;

      RSS( ref_smooth_local_pris_about( ref_grid, node, &allowed ), "para" );
      if ( !allowed )
	{
	  ref_node_age(ref_node,node)++;
	  continue;
	}

      ref_node_age(ref_node,node) = 0;
      RSS( ref_smooth_tri_improve( ref_grid, node ), "improve" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_quality_around( REF_GRID ref_grid,
                                          REF_INT node,
                                          REF_DBL *min_quality )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL none_found = REF_TRUE;
  REF_DBL quality;

  *min_quality = 1.0;
  each_ref_cell_having_node( ref_cell, node, item, cell )
  {
    RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
    none_found = REF_FALSE;
    RSS( ref_node_tet_quality( ref_node,
                               nodes,
                               &quality ), "qual" );
    *min_quality = MIN( *min_quality, quality );
  }

  if ( none_found )
    {
      *min_quality = -2.0;
      return REF_NOT_FOUND;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_ideal( REF_GRID ref_grid,
                                 REF_INT node,
                                 REF_INT tet,
                                 REF_DBL *ideal_location )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_node, tri_node;
  REF_INT ixyz;
  REF_DBL dn[3];
  REF_DBL scale, length_in_metric;

  RSS(ref_cell_nodes(ref_cell, tet, nodes ), "get tri");
  tri_nodes[0] = REF_EMPTY; tri_nodes[1] = REF_EMPTY; tri_nodes[2] = REF_EMPTY;
  for ( tet_node=0 ; tet_node<4; tet_node++ )
    if ( node == nodes[tet_node])
      {
	for ( tri_node=0 ; tri_node<3; tri_node++ )
	  tri_nodes[tri_node] = 
	    nodes[ref_cell_f2n_gen(ref_cell,tri_node,tet_node)];
      }
  if ( tri_nodes[0] == REF_EMPTY ||
       tri_nodes[1] == REF_EMPTY ||
       tri_nodes[2] == REF_EMPTY )
    THROW("empty tetrahedra face");

  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] = ( ref_node_xyz(ref_node,ixyz,tri_nodes[0]) +
			     ref_node_xyz(ref_node,ixyz,tri_nodes[1]) +
			     ref_node_xyz(ref_node,ixyz,tri_nodes[2]) ) / 3.0;

  RSS( ref_node_tri_normal( ref_node, tri_nodes, dn ), "tri normal" );

  RSS( ref_math_normalize( dn ), "normalize direction" );
  /* would an averaged metric be more appropriate? */
  length_in_metric =
    ref_matrix_sqrt_vt_m_v( ref_node_metric_ptr(ref_node,node), dn );

  scale = sqrt(6.0)/3.0; /* altitude of regular tetrahedra */
  if ( ref_math_divisible(scale,length_in_metric) )
    {
      scale = scale/length_in_metric;
    }
  else
    {
      printf(" length_in_metric = %e, not invertable\n", length_in_metric);
      return REF_DIV_ZERO;
    }

  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] += scale*dn[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_weighted_ideal( REF_GRID ref_grid,
					  REF_INT node,
					  REF_DBL *ideal_location )
{
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ixyz;
  REF_DBL tet_ideal[3];
  REF_DBL quality, weight, normalization;

  normalization = 0.0;
  for (ixyz = 0; ixyz<3; ixyz++)
    ideal_location[ixyz] = 0.0;

  each_ref_cell_having_node( ref_grid_tet(ref_grid), node, item, cell )
    {
      RSS( ref_smooth_tet_ideal( ref_grid, node, cell, 
				 tet_ideal ), "tet ideal");
      RSS( ref_cell_nodes( ref_grid_tet(ref_grid), cell, nodes ), "nodes" );
      RSS( ref_node_tet_quality( ref_grid_node(ref_grid), 
				 nodes,  
				 &quality ), "tet qual");
      quality = MAX(quality,ref_adapt_smooth_min_quality);
      weight = 1.0/quality;
      normalization += weight;
      for (ixyz = 0; ixyz<3; ixyz++)
	ideal_location[ixyz] += weight*tet_ideal[ixyz];
    }

  if ( ref_math_divisible(1.0,normalization) )
    {
      for (ixyz = 0; ixyz<3; ixyz++)
        ideal_location[ixyz] =  (1.0/normalization) * ideal_location[ixyz];
    }
  else
    {
      printf("normalization = %e\n",normalization);
      return REF_DIV_ZERO;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_tet_improve( REF_GRID ref_grid,
				   REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT tries;
  REF_DBL ideal[3], original[3];
  REF_DBL backoff, quality0, quality;
  REF_INT ixyz;

  /* can't handle boundaries yet */
  if ( !ref_cell_node_empty( ref_grid_tri( ref_grid ), node ) ||
       !ref_cell_node_empty( ref_grid_qua( ref_grid ), node ) )
    return REF_SUCCESS;

  for (ixyz = 0; ixyz<3; ixyz++)
    original[ixyz] = ref_node_xyz(ref_node,ixyz,node);

  RSS( ref_smooth_tet_weighted_ideal( ref_grid, node, ideal ), "ideal" );

  RSS( ref_smooth_tet_quality_around( ref_grid, node, &quality0),"q");

  backoff = 1.0;
  for (tries = 0; tries < 8; tries++)
    {
      for (ixyz = 0; ixyz<3; ixyz++)
	ref_node_xyz(ref_node,ixyz,node) = backoff*ideal[ixyz] +
	  (1.0 - backoff) * original[ixyz];
      RSS( ref_smooth_tet_quality_around( ref_grid, node, &quality),"q");
      if ( quality > quality0 )
	{
	  return REF_SUCCESS;
	}
      backoff *= 0.5;
    }

  for (ixyz = 0; ixyz<3; ixyz++)
    ref_node_xyz(ref_node,ixyz,node) = original[ixyz];

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_local_tet_about( REF_GRID ref_grid, 
				       REF_INT about_node, 
				       REF_BOOL *allowed )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT item, cell, node;

  *allowed =  REF_FALSE;

  ref_cell = ref_grid_tet(ref_grid);

  each_ref_cell_having_node( ref_cell, about_node, item, cell )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
      if ( ref_mpi_id != ref_node_part(ref_node,
				       ref_cell_c2n(ref_cell,node,cell)) )
	return REF_SUCCESS;

  *allowed =  REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_smooth_threed_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT node, i;
  REF_BOOL allowed, it_is;

  each_ref_node_valid_node( ref_node, node )
    {
      /* don't move geom nodes */
      RSS( ref_geom_is_a(ref_grid_geom(ref_grid), node, REF_GEOM_NODE, &it_is),
	   "node check");
      if ( it_is ) continue;
      
      /* can't handle boundaries yet */
      allowed = ref_cell_node_empty( ref_grid_tri( ref_grid ), node ) ||
                ref_cell_node_empty( ref_grid_qua( ref_grid ), node );
      if ( !allowed ) continue;

      RSS( ref_smooth_local_tet_about( ref_grid, node, &allowed ), "para" );
      if ( !allowed )
	{
	  ref_node_age(ref_node,node)++;
	  continue;
	}

      RSS( ref_smooth_tet_improve( ref_grid, node ), "ideal node for tet" );
      
      ref_node_age(ref_node,node) = 0;
      for ( i = 0 ; i < 3 ; i++ )
	RSS( ref_smooth_tet_improve( ref_grid, node ), "improve" );
    }

  return REF_SUCCESS;
}

