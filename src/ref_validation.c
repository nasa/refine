
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_validation.h"

#include "ref_face.h"
#include "ref_matrix.h"
#include "ref_export.h"
#include "ref_mpi.h"

REF_STATUS ref_validation_all( REF_GRID ref_grid )
{
  RSS( ref_validation_cell_face( ref_grid ), "cell face");
  RSS( ref_validation_cell_node( ref_grid ), "cell node valid");
  RSS( ref_validation_unused_node( ref_grid ), "unused node");

  return REF_SUCCESS;
}

REF_STATUS ref_validation_unused_node( REF_GRID ref_grid )
{
  REF_INT node;
  REF_BOOL problem;
  REF_ADJ ref_adj;
  REF_CELL ref_cell;
  REF_INT group, cell, nodes[REF_CELL_MAX_SIZE_PER];

  problem = REF_FALSE;
  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      if ( ref_adj_empty( ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
	   ref_adj_empty( ref_cell_adj(ref_grid_hex(ref_grid)), node) )
	{
	  problem = REF_TRUE;
	  printf(" unused node %d: %e %e %e\n",
		 node,
		 ref_node_xyz(ref_grid_node(ref_grid),0,node),
		 ref_node_xyz(ref_grid_node(ref_grid),1,node),
		 ref_node_xyz(ref_grid_node(ref_grid),2,node));
	}
    }

  ref_adj_create( &ref_adj );

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  for ( node = 0; node < ref_cell_node_per(ref_cell); node++ )
	    RSS( ref_adj_add( ref_adj, nodes[node], group+4*cell ), "add");
	}
    }

  each_ref_node_valid_node( ref_grid_node(ref_grid), node )
    {
      if ( ref_adj_empty( ref_adj, node ) )
	{
	  problem = REF_TRUE;
	  printf(" unused node %d\n",node);	  
	}
    }

  ref_adj_free(ref_adj);

  return (problem?REF_FAILURE:REF_SUCCESS);
}

REF_STATUS ref_validation_boundary_face( REF_GRID ref_grid )
{
  REF_CELL ref_cell;
  REF_BOOL has_face;
  REF_INT cell;
  REF_INT node;
  REF_INT nodes[4];
  REF_BOOL problem;

  problem = REF_FALSE;
 
  ref_cell = ref_grid_tri( ref_grid );
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      for(node=0;node<3;node++)
	nodes[node]=ref_cell_c2n(ref_cell,node,cell);
      nodes[3]=nodes[0];
      RSS( ref_grid_cell_has_face( ref_grid, nodes, &has_face ), "has_face");
      if ( !has_face ) problem = REF_TRUE;
    }
 
  ref_cell = ref_grid_qua( ref_grid );
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      for(node=0;node<4;node++)
	nodes[node]=ref_cell_c2n(ref_cell,node,cell);
      RSS( ref_grid_cell_has_face( ref_grid, nodes, &has_face ), "has_face");
      if ( !has_face ) problem = REF_TRUE;
    }
 
  return (problem?REF_FAILURE:REF_SUCCESS);
}

REF_STATUS ref_validation_cell_face( REF_GRID ref_grid )
{
  REF_FACE ref_face;
  REF_CELL ref_cell;
  REF_INT *hits;
  REF_INT face;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[4];
  REF_BOOL problem;
  REF_STATUS code;

  problem = REF_FALSE;

  RSS( ref_face_create( &ref_face, ref_grid ), "face");

  hits = (REF_INT *)malloc( ref_face_n(ref_face) * sizeof(REF_INT) );
  RNS(hits,"malloc hits NULL");

  for ( face=0; face< ref_face_n(ref_face) ; face++ )
    hits[face]=0;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_face( ref_cell, cell_face )
        {
	  for(node=0;node<4;node++)
	    nodes[node]=ref_cell_f2n(ref_cell,node,cell_face,cell);
	  RSS( ref_face_with( ref_face, nodes, &face ), "find cell face");
	  hits[face]++;
	}
 
  ref_cell = ref_grid_tri( ref_grid );
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      for(node=0;node<3;node++)
	nodes[node]=ref_cell_c2n(ref_cell,node,cell);
      nodes[3]=nodes[0];
      code = ref_face_with( ref_face, nodes, &face );
      if ( REF_SUCCESS != code)
	{
	  ref_node_location( ref_grid_node(ref_grid), nodes[0] );
	  ref_node_location( ref_grid_node(ref_grid), nodes[1] );
	  ref_node_location( ref_grid_node(ref_grid), nodes[2] );
	  ref_node_location( ref_grid_node(ref_grid), nodes[3] );
	}
      RSS( code, "find tri");
      hits[face]++;
    }
 
  ref_cell = ref_grid_qua( ref_grid );
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      for(node=0;node<4;node++)
	nodes[node]=ref_cell_c2n(ref_cell,node,cell);
      RSS( ref_face_with( ref_face, nodes, &face ), "find qua");
      hits[face]++;
    }
 
  for ( face=0; face< ref_face_n(ref_face) ; face++ )
    if ( 2 != hits[face] )
      {
	problem = REF_TRUE;	
	printf(" hits %d\n",hits[face]);	  
      }

  free(hits);

  RSS( ref_face_free( ref_face ), "face free");

  return (problem?REF_FAILURE:REF_SUCCESS);
}

REF_STATUS ref_validation_cell_node( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT cell, node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL has_local;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      {
	has_local = REF_FALSE;
	for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
	  {
	    if ( ! ref_node_valid(ref_grid_node(ref_grid),nodes[node]))
	      {
		RSS( REF_FAILURE, "cell with invalid node" );
	      }
	    has_local = has_local || 
	      ( ref_mpi_id == ref_node_part(ref_node,nodes[node]) );
	  }
	if ( !has_local )
	  {
	    RSS( REF_FAILURE, "cell with all ghost nodes" );
	  }
      }

  ref_cell = ref_grid_tri(ref_grid);
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      {
	has_local = REF_FALSE;
	for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
	  {
	    if ( ! ref_node_valid(ref_grid_node(ref_grid),nodes[node]))
	      {
		RSS( REF_FAILURE, "cell with invalid node" );
	      }
	    has_local = has_local || 
	      ( ref_mpi_id == ref_node_part(ref_node,nodes[node]) );
	  }
	if ( !has_local )
	  {
	    RSS( REF_FAILURE, "cell with all ghost nodes" );
	  }
      }

  ref_cell = ref_grid_qua(ref_grid);
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
      {
	has_local = REF_FALSE;
	for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
	  {
	    if ( ! ref_node_valid(ref_grid_node(ref_grid),nodes[node]))
	      {
		RSS( REF_FAILURE, "cell with invalid node" );
	      }
	    has_local = has_local || 
	      ( ref_mpi_id == ref_node_part(ref_node,nodes[node]) );
	  }
	if ( !has_local )
	  {
	    RSS( REF_FAILURE, "cell with all ghost nodes" );
	  }
      }

  return REF_SUCCESS;
}

REF_STATUS ref_validation_cell_volume( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL volume;
  REF_DBL min_volume, max_volume;
  REF_INT cell_node;
  REF_DBL part_complexity, complexity, det;
  REF_BOOL first_volume;
  REF_INT part_nnode, total_nnode, node;

  ref_cell = ref_grid_tet(ref_grid);
  if (ref_grid_twod(ref_grid) ) ref_cell = ref_grid_tri(ref_grid);

  min_volume = 1.0e100; max_volume = -1.0e100;
  first_volume = REF_TRUE;
  complexity = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      if ( ref_grid_twod(ref_grid) )
	{
	  RSS( ref_node_tri_area( ref_node, nodes, &volume ), "area" );
	}
      else
	{
	  RSS( ref_node_tet_vol( ref_node, nodes, &volume ), "vol" );
	}
      RAS ( volume>0.0, "negative volume tet");
      if ( first_volume )
	{
	  min_volume = volume;
	  max_volume = volume;
	  first_volume = REF_FALSE;
	}
      else
	{
	  min_volume = MIN( min_volume, volume);
	  max_volume = MAX( max_volume, volume);
	}
      for ( cell_node = 0 ; 
	    cell_node < ref_cell_node_per( ref_cell ) ;
	    cell_node++ )
	{
	  RSS( ref_matrix_det_m( ref_node_metric_ptr(ref_node, 
						     nodes[cell_node]), 
				 &det),"det");
	  complexity += sqrt(det)*volume/((REF_DBL)ref_cell_node_per(ref_cell));
	}
    }

  volume = min_volume;
  RSS( ref_mpi_min( &volume, &min_volume, REF_DBL_TYPE ), "mpi min");
  volume = max_volume;
  RSS( ref_mpi_max( &volume, &max_volume, REF_DBL_TYPE ), "mpi max");

  part_nnode=0;
  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id == ref_node_part(ref_node,node) ) part_nnode++;
  RSS( ref_mpi_sum( &part_nnode, &total_nnode, 1, REF_INT_TYPE ), "int sum");

  part_complexity=complexity;
  RSS( ref_mpi_sum( &part_complexity, &complexity, 1, REF_DBL_TYPE ),"dbl sum");

  if ( ref_mpi_master )
    {
      if (ref_grid_twod(ref_grid) ) total_nnode = total_nnode / 2;
      printf("nnode %10d complexity %12.1f ratio %5.2f\nvolume range %e %e\n",
	     total_nnode, complexity, (REF_DBL)total_nnode/complexity,
	     max_volume, min_volume);
    }

  return REF_SUCCESS;
}
