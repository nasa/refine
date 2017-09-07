
#include <stdlib.h>
#include <stdio.h>

#include "ref_layer.h"

#include "ref_math.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

REF_STATUS ref_layer_create( REF_LAYER *ref_layer_ptr )
{
  REF_LAYER ref_layer;

  ref_malloc( *ref_layer_ptr, 1, REF_LAYER_STRUCT );

  ref_layer = *ref_layer_ptr;

  RSS(ref_list_create(&(ref_layer_list(ref_layer))),"create list");
  RSS(ref_grid_create(&(ref_layer_grid(ref_layer))),"create grid");

  return REF_SUCCESS;
}

REF_STATUS ref_layer_free( REF_LAYER ref_layer )
{
  if ( NULL == (void *)ref_layer ) return REF_NULL;

  ref_grid_free( ref_layer_grid(ref_layer) );
  ref_list_free( ref_layer_list(ref_layer) );
  ref_free( ref_layer );

  return REF_SUCCESS;
}

REF_STATUS ref_layer_attach( REF_LAYER ref_layer,
			     REF_GRID ref_grid, REF_INT faceid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  /* copy nodes into local copy that provides compact index */
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( faceid == nodes[ref_cell_node_per(ref_cell)] )
      RSS( ref_list_add( ref_layer_list(ref_layer), cell ), "parent" );
  
  return REF_SUCCESS;
}

REF_STATUS ref_layer_normal( REF_LAYER ref_layer, REF_GRID ref_grid,
			     REF_INT node, REF_DBL *norm )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT i, item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL contains;
  REF_DBL angle, total, triangle_norm[3];

  total = 0.0;
  norm[0] = 0.0;
  norm[1] = 0.0;
  norm[2] = 0.0;

  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_list_contains( ref_layer_list(ref_layer), cell,
			      &contains ), "in layer" );
      if ( ! contains ) 
	continue;
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "tri nodes");
      RSS( ref_node_tri_node_angle( ref_grid_node(ref_grid), nodes, node,
				    &angle ), "angle" );
      RSS( ref_node_tri_normal( ref_grid_node(ref_grid), nodes,
				triangle_norm ), "norm" );
      RSS( ref_math_normalize( triangle_norm ), "normalize tri norm" );
      total += angle;
      for(i=0;i<3;i++)
	norm[i] += angle*triangle_norm[i];
    }
  
  if ( !ref_math_divisible(norm[0],total) ||
       !ref_math_divisible(norm[1],total) ||
       !ref_math_divisible(norm[2],total) ) 
    return REF_DIV_ZERO;

  for(i=0;i<3;i++)
    norm[i] /= angle; 

  RSS( ref_math_normalize( norm ), "normalize average norm" );

  return REF_SUCCESS;
}

REF_STATUS ref_layer_puff( REF_LAYER ref_layer, REF_GRID ref_grid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE layer_node = ref_grid_node(ref_layer_grid(ref_layer));
  REF_CELL layer_prism = ref_grid_pri(ref_layer_grid(ref_layer));
  REF_CELL layer_edge = ref_grid_edg(ref_layer_grid(ref_layer));
  REF_INT item, cell, cell_node, cell_edge, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT prism[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT node, local, global, i, nnode;
  REF_DBL norm[3];

  /* first layer of nodes */
  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	{
	  RSS( ref_node_add(layer_node,
			    nodes[cell_node], &node), "add");
	  for (i=0;i<3;i++)
	    ref_node_xyz(layer_node, i, node) =
	      ref_node_xyz(ref_grid_node(ref_grid), i, nodes[cell_node]);
	}
    }
  nnode = ref_node_n(layer_node);

  /* second layer of nodes */
  for (local = 0;local<nnode;local++)
    {
      global = local+ref_node_n_global(ref_grid_node(ref_grid));
      RSS( ref_node_add(layer_node, global, &node), "add");
      RSS( ref_layer_normal(ref_layer,ref_grid,
			    ref_node_global(layer_node,
					    local), norm ), "normal");
      for (i=0;i<3;i++)
	ref_node_xyz(layer_node, i, node) =
	  0.1*norm[i] + ref_node_xyz(layer_node, i, local);
    }

  /* layer of prisms */
  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	{
	  RSS( ref_node_local(layer_node,
			      nodes[cell_node], &local), "local");
	  prism[cell_node] = local;
	  prism[3+cell_node] = local+nnode;
	}
      RSS(ref_cell_add(layer_prism, prism, &new_cell ), "add");
    }

  /* constrain faces */
  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_edge( ref_cell, cell_edge )
	{
	  REF_INT node0;
	  REF_INT node1;
	  REF_INT ncell, cell_list[2];
	  REF_BOOL contains0, contains1;
	  REF_INT edge_nodes[REF_CELL_MAX_SIZE_PER];
	  node0 = nodes[ref_cell_e2n_gen(ref_cell,0,cell_edge)];
	  node1 = nodes[ref_cell_e2n_gen(ref_cell,1,cell_edge)];
	  RSS( ref_cell_list_with2(ref_cell,node0,node1,
				   2, &ncell, cell_list ), "find with 2");
	  REIS(2,ncell, "expected two tri for tri side");
	  RSS( ref_list_contains( ref_layer_list(ref_layer), cell_list[0],
				  &contains0 ), "0 in layer" );
	  RSS( ref_list_contains( ref_layer_list(ref_layer), cell_list[1],
				  &contains1 ), "1 in layer" );
	  if ( contains0 && contains1 )
	    continue; /* tri side interior to layer */
	  if ( !contains0 && !contains1 )
	    THROW("tri side is not in layer");
	  RSS( ref_node_local(layer_node,
			      node0, &local), "local");
	  edge_nodes[0] = local+nnode;
	  RSS( ref_node_local(layer_node,
			      node1, &local), "local");
	  edge_nodes[1] = local+nnode;
	  if ( contains0 )
	    {
	      REIS(cell,cell_list[0], "cell should be in layer");
	      edge_nodes[2] = ref_cell_c2n(ref_cell,3,cell_list[1]);
	    }
	  if ( contains1 )
	    {
	      REIS(cell,cell_list[1], "cell should be in layer");
	      edge_nodes[2] = ref_cell_c2n(ref_cell,3,cell_list[0]);
	    }
	  RSS(ref_cell_add(layer_edge, edge_nodes, &new_cell ), "add");
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_layer_insert( REF_LAYER ref_layer, REF_GRID ref_grid )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE layer_node = ref_grid_node(ref_layer_grid(ref_layer));
  REF_INT nnode, node, local, global;
  REF_INT tet, i;
  REF_DBL bary[4];
  REF_INT zeros;
  REF_DBL zero_tol = 1.0e-10;
  nnode = ref_node_n(layer_node)/2; /* should persist in ref_layer */

  for ( node = 0 ; node < nnode ; node++ )
    {
      local = node+nnode;
      global = ref_node_global(layer_node,node);
      tet = ref_adj_first(ref_cell_adj(ref_cell),global);
      RSS( ref_grid_enclosing_tet( ref_grid,
				   ref_node_xyz_ptr(layer_node,local),
				   &tet, bary ), "enclosing tet" );
      zeros = 0;
      for (i=0;i<4;i++)
	if (ABS(bary[i]) < zero_tol)
	  zeros++;
      switch ( zeros )
	{
	case 2:
	  printf("%d bary %f %f %f %f\n",zeros,bary[0],bary[1],bary[2],bary[3]);
	  break;
	}
    }
  return REF_SUCCESS;
}
