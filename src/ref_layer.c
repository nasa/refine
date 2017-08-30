
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
  RSS(ref_node_create(&(ref_layer_node(ref_layer))),"create node");
  RSS(ref_cell_create(&(ref_layer_cell(ref_layer)),6,REF_FALSE),"create cell");

  return REF_SUCCESS;
}

REF_STATUS ref_layer_free( REF_LAYER ref_layer )
{
  if ( NULL == (void *)ref_layer ) return REF_NULL;

  ref_cell_free( ref_layer_cell(ref_layer) );
  ref_node_free( ref_layer_node(ref_layer) );
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
  REF_INT item, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT prism[REF_CELL_MAX_SIZE_PER];
  REF_INT node, local, global, i, nnode;
  REF_DBL norm[3];

  /* first layer of nodes */
  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	{
	  RSS( ref_node_add(ref_layer_node(ref_layer),
			    nodes[cell_node], &node), "add");
	  for (i=0;i<3;i++)
	    ref_node_xyz(ref_layer_node(ref_layer), i, node) =
	      ref_node_xyz(ref_grid_node(ref_grid), i, nodes[cell_node]);
	}
    }
  nnode = ref_node_n(ref_layer_node(ref_layer));

  /* second layer of nodes */
  for (local = 0;local<nnode;local++)
    {
      global = local+ref_node_n_global(ref_grid_node(ref_grid));
      RSS( ref_node_add(ref_layer_node(ref_layer), global, &node), "add");
      RSS( ref_layer_normal(ref_layer,ref_grid,
			    ref_node_global(ref_layer_node(ref_layer),
					    local), norm ), "normal");
      for (i=0;i<3;i++)
	ref_node_xyz(ref_layer_node(ref_layer), i, node) =
	  0.1*norm[i] + ref_node_xyz(ref_layer_node(ref_layer), i, local);
    }

  /* layer of prisms */
  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	{
	  RSS( ref_node_local(ref_layer_node(ref_layer),
			      nodes[cell_node], &local), "local");
	  prism[cell_node] = local;
	  prism[3+cell_node] = local+nnode;
	}
      RSS(ref_cell_add(ref_layer_cell(ref_layer), prism, &item ), "add");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_layer_tec( REF_LAYER ref_layer, const char *filename )
{
  REF_NODE ref_node = ref_layer_node(ref_layer);
  REF_CELL ref_cell = ref_layer_cell(ref_layer);
  FILE *file;
  REF_INT node, i, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT brick[REF_CELL_MAX_SIZE_PER];
  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"layer\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");
  
  fprintf(file,
	  "zone t=layer, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  ref_node_n(ref_node), ref_cell_n(ref_cell), "point", "febrick" );

  for ( node = 0; node < ref_node_n(ref_node); node++ )
    fprintf(file, " %.16e %.16e %.16e\n",
	    ref_node_xyz(ref_node,0,node),
	    ref_node_xyz(ref_node,1,node),
	    ref_node_xyz(ref_node,2,node) ) ;

#define TEC_BRICK_PRI(brick,nodes)					\
  {									\
    brick[0] = nodes[0]; brick[1] = nodes[1];				\
    brick[2] = nodes[2]; brick[3] = nodes[2];				\
    brick[4] = nodes[3]; brick[5] = nodes[4];				\
    brick[6] = nodes[5]; brick[7] = nodes[5];				\
  }
  
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    {
      TEC_BRICK_PRI(brick,nodes);
      for ( i = 0; i < 8; i++ )
	{
	  fprintf(file," %d",brick[i] + 1);
	}
      fprintf(file,"\n");
    }
  fclose(file);
  return REF_SUCCESS;
}
