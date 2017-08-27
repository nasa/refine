
#include <stdlib.h>
#include <stdio.h>

#include "ref_layer.h"

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

REF_STATUS ref_layer_puff( REF_LAYER ref_layer, REF_GRID ref_grid )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT prism[REF_CELL_MAX_SIZE_PER];
  REF_INT node, local, global, i, nnode;
  
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
  printf(" layer ntri %d nnode %d\n",
	 ref_list_n(ref_layer_list(ref_layer)),
	 ref_node_n(ref_layer_node(ref_layer)));

  nnode = ref_node_n(ref_layer_node(ref_layer));
  for (local = 0;local<node;local++)
    {
      global = local+ref_node_n_global(ref_grid_node(ref_grid));
      RSS( ref_node_add(ref_layer_node(ref_layer), global, &node), "add");
      for (i=0;i<3;i++)
	ref_node_xyz(ref_layer_node(ref_layer), i, node) =
	  ref_node_xyz(ref_grid_node(ref_grid), i, local);
    }

  each_ref_list_item( ref_layer_list(ref_layer), item )
    {
      cell = ref_list_value( ref_layer_list(ref_layer), item );
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
	{
	  RSS( ref_node_local(ref_layer_node(ref_layer),
			      nodes[cell_node], &local), "local");
	  prism[cell_node] = local;
	  global = local+ref_node_n_global(ref_grid_node(ref_grid));
	  RSS( ref_node_local(ref_layer_node(ref_layer),
			      global, &local), "local");
	  prism[3+cell_node] = local;
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
