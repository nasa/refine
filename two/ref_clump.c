
#include <stdlib.h>
#include <stdio.h>

#include "ref_clump.h"

#include "ref_dict.h"
#include "ref_cell.h"

REF_STATUS ref_clump_around( REF_GRID ref_grid, REF_INT node, char *filename )
{
  REF_DICT node_dict, tri_dict;
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT local;

  FILE *f;

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&tri_dict),"create tris");
    
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      RSS( ref_dict_store( tri_dict, cell, 0 ), "store");
      each_ref_cell_cell_node(ref_cell,cell_node)
	RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
    }

  f = fopen(filename,"w");
  if (NULL == (void *)f) printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  fprintf(f, "title=\"tecplot refine clump file\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(f,
	  "zone t=clump, nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
	  ref_dict_n(node_dict), ref_dict_n(tri_dict), 
	  "point", "fequadrilateral" );

  for ( item = 0; item < ref_dict_n(node_dict); item++ )
    {
      node = ref_dict_key(node_dict,item);
      fprintf(f, " %.16e %.16e %.16e\n",
	      ref_node_xyz(ref_grid_node(ref_grid),0,node),
	      ref_node_xyz(ref_grid_node(ref_grid),1,node),
	      ref_node_xyz(ref_grid_node(ref_grid),2,node) );
    }

  for ( item = 0; item < ref_dict_n(tri_dict); item++ )
    {
      cell = ref_dict_key(tri_dict,item);
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      each_ref_cell_cell_node(ref_cell,cell_node)
	{
	  RSS( ref_dict_location( node_dict, 
				  nodes[cell_node], &local), "ret");
	  fprintf(f," %d",local + 1);
	}
      RSS( ref_dict_location( node_dict, 
			      nodes[0], &local), "ret");
      fprintf(f," %d",local + 1);
      fprintf(f,"\n");
    }

  fclose(f);

  RSS(ref_dict_free(tri_dict),"free tris");
  RSS(ref_dict_free(node_dict),"free nodes");

  return REF_SUCCESS;
}
