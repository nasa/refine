
#include <stdlib.h>
#include <stdio.h>

#include "ref_clump.h"

#include "ref_dict.h"
#include "ref_cell.h"

REF_STATUS ref_clump_around( REF_GRID ref_grid, REF_INT node )
{
  REF_DICT node_dict, tri_dict;
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT item, cell, cell_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&tri_dict),"create tris");
    
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_nodes(ref_cell,cell,nodes), "n");
      each_ref_cell_cell_node(ref_cell,cell_node)
	RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
    }

  RSS(ref_dict_free(tri_dict),"free tris");
  RSS(ref_dict_free(node_dict),"free nodes");

  return REF_SUCCESS;
}
