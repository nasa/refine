#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"

#include "ref_cell.h"

#include "ref_edge.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_dict.h"

#include "ref_quality.h"
#include "ref_shard.h"
#include "ref_subdiv.h"

#include "ref_validation.h"

#include "ref_math.h"
#include "ref_swap.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_DICT orig_dict, o2n_dict;

  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT bc;

  REF_INT bc_to_lump, lump_to;
  REF_INT key_index, old_face, new_face, face;

  if (argc<4) 
    {
      printf("usage: %s filename.extension filename.mapbc bc_to_lump\n",
	     argv[0]);
      return 0;
    }

  printf("reading %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"by extension");
  printf("complete.\n");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("reading %s\n",argv[2]);
  RSS(ref_import_mapbc( &orig_dict, argv[2] ),"mapbc");
  printf("complete.\n");

  RSS(ref_dict_inspect( orig_dict ), "inspection");

  RSS( ref_dict_create( &o2n_dict ), "o2n_dict" );

  bc_to_lump = atoi(argv[3]);
  printf("lumping bc %d\n",bc_to_lump);

  lump_to = REF_EMPTY;
  new_face = 0;
  each_ref_dict_key( orig_dict, key_index, old_face )
    {
      RSS( ref_dict_value(orig_dict, old_face, &bc),"val");
      if ( bc == bc_to_lump )
	{
	  if ( REF_EMPTY == lump_to ) 
	    {
	      new_face++;
	      lump_to = new_face;
	      printf("%d %d\n",lump_to,bc);
	    }
	  face = lump_to;
	}
      else
	{
	  new_face++;
	  face = new_face;
	  printf("%d %d\n",face,bc);
	}
      RSS( ref_dict_store( o2n_dict, old_face, face ), "store o2n");
    }

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_dict_value(o2n_dict, 
			  nodes[ref_cell_node_per(ref_cell)], &face),"val");
      nodes[ref_cell_node_per(ref_cell)] = face;
      RSS( ref_cell_replace_whole( ref_cell, cell, nodes), "replace" );
    }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_dict_value(o2n_dict, 
			  nodes[ref_cell_node_per(ref_cell)], &face),"val");
      nodes[ref_cell_node_per(ref_cell)] = face;
      RSS( ref_cell_replace_whole( ref_cell, cell, nodes), "replace" );
    }

  printf("ugrid.\n");
  RSS(ref_export_b8_ugrid(ref_grid, "lumped.b8.ugrid"),"to ugrid");

  printf("free.\n");
  RSS(ref_dict_free(orig_dict),"free");
  RSS(ref_dict_free(o2n_dict),"free");
  RSS(ref_grid_free(ref_grid),"free");

  printf("done.\n");
  return 0;
}
