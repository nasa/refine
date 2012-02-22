#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_metric.h"
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
  REF_DICT ref_dict;

  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT faceid, boundary_condition;
  REF_INT open_node0, open_node1;

  REF_SHARD ref_shard;
  REF_INT face_marks, hex_marks;

  if (argc<3) 
    {
      printf("usage: %s filename.extension filename.mapbc\n",argv[0]);
      return 0;
    }

  printf("reading %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"by extension");
  printf("complete.\n");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("reading %s\n",argv[2]);
  RSS(ref_import_mapbc( &ref_dict, argv[2] ),"mapbc");
  printf("complete.\n");

  RSS(ref_dict_inspect( ref_dict ), "inspection");

  printf("create shard.\n");
  RSS( ref_shard_create( &ref_shard, ref_grid ), "make shard");

  ref_node = ref_grid_node(ref_grid);
  ref_cell = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      faceid = nodes[4];
      RSS( ref_dict_value(ref_dict,faceid,&boundary_condition), "get bc");
      if ( boundary_condition == 4000 )
	{
	  RSS( ref_face_open_node( &ref_node_xyz(ref_node,0,nodes[0]),
				   &ref_node_xyz(ref_node,0,nodes[1]),
				   &ref_node_xyz(ref_node,0,nodes[2]),
				   &ref_node_xyz(ref_node,0,nodes[3]),
				   &open_node0 ), "find open" );
	  open_node1 = open_node0+2;
	  if ( open_node1 >= 4 ) open_node1 -= 4;
	  RSS( ref_shard_mark_to_split( ref_shard, 
					nodes[open_node0], nodes[open_node1]), 
	       "mark");
	}
    }

  RSS( ref_shard_mark_n( ref_shard, &face_marks, &hex_marks ), "count marks");
  printf("marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS( ref_shard_mark_relax( ref_shard ), "relax" );

  RSS( ref_shard_mark_n( ref_shard, &face_marks, &hex_marks ), "count marks");
  printf("relaxed marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS( ref_shard_split( ref_shard ), "split hex to prism" );

  RSS( ref_shard_free( ref_shard ), "free shard" );

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("ugrid.\n");
  RSS(ref_export_b8_ugrid(ref_grid, "prism.b8.ugrid"),"to ugrid");

  printf("free.\n");
  RSS(ref_grid_free(ref_grid),"free");

  printf("done.\n");
  return 0;
}
