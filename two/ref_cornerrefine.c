#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_inflate.h"

#include "ref_math.h"

#include "ref_import.h"
#include "ref_export.h"

#include "ref_cell.h"
#include "ref_grid.h"
#include "ref_sort.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_matrix.h"
#include "ref_mpi.h"
#include "ref_dict.h"
#include "ref_list.h"

#include "ref_edge.h"

#include "ref_fixture.h"
#include "ref_metric.h"
#include "ref_gather.h"

#include "ref_subdiv.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_CELL ref_cell;
  REF_SUBDIV ref_subdiv;
  REF_BOOL valid_inputs;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node, item;
  REF_INT node_refinements, extra_at_node;

  valid_inputs = (4 == argc);  

  if ( !valid_inputs )
    {
      printf("usage: %s input_grid.extension output_grid.extension extra_at_node\n",argv[0]);
      return 1;
    }

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "in");
  RSS( ref_metric_unit_node( ref_grid_node(ref_grid) ), "met");

  ref_cell = ref_grid_pri(ref_grid);

  RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "init" );
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[1], nodes[2] ), "o0" );
      RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[2], nodes[0] ), "o1" );
      RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[0], nodes[1] ), "o2" );
    }
  RSS(ref_subdiv_split(ref_subdiv),"split");
  RSS(ref_subdiv_free(ref_subdiv),"free");

  extra_at_node = atoi(argv[3]);

  for (node_refinements = 0 ; 
       node_refinements < extra_at_node ;
       node_refinements++ )
    {
      RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "init" );
      node = 0;
      each_ref_cell_having_node( ref_cell, node, item, cell )
	{
	  RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes" );
	  RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[1], nodes[2] ), "o0" );
	  RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[2], nodes[0] ), "o1" );
	  RSS( ref_subdiv_mark_to_split( ref_subdiv, nodes[0], nodes[1] ), "o2" );
	}
      RSS(ref_subdiv_split(ref_subdiv),"split");
      RSS(ref_subdiv_free(ref_subdiv),"free");
    }

  RSS( ref_export_by_extension( ref_grid, argv[2] ), "out");

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
