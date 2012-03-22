
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_collapse.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_mpi.h"
#include "ref_sort.h"
#include "ref_malloc.h"

#define MAX_CELL_COLLAPSE (100)

REF_STATUS ref_collapse_edge( REF_GRID ref_grid, 
			      REF_INT node0, REF_INT node1 )
{
  REF_CELL ref_cell;
  REF_INT cell;
  REF_INT ncell, cell_in_list;
  REF_INT cell_to_collapse[MAX_CELL_COLLAPSE];

  ref_cell = ref_grid_tet(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_COLLAPSE, &ncell, cell_to_collapse ),"list");

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_collapse[cell_in_list];
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );
    }
  RSS( ref_cell_replace_node( ref_cell, node1, node0 ), "replace node" );

  ref_cell = ref_grid_tri(ref_grid);
  RSS( ref_cell_list_with(ref_cell,node0,node1,
			  MAX_CELL_COLLAPSE, &ncell, cell_to_collapse ),"list");

  for ( cell_in_list = 0; cell_in_list < ncell ; cell_in_list++ )
    {
      cell = cell_to_collapse[cell_in_list];
      RSS( ref_cell_remove( ref_cell, cell ), "remove" );
    }
  RSS( ref_cell_replace_node( ref_cell, node1, node0 ), "replace node" );

  RSS( ref_node_remove(ref_grid_node(ref_grid),node1), "rm");

  return REF_SUCCESS;
}

