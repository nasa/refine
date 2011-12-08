
#include <stdlib.h>
#include <stdio.h>

#include "ref_quality.h"

#include "ref_cell.h"


REF_STATUS ref_quality_hex( REF_GRID ref_grid )
{
  REF_CELL ref_cell;
  REF_INT cell, cell_edge;
  REF_INT nodes[8];

  ref_cell = ref_grid_hex(ref_grid);

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    each_ref_cell_cell_edge( ref_cell, cell_edge )      
      {
	
      }


  return REF_SUCCESS;
}
