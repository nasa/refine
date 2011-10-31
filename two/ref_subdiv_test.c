#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"

#include "ref_test.h"

int main( void )
{
  REF_GRID ref_grid;
  REF_CELL ref_cell;
  REF_INT nodes[6] = {0,1,2,3,4,5};
  REF_INT cell;

  TSS(ref_grid_create(&ref_grid),"create");

  ref_cell = ref_grid_pri(ref_grid);

  TSS(ref_cell_add(ref_cell,nodes,&cell),"add prism");

  TSS( ref_grid_free(ref_grid),"free" );

  return 0;
}
