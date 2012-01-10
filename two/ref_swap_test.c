#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_metric.h"

#include "ref_swap.h"

#include "ref_fixture.h"

int main( void )
{

  { /* leave single faces alone */
    REF_GRID ref_grid;

    TSS(ref_fixture_tet_grid(&ref_grid),"set up");

    TEIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");

    TSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
