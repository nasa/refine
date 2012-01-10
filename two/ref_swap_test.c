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
    TEIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    TSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* leave different faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,50}, cell;

    TSS(ref_fixture_tet_grid(&ref_grid),"set up");
    TSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    TEIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");
    TEIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    TSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* same faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,3}, cell;

    TSS(ref_fixture_tet_grid(&ref_grid),"set up");
    TSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    TSS(ref_swap_remove_two_face_cell(ref_grid,0),"removal failed");
    TEIS(0, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    TSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
