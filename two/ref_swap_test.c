#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_mpi.h"

#include "ref_swap.h"

#include "ref_fixture.h"

int main( void )
{

  { /* leave single faces alone */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    REIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* leave different 2 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,50}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    REIS(REF_INVALID,ref_swap_remove_two_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* remove same 2 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4] = {0,3,1,3}, cell;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    RSS(ref_swap_remove_two_face_cell(ref_grid,0),"removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* leave different 3 faces alone */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    nodes[0]=0;nodes[1]=3;nodes[2]=1;nodes[3]=50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    nodes[0]=0;nodes[1]=2;nodes[2]=3;nodes[3]=50;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");
    
    REIS(REF_INVALID,ref_swap_remove_three_face_cell(ref_grid,0),"cell 0");
    REIS(1, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* remove same 3 faces */
    REF_GRID ref_grid;
    REF_INT nodes[4], cell;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    nodes[0]=0;nodes[1]=3;nodes[2]=1;nodes[3]=3;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    nodes[0]=0;nodes[1]=2;nodes[2]=3;nodes[3]=3;
    RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"other tri");

    RSS(ref_swap_remove_three_face_cell(ref_grid,0),"removal failed");
    REIS(0, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    REIS(3, ref_node_n(ref_grid_node(ref_grid)),"nodes");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
