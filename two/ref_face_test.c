#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"
#include "ref_metric.h"

#include "ref_face.h"

#include "ref_test.h"

int main( void )
{

  {  /* make faces shared by two elements */
    REF_FACE ref_face;
    REF_GRID ref_grid;
    REF_INT nodes[6];
    REF_INT cell;

    TSS(ref_grid_create(&ref_grid),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2;
    nodes[3] = 3; nodes[4] = 4; nodes[5] = 5;
    TSS(ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add pri");

    nodes[0] = 3; nodes[1] = 4; nodes[2] = 5; nodes[3] = 6;
    TSS(ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add tet");

    TSS(ref_face_create(&ref_face,ref_grid),"create");

    TEIS( 7, ref_face_n(ref_face), "check total faces");

    TSS(ref_face_free(ref_face),"face");
    TSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
