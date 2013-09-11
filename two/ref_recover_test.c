#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_recover.h"
#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"

int main( void )
{

  { /* init */
    REF_GRID ref_grid;
    REF_RECOVER ref_recover;
    REIS(REF_NULL, ref_recover_free(NULL),"dont free NULL");
    RSS(ref_grid_create(&ref_grid),"create");
    RSS(ref_recover_create(&ref_recover,ref_grid),"create");
    REIS( 0, ref_recover_n(ref_recover), "init no recover");
    RSS(ref_recover_free(ref_recover),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
