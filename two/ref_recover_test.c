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
#include "ref_fixture.h"

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

  SKIP_BLOCK("implement")  
  { /* insert a node */
    REF_GRID ref_grid;
    REF_RECOVER ref_recover;
    REF_DBL xz[2];
    REF_INT node;

    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fixture" );
    ref_grid_twod(ref_grid) = REF_TRUE;

    RSS(ref_recover_create(&ref_recover,ref_grid),"create");

    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    xz[0] = 0.3; xz[1] = 0.3;
    RSS(ref_recover_insert_twod(ref_recover,xz,&node),"create");
 
    REIS(8, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS(ref_recover_free(ref_recover),"free");
    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
