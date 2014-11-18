
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_twod.h"

#include "ref_fixture.h"
#include  "ref_adj.h"
#include  "ref_node.h"
#include   "ref_list.h"
#include   "ref_sort.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include  "ref_grid.h"

int main( void )
{

  {   
    REF_GRID ref_grid;
    REF_INT node, opposite;
    
    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fix" );

    node = 0;
    RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid), 
				 node, &opposite ), "opp node" );
    REIS( 3, opposite, "wrong pair" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
