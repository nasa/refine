
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

  { /* opposite prism node */
    REF_GRID ref_grid;
    REF_INT node, opposite;
    
    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fix" );

    node = 0;
    RSS( ref_twod_opposite_node( ref_grid_pri(ref_grid), 
				 node, &opposite ), "opp node" );
    REIS( 3, opposite, "wrong pair" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* opposite prism edge */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, node3;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");

    node0=0;node1=1;
    RSS(ref_twod_opposite_edge(ref_grid_pri(ref_grid),
			       node0,node1,&node2,&node3),"opp");
    REIS(3,node2,"n2");
    REIS(4,node3,"n3");

    node0=1;node1=0;
    RSS(ref_twod_opposite_edge(ref_grid_pri(ref_grid),
			       node0,node1,&node2,&node3),"opp");
    REIS(4,node2,"n2");
    REIS(3,node3,"n3");

    node0=4;node1=5;
    RSS(ref_twod_opposite_edge(ref_grid_pri(ref_grid),
			       node0,node1,&node2,&node3),"opp");
    REIS(1,node2,"n2");
    REIS(2,node3,"n3");

    node0=5;node1=4;
    RSS(ref_twod_opposite_edge(ref_grid_pri(ref_grid),
			       node0,node1,&node2,&node3),"opp");
    REIS(2,node2,"n2");
    REIS(1,node3,"n3");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
