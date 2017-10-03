
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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

  { /* pri tri stack */
    REF_GRID ref_grid;
    REF_INT cell, tri, pri;
    
    RSS( ref_fixture_pri_grid( &ref_grid ), "pri fix" );

    cell = 0;
    RSS( ref_twod_tri_pri_tri( ref_grid_tri(ref_grid), ref_grid_pri(ref_grid), 
			       cell, &pri, &tri ), "stack" );
    REIS( 0, pri, "only pri" );
    REIS( 1, tri, "other tri" );

    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
