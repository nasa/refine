
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

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_matrix.h"
#include "ref_math.h"

#include "ref_mpi.h"

#include "ref_split.h"
#include "ref_adapt.h"
#include "ref_collapse.h"
#include "ref_smooth.h"
#include  "ref_twod.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_metric.h"

#include "ref_gather.h"
#include "ref_geom.h"
#include "ref_clump.h"

int main( void )
{

  { /* split tet in two */
    REF_GRID ref_grid;
    REF_INT node0, node1, new_node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    node0 = 0; node1 = 3;

    RSS( ref_node_add(ref_grid_node(ref_grid),4,&new_node), "new");
    RSS(ref_split_edge(ref_grid,node0,node1,new_node),"split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)),"tri");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split tet and tri in two */
    REF_GRID ref_grid;
    REF_INT node0, node1, new_node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    node0 = 0; node1 = 1;

    RSS( ref_node_add(ref_grid_node(ref_grid),4,&new_node), "new");
    RSS(ref_split_edge(ref_grid,node0,node1,new_node),"split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)),"edg");
    
    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split tet, tri, edg in two */
    REF_GRID ref_grid;
    REF_INT node0, node1, new_node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    node0 = 1; node1 = 2;

    RSS( ref_node_add(ref_grid_node(ref_grid),4,&new_node), "new");
    RSS(ref_split_edge(ref_grid,node0,node1,new_node),"split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(2, ref_cell_n(ref_grid_edg(ref_grid)),"edg");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split tet in three */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, new_node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    node0 = 1; node1 = 2; node2 = 3;

    RSS( ref_node_add(ref_grid_node(ref_grid),4,&new_node), "new");
    RSS(ref_split_face(ref_grid,node0,node1,node2,new_node),"split");

    REIS(3, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)),"tri");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split tet and tri in three */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, new_node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    node0 = 0; node1 = 1; node2 = 2;

    RSS( ref_node_add(ref_grid_node(ref_grid),4,&new_node), "new");
    RSS(ref_split_face(ref_grid,node0,node1,node2,new_node),"split");

    REIS(3, ref_cell_n(ref_grid_tet(ref_grid)),"tet");
    REIS(3, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(1, ref_cell_n(ref_grid_edg(ref_grid)),"tri");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split two tet in six */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, new_node;

    RSS(ref_fixture_tet2_grid(&ref_grid),"set up");
    node0 = 1; node1 = 2; node2 = 3;

    RSS( ref_node_add(ref_grid_node(ref_grid),5,&new_node), "new");
    RSS(ref_split_face(ref_grid,node0,node1,node2,new_node),"split");

    REIS(6, ref_cell_n(ref_grid_tet(ref_grid)),"tet");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split tet allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    node0 = 0; node1 = 1;
    RSS(ref_split_edge_mixed(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_TRUE,allowed,"pure tet allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split near mixed allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_tet_cap_grid(&ref_grid),"set up");

    node0 = 5; node1 = 6;
    RSS(ref_split_edge_mixed(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_TRUE,allowed,"tet near mixed allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split mixed allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_tet_cap_grid(&ref_grid),"set up");

    node0 = 3; node1 = 4;
    RSS(ref_split_edge_mixed(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_FALSE,allowed,"mixed split allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split local allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    node0 = 0; node1 = 1;
    RSS(ref_split_edge_local_tets(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_TRUE,allowed,"local split allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split partition border allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    ref_node_part(ref_grid_node(ref_grid),2) =
      ref_node_part(ref_grid_node(ref_grid),2) + 1;

    node0 = 0; node1 = 1;
    RSS(ref_split_edge_local_tets(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_FALSE,allowed,"ghost split allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* no split, close enough */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid)), "id metric" );

    RSS(ref_split_pass(ref_grid),"pass");

    REIS( 4, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS( 1, ref_cell_n(ref_grid_tet(ref_grid)), "tets");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* top small */
    REF_GRID ref_grid;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid)), "id metric" );

    ref_node_metric(ref_grid_node(ref_grid),5,3) = 1.0/(0.25*0.25);

    RSS(ref_split_pass(ref_grid),"pass");

    REIS( 7, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS( 4, ref_cell_n(ref_grid_tet(ref_grid)), "tets");

    /* ref_export_by_extension(ref_grid,"ref_split_test.tec"); */

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split twod prism in two */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, node3, new_node0, new_node1;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");
    node0 = 1; node1 = 2;
    node2 = 4; node3 = 5;

    RSS( ref_node_add(ref_grid_node(ref_grid),6,&new_node0), "new");
    RSS( ref_node_add(ref_grid_node(ref_grid),7,&new_node1), "new");

    RSS(ref_split_twod_edge(ref_grid,node0,node1,new_node0,
			    node2,node3,new_node1),"split");

    REIS(2, ref_cell_n(ref_grid_pri(ref_grid)),"tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(1, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split twod prism in two with quad */
    REF_GRID ref_grid;
    REF_INT node0, node1, node2, node3, new_node0, new_node1;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");
    node0 = 0; node1 = 1;
    node2 = 3; node3 = 4;

    RSS( ref_node_add(ref_grid_node(ref_grid),6,&new_node0), "new");
    RSS( ref_node_add(ref_grid_node(ref_grid),7,&new_node1), "new");

    RSS(ref_split_twod_edge(ref_grid,node0,node1,new_node0,
			    node2,node3,new_node1),"split");

    REIS(2, ref_cell_n(ref_grid_pri(ref_grid)),"tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* twod prism, no split, close enough */
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid)), "id metric" );

    RSS(ref_split_twod_pass(ref_grid),"pass");

    REIS( 6, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS( 1, ref_cell_n(ref_grid_pri(ref_grid)), "tets");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* twod prism, splits */
    REF_GRID ref_grid;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");
    RSS( ref_metric_unit_node( ref_grid_node(ref_grid)), "id metric" );

    ref_node_metric(ref_grid_node(ref_grid),5,1) = 1.0/(0.25*0.25);
    ref_node_metric(ref_grid_node(ref_grid),5,4) = 1.0/(0.25*0.25);

    RSS(ref_split_twod_pass(ref_grid),"pass");

    REIS( 10, ref_node_n(ref_grid_node(ref_grid)), "nodes");
    REIS( 3, ref_cell_n(ref_grid_pri(ref_grid)), "tets");
    REIS( 6, ref_cell_n(ref_grid_tri(ref_grid)), "tri");
    REIS( 2, ref_cell_n(ref_grid_qua(ref_grid)), "qua");

    /*
    ref_export_by_extension( ref_grid, "splitpri.tec" );
    */

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split local allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");

    node0 = 0; node1 = 1;
    RSS(ref_split_edge_local_prisms(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_TRUE,allowed,"local split allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* split border allowed? */
    REF_GRID ref_grid;
    REF_INT node0, node1;
    REF_BOOL allowed;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");

    ref_node_part(ref_grid_node(ref_grid),2) =
      ref_node_part(ref_grid_node(ref_grid),2) + 1;

    node0 = 0; node1 = 1;
    RSS(ref_split_edge_local_prisms(ref_grid,node0,node1,&allowed),"split");

    REIS(REF_FALSE,allowed,"local split allowed?");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
