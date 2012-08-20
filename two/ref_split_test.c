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

#include "ref_mpi.h"

#include "ref_split.h"

#include "ref_fixture.h"
#include "ref_export.h"

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
    REF_NODE ref_node;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");


    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node( ref_node, node )
      {
	ref_node_metric(ref_node,0,node) = 1.0;
	ref_node_metric(ref_node,1,node) = 0.0;
	ref_node_metric(ref_node,2,node) = 0.0;
	ref_node_metric(ref_node,3,node) = 1.0;
	ref_node_metric(ref_node,4,node) = 0.0;
	ref_node_metric(ref_node,5,node) = 1.0;
      }

    RSS(ref_split_pass(ref_grid),"pass");

    REIS( 4, ref_node_n(ref_node), "nodes");
    REIS( 1, ref_cell_n(ref_grid_tet(ref_grid)), "tets");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* top small */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;

    RSS(ref_fixture_tet_grid(&ref_grid),"set up");

    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node( ref_node, node )
      {
	ref_node_metric(ref_node,0,node) = 1.0;
	ref_node_metric(ref_node,1,node) = 0.0;
	ref_node_metric(ref_node,2,node) = 0.0;
	ref_node_metric(ref_node,3,node) = 1.0;
	ref_node_metric(ref_node,4,node) = 0.0;
	ref_node_metric(ref_node,5,node) = 1.0;
      }

    ref_node_metric(ref_node,5,3) = 1/(0.25*0.25);

    RSS(ref_split_pass(ref_grid),"pass");

    REIS( 7, ref_node_n(ref_node), "nodes");
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

    RSS(ref_split_face(ref_grid,node0,node1,new_node0,
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

    RSS(ref_split_face(ref_grid,node0,node1,new_node0,
		                node2,node3,new_node1),"split");

    REIS(2, ref_cell_n(ref_grid_pri(ref_grid)),"tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  { /* twod prism, no split, close enough */
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node;

    RSS(ref_fixture_pri_grid(&ref_grid),"set up");

    ref_node = ref_grid_node(ref_grid);
    each_ref_node_valid_node( ref_node, node )
      {
	ref_node_metric(ref_node,0,node) = 1.0;
	ref_node_metric(ref_node,1,node) = 0.0;
	ref_node_metric(ref_node,2,node) = 0.0;
	ref_node_metric(ref_node,3,node) = 1.0;
	ref_node_metric(ref_node,4,node) = 0.0;
	ref_node_metric(ref_node,5,node) = 1.0;
      }

    RSS(ref_split_twod_pass(ref_grid),"pass");

    REIS( 6, ref_node_n(ref_node), "nodes");
    REIS( 1, ref_cell_n(ref_grid_pri(ref_grid)), "tets");

    RSS( ref_grid_free( ref_grid ), "free grid");
  }

  return 0;
}
