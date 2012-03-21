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

#include "ref_mpi.h"

#include "ref_split.h"

#include "ref_fixture.h"

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

  return 0;
}
