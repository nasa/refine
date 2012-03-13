#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"
#include "ref_face.h"

#include "ref_subdiv.h"

#include "ref_fixture.h"
#include "ref_export.h"
#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_validation.h"

#include "ref_part.h"
#include "ref_migrate.h"
#include "ref_gather.h"

static REF_STATUS set_up_tet_for_subdiv( REF_SUBDIV *ref_subdiv_ptr )
{
  REF_GRID ref_grid;

  RSS(ref_fixture_tet_grid( &ref_grid ), "tet");

  RSS(ref_subdiv_create(ref_subdiv_ptr,ref_grid),"create");

  return REF_SUCCESS;
}

static REF_STATUS set_up_pyramid_for_subdiv( REF_SUBDIV *ref_subdiv_ptr )
{
  REF_GRID ref_grid;

  RSS(ref_fixture_pyr_grid( &ref_grid ), "pri");

  RSS(ref_subdiv_create(ref_subdiv_ptr,ref_grid),"create");

  return REF_SUCCESS;
}

static REF_STATUS set_up_prism_for_subdiv( REF_SUBDIV *ref_subdiv_ptr )
{
  REF_GRID ref_grid;

  RSS(ref_fixture_pri_grid( &ref_grid ), "pri");

  RSS(ref_subdiv_create(ref_subdiv_ptr,ref_grid),"create");

  return REF_SUCCESS;
}

static REF_STATUS tear_down( REF_SUBDIV ref_subdiv )
{
  REF_GRID ref_grid;

  ref_grid = ref_subdiv_grid(ref_subdiv);

  RSS(ref_subdiv_free(ref_subdiv),"free");

  RSS( ref_grid_free(ref_grid),"free" );

  return REF_SUCCESS;
}

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( 1 < argc )
    {
      REF_SUBDIV ref_subdiv;
      REF_GRID ref_grid;
      REF_NODE ref_node;
      REF_EDGE ref_edge;
      REF_INT edge;
      
      RSS(ref_part_b8_ugrid( &ref_grid, argv[1] ), "import" );
      RSS(ref_migrate_new_part(ref_grid),"new part");
      RSS( ref_migrate_shufflin( ref_grid ), "shufflin");

      RSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

      RSS(ref_export_tec_part(ref_grid,"ref_subdiv_orig"),"orig part");

      ref_node = ref_grid_node(ref_grid);
      ref_edge = ref_subdiv_edge(ref_subdiv);

      for ( edge=0; edge<ref_edge_n(ref_edge);edge++ )
	if ( ref_node_xyz(ref_node, 2, ref_edge_e2n(ref_edge,edge,0))<0.001 &&
	     ref_node_xyz(ref_node, 2, ref_edge_e2n(ref_edge,edge,1))<0.001 )
	  RSS(ref_subdiv_mark_to_split(ref_subdiv,
				       ref_edge_e2n(ref_edge,edge,0),
				       ref_edge_e2n(ref_edge,edge,1)
				       ),"mark edge");

      RSS(ref_subdiv_split(ref_subdiv),"split");

      RSS(ref_export_tec_part(ref_grid,"ref_subdiv_splt"),"split part");

      RSS(ref_gather_b8_ugrid(ref_grid,"ref_subdiv_test.b8.ugrid"),"gather");

      RSS( tear_down( ref_subdiv ), "tear down");
    }

  { /* split stack in two */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    REF_NODE ref_node;
    REF_INT node0, node1;

    RSS( ref_fixture_pri_stack_grid( &ref_grid ), "stack" );
    ref_node = ref_grid_node(ref_grid);
    RSS(ref_subdiv_create(&ref_subdiv,ref_grid),"create");

    if ( 1 < ref_mpi_n )
      RSS(ref_export_tec_part(ref_grid,"stack_orig"),"stack part");

    REIS( 12, ref_node_n_global(ref_node), "start with 12" );

    if ( REF_SUCCESS == ref_node_local(ref_node,0,&node0) &&
	 REF_SUCCESS == ref_node_local(ref_node,1,&node1) ) 
      if ( ref_mpi_id == ref_node_part(ref_node,node0) )
	RSS(ref_subdiv_mark_to_split(ref_subdiv,node0,node1),"mark edge");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    if ( 1 < ref_mpi_n )
      RSS(ref_export_tec_part(ref_grid,"stack_split"),"stack part");

    REIS( 16, ref_node_n_global(ref_node), "where my nodes?" );

    RSS(ref_validation_cell_node(ref_grid),"validate");
    RSS(ref_validation_hanging_node(ref_grid),"validate");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  RSS( ref_mpi_stop( ), "stop" );

  { /* mark and relax prism*/
    REF_SUBDIV ref_subdiv;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");

    RES(0,ref_subdiv_mark(ref_subdiv,0),"init mark");
    RES(0,ref_subdiv_mark(ref_subdiv,6),"init mark");

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    RES(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    RES(0,ref_subdiv_mark(ref_subdiv,6),"not been marked");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    RES(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    RES(1,ref_subdiv_mark(ref_subdiv,6),"been relaxed");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* mark and relax prism*/
    REF_SUBDIV ref_subdiv;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1-2");

    RES(0,ref_subdiv_mark(ref_subdiv,1),"no yet");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    RES(1,ref_subdiv_mark(ref_subdiv,1),"yet");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax tet */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1-2");

    REIS(1,ref_subdiv_mark(ref_subdiv,0),"yet");
    REIS(1,ref_subdiv_mark(ref_subdiv,3),"yet");

    REIS(0,ref_subdiv_mark(ref_subdiv,1),"no yet");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,1),"yet");

    REIS(0,ref_subdiv_mark(ref_subdiv,2),"no yet");
    REIS(0,ref_subdiv_mark(ref_subdiv,4),"no yet");
    REIS(0,ref_subdiv_mark(ref_subdiv,5),"no yet");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax tet oppisite edges */
    REF_SUBDIV ref_subdiv;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,2,3),"mark edge 5");

    REIS(1,ref_subdiv_mark(ref_subdiv,0),"yet");
    REIS(1,ref_subdiv_mark(ref_subdiv,5),"yet");

    REIS(0,ref_subdiv_mark(ref_subdiv,1),"no yet");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,1),"yet");

    REIS(1,ref_subdiv_mark(ref_subdiv,2),"yet");
    REIS(1,ref_subdiv_mark(ref_subdiv,4),"yet");
    REIS(1,ref_subdiv_mark(ref_subdiv,5),"yet");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* new nodes */
    REF_SUBDIV ref_subdiv;
    REF_NODE ref_node;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");
    ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,3,4),"mark edge 6 3-4");

    REIS(6, ref_node_n(ref_node), "6 node prism" );
    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(8, ref_node_n(ref_node), "two new nodes" );

    REIS(6, ref_subdiv_node(ref_subdiv,0), "new 6" );
    REIS(7, ref_subdiv_node(ref_subdiv,6), "new 6" );

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split prism in two */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(2, ref_cell_n(ref_grid_pri(ref_grid)),"two");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split prism in two with bcs */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(2, ref_cell_n(ref_grid_pri(ref_grid)),"two pri");

    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"two qua");

    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"four tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split prism in four with bcs */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_prism_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1-2");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,2,0),"mark edge 2-0");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(4, ref_cell_n(ref_grid_pri(ref_grid)),"two pri");

    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"two qua");

    REIS(8, ref_cell_n(ref_grid_tri(ref_grid)),"four tri");

    /*    ref_export_tec(ref_grid,"pri4.tec");  */

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in two, map 1 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)),"two tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"two tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in two, map 4 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,3),"mark edge 0-3");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(2, ref_cell_n(ref_grid_tet(ref_grid)),"two tet");
    REIS(1, ref_cell_n(ref_grid_tri(ref_grid)),"still one tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in 4 around node 3 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0-1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1-2");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)),"four tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"four tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in 4 around node 0 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,3,2),"mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)),"four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in 4 around node 1 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,3,2),"mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,2),"mark edge");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)),"four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in 4 around node 2 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,3,1),"mark edge");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(4, ref_cell_n(ref_grid_tet(ref_grid)),"four tet");
    REIS(2, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* split tet in 8 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_tet_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,2,3),"mark edge 5");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(8, ref_cell_n(ref_grid_tet(ref_grid)),"eight tet");
    REIS(4, ref_cell_n(ref_grid_tri(ref_grid)),"tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* unsplit pyramid */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)),"tri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* unsplit pyramid */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)),"pri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax and split pyramid in two 0-1 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 0");

    REIS(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    REIS(0,ref_subdiv_mark(ref_subdiv,7),"been marked");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    REIS(1,ref_subdiv_mark(ref_subdiv,7),"been marked");

    REIS(2, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(0, ref_cell_n(ref_grid_pri(ref_grid)),"pri");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax and split pyramid in to pyr and pri 1-2 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,1,2),"mark edge 1");

    REIS(1,ref_subdiv_mark(ref_subdiv,3),"been marked");
    REIS(0,ref_subdiv_mark(ref_subdiv,6),"been marked");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,3),"been marked");
    REIS(1,ref_subdiv_mark(ref_subdiv,6),"been marked");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(1, ref_cell_n(ref_grid_pri(ref_grid)),"pri");
    REIS(6, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax and split pyramid in to pyr and pri 0-2 */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,2),"mark edge 1");

    REIS(1,ref_subdiv_mark(ref_subdiv,1),"been marked");
    REIS(0,ref_subdiv_mark(ref_subdiv,5),"been marked");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,1),"been marked");
    REIS(1,ref_subdiv_mark(ref_subdiv,5),"been marked");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(1, ref_cell_n(ref_grid_pri(ref_grid)),"pri");
    REIS(6, ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(2, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  { /* relax and split pyramid in to pyr and 3 pri rpomoting two edges */
    REF_SUBDIV ref_subdiv;
    REF_GRID ref_grid;
    RSS(set_up_pyramid_for_subdiv(&ref_subdiv),"set up");
    ref_grid = ref_subdiv_grid(ref_subdiv);

    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,1),"mark edge 1");
    RSS(ref_subdiv_mark_to_split(ref_subdiv,0,2),"mark edge 1");

    RSS(ref_subdiv_split(ref_subdiv),"split");

    REIS(1,ref_subdiv_mark(ref_subdiv,0),"been marked");
    REIS(1,ref_subdiv_mark(ref_subdiv,1),"been marked");
    REIS(1,ref_subdiv_mark(ref_subdiv,3),"been marked");

    REIS(1, ref_cell_n(ref_grid_pyr(ref_grid)),"pyr");
    REIS(3, ref_cell_n(ref_grid_pri(ref_grid)),"pri");
    REIS(10,ref_cell_n(ref_grid_tri(ref_grid)),"tri");
    REIS(4, ref_cell_n(ref_grid_qua(ref_grid)),"qua");

    RSS( tear_down( ref_subdiv ), "tear down");
  }

  return 0;
}
