#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_metric.h"
#include "ref_cell.h"
#include "ref_sort.h"

#include "ref_test.h"
#include "ref_axi.h"

REF_STATUS ref_quad_grid( REF_GRID *ref_grid_ptr, REF_DBL z0, REF_DBL z1 )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[5] = {0,1,2,3,10};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = z0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = z1;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = z1;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = z0;

  ref_node_n_global(ref_node) = ref_node_n(ref_node);

  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  return REF_SUCCESS;
}

int main( void )
{

  { /* collapse quad 0 */
    REF_GRID ref_grid;

    RSS( ref_quad_grid( &ref_grid, 0.0, 0.0 ), "quad fixture" );
    RSS( ref_axi_wedge( ref_grid ), "wedge");

    REIS( 2, ref_node_n(ref_grid_node(ref_grid)), "node n");
    REIS( 2, ref_node_n_global(ref_grid_node(ref_grid)), "node n global");

    REIS( 0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS( 0, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS( ref_grid_free( ref_grid ), "free" );
  }

  { /* collapse quad 1a */
    REF_GRID ref_grid;

    RSS( ref_quad_grid( &ref_grid, 0.0, 0.1 ), "quad fixture" );
    RSS( ref_axi_wedge( ref_grid ), "wedge");

    REIS( 3, ref_node_n(ref_grid_node(ref_grid)), "node n");
    REIS( 3, ref_node_n_global(ref_grid_node(ref_grid)), "node n global");

    REIS( 0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS( 1, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS( ref_grid_free( ref_grid ), "free" );
  }

  { /* collapse quad 1b */
    REF_GRID ref_grid;

    RSS( ref_quad_grid( &ref_grid, 0.1, 0.0 ), "quad fixture" );
    RSS( ref_axi_wedge( ref_grid ), "wedge");

    REIS( 3, ref_node_n(ref_grid_node(ref_grid)), "node n");
    REIS( 3, ref_node_n_global(ref_grid_node(ref_grid)), "node n global");

    REIS( 0, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS( 1, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS( ref_grid_free( ref_grid ), "free" );
  }

  { /* collapse quad 2 */
    REF_GRID ref_grid;

    RSS( ref_quad_grid( &ref_grid, 0.1, 0.1 ), "quad fixture" );
    RSS( ref_axi_wedge( ref_grid ), "wedge");

    REIS( 4, ref_node_n(ref_grid_node(ref_grid)), "node n");
    REIS( 4, ref_node_n_global(ref_grid_node(ref_grid)), "node n global");

    REIS( 1, ref_cell_n(ref_grid_qua(ref_grid)), "qua n");
    REIS( 0, ref_cell_n(ref_grid_tri(ref_grid)), "tri n");

    RSS( ref_grid_free( ref_grid ), "free" );
  }

  return 0;
}
