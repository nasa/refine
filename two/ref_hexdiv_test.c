#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_adj.h"
#include "ref_metric.h"
#include "ref_sort.h"

#include "ref_hexdiv.h"

#include "ref_test.h"

static REF_STATUS set_up_hex_for_hexdiv( REF_HEXDIV *ref_hexdiv_ptr )
{
  REF_GRID ref_grid;
  REF_INT nodes[8] = {0,1,2,3,4,5,6,7};
  REF_INT cell, node;

  TSS(ref_grid_create(&ref_grid),"create");

  TSS(ref_node_add( ref_grid_node(ref_grid), 0, &node),"n0");
  TSS(ref_node_add( ref_grid_node(ref_grid), 1, &node),"n1");
  TSS(ref_node_add( ref_grid_node(ref_grid), 2, &node),"n2");
  TSS(ref_node_add( ref_grid_node(ref_grid), 3, &node),"n3");
  TSS(ref_node_add( ref_grid_node(ref_grid), 4, &node),"n4");
  TSS(ref_node_add( ref_grid_node(ref_grid), 5, &node),"n5");
  TSS(ref_node_add( ref_grid_node(ref_grid), 5, &node),"n6");
  TSS(ref_node_add( ref_grid_node(ref_grid), 5, &node),"n7");

  TSS(ref_cell_add(ref_grid_hex(ref_grid),nodes,&cell),"add hex");
  TSS(ref_hexdiv_create(ref_hexdiv_ptr,ref_grid),"create");

  return REF_SUCCESS;
}

static REF_STATUS tear_down( REF_HEXDIV ref_hexdiv )
{
  REF_GRID ref_grid;

  ref_grid = ref_hexdiv_grid(ref_hexdiv);

  TSS(ref_hexdiv_free(ref_hexdiv),"free");

  TSS( ref_grid_free(ref_grid),"free" );

  return REF_SUCCESS;
}

int main( void )
{

  { /* mark and relax */
    REF_HEXDIV ref_hexdiv;
    TSS(set_up_hex_for_hexdiv(&ref_hexdiv),"set up");

    TSS( tear_down( ref_hexdiv ), "tear down");
  }

  return 0;
}
