
#include <stdlib.h>
#include <stdio.h>

#include "ref_fixture.h"

REF_STATUS ref_fixture_tet_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[4] = {0,1,2,3};
  REF_INT cell, node;

  RSS(ref_grid_create(ref_grid_ptr),"create");
  ref_grid =  *ref_grid_ptr;

  ref_node = ref_grid_node(ref_grid);

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");

  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}
