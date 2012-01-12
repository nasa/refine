
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

REF_STATUS ref_fixture_pyr_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[5] = {0,3,4,1,2};
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
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.5;
  ref_node_xyz(ref_node,1,node) = 0.5;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_cell_add(ref_grid_pyr(ref_grid),nodes,&cell),"add pyr");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add qua");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 4; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 1; nodes[1] = 2; nodes[2] = 4; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 2; nodes[1] = 3; nodes[2] = 4; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 3; nodes[1] = 0; nodes[2] = 4; nodes[3] = 20;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  RSS(ref_node_add(ref_node,5,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 2.0;

  RSS(ref_node_add(ref_node,6,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 2.0;

  RSS(ref_node_add(ref_node,7,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 2.0;

  RSS(ref_node_add(ref_node,8,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 3.0;

  /*
  8    
    7
  5    6
  */
  nodes[0] = 5; nodes[1] = 6; nodes[2] = 7; nodes[3] = 8;
  RSS(ref_cell_add(ref_grid_tet(ref_grid),nodes,&cell),"add tet");

  nodes[0] = 6; nodes[1] = 8; nodes[2] = 7; nodes[3] = 30;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 5; nodes[1] = 7; nodes[2] = 8; nodes[3] = 30;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 5; nodes[1] = 8; nodes[2] = 6; nodes[3] = 30;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 5; nodes[1] = 6; nodes[2] = 7; nodes[3] = 30;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_pri_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[6] = {0,1,2,3,4,5};
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

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,4,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add prism");

  nodes[0] = 0; nodes[1] = 3; nodes[2] = 4; nodes[3] = 1; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 3; nodes[1] = 5; nodes[2] = 4; nodes[3] = 100;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 101;
  RSS(ref_cell_add(ref_grid_tri(ref_grid),nodes,&cell),"add tri");

  return REF_SUCCESS;
}

REF_STATUS ref_fixture_hex_grid( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_INT nodes[8] = {0,1,2,3,4,5,6,7};
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
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 0.0;

  RSS(ref_node_add(ref_node,0,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,1,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 0.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,2,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 1.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_node_add(ref_node,3,&node),"add node");
  ref_node_xyz(ref_node,0,node) = 0.0;
  ref_node_xyz(ref_node,1,node) = 1.0;
  ref_node_xyz(ref_node,2,node) = 1.0;

  RSS(ref_cell_add(ref_grid_pri(ref_grid),nodes,&cell),"add prism");

  nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 4; nodes[1] = 7; nodes[2] = 6; nodes[3] = 5; nodes[4] = 10;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 1; nodes[1] = 5; nodes[2] = 6; nodes[3] = 2; nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 0; nodes[1] = 3; nodes[2] = 7; nodes[3] = 4; nodes[4] = 20;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 0; nodes[1] = 4; nodes[2] = 5; nodes[3] = 1; nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  nodes[0] = 2; nodes[1] = 6; nodes[2] = 7; nodes[3] = 3; nodes[4] = 30;
  RSS(ref_cell_add(ref_grid_qua(ref_grid),nodes,&cell),"add quad");

  return REF_SUCCESS;
}
