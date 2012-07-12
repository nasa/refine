#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_cell.h"
#include "ref_sort.h"
#include "ref_matrix.h"

#include "ref_node.h"
#include "ref_list.h"

#include "ref_mpi.h"

static REF_STATUS ref_tri(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,3,REF_TRUE);
}
static REF_STATUS ref_qua(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,4,REF_TRUE);
}

static REF_STATUS ref_tet(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,4,REF_FALSE);
}
static REF_STATUS ref_pyr(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,5,REF_FALSE);
}
static REF_STATUS ref_pri(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,6,REF_FALSE);
}
static REF_STATUS ref_hex(REF_CELL *ref_cell_ptr)
{
  return ref_cell_create(ref_cell_ptr,8,REF_FALSE);
}

int main( void )
{

  REIS(REF_NULL,ref_cell_free(NULL),"dont free NULL");

  { /* add */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tet(&ref_cell),"create");
    RES(0,ref_cell_n(ref_cell),"init zero cells");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RES(0,cell,"first cell is zero");
    RES(1,ref_cell_n(ref_cell),"first cell incements n");
    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RES(1,cell,"second cell is one");
    RES(2,ref_cell_n(ref_cell),"second cell incements n");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* add many global */
    REF_CELL ref_cell;
    REF_NODE ref_node;
    REF_INT nodes[4];
    REF_INT parts[4];
    REF_INT retrieved[4];

    RSS(ref_node_create(&ref_node),"create node");

    RSS(ref_tet(&ref_cell),"create");
    RES(0,ref_cell_n(ref_cell),"init zero cells");

    nodes[0]= 0; nodes[1]= 10; nodes[2]= 20; nodes[3]= 30;
    parts[0]= 0; parts[1]= 1; parts[2]= 2; parts[3]= 3;
    RSS(ref_cell_add_many_global(ref_cell,ref_node,
				 1,nodes,parts,REF_EMPTY),"add many");

    RSS(ref_cell_nodes(ref_cell,0,retrieved),"cell should exist");
    RES(0,retrieved[0],"node 0");
    RES(1,retrieved[1],"node 1");
    RES(2,retrieved[2],"node 2");
    RES(3,retrieved[3],"node 3");

    RES(20,ref_node_global(ref_node,2),"global mapped");
    RES(2,ref_node_part(ref_node,2),"part rec");

    RSS(ref_cell_free(ref_cell),"cleanup");
    RSS(ref_node_free(ref_node),"cleanup");
  }

  { /* remove */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT item;
    REF_INT cell;

    RSS(ref_tet(&ref_cell),"create");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    nodes[0]= 0; nodes[1]= 4; nodes[2]= 5; nodes[3]= 6; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    REIS(REF_INVALID,ref_cell_remove(ref_cell,3),"remove cell missing cell");
    RES(2,ref_cell_n(ref_cell),"still there");

    RSS(ref_cell_remove(ref_cell,0),"remove cell");
    RES(1,ref_cell_n(ref_cell),"reduced count")
    RAS(!ref_cell_valid(ref_cell,0),"cell is still here");

    each_ref_adj_node_item_with_ref( (ref_cell)->ref_adj, 0, item, cell)
      RAS( !(cell == 0), "removed cell still in adj");

    each_ref_adj_node_item_with_ref( (ref_cell)->ref_adj, 2, item, cell)
      RAS( !(cell == 0), "removed cell still in adj");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* remove tri*/
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell),"create");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 1;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_remove(ref_cell,cell),"remove cell");

    RAS( ref_adj_empty( ref_cell_adj(ref_cell), 1 ), "id");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* replace whole tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell),"create");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    nodes[0]= 2; nodes[1]= 1; nodes[2]= 3; nodes[3]= 20;
    RSS(ref_cell_replace_whole(ref_cell,cell,nodes),"replace cell");

    RAS(  ref_adj_empty( ref_cell_adj(ref_cell), 0 ), "old node");
    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 3 ), "new node");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* replace node of tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT retrieved[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell),"create");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_replace_node(ref_cell,2,3),"replace node");

    RSS(ref_cell_nodes(ref_cell,cell,retrieved),"cell should exist");
    REIS(0,retrieved[0],"node 0");
    REIS(1,retrieved[1],"node 1");
    REIS(3,retrieved[2],"node 2");
    REIS(10,retrieved[3],"id");

    RAS(  ref_adj_empty( ref_cell_adj(ref_cell), 2 ), "old node");
    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 3 ), "new node");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* replace node of 2 tri */
    REF_CELL ref_cell;
    REF_INT nodes[4];
    REF_INT cell;

    RSS(ref_tri(&ref_cell),"create");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    nodes[0]= 2; nodes[1]= 1; nodes[2]= 3; nodes[3]= 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_replace_node(ref_cell,2,4),"replace node");

    RAS(  ref_adj_empty( ref_cell_adj(ref_cell), 2 ), "old node");
    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 4 ), "new node");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* force realloc twice */
    REF_CELL ref_cell;
    REF_INT max, i;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell),"create");

    max = ref_cell_max(ref_cell);
    for (i = 0; i < max+1; i++ ){
      nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
      RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
      RES(i,cell,"expected ordered new cells first block");
    }
    RAS(ref_cell_max(ref_cell)>max,"realloc max");

    max = ref_cell_max(ref_cell);
    for (i = ref_cell_n(ref_cell); i < max+1; i++ ){
      nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
      RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
      RES(i,cell,"expected ordered new cells second block");
    }
    RAS(ref_cell_max(ref_cell)>max,"realloc max");

    RSS(ref_cell_free(ref_cell),"free");
  }

  { /* get nodes */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT retrieved[4];

    RSS(ref_tet(&ref_cell),"create new");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3;
    REIS(REF_INVALID,ref_cell_nodes(ref_cell,0,nodes),"missing cell");
    REIS(REF_INVALID,ref_cell_nodes(ref_cell,-1,nodes),"-1 cell should fail");
    REIS(REF_INVALID,ref_cell_nodes(ref_cell,1000000000,nodes),"large cell");

    nodes[0]= 10; nodes[1]= 20; nodes[2]= 30; nodes[3]= 40; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_nodes(ref_cell,0,retrieved),"cell should exist");
    RES(10,retrieved[0],"node 0");
    RES(20,retrieved[1],"node 1");
    RES(30,retrieved[2],"node 2");
    RES(40,retrieved[3],"node 3");

    RSS(ref_cell_free(ref_cell),"free");
  }

  { /* valid? */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    RSS(ref_tet(&ref_cell),"create new");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
    RAS(!ref_cell_valid(ref_cell,-1),"invlid -1");
    RAS(!ref_cell_valid(ref_cell,-1),"invlid 0");
    RAS(!ref_cell_valid(ref_cell,-1),"invlid 1");
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RAS(ref_cell_valid(ref_cell,0),"valid 0");
  
    RSS(ref_cell_free(ref_cell),"free");
  }

  { /* adj */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell),"create new");

    RAS( ref_adj_empty( ref_cell_adj(ref_cell), 0 ), "first node");
    RAS( ref_adj_empty( ref_cell_adj(ref_cell), 3 ), "last node");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 0 ), "first node");
    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 3 ), "last node");

    RSS(ref_cell_free(ref_cell),"free");
  }

  {/* adj with id */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];

    RSS(ref_tri(&ref_cell),"create new");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 0 ), "first node");
    RAS( !ref_adj_empty( ref_cell_adj(ref_cell), 2 ), "last node");
    RAS(  ref_adj_empty( ref_cell_adj(ref_cell), 3 ), "id");

    RSS(ref_cell_free(ref_cell),"free");
  }

  { /* loop cells */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell;

    RSS(ref_tet(&ref_cell),"create new");

    ncell = 0;
    each_ref_cell_valid_cell( ref_cell, cell )
      ncell++;

    RES(0,ncell,"start zero cells");

    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    nodes[0]= -1; nodes[1]= -1; nodes[2]= -1; nodes[3]= -1; 

    ncell = 0;
    each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
      ncell++;

    RES(1,ncell,"found cell");
    RES(0,nodes[0],"got node 0");
    RES(1,nodes[1],"got node 1");
    RES(2,nodes[2],"got node 2");
    RES(3,nodes[3],"got node 3");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* edge_per */
    REF_CELL ref_cell;

    RSS(ref_tet(&ref_cell),"create");
    REIS(6,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_pyr(&ref_cell),"create");
    REIS(8,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_pri(&ref_cell),"create");
    REIS(9,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_hex(&ref_cell),"create");
    REIS(12,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_tri(&ref_cell),"create");
    REIS(3,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_qua(&ref_cell),"create");
    REIS(4,ref_cell_edge_per(ref_cell),"edge_per");
    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* empty edges */
    REF_CELL ref_cell;

    RSS(ref_tet(&ref_cell),"create");
    RSS(ref_cell_empty_edges(ref_cell),"empty edges");
    RES(REF_EMPTY,ref_cell_c2e(ref_cell,0,0),"edge");
    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* set tet edges */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];

    RSS(ref_tet(&ref_cell),"create");
    nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RSS(ref_cell_empty_edges(ref_cell),"empty edges");

    RSS(ref_cell_set_edge(ref_cell,0,1,0),"set edge");
    RSS(ref_cell_set_edge(ref_cell,0,2,1),"set edge");
    RSS(ref_cell_set_edge(ref_cell,0,3,2),"set edge");
    RSS(ref_cell_set_edge(ref_cell,1,2,3),"set edge");
    RSS(ref_cell_set_edge(ref_cell,1,3,4),"set edge");
    RSS(ref_cell_set_edge(ref_cell,2,3,5),"set edge");

    RES(0,ref_cell_c2e(ref_cell,0,0),"edge");
    RES(1,ref_cell_c2e(ref_cell,1,0),"edge");
    RES(2,ref_cell_c2e(ref_cell,2,0),"edge");
    RES(3,ref_cell_c2e(ref_cell,3,0),"edge");
    RES(4,ref_cell_c2e(ref_cell,4,0),"edge");
    RES(5,ref_cell_c2e(ref_cell,5,0),"edge");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* face_per */
    REF_CELL ref_cell;

    RSS(ref_tet(&ref_cell),"create");
    RES(4,ref_cell_face_per(ref_cell),"face_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_pyr(&ref_cell),"create");
    RES(5,ref_cell_face_per(ref_cell),"face_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_pri(&ref_cell),"create");
    RES(5,ref_cell_face_per(ref_cell),"face_per");
    RSS(ref_cell_free(ref_cell),"cleanup");

    RSS(ref_hex(&ref_cell),"create");
    RES(6,ref_cell_face_per(ref_cell),"face_per");
    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* tet face */
    REF_CELL ref_cell;
    RSS(ref_tet(&ref_cell),"create");
    REIS(1,ref_cell_f2n_gen(ref_cell,0,0),"tri face nodes");
    REIS(3,ref_cell_f2n_gen(ref_cell,1,0),"tri face nodes");
    REIS(2,ref_cell_f2n_gen(ref_cell,2,0),"tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell,0,0),ref_cell_f2n_gen(ref_cell,3,0),
	 "quad not tri");
    REIS(0,ref_cell_f2n_gen(ref_cell,0,1),"tri face nodes");
    REIS(2,ref_cell_f2n_gen(ref_cell,1,1),"tri face nodes");
    REIS(3,ref_cell_f2n_gen(ref_cell,2,1),"tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell,0,1),ref_cell_f2n_gen(ref_cell,3,1),
	 "quad not tri");
    REIS(0,ref_cell_f2n_gen(ref_cell,0,2),"tri face nodes");
    REIS(3,ref_cell_f2n_gen(ref_cell,1,2),"tri face nodes");
    REIS(1,ref_cell_f2n_gen(ref_cell,2,2),"tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell,0,2),ref_cell_f2n_gen(ref_cell,3,2),
	 "quad not tri");
    REIS(0,ref_cell_f2n_gen(ref_cell,0,3),"tri face nodes");
    REIS(1,ref_cell_f2n_gen(ref_cell,1,3),"tri face nodes");
    REIS(2,ref_cell_f2n_gen(ref_cell,2,3),"tri face nodes");
    REIS(ref_cell_f2n_gen(ref_cell,0,3),ref_cell_f2n_gen(ref_cell,3,3),
	 "quad not tri");
    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* tet face nodes*/
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT cell_face;
    RSS(ref_tet(&ref_cell),"create");

    nodes[0]= 10; nodes[1]= 20; nodes[2]= 30; nodes[3]= 40; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    cell = 0;
    cell_face = 3;
    REIS(10,ref_cell_f2n(ref_cell,0,cell_face,cell),"tri face nodes");
    REIS(20,ref_cell_f2n(ref_cell,1,cell_face,cell),"tri face nodes");
    REIS(30,ref_cell_f2n(ref_cell,2,cell_face,cell),"tri face nodes");
    REIS(ref_cell_f2n(ref_cell,0,cell_face,cell),
	 ref_cell_f2n(ref_cell,3,cell_face,cell),
	 "tri face nodes");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* hex edge faces */
    REF_CELL ref_cell;
    REF_INT edge, face0, face1;
    RSS(ref_hex(&ref_cell),"create");

    edge = 0;
    RSS(ref_cell_gen_edge_face(ref_cell,edge,&face0,&face1),"edge face");
    REIS(0,face0,"face0");
    REIS(4,face1,"face1");

    edge = 1;
    RSS(ref_cell_gen_edge_face(ref_cell,edge,&face0,&face1),"edge face");
    REIS(3,face0,"face0");
    REIS(4,face1,"face1");

    edge = 2;
    RSS(ref_cell_gen_edge_face(ref_cell,edge,&face0,&face1),"edge face");
    REIS(0,face0,"face0");
    REIS(3,face1,"face1");

    edge = 6;
    RSS(ref_cell_gen_edge_face(ref_cell,edge,&face0,&face1),"edge face");
    REIS(1,face0,"face0");
    REIS(2,face1,"face1");

    edge = 11;
    RSS(ref_cell_gen_edge_face(ref_cell,edge,&face0,&face1),"edge face");
    REIS(2,face0,"face0");
    REIS(5,face1,"face1");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* tet with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT found;

    RSS(ref_tet(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_with(ref_cell,nodes,&found),"not found");
    REIS(cell,found, "not same");

    nodes[0]=5;
    REIS(REF_NOT_FOUND,ref_cell_with(ref_cell,nodes,&found),"found");
    REIS(REF_EMPTY,found, "not same")

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* list of one tet with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell,list[5];

    RSS(ref_tet(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_list_with(ref_cell,0,1,5,&ncell,list),"no list");
    REIS(1,ncell, "mis count");
    REIS(0,list[0], "not in list");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* list of two tet with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT ncell,list[5];

    RSS(ref_tet(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    nodes[0] = 0; nodes[1] = 2; nodes[2] = 1; nodes[3] = 4;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_list_with(ref_cell,0,1,5,&ncell,list),"no list");
    REIS(2,ncell, "mis count");
    REIS(1,list[0], "not in list");
    REIS(0,list[1], "not in list");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* list of node */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT nnode,list[6];

    RSS(ref_tet(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 3;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_node_list_around(ref_cell,0,6,&nnode,list),"no list");
    REIS(3,nnode, "mis count");
    REIS(1,list[0], "not in list");
    REIS(2,list[1], "not in list");
    REIS(3,list[2], "not in list");

    nodes[0] = 0; nodes[1] = 2; nodes[2] = 1; nodes[3] = 4;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_node_list_around(ref_cell,0,6,&nnode,list),"no list");
    REIS(4,nnode, "mis count");
    REIS(2,list[0], "not in list");
    REIS(1,list[1], "not in list");
    REIS(4,list[2], "not in list");
    REIS(3,list[3], "not in list");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* tri with */
    REF_CELL ref_cell;
    REF_INT cell, nodes[4];
    REF_INT found;

    RSS(ref_tri(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 10;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    RSS(ref_cell_with(ref_cell,nodes,&found),"not found");
    REIS(cell,found, "not same");

    nodes[3]=5;
    RSS(ref_cell_with(ref_cell,nodes,&found),"not found");
    REIS(cell,found, "not same");

    nodes[0]=5;
    REIS(REF_NOT_FOUND,ref_cell_with(ref_cell,nodes,&found),"found");
    REIS(REF_EMPTY,found, "not same")

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* cell has face */
    REF_CELL ref_cell;
    REF_INT cell;
    REF_INT nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT face_nodes[REF_CELL_MAX_SIZE_PER];
    REF_INT found;

    RSS(ref_pri(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; 
    nodes[3] = 3; nodes[4] = 4; nodes[5] = 5; 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    face_nodes[0]=0; face_nodes[1]=1; face_nodes[2]=5; face_nodes[3]=0; 
    REIS(REF_NOT_FOUND,ref_cell_with_face(ref_cell,face_nodes,&found),"with");
    REIS(REF_EMPTY,found,"false positive");

    face_nodes[0]=1; face_nodes[1]=0; face_nodes[2]=3; face_nodes[3]=4; 
    RSS(ref_cell_with_face(ref_cell,face_nodes,&found),"with");
    REIS(0,found,"false negative");

    face_nodes[0]=0; face_nodes[1]=1; face_nodes[2]=2; face_nodes[3]=0; 
    RSS(ref_cell_with_face(ref_cell,face_nodes,&found),"with");
    REIS(0,found,"false negative");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }

  { /* has side */
    REF_CELL ref_cell;
    REF_INT cell, nodes[8];
    REF_INT node0, node1;
    REF_BOOL has_side;

    RSS(ref_hex(&ref_cell),"create");

    nodes[0] = 0; nodes[1] = 1; nodes[2] = 2; nodes[3] = 4;
    nodes[4] = 4; nodes[5] = 5; nodes[6] = 6; nodes[7] = 7;
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    node0 = 0; node1 = 1;
    has_side = REF_FALSE;
    RSS(ref_cell_has_side(ref_cell,node0,node1,&has_side),"side");
    REIS(REF_TRUE,has_side, "side expected");

    node0 = 0; node1 = 2;
    RSS(ref_cell_has_side(ref_cell,node0,node1,&has_side),"side");
    REIS(REF_FALSE,has_side, "diagonal, not side");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }
  
  { /* triangle has face */
    REF_CELL ref_cell;
    REF_INT cell, nodes[8];

    RSS(ref_tri(&ref_cell),"create");
    nodes[0] = 0;
    nodes[1] = 1;
    nodes[2] = 2;
    nodes[3] = 10;
 
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

    nodes[3] = 0;
    cell = REF_EMPTY;
    RSS( ref_cell_with_face(ref_cell, nodes, &cell),"has");
    REIS( 0, cell, "wrong cell");

    nodes[0] = 1;
    nodes[1] = 2;
    nodes[2] = 0;
    nodes[3] = 1;
    cell = REF_EMPTY;
    RSS( ref_cell_with_face(ref_cell, nodes, &cell),"has");
    REIS( 0, cell, "wrong cell");

    RSS(ref_cell_free(ref_cell),"cleanup");
  }


  return 0;
}
