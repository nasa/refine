#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_cell.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_CELL ref_cell;
  REF_INT nodes[4] = {0,1,2,3};
  REF_INT retrieved[4];
  REF_INT cell;
  REF_INT max, i;
  REF_INT ncell;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TFS(ref_cell_free(NULL),"dont free NULL");

  TSS(ref_cell_create(4,&ref_cell),"create");
  TES(0,ref_cell_n(ref_cell),"init zero cells");

  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TES(0,cell,"first cell is zero");
  TES(1,ref_cell_n(ref_cell),"first cell incements n");
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TES(1,cell,"second cell is one");
  TES(2,ref_cell_n(ref_cell),"second cell incements n");

  TSS(ref_cell_free(ref_cell),"cleanup");
  /* force realloc twice */
  TSS(ref_cell_create(4,&ref_cell),"create");

  max = ref_cell_max(ref_cell);
  for (i = 0; i < max+1; i++ ){
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RES(i,cell,"expected ordered new cells first block");
  }
  TAS(ref_cell_max(ref_cell)>max,"realloc max");

  max = ref_cell_max(ref_cell);
  for (i = ref_cell_n(ref_cell); i < max+1; i++ ){
    RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
    RES(i,cell,"expected ordered new cells second block");
  }
  TAS(ref_cell_max(ref_cell)>max,"realloc max");

  TSS(ref_cell_free(ref_cell),"free");
  /* get nodes */
  TSS(ref_cell_create(4,&ref_cell),"create new");

  TFS(ref_cell_nodes(ref_cell,0,nodes),"missing cell should fail");
  TFS(ref_cell_nodes(ref_cell,-1,nodes),"-1 cell should fail");
  TFS(ref_cell_nodes(ref_cell,1000000000,nodes),"1000000000 cell should fail");

  nodes[0]= 10; nodes[1]= 20; nodes[2]= 30; nodes[3]= 40; 
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

  TSS(ref_cell_nodes(ref_cell,0,retrieved),"cell should exist");
  TES(10,retrieved[0],"node 0");
  TES(20,retrieved[1],"node 1");
  TES(30,retrieved[2],"node 2");
  TES(40,retrieved[3],"node 3");

  TSS(ref_cell_free(ref_cell),"free");
  /* get nodes */
  TSS(ref_cell_create(4,&ref_cell),"create new");

  nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
  TAS(!ref_cell_valid(ref_cell,-1),"invlid -1");
  TAS(!ref_cell_valid(ref_cell,-1),"invlid 0");
  TAS(!ref_cell_valid(ref_cell,-1),"invlid 1");
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TAS(ref_cell_valid(ref_cell,0),"valid 0");
  
  TSS(ref_cell_free(ref_cell),"free");
  /* loop cells */
  TSS(ref_cell_create(4,&ref_cell),"create new");

  ncell = 0;
  each_ref_cell_valid_cell( ref_cell, cell )
    ncell++;

  TES(0,ncell,"start zero cells");

  nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

  nodes[0]= -1; nodes[1]= -1; nodes[2]= -1; nodes[3]= -1; 

  ncell = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    ncell++;

  TES(1,ncell,"found cell");
  TES(0,nodes[0],"got node 0");
  TES(1,nodes[1],"got node 1");
  TES(2,nodes[2],"got node 2");
  TES(3,nodes[3],"got node 3");

  TSS(ref_cell_free(ref_cell),"cleanup");

  /* edge_per */
  TSS(ref_cell_create(4,&ref_cell),"create");
  TES(6,ref_cell_edge_per(ref_cell),"edge_per");
  TSS(ref_cell_free(ref_cell),"cleanup");

  TSS(ref_cell_create(5,&ref_cell),"create");
  TES(8,ref_cell_edge_per(ref_cell),"edge_per");
  TSS(ref_cell_free(ref_cell),"cleanup");

  TSS(ref_cell_create(6,&ref_cell),"create");
  TES(9,ref_cell_edge_per(ref_cell),"edge_per");
  TSS(ref_cell_free(ref_cell),"cleanup");

  TSS(ref_cell_create(8,&ref_cell),"create");
  TES(12,ref_cell_edge_per(ref_cell),"edge_per");
  TSS(ref_cell_free(ref_cell),"cleanup");

  /* empty edges */

  TSS(ref_cell_create(4,&ref_cell),"create");
  TSS(ref_cell_empty_edges(ref_cell),"empty edges");
  TES(REF_EMPTY,ref_cell_c2e(ref_cell,0,0),"edge");
  TSS(ref_cell_free(ref_cell),"cleanup");

  /* set tet edges */

  TSS(ref_cell_create(4,&ref_cell),"create");
  nodes[0]= 0; nodes[1]= 1; nodes[2]= 2; nodes[3]= 3; 
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TSS(ref_cell_empty_edges(ref_cell),"empty edges");

  TSS(ref_cell_set_edge(ref_cell,0,1,0),"set edge");
  TSS(ref_cell_set_edge(ref_cell,0,2,1),"set edge");
  TSS(ref_cell_set_edge(ref_cell,0,3,2),"set edge");
  TSS(ref_cell_set_edge(ref_cell,1,2,3),"set edge");
  TSS(ref_cell_set_edge(ref_cell,1,3,4),"set edge");
  TSS(ref_cell_set_edge(ref_cell,2,3,5),"set edge");

  TES(0,ref_cell_c2e(ref_cell,0,0),"edge");
  TES(1,ref_cell_c2e(ref_cell,1,0),"edge");
  TES(2,ref_cell_c2e(ref_cell,2,0),"edge");
  TES(3,ref_cell_c2e(ref_cell,3,0),"edge");
  TES(4,ref_cell_c2e(ref_cell,4,0),"edge");
  TES(5,ref_cell_c2e(ref_cell,5,0),"edge");

  TSS(ref_cell_free(ref_cell),"cleanup");

  return 0;
}
