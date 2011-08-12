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

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TFS(ref_cell_free(NULL),"dont free NULL");

  TSS(ref_cell_create(4,&ref_cell),"create");
  TES(0,ref_cell_n(ref_cell),"init zero cells");

  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TES(0,cell,"first cell is zero");
  TSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");
  TES(1,cell,"second cell is one");

  max = ref_cell_max(ref_cell);
  for (i = 0; i < max; i++ ) RSS(ref_cell_add(ref_cell,nodes,&cell),"add cell");

  TAS(ref_cell_max(ref_cell)>max,"realloc max");

  TSS(ref_cell_free(ref_cell),"free");
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

  return 0;
}
