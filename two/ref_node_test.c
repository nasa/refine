#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_node.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_NODE ref_node;
  REF_INT global, node;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  TFS(ref_node_free(NULL),"dont free NULL");

  TSS(ref_node_create(&ref_node),"create");
  TES(0,ref_node_n(ref_node),"init zero nodes");
  TES(10,ref_node_max(ref_node),"init 10 max");

  TSS(ref_node_free(ref_node),"free");
  TSS(ref_node_create(&ref_node),"create");

  TES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for missing node");

  global = 10;
  TSS(ref_node_add(ref_node,global,&node),"firat add");
  TES(0,node,"first node is zero");
  TES(global,ref_node_global(ref_node,0),"global match for first node");

  global = 20;
  TSS(ref_node_add(ref_node,global,&node),"second add");
  TES(1,node,"second node is one");
  TES(global,ref_node_global(ref_node,1),"global match for second node");

  TSS(ref_node_free(ref_node),"free");

  return 0;
}
