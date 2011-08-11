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

  TSS(ref_node_create(&ref_node),"create");
  TES(0,ref_node_n(ref_node),"init zero nodes");
  TSS(ref_node_free(ref_node),"free");

  global = 10;
  TSS(ref_node_add(ref_node,global,&node),"add");
  TES(0,node,"first node is zero");


  return 0;
}
