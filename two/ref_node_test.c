#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_node.h"
#include "ref_test.h"

int main( int argc, char *argv[] )
{
  REF_NODE ref_node;
  REF_INT global, node, max;

  if (argc>1) {printf("%s ignored\n",argv[0]);}

  /* init */

  TFS(ref_node_free(NULL),"dont free NULL");

  TSS(ref_node_create(&ref_node),"create");
  TES(0,ref_node_n(ref_node),"init zero nodes");
  TES(10,ref_node_max(ref_node),"init 10 max");

  TSS(ref_node_free(ref_node),"free");
  TSS(ref_node_create(&ref_node),"create");

  TES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for missing node");

  /* dont allow neg globals */
  TFS(ref_node_add(ref_node,-3,&node),"negative global no allowed");
  TES(0,ref_node_n(ref_node),"count not changed");

  /* first add in order */

  global = 10;
  TSS(ref_node_add(ref_node,global,&node),"first add");
  TES(0,node,"first node is zero");
  TES(1,ref_node_n(ref_node),"count incremented");
  TES(global,ref_node_global(ref_node,0),"global match for first node");

  global = 20;
  TSS(ref_node_add(ref_node,global,&node),"second add");
  TES(1,node,"second node is one");
  TES(2,ref_node_n(ref_node),"count incremented");
  TES(global,ref_node_global(ref_node,1),"global match for second node");

  /* removed node invalid */

  TFS(ref_node_remove(ref_node,-1),"remove invalid node");
  TFS(ref_node_remove(ref_node,2),"remove invalid node");

  TSS(ref_node_remove(ref_node,0),"remove first node");
  TES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for removed node");
  TES(1,ref_node_n(ref_node),"count decremented");

  global = 30;
  TSS(ref_node_add(ref_node,global,&node),"replace");
  TES(0,node,"reuse removed node");
  TES(global,ref_node_global(ref_node,node),"global match for replaced node");
  TES(2,ref_node_n(ref_node),"count incremented");

  TES(20,ref_node_global(ref_node,1),"global match for second node");

  /* add bunch testing realloc */

  TSS(ref_node_free(ref_node),"free");
  TSS(ref_node_create(&ref_node),"create");

  max = ref_node_max(ref_node);
  for ( global = 10; global < 10*(max+2) ; global += 10 )
    RSS(ref_node_add(ref_node,global,&node),"realloc");

  TAS(max < ref_node_max(ref_node),"grow max");

  /* lookup local from global */

  TSS(ref_node_free(ref_node),"free");

  return 0;
}
