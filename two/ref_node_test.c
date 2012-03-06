#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_node.h"
#include "ref_sort.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{

  TSS( ref_mpi_start( argc, argv ), "start" );

  TFS(ref_node_free(NULL),"dont free NULL");

  { /* init */
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");
    TEIS(0,ref_node_n(ref_node),"init zero nodes");
    TEIS(10,ref_node_max(ref_node),"init 10 max");

    TSS(ref_node_free(ref_node),"free");
  }

  {
    REF_INT global, node;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    TES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for missing node");

    /* dont allow neg globals */
    /* skipped, it is an assert 
    TFS(ref_node_add(ref_node,-3,&node),"negative global no allowed");
    TES(0,ref_node_n(ref_node),"count not changed");
    */

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

    TSS(ref_node_free(ref_node),"free");
  }

  { /* add bunch testing realloc */
    REF_INT global, node, max;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    max = ref_node_max(ref_node);
    for ( global = 10; global < 10*(max+2) ; global += 10 )
      RSS(ref_node_add(ref_node,global,&node),"realloc");

    TAS(max < ref_node_max(ref_node),"grow max");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* lookup local from global */
    REF_INT global, node;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    node = 0;
    TFS(ref_node_local(ref_node,-1,&node),"returned invalid global");
    TES(REF_EMPTY,node,"expect node empty for invalid global");
    TFS(ref_node_local(ref_node,5,&node),"returned invalid global");
    TFS(ref_node_local(ref_node,200,&node),"returned invalid global");

    TSS(ref_node_local(ref_node,10,&node),"return global");
    TEIS(0,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* lookup local from global after remove */
    REF_INT global, node;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"add");

    TSS(ref_node_remove(ref_node,0),"remove");

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(1,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  {  /* compact nodes */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT *o2n;
    TSS(ref_node_create(&ref_node),"create");

    TSS(ref_node_add(ref_node,1,&node),"add");
    TSS(ref_node_add(ref_node,3,&node),"add");
    TSS(ref_node_add(ref_node,2,&node),"add");
    TSS(ref_node_remove(ref_node,1),"remove");

    TSS(ref_node_compact(ref_node,&o2n),"compact");
 
    TES(0,o2n[0],"o2n");
    TES(REF_EMPTY,o2n[1],"o2n");
    TES(1,o2n[2],"o2n");
    free(o2n);
    
    TSS(ref_node_free(ref_node),"free");
  }

  {  /* valid */
    REF_INT node;
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    TAS(!ref_node_valid(ref_node,0),"empty invalid");
    TSS(ref_node_add(ref_node,0,&node),"add 0 global");
    TAS(ref_node_valid(ref_node,0),"zero is valid global");
    TES(0,ref_node_global(ref_node,0),"zero global");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* unique */
    REF_INT global, node;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    global = 10;
    TSS(ref_node_add(ref_node,global,&node),"first");
    global = 20;
    TSS(ref_node_add(ref_node,global,&node),"second");

    global = 10;
    TSS(ref_node_add(ref_node,global,&node),"first");
    TEIS(0,node,"return first");
    global = 20;
    TSS(ref_node_add(ref_node,global,&node),"second");
    TEIS(1,node,"return second");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* sorted_global rebuild */
    REF_INT global, node;
    REF_NODE ref_node;
    TSS(ref_node_create(&ref_node),"create");

    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    global = 30;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(0,node,"wrong local");
    TSS(ref_node_local(ref_node,10,&node),"return global");
    TEIS(1,node,"wrong local");
    TSS(ref_node_local(ref_node,30,&node),"return global");
    TEIS(2,node,"wrong local");

    RSS( ref_node_rebuild_sorted_global( ref_node ), "rebuild" );

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(0,node,"wrong local");
    TSS(ref_node_local(ref_node,10,&node),"return global");
    TEIS(1,node,"wrong local");
    TSS(ref_node_local(ref_node,30,&node),"return global");
    TEIS(2,node,"wrong local");


    TSS(ref_node_free(ref_node),"free");
  }

  { /* add many to empty */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    TEIS(2,ref_node_n(ref_node),"init zero nodes");

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(0,node,"wrong local");
    TSS(ref_node_local(ref_node,10,&node),"return global");
    TEIS(1,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* add many to existing */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,10,&node),"many");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    TEIS(2,ref_node_n(ref_node),"init zero nodes");

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(1,node,"wrong local");
    TSS(ref_node_local(ref_node,10,&node),"return global");
    TEIS(0,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* add many duplicates */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    global[0] = 20;
    global[1] = 20;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    TEIS(1,ref_node_n(ref_node),"init zero nodes");

    TSS(ref_node_local(ref_node,20,&node),"return global");
    TEIS(0,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* reuse removed global */
    REF_NODE ref_node;
    REF_INT node, global;

    TSS(ref_node_create(&ref_node),"create");

    global = 3542;
    RSS(ref_node_add(ref_node,global,&node),"add orig");

    TSS(ref_node_remove(ref_node,node),"remove node");

    TSS( ref_node_next_global(ref_node,&node), "next gloabal");
    TEIS( global, node, "not reused");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* shift globals */
    REF_INT local, global, node;
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS( ref_node_initialize_n_global( ref_node, 30 ), "init n glob" );

    RSS( ref_node_next_global( ref_node, &global ), "next");
    REIS( 30, global, "expected n global");
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS(ref_node_sync_new_globals(ref_node),"sync");

    TSS(ref_node_local(ref_node,30,&node),"return global");
    TEIS(2,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  { /* eliminate unused globals */
    REF_INT local, global, node;
    REF_NODE ref_node;

    TSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 30;
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS(ref_node_remove(ref_node,1),"rm");

    RSS( ref_node_initialize_n_global( ref_node, 30 ), "init n glob" );

    RSS(ref_node_eliminate_unused_globals(ref_node),"sync");

    TSS(ref_node_local(ref_node,29,&node),"return global");
    TEIS(2,node,"wrong local");

    TSS(ref_node_free(ref_node),"free");
  }

  TSS( ref_mpi_stop( ), "stop" );

  return 0;
}
