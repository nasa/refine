#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_node.h"
#include "ref_sort.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  REIS(REF_NULL,ref_node_free(NULL),"dont free NULL");

  { /* init */
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");
    REIS(0,ref_node_n(ref_node),"init zero nodes");

    RSS(ref_node_free(ref_node),"free");
  }

  {
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    RES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for missing node");

    /* first add in order */
    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"first add");
    RES(0,node,"first node is zero");
    RES(1,ref_node_n(ref_node),"count incremented");
    RES(global,ref_node_global(ref_node,0),"global match for first node");

    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"second add");
    RES(1,node,"second node is one");
    RES(2,ref_node_n(ref_node),"count incremented");
    RES(global,ref_node_global(ref_node,1),"global match for second node");

    /* removed node invalid */
    REIS(REF_INVALID,ref_node_remove(ref_node,-1),"remove invalid node");
    REIS(REF_INVALID,ref_node_remove(ref_node,2),"remove invalid node");

    RSS(ref_node_remove(ref_node,0),"remove first node");
    RES(REF_EMPTY,ref_node_global(ref_node,0),"global empty for removed node");
    RES(1,ref_node_n(ref_node),"count decremented");

    global = 30;
    RSS(ref_node_add(ref_node,global,&node),"replace");
    RES(0,node,"reuse removed node");
    RES(global,ref_node_global(ref_node,node),"global match for replaced node");
    RES(2,ref_node_n(ref_node),"count incremented");

    RES(20,ref_node_global(ref_node,1),"global match for second node");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* add bunch testing realloc */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    max = ref_node_max(ref_node);
    for ( global = 10; global < 10*(max+2) ; global += 10 )
      RSS(ref_node_add(ref_node,global,&node),"realloc");

    RAS(max < ref_node_max(ref_node),"grow max");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* lookup local from global */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    node = 0;
    REIS(REF_NOT_FOUND,ref_node_local(ref_node,-1,&node),"invalid global");
    RES(REF_EMPTY,node,"expect node empty for invalid global");
    REIS(REF_NOT_FOUND,ref_node_local(ref_node,5,&node),"invalid global");
    REIS(REF_NOT_FOUND,ref_node_local(ref_node,200,&node),"invalid global");

    RSS(ref_node_local(ref_node,10,&node),"return global");
    REIS(0,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* lookup local from global after remove */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"add");

    RSS(ref_node_remove(ref_node,0),"remove");

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(1,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  {  /* compact nodes */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT *o2n;
    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,1,&node),"add");
    RSS(ref_node_add(ref_node,3,&node),"add");
    RSS(ref_node_add(ref_node,2,&node),"add");
    RSS(ref_node_remove(ref_node,1),"remove");

    RSS(ref_node_compact(ref_node,&o2n),"compact");
 
    RES(0,o2n[0],"o2n");
    RES(REF_EMPTY,o2n[1],"o2n");
    RES(1,o2n[2],"o2n");
    free(o2n);
    
    RSS(ref_node_free(ref_node),"free");
  }

  {  /* valid */
    REF_INT node;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    RAS(!ref_node_valid(ref_node,0),"empty invalid");
    RSS(ref_node_add(ref_node,0,&node),"add 0 global");
    RAS(ref_node_valid(ref_node,0),"zero is valid global");
    RES(0,ref_node_global(ref_node,0),"zero global");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* unique */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"first");
    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"second");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"first");
    REIS(0,node,"return first");
    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"second");
    REIS(1,node,"return second");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* sorted_global rebuild */
    REF_INT global, node;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    global = 20;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    global = 10;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    global = 30;
    RSS(ref_node_add(ref_node,global,&node),"realloc");

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(0,node,"wrong local");
    RSS(ref_node_local(ref_node,10,&node),"return global");
    REIS(1,node,"wrong local");
    RSS(ref_node_local(ref_node,30,&node),"return global");
    REIS(2,node,"wrong local");

    RSS( ref_node_rebuild_sorted_global( ref_node ), "rebuild" );

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(0,node,"wrong local");
    RSS(ref_node_local(ref_node,10,&node),"return global");
    REIS(1,node,"wrong local");
    RSS(ref_node_local(ref_node,30,&node),"return global");
    REIS(2,node,"wrong local");


    RSS(ref_node_free(ref_node),"free");
  }

  { /* add many to empty */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    REIS(2,ref_node_n(ref_node),"init zero nodes");

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(0,node,"wrong local");
    RSS(ref_node_local(ref_node,10,&node),"return global");
    REIS(1,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* add many to existing */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,10,&node),"many");

    global[0] = 20;
    global[1] = 10;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    REIS(2,ref_node_n(ref_node),"init zero nodes");

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(1,node,"wrong local");
    RSS(ref_node_local(ref_node,10,&node),"return global");
    REIS(0,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* add many duplicates */
    REF_INT n = 2, node;
    REF_INT global[2];
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    global[0] = 20;
    global[1] = 20;

    RSS(ref_node_add_many(ref_node,n,global),"many");

    REIS(1,ref_node_n(ref_node),"init zero nodes");

    RSS(ref_node_local(ref_node,20,&node),"return global");
    REIS(0,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* reuse removed global */
    REF_NODE ref_node;
    REF_INT node, global;

    RSS(ref_node_create(&ref_node),"create");

    global = 3542;
    RSS(ref_node_add(ref_node,global,&node),"add orig");

    RSS(ref_node_remove(ref_node,node),"remove node");

    RSS( ref_node_next_global(ref_node,&node), "next gloabal");
    REIS( global, node, "not reused");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* shift globals */
    REF_INT local, global, node;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS( ref_node_initialize_n_global( ref_node, 30 ), "init n glob" );

    RSS( ref_node_next_global( ref_node, &global ), "next");
    REIS( 30, global, "expected n global");
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS(ref_node_shift_new_globals(ref_node),"shift");

    RSS(ref_node_local(ref_node,30,&node),"return global");
    REIS(2,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* eliminate unused globals */
    REF_INT local, global, node;
    REF_NODE ref_node;

    RSS(ref_node_create(&ref_node),"create");

    global = 10;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 20;
    RSS(ref_node_add(ref_node,global,&local),"add");
    global = 30;
    RSS(ref_node_add(ref_node,global,&local),"add");

    RSS(ref_node_remove(ref_node,1),"rm");

    RSS( ref_node_initialize_n_global( ref_node, 30 ), "init n glob" );

    RSS(ref_node_eliminate_unused_globals(ref_node),"unused");

    RSS(ref_node_local(ref_node,29,&node),"return global");
    REIS(2,node,"wrong local");

    RSS(ref_node_free(ref_node),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
