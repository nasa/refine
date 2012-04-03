#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_node.h"
#include "ref_sort.h"
#include "ref_mpi.h"
#include "ref_matrix.h"

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

  { /* distance in metric */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio, h;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&node0),"add");
    ref_node_xyz(ref_node,0,node0) = 0.0;
    ref_node_xyz(ref_node,1,node0) = 0.0;
    ref_node_xyz(ref_node,2,node0) = 0.0;
    ref_node_metric(ref_node,0,node0) = 1.0;
    ref_node_metric(ref_node,1,node0) = 0.0;
    ref_node_metric(ref_node,2,node0) = 0.0;
    ref_node_metric(ref_node,3,node0) = 1.0;
    ref_node_metric(ref_node,4,node0) = 0.0;
    ref_node_metric(ref_node,5,node0) = 1.0;

    global = 1;
    RSS(ref_node_add(ref_node,global,&node1),"add");
    ref_node_xyz(ref_node,0,node1) = 0.0;
    ref_node_xyz(ref_node,1,node1) = 0.0;
    ref_node_xyz(ref_node,2,node1) = 0.0;
    ref_node_metric(ref_node,0,node1) = 1.0;
    ref_node_metric(ref_node,1,node1) = 0.0;
    ref_node_metric(ref_node,2,node1) = 0.0;
    ref_node_metric(ref_node,3,node1) = 1.0;
    ref_node_metric(ref_node,4,node1) = 0.0;
    ref_node_metric(ref_node,5,node1) = 1.0;

    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RWDS( 0.0, ratio, -1.0, "ratio expected" );

    ref_node_xyz(ref_node,0,node1) = 1.0;
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RWDS( 1.0, ratio, -1.0, "ratio expected" );

    h = 0.5;
    ref_node_metric(ref_node,0,node0) = 1.0/(h*h);
    ref_node_metric(ref_node,0,node1) = 1.0/(h*h);
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RWDS( 2.0, ratio, -1.0, "ratio expected" );

    h = 0.1;
    ref_node_metric(ref_node,0,node0) = 1.0/(h*h);
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RWDS( 4.970679, ratio, 0.00001, "ratio expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* tet volume */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL vol;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");
    global = 3;
    RSS(ref_node_add(ref_node,global,&(nodes[3])),"add");

    for ( global=0;global<4;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
      }

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,1,nodes[2]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[3]) = 1.0;

    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS( 1.0/6.0, vol, -1.0, "vol expected" );

    /* inverted tet is negative volume */
    ref_node_xyz(ref_node,2,nodes[3]) = -1.0;
    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS( -1.0/6.0, vol, -1.0, "vol expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* right tet quality */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");
    global = 3;
    RSS(ref_node_add(ref_node,global,&(nodes[3])),"add");

    for ( global=0;global<4;global++)
      {
	ref_node_xyz(ref_node,0,global) = 0.0;
	ref_node_xyz(ref_node,1,global) = 0.0;
	ref_node_xyz(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
      }

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,1,nodes[2]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[3]) = 1.0;

    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "qual expected" );

    /* inverted tet is negative volume */
    ref_node_xyz(ref_node,2,nodes[3]) = -1.0;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( -1.0/6.0, qual, -1.0, "qual expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* Regular Tetrahedron vol, quality, ratio */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL qual, vol, ratio;
    REF_DBL a;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");
    global = 3;
    RSS(ref_node_add(ref_node,global,&(nodes[3])),"add");

    for ( global=0;global<4;global++)
      {
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
      }

    a = 1.0;

    ref_node_xyz(ref_node,0,nodes[0]) = 1.0/3.0*sqrt(3.0)*a;
    ref_node_xyz(ref_node,1,nodes[0]) = 0.0;
    ref_node_xyz(ref_node,2,nodes[0]) = 0.0;

    ref_node_xyz(ref_node,0,nodes[1]) = -1.0/6.0*sqrt(3.0)*a;
    ref_node_xyz(ref_node,1,nodes[1]) = 0.5*a;
    ref_node_xyz(ref_node,2,nodes[1]) = 0.0;

    ref_node_xyz(ref_node,0,nodes[2]) = -1.0/6.0*sqrt(3.0)*a;
    ref_node_xyz(ref_node,1,nodes[2]) = -0.5*a;
    ref_node_xyz(ref_node,2,nodes[2]) = 0.0;

    ref_node_xyz(ref_node,0,nodes[3]) = 0.0;
    ref_node_xyz(ref_node,1,nodes[3]) = 0.0;
    ref_node_xyz(ref_node,2,nodes[3]) = 1.0/3.0*sqrt(6.0)*a;

    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 1.0, qual, -1.0, "qual expected" );

    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS( 1.0/12.0*sqrt(2.0), vol, -1.0, "vol expected" );

    RSS( ref_node_ratio(ref_node, nodes[2], nodes[3], &ratio), "ratio" );
    RWDS( 1.0, ratio, -1.0, "ratio expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* interpolate */
    REF_NODE ref_node;
    REF_INT node0, node1, new_node, global;

    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_next_global( ref_node, &global ), "next_global");
    RSS(ref_node_add(ref_node,global,&node0),"add");
    ref_node_xyz(ref_node,0,node0) = 0.0;
    ref_node_xyz(ref_node,1,node0) = 0.0;
    ref_node_xyz(ref_node,2,node0) = 0.0;
    ref_node_metric(ref_node,0,node0) = 1.0;
    ref_node_metric(ref_node,1,node0) = 0.0;
    ref_node_metric(ref_node,2,node0) = 0.0;
    ref_node_metric(ref_node,3,node0) = 1.0;
    ref_node_metric(ref_node,4,node0) = 0.0;
    ref_node_metric(ref_node,5,node0) = 1.0;

    RSS(ref_node_next_global( ref_node, &global ), "next_global");
    RSS(ref_node_add(ref_node,global,&node1),"add");
    ref_node_xyz(ref_node,0,node1) = 1.0;
    ref_node_xyz(ref_node,1,node1) = 0.0;
    ref_node_xyz(ref_node,2,node1) = 0.0;
    ref_node_metric(ref_node,0,node1) = 1.0;
    ref_node_metric(ref_node,1,node1) = 0.0;
    ref_node_metric(ref_node,2,node1) = 0.0;
    ref_node_metric(ref_node,3,node1) = 1.0;
    ref_node_metric(ref_node,4,node1) = 0.0;
    ref_node_metric(ref_node,5,node1) = 1.0;

    RSS(ref_node_next_global( ref_node, &global ), "next_global");
    RSS(ref_node_add(ref_node,global,&new_node),"add");
    RSS(ref_node_interpolate_edge(ref_node, node0, node1, new_node),"interp");

    RWDS( 0.5, ref_node_xyz(ref_node,0,new_node), -1.0, "x");
    RWDS( 0.0, ref_node_xyz(ref_node,1,new_node), -1.0, "y");
    RWDS( 0.0, ref_node_xyz(ref_node,2,new_node), -1.0, "z");

    RWDS( 1.0, ref_node_metric(ref_node,0,new_node), -1.0, "m0");
    RWDS( 0.0, ref_node_metric(ref_node,1,new_node), -1.0, "m1");
    RWDS( 0.0, ref_node_metric(ref_node,2,new_node), -1.0, "m2");
    RWDS( 1.0, ref_node_metric(ref_node,3,new_node), -1.0, "m3");
    RWDS( 0.0, ref_node_metric(ref_node,4,new_node), -1.0, "m4");
    RWDS( 1.0, ref_node_metric(ref_node,5,new_node), -1.0, "m5");

    ref_node_metric(ref_node,0,node1) = 1.0/(0.1*0.1);

    RSS(ref_node_interpolate_edge(ref_node, node0, node1, new_node),"interp");
    RWDS( 1.0/0.1, ref_node_metric(ref_node,0,new_node), -1.0, "m0");
    RWDS( 0.0, ref_node_metric(ref_node,1,new_node), -1.0, "m1");
    RWDS( 0.0, ref_node_metric(ref_node,2,new_node), -1.0, "m2");
    RWDS( 1.0, ref_node_metric(ref_node,3,new_node), -1.0, "m3");
    RWDS( 0.0, ref_node_metric(ref_node,4,new_node), -1.0, "m4");
    RWDS( 1.0, ref_node_metric(ref_node,5,new_node), -1.0, "m5");

    RSS(ref_node_free(ref_node),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
