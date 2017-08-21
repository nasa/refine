
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_node.h"
#include  "ref_mpi.h"
#include  "ref_matrix.h"
#include  "ref_sort.h"
#include  "ref_list.h"

#include "ref_malloc.h"

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

  { /* deep copy empty */
    REF_NODE original, copy;
    RSS(ref_node_create(&original),"create");
    RSS(ref_node_deep_copy(&copy,original),"deep copy");

    RSS(ref_node_free(original),"free");
    RSS(ref_node_free(copy),"free");
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

  { /* remove max node */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    max = ref_node_max(ref_node);
    for ( global = 0; global < max ; global += 1 )
      RSS(ref_node_add(ref_node,global,&node),"realloc");

    RSS(ref_node_remove(ref_node,max-1),"remove last node");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* remove max node without global */
    REF_INT global, node, max;
    REF_NODE ref_node;
    RSS(ref_node_create(&ref_node),"create");

    max = ref_node_max(ref_node);
    for ( global = 0; global < max ; global += 1 )
      RSS(ref_node_add(ref_node,global,&node),"realloc");

    RSS(ref_node_remove_without_global(ref_node,max-1),"remove last node");

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
    REF_INT *o2n, *n2o;
    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,1,&node),"add");
    RSS(ref_node_add(ref_node,3,&node),"add");
    RSS(ref_node_add(ref_node,2,&node),"add");
    RSS(ref_node_remove(ref_node,1),"remove");

    RSS(ref_node_compact(ref_node,&o2n,&n2o),"compact");
 
    REIS(0,o2n[0],"o2n");
    REIS(REF_EMPTY,o2n[1],"o2n");
    REIS(1,o2n[2],"o2n");

    REIS(0,n2o[0],"n2o");
    REIS(2,n2o[1],"n2o");

    ref_free(n2o);
    ref_free(o2n);
    
    RSS(ref_node_free(ref_node),"free");
  }

  {  /* compact local nodes first */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT *o2n, *n2o;
    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,1,&node),"add");
    ref_node_part(ref_node,node) = ref_mpi_id+1;
    RSS(ref_node_add(ref_node,3,&node),"add");
    RSS(ref_node_add(ref_node,2,&node),"add");
    RSS(ref_node_remove(ref_node,1),"remove");

    RSS(ref_node_compact(ref_node,&o2n,&n2o),"compact");
 
    REIS(1,o2n[0],"o2n");
    REIS(REF_EMPTY,o2n[1],"o2n");
    REIS(0,o2n[2],"o2n");

    REIS(2,n2o[0],"n2o");
    REIS(0,n2o[1],"n2o");

    ref_free(n2o);
    ref_free(o2n);
    
    RSS(ref_node_free(ref_node),"free");
  }

  {  /* in bounding box */
    REF_INT node;
    REF_NODE ref_node;
    REF_INT n, *o2n, *n2o;
    REF_DBL bounding_box[6];
    RSS(ref_node_create(&ref_node),"create");

    RSS(ref_node_add(ref_node,1,&node),"add");
    ref_node_xyz(ref_node,0,node) = 0.0;
    ref_node_xyz(ref_node,1,node) = 0.0;
    ref_node_xyz(ref_node,2,node) = 0.0;
    RSS(ref_node_add(ref_node,2,&node),"add");
    ref_node_xyz(ref_node,0,node) = 0.0;
    ref_node_xyz(ref_node,1,node) = 1.0;
    ref_node_xyz(ref_node,2,node) = 0.0;
    RSS(ref_node_add(ref_node,3,&node),"add");
    ref_node_xyz(ref_node,0,node) = 0.0;
    ref_node_xyz(ref_node,1,node) = 0.0;
    ref_node_xyz(ref_node,2,node) = 1.0;

    bounding_box[0] = -0.5; bounding_box[1] =  0.5; 
    bounding_box[2] = -0.0; bounding_box[3] =  0.5; 
    bounding_box[4] = -0.5; bounding_box[5] =  1.5; 

    RSS(ref_node_in_bounding_box(ref_node,bounding_box,&n,&o2n,&n2o),"bbox");
 
    REIS(2,n,"o2n");

    REIS(0,o2n[0],"o2n");
    REIS(REF_EMPTY,o2n[1],"o2n");
    REIS(1,o2n[2],"o2n");

    REIS(0,n2o[0],"n2o");
    REIS(2,n2o[1],"n2o");

    ref_free(n2o);
    ref_free(o2n);
    
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

  { /* twod edge */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_BOOL twod;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&node0),"add");
    ref_node_xyz(ref_node,0,node0) = 0.0;
    ref_node_xyz(ref_node,1,node0) = 0.0;
    ref_node_xyz(ref_node,2,node0) = 0.0;
    global = 1;
    RSS(ref_node_add(ref_node,global,&node1),"add");
    ref_node_xyz(ref_node,0,node1) = 1.0;
    ref_node_xyz(ref_node,1,node1) = 0.0;
    ref_node_xyz(ref_node,2,node1) = 0.0;

    RSS( ref_node_edge_twod( ref_node, node0, node1, &twod ), "twod" ); 
    REIS( REF_TRUE, twod, "twod" );

    ref_node_xyz(ref_node,1,node1) = 1.0;
     
    RSS( ref_node_edge_twod( ref_node, node0, node1, &twod ), "twod" ); 
    REIS( REF_FALSE, twod, "twod" );

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

 #define FD_NODE0( xfuncx )				\
    {							\
      REF_DBL f, d[3];					\
      REF_DBL fd[3], x0, step = 1.0e-7, tol = 1.0e-6;	\
      REF_INT dir;					\
      RSS(xfuncx(ref_node,node0,node1,&f,d), "fd0");	\
      for ( dir=0;dir<3;dir++)				\
	{						\
	  x0 = ref_node_xyz(ref_node,dir,node0);	\
	  ref_node_xyz(ref_node,dir,node0)=		\
	    x0+step;					\
	  RSS(xfuncx(ref_node,node0,node1,		\
		     &(fd[dir]),d), "fd+");		\
	  fd[dir] = (fd[dir]-f)/step;			\
	  ref_node_xyz(ref_node,dir,node0) = x0;	\
	}						\
      RSS(xfuncx(ref_node,node0,node1,&f,d), "exact");	\
      RWDS( fd[0], d[0], tol, "dx expected" );		\
      RWDS( fd[1], d[1], tol, "dy expected" );		\
      RWDS( fd[2], d[2], tol, "dz expected" );		\
    }

  { /* derivative of node0 distance in metric */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio;
    REF_DBL f_ratio, d_ratio[3];

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

    /* same node */
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RSS( ref_node_dratio_dnode0(ref_node, node0, node1, 
			      &f_ratio, d_ratio), "ratio deriv" );
    RWDS( ratio, f_ratio, -1.0, "ratio expected" );
    RWDS( 0.0, d_ratio[0], -1.0, "dx expected" );
    RWDS( 0.0, d_ratio[1], -1.0, "dy expected" );
    RWDS( 0.0, d_ratio[2], -1.0, "dz expected" );

    /* length one in x */
    ref_node_xyz(ref_node,0,node1) = 1.0;

    FD_NODE0( ref_node_dratio_dnode0 );
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RSS( ref_node_dratio_dnode0(ref_node, node0, node1, 
				&f_ratio, d_ratio), "ratio deriv" );
    RWDS( ratio, f_ratio, -1.0, "ratio expected" );

    /* length one in xyz */
    ref_node_xyz(ref_node,0,node1) = 1.0;
    ref_node_xyz(ref_node,1,node1) = 1.0;
    ref_node_xyz(ref_node,2,node1) = 1.0;

    FD_NODE0( ref_node_dratio_dnode0 );
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RSS( ref_node_dratio_dnode0(ref_node, node0, node1, 
				&f_ratio, d_ratio), "ratio deriv" );
    RWDS( ratio, f_ratio, -1.0, "ratio expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* derivative of node0 distance in metric  gen */
    REF_NODE ref_node;
    REF_INT node0, node1, global;
    REF_DBL ratio, f_ratio, d_ratio[3];

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&node0),"add");
    ref_node_xyz(ref_node,0,node0) = 0.0;
    ref_node_xyz(ref_node,1,node0) = 0.0;
    ref_node_xyz(ref_node,2,node0) = 0.0;
    ref_node_metric(ref_node,0,node0) = 1.0;
    ref_node_metric(ref_node,1,node0) = 1.3;
    ref_node_metric(ref_node,2,node0) = 0.4;
    ref_node_metric(ref_node,3,node0) = 1.8;
    ref_node_metric(ref_node,4,node0) = 0.5;
    ref_node_metric(ref_node,5,node0) = 0.5;

    global = 1;
    RSS(ref_node_add(ref_node,global,&node1),"add");
    ref_node_xyz(ref_node,0,node1) = 0.6;
    ref_node_xyz(ref_node,1,node1) = 0.7;
    ref_node_xyz(ref_node,2,node1) = 0.8;
    ref_node_metric(ref_node,0,node1) = 1.0;
    ref_node_metric(ref_node,1,node1) = 0.0;
    ref_node_metric(ref_node,2,node1) = 0.0;
    ref_node_metric(ref_node,3,node1) = 1.0;
    ref_node_metric(ref_node,4,node1) = 0.0;
    ref_node_metric(ref_node,5,node1) = 1.0;

    FD_NODE0( ref_node_dratio_dnode0 );

    RSS( ref_node_dratio_dnode0(ref_node, node0, node1, 
			      &f_ratio, d_ratio), "ratio deriv" );
    RSS( ref_node_ratio(ref_node, node0, node1, &ratio), "ratio" );
    RWDS( ratio, f_ratio, -1.0, "ratio expected" );

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

#define FD_NODES0( xfuncx )				\
    {							\
      REF_DBL f, d[3];					\
      REF_DBL fd[3], x0, step = 1.0e-8, tol = 5.0e-6;	\
      REF_INT dir;					\
      RSS(xfuncx(ref_node,nodes,&f,d), "fd0");		\
      for ( dir=0;dir<3;dir++)				\
	{						\
	  x0 = ref_node_xyz(ref_node,dir,nodes[0]);	\
	  ref_node_xyz(ref_node,dir,nodes[0])= x0+step;	\
	  RSS(xfuncx(ref_node,nodes,&(fd[dir]),d),	\
	      "fd+");					\
	  fd[dir] = (fd[dir]-f)/step;			\
	  ref_node_xyz(ref_node,dir,nodes[0]) = x0;	\
	}						\
      RSS(xfuncx(ref_node,nodes,&f,d), "exact");	\
      RWDS( fd[0], d[0], tol, "dx expected" );		\
      RWDS( fd[1], d[1], tol, "dy expected" );		\
      RWDS( fd[2], d[2], tol, "dz expected" );		\
    }

  { /* tet volume quality deriv */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL f_vol, d_vol[3], vol;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");
    global = 3;
    RSS(ref_node_add(ref_node,global,&(nodes[3])),"add");

    ref_node_xyz(ref_node,0,nodes[0]) = 0.1;
    ref_node_xyz(ref_node,1,nodes[0]) = 0.2;
    ref_node_xyz(ref_node,2,nodes[0]) = 0.3;

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,1,nodes[1]) = 0.5;
    ref_node_xyz(ref_node,2,nodes[1]) = 0.7;

    ref_node_xyz(ref_node,0,nodes[2]) = 0.3;
    ref_node_xyz(ref_node,1,nodes[2]) = 1.4;
    ref_node_xyz(ref_node,2,nodes[2]) = 0.6;

    ref_node_xyz(ref_node,0,nodes[3]) = 0.2;
    ref_node_xyz(ref_node,1,nodes[3]) = 0.7;
    ref_node_xyz(ref_node,2,nodes[3]) = 1.9;

    FD_NODES0( ref_node_tet_dvol_dnode0 );

    RSS( ref_node_tet_dvol_dnode0(ref_node, nodes, 
				  &f_vol, d_vol), "ratio deriv" );
    RSS(ref_node_tet_vol(ref_node, nodes, &vol), "vol");
    RWDS( vol, f_vol, -1.0, "vol expected" );

    for ( global=0;global<4;global++)
      {
	ref_node_metric(ref_node,0,nodes[global]) = 20.0;
	ref_node_metric(ref_node,1,nodes[global]) = 0.0;
	ref_node_metric(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,3,nodes[global]) = 20.0;
	ref_node_metric(ref_node,4,nodes[global]) = 0.0;
	ref_node_metric(ref_node,5,nodes[global]) = 20.0;
       }

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;

    FD_NODES0( ref_node_tet_dquality_dnode0 );

    RSS( ref_node_tet_dquality_dnode0(ref_node, nodes, 
				      &f_quality, d_quality), "deriv" );
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    RWDS( quality, f_quality, -1.0, "vol expected" );

    /* test negative tet */
    ref_node_xyz(ref_node,2,nodes[3]) = -1.9;

    FD_NODES0( ref_node_tet_dquality_dnode0 );

    RSS( ref_node_tet_dquality_dnode0(ref_node, nodes, 
				      &f_quality, d_quality), "deriv" );
    RSS(ref_node_tet_quality(ref_node, nodes, &quality), "qual");
    RWDS( quality, f_quality, -1.0, "vol expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  /* FIXME, break this test up into pieces */

  { /* right tri normal, area */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL norm[3];
    REF_DBL area;
    REF_BOOL valid;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
       }

    ref_node_xyz(ref_node,0,nodes[1]) = 2.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 2.0;

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS( 2.0, area, -1.0, "expected area" );
    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(-2.0, area, -1.0, "expected area" );

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS( 0.5, area, -1.0, "expected area" );
    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS(-0.5, area, -1.0, "expected area" );
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS( !valid, "expected invalid" );

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS( 0.0, norm[0], -1.0, "nx" );
    RWDS(-1.0, norm[1], -1.0, "ny" );
    RWDS( 0.0, norm[2], -1.0, "nz" );

    global=nodes[2];
    nodes[2]=nodes[1];
    nodes[1]=global;

    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS( 0.5, area, -1.0, "expected area" );
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS( valid, "expected valid" );
    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS( 0.0, norm[0], -1.0, "nx" );
    RWDS( 1.0, norm[1], -1.0, "ny" );
    RWDS( 0.0, norm[2], -1.0, "nz" );

    ref_node_xyz(ref_node,0,nodes[1]) = 0.5;
    ref_node_xyz(ref_node,2,nodes[1]) = 0.0;

    RSS(ref_node_tri_y_projection(ref_node, nodes, &area), "area");
    RWDS( 0.0, area, -1.0, "expected area" );
    RSS(ref_node_tri_twod_orientation(ref_node, nodes, &valid), "valid");
    RAS( !valid, "expected zero area is invalid" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* right tri qual epic */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node),"create");
    ref_node->tri_quality = REF_NODE_EPIC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
       }
    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.5*sqrt(3.0), qual, -1.0, "qual expected" );

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 100.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 100.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 100.0;
       }
    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.5*sqrt(3.0), qual, -1.0, "qual expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* right tri qual jac */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL qual;

    RSS(ref_node_create(&ref_node),"create");
    ref_node->tri_quality = REF_NODE_JAC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
       }
    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,0,nodes[2]) = 0.5;
    ref_node_xyz(ref_node,1,nodes[2]) = 0.5*sqrt(3.0);

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 1.0, qual, -1.0, "qual expected" );

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
       }
    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.5*sqrt(3.0), qual, -1.0, "qual expected" );

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 100.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 100.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 100.0;
       }
    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.5*sqrt(3.0), qual, -1.0, "qual expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* tri area quality deriv */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL f_area, d_area[3], area;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_metric(ref_node,0,global) = 10.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 10.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 10.0;
       }

    ref_node_xyz(ref_node,0,nodes[0]) = 0.1;
    ref_node_xyz(ref_node,1,nodes[0]) = 0.2;
    ref_node_xyz(ref_node,2,nodes[0]) = 0.3;

    ref_node_xyz(ref_node,0,nodes[1]) = 1.1;
    ref_node_xyz(ref_node,1,nodes[1]) = 0.5;
    ref_node_xyz(ref_node,2,nodes[1]) = 0.7;

    ref_node_xyz(ref_node,0,nodes[2]) = 0.4;
    ref_node_xyz(ref_node,1,nodes[2]) = 2.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 0.6;

    FD_NODES0( ref_node_tri_darea_dnode0 );
    RSS( ref_node_tri_darea_dnode0(ref_node, nodes, 
				 &f_area, d_area), "area deriv" );
    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS( area, f_area, -1.0, "expected area" );

    /* quality */

    FD_NODES0( ref_node_tri_dquality_dnode0 );
    RSS( ref_node_tri_dquality_dnode0(ref_node, nodes, 
				      &f_quality, d_quality), "qual deriv" );
    RSS( ref_node_tri_quality(ref_node, nodes, 
			      &quality), "qual deriv" );
    RWDS( quality, f_quality, -1.0, "expected quality" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* tri quality jac deriv */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL f_quality, d_quality[3], quality;

    RSS(ref_node_create(&ref_node),"create");
    ref_node->tri_quality = REF_NODE_JAC_QUALITY;

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_metric(ref_node,0,global) = 14.0;
	ref_node_metric(ref_node,1,global) = 25.0;
	ref_node_metric(ref_node,2,global) = 40.0;
	ref_node_metric(ref_node,3,global) = 45.0;
	ref_node_metric(ref_node,4,global) = 71.0;
	ref_node_metric(ref_node,5,global) = 115.0;
       }

    ref_node_xyz(ref_node,0,nodes[0]) = 0.1;
    ref_node_xyz(ref_node,1,nodes[0]) = 0.2;
    ref_node_xyz(ref_node,2,nodes[0]) = 0.3;

    ref_node_xyz(ref_node,0,nodes[1]) = 1.1;
    ref_node_xyz(ref_node,1,nodes[1]) = 0.5;
    ref_node_xyz(ref_node,2,nodes[1]) = 0.7;

    ref_node_xyz(ref_node,0,nodes[2]) = 0.4;
    ref_node_xyz(ref_node,1,nodes[2]) = 2.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 0.6;

    /* quality */

    FD_NODES0( ref_node_tri_dquality_dnode0 );
    RSS( ref_node_tri_dquality_dnode0(ref_node, nodes, 
				      &f_quality, d_quality), "qual deriv" );
    RSS( ref_node_tri_quality(ref_node, nodes, 
			      &quality), "qual deriv" );
    RWDS( quality, f_quality, -1.0, "expected quality" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* equilateral tri normal, area, qual */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL norm[3];
    REF_DBL area, qual;

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
	ref_node_metric(ref_node,0,global) = 1.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 1.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 1.0;
       }

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,0,nodes[2]) = 0.5;
    ref_node_xyz(ref_node,1,nodes[2]) = 0.5*sqrt(3.0);

    RSS(ref_node_tri_area(ref_node, nodes, &area), "area");
    RWDS( 0.25*sqrt(3.0), area, -1.0, "expected area" );

    RSS(ref_node_tri_quality(ref_node, nodes, &qual), "q");
    RWDS( 1.0, qual, -1.0, "qual expected" );

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS( 0.0, norm[0], -1.0, "nx" );
    RWDS( 0.0, norm[1], -1.0, "ny" );
    RWDS( 0.5*sqrt(3.0), norm[2], -1.0, "nz" );

    global=nodes[2];
    nodes[2]=nodes[1];
    nodes[1]=global;

    RSS(ref_node_tri_normal(ref_node, nodes, norm), "norm");
    RWDS( 0.0, norm[0], -1.0, "nx" );
    RWDS( 0.0, norm[1], -1.0, "ny" );
    RWDS(-0.5*sqrt(3.0), norm[2], -1.0, "nz" );

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

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "epic qual expected" );

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "jac qual expected" );

    for ( global=0;global<4;global++)
      {
	ref_node_metric(ref_node,0,global) = 100.0;
	ref_node_metric(ref_node,1,global) = 0.0;
	ref_node_metric(ref_node,2,global) = 0.0;
	ref_node_metric(ref_node,3,global) = 100.0;
	ref_node_metric(ref_node,4,global) = 0.0;
	ref_node_metric(ref_node,5,global) = 100.0;
      }

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "qual expected not metric dep" );

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "jac qual expected" );

    /* inverted tet is negative volume */
    ref_node_xyz(ref_node,2,nodes[3]) = -1.0;
    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( -1.0/6.0, qual, -1.0, "qual expected" );
    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    /* mapped volume */
    RWDS( -1000.0/6.0, qual, -1.0, "qual expected" );

    RSS(ref_node_free(ref_node),"free");
  }

  if (REF_FALSE)
  { /* right tet jac quality sign */
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
	ref_node_metric(ref_node,0,global) = 5.569680186165702e+00;
	ref_node_metric(ref_node,1,global) = -1.669546513836191e-01;
	ref_node_metric(ref_node,2,global) = -5.647583642671333e-02;
	ref_node_metric(ref_node,3,global) = 5.645045675922971e+00;
	ref_node_metric(ref_node,4,global) = 1.230436787883563e+00;
	ref_node_metric(ref_node,5,global) = 2.881886210688973e+0;
      }

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,1,nodes[2]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[3]) = 1.0;

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.81577123, qual, 0.00001, "epic qual expected" );

    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 0.839947, qual, 0.00001, "jac qual expected" );

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

    ref_node->tet_quality = REF_NODE_EPIC_QUALITY;
    RSS(ref_node_tet_quality(ref_node, nodes, &qual), "q");
    RWDS( 1.0, qual, -1.0, "qual expected" );
    ref_node->tet_quality = REF_NODE_JAC_QUALITY;
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

  { /* clone a twod node  */
    REF_NODE ref_node;
    REF_INT node0, node1, global;

    RSS(ref_node_create(&ref_node),"create");

    ref_node_naux(ref_node) = 2;
    RSS(ref_node_resize_aux(ref_node),"resize aux");

    RSS(ref_node_next_global( ref_node, &global ), "next_global");
    RSS(ref_node_add(ref_node,global,&node0),"add");
    ref_node_xyz(ref_node,0,node0) = 0.1;
    ref_node_xyz(ref_node,1,node0) = 0.0;
    ref_node_xyz(ref_node,2,node0) = 0.3;
    ref_node_metric(ref_node,0,node0) = 1.1;
    ref_node_metric(ref_node,1,node0) = 0.2;
    ref_node_metric(ref_node,2,node0) = 0.3;
    ref_node_metric(ref_node,3,node0) = 1.4;
    ref_node_metric(ref_node,4,node0) = 0.5;
    ref_node_metric(ref_node,5,node0) = 1.6;

    ref_node_aux(ref_node,0,node0) =  1.0;
    ref_node_aux(ref_node,1,node0) = 20.0;

    RSS(ref_node_twod_clone(ref_node, node0, &node1),"clone");

    RWDS(  1.0, ref_node_aux(ref_node,0,node1), -1.0, "a");
    RWDS( 20.0, ref_node_aux(ref_node,1,node1), -1.0, "a");

    RWDS(  0.1, ref_node_xyz(ref_node,0,node1), -1.0, "x");
    RWDS(  1.0, ref_node_xyz(ref_node,1,node1), -1.0, "y");
    RWDS(  0.3, ref_node_xyz(ref_node,2,node1), -1.0, "z");

    RWDS(  1.4, ref_node_metric(ref_node,3,node1), -1.0, "m");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* interpolate aux */
    REF_NODE ref_node;
    REF_INT node0, node1, new_node, global;

    RSS(ref_node_create(&ref_node),"create");

    ref_node_naux(ref_node) = 2;
    RSS(ref_node_resize_aux(ref_node),"resize aux");

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

    ref_node_aux(ref_node,0,node0) =  1.0;
    ref_node_aux(ref_node,1,node0) = 20.0;


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

    ref_node_aux(ref_node,0,node1) =  2.0;
    ref_node_aux(ref_node,1,node1) = 40.0;

    RSS(ref_node_next_global( ref_node, &global ), "next_global");
    RSS(ref_node_add(ref_node,global,&new_node),"add");
    RSS(ref_node_interpolate_edge(ref_node, node0, node1, new_node),"interp");

    RWDS(  1.5, ref_node_aux(ref_node,0,new_node), -1.0, "a");
    RWDS( 30.0, ref_node_aux(ref_node,1,new_node), -1.0, "a");

    RSS(ref_node_free(ref_node),"free");
  }

  { /* twod bary */
    REF_NODE ref_node;
    REF_INT nodes[3], global;
    REF_DBL xyz[3];
    REF_DBL bary[3];

    RSS(ref_node_create(&ref_node),"create");

    global = 0;
    RSS(ref_node_add(ref_node,global,&(nodes[0])),"add");
    global = 1;
    RSS(ref_node_add(ref_node,global,&(nodes[1])),"add");
    global = 2;
    RSS(ref_node_add(ref_node,global,&(nodes[2])),"add");

    for ( global=0;global<3;global++)
      {
	ref_node_xyz(ref_node,0,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,1,nodes[global]) = 0.0;
	ref_node_xyz(ref_node,2,nodes[global]) = 0.0;
       }

    ref_node_xyz(ref_node,0,nodes[1]) = 1.0;
    ref_node_xyz(ref_node,2,nodes[2]) = 1.0;

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 1.0, bary[0], -1.0, "b0" );
    RWDS( 0.0, bary[1], -1.0, "b1" );
    RWDS( 0.0, bary[2], -1.0, "b2" );

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 1.0, bary[1], -1.0, "b1" );
    RWDS( 0.0, bary[2], -1.0, "b2" );

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 0.0, bary[1], -1.0, "b1" );
    RWDS( 1.0, bary[2], -1.0, "b2" );

    xyz[0] = 0.5;
    xyz[1] = 0.0;
    xyz[2] = 0.5;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 0.5, bary[1], -1.0, "b1" );
    RWDS( 0.5, bary[2], -1.0, "b2" );

    xyz[0] = 1.0/3.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0/3.0;

    RSS(ref_node_bary3(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 1.0/3.0, bary[0], -1.0, "b0" );
    RWDS( 1.0/3.0, bary[1], -1.0, "b1" );
    RWDS( 1.0/3.0, bary[2], -1.0, "b2" );

    RSS(ref_node_free(ref_node),"free");
  }

  { /* threed bary */
    REF_NODE ref_node;
    REF_INT nodes[4], global;
    REF_DBL xyz[3];
    REF_DBL bary[4];

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

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 1.0, bary[0], -1.0, "b0" );
    RWDS( 0.0, bary[1], -1.0, "b1" );
    RWDS( 0.0, bary[2], -1.0, "b2" );
    RWDS( 0.0, bary[3], -1.0, "b3" );

    xyz[0] = 1.0;
    xyz[1] = 0.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 1.0, bary[1], -1.0, "b1" );
    RWDS( 0.0, bary[2], -1.0, "b2" );
    RWDS( 0.0, bary[3], -1.0, "b3" );

    xyz[0] = 0.0;
    xyz[1] = 1.0;
    xyz[2] = 0.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 0.0, bary[1], -1.0, "b1" );
    RWDS( 1.0, bary[2], -1.0, "b2" );
    RWDS( 0.0, bary[3], -1.0, "b3" );

    xyz[0] = 0.0;
    xyz[1] = 0.0;
    xyz[2] = 1.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 0.0, bary[1], -1.0, "b1" );
    RWDS( 0.0, bary[2], -1.0, "b2" );
    RWDS( 1.0, bary[3], -1.0, "b3" );

    xyz[0] = 1.0/3.0;
    xyz[1] = 1.0/3.0;
    xyz[2] = 1.0/3.0;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.0, bary[0], -1.0, "b0" );
    RWDS( 1.0/3.0, bary[1], -1.0, "b1" );
    RWDS( 1.0/3.0, bary[2], -1.0, "b2" );
    RWDS( 1.0/3.0, bary[3], -1.0, "b3" );

    xyz[0] = 0.25;
    xyz[1] = 0.25;
    xyz[2] = 0.25;

    RSS(ref_node_bary4(ref_node, nodes, xyz, bary ), "bary");
    RWDS( 0.25, bary[0], -1.0, "b0" );
    RWDS( 0.25, bary[1], -1.0, "b1" );
    RWDS( 0.25, bary[2], -1.0, "b2" );
    RWDS( 0.25, bary[3], -1.0, "b3" );

    RSS(ref_node_free(ref_node),"free");
  }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
