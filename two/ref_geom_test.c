#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "ref_geom.h"

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_sort.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include "ref_fixture.h"
#include "ref_export.h"
#include  "ref_dict.h"
#include  "ref_edge.h"

int main( int argc, char *argv[] )
{

  if ( 2 == argc )
    { /* fixture */
      char *filename = argv[1];
      if ( 0 == access( filename, R_OK ) )
	{
	  printf("EGADS project %s exisits, deleting\n",filename);
	  REIS(0, remove( filename ), "test clean up");
	}
      RSS( ref_geom_egads_fixture( filename ), "egads fixture" );
      printf("wrote EGADS project %s\n",filename);
    }

  if ( 3 == argc )
    { /* egads to grid */
      REF_GRID ref_grid;
      RSS( ref_geom_grid_from_egads( &ref_grid, argv[1] ), "from egads" );
      RSS( ref_export_by_extension( ref_grid, argv[2] ), "argv export" );
      RSS( ref_geom_tec( ref_grid, "ref_geom_test.tec" ), "geom export" );
      RSS( ref_grid_free(ref_grid),"free");
    }

  REIS(REF_NULL,ref_geom_free(NULL),"dont free NULL");

  { /* create and destroy */
    REF_GEOM ref_geom;
    RSS(ref_geom_create(&ref_geom),"create");
    REIS( 0, ref_geom_n(ref_geom), "init no nodes" );
    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* deep copy empty */
    REF_GEOM ref_geom;
    REF_GEOM original;

    RSS(ref_geom_create(&original),"create");
    RSS(ref_geom_deep_copy(&ref_geom,original),"deep copy");
    REIS( ref_geom_n(original), ref_geom_n(ref_geom), "items" );
    REIS( ref_geom_max(original), ref_geom_max(ref_geom), "items" );	
    RSS(ref_geom_free(original),"cleanup");
    RSS(ref_geom_free(ref_geom),"cleanup");
  }

  { /* add geom node */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL *params;
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2;
    type = REF_GEOM_NODE;
    id = 5;
    params = NULL;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add node" );
    REIS( 1, ref_geom_n(ref_geom), "items" );
    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* add geom edge */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL params[1];
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2;
    type = REF_GEOM_EDGE;
    id = 5;
    params[0] = 11.0;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add edge" );
    REIS( 1, ref_geom_n(ref_geom), "items" );
    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* add geom face */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL params[2];
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2;
    type = REF_GEOM_FACE;
    id = 5;
    params[0] = 11.0;
    params[1] = 21.0;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
    REIS( 1, ref_geom_n(ref_geom), "items" );
    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* add updates face uv*/
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL params[2], uv[2];
    REF_DBL tol = 1.0e-14;
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2;
    type = REF_GEOM_FACE;
    id = 5;
    params[0] = 11.0; params[1] = 21.0;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
    REIS( 0, ref_geom_tuv(ref_geom,node,type,id,uv), "face uv" );
    RWDS( params[0], uv[0], tol, "u" );
    RWDS( params[1], uv[1], tol, "v" );
    params[0] = 12.0; params[1] = 22.0;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
    REIS( 0, ref_geom_tuv(ref_geom,node,type,id,uv), "face uv" );
    RWDS( params[0], uv[0], tol, "u" );
    RWDS( params[1], uv[1], tol, "v" );
    REIS( 1, ref_geom_n(ref_geom), "items" );
    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* add and remove */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL params[2];
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2; type = REF_GEOM_FACE; id = 5; params[0] = 11.0; params[1] = 21.0;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
    node = 4; type = REF_GEOM_NODE; id = 2;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add node" );
    REIS( 2, ref_geom_n(ref_geom), "items" );
    node = 4; type = REF_GEOM_EDGE; id = 2;
    REIS( REF_NOT_FOUND, ref_geom_remove(ref_geom,node,type,id),
	  "should not remove missing edge" );
    REIS( 2, ref_geom_n(ref_geom), "items" );
    node = 4; type = REF_GEOM_NODE; id = 2;
    REIS( 0, ref_geom_remove(ref_geom,node,type,id),
	  "remove node" );
    REIS( 1, ref_geom_n(ref_geom), "items" );
    REIS( REF_NOT_FOUND, ref_geom_remove(ref_geom,node,type,id),
	  "not really gone" );
    RSS(ref_geom_free(ref_geom),"free");
  }

  { /* reuse without reallocation */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_INT geom, max;
    REF_DBL params[2];
    RSS(ref_geom_create(&ref_geom),"create");
    max = ref_geom_max(ref_geom);
    
    for ( geom = 0 ; geom < max+10 ; geom++ )
      {
	node = 2; type = REF_GEOM_FACE; id = 5; params[0] =1.0; params[1] =2.0;
	REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
	REIS( 0, ref_geom_remove(ref_geom,node,type,id), "rm" );
      }
    REIS( max, ref_geom_max(ref_geom), "items" );
    RSS(ref_geom_free(ref_geom),"free");
  }

  { /* force realloc twice */
    REF_GEOM ref_geom;
    REF_INT max, i;
    REF_INT node, type, id;
    REF_DBL params[2];

    RSS(ref_geom_create(&ref_geom),"create");
    
    max = ref_geom_max(ref_geom);
    for (i = 0; i < max+1; i++ )
      {
	node = i; type = REF_GEOM_FACE; id = 5; params[0] =1.0; params[1] =2.0;
	REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
      }
    RAS(ref_geom_max(ref_geom)>max,"realloc max");

    max = ref_geom_max(ref_geom);
    for (i = ref_geom_n(ref_geom); i < max+1; i++ )
      {
	node = i; type = REF_GEOM_FACE; id = 5; params[0] =1.0; params[1] =2.0;
	REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add face" );
      }
    RAS(ref_geom_max(ref_geom)>max,"realloc max");

    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* evaluate common edge */
    REF_GEOM ref_geom;
    REF_INT node0, node1, new_node;
    REF_INT type, id;
    REF_DBL params[2];
    REF_DBL tol =1.0e-12;

    RSS(ref_geom_create(&ref_geom),"create");

    node0 = 0; node1 = 1; new_node = 10;

    type = REF_GEOM_EDGE; id = 5; params[0] = 11.0;
    RSS( ref_geom_add(ref_geom,node0,type,id,params), "node0 edge" );
    type = REF_GEOM_EDGE; id = 5; params[0] = 13.0;
    RSS( ref_geom_add(ref_geom,node1,type,id,params), "node1 edge" );
    params[0] = 0.0;
    
    RSS( ref_geom_add_between(ref_geom,node0,node1,new_node), "eval" );
    RSS( ref_geom_tuv(ref_geom,new_node,type,id,params), "node1 edge" );
    RWDS( 12.0, params[0], tol, "v" );

    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* evaluate skip different edge */
    REF_GEOM ref_geom;
    REF_INT node0, node1, new_node;
    REF_INT type, id, geom;
    REF_DBL params[2];

    RSS(ref_geom_create(&ref_geom),"create");

    node0 = 0; node1 = 1; new_node = 10;

    type = REF_GEOM_EDGE; id = 5; params[0] = 11.0;
    RSS( ref_geom_add(ref_geom,node0,type,id,params), "node0 edge" );
    type = REF_GEOM_EDGE; id = 7; params[0] = 13.0;
    RSS( ref_geom_add(ref_geom,node1,type,id,params), "node1 edge" );
    params[0] = 0.0;
    
    RSS( ref_geom_add_between(ref_geom,node0,node1,new_node), "eval" );
    REIS( REF_NOT_FOUND,
	  ref_geom_find(ref_geom,new_node,type,id,&geom), "what edge" );
    id = 5;
    REIS( REF_NOT_FOUND,
	  ref_geom_find(ref_geom,new_node,type,id,&geom), "what edge" );

    RSS(ref_geom_free(ref_geom),"free");
  }
  
  return 0;
}
