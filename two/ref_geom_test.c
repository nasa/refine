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
      RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );
      RSS( ref_grid_free(ref_grid),"free");
    }

  REIS(REF_NULL,ref_geom_free(NULL),"dont free NULL");

  { /* create and destroy */
    REF_GEOM ref_geom;
    RSS(ref_geom_create(&ref_geom),"create");
    REIS( 0, ref_geom_n(ref_geom), "init no nodes" );
    RSS(ref_geom_free(ref_geom),"free");
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
  
  return 0;
}