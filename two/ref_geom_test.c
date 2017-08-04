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
#include "ref_import.h"
#include "ref_export.h"
#include  "ref_dict.h"
#include  "ref_edge.h"

#include  "ref_math.h"
#include  "ref_args.h"

int main( int argc, char *argv[] )
{
  REF_INT assoc_pos = REF_EMPTY;
  RXS( ref_args_find( argc, argv, "--assoc", &assoc_pos ),
       REF_NOT_FOUND, "arg search" );

  if ( assoc_pos != REF_EMPTY )
    {
      REF_GRID ref_grid;
      REIS( 5, argc,
	    "required args: --assoc grid.ext input.gas grid_geom_assoc.meshb");
      REIS( 1, assoc_pos,
	    "required args: --assoc grid.ext input.gas grid_geom_assoc.meshb");
      printf("merge geometry association into meshb\n");
      printf("grid source %s\n",argv[2]);
      printf("geometry association source %s\n",argv[3]);
      printf("output %s\n",argv[4]);
      RSS( ref_import_by_extension( &ref_grid, argv[2] ), "argv import" );
      RSS( ref_geom_load( ref_grid, argv[3] ), "geom gas import" );
      RSS( ref_export_by_extension( ref_grid, argv[4] ), "argv export" );
      return 0;
    }
  
  if ( 2 == argc )
    { /* fixture */
      char *filename = argv[1];
      if ( 0 == access( filename, R_OK ) )
	{
	  printf("EGADS project %s exists, deleting\n",filename);
	  REIS(0, remove( filename ), "test clean up");
	}
      RSS( ref_geom_egads_export( filename ), "egads fixture" );
      printf("wrote EGADS project %s\n",filename);
    }

  if ( 3 == argc || 4 == argc )
    { /* egads to grid */
      REF_GRID ref_grid;
      REF_INT node;
      REF_INT nedge;
      
      REF_DBL max_edge = -0.25;
      if ( 4 == argc ) max_edge = atof( argv[3] );

      RSS(ref_grid_create(&ref_grid),"create");

      RSS(ref_geom_egads_load(ref_grid_geom(ref_grid), argv[1] ), "ld egads" );
      RSS(ref_geom_egads_tess( ref_grid, max_edge ), "tess egads" );
      RSS(ref_geom_tetgen_volume(ref_grid ), "tetgen surface to volume " );

      RSS( ref_export_by_extension( ref_grid, argv[2] ), "argv export" );
      RSS( ref_geom_tec( ref_grid, "ref_geom_test.tec" ), "geom export" );
      RSS( ref_geom_verify_param( ref_grid ), "original params" );
      printf("constrain\n");
      each_ref_node_valid_node( ref_grid_node(ref_grid), node )
	RSS( ref_geom_constrain( ref_grid, node ), "original params" );
      printf("verify\n");
      RSS( ref_geom_verify_param( ref_grid ), "constrained params" );
      nedge = ref_cell_n(ref_grid_edg(ref_grid));
      printf("save %d edge\n",nedge);
      RSS( ref_geom_save( ref_grid, "ref_geom_test.gas" ), "save" );
      printf("clear\n");
      each_ref_node_valid_node( ref_grid_node(ref_grid), node )
	RSS( ref_geom_remove_all( ref_grid_geom(ref_grid), node ), "erasure" );
      RSS( ref_cell_free( ref_grid_edg(ref_grid) ), "free edge");
      RSS( ref_cell_create( &ref_grid_edg(ref_grid), 2, REF_TRUE ), "new edg" );
      printf("load\n");
      RSS( ref_geom_load( ref_grid, "ref_geom_test.gas" ), "load" );
      REIS( nedge, ref_cell_n(ref_grid_edg(ref_grid)), "nedge");
      printf("verify\n");
      RSS( ref_geom_verify_param( ref_grid ), "loaded params" );
      printf("constrain\n");
      each_ref_node_valid_node( ref_grid_node(ref_grid), node )
	RSS( ref_geom_constrain( ref_grid, node ), "original params" );
      printf("verify\n");
      RSS( ref_geom_verify_param( ref_grid ), "constrained params" );
      RSS( ref_grid_free(ref_grid),"free");
      /* REIS(0, remove( "ref_geom_test.gas" ), "test clean up"); */
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

  { /* add and remove all */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL params[2];
    RSS(ref_geom_create(&ref_geom),"create");
    node = 2; type = REF_GEOM_FACE; id = 5; params[0] = 11.0; params[1] = 21.0;
    RSS( ref_geom_add(ref_geom,node,type,id,params), "add face" );
    node = 2; type = REF_GEOM_EDGE; id = 2; params[0] = 5.0;
    RSS( ref_geom_add(ref_geom,node,type,id,params), "add edge" );
    REIS( 2, ref_geom_n(ref_geom), "items" );
    node = 4; type = REF_GEOM_EDGE; id = 2; params[0] = 5.0;
    RSS( ref_geom_remove_all(ref_geom,node), "ok with nothing there" );
    REIS( 2, ref_geom_n(ref_geom), "items" );
    node = 2;
    RSS( ref_geom_remove_all(ref_geom,node), "remove all at node" );
    REIS( 0, ref_geom_n(ref_geom), "items" );
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
  
  { /* add between common edge */
    REF_GRID ref_grid;
    REF_GEOM ref_geom;
    REF_INT node0, node1, new_node;
    REF_INT type, id, geom;
    REF_DBL params[2];
    REF_DBL tol =1.0e-12;
    REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
    RSS(ref_grid_create(&ref_grid),"create");
    ref_geom = ref_grid_geom(ref_grid);
 
    node0 = 0; node1 = 1; new_node = 10;

    type = REF_GEOM_EDGE; id = 5; params[0] = 11.0;
    RSS( ref_geom_add(ref_geom,node0,type,id,params), "node0 edge" );
    type = REF_GEOM_EDGE; id = 5; params[0] = 13.0;
    RSS( ref_geom_add(ref_geom,node1,type,id,params), "node1 edge" );
    params[0] = 0.0;
    
    RSS( ref_geom_add_between(ref_grid,node0,node1,new_node), "between" );
    id = 5;
    REIS( REF_NOT_FOUND,
	  ref_geom_find(ref_geom,new_node,type,id,&geom), "what edge" );

    nodes[0]=node0;nodes[1]=node1;nodes[2]=id;
    RSS(ref_cell_add(ref_grid_edg(ref_grid),nodes,&cell),"add edge");
    RSS( ref_geom_add_between(ref_grid,node0,node1,new_node), "between" );
    RSS( ref_geom_tuv(ref_geom,new_node,type,id,params), "node1 edge" );
    RWDS( 12.0, params[0], tol, "v" );

    RSS(ref_grid_free(ref_grid),"free");
  }
  
  { /* add between skip different edge */
    REF_GRID ref_grid;
    REF_GEOM ref_geom;
    REF_INT node0, node1, new_node;
    REF_INT type, id, geom;
    REF_DBL params[2];

    RSS(ref_grid_create(&ref_grid),"create");
    ref_geom = ref_grid_geom(ref_grid);

    node0 = 0; node1 = 1; new_node = 10;

    type = REF_GEOM_EDGE; id = 5; params[0] = 11.0;
    RSS( ref_geom_add(ref_geom,node0,type,id,params), "node0 edge" );
    type = REF_GEOM_EDGE; id = 7; params[0] = 13.0;
    RSS( ref_geom_add(ref_geom,node1,type,id,params), "node1 edge" );
    params[0] = 0.0;
    
    RSS( ref_geom_add_between(ref_grid,node0,node1,new_node), "between" );
    id = 7;
    REIS( REF_NOT_FOUND,
	  ref_geom_find(ref_geom,new_node,type,id,&geom), "what edge" );
    id = 5;
    REIS( REF_NOT_FOUND,
	  ref_geom_find(ref_geom,new_node,type,id,&geom), "what edge" );

    RSS(ref_grid_free(ref_grid),"free");
  }
  
  { /* eval out of range */
    REF_GEOM ref_geom;
    REF_INT geom;
    REF_DBL xyz[3];

    RSS(ref_geom_create(&ref_geom),"create");

    geom = -1;
    REIS( REF_INVALID,
	  ref_geom_eval(ref_geom,geom,xyz,NULL), "invalid geom" );
    geom = 1+ref_geom_max(ref_geom);
    REIS( REF_INVALID,
	  ref_geom_eval(ref_geom,geom,xyz,NULL), "invalid geom" );

    RSS(ref_geom_free(ref_geom),"free");
  }

  { /* constrain without geom or on node */
    REF_GRID ref_grid;
    REF_GEOM ref_geom;
    REF_INT node;
    REF_INT type, id;
    REF_DBL params[2];

    RSS(ref_grid_create(&ref_grid),"create");
    ref_geom = ref_grid_geom(ref_grid);
    
    node = 4;
    RSS( ref_geom_constrain(ref_grid,node), "no geom" );
    node = 4; type = REF_GEOM_NODE; id = 2;
    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add node" );
    RSS( ref_geom_constrain(ref_grid,node), "no geom" );

    RSS(ref_grid_free(ref_grid),"free");
  }
  
  { /* has_support */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL *params;
    REF_BOOL has_support;

    RSS(ref_geom_create(&ref_geom),"create");

    node = 2;
    type = REF_GEOM_NODE;
    id = 5;
    params = NULL;

    RSS( ref_geom_supported( ref_geom, node, &has_support ), "empty" );
    RAS( REF_FALSE == has_support, "phantom support" );

    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add node" );

    RSS( ref_geom_supported( ref_geom, node, &has_support ), "empty" );
    RAS( REF_TRUE == has_support, "missing support" );

    RSS(ref_geom_free(ref_geom),"free");
  }

  { /* is a */
    REF_GEOM ref_geom;
    REF_INT node, type, id;
    REF_DBL *params;
    REF_BOOL it_is;

    RSS(ref_geom_create(&ref_geom),"create");

    node = 2;
    type = REF_GEOM_NODE;
    id = 5;
    params = NULL;

    RSS( ref_geom_is_a( ref_geom, node, REF_GEOM_NODE, &it_is ), "empty" );
    RAS( REF_FALSE == it_is, "expected nothing" );

    REIS( 0, ref_geom_add(ref_geom,node,type,id,params), "add node" );

    RSS( ref_geom_is_a( ref_geom, node, REF_GEOM_NODE, &it_is ), "has node" );
    RAS( REF_TRUE == it_is, "expected node" );
    RSS( ref_geom_is_a( ref_geom, node, REF_GEOM_FACE, &it_is ), "empty face" );
    RAS( REF_FALSE == it_is, "expected no face" );

    RSS(ref_geom_free(ref_geom),"free");
  }

  { /* determine unique id */
    REF_GEOM ref_geom;
    REF_INT node;
    REF_INT type, id;
    REF_DBL params[2];

    RSS(ref_geom_create(&ref_geom),"create");

    type = REF_GEOM_EDGE;
    node = 0;
    
    REIS( REF_NOT_FOUND, ref_geom_unique_id(ref_geom,node,type,&id),
	  "found nothing" );

    id = 5; params[0] = 15.0;
    RSS( ref_geom_add(ref_geom,node,type,id,params), "node edge 5" );

    id = REF_EMPTY;
    RSS( ref_geom_unique_id(ref_geom,node,type,&id),
	 "found it" );
    REIS( 5, id, "expected 5" );
    
    id = 7; params[0] = 17.0;
    RSS( ref_geom_add(ref_geom,node,type,id,params), "node edge 7" );
    REIS( REF_INVALID, ref_geom_unique_id(ref_geom,node,type,&id),
	  "two found" );

    RSS(ref_geom_free(ref_geom),"free");
  }
  
  { /* system transform */
    REF_DBL tol = 1.0e-12;
    REF_DBL duv[6];
    REF_DBL r[3], s[3], n[3], drsduv[4];
    duv[0] = 2.0;
    duv[1] = 0.0;
    duv[2] = 0.0;
    duv[3] = 3.0;
    duv[4] = 4.0;
    duv[5] = 0.0;
    
    RSS( ref_geom_uv_rsn( duv, r, s, n, drsduv ), "make orthog" );
    RWDS( 1.0, r[0], tol, "r[0]" );
    RWDS( 0.0, r[1], tol, "r[1]" );
    RWDS( 0.0, r[2], tol, "r[2]" );
    RWDS( 0.0, s[0], tol, "s[0]" );
    RWDS( 1.0, s[1], tol, "s[1]" );
    RWDS( 0.0, s[2], tol, "s[2]" );
    RWDS( 0.0, n[0], tol, "n[0]" );
    RWDS( 0.0, n[1], tol, "n[1]" );
    RWDS( 1.0, n[2], tol, "n[2]" );

    RWDS( 0.5,   drsduv[0], tol, "drdu" );
    RWDS( 0.0,   drsduv[1], tol, "drdv" );
    RWDS(-0.375, drsduv[2], tol, "dsdu" );
    RWDS( 0.25,  drsduv[3], tol, "dsdv" );
  }

  return 0;
}

