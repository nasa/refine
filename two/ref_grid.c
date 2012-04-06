
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid.h"
#include "ref_adj.h"

#include "ref_malloc.h"
#include "ref_matrix.h"

REF_STATUS ref_grid_create( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  RSS( ref_node_create( &ref_grid_node(ref_grid) ), "node create" );

  RSS( ref_cell_create( &ref_grid_tet(ref_grid), 4, REF_FALSE ), "tet create" );
  RSS( ref_cell_create( &ref_grid_pyr(ref_grid), 5, REF_FALSE ), "pyr create" );
  RSS( ref_cell_create( &ref_grid_pri(ref_grid), 6, REF_FALSE ), "pri create" );
  RSS( ref_cell_create( &ref_grid_hex(ref_grid), 8, REF_FALSE ), "hex create" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_create( &ref_grid_tri(ref_grid), 3, REF_TRUE ), "tri create" );
  RSS( ref_cell_create( &ref_grid_qua(ref_grid), 4, REF_TRUE ), "qua create" );

  return REF_SUCCESS;
}

REF_STATUS ref_grid_empty_cell_clone( REF_GRID *ref_grid_ptr, REF_GRID parent )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  ref_grid_node(ref_grid) = ref_grid_node(parent);

  RSS( ref_cell_create( &ref_grid_tet(ref_grid), 4, REF_FALSE ), "tet create" );
  RSS( ref_cell_create( &ref_grid_pyr(ref_grid), 5, REF_FALSE ), "pyr create" );
  RSS( ref_cell_create( &ref_grid_pri(ref_grid), 6, REF_FALSE ), "pri create" );
  RSS( ref_cell_create( &ref_grid_hex(ref_grid), 8, REF_FALSE ), "hex create" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_create( &ref_grid_tri(ref_grid), 3, REF_TRUE ), "tri create" );
  RSS( ref_cell_create( &ref_grid_qua(ref_grid), 4, REF_TRUE ), "qua create" );

  return REF_SUCCESS;
}

REF_STATUS ref_grid_free( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_node_free( ref_grid_node(ref_grid) ), "node free");

  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");

  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");

  ref_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_free_cell_clone( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");

  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");

  ref_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_inspect( REF_GRID ref_grid )
{
  printf("ref_grid = %p\n",(void *)ref_grid);
  printf(" %d node\n",ref_node_n(ref_grid_node(ref_grid)));
  printf(" %d tet\n",ref_cell_n(ref_grid_tet(ref_grid)));
  printf(" %d pyr\n",ref_cell_n(ref_grid_pyr(ref_grid)));
  printf(" %d pri\n",ref_cell_n(ref_grid_pri(ref_grid)));
  printf(" %d hex\n",ref_cell_n(ref_grid_hex(ref_grid)));
  printf(" %d tri\n",ref_cell_n(ref_grid_tri(ref_grid)));
  printf(" %d qua\n",ref_cell_n(ref_grid_qua(ref_grid)));

  return REF_SUCCESS;
}

REF_STATUS ref_grid_imply_metric( REF_GRID ref_grid, REF_DBL *metric )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT *hits;
  REF_DBL m[6], log_m[6];
  REF_INT node, im;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  REF_CELL ref_cell = ref_grid_tet(ref_grid);

  ref_malloc_init( hits, ref_node_max(ref_node),
		   REF_INT, 0 );

  for( node=0; node<ref_node_max(ref_node); node++ )
    for( im=0; im<6; im++ )
      metric[im+6*node] = 0.0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_matrix_imply_m( m, 
			       ref_node_xyz_ptr(ref_node,nodes[0]), 
			       ref_node_xyz_ptr(ref_node,nodes[1]), 
			       ref_node_xyz_ptr(ref_node,nodes[2]), 
			       ref_node_xyz_ptr(ref_node,nodes[3])), "imply" );
      RSS( ref_matrix_log_m( m, log_m ), "log" );

      for( node=0; node<4; node++ )
	{
	  hits[nodes[node]]++;
	  for( im=0; im<6; im++ )
	    metric[im+6*nodes[node]] += log_m[im];
	}
    }

  each_ref_node_valid_node( ref_node, node )
    {
      RAS( 0 != hits[node], "zero metric contributions" );
      for( im=0; im<6; im++ )
	log_m[im] = metric[im+6*node] / ((REF_DBL)hits[node]);
      RSS( ref_matrix_exp_m( log_m, m ), "exp" );
      for( im=0; im<6; im++ )
	metric[im+6*node] = m[im];
    }

  ref_free( hits );

  return REF_SUCCESS;
}
