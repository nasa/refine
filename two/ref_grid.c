
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

REF_STATUS ref_grid_cell_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 4:
      *ref_cell = ref_grid_tet(ref_grid);
      break;
    case 5:
      *ref_cell = ref_grid_pyr(ref_grid);
      break;
    case 6:
      *ref_cell = ref_grid_pri(ref_grid);
      break;
    case 8:
      *ref_cell = ref_grid_hex(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_face_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 3:
      *ref_cell = ref_grid_tri(ref_grid);
      break;
    case 4:
      *ref_cell = ref_grid_qua(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_cell_has_face( REF_GRID ref_grid, 
				   REF_INT *face_nodes,
				   REF_BOOL *has_face )
{
  REF_CELL ref_cell;
  REF_INT group;

  *has_face = REF_FALSE;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      RSS( ref_cell_has_face( ref_cell, face_nodes, has_face),
	   "cell has face" );
      if ( *has_face ) return REF_SUCCESS;
    }

 return REF_SUCCESS;
}
