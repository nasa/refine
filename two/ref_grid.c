
#include <stdlib.h>
#include <stdio.h>

#include "ref_grid.h"

REF_STATUS ref_grid_create( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;
  (*ref_grid_ptr) = NULL;
  (*ref_grid_ptr) = (REF_GRID)malloc( sizeof(REF_GRID_STRUCT) );
  RNS(*ref_grid_ptr,"malloc ref_grid NULL");

  ref_grid = *ref_grid_ptr;

  RSS( ref_node_create( &ref_grid_node(ref_grid) ), "node create" );

  RSS( ref_cell_create( 4, &ref_grid_tet(ref_grid) ), "tet create" );
  RSS( ref_cell_create( 5, &ref_grid_pyr(ref_grid) ), "pyr create" );
  RSS( ref_cell_create( 6, &ref_grid_pri(ref_grid) ), "pri create" );
  RSS( ref_cell_create( 8, &ref_grid_hex(ref_grid) ), "hex create" );

  RSS( ref_cell_create( 4, &ref_grid_tri(ref_grid) ), "tri create" );
  RSS( ref_cell_create( 5, &ref_grid_qua(ref_grid) ), "qua create" );

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

  ref_cond_free( ref_grid );
  return REF_SUCCESS;
}
