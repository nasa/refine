
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
  RSS( ref_metric_create( &ref_grid_metric(ref_grid) ), "metric create" );

  RSS( ref_cell_create( 4, &ref_grid_tet(ref_grid) ), "tet create" );
  RSS( ref_cell_create( 5, &ref_grid_pyr(ref_grid) ), "pyr create" );
  RSS( ref_cell_create( 6, &ref_grid_pri(ref_grid) ), "pri create" );
  RSS( ref_cell_create( 8, &ref_grid_hex(ref_grid) ), "hex create" );

  ref_grid->cell[4] = NULL;

  RSS( ref_cell_create( 4, &ref_grid_tri(ref_grid) ), "tri create" );
  RSS( ref_cell_create( 5, &ref_grid_qua(ref_grid) ), "qua create" );

  ref_grid_nedge(ref_grid) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_grid_free( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_node_free( ref_grid_node(ref_grid) ), "node free");
  RSS( ref_metric_free( ref_grid_metric(ref_grid) ), "metric free");

  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");

  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");

  ref_cond_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_inspect( REF_GRID ref_grid )
{
  printf("ref_grid = %p\n",(void *)ref_grid);
  printf(" tet = %p\n",(void *)ref_grid_tet(ref_grid));
  printf(" pyr = %p\n",(void *)ref_grid_pyr(ref_grid));
  printf(" pri = %p\n",(void *)ref_grid_pri(ref_grid));
  printf(" hex = %p\n",(void *)ref_grid_hex(ref_grid));

  return REF_SUCCESS;
}

REF_STATUS ref_grid_make_edges( REF_GRID ref_grid )
{
  REF_CELL ref_cell, ref_cell2;
  REF_INT group, group2;
  REF_INT nedge;
  REF_INT cell, cell_edge;
  REF_INT n0, n1;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    ref_cell_empty_edges( ref_cell);

  nedge = 0;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_edge( ref_cell, cell_edge )
        if ( REF_EMPTY == ref_cell_c2e(ref_cell,cell_edge,cell) )
	    {
	      n0 = ref_cell_e2n(ref_cell,0,cell,cell_edge);
	      n1 = ref_cell_e2n(ref_cell,1,cell,cell_edge);
	      each_ref_grid_ref_cell( ref_grid, group2, ref_cell2 )
		RSS( ref_cell_set_edge( ref_cell2, n0, n1, nedge ), "set edge");
	      nedge++;
	    }

  ref_grid_nedge(ref_grid) = nedge;

  return REF_SUCCESS;
}
