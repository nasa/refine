
#include <stdlib.h>
#include <stdio.h>

#include "ref_edge.h"

#include "ref_malloc.h"

REF_STATUS ref_edge_create( REF_EDGE *ref_edge_ptr, REF_GRID ref_grid )
{
  REF_EDGE ref_edge;
  REF_INT edge;
  REF_INT group, group2, cell, cell_edge;
  REF_INT n0, n1;
  REF_CELL ref_cell, ref_cell2;

  ref_malloc( *ref_edge_ptr, 1, REF_EDGE_STRUCT );

  ref_edge = *ref_edge_ptr;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    ref_cell_empty_edges( ref_cell);

  ref_edge_n(ref_edge) = 0;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_edge( ref_cell, cell_edge )
        if ( REF_EMPTY == ref_cell_c2e(ref_cell,cell_edge,cell) )
	    {
	      n0 = ref_cell_e2n(ref_cell,0,cell,cell_edge);
	      n1 = ref_cell_e2n(ref_cell,1,cell,cell_edge);
	      each_ref_grid_ref_cell( ref_grid, group2, ref_cell2 )
		RSS( ref_cell_set_edge( ref_cell2, n0, n1, 
					ref_edge_n(ref_edge) ), "set edge");
	      ref_edge_n(ref_edge)++;
	    }
  ref_malloc( ref_edge->e2n, 2*ref_edge_n(ref_edge), REF_INT );

  for ( edge=0 ; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      ref_edge_e2n( ref_edge, edge, 0 ) = REF_EMPTY;
      ref_edge_e2n( ref_edge, edge, 1 ) = REF_EMPTY;
    }

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_edge( ref_cell, cell_edge )
        {
	  edge = ref_cell_c2e(ref_cell,cell_edge,cell);
	  ref_edge_e2n( ref_edge, edge, 0 ) = 
	    ref_cell_e2n(ref_cell,0,cell,cell_edge);
	  ref_edge_e2n( ref_edge, edge, 1 ) = 
	    ref_cell_e2n(ref_cell,1,cell,cell_edge);
	}

  RSS( ref_adj_create( &(ref_edge_adj( ref_edge )) ), "create adj" );

  for ( edge=0 ; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      RUS(REF_EMPTY,ref_edge_e2n( ref_edge, edge, 0 ),"edge n0 empty");
      RUS(REF_EMPTY,ref_edge_e2n( ref_edge, edge, 1 ),"edge n1 empty");
      RSS( ref_adj_add( ref_edge_adj( ref_edge ), 
			ref_edge_e2n( ref_edge, edge, 0 ), 
			edge ), "adj n0");
      RSS( ref_adj_add( ref_edge_adj( ref_edge ), 
			ref_edge_e2n( ref_edge, edge, 1 ), 
			edge ), "adj n1");
    }

  ref_edge_node(ref_edge) = ref_grid_node(ref_grid);

  return REF_SUCCESS;
}

REF_STATUS ref_edge_free( REF_EDGE ref_edge )
{
  if ( NULL == (void *)ref_edge ) return REF_NULL;

  RSS( ref_adj_free( ref_edge_adj( ref_edge ) ), "free adj" );
  ref_free( ref_edge->e2n );

  ref_free( ref_edge );

  return REF_SUCCESS;
}

REF_STATUS ref_edge_with( REF_EDGE ref_edge, 
			  REF_INT node0, REF_INT node1,
			  REF_INT *edge )
{
  REF_INT item, ref;
  REF_INT n0, n1;

  *edge=REF_EMPTY;

  each_ref_adj_node_item_with_ref( ref_edge_adj(ref_edge), node0, item, ref)
    {
      n0 = ref_edge_e2n(ref_edge,ref,0);
      n1 = ref_edge_e2n(ref_edge,ref,1);
      if ( ( n0==node0 && n1==node1 ) ||
	   ( n0==node1 && n1==node0 ) )
	{
	  *edge=ref;
	  return REF_SUCCESS;
	}
    }

  return REF_NOT_FOUND;
}
