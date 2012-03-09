
#include <stdlib.h>
#include <stdio.h>

#include "ref_edge.h"

#include "ref_malloc.h"
#include "ref_mpi.h"

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

static REF_STATUS ref_edge_part( REF_EDGE ref_edge, REF_INT edge, 
				 REF_INT *part )
{
  REF_NODE ref_node = ref_edge_node(ref_edge);

  if ( ref_node_global(ref_node,ref_edge_e2n( ref_edge, edge, 0 ) ) <
       ref_node_global(ref_node,ref_edge_e2n( ref_edge, edge, 1 ) ) )
    {
      *part = ref_node_part(ref_node,ref_edge_e2n( ref_edge, edge, 0 ));
    }
  else
    {
      *part = ref_node_part(ref_node,ref_edge_e2n( ref_edge, edge, 1 ));
    }

  return REF_SUCCESS;
}

REF_STATUS ref_edge_ghost_int( REF_EDGE ref_edge, REF_INT *data )
{
  REF_NODE ref_node = ref_edge_node(ref_edge);
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT edge;
  REF_INT part;
  REF_INT *a_next, *a_edge;
  REF_INT *a_nodes, *b_nodes;
  REF_INT *a_data, *b_data;

  REF_INT node0, node1;
  REF_INT request;

  ref_malloc_init( a_size, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_size, ref_mpi_n, REF_INT, 0 );

  for ( edge=0 ; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      RSS( ref_edge_part( ref_edge, edge, &part ), "edge part" );
      if ( part != ref_mpi_id ) a_size[part]++;
    }

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  ref_malloc( a_nodes, 2*a_total, REF_INT );
  ref_malloc( a_data, a_total, REF_INT );
  ref_malloc( a_edge, a_total, REF_INT );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_nodes, 2*b_total, REF_INT );
  ref_malloc( b_data, b_total, REF_INT );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  for ( edge=0 ; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      RSS( ref_edge_part( ref_edge, edge, &part ), "edge part" );
      if ( part != ref_mpi_id )
	{
	  a_edge[a_next[part]] = edge;
	  a_nodes[0+2*a_next[part]] = 
	    ref_node_global(ref_node,ref_edge_e2n( ref_edge, edge, 0 ) );
	  a_nodes[1+2*a_next[part]] = 
	    ref_node_global(ref_node,ref_edge_e2n( ref_edge, edge, 1 ) );
	  a_next[part]++;
	}
    }

  RSS( ref_mpi_alltoallv( a_nodes, a_size, b_nodes, b_size, 
			  2, REF_INT_TYPE ), 
       "alltoallv requested nodes");

  for ( request=0 ; request < b_total ; request++ )
    {
      RSS( ref_node_local( ref_node, b_nodes[0+2*request], &node0 ), "loc 0");
      RSS( ref_node_local( ref_node, b_nodes[1+2*request], &node1 ), "loc 1");
      RSS( ref_edge_with( ref_edge, node0, node1, &edge ), "find edge" );
      b_data[request] = data[edge];
    }

  RSS( ref_mpi_alltoallv( b_data, b_size, a_data, a_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv return data");

  for ( request=0 ; request < a_total ; request++ )
    {
      data[a_edge[request]] = a_data[request];
    }

  ref_free( a_next );

  ref_free( b_data );
  ref_free( b_nodes );

  ref_free( a_edge );

  ref_free( a_data );
  ref_free( a_nodes );

  ref_free( b_size );
  ref_free( a_size );

  return REF_SUCCESS;
}
