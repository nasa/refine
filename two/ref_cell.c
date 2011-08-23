
#include <stdlib.h>
#include <stdio.h>

#include "ref_cell.h"

REF_STATUS ref_cell_create( REF_INT node_per, REF_CELL *ref_cell_ptr )
{
  REF_INT cell;
  REF_INT max;
  REF_CELL ref_cell;

  (*ref_cell_ptr) = NULL;

  if ( node_per > REF_CELL_MAX_NODE_PER)
    {
      RSS( REF_FAILURE, "node_per limited to REF_CELL_MAX_NODE_PER");
    }

  (*ref_cell_ptr) = (REF_CELL)malloc( sizeof(REF_CELL_STRUCT) );
  RNS(*ref_cell_ptr,"malloc ref_cell NULL");

  ref_cell = (*ref_cell_ptr);

  max = 100;

  ref_cell_node_per(ref_cell) = node_per;
  switch ( node_per )
    {
    case 4:
      ref_cell_edge_per(ref_cell) = 6;
      break;
    case 5:
      ref_cell_edge_per(ref_cell) = 8;
      break;
    case 6:
      ref_cell_edge_per(ref_cell) = 9;
      break;
    case 8:
      ref_cell_edge_per(ref_cell) = 12;
      break;
    default:
      ref_cell_edge_per(ref_cell) = 0;
      break;    
    }

  ref_cell_n(ref_cell) = 0;
  ref_cell_max(ref_cell) = max;

  ref_cell->c2n = (REF_INT *)malloc(ref_cell_max(ref_cell) *
				    ref_cell_node_per(ref_cell) *
				    sizeof(REF_INT));
  RNS(ref_cell->c2n,"malloc c2n NULL");
  ref_cell->c2e = (REF_INT *)malloc(ref_cell_max(ref_cell) *
				    ref_cell_edge_per(ref_cell) *
				    sizeof(REF_INT));
  RNS(ref_cell->c2e,"malloc c2e NULL");
  for ( cell = 0 ; cell < max ; cell++ ) 
    {
      ref_cell_c2n(ref_cell,0,cell) = REF_EMPTY;
      ref_cell_c2n(ref_cell,1,cell) = cell+1;
    }
  ref_cell_c2n(ref_cell,1,max-1) = REF_EMPTY;
  ref_cell_blank(ref_cell) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_free( REF_CELL ref_cell )
{
  if ( NULL == (void *)ref_cell ) return REF_NULL;
  ref_cond_free( ref_cell->c2n );
  ref_cond_free( ref_cell->c2e );
  ref_cond_free( ref_cell );
  return REF_SUCCESS;
}

REF_STATUS ref_cell_inspect( REF_CELL ref_cell )
{
  REF_INT cell, node;
  printf("ref_cell = %p\n",(void *)ref_cell);
  printf(" node_per = %d\n",ref_cell_node_per(ref_cell));
  printf(" n = %d\n",ref_cell_n(ref_cell));
  printf(" max = %d\n",ref_cell_max(ref_cell));
  printf(" blank = %d\n",ref_cell->blank);
  for (cell=0;cell<ref_cell_max(ref_cell);cell++)
    {
      printf(" %d:",cell);
      if ( ref_cell_valid(ref_cell,cell) )
	{
	  for (node=0;node<ref_cell_node_per(ref_cell);node++)
	    printf(" %d",ref_cell_c2n(ref_cell,node,cell));
	}
      else
	{
	  for (node=0;node<2;node++)
	    printf(" %d",ref_cell_c2n(ref_cell,node,cell));

	}
      printf("\n");

    }
  return REF_SUCCESS;
}

REF_STATUS ref_cell_add( REF_CELL ref_cell, REF_INT *nodes, REF_INT *new_cell )
{
  REF_INT node, cell;
  REF_INT orig, chunk;

  (*new_cell) = REF_EMPTY;

  if ( REF_EMPTY == ref_cell_blank(ref_cell) ) 
    {
      orig = ref_cell_max(ref_cell);
      chunk = 5000;
      ref_cell->max = orig + chunk;
      ref_cell->c2n = (REF_INT *)realloc( ref_cell->c2n,
					  ref_cell_node_per(ref_cell) *
					  ref_cell_max(ref_cell) *
					  sizeof(REF_INT) );
      RNS(ref_cell->c2n,"remalloc c2n NULL");
      ref_cell->c2e = (REF_INT *)realloc( ref_cell->c2e,
					  ref_cell_edge_per(ref_cell) *
					  ref_cell_max(ref_cell) *
					  sizeof(REF_INT) );
      RNS(ref_cell->c2e,"remalloc c2e NULL");
      for (cell=orig;cell < ref_cell_max(ref_cell); cell++ ) 
	{
	  ref_cell_c2n(ref_cell,0,cell)= REF_EMPTY; 
	  ref_cell_c2n(ref_cell,1,cell) = cell+1; 
	}
      ref_cell_c2n(ref_cell,1,(ref_cell->max)-1) = REF_EMPTY; 
      ref_cell_blank(ref_cell) = orig;
    }

  cell = ref_cell_blank(ref_cell);
  ref_cell_blank(ref_cell) = ref_cell_c2n(ref_cell,1,cell);
  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    ref_cell_c2n(ref_cell,node,cell) = nodes[node];

  ref_cell_n(ref_cell)++;

  (*new_cell) = cell;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_nodes( REF_CELL ref_cell, REF_INT cell, REF_INT *nodes )
{
  REF_INT node;
  if ( cell < 0 || cell > ref_cell_max(ref_cell) ) return REF_INVALID;
  if ( REF_EMPTY == ref_cell_c2n(ref_cell,0,cell) ) 
    return REF_INVALID;
  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    nodes[node] = ref_cell_c2n(ref_cell,node,cell);
  return REF_SUCCESS;
}

REF_STATUS ref_cell_empty_edges( REF_CELL ref_cell)
{
  REF_INT cell, edge;
  for ( cell = 0 ; cell < ref_cell_max(ref_cell) ; cell++ )
    for ( edge = 0 ; edge < ref_cell_edge_per(ref_cell) ; edge++ )
      ref_cell_c2e(ref_cell,edge,cell) = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_set_edge( REF_CELL ref_cell, 
			      REF_INT n0, REF_INT n1, REF_INT edge)
{
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT e2n0[6] = { 0, 0, 0, 1, 1, 2 };
  REF_INT e2n1[6] = { 1, 2, 3, 2, 3, 3 };
  REF_INT cell, cell_edge;
  REF_INT e0, e1;

  ref_cell_for_with_nodes( ref_cell, cell, nodes)
    for (cell_edge = 0; cell_edge < ref_cell_edge_per(ref_cell); cell_edge++)
      {
	e0 = nodes[e2n0[cell_edge]];
	e1 = nodes[e2n1[cell_edge]];
	if ( MAX(e0,e1) == MAX(n0,n1) && MIN(e0,e1) == MIN(n0,n1) )
	  ref_cell_c2e(ref_cell,cell_edge,cell) = edge;
      }

  return REF_SUCCESS;
}
