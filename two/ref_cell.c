
#include <stdlib.h>
#include <stdio.h>

#include "ref_cell.h"

REF_STATUS ref_cell_create( REF_INT node_per, REF_CELL *ref_cell_ptr )
{
  REF_INT i;
  REF_INT max;
  (*ref_cell_ptr) = NULL;
  (*ref_cell_ptr) = (REF_CELL)malloc( sizeof(REF_CELL_STRUCT) );
  RNS(*ref_cell_ptr,"malloc ref_cell NULL");

  max = 100;

  (*ref_cell_ptr)->node_per = node_per;
  (*ref_cell_ptr)->n = 0;
  (*ref_cell_ptr)->max = max;
  (*ref_cell_ptr)->c2n = (REF_INT *)malloc(node_per*max*sizeof(REF_INT));
  RNS((*ref_cell_ptr)->c2n,"malloc c2n NULL");
  for ( i=0 ; i < max ; i++ ) 
    {
      (*ref_cell_ptr)->c2n[0+node_per*i] = REF_EMPTY; 
      (*ref_cell_ptr)->c2n[1+node_per*i] = i+1; 
    }
  (*ref_cell_ptr)->c2n[1+node_per*(max-1)] = REF_EMPTY;
  (*ref_cell_ptr)->blank = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_free( REF_CELL ref_cell )
{
  if ( NULL == (void *)ref_cell ) return REF_NULL;
  ref_cond_free( ref_cell->c2n );
  ref_cond_free( ref_cell );
  return REF_SUCCESS;
}

REF_STATUS ref_cell_add( REF_CELL ref_cell, REF_INT *nodes, REF_INT *new_cell )
{
  int node, cell;
  int orig, chunk;

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
      for (cell=orig;cell < ref_cell_max(ref_cell); cell++ ) 
	{
	  ref_cell->c2n[0+ref_cell_node_per(ref_cell)*cell] = REF_EMPTY; 
	  ref_cell->c2n[1+ref_cell_node_per(ref_cell)*cell] = cell+1; 
	}
      ref_cell->c2n[1+4*(ref_cell_max(ref_cell)-1)] = REF_EMPTY; 
      ref_cell->blank = orig;
    }

  cell = ref_cell_blank(ref_cell);
  ref_cell_blank(ref_cell) = ref_cell_c2n(ref_cell,1,cell);
  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    ref_cell_c2n(ref_cell,node,cell) = nodes[node];

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
