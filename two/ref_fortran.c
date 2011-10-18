
#include <stdlib.h>
#include <stdio.h>

#include "ref_fortran.h"
#include "ref_grid.h"

static REF_GRID ref_grid = NULL;

REF_STATUS ref_init_node_(REF_INT *nnodes, REF_INT *nnodesg,
			  REF_INT *l2g, REF_INT *part, REF_INT *partition, 
			  REF_DBL *x, REF_DBL *y, REF_DBL *z )
{
  REF_NODE ref_node;
  int node, pos;
  RSS( ref_grid_create( &ref_grid ), "create grid" );
  ref_node = ref_grid_node(ref_grid);

  ref_node_total(ref_node) = *nnodesg;
  ref_node_partition(ref_node) = *partition;

  for (node=0;node<(*nnodes);node++)
    {
      RSS( ref_node_add( ref_node, l2g[node]-1, &pos ), 
	   "add node");
      ref_node_xyz(ref_node,0,pos) = x[node];
      ref_node_xyz(ref_node,1,pos) = y[node];
      ref_node_xyz(ref_node,2,pos) = z[node];
      ref_node_part(ref_node,pos) = part[node]-1;
    }
  return REF_SUCCESS;
}

REF_STATUS ref_free_( void )
{
  RSS( ref_grid_free( ref_grid ), "free grid");
  ref_grid = NULL;
  return REF_SUCCESS;
}

