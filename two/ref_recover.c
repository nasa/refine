
#include <stdlib.h>
#include <stdio.h>

#include "ref_recover.h"
#include "ref_malloc.h"

REF_STATUS ref_recover_create( REF_RECOVER *ref_recover_ptr, REF_GRID ref_grid )
{
  REF_RECOVER ref_recover;

  ref_malloc( *ref_recover_ptr, 1, REF_RECOVER_STRUCT );
  ref_recover = (*ref_recover_ptr);

  ref_recover_grid(ref_recover) = ref_grid;
  ref_recover_n(ref_recover) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_recover_free( REF_RECOVER ref_recover )
{
  if ( NULL == (void *)ref_recover ) return REF_NULL;
  ref_free( ref_recover );
  return REF_SUCCESS;
}

REF_STATUS ref_recover_insert_twod( REF_RECOVER ref_recover, REF_DBL *xz,
				    REF_INT *node_ptr )
{
  REF_GRID ref_grid = ref_recover_grid(ref_recover);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT global;
  REF_INT node0, node1;

  *node_ptr = REF_EMPTY;

  RSS( ref_node_next_global( ref_node, &global ), "next global");
  RSS( ref_node_add( ref_node, global, &node0 ), "add node");
  ref_node_xyz(ref_node,0,node0) = xz[0];
  ref_node_xyz(ref_node,1,node0) = 0.0;
  ref_node_xyz(ref_node,2,node0) = xz[1];

  RSS( ref_node_next_global( ref_node, &global ), "next global");
  RSS( ref_node_add( ref_node, global, &node1 ), "add node");
  ref_node_xyz(ref_node,0,node1) = xz[0];
  ref_node_xyz(ref_node,1,node1) = 1.0;
  ref_node_xyz(ref_node,2,node1) = xz[1];

  return REF_SUCCESS;
}
