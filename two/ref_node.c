
#include <stdlib.h>
#include <stdio.h>

#include "ref_node.h"

REF_STATUS ref_node_create( REF_NODE *ref_node_ptr )
{
  REF_INT max;
  *ref_node_ptr = (REF_NODE)malloc( sizeof(REF_NODE_STRUCT) );
  RNS(*ref_node_ptr,"malloc ref_node NULL");

  max = 100;

  (*ref_node_ptr)->n = 0;
  (*ref_node_ptr)->max = max;

  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  free( ref_node );
  return REF_SUCCESS;
}
