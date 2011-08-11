
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
  (*ref_node_ptr)->global       = (REF_INT *)malloc(max*sizeof(REF_INT));
  RNS((*ref_node_ptr)->global,"malloc global NULL");

  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  free( ref_node );
  return REF_SUCCESS;
}

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node )
{
  *node = 0;
  ref_node->global[*node] = global;
  return REF_SUCCESS;
}
