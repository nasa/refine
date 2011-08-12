
#include <stdlib.h>
#include <stdio.h>

#include "ref_node.h"

REF_STATUS ref_node_create( REF_NODE *ref_node_ptr )
{
  REF_INT max, node;
  *ref_node_ptr = (REF_NODE)malloc( sizeof(REF_NODE_STRUCT) );
  RNS(*ref_node_ptr,"malloc ref_node NULL");

  max = 10;

  (*ref_node_ptr)->n = 0;
  (*ref_node_ptr)->max = max;
  (*ref_node_ptr)->global       = (REF_INT *)malloc(max*sizeof(REF_INT));
  RNS((*ref_node_ptr)->global,"malloc global NULL");

  for (node=0;node<ref_node_max(*ref_node_ptr);node++)
    (*ref_node_ptr)->global[node] = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  if ( NULL == (void *)ref_node ) return REF_NULL;
  ref_cond_free( ref_node->global );
  ref_cond_free( ref_node );
  return REF_SUCCESS;
}

REF_STATUS ref_node_inspect( REF_NODE ref_node )
{
  REF_INT node;
  printf("ref_node = %p\n",(void *)ref_node);
  printf(" n = %d\n",ref_node_n(ref_node));
  printf(" max = %d\n",ref_node_max(ref_node));
  for (node=0;node<ref_node_max(ref_node);node++)
    printf(" global[%d] = %d\n",node,ref_node_global(ref_node,node));
  return REF_SUCCESS;
}

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node )
{
  *node = ref_node_n(ref_node);
  ref_node->global[*node] = global;
  (ref_node->n)++;
  return REF_SUCCESS;
}

