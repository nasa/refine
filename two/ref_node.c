
#include <stdlib.h>
#include <stdio.h>

#include "ref_node.h"

/* REF_EMPTY is terminatior, next avalable is shifted by 2*/
#define next2index(next) (-(next)-2)
#define index2next(index) (-2-(index))

REF_STATUS ref_node_create( REF_NODE *ref_node_ptr )
{
  REF_INT max, node;
  REF_NODE ref_node;

  (*ref_node_ptr) = NULL;
  (*ref_node_ptr) = (REF_NODE)malloc( sizeof(REF_NODE_STRUCT) );
  RNS(*ref_node_ptr,"malloc ref_node NULL");

  ref_node = *ref_node_ptr;

  max = 10;

  ref_node_n(ref_node) = 0;
  ref_node_max(ref_node) = max;

  ref_node_partition(ref_node) = REF_EMPTY;

  ref_node_n_global(ref_node) = 0;

  ref_node->global = (REF_INT *)malloc(max*sizeof(REF_INT));
  RNS(ref_node->global,"malloc global NULL");

  for (node=0;node<ref_node_max(ref_node);node++)
    ref_node->global[node] = index2next(node+1);
  ref_node->global[(ref_node->max)-1] = REF_EMPTY;
  ref_node->blank = index2next(0);

  ref_node->part = (REF_INT *)malloc(max*sizeof(REF_INT));
  RNS(ref_node->part,"malloc part NULL");

  ref_node->xyz = (REF_DBL *)malloc(max*3*sizeof(REF_DBL));
  RNS(ref_node->xyz,"malloc xyz NULL");


  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  if ( NULL == (void *)ref_node ) return REF_NULL;
  ref_cond_free( ref_node->xyz );
  ref_cond_free( ref_node->part );
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
  printf(" blank = %d\n",ref_node->blank);
  for (node=0;node<ref_node_max(ref_node);node++)
    printf(" raw global[%d] = %d\n",node,ref_node->global[node]);
  return REF_SUCCESS;
}

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node )
{
  REF_INT extra;
  REF_INT orig, chunk;

  if ( global < 0 ) RSS( REF_INVALID, "invalid global node");

  if ( REF_EMPTY == ref_node->blank )
    {
      orig = ref_node_max(ref_node);
      chunk = 5000;
      ref_node->max = orig + chunk;
      ref_node->global = (REF_INT *)realloc( ref_node->global,
					     ref_node_max(ref_node) *
					     sizeof(REF_INT) );
      RNS(ref_node->global,"remalloc global NULL");
      for (extra=orig;extra < ref_node_max(ref_node); extra++ ) 
	  ref_node->global[extra] = index2next(extra+1);
      ref_node->global[ref_node_max(ref_node)-1] = REF_EMPTY; 
      ref_node->blank = index2next(orig);

      ref_node->part = (REF_INT *)realloc( ref_node->part,
					   ref_node_max(ref_node) *
					   sizeof(REF_INT) );
      RNS(ref_node->part,"remalloc part NULL");

      ref_node->xyz = (REF_DBL *)realloc( ref_node->xyz,
					  ref_node_max(ref_node) *
					  3 * sizeof(REF_DBL) );
      RNS(ref_node->xyz,"remalloc xyz NULL");
    }

  *node = next2index(ref_node->blank);
  ref_node->blank = ref_node->global[*node];

  ref_node->global[*node] = global;
  (ref_node->n)++;
  return REF_SUCCESS;
}

REF_STATUS ref_node_remove( REF_NODE ref_node, REF_INT node )
{
  if ( ! ref_node_valid(ref_node,node) ) return REF_INVALID;
 
  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_next_global( REF_NODE ref_node, REF_INT *global )
{
  (*global) = ref_node_n_global(ref_node);
  ref_node_n_global(ref_node)++;

  return REF_SUCCESS;
}

REF_STATUS ref_node_local( REF_NODE ref_node, REF_INT global, REF_INT *local )
{
  REF_INT node;

  (*local) = REF_EMPTY;
  if ( global < 0 ) return REF_INVALID;

  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( ref_node->global[node] == global ) (*local) = node;

  if ( (*local) == REF_EMPTY ) return REF_FAILURE;
  return REF_SUCCESS;
}

REF_STATUS ref_node_compact( REF_NODE ref_node, REF_INT **o2n_ptr )
{
  REF_INT node;
  REF_INT nnode;
  REF_INT *o2n;
  
  *o2n_ptr = (REF_INT *)malloc( ref_node_max(ref_node) * sizeof(REF_INT) );
  RNS(*o2n_ptr,"malloc o2n NULL");
  o2n = *o2n_ptr;

  nnode = 0;
  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( ref_node_valid(ref_node,node) )
      {
	o2n[node] = nnode;
	nnode++;
      }
    else
      {
	o2n[node] = REF_EMPTY;
      }
  RES( nnode, ref_node_n(ref_node), "nnode miscount" );

  return REF_SUCCESS;
}

