
#include <stdlib.h>
#include <stdio.h>

#include "ref_adj.h"

REF_STATUS ref_adj_create( REF_ADJ *ref_adj_ptr )
{
  REF_ADJ ref_adj;
  REF_INT i;

  (*ref_adj_ptr) = NULL;
  (*ref_adj_ptr) = (REF_ADJ)malloc( sizeof(REF_ADJ_STRUCT) );
  RNS(*ref_adj_ptr,"malloc ref_adj NULL");
  ref_adj = (*ref_adj_ptr);

  ref_adj_nnode(ref_adj) = 10;
  ref_adj_nitem(ref_adj) = 20;

  ref_adj->first = (REF_INT *)malloc( ref_adj_nnode(ref_adj) * 
				      sizeof(REF_INT) );
  RNS(ref_adj->first,"malloc ref_adj->first NULL");
  for (i = 0; i < ref_adj_nnode(ref_adj) ; i++ )
    ref_adj->first[i] = REF_EMPTY;

  ref_adj->item = (REF_ADJ_ITEM)malloc( ref_adj_nitem(ref_adj) * 
					sizeof(REF_ADJ_ITEM_STRUCT) );
  RNS(ref_adj->item,"malloc ref_adj->item NULL");
  for (i = 0; i < ref_adj_nitem(ref_adj) ; i++ )
    {
      ref_adj->item[i].ref = REF_EMPTY;
      ref_adj->item[i].next = i+1;
    }
  ref_adj->item[ref_adj_nitem(ref_adj)-1].next = REF_EMPTY;
  ref_adj->blank = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_free( REF_ADJ ref_adj )
{
  if ( NULL == (void *)ref_adj ) return REF_NULL;
  ref_cond_free( ref_adj->first );
  ref_cond_free( ref_adj->item );
  ref_cond_free( ref_adj );
  return REF_SUCCESS;
}

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj )
{
  REF_INT node, item;
  printf("ref_adj = %p\n",(void *)ref_adj);
  printf(" blank = %d\n",ref_adj->blank);
  for ( node = 0 ; node < ref_adj_nnode( ref_adj ) ; node++ )
    printf(" first[%d] = %d\n",node,ref_adj->first[node]);
  for ( item = 0 ; item < ref_adj_nitem( ref_adj ) ; item++ )
    printf(" item[%d] = %d : %d\n",item,
	   ref_adj->item[item].next,ref_adj->item[item].ref);

  return REF_SUCCESS;
}

REF_STATUS ref_adj_add( REF_ADJ ref_adj, REF_INT node, REF_INT reference )
{
  REF_INT item;
  REF_INT orig, chunk, i;

  if ( node < 0 ) return REF_INVALID;

  if ( node >= ref_adj_nnode(ref_adj) )
    {
      orig = ref_adj_nnode( ref_adj );
      chunk = 100 + MAX( 0, node-orig );
      ref_adj_nnode(ref_adj) = orig + chunk;
      ref_adj->first = (REF_INT *)realloc( ref_adj->first,
					   ref_adj_nnode(ref_adj) * 
					   sizeof(REF_INT) );
      RNS(ref_adj->first,"realloc ref_adj->first NULL");
      for (i = orig; i < ref_adj_nnode(ref_adj) ; i++ )
	ref_adj->first[i] = REF_EMPTY;
    }

  if ( REF_EMPTY == ref_adj_blank(ref_adj) )
    {
      orig = ref_adj_nitem( ref_adj );
      chunk = 100;
      ref_adj_nitem( ref_adj ) =  orig + chunk;
      ref_adj->item = (REF_ADJ_ITEM)realloc( ref_adj->item,
					     ref_adj_nitem(ref_adj) * 
					     sizeof(REF_ADJ_ITEM_STRUCT) );
      RNS(ref_adj->item,"realloc ref_adj->item NULL");
      for (i = orig; i < ref_adj_nitem(ref_adj) ; i++ )
	{
	  ref_adj->item[i].ref = REF_EMPTY;
	  ref_adj->item[i].next = i+1;
	}
      ref_adj->item[ref_adj_nitem(ref_adj)-1].next = REF_EMPTY;
      ref_adj->blank = orig;
    }

  item = ref_adj_blank(ref_adj);
  ref_adj_blank(ref_adj) = ref_adj_item_next(ref_adj,item);

  ref_adj_item_ref(ref_adj,item) = reference;
  ref_adj_item_next(ref_adj,item) = ref_adj_first(ref_adj,node);

  ref_adj->first[node] = item;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_remove( REF_ADJ ref_adj, REF_INT node, REF_INT reference )
{
  REF_INT item, ref;
  REF_INT target, parent;

  item = ref_adj_first( ref_adj, node );

  if ( !ref_adj_valid( item ) ) return REF_INVALID;

  if ( reference == ref_adj_item_ref( ref_adj, item ) )
    {
      ref_adj->first[node] = ref_adj_item_next( ref_adj, item );
      ref_adj_item_next( ref_adj, item ) = ref_adj_blank( ref_adj );
      ref_adj_blank( ref_adj ) = item;
      ref_adj_item_ref( ref_adj, item ) = REF_EMPTY;
      return REF_SUCCESS;
    }
  
  target = REF_EMPTY;
  parent = REF_EMPTY;
  ref_adj_for( ref_adj, node, item, ref)
    {
      if ( ref == reference )
	{
	  target = item;
	  break;
	}
      else
	{
	  parent = item;
	}
    }

  if ( REF_EMPTY == target )
    return REF_INVALID;

  if ( REF_EMPTY == parent )
    RSS( REF_FAILURE, "parent empty");
  
  ref_adj_item_next( ref_adj, parent ) = ref_adj_item_next( ref_adj, item );
  ref_adj_item_next( ref_adj, item ) = ref_adj_blank( ref_adj );
  ref_adj_blank( ref_adj ) = item;
  ref_adj_item_ref( ref_adj, item ) = REF_EMPTY;

  return REF_SUCCESS;
}
