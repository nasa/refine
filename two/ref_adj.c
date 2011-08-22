
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
      ref_adj->item[i].reference = REF_EMPTY;
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
  printf("ref_adj = %p\n",(void *)ref_adj);
  printf(" blank = %d\n",ref_adj->blank);

  return REF_SUCCESS;
}
