
#include <stdlib.h>
#include <stdio.h>

#include "ref_adj.h"

REF_STATUS ref_adj_create( REF_ADJ *ref_adj_ptr )
{
  REF_ADJ ref_adj;

  (*ref_adj_ptr) = NULL;
  (*ref_adj_ptr) = (REF_ADJ)malloc( sizeof(REF_ADJ_STRUCT) );
  RNS(*ref_adj_ptr,"malloc ref_adj NULL");
  ref_adj = (*ref_adj_ptr);

  ref_adj->blank = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_adj_free( REF_ADJ ref_adj )
{
  if ( NULL == (void *)ref_adj ) return REF_NULL;
  ref_cond_free( ref_adj );
  return REF_SUCCESS;
}

REF_STATUS ref_adj_inspect( REF_ADJ ref_adj )
{
  printf("ref_adj = %p\n",(void *)ref_adj);
  printf(" blank = %d\n",ref_adj->blank);

  return REF_SUCCESS;
}
