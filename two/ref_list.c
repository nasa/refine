
#include <stdlib.h>
#include <stdio.h>

#include "ref_list.h"

REF_STATUS ref_list_create( REF_LIST *ref_list_ptr )
{
  REF_LIST ref_list;

  (*ref_list_ptr) = NULL;
  (*ref_list_ptr) = (REF_LIST)malloc( sizeof(REF_LIST_STRUCT) );
  RNS(*ref_list_ptr,"malloc ref_list NULL");
  ref_list = (*ref_list_ptr);

  ref_list_n(ref_list) = 0;
  ref_list_max(ref_list) = 10;

  ref_list->value = (REF_INT *)malloc( ref_list_max(ref_list) * 
				       sizeof(REF_INT) );
  RNS(ref_list->value,"malloc ref_list->value NULL");

  return REF_SUCCESS;
}

REF_STATUS ref_list_free( REF_LIST ref_list )
{
  if ( NULL == (void *)ref_list ) return REF_NULL;
  ref_cond_free( ref_list->value );
  ref_cond_free( ref_list );
  return REF_SUCCESS;
}

REF_STATUS ref_list_add( REF_LIST ref_list, REF_INT value )
{

  if ( ref_list_max(ref_list) == ref_list_n(ref_list) )
    {
      ref_list_max(ref_list) += 1000;
      ref_list->value = (REF_INT *)realloc( (void *)(ref_list->value),
					    ref_list_max(ref_list) * 
					    sizeof(REF_INT) );
      RNS(ref_list->value,"realloc ref_list->value NULL");
    }

  ref_list->value[ref_list_n( ref_list )] = value;

  ref_list_n( ref_list )++;

  return REF_SUCCESS;
}

REF_STATUS ref_list_remove( REF_LIST ref_list, REF_INT *value )
{

  if ( 0 == ref_list_n(ref_list) )
    {
      *value = REF_EMPTY;
      return REF_FAILURE;
    }

  ref_list_n( ref_list )--;
  *value = ref_list->value[ref_list_n( ref_list )];

  return REF_SUCCESS;
}

