
#include <stdlib.h>
#include <stdio.h>

#include "ref_list.h"
#include "ref_malloc.h"
#include "ref_sort.h"

REF_STATUS ref_list_create( REF_LIST *ref_list_ptr )
{
  REF_LIST ref_list;

  ref_malloc( *ref_list_ptr, 1, REF_LIST_STRUCT );

  ref_list = (*ref_list_ptr);

  ref_list_n(ref_list) = 0;
  ref_list_max(ref_list) = 10;

  ref_malloc( ref_list->value, ref_list_max(ref_list), REF_INT);

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
      ref_realloc( ref_list->value, ref_list_max(ref_list), REF_INT );
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

REF_STATUS ref_list_shift( REF_LIST ref_list, 
			   REF_INT equal_and_above, REF_INT offset )
{
  REF_INT i;

  for(i=0;i< ref_list_n(ref_list);i++)
    if ( ref_list->value[i] >= equal_and_above )
      ref_list->value[i] += offset;

  return REF_SUCCESS;
}

REF_STATUS ref_list_sort( REF_LIST ref_list )
{
  REF_INT *order;
  REF_INT i;

  /* see if it is too short to require sorting */
  if ( 2 > ref_list_n( ref_list ) ) return REF_SUCCESS;

  ref_malloc( order, ref_list_n( ref_list ), REF_INT );

  RSS( ref_sort_heap( ref_list_n(ref_list), ref_list->value, order ), "heap" );

  for ( i=0; i<  ref_list_n(ref_list); i++ )
    order[i] = ref_list->value[order[i]];

  for ( i=0; i<  ref_list_n(ref_list); i++ )
    ref_list->value[i] = order[i];

  ref_free( order );

  return REF_SUCCESS;
}
