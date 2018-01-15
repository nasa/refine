
/* Copyright 2014 United States Government as represented by the
 * Administrator of the National Aeronautics and Space
 * Administration. No copyright is claimed in the United States under
 * Title 17, U.S. Code.  All Other Rights Reserved.
 *
 * The refine platform is licensed under the Apache License, Version
 * 2.0 (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

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
  ref_free( ref_list->value );
  ref_free( ref_list );
  return REF_SUCCESS;
}

REF_STATUS ref_list_deep_copy( REF_LIST *ref_list_ptr, REF_LIST original )
{
  REF_LIST ref_list;
  REF_INT i;

  ref_malloc( *ref_list_ptr, 1, REF_LIST_STRUCT );
  ref_list = (*ref_list_ptr);

  ref_list_n(ref_list)   = ref_list_n(original);
  ref_list_max(ref_list) = ref_list_max(original);

  ref_malloc( ref_list->value, ref_list_max(ref_list), REF_INT);
  for(i=0;i< ref_list_n(ref_list);i++)
    ref_list->value[i] = original->value[i];
    
  return REF_SUCCESS;
}

REF_STATUS ref_list_inspect( REF_LIST ref_list )
{
  REF_INT i;

  if ( NULL == (void *)ref_list ) return REF_NULL;

  printf("list has %d items\n",ref_list_n(ref_list));
  for(i=0;i< ref_list_n(ref_list);i++)
    printf( "value[%d] = %d\n", i, ref_list->value[i] );

  return REF_SUCCESS;
}

REF_STATUS ref_list_push( REF_LIST ref_list, REF_INT last )
{

  if ( ref_list_max(ref_list) == ref_list_n(ref_list) )
    {
      ref_list_max(ref_list) += 1000;
      ref_realloc( ref_list->value, ref_list_max(ref_list), REF_INT );
    }

  ref_list->value[ref_list_n( ref_list )] = last;

  ref_list_n( ref_list )++;

  return REF_SUCCESS;
}

REF_STATUS ref_list_pop( REF_LIST ref_list, REF_INT *last )
{

  if ( 0 == ref_list_n(ref_list) )
    {
      *last = REF_EMPTY;
      return REF_FAILURE;
    }

  ref_list_n( ref_list )--;
  *last = ref_list->value[ref_list_n( ref_list )];

  return REF_SUCCESS;
}

REF_STATUS ref_list_shift( REF_LIST ref_list, REF_INT *first )
{
  REF_INT i;
  if ( 0 == ref_list_n(ref_list) )
    {
      *first = REF_EMPTY;
      return REF_FAILURE;
    }

  *first = ref_list->value[0];

  ref_list_n( ref_list )--;
  for ( i=0;i<ref_list_n( ref_list );i++)
    ref_list->value[i] = ref_list->value[i+1];

  return REF_SUCCESS;
}

REF_STATUS ref_list_delete( REF_LIST ref_list, REF_INT item )
{
  REF_INT to, from;

  to = 0;
  for(from=0;from< ref_list_n(ref_list);from++)
    {
      if ( item != ref_list->value[from] ) 
	{
	  ref_list->value[to] = ref_list->value[from];
	  to++;
	}
    }

  if ( to == ref_list_n(ref_list) ) return REF_NOT_FOUND;

  ref_list_n(ref_list) = to;

  return REF_SUCCESS;
}

REF_STATUS ref_list_apply_offset( REF_LIST ref_list,
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

  RSS( ref_sort_heap_int( ref_list_n(ref_list), ref_list->value, order ), 
       "heap" );

  for ( i=0; i < ref_list_n(ref_list); i++ )
    order[i] = ref_list->value[order[i]];

  for ( i=0; i < ref_list_n(ref_list); i++ )
    ref_list->value[i] = order[i];

  ref_free( order );

  return REF_SUCCESS;
}

REF_STATUS ref_list_erase( REF_LIST ref_list )
{
  ref_list_n( ref_list ) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_list_allgather( REF_LIST ref_list, REF_MPI ref_mpi )
{
  REF_INT i;
  REF_INT *local_copy;
  REF_INT proc;
  REF_INT *counts;
  REF_INT total_count;

  ref_malloc( counts, ref_mpi_m(ref_mpi), REF_INT );

  RSS( ref_mpi_allgather( ref_mpi,
			  &(ref_list_n(ref_list)), counts, REF_INT_TYPE ), 
       "gather size");

  total_count = 0;
  each_ref_mpi_part( ref_mpi, proc )
    total_count += counts[proc];

  ref_malloc( local_copy, ref_list_n( ref_list ), REF_INT );
  for ( i=0; i < ref_list_n(ref_list); i++ )
    local_copy[i] = ref_list->value[i];

  if ( total_count > ref_list_max(ref_list) )
    { 
      ref_list_max(ref_list) = total_count;
      ref_free( ref_list->value );
      ref_malloc( ref_list->value, ref_list_max(ref_list), REF_INT );
    }

  RSS( ref_mpi_allgatherv( ref_mpi,
			   local_copy, counts, ref_list->value, REF_INT_TYPE ), 
       "gather values");

  ref_list_n( ref_list ) = total_count;

  ref_free( local_copy );
  ref_free( counts );

  return REF_SUCCESS;
}

REF_STATUS ref_list_contains( REF_LIST ref_list, REF_INT item, 
			      REF_BOOL *contains )
{
  REF_INT i;

  *contains = REF_FALSE;

  for(i=0;i< ref_list_n(ref_list);i++)
    if (ref_list->value[i] == item)
      {
	*contains = REF_TRUE;
	return REF_SUCCESS;
      }

  return REF_SUCCESS;
}
