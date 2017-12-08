
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
#include <math.h>

#include "ref_search.h"

#include "ref_malloc.h"

#define MAX_NODE_LIST ( 100 )

REF_STATUS ref_search_create( REF_SEARCH *ref_search_ptr, REF_INT n )
{
  REF_SEARCH ref_search;
  
  ref_malloc( *ref_search_ptr, 1, REF_SEARCH_STRUCT );
  ref_search = ( *ref_search_ptr );

  ref_search->d = 3;
  ref_search->n = n;
  ref_search->empty = 0;

  ref_malloc_init( ref_search->item, ref_search->n, REF_INT, REF_EMPTY );
  ref_malloc_init( ref_search->left, ref_search->n, REF_INT, REF_EMPTY );
  ref_malloc_init( ref_search->right, ref_search->n, REF_INT, REF_EMPTY );

  ref_malloc( ref_search->pos, ref_search->d*ref_search->n, REF_DBL );
  ref_malloc( ref_search->radius, ref_search->n, REF_DBL );
  ref_malloc( ref_search->left_radius, ref_search->n, REF_DBL );
  ref_malloc( ref_search->right_radius, ref_search->n, REF_DBL );

  return REF_SUCCESS;
}

REF_STATUS ref_search_free( REF_SEARCH ref_search )
{
  if ( NULL == (void *)ref_search )
    return REF_NULL;
  ref_free( ref_search->right_radius );
  ref_free( ref_search->left_radius );
  ref_free( ref_search->radius );
  ref_free( ref_search->pos );
  ref_free( ref_search->right );
  ref_free( ref_search->left );
  ref_free( ref_search->item );
  ref_free( ref_search );
  return REF_SUCCESS;
}

REF_STATUS ref_search_insert( REF_SEARCH ref_search,
			      REF_INT item, REF_DBL *position, REF_DBL radius )
{
  REF_INT i, location;
  if ( ref_search->empty >= ref_search->n )
    RSS( REF_INCREASE_LIMIT, "need larger tree for more items" );

  location = ref_search->empty;
  (ref_search->empty)++;

  ref_search->item[location] = item;
  for (i=0;i<ref_search->d;i++)
    ref_search->pos[i+ref_search->d*location] = position[i];
  ref_search->radius[location] = radius;  

  return REF_SUCCESS;
}
