
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

#include "ref_comprow.h"
#include "ref_grid.h"
#include "ref_node.h"

#include "ref_malloc.h"

REF_STATUS ref_comprow_create( REF_COMPROW *ref_comprow_ptr, REF_GRID ref_grid )
{
  REF_COMPROW ref_comprow;
  REF_NODE ref_node = ref_grid_node(ref_grid);
  ref_malloc( *ref_comprow_ptr, 1, REF_COMPROW_STRUCT );

  ref_comprow = *ref_comprow_ptr;

  ref_comprow_nnz(ref_comprow) = 0;
  ref_malloc( ref_comprow->first, 1+ref_node_max(ref_node), REF_INT );
  ref_comprow->col = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_comprow_free( REF_COMPROW ref_comprow )
{
  if ( NULL == (void *)ref_comprow ) return REF_NULL;

  ref_free( ref_comprow->col );
  ref_free( ref_comprow->first );

  ref_free( ref_comprow );
  
  return REF_SUCCESS;
}
