
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

#include "ref_elast.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_edge.h"

#include "ref_cell.h"

#include "ref_malloc.h"

REF_STATUS ref_elast_create( REF_ELAST *ref_elast_ptr, REF_GRID ref_grid )
{
  REF_ELAST ref_elast;
  
  ref_malloc( *ref_elast_ptr, 1, REF_ELAST_STRUCT );
  
  ref_elast = *ref_elast_ptr;

  ref_elast_grid(ref_elast) = ref_grid;
  RSS( ref_comprow_create( &(ref_elast->ref_comprow),
                           ref_elast_grid(ref_elast) ), "comprow" );

  ref_malloc_init( ref_elast->a,
                   3*3*ref_comprow_nnz(ref_elast_comprow(ref_elast)),
                   REF_DBL, 0.0 );
  
  return REF_SUCCESS;
}

REF_STATUS ref_elast_free( REF_ELAST ref_elast )
{
  if ( NULL == (void *)ref_elast ) return REF_NULL;

  ref_free( ref_elast->a );
  ref_comprow_free( ref_elast->ref_comprow );
  
  ref_free( ref_elast );
  
  return REF_SUCCESS;
}

REF_STATUS ref_elast_assemble( REF_ELAST ref_elast )
{
  REF_COMPROW ref_comprow = ref_elast_comprow(ref_elast);
  REF_GRID ref_grid = ref_elast_grid(ref_elast);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT i;
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];

  for(i=0;i<9*ref_comprow_nnz(ref_comprow);i++)
    ref_elast->a[i]=0.0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
    }
  
  return REF_SUCCESS;
}
