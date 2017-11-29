
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
#include <time.h>

#include "ref_interp.h"

#include "ref_malloc.h"

REF_STATUS ref_interp_create( REF_INTERP *ref_interp_ptr )
{
  REF_INTERP ref_interp;
  
  ref_malloc( *ref_interp_ptr, 1, REF_INTERP_STRUCT );
  ref_interp = ( *ref_interp_ptr );

  ref_interp->steps = 0;
  ref_interp->node = NULL;
  ref_interp->bary = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_interp_free( REF_INTERP ref_interp )
{
  if ( NULL == (void *)ref_interp )
    return REF_NULL;
  ref_free( ref_interp->bary );
  ref_free( ref_interp->node );
  ref_free( ref_interp );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_exhaustive_enclosing_tet( REF_GRID ref_grid, REF_DBL *xyz,
						REF_INT *node, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT guess, best_guess;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
 
  best_guess = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_valid_cell( ref_cell, guess)
    {
      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RXS( ref_node_bary4( ref_node, nodes, xyz, current_bary ), 
	   REF_DIV_ZERO, "bary");
      min_bary = MIN( MIN(current_bary[0],current_bary[1]),
		      MIN(current_bary[2],current_bary[3]));
      if ( REF_EMPTY == best_guess || min_bary > best_bary )
	{
	  best_guess = guess;
	  best_bary = min_bary;
	}
    }
  
  RUS( REF_EMPTY, best_guess, "failed to find cell");

  RSS( ref_cell_nodes( ref_cell, best_guess, node), "cell" );
  RSS( ref_node_bary4( ref_node, node, xyz, bary ), "bary");
  
  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate( REF_INTERP ref_interp, 
			      REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT node;

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  ref_malloc( ref_interp->node, 
	      4*ref_node_max(to_node), 
	      REF_INT );
  ref_malloc( ref_interp->bary, 
	      4*ref_node_max(to_node), 
	      REF_DBL );

  each_ref_node_valid_node( to_node, node )
    {
      RSS(ref_interp_exhaustive_enclosing_tet( from_grid,
					       ref_node_xyz_ptr(to_node,node),
					       &(ref_interp->node[4*node]), 
					       &(ref_interp->bary[4*node]) ), 
	  "exhast");
    }

  return REF_SUCCESS;
}
