
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

#include "ref_interp.h"

#include "ref_malloc.h"

REF_STATUS ref_interp_create( REF_INTERP *ref_interp_ptr )
{
  REF_INTERP ref_interp;
  
  ref_malloc( *ref_interp_ptr, 1, REF_INTERP_STRUCT );
  ref_interp = ( *ref_interp_ptr );

  ref_interp->nexhaustive = 0;
  ref_interp->cell = NULL;
  ref_interp->bary = NULL;

  return REF_SUCCESS;
}

REF_STATUS ref_interp_free( REF_INTERP ref_interp )
{
  if ( NULL == (void *)ref_interp )
    return REF_NULL;
  ref_free( ref_interp->bary );
  ref_free( ref_interp->cell );
  ref_free( ref_interp );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_exhaustive_enclosing_tet( REF_GRID ref_grid, REF_DBL *xyz,
						REF_INT *cell, REF_DBL *bary )
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

  *cell = best_guess;
  RSS( ref_cell_nodes( ref_cell, best_guess, nodes), "cell" );
  RSS( ref_node_bary4( ref_node, nodes, xyz, bary ), "bary");
  
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

  ref_malloc_init( ref_interp->cell, 
		   ref_node_max(to_node), 
		   REF_INT, REF_EMPTY );
  ref_malloc( ref_interp->bary, 
	      4*ref_node_max(to_node), 
	      REF_DBL );

  each_ref_node_valid_node( to_node, node )
    {
      if ( REF_EMPTY != ref_interp->cell[node] )
	continue;
      RSS(ref_interp_exhaustive_enclosing_tet( from_grid,
					       ref_node_xyz_ptr(to_node,node),
					       &(ref_interp->cell[node]), 
					       &(ref_interp->bary[4*node]) ), 
	  "exhast");
      (ref_interp->nexhaustive)++;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_max_error( REF_INTERP ref_interp, 
				 REF_GRID from_grid, REF_GRID to_grid,
				 REF_DBL *max_error)
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL xyz[3], error;
  REF_INT i;

  *max_error = 0.0;

  RNS( ref_interp->cell, "locate first" );
  RNS( ref_interp->bary, "locate first" );

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, node )
    {
      RSS( ref_cell_nodes( from_cell, ref_interp->cell[node], nodes), "cn" );  
      for(i=0;i<3;i++)
	xyz[i] = 
	  ref_interp->bary[0+4*node] *
	  ref_node_xyz(from_node,i,nodes[0]) +
	  ref_interp->bary[1+4*node] *
	  ref_node_xyz(from_node,i,nodes[1]) +
	  ref_interp->bary[2+4*node] *
	  ref_node_xyz(from_node,i,nodes[2]) +
	  ref_interp->bary[3+4*node] *
	  ref_node_xyz(from_node,i,nodes[3]);
      error = 
	pow(xyz[0]-ref_node_xyz(to_node,0,node),2) + 
	pow(xyz[1]-ref_node_xyz(to_node,1,node),2) + 
	pow(xyz[2]-ref_node_xyz(to_node,2,node),2) ;
      *max_error = MAX( *max_error, sqrt(error) );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_stats( REF_INTERP ref_interp, 
			     REF_GRID from_grid, REF_GRID to_grid)
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_CELL from_cell = ref_grid_tet(from_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT node;
  REF_DBL xyz[3], error;
  REF_INT i;
  REF_DBL max_error = 0.0;
  REF_DBL min_bary = 1.0;

  RNS( ref_interp->cell, "locate first" );
  RNS( ref_interp->bary, "locate first" );

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, node )
    {
      RSS( ref_cell_nodes( from_cell, ref_interp->cell[node], nodes), "cn" );  
      for(i=0;i<3;i++)
	xyz[i] = 
	  ref_interp->bary[0+4*node] *
	  ref_node_xyz(from_node,i,nodes[0]) +
	  ref_interp->bary[1+4*node] *
	  ref_node_xyz(from_node,i,nodes[1]) +
	  ref_interp->bary[2+4*node] *
	  ref_node_xyz(from_node,i,nodes[2]) +
	  ref_interp->bary[3+4*node] *
	  ref_node_xyz(from_node,i,nodes[3]);
      error = 
	pow(xyz[0]-ref_node_xyz(to_node,0,node),2) + 
	pow(xyz[1]-ref_node_xyz(to_node,1,node),2) + 
	pow(xyz[2]-ref_node_xyz(to_node,2,node),2) ;
      max_error = MAX( max_error, sqrt(error) );
      for(i=0;i<4;i++)
	min_bary= MIN( min_bary, ref_interp->bary[i+4*node] );
    }

  if ( ref_mpi_once(ref_mpi) )
    {
      printf("interp min bary %e max error %e\n", min_bary, max_error);
      printf("exhastive %5.1f%% %d of %d\n", 
	     (REF_DBL)ref_interp->nexhaustive / 
	     (REF_DBL)ref_node_n(to_node) * 100.0,
	     ref_interp->nexhaustive, ref_node_n(to_node));
    }

  return REF_SUCCESS;
}
