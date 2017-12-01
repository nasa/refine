
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

#define MAX_NODE_LIST ( 100 )

REF_STATUS ref_interp_create( REF_INTERP *ref_interp_ptr )
{
  REF_INTERP ref_interp;
  
  ref_malloc( *ref_interp_ptr, 1, REF_INTERP_STRUCT );
  ref_interp = ( *ref_interp_ptr );

  ref_interp->nexhaustive = 0;
  ref_interp->nwalk = 0;
  ref_interp->nfail = 0;
  ref_interp->steps = 0;
  ref_interp->wasted = 0;
  ref_interp->ngeom = 0;
  ref_interp->ngeomfail = 0;
  ref_interp->guess = NULL;
  ref_interp->cell = NULL;
  ref_interp->bary = NULL;
  ref_interp->inside = -1.0e-12; /* inside tolerence */

  RSS( ref_list_create( &( ref_interp->ref_list ) ), "add list" );

  return REF_SUCCESS;
}

REF_STATUS ref_interp_free( REF_INTERP ref_interp )
{
  if ( NULL == (void *)ref_interp )
    return REF_NULL;
  ref_list_free( ref_interp->ref_list );
  ref_free( ref_interp->bary );
  ref_free( ref_interp->cell );
  ref_free( ref_interp->guess );
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

REF_STATUS ref_interp_exhaustive_tet_around_node( REF_GRID ref_grid,
						  REF_INT node, REF_DBL *xyz,
						  REF_INT *cell, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, guess, best_guess;
  REF_DBL current_bary[4];
  REF_DBL best_bary, min_bary;
 
  best_guess = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_having_node( ref_cell, node, item, guess )
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

static REF_STATUS ref_update_tet_guess( REF_CELL ref_cell,
					REF_INT node0,
					REF_INT node1,
					REF_INT node2,
					REF_INT *guess)
{
  REF_INT face_nodes[4], cell0, cell1;

  face_nodes[0]=node0;
  face_nodes[1]=node1;
  face_nodes[2]=node2;
  face_nodes[3]=node0;

  RSS( ref_cell_with_face( ref_cell, face_nodes, &cell0, &cell1 ), "next" );
  if ( REF_EMPTY == cell0 )
    THROW("bary update missing first");
  if ( REF_EMPTY == cell1 )
    { /* hit boundary */
      *guess = REF_EMPTY;
      return REF_SUCCESS;
    }

  if ( *guess == cell0 )
    {
      *guess = cell1;
      return REF_SUCCESS;
    }
  if ( *guess == cell1 )
    {
      *guess = cell0;
      return REF_SUCCESS;
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_interp_extrap( REF_INTERP ref_interp, REF_GRID ref_grid, 
			      REF_DBL *xyz, REF_INT guess,
			      REF_INT *tet, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT step, limit;

  *tet = REF_EMPTY;

  limit = 1000; /* was 10e6^(1/3), required 108 for twod testcase  */

  for ( step=0; step < limit; step++)
    {
      /* give up if cell is invalid */
      if ( !ref_cell_valid(ref_cell,guess) )
	{
	  printf("out %d, tet %d, bary %e %e %e %e inside %e\n",
		 step,guess,bary[0],bary[1],bary[2],bary[3],ref_interp->inside);
	  return REF_SUCCESS;
	}

      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RXS( ref_node_bary4( ref_node, nodes, xyz, bary ), REF_DIV_ZERO, "bary");

      if ( step > 990 )
	{
	  printf("step %d, tet %d, bary %e %e %e %e inside %e\n",
		 step,guess,bary[0],bary[1],bary[2],bary[3],ref_interp->inside);
	}
      
      if ( bary[0] >= ref_interp->inside &&
	   bary[1] >= ref_interp->inside &&
	   bary[2] >= ref_interp->inside &&
	   bary[3] >= ref_interp->inside )
	{
	  (ref_interp->steps) += (step+1);
	  (ref_interp->nwalk)++;
	  *tet = guess;
	  return REF_SUCCESS;
	}
      
      /* less than */
      if ( bary[0] < bary[1] && bary[0] < bary[2] && bary[0] < bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[1], nodes[2], nodes[3],
				     &guess ), "update 1 2 3");
	  continue;
	}

      if ( bary[1] < bary[0] && bary[1] < bary[3] && bary[1] < bary[2] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[3], nodes[2],
				     &guess ), "update 0 3 2");
	  continue;
	}

      if ( bary[2] < bary[0] && bary[2] < bary[1] && bary[2] < bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[1], nodes[3],
				     &guess ), "update 0 1 3");
	  continue;
	}

      if ( bary[3] < bary[0] && bary[3] < bary[2] && bary[3] < bary[1] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[2], nodes[1],
				     &guess ), "update 0 2 1");
	  continue;
	}
      
      /* less than or equal */
      if ( bary[0] <= bary[1] && bary[0] <= bary[2] && bary[0] <= bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[1], nodes[2], nodes[3],
				     &guess ), "update 1 2 3");
	  continue;
	}

      if ( bary[1] <= bary[0] && bary[1] <= bary[3] && bary[1] <= bary[2] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[3], nodes[2],
				     &guess ), "update 0 3 2");
	  continue;
	}

      if ( bary[2] <= bary[0] && bary[2] <= bary[1] && bary[2] <= bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[1], nodes[3],
				     &guess ), "update 0 1 3");
	  continue;
	}

      if ( bary[3] <= bary[0] && bary[3] <= bary[2] && bary[3] <= bary[1] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[2], nodes[1],
				     &guess ), "update 0 2 1");
	  continue;
	}

      THROW("unable to find the next step");
    }
  
  THROW("out of iterations");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_enclosing_tet( REF_INTERP ref_interp, REF_GRID ref_grid, 
				     REF_DBL *xyz, REF_INT guess,
				     REF_INT *tet, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT step, limit;

  *tet = REF_EMPTY;

  limit = 1000; /* was 10e6^(1/3), required 108 for twod testcase  */

  for ( step=0; step < limit; step++)
    {
      /* give up if cell is invalid */
      if ( !ref_cell_valid(ref_cell,guess) )
	{
	  (ref_interp->nfail)++;
	  (ref_interp->wasted)+=(step+1);
	  return REF_SUCCESS;
	}

      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RXS( ref_node_bary4( ref_node, nodes, xyz, bary ), REF_DIV_ZERO, "bary");

      if ( step > 990 )
	{
	  printf("step %d, tet %d, bary %e %e %e %e inside %e\n",
		 step,guess,bary[0],bary[1],bary[2],bary[3],ref_interp->inside);
	}
      
      if ( bary[0] >= ref_interp->inside &&
	   bary[1] >= ref_interp->inside &&
	   bary[2] >= ref_interp->inside &&
	   bary[3] >= ref_interp->inside )
	{
	  (ref_interp->steps) += (step+1);
	  (ref_interp->nwalk)++;
	  *tet = guess;
	  return REF_SUCCESS;
	}
      
      /* less than */
      if ( bary[0] < bary[1] && bary[0] < bary[2] && bary[0] < bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[1], nodes[2], nodes[3],
				     &guess ), "update 1 2 3");
	  continue;
	}

      if ( bary[1] < bary[0] && bary[1] < bary[3] && bary[1] < bary[2] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[3], nodes[2],
				     &guess ), "update 0 3 2");
	  continue;
	}

      if ( bary[2] < bary[0] && bary[2] < bary[1] && bary[2] < bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[1], nodes[3],
				     &guess ), "update 0 1 3");
	  continue;
	}

      if ( bary[3] < bary[0] && bary[3] < bary[2] && bary[3] < bary[1] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[2], nodes[1],
				     &guess ), "update 0 2 1");
	  continue;
	}
      
      /* less than or equal */
      if ( bary[0] <= bary[1] && bary[0] <= bary[2] && bary[0] <= bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[1], nodes[2], nodes[3],
				     &guess ), "update 1 2 3");
	  continue;
	}

      if ( bary[1] <= bary[0] && bary[1] <= bary[3] && bary[1] <= bary[2] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[3], nodes[2],
				     &guess ), "update 0 3 2");
	  continue;
	}

      if ( bary[2] <= bary[0] && bary[2] <= bary[1] && bary[2] <= bary[3] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[1], nodes[3],
				     &guess ), "update 0 1 3");
	  continue;
	}

      if ( bary[3] <= bary[0] && bary[3] <= bary[2] && bary[3] <= bary[1] )
	{
	  RSS( ref_update_tet_guess( ref_cell, nodes[0], nodes[2], nodes[1],
				     &guess ), "update 0 2 1");
	  continue;
	}

      THROW("unable to find the next step");
    }
  
  THROW("out of iterations");

  return REF_SUCCESS;
}

REF_STATUS ref_interp_push_onto_queue( REF_INTERP ref_interp, 
				       REF_GRID ref_grid, REF_INT node )
{ 
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT neighbor, nneighbor, neighbors[MAX_NODE_LIST];
  REF_INT other;

  RUS( REF_EMPTY, ref_interp->cell[node], "no cell for guess" );

  RSS( ref_cell_node_list_around( ref_cell, node, MAX_NODE_LIST,
                                  &nneighbor, neighbors ), "list too small");
  for ( neighbor = 0; neighbor < nneighbor; neighbor++ )
    {
      other = neighbors[neighbor];
      if ( ref_interp->cell[other] == REF_EMPTY &&
	   ref_interp->guess[other] == REF_EMPTY )
	{
	  RSS(ref_list_add(ref_interp->ref_list, other), "enqueue" );
	  ref_interp->guess[other] = ref_interp->cell[node];
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_drain_queue( REF_INTERP ref_interp, 
				   REF_GRID from_grid, REF_GRID to_grid )
{
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT node;

  while ( 0 < ref_list_n( ref_interp->ref_list ) )
    {
      RSS( ref_list_pop( ref_interp->ref_list, &node ), "pop queue");
      RUS( REF_EMPTY, ref_interp->guess[node], "no guess" );
      REIS( REF_EMPTY, ref_interp->cell[node], "queued to node already found?");
      RSS( ref_interp_enclosing_tet( ref_interp,
				     from_grid,
				     ref_node_xyz_ptr(to_node,node),
				     ref_interp->guess[node],
				     &(ref_interp->cell[node]),
				     &(ref_interp->bary[4*node]) ), 
	   "walk");
      if ( REF_EMPTY == ref_interp->cell[node] )
	{
	  ref_interp->guess[node] = REF_EMPTY;
	}
      else
	{
	  RSS( ref_interp_push_onto_queue(ref_interp,to_grid,node), "push" ); 
	}
    }
  return REF_SUCCESS;
}

REF_STATUS ref_interp_geom_nodes( REF_INTERP ref_interp,
				  REF_GRID from_grid, REF_GRID to_grid );
REF_STATUS ref_interp_try_adj( REF_INTERP ref_interp,
			       REF_GRID from_grid, REF_GRID to_grid );

REF_STATUS ref_interp_locate( REF_INTERP ref_interp, 
			      REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT node;

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  ref_malloc_init( ref_interp->guess, 
		   ref_node_max(to_node), 
		   REF_INT, REF_EMPTY );
  ref_malloc_init( ref_interp->cell, 
		   ref_node_max(to_node), 
		   REF_INT, REF_EMPTY );
  ref_malloc( ref_interp->bary, 
	      4*ref_node_max(to_node), 
	      REF_DBL );

  RSS( ref_interp_geom_nodes( ref_interp, from_grid, to_grid ), "geom nodes");
  
  RSS( ref_interp_drain_queue( ref_interp, from_grid, to_grid), "drain" );

  RSS( ref_interp_try_adj( ref_interp, from_grid, to_grid), "adj" );

  each_ref_node_valid_node( to_node, node )
    {
      if ( REF_EMPTY != ref_interp->cell[node] )
	continue;
      RSS(ref_interp_exhaustive_enclosing_tet( from_grid,
					       ref_node_xyz_ptr(to_node,node),
					       &(ref_interp->cell[node]), 
					       &(ref_interp->bary[4*node]) ), 
	  "exhaust");
      (ref_interp->nexhaustive)++;
      RSS( ref_interp_push_onto_queue(ref_interp,to_grid,node), "push" ); 
      RSS( ref_interp_drain_queue( ref_interp, from_grid, to_grid), "drain" );
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
      RSS( ref_cell_nodes( from_cell, ref_interp->cell[node], nodes),
	   "node needs to be localized" );
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
  REF_DBL this_bary;
  REF_DBL min_bary = 1.0;
  REF_INT extrapolate = 0;

  RNS( ref_interp->cell, "locate first" );
  RNS( ref_interp->bary, "locate first" );

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, node )
    {
      RSS( ref_cell_nodes( from_cell, ref_interp->cell[node], nodes),
	   "node needs to be localized" );
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
      this_bary = MIN( MIN( ref_interp->bary[0+4*node],
			    ref_interp->bary[1+4*node] ),
		       MIN( ref_interp->bary[2+4*node],
			    ref_interp->bary[3+4*node] ) );
      min_bary= MIN( min_bary, this_bary );
      if ( this_bary < ref_interp->inside ) 
	extrapolate++;
    }

  if ( ref_mpi_once(ref_mpi) )
    {
      printf("interp min bary %e max error %e extrap %d\n", 
	     min_bary, max_error, extrapolate );
      printf("exhastive %5.1f%% %d of %d\n", 
	     (REF_DBL)ref_interp->nexhaustive / 
	     (REF_DBL)ref_node_n(to_node) * 100.0,
	     ref_interp->nexhaustive, ref_node_n(to_node));
      printf("walks: failed %d of %f, successful %d of %f\n",
	     ref_interp->nfail, 
	     (REF_DBL)ref_interp->wasted / (REF_DBL)ref_interp->nfail,
	     ref_interp->nwalk, 
	     (REF_DBL)ref_interp->steps / (REF_DBL)ref_interp->nwalk );
      printf("geom nodes: failed %d, successful %d, of %d\n",
	     ref_interp->ngeomfail,  ref_interp->ngeom,
	     ref_interp->ngeom+ref_interp->ngeomfail);
    }

  return REF_SUCCESS;
}

REF_STATUS ref_interp_geom_nodes( REF_INTERP ref_interp,
				  REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_CELL to_tri = ref_grid_tri(to_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tri = ref_grid_tri(from_grid);
  REF_LIST ref_list;
  REF_INT to_n, from_n;
  REF_INT nfaceid, faceids[3];
  REF_DBL *xyz;
  REF_INT item;
  REF_INT i, best_item, cell;
  REF_DBL dist, best_dist, bary[4];
  RSS( ref_list_create( &ref_list ), "create list" );
  
  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, to_n )
    {
      if ( !ref_cell_node_empty( to_tri, to_n ) )
	{
	  RXS( ref_cell_faceid_list_around( to_tri,
					    to_n,
					    3,
					    &nfaceid, faceids ),
	       REF_INCREASE_LIMIT, "count faceids" );
	  if ( nfaceid >= 3 )
	    RSS( ref_list_add( ref_list, to_n ), "add geom node" );
	}
    }
	      
  each_ref_node_valid_node( from_node, from_n )
    {
      if ( !ref_cell_node_empty( from_tri, from_n ) )
	{
	  RXS( ref_cell_faceid_list_around( from_tri,
					    from_n,
					    3,
					    &nfaceid, faceids ),
	       REF_INCREASE_LIMIT, "count faceids" );
	  if ( nfaceid >= 3 )
	    {
	      best_dist = 1.0e20;
	      best_item = REF_EMPTY;
	      each_ref_list_item( ref_list, item )
		{
		  xyz = ref_node_xyz_ptr(to_node,ref_list_value(ref_list,item));
		  dist =
		    pow(xyz[0]-ref_node_xyz(from_node,0,from_n),2) +
		    pow(xyz[1]-ref_node_xyz(from_node,1,from_n),2) +
		    pow(xyz[2]-ref_node_xyz(from_node,2,from_n),2) ;
		  dist = sqrt(dist);
		  if ( dist < best_dist )
		    {
		      best_dist = dist;
		      best_item = item;
		    }
		}
	      to_n = ref_list_value(ref_list,best_item);
	      xyz = ref_node_xyz_ptr( to_node, to_n );
	      RSS( ref_interp_exhaustive_tet_around_node( from_grid, from_n,
							  xyz, &cell, bary),
		   "tet around node");
	      if ( bary[0] > ref_interp->inside &&
		   bary[1] > ref_interp->inside &&
		   bary[2] > ref_interp->inside &&
		   bary[3] > ref_interp->inside )
		{
		  ref_interp->ngeom++;
		  REIS( REF_EMPTY, ref_interp->cell[to_n],
			"geom already found?" );
		  if ( REF_EMPTY != ref_interp->guess[to_n] )
		    { /* need to dequeue */
		      ref_interp->guess[to_n] = REF_EMPTY;
		      RSS( ref_list_delete( ref_interp->ref_list, to_n ),"deq");
		    }
		  ref_interp->cell[to_n] = cell;
		  for(i=0;i<4;i++)
		    ref_interp->bary[i+4*to_n] = bary[i];
		  RSS( ref_interp_push_onto_queue(ref_interp,to_grid,to_n),
		       "push" );
		}
	      else
		{
		  ref_interp->ngeomfail++;
		}
	    }
	}
    }
  ref_list_free( ref_list );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_try_adj( REF_INTERP ref_interp, 
			       REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_CELL to_tet = ref_grid_tet(to_grid);
  REF_INT node;
  REF_INT neighbor, nneighbor, neighbors[MAX_NODE_LIST];
  REF_INT other;
  REF_INT have, havenot;

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, node )
    {
      if ( REF_EMPTY != ref_interp->cell[node] )
	continue;
      printf("for node %d:\n",node);
      RSS( ref_cell_node_list_around( to_tet, node, MAX_NODE_LIST,
				      &nneighbor, neighbors ), "list too small");
      have = 0;
      havenot=0;
      for ( neighbor = 0; neighbor < nneighbor; neighbor++ )
	{
	  other = neighbors[neighbor];
	  if ( ref_interp->cell[other] != REF_EMPTY )
	    {
	      have++;
	      RSS( ref_interp_extrap( ref_interp,
				      from_grid,
				      ref_node_xyz_ptr(to_node,node),
				      ref_interp->cell[other],
				      &(ref_interp->cell[node]),
				      &(ref_interp->bary[4*node]) ), 
		   "walk");
	      if ( ref_interp->cell[node] != REF_EMPTY )
		{
	  printf("tet %d, bary %e %e %e %e\n",ref_interp->cell[node],
		 ref_interp->bary[0+4*node],
		 ref_interp->bary[1+4*node],
		 ref_interp->bary[2+4*node],
		 ref_interp->bary[3+4*node]);
		  printf("got it!\n");
		  break;
		}
	    }
	  else
	    {
	      havenot++;
	    }
	}
      printf("have %d havenot %d\n",have,havenot);
      if ( REF_EMPTY == ref_interp->cell[node] )
	{
	  RSS(ref_interp_exhaustive_enclosing_tet( from_grid,
						   ref_node_xyz_ptr(to_node,node),
						   &(ref_interp->cell[node]), 
						   &(ref_interp->bary[4*node]) ), 
	      "exhaust");
	  printf("tet %d, bary %e %e %e %e\n",ref_interp->cell[node],
		 ref_interp->bary[0+4*node],
		 ref_interp->bary[1+4*node],
		 ref_interp->bary[2+4*node],
		 ref_interp->bary[3+4*node]);
	  (ref_interp->nexhaustive)++;
	  RSS( ref_interp_push_onto_queue(ref_interp,to_grid,node), "push" ); 
	  RSS( ref_interp_drain_queue( ref_interp, from_grid, to_grid), "drain" );
	}
      if ( REF_EMPTY != ref_interp->cell[node] )
	{
	  RSS( ref_interp_push_onto_queue(ref_interp,to_grid,node), "push" ); 
	  RSS( ref_interp_drain_queue( ref_interp, from_grid, to_grid),
	       "drain" );
	}
    }

  return REF_SUCCESS;
}
