
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

#include "ref_search.h"

#include "ref_malloc.h"

#define MAX_NODE_LIST ( 100 )

REF_STATUS ref_interp_create( REF_INTERP *ref_interp_ptr )
{
  REF_INTERP ref_interp;
  
  ref_malloc( *ref_interp_ptr, 1, REF_INTERP_STRUCT );
  ref_interp = ( *ref_interp_ptr );

  ref_interp->instrument = REF_FALSE;
  ref_interp->n_walk = 0;
  ref_interp->walk_steps = 0;
  ref_interp->n_geom = 0;
  ref_interp->n_geom_fail = 0;
  ref_interp->n_tree = 0;
  ref_interp->tree_cells = 0;
  ref_interp->guess = NULL;
  ref_interp->cell = NULL;
  ref_interp->bary = NULL;
  ref_interp->inside = -1.0e-12; /* inside tolerence */
  ref_interp->bound = -0.1; /* bound tolerence */

  RSS( ref_list_create( &( ref_interp->ref_list ) ), "add list" );
  RSS( ref_list_create( &( ref_interp->visualize ) ), "add list" );

  return REF_SUCCESS;
}

REF_STATUS ref_interp_free( REF_INTERP ref_interp )
{
  if ( NULL == (void *)ref_interp )
    return REF_NULL;
  ref_list_free( ref_interp->visualize );
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

REF_STATUS ref_interp_enclosing_tet_in_list( REF_GRID ref_grid,
					     REF_LIST ref_list,
					     REF_DBL *xyz,
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
  each_ref_list_item( ref_list, item )
    {
      guess = ref_list_value( ref_list, item );
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
	  (ref_interp->walk_steps) += (step+1);
	  (ref_interp->n_walk)++;
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

REF_STATUS ref_interp_geom_node_list( REF_GRID ref_grid, REF_LIST ref_list )
{
  REF_NODE ref_node = ref_grid_node( ref_grid );
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL edg = ref_grid_edg(ref_grid);
  REF_INT nfaceid, faceids[3];
  REF_INT nedgeid, edgeids[2];
  REF_INT node;
  each_ref_node_valid_node( ref_node, node )
    {
      if ( ref_node_owned(ref_node,node) )
	{
	  RXS( ref_cell_faceid_list_around( edg,
					    node,
					    2,
					    &nedgeid, edgeids ),
	       REF_INCREASE_LIMIT, "count faceids" );
	  RXS( ref_cell_faceid_list_around( tri,
					    node,
					    3,
					    &nfaceid, faceids ),
	       REF_INCREASE_LIMIT, "count faceids" );
	  if ( nfaceid >= 3 || nedgeid >= 2 )
	    RSS( ref_list_add( ref_list, node ), "add geom node" );
	}
    }
  return REF_SUCCESS;
}

REF_STATUS ref_interp_geom_nodes( REF_INTERP ref_interp,
				  REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_LIST to_geom_list, from_geom_list;
  REF_INT to_geom_node, from_geom_node;
  REF_INT to_item, from_item;
  REF_DBL *xyz;
  REF_DBL *local_xyz, *global_xyz;
  REF_INT total_xyz, i, *best_node, *best_proc, cell;
  REF_DBL dist, *best_dist, bary[4];

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  RSS( ref_list_create( &to_geom_list ), "create list" );
  RSS( ref_list_create( &from_geom_list ), "create list" );
  RSS( ref_interp_geom_node_list( to_grid, to_geom_list ), "to list" );
  RSS( ref_interp_geom_node_list( from_grid, from_geom_list ), "from list" );

  ref_malloc( local_xyz, 3*ref_list_n(to_geom_list), REF_DBL );
  each_ref_list_item_value( to_geom_list, to_item, to_geom_node )
    {
      local_xyz[0+3*to_item] = ref_node_xyz(to_node,0,to_geom_node);
      local_xyz[1+3*to_item] = ref_node_xyz(to_node,1,to_geom_node);
      local_xyz[2+3*to_item] = ref_node_xyz(to_node,2,to_geom_node);
    }
  RSS( ref_mpi_allconcat( ref_mpi, 3*ref_list_n(to_geom_list), 
			  (void *)local_xyz,
			  &total_xyz, (void **)&global_xyz, 
			  REF_DBL_TYPE ), "cat");
  
  ref_malloc( best_dist, total_xyz/3, REF_DBL );
  ref_malloc( best_node, total_xyz/3, REF_INT );
  ref_malloc( best_proc, total_xyz/3, REF_INT );
  for ( to_item = 0; to_item < total_xyz/3; to_item++ )
    {
      xyz = &(global_xyz[3*to_item]);
      best_dist[to_item] = 1.0e20;
      best_node[to_item] = REF_EMPTY;
      each_ref_list_item_value( from_geom_list, from_item, from_geom_node )
	{
	  dist =
	    pow(xyz[0]-ref_node_xyz(from_node,0,from_geom_node),2) +
	    pow(xyz[1]-ref_node_xyz(from_node,1,from_geom_node),2) +
	    pow(xyz[2]-ref_node_xyz(from_node,2,from_geom_node),2) ;
	  dist = sqrt(dist);
	  if ( dist < best_dist[to_item] || 0 == from_item )
	    {
	      best_dist[to_item] = dist;
	      best_node[to_item] = from_geom_node;
	    }
	}
    }

  RSS( ref_mpi_allminwho( ref_mpi, best_dist, best_proc, total_xyz/3), "who" );

  for ( to_item = 0; to_item < total_xyz/3; to_item++ )
    if ( ref_mpi_rank(ref_mpi) == best_proc[to_item] )
      {
	RUS( REF_EMPTY, best_node[to_item], "no geom node" );
	xyz = &(global_xyz[3*to_item]);
	RSS( ref_interp_exhaustive_tet_around_node( from_grid, 
						    best_node[to_item],
						    xyz, &cell, bary),
	     "tet around node");
	if ( bary[0] > ref_interp->inside &&
	     bary[1] > ref_interp->inside &&
	     bary[2] > ref_interp->inside &&
	     bary[3] > ref_interp->inside )
	  {
	    ref_interp->n_geom++;
	    to_geom_node = ref_list_value(to_geom_list,to_item);
	    REIS( REF_EMPTY, ref_interp->cell[to_geom_node],
		  "geom already found?" );
	    if ( REF_EMPTY != ref_interp->guess[to_geom_node] )
	      { /* need to dequeue */
		ref_interp->guess[to_geom_node] = REF_EMPTY;
		RSS( ref_list_delete( ref_interp->ref_list, to_geom_node ),"deq");
	      }
	    ref_interp->cell[to_geom_node] = cell;
	    for(i=0;i<4;i++)
	      ref_interp->bary[i+4*to_geom_node] = bary[i];
	    RSS( ref_interp_push_onto_queue(ref_interp,to_grid,to_geom_node),
		 "push" );
	  }
	else
	  {
	    ref_interp->n_geom_fail++;
	  }
      }

  ref_free( best_proc );
  ref_free( best_node );
  ref_free( best_dist );

  ref_free( global_xyz );
  ref_free( local_xyz );
  ref_list_free( from_geom_list );
  ref_list_free( to_geom_list );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_bounding_sphere( REF_NODE ref_node, REF_INT *nodes,
				       REF_DBL *center, REF_DBL *radius )
{
  REF_INT i;
  for(i=0;i<3;i++)
    center[i] = 0.25 * ( ref_node_xyz(ref_node,i,nodes[0]) +
			 ref_node_xyz(ref_node,i,nodes[1]) +
			 ref_node_xyz(ref_node,i,nodes[2]) +
			 ref_node_xyz(ref_node,i,nodes[3]) );
  *radius = 0.0;
  for(i=0;i<4;i++)
    *radius = MAX( *radius, sqrt( pow( ref_node_xyz(ref_node,0,nodes[i]) -
				       center[0], 2) +
				  pow( ref_node_xyz(ref_node,1,nodes[i]) -
				       center[1], 2) +
				  pow( ref_node_xyz(ref_node,2,nodes[i]) -
				       center[2], 2) ) );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_tree( REF_INTERP ref_interp, 
			    REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE from_node = ref_grid_node(from_grid);
  REF_CELL from_tet = ref_grid_tet(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_SEARCH ref_search;
  REF_INT node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL center[3], radius;
  REF_LIST ref_list;
  REF_DBL fuzz = 1.0e-12;
  REF_BOOL in_search;
  REF_BOOL verify = REF_FALSE;
  
  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  RSS( ref_list_create( &ref_list ), "create list" );
  RSS( ref_search_create( &ref_search, ref_cell_n(from_tet) ), "mk sr" );
  each_ref_cell_valid_cell_with_nodes( from_tet, cell, nodes )
    {
      RSS( ref_interp_bounding_sphere( from_node, nodes, center,&radius ), "b");
      RSS( ref_search_insert( ref_search, cell, center, 2.0*radius ), "ins" );
    }

  each_ref_node_valid_node( to_node, node )
    {
      if ( REF_EMPTY != ref_interp->cell[node] )
	continue;
      RSS( ref_search_touching( ref_search, ref_list,
				ref_node_xyz_ptr(to_node,node), fuzz ), "tch" );
      RSS( ref_interp_enclosing_tet_in_list( from_grid, ref_list,
					     ref_node_xyz_ptr(to_node,node),
					     &(ref_interp->cell[node]), 
					     &(ref_interp->bary[4*node])),
	   "best in list");
      (ref_interp->n_tree)++;
      (ref_interp->tree_cells)+=ref_list_n(ref_list);
      if ( verify )
	{
	  REF_INT verify_cell;
	  REF_DBL verify_bary[4];
	  RSS( ref_interp_exhaustive_enclosing_tet( from_grid,
						    ref_node_xyz_ptr(to_node,
								     node),
						    &verify_cell, 
						    verify_bary),
	       "exhaust verification");
	  RSS( ref_list_contains( ref_list, verify_cell, &in_search ), "hav");
	  if ( !in_search ) printf("didn't get it :(");
	  REIS( verify_cell, ref_interp->cell[node], "bad validation" );
	}
      RSS( ref_list_erase(ref_list), "reset list" );
    }
  RSS( ref_search_free(ref_search), "free list" );
  RSS( ref_list_free(ref_list), "free list" );
  return REF_SUCCESS;
}

REF_STATUS ref_interp_locate( REF_INTERP ref_interp, 
			      REF_GRID from_grid, REF_GRID to_grid )
{
  REF_MPI ref_mpi = ref_grid_mpi(from_grid);
  REF_NODE to_node = ref_grid_node(to_grid);

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

  if ( ref_interp->instrument)
    RSS( ref_mpi_stopwatch_start( ref_mpi ), "locate clock");

  RSS( ref_interp_geom_nodes( ref_interp, from_grid, to_grid ), "geom nodes");
  if ( ref_interp->instrument)
    RSS( ref_mpi_stopwatch_stop( ref_mpi, "geom" ), "locate clock");
  
  RSS( ref_interp_drain_queue( ref_interp, from_grid, to_grid), "drain" );
  if ( ref_interp->instrument)
    RSS( ref_mpi_stopwatch_stop( ref_mpi, "drain" ), "locate clock");

  RSS( ref_interp_tree( ref_interp, from_grid, to_grid), "tree" );
  if ( ref_interp->instrument)
    RSS( ref_mpi_stopwatch_stop( ref_mpi, "tree" ), "locate clock");
  
  return REF_SUCCESS;
}

REF_STATUS ref_interp_min_bary( REF_INTERP ref_interp, 
				REF_GRID to_grid,
				REF_DBL *min_bary )
{
  REF_MPI ref_mpi = ref_grid_mpi(to_grid);
  REF_NODE to_node = ref_grid_node(to_grid);
  REF_INT node;
  REF_DBL this_bary;
  
  *min_bary = 1.0;

  RNS( ref_interp->cell, "locate first" );
  RNS( ref_interp->bary, "locate first" );

  if ( ref_mpi_para(ref_mpi) )
    RSS( REF_IMPLEMENT, "not para" );

  each_ref_node_valid_node( to_node, node )
    {
      RUS( REF_EMPTY, ref_interp->cell[node], "node needs to be localized");
      this_bary = MIN( MIN( ref_interp->bary[0+4*node],
			    ref_interp->bary[1+4*node] ),
		       MIN( ref_interp->bary[2+4*node],
			    ref_interp->bary[3+4*node] ) );
      *min_bary= MIN( *min_bary, this_bary );
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
      printf("tree search: %d found, %.2f avg cells\n",
	     ref_interp->n_tree,
	     (REF_DBL)ref_interp->tree_cells / (REF_DBL)ref_interp->n_tree);
      printf("walks: %d successful, %.2f avg cells\n",
	     ref_interp->n_walk, 
	     (REF_DBL)ref_interp->walk_steps / (REF_DBL)ref_interp->n_walk );
      printf("geom nodes: %d failed, %d successful\n",
	     ref_interp->n_geom_fail, ref_interp->n_geom);
      printf("total: %d nodes, from %d tet\n",
	     ref_node_n(to_node), ref_cell_n(ref_grid_tet(from_grid)) );
    }

  return REF_SUCCESS;
}


REF_STATUS ref_interp_tec( REF_INTERP ref_interp, 
			   REF_GRID to_grid, const char *filename )
{
  REF_NODE ref_node = ref_grid_node(to_grid); 
  FILE *file;
  REF_INT item;

  /* skip if noting to show */
  if ( 0 == ref_list_n(ref_interp->visualize) )
    return REF_SUCCESS;
  
  file = fopen(filename,"w");
  if (NULL == (void *)file) printf("unable to open %s\n",filename);
  RNS(file, "unable to open file" );

  fprintf(file, "title=\"refine interp\"\n");
  fprintf(file, "variables = \"x\" \"y\" \"z\"\n");

  fprintf(file,
	  "zone t=exhaust, i=%d, datapacking=%s\n",
	  ref_list_n(ref_interp->visualize), "point");

  each_ref_list_item( ref_interp->visualize, item )
    fprintf(file,"%.15e %.15e %.15e\n",
	    ref_node_xyz(ref_node,0,ref_list_value(ref_interp->visualize,item)),
	    ref_node_xyz(ref_node,1,ref_list_value(ref_interp->visualize,item)),
	    ref_node_xyz(ref_node,2,ref_list_value(ref_interp->visualize,item))
	    );

  fclose(file);
  return REF_SUCCESS;
}
