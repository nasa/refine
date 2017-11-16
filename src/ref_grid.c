
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

#include "ref_grid.h"
#include "ref_adj.h"

#include "ref_malloc.h"
#include "ref_matrix.h"

REF_STATUS ref_grid_create( REF_GRID *ref_grid_ptr )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  RSS( ref_mpi_create( &ref_grid_mpi(ref_grid) ), "mpi create" );

  RSS( ref_node_create( &ref_grid_node(ref_grid), 
			ref_grid_mpi(ref_grid) ), "node create" );

  RSS( ref_cell_create( &ref_grid_tet(ref_grid), 4, REF_FALSE ), "tet create" );
  RSS( ref_cell_create( &ref_grid_pyr(ref_grid), 5, REF_FALSE ), "pyr create" );
  RSS( ref_cell_create( &ref_grid_pri(ref_grid), 6, REF_FALSE ), "pri create" );
  RSS( ref_cell_create( &ref_grid_hex(ref_grid), 8, REF_FALSE ), "hex create" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_create( &ref_grid_edg(ref_grid), 2, REF_TRUE ), "edg create" );
  RSS( ref_cell_create( &ref_grid_tri(ref_grid), 3, REF_TRUE ), "tri create" );
  RSS( ref_cell_create( &ref_grid_qua(ref_grid), 4, REF_TRUE ), "qua create" );

  RSS( ref_geom_create( &ref_grid_geom(ref_grid) ), "geom create" );
  RSS( ref_gather_create( &ref_grid_gather(ref_grid) ), "gather create" );
  RSS( ref_adapt_create( &(ref_grid->adapt) ), "adapt create" );

  ref_grid_twod(ref_grid) = REF_FALSE;

  return REF_SUCCESS;
}

REF_STATUS ref_grid_deep_copy( REF_GRID *ref_grid_ptr, REF_GRID original )
{
  REF_GRID ref_grid;

  ref_malloc( *ref_grid_ptr, 1, REF_GRID_STRUCT );

  ref_grid = *ref_grid_ptr;

  RSS( ref_node_deep_copy( &ref_grid_node(ref_grid),
			   ref_grid_node(original) ), "node deep copy" );

  RSS( ref_cell_deep_copy( &ref_grid_tet(ref_grid),
			   ref_grid_tet(original) ), "tet deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_pyr(ref_grid),
			   ref_grid_pyr(original) ), "pyr deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_pri(ref_grid),
			   ref_grid_pri(original) ), "pri deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_hex(ref_grid),
			   ref_grid_hex(original) ), "hex deep copy" );

  ref_grid_cell(ref_grid,4) = NULL;

  RSS( ref_cell_deep_copy( &ref_grid_edg(ref_grid),
			   ref_grid_edg(original) ), "edg deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_tri(ref_grid),
			   ref_grid_tri(original) ), "tri deep copy" );
  RSS( ref_cell_deep_copy( &ref_grid_qua(ref_grid),
			   ref_grid_qua(original) ), "qua deep copy" );

  RSS( ref_mpi_deep_copy( &ref_grid_mpi(ref_grid),
			   ref_grid_mpi(original) ), "mpi deep copy" );
  RSS( ref_geom_deep_copy( &ref_grid_geom(ref_grid),
			   ref_grid_geom(original) ), "geom deep copy" );
  RSS( ref_gather_create( &ref_grid_gather(ref_grid) ), "gather create" );
  RSS( ref_adapt_deep_copy( &(ref_grid->adapt),
			    original->adapt ), "adapt deep copy" );
  
  ref_grid_twod(ref_grid) = ref_grid_twod(original);

  return REF_SUCCESS;
}

REF_STATUS ref_grid_free( REF_GRID ref_grid )
{
  if ( NULL == (void *)ref_grid ) return REF_NULL;

  RSS( ref_adapt_free( ref_grid->adapt ), "adapt free");
  RSS( ref_gather_free( ref_grid_gather(ref_grid) ), "gather free");
  RSS( ref_geom_free( ref_grid_geom(ref_grid) ), "geom free");
  RSS( ref_mpi_free( ref_grid_mpi(ref_grid) ), "mpi free");
  
  RSS( ref_cell_free( ref_grid_qua(ref_grid) ), "qua free");
  RSS( ref_cell_free( ref_grid_tri(ref_grid) ), "tri free");
  RSS( ref_cell_free( ref_grid_edg(ref_grid) ), "edg free");

  RSS( ref_cell_free( ref_grid_hex(ref_grid) ), "hex free");
  RSS( ref_cell_free( ref_grid_pri(ref_grid) ), "pri free");
  RSS( ref_cell_free( ref_grid_pyr(ref_grid) ), "pyr free");
  RSS( ref_cell_free( ref_grid_tet(ref_grid) ), "tet free");

  RSS( ref_node_free( ref_grid_node(ref_grid) ), "node free");

  ref_free( ref_grid );
  return REF_SUCCESS;
}

REF_STATUS ref_grid_inspect( REF_GRID ref_grid )
{
  printf("ref_grid = %p\n",(void *)ref_grid);
  printf(" %d node\n",ref_node_n(ref_grid_node(ref_grid)));
  printf(" %d tet\n",ref_cell_n(ref_grid_tet(ref_grid)));
  printf(" %d pyr\n",ref_cell_n(ref_grid_pyr(ref_grid)));
  printf(" %d pri\n",ref_cell_n(ref_grid_pri(ref_grid)));
  printf(" %d hex\n",ref_cell_n(ref_grid_hex(ref_grid)));
  printf(" %d edg\n",ref_cell_n(ref_grid_edg(ref_grid)));
  printf(" %d tri\n",ref_cell_n(ref_grid_tri(ref_grid)));
  printf(" %d qua\n",ref_cell_n(ref_grid_qua(ref_grid)));
  printf(" %d geom\n",ref_geom_n(ref_grid_geom(ref_grid)));
  printf(" %p gather\n",(void *)(ref_grid_gather(ref_grid)->file));
  printf(" %p adapt\n",(void *)(ref_grid->adapt));

  return REF_SUCCESS;
}

REF_STATUS ref_grid_tattle( REF_GRID ref_grid, REF_INT node )
{
  REF_CELL ref_cell;
  REF_INT item, cell;
  
  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_tattle( ref_cell, cell ),"cell tat");
    }
  ref_cell = ref_grid_tet(ref_grid);
  each_ref_cell_having_node( ref_cell, node, item, cell )
    {
      RSS( ref_cell_tattle( ref_cell, cell ),"cell tat");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_grid_cell_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 4:
      *ref_cell = ref_grid_tet(ref_grid);
      break;
    case 5:
      *ref_cell = ref_grid_pyr(ref_grid);
      break;
    case 6:
      *ref_cell = ref_grid_pri(ref_grid);
      break;
    case 8:
      *ref_cell = ref_grid_hex(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_face_with( REF_GRID ref_grid, REF_INT node_per,
			       REF_CELL *ref_cell )
{

  *ref_cell = NULL;

  switch ( node_per )
    {
    case 3:
      *ref_cell = ref_grid_tri(ref_grid);
      break;
    case 4:
      *ref_cell = ref_grid_qua(ref_grid);
      break;
    default:
      printf("node_per %d\n",node_per);
      RSS(REF_FAILURE, "unexpected node_per");
      break;    
    }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_cell_has_face( REF_GRID ref_grid, 
				   REF_INT *face_nodes,
				   REF_BOOL *has_face )
{
  REF_CELL ref_cell;
  REF_INT group;
  REF_INT cell0, cell1;

  *has_face = REF_FALSE;
  
  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    {
      RSS( ref_cell_with_face( ref_cell, face_nodes, &cell0, &cell1),
	   "face search failed" );
      *has_face = ( REF_EMPTY != cell0 );
      if ( *has_face ) return REF_SUCCESS;
    }

 return REF_SUCCESS;
}

REF_STATUS ref_grid_boundary_nodes( REF_GRID ref_grid, 
				    REF_INT boundary_tag, 
				    REF_INT *nnode, REF_INT *nface, 
				    REF_INT **g2l, REF_INT **l2g )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_node = ref_grid_node(ref_grid);

  ref_malloc_init( *g2l, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  (*nnode) = 0;
  (*nface) = 0;

  ref_cell = ref_grid_tri(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( boundary_tag == nodes[3] )
      {
	(*nface)++;
	for ( node = 0; node < 3; node++ )
	  if ( REF_EMPTY == (*g2l)[nodes[node]] )
	    { (*g2l)[nodes[node]] = (*nnode); (*nnode)++; }
      }

  ref_cell = ref_grid_qua(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( boundary_tag == nodes[4] )
      {
	(*nface)++;
	for ( node = 0; node < 4; node++ )
	  if ( REF_EMPTY == (*g2l)[nodes[node]] )
	    { (*g2l)[nodes[node]] = (*nnode); (*nnode)++; }
      }

  ref_malloc( *l2g, *nnode, REF_INT );

  (*nnode) = 0;
  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( REF_EMPTY != (*g2l)[node] ) 
      {
	(*g2l)[node] = (*nnode);
	(*l2g)[(*nnode)] = node;
	(*nnode)++;
      }

  return REF_SUCCESS;
}

REF_STATUS ref_grid_edge_nodes( REF_GRID ref_grid, 
				REF_INT edge_tag, 
				REF_INT *nnode, REF_INT *nedge, 
				REF_INT **g2l, REF_INT **l2g )
{
  REF_NODE ref_node;
  REF_CELL ref_cell;
  REF_INT cell, node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

  ref_node = ref_grid_node(ref_grid);

  ref_malloc_init( *g2l, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  (*nnode) = 0;
  (*nedge) = 0;

  ref_cell = ref_grid_edg(ref_grid);
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
    if ( edge_tag == nodes[2] )
      {
	(*nedge)++;
	for ( node = 0; node < 2; node++ )
	  if ( REF_EMPTY == (*g2l)[nodes[node]] )
	    { (*g2l)[nodes[node]] = (*nnode); (*nnode)++; }
      }

  ref_malloc( *l2g, *nnode, REF_INT );

  (*nnode) = 0;
  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( REF_EMPTY != (*g2l)[node] ) 
      {
	(*g2l)[node] = (*nnode);
	(*l2g)[(*nnode)] = node;
	(*nnode)++;
      }

  return REF_SUCCESS;
}

static REF_STATUS ref_update_tri_guess( REF_CELL ref_cell,
				    REF_INT node0, REF_INT node1,
				    REF_INT *guess)
{
  REF_INT ncell, max_cell = 2, cell_list[2];

  RSS( ref_cell_list_with2( ref_cell, node0, node1,
			    max_cell, &ncell, cell_list ), "next" );
  if ( ncell == 1 )
    THROW("bary update hit boundary");

  if ( *guess == cell_list[0] )
    {
      *guess = cell_list[1];
      return REF_SUCCESS;
    }
  if ( *guess == cell_list[1] )
    {
      *guess = cell_list[0];
      return REF_SUCCESS;
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_grid_identity_interp_guess( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell;
  REF_ADJ ref_adj;
  REF_INT node;
  
  RSS( ref_node_allocate_guess( ref_node ), "alloc");

  if (ref_grid_twod(ref_grid))
    {
      ref_cell = ref_grid_tri(ref_grid);
    }
  else
    {
      ref_cell = ref_grid_tet(ref_grid);
    }
  ref_adj = ref_cell_adj(ref_cell);
  
  each_ref_node_valid_node( ref_node, node )
    {
      ref_node_raw_guess(ref_node,node) =
	ref_adj_item_ref(ref_adj,ref_adj_first(ref_adj,node));
    }
  
  return REF_SUCCESS;
}

REF_STATUS ref_grid_enclosing_tri( REF_GRID ref_grid, REF_DBL *xyz,
                                  REF_INT *tri, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT guess;
  REF_DBL inside;
  REF_INT step, limit;
  
  guess = *tri;
  if ( ! ref_cell_valid(ref_cell,guess) )
    {
      each_ref_cell_valid_cell( ref_cell, guess)
	{
	  break;
	}
    }
  RAS( ref_cell_valid(ref_cell,guess), "unable to find start");

  limit = 1000; /* was 10e6^(1/3), required 108 for twod testcase  */
  inside = -1.0e-12; /* initial inside tolerence  */

  for ( step=0; step < limit; step++)
    {
      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RSS( ref_node_bary3( ref_node, nodes, xyz, bary ), "bary");

      if ( step > 990 )
	{
	  printf("step %d, tri %d, bary %e %e %e inside %e\n",
		 step,guess,bary[0],bary[1],bary[2],inside);
	}
      
      if ( bary[0] >= inside &&
	   bary[1] >= inside &&
	   bary[2] >= inside )
	{
	  *tri = guess;
	  return REF_SUCCESS;
	}
      
      if ( bary[0] < bary[1] && bary[0] < bary[2] )
	{
	  RSS( ref_update_tri_guess( ref_cell, nodes[1], nodes[2], &guess ),
	       "update 1 2");
	  continue;
	}

      if ( bary[1] < bary[0] && bary[1] < bary[2] )
	{
	  RSS( ref_update_tri_guess( ref_cell, nodes[2], nodes[0], &guess ),
	       "update 2 0");
	  continue;
	}

      if ( bary[2] < bary[0] && bary[2] < bary[1] )
	{
	  RSS( ref_update_tri_guess( ref_cell, nodes[0], nodes[1], &guess ),
	       "update 0 1");
	  continue;
	}

      THROW("unable to find the next step");
    }

  THROW("max steps exceeded");

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

REF_STATUS ref_grid_exhaustive_enclosing_tet( REF_GRID ref_grid, REF_DBL *xyz,
					     REF_INT *tet, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT guess, best_guess;
  REF_DBL best_bary, min_bary;
 
  best_guess = REF_EMPTY;
  best_bary = -999.0;
  each_ref_cell_valid_cell( ref_cell, guess)
    {
      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RXS( ref_node_bary4( ref_node, nodes, xyz, bary ), REF_DIV_ZERO, "bary");
      min_bary = MIN( MIN(bary[0],bary[1]),MIN(bary[2],bary[3]));
      if ( REF_EMPTY == best_guess || min_bary > best_bary )
	{
	  best_guess = guess;
	  best_bary = min_bary;
	}
    }
  
  RUS( REF_EMPTY, best_guess, "failed to find cell");

  *tet = best_guess;
  RSS( ref_cell_nodes( ref_cell, best_guess, nodes), "cell" );
  RSS( ref_node_bary4( ref_node, nodes, xyz, bary ), "bary");
  
  return REF_SUCCESS;
}

REF_STATUS ref_grid_enclosing_tet( REF_GRID ref_grid, REF_DBL *xyz,
				   REF_INT *tet, REF_DBL *bary )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT guess;
  REF_DBL inside;
  REF_INT step, limit;
  
  guess = *tet;
  if ( ! ref_cell_valid(ref_cell,guess) )
    {
      each_ref_cell_valid_cell( ref_cell, guess)
	{
	  break;
	}
    }

  limit = 1000; /* was 10e6^(1/3), required 108 for twod testcase  */
  inside = -1.0e-12; /* initial inside tolerence  */

  for ( step=0; step < limit; step++)
    {
      /* exhastive serach if cell is invalid */
      if ( !ref_cell_valid(ref_cell,guess) )
	{
	  RSS( ref_grid_exhaustive_enclosing_tet( ref_grid, xyz,
						  tet, bary ), "enclose");
	  return REF_SUCCESS;
	}

      RSS( ref_cell_nodes( ref_cell, guess, nodes), "cell" );
      RXS( ref_node_bary4( ref_node, nodes, xyz, bary ), REF_DIV_ZERO, "bary");

      if ( step > 990 )
	{
	  printf("step %d, tet %d, bary %e %e %e %e inside %e\n",
		 step,guess,bary[0],bary[1],bary[2],bary[3],inside);
	}
      
      if ( bary[0] >= inside &&
	   bary[1] >= inside &&
	   bary[2] >= inside &&
	   bary[3] >= inside )
	{
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
  
  /* exhaustive search when walk runs out of iterations */
  RSS( ref_grid_exhaustive_enclosing_tet( ref_grid, xyz,
					  tet, bary ), "enclose");

  return REF_SUCCESS;
}
