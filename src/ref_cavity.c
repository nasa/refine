
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

#include "ref_cavity.h"

#include "ref_malloc.h"
#include "ref_list.h"

#include "ref_twod.h"

#include "ref_edge.h"
#include "ref_adapt.h"
#include "ref_sort.h"

#include "ref_dict.h"

REF_STATUS ref_cavity_create( REF_CAVITY *ref_cavity_ptr, REF_INT node_per )
{
  REF_CAVITY ref_cavity;
  REF_INT face;

  ref_malloc( *ref_cavity_ptr, 1, REF_CAVITY_STRUCT );
  ref_cavity = ( *ref_cavity_ptr );

  ref_cavity_n(ref_cavity) = 0;
  ref_cavity_node_per(ref_cavity) = node_per;

  ref_cavity_max(ref_cavity) = 10;

  ref_malloc_init( ref_cavity->f2n, ref_cavity_max(ref_cavity) *
                   ref_cavity_node_per(ref_cavity), REF_INT, 0);
  for ( face = 0; face < ref_cavity_max(ref_cavity); face++ )
    {
      ref_cavity_f2n(ref_cavity,0,face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity,1,face) = face+1;
    }
  ref_cavity_f2n(ref_cavity,1,ref_cavity_max(ref_cavity)-1) = REF_EMPTY;
  ref_cavity_blank(ref_cavity) = 0;

  RSS( ref_list_create( &( ref_cavity->ref_list ) ), "add list" );

  ref_cavity->debug = REF_FALSE;
  
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_free( REF_CAVITY ref_cavity )
{
  if ( NULL == (void *)ref_cavity )
    return REF_NULL;
  ref_list_free( ref_cavity->ref_list );
  ref_free( ref_cavity->f2n );
  ref_free( ref_cavity );
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_inspect( REF_CAVITY ref_cavity )
{
  REF_INT face, node;
  if ( NULL == (void *)ref_cavity )
    return REF_NULL;
  printf("node_per = %d n = %d max = %d blank = %d\n",
         ref_cavity_node_per( ref_cavity ),
         ref_cavity_n( ref_cavity ),
         ref_cavity_max( ref_cavity ),
         ref_cavity_blank( ref_cavity ));
  for ( face = 0; face < ref_cavity_max(ref_cavity); face++ )
    {
      printf(" f2n[%d] = ",face);
      for ( node = 0; node < ref_cavity_node_per(ref_cavity); node++ )
        printf(" %d ",ref_cavity_f2n(ref_cavity,node,face));
      printf("\n");
    }
  RSS( ref_list_inspect( ref_cavity_list(ref_cavity) ), "insp" );
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_insert( REF_CAVITY ref_cavity, REF_INT *nodes )
{
  REF_INT node, face;
  REF_INT orig, chunk;
  REF_BOOL reversed;

  RXS( ref_cavity_find( ref_cavity, nodes, &face, &reversed),
       REF_NOT_FOUND, "find existing" );

  if ( REF_EMPTY != face )
    {
      if ( reversed )
        { /* two faces with opposite orientation destroy each other */
          ref_cavity_f2n(ref_cavity,0,face) = REF_EMPTY;
          ref_cavity_f2n(ref_cavity,1,face) = ref_cavity_blank(ref_cavity);
          ref_cavity_blank(ref_cavity) = face;
          ref_cavity_n(ref_cavity)--;
          return REF_SUCCESS;
        }
      else
        { /* can't happen, added same face twice */
          return REF_INVALID;
        }
    }

  /* if I need to grow my array of faces */
  if ( REF_EMPTY == ref_cavity_blank(ref_cavity) )
    {
      orig = ref_cavity_max(ref_cavity);
      chunk = MAX(100,(REF_INT)( 1.5*(REF_DBL)orig ));
      ref_cavity_max(ref_cavity) = orig + chunk;

      ref_realloc( ref_cavity->f2n, ref_cavity_node_per(ref_cavity) *
                   ref_cavity_max(ref_cavity), REF_INT );

      for (face = orig; face < ref_cavity_max(ref_cavity); face++ )
        {
          ref_cavity_f2n(ref_cavity,0,face) = REF_EMPTY;
          ref_cavity_f2n(ref_cavity,1,face) = face+1;
        }
      ref_cavity_f2n(ref_cavity,1,( ref_cavity->max )-1) = REF_EMPTY;
      ref_cavity_blank(ref_cavity) = orig;
    }

  face = ref_cavity_blank(ref_cavity);
  ref_cavity_blank(ref_cavity) = ref_cavity_f2n(ref_cavity,1,face);
  for ( node = 0; node < ref_cavity_node_per(ref_cavity); node++ )
    ref_cavity_f2n(ref_cavity,node,face) = nodes[node];

  ref_cavity_n(ref_cavity)++;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_find( REF_CAVITY ref_cavity, REF_INT *nodes,
                            REF_INT *found_face, REF_BOOL *reversed)
{
  REF_INT face;

  *found_face = REF_EMPTY;

  if ( 2 == ref_cavity_node_per( ref_cavity ) )
    {
      each_ref_cavity_valid_face( ref_cavity, face )
      {
        if ( nodes[0] == ref_cavity_f2n(ref_cavity,0,face) &&
             nodes[1] == ref_cavity_f2n(ref_cavity,1,face) )
          {
            *found_face = face;
            *reversed = REF_FALSE;
            return REF_SUCCESS;
          }
        if ( nodes[1] == ref_cavity_f2n(ref_cavity,0,face) &&
             nodes[0] == ref_cavity_f2n(ref_cavity,1,face) )
          {
            *found_face = face;
            *reversed = REF_TRUE;
            return REF_SUCCESS;
          }
      }
    }
  else
    {
      each_ref_cavity_valid_face( ref_cavity, face )
      {
        if ( ( nodes[0] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[1] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[2] == ref_cavity_f2n(ref_cavity,2,face) ) ||
             ( nodes[1] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[2] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[0] == ref_cavity_f2n(ref_cavity,2,face) ) ||
             ( nodes[2] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[0] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[1] == ref_cavity_f2n(ref_cavity,2,face) ) )
          {
            *found_face = face;
            *reversed = REF_FALSE;
            return REF_SUCCESS;
          }
        if ( ( nodes[2] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[1] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[0] == ref_cavity_f2n(ref_cavity,2,face) ) ||
             ( nodes[1] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[0] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[2] == ref_cavity_f2n(ref_cavity,2,face) ) ||
             ( nodes[0] == ref_cavity_f2n(ref_cavity,0,face) &&
               nodes[2] == ref_cavity_f2n(ref_cavity,1,face) &&
               nodes[1] == ref_cavity_f2n(ref_cavity,2,face) ) )
          {
            *found_face = face;
            *reversed = REF_TRUE;
            return REF_SUCCESS;
          }
      }
    }

  return REF_NOT_FOUND;
}

REF_STATUS ref_cavity_add_tet( REF_CAVITY ref_cavity,
                               REF_GRID ref_grid, REF_INT tet )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell_face, node;
  REF_INT face_nodes[4];
  REF_INT already_have_it;

  RSS( ref_list_contains( ref_cavity_list(ref_cavity), tet, &already_have_it ),
       "have tet?");
  if ( already_have_it )
    return REF_SUCCESS;

  RSS( ref_list_push( ref_cavity_list(ref_cavity), tet ),
       "save tet");

  each_ref_cell_cell_face( ref_cell, cell_face )
  {
    for (node = 0; node<3; node++)
      face_nodes[node] = ref_cell_f2n(ref_cell,node,cell_face,tet);
    RSS( ref_cavity_insert( ref_cavity, face_nodes ), "tet side" );
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_rm_tet( REF_CAVITY ref_cavity,
                              REF_GRID ref_grid, REF_INT tet )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT cell_face;
  REF_INT face_nodes[4];

  RSS( ref_list_delete( ref_cavity_list(ref_cavity), tet ),
       "dump tet");

  each_ref_cell_cell_face( ref_cell, cell_face )
  {
    /* reverse face nodes orientation */
    face_nodes[0] = ref_cell_f2n(ref_cell,1,cell_face,tet);
    face_nodes[1] =  ref_cell_f2n(ref_cell,0,cell_face,tet);
    face_nodes[2] =  ref_cell_f2n(ref_cell,2,cell_face,tet);
    RSS( ref_cavity_insert( ref_cavity, face_nodes ), "tet side" );
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace_tet( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT node )
{
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face;
  REF_INT cell;
  REF_INT i;
  REF_DBL volume;

  if ( ref_cavity_debug(ref_cavity) )
    RSS( ref_cavity_tec( ref_cavity, ref_grid, node,
			 "ref_cavity_tet_replace.tec" ),
	 "tec for enlarge_face fail" );
  
  each_ref_cavity_valid_face( ref_cavity, face )
  {
    nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
    nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
    nodes[2] = ref_cavity_f2n(ref_cavity,2,face);
    nodes[3] = node;
    if ( node == nodes[0] || node == nodes[1] || node == nodes[2] )
      continue; /* attached face */
    RSS( ref_cell_add( ref_cell, nodes, &cell ), "add" );
    RSS( ref_node_tet_vol( ref_grid_node(ref_grid), nodes, &volume ), "norm");
    if ( volume <= ref_node_min_volume(ref_grid_node(ref_grid)) )
      printf("%d %d %d %d %e\n",nodes[0],nodes[1],nodes[2],nodes[3],volume);
  }

  while ( ref_list_n( ref_cavity_list(ref_cavity) ) > 0 )
    {
      RSS( ref_list_pop( ref_cavity_list(ref_cavity), &cell ), "list" );
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "rm" );
      RSS( ref_cell_remove( ref_cell, cell ), "rm" );
      for (i=0;i<4;i++)
	if ( ref_adj_empty( ref_cell_adj(ref_cell), nodes[i] ) )
	  RSS( ref_node_remove( ref_grid_node(ref_grid), nodes[i] ), "remove" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_tri( REF_CAVITY ref_cavity,
                               REF_GRID ref_grid, REF_INT tri )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face[2];

  RUS( REF_EMPTY, tri, "index marked empty" ); 

  RSS( ref_list_push( ref_cavity_list(ref_cavity), tri ),
       "save tri");

  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ),
       "grab faceid");

  face[0] = nodes[1]; face[1] = nodes[2];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 0" );
  face[0] = nodes[2]; face[1] = nodes[0];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 1" );
  face[0] = nodes[0]; face[1] = nodes[1];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 2" );

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_rm_tri( REF_CAVITY ref_cavity,
                              REF_GRID ref_grid, REF_INT tri )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face[2];

  RSS( ref_list_delete( ref_cavity_list(ref_cavity), tri ),
       "dump tri");

  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ),
       "grab faceid");

  /* add faces backwards */
  face[1] = nodes[1]; face[0] = nodes[2];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 0" );
  face[1] = nodes[2]; face[0] = nodes[0];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 1" );
  face[1] = nodes[0]; face[0] = nodes[1];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 2" );

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_ball( REF_CAVITY ref_cavity,
                                REF_GRID ref_grid, REF_INT node )
{
  REF_INT item, cell;

  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      each_ref_cell_having_node( ref_grid_tri(ref_grid), node, item, cell )
      {
        RSS( ref_cavity_add_tri( ref_cavity, ref_grid, cell ), "insert");
      }
      break;
    case ( 3 ):
      each_ref_cell_having_node( ref_grid_tet(ref_grid), node, item, cell )
      {
        RSS( ref_cavity_add_tet( ref_cavity, ref_grid, cell ), "insert");
      }
      break;
    default:
      THROW("add_ball unknown node_per");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_edge( REF_CAVITY ref_cavity,
                                REF_GRID ref_grid,
				REF_INT node0, REF_INT node1 )
{
  REF_INT cell, ncell;
  REF_INT cell_to_add[50];
  
  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      RSS( REF_IMPLEMENT, "twod" );
      break;
    case ( 3 ):
      RSS( ref_cell_list_with2(ref_grid_tet(ref_grid),node0,node1,
			       50, &ncell, cell_to_add ), "get list" );
      for ( cell = 0 ; cell < ncell ; cell++ )
	{
	  RSS( ref_cavity_add_tet( ref_cavity, ref_grid,
				   cell_to_add[cell] ), "insert");
	}
      break;
    default:
      THROW("add_ball unknown node_per");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_split_edge( REF_CAVITY ref_cavity,
				  REF_INT node0, REF_INT node1,
				  REF_INT new_node )
{
  REF_INT face;
  REF_INT node, nodes[3];
  
  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      RSS( REF_IMPLEMENT, "twod" );
      break;
    case ( 3 ):
      each_ref_cavity_valid_face( ref_cavity, face )
	{
	  each_ref_cavity_face_node( ref_cavity, node )
	    nodes[node] = ref_cavity_f2n(ref_cavity,node,face);
	  if ( ( node0 == nodes[1] && node1 == nodes[2] ) ||
	       ( node1 == nodes[1] && node0 == nodes[2] ) )
	    {
	      ref_cavity_f2n(ref_cavity,1,face) = new_node;
	      nodes[2] = new_node;
	      RSS( ref_cavity_insert( ref_cavity, nodes ), "insert edge 0" );
	      continue;
	    }
	  if ( ( node0 == nodes[0] && node1 == nodes[2] ) ||
	       ( node1 == nodes[0] && node0 == nodes[2] ) )
	    {
	      ref_cavity_f2n(ref_cavity,0,face) = new_node;
	      nodes[2] = new_node;
	      RSS( ref_cavity_insert( ref_cavity, nodes ), "insert edge 1" );
	      continue;
	    }
	  if ( ( node0 == nodes[0] && node1 == nodes[1] ) ||
	       ( node1 == nodes[0] && node0 == nodes[1] ) )
	    {
	      ref_cavity_f2n(ref_cavity,0,face) = new_node;
	      nodes[1] = new_node;
	      RSS( ref_cavity_insert( ref_cavity, nodes ), "insert edge 2" );
	      continue;
	    }
	}
      break;
    default:
      THROW("add_ball unknown node_per");
    }
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace_tri( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid,
                                   REF_INT node, REF_INT clone )
{
  REF_INT cell, pri, tri;
  REF_INT face;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT faceid0, faceid1;
  REF_INT node2, node3;
  REF_INT cell_node;
  REF_INT existing_cell;

  if ( 0 == ref_list_n( ref_cavity_list(ref_cavity) ) )
    return REF_INVALID;

  cell = ref_list_value( ref_cavity_list(ref_cavity), 0 );
  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes ),
       "grab faceid");
  faceid0 = nodes[3];
  RSS( ref_twod_tri_pri_tri( ref_grid_tri(ref_grid), ref_grid_pri(ref_grid),
                             cell, &pri, &tri ), "tpt");
  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ),
       "grab faceid");
  faceid1 = nodes[3];

  each_ref_cavity_valid_face( ref_cavity, face )
  {

    /* skip a collapsed triangle that in on the boundary of cavity */
    if ( node == ref_cavity_f2n(ref_cavity,0,face) ||
         node == ref_cavity_f2n(ref_cavity,1,face) )
      continue;

    nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
    nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
    nodes[2] = node;
    nodes[3] = faceid0;
    /* skip exisiting triangle */
    RXS( ref_cell_with( ref_grid_tri(ref_grid), nodes, &existing_cell),
         REF_NOT_FOUND, "with failed");
    if ( REF_EMPTY != existing_cell )
      {
        RSS( ref_list_delete( ref_cavity_list(ref_cavity), existing_cell ),
             "existing tri was not marked for removal");
        continue;
      }
    RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &cell ), "add" );

    RSS( ref_twod_opposite_edge(ref_grid_pri(ref_grid),
                                ref_cavity_f2n(ref_cavity,0,face),
                                ref_cavity_f2n(ref_cavity,1,face),
                                &node2,&node3),"opp");
    nodes[0] = node3;
    nodes[1] = node2;
    nodes[2] = clone;
    nodes[3] = faceid1;
    RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &cell ), "add" );

    nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
    nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
    nodes[2] = node;
    nodes[3] = node2;
    nodes[4] = node3;
    nodes[5] = clone;
    RSS( ref_cell_add( ref_grid_pri(ref_grid), nodes, &cell ), "add" );
  }

  while ( ref_list_n( ref_cavity_list(ref_cavity) ) > 0 )
    {
      RSS( ref_list_pop( ref_cavity_list(ref_cavity), &cell ), "list" );
      RSS( ref_twod_tri_pri_tri( ref_grid_tri(ref_grid),
                                 ref_grid_pri(ref_grid),
                                 cell, &pri, &tri ), "tpt");

      RSS( ref_cell_remove( ref_grid_pri(ref_grid), pri ), "rm" );

      RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes ),
           "save nodes");
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), cell ), "rm" );
      each_ref_cell_cell_node( ref_grid_tri(ref_grid), cell_node )
      if ( ref_cell_node_empty( ref_grid_tri(ref_grid), nodes[cell_node] ) )
        {
          RSS( ref_node_remove( ref_grid_node(ref_grid),
                                nodes[cell_node] ),"rm");
          RSS( ref_geom_remove_all( ref_grid_geom(ref_grid),
                                    nodes[cell_node] ),"rm");
        }

      RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ),
           "save nodes");
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), tri ), "rm" );
      each_ref_cell_cell_node( ref_grid_tri(ref_grid), cell_node )
      if ( ref_cell_node_empty( ref_grid_tri(ref_grid), nodes[cell_node] ) )
        {
          RSS( ref_node_remove( ref_grid_node(ref_grid),
                                nodes[cell_node] ),"rm");
          RSS( ref_geom_remove_all( ref_grid_geom(ref_grid),
                                    nodes[cell_node] ),"rm");
        }

    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_visible( REF_CAVITY ref_cavity,
                               REF_NODE ref_node, REF_INT node, REF_INT face,
                               REF_BOOL *visible )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL normal[3];
  REF_DBL volume;

  *visible = REF_FALSE;

  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
      nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
      nodes[2] = node;

      RSS( ref_node_tri_normal( ref_node,nodes,normal ), "norm");

      if ( ( ref_node_xyz(ref_node,1,nodes[0]) > 0.5 &&
             normal[1] >= 0.0 ) ||
           ( ref_node_xyz(ref_node,1,nodes[0]) < 0.5 &&
             normal[1] <= 0.0 ) )
        return REF_SUCCESS;
      break;
    case ( 3 ):
      nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
      nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
      nodes[2] = ref_cavity_f2n(ref_cavity,2,face);
      nodes[3] = node;

      RSS( ref_node_tet_vol( ref_node, nodes, &volume ), "norm");

      if ( volume <= ref_node_min_volume(ref_node) )
        return REF_SUCCESS;

      break;
    default:
      THROW("unknown node_per");
    }

  *visible = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_visible( REF_CAVITY ref_cavity,
                                       REF_GRID ref_grid, REF_INT node )
{
  REF_INT face;
  REF_BOOL visible;
  REF_BOOL keep_growing;
  REF_STATUS status;

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge start %d tets %d faces\n",
	   ref_list_n(ref_cavity_list(ref_cavity)),
	   ref_cavity_n(ref_cavity));
  
  keep_growing = REF_TRUE;
  while (keep_growing)
    {
      keep_growing = REF_FALSE;
      each_ref_cavity_valid_face( ref_cavity, face )
      {
        /* skip a face attached to node */
        if ( node == ref_cavity_f2n(ref_cavity,0,face) ||
             node == ref_cavity_f2n(ref_cavity,1,face) )
          continue;
        if ( 3 == ref_cavity_node_per(ref_cavity) &&
             node == ref_cavity_f2n(ref_cavity,2,face) )
          continue;

        RSS(ref_cavity_visible(ref_cavity, ref_grid_node(ref_grid),
                               node, face, &visible ),"free");
        if ( !visible )
          {
            status =  ref_cavity_enlarge_face( ref_cavity, ref_grid, face );
            RXS( status, REF_INVALID, "enlarge face" );
            if ( REF_SUCCESS == status )
              {
                keep_growing = REF_TRUE;
              }
            else
              {
		if ( ref_cavity_debug(ref_cavity) )
		  RSS( ref_cavity_tec( ref_cavity, ref_grid, node,
				       "enlarge.tec" ),
		       "tec for enlarge_face fail" );
                return status;
              }
          }
      }
    }

  if (ref_cavity_debug(ref_cavity))
    printf(" enlarge final %d tets %d faces\n",
	   ref_list_n(ref_cavity_list(ref_cavity)),
	   ref_cavity_n(ref_cavity));

  if (ref_cavity_debug(ref_cavity))
    RSS( ref_cavity_topo( ref_cavity, ref_grid, node ), "topo");
    
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_shrink_visible( REF_CAVITY ref_cavity,
                                      REF_GRID ref_grid, REF_INT node )
{
  REF_INT face;
  REF_BOOL visible;
  REF_BOOL keep_growing;
  REF_STATUS status;

  keep_growing = REF_TRUE;
  while (keep_growing)
    {
      keep_growing = REF_FALSE;
      each_ref_cavity_valid_face( ref_cavity, face )
      {
        /* skip a face attached to node */
        if ( node == ref_cavity_f2n(ref_cavity,0,face) ||
             node == ref_cavity_f2n(ref_cavity,1,face) )
          continue;
        if ( 3 == ref_cavity_node_per(ref_cavity) &&
             node == ref_cavity_f2n(ref_cavity,2,face) )
          continue;

        RSS(ref_cavity_visible(ref_cavity, ref_grid_node(ref_grid),
                               node, face, &visible ),"free");
        if ( !visible )
          {
            status =  ref_cavity_shrink_face( ref_cavity, ref_grid, face );
            RXS( status, REF_INVALID, "shrink face" );
            if ( REF_SUCCESS == status )
              {
                keep_growing = REF_TRUE;
              }
            else
              {
                RSS( ref_cavity_tec( ref_cavity, ref_grid, node,
                                     "ref_cavity_debug_shrink.tec" ), "tec");
                THROW("boundary, see debug");
              }
          }
      }
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_make_visible( REF_CAVITY ref_cavity,
                                    REF_GRID ref_grid, REF_INT node )
{
  RSB( ref_cavity_shrink_visible(ref_cavity,ref_grid,node),
       "shrink failed",
       ref_cavity_tec(ref_cavity,ref_grid,node,"ref_cavity_fail.tec"));
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_metric( REF_CAVITY ref_cavity,
                                      REF_GRID ref_grid, REF_INT node )
{
  REF_INT face;
  REF_BOOL keep_growing;
  REF_DBL ratio, largest_ratio;
  REF_INT face_node;
  REF_STATUS status;

  keep_growing = REF_TRUE;
  while (keep_growing)
    {
      keep_growing = REF_FALSE;
      each_ref_cavity_valid_face( ref_cavity, face )
      {
        /* skip a face attached to node */
        if ( node == ref_cavity_f2n(ref_cavity,0,face) ||
             node == ref_cavity_f2n(ref_cavity,1,face) )
          continue;
        if ( 3 == ref_cavity_node_per(ref_cavity) &&
             node == ref_cavity_f2n(ref_cavity,2,face) )
          continue;

        largest_ratio = 0.0;
        each_ref_cavity_face_node( ref_cavity, face_node )
        {
          RSS( ref_node_ratio( ref_grid_node(ref_grid),
                               node,
                               ref_cavity_f2n(ref_cavity,face_node,face),
                               &ratio ), "ratio" );
          largest_ratio = MAX( largest_ratio, ratio );
        }
        if ( largest_ratio < 1.0 )
          {
            status = ref_cavity_enlarge_face( ref_cavity, ref_grid, face );
            RXS( status, REF_INVALID, "enlarge face" );
            if ( REF_SUCCESS == status )
              keep_growing = REF_TRUE;
          }
      }
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_enlarge_face( REF_CAVITY ref_cavity,
                                    REF_GRID ref_grid, REF_INT face )
{
  REF_INT ncell, cells[2];
  REF_INT face_nodes[4];
  REF_BOOL have_cell0, have_cell1;
  REF_INT tet0, tet1;

  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      RSS( ref_cell_list_with2( ref_grid_tri(ref_grid),
                                ref_cavity_f2n(ref_cavity,0,face),
                                ref_cavity_f2n(ref_cavity,1,face),
                                2, &ncell, cells), "more than two" );
      if ( 0 == ncell )
        THROW("cavity triangle missing");
      if ( 1 == ncell )
        return REF_INVALID; /* boundary */
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), cells[0],
                              &have_cell0 ), "cell0" );
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), cells[1],
                              &have_cell1 ), "cell1" );
      if ( have_cell0 == have_cell1 )
        THROW("cavity same state");
      if ( have_cell0 )
        RSS( ref_cavity_add_tri( ref_cavity, ref_grid,
                                 cells[1] ), "add c1" );
      if ( have_cell1 )
        RSS( ref_cavity_add_tri( ref_cavity, ref_grid,
                                 cells[0] ), "add c0" );
      break;
    case ( 3 ):
      face_nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
      face_nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
      face_nodes[2] = ref_cavity_f2n(ref_cavity,2,face);
      face_nodes[3] = face_nodes[0];
      RSS( ref_cell_with_face( ref_grid_tet(ref_grid), face_nodes,
                               &tet0, &tet1 ),
           "unable to find tets with face");
      if ( REF_EMPTY == tet0 )
        THROW("cavity tets missing");
      if ( REF_EMPTY == tet1 )
        return REF_INVALID; /* boundary */

      RSS( ref_list_contains( ref_cavity_list(ref_cavity), tet0,
                              &have_cell0 ), "cell0" );
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), tet1,
                              &have_cell1 ), "cell1" );
      if ( have_cell0 == have_cell1 )
        THROW("cavity same state");
      if ( have_cell0 )
        RSS( ref_cavity_add_tet( ref_cavity, ref_grid,
                                 tet1 ), "add c1" );
      if ( have_cell1 )
        RSS( ref_cavity_add_tet( ref_cavity, ref_grid,
                                 tet0 ), "add c0" );
      break;
    default:
      THROW("enlarge unknown node_per");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_shrink_face( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT face )
{
  REF_INT ncell, cells[2];
  REF_INT face_nodes[4];
  REF_BOOL have_cell0, have_cell1;
  REF_INT tet0, tet1;


  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      RSS( ref_cell_list_with2( ref_grid_tri(ref_grid),
                                ref_cavity_f2n(ref_cavity,0,face),
                                ref_cavity_f2n(ref_cavity,1,face),
                                2, &ncell, cells), "more than two" );
      if ( 0 == ncell )
        THROW("cavity triangle missing");
      /* boundary is allowed, use the interior tri */
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), cells[0],
                              &have_cell0 ), "cell0" );
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), cells[1],
                              &have_cell1 ), "cell1" );
      if ( have_cell0 == have_cell1 )
        THROW("cavity same state");
      if ( !have_cell0 )
        RSS( ref_cavity_rm_tri( ref_cavity, ref_grid,
                                cells[1] ), "rm c1" );
      if ( !have_cell1 )
        RSS( ref_cavity_rm_tri( ref_cavity, ref_grid,
                                cells[0] ), "rm c0" );
      break;
    case ( 3 ):
      face_nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
      face_nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
      face_nodes[2] = ref_cavity_f2n(ref_cavity,2,face);
      face_nodes[3] = face_nodes[0];
      RSS( ref_cell_with_face( ref_grid_tet(ref_grid), face_nodes,
                               &tet0, &tet1 ),
           "unable to find tets with face");
      if ( REF_EMPTY == tet0 )
        THROW("cavity tets missing");
      /* boundary is allowed, use the interior tet */

      RSS( ref_list_contains( ref_cavity_list(ref_cavity), tet0,
                              &have_cell0 ), "cell0" );
      RSS( ref_list_contains( ref_cavity_list(ref_cavity), tet1,
                              &have_cell1 ), "cell1" );
      if ( have_cell0 == have_cell1 )
        THROW("cavity same state");
      if ( !have_cell0 )
        RSS( ref_cavity_rm_tet( ref_cavity, ref_grid,
                                tet1 ), "add c1" );
      if ( !have_cell1 )
        RSS( ref_cavity_rm_tet( ref_cavity, ref_grid,
                                tet0 ), "add c0" );
      break;
    default:
      THROW("shrink unknown node_per");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_tet_quality( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tet(ref_grid);
  REF_CAVITY ref_cavity;
  REF_INT node, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL min_quality = 0.10;
  REF_DBL quality;
  REF_INT count;
  char filename[1024];
  REF_BOOL improved;
  REF_BOOL debug = REF_FALSE;
  count = 0;
  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
  {
    node = nodes[0];
    RSS( ref_node_tet_quality( ref_node, nodes, &quality ), "qual");
    if ( quality < min_quality )
      {
        snprintf( filename, 1024, "cavity%04d.tec", count ); count++;
        RSS(ref_cavity_create(&ref_cavity,3),"create");
        RSS(ref_cavity_add_tet(ref_cavity,ref_grid,cell),"insert first");
        if (debug)
          RSS(ref_cavity_tec(ref_cavity, ref_grid, node, filename ),"tec");

        RXS(ref_cavity_enlarge_face(ref_cavity,ref_grid,0),
            REF_INVALID,"remove problem");
        if (debug)
          RSS(ref_cavity_tec(ref_cavity, ref_grid, node, filename ),"tec");

        RXS(ref_cavity_enlarge_visible(ref_cavity,ref_grid,node),
            REF_INVALID,"enlarge viz");
        RSS(ref_cavity_change(ref_cavity, ref_grid, node,
                              &improved), "change" );
        if (debug)
          RSS(ref_cavity_tec(ref_cavity, ref_grid, node, filename ),"tec");
        if (REF_FALSE && improved)
          RSS(ref_cavity_replace_tet(ref_cavity, ref_grid, node ),"replace");
        RSS(ref_cavity_free(ref_cavity),"free");
      }
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_twod_pass( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_EDGE ref_edge;
  REF_DBL *ratio;
  REF_INT *order;
  REF_INT ntarget, *target, *node2target;
  REF_INT node, node0, node1, opp;
  REF_INT i, edge;
  REF_DBL edge_ratio;
  REF_BOOL active;
  REF_CAVITY ref_cavity;
  REF_BOOL improved;

  RSS( ref_edge_create( &ref_edge, ref_grid ), "orig edges" );

  ref_malloc_init( ratio, ref_node_max(ref_node),
                   REF_DBL, 2.0*ref_grid_adapt(ref_grid,collapse_ratio) );

  for (edge = 0; edge<ref_edge_n(ref_edge); edge++)
    {
      node0 = ref_edge_e2n( ref_edge, 0, edge );
      node1 = ref_edge_e2n( ref_edge, 1, edge );
      RSS( ref_node_edge_twod( ref_node, node0, node1,
                               &active ), "act" );
      if ( !active )
        continue;

      RSS( ref_node_ratio( ref_node, node0, node1,
                           &edge_ratio ), "ratio");
      ratio[node0] = MIN( ratio[node0], edge_ratio );
      ratio[node1] = MIN( ratio[node1], edge_ratio );
    }

  ref_malloc( target, ref_node_n(ref_node), REF_INT );
  ref_malloc_init( node2target, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  ntarget = 0;
  for ( node = 0; node < ref_node_max(ref_node); node++ )
    if ( ratio[node] < ref_grid_adapt(ref_grid,collapse_ratio) )
      {
        node2target[node] = ntarget;
        target[ntarget] = node;
        ratio[ntarget] = ratio[node];
        ntarget++;
      }

  ref_malloc( order, ntarget, REF_INT );

  RSS( ref_sort_heap_dbl( ntarget, ratio, order), "sort lengths" );

  for ( i = 0; i < ntarget; i++ )
    {
      if ( ratio[order[i]] > ref_grid_adapt(ref_grid,collapse_ratio) )
        continue;
      node = target[order[i]];
      if ( ref_node_valid(ref_node,node) )
        {
          RSS(ref_cavity_create(&ref_cavity,2),"create");
          RSS(ref_cavity_add_ball(ref_cavity,ref_grid,node),"insert first");
          RSS(ref_cavity_enlarge_metric(ref_cavity,ref_grid,node),"enlarge short");
          RSS(ref_cavity_make_visible(ref_cavity,ref_grid,node),"make valid");
          RSS(ref_twod_opposite_node(ref_grid_pri(ref_grid), node, &opp), "opp");
          RSS(ref_cavity_change(ref_cavity, ref_grid, node,
                                &improved), "change" );
          if ( improved )
            RSS(ref_cavity_replace_tri(ref_cavity, ref_grid,
                                       node, opp ),"free");
          RSS(ref_cavity_free(ref_cavity),"free");
        }
    }

  ref_free( order );
  ref_free( node2target );
  ref_free( target );
  ref_free( ratio );

  ref_edge_free( ref_edge );

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_tec( REF_CAVITY ref_cavity, REF_GRID ref_grid,
                           REF_INT node, const char *filename )
{
  REF_DICT node_dict, face_dict;
  REF_CELL ref_cell;
  REF_INT face, face_node;
  REF_INT cell, cell_node, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, local;
  REF_DBL xyz_phys[3];
  char *zonetype;
  FILE *f;

  switch ( ref_cavity_node_per( ref_cavity ) )
    {
    case ( 2 ):
      ref_cell = ref_grid_tri(ref_grid);
      zonetype = "fetriangle";
      break;
    case ( 3 ):
      ref_cell = ref_grid_tet(ref_grid);
      zonetype = "fetetrahedron";
      break;
    default:
      THROW("tec unknown node_per");
    }

  f = fopen(filename,"w");
  if (NULL == (void *)f)
    printf("unable to open %s\n",filename);
  RNS(f, "unable to open file" );

  fprintf(f, "title=\"tecplot refine cavity\"\n");
  fprintf(f, "variables = \"x\" \"y\" \"z\"\n");

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&face_dict),"create faces");

  each_ref_list_item( ref_cavity_list(ref_cavity), item )
  {
    cell = ref_list_value( ref_cavity_list(ref_cavity), item );
    RSS( ref_dict_store( face_dict, cell, 0 ), "store");
    RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
    each_ref_cell_cell_node( ref_cell, cell_node )
    RSS( ref_dict_store( node_dict, nodes[cell_node], 0 ), "store");
  }

  fprintf(f,
          "zone t=\"old\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
          ref_dict_n(node_dict), ref_dict_n(face_dict),
          "point", zonetype );
  for ( item = 0; item < ref_dict_n(node_dict); item++ )
    {
      local = ref_dict_key(node_dict,item);
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid),0,local);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid),1,local);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid),2,local);
      fprintf(f, " %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
    }

  for ( item = 0; item < ref_dict_n(face_dict); item++ )
    {
      cell = ref_dict_key(face_dict,item);
      RSS( ref_cell_nodes( ref_cell, cell, nodes), "nodes");
      each_ref_cell_cell_node( ref_cell, cell_node )
      {
        RSS( ref_dict_location( node_dict,
                                nodes[cell_node],
                                &local), "ret");
        fprintf(f," %d",local + 1);
      }
      fprintf(f,"\n");
    }

  RSS(ref_dict_free(face_dict),"free tris");
  RSS(ref_dict_free(node_dict),"free nodes");

  RSS(ref_dict_create(&node_dict),"create nodes");
  RSS(ref_dict_create(&face_dict),"create faces");

  RSS( ref_dict_store( node_dict, node, 0 ), "store");
  each_ref_cavity_valid_face( ref_cavity, face )
  {
    RSS( ref_dict_store( face_dict, face, 0 ), "store");
    each_ref_cavity_face_node( ref_cavity, face_node )
    RSS( ref_dict_store( node_dict,
                         ref_cavity_f2n(ref_cavity,face_node,face), 0 ),
         "store");
  }

  fprintf(f,
          "zone t=\"new\", nodes=%d, elements=%d, datapacking=%s, zonetype=%s\n",
          ref_dict_n(node_dict), ref_dict_n(face_dict),
          "point", zonetype );
  for ( item = 0; item < ref_dict_n(node_dict); item++ )
    {
      local = ref_dict_key(node_dict,item);
      xyz_phys[0] = ref_node_xyz(ref_grid_node(ref_grid),0,local);
      xyz_phys[1] = ref_node_xyz(ref_grid_node(ref_grid),1,local);
      xyz_phys[2] = ref_node_xyz(ref_grid_node(ref_grid),2,local);
      fprintf(f, " %.16e %.16e %.16e\n", xyz_phys[0], xyz_phys[1], xyz_phys[2]);
    }

  for ( item = 0; item < ref_dict_n(face_dict); item++ )
    {
      face = ref_dict_key(face_dict,item);
      RSS( ref_dict_location( node_dict, node, &local), "center node");
      fprintf(f," %d",local + 1);
      each_ref_cavity_face_node( ref_cavity, face_node )
      {
        RSS( ref_dict_location( node_dict,
                                ref_cavity_f2n(ref_cavity,face_node,face),
                                &local), "ret");
        fprintf(f," %d",local + 1);
      }
      fprintf(f,"\n");
    }

  RSS(ref_dict_free(face_dict),"free face");
  RSS(ref_dict_free(node_dict),"free nodes");

  fclose(f);

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_local( REF_CAVITY ref_cavity, REF_GRID ref_grid,
                              REF_INT node, REF_BOOL *local )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face, face_node;
  REF_INT rank = ref_mpi_rank(ref_grid_mpi(ref_grid));
  
  REIS(3,ref_cavity_node_per( ref_cavity ), "only implemented for tets" );

  *local = REF_FALSE;
  if ( rank != ref_node_part( ref_node, node ) )
    return REF_SUCCESS;

  each_ref_list_item( ref_cavity_list(ref_cavity), item )
    {
      cell = ref_list_value( ref_cavity_list(ref_cavity), item );
      RSS( ref_cell_nodes( ref_grid_tet(ref_grid), cell, nodes), "cell");
      if ( rank != ref_node_part( ref_node, nodes[0] ) ||
	   rank != ref_node_part( ref_node, nodes[1] ) ||
	   rank != ref_node_part( ref_node, nodes[2] ) ||
	   rank != ref_node_part( ref_node, nodes[3] ) )
	{
	  *local = REF_FALSE;
	  return REF_SUCCESS;
	}
    }

  each_ref_cavity_valid_face( ref_cavity, face )
    {
      each_ref_cavity_face_node( ref_cavity, face_node )
	nodes[face_node] = ref_cavity_f2n(ref_cavity,face_node,face);
      nodes[ref_cavity_node_per(ref_cavity)] = node;
      if ( rank != ref_node_part( ref_node, nodes[0] ) ||
	   rank != ref_node_part( ref_node, nodes[1] ) ||
	   rank != ref_node_part( ref_node, nodes[2] ) ||
	   rank != ref_node_part( ref_node, nodes[3] ) )
	{
	  *local = REF_FALSE;
	  return REF_SUCCESS;
	}
  }

  *local = REF_TRUE;
  return REF_SUCCESS;
}

REF_STATUS ref_cavity_change( REF_CAVITY ref_cavity, REF_GRID ref_grid,
                              REF_INT node, REF_BOOL *improved )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT item, cell;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL quality, min_quality, total_quality;
  REF_INT face, face_node, n, n_del, n_add;
  REF_BOOL skip;
  REF_DBL min_del, min_add;
  *improved = REF_FALSE;
  n = 0;
  min_quality = 1.0;
  total_quality = 0.0;
  each_ref_list_item( ref_cavity_list(ref_cavity), item )
  {
    cell = ref_list_value( ref_cavity_list(ref_cavity), item );
    switch ( ref_cavity_node_per( ref_cavity ) )
      {
      case ( 2 ):
        RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes), "cell");
        RSS( ref_node_tri_quality( ref_node, nodes, &quality ), "new qual");
        break;
      case ( 3 ):
        RSS( ref_cell_nodes( ref_grid_tet(ref_grid), cell, nodes), "cell");
        RSS( ref_node_tet_quality( ref_node, nodes, &quality ), "new qual");
        break;
      default:
        THROW("change unknown node_per");
      }
    n++;
    total_quality += quality;
    min_quality = MIN( min_quality, quality );
  }
  if ( REF_FALSE && n > 0 )
    printf("- min %12.8f avg %12.8f n %d\n",
           min_quality, total_quality/((REF_DBL)n ), n);
  min_del = min_quality;
  n_del = n;

  n = 0;
  min_quality = 1.0;
  total_quality = 0.0;
  each_ref_cavity_valid_face( ref_cavity, face )
  {
    skip = REF_FALSE;
    /* skip a collapsed triangle that in on the boundary of cavity */
    each_ref_cavity_face_node( ref_cavity, face_node )
    if ( node == ref_cavity_f2n(ref_cavity,face_node,face) )
      skip = REF_TRUE;
    if (skip)
      continue;
    each_ref_cavity_face_node( ref_cavity, face_node )
    nodes[face_node] = ref_cavity_f2n(ref_cavity,face_node,face);
    nodes[ref_cavity_node_per(ref_cavity)] = node;
    switch ( ref_cavity_node_per( ref_cavity ) )
      {
      case ( 2 ):
        RSS( ref_node_tri_quality( ref_node, nodes, &quality ), "new qual");
        break;
      case ( 3 ):
        RSS( ref_node_tet_quality( ref_node, nodes, &quality ), "new qual");
        break;
      default:
        THROW("change unknown node_per");
      }
    n++;
    total_quality += quality;
    min_quality = MIN( min_quality, quality );
  }
  if ( REF_FALSE && n > 0 )
    printf("+ min %12.8f avg %12.8f n %d\n",
           min_quality, total_quality/((REF_DBL)n ), n);
  min_add = min_quality;
  n_add = n;

  printf(" min %12.8f <- %12.8f diff %12.8f n %d <- %d\n",
         min_add, min_del, min_add-min_del, n_add, n_del);

  if ( min_add > min_del )
    *improved = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_topo( REF_CAVITY ref_cavity, REF_GRID ref_grid,
			    REF_INT node )
{
  REF_INT item, cell, face, face_node;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
 
  each_ref_list_item( ref_cavity_list(ref_cavity), item )
    {
      cell = ref_list_value( ref_cavity_list(ref_cavity), item );
      switch ( ref_cavity_node_per( ref_cavity ) )
	{
	case ( 2 ):
	  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes), "cell");
	  break;
	case ( 3 ):
	  RSS( ref_cell_nodes( ref_grid_tet(ref_grid), cell, nodes), "cell");
	  break;
	default:
	  THROW("change unknown node_per");
	}
      printf("old %d %d %d %d\n",nodes[0],nodes[1],nodes[2],nodes[3]);
    }

  each_ref_cavity_valid_face( ref_cavity, face )
    {
      if ( node == ref_cavity_f2n(ref_cavity,0,face) )
	continue;
      if ( node == ref_cavity_f2n(ref_cavity,1,face) )
	continue;
      if ( 2 < ref_cavity_node_per( ref_cavity ) &&
	   node == ref_cavity_f2n(ref_cavity,2,face) )
	continue;
      printf("new ");
      for ( face_node = 0;
	    face_node < ref_cavity_node_per(ref_cavity);
	    face_node++ )
        printf(" %d ",ref_cavity_f2n(ref_cavity,face_node,face));
      printf(" %d ",node);
      printf("\n");
  }

  return REF_SUCCESS;
}
