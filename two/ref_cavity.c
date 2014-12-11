
#include <stdlib.h>
#include <stdio.h>

#include "ref_cavity.h"

#include "ref_malloc.h"
#include "ref_list.h"

#include "ref_twod.h"

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
        {
          ref_cavity_f2n(ref_cavity,0,face) = REF_EMPTY;
          ref_cavity_f2n(ref_cavity,1,face) = ref_cavity_blank(ref_cavity);
          ref_cavity_blank(ref_cavity) = face;
          ref_cavity_n(ref_cavity)--;
          return REF_SUCCESS;
        }
      else
        {
          return REF_INVALID;
        }
    }

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

  RSS( ref_list_add( ref_cavity_list(ref_cavity), tet ),
       "save tet");

  each_ref_cell_cell_face( ref_cell, cell_face )
  {
    for (node = 0; node<3; node++)
      face_nodes[node] = ref_cell_f2n(ref_cell,node,cell_face,tet);
    RSS( ref_cavity_insert( ref_cavity, face_nodes ), "tet side" );
  }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace_tet( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT node )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face;
  REF_INT cell;

  each_ref_cavity_valid_face( ref_cavity, face )
  {
    nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
    nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
    nodes[2] = ref_cavity_f2n(ref_cavity,2,face);
    nodes[3] = node;
    RSS( ref_cell_add( ref_grid_tet(ref_grid), nodes, &cell ), "add" );
  }

  while ( ref_list_n( ref_cavity_list(ref_cavity) ) > 0 )
    {
      RSS( ref_list_pop( ref_cavity_list(ref_cavity), &cell ), "list" );
      RSS( ref_cell_remove( ref_grid_tet(ref_grid), cell ), "rm" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_add_tri( REF_CAVITY ref_cavity,
                               REF_GRID ref_grid, REF_INT tri )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face[2];

  RSS( ref_list_add( ref_cavity_list(ref_cavity), tri ),
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
       "save tri");

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

REF_STATUS ref_cavity_add_disk( REF_CAVITY ref_cavity,
                                REF_GRID ref_grid, REF_INT node )
{
  REF_INT item, cell;

  each_ref_cell_having_node( ref_grid_tri(ref_grid), node, item, cell )
  {
    RSS( ref_cavity_add_tri( ref_cavity, ref_grid, cell ), "insert");
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
        RSS( ref_node_remove( ref_grid_node(ref_grid),
                              nodes[cell_node] ),"rm");

      RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ),
           "save nodes");
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), tri ), "rm" );
      each_ref_cell_cell_node( ref_grid_tri(ref_grid), cell_node )
      if ( ref_cell_node_empty( ref_grid_tri(ref_grid), nodes[cell_node] ) )
        RSS( ref_node_remove( ref_grid_node(ref_grid),
                              nodes[cell_node] ),"rm");

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

      if ( volume <= 0.0 )
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
            keep_growing = REF_TRUE;
            RSS( ref_cavity_enlarge_face( ref_cavity, ref_grid, face ),"ef" );
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
	THROW("boundary");
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
      face_nodes[0]=ref_cavity_f2n(ref_cavity,0,face);
      face_nodes[1]=ref_cavity_f2n(ref_cavity,1,face);
      face_nodes[2]=ref_cavity_f2n(ref_cavity,2,face);
      face_nodes[3]=face_nodes[0];
      break;
    default:
      THROW("unknown node_per");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_shrink_face( REF_CAVITY ref_cavity,
                                   REF_GRID ref_grid, REF_INT face )
{
  REF_INT ncell, cells[2];
  REF_BOOL have_cell0, have_cell1;


  RSS( ref_cell_list_with2( ref_grid_tri(ref_grid),
			    ref_cavity_f2n(ref_cavity,0,face),
			    ref_cavity_f2n(ref_cavity,1,face),
			    2, &ncell, cells), "more than two" );
  if ( 0 == ncell )
    THROW("cavity triangle missing");
  if ( 1 == ncell )
    THROW("boundary");
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

  return REF_SUCCESS;
}
