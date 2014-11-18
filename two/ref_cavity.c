
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
  ref_cavity = (*ref_cavity_ptr);

  ref_cavity_n(ref_cavity) = 0;
  ref_cavity_node_per(ref_cavity) = node_per;

  ref_cavity_max(ref_cavity) = 10;

  ref_malloc( ref_cavity->f2n, ref_cavity_max(ref_cavity) *
	      ref_cavity_node_per(ref_cavity), REF_INT);
  for ( face = 0 ; face < ref_cavity_max(ref_cavity) ; face++ ) 
    {
      ref_cavity_f2n(ref_cavity,0,face) = REF_EMPTY;
      ref_cavity_f2n(ref_cavity,1,face) = face+1;
    }
  ref_cavity_f2n(ref_cavity,1,ref_cavity_max(ref_cavity)-1) = REF_EMPTY;
  ref_cavity_blank(ref_cavity) = 0;

  RSS( ref_list_create( &(ref_cavity->ref_list) ), "add list" );

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_free( REF_CAVITY ref_cavity )
{
  if ( NULL == (void *)ref_cavity ) return REF_NULL;
  ref_list_free( ref_cavity->ref_list );
  ref_free( ref_cavity->f2n );
  ref_free( ref_cavity );
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
      chunk = MAX(100,(REF_INT)(1.5*(REF_DBL)orig));
      ref_cavity_max(ref_cavity) = orig + chunk;

      ref_realloc( ref_cavity->f2n, ref_cavity_node_per(ref_cavity) *
		   ref_cavity_max(ref_cavity), REF_INT );

      for (face=orig;face < ref_cavity_max(ref_cavity); face++ ) 
	{
	  ref_cavity_f2n(ref_cavity,0,face)= REF_EMPTY; 
	  ref_cavity_f2n(ref_cavity,1,face) = face+1; 
	}
      ref_cavity_f2n(ref_cavity,1,(ref_cavity->max)-1) = REF_EMPTY; 
      ref_cavity_blank(ref_cavity) = orig;
    }

  face = ref_cavity_blank(ref_cavity);
  ref_cavity_blank(ref_cavity) = ref_cavity_f2n(ref_cavity,1,face);
  for ( node = 0 ; node < ref_cavity_node_per(ref_cavity) ; node++ )
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
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT cell_face, node;
  REF_INT face_nodes[4];

  RSS( ref_list_add( ref_cavity_list(ref_cavity), tet ), 
       "save tet");

  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tet, nodes ), 
       "grab faceid");

  each_ref_cell_cell_face( ref_cell, cell_face )
    {
      for(node=0;node<3;node++)
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
      RSS( ref_list_remove( ref_cavity_list(ref_cavity), &cell ), "list" );
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

  face[0]=nodes[1];face[1]=nodes[2];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 0" ); 
  face[0]=nodes[2];face[1]=nodes[0];
  RSS( ref_cavity_insert( ref_cavity, face ), "side 1" ); 
  face[0]=nodes[0];face[1]=nodes[1];
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

REF_STATUS ref_cavity_tri_pri_tri( REF_GRID ref_grid, REF_INT cell,
				   REF_INT *pri, REF_INT *tri )
{
  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT pri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face[4];
  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, tri_nodes ), 
       "grab tri");
  face[0]=tri_nodes[0];face[1]=tri_nodes[1];face[2]=tri_nodes[2];
  face[3]=face[0];
  RSS( ref_cell_with_face( ref_grid_pri( ref_grid) , face, pri ), "pri" );

  RSS( ref_cell_nodes( ref_grid_pri(ref_grid), *pri, pri_nodes ), 
       "grab pri");
  if ( tri_nodes[0] == pri_nodes[0] ||
       tri_nodes[0] == pri_nodes[1] ||
       tri_nodes[0] == pri_nodes[2] )
    {
      tri_nodes[0] = pri_nodes[3];
      tri_nodes[1] = pri_nodes[5];
      tri_nodes[2] = pri_nodes[4];
    }
  else
    {
      tri_nodes[0] = pri_nodes[0];
      tri_nodes[1] = pri_nodes[1];
      tri_nodes[2] = pri_nodes[2];
    }
    
  RSS( ref_cell_with( ref_grid_tri(ref_grid), tri_nodes, tri ), "tri" );

  return REF_SUCCESS;
}

REF_STATUS ref_cavity_replace_tri( REF_CAVITY ref_cavity, 
				   REF_GRID ref_grid, REF_INT node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT cell, pri, tri;
  REF_INT face;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT faceid0, faceid1;
  REF_INT node2, node3;
  REF_INT clone;
  
  if ( 0 == ref_list_n( ref_cavity_list(ref_cavity) ) )
    return REF_INVALID;

  cell = ref_list_value( ref_cavity_list(ref_cavity), 0 );
  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes ), 
       "grab faceid");
  faceid0 = nodes[3];
  RSS( ref_cavity_tri_pri_tri( ref_grid, cell, &pri, &tri ), "tpt");
  RSS( ref_cell_nodes( ref_grid_tri(ref_grid), tri, nodes ), 
       "grab faceid");
  faceid1 = nodes[3];

  RSS( ref_node_twod_clone( ref_node, node, &clone ), "clone" );

  each_ref_cavity_valid_face( ref_cavity, face )
    {
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
      RSS( ref_list_remove( ref_cavity_list(ref_cavity), &cell ), "list" );
      RSS( ref_cavity_tri_pri_tri( ref_grid, cell, &pri, &tri ), "tpt");
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), cell ), "rm" );
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), tri ), "rm" );
      RSS( ref_cell_remove( ref_grid_pri(ref_grid), pri ), "rm" );
    }

  return REF_SUCCESS;
}

