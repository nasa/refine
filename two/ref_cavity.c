
#include <stdlib.h>
#include <stdio.h>

#include "ref_cavity.h"

#include "ref_malloc.h"
#include "ref_list.h"

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

  return REF_NOT_FOUND;
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

REF_STATUS ref_cavity_replace( REF_CAVITY ref_cavity, 
			       REF_GRID ref_grid, REF_INT node )
{
  REF_INT cell;
  REF_INT face;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  
  nodes[3] = REF_EMPTY; /* faceid */
  while ( ref_list_n( ref_cavity_list(ref_cavity) ) > 0 )
    {
      RSS( ref_list_remove( ref_cavity_list(ref_cavity), &cell ), "list" );
      RSS( ref_cell_nodes( ref_grid_tri(ref_grid), cell, nodes ), 
	   "grab faceid");
      RSS( ref_cell_remove( ref_grid_tri(ref_grid), cell ), "rm" );
    }

  each_ref_cavity_valid_face( ref_cavity, face )
    {
      nodes[0] = ref_cavity_f2n(ref_cavity,0,face);
      nodes[1] = ref_cavity_f2n(ref_cavity,1,face);
      nodes[2] = node;
      RSS( ref_cell_add( ref_grid_tri(ref_grid), nodes, &cell ), "add" );
    }

  return REF_SUCCESS;
}

