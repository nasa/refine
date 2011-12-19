
#include <stdlib.h>
#include <stdio.h>

#include "ref_hexdiv.h"

REF_STATUS ref_hexdiv_create( REF_HEXDIV *ref_hexdiv_ptr, REF_GRID ref_grid )
{
  REF_HEXDIV ref_hexdiv;
  REF_INT face;

  (*ref_hexdiv_ptr) = NULL;
  (*ref_hexdiv_ptr) = (REF_HEXDIV)malloc( sizeof(REF_HEXDIV_STRUCT) );
  RNS(*ref_hexdiv_ptr,"malloc ref_hexdiv NULL");

  ref_hexdiv = *ref_hexdiv_ptr;

  ref_hexdiv_grid(ref_hexdiv) = ref_grid;

  RSS( ref_face_create( &(ref_hexdiv_face( ref_hexdiv )), 
			ref_hexdiv_grid(ref_hexdiv) ), "create face" );

  ref_hexdiv->mark = (REF_INT *)malloc( ref_face_n(ref_hexdiv_face(ref_hexdiv)) 
					* sizeof(REF_INT));
  RNS(ref_hexdiv->mark,"malloc mark NULL");

  for ( face=0 ; face < ref_face_n(ref_hexdiv_face(ref_hexdiv)) ; face++ )
    ref_hexdiv_mark( ref_hexdiv, face ) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_free( REF_HEXDIV ref_hexdiv )
{
  if ( NULL == (void *)ref_hexdiv ) return REF_NULL;

  free( ref_hexdiv->mark );
  RSS( ref_face_free( ref_hexdiv_face( ref_hexdiv ) ), "free face" );

  ref_cond_free( ref_hexdiv );

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_mark_to_split( REF_HEXDIV ref_hexdiv, 
				     REF_INT node0, REF_INT node1 )
{
  REF_INT face;

  RSS( ref_face_spanning( ref_hexdiv_face( ref_hexdiv ), 
			  node0, node1,
			  &face ), "missing face");

  if ( node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),2,face) ||
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),2,face) )
    {
      if ( 3 == ref_hexdiv_mark(ref_hexdiv,face) )
	RSS( REF_FAILURE, "2-3 mark mismatch");
      ref_hexdiv_mark(ref_hexdiv,face) = 2;
      return REF_SUCCESS;
    }

  if ( node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) ||
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) )
    {
      if ( 2 == ref_hexdiv_mark(ref_hexdiv,face) )
	RSS( REF_FAILURE, "3-2 mark mismatch");
      ref_hexdiv_mark(ref_hexdiv,face) = 3;
      return REF_SUCCESS;
    }

  return REF_FAILURE;
}

REF_STATUS ref_hexdiv_marked( REF_HEXDIV ref_hexdiv, 
			      REF_INT node0, REF_INT node1,
			      REF_BOOL *marked )
{
  REF_INT face;

  *marked = REF_FALSE;

  RSS( ref_face_spanning( ref_hexdiv_face( ref_hexdiv ), 
			  node0, node1,
			  &face ), "missing face");


  if ( 0 == ref_hexdiv_mark(ref_hexdiv,face) ) return REF_SUCCESS;

  if ( 2 == ref_hexdiv_mark(ref_hexdiv,face) &&
       node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),2,face) &&
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),0,face) ) 
    *marked = REF_TRUE;
  
  if ( 2 == ref_hexdiv_mark(ref_hexdiv,face) &&
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),2,face) &&
       node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),0,face) ) 
    *marked = REF_TRUE;

  if ( 3 == ref_hexdiv_mark(ref_hexdiv,face) &&
       node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) &&
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),1,face) ) 
    *marked = REF_TRUE;
  
  if ( 3 == ref_hexdiv_mark(ref_hexdiv,face) &&
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) &&
       node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),1,face) ) 
    *marked = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_mark_n( REF_HEXDIV ref_hexdiv, REF_INT *marks )
{
  REF_INT face;

  (*marks) = 0;

  for (face=0; face < ref_face_n(ref_hexdiv_face(ref_hexdiv)) ; face++ )
    if ( 0 != ref_hexdiv_mark( ref_hexdiv, face ) ) (*marks)++;
  
  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_mark_cell_edge_split( REF_HEXDIV ref_hexdiv, 
					    REF_INT cell, REF_INT cell_edge )
{
  REF_CELL ref_cell;
  REF_INT nodes[8];

  ref_cell = ref_grid_hex( ref_hexdiv_grid(ref_hexdiv) );

  RSS( ref_cell_nodes( ref_cell, cell, nodes ), "" );

  switch ( cell_edge )
    {
    case 0:  case 11:
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv, nodes[1], nodes[6] ), "mark1" );
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv, nodes[0], nodes[7] ), "mark2" );
      break;  
    case 5:  case 8:
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv, nodes[2], nodes[5] ), "mark1" );
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv, nodes[3], nodes[4] ), "mark2" );
      break;
    default:
      return REF_IMPLEMENT;
      break;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_pair( REF_HEXDIV ref_hexdiv, REF_BOOL *again,
			    REF_INT a0, REF_INT a1, 
			    REF_INT b0, REF_INT b1 )
{
  REF_BOOL a_marked,b_marked;
  
  RSS( ref_hexdiv_marked(ref_hexdiv,a0,a1,&a_marked), "marked? a0-a1" );
  RSS( ref_hexdiv_marked(ref_hexdiv,b0,b1,&b_marked), "marked? b0-b1" );
  if ( a_marked != b_marked )
    {
      *again = REF_TRUE;
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv,a0,a1), "mark a0-a1" );
      RSS( ref_hexdiv_mark_to_split(ref_hexdiv,b0,b1), "mark b0-b1" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_mark_relax( REF_HEXDIV ref_hexdiv )
{
  REF_CELL ref_cell;
  REF_BOOL again;
  REF_INT cell;
  REF_BOOL nodes[8];

  ref_cell = ref_grid_hex(ref_hexdiv_grid(ref_hexdiv));

  again = REF_TRUE;

  while (again)
    {

      again = REF_FALSE;

      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  RSS( ref_hexdiv_pair( ref_hexdiv, &again, 
				nodes[1], nodes[6],
				nodes[0], nodes[7] ), "not consist");
	  RSS( ref_hexdiv_pair( ref_hexdiv, &again, 
				nodes[5], nodes[2],
				nodes[4], nodes[3] ), "not consist");
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_split( REF_HEXDIV ref_hexdiv )
{
  REF_GRID ref_grid;
  REF_CELL hex, pri, tri, qua;
  REF_INT cell, hex_nodes[8];
  REF_INT pri_nodes[6], new_cell;
  REF_INT tri_nodes[4], qua_nodes[5];
  REF_BOOL marked;

  ref_grid = ref_hexdiv_grid(ref_hexdiv);
  hex = ref_grid_hex(ref_grid);
  pri = ref_grid_pri(ref_grid);
  tri = ref_grid_tri(ref_grid);
  qua = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes( qua, cell, qua_nodes)
    {
      RSS( ref_hexdiv_marked( ref_hexdiv, qua_nodes[0], qua_nodes[2], 
			      &marked ), "0-2"); 
      if ( marked )
	{
	  RSS( ref_cell_remove( qua, cell ), "remove qua");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[2];
	  tri_nodes[3] = qua_nodes[4]; /* bound id */
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[2];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = qua_nodes[4]; /* bound id */
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	}
      RSS( ref_hexdiv_marked( ref_hexdiv, qua_nodes[1], qua_nodes[3], 
			      &marked ), "1-3"); 
      if ( marked )
	{
	  RSS( ref_cell_remove( qua, cell ), "remove qua");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = qua_nodes[4]; /* bound id */
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  tri_nodes[0] = qua_nodes[1];
	  tri_nodes[1] = qua_nodes[2];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = qua_nodes[4]; /* bound id */
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	}
    }

  each_ref_cell_valid_cell_with_nodes( hex, cell, hex_nodes)
    {
      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[0], hex_nodes[5], 
			      &marked ), "0-5"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[1], hex_nodes[4], 
			      &marked ), "1-4"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")

      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[1], hex_nodes[6], 
			      &marked ), "1-6"); 
      if ( marked )
	{
	  RSS( ref_cell_remove( hex, cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[6];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[7];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "add hex pri 1");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[6];
	  pri_nodes[2] = hex_nodes[2];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[7];
	  pri_nodes[5] = hex_nodes[3];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "add hex_pri 2");
	  break;
	}
      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[2], hex_nodes[5], 
			      &marked ), "2-5"); 
      if ( marked )
	{
	  RSS( ref_cell_remove( hex, cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[2];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[3];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "add hex pri 1");
	  pri_nodes[0] = hex_nodes[2];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[6];
	  pri_nodes[3] = hex_nodes[3];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[7];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "add hex pri 2");
	  break;
	}

      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[0], hex_nodes[2], 
			      &marked ), "0-2"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
      RSS( ref_hexdiv_marked( ref_hexdiv, hex_nodes[1], hex_nodes[3], 
			      &marked ), "1-3"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
    }

  return REF_SUCCESS;
}
