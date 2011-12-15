
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
      ref_hexdiv_mark(ref_hexdiv,face) = 2;
      return REF_SUCCESS;
    }

  if ( node0 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) ||
       node1 == ref_face_f2n(ref_hexdiv_face( ref_hexdiv ),3,face) )
    {
      ref_hexdiv_mark(ref_hexdiv,face) = 3;
      return REF_SUCCESS;
    }

  return REF_FAILURE;
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

REF_STATUS ref_hexdiv_pair( REF_HEXDIV ref_hexdiv, 
			    REF_CELL ref_cell, REF_INT cell,
			    REF_INT cell_face0, REF_INT mark0, 
			    REF_INT cell_face1, REF_INT mark1 )
{
  REF_INT face0, face1;
  REF_INT nodes[4], node;
  
  for(node=0;node<4;node++)
    nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face0);
  RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), nodes, &face0 ), "face0");
  for(node=0;node<4;node++)
    nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face1);
  RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), nodes, &face1 ), "face1");

  if ( ref_hexdiv_mark( ref_hexdiv, face0 ) == mark0 &&
       ref_hexdiv_mark( ref_hexdiv, face1 ) == mark1 )
    return REF_SUCCESS;

  if ( ref_hexdiv_mark( ref_hexdiv, face0 ) == mark0 )
    {
      if ( ref_hexdiv_mark( ref_hexdiv, face1 ) != 0 ) return REF_FAILURE;
      ref_hexdiv_mark( ref_hexdiv, face1 ) = mark1;  
      return REF_SUCCESS;
    }

  if ( ref_hexdiv_mark( ref_hexdiv, face1 ) == mark1 )
    {
      if ( ref_hexdiv_mark( ref_hexdiv, face0 ) != 0 ) return REF_FAILURE;
      ref_hexdiv_mark( ref_hexdiv, face0 ) = mark0;  
      return REF_SUCCESS;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_mark_relax( REF_HEXDIV ref_hexdiv )
{
  REF_INT group, cell;
  REF_CELL ref_cell;
  REF_BOOL again;

  again = REF_TRUE;

  while (again)
    {

      again = REF_FALSE;

      each_ref_grid_ref_cell( ref_hexdiv_grid(ref_hexdiv), group, ref_cell )
	each_ref_cell_valid_cell( ref_cell, cell )
	  {
	    RSS( ref_hexdiv_pair( ref_hexdiv, ref_cell, cell,
				  1, 2, 3, 2 ), "not consist");
	  }
    }

  return REF_SUCCESS;
}

REF_STATUS ref_hexdiv_split( REF_HEXDIV ref_hexdiv )
{
  REF_GRID ref_grid;
  REF_CELL hex, pri, tri, qua;
  REF_INT cell, hex_nodes[8], face_nodes[4];
  REF_INT node, face;
  REF_INT pri_nodes[6], new_cell;
  REF_INT tri_nodes[4], qua_nodes[5];

  ref_grid = ref_hexdiv_grid(ref_hexdiv);
  hex = ref_grid_hex(ref_grid);
  pri = ref_grid_pri(ref_grid);
  tri = ref_grid_tri(ref_grid);
  qua = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes( qua, cell, qua_nodes)
    {
      RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), qua_nodes, &face ), 
	  "face");
      if ( 2 == ref_hexdiv_mark( ref_hexdiv, face ) )
	{
	  RSS( ref_cell_remove( qua, cell ), "remove qua");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[2];
	  tri_nodes[3] = REF_EMPTY;
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  ref_cell_c2n(tri,3,new_cell) = qua_nodes[4];
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[2];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = REF_EMPTY;
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  ref_cell_c2n(tri,3,new_cell) = qua_nodes[4];
	}
      if ( 3 == ref_hexdiv_mark( ref_hexdiv, face ) )
	{
	  RSS( ref_cell_remove( qua, cell ), "remove qua");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = REF_EMPTY;
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  ref_cell_c2n(tri,3,new_cell) = qua_nodes[4];
	  tri_nodes[0] = qua_nodes[1];
	  tri_nodes[1] = qua_nodes[2];
	  tri_nodes[2] = qua_nodes[3];
	  tri_nodes[3] = REF_EMPTY;
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  ref_cell_c2n(tri,3,new_cell) = qua_nodes[4];
	}
    }

  each_ref_cell_valid_cell_with_nodes( hex, cell, hex_nodes)
    {
      for(node=0;node<4;node++)
	face_nodes[node]=ref_cell_f2n(hex,node,cell,0);
      RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), face_nodes, &face ), 
	  "face 0");
      if ( 2 == ref_hexdiv_mark( ref_hexdiv, face ) )
	return REF_IMPLEMENT;
      if ( 3 == ref_hexdiv_mark( ref_hexdiv, face ) )
	return REF_IMPLEMENT;

      for(node=0;node<4;node++)
	face_nodes[node]=ref_cell_f2n(hex,node,cell,1);
      RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), face_nodes, &face ), 
	  "face 1");
      if ( 2 == ref_hexdiv_mark( ref_hexdiv, face ) )
	{
	  RSS( ref_cell_remove( hex, cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[6];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[7];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[6];
	  pri_nodes[2] = hex_nodes[3];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[7];
	  pri_nodes[5] = hex_nodes[3];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "remove hex");
	  break;
	}
      if ( 3 == ref_hexdiv_mark( ref_hexdiv, face ) )
	{
	  RSS( ref_cell_remove( hex, cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[1];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[2];
	  pri_nodes[3] = hex_nodes[0];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[3];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "remove hex");
	  pri_nodes[0] = hex_nodes[2];
	  pri_nodes[1] = hex_nodes[5];
	  pri_nodes[2] = hex_nodes[6];
	  pri_nodes[3] = hex_nodes[3];
	  pri_nodes[4] = hex_nodes[4];
	  pri_nodes[5] = hex_nodes[7];
	  RSS( ref_cell_add( pri, pri_nodes, &new_cell ), "remove hex");
	  break;
	}

      for(node=0;node<4;node++)
	face_nodes[node]=ref_cell_f2n(hex,node,cell,4);
      RSS(ref_face_with( ref_hexdiv_face(ref_hexdiv), face_nodes, &face ), 
	  "face 4");
      if ( 2 == ref_hexdiv_mark( ref_hexdiv, face ) )
	return REF_IMPLEMENT;
      if ( 3 == ref_hexdiv_mark( ref_hexdiv, face ) )
	return REF_IMPLEMENT;
    }

  return REF_SUCCESS;
}
