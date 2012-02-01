
#include <stdlib.h>
#include <stdio.h>

#include "ref_shard.h"

REF_STATUS ref_shard_create( REF_SHARD *ref_shard_ptr, REF_GRID ref_grid )
{
  REF_SHARD ref_shard;
  REF_INT face;

  (*ref_shard_ptr) = NULL;
  (*ref_shard_ptr) = (REF_SHARD)malloc( sizeof(REF_SHARD_STRUCT) );
  RNS(*ref_shard_ptr,"malloc ref_shard NULL");

  ref_shard = *ref_shard_ptr;

  ref_shard_grid(ref_shard) = ref_grid;

  RSS( ref_face_create( &(ref_shard_face( ref_shard )), 
			ref_shard_grid(ref_shard) ), "create face" );

  ref_shard->mark = (REF_INT *)malloc( ref_face_n(ref_shard_face(ref_shard)) 
					* sizeof(REF_INT));
  RNS(ref_shard->mark,"malloc mark NULL");

  for ( face=0 ; face < ref_face_n(ref_shard_face(ref_shard)) ; face++ )
    ref_shard_mark( ref_shard, face ) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_shard_free( REF_SHARD ref_shard )
{
  if ( NULL == (void *)ref_shard ) return REF_NULL;

  free( ref_shard->mark );
  RSS( ref_face_free( ref_shard_face( ref_shard ) ), "free face" );

  ref_cond_free( ref_shard );

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_to_split( REF_SHARD ref_shard, 
				     REF_INT node0, REF_INT node1 )
{
  REF_INT face;

  RSS( ref_face_spanning( ref_shard_face( ref_shard ), 
			  node0, node1,
			  &face ), "missing face");

  if ( node0 == ref_face_f2n(ref_shard_face( ref_shard ),2,face) ||
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),2,face) )
    {
      if ( 3 == ref_shard_mark(ref_shard,face) )
	RSS( REF_FAILURE, "2-3 mark mismatch");
      ref_shard_mark(ref_shard,face) = 2;
      return REF_SUCCESS;
    }

  if ( node0 == ref_face_f2n(ref_shard_face( ref_shard ),3,face) ||
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),3,face) )
    {
      if ( 2 == ref_shard_mark(ref_shard,face) )
	RSS( REF_FAILURE, "3-2 mark mismatch");
      ref_shard_mark(ref_shard,face) = 3;
      return REF_SUCCESS;
    }

  return REF_FAILURE;
}

REF_STATUS ref_shard_marked( REF_SHARD ref_shard, 
			      REF_INT node0, REF_INT node1,
			      REF_BOOL *marked )
{
  REF_INT face;

  *marked = REF_FALSE;

  RSS( ref_face_spanning( ref_shard_face( ref_shard ), 
			  node0, node1,
			  &face ), "missing face");


  if ( 0 == ref_shard_mark(ref_shard,face) ) return REF_SUCCESS;

  if ( 2 == ref_shard_mark(ref_shard,face) &&
       node0 == ref_face_f2n(ref_shard_face( ref_shard ),2,face) &&
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),0,face) ) 
    *marked = REF_TRUE;
  
  if ( 2 == ref_shard_mark(ref_shard,face) &&
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),2,face) &&
       node0 == ref_face_f2n(ref_shard_face( ref_shard ),0,face) ) 
    *marked = REF_TRUE;

  if ( 3 == ref_shard_mark(ref_shard,face) &&
       node0 == ref_face_f2n(ref_shard_face( ref_shard ),3,face) &&
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),1,face) ) 
    *marked = REF_TRUE;
  
  if ( 3 == ref_shard_mark(ref_shard,face) &&
       node1 == ref_face_f2n(ref_shard_face( ref_shard ),3,face) &&
       node0 == ref_face_f2n(ref_shard_face( ref_shard ),1,face) ) 
    *marked = REF_TRUE;

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_n( REF_SHARD ref_shard, 
			      REF_INT *face_marks, REF_INT *hex_marks )
{
  REF_INT face;
  REF_INT cell, hex_nodes[8];
  REF_BOOL marked16, marked25;

  (*face_marks) = 0;

  for (face=0; face < ref_face_n(ref_shard_face(ref_shard)) ; face++ )
    if ( 0 != ref_shard_mark( ref_shard, face ) ) (*face_marks)++;
  
  (*hex_marks) = 0;
  each_ref_cell_valid_cell_with_nodes( ref_grid_hex(ref_shard_grid(ref_shard)), 
				       cell, hex_nodes)
    {
      
      RSS( ref_shard_marked( ref_shard, hex_nodes[1], hex_nodes[6], 
			      &marked16 ), "1-6");
      RSS( ref_shard_marked( ref_shard, hex_nodes[2], hex_nodes[5], 
			      &marked25 ), "2-5");
      if ( marked16 || marked25 ) (*hex_marks)++;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_cell_edge_split( REF_SHARD ref_shard, 
					    REF_INT cell, REF_INT cell_edge )
{
  REF_CELL ref_cell;
  REF_INT nodes[8];

  ref_cell = ref_grid_hex( ref_shard_grid(ref_shard) );

  RSS( ref_cell_nodes( ref_cell, cell, nodes ), "cell nodes" );

  switch ( cell_edge )
    {
    case 0:  case 11:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[1], nodes[6] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[0], nodes[7] ), "mark2" );
      break;  
    case 5:  case 8:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[2], nodes[5] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[3], nodes[4] ), "mark2" );
      break;

      /*
    case 2:  case 6:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[4], nodes[6] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[0], nodes[2] ), "mark2" );
      break;  
    case 4:  case 7:
      RSS( ref_shard_mark_to_split(ref_shard, nodes[7], nodes[5] ), "mark1" );
      RSS( ref_shard_mark_to_split(ref_shard, nodes[3], nodes[1] ), "mark2" );
      break;
      */

    default:
      printf("%s: %d: %s: cell edge %d not implemented, skipping.\n",
	     __FILE__,__LINE__,__func__,cell_edge);
      /*
      RSB( REF_IMPLEMENT, "can not hadle cell edge",
	   printf("cell edge %d\n",cell_edge););
      */
      break;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_pair( REF_SHARD ref_shard, REF_BOOL *again,
			    REF_INT a0, REF_INT a1, 
			    REF_INT b0, REF_INT b1 )
{
  REF_BOOL a_marked,b_marked;
  
  RSS( ref_shard_marked(ref_shard,a0,a1,&a_marked), "marked? a0-a1" );
  RSS( ref_shard_marked(ref_shard,b0,b1,&b_marked), "marked? b0-b1" );
  if ( a_marked != b_marked )
    {
      *again = REF_TRUE;
      RSS( ref_shard_mark_to_split(ref_shard,a0,a1), "mark a0-a1" );
      RSS( ref_shard_mark_to_split(ref_shard,b0,b1), "mark b0-b1" );
    }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_mark_relax( REF_SHARD ref_shard )
{
  REF_CELL ref_cell;
  REF_BOOL again;
  REF_INT cell;
  REF_BOOL nodes[8];

  ref_cell = ref_grid_hex(ref_shard_grid(ref_shard));

  again = REF_TRUE;

  while (again)
    {

      again = REF_FALSE;

      each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes )
	{
	  RSS( ref_shard_pair( ref_shard, &again, 
				nodes[1], nodes[6],
				nodes[0], nodes[7] ), "not consist");
	  RSS( ref_shard_pair( ref_shard, &again, 
				nodes[5], nodes[2],
				nodes[4], nodes[3] ), "not consist");
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_split( REF_SHARD ref_shard )
{
  REF_GRID ref_grid;
  REF_CELL hex, pri, tri, qua;
  REF_INT cell, hex_nodes[8];
  REF_INT pri_nodes[6], new_cell;
  REF_INT tri_nodes[4], qua_nodes[5];
  REF_BOOL marked;

  ref_grid = ref_shard_grid(ref_shard);
  hex = ref_grid_hex(ref_grid);
  pri = ref_grid_pri(ref_grid);
  tri = ref_grid_tri(ref_grid);
  qua = ref_grid_qua(ref_grid);

  each_ref_cell_valid_cell_with_nodes( qua, cell, qua_nodes)
    {
      RSS( ref_shard_marked( ref_shard, qua_nodes[0], qua_nodes[2], 
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
      RSS( ref_shard_marked( ref_shard, qua_nodes[1], qua_nodes[3], 
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
      RSS( ref_shard_marked( ref_shard, hex_nodes[0], hex_nodes[5], 
			      &marked ), "0-5"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
      RSS( ref_shard_marked( ref_shard, hex_nodes[1], hex_nodes[4], 
			      &marked ), "1-4"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")

      RSS( ref_shard_marked( ref_shard, hex_nodes[1], hex_nodes[6], 
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
	  continue;
	}
      RSS( ref_shard_marked( ref_shard, hex_nodes[2], hex_nodes[5], 
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
	  continue;
	}

      RSS( ref_shard_marked( ref_shard, hex_nodes[0], hex_nodes[2], 
			      &marked ), "0-2"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
      RSS( ref_shard_marked( ref_shard, hex_nodes[1], hex_nodes[3], 
			      &marked ), "1-3"); 
      if ( marked ) RSS( REF_IMPLEMENT, "add split")
    }

  return REF_SUCCESS;
}

REF_STATUS ref_shard_prism_into_tet( REF_GRID ref_grid )
{
  REF_INT cell, temp, new_cell;
  REF_INT pri_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT tet_nodes[REF_CELL_MAX_NODE_PER];
  REF_CELL pri, tet;

  pri = ref_grid_pri(ref_grid);
  tet = ref_grid_tet(ref_grid);
  
  each_ref_cell_valid_cell_with_nodes( pri, cell, pri_nodes )
    {
      RSS( ref_cell_remove( pri, cell ), "remove pri");
      if ( pri_nodes[1] < pri_nodes[2] && pri_nodes[1] < pri_nodes[0] )
	{ /* rotate node 1 into 0 position */
	  temp = pri_nodes[0];
	  pri_nodes[0] = pri_nodes[1];
	  pri_nodes[1] = pri_nodes[2];
	  pri_nodes[2] = temp;
	  temp = pri_nodes[3];
	  pri_nodes[3] = pri_nodes[4];
	  pri_nodes[4] = pri_nodes[5];
	  pri_nodes[5] = temp;
	}
      if ( pri_nodes[2] < pri_nodes[1] && pri_nodes[2] < pri_nodes[0] )
	{ /* rotate node 2 into 0 position */
	  temp = pri_nodes[2];
	  pri_nodes[2] = pri_nodes[1];
	  pri_nodes[1] = pri_nodes[0];
	  pri_nodes[0] = temp;
	  temp = pri_nodes[5];
	  pri_nodes[5] = pri_nodes[4];
	  pri_nodes[4] = pri_nodes[3];
	  pri_nodes[3] = temp;
	}
      /* node 0 is now the smallest index of 0,1,2 */
      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[4];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[1];
      RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");

      if ( pri_nodes[1] < pri_nodes[2] )
	{
	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[5];
	  tet_nodes[3] = pri_nodes[4];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[2];
	  tet_nodes[3] = pri_nodes[5];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	}
      else
	{
	  tet_nodes[0] = pri_nodes[2];
	  tet_nodes[1] = pri_nodes[0];
	  tet_nodes[2] = pri_nodes[4];
	  tet_nodes[3] = pri_nodes[5];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[2];
	  tet_nodes[3] = pri_nodes[4];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	}

    }

  return REF_SUCCESS;
}
