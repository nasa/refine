
#include <stdlib.h>
#include <stdio.h>

#include "ref_shard.h"
#include "ref_quality.h"
#include "ref_export.h"

#include "ref_malloc.h"

REF_STATUS ref_shard_create( REF_SHARD *ref_shard_ptr, REF_GRID ref_grid )
{
  REF_SHARD ref_shard;
  REF_INT face;

  ref_malloc( *ref_shard_ptr, 1, REF_SHARD_STRUCT );

  ref_shard = *ref_shard_ptr;

  ref_shard_grid(ref_shard) = ref_grid;

  RSS( ref_face_create( &(ref_shard_face( ref_shard )), 
			ref_shard_grid(ref_shard) ), "create face" );

  ref_malloc(ref_shard->mark, ref_face_n(ref_shard_face(ref_shard)), REF_INT );

  for ( face=0 ; face < ref_face_n(ref_shard_face(ref_shard)) ; face++ )
    ref_shard_mark( ref_shard, face ) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_shard_free( REF_SHARD ref_shard )
{
  if ( NULL == (void *)ref_shard ) return REF_NULL;

  ref_free( ref_shard->mark );
  RSS( ref_face_free( ref_shard_face( ref_shard ) ), "free face" );

  ref_free( ref_shard );

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
  REF_INT cell, hex_nodes[REF_CELL_MAX_SIZE_PER];
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
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];

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

static REF_STATUS ref_shard_pair( REF_SHARD ref_shard, REF_BOOL *again,
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
  REF_BOOL nodes[REF_CELL_MAX_SIZE_PER];

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
  REF_INT cell, hex_nodes[REF_CELL_MAX_SIZE_PER];
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

#define check_tet_volume()						\
  {									\
    REF_DBL vol;							\
    RSS( ref_node_tet_vol( ref_node, tet_nodes, &vol ), "tet vol");	\
    if( vol<=0.0 )							\
      {									\
	REF_GRID viz;							\
	REF_INT newnew;							\
	printf("tet vol %e\n",vol);					\
	printf("minnode %d\n",minnode);					\
	printf("orig %d %d %d %d %d %d\n",				\
	       orig[0],orig[1],orig[2],					\
	       orig[3],orig[4],orig[5]);				\
	printf("prism %d %d %d %d %d %d\n",				\
	       pri_nodes[0],pri_nodes[1],pri_nodes[2],			\
	       pri_nodes[3],pri_nodes[4],pri_nodes[5]);			\
	printf("tet %d %d %d %d\n",					\
	       tet_nodes[0],tet_nodes[1],tet_nodes[2],			\
	       tet_nodes[3]);						\
	RSS( ref_grid_empty_cell_clone(&viz,ref_grid),"viz");		\
	RSS( ref_cell_add( ref_grid_pri(viz), orig, &newnew ), "o");	\
	RSS( ref_cell_add( ref_grid_pri(viz), pri_nodes, &newnew ), "p"); \
	RSS( ref_cell_add( ref_grid_tet(viz), tet_nodes, &newnew ), "t"); \
	RSS( ref_grid_free_cell_clone(viz),"free temp grid");		\
      }									\
  }

/*
	RSS(ref_export_by_extension(viz, "neg.tec"),"to tec");		\
	RSS( REF_FAILURE, "neg vol tet");                               \
 
*/

REF_STATUS ref_shard_prism_into_tet( REF_GRID ref_grid, 
				     REF_INT keeping_n_layers, 
				     REF_INT of_faceid )
{
  REF_INT cell, new_cell, minnode;

  REF_INT orig[REF_CELL_MAX_SIZE_PER];
  REF_INT pri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT tet_nodes[REF_CELL_MAX_SIZE_PER];
  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_CELL pyr = ref_grid_pyr(ref_grid);
  REF_CELL tet = ref_grid_tet(ref_grid);

  REF_INT tri_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT qua_nodes[REF_CELL_MAX_SIZE_PER];
  REF_CELL qua = ref_grid_qua(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);

  REF_NODE ref_node = ref_grid_node(ref_grid);

  REF_INT relaxation;
  REF_INT node;
  REF_INT *mark, *mark_copy;

  ref_malloc_init( mark,      ref_node_max(ref_node), REF_INT, REF_EMPTY );
  ref_malloc_init( mark_copy, ref_node_max(ref_node), REF_INT, REF_EMPTY );

  /* mark nodes on prism tris */
  each_ref_cell_valid_cell_with_nodes( pri, cell, orig )
    {
      tri_nodes[0] = orig[0];
      tri_nodes[1] = orig[1];
      tri_nodes[2] = orig[2];
      RXS( ref_cell_with( tri, tri_nodes, &new_cell ), REF_NOT_FOUND, "with");
      if ( REF_EMPTY != new_cell && of_faceid == ref_cell_c2n(tri,3,new_cell) )
	{ 
	  mark[ tri_nodes[0] ] = 0;
	  mark[ tri_nodes[1] ] = 0;
	  mark[ tri_nodes[2] ] = 0;
	}
      tri_nodes[0] = orig[3];
      tri_nodes[1] = orig[5];
      tri_nodes[2] = orig[4];
      RXS( ref_cell_with( tri, tri_nodes, &new_cell ), REF_NOT_FOUND, "with");
      if ( REF_EMPTY != new_cell && of_faceid == ref_cell_c2n(tri,3,new_cell) )
	{ 
	  mark[ tri_nodes[0] ] = 0;
	  mark[ tri_nodes[1] ] = 0;
	  mark[ tri_nodes[2] ] = 0;
	}
    }

  for (relaxation=0;relaxation<keeping_n_layers;relaxation++)
    {
      for (node=0;node<ref_node_max(ref_node);node++)
	mark_copy[node]=mark[node];
      each_ref_cell_valid_cell_with_nodes( pri, cell, orig )
	{
	  if ( mark_copy[orig[0]] == REF_EMPTY &&
	       mark_copy[orig[3]] != REF_EMPTY )
	    mark[orig[0]] = mark_copy[orig[3]] + 1;
	  if ( mark_copy[orig[3]] == REF_EMPTY &&
	       mark_copy[orig[0]] != REF_EMPTY )
	    mark[orig[3]] = mark_copy[orig[0]] + 1;

	  if ( mark_copy[orig[1]] == REF_EMPTY &&
	       mark_copy[orig[4]] != REF_EMPTY )
	    mark[orig[1]] = mark_copy[orig[4]] + 1;
	  if ( mark_copy[orig[4]] == REF_EMPTY &&
	       mark_copy[orig[1]] != REF_EMPTY )
	    mark[orig[4]] = mark_copy[orig[1]] + 1;

	  if ( mark_copy[orig[2]] == REF_EMPTY &&
	       mark_copy[orig[5]] != REF_EMPTY )
	    mark[orig[2]] = mark_copy[orig[5]] + 1;
	  if ( mark_copy[orig[5]] == REF_EMPTY &&
	       mark_copy[orig[2]] != REF_EMPTY )
	    mark[orig[5]] = mark_copy[orig[2]] + 1;
	}
    }

  each_ref_cell_valid_cell_with_nodes( pri, cell, orig )
    {

      if ( mark[orig[0]] != REF_EMPTY &&
	   mark[orig[3]] != REF_EMPTY ) continue;

      RSS( ref_cell_remove( pri, cell ), "remove pri");

      minnode = MIN( MIN( orig[0], orig[1] ), MIN( orig[2], orig[3] ) );
      minnode = MIN( MIN( orig[4], orig[5] ), minnode );

      pri_nodes[0] = orig[0];
      pri_nodes[1] = orig[1];
      pri_nodes[2] = orig[2];
      pri_nodes[3] = orig[3];
      pri_nodes[4] = orig[4];
      pri_nodes[5] = orig[5];
 
      if ( orig[1] == minnode )
	{
	  pri_nodes[0] = orig[1];
	  pri_nodes[1] = orig[2];
	  pri_nodes[2] = orig[0];
	  pri_nodes[3] = orig[4];
	  pri_nodes[4] = orig[5];
	  pri_nodes[5] = orig[3];
	}
 
      if ( orig[2] == minnode )
	{
	  pri_nodes[0] = orig[2];
	  pri_nodes[1] = orig[0];
	  pri_nodes[2] = orig[1];
	  pri_nodes[3] = orig[5];
	  pri_nodes[4] = orig[3];
	  pri_nodes[5] = orig[4];
	}
 
      if ( orig[3] == minnode )
	{
	  pri_nodes[0] = orig[3];
	  pri_nodes[1] = orig[5];
	  pri_nodes[2] = orig[4];
	  pri_nodes[3] = orig[0];
	  pri_nodes[4] = orig[2];
	  pri_nodes[5] = orig[1];
	}
 
      if ( orig[4] == minnode )
	{
	  pri_nodes[0] = orig[4];
	  pri_nodes[1] = orig[3];
	  pri_nodes[2] = orig[5];
	  pri_nodes[3] = orig[1];
	  pri_nodes[4] = orig[0];
	  pri_nodes[5] = orig[2];
	}
 
      if ( orig[5] == minnode )
	{
	  pri_nodes[0] = orig[5];
	  pri_nodes[1] = orig[4];
	  pri_nodes[2] = orig[3];
	  pri_nodes[3] = orig[2];
	  pri_nodes[4] = orig[1];
	  pri_nodes[5] = orig[0];
	}
 
      /* node 0 is now the smallest index of prism */
 
      tet_nodes[0] = pri_nodes[0];
      tet_nodes[1] = pri_nodes[4];
      tet_nodes[2] = pri_nodes[5];
      tet_nodes[3] = pri_nodes[3];
      RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
      check_tet_volume( );

      if ( ( pri_nodes[1] < pri_nodes[2] && pri_nodes[1] < pri_nodes[4] ) ||
	   ( pri_nodes[5] < pri_nodes[2] && pri_nodes[5] < pri_nodes[4] ) )
	{
	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[5];
	  tet_nodes[3] = pri_nodes[4];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  check_tet_volume( );

	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[2];
	  tet_nodes[3] = pri_nodes[5];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  check_tet_volume(  );
	}
      else
	{
	  tet_nodes[0] = pri_nodes[2];
	  tet_nodes[1] = pri_nodes[0];
	  tet_nodes[2] = pri_nodes[4];
	  tet_nodes[3] = pri_nodes[5];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  check_tet_volume(  );

	  tet_nodes[0] = pri_nodes[0];
	  tet_nodes[1] = pri_nodes[1];
	  tet_nodes[2] = pri_nodes[2];
	  tet_nodes[3] = pri_nodes[4];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  check_tet_volume(  );
	}

    }

  each_ref_cell_valid_cell_with_nodes( pyr, cell, orig )
    {
      if ( mark[orig[0]] != REF_EMPTY &&
	   mark[orig[1]] != REF_EMPTY &&
	   mark[orig[3]] != REF_EMPTY &&
	   mark[orig[4]] != REF_EMPTY ) continue;

      RSS( ref_cell_remove( pyr, cell ), "remove qua");
      if ( ( orig[0] < orig[1] && orig[0] < orig[3] ) ||
	   ( orig[4] < orig[1] && orig[4] < orig[3] ) )
	{ /* 0-4 diag split of quad */
  /* 4-1\
     |\| 2
     3-0/ */
	  tet_nodes[0] = orig[0];
	  tet_nodes[1] = orig[4];
	  tet_nodes[2] = orig[1];
	  tet_nodes[3] = orig[2];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  tet_nodes[0] = orig[0];
	  tet_nodes[1] = orig[3];
	  tet_nodes[2] = orig[4];
	  tet_nodes[3] = orig[2];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	}
      else
	{ /* 3-1 diag split of quad */
  /* 4-1\
     |/| 2
     3-0/ */
	  tet_nodes[0] = orig[0];
	  tet_nodes[1] = orig[3];
	  tet_nodes[2] = orig[1];
	  tet_nodes[3] = orig[2];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	  tet_nodes[0] = orig[1];
	  tet_nodes[1] = orig[3];
	  tet_nodes[2] = orig[4];
	  tet_nodes[3] = orig[2];
	  RSS( ref_cell_add( tet, tet_nodes, &new_cell ), "add tet");
	}
    }

  each_ref_cell_valid_cell_with_nodes( qua, cell, qua_nodes )
    {
      if ( mark[qua_nodes[0]] != REF_EMPTY &&
	   mark[qua_nodes[1]] != REF_EMPTY &&
	   mark[qua_nodes[2]] != REF_EMPTY &&
	   mark[qua_nodes[3]] != REF_EMPTY ) continue;

      RSS( ref_cell_remove( qua, cell ), "remove qua");
      tri_nodes[3] = qua_nodes[4]; /* patch id */
      if ( ( qua_nodes[0] < qua_nodes[1] && qua_nodes[0] < qua_nodes[3] ) ||
	   ( qua_nodes[2] < qua_nodes[1] && qua_nodes[2] < qua_nodes[3] ) )
	{ /* 0-2 diag split of quad */
  /* 2-1
     |\|
     3-0 */
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[2];
	  tri_nodes[2] = qua_nodes[3];
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[2];
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	}
      else
	{ /* 3-1 diag split of quad */
  /* 2-1
     |/|
     3-0 */
	  tri_nodes[0] = qua_nodes[0];
	  tri_nodes[1] = qua_nodes[1];
	  tri_nodes[2] = qua_nodes[3];
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	  tri_nodes[0] = qua_nodes[2];
	  tri_nodes[1] = qua_nodes[3];
	  tri_nodes[2] = qua_nodes[1];
	  RSS( ref_cell_add( tri, tri_nodes, &new_cell ), "add tri");
	}
    }

  ref_free(mark_copy);
  ref_free(mark);

  return REF_SUCCESS;
}
