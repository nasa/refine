
#include <stdlib.h>
#include <stdio.h>

#include "ref_cell.h"
#include "ref_sort.h"
#include "ref_malloc.h"

REF_STATUS ref_cell_create( REF_CELL *ref_cell_ptr, 
			    REF_INT node_per, REF_BOOL last_node_is_an_id )
{
  REF_INT cell;
  REF_INT max;
  REF_CELL ref_cell;

  (*ref_cell_ptr) = NULL;
  
  if ( node_per+(last_node_is_an_id?1:0) > REF_CELL_MAX_SIZE_PER)
    {
      RSS( REF_FAILURE, "node_per limited to REF_CELL_MAX_SIZE_PER");
    }

  ref_malloc( *ref_cell_ptr, 1, REF_CELL_STRUCT );

  ref_cell = (*ref_cell_ptr);

  ref_cell_last_node_is_an_id(ref_cell) = last_node_is_an_id;

  ref_cell_node_per(ref_cell) = node_per;
  ref_cell_size_per(ref_cell) = node_per+(last_node_is_an_id?1:0);
  switch ( ref_cell_node_per(ref_cell) )
    {
    case 4:
      ref_cell_edge_per(ref_cell) = 6;
      break;
    case 5:
      ref_cell_edge_per(ref_cell) = 8;
      break;
    case 6:
      ref_cell_edge_per(ref_cell) = 9;
      break;
    case 8:
      ref_cell_edge_per(ref_cell) = 12;
      break;
    default:
      ref_cell_edge_per(ref_cell) = 0;
      break;    
    }

  ref_cell->e2n = NULL;
  if ( ref_cell_edge_per(ref_cell) > 0 )
    ref_malloc( ref_cell->e2n, 2 * ref_cell_edge_per(ref_cell), REF_INT);

  switch ( ref_cell_edge_per(ref_cell) )
    {
    case 6:
      ref_cell_e2n_gen(ref_cell,0,0)  = 0; ref_cell_e2n_gen(ref_cell,1,0)  = 1;
      ref_cell_e2n_gen(ref_cell,0,1)  = 0; ref_cell_e2n_gen(ref_cell,1,1)  = 2;
      ref_cell_e2n_gen(ref_cell,0,2)  = 0; ref_cell_e2n_gen(ref_cell,1,2)  = 3;
      ref_cell_e2n_gen(ref_cell,0,3)  = 1; ref_cell_e2n_gen(ref_cell,1,3)  = 2;
      ref_cell_e2n_gen(ref_cell,0,4)  = 1; ref_cell_e2n_gen(ref_cell,1,4)  = 3;
      ref_cell_e2n_gen(ref_cell,0,5)  = 2; ref_cell_e2n_gen(ref_cell,1,5)  = 3;
      break;
    case 8:
      ref_cell_e2n_gen(ref_cell,0,0)  = 0; ref_cell_e2n_gen(ref_cell,1,0)  = 1;
      ref_cell_e2n_gen(ref_cell,0,1)  = 0; ref_cell_e2n_gen(ref_cell,1,1)  = 2;
      ref_cell_e2n_gen(ref_cell,0,2)  = 0; ref_cell_e2n_gen(ref_cell,1,2)  = 3;
      ref_cell_e2n_gen(ref_cell,0,3)  = 1; ref_cell_e2n_gen(ref_cell,1,3)  = 2;
      ref_cell_e2n_gen(ref_cell,0,4)  = 1; ref_cell_e2n_gen(ref_cell,1,4)  = 4;
      ref_cell_e2n_gen(ref_cell,0,5)  = 2; ref_cell_e2n_gen(ref_cell,1,5)  = 3;
      ref_cell_e2n_gen(ref_cell,0,6)  = 2; ref_cell_e2n_gen(ref_cell,1,6)  = 4;
      ref_cell_e2n_gen(ref_cell,0,7)  = 3; ref_cell_e2n_gen(ref_cell,1,7)  = 4;
      break;
    case 9:
      ref_cell_e2n_gen(ref_cell,0,0)  = 0; ref_cell_e2n_gen(ref_cell,1,0)  = 1;
      ref_cell_e2n_gen(ref_cell,0,1)  = 0; ref_cell_e2n_gen(ref_cell,1,1)  = 2;
      ref_cell_e2n_gen(ref_cell,0,2)  = 0; ref_cell_e2n_gen(ref_cell,1,2)  = 3;
      ref_cell_e2n_gen(ref_cell,0,3)  = 1; ref_cell_e2n_gen(ref_cell,1,3)  = 2;
      ref_cell_e2n_gen(ref_cell,0,4)  = 1; ref_cell_e2n_gen(ref_cell,1,4)  = 4;
      ref_cell_e2n_gen(ref_cell,0,5)  = 2; ref_cell_e2n_gen(ref_cell,1,5)  = 5;
      ref_cell_e2n_gen(ref_cell,0,6)  = 3; ref_cell_e2n_gen(ref_cell,1,6)  = 4;
      ref_cell_e2n_gen(ref_cell,0,7)  = 3; ref_cell_e2n_gen(ref_cell,1,7)  = 5;
      ref_cell_e2n_gen(ref_cell,0,8)  = 4; ref_cell_e2n_gen(ref_cell,1,8)  = 5;
      break;
    case 12:
      ref_cell_e2n_gen(ref_cell,0,0)  = 0; ref_cell_e2n_gen(ref_cell,1,0)  = 1;
      ref_cell_e2n_gen(ref_cell,0,1)  = 0; ref_cell_e2n_gen(ref_cell,1,1)  = 3;
      ref_cell_e2n_gen(ref_cell,0,2)  = 0; ref_cell_e2n_gen(ref_cell,1,2)  = 4;
      ref_cell_e2n_gen(ref_cell,0,3)  = 1; ref_cell_e2n_gen(ref_cell,1,3)  = 2;
      ref_cell_e2n_gen(ref_cell,0,4)  = 1; ref_cell_e2n_gen(ref_cell,1,4)  = 5;
      ref_cell_e2n_gen(ref_cell,0,5)  = 2; ref_cell_e2n_gen(ref_cell,1,5)  = 3;
      ref_cell_e2n_gen(ref_cell,0,6)  = 2; ref_cell_e2n_gen(ref_cell,1,6)  = 6;
      ref_cell_e2n_gen(ref_cell,0,7)  = 3; ref_cell_e2n_gen(ref_cell,1,7)  = 7;
      ref_cell_e2n_gen(ref_cell,0,8)  = 4; ref_cell_e2n_gen(ref_cell,1,8)  = 5;
      ref_cell_e2n_gen(ref_cell,0,9)  = 4; ref_cell_e2n_gen(ref_cell,1,9)  = 7;
      ref_cell_e2n_gen(ref_cell,0,10) = 5; ref_cell_e2n_gen(ref_cell,1,10) = 6;
      ref_cell_e2n_gen(ref_cell,0,11) = 6; ref_cell_e2n_gen(ref_cell,1,11) = 7;
      break;
    }

  switch ( ref_cell_node_per(ref_cell) )
    {
    case 4:
      ref_cell_face_per(ref_cell) = 4;
      break;
    case 5:
      ref_cell_face_per(ref_cell) = 5;
      break;
    case 6:
      ref_cell_face_per(ref_cell) = 5;
      break;
    case 8:
      ref_cell_face_per(ref_cell) = 6;
      break;
    default:
      ref_cell_face_per(ref_cell) = 0;
      break;    
    }

  ref_cell->f2n = NULL;
  if ( ref_cell_face_per(ref_cell) > 0 )
    ref_malloc( ref_cell->f2n, 4 * ref_cell_face_per(ref_cell), REF_INT);

  switch ( ref_cell_node_per(ref_cell) )
    {
    case 4:
      ref_cell_f2n_gen(ref_cell,0,0) = 1; 
      ref_cell_f2n_gen(ref_cell,1,0) = 3;
      ref_cell_f2n_gen(ref_cell,2,0) = 2;
      ref_cell_f2n_gen(ref_cell,3,0) = ref_cell_f2n_gen(ref_cell,0,0);
      ref_cell_f2n_gen(ref_cell,0,1) = 0; 
      ref_cell_f2n_gen(ref_cell,1,1) = 2;
      ref_cell_f2n_gen(ref_cell,2,1) = 3;
      ref_cell_f2n_gen(ref_cell,3,1) = ref_cell_f2n_gen(ref_cell,0,1);
      ref_cell_f2n_gen(ref_cell,0,2) = 0; 
      ref_cell_f2n_gen(ref_cell,1,2) = 3;
      ref_cell_f2n_gen(ref_cell,2,2) = 1;
      ref_cell_f2n_gen(ref_cell,3,2) = ref_cell_f2n_gen(ref_cell,0,2);
      ref_cell_f2n_gen(ref_cell,0,3) = 0; 
      ref_cell_f2n_gen(ref_cell,1,3) = 1;
      ref_cell_f2n_gen(ref_cell,2,3) = 2;
      ref_cell_f2n_gen(ref_cell,3,3) = ref_cell_f2n_gen(ref_cell,0,3);
      break;
    case 5:
      ref_cell_f2n_gen(ref_cell,0,0) = 0; 
      ref_cell_f2n_gen(ref_cell,1,0) = 1;
      ref_cell_f2n_gen(ref_cell,2,0) = 2;
      ref_cell_f2n_gen(ref_cell,3,0) = ref_cell_f2n_gen(ref_cell,0,0);
      ref_cell_f2n_gen(ref_cell,0,1) = 1; 
      ref_cell_f2n_gen(ref_cell,1,1) = 4;
      ref_cell_f2n_gen(ref_cell,2,1) = 2;
      ref_cell_f2n_gen(ref_cell,3,1) = ref_cell_f2n_gen(ref_cell,0,1);
      ref_cell_f2n_gen(ref_cell,0,2) = 2; 
      ref_cell_f2n_gen(ref_cell,1,2) = 4;
      ref_cell_f2n_gen(ref_cell,2,2) = 3;
      ref_cell_f2n_gen(ref_cell,3,2) = ref_cell_f2n_gen(ref_cell,0,2);
      ref_cell_f2n_gen(ref_cell,0,3) = 0; 
      ref_cell_f2n_gen(ref_cell,1,3) = 2;
      ref_cell_f2n_gen(ref_cell,2,3) = 3;
      ref_cell_f2n_gen(ref_cell,3,3) = ref_cell_f2n_gen(ref_cell,0,3);
      ref_cell_f2n_gen(ref_cell,0,4) = 0; 
      ref_cell_f2n_gen(ref_cell,1,4) = 3;
      ref_cell_f2n_gen(ref_cell,2,4) = 4;
      ref_cell_f2n_gen(ref_cell,3,4) = 1;
      break;
    case 6:
      ref_cell_f2n_gen(ref_cell,0,0) = 0; 
      ref_cell_f2n_gen(ref_cell,1,0) = 3;
      ref_cell_f2n_gen(ref_cell,2,0) = 4;
      ref_cell_f2n_gen(ref_cell,3,0) = 1;

      ref_cell_f2n_gen(ref_cell,0,1) = 1; 
      ref_cell_f2n_gen(ref_cell,1,1) = 4;
      ref_cell_f2n_gen(ref_cell,2,1) = 5;
      ref_cell_f2n_gen(ref_cell,3,1) = 2;

      ref_cell_f2n_gen(ref_cell,0,2) = 0; 
      ref_cell_f2n_gen(ref_cell,1,2) = 2;
      ref_cell_f2n_gen(ref_cell,2,2) = 5;
      ref_cell_f2n_gen(ref_cell,3,2) = 3;

      ref_cell_f2n_gen(ref_cell,0,3) = 0; 
      ref_cell_f2n_gen(ref_cell,1,3) = 1;
      ref_cell_f2n_gen(ref_cell,2,3) = 2;
      ref_cell_f2n_gen(ref_cell,3,3) = ref_cell_f2n_gen(ref_cell,0,3);

      ref_cell_f2n_gen(ref_cell,0,4) = 3; 
      ref_cell_f2n_gen(ref_cell,1,4) = 5;
      ref_cell_f2n_gen(ref_cell,2,4) = 4;
      ref_cell_f2n_gen(ref_cell,3,4) = ref_cell_f2n_gen(ref_cell,0,4);
      break;
    case 8:
      ref_cell_f2n_gen(ref_cell,0,0) = 0; 
      ref_cell_f2n_gen(ref_cell,1,0) = 4;
      ref_cell_f2n_gen(ref_cell,2,0) = 5;
      ref_cell_f2n_gen(ref_cell,3,0) = 1;

      ref_cell_f2n_gen(ref_cell,0,1) = 1; 
      ref_cell_f2n_gen(ref_cell,1,1) = 5;
      ref_cell_f2n_gen(ref_cell,2,1) = 6;
      ref_cell_f2n_gen(ref_cell,3,1) = 2;

      ref_cell_f2n_gen(ref_cell,0,2) = 2; 
      ref_cell_f2n_gen(ref_cell,1,2) = 6;
      ref_cell_f2n_gen(ref_cell,2,2) = 7;
      ref_cell_f2n_gen(ref_cell,3,2) = 3;

      ref_cell_f2n_gen(ref_cell,0,3) = 0; 
      ref_cell_f2n_gen(ref_cell,1,3) = 3;
      ref_cell_f2n_gen(ref_cell,2,3) = 7;
      ref_cell_f2n_gen(ref_cell,3,3) = 4;

      ref_cell_f2n_gen(ref_cell,0,4) = 0; 
      ref_cell_f2n_gen(ref_cell,1,4) = 1;
      ref_cell_f2n_gen(ref_cell,2,4) = 2;
      ref_cell_f2n_gen(ref_cell,3,4) = 3;

      ref_cell_f2n_gen(ref_cell,0,5) = 4; 
      ref_cell_f2n_gen(ref_cell,1,5) = 7;
      ref_cell_f2n_gen(ref_cell,2,5) = 6;
      ref_cell_f2n_gen(ref_cell,3,5) = 5;

      break;
    }

  max = 100;

  ref_cell_n(ref_cell) = 0;
  ref_cell_max(ref_cell) = max;

  ref_malloc( ref_cell->c2n, ref_cell_max(ref_cell) *
	      ref_cell_size_per(ref_cell), REF_INT);

  ref_malloc( ref_cell->c2e, ref_cell_max(ref_cell) *
	      ref_cell_edge_per(ref_cell), REF_INT);

  for ( cell = 0 ; cell < max ; cell++ ) 
    {
      ref_cell_c2n(ref_cell,0,cell) = REF_EMPTY;
      ref_cell_c2n(ref_cell,1,cell) = cell+1;
    }
  ref_cell_c2n(ref_cell,1,max-1) = REF_EMPTY;
  ref_cell_blank(ref_cell) = 0;

  RSS( ref_adj_create( &(ref_cell->ref_adj) ), "create ref_adj for ref_cell" );

  return REF_SUCCESS;
}

REF_STATUS ref_cell_free( REF_CELL ref_cell )
{
  if ( NULL == (void *)ref_cell ) return REF_NULL;
  ref_adj_free( ref_cell->ref_adj );
  ref_free( ref_cell->c2n );
  ref_free( ref_cell->c2e );
  ref_free( ref_cell->f2n );
  ref_free( ref_cell->e2n );
  ref_free( ref_cell );
  return REF_SUCCESS;
}

REF_STATUS ref_cell_inspect( REF_CELL ref_cell )
{
  REF_INT cell, node;
  REF_INT cell_edge;
  printf("ref_cell = %p\n",(void *)ref_cell);
  printf(" size_per = %d\n",ref_cell_size_per(ref_cell));
  printf(" node_per = %d\n",ref_cell_node_per(ref_cell));
  printf(" n = %d\n",ref_cell_n(ref_cell));
  printf(" max = %d\n",ref_cell_max(ref_cell));
  printf(" blank = %d\n",ref_cell->blank);
  for (cell=0;cell<ref_cell_max(ref_cell);cell++)
    {
      if ( ref_cell_valid(ref_cell,cell) )
	{
	  printf(" %d:",cell);
	  for (node=0;node<ref_cell_size_per(ref_cell);node++)
	    printf(" %d",ref_cell_c2n(ref_cell,node,cell));
	  printf("\n");
	}
    }
  
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      printf(" cell %d\n",cell);
      each_ref_cell_cell_edge( ref_cell, cell_edge )
	printf("  edge %d nodes %d %d\n",
	       ref_cell_c2e(ref_cell,cell_edge,cell),
	       ref_cell_e2n(ref_cell,0,cell,cell_edge),
	       ref_cell_e2n(ref_cell,1,cell,cell_edge));
      printf("\n");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_taddle( REF_CELL ref_cell, REF_INT cell )
{
  REF_INT node;
  printf("cell %d:",cell);
  for (node=0;node<ref_cell_size_per(ref_cell);node++)
    {
      printf(" %d",ref_cell_c2n(ref_cell,node,cell));
    }
  printf("\n");
  return REF_SUCCESS;
}

REF_STATUS ref_cell_add( REF_CELL ref_cell, REF_INT *nodes, REF_INT *new_cell )
{
  REF_INT node, cell;
  REF_INT orig, chunk;

  (*new_cell) = REF_EMPTY;

  if ( REF_EMPTY == ref_cell_blank(ref_cell) ) 
    {
      orig = ref_cell_max(ref_cell);
      chunk = 5000;
      ref_cell->max = orig + chunk;

      ref_realloc( ref_cell->c2n, ref_cell_size_per(ref_cell) *
		   ref_cell_max(ref_cell), REF_INT );
      ref_realloc( ref_cell->c2e, ref_cell_edge_per(ref_cell) *
		   ref_cell_max(ref_cell), REF_INT );

      for (cell=orig;cell < ref_cell_max(ref_cell); cell++ ) 
	{
	  ref_cell_c2n(ref_cell,0,cell)= REF_EMPTY; 
	  ref_cell_c2n(ref_cell,1,cell) = cell+1; 
	}
      ref_cell_c2n(ref_cell,1,(ref_cell->max)-1) = REF_EMPTY; 
      ref_cell_blank(ref_cell) = orig;
    }

  cell = ref_cell_blank(ref_cell);
  ref_cell_blank(ref_cell) = ref_cell_c2n(ref_cell,1,cell);
  for ( node = 0 ; node < ref_cell_size_per(ref_cell) ; node++ )
    ref_cell_c2n(ref_cell,node,cell) = nodes[node];

  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    RSS( ref_adj_add(ref_cell->ref_adj, nodes[node], cell), 
	 "register cell" );

  ref_cell_n(ref_cell)++;

  (*new_cell) = cell;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_add_many_global( REF_CELL ref_cell, REF_NODE ref_node,
				     REF_INT n, REF_INT *c2n, REF_INT *part,
				     REF_INT exclude_part_id )
{
  REF_INT *global;
  REF_INT nnode;
  REF_INT node, cell;
  REF_INT local, local_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;

  ref_malloc( global, ref_cell_node_per(ref_cell)*n, REF_INT );

  nnode = 0;
  for ( cell=0; cell<n; cell++ )
    for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
      if ( exclude_part_id != part[node+ref_cell_size_per(ref_cell)*cell] )
	{
	global[nnode] =	c2n[node+ref_cell_size_per(ref_cell)*cell];
	nnode++;
      }
  RSS( ref_node_add_many( ref_node, nnode, global ), "many nodes" );

  ref_free( global );

  /* set parts */
  for ( cell=0; cell<n; cell++ )
    {
      for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
	{
	  RSS( ref_node_local( ref_node,
			       c2n[node+ref_cell_size_per(ref_cell)*cell], 
			       &local ), "local" );
	  ref_node_part(ref_node,local) = 
	    part[node+ref_cell_size_per(ref_cell)*cell];
	  local_nodes[node] = local;
	}
      if ( ref_cell_last_node_is_an_id(ref_cell) )
	local_nodes[ref_cell_size_per(ref_cell)-1] = 
	  c2n[(ref_cell_size_per(ref_cell)-1)+ref_cell_size_per(ref_cell)*cell];
      
      RXS( ref_cell_with( ref_cell, local_nodes, &new_cell),
	   REF_NOT_FOUND, "with failed");
      
      if ( REF_EMPTY == new_cell ) 
	RSS( ref_cell_add( ref_cell, local_nodes, &new_cell ), "add cell") ;
    }

  return REF_SUCCESS;
}


REF_STATUS ref_cell_remove( REF_CELL ref_cell, REF_INT cell )
{
  REF_INT node;
  if ( !ref_cell_valid(ref_cell,cell) ) return REF_INVALID;
  ref_cell_n(ref_cell)--;

  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    RSS( ref_adj_remove(ref_cell->ref_adj, 
			ref_cell_c2n(ref_cell,node,cell), cell), 
	 "unregister cell" );

  ref_cell_c2n(ref_cell,0,cell) = REF_EMPTY;
  ref_cell_c2n(ref_cell,1,cell) = ref_cell_blank(ref_cell); 

  return REF_SUCCESS;
}

REF_STATUS ref_cell_replace( REF_CELL ref_cell, REF_INT cell, REF_INT *nodes )
{
  REF_INT node;
  if ( !ref_cell_valid(ref_cell,cell) ) return REF_FAILURE;

  for ( node = 0 ; node < ref_cell_node_per(ref_cell) ; node++ )
    {
      RSS( ref_adj_remove(ref_cell->ref_adj, 
			  ref_cell_c2n(ref_cell,node,cell), cell), 
	   "unregister cell" );
      ref_cell_c2n(ref_cell,node,cell) = nodes[node];
      RSS( ref_adj_add(ref_cell->ref_adj, nodes[node], cell), 
	   "register cell with id" );
    }

  if ( ref_cell_last_node_is_an_id(ref_cell) )
    {
      node = ref_cell_size_per(ref_cell)-1;
      ref_cell_c2n(ref_cell,node,cell) = nodes[node];
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_nodes( REF_CELL ref_cell, REF_INT cell, REF_INT *nodes )
{
  REF_INT node;
  if ( cell < 0 || cell > ref_cell_max(ref_cell) ) return REF_INVALID;
  if ( REF_EMPTY == ref_cell_c2n(ref_cell,0,cell) ) 
    return REF_INVALID;
  for ( node = 0 ; node < ref_cell_size_per(ref_cell) ; node++ )
    nodes[node] = ref_cell_c2n(ref_cell,node,cell);
  return REF_SUCCESS;
}

static REF_STATUS ref_cell_make_canonical( REF_INT n, 
					   REF_INT *original, 
					   REF_INT *canonical )
{
  RSS( ref_sort_insertion_int( n, original, canonical ), "sort" );

  return REF_SUCCESS;
}

REF_STATUS ref_cell_has_side( REF_CELL ref_cell, 
			      REF_INT node0, REF_INT node1, 
			      REF_BOOL *has_side)
{
  REF_INT item, cell;
  REF_INT cell_edge;

  *has_side = REF_FALSE;

  each_ref_adj_node_item_with_ref( ref_cell_adj(ref_cell), node0, item, cell)
    each_ref_cell_cell_edge( ref_cell, cell_edge )
      if ( ( node0 == ref_cell_e2n(ref_cell,0,cell,cell_edge) &&
	     node1 == ref_cell_e2n(ref_cell,1,cell,cell_edge) ) ||
	   ( node0 == ref_cell_e2n(ref_cell,1,cell,cell_edge) &&
	     node1 == ref_cell_e2n(ref_cell,0,cell,cell_edge) ) )
	{	   
	  *has_side = REF_TRUE;
	  return REF_SUCCESS;
	}

  return REF_SUCCESS;
}

REF_STATUS ref_cell_with( REF_CELL ref_cell, REF_INT *nodes, REF_INT *cell )
{
  REF_INT item, ref, node, same;
  REF_INT target[REF_CELL_MAX_SIZE_PER];
  REF_INT canidate[REF_CELL_MAX_SIZE_PER], orig[REF_CELL_MAX_SIZE_PER];

  (*cell) = REF_EMPTY;

  RSS( ref_cell_make_canonical( ref_cell_node_per(ref_cell), 
				nodes, target ), "canonical" );

  each_ref_adj_node_item_with_ref( ref_cell_adj(ref_cell), nodes[0], item, ref)
    {
      RSS( ref_cell_nodes( ref_cell, ref, orig ), "get orig");
      RSS( ref_cell_make_canonical( ref_cell_node_per(ref_cell), 
				    orig, canidate ), "canonical" );
      same = 0;
      for (node=0;node<ref_cell_node_per(ref_cell); node++)
	if ( target[node] == canidate[node]) same++;
      if ( ref_cell_node_per(ref_cell) == same )
	{
	  (*cell) = ref;
	  return REF_SUCCESS;
	}
    }
  
  return REF_NOT_FOUND;
}

REF_STATUS ref_cell_list_with( REF_CELL ref_cell, 
			       REF_INT node0, REF_INT node1,
			       REF_INT max_cell, REF_INT *ncell,
			       REF_INT *cell_list ) 
{
  REF_INT cell, item, node;

  *ncell = 0;
  each_ref_cell_having_node( ref_cell, node0, item, cell )
    for ( node = 0 ; node < ref_cell_node_per(ref_cell); node++ )
      if ( node1 == ref_cell_c2n(ref_cell,node,cell) )
        {
          if ( *ncell >= max_cell )
            RSS( REF_INCREASE_LIMIT, "max_cell too small" );
          cell_list[*ncell] = cell;
	  (*ncell)++;
          continue; /* node loop */
        }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_empty_edges( REF_CELL ref_cell)
{
  REF_INT cell, edge;
  for ( cell = 0 ; cell < ref_cell_max(ref_cell) ; cell++ )
    for ( edge = 0 ; edge < ref_cell_edge_per(ref_cell) ; edge++ )
      ref_cell_c2e(ref_cell,edge,cell) = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_cell_set_edge( REF_CELL ref_cell, 
			      REF_INT n0, REF_INT n1, REF_INT edge)
{
  REF_INT item, cell, cell_edge;
  REF_INT e0, e1;

  each_ref_cell_having_node( ref_cell, n0, item, cell)
    {
      for (cell_edge = 0; cell_edge < ref_cell_edge_per(ref_cell); cell_edge++)
	{
	  e0 = ref_cell_e2n(ref_cell,0,cell,cell_edge);
	  e1 = ref_cell_e2n(ref_cell,1,cell,cell_edge);
	  if ( ( e0 == n0 && e1 == n1 ) ||
	       ( e0 == n1 && e1 == n0 )  )
	    ref_cell_c2e(ref_cell,cell_edge,cell) = edge;
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_cell_gen_edge_face( REF_CELL ref_cell, REF_INT edge, 
				   REF_INT *face0, REF_INT *face1 )
{
  REF_INT face, node0, node1;
  REF_BOOL have_node0;
  REF_BOOL have_node1;

  *face0 = REF_EMPTY;
  *face1 = REF_EMPTY;

  node0 = ref_cell_e2n_gen(ref_cell,0,edge);
  node1 = ref_cell_e2n_gen(ref_cell,1,edge);

  for(face=0;face<ref_cell_face_per(ref_cell);face++)
    {
      have_node0 = ( node0 == ref_cell_f2n_gen(ref_cell,0,face) ||
		     node0 == ref_cell_f2n_gen(ref_cell,1,face) ||
		     node0 == ref_cell_f2n_gen(ref_cell,2,face) ||
		     node0 == ref_cell_f2n_gen(ref_cell,3,face) );
      have_node1 = ( node1 == ref_cell_f2n_gen(ref_cell,0,face) ||
		     node1 == ref_cell_f2n_gen(ref_cell,1,face) ||
		     node1 == ref_cell_f2n_gen(ref_cell,2,face) ||
		     node1 == ref_cell_f2n_gen(ref_cell,3,face) );
      if ( have_node0 && have_node1 )
	{
	  if ( (*face0) == REF_EMPTY )
	    {
	      (*face0) = face;
	    }
	  else
	    {
	      RAS( REF_EMPTY == (*face1), "face1 set twice" );
	      (*face1) = face;
	    }
	}
    }

  RAS( REF_EMPTY != (*face0), "face0 not set" );
  RAS( REF_EMPTY != (*face1), "face1 not set" );

  return REF_SUCCESS;
}
