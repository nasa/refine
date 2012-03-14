
#include <stdlib.h>
#include <stdio.h>

#include "ref_node.h"
#include "ref_sort.h"
#include "ref_malloc.h"
#include "ref_mpi.h"

/* REF_EMPTY is terminatior, next avalable is shifted by 2*/
#define next2index(next) (-(next)-2)
#define index2next(index) (-2-(index))

REF_STATUS ref_node_create( REF_NODE *ref_node_ptr )
{
  REF_INT max, node;
  REF_NODE ref_node;

  ref_malloc( *ref_node_ptr, 1, REF_NODE_STRUCT );

  ref_node = *ref_node_ptr;

  max = 20;

  ref_node_n(ref_node) = 0;
  ref_node_max(ref_node) = max;

  ref_malloc( ref_node->global, max, REF_INT );

  for (node=0;node<ref_node_max(ref_node);node++)
    ref_node->global[node] = index2next(node+1);
  ref_node->global[(ref_node->max)-1] = REF_EMPTY;
  ref_node->blank = index2next(0);

  ref_malloc( ref_node->sorted_global, max, REF_INT );
  ref_malloc( ref_node->sorted_local, max, REF_INT );

  ref_malloc( ref_node->part, max, REF_INT );

  ref_malloc( ref_node->xyz, 3*max, REF_DBL );

  RSS( ref_list_create( &(ref_node->unused_global_list) ), "create list");

  ref_node->old_n_global = REF_EMPTY;
  ref_node->new_n_global = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  if ( NULL == (void *)ref_node ) return REF_NULL;
  ref_list_free( ref_node->unused_global_list );
  ref_free( ref_node->xyz );
  ref_free( ref_node->part );
  ref_free( ref_node->sorted_local );
  ref_free( ref_node->sorted_global );
  ref_free( ref_node->global );
  ref_free( ref_node );
  return REF_SUCCESS;
}

REF_STATUS ref_node_inspect( REF_NODE ref_node )
{
  REF_INT node;
  printf("ref_node = %p\n",(void *)ref_node);
  printf(" n = %d\n",ref_node_n(ref_node));
  printf(" max = %d\n",ref_node_max(ref_node));
  printf(" blank = %d\n",ref_node->blank);
  for (node=0;node<ref_node_max(ref_node);node++)
    if ( 0 <= ref_node->global[node] ) 
      printf(" global[%d] = %3d; part[%d] = %3d;\n",
	     node,ref_node->global[node],
	     node,ref_node->part[node]	     );
  for (node=0;node<ref_node_n(ref_node);node++)
    printf(" sorted_global[%d] = %d sorted_local[%d] = %d\n",
	   node,ref_node->sorted_global[node],
	   node,ref_node->sorted_local[node]);
  printf(" old_n_global = %d\n",ref_node->old_n_global);
  printf(" new_n_global = %d\n",ref_node->new_n_global);
  return REF_SUCCESS;
}

REF_STATUS ref_node_location( REF_NODE ref_node, REF_INT node )
{
  printf("ref_node %d\n",node);
  if ( ref_node_valid(ref_node,node) )
    {
      printf ("(%.15e,%.15e,%.15e)\n",
	      ref_node_xyz(ref_node,0,node),
	      ref_node_xyz(ref_node,1,node),
	      ref_node_xyz(ref_node,2,node));
      
    }
  else
    {
      printf ("invalid\n");
    }

  return REF_SUCCESS;
}

REF_STATUS ref_node_tattle_global( REF_NODE ref_node, REF_INT global )
{
  REF_INT local_from_sorted;
  REF_INT local_from_exhaustive;
  REF_BOOL found_from_sorted;
  REF_BOOL found_from_exhaustive;
  REF_INT node;

  REF_STATUS ref_status;
  
  ref_status = ref_node_local(ref_node, global, &local_from_sorted );
  if ( REF_NOT_FOUND == ref_status )
    {
      found_from_sorted = REF_FALSE;
    }
  else
    {
      RSS( ref_status, "local search" );
      found_from_sorted = REF_TRUE;
    }

  found_from_exhaustive = REF_FALSE;
  local_from_exhaustive = REF_EMPTY;
  each_ref_node_valid_node( ref_node, node )
    {
      if ( global == ref_node_global(ref_node,node) )
	{
	  if ( found_from_exhaustive ) RSS(REF_FAILURE, "twice");
	  local_from_exhaustive = node;
	  found_from_exhaustive = REF_TRUE;
	}
    }

  printf("%d: global %d: search%d %d exhast%d %d\n",
	 ref_mpi_id, global,
	 found_from_sorted, local_from_sorted,
	 found_from_exhaustive, local_from_exhaustive );
      
  return REF_SUCCESS;
}

static REF_STATUS ref_node_add_core( REF_NODE ref_node, 
				     REF_INT global, REF_INT *node )
{
  REF_INT orig, chunk, extra;

  if ( global < 0 ) RSS( REF_INVALID, "invalid global node");

  if ( REF_EMPTY == ref_node->blank )
    {
      orig = ref_node_max(ref_node);
      chunk = 5000;
      ref_node->max = orig + chunk;
      ref_realloc( ref_node->global, ref_node_max(ref_node), REF_INT);
      for (extra=orig;extra < ref_node_max(ref_node); extra++ ) 
	  ref_node->global[extra] = index2next(extra+1);
      ref_node->global[ref_node_max(ref_node)-1] = REF_EMPTY; 
      ref_node->blank = index2next(orig);

      ref_realloc( ref_node->sorted_global, ref_node_max(ref_node), REF_INT);
      ref_realloc( ref_node->sorted_local, ref_node_max(ref_node), REF_INT);

      ref_realloc( ref_node->part, ref_node_max(ref_node), REF_INT);

      ref_realloc( ref_node->xyz, 3*ref_node_max(ref_node), REF_DBL);
    }

  *node = next2index(ref_node->blank);
  ref_node->blank = ref_node->global[*node];

  ref_node->global[*node] = global;
  ref_node->part[*node] = ref_mpi_id; /* default local node */

  (ref_node->n)++;
  return REF_SUCCESS;
}

REF_STATUS ref_node_add( REF_NODE ref_node, REF_INT global, REF_INT *node )
{
  REF_INT location, insert_point;
  REF_STATUS status;

  if ( global < 0 ) RSS( REF_INVALID, "invalid global node");

  status = ref_node_local( ref_node, global, node );
  if ( REF_SUCCESS == status ) return REF_SUCCESS;

  RSS( ref_node_add_core( ref_node,global, node ), "core");

  /* general case of non-ascending global node, requires: 
     search and shift (but looks to see if bigger than last early) */
  insert_point = 0;
  for (location=ref_node_n(ref_node)-2; location>=0; location--) {
    if (ref_node->sorted_global[location] < global) {
      insert_point = location+1;
      break;
    }
  }

  /* shift down to clear insert_point */
  for(location=ref_node_n(ref_node)-1;location>insert_point;location--)
    ref_node->sorted_global[location] = ref_node->sorted_global[location-1];
  for(location=ref_node_n(ref_node)-1;location>insert_point;location--)
    ref_node->sorted_local[location] = ref_node->sorted_local[location-1];

  /* insert in empty location */
  ref_node->sorted_global[insert_point] = global;
  ref_node->sorted_local[insert_point] = *node;

  return REF_SUCCESS;
}

REF_STATUS ref_node_add_many( REF_NODE ref_node, 
			      REF_INT n, REF_INT *global_orig )
{
  REF_STATUS status;
  REF_INT i, j, local, new;

  REF_INT *global;
  REF_INT *sorted;

  /* copy, removing existing nodes from list */

  ref_malloc( global, n, REF_INT );

  new = 0;
  for (i=0;i<n;i++)
    {
      status = ref_node_local( ref_node, global_orig[i], &local );
      if ( REF_NOT_FOUND == status ) 
	{
	  global[new] = global_orig[i];
	  new++;
	}
    }

  /* remove duplicates from list so core add can be used with existing check */

  ref_malloc( sorted, new, REF_INT );

  RSS( ref_sort_heap( new, global, sorted ), "heap" );

  j = 0;
  for (i=1;i<new;i++)
    {
      if ( global[sorted[i]] != global[sorted[j]] )
	{
	  j = i;
	  continue;
	}
      global[sorted[i]] = REF_EMPTY;
    }

  /* add remaining via core */

  for (i=0;i<new;i++)
    if ( REF_EMPTY != global[i] )
      {
	RSS( ref_node_add_core( ref_node, global[i], &local ), "add core" );
      }

  RSS( ref_node_rebuild_sorted_global( ref_node ), "rebuild globals" );

  ref_free( sorted );
  ref_free( global );

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove( REF_NODE ref_node, REF_INT node )
{
  REF_INT location, sorted_node;
  if ( ! ref_node_valid(ref_node,node) ) return REF_INVALID;
 
  RSS( ref_sort_search( ref_node_n(ref_node), ref_node->sorted_global, 
			ref_node->global[node], &location ), 
       "find global in sort list" );

  for(sorted_node=location;sorted_node<ref_node_n(ref_node);sorted_node++)
    ref_node->sorted_global[sorted_node]=ref_node->sorted_global[sorted_node+1];
  for(sorted_node=location;sorted_node<ref_node_n(ref_node);sorted_node++)
    ref_node->sorted_local[sorted_node]=ref_node->sorted_local[sorted_node+1];

  RSS( ref_list_add( ref_node->unused_global_list, ref_node->global[node] ),
       "store unused global" );
  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_remove_without_global( REF_NODE ref_node, REF_INT node )
{
  if ( ! ref_node_valid(ref_node,node) ) return REF_INVALID;
 
  ref_node->global[node] = ref_node->blank;
  ref_node->blank = index2next(node);

  (ref_node->n)--;

  return REF_SUCCESS;
}

REF_STATUS ref_node_rebuild_sorted_global( REF_NODE ref_node )
{
  REF_INT node, nnode, *pack;

  ref_malloc( pack, ref_node_n(ref_node), REF_INT );

  nnode = 0;
  each_ref_node_valid_node( ref_node, node )
    {
      ref_node->sorted_global[nnode] = ref_node->global[node];
      pack[nnode] = node;
      nnode++;
    }

  RSS( ref_sort_heap( ref_node_n(ref_node),
		      ref_node->sorted_global,
		      ref_node->sorted_local ), "heap" );
  
  for(node=0;node<ref_node_n(ref_node);node++)
    {
      ref_node->sorted_local[node]=pack[ref_node->sorted_local[node]];
      ref_node->sorted_global[node] = 
	ref_node->global[ref_node->sorted_local[node]];
    }

  free(pack);
  return REF_SUCCESS;
}

REF_STATUS ref_node_initialize_n_global(  REF_NODE ref_node, REF_INT n_global )
{
  ref_node->old_n_global = n_global;
  ref_node->new_n_global = n_global;

  return REF_SUCCESS;
}

REF_STATUS ref_node_next_global( REF_NODE ref_node, REF_INT *global )
{
  if ( 0 < ref_list_n( ref_node->unused_global_list ) )
    {
      RSS( ref_list_remove( ref_node->unused_global_list, global ), 
	   "grab an unused global from list");
    }
  else
    {
      if ( REF_EMPTY == ref_node->new_n_global )
	RSS( ref_node_initialize_n_global( ref_node, ref_node_n(ref_node) ), 
	     "init with n");
      (*global) = ref_node->new_n_global;
      (ref_node->new_n_global)++;
    }

  return REF_SUCCESS;
}
REF_STATUS ref_node_synchronize_globals( REF_NODE ref_node )
{
  RSS( ref_node_shift_new_globals( ref_node ), "shift" );
  RSS( ref_node_eliminate_unused_globals( ref_node ), "shift" );

  return REF_SUCCESS;
}

REF_STATUS ref_node_shift_new_globals( REF_NODE ref_node )
{
  REF_INT new_nodes;
  REF_INT *everyones_new_nodes;
  REF_INT offset, proc, total_new_nodes, node;

  ref_malloc( everyones_new_nodes, ref_mpi_n, REF_INT );

  new_nodes = ref_node->new_n_global - ref_node->old_n_global;

  RSS( ref_mpi_allgather( &new_nodes, everyones_new_nodes, REF_INT_TYPE ),
       "allgather");

  offset = 0;
  for( proc=0;proc<ref_mpi_id;proc++)
    offset += everyones_new_nodes[proc];

  total_new_nodes = 0;
  for( proc=0;proc<ref_mpi_n;proc++)
    total_new_nodes += everyones_new_nodes[proc];

  ref_free( everyones_new_nodes );

  if ( 0 != offset ) 
    {

      each_ref_node_valid_node( ref_node, node )
	if ( ref_node_global(ref_node,node) >= ref_node->old_n_global )
	  (ref_node->global[node]) += offset;

      for ( node = ref_node_n(ref_node)-1 ; 
	    node>=0 && ref_node->sorted_global[node] >= ref_node->old_n_global; 
	    node-- )
	ref_node->sorted_global[node] += offset;

      RSS( ref_list_shift( ref_node->unused_global_list, 
			   ref_node->old_n_global, offset ), "shift" );
    }

  RSS( ref_node_initialize_n_global( ref_node, 
				     total_new_nodes + ref_node->old_n_global),
       "re-init" );

  return REF_SUCCESS;
}

REF_STATUS ref_node_eliminate_unused_globals( REF_NODE ref_node )
{
  REF_LIST ref_list = ref_node->unused_global_list;
  REF_INT sort, offset, local;

  RSS( ref_list_sort( ref_list ), "sort unused global" );

  offset = 0;
  for (sort=0;sort<ref_node_n(ref_node);sort++) {
    while ( (offset < ref_list_n(ref_list) ) &&
	    ( ref_list_value(ref_list,offset) < 
	      ref_node->sorted_global[sort]  ) ) {
      offset++;
    }
    local = ref_node->sorted_local[sort];
    ref_node->global[local] -= offset; /* move to separate loop for cashe? */
    ref_node->sorted_global[sort] -= offset;
  }

  RSS( ref_node_initialize_n_global( ref_node, 
				     ref_node->old_n_global - 
				     ref_list_n(ref_list) ),
       "re-init" );

  RSS( ref_list_erase( ref_list ), "erase unused list");

  return REF_SUCCESS;
}

REF_STATUS ref_node_local( REF_NODE ref_node, REF_INT global, REF_INT *local )
{
  REF_INT location;

  (*local) = REF_EMPTY;

  RAISE( ref_sort_search( ref_node_n(ref_node), ref_node->sorted_global, 
			  global, &location ) );

  if ( (location) == REF_EMPTY ) return REF_NOT_FOUND;

  (*local) = ref_node->sorted_local[location];

  return REF_SUCCESS;
}

REF_STATUS ref_node_compact( REF_NODE ref_node, REF_INT **o2n_ptr )
{
  REF_INT node;
  REF_INT nnode;
  REF_INT *o2n;
  
  ref_malloc( *o2n_ptr, ref_node_max(ref_node), REF_INT );
  o2n = *o2n_ptr;

  nnode = 0;
  for ( node = 0 ; node < ref_node_max(ref_node) ; node++ )
    if ( ref_node_valid(ref_node,node) )
      {
	o2n[node] = nnode;
	nnode++;
      }
    else
      {
	o2n[node] = REF_EMPTY;
      }
  RES( nnode, ref_node_n(ref_node), "nnode miscount" );

  return REF_SUCCESS;
}


REF_STATUS ref_node_ghost_xyz( REF_NODE ref_node )
{
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_xyz, *b_xyz;
  REF_INT local;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc_init( a_size, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_size, ref_mpi_n, REF_INT, 0 );

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      a_size[ref_node_part(ref_node,node)]++;

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  ref_malloc( a_global, a_total, REF_INT );
  ref_malloc( a_xyz, 3*a_total, REF_DBL );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_global, b_total, REF_INT );
  ref_malloc( b_xyz, 3*b_total, REF_DBL );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	part = ref_node_part(ref_node,node);
	a_global[a_next[part]] = ref_node_global(ref_node,node);
	a_next[ref_node_part(ref_node,node)]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_size, b_global, b_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<b_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      b_xyz[0+3*node] = ref_node_xyz(ref_node,0,local);
      b_xyz[1+3*node] = ref_node_xyz(ref_node,1,local);
      b_xyz[2+3*node] = ref_node_xyz(ref_node,2,local);
    }

  RSS( ref_mpi_alltoallv( b_xyz, b_size, a_xyz, a_size, 
			  3, REF_DBL_TYPE ), 
       "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      ref_node_xyz(ref_node,0,local) = a_xyz[0+3*node];
      ref_node_xyz(ref_node,1,local) = a_xyz[1+3*node];
      ref_node_xyz(ref_node,2,local) = a_xyz[2+3*node];
    }

  free(a_next);
  free(b_xyz);
  free(b_global);
  free(a_xyz);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;  
}

REF_STATUS ref_node_ghost_int( REF_NODE ref_node, REF_INT *scalar )
{
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_INT *a_scalar, *b_scalar;
  REF_INT local;

  if ( 1 == ref_mpi_n ) return REF_SUCCESS;

  ref_malloc_init( a_size, ref_mpi_n, REF_INT, 0 );
  ref_malloc_init( b_size, ref_mpi_n, REF_INT, 0 );

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      a_size[ref_node_part(ref_node,node)]++;

  RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");

  a_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    a_total += a_size[part];
  ref_malloc( a_global, a_total, REF_INT );
  ref_malloc( a_scalar, a_total, REF_INT );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_global, b_total, REF_INT );
  ref_malloc( b_scalar, b_total, REF_INT );

  ref_malloc( a_next, ref_mpi_n, REF_INT );
  a_next[0] = 0;
  for ( part = 1; part<ref_mpi_n ; part++ )
    a_next[part] = a_next[part-1]+a_size[part-1];

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	part = ref_node_part(ref_node,node);
	a_global[a_next[part]] = ref_node_global(ref_node,node);
	a_next[ref_node_part(ref_node,node)]++;
      }

  RSS( ref_mpi_alltoallv( a_global, a_size, b_global, b_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<b_total;node++)
    {
      RSS( ref_node_local( ref_node, b_global[node], &local ), "g2l");
      b_scalar[node] = scalar[local];
    }

  RSS( ref_mpi_alltoallv( b_scalar, b_size, a_scalar, a_size, 
			  1, REF_INT_TYPE ), 
       "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      scalar[local] = a_scalar[node];
    }

  free(a_next);
  free(b_scalar);
  free(b_global);
  free(a_scalar);
  free(a_global);
  free(b_size);
  free(a_size);

  return REF_SUCCESS;  
}
