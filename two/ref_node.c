
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_node.h"
#include "ref_sort.h"
#include "ref_malloc.h"
#include "ref_mpi.h"
#include "ref_math.h"
#include "ref_matrix.h"

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

  ref_malloc( ref_node->real, REF_NODE_REAL_PER*max, REF_DBL );

  ref_node_naux(ref_node) = 0;
  ref_node->aux = NULL;

  RSS( ref_list_create( &(ref_node->unused_global_list) ), "create list");

  ref_node->old_n_global = REF_EMPTY;
  ref_node->new_n_global = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_node_free( REF_NODE ref_node )
{
  if ( NULL == (void *)ref_node ) return REF_NULL;
  ref_list_free( ref_node->unused_global_list );
  ref_free( ref_node->aux );
  ref_free( ref_node->real );
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

      ref_realloc( ref_node->real, 
		   REF_NODE_REAL_PER*ref_node_max(ref_node), REF_DBL);

      if ( ref_node_naux(ref_node) > 0 )
	ref_realloc( ref_node->aux, 
		     ref_node_naux(ref_node)*ref_node_max(ref_node), REF_DBL);	
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

  RSS( ref_sort_heap_int( new, global, sorted ), "heap" );

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
  REF_INT location, sorted_node;
  if ( ! ref_node_valid(ref_node,node) ) return REF_INVALID;
 
  RSS( ref_sort_search( ref_node_n(ref_node), ref_node->sorted_global, 
			ref_node->global[node], &location ), 
       "find global in sort list" );

  for(sorted_node=location;sorted_node<ref_node_n(ref_node);sorted_node++)
    ref_node->sorted_global[sorted_node]=ref_node->sorted_global[sorted_node+1];
  for(sorted_node=location;sorted_node<ref_node_n(ref_node);sorted_node++)
    ref_node->sorted_local[sorted_node]=ref_node->sorted_local[sorted_node+1];

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

  RSS( ref_sort_heap_int( ref_node_n(ref_node),
			  ref_node->sorted_global,
			  ref_node->sorted_local ), "heap" );
  
  for(node=0;node<ref_node_n(ref_node);node++)
    {
      ref_node->sorted_local[node]=pack[ref_node->sorted_local[node]];
      ref_node->sorted_global[node] = 
	ref_node->global[ref_node->sorted_local[node]];
    }

  ref_free(pack);
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

  RSS( ref_list_allgather( ref_list ), "gather unused global" );
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

REF_STATUS ref_node_compact( REF_NODE ref_node,
			     REF_INT **o2n_ptr, REF_INT **n2o_ptr )
{
  REF_INT node;
  REF_INT nnode;
  REF_INT *o2n, *n2o;
  
  ref_malloc_init( *o2n_ptr, ref_node_max(ref_node), REF_INT, REF_EMPTY );
  o2n = *o2n_ptr;
  ref_malloc( *n2o_ptr, ref_node_n(ref_node), REF_INT );
  n2o = *n2o_ptr;

  nnode = 0;    
  
  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id == ref_node_part(ref_node,node) )
      {
	o2n[node] = nnode;
	nnode++;
      }

  each_ref_node_valid_node( ref_node, node )
    if ( ref_mpi_id != ref_node_part(ref_node,node) )
      {
	o2n[node] = nnode;
	nnode++;
      }

  RES( nnode, ref_node_n(ref_node), "nnode miscount" );

  each_ref_node_valid_node( ref_node, node )
    n2o[o2n[node]] = node;

  return REF_SUCCESS;
}


REF_STATUS ref_node_ghost_real( REF_NODE ref_node )
{
  REF_INT *a_size, *b_size;
  REF_INT a_total, b_total;
  REF_INT *a_global, *b_global;
  REF_INT part, node;
  REF_INT *a_next;
  REF_DBL *a_real, *b_real;
  REF_DBL *a_aux, *b_aux;
  REF_INT local;
  REF_INT i;

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
  ref_malloc( a_real, REF_NODE_REAL_PER*a_total, REF_DBL );
  a_aux = NULL;
  if ( ref_node_naux(ref_node) > 0 )
    ref_malloc( a_aux, ref_node_naux(ref_node)*a_total, REF_DBL );

  b_total = 0;
  for ( part = 0; part<ref_mpi_n ; part++ )
    b_total += b_size[part];
  ref_malloc( b_global, b_total, REF_INT );
  ref_malloc( b_real, REF_NODE_REAL_PER*b_total, REF_DBL );
  b_aux = NULL;
  if ( ref_node_naux(ref_node) > 0 )
    ref_malloc( b_aux, ref_node_naux(ref_node)*b_total, REF_DBL );

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
      for ( i=0; i < REF_NODE_REAL_PER ; i++ )
	b_real[i+REF_NODE_REAL_PER*node] = ref_node_real(ref_node,i,local);
      for ( i=0; i < ref_node_naux(ref_node) ; i++ )
	b_aux[i+ref_node_naux(ref_node)*node] = ref_node_aux(ref_node,i,local);
    }

  RSS( ref_mpi_alltoallv( b_real, b_size, a_real, a_size, 
			  REF_NODE_REAL_PER, REF_DBL_TYPE ), 
       "alltoallv global");

  if ( ref_node_naux(ref_node) > 0 )
    RSS( ref_mpi_alltoallv( b_aux, b_size, a_aux, a_size, 
			    ref_node_naux(ref_node), REF_DBL_TYPE ), 
	 "alltoallv global");

  for (node=0;node<a_total;node++)
    {
      RSS( ref_node_local( ref_node, a_global[node], &local ), "g2l");
      for ( i=0; i < REF_NODE_REAL_PER ; i++ )
	ref_node_real(ref_node,i,local) = a_real[i+REF_NODE_REAL_PER*node];
      for ( i=0; i < ref_node_naux(ref_node) ; i++ )
	ref_node_aux(ref_node,i,local) = a_aux[i+ref_node_naux(ref_node)*node];
    }

  ref_free(a_next);
  ref_free(b_aux);
  ref_free(b_real);
  ref_free(b_global);
  ref_free(a_aux);
  ref_free(a_real);
  ref_free(a_global);
  ref_free(b_size);
  ref_free(a_size);

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

REF_STATUS ref_node_ratio( REF_NODE ref_node, REF_INT node0, REF_INT node1, 
			   REF_DBL *ratio )
{
  REF_DBL direction[3], length;
  REF_DBL ratio0, ratio1;
  REF_DBL r, r_min, r_max;

  if ( !ref_node_valid(ref_node,node0) ||
       !ref_node_valid(ref_node,node1) ) 
    RSS( REF_INVALID, "node invalid" );

  direction[0] = ( ref_node_xyz(ref_node,0,node1) -
		   ref_node_xyz(ref_node,0,node0) );
  direction[1] = ( ref_node_xyz(ref_node,1,node1) -
		   ref_node_xyz(ref_node,1,node0) );
  direction[2] = ( ref_node_xyz(ref_node,2,node1) -
		   ref_node_xyz(ref_node,2,node0) );

  length = ref_math_dot(direction,direction);
  length = sqrt(length);
	  
  if ( !ref_math_divisible(direction[0],length) ||
       !ref_math_divisible(direction[1],length) ||
       !ref_math_divisible(direction[2],length) ) 
    {
      *ratio = 0.0;
      return REF_SUCCESS;  
    }

  ratio0 = ref_matrix_sqrt_vt_m_v( ref_node_metric_ptr(ref_node,node0), 
				   direction );
  ratio1 = ref_matrix_sqrt_vt_m_v( ref_node_metric_ptr(ref_node,node1), 
				   direction );

  /* Loseille Lohner IMR 18 (2009) pg 613 */
  /* Alauzet Finite Elements in Analysis and Design 46 (2010) pg 185 */
  
  if ( ratio0 < 1.0e-12 || ratio1 < 1.0e-12 )
    {
      *ratio = MIN(ratio0,ratio1);
      return REF_SUCCESS;  
    }

  r_min = MIN( ratio0, ratio1 );
  r_max = MAX( ratio0, ratio1 );

  r = r_min/r_max;

  if ( ABS(r-1.0) < 1.0e-12 )
    {
      *ratio = 0.5*(ratio0+ratio1);
      return REF_SUCCESS;  
    }    
 
  *ratio = r_min * (r-1.0) / ( r * log(r) );

  return REF_SUCCESS;  
}

REF_STATUS ref_node_tet_quality( REF_NODE ref_node, 
				 REF_INT *nodes, 
				 REF_DBL *quality )
{
  REF_DBL l0,l1,l2,l3,l4,l5;

  REF_DBL min_det, volume;
  REF_DBL volume_in_metric;
  REF_DBL num, denom;

  RSS( ref_node_ratio( ref_node, nodes[0], nodes[1], &l0 ), "l0" );
  RSS( ref_node_ratio( ref_node, nodes[0], nodes[2], &l1 ), "l1" );
  RSS( ref_node_ratio( ref_node, nodes[0], nodes[3], &l2 ), "l2" );
  RSS( ref_node_ratio( ref_node, nodes[1], nodes[2], &l3 ), "l3" );
  RSS( ref_node_ratio( ref_node, nodes[1], nodes[3], &l4 ), "l4" );
  RSS( ref_node_ratio( ref_node, nodes[2], nodes[3], &l5 ), "l5" );
  
  RSS( ref_node_tet_vol( ref_node, nodes, &volume ), "vol");

  if ( volume <= 0.0 )
    {
      *quality = volume;
       return REF_SUCCESS;
    }

  min_det = MIN( 
		MIN( 
		    ref_matrix_det_m(ref_node_metric_ptr(ref_node, nodes[0])),
		    ref_matrix_det_m(ref_node_metric_ptr(ref_node, nodes[1]))
		    ),
		MIN( 
		    ref_matrix_det_m(ref_node_metric_ptr(ref_node, nodes[2])),
		    ref_matrix_det_m(ref_node_metric_ptr(ref_node, nodes[3]))
		    )
		);

  volume_in_metric = sqrt( min_det ) * volume;

  num = pow(volume_in_metric,2.0/3.0);
  denom = l0*l0 + l1*l1 + l2*l2 + l3*l3 + l4*l4 + l5*l5;

  if ( ref_math_divisible(num,denom) )
    {
      /* 36/3^(1/3) */
      *quality = 24.9610058766228 * num / denom;
    }
  else
    {
      *quality = -1.0;
      RSS( REF_DIV_ZERO, "in quality");
    }

  return REF_SUCCESS;  
}

REF_STATUS ref_node_tri_normal( REF_NODE ref_node, 
				REF_INT *nodes, 
				REF_DBL *normal )
{
  REF_DBL *xyz0, *xyz1, *xyz2;
  REF_DBL edge10[3], edge20[3];

  if ( !ref_node_valid(ref_node,nodes[0]) ||
       !ref_node_valid(ref_node,nodes[1]) ||
       !ref_node_valid(ref_node,nodes[2]) ) 
    RSS( REF_INVALID, "node invalid" );

  xyz0 = ref_node_xyz_ptr(ref_node,nodes[0]);
  xyz1 = ref_node_xyz_ptr(ref_node,nodes[1]);
  xyz2 = ref_node_xyz_ptr(ref_node,nodes[2]);

  edge10[0] = xyz1[0] - xyz0[0];
  edge10[1] = xyz1[1] - xyz0[1];
  edge10[2] = xyz1[2] - xyz0[2];
  
  edge20[0] = xyz2[0] - xyz0[0];
  edge20[1] = xyz2[1] - xyz0[1];
  edge20[2] = xyz2[2] - xyz0[2];

  ref_math_cross_product(edge10,edge20,normal);

  normal[0] *= 0.5;
  normal[1] *= 0.5;
  normal[2] *= 0.5;

  return REF_SUCCESS;  
}

REF_STATUS ref_node_tet_vol( REF_NODE ref_node, 
			     REF_INT *nodes, 
			     REF_DBL *volume )
{
  REF_DBL *a, *b, *c, *d;
  REF_DBL m11, m12, m13;
  REF_DBL det;

  if ( !ref_node_valid(ref_node,nodes[0]) ||
       !ref_node_valid(ref_node,nodes[1]) ||
       !ref_node_valid(ref_node,nodes[2]) ||
       !ref_node_valid(ref_node,nodes[3]) ) 
    RSS( REF_INVALID, "node invalid" );

  a = ref_node_xyz_ptr(ref_node,nodes[0]);
  b = ref_node_xyz_ptr(ref_node,nodes[1]);
  c = ref_node_xyz_ptr(ref_node,nodes[2]);
  d = ref_node_xyz_ptr(ref_node,nodes[3]);
  
  m11 = (a[0]-d[0])*((b[1]-d[1])*(c[2]-d[2])-(c[1]-d[1])*(b[2]-d[2]));
  m12 = (a[1]-d[1])*((b[0]-d[0])*(c[2]-d[2])-(c[0]-d[0])*(b[2]-d[2]));
  m13 = (a[2]-d[2])*((b[0]-d[0])*(c[1]-d[1])-(c[0]-d[0])*(b[1]-d[1]));
  det = ( m11 - m12 + m13 );

  *volume = -det/6.0;

  return REF_SUCCESS;  
}

REF_STATUS ref_node_interpolate_edge( REF_NODE ref_node, 
				      REF_INT node0, REF_INT node1, 
				      REF_INT new_node )
{
  REF_DBL log_m0[6], log_m1[6], m[6];
  REF_INT i;

  if ( !ref_node_valid(ref_node,node0) ||
       !ref_node_valid(ref_node,node1) ) 
    RSS( REF_INVALID, "node invalid" );

  for ( i = 0; i < 3 ; i++ )
    ref_node_xyz(ref_node,i,new_node) = 
      0.5 * (ref_node_xyz(ref_node,i,node0) + ref_node_xyz(ref_node,i,node1));

  for ( i = 0; i < ref_node_naux(ref_node) ; i++ )
    ref_node_aux(ref_node,i,new_node) = 
      0.5 * (ref_node_aux(ref_node,i,node0) + ref_node_aux(ref_node,i,node1));

  RSS( ref_matrix_log_m( ref_node_metric_ptr(ref_node,node0), log_m0 ),"log 0");
  RSS( ref_matrix_log_m( ref_node_metric_ptr(ref_node,node1), log_m1 ),"log 1");

  RSS( ref_matrix_average_m( log_m0, log_m1, m ),"log 1");
  
  RSS( ref_matrix_exp_m( m, ref_node_metric_ptr(ref_node,new_node) ),"exp m");

  return REF_SUCCESS;
}

REF_STATUS ref_node_resize_aux( REF_NODE ref_node )
{
  if ( NULL == ref_node->aux )
    {
      ref_malloc( ref_node->aux, 
		  ref_node_naux(ref_node)*ref_node_max(ref_node), REF_DBL );
    }
  else
    {
      ref_realloc( ref_node->aux, 
		   ref_node_naux(ref_node)*ref_node_max(ref_node), REF_DBL );
    }
  return REF_SUCCESS;
}
