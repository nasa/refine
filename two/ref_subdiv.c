
#include <stdlib.h>
#include <stdio.h>

#include "ref_subdiv.h"

#include "ref_malloc.h"
#include "ref_mpi.h"
#include "ref_adj.h"

static REF_INT ref_subdiv_map( REF_SUBDIV ref_subdiv, 
			       REF_CELL ref_cell, REF_INT cell )
{
  REF_INT edge, map, bit;

  map = 0;
  bit = 1;
  for ( edge = 0; edge < ref_cell_edge_per(ref_cell) ; edge++ )
    {
      map += bit*ref_subdiv_mark(ref_subdiv,ref_cell_c2e(ref_cell, edge, cell));
      bit *= 2;
      /*
      printf("edge %d bit %d mark %d map %d\n",edge,bit,
	     ref_subdiv_mark(ref_subdiv,ref_cell_c2e(ref_cell, edge, cell)),
	     map);
      */
    }

  return map;
}

static REF_STATUS ref_subdiv_map_to_edge( REF_INT map )
{
  REF_INT edge, bit;
  
  bit = 2048;
  for ( edge = 11; edge >= 0; edge--)
    {
      if ( map >= bit )
	{ 
	  map -= bit;
	  printf("edge %d bit %d\n",edge,bit);
	}
      bit /= 2;
    }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_create( REF_SUBDIV *ref_subdiv_ptr, REF_GRID ref_grid )
{
  REF_SUBDIV ref_subdiv;

  ref_malloc( *ref_subdiv_ptr, 1, REF_SUBDIV_STRUCT );

  ref_subdiv = *ref_subdiv_ptr;

  ref_subdiv_grid(ref_subdiv) = ref_grid;

  RSS( ref_edge_create( &(ref_subdiv_edge( ref_subdiv )), 
			ref_subdiv_grid(ref_subdiv) ), "create edge" );

  ref_malloc_init( ref_subdiv->mark, ref_edge_n(ref_subdiv_edge(ref_subdiv)),
		   REF_INT, 0 );
  ref_malloc_init( ref_subdiv->node, ref_edge_n(ref_subdiv_edge(ref_subdiv)),
		   REF_INT, REF_EMPTY );

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv )
{
  if ( NULL == (void *)ref_subdiv ) return REF_NULL;

  ref_free( ref_subdiv->node );
  ref_free( ref_subdiv->mark );
  RSS( ref_edge_free( ref_subdiv_edge( ref_subdiv ) ), "free edge" );

  ref_free( ref_subdiv );

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_inspect( REF_SUBDIV ref_subdiv )
{
  REF_INT group, cell, cell_edge, edge, map;
  REF_CELL ref_cell;

  each_ref_grid_ref_cell( ref_subdiv_grid(ref_subdiv), group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      {
	map = ref_subdiv_map( ref_subdiv, ref_cell, cell );
	printf(" group %d cell %d map %d\n",group,cell, map);
	each_ref_cell_cell_edge( ref_cell, cell_edge )
	  {
	    edge = ref_cell_c2e(ref_cell,cell_edge,cell);
	    printf("  edge %d nodes %d %d mark %d\n",
		   edge,
		   ref_cell_e2n(ref_cell,0,cell,cell_edge),
		   ref_cell_e2n(ref_cell,1,cell,cell_edge),
		   ref_subdiv_mark(ref_subdiv,edge) );
	  }

      }
  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_n( REF_SUBDIV ref_subdiv, REF_INT *n )
{
  REF_INT edge;

  *n = 0;
  for ( edge=0 ; edge < ref_edge_n(ref_subdiv_edge(ref_subdiv)) ; edge++ )
    if ( 0 != ref_subdiv_mark( ref_subdiv, edge ) ) (*n)++;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_mark_to_split( REF_SUBDIV ref_subdiv, 
				     REF_INT node0, REF_INT node1 )
{
  REF_INT edge;

  RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ), 
		      node0, node1,
		      &edge ), "missing edge");

  ref_subdiv_mark(ref_subdiv,edge) = 1;

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_node_between( REF_SUBDIV ref_subdiv, 
					   REF_INT node0, REF_INT node1,
					   REF_INT *new_node )
{
  REF_INT edge;

  (*new_node) = REF_EMPTY;

  RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ), 
		      node0, node1,
		      &edge ), "missing edge");

  (*new_node) = ref_subdiv_node(ref_subdiv,edge);

  return REF_SUCCESS;
}

#define edge_or(ce0,ce1)			\
  {						\
    REF_INT ge0, ge1;				\
    ge0 = ref_cell_c2e( ref_cell, ce0, cell );	\
    ge1 = ref_cell_c2e( ref_cell, ce1, cell );	\
    if ( ref_subdiv_mark( ref_subdiv, ge0 ) !=	\
	 ref_subdiv_mark( ref_subdiv, ge1 ) )	\
      {						\
	again = REF_TRUE;			\
	ref_subdiv_mark( ref_subdiv, ge0 ) = 1;	\
	ref_subdiv_mark( ref_subdiv, ge1 ) = 1;	\
      }						\
  }

#define promote_2_3(ce0,ce1,ce2)		\
  {						\
    REF_INT ge0, ge1, ge2, sum;			\
    ge0 = ref_cell_c2e( ref_cell, ce0, cell );	\
    ge1 = ref_cell_c2e( ref_cell, ce1, cell );	\
    ge2 = ref_cell_c2e( ref_cell, ce2, cell );	\
    sum = ref_subdiv_mark( ref_subdiv, ge0 )	\
        + ref_subdiv_mark( ref_subdiv, ge1 )	\
        + ref_subdiv_mark( ref_subdiv, ge2 );	\
    if ( 2 == sum )				\
      {						\
	again = REF_TRUE;			\
	ref_subdiv_mark( ref_subdiv, ge0 ) = 1;	\
	ref_subdiv_mark( ref_subdiv, ge1 ) = 1;	\
	ref_subdiv_mark( ref_subdiv, ge2 ) = 1;	\
      }						\
  }

#define promote_2_all()							\
  {									\
    REF_INT ge0, ge1, ge2, ge3, ge4, ge5, sum;				\
    ge0 = ref_cell_c2e( ref_cell, 0, cell );				\
    ge1 = ref_cell_c2e( ref_cell, 1, cell );				\
    ge2 = ref_cell_c2e( ref_cell, 2, cell );				\
    ge3 = ref_cell_c2e( ref_cell, 3, cell );				\
    ge4 = ref_cell_c2e( ref_cell, 4, cell );				\
    ge5 = ref_cell_c2e( ref_cell, 5, cell );				\
    sum = ref_subdiv_mark( ref_subdiv, ge0 )				\
        + ref_subdiv_mark( ref_subdiv, ge1 )				\
        + ref_subdiv_mark( ref_subdiv, ge2 )				\
        + ref_subdiv_mark( ref_subdiv, ge3 )				\
        + ref_subdiv_mark( ref_subdiv, ge4 )				\
        + ref_subdiv_mark( ref_subdiv, ge5 );				\
    if ( 2 == sum )							\
      if ( ( ge0 > 0 && ge5 > 0) ||					\
	   ( ge1 > 0 && ge4 > 0) ||					\
	   ( ge2 > 0 && ge1 > 3) )					\
	{								\
	  again = REF_TRUE;						\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 0, cell ) ) = 1;	\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 1, cell ) ) = 1;	\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 2, cell ) ) = 1;	\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 3, cell ) ) = 1;	\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 4, cell ) ) = 1;	\
	  ref_subdiv_mark( ref_subdiv,					\
			   ref_cell_c2e( ref_cell, 5, cell ) ) = 1;	\
	}								\
  }

static REF_STATUS ref_subdiv_mark_relax( REF_SUBDIV ref_subdiv )
{
  REF_INT group, cell;
  REF_CELL ref_cell;
  REF_BOOL again;

  again = REF_TRUE;
  while (again)
    {
      again = REF_FALSE;

      RSS( ref_edge_ghost_int( ref_subdiv_edge(ref_subdiv),
			       ref_subdiv->mark), "ghost mark" );

      each_ref_grid_ref_cell( ref_subdiv_grid(ref_subdiv), group, ref_cell )
	each_ref_cell_valid_cell( ref_cell, cell )
	  {
	    switch ( ref_cell_node_per(ref_cell) )
	      {
	      case 4:
		promote_2_3(3,4,5);
		promote_2_3(1,2,5);
		promote_2_3(0,2,4);
		promote_2_3(0,1,3);
		promote_2_all();
		break;
	      case 5:
		edge_or(0,7);
		edge_or(1,5);
		edge_or(3,6);
		promote_2_3(0,1,3);
		promote_2_3(5,6,7);
		break;
	      case 6:
		edge_or(0,6);
		edge_or(3,8);
		edge_or(1,7);
		promote_2_3(0,1,3);
		promote_2_3(6,7,8);
		break;
	      default:
		RSS(REF_IMPLEMENT,"implement cell type");
		break;    
	      }
	  }

      RSS( ref_mpi_all_or( &again ), "mpi all or" );
    }

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_new_node( REF_SUBDIV ref_subdiv )
{
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_EDGE ref_edge = ref_subdiv_edge(ref_subdiv);
  REF_INT edge, global, node, node0, node1, ixyz;
  REF_INT part;

  REF_INT *edge_global, *edge_part;
  REF_DBL *edge_xyz;

  RSS( ref_node_synchronize_globals( ref_node ), "sync glob" );

  for ( edge = 0; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      if ( ref_subdiv_mark( ref_subdiv, edge ) )
	{
	  RSS( ref_edge_part(ref_edge,edge,&part), "edge part" );
	  if ( ref_mpi_id == part )
	    {
	      RSS( ref_node_next_global( ref_node, &global ),
		   "next global");
	      RSS( ref_node_add( ref_node, global, &node), 
		   "add node");
	      ref_subdiv_node( ref_subdiv, edge ) = node;

	      node0 = ref_edge_e2n(ref_edge, edge, 0 );
	      node1 = ref_edge_e2n(ref_edge, edge, 1 );
	      for (ixyz=0;ixyz<3;ixyz++)
		ref_node_xyz(ref_node,ixyz,node) = 
		  0.5 * ( ref_node_xyz(ref_node,ixyz,node0) +
			  ref_node_xyz(ref_node,ixyz,node1) );
	    }
	}
    }

  RSS( ref_node_shift_new_globals( ref_node ), "shift glob" );
  
  ref_malloc_init( edge_global, ref_edge_n(ref_edge), REF_INT, REF_EMPTY );
  ref_malloc_init( edge_part,   ref_edge_n(ref_edge), REF_INT, REF_EMPTY );
  ref_malloc_init( edge_xyz,  3*ref_edge_n(ref_edge), REF_DBL, -999.0 );
  
  for ( edge = 0; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      node = ref_subdiv_node( ref_subdiv, edge );
      if ( REF_EMPTY != node )
	{
	  edge_global[edge] = ref_node_global(ref_node,node);
	  edge_part[edge] = ref_node_part(ref_node,node);
	  for (ixyz=0;ixyz<3;ixyz++)
	    edge_xyz[ixyz+3*edge] = ref_node_xyz(ref_node,ixyz,node);
	}
    }

  RSS( ref_edge_ghost_int( ref_edge, edge_global ), "global ghost" );
  RSS( ref_edge_ghost_int( ref_edge, edge_part ), "part ghost" );
  RSS( ref_edge_ghost_dbl( ref_edge, edge_xyz, 3 ), "xyz ghost" );

  for ( edge = 0; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      node = ref_subdiv_node( ref_subdiv, edge );
      global = edge_global[edge];
      if ( REF_EMPTY == node && REF_EMPTY != global )
	{
	  RSS( ref_node_add( ref_node, global, &node), 
	       "add node");
	  ref_subdiv_node( ref_subdiv, edge ) = node;
	  ref_node_part(ref_node,node) = edge_part[edge];
	  for (ixyz=0;ixyz<3;ixyz++)
	    ref_node_xyz(ref_node,ixyz,node) = edge_xyz[ixyz+3*edge];

	}
    }

  ref_free( edge_xyz );
  ref_free( edge_part );
  ref_free( edge_global );

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_add_local_cell( REF_SUBDIV ref_subdiv, 
					     REF_CELL ref_cell, 
					     REF_INT *nodes)
{
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_BOOL has_local;
  REF_INT node, new_cell;

  has_local = REF_FALSE;

  for ( node=0; node<ref_cell_node_per(ref_cell); node++ )
    has_local = has_local || 
      ( ref_mpi_id == ref_node_part(ref_node,nodes[node]) );
  
  if ( has_local )
    RSS(ref_cell_add(ref_cell,nodes,&new_cell),"add");
	
  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_qua( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01, edge12, edge23, edge30;

  ref_cell = ref_grid_qua(ref_subdiv_grid(ref_subdiv));
  marked_for_removal = 
    (REF_INT *)malloc(ref_cell_max(ref_cell)*sizeof(REF_INT));
  RNS(marked_for_removal,"malloc failed");
  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    marked_for_removal[cell]=0;
  RSS( ref_cell_create( &ref_cell_split, 
			ref_cell_node_per(ref_cell), 
			ref_cell_last_node_is_an_id(ref_cell)), 
       "temp cell");
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes");
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[0], nodes[1], &edge01 ), "e01" );
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[1], nodes[2], &edge12 ), "e12" );
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[2], nodes[3], &edge23 ), "e23" );
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[3], nodes[0], &edge30 ), "e30" );
      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge23 ) &&
	  ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge30 ) )
	RSS( REF_IMPLEMENT, "all quad edges" );
      if( ref_subdiv_mark( ref_subdiv, edge01 ) !=
	  ref_subdiv_mark( ref_subdiv, edge23 ) ||
	  ref_subdiv_mark( ref_subdiv, edge12 ) !=
	  ref_subdiv_mark( ref_subdiv, edge30 ) )
	RSS( REF_IMPLEMENT, "quad edges not paired" );
	  
      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge23 ) )
	{
	  marked_for_removal[cell]=1;

	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[3], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[3], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	}
	  
      if( ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge30 ) )
	{
	  marked_for_removal[cell]=1;

	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[0], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	}
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");
      
  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes),"add local");
      
  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_tri( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL tri, qua;
  REF_CELL tri_split, qua_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01, edge12, edge20;

  tri = ref_grid_tri(ref_subdiv_grid(ref_subdiv));
  marked_for_removal = 
    (REF_INT *)malloc(ref_cell_max(tri)*sizeof(REF_INT));
  RNS(marked_for_removal,"malloc failed");
  for(cell=0;cell<ref_cell_max(tri);cell++)
    marked_for_removal[cell]=0;

  RSS( ref_cell_create( &tri_split, 
			ref_cell_node_per(tri), 
			ref_cell_last_node_is_an_id(tri)), 
       "temp tri");

  qua = ref_grid_qua(ref_subdiv_grid(ref_subdiv));
  RSS( ref_cell_create( &qua_split, 
			ref_cell_node_per(qua), 
			ref_cell_last_node_is_an_id(qua)), 
       "temp qua");

  each_ref_cell_valid_cell( tri, cell )
    {
      RSS( ref_cell_nodes( tri, cell, nodes ), "nodes");
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[0], nodes[1], &edge01 ), "e01" );
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[1], nodes[2], &edge12 ), "e12" );
      RSS( ref_edge_with( ref_subdiv_edge( ref_subdiv ),
			  nodes[2], nodes[0], &edge20 ), "e20" );
	  
      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge20 ) )
	{
	  marked_for_removal[cell]=1;
	  /* near node 0 */
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  /* near node 1 */
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  /* near node 2 */
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  /* center */
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	  continue;
	}

      /*
	2   3-2
	|\  | |
	0-1 0-1
       */

      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge20 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  new_nodes[4] =  new_nodes[3]; /* faceid */
	  new_nodes[3] =  new_nodes[2]; /* last node */
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(qua_split,new_nodes,&new_cell),"add");
	  continue;
	}

      if( ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge20 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge01 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  new_nodes[4] =  new_nodes[3]; /* faceid */
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(qua_split,new_nodes,&new_cell),"add");
	  continue;
	}

      if( ref_subdiv_mark( ref_subdiv, edge20 ) &&
	  ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge12 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  new_nodes[4] =  new_nodes[3]; /* faceid */
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(qua_split,new_nodes,&new_cell),"add");
	  continue;
	}


      if( ref_subdiv_mark( ref_subdiv, edge01 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	}
	  
      if( ref_subdiv_mark( ref_subdiv, edge12 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	}
	  
      if( ref_subdiv_mark( ref_subdiv, edge20 ) )
	{
	  marked_for_removal[cell]=1;
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	  RSS( ref_cell_nodes( tri, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(tri_split,new_nodes,&new_cell),"add");
	}
	  
    }

  for(cell=0;cell<ref_cell_max(tri);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(tri,cell),"remove");
      
  each_ref_cell_valid_cell_with_nodes( tri_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, tri, nodes),"add local");

  each_ref_cell_valid_cell_with_nodes( qua_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, qua, nodes),"add local");
      
  RSS( ref_cell_free( tri_split ), "temp tri free");
  RSS( ref_cell_free( qua_split ), "temp qua free");

  free(marked_for_removal);

  return REF_SUCCESS;
}

static REF_STATUS ref_subdiv_split_pri( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;

  ref_cell = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  marked_for_removal = 
    (REF_INT *)malloc(ref_cell_max(ref_cell)*sizeof(REF_INT));
  RNS(marked_for_removal,"malloc failed");
  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    marked_for_removal[cell]=0;
  RSS( ref_cell_create( &ref_cell_split, 
			ref_cell_node_per(ref_cell), 
			ref_cell_last_node_is_an_id(ref_cell)), 
       "temp cell");
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      map = ref_subdiv_map( ref_subdiv, ref_cell, cell );
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes");
      switch ( map )
	{
	case 0: /* don't split */
	  break;
	case 65: /* prism split */
	  marked_for_removal[cell]=1;
	  
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  break;
	case 459: /* prism split */
	  marked_for_removal[cell]=1;
	  
	  /* near edge 0-3 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[5], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near edge 1-4 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[3], 
				       &(new_nodes[3])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[5], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near edge 2-5 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[5],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[5],nodes[3], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* center */
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");

	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[5], 
				       &(new_nodes[4])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[5],nodes[3], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");


	  break;
	default:
	  RSS( ref_subdiv_map_to_edge( map ), "map2edge");
	  printf("pri %d, map %d\n",cell,map);
	  RSS( REF_IMPLEMENT, "map not implemented yet" )
	}
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");

  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes),"add local");

  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

#define node_swap(nodes,a,b) \
  {REF_INT nst;nst=(nodes)[(a)];(nodes)[(a)]=(nodes)[(b)];(nodes)[(b)]=nst;}

#define add_cell_with(fnnw0, fnnw1, fnnw2, fnnw3) \
  new_nodes[0] = (fnnw0); new_nodes[1] = (fnnw1); \
  new_nodes[2] = (fnnw2); new_nodes[3] = (fnnw3); \
  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

static REF_STATUS ref_subdiv_split_tet( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;

  REF_INT edge,split_edge, global_edge;

  ref_cell = ref_grid_tet(ref_subdiv_grid(ref_subdiv));
  marked_for_removal = 
    (REF_INT *)malloc(ref_cell_max(ref_cell)*sizeof(REF_INT));
  RNS(marked_for_removal,"malloc failed");
  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    marked_for_removal[cell]=0;
  RSS( ref_cell_create( &ref_cell_split, 
			ref_cell_node_per(ref_cell), 
			ref_cell_last_node_is_an_id(ref_cell)), 
       "temp cell");
  each_ref_cell_valid_cell( ref_cell, cell )
    {
      map = ref_subdiv_map( ref_subdiv, ref_cell, cell );
      RSS( ref_cell_nodes( ref_cell, cell, nodes ), "nodes");
      switch ( map )
	{
	case 0: /* don't split */
	  break;
	case  1: 
	case  2: 
	case  4:
	case  8:
	case 16:
	case 32:
	  split_edge=REF_EMPTY;
	  for ( edge = 0; edge < ref_cell_edge_per(ref_cell) ; edge++ )
	    if (ref_subdiv_mark(ref_subdiv,
				ref_cell_c2e(ref_cell, edge, cell)))
		split_edge = edge;
	  RAS( REF_EMPTY != split_edge, "edge not found");
	  global_edge = ref_cell_c2e(ref_cell, split_edge, cell);
	  
	  marked_for_removal[cell]=1;
	  
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  new_nodes[ref_cell_e2n_gen(ref_cell,0,split_edge)] = 
	    ref_subdiv_node(ref_subdiv, global_edge);
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  new_nodes[ref_cell_e2n_gen(ref_cell,1,split_edge)] = 
	    ref_subdiv_node(ref_subdiv, global_edge);
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  break;
	case 11:
	case 56:
	case 38:
	case 21:
	  /* orient cell for other cases */

	  marked_for_removal[cell]=1;

	  if ( 56 == map ) { node_swap(nodes,0,3); node_swap(nodes,1,2); }
	  if ( 38 == map ) { node_swap(nodes,1,3); node_swap(nodes,0,2); }
	  if ( 21 == map ) { node_swap(nodes,2,3); node_swap(nodes,0,1); }

	  /* near node 0 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 1 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 2 */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");


	  /* center */
	  RSS( ref_cell_nodes( ref_cell, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  
	  break;
	  /*
                              inode3------5------inode2
                                 / \              . /
                                /   \          .   /
                               /     \      .     /
                              /       \  .       /
                             /        .\        /
                            2      1    4      3
                           /    .        \    /
                          /  .            \  /
                         /.                \/
                      inode0------0------inode1
	   */
	case 63:
	  marked_for_removal[cell]=1;
	  { 
	    REF_INT n0, n1, n2, n3;
	    REF_INT e0, e1, e2, e3,e4, e5;
	    n0 = nodes[0]; n1 = nodes[1]; n2 = nodes[2]; n3 = nodes[3];
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1],&e0),"e");
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2],&e1),"e");
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[3],&e2),"e");
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2],&e3),"e");
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[3],&e4),"e");
	    RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[3],&e5),"e");
	    add_cell_with(e0, e2, e1, n0);
	    add_cell_with(e0, e3, e4, n1);
	    add_cell_with(e1, e5, e3, n2);
	    add_cell_with(e2, e4, e5, n3);
	    add_cell_with(e0, e5, e1, e2);
	    add_cell_with(e0, e5, e2, e4);
	    add_cell_with(e0, e5, e4, e3);
	    add_cell_with(e0, e5, e3, e1);
	  }
	  break;
	default:
	  RSS( ref_subdiv_map_to_edge( map ), "map2edge");
	  printf("tet %d, map %d\n",cell,map);
	  RSS( REF_IMPLEMENT, "map not implemented yet" )
	}
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");

  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, ref_cell, nodes),"add local");

  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}
static REF_STATUS ref_subdiv_split_pyr( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL pyr, pri;
  REF_CELL pyr_split, pri_split;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_cell;
  REF_INT *marked_for_removal;

  REF_INT map;

  pyr = ref_grid_pyr(ref_subdiv_grid(ref_subdiv));
  marked_for_removal = 
    (REF_INT *)malloc(ref_cell_max(pyr)*sizeof(REF_INT));
  RNS(marked_for_removal,"malloc failed");
  for(cell=0;cell<ref_cell_max(pyr);cell++)
    marked_for_removal[cell]=0;

  RSS( ref_cell_create( &pyr_split, 
			ref_cell_node_per(pyr), 
			ref_cell_last_node_is_an_id(pyr)), 
       "temp pyr");

  pri = ref_grid_pri(ref_subdiv_grid(ref_subdiv));
  RSS( ref_cell_create( &pri_split, 
			ref_cell_node_per(pri), 
			ref_cell_last_node_is_an_id(pri)), 
       "temp pri");

  each_ref_cell_valid_cell( pyr, cell )
    {
      map = ref_subdiv_map( ref_subdiv, pyr, cell );
      RSS( ref_cell_nodes( pyr, cell, nodes ), "nodes");
      switch ( map )
	{
	case 0: /* don't split */
	  break;
	case 129: /* split into two pyr*/
	  marked_for_removal[cell]=1;
	  
	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(pyr_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS(ref_cell_add(pyr_split,new_nodes,&new_cell),"add");
	  break;
	case 72: /* split into pyr and pri edge 3-6 */
	  marked_for_removal[cell]=1;
	  
	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[2], 
				       &(new_nodes[4])), "mis");
	  RSS(ref_cell_add(pyr_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[2], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(pri_split,new_nodes,&new_cell),"add");
	  break;
	case 34: /* split into pyr and pri edge 1-5 */
	  marked_for_removal[cell]=1;
	  
	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[2], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(pyr_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[2], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(pri_split,new_nodes,&new_cell),"add");
	  break;
	case 235: /* split into 1 pyr, 2 pri*/
	  marked_for_removal[cell]=1;

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[2], 
				       &(new_nodes[3])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[2], 
				       &(new_nodes[4])), "mis");
	  RSS(ref_cell_add(pyr_split,new_nodes,&new_cell),"add");

	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[2], 
				       &(new_nodes[4])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[3], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(pri_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[2], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(pri_split,new_nodes,&new_cell),"add");

	  RSS( ref_cell_nodes( pyr, cell, new_nodes ), "nodes");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[4],nodes[2], 
				       &(new_nodes[5])), "mis");
	  RSS(ref_cell_add(pri_split,new_nodes,&new_cell),"add");
	  
	  break;
	default:
	  RSS( ref_subdiv_map_to_edge( map ), "map2edge");
	  printf("pyr %d, map %d\n",cell,map);
	  RSS( REF_IMPLEMENT, "map not implemented yet" )
	}
    }

  for(cell=0;cell<ref_cell_max(pyr);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(pyr,cell),"remove");

  each_ref_cell_valid_cell_with_nodes( pyr_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, pyr, nodes),"add local");

  each_ref_cell_valid_cell_with_nodes( pri_split, cell, nodes)
    RSS(ref_subdiv_add_local_cell(ref_subdiv, pri, nodes),"add local");

  RSS( ref_cell_free( pyr_split ), "temp pyr free");
  RSS( ref_cell_free( pri_split ), "temp pri free");

  free(marked_for_removal);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split( REF_SUBDIV ref_subdiv )
{
  REF_GRID ref_grid = ref_subdiv_grid(ref_subdiv);
  REF_NODE ref_node = ref_grid_node(ref_subdiv_grid(ref_subdiv));
  REF_INT node;

  RSS(ref_subdiv_mark_relax(ref_subdiv),"relax marks");
  RSS(ref_subdiv_new_node(ref_subdiv),"new nodes");

  RSS( ref_subdiv_split_tet( ref_subdiv ), "split tet" );
  RSS( ref_subdiv_split_pri( ref_subdiv ), "split pri" );
  /* pyr comes last, it can make other elements too */
  RSS( ref_subdiv_split_pyr( ref_subdiv ), "split pyr" );

  RSS( ref_subdiv_split_qua( ref_subdiv ), "split qua" );
  /* tri comes last, it can make qua elements too */
  RSS( ref_subdiv_split_tri( ref_subdiv ), "split tri" );

  /* remove unused nods on partition boundaries */
  each_ref_node_valid_node( ref_node, node )
    if ( ref_adj_empty( ref_cell_adj(ref_grid_tet(ref_grid)), node) &&
	 ref_adj_empty( ref_cell_adj(ref_grid_pyr(ref_grid)), node) &&
	 ref_adj_empty( ref_cell_adj(ref_grid_pri(ref_grid)), node) &&
	 ref_adj_empty( ref_cell_adj(ref_grid_hex(ref_grid)), node) )
      {
	if ( ref_mpi_id == ref_node_part(ref_node,node) )
	  RSS( REF_FAILURE, "unused local node");
	if ( !ref_adj_empty( ref_cell_adj(ref_grid_tri(ref_grid)), node) ||
	     !ref_adj_empty( ref_cell_adj(ref_grid_qua(ref_grid)), node) )
	  RSS( REF_FAILURE, "boundary face node not in vol cells");
	RSS( ref_node_remove_without_global( ref_node, node ), "rm");
      }

  return REF_SUCCESS;
}

