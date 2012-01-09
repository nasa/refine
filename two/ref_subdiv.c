
#include <stdlib.h>
#include <stdio.h>

#include "ref_subdiv.h"

static int ref_subdiv_map( REF_SUBDIV ref_subdiv, 
			   REF_CELL ref_cell, REF_INT cell )
{
  REF_INT edge, map, bit;

  map = 0;
  bit = 1;
  for ( edge = 0; edge < ref_cell_edge_per(ref_cell) ; edge++ )
    {
      map += bit*ref_subdiv_mark(ref_subdiv,ref_cell_c2e(ref_cell, edge, cell));
      bit *= 2;
    }

  return map;
}

REF_STATUS ref_subdiv_create( REF_SUBDIV *ref_subdiv_ptr, REF_GRID ref_grid )
{
  REF_SUBDIV ref_subdiv;
  REF_INT edge;

  (*ref_subdiv_ptr) = NULL;
  (*ref_subdiv_ptr) = (REF_SUBDIV)malloc( sizeof(REF_SUBDIV_STRUCT) );
  RNS(*ref_subdiv_ptr,"malloc ref_subdiv NULL");

  ref_subdiv = *ref_subdiv_ptr;

  ref_subdiv_grid(ref_subdiv) = ref_grid;

  RSS( ref_edge_create( &(ref_subdiv_edge( ref_subdiv )), 
			ref_subdiv_grid(ref_subdiv) ), "create edge" );

  ref_subdiv->mark = (REF_INT *)malloc( ref_edge_n(ref_subdiv_edge(ref_subdiv)) 
					* sizeof(REF_INT));
  RNS(ref_subdiv->mark,"malloc mark NULL");

  for ( edge=0 ; edge < ref_edge_n(ref_subdiv_edge(ref_subdiv)) ; edge++ )
    ref_subdiv_mark( ref_subdiv, edge ) = 0;

  ref_subdiv->node = (REF_INT *)malloc( ref_edge_n(ref_subdiv_edge(ref_subdiv))
					* sizeof(REF_INT));
  RNS(ref_subdiv->node,"malloc node NULL");

  for ( edge=0 ; edge < ref_edge_n(ref_subdiv_edge(ref_subdiv)) ; edge++ )
    ref_subdiv_node( ref_subdiv, edge ) = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv )
{
  if ( NULL == (void *)ref_subdiv ) return REF_NULL;

  free( ref_subdiv->node );
  free( ref_subdiv->mark );
  RSS( ref_edge_free( ref_subdiv_edge( ref_subdiv ) ), "free edge" );

  ref_cond_free( ref_subdiv );

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

REF_STATUS ref_subdiv_node_between( REF_SUBDIV ref_subdiv, 
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

#define edge_or(ce0,ce1) \
  { \
    REF_INT ge0, ge1; \
    ge0 = ref_cell_c2e( ref_cell, ce0, cell );\
    ge1 = ref_cell_c2e( ref_cell, ce1, cell );\
    if ( ref_subdiv_mark( ref_subdiv, ge0 ) != \
	 ref_subdiv_mark( ref_subdiv, ge1 ) ) \
      {\
	again = REF_TRUE;\
	ref_subdiv_mark( ref_subdiv, ge0 ) = 1;\
	ref_subdiv_mark( ref_subdiv, ge1 ) = 1;\
      }\
  }

#define promote_2_3(ce0,ce1,ce2)			\
  { \
    REF_INT ge0, ge1, ge2, sum;			\
    ge0 = ref_cell_c2e( ref_cell, ce0, cell );\
    ge1 = ref_cell_c2e( ref_cell, ce1, cell );\
    ge2 = ref_cell_c2e( ref_cell, ce2, cell );\
    sum = ref_subdiv_mark( ref_subdiv, ge0 ) \
        + ref_subdiv_mark( ref_subdiv, ge1 ) \
        + ref_subdiv_mark( ref_subdiv, ge2 );    \
    if ( 2 == sum ) \
      {\
	again = REF_TRUE;\
	ref_subdiv_mark( ref_subdiv, ge0 ) = 1;\
	ref_subdiv_mark( ref_subdiv, ge1 ) = 1;\
	ref_subdiv_mark( ref_subdiv, ge2 ) = 1;\
      }\
  }

REF_STATUS ref_subdiv_mark_relax( REF_SUBDIV ref_subdiv )
{
  REF_INT group, cell;
  REF_CELL ref_cell;
  REF_BOOL again;

  again = REF_TRUE;

  while (again)
    {

      again = REF_FALSE;

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
    }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_new_node( REF_SUBDIV ref_subdiv )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;
  REF_EDGE ref_edge;
  REF_INT edge, global, node, node0, node1, ixyz;

  ref_grid = ref_subdiv_grid(ref_subdiv);
  ref_node = ref_grid_node(ref_grid);
  ref_edge = ref_subdiv_edge(ref_subdiv);

  for ( edge = 0; edge < ref_edge_n(ref_edge) ; edge++ )
    {
      if ( ref_subdiv_mark( ref_subdiv, edge ) )
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

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split_qua( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT new_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, new_cell;
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
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[3], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[0], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
    RSS(ref_cell_add(ref_cell,nodes,&new_cell),"add");
      
  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split_tri( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT new_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, new_cell;
  REF_INT *marked_for_removal;

  REF_INT edge01, edge12, edge20;

  ref_cell = ref_grid_tri(ref_subdiv_grid(ref_subdiv));
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
			  nodes[2], nodes[0], &edge20 ), "e20" );

      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge20 ) )
	RSS( REF_IMPLEMENT, "code" );
      if( ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge20 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge01 ) )
	RSS( REF_IMPLEMENT, "code" );
      if( ref_subdiv_mark( ref_subdiv, edge20 ) &&
	  ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  !ref_subdiv_mark( ref_subdiv, edge12 ) )
	RSS( REF_IMPLEMENT, "code" );
	  
      if( ref_subdiv_mark( ref_subdiv, edge01 ) &&
	  ref_subdiv_mark( ref_subdiv, edge12 ) &&
	  ref_subdiv_mark( ref_subdiv, edge20 ) )
	{
	  marked_for_removal[cell]=1;
	  /* near node 0 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 1 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 2 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* center */
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  continue;
	}

      if( ref_subdiv_mark( ref_subdiv, edge01 ) )
	{
	  marked_for_removal[cell]=1;
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	}
	  
      if( ref_subdiv_mark( ref_subdiv, edge12 ) )
	{
	  marked_for_removal[cell]=1;
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	}
	  
      if( ref_subdiv_mark( ref_subdiv, edge20 ) )
	{
	  marked_for_removal[cell]=1;
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	}
	  
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");
      
  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_cell_add(ref_cell,nodes,&new_cell),"add");
      
  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split_pri( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT new_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, new_cell;
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
	  
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[3])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[3],nodes[4], 
				       &(new_nodes[4])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  break;
	case 459: /* prism split */
	  marked_for_removal[cell]=1;
	  
	  /* near edge 0-3 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
	  printf("cell %d, map %d\n",cell,map);
	  RSS( REF_IMPLEMENT, "map not implemented yet" )
	}
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");

  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_cell_add(ref_cell,nodes,&new_cell),"add");

  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

#define node_swap(nodes,a,b) \
  {REF_INT nst;nst=(nodes)[(a)];(nodes)[(a)]=(nodes)[(b)];(nodes)[(b)]=nst;}

REF_STATUS ref_subdiv_split_tet( REF_SUBDIV ref_subdiv )
{
  REF_INT cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT new_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, new_cell;
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
	  global_edge = ref_cell_c2e(ref_cell, split_edge, cell);
	  
	  marked_for_removal[cell]=1;
	  
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  new_nodes[ref_cell_e2n_gen(ref_cell,0,split_edge)] = 
	    ref_subdiv_node(ref_subdiv, global_edge);
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
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
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 1 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	  /* near node 2 */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[1], 
				       &(new_nodes[1])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");


	  /* center */
	  for(node=0;node<ref_cell_node_per(ref_cell);node++)
	    new_nodes[node] = nodes[node];
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[0],nodes[1], 
				       &(new_nodes[0])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[1],nodes[2], 
				       &(new_nodes[1])), "mis");
	  RSS( ref_subdiv_node_between(ref_subdiv,nodes[2],nodes[0], 
				       &(new_nodes[2])), "mis");
	  RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");
	  
	  break;
	default:
	  printf("cell %d, map %d\n",cell,map);
	  RSS( REF_IMPLEMENT, "map not implemented yet" )
	}
    }

  for(cell=0;cell<ref_cell_max(ref_cell);cell++)
    if ( 1 == marked_for_removal[cell] )
      RSS(ref_cell_remove(ref_cell,cell),"remove");

  each_ref_cell_valid_cell_with_nodes( ref_cell_split, cell, nodes)
    RSS(ref_cell_add(ref_cell,nodes,&new_cell),"add");

  RSS( ref_cell_free( ref_cell_split ), "temp ref_cell free");
  free(marked_for_removal);

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split( REF_SUBDIV ref_subdiv )
{

  RSS( ref_subdiv_split_tet( ref_subdiv ), "split tet" );
  RSS( ref_subdiv_split_pri( ref_subdiv ), "split pri" );

  RSS( ref_subdiv_split_qua( ref_subdiv ), "split qua" );
  RSS( ref_subdiv_split_tri( ref_subdiv ), "split tri" );

  return REF_SUCCESS;
}

