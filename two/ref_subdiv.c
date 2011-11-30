
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
  REF_INT group, cell, cell_edge;
  REF_CELL ref_cell;

  (*ref_subdiv_ptr) = NULL;
  (*ref_subdiv_ptr) = (REF_SUBDIV)malloc( sizeof(REF_SUBDIV_STRUCT) );
  RNS(*ref_subdiv_ptr,"malloc ref_subdiv NULL");

  ref_subdiv = *ref_subdiv_ptr;

  ref_subdiv_grid(ref_subdiv) = ref_grid;

  RSS( ref_grid_make_edges( ref_grid ), "edges");

  RSS( ref_adj_create( &(ref_subdiv_adj( ref_subdiv )) ), "create adj" );

  ref_subdiv->e2n = (REF_INT *)malloc( ref_grid_nedge(ref_grid) 
				       * 2 * sizeof(REF_INT));
  RNS(ref_subdiv->e2n,"malloc global NULL");

  for ( edge=0 ; edge < ref_grid_nedge(ref_grid) ; edge++ )
    {
      ref_subdiv_e2n( ref_subdiv, edge, 0 ) = REF_EMPTY;
      ref_subdiv_e2n( ref_subdiv, edge, 1 ) = REF_EMPTY;
    }

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
      each_ref_cell_cell_edge( ref_cell, cell_edge )
        {
	  edge = ref_cell_c2e(ref_cell,cell_edge,cell);
	  ref_subdiv_e2n( ref_subdiv, edge, 0 ) = 
	    ref_cell_e2n(ref_cell,0,cell,cell_edge);
	  ref_subdiv_e2n( ref_subdiv, edge, 1 ) = 
	    ref_cell_e2n(ref_cell,1,cell,cell_edge);
	}

  for ( edge=0 ; edge < ref_grid_nedge(ref_grid) ; edge++ )
    {
      RUS(REF_EMPTY,ref_subdiv_e2n( ref_subdiv, edge, 0 ),"edge n0 empty");
      RUS(REF_EMPTY,ref_subdiv_e2n( ref_subdiv, edge, 1 ),"edge n1 empty");
      RSS( ref_adj_add( ref_subdiv_adj( ref_subdiv ), 
			ref_subdiv_e2n( ref_subdiv, edge, 0 ), 
			edge ), "adj n0");
      RSS( ref_adj_add( ref_subdiv_adj( ref_subdiv ), 
			ref_subdiv_e2n( ref_subdiv, edge, 1 ), 
			edge ), "adj n1");
    }

  ref_subdiv->mark = (REF_INT *)malloc( ref_grid_nedge(ref_grid) 
					* sizeof(REF_INT));
  RNS(ref_subdiv->mark,"malloc mark NULL");

  for ( edge=0 ; edge < ref_grid_nedge(ref_grid) ; edge++ )
    ref_subdiv_mark( ref_subdiv, edge ) = 0;

  ref_subdiv->node = (REF_INT *)malloc( ref_grid_nedge(ref_grid) 
					* sizeof(REF_INT));
  RNS(ref_subdiv->node,"malloc node NULL");

  for ( edge=0 ; edge < ref_grid_nedge(ref_grid) ; edge++ )
    ref_subdiv_node( ref_subdiv, edge ) = REF_EMPTY;

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_free( REF_SUBDIV ref_subdiv )
{
  if ( NULL == (void *)ref_subdiv ) return REF_NULL;

  free( ref_subdiv->node );
  free( ref_subdiv->mark );
  free( ref_subdiv->e2n );
  RSS( ref_adj_free( ref_subdiv_adj( ref_subdiv ) ), "free adj" );

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


REF_STATUS ref_subdiv_edge_with( REF_SUBDIV ref_subdiv, 
				 REF_INT node0, REF_INT node1,
				 REF_INT *edge )
{
  REF_INT item, ref;
  REF_INT n0, n1;

  each_ref_adj_node_item_with_ref( ref_subdiv_adj(ref_subdiv), node0, item, ref)
    {
      n0 = ref_subdiv_e2n(ref_subdiv,ref,0);
      n1 = ref_subdiv_e2n(ref_subdiv,ref,1);
      if ( ( n0==node0 && n1==node1 ) ||
	   ( n0==node1 && n1==node0 ) )
	{
	  *edge=ref;
	  return REF_SUCCESS;
	}
    }
  return REF_FAILURE;
}


REF_STATUS ref_subdiv_mark_to_split( REF_SUBDIV ref_subdiv, 
				     REF_INT node0, REF_INT node1 )
{
  REF_INT edge;

  RSS( ref_subdiv_edge_with(ref_subdiv,node0, node1, &edge), "missing edge");

  ref_subdiv_mark(ref_subdiv,edge) = 1;

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
  REF_INT edge, global, node;

  ref_grid = ref_subdiv_grid(ref_subdiv);

  for ( edge = 0; edge < ref_grid_nedge(ref_grid) ; edge++ )
    {
      if ( ref_subdiv_mark( ref_subdiv, edge ) )
	{
	  RSS( ref_node_next_global( ref_grid_node(ref_grid), &global ),
	       "next global");
	  RSS( ref_node_add( ref_grid_node(ref_grid), global, &node), 
	       "add node");
	  ref_subdiv_node( ref_subdiv, edge ) = node;
	}
    }

  return REF_SUCCESS;
}

REF_STATUS ref_subdiv_split( REF_SUBDIV ref_subdiv )
{
  REF_GRID ref_grid;
  REF_INT group, cell;
  REF_CELL ref_cell;
  REF_CELL ref_cell_split;
  REF_INT map;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT new_nodes[REF_CELL_MAX_NODE_PER];
  REF_INT node, edge, new_cell;
  REF_INT *marked_for_removal;

  ref_grid = ref_subdiv_grid(ref_subdiv);

  each_ref_grid_ref_cell( ref_subdiv_grid(ref_subdiv), group, ref_cell )
    {
      marked_for_removal = 
	(REF_INT *)malloc(ref_cell_max(ref_cell)*sizeof(REF_INT));
      RNS(marked_for_removal,"malloc failed");
      for(cell=0;cell<ref_cell_max(ref_cell);cell++)
	marked_for_removal[cell]=0;
      RSS( ref_cell_create( ref_cell_node_per(ref_cell), &ref_cell_split), "");
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
	      RSS( ref_subdiv_edge_with(ref_subdiv,nodes[0],nodes[1], &edge),
		   "mis");
	      new_nodes[0] = ref_subdiv_node(ref_subdiv,edge);
	      RSS( ref_subdiv_edge_with(ref_subdiv,nodes[3],nodes[4], &edge),
		   "mis");
	      new_nodes[3] = ref_subdiv_node(ref_subdiv,edge);
	      RSS(ref_cell_add(ref_cell_split,new_nodes,&new_cell),"add");

	      for(node=0;node<ref_cell_node_per(ref_cell);node++)
		new_nodes[node] = nodes[node];
	      RSS( ref_subdiv_edge_with(ref_subdiv,nodes[0],nodes[1], &edge),
		   "mis");
	      new_nodes[1] = ref_subdiv_node(ref_subdiv,edge);
	      RSS( ref_subdiv_edge_with(ref_subdiv,nodes[3],nodes[4], &edge),
		   "mis");
	      new_nodes[4] = ref_subdiv_node(ref_subdiv,edge);
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
    }

  return REF_SUCCESS;
}

