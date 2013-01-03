
#include <stdlib.h>
#include <stdio.h>

#include "ref_inflate.h"

#include "ref_cell.h"
#include "ref_malloc.h"
#include "ref_math.h"

REF_STATUS ref_inflate_pri_min_dot( REF_NODE ref_node, 
				    REF_INT *nodes,  
				    REF_DBL *min_dot )
{
  REF_INT tri_nodes[3];
  REF_DBL top_normal[3];
  REF_DBL bot_normal[3];
  REF_DBL edge[3];
  REF_INT node, i;

  tri_nodes[0]= nodes[0];
  tri_nodes[1]= nodes[1];
  tri_nodes[2]= nodes[2];
  RSS( ref_node_tri_normal( ref_node, tri_nodes, bot_normal ), "bot"); 
  RSS( ref_math_normalize( bot_normal ), "norm bot");

  tri_nodes[0]= nodes[3];
  tri_nodes[1]= nodes[4];
  tri_nodes[2]= nodes[5];
  RSS( ref_node_tri_normal( ref_node, tri_nodes, top_normal ), "top"); 
  RSS( ref_math_normalize( top_normal ), "norm top");

  *min_dot = 1.0;
  for ( node=0;node<3;node++)
    {
      for ( i = 0; i < 3 ; i++ )
	edge[i] = ref_node_xyz(ref_node,i,nodes[node+3]) 
	  - ref_node_xyz(ref_node,i,nodes[node]); 
      RSS( ref_math_normalize( edge ), "norm edge0");
      *min_dot = MIN(*min_dot,ref_math_dot(edge,bot_normal));
      *min_dot = MIN(*min_dot,ref_math_dot(edge,top_normal));
    }

  return REF_SUCCESS;
}

REF_STATUS ref_inflate_face( REF_GRID ref_grid, 
			     REF_DICT faceids, 
			     REF_DBL thickness, REF_DBL xshift )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_CELL qua = ref_grid_qua(ref_grid);
  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_INT cell, tri_side, node0, node1;
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT ntri, tris[2], nquad, quads[2];
  REF_INT tri_node;
  REF_INT *o2n;
  REF_INT global, new_node;
  REF_INT new_cell;
  REF_DBL min_dot;

  REF_DBL normal[3], ref_normal[3], dot, len;
  REF_INT ref_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT item, ref;
  REF_DBL radius, scale;

#define ref_inflate_hyperbolically 0
#define ref_inflate_parallel 1
#define ref_inflate_cylindrically 2

  REF_INT ref_inflation_style = ref_inflate_parallel;

  ref_malloc_init( o2n, ref_node_max(ref_node), 
		   REF_INT, REF_EMPTY );

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      for(tri_node=0;tri_node<3;tri_node++)
	{
	  node0 = nodes[tri_node];
	  if ( REF_EMPTY == o2n[node0] )
	    {
	      RSS( ref_node_next_global( ref_node, &global ), "global" );
	      RSS( ref_node_add( ref_node, global, &new_node ), "add node" );
	      o2n[node0] = new_node;
	      switch (ref_inflation_style) 
		{
		case ref_inflate_hyperbolically:
		case ref_inflate_parallel:
		  RAS( ref_node_valid(ref_node,node0),"inlvalid tri node");
		  normal[0]=0.0;
		  normal[1]=ref_node_xyz(ref_node,1,node0);
		  normal[2]=ref_node_xyz(ref_node,2,node0);
		  RSS( ref_math_normalize( normal ), "make norm" );
		  each_ref_cell_having_node( tri, node0, item, ref )
		    {
		      ref_cell_nodes( tri, ref, ref_nodes );
		      if ( !ref_dict_has_key( faceids, ref_nodes[3] ) )
			continue;
		      RSS( ref_node_tri_normal( ref_node, 
						ref_nodes, ref_normal ), "n" );
		      RSS( ref_math_normalize( ref_normal ), "make norm" );
		      dot = -ref_math_dot(normal, ref_normal);
		      if ( dot < 0.5 || dot > 2.0 ) printf("dot %f\n",dot);
		      normal[1] /= dot;
		      normal[2] /= dot;
		    }
		  len = sqrt( ref_math_dot(normal,normal) );
		  if ( ref_inflate_parallel ) len = 1.0;
		  ref_node_xyz(ref_node,0,new_node) = 
		    len*xshift + ref_node_xyz(ref_node,0,node0);
		  ref_node_xyz(ref_node,1,new_node) = 
		    thickness*normal[1] + ref_node_xyz(ref_node,1,node0);
		  ref_node_xyz(ref_node,2,new_node) = 
		    thickness*normal[2] + ref_node_xyz(ref_node,2,node0);
		  break;
		case ref_inflate_cylindrically:
		  radius = sqrt( ref_node_xyz(ref_node,1,node0) *
				 ref_node_xyz(ref_node,1,node0) +
				 ref_node_xyz(ref_node,2,node0) *
				 ref_node_xyz(ref_node,2,node0) );
		  if ( !ref_math_divisible(radius+thickness,radius) )
		    THROW("div zero");
		  scale = (radius+thickness)/radius;
		  ref_node_xyz(ref_node,0,new_node) = 
		    xshift + ref_node_xyz(ref_node,0,node0);
		  ref_node_xyz(ref_node,1,new_node) = 
		    scale * ref_node_xyz(ref_node,1,node0);
		  ref_node_xyz(ref_node,2,new_node) = 
		    scale * ref_node_xyz(ref_node,2,node0);
		  break;
		}
	    }
	}

  if ( ref_inflation_style != ref_inflate_parallel )
    RSS( ref_inflate_fix(ref_grid,faceids, o2n), "fix" );

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	for(tri_side=0;tri_side<3;tri_side++)
	  {
	    node0 = ref_cell_e2n(tri,0,tri_side,cell);
	    node1 = ref_cell_e2n(tri,1,tri_side,cell);
	    RSS( ref_cell_list_with( tri, 
				     node0, node1,
				     2, &ntri,
				     tris ),"bad tri count");
	    if ( 1 == ntri ) 
	      {
		RSS( ref_cell_list_with( qua, 
					 node0, node1,
					 2, &nquad,
					 quads ),"bad quad count");
		if ( 1 != nquad ) THROW("tri without quad");
		new_nodes[4] = ref_cell_c2n(qua,4,quads[0]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	    if (  ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[0]) ) &&
		 !ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[1]) ) )
	      {
		new_nodes[4] = ref_cell_c2n(tri,3,tris[1]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	    if ( !ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[0]) ) &&
		  ref_dict_has_key( faceids, ref_cell_c2n(tri,3,tris[1]) ) )
	      {
		new_nodes[4] = ref_cell_c2n(tri,3,tris[0]);
		new_nodes[0] = node0;
		new_nodes[1] = node1;
		new_nodes[2] = o2n[node1];
		new_nodes[3] = o2n[node0];
		RSS( ref_cell_add( qua, new_nodes, &new_cell ), "qua tri1");
		continue;
	      }
	  }
      }

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	new_nodes[0] = nodes[0];
	new_nodes[1] = nodes[2];
	new_nodes[2] = nodes[1];
	new_nodes[3] = o2n[nodes[0]];
	new_nodes[4] = o2n[nodes[2]];
	new_nodes[5] = o2n[nodes[1]];
	
	RSS( ref_inflate_pri_min_dot( ref_node, new_nodes, &min_dot ), "md");
	if ( min_dot <= 0.0 ) 
	  {
	    printf("min_dot %f\n",min_dot);
	    THROW("malformed prism");
	  }
	
	RSS( ref_cell_add( pri, new_nodes, &new_cell ), "pri");
      }

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	nodes[0] = o2n[nodes[0]];
	nodes[1] = o2n[nodes[1]];
	nodes[2] = o2n[nodes[2]];
	RSS( ref_cell_replace_whole( tri, cell, nodes ), "repl");
      }

  ref_free( o2n );

  return REF_SUCCESS;
}

REF_STATUS ref_inflate_fix( REF_GRID ref_grid, 
			    REF_DICT faceids, 
			    REF_INT *o2n )
{
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];

  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);

  REF_INT cell;

  REF_DBL theta[3];
  REF_DBL area0, area1;
  REF_INT i;
  REF_INT node, node0, node1;
  REF_DBL projection[3], edge[3];
  REF_DBL dot;

  each_ref_cell_valid_cell_with_nodes( tri, cell, nodes)
    if ( ref_dict_has_key( faceids, nodes[3] ) )
      {
	RSS( ref_node_tri_area( ref_node, nodes, &area0 ), "a0");
	new_nodes[0] = o2n[nodes[0]];
	new_nodes[1] = o2n[nodes[1]];
	new_nodes[2] = o2n[nodes[2]];
	RSS( ref_node_tri_area( ref_node, new_nodes, &area1 ), "a1");
	
	if ( area1 / area0 < 1.00 ) printf("ratio %f\n",area1/area0);

	if ( area1 / area0 > 0.95 ) continue;

	for(i=0;i<3;i++)
	  theta[i] = 180.0/ref_math_pi*atan2(ref_node_xyz(ref_node,1,nodes[i]),
					     ref_node_xyz(ref_node,2,nodes[i]));

	node = REF_EMPTY;
	if ( theta[1] < theta[0] && theta[0] < theta[2] ) node = 0;
	if ( theta[2] < theta[0] && theta[0] < theta[1] ) node = 0;

	if ( theta[0] < theta[1] && theta[1] < theta[2] ) node = 1;
	if ( theta[2] < theta[1] && theta[1] < theta[0] ) node = 1;

	if ( theta[0] < theta[2] && theta[2] < theta[1] ) node = 2;
	if ( theta[1] < theta[2] && theta[2] < theta[0] ) node = 2;

	RUS(REF_EMPTY,node,"middle node not found, symmetry plane included?");

	node0 = node + 1;
	if ( 2 < node0 ) node0 -= 3;
	node1 = node0 + 1;
	if ( 2 < node1 ) node1 -= 3;

	for ( i = 0; i < 3 ; i++ )
	  projection[i] = ref_node_xyz(ref_node,i,new_nodes[node]) 
	    - ref_node_xyz(ref_node,i,new_nodes[node0]); 

	for ( i = 0; i < 3 ; i++ )
	  edge[i] = ref_node_xyz(ref_node,i,new_nodes[node1]) 
	    - ref_node_xyz(ref_node,i,new_nodes[node0]); 

	RSS( ref_math_normalize( edge ), "edge");

	dot = ref_math_dot( edge, projection );

	for ( i = 0; i < 3 ; i++ )
	  projection[i] -= dot*edge[i];
  
	for ( i = 0; i < 3 ; i++ )
	  ref_node_xyz(ref_node,i,new_nodes[node]) 
	    += (area0/area1-1.0)*projection[i];

	RSS( ref_node_tri_area( ref_node, nodes, &area0 ), "a0");
	RSS( ref_node_tri_area( ref_node, new_nodes, &area1 ), "a1");

      }    

  return REF_SUCCESS;
}
