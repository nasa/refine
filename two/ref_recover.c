
#include <stdlib.h>
#include <stdio.h>

#include "ref_recover.h"
#include "ref_malloc.h"
#include "ref_math.h"

REF_STATUS ref_recover_create( REF_RECOVER *ref_recover_ptr, REF_GRID ref_grid )
{
  REF_RECOVER ref_recover;

  ref_malloc( *ref_recover_ptr, 1, REF_RECOVER_STRUCT );
  ref_recover = (*ref_recover_ptr);

  ref_recover_grid(ref_recover) = ref_grid;
  ref_recover_n(ref_recover) = 0;

  return REF_SUCCESS;
}

REF_STATUS ref_recover_free( REF_RECOVER ref_recover )
{
  if ( NULL == (void *)ref_recover ) return REF_NULL;
  ref_free( ref_recover );
  return REF_SUCCESS;
}

REF_STATUS ref_recover_enclosing_triangle( REF_RECOVER ref_recover, 
					   REF_INT node,
					   REF_INT *enclosing_cell,
					   REF_DBL *enclosing_bary)
{
  REF_GRID ref_grid = ref_recover_grid(ref_recover);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL ref_cell = ref_grid_tri(ref_grid);
  REF_INT cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT subnodes[REF_CELL_MAX_SIZE_PER];
  REF_BOOL twod;
  REF_DBL bary[3], best_min, bary_min, total;
  REF_INT i;

  *enclosing_cell = REF_EMPTY;
  enclosing_bary[0]=0.0;enclosing_bary[1]=0.0;enclosing_bary[2]=0.0;
  best_min = -10.0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, nodes)
    {
      RSS( ref_node_edge_twod( ref_node, nodes[0], nodes[1], &twod ), "2d");
      if (twod)
	{
	  for (i=0;i<3;i++)
	    {
	      subnodes[0]=nodes[0];subnodes[1]=nodes[1];subnodes[2]=nodes[2];
	      subnodes[i]=node;
	      RSS( ref_node_tri_y_projection( ref_node, subnodes, 
					      &(bary[i]) ),"bar");
	    }
	  total = bary[0]+bary[1]+bary[2];
	  if ( ref_math_divisible(bary[0],total) &&
	       ref_math_divisible(bary[1],total) &&
	       ref_math_divisible(bary[2],total) )
	    {
	      bary[0] /= total;
	      bary[1] /= total;
	      bary[2] /= total;
	    }
	  else
	    {
	      return REF_DIV_ZERO;
	    }
	  bary_min = MIN(MIN(bary[0],bary[1]),bary[2]);
	  if ( bary_min > best_min )
	    {
	      best_min = bary_min;
	      for (i=0;i<3;i++)
		enclosing_bary[i]=bary[i];
	      *enclosing_cell = cell;
	    }
	}
    }

  return (REF_EMPTY==(*enclosing_cell)?REF_FAILURE:REF_SUCCESS);
}

REF_STATUS ref_recover_insert_twod( REF_RECOVER ref_recover, REF_DBL *xz,
				    REF_INT *node_ptr )
{
  REF_GRID ref_grid = ref_recover_grid(ref_recover);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_INT global;
  REF_INT node0, node1;
  REF_INT tri0, tri1, pri;
  REF_DBL bary[3];
  REF_INT nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT face_nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT new, cell;
  REF_INT node0_offset, node1_offset;

  *node_ptr = REF_EMPTY;

  RSS( ref_node_next_global( ref_node, &global ), "next global");
  RSS( ref_node_add( ref_node, global, &node0 ), "add node");
  ref_node_xyz(ref_node,0,node0) = xz[0];
  ref_node_xyz(ref_node,1,node0) = 0.0;
  ref_node_xyz(ref_node,2,node0) = xz[1];

  RSS(ref_recover_enclosing_triangle(ref_recover,node0,&tri0,bary),"create");
  RSS(ref_cell_nodes(ref_grid_tri(ref_grid),tri0,nodes),"tri0 nodes");
  face_nodes[0]=nodes[0];face_nodes[1]=nodes[1];face_nodes[2]=nodes[2];
  face_nodes[3]=nodes[0];
  RSS(ref_cell_with_face(ref_grid_pri(ref_grid),face_nodes,&pri),"pri");

  RSS(ref_cell_nodes(ref_grid_pri(ref_grid),pri,nodes),"pri nodes");
  node0_offset = 3;
  if ( face_nodes[0]==nodes[0] ||
       face_nodes[0]==nodes[1] ||
       face_nodes[0]==nodes[2] ) node0_offset = 0;
  node1_offset = 3-node0_offset;

  RSS( ref_recover_opposite_node( ref_grid, face_nodes[0], 
				  &(face_nodes[0]) ),"n0");
  RSS( ref_recover_opposite_node( ref_grid, face_nodes[1], 
				  &(face_nodes[1]) ),"n1");
  RSS( ref_recover_opposite_node( ref_grid, face_nodes[2], 
				  &(face_nodes[2]) ),"n2");
  RSS( ref_recover_opposite_node( ref_grid, face_nodes[3], 
				  &(face_nodes[3]) ),"n3");
  RSS(ref_cell_with_face(ref_grid_tri(ref_grid),face_nodes,&tri1),"tri1");

  RSS( ref_node_next_global( ref_node, &global ), "next global");
  RSS( ref_node_add( ref_node, global, &node1 ), "add node");
  ref_node_xyz(ref_node,0,node1) = xz[0];
  ref_node_xyz(ref_node,1,node1) = 1.0;
  ref_node_xyz(ref_node,2,node1) = xz[1];

  for (new=0;new<3;new++)
    {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid),tri0,new_nodes),"tri0 nodes");
      new_nodes[new]=node0;
      RSS(ref_cell_add(ref_grid_tri(ref_grid),new_nodes,&cell),"add");
    }
  RSS(ref_cell_remove(ref_grid_tri(ref_grid),tri0),"rm");

  for (new=0;new<3;new++)
    {
      RSS(ref_cell_nodes(ref_grid_tri(ref_grid),tri1,new_nodes),"tri0 nodes");
      new_nodes[new]=node1;
      RSS(ref_cell_add(ref_grid_tri(ref_grid),new_nodes,&cell),"add");
    }
  RSS(ref_cell_remove(ref_grid_tri(ref_grid),tri1),"rm");

  for (new=0;new<3;new++)
    {
      RSS(ref_cell_nodes(ref_grid_pri(ref_grid),pri,new_nodes),"tri0 nodes");
      new_nodes[new+node0_offset]=node0;
      new_nodes[new+node1_offset]=node1;
      RSS(ref_cell_add(ref_grid_pri(ref_grid),new_nodes,&cell),"add");
    }
  RSS(ref_cell_remove(ref_grid_pri(ref_grid),pri),"rm");

  *node_ptr = node0;

  return REF_SUCCESS;
}

REF_STATUS ref_recover_opposite_node( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT *node1 )
{
  REF_CELL pri = ref_grid_pri(ref_grid);
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];

  *node1 = REF_EMPTY;

  each_ref_cell_having_node( pri, node0, item, cell )
    {
      RSS( ref_cell_nodes(pri, cell, nodes), "nodes" );

      if ( node0 == nodes[0]  )
	{ *node1 = nodes[3]; return REF_SUCCESS; }
      if ( node0 == nodes[1]  )
	{ *node1 = nodes[4]; return REF_SUCCESS; }
      if ( node0 == nodes[2]  )
	{ *node1 = nodes[5]; return REF_SUCCESS; }

      if ( node0 == nodes[3]  )
	{ *node1 = nodes[0]; return REF_SUCCESS; }
      if ( node0 == nodes[4]  )
	{ *node1 = nodes[1]; return REF_SUCCESS; }
      if ( node0 == nodes[5]  )
	{ *node1 = nodes[2]; return REF_SUCCESS; }
    }

   return REF_FAILURE;
}

REF_STATUS ref_recover_edge_twod( REF_RECOVER ref_recover, 
				  REF_INT node0, REF_INT node1 )
{
  REF_GRID ref_grid = ref_recover_grid(ref_recover);
  REF_BOOL has_side;
  
  RSS( ref_cell_has_side( ref_grid_tri(ref_grid), node0, node1, 
			  &has_side),"has");

  return (has_side?REF_SUCCESS:REF_FAILURE);
}

REF_STATUS ref_recover_first_block_twod( REF_RECOVER ref_recover, 
					 REF_INT node0, REF_INT node1,
					 REF_INT *block1, REF_INT *block2 )
{
  REF_GRID ref_grid = ref_recover_grid(ref_recover);
  REF_NODE ref_node = ref_grid_node(ref_grid);
  REF_CELL tri = ref_grid_tri(ref_grid);
  REF_INT item, cell, nodes[REF_CELL_MAX_SIZE_PER];
  REF_INT sub_nodes[REF_CELL_MAX_SIZE_PER];
  REF_DBL proj, proj1, proj2;
  REF_DBL best_proj, min_proj;
  REF_INT best_block1, best_block2;

  *block1 = REF_EMPTY;
  *block2 = REF_EMPTY;

  best_proj = -1.0;
  best_block1 = REF_EMPTY;
  best_block2 = REF_EMPTY;

  each_ref_cell_having_node( tri, node0, item, cell )
    {
      /* get the nodes rotated so node0 is nodes[0] */
      RSS( ref_cell_nodes(tri, cell, nodes), "nodes" );
      if ( node0 == nodes[1] )
	{
	  nodes[3] = nodes[0];
	  nodes[0] = nodes[1];
	  nodes[1] = nodes[2];
	  nodes[2] = nodes[3];
	}
      if ( node0 == nodes[2] )
	{
	  nodes[3] = nodes[0];
	  nodes[0] = nodes[2];
	  nodes[2] = nodes[1];
	  nodes[1] = nodes[3];
	}
      /*
              nodes[2] (block2)
             / |  
	node0- | -  -  - node1
             \ |
              nodes[1] (block1)
      */
      sub_nodes[0] = nodes[1];
      sub_nodes[1] = node1;
      sub_nodes[2] = node0;
      RSS( ref_node_tri_y_projection( ref_node, sub_nodes, &proj1 ),"p1");
      sub_nodes[0] = nodes[2];
      sub_nodes[1] = node0;
      sub_nodes[2] = node1;
      RSS( ref_node_tri_y_projection( ref_node, sub_nodes, &proj2 ),"p2");
      RSS( ref_node_tri_y_projection( ref_node, nodes, &proj ),"p");
      if ( !ref_math_divisible(proj1,proj) ||
	   !ref_math_divisible(proj1,proj) ) return REF_DIV_ZERO;
      proj1 /= proj;
      proj2 /= proj;
      min_proj = MIN(proj1, proj2);
      if ( min_proj > best_proj )
	{ 
	  best_proj = min_proj;
	  best_block1 = nodes[1];
	  best_block2 = nodes[2];
	}
    }

  if ( best_proj < 1.0e-14 )
    {
      printf("%s: %d: %s: best proj %e\n",__FILE__,__LINE__,__func__,best_proj);
      return REF_FAILURE;
    }

  *block1 = best_block1;
  *block2 = best_block2;

  return REF_SUCCESS;
}
