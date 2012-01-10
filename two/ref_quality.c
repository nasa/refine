
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_quality.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_face.h"
#include "ref_shard.h"
#include "ref_subdiv.h"
#include "ref_export.h"

#include "ref_math.h"

REF_STATUS ref_quality_hex( REF_GRID ref_grid )
{
  REF_CELL ref_cell;

  REF_SHARD ref_shard;

  REF_NODE ref_node;
  REF_INT cell, cell_edge;
  REF_INT cell_nodes[8];
  REF_INT face_nodes[4];
  REF_INT node,face0,face1;
  REF_DBL normal0[3], normal1[3];
  REF_DBL angle;

  REF_GRID viz;

  REF_INT marks, new_cell;
  REF_CELL marked_cell;

  REF_INT face_marks, hex_marks;

  RSS( ref_shard_create( &ref_shard, ref_grid ), "make shard");

  ref_cell = ref_grid_hex(ref_grid);
  ref_node = ref_grid_node(ref_grid);

  RSS( ref_grid_empty_cell_clone(&viz,ref_grid),"temp grid for marks");
  marked_cell =  ref_grid_hex(viz);
  marks = 0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, cell_nodes)
    each_ref_cell_cell_edge( ref_cell, cell_edge )      
      {
	RSS( ref_cell_gen_edge_face( ref_cell, cell_edge, 
				     &face0, &face1 ), "edge faces" );
	
	for(node=0;node<4;node++)
	  face_nodes[node] = ref_cell_f2n(ref_cell,node,cell,face0);
	RSS( ref_face_normal( &ref_node_xyz(ref_node,0,face_nodes[0]),
			      &ref_node_xyz(ref_node,0,face_nodes[1]),
			      &ref_node_xyz(ref_node,0,face_nodes[2]),
			      &ref_node_xyz(ref_node,0,face_nodes[3]),
			      normal0 ), "normal" );
	for(node=0;node<4;node++)
	  face_nodes[node] = ref_cell_f2n(ref_cell,node,cell,face1);
	RSS( ref_face_normal( &ref_node_xyz(ref_node,0,face_nodes[0]),
			      &ref_node_xyz(ref_node,0,face_nodes[1]),
			      &ref_node_xyz(ref_node,0,face_nodes[2]),
			      &ref_node_xyz(ref_node,0,face_nodes[3]),
			      normal1 ), "normal" );
			     
	RSS( ref_math_normalize( normal0 ), "norm0");
	RSS( ref_math_normalize( normal1 ), "norm1");
	angle = ref_math_dot(normal0,normal1);
	if ( ABS(angle) > (1-1.0e-8) ) 
	  {
	    marks++;
	    RSS( ref_cell_add( marked_cell, cell_nodes, &new_cell ), 
		 "add marked");
	    RSB( ref_shard_mark_cell_edge_split( ref_shard, 
						  cell, cell_edge ), 
		 "mark cell edge",
		 ref_node_location(ref_node,face_nodes[0]););

	  }
      }

  printf("marks %d\n",marks);

  RSS( ref_shard_mark_n( ref_shard, &face_marks, &hex_marks ), "count marks");
  printf("marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS(ref_export_vtk(viz, "pole.vtk"),"to vtk");

  RSS( ref_grid_free_cell_clone(viz),"free temp grid");

  RSS( ref_shard_mark_relax( ref_shard ), "relax" );

  RSS( ref_shard_mark_n( ref_shard, &face_marks, &hex_marks ), "count marks");
  printf("relaxed marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS( ref_shard_split( ref_shard ), "split hex to prism" );

  return REF_SUCCESS;
}

REF_STATUS ref_quality_multiple_face_cell( REF_GRID ref_grid )
{
  REF_CELL ref_cell;
  REF_INT group, cell, cell_face;
  REF_INT node;
  REF_INT nodes[REF_CELL_MAX_NODE_PER];
  REF_INT face_nodes[REF_CELL_MAX_NODE_PER];
  REF_BOOL problem;
  REF_INT boundary_faces, found;
  REF_INT bcface[REF_CELL_MAX_FACE_PER];
  REF_INT targeted[5], marks;

  REF_GRID viz;
  REF_INT viz_cell;
  REF_INT edge, face0, face1;

  REF_SUBDIV ref_subdiv;

  RSS( ref_subdiv_create( &ref_subdiv, ref_grid ), "make subdiv");

  RSS( ref_grid_empty_cell_clone( &viz, ref_grid ), "viz grid" );
  
  problem = REF_FALSE;
  for ( boundary_faces = 0 ; boundary_faces<= 4; boundary_faces++)
    targeted[boundary_faces] = 0;

  each_ref_grid_ref_cell( ref_grid, group, ref_cell )
    each_ref_cell_valid_cell( ref_cell, cell )
    {
      boundary_faces = 0;
      each_ref_cell_cell_face( ref_cell, cell_face )
        {
	  bcface[cell_face] = REF_EMPTY;
	  for(node=0;node<4;node++)
	    nodes[node]=ref_cell_f2n(ref_cell,node,cell,cell_face);
	  
	  if ( nodes[0] == nodes[3] )
	    {
	      if ( REF_SUCCESS == ref_cell_with( ref_grid_tri( ref_grid ), 
						 nodes, &found ) )
		{
		  RSS( ref_cell_nodes( ref_grid_tri( ref_grid ), 
				       found, face_nodes), "tri");
		  bcface[cell_face] = face_nodes[3];
		  boundary_faces++;
		}
	    }
	  else
	    {
	      if ( REF_SUCCESS == ref_cell_with( ref_grid_qua( ref_grid ), 
						 nodes, &found ) )
		{
		  RSS( ref_cell_nodes( ref_grid_qua( ref_grid ), 
				       found, face_nodes), "qua");
		  bcface[cell_face] = face_nodes[4];
		  boundary_faces++;
		}
	    }
	}
      targeted[boundary_faces]++;
      if ( boundary_faces > 1 )
	{
	  problem = REF_TRUE;
	  RSS( ref_cell_nodes( ref_cell, cell, nodes), "cell nodes");
	  RSS( ref_cell_add( ref_grid_cell(viz,group), nodes, &viz_cell), 
	       "add viz cell");	  
	}
      if ( 4 == ref_cell_node_per(ref_cell) && 2 == boundary_faces )
	{
	  for(edge=0;edge<ref_cell_edge_per(ref_cell);edge++)
	    {
	      face0 = ref_cell_e2n_gen(ref_cell,0,edge);
	      face1 = ref_cell_e2n_gen(ref_cell,1,edge);
	      if ( REF_EMPTY != bcface[face0] &&
		   REF_EMPTY != bcface[face1] &&
		   bcface[face0] != bcface[face1] )
		{
		  ref_subdiv_mark( ref_subdiv, 
				   ref_cell_c2e(ref_cell,edge,cell) ) = 1; 
		}
	      if ( REF_EMPTY != bcface[face0] &&
		   REF_EMPTY != bcface[face1] &&
		   bcface[face0] == bcface[face1] )
		{
		  printf("same bcs %d %d\n",bcface[face0],bcface[face1]);
		}
	    }
	}
    }

  for ( boundary_faces = 0 ; boundary_faces <= 4; boundary_faces++ )
    printf(" %d : %d\n", boundary_faces, targeted[boundary_faces]);

  if ( problem )
    {
      ref_export_tec( viz, "ref_quality_multiple_face_cell.tec" );
    }

  RSS( ref_grid_free_cell_clone( viz ), "free viz");

  RSS( ref_subdiv_mark_n( ref_subdiv, &marks ), "n mark");
  printf("original marks %d\n",marks);

  RSS( ref_subdiv_mark_relax( ref_subdiv ), "relax subdiv");

  RSS( ref_subdiv_mark_n( ref_subdiv, &marks ), "n mark");
  printf("relaxed marks %d\n",marks);

  RSS( ref_subdiv_new_node( ref_subdiv ), "new node subdiv");
  RSS( ref_subdiv_split( ref_subdiv ), "split subdiv");
  RSS( ref_subdiv_free( ref_subdiv ), "free subdiv");

  return REF_SUCCESS;
}
