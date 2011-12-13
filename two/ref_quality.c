
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_quality.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_grid_export.h"

REF_STATUS ref_tri_normal(REF_DBL *xyz0, REF_DBL *xyz1, REF_DBL *xyz2, 
			  REF_DBL *normal )
{
  REF_DBL edge1[3], edge2[3];

  edge1[0] = xyz1[0] - xyz0[0];
  edge1[1] = xyz1[1] - xyz0[1];
  edge1[2] = xyz1[2] - xyz0[2];

  edge2[0] = xyz2[0] - xyz0[0];
  edge2[1] = xyz2[1] - xyz0[1];
  edge2[2] = xyz2[2] - xyz0[2];

  normal[0] = edge1[1]*edge2[2] - edge1[2]*edge2[1];
  normal[1] = edge1[2]*edge2[0] - edge1[0]*edge2[2];
  normal[2] = edge1[0]*edge2[1] - edge1[1]*edge2[0];

  return REF_SUCCESS;
}

REF_STATUS ref_quad_normal( REF_DBL *xyz0, REF_DBL *xyz1, 
			    REF_DBL *xyz2, REF_DBL *xyz3, 
			    REF_DBL *normal )
{
  REF_DBL normal0[3], normal1[3];

  RSS( ref_tri_normal( xyz0, xyz1, xyz2, normal0 ), "normal0" );
  RSS( ref_tri_normal( xyz2, xyz3, xyz0, normal1 ), "normal1" );

  normal[0] = normal0[0] + normal1[0];
  normal[1] = normal0[1] + normal1[1];
  normal[2] = normal0[2] + normal1[2];

  return REF_SUCCESS;
}

#define ref_dot(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1]+(a)[2]*(b)[2])

REF_STATUS ref_normalize( REF_DBL *normal )
{
  REF_DBL length;

  length = ref_dot(normal,normal);
  if (ABS(length) < 1.0e-15) return REF_DIV_ZERO;

  length = sqrt(length);
  
  normal[0] /= length;
  normal[1] /= length;
  normal[2] /= length;

  length = ref_dot(normal,normal);
  RAS( (ABS(length-1.0) < 1.0e-13), "vector length not unity"); 

  return REF_SUCCESS;
}

REF_STATUS ref_quality_hex( REF_GRID ref_grid )
{
  REF_CELL ref_cell;

  REF_NODE ref_node;
  REF_INT cell, cell_edge;
  REF_INT cell_nodes[8];
  REF_INT face_nodes[4];
  REF_INT node,face0,face1;
  REF_DBL normal0[3], normal1[3];
  REF_DBL angle;

  REF_GRID viz;
  REF_NODE orig_viz_node;

  REF_INT marks, new_cell;
  REF_CELL marked_cell;


  ref_cell = ref_grid_hex(ref_grid);
  ref_node = ref_grid_node(ref_grid);

  RSS( ref_grid_create(&viz),"temp grid for marks");
  orig_viz_node = ref_grid_node(viz); /*save before replacement */
  ref_grid_node(viz) = ref_grid_node(ref_grid);
  marked_cell =  ref_grid_hex(viz);
  marks = 0;

  each_ref_cell_valid_cell_with_nodes( ref_cell, cell, cell_nodes)
    each_ref_cell_cell_edge( ref_cell, cell_edge )      
      {
	RSS( ref_cell_gen_edge_face( ref_cell, cell_edge, 
				     &face0, &face1 ), "edge faces" );
	
	for(node=0;node<4;node++)
	  face_nodes[node] = ref_cell_f2n(ref_cell,node,cell,face0);
	RSS( ref_quad_normal( &ref_node_xyz(ref_node,0,face_nodes[0]),
			      &ref_node_xyz(ref_node,0,face_nodes[1]),
			      &ref_node_xyz(ref_node,0,face_nodes[2]),
			      &ref_node_xyz(ref_node,0,face_nodes[3]),
			      normal0 ), "normal" );
	for(node=0;node<4;node++)
	  face_nodes[node] = ref_cell_f2n(ref_cell,node,cell,face1);
	RSS( ref_quad_normal( &ref_node_xyz(ref_node,0,face_nodes[0]),
			      &ref_node_xyz(ref_node,0,face_nodes[1]),
			      &ref_node_xyz(ref_node,0,face_nodes[2]),
			      &ref_node_xyz(ref_node,0,face_nodes[3]),
			      normal1 ), "normal" );
			     
	RSS( ref_normalize( normal0 ), "norm0");
	RSS( ref_normalize( normal1 ), "norm1");
	angle = ref_dot(normal0,normal1);
	if ( ABS(angle) > (1-1.0e-8) ) 
	  {
	    marks++;
	    printf("angle %d, cell %d, edge %d : %f\n",
		   marks,cell,cell_edge,angle);
	    RSS( ref_cell_add( marked_cell, cell_nodes, &new_cell ), 
		 "add marked");
	  }
      }

  RSS(ref_grid_export_vtk(viz, "pole.vtk"),"to vtk");

  ref_grid_node(viz) = orig_viz_node;/* replace before free */
  RXS( ref_grid_free(viz),REF_NULL,"free temp grid");

  return REF_SUCCESS;
}
