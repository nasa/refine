
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_quality.h"

#include "ref_cell.h"
#include "ref_node.h"
#include "ref_face.h"
#include "ref_hexdiv.h"
#include "ref_export.h"

#include "ref_math.h"

REF_STATUS ref_quality_hex( REF_GRID ref_grid )
{
  REF_CELL ref_cell;

  REF_HEXDIV ref_hexdiv;

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

  REF_INT face_marks, hex_marks;

  RSS( ref_hexdiv_create( &ref_hexdiv, ref_grid ), "make hexdiv");

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
	    RSS( ref_hexdiv_mark_cell_edge_split( ref_hexdiv, 
						  cell, cell_edge ), 
		 "mark cell edge");

	  }
      }

  printf("marks %d\n",marks);

  RSS( ref_hexdiv_mark_n( ref_hexdiv, &face_marks, &hex_marks ), "count marks");
  printf("marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS(ref_export_vtk(viz, "pole.vtk"),"to vtk");

  ref_grid_node(viz) = orig_viz_node;/* replace before free */
  RXS( ref_grid_free(viz),REF_NULL,"free temp grid");

  RSS( ref_hexdiv_mark_relax( ref_hexdiv ), "relax" );

  RSS( ref_hexdiv_mark_n( ref_hexdiv, &face_marks, &hex_marks ), "count marks");
  printf("relaxed marked faces %d hexes %d\n",face_marks,hex_marks);

  RSS( ref_hexdiv_split( ref_hexdiv ), "split hex to prism" );

  return REF_SUCCESS;
}
