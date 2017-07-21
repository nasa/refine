
#ifndef REF_COLLAPSE_H
#define REF_COLLAPSE_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

/* node1 is removed */

REF_STATUS ref_collapse_pass( REF_GRID ref_grid );

REF_STATUS ref_collapse_to_remove_node1( REF_GRID ref_grid, 
					 REF_INT *node0, REF_INT node1 );

REF_STATUS ref_collapse_edge( REF_GRID ref_grid, 
			      REF_INT node0, REF_INT node1 );
/*                               keep node0,  remove node1 */

REF_STATUS ref_collapse_edge_geometry( REF_GRID ref_grid, 
				       REF_INT node0, REF_INT node1,
				       REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_same_normal( REF_GRID ref_grid, 
					  REF_INT node0, REF_INT node1,
					  REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_mixed( REF_GRID ref_grid, 
				    REF_INT node0, REF_INT node1,
				    REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_local_tets( REF_GRID ref_grid, 
					 REF_INT node0, REF_INT node1,
					 REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_cad_constrained( REF_GRID ref_grid, 
					      REF_INT node0, REF_INT node1,
					      REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_quality( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT node1,
				      REF_BOOL *allowed );

REF_STATUS ref_collapse_face( REF_GRID ref_grid,
			      REF_INT keep0, REF_INT remove0,
			      REF_INT keep1, REF_INT remove1);


REF_STATUS ref_collapse_face_local_pris( REF_GRID ref_grid, 
					 REF_INT keep, REF_INT remove,
					 REF_BOOL *allowed );

REF_STATUS ref_collapse_face_quality( REF_GRID ref_grid, 
				      REF_INT keep, REF_INT remove,
				      REF_BOOL *allowed );

REF_STATUS ref_collapse_face_outward_norm( REF_GRID ref_grid, 
					   REF_INT keep, REF_INT remove,
					   REF_BOOL *allowed );

REF_STATUS ref_collapse_face_geometry( REF_GRID ref_grid, 
				       REF_INT keep, REF_INT remove,
				       REF_BOOL *allowed );

REF_STATUS ref_collapse_face_same_tangent( REF_GRID ref_grid, 
					   REF_INT keep, REF_INT remove,
					   REF_BOOL *allowed );

REF_STATUS ref_collapse_twod_pass( REF_GRID ref_grid );

REF_STATUS ref_collapse_face_remove_node1( REF_GRID ref_grid, 
					   REF_INT *actual_node0, 
					   REF_INT node1 );

END_C_DECLORATION

#endif /* REF_COLLAPSE_H */
