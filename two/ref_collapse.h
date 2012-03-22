
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

REF_STATUS ref_collapse_edge_geometry( REF_GRID ref_grid, 
				       REF_INT node0, REF_INT node1,
				       REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_mixed( REF_GRID ref_grid, 
				    REF_INT node0, REF_INT node1,
				    REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_local_tets( REF_GRID ref_grid, 
					 REF_INT node0, REF_INT node1,
					 REF_BOOL *allowed );

REF_STATUS ref_collapse_edge_quality( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT node1,
				      REF_BOOL *allowed );

END_C_DECLORATION

#endif /* REF_COLLAPSE_H */
