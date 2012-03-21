
#ifndef REF_SPLIT_H
#define REF_SPLIT_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_split_pass( REF_GRID ref_grid );

REF_STATUS ref_split_edge( REF_GRID ref_grid, 
			   REF_INT node0, REF_INT node1,
			   REF_INT new_node );

REF_STATUS ref_split_edge_mixed( REF_GRID ref_grid, 
				 REF_INT node0, REF_INT node1,
				 REF_BOOL *allowed );

REF_STATUS ref_split_edge_local_tets( REF_GRID ref_grid, 
				      REF_INT node0, REF_INT node1,
				      REF_BOOL *allowed );

REF_STATUS ref_split_edge_quality( REF_GRID ref_grid, 
				   REF_INT node0, REF_INT node1,
				   REF_INT new_node,
				   REF_BOOL *allowed );

END_C_DECLORATION

#endif /* REF_SPLIT_H */
