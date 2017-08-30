
#ifndef REF_INFLATE_H
#define REF_INFLATE_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_node.h"
#include "ref_dict.h"

BEGIN_C_DECLORATION

REF_STATUS ref_inflate_pri_min_dot( REF_NODE ref_node, 
				    REF_INT *nodes,  
				    REF_DBL *min_dot );

REF_STATUS ref_inflate_face( REF_GRID ref_grid, 
			     REF_DICT faceids, 
			     REF_DBL *origin, 
			     REF_DBL thickness, REF_DBL xshift );

REF_STATUS ref_inflate_radially( REF_GRID ref_grid, 
				 REF_DICT faceids,
				 REF_DBL *origin, 
				 REF_DBL thickness, 
				 REF_DBL mach_angle_rad,
				 REF_DBL alpha_rad );

REF_STATUS ref_inflate_rate( REF_INT nlayers,
			     REF_DBL first_thickness,
			     REF_DBL total_thickness,
			     REF_DBL *rate );
REF_STATUS ref_inflate_total_thickness( REF_INT nlayers,
					REF_DBL first_thickness,
					REF_DBL rate,
					REF_DBL *total_thickness );
REF_STATUS ref_inflate_dthickness( REF_INT nlayers,
				   REF_DBL first_thickness,
				   REF_DBL rate,
				   REF_DBL *dHdr );

REF_STATUS ref_inflate_origin( REF_GRID ref_grid,
			       REF_DICT faceids,
			       REF_DBL *origin );

END_C_DECLORATION

#endif /* REF_INFLATE_H */
