
#ifndef REF_SMOOTH_H
#define REF_SMOOTH_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_smooth_tri_steepest_descent( REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_smooth_tri_quality_around( REF_GRID ref_grid,
                                          REF_INT node,
                                          REF_DBL *min_quality );

REF_STATUS ref_smooth_outward_norm( REF_GRID ref_grid, 
				    REF_INT node,
				    REF_BOOL *allowed );

REF_STATUS ref_smooth_tri_ideal( REF_GRID ref_grid,
                                 REF_INT node,
                                 REF_INT tri,
                                 REF_DBL *ideal_location );

REF_STATUS ref_smooth_tri_weighted_ideal( REF_GRID ref_grid,
					  REF_INT node,
					  REF_DBL *ideal_location );

REF_STATUS ref_smooth_tri_improve( REF_GRID ref_grid,
				   REF_INT node );

REF_STATUS ref_smooth_tri_pass( REF_GRID ref_grid );

REF_STATUS ref_smooth_tet_ideal( REF_GRID ref_grid,
                                 REF_INT node,
                                 REF_INT tet,
                                 REF_DBL *ideal_location );

REF_STATUS ref_smooth_tet_weighted_ideal( REF_GRID ref_grid,
					  REF_INT node,
					  REF_DBL *ideal_location );


END_C_DECLORATION

#endif /* REF_SMOOTH_H */
