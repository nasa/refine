
#ifndef REF_SMOOTH_H
#define REF_SMOOTH_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_smooth_twod( REF_GRID ref_grid, REF_INT node );

REF_STATUS ref_smooth_tri_quality_around( REF_GRID ref_grid, 
					  REF_INT node,
					  REF_DBL *min_quality );

END_C_DECLORATION

#endif /* REF_SMOOTH_H */
