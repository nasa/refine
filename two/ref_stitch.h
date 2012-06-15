
#ifndef REF_STITCH_H
#define REF_STITCH_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_stitch_together( REF_GRID ref_grid, 
				REF_INT tri_boundary, REF_INT qua_boundary );

END_C_DECLORATION

#endif /* REF_STITCH_H */
