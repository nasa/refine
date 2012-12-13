
#ifndef REF_INFLATE_H
#define REF_INFLATE_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_dict.h"

BEGIN_C_DECLORATION

REF_STATUS ref_inflate_face( REF_GRID ref_grid, 
			     REF_DICT faceids, 
			     REF_DBL thickness, REF_DBL xshift );

REF_STATUS ref_inflate_normal( REF_GRID ref_grid, 
			       REF_DICT faceids, 
			       REF_DBL thickness, REF_DBL xshift );

END_C_DECLORATION

#endif /* REF_INFLATE_H */
