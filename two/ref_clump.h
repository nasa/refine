
#ifndef REF_CLUMP_H
#define REF_CLUMP_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

  END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_clump_around( REF_GRID ref_grid, REF_INT node,
			     char *filename );
REF_STATUS ref_clump_tri_around( REF_GRID ref_grid, REF_INT node,
				 char *filename );
REF_STATUS ref_clump_stuck_edges( REF_GRID ref_grid );
REF_STATUS ref_clump_stuck_edges_twod( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_CLUMP_H */
