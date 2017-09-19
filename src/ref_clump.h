
#ifndef REF_CLUMP_H
#define REF_CLUMP_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

  END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_clump_around( REF_GRID ref_grid, REF_INT node,
                             const char *filename );
REF_STATUS ref_clump_between( REF_GRID ref_grid, REF_INT node0, REF_INT node1,
			      const char *filename );
REF_STATUS ref_clump_tri_around( REF_GRID ref_grid, REF_INT node,
                                 const char *filename );
REF_STATUS ref_clump_stuck_edges( REF_GRID ref_grid, REF_DBL ratio_tol );
REF_STATUS ref_clump_stuck_edges_twod( REF_GRID ref_grid );

REF_STATUS ref_clump_tet_quality( REF_GRID ref_grid, REF_DBL min_quality,
                                  const char *filename );

END_C_DECLORATION

#endif /* REF_CLUMP_H */
