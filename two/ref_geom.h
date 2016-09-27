
#ifndef REF_GEOM_H
#define REF_GEOM_H

#include "ref_defs.h"

BEGIN_C_DECLORATION

  END_C_DECLORATION

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_geom_egads_fixture( char *filename );

REF_STATUS ref_geom_brep_from_egads( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_geom_tetgen_volume( REF_GRID ref_grid );
REF_STATUS ref_geom_grid_from_egads( REF_GRID *ref_grid, char *filename );

END_C_DECLORATION

#endif /* REF_GEOM_H */
