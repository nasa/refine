
#ifndef REF_VALIDATION_H
#define REF_VALIDATION_H

#include "ref_defs.h"

#include "ref_grid.h"

REF_STATUS ref_validation_all( REF_GRID ref_grid );
REF_STATUS ref_validation_hanging_node( REF_GRID ref_grid );
REF_STATUS ref_validation_cell_face( REF_GRID ref_grid );
REF_STATUS ref_validation_cell_node( REF_GRID ref_grid );
REF_STATUS ref_validation_cell_volume( REF_GRID ref_grid );

END_C_DECLORATION

#endif /* REF_VALIDATION_H */
