
#ifndef REF_FIXTURE_H
#define REF_FIXTURE_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_fixture_tet_grid( REF_GRID *ref_grid );

REF_STATUS ref_fixture_pyr_grid( REF_GRID *ref_grid );

REF_STATUS ref_fixture_pri_grid( REF_GRID *ref_grid );
REF_STATUS ref_fixture_pri_tet_cap_grid( REF_GRID *ref_grid_ptr );
REF_STATUS ref_fixture_pri_stack_grid( REF_GRID *ref_grid );

REF_STATUS ref_fixture_hex_grid( REF_GRID *ref_grid );

REF_STATUS ref_fixture_hex_brick_grid( REF_GRID *ref_grid );
REF_STATUS ref_fixture_pri_brick_grid( REF_GRID *ref_grid );

END_C_DECLORATION

#endif /* REF_FIXTURE_H */
