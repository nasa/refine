
#ifndef REF_PROJECT_H
#define REF_PROJECT_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

REF_STATUS ref_project_edge( REF_GRID ref_grid,
			     REF_INT node0, REF_INT node1,
			     REF_INT new_node );

END_C_DECLORATION

#endif /* REF_PROJECT_H */
