
#ifndef REF_TWOD_H
#define REF_TWOD_H

#include "ref_defs.h"

#include "ref_cell.h"

BEGIN_C_DECLORATION

REF_STATUS ref_twod_opposite_node( REF_CELL pri,
				   REF_INT node, REF_INT *opposite);

REF_STATUS ref_twod_opposite_edge( REF_CELL pri,
				   REF_INT node0, REF_INT node1, 
				   REF_INT *node2, REF_INT *node3);

END_C_DECLORATION

#endif /* REF_TWOD_H */
