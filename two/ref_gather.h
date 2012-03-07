
#ifndef REF_GATHER_H
#define REF_GATHER_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

REF_STATUS ref_gather_b8_ugrid( REF_GRID ref_grid, char *filename );

REF_STATUS ref_gather_ncell( REF_NODE ref_node, REF_CELL ref_cell, 
			     REF_INT *ncell );

REF_STATUS ref_gather_node( REF_NODE ref_node, FILE *file );
REF_STATUS ref_gather_cell( REF_NODE ref_node, REF_CELL ref_cell, 
			    REF_BOOL faceid_insted_of_c2n, FILE *file );

END_C_DECLORATION

#endif /* REF_GATHER_H */
