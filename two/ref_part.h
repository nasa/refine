
#ifndef REF_PART_H
#define REF_PART_H

#include "ref_defs.h"

#include "ref_grid.h"
#include "ref_node.h"

BEGIN_C_DECLORATION

#define ref_part_first( total_things, total_parts, part ) \
  (MIN(((((total_things)-1)/(total_parts))+1)*(part),(total_things)))

#define ref_part_implicit( total_things, total_parts, thing ) \
  ((thing) / ((((total_things)-1)/(total_parts))+1) )

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid, char *filename );
REF_STATUS ref_part_b8_ugrid_cell( REF_CELL ref_cell, REF_INT ncell,
				   REF_NODE ref_node, REF_INT nnode,
				   FILE *file, 
				   long conn_offset,
				   long faceid_offset );
REF_STATUS ref_part_ghost_xyz( REF_GRID ref_grid );
REF_STATUS ref_part_ghost_int( REF_GRID ref_grid, REF_INT *scalar );

REF_STATUS ref_part_metric( REF_NODE ref_node, char *filename );
REF_STATUS ref_part_ratio( REF_NODE ref_node, REF_DBL *ratio, char *filename );

END_C_DECLORATION

#endif /* REF_PART_H */
