
#ifndef REF_PART_H
#define REF_PART_H

#include "ref_defs.h"

#include "ref_grid.h"

BEGIN_C_DECLORATION

#define ref_part_first( total_things, total_parts, part ) \
  (MIN(((((total_things)-1)/(total_parts))+1)*(part),(total_things)))

#define ref_part_implicit( total_things, total_parts, thing ) \
  ((thing) / ((((total_things)-1)/(total_parts))+1) )

REF_STATUS ref_part_b8_ugrid( REF_GRID *ref_grid, char *filename );

END_C_DECLORATION

#endif /* REF_PART_H */
