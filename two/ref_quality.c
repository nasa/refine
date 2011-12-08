
#include <stdlib.h>
#include <stdio.h>

#include "ref_quality.h"

#include "ref_cell.h"

REF_STATUS ref_quality_hex( REF_GRID ref_grid )
{
  REF_CELL ref_cell;

  ref_cell = ref_grid_hex(ref_grid);

  return REF_SUCCESS;
}
