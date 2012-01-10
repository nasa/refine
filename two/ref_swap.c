
#include <stdlib.h>
#include <stdio.h>

#include "ref_swap.h"

REF_STATUS ref_swap_remove_two_face_cell( REF_GRID ref_grid, REF_INT cell )
{
  SUPRESS_UNUSED_COMPILER_WARNING(*ref_grid);
  printf("cell %d\n",cell);
  return REF_INVALID;
}
