
#include <stdlib.h>
#include <stdio.h>
#include "ref_cell.h"

REF_STATUS ref_cell_create( REF_INT nodes, REF_CELL *ref_cell )
{
  *cell = (REF_CELL)malloc( sizeof(REF_CELL_STRUCT) );
  RNS(*cell,'malloc NULL');

  (*cell)->nodes = nodes;
  (*cell)->n = 0;
  (*cell)->max = 0;
  (*cell)->c2n = NULL;

  return REF_SUCCESS;
}
