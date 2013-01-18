
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_gen.h"

#include "ref_node.h"

REF_STATUS ref_gen_fill( REF_GRID ref_grid )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);

  printf("nnodes %d",ref_node_n(ref_node));
  
  return REF_SUCCESS;
}

