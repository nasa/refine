
#include <stdlib.h>
#include <stdio.h>

#include "ref_geom.h"

#include "ref_dict.h"
#include "ref_cell.h"
#include "ref_edge.h"
#include "ref_malloc.h"
#include "ref_adapt.h"
#include "ref_matrix.h"

REF_STATUS ref_geom_from_egads( REF_GRID *ref_grid_ptr, char *filename )
{
  REF_GRID ref_grid;
  REF_NODE ref_node;

  printf("EGAGS project %s\n",filename);
  RSS( ref_grid_create( ref_grid_ptr ), "create grid");
  ref_grid = (*ref_grid_ptr);
  ref_node = ref_grid_node(ref_grid);

  return REF_SUCCESS;
}
