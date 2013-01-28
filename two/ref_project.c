
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ref_project.h"

#include "ref_node.h"
#include "ref_cell.h"

#include "ref_math.h"

REF_STATUS ref_project_edge( REF_GRID ref_grid,
			     REF_INT node0, REF_INT node1,
			     REF_INT new_node )
{
  REF_NODE ref_node = ref_grid_node(ref_grid);

  printf("%d %d %d %d\n",node0,node1,new_node,ref_node_n(ref_node));

  return REF_SUCCESS;
}

