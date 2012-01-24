#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_import.h"
#include "ref_export.h"

#include "ref_adj.h"
#include "ref_grid.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

#include "ref_edge.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_dict.h"

#include "ref_quality.h"
#include "ref_shard.h"
#include "ref_subdiv.h"

#include "ref_validation.h"

#include "ref_math.h"
#include "ref_swap.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc<2) 
    {
      printf("usage: %s filename.extension\n",argv[0]);
      return 0;
    }

  printf("reading %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"by extension");
  printf("complete.\n");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("hex quality.\n");
  RSS(ref_quality_hex( ref_grid ),"quality");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("ugrid.\n");
  RSS(ref_export_b8_ugrid(ref_grid, "ref2.b8.ugrid"),"to ugrid");

  printf("free.\n");
  RSS(ref_grid_free(ref_grid),"free");

  printf("done.\n");
  return 0;
}
