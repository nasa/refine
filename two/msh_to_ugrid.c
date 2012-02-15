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

#include "ref_axi.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc<3) 
    {
      printf("usage: %s filename.msh filename.ugrid [-axi]\n",argv[0]);
      return 0;
    }

  printf("importing %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from msh");
  printf("complete.\n");
      
  RSS(ref_grid_inspect( ref_grid ), "inspection");

  if ( 4 == argc )
    {
      RSS(ref_axi_wedge(ref_grid),"axi wedge");
      RSS(ref_grid_inspect( ref_grid ), "inspection");
    }

  printf("exporting %s\n",argv[2]);
  RSS(ref_export_ugrid( ref_grid, argv[2] ),"to ugrid");
  printf("done.\n");

  return 0;
}
