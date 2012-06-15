#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



#include "ref_grid.h"
#include "ref_adj.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_matrix.h"

#include "ref_cell.h"
#include "ref_sort.h"

#include "ref_stitch.h"

#include "ref_import.c"
#include "ref_export.c"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

int main( int argc, char *argv[] )
{

  if (argc>1)
    {      
      REF_GRID ref_grid;
      printf("importing %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"imp");
      printf("complete.\n");

      RSS(ref_stitch_together( ref_grid ),"stitch");

      RSS(ref_grid_inspect( ref_grid ), "inspection");
      return 0;
    }

  return 0;
}
