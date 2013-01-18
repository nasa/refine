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


#include "ref_gen.h"

#include "ref_import.c"
#include "ref_export.c"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

#include "ref_fixture.h"

int main( int argc, char *argv[] )
{

  if (argc==2) 
    {
      REF_GRID ref_grid;

      RSS( ref_fixture_clock( &ref_grid ), "clock" );
      printf("export to %s\n",argv[1]);
      RSS(ref_export_by_extension( ref_grid, argv[1] ),"export");
      printf("complete.\n");

      RSS( ref_grid_free( ref_grid ), "free" );
      printf("done.\n");

      return 0;
    }

  return 0;
}
