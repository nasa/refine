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


#include "ref_project.h"

#include "ref_import.c"
#include "ref_export.c"

#include "ref_dict.h"
#include "ref_mpi.h"
#include "ref_edge.h"

#include "ref_math.h"

#include "ref_fixture.h"

int main( int argc, char *argv[] )
{

  if (argc==2) 
    {
      REF_GRID ref_grid;

      printf("import from %s\n",argv[1]);
      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );
      printf(" complete.\n");

      printf("export ref_project_test_orig.tec\n");
      RSS(ref_export_by_extension( ref_grid, "ref_project_test_orig" ),"exp");
      printf(" complete.\n");

      RSS( ref_grid_free( ref_grid ), "free" );
      printf("done.\n");

      return 0;
    }

  return 0;
}
