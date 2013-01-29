#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_sort.h"
#include "ref_dict.h"
#include "ref_list.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (3 != argc) 
    {
      printf("usage: %s input_grid.extension output_grid.extension\n",argv[0]);
      return 0;
    }

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "import" );
  RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
