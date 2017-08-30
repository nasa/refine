#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_grid.h"
#include "ref_axi.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (3 != argc) 
    {
      printf("usage: %s input_grid.extension output_grid.extension\n",argv[0]);
      return 0;
    }

  printf("importing %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from msh");
  printf("complete.\n");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  RSS(ref_axi_wedge(ref_grid),"axi wedge");

  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("exporting %s\n",argv[2]);
  RSS(ref_export_by_extension( ref_grid, argv[2] ),"to file");
  printf("done.\n");

  return 0;
}
