#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_part.h"
#include "ref_grid.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_BOOL valid_inputs;

  valid_inputs = (3 == argc);

  if ( !valid_inputs )
    {
      printf("usage: %s project.extension project.metric\n",argv[0]);
      return 0;
    }

  RSS( ref_import_by_extension( &ref_grid, argv[1] ), "import" );
  RSS( ref_export_by_extension( ref_grid, "bamg.msh" ), "export" );
  RSS( ref_part_metric( ref_grid_node(ref_grid), argv[2] ), "part metric" );
  RSS( ref_export_metric2d( ref_grid, "bamg.metric2d" ), "export" );

  RSS(ref_grid_free(ref_grid),"free");

  REIS( 0, system("bamg -M bamg.metric2d -b bamg.msh -o bamg-out.msh"),
	"bamg failed" );

  RSS( ref_import_by_extension( &ref_grid, "bamg-out.msh" ), "import" );
  RSS( ref_export_by_extension( ref_grid, "ref_bamg_test.b8.ugrid" ), "export");
  RSS(ref_export_tec_surf( ref_grid, "ref_bamg_test.tec" ),"ex" );

  RSS(ref_grid_free(ref_grid),"free");

  return 0;
}
