#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "ref_geom.h"

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_sort.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"
#include "ref_fixture.h"
#include "ref_export.h"
#include  "ref_dict.h"
#include  "ref_edge.h"

int main( int argc, char *argv[] )
{

  if ( 2 == argc )
    { /* fixture */
      char *filename = argv[1];
      if ( 0 == access( filename, R_OK ) )
	{
	  printf("EGADS project %s exisits, deleting\n",filename);
	  REIS(0, remove( filename ), "test clean up");
	}
      RSS( ref_geom_egads_fixture( filename ), "egads fixture" );
      printf("wrote EGADS project %s\n",filename);
    }

  if ( 3 == argc )
    { /* egads to grid */
      REF_GRID ref_grid;
      RSS( ref_geom_grid_from_egads( &ref_grid, argv[1] ), "from egads" );
      RSS( ref_export_by_extension( ref_grid, argv[2] ), "export" );
      RSS( ref_grid_free(ref_grid),"free");
    }

  REIS(REF_NULL,ref_geom_free(NULL),"dont free NULL");

  
  return 0;
}
