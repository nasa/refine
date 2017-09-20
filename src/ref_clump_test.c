#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_clump.h"

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
#include "ref_split.h"
#include  "ref_adapt.h"
#include   "ref_collapse.h"
#include    "ref_math.h"
#include   "ref_smooth.h"
#include    "ref_twod.h"
#include   "ref_metric.h"

#include   "ref_import.h"
#include   "ref_part.h"
#include   "ref_export.h"
#include   "ref_gather.h"
#include   "ref_swap.h"

int main( int argc, char *argv[] )
{

  if ( 1 < argc )
    {
      REF_GRID ref_grid;
      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "import" );
      if ( 2 < argc )
	{
	  RSS(ref_part_metric( ref_grid_node(ref_grid), argv[2] ), "part m");
	}
      else
	{
	  RSS(ref_metric_unit_node( ref_grid_node(ref_grid) ), "unit m");
	}

      RSS( ref_swap_triage( ref_grid ), "flip traige" );

      RSS( ref_clump_stuck_edges( ref_grid, 0.5 ), "stuck edge" );
      RSS( ref_export_tec_surf( ref_grid, "clump_surf.tec" ), "surf" );
      RSS( ref_gather_tec_movie_record_button( ref_grid_gather(ref_grid),
					       REF_TRUE ), "gather on" );      
      RSS( ref_gather_tec_movie_frame( ref_grid, "clump" ), "gather" );      
      RSS( ref_grid_free(ref_grid),"free");
    }

  return 0;
}
