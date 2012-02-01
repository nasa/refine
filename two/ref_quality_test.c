#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_test.h"

#include "ref_grid.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_quality.h"

#include "ref_adj.h"
#include "ref_node.h"
#include "ref_metric.h"
#include "ref_cell.h"

#include "ref_face.h"
#include "ref_sort.h"
#include "ref_shard.h"
#include "ref_subdiv.h"
#include "ref_edge.h"
#include "ref_math.h"
#include "ref_dict.h"

#include "ref_swap.h"

#include "ref_test.h"
#include "ref_fixture.h"

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;

  if (argc>1) 
    {

      printf("reading %s\n",argv[1]);
      RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from ugrid");
      printf("complete.\n");
      
      RSS(ref_grid_inspect( ref_grid ), "inspection");

      printf("tec surf.\n");
      RSS(ref_export_tec_surf( ref_grid, "qual-surf-orig.tec" ), "tec surf");

      printf("check for multiple_face_cell.\n");
      RSS( ref_quality_report_multiple_face_cell( ref_grid, "qual-orig.tec"), 
	   "report" );

      printf("try swapping.\n");
      RSS( ref_quality_swap_multiple_face_cell( ref_grid), "swap" );

      printf("tec surf.\n");
      RSS(ref_export_tec_surf( ref_grid, "qual-surf-swap.tec" ), "tec surf");

      printf("check for multiple_face_cell.\n");
      RSS( ref_quality_report_multiple_face_cell( ref_grid, "qual-swap.tec"), 
	   "report" );

      printf("try split.\n");
      RSS( ref_quality_split_multiple_face_cell( ref_grid ), "split" );

      printf("check for multiple_face_cell.\n");
      RSS( ref_quality_report_multiple_face_cell( ref_grid, "qual-split.tec"), 
	   "report" );

      printf("export\n");
      RSS(ref_export_b8_ugrid( ref_grid, "quality.b8.ugrid" ),"from ugrid");

      printf("done.\n");
      return 0;
    }

  { /* find mark */
    REF_GRID ref_grid;
    REF_DBL vol;

    TSS( ref_fixture_tet_grid( &ref_grid ), "tet fixture" );
    TSS( ref_quality_tet_vol( ref_grid, 0, &vol), "get vol");

    TWDS(1.0/6.0,vol,-1.0,"expected vol");

    TSS( ref_grid_free( ref_grid ), "free" );
  }


  return 0;
}
