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

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid;
  REF_INT count;

  if (argc<2) 
    {
      printf("usage: %s filename.extension\n",argv[0]);
      return 0;
    }

  printf("reading %s\n",argv[1]);
  RSS(ref_import_by_extension( &ref_grid, argv[1] ),"from ugrid");
  printf("complete.\n");
      
  RSS(ref_grid_inspect( ref_grid ), "inspection");

  printf("check for multiple_face_cell.\n");
  RSS( ref_quality_report_multiple_face_cell(ref_grid,&count,"qual-orig.tec"), 
       "report" );
  if ( 0 == count )
    {
      printf("zero mulitple face tets, export\n");
      RSS(ref_export_b8_ugrid( ref_grid, "quality.b8.ugrid" ),"from ugrid");
      printf("done.\n");
      return 0;
    }

  printf("try swapping.\n");
  RSS( ref_quality_swap_multiple_face_cell( ref_grid), "swap" );

  printf("check for multiple_face_cell.\n");
  RSS( ref_quality_report_multiple_face_cell(ref_grid,&count,"qual-swap.tec"), 
       "report" );
  if ( 0 == count )
    {
      printf("zero mulitple face tets, export\n");
      RSS(ref_export_b8_ugrid( ref_grid, "quality.b8.ugrid" ),"from ugrid");
      printf("done.\n");
      return 0;
    }

  printf("try split.\n");
  RSS( ref_quality_split_multiple_face_cell( ref_grid ), "split" );

  printf("check for multiple_face_cell.\n");
  RSS( ref_quality_report_multiple_face_cell(ref_grid,&count,"qual-split.tec"), 
       "report" );
  if ( 0 == count )
    {
      printf("zero mulitple face tets, export\n");
      RSS(ref_export_b8_ugrid( ref_grid, "quality.b8.ugrid" ),"from ugrid");
      printf("done.\n");
      return 0;
    }

  printf("try split.\n");
  RSS( ref_quality_split_multiple_face_cell( ref_grid ), "split" );

  printf("check for multiple_face_cell.\n");
  RSS( ref_quality_report_multiple_face_cell(ref_grid,&count,"qual-split2.tec"), 
       "report" );

  if ( 0 == count )
    {
      printf("zero mulitple face tets, export\n");
      RSS(ref_export_b8_ugrid( ref_grid, "quality.b8.ugrid" ),"from ugrid");
      printf("done.\n");
      return 0;
    }

  return 1;
}
