#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include  "ref_dict.h"
#include  "ref_edge.h"

#include "ref_grid.h"
#include  "ref_node.h"
#include   "ref_mpi.h"
#include   "ref_matrix.h"
#include   "ref_sort.h"
#include   "ref_list.h"
#include  "ref_cell.h"
#include   "ref_adj.h"

#include "ref_fixture.h"

#include "ref_import.h"

int main( int argc, char *argv[] )
{

  if (2 == argc) 
    {
      REF_GRID ref_grid;
      char file[] = "ref_export_test.tec";
      RSS( ref_mpi_start( argc, argv ), "start" );
      ref_mpi_stopwatch_start();
      RSS( ref_import_by_extension( &ref_grid, argv[1] ), "examine header" );
      ref_mpi_stopwatch_stop("import");
      RSS(ref_export_tec_surf( ref_grid, file ),"export" );
      ref_mpi_stopwatch_stop("export");
      RSS(ref_grid_free(ref_grid),"free");
      RSS( ref_mpi_stop(  ), "stop" );
      return 0;
    }

  if (3 == argc) 
    {
      REF_GRID ref_grid;
      char grid[] = "ref_export_test.msh";
      char metric[] = "ref_export_test.metric2d";
      RSS(ref_fixture_twod_brick_grid( &ref_grid ), "set up brick" );
      RSS(ref_export_twod_msh( ref_grid, grid ),"export" );
      RSS(ref_export_metric2d( ref_grid, metric ),"export" );
      RSS(ref_grid_free(ref_grid),"free");
      return 0;
    }

  { /* export .vtk tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.vtk";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_vtk( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* export .tec tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_tec( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* export .tec hex */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.tec";
    RSS(ref_fixture_hex_grid( &ref_grid ), "set up hex" );
    RSS(ref_export_tec( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* export .fgrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.fgrid";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_fgrid( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* export .ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.ugrid";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_ugrid( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* export .b8.ugrid tet */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.b8.ugrid";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_b8_ugrid( ref_grid, file ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( file ), "test clean up");
  }

  { /* face id flag range */
    REF_GRID ref_grid;
    REF_INT min_faceid, max_faceid;
    RSS(ref_fixture_pri_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_faceid_range( ref_grid, &min_faceid, &max_faceid ),"range" );
    REIS( 10, min_faceid, "min");
    REIS(101, max_faceid, "max");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* export .cogsg tet */
    REF_GRID ref_grid;
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_cogsg( ref_grid, "ref_export_test.cogsg" ),"export" );
    RSS(ref_grid_free(ref_grid),"free");
    REIS(0, remove( "ref_export_test.cogsg" ), "test clean up");
    REIS(0, remove( "ref_export_test.bc" ), "test clean up");
  }

  if ( REF_FALSE ) /* removes gnuplot dependency from unit tests */
  { /* export .eps pri */
    REF_GRID ref_grid;
    RSS(ref_fixture_pri_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_eps( ref_grid, "ref_export_test.eps" ),"export" );
    REIS(0, remove( "ref_export_test.eps" ), "test clean up");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* export .html */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.html";
    RSS(ref_fixture_pri_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_html( ref_grid, file ),"export" );
    REIS(0, remove( file ), "test clean up");
    RSS(ref_grid_free(ref_grid),"free");
  }

  { /* export .meshb */
    REF_GRID ref_grid;
    char file[] = "ref_export_test.meshb";
    RSS(ref_fixture_tet_grid( &ref_grid ), "set up tet" );
    RSS(ref_export_meshb( ref_grid, file ),"export" );
    REIS(0, remove( file ), "test clean up");
    RSS(ref_grid_free(ref_grid),"free");
  }

  return 0;
}
