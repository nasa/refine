#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "ref_grid.h"
#include "ref_cell.h"
#include "ref_node.h"
#include "ref_list.h"
#include "ref_adj.h"
#include "ref_matrix.h"

#include "ref_sort.h"

#include "ref_migrate.h"

#include "ref_fixture.h"
#include "ref_import.h"
#include "ref_export.h"
#include "ref_dict.h"

#include "ref_mpi.h"
#include "ref_part.h"

#include "ref_gather.h"
#include "ref_adapt.h"

#include "ref_smooth.h"
#include  "ref_twod.h"
#include "ref_collapse.h"
#include "ref_split.h"
#include "ref_edge.h"

#include "ref_subdiv.h"
#include "ref_validation.h"
#include "ref_face.h"

#include "ref_malloc.h"
#include "ref_metric.h"
#include "ref_math.h"

#include "ref_histogram.h"

#include "ref_cavity.h"
#include "ref_clump.h"

static void echo_argv( int argc, char *argv[] )
{
  int pos;
  printf("\n");
  for ( pos = 0 ; pos < argc ; pos++ ) 
    printf(" %s",argv[pos]);
  printf("\n\n");
}

int main( int argc, char *argv[] )
{
  REF_GRID ref_grid = NULL;
  int opt;
  int passes, pass;
  
  echo_argv( argc, argv );

  while ((opt = getopt(argc, argv, "i:m:g:p:")) != -1)
    {
      switch (opt)
	{
	case 'i':
	  RSS( ref_import_by_extension( &ref_grid, optarg ), "import" );
	  break;
	case 'g':
	  RSS( ref_geom_egads_load( ref_grid_geom(ref_grid), optarg ), "ld e" );
	  break;
	case 'p':
	  RSS( ref_geom_load( ref_grid, optarg ), "load geom" );
	  break;
	case 'm':
	  RSS(ref_part_metric( ref_grid_node(ref_grid), optarg ), "part metric" );
	  break;
	case '?':
	default:
	  printf("parse error -%c\n",optopt);
	  printf("usage: \n %s\n",argv[0]);
	  printf("       [-i input_grid.ext]\n");
	  printf("./ref_geom_test ega.egads \n");
	  printf("./ref_geom_test ega.egads ega.ugrid\n");
	  printf("./ref_acceptance ega.ugrid ega.metric 0.1\n");
	  printf("./ref_driver -i ega.ugrid -g ega.egads -p ref_geom_test.gas -m ega.metric\n");
	  return 1;
	}
    }
  
  ref_mpi_stopwatch_stop("read grid");
  RSS(ref_validation_cell_volume(ref_grid),"vol");
  RSS( ref_histogram_quality( ref_grid ), "gram");
  RSS( ref_histogram_ratio( ref_grid ), "gram");
  RSS( ref_grid_inspect( ref_grid ), "inspect");

  passes = 5;
  for (pass = 0; pass<passes; pass++ )
    {
      printf(" pass %d of %d\n",pass,passes);
      RSS( ref_adapt_pass( ref_grid ), "pass");
      ref_mpi_stopwatch_stop("pass");
      RSS(ref_validation_cell_volume(ref_grid),"vol");
      RSS( ref_histogram_quality( ref_grid ), "gram");
      RSS( ref_histogram_ratio( ref_grid ), "gram");
      RSS( ref_grid_inspect( ref_grid ), "inspect");
      RSS(ref_migrate_to_balance(ref_grid),"balance");
      ref_mpi_stopwatch_stop("balance");
    }

  RSS( ref_gather_b8_ugrid( ref_grid, "ref_driver.b8.ugrid" ),
       "gather");
  ref_mpi_stopwatch_stop("gather");
  RSS(ref_export_tec_surf( ref_grid, "ref_driver_surf.tec" ),"surf tec" );
  ref_mpi_stopwatch_stop("surf tec");
  RSS(ref_geom_tec( ref_grid, "ref_driver_geom.tec" ),"geom tec" );
  ref_mpi_stopwatch_stop("geom tec");
  RSS(ref_geom_save( ref_grid, "ref_driver.gas" ),"geom tec" );
  ref_mpi_stopwatch_stop("geom assoc");
  RSS(ref_clump_stuck_edges( ref_grid, 0.5 ), "clump" );
  ref_mpi_stopwatch_stop("clump stuck");
 
  if ( NULL != ref_grid ) RSS(ref_grid_free( ref_grid ), "free");

  return 0;
}

