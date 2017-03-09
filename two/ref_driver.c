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
  REF_GRID background_grid = NULL;
  int opt;
  int passes = 15, pass;
  REF_BOOL output_clumps = REF_FALSE;
  char output_project[1024];
  char output_filename[1024];

  snprintf( output_project, 1024, "ref_driver" );
	  
  echo_argv( argc, argv );

  while ((opt = getopt(argc, argv, "i:m:g:p:o:s:c")) != -1)
    {
      switch (opt)
	{
	case 'i':
	  RSS( ref_import_by_extension( &ref_grid, optarg ), "import" );
	  RSS( ref_import_by_extension( &background_grid, optarg ), "import" );
	  break;
	case 'g':
	  RSS( ref_geom_egads_load( ref_grid_geom(ref_grid), optarg ), "ld e" );
	  break;
	case 'p':
	  RSS( ref_geom_load( ref_grid, optarg ), "load geom" );
	  break;
	case 'm':
	  RSS(ref_part_metric( ref_grid_node(ref_grid), optarg ), "part m");
	  if ( NULL != background_grid )
	    RSS(ref_part_metric( ref_grid_node(background_grid), optarg ),
		"part m back" );
	  break;
	case 'o':
	  snprintf( output_project, 1024, "%s", optarg );
	  break;
	case 's':
	  passes = atoi(optarg);
	  break;
	case 'c':
	  output_clumps = REF_TRUE;
	  break;
	case '?':
	default:
	  printf("parse error -%c\n",optopt);
	  printf("usage: \n %s\n",argv[0]);
	  printf("       [-i input_grid.ext]\n");
	  printf("       [-g geometry.egads]\n");
	  printf("       [-p parameterization-restart.gas]\n");
	  printf("       [-m input_project.metric]\n");
	  printf("       [-s adapt_cycles] default is 15\n");
	  printf("       [-o output_project]\n");
	  printf("       [-c] output clumps\n");
	  printf("./ref_geom_test ega.egads \n");
	  printf("./ref_geom_test ega.egads ega.ugrid\n");
	  printf("./ref_acceptance ega.ugrid ega.metric 0.1\n");
	  printf("./ref_driver -i ega.ugrid -g ega.egads -p ref_geom_test.gas -m ega.metric\n");
	  printf("cp ref_driver.b8.ugrid ref_driver1.b8.ugrid\n");
	  printf("cp ref_driver.gas ref_driver1.gas\n");
	  printf("./ref_acceptance ref_driver1.b8.ugrid ref_driver1.metric 0.1\n");
	  printf("./ref_driver -i ref_driver1.b8.ugrid -g ega.egads -p ref_driver1.gas -m ref_driver1.metric\n");
	  return 1;
	}
    }
  
  ref_mpi_stopwatch_stop("read grid");
  RSS(ref_validation_cell_volume(ref_grid),"vol");
  RSS( ref_histogram_quality( ref_grid ), "gram");
  RSS( ref_histogram_ratio( ref_grid ), "gram");

  for (pass = 0; pass<passes; pass++ )
    {
      printf(" pass %d of %d\n",pass,passes);
      RSS( ref_adapt_pass( ref_grid ), "pass");
      ref_mpi_stopwatch_stop("pass");
      if ( NULL != background_grid )
	{
	  RSS( ref_metric_interpolate( ref_grid, background_grid ), "interp" );
	  ref_mpi_stopwatch_stop("interp");
	}
      RSS(ref_validation_cell_volume(ref_grid),"vol");
      RSS( ref_histogram_quality( ref_grid ), "gram");
      RSS( ref_histogram_ratio( ref_grid ), "gram");
      RSS(ref_migrate_to_balance(ref_grid),"balance");
      ref_mpi_stopwatch_stop("balance");
    }

  RSS( ref_geom_verify_param( ref_grid ), "loaded params" );
  ref_mpi_stopwatch_stop("verify params");
  snprintf( output_filename, 1024, "%s.b8.ugrid", output_project );
  RSS( ref_gather_b8_ugrid( ref_grid, output_filename ),
       "gather");
  ref_mpi_stopwatch_stop("gather");
  snprintf( output_filename, 1024, "%s_surf.tec", output_project );
  RSS(ref_export_tec_surf( ref_grid, output_filename ),"surf tec" );
  snprintf( output_filename, 1024, "%s_geom.tec", output_project );
  RSS(ref_geom_tec( ref_grid, output_filename ),"geom tec" );
  snprintf( output_filename, 1024, "%s.gas", output_project );
  RSS(ref_geom_save( ref_grid, output_filename ),"geom tec" );
  ref_mpi_stopwatch_stop("tec");
  if ( output_clumps)
    {
      RSS(ref_clump_stuck_edges( ref_grid, 0.5 ), "clump" );
      ref_mpi_stopwatch_stop("clump stuck");
    }
  
  if ( NULL != background_grid ) RSS(ref_grid_free( background_grid ), "free");
  if ( NULL != ref_grid ) RSS(ref_grid_free( ref_grid ), "free");

  return 0;
}

