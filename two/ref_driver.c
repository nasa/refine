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
  
  echo_argv( argc, argv );

  while ((opt = getopt(argc, argv, "i:m:g:p:")) != -1)
    {
      switch (opt)
	{
	case 'i':
	  RSS( ref_import_by_extension( &ref_grid, optarg ), "import" );
	  break;
	case 'g':
	  RSS( REF_IMPLEMENT, "load geom" );
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
	  printf("./ref_driver -i ega.ugrid -g ega.egads -p ref_geom_test.gas -m ega.metric\n");
	  return 1;
	}
    }
  
  if ( NULL != ref_grid ) RSS(ref_grid_free( ref_grid ), "free");

  return 0;
}

