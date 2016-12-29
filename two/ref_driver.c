#include <stdlib.h>
#include <stdio.h>
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

static int print_usage( char *name )
{
  if ( ref_mpi_master )
    {
      printf("usage: \n %s\n",name);
      printf("       [-ig input_grid.ext]\n");
      printf("       [-og output_grid.ext]\n");
    }
  RSS( ref_mpi_stop(  ), "stop" );
  return 1;
}

int main( int argc, char *argv[] )
{
  REF_INT pos;
  REF_GRID ref_grid = NULL;
  REF_BOOL parse_error;

  RSS( ref_mpi_start( argc, argv ), "start" );
  
  if ( ref_mpi_master )
    {
      printf("\n");
      for ( pos = 0 ; pos < argc ; pos++ ) 
	printf(" %s",argv[pos]);
      printf("\n\n");
    }

  parse_error = REF_TRUE;
  pos = 1;
  while( pos < argc ) {
    if( strcmp(argv[pos],"-ic") == 0 ) {
      if ( ref_mpi_master )
	printf("%d: -ic\n",pos);
      if ( pos+2 > argc )
	return(print_usage(argv[0]));
      pos++;
      if ( ref_mpi_master )
	printf("%d: %s\n",pos,argv[pos]);
      if ( NULL != ref_grid ) RSS(ref_grid_free( ref_grid ), "free");
      RSS( ref_import_by_extension( &ref_grid, argv[pos] ), "import" );
      parse_error = REF_FALSE;
    }
    if ( parse_error ) return(print_usage(argv[0]));
   pos++; 
  }
  
  if ( NULL != ref_grid ) RSS(ref_grid_free( ref_grid ), "free");

  RSS( ref_mpi_stop(  ), "stop" );

  return 0;
}
