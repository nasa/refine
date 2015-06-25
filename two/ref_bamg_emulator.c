#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_export.h"
#include "ref_import.h"
#include "ref_part.h"
#include "ref_grid.h"
#include "ref_validation.h"
#include "ref_histogram.h"
#include "ref_args.h"

int main( int argc, char *argv[] )
{
  char *input_filename = NULL;
  char *output_filename = NULL;
  char *metric_filename = NULL;
  REF_INT location;

  REF_GRID ref_grid;
  
  RSS( ref_args_find( argc, argv, "-b", &location), "-b argument missing" );
  printf(" %s ",argv[location]);
  RAS( location<argc-1, "-b missing");
  input_filename = argv[1+location];
  printf("'%s'\n",input_filename);
  
  RSS( ref_args_find( argc, argv, "-M", &location), "-M argument missing" );
  printf(" %s ",argv[location]);
  RAS( location<argc-1, "-M missing");
  metric_filename = argv[1+location];
  printf("'%s'\n",metric_filename);
  
  RSS( ref_args_find( argc, argv, "-o", &location), "-o argument missing" );
  printf(" %s ",argv[location]);
  RAS( location<argc-1, "-o missing");
  output_filename = argv[1+location];
  printf("'%s'\n",output_filename);

  RSS( ref_import_by_extension( &ref_grid, input_filename ), "in" );


  
  RSS( ref_export_by_extension( ref_grid, output_filename ), "out" );

  return 0;
}
