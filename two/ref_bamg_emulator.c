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
  
  RSS( ref_args_find( argc, argv, "-b", &location), "-b argument missing" );
  printf(" %s : ",argv[location]);
  RAS( location<argc-1, "-b missing");
  input_filename = argv[1+location];
  printf(" %s : ",input_filename);
  
  RSS( ref_args_find( argc, argv, "-M", &location), "-M argument missing" );
  printf(" %s : ",argv[location]);
  RAS( location<argc-1, "-M missing");
  metric_filename = argv[1+location];
  printf(" %s : ",metric_filename);
  
  RSS( ref_args_find( argc, argv, "-o", &location), "-o argument missing" );
  printf(" %s : ",argv[location]);
  RAS( location<argc-1, "-o missing");
  output_filename = argv[1+location];
  printf(" %s : ",output_filename);
  
  return 0;
}
