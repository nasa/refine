#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ref_layer.h"
#include "ref_mpi.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  {  /* make layers */
    REF_LAYER ref_layer;

    RSS(ref_layer_create(&ref_layer),"create");

    REIS( 0, ref_layer_n(ref_layer), "check total layers");

    RSS(ref_layer_free(ref_layer),"layer");
  }

  RSS( ref_mpi_stop(  ), "stop" );
  
  return 0;
}
