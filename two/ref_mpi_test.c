#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_mpi.h"

#include "ref_test.h"

int main( int argc, char *argv[] )
{

  TSS( ref_mpi_start( argc, argv ), "start" );

  if ( ref_mpi_n == 1 )
    { /* start */
      TEIS( 1, ref_mpi_n, "n" );
      TEIS( 0, ref_mpi_id, "n" );
      TAS( ref_mpi_master, "master" );
    }
  else
    {
      REF_INT bc = 5;
      if ( ref_mpi_master ) printf("number of processors %d \n",ref_mpi_n);
      RSS( ref_mpi_stopwatch_start(), "sw start");
      RSS( ref_mpi_bcast( &bc, 1, REF_INT_TYPE ), "bcast" );
      RSS( ref_mpi_stopwatch_stop( "integer broadcast" ), "sw start");
    }

  TSS( ref_mpi_stop( ), "stop" );

  return 0;
}
