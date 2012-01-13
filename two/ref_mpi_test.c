#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_mpi.h"

#include "ref_test.h"

int main( int argc, char *argv[] )
{

  { /* start */
    TSS( ref_mpi_start( argc, argv ), "start" );
    TEIS( 1, ref_mpi_n, "n" );
    TEIS( 0, ref_mpi_id, "n" );
  }

  return 0;
}
