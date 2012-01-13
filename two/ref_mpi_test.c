#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "ref_mpi.h"

#include "ref_test.h"

int main( int argc, char *argv[] )
{

  { /* start */
    TSS( ref_mpi_start( argc, argv ), "start" );
  }

  return 0;
}
