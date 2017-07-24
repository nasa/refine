#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "ref_mpi.h"

#include "ref_malloc.h"

int main( int argc, char *argv[] )
{

  RSS( ref_mpi_start( argc, argv ), "start" );

  if ( ref_mpi_n == 1 )
    { /* start */
      REIS( 1, ref_mpi_n, "n" );
      REIS( 0, ref_mpi_id, "n" );
      RAS( ref_mpi_master, "master" );
    }
  else
    {
      REF_INT part;
      REF_INT *a_size, *b_size;
      REF_INT bc = 5;

      if ( ref_mpi_master ) printf("number of processors %d \n",ref_mpi_n);

      RSS( ref_mpi_stopwatch_start(), "sw start");
      RSS( ref_mpi_bcast( &bc, 1, REF_INT_TYPE ), "bcast" );
      RSS( ref_mpi_stopwatch_stop( "integer broadcast" ), "sw start");

      ref_malloc_init( a_size, ref_mpi_n, REF_INT, REF_EMPTY );
      ref_malloc_init( b_size, ref_mpi_n, REF_INT, REF_EMPTY );

      for ( part = 0; part<ref_mpi_n ; part++ )
	a_size[part] = part;

      RSS( ref_mpi_stopwatch_start(), "sw start");
      RSS( ref_mpi_alltoall( a_size, b_size, REF_INT_TYPE ), "alltoall sizes");
      RSS( ref_mpi_stopwatch_stop( "integer alltoall" ), "sw start");

      for ( part = 0; part<ref_mpi_n ; part++ )
	REIS(part, a_size[part], "a_size changed" );
      for ( part = 0; part<ref_mpi_n ; part++ )
	REIS(ref_mpi_id, b_size[part], "b_size wrong" );

      ref_free( b_size );
      ref_free( a_size );
    }

  RSS( ref_mpi_stop( ), "stop" );

  return 0;
}
