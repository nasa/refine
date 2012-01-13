
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "ref_mpi.h"

REF_INT ref_mpi_n = 1;
REF_INT ref_mpi_id = 0;

REF_STATUS ref_mpi_start( int argc, char *argv[] )
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&ref_mpi_n);
  MPI_Comm_rank(MPI_COMM_WORLD,&ref_mpi_id);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(argc);
  SUPRESS_UNUSED_COMPILER_WARNING(argv[0]);
  ref_mpi_n = 1;
  ref_mpi_id = 0;
#endif

  return REF_SUCCESS;
}

REF_STATUS ref_mpi_stop( )
{

#ifdef HAVE_MPI
  MPI_Finalize( );
#endif

  return REF_SUCCESS;
}
