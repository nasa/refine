
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_MPI
#include "mpi.h"
#endif

#include "ref_mpi.h"

REF_STATUS ref_mpi_start( int argc, char *argv[] )
{
  int n, id;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);

  MPI_Comm_size(MPI_COMM_WORLD,&n);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
#else
  SUPRESS_UNUSED_COMPILER_WARNING(argc);
  SUPRESS_UNUSED_COMPILER_WARNING(argv[0]);
  n = 1;
  id = 0;
#endif

  printf("proc %d of %d\n",id,n);

  return REF_SUCCESS;
}
