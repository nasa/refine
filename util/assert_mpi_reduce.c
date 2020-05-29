/* mpicc -g -o assert_mpi_reduce assert_mpi_reduce.c */
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"
int main( int argc, char *argv[] )
{
  int n = 10;
  int nproc, rank;
  double *in, *out;

  if ( argc > 1 ) n = atoi(argv[1]);
    	
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  in = (double *)malloc(n*sizeof(double));
  out = (double *)malloc(n*sizeof(double));

  for (i=0; i<n; i++) in[i] = 0;
  if ( 0 == rank ){
    printf("calling MPI_Reduce( %d, MPI_DOUBLE, MPI_SUM) with %d\n",
	   n, nproc);
    fflush(stdout);
  }
  MPI_Allreduce( in, out, n, MPI_DOUBLE, MPI_SUM, 
		 MPI_COMM_WORLD );
  if ( 0 == rank ){
    printf("MPI_Allreduce returned\n");
    fflush(stdout);
  }

  free( out );
  free( in );
  
  MPI_Finalize( );
  return 0;
}
