/* mpicc -g -o assert_mpi_min_loc assert_mpi_min_loc.c */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
int main(int argc, char *argv[]) {
  int i, n = 1403;
  int nproc, rank;
  typedef struct MY_MPI_DBLWHO MY_MPI_DBLWHO;
  struct MY_MPI_DBLWHO {
    double val;
    int rank;
  };
  MY_MPI_DBLWHO *in, *out;

  if (argc > 1) {
    n = atoi(argv[1]);
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  in = (MY_MPI_DBLWHO *)malloc(n * sizeof(MY_MPI_DBLWHO));
  out = (MY_MPI_DBLWHO *)malloc(n * sizeof(MY_MPI_DBLWHO));

  for (i = 0; i < n; i++) {
    in[i].val = (double)rank;
    in[i].rank = rank;
  }
  if (0 == rank)
    printf("calling MPI_Allreduce( %d, MPI_DOUBLE_INT, MPI_MINLOC) with %d\n",
           n, nproc);
  MPI_Allreduce(in, out, n, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
  if (0 == rank) printf("MPI_Allreduce returned\n");

  free(out);
  free(in);

  MPI_Finalize();
  return 0;
}
