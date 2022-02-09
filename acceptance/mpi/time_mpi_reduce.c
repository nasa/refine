/* mpicc -g -o time_mpi_reduce time_mpi_reduce.c */
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
int main(int argc, char *argv[]) {
  int i, nnode = 1000000, ldim = 6, chunk, n;
  int nproc, rank;
  double *in, *out;
  double start_time = 0, end_time = 0, delta_time = 0, total_time = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (0 == rank) {
    printf("usage: \n %s <nnode=1M> <ldim=6>\n", argv[0]);
  }
  if (argc > 1) nnode = atoi(argv[1]);
  if (argc > 2) ldim = atoi(argv[2]);

  chunk = nnode / nproc + 1;

  n = (ldim + 1) * chunk;
  in = (double *)malloc(n * sizeof(double));
  out = (double *)malloc(n * sizeof(double));

  if (0 == rank) {
    printf("chunk size %d doubles %lu bytes\n", n, n * sizeof(double));
    fflush(stdout);
  }

  for (i = 0; i < n; i++) in[i] = 0.0;
  for (i = 0; i < nproc; i++) {
    if (0 == rank) {
      printf("calling MPI_Reduce( %d, MPI_DOUBLE, MPI_SUM) on %d\n", n, nproc);
      fflush(stdout);
      start_time = MPI_Wtime();
    }
    MPI_Reduce(in, out, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (0 == rank) {
      end_time = MPI_Wtime();
      delta_time = end_time - start_time;
      total_time += delta_time;
      printf("MPI_Allreduce returned in %f sec\n", delta_time);
      fflush(stdout);
    }
  }
  if (0 == rank) {
    printf("total %f sec\n", total_time);
    fflush(stdout);
  }
  free(out);
  free(in);

  MPI_Finalize();
  return 0;
}
