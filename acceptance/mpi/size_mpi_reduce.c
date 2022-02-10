/* mpicc -g -o size_mpi_reduce size_mpi_reduce.c */
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
int main(int argc, char *argv[]) {
  int repeat, repeats = 1000, step, steps = 10, i, n, target, increment,
              byte_target;
  int nproc, rank;
  double *in, *out;
  double start_time = 0, end_time = 0, delta_time = 0, total_time = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  target = 1000000;
  increment = 10000;
  for (step = 0; step < steps; step++) {
    total_time = 0;
    n = (target + increment * (step - steps / 2)) / sizeof(double);
    in = (double *)malloc(n * sizeof(double));
    out = (double *)malloc(n * sizeof(double));

    for (i = 0; i < n; i++) in[i] = 0.0;
    for (repeat = 0; i < repeats; repeat++) {
      if (0 == rank) start_time = MPI_Wtime();
      MPI_Reduce(in, out, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (0 == rank) {
        end_time = MPI_Wtime();
        delta_time = end_time - start_time;
        total_time += delta_time;
        if (delta_time > 0.5) {
          printf("MPI_Allreduce returned in %f sec. on %d of %d %lu bytes\n",
                 delta_time, repeat, repeats, n * sizeof(double));
          fflush(stdout);
        }
        if (delta_time > 1.0) {
          printf("MPI_Allreduce took longer than 1 sec. abort %lu bytes\n",
                 n * sizeof(double));
          fflush(stdout);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
    }
    if (0 == rank) {
      printf("total %f sec for %d repeats at %lu bytes\n", total_time, repeats,
             n * sizeof(double));
      fflush(stdout);
    }
    free(out);
    free(in);
  }
  MPI_Finalize();
  return 0;
}
