/* mpicc -g -o size_mpi_reduce size_mpi_reduce.c */
#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"
int main(int argc, char *argv[]) {
  int repeat, repeats = 1000, step, steps = 10, i, n, target, increment;
  int nproc, rank;
  double *in, *out;
  double start_time = 0, end_time = 0, delta_time = 0, total_time = 0;
  double time_limit = 10.0;
  int slow_down = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc > 1) slow_down = 0;

  if (0 == rank) {
    if (slow_down) {
      printf("slow version, only MPI_Reduce\n");
    } else {
      printf("fast version, MPI_Reduce and MPI_Bcast interleave\n");
    }
  }

  target = 1000000;
  increment = 100000;
  for (step = 0; step < steps; step++) {
    total_time = 0;
    n = (target + increment * (step - steps / 2)) / (int)sizeof(double);
    in = (double *)malloc((size_t)n * sizeof(double));
    out = (double *)malloc((size_t)n * sizeof(double));
    for (i = 0; i < n; i++) in[i] = 0.0;

    for (repeat = 0; repeat < repeats; repeat++) {
      if (0 == rank) {
        start_time = MPI_Wtime();
      }
      MPI_Reduce(in, out, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (!slow_down) {
        i = 0;
        MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
      if (0 == rank) {
        end_time = MPI_Wtime();
        delta_time = end_time - start_time;
        total_time += delta_time;
        if (total_time > time_limit) {
          printf("total %f sec for %d repeats at %lu bytes [GAVE UP]\n",
                 total_time, repeat + 1, (size_t)n * sizeof(double));
          fflush(stdout);
        }
      }
    }
    if (0 == rank) {
      printf("total %f sec for %d repeats at %lu bytes\n", total_time, repeats,
             (size_t)n * sizeof(double));
      fflush(stdout);
    }
    free(out);
    free(in);
  }
  MPI_Finalize();
  return 0;
}
