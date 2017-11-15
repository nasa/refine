#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

mpiexec -np 2 ./ref_mpi_test
mpiexec -np 2 ./ref_part_test
mpiexec -np 2 ./ref_gather_test
mpiexec -np 8 ./ref_mpi_test
mpiexec -np 8 ./ref_part_test
mpiexec -np 8 ./ref_gather_test
