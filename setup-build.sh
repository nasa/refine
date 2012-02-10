#!/usr/bin/env bash

set -x

mkdir -p strict
( cd strict && ../configure --prefix=`pwd` CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused' TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full' )

mkdir -p mpi
( cd mpi && ../configure --prefix=`pwd` CC=mpicc CFLAGS='-DHAVE_MPI -g -O2 -traceback' TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full' )

mkdir -p profile
( cd profile && ../configure --prefix=`pwd` CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -pg' )

