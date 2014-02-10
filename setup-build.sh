#!/usr/bin/env bash

set -x

# TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full'

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-parmetis=/usr/local/pkgs-modules/ParMETIS_4.0.3-openmpi_1.6.5_intel_2013-intel.2013.4.183_64 \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    ) \
    || exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${HOME}/local/pkgs/zoltan_distrib_v3.6/Zoltan_v3.6/build \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    ) \
    || exit

mkdir -p profile
( cd profile && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -pg' \
    FC=gfortran \
    ) \
    || exit

