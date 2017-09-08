#!/usr/bin/env bash

set -x

./bootstrap

zoltan_path="/Users/mpark/spack/opt/spack/darwin-elcapitan-x86_64/gcc-4.9.2/zoltan-3.83-xvmrwrify56h4ceh4oacyns2t36wlx25"
egads_path="/Users/mpark/esp/EngSketchPad"

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran \
    ) \
    || exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS='-DHAVE_MPI -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized -Wno-long-long' \
    ) \
    || exit
