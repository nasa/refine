#!/usr/bin/env bash

set -x

./bootstrap

zoltan_path="/Users/mpark/spack/opt/spack/darwin-elcapitan-x86_64/gcc-7.2.0/zoltan-3.83-xxyq2a3qxwtcrdt3otmiflhgqaz7gial"
egads_path="/Users/mpark/esp/EngSketchPad"
gccflags='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized'

# set DYLD_LIBRARY_PATH for libzoltan.dylib

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS="${gccflags}" \
    CC=gcc-7 \
    FC=gfortran-7 \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    CFLAGS="${gccflags}" \
    CC=gcc-7 \
    FC=gfortran-7 \
    ) \
    || exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS="-DHAVE_MPI ${gccflags} -Wno-long-long" \
    ) \
    || exit
