#!/usr/bin/env bash

set -x

./bootstrap

zoltan_path="/Users/mpark/spack/opt/spack/darwin-elcapitan-x86_64/gcc-7.2.0/zoltan-3.83-xxyq2a3qxwtcrdt3otmiflhgqaz7gial"
egads_path="/Users/mpark/local/pkgs/EngSketchPad"
gccflags='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized'
parmetis_path="/Users/mpark/spack/opt/spack/darwin-elcapitan-x86_64/gcc-7.2.0/parmetis-4.0.3-lxz27qnvcp7blhtqltxvngthbnabs24h"
metis_path="/Users/mpark/spack/opt/spack/darwin-elcapitan-x86_64/gcc-7.2.0/metis-5.1.0-x5ptqjut6c5wzhzs74pcje3lpwfl5lct"

# set DYLD_LIBRARY_PATH for libzoltan.dylib

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS="${gccflags}" \
    CC=gcc-8 \
    FC=gfortran-8 \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    CFLAGS="${gccflags}" \
    CC=gcc-8 \
    FC=gfortran-8 \
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

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-metis=${metis_path} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS="-DHAVE_MPI ${gccflags} -Wno-long-long" \
    ) \
    || exit
