#!/usr/bin/env bash

set -x

./bootstrap

zoltan_path="/Users/mpark/spack/opt/spack/darwin-highsierra-x86_64/gcc-8.2.0/zoltan-3.83-6gpgaz4xpcw73nqwhup5wkq45m5xehb2"
egads_path="/Users/mpark/local/pkgs/EngSketchPad"
opencascade_path="/Users/mpark/local/pkgs/OpenCASCADE-6.8.1"
gccflags='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized'
parmetis_path="/Users/mpark/spack/opt/spack/darwin-highsierra-x86_64/gcc-8.2.0/parmetis-4.0.3-xasrm372i7dgoygsou7f2m7pqppyubmd"
metis_path="/Users/mpark/spack/opt/spack/darwin-highsierra-x86_64/gcc-8.2.0/metis-5.1.0-nugla3v2bdfmvp5ajz67ypecdymm6wmv"

# set DYLD_LIBRARY_PATH for libzoltan.dylib

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS="${gccflags}" \
    CC=gcc-8 \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS="${gccflags}" \
    CC=gcc-8 \
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

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS="-DHAVE_MPI ${gccflags} -Wno-long-long" \
    ) \
    || exit

