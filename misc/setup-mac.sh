#!/usr/bin/env bash

set -x

./bootstrap

gccflags='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized -Wconversion'

zoltan_path="/Users/mpark/spack/opt/spack/darwin-mojave-x86_64/gcc-9.1.0/zoltan-3.83-5uh3ojfi7bp5ge7aavovf6lldduugwep"
egads_path="/Users/mpark/local/pkgs/EngSketchPad"
opencascade_path="/Users/mpark/local/pkgs/OpenCASCADE"

# production spack packages
parmetis_path="/Users/mpark/spack/opt/spack/darwin-mojave-x86_64/gcc-9.1.0/parmetis-4.0.3-jwaxhhbilbtvsmt2tskek4k72nel7gtc"
metis_path="/Users/mpark/spack/opt/spack/darwin-mojave-x86_64/gcc-9.1.0/metis-5.1.0-czqd5zteq5zfccffypwbjrlb5joqeoyw"

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS="${gccflags}" \
    CC=gcc-9 \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS="${gccflags}" \
    CC=gcc-9 \
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

mkdir -p both
( cd both && \
    ../configure \
    --prefix=`pwd` \
    --with-metis=${metis_path} \
    --with-parmetis=${parmetis_path} \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS="-DHAVE_MPI ${gccflags} -Wno-long-long" \
    ) \
    || exit

