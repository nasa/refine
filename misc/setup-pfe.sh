#!/usr/bin/env bash

set -x

./bootstrap

module_path="/swbuild/fun3d/shared/fun3d_users/modules"
zoltan_path="${module_path}/Zoltan/3.82_mpt-2.17r13_ifort-2018.3.222"
parmetis_path="${module_path}/ParMETIS/4.0.3_mpt-2.17r13_ifort-2018.3.222"
egads_path="${module_path}/ESP/118/EngSketchPad"
opencascade_path="${module_path}/ESP/118/OpenCASCADE-7.3.1"

mpi_path="/nasa/hpe/mpt/2.17r13"

gcc_flags="-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized"

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS="${gcc_flags}" \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${occ_path} \
    --with-mpi=${mpi_path} \
    --with-parmetis=${parmetis_path} \
    CFLAGS="${gcc_flags}" \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS="-DHAVE_MPI ${gcc_flags}" \
    LIBS=-lmpi \
    ) \
    || exit

