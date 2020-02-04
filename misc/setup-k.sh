#!/usr/bin/env bash

set -x

./bootstrap

module_path="/u/shared/fun3d/fun3d_users/modules"
zoltan_path="${module_path}/Zoltan/3.82-mpt-2.19-intel_2018.3.222"
parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.19-intel_2018.3.222"
egads_path="${module_path}/ESP/116/EngSketchPad"
occ_path="${module_path}/ESP/116/OpenCASCADE-7.3.1"

mpi_path="/opt/hpe/hpc/mpt/mpt-2.19"

gcc_flags="-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized"
icc_flags="-g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,49,1419 -ftrapuv"

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
    CC=icc \
    CFLAGS="${icc_flags}" \
    ) \
    || exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=icc \
    CFLAGS="-DHAVE_MPI ${icc_flags}" \
    LIBS=-lmpi \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=icc \
    CFLAGS="-DHAVE_MPI ${icc_flags}" \
    LIBS=-lmpi \
    ) \
    || exit


