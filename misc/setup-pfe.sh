#!/usr/bin/env bash

set -x

./bootstrap

module_path="/u/fun3d/shared/n1337/modules"
zoltan_path="${module_path}/Zoltan/3.82_mpt-2.17r13_ifort-2018.3.222"
parmetis_path="${module_path}/ParMETIS/4.0.3_mpt-2.17r13_ifort-2018.3.222"
egads_path="${module_path}/ESP/114/EngSketchPad"

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


