#!/usr/bin/env bash

set -x

./bootstrap

module_path="/u/shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS-64/4.0.3-mpt-2.23-intel_2018.3.222"
gcc_parmetis_path="${module_path}/ParMETIS/4.0.3-mpt-2.23-gcc_6.2.0"
egads_path="${module_path}/ESP/120-beta.2021.09.20.1202/EngSketchPad"
egads_svn_path="/u/mpark/local/pkgs/EGADS/trunk"
occ_path="${module_path}/ESP/120-beta.2021.09.20.1202/OpenCASCADE"
meshlink_path="/u/mpark/local/pkgs/MeshLink"

mpi_path="/opt/hpe/hpc/mpt/mpt-2.23"

export openmpi_path="/usr/local/pkgs-modules/openmpi_2.1.1_intel_2017"

export openmpi_parmetis_path="${module_path}/ParMETIS/4.0.3-openmpi-2.1.1-intel_2017.2.174"


gcc_flags="-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized"
icc_flags="-g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,49,1419 -ftrapuv"

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-MeshLink=${meshlink_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${occ_path} \
    --with-mpi=${mpi_path} \
    --with-parmetis=${parmetis_path} \
    CC=icc \
    CFLAGS="${icc_flags}" \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-metis=${gcc_parmetis_path} \
    --with-parmetis=${gcc_parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS="-DHAVE_MPI ${gcc_flags}" \
    ) \
    || exit

mkdir -p openmpi
( cd openmpi && \
    ../configure \
    --prefix=`pwd` \
    --with-metis=${openmpi_parmetis_path} \
    --with-parmetis=${openmpi_parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=${openmpi_path}/bin/mpicc \
    CFLAGS="-DHAVE_MPI ${icc_flags}" \
    ) \
    || exit

