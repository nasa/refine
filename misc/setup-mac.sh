#!/usr/bin/env bash

set -u
set -e
set -x

./bootstrap

clangflags='-g -O2  -Werror -Wall -Wextra -Wpedantic -Weverything -Wno-unused-macros -Wno-unreachable-code-return -Wno-padded -Wno-covered-switch-default -Wno-reserved-id-macro -Wno-documentation-unknown-command'
# -Wno-padded ref_mpi struct
# -Wno-covered-switch-default allow default on ENUM case
# -Wno-reserved-id-macro egads.h uses reserved macro names
# -Wdocumentation-unknown-command meshlink Types.h

gcc9flags='-g -O2 -Werror -pedantic-errors -Wall -Wextra -Wunused -Wuninitialized -Wconversion'

zoltan_path="/Users/mpark/spack/opt/spack/darwin-mojave-x86_64/gcc-9.1.0/zoltan-3.83-5uh3ojfi7bp5ge7aavovf6lldduugwep"
egads_path="/Users/mpark/local/pkgs/EngSketchPad"
opencascade_path="/Users/mpark/local/pkgs/OpenCASCADE"
meshlink_path="/Users/mpark/local/pkgs/MeshLink"

# production spack packages
parmetis_path="/Users/mpark/local/pkgs/parmetis-4.0.3"
metis_path="/Users/mpark/local/pkgs/parmetis-4.0.3"
mpi_path="/Users/mpark/local/pkgs/openmpi-4.0.5/build"

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-mpi=${mpi_path} \
    --with-metis=${metis_path} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    --with-MeshLink=${meshlink_path} \
    CFLAGS="${clangflags}" \
    ) \
    || exit

mkdir -p all
( cd all && \
      cmake .. \
	    -DCMAKE_INSTALL_PREFIX=`pwd` \
	    -DCMAKE_PREFIX_PATH="${mpi_path};${metis_path};${parmetis_path};${egads_path};${opencascade_path}/lib/cmake/opencascade"
) \
    || exit
