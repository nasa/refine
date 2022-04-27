#!/usr/bin/env bash

set -u
set -e
set -x

./bootstrap

clangflags='-g -O2  -Werror -Wall -Wextra -Wpedantic -Weverything -Wno-unused-macros -Wno-unreachable-code-return -Wno-padded -Wno-covered-switch-default -Wno-reserved-id-macro -Wno-documentation-unknown-command -Wno-poison-system-directories'
# -Wno-padded ref_mpi struct
# -Wno-covered-switch-default allow default on ENUM case
# -Wno-reserved-id-macro egads.h uses reserved macro names
# -Wdocumentation-unknown-command meshlink Types.h
# -Wno-poison-system-directories '/usr/local/include'
#    is unsafe for cross-compilation and I'm not cross-compiling

gcc9flags='-g -O2 -Werror -pedantic-errors -Wall -Wextra -Wunused -Wuninitialized -Wconversion'

zoltan_path="/Users/mpark/spack/opt/spack/darwin-mojave-x86_64/gcc-9.1.0/zoltan-3.83-5uh3ojfi7bp5ge7aavovf6lldduugwep"
egads_path="/Users/mpark/local/pkgs/EngSketchPad"
egads_svn_path="/Users/mpark/local/pkgs/EGADS/trunk"
egads_path="/Users/mpark/local/pkgs/ESPbeta-2022-02-02"
opencascade_path="/Users/mpark/local/pkgs/OpenCASCADE"
meshlink_path="/Users/mpark/local/pkgs/MeshLink"

# production spack packages
parmetis_path="/Users/mpark/local/pkgs/parmetis-4.0.3"
metis_path="/Users/mpark/local/pkgs/parmetis-4.0.3/metis"
mpi_path="/Users/mpark/homebrew"

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-mpi=${mpi_path} \
    --with-metis=${metis_path} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS="${clangflags}" \
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
    CFLAGS="-DHAVE_MPI ${gcc9flags}" \
    CC=mpicc \
    ) \
    || exit

mkdir -p clang
( cd clang && \
    ../configure \
    --prefix=`pwd` \
    --with-mpi=${mpi_path} \
    --with-metis=${metis_path} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS="${clangflags} -fsanitize=address -fsanitize=alignment -fsanitize=bounds -fsanitize=enum -fsanitize=vptr -fsanitize=integer-divide-by-zero -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fsanitize=nonnull-attribute -fsanitize=nullability-arg -fsanitize=nullability-assign -fsanitize=returns-nonnull-attribute -fsanitize=null -fsanitize=object-size -fsanitize=shift -fsanitize=signed-integer-overflow -fsanitize=unreachable -fsanitize=vla-bound" \
    CC=clang \
    ) \
    || exit


mkdir -p all
( cd all && \
      cmake .. \
	    -DCMAKE_INSTALL_PREFIX=`pwd` \
	    -DCMAKE_PREFIX_PATH="${mpi_path};${metis_path};${parmetis_path};${egads_path};${opencascade_path}/lib/cmake/opencascade"
) \
    || exit

