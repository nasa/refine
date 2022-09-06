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

gcc12flags='-g -O2 -Werror -pedantic-errors -Wall -Wextra -Wunused -Wuninitialized -Wconversion'

egads_path="${HOME}/local/pkgs/EngSketchPad"
egads_svn_path="${HOME}/local/pkgs/EGADS"
opencascade_path="${HOME}/local/pkgs/OpenCASCADE"
meshlink_path="${HOME}/local/pkgs/MeshLink"

parmetis_path="${HOME}/local/pkgs/parmetis-4.0.3-gcc-12-mpich"
metis_path="${HOME}/local/pkgs/parmetis-4.0.3-gcc-12-mpich"
mpi_path="${HOME}/local/pkgs/mpich-4.0.2/gcc-12-install/bin/mpiexec"

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

