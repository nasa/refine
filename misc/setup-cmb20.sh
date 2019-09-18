#!/usr/bin/env bash

set -x

./bootstrap

module_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules"
mpi_path="/usr/local/pkgs-modules/openmpi_1.10.2_intel_2017"
zoltan_path="${module_path}/Zoltan/3.82-1.10.2_intel_2017-2017.2.174"
egads_path="${module_path}/ESP/114/EngSketchPad"
opencascade_path="${module_path}/ESP/114/OpenCASCADE-6.8.1/Linux"
parmetis_path="${module_path}/ParMETIS/4.0.3-1.10.2_intel_2017-2017.2.174"
parmetis_path="/ump/fldmd/home/mpark/local/pkgs/parmetis-4.0.3/build/Linux-x86_64"

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-mpi=${mpi_path} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CC=icc \
    CFLAGS='-g -O3 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,1419 -fp-stack-check -fstack-security-check' \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-parmetis=${parmetis_path} \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS='-DHAVE_MPI -g -O3 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,1419 -fp-stack-check -fstack-security-check' \
    ) \
    || exit

exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CC=mpicc \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,49,1419 -fp-stack-check -fstack-protector-all -fstack-security-check' \
    ) \
    || exit

mkdir -p sanitize
( cd sanitize && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-fsanitize=address -ggdb -g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    LDFLAGS='-fsanitize=address' \
    ) \
    || exit

