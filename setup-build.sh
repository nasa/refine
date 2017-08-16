#!/usr/bin/env bash

set -x

./bootstrap

# TESTS_ENVIRONMENT='valgrind --quiet --leak-check=full'

module_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS/4.0.3-1.10.2_intel_2013-2013.4.183_64"
zoltan_path="${module_path}/Zoltan/3.82-1.10.2_intel_2013-2013.4.183_64"
#egads_path="/ump/fldmd/home/mpark/esp/LINUX64"
egads_path="/ump/fldmd/home/mpark/local/pkgs/ESP110/EngSketchPad"

mkdir -p strict
( cd strict && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran \
    ) \
    || exit

mkdir -p egads
( cd egads && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran \
    ) \
    || exit

mkdir -p parmetis
( cd parmetis && \
    ../configure \
    --prefix=`pwd` \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    ) \
    || exit

mkdir -p zoltan
( cd zoltan && \
    ../configure \
    --prefix=`pwd` \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    CC=mpicc \
    FC=mpif90 \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    ) \
    || exit

mkdir -p profile
( cd profile && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -pg' \
    FC=gfortran \
    ) \
    || exit

