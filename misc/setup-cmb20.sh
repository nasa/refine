#!/usr/bin/env bash

set -x

./bootstrap

module_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS/4.0.3-1.10.2_intel_2017-2017.2.174"
zoltan_path="${module_path}/Zoltan/3.82-1.10.2_intel_2017-2017.2.174"
egads_path="${module_path}/ESP/114/EngSketchPad"
opencascade_path="${module_path}/ESP/114/OpenCASCADE-6.8.1/Linux"

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
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
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
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,1419' \
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
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -w3 -wd1418,2259,2547,981,11074,11076,1572,49,1419' \
    ) \
    || exit

mkdir -p profile
( cd profile && \
    ../configure \
    --prefix=`pwd` \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -pg' \
    ) \
    || exit

