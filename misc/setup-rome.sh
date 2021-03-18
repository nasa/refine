#!/usr/bin/env bash

set -x

./bootstrap

module_path="/swbuild/fun3d/shared/fun3d_users/modules"
egads_path="${module_path}/ESP/119-beta.2021.03.05.1211/EngSketchPad"
occ_path="${module_path}/ESP/119-beta.2021.03.05.1211/OpenCASCADE-7.3.1"
parmetis_path="/swbuild/fun3d/shared/fun3d_users/sles15/test_modules/ParMETIS-64/4.0.3_mpt-2.23_ifort-2018.3.222"

mpi_path="/nasa/hpe/mpt/2.23_sles15_patch11654"

gcc_flags="-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized"

mkdir -p rome
(cd rome && \
    ../configure \
    --prefix=`pwd` \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${occ_path} \
    --with-mpi=${mpi_path} \
    --with-parmetis=${parmetis_path} \
    --with-metis=${parmetis_path} \
    CFLAGS="${gcc_flags}" \
    ) \
    || exit
