#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

module_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS/4.0.3-1.10.2_intel_2013-2013.4.183_64"
zoltan_path="${module_path}/Zoltan/3.82-1.10.2_intel_2013-2013.4.183_64"

egads_path=/ump/fldmd/home/wtjones1/local/pkgs-modules/ESP/svn

root_dir=$(dirname $PWD)
source_dir=${root_dir}/refine
build_dir=${root_dir}/_refine

cd ${source_dir}
./bootstrap

cd ${build_dir}
${source_dir}/configure \
    --prefix=${build_dir} \
    --with-EGADS=${egads_path} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran

env TMPDIR=${PWD} make -j 8
make check
make install





