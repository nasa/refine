#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

set +x # echo commands
. /usr/local/pkgs/modules/init/bash

date

module load gcc_4.9.1_64
module load openmpi_1.10.2_intel_2017
module load intel.2017.2.174
module load git # for git describe

module use --append /ump/fldmd/home/wtjones1/Modules/modulefiles
module load OpenCASCADE/6.6.0
module load ESP/svn

module list
set -x # echo commands

module_path="/ump/fldmd/home/casb-shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS/4.0.3-1.10.2_intel_2017-2017.2.174"
zoltan_path="${module_path}/Zoltan/3.82-1.10.2_intel_2017-2017.2.174"

egads_path=/ump/fldmd/home/wtjones1/local/pkgs-modules/ESP/svn

root_dir=$(dirname $PWD)
source_dir=${root_dir}/refine
build_dir=${root_dir}/_refine
strict_dir=${root_dir}/_strict

cd ${source_dir}
LOG=${root_dir}/log.bootstrap
trap "cat $LOG" EXIT
./bootstrap > $LOG 2>&1
trap - EXIT

date

mkdir -p ${strict_dir}
cd ${strict_dir}

LOG=${root_dir}/log.strict-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${strict_dir} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    FC=gfortran  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.make-valgrind
trap "cat $LOG" EXIT
make -j 8 > $LOG 2>&1
make check TESTS_ENVIRONMENT='valgrind --quiet  --error-exitcode=1 --leak-check=full' >> $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.make-distcheck
trap "cat $LOG" EXIT
make distcheck > $LOG 2>&1
cp refine-*.tar.gz ${root_dir}
trap - EXIT

date

mkdir -p ${build_dir}
cd ${build_dir}

LOG=${root_dir}/log.configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${build_dir} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    CC=mpicc \
    FC=mpif90  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.make-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.unit-para
trap "cat $LOG" EXIT
cd ${build_dir}/src
echo para-unit > $LOG 2>&1
mpiexec -np 2 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 2 ./ref_part_test >> $LOG 2>&1
mpiexec -np 2 ./ref_gather_test >> $LOG 2>&1
mpiexec -np 8 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 8 ./ref_part_test >> $LOG 2>&1
mpiexec -np 8 ./ref_gather_test >> $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-2d-linear-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/linear/two
./accept-2d-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-2d-polar-2-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/polar-2/two
./accept-2d-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-3d-linear-one
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear/one
./accept-3d-one.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-3d-linear-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear/two
./accept-3d-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-3d-linear-two-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear/two
./accept-3d-two-para.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-cube-cylinder-uniform-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/uniform/two
./accept-cube-cylinder-uniform-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-cube-cylinder-linear010-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/linear010/two
./accept-cube-cylinder-linear010-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-cube-cylinder-polar-2-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2/two
./accept-cube-cylinder-polar-2-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-3d-polar-1-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/polar-1/two
./accept-3d-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-cube-sphere-uniform-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-sphere/uniform/two
./accept-cube-sphere-uniform-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-cube-sphere-ring-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-sphere/ring/two
./accept-cube-sphere-ring-two.sh ${build_dir} > $LOG 2>&1
trap - EXIT

date
