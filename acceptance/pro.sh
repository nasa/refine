#!/usr/bin/env bash

set -x # echo commands

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

module purge
source acceptance/pro-modules.sh

root_dir=$(dirname $PWD)
source_dir=${root_dir}/refine

build32=${root_dir}/_refine-32bit
build64=${root_dir}/_refine-64bit

cd ${source_dir}
LOG=${root_dir}/log.bootstrap
trap "cat $LOG" EXIT
./bootstrap > $LOG 2>&1
trap - EXIT

date

mkdir -p ${build32}
cd ${build32}

LOG=${root_dir}/log.build32-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${build32} \
    --with-parmetis=${parmetis32_path} \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv -fp-stack-check -fstack-protector-all -fstack-security-check' \
    CC=mpicc \
    > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.build32-make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.build32-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.build32-unit
trap "cat $LOG" EXIT
cd ${build32}/src
echo para-unit > $LOG 2>&1
mpiexec -np 2 ./ref_agents_test >> $LOG 2>&1
mpiexec -np 2 ./ref_edge_test >> $LOG 2>&1
mpiexec -np 2 ./ref_gather_test >> $LOG 2>&1
mpiexec -np 2 ./ref_interp_test >> $LOG 2>&1
mpiexec -np 2 ./ref_metric_test >> $LOG 2>&1
mpiexec -np 2 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 2 ./ref_node_test >> $LOG 2>&1
mpiexec -np 2 ./ref_part_test >> $LOG 2>&1
mpiexec -np 2 ./ref_migrate_test >> $LOG 2>&1
mpiexec -np 2 ./ref_cavity_test >> $LOG 2>&1
mpiexec -np 2 ./ref_elast_test >> $LOG 2>&1
mpiexec -np 2 ./ref_recon_test >> $LOG 2>&1
mpiexec -np 2 ./ref_search_test >> $LOG 2>&1
mpiexec -np 2 ./ref_subdiv_test >> $LOG 2>&1
mpiexec -np 8 ./ref_agents_test >> $LOG 2>&1
mpiexec -np 8 ./ref_edge_test >> $LOG 2>&1
mpiexec -np 8 ./ref_gather_test >> $LOG 2>&1
mpiexec -np 8 ./ref_interp_test >> $LOG 2>&1
mpiexec -np 8 ./ref_metric_test >> $LOG 2>&1
mpiexec -np 8 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 8 ./ref_node_test >> $LOG 2>&1
mpiexec -np 8 ./ref_part_test >> $LOG 2>&1
mpiexec -np 8 ./ref_migrate_test >> $LOG 2>&1
mpiexec -np 8 ./ref_cavity_test >> $LOG 2>&1
mpiexec -np 8 ./ref_elast_test >> $LOG 2>&1
mpiexec -np 8 ./ref_recon_test >> $LOG 2>&1
mpiexec -np 8 ./ref_search_test >> $LOG 2>&1
mpiexec -np 8 ./ref_subdiv_test >> $LOG 2>&1
trap - EXIT

date

mkdir -p ${build64}
cd ${build64}

LOG=${root_dir}/log.build64-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${build64} \
    --with-parmetis=${parmetis64_path} \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv -fp-stack-check -fstack-protector-all -fstack-security-check' \
    CC=mpicc \
    > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.build64-make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.build64-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.build64-unit
trap "cat $LOG" EXIT
cd ${build64}/src
echo para-unit > $LOG 2>&1
mpiexec -np 2 ./ref_agents_test >> $LOG 2>&1
mpiexec -np 2 ./ref_edge_test >> $LOG 2>&1
mpiexec -np 2 ./ref_gather_test >> $LOG 2>&1
mpiexec -np 2 ./ref_interp_test >> $LOG 2>&1
mpiexec -np 2 ./ref_metric_test >> $LOG 2>&1
mpiexec -np 2 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 2 ./ref_node_test >> $LOG 2>&1
mpiexec -np 2 ./ref_part_test >> $LOG 2>&1
mpiexec -np 2 ./ref_migrate_test >> $LOG 2>&1
mpiexec -np 2 ./ref_cavity_test >> $LOG 2>&1
mpiexec -np 2 ./ref_elast_test >> $LOG 2>&1
mpiexec -np 2 ./ref_recon_test >> $LOG 2>&1
mpiexec -np 2 ./ref_search_test >> $LOG 2>&1
mpiexec -np 2 ./ref_subdiv_test >> $LOG 2>&1
mpiexec -np 8 ./ref_agents_test >> $LOG 2>&1
mpiexec -np 8 ./ref_edge_test >> $LOG 2>&1
mpiexec -np 8 ./ref_gather_test >> $LOG 2>&1
mpiexec -np 8 ./ref_interp_test >> $LOG 2>&1
mpiexec -np 8 ./ref_metric_test >> $LOG 2>&1
mpiexec -np 8 ./ref_mpi_test >> $LOG 2>&1
mpiexec -np 8 ./ref_node_test >> $LOG 2>&1
mpiexec -np 8 ./ref_part_test >> $LOG 2>&1
mpiexec -np 8 ./ref_migrate_test >> $LOG 2>&1
mpiexec -np 8 ./ref_cavity_test >> $LOG 2>&1
mpiexec -np 8 ./ref_elast_test >> $LOG 2>&1
mpiexec -np 8 ./ref_recon_test >> $LOG 2>&1
mpiexec -np 8 ./ref_search_test >> $LOG 2>&1
mpiexec -np 8 ./ref_subdiv_test >> $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.accept-3d-linear-mpt-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear
( ./accept-3d-linear-mpt.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-polar-2-mpt-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2
( ./accept-cube-cylinder-polar-2-mpt.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-3d-linear-mpt-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear
( ./accept-3d-linear-mpt.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-polar-2-mpt-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2
( ./accept-cube-cylinder-polar-2-mpt.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

grep RAC ${root_dir}/log.accept-* > ${root_dir}/log.summary

find ${source_dir} -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find ${source_dir} -name FAILED | wc -l`

