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

wait
sleep 2

# 2 procs
LOG=${root_dir}/log.accept-2d-linear-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/linear
( ./accept-2d-linear-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-2d-masabl-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/masabl
( ./accept-2d-masabl-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-2d-polar-2-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/polar-2
( ./accept-2d-polar-2-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 2 procs
LOG=${root_dir}/log.accept-2d-linear-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/linear
( ./accept-2d-linear-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

 # 4 procs
LOG=${root_dir}/log.accept-2d-masabl-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/masabl
( ./accept-2d-masabl-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-2d-polar-2-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/polar-2
( ./accept-2d-polar-2-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 2 procs
LOG=${root_dir}/log.accept-facebody-side-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/facebody/side
( ./accept-facebody-side-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 2 procs
LOG=${root_dir}/log.accept-3d-linear-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear
( ./accept-3d-linear-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 1-1-4-4 procs
LOG=${root_dir}/log.accept-cube-cylinder-polar-2-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2
( ./accept-cube-cylinder-polar-2-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-hemisphere-uniform-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/hemisphere/uniform
( ./accept-hemisphere-uniform-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 2 procs
LOG=${root_dir}/log.accept-facebody-side-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/facebody/side
( ./accept-facebody-side-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 2 procs
LOG=${root_dir}/log.accept-3d-linear-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear
( ./accept-3d-linear-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 1-1-4-4 procs
LOG=${root_dir}/log.accept-cube-cylinder-polar-2-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2
( ./accept-cube-cylinder-polar-2-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-hemisphere-uniform-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/hemisphere/uniform
( ./accept-hemisphere-uniform-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-normal-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/normal
( ./inflate-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-3d-subset-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/subset
( ./accept-3d-subset-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-normal-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/normal
( ./inflate-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-3d-subset-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/subset
( ./accept-3d-subset-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-radial-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/radial
( ./inflate-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-cube-subset-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube/subset
( ./accept-cube-subset-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-radial-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/radial
( ./inflate-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-cube-subset-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube/subset
( ./accept-cube-subset-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-mapbc-para-32
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/mapbc
( ./inflate-para.sh ${build32} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

# 8 procs
LOG=${root_dir}/log.accept-inflate-mapbc-para-64
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/mapbc
( ./inflate-para.sh ${build64} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait
sleep 2

grep RAC ${root_dir}/log.accept-* > ${root_dir}/log.summary

find ${source_dir} -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find ${source_dir} -name FAILED | wc -l`

