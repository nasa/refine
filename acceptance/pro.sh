#!/usr/bin/env bash

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
set +x # echo commands off for module
source /usr/local/pkgs/modules/init/bash
module purge
source acceptance/pro-modules.sh
module list
set -x # echo commands

root_dir=$(dirname $PWD)
source_dir=${root_dir}/refine

cov_dir=${root_dir}/_coverage
build32=${root_dir}/_refine-32bit
build64=${root_dir}/_refine-64bit

cd ${source_dir}
LOG=${root_dir}/log.bootstrap
trap "cat $LOG" EXIT
./bootstrap > $LOG 2>&1
trap - EXIT

date

mkdir -p ${cov_dir}
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
    CFLAGS="-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv -fp-stack-check -fstack-protector-all -fstack-security-check -prof-gen=srcpos -prof-dir=${cov_dir}" \
    CC=mpicc \
    > $LOG 2>&1
trap - EXIT

# remove spi to exclude conftest.c instead of using -comp
rm -f ${cov_dir}/*.spi

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
which mpiexec
which mpiexec >> $LOG 2>&1
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
which mpiexec
which mpiexec >> $LOG 2>&1
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

wait

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

wait

LOG=${root_dir}/log.build32-cov
trap "cat $LOG" EXIT
cd ${cov_dir}
echo "Create code coverage report" > $LOG 2>&1
profmerge                          >> $LOG 2>&1
codecov -prj unit_tests_cov -spi pgopti.spi \
        -dpi pgopti.dpi >> $LOG 2>&1
cov_report=CodeCoverageReport
mkdir -p ${cov_report}
mv CODE_COVERAGE.HTML ${cov_report}
mv CodeCoverage ${cov_report}
tar zcf ${cov_report}.tar.gz ${cov_report} && rm -rf ${cov_report}
echo "${cov_report}.tar.gz"  >> $LOG 2>&1
#get the summary in log file
codecov -prj unit_tests_cov -spi pgopti.spi \
        -dpi pgopti.dpi -txtlcov >> $LOG 2>&1
head -24 CODE_COVERAGE.TXT >> $LOG 2>&1
echo "done"  >> $LOG 2>&1
trap - EXIT

find ${source_dir} -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find ${source_dir} -name FAILED | wc -l`

