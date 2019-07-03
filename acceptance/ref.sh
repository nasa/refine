#!/usr/bin/env bash

set -x # echo commands

set -e # exit on first error
set -u # Treat unset variables as error

# Setup bash module environment
. /usr/local/pkgs/modules/init/bash

module purge
source acceptance/ref-modules.sh

module_path="/u/shared/fun3d/fun3d_users/modules"
parmetis_path="${module_path}/ParMETIS/4.0.3-openmpi-2.1.1-intel_2017.2.174"
zoltan_path="${module_path}/Zoltan/3.82-openmpi-1.10.7-intel_2017.2.174"

egads_path="${module_path}/ESP/114/EngSketchPad"
opencascade_path="${module_path}/ESP/114/OpenCASCADE-6.8.1/lin64/gcc"

root_dir=$(dirname $PWD)
source_dir=${root_dir}/refine

parmetis_dir=${root_dir}/_refine-parmetis-egadslite
zoltan_dir=${root_dir}/_refine-zoltan-egadslite
egads_dir=${root_dir}/_refine-egads-full
strict_dir=${root_dir}/_refine-strict

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
    CC=gcc  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.strict-make-check-valgrind
trap "cat $LOG" EXIT
make -j 8 > $LOG 2>&1
make check TESTS_ENVIRONMENT='valgrind --quiet  --error-exitcode=1 --leak-check=full' >> $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.strict-make-distcheck
trap "cat $LOG" EXIT
( make distcheck > $LOG 2>&1 && cp refine-*.tar.gz ${root_dir} || touch FAILED ) &
trap - EXIT

date

mkdir -p ${egads_dir}
cd ${egads_dir}

LOG=${root_dir}/log.egads-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${egads_dir} \
    --with-EGADS=${egads_path} \
    --with-OpenCASCADE=${opencascade_path} \
    CFLAGS='-g -O2 -pedantic-errors -Wall -Wextra -Werror -Wunused -Wuninitialized' \
    CC=gcc  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.egads-make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.egads-make-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

mkdir -p ${zoltan_dir}
cd ${zoltan_dir}

LOG=${root_dir}/log.zoltan-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${zoltan_dir} \
    --with-zoltan=${zoltan_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    CC=mpicc > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.zoltan-make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.zoltan-make-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

mkdir -p ${parmetis_dir}
cd ${parmetis_dir}

LOG=${root_dir}/log.parmetis-configure
trap "cat $LOG" EXIT
${source_dir}/configure \
    --prefix=${parmetis_dir} \
    --with-parmetis=${parmetis_path} \
    --with-EGADS=${egads_path} \
    --enable-lite \
    CFLAGS='-DHAVE_MPI -g -O2 -traceback -Wall -ftrapuv' \
    CC=mpicc > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.parmetis-make
trap "cat $LOG" EXIT
env TMPDIR=${PWD} make -j 8  > $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.parmetis-make-install
trap "cat $LOG" EXIT
make install > $LOG 2>&1
trap - EXIT

date

LOG=${root_dir}/log.zoltan-unit
trap "cat $LOG" EXIT
cd ${zoltan_dir}/src
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
mpiexec -np 2 ./ref_subdiv_test >> $LOG 2>&1
mpiexec -np 2 ./ref_subdiv_test >> $LOG 2>&1
mpiexec -np 4 ./ref_agents_test >> $LOG 2>&1
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
mpiexec -np 8 ./ref_subdiv_test >> $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.parmetis-unit
trap "cat $LOG" EXIT
cd ${parmetis_dir}/src
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
mpiexec -np 2 ./ref_subdiv_test >> $LOG 2>&1
mpiexec -np 4 ./ref_subdiv_test >> $LOG 2>&1
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
mpiexec -np 8 ./ref_subdiv_test >> $LOG 2>&1
trap - EXIT

LOG=${root_dir}/log.accept-2d-linear-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/linear/two
( ./accept-2d-two.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-2d-polar-2-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/polar-2/two
( ./accept-2d-two.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-2d-masabl-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/masabl/two
( ./accept-2d-two.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-2d-mixed
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/2d/mixed
( ./accept-2d-mixed.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-3d-linear-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear/two
( ./accept-3d-two.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

sleep 10 # allow some tests to complete before making more

LOG=${root_dir}/log.accept-cube-cylinder-gen-aflr3
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/gen/aflr3
( ./accept-cube-cylinder-aflr3.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-uniform-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/uniform/two
( ./accept-cube-cylinder-uniform-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-linear010-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/linear010/two
( ./accept-cube-cylinder-linear010-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-polar-2-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2/two
( ./accept-cube-cylinder-polar-2-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-3d-polar-1-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/polar-1/two
( ./accept-3d-two.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

sleep 10 # allow some tests to complete before making more

LOG=${root_dir}/log.accept-cube-sphere-uniform-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-sphere/uniform/two
( ./accept-cube-sphere-uniform-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cube-sphere-ring-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-sphere/ring/two
( ./accept-cube-sphere-ring-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-annulus-uniform
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/annulus/uniform
( ./accept-annulus-uniform.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cone-cone-uniform-two
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cone-cone/uniform/two
( ./accept-cone-cone-uniform-two.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-cone-cone-recon
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cone-cone/recon
( ./accept-cone-cone-recon.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

sleep 10 # allow some tests to complete before making more

LOG=${root_dir}/log.accept-om6-recon
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/om6/recon
( ./accept-om6-recon.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-revolve-pencil-curve
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/revolve-pencil/curve
( ./revolve-pencil-curve.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-hemisphere-uniform
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/hemisphere/uniform
( ./accept-hemisphere-uniform.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-sphere-cube-tetgen
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/sphere-cube/tetgen
( ./accept-sphere-cube-init.sh ${egads_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-3d-u5
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/u5
( ./accept-3d-u5.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

sleep 10 # allow some tests to complete before making more

LOG=${root_dir}/log.accept-inflate-normal
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/normal
( ./inflate.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-inflate-radial
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/radial
( ./inflate.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

LOG=${root_dir}/log.accept-inflate-mapbc
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/mapbc
( ./inflate.sh ${strict_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait

# 2 procs
LOG=${root_dir}/log.accept-3d-linear-two-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/3d/linear/two
( ./accept-3d-two-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 1-1-4-4 procs
LOG=${root_dir}/log.accept-cube-cylinder-polar-2-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/polar-2/two
( ./accept-cube-cylinder-polar-2-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

# 4 procs
LOG=${root_dir}/log.accept-hemisphere-uniform-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/hemisphere/uniform
( ./accept-hemisphere-uniform-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait

# 8 procs
LOG=${root_dir}/log.accept-inflate-normal-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/normal
( ./inflate-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait

# 8 procs
LOG=${root_dir}/log.accept-inflate-radial-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/radial
( ./inflate-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait

# 8 procs
LOG=${root_dir}/log.accept-inflate-mapbc-para
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/inflate/mapbc
( ./inflate-para.sh ${parmetis_dir} > $LOG 2>&1 || touch FAILED ) &
trap - EXIT

wait

LOG=${root_dir}/log.accept-cube-cylinder-uniform-valgrind
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/uniform/two
( ./accept-cube-cylinder-uniform-two-valgrind.sh ${egads_dir} > $LOG 2>&1 || touch FAILED )
trap - EXIT

LOG=${root_dir}/log.accept-cube-cylinder-uniform-valgrind-mpi
trap "cat $LOG" EXIT
cd ${source_dir}/acceptance/cube-cylinder/uniform/two
( ./accept-cube-cylinder-uniform-two-valgrind-mpi.sh ${zoltan_dir} > $LOG 2>&1 || touch FAILED )
trap - EXIT

grep RAC ${root_dir}/log.accept-* > ${root_dir}/log.summary

find ${source_dir} -name FAILED

echo -e \\n\
# Build has failed if any failed cases have been reported
exit `find ${source_dir} -name FAILED | wc -l`

