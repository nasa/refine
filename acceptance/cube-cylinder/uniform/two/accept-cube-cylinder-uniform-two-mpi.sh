#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/src
    two=$1/two
else
    one=${HOME}/refine/zoltan/src
    two=${HOME}/refine/zoltan/two
fi

nproc=16

# ${two}/ref_geom_test ega.egads
# ${two}/ref_geom_test ega.egads ega.ugrid
# mv ref_geom_test.gas ega.gas
${two}/ref_acceptance ega.meshb ega.metric 0.1
mpiexec -np ${nproc} ${two}/ref_driver -i ega.meshb -g ega.egads -m ega.metric -o ref_driver1 -t
cp ref_gather_movie.tec ref_driver1_movie.tec
${two}/ref_acceptance ref_driver1.meshb ref_driver1.metric 0.1
${two}/ref_metric_test ref_driver1.meshb ref_driver1.metric > accept-cube-cylinder-uniform-two-mpi-01.status

mpiexec -np ${nproc} ${two}/ref_driver -i ref_driver1.meshb -g ega.egads -m ref_driver1.metric -o ref_driver2 -t
cp ref_gather_movie.tec ref_driver2_movie.tec
${two}/ref_acceptance ref_driver2.meshb ref_driver2.metric 0.1
${two}/ref_metric_test ref_driver2.meshb ref_driver2.metric > accept-cube-cylinder-uniform-two-mpi-02.status

cat accept-cube-cylinder-uniform-two-mpi-02.status
../../../check.rb accept-cube-cylinder-uniform-two-mpi-02.status 0.15 1.8



