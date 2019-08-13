#!/usr/bin/env bash

set -x # echo commands
set -e # exit on first error
set -u # Treat unset variables as error

if [ $# -gt 0 ] ; then
    one=$1/one
    two=$1/src
else
    one=${HOME}/refine/egads/one
    two=${HOME}/refine/egads/src
fi

geomfile=ega.egads

# ${two}/ref_geom_test ${geomfile} ega.meshb

${two}/ref_acceptance ega.meshb ega.metric 0.1
${two}/ref_driver -i ega.meshb -g ${geomfile} -m ega.metric -o ref_driver1 -t
mv ref_gather_movie.tec ref_driver1_movie.tec
${two}/ref_acceptance ref_driver1.meshb ref_driver1.metric 0.1
${two}/ref_metric_test ref_driver1.meshb ref_driver1.metric > accept-cube-cylinder-uniform-01.status

${two}/ref_driver -i ref_driver1.meshb -g ${geomfile} -m ref_driver1.metric -o ref_driver2 -t
mv ref_gather_movie.tec ref_driver2_movie.tec
${two}/ref_acceptance ref_driver2.meshb ref_driver2.metric 0.1
${two}/ref_metric_test ref_driver2.meshb ref_driver2.metric > accept-cube-cylinder-uniform-02.status

cat accept-cube-cylinder-uniform-02.status
../../check.rb accept-cube-cylinder-uniform-02.status 0.3 3.0



