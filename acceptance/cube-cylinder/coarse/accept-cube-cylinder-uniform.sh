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

geomfile=cube-cylinder.egads

# ${two}/ref_geom_test ${geomfile} cube-cylinder.meshb.meshb

${two}/ref_acceptance cube-cylinder.meshb cube-cylinder.metric 0.1
${two}/ref_driver -i cube-cylinder.meshb -g ${geomfile} -m cube-cylinder.metric -o ref_driver1
${two}/ref_acceptance ref_driver1.meshb ref_driver1.metric 0.1
${two}/ref_metric_test ref_driver1.meshb ref_driver1.metric > accept-cube-cylinder-uniform-01.status

${two}/ref_driver -i ref_driver1.meshb -g ${geomfile} -m ref_driver1.metric -o ref_driver2
${two}/ref_acceptance ref_driver2.meshb ref_driver2.metric 0.1
${two}/ref_metric_test ref_driver2.meshb ref_driver2.metric > accept-cube-cylinder-uniform-02.status

cat accept-cube-cylinder-uniform-02.status
../../check.rb accept-cube-cylinder-uniform-02.status 0.16 1.6



